.BCH.XZ.multinomial <- function(design, weight, ref.class,
                                method.regression = "Analytic",
                                method.SE = "Analytic",
                                lower = -10, upper = 10,
                                tol = 1e-4, maxiter = 5000,
                                cluster = NULL) {
  design <- as.matrix(design)
  weight <- as.matrix(weight)
  L <- ncol(weight)
  p <- ncol(design)
  if(nrow(design) != nrow(weight)){
    stop("design and BCH weights must have the same number of rows")
  }
  evaluate <- function(params, gradient = FALSE){
    value <- bch_multinomial_cpp(
      params, design, weight, ref.class, gradient, FALSE
    )
    if(gradient) value else value$objective
  }
  optimization <- stats::optim(
    rep(0, p * (L - 1L)),
    fn = function(x) evaluate(x),
    gr = if(method.regression == "Analytic"){
      function(x) evaluate(x, gradient = TRUE)$gradient
    }else NULL,
    method = "L-BFGS-B",
    lower = rep(lower, p * (L - 1L)),
    upper = rep(upper, p * (L - 1L)),
    control = list(maxit = maxiter, factr = max(tol, 1e-12) / .Machine$double.eps)
  )
  fit <- bch_multinomial_cpp(
    optimization$par, design, weight, ref.class, TRUE, TRUE
  )
  information <- if(method.SE == "Numeric"){
    numDeriv::hessian(function(x) evaluate(x), optimization$par)
  }else{
    fit$information
  }
  score <- .ZY.cluster.score(fit$score, cluster)
  inverse.information <- MASS::ginv(information)
  vcov <- inverse.information %*% crossprod(score) %*% inverse.information
  list(
    params = optimization$par,
    vcov = vcov,
    information = information,
    score = score,
    objective = fit$objective,
    probability = fit$probability,
    converged = optimization$convergence == 0L,
    iterations = unname(optimization$counts["function"]),
    bound.diagnostics = .three.step.bound.diagnostics(
      optimization$par,
      rep(lower, length(optimization$par)),
      rep(upper, length(optimization$par))
    )
  )
}

.BCH.XZ.block.diagonal <- function(matrices) {
  dimensions <- vapply(matrices, nrow, integer(1))
  value <- matrix(0, sum(dimensions), sum(dimensions))
  start <- 0L
  for(j in seq_along(matrices)){
    index <- start + seq_len(dimensions[j])
    value[index, index] <- matrices[[j]]
    start <- start + dimensions[j]
  }
  value
}

.BCH.XZ.statistics <- function(params, se, covariates, L, ref.class,
                               covariates.time.cross = FALSE) {
  estimate <- LTA.vector.to.parameters(
    params, covariates, L, ref.class, covariates.time.cross
  )
  se.object <- .three.step.reference.na(
    LTA.vector.to.parameters(
      se, covariates, L, ref.class, covariates.time.cross
    ), ref.class
  )
  z <- params / se
  tail1 <- stats::pnorm(-abs(z))
  tail2 <- 2 * tail1
  list(
    estimate = estimate,
    se = se.object,
    z = .three.step.reference.na(
      LTA.vector.to.parameters(
        z, covariates, L, ref.class, covariates.time.cross
      ), ref.class
    ),
    tail1 = .three.step.reference.na(
      LTA.vector.to.parameters(
        tail1, covariates, L, ref.class, covariates.time.cross
      ), ref.class
    ),
    tail2 = .three.step.reference.na(
      LTA.vector.to.parameters(
        tail2, covariates, L, ref.class, covariates.time.cross
      ), ref.class
    )
  )
}
