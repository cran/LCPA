flexmix.LPA <- function(response, L, constraint = "V0", nrep = 20,
                        starts = 100, maxiter.warmup = 20, vis = TRUE,
                        control.flexmix){
  if(!requireNamespace("flexmix", quietly = TRUE)){
    stop("method = 'flexmix' requires the optional package 'flexmix'.", call. = FALSE)
  }

  response <- as.matrix(response)
  N <- nrow(response)
  I <- ncol(response)
  if(anyNA(response) || any(!is.finite(response))){
    stop("flexmix SEM requires complete finite LPA responses.")
  }
  constraint <- .validate.LPA.constraint(constraint, I)
  if(N < L){
    stop("flexmix requires N >= L so every latent profile can be non-empty.")
  }

  data.flexmix <- data.frame(.response = I(response))
  model.flexmix <- function(){
    .flexmix.LPA.model(response, L, constraint)
  }
  fit.result <- .fit.flexmix.SEM(
    .response ~ 1, data.flexmix, model.flexmix, L,
    nrep, starts, maxiter.warmup, vis, control.flexmix, "LPA"
  )
  fit <- fit.result$fit
  components <- methods::slot(fit, "components")

  means <- matrix(NA_real_, nrow = L, ncol = I)
  covs <- array(NA_real_, dim = c(I, I, L))
  for(l in seq_len(L)){
    parameters.l <- methods::slot(components[[l]][[1L]], "parameters")
    means[l, ] <- as.numeric(parameters.l$center)
    covs[, , l] <- as.matrix(parameters.l$cov)
  }
  if(any(!is.finite(means)) || any(!is.finite(covs))){
    stop("flexmix SEM returned invalid Gaussian parameters.")
  }
  covariance.global <- matrix(stats::cov(response), I, I)
  if(any(!is.finite(covariance.global))) covariance.global <- diag(I)
  covariance.global <- .repair.flexmix.covariance(covariance.global, diag(I))
  for(l in seq_len(L)){
    covs[, , l] <- .repair.flexmix.covariance(covs[, , l], covariance.global)
  }

  P.Z <- as.numeric(methods::slot(fit, "prior"))
  P.Z <- P.Z / sum(P.Z)
  expectation <- lpa_expectation_cpp(
    response, means, covs, P.Z, repair = FALSE
  )
  if(!isTRUE(expectation$valid)){
    stop("flexmix SEM returned an invalid stabilized covariance matrix.")
  }
  P.Z.Xn <- expectation$posterior
  Z <- max.col(P.Z.Xn, ties.method = "first")
  Log.Lik <- expectation$Log.Lik
  npar <- get.npar.LPA(I, L, constraint)
  AIC <- -2 * Log.Lik + 2 * npar
  BIC <- -2 * Log.Lik + npar * log(N)

  if(vis){
    .print.estimation.summary("flexmix SEM", Log.Lik, BIC)
  }

  list(
    params = list(means = means, covs = covs, P.Z = P.Z),
    model = fit,
    npar = npar,
    Log.Lik = Log.Lik,
    AIC = AIC,
    BIC = BIC,
    best_BIC = BIC,
    P.Z.Xn = P.Z.Xn,
    P.Z = P.Z,
    Z = Z,
    Log.Lik.nrep = fit.result$Log.Lik.nrep
  )
}
