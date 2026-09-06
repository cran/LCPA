flexmix.LCA <- function(response, L, poly.value, nrep = 20,
                        starts = 100, maxiter.warmup = 20, vis = TRUE,
                        control.flexmix){
  if(!requireNamespace("flexmix", quietly = TRUE)){
    stop("method = 'flexmix' requires the optional package 'flexmix'.", call. = FALSE)
  }

  response <- as.matrix(response)
  N <- nrow(response)
  I <- ncol(response)
  poly.value <- as.integer(poly.value)
  invalid <- vapply(seq_len(I), function(i){
    any(!is.finite(response[, i])) ||
      any(!(response[, i] %in% 0:(poly.value[i] - 1L)))
  }, logical(1))
  if(any(invalid)){
    stop("flexmix SEM requires complete LCA responses coded from 0 to C_i - 1.")
  }
  if(N < L){
    stop("flexmix requires N >= L so every latent class can be non-empty.")
  }

  if(all(poly.value == 2L)){
    data.flexmix <- data.frame(.response = I(response))
    formula.flexmix <- .response ~ 1
    model.flexmix <- flexmix::FLXMCmvbinary()
    binary.driver <- TRUE
  }else{
    data.flexmix <- data.frame(.response = I(response))
    formula.flexmix <- .response ~ 1
    model.flexmix <- .flexmix.LCA.model(response, poly.value)
    binary.driver <- FALSE
  }

  fit.result <- .fit.flexmix.SEM(
    formula.flexmix, data.flexmix, model.flexmix, L,
    nrep, starts, maxiter.warmup, vis, control.flexmix, "LCA"
  )
  fit <- fit.result$fit
  components <- methods::slot(fit, "components")

  poly.max <- max(poly.value)
  par <- array(NA_real_, dim = c(L, I, poly.max))
  for(l in seq_len(L)){
    if(binary.driver){
      prob <- as.numeric(methods::slot(components[[l]][[1L]], "parameters")$center)
      par[l, , 1L] <- 1 - prob
      par[l, , 2L] <- prob
    }else{
      probability <- methods::slot(
        components[[l]][[1L]], "parameters"
      )$probability
      par[l, , ] <- probability
    }
  }
  if(any(!is.finite(par[!is.na(par)])) || any(par[!is.na(par)] < 0)){
    stop("flexmix SEM returned invalid multinomial parameters.")
  }

  P.Z <- as.numeric(methods::slot(fit, "prior"))
  P.Z <- P.Z / sum(P.Z)
  P.Z.Xn <- as.matrix(flexmix::posterior(fit))
  P.Z.Xn <- P.Z.Xn / rowSums(P.Z.Xn)
  Z <- max.col(P.Z.Xn, ties.method = "first")
  Log.Lik <- get.Log.Lik.LCA(response, par, P.Z)
  npar <- get.npar.LCA(poly.value, L)
  AIC <- -2 * Log.Lik + 2 * npar
  BIC <- -2 * Log.Lik + npar * log(N)

  if(vis){
    .print.estimation.summary("flexmix SEM", Log.Lik, BIC)
  }

  list(
    params = list(par = par, P.Z = P.Z),
    model = fit,
    npar = npar,
    Log.Lik = Log.Lik,
    AIC = AIC,
    BIC = BIC,
    best_BIC = BIC,
    P.Z.Xn = P.Z.Xn,
    P.Z = P.Z,
    Z = Z,
    probability = NULL,
    Log.Lik.nrep = fit.result$Log.Lik.nrep
  )
}
