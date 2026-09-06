.extract.RMixtComp.LPA <- function(fit.cur, response, L, variable.names){
  N <- nrow(response)
  I <- ncol(response)
  P.Z <- as.numeric(RMixtCompUtilities::getProportion(fit.cur))
  if(length(P.Z) != L || any(!is.finite(P.Z)) || any(P.Z <= 0)){
    stop("RMixtComp SEM returned invalid latent profile proportions.")
  }
  P.Z <- P.Z / sum(P.Z)

  P.Z.Xn <- as.matrix(RMixtCompUtilities::getTik(fit.cur, log = FALSE))
  if(!identical(dim(P.Z.Xn), c(N, as.integer(L))) ||
     any(!is.finite(P.Z.Xn)) || any(rowSums(P.Z.Xn) <= 0)){
    stop("RMixtComp SEM returned invalid posterior profile probabilities.")
  }
  P.Z.Xn <- P.Z.Xn / rowSums(P.Z.Xn)

  means <- matrix(NA_real_, nrow = L, ncol = I)
  covs <- array(0, dim = c(I, I, L))
  for(i in seq_len(I)){
    param.i <- as.matrix(RMixtCompUtilities::getParam(fit.cur, variable.names[i]))
    if(nrow(param.i) != L || !all(c("mean", "sd") %in% colnames(param.i)) ||
       any(!is.finite(param.i[, c("mean", "sd"), drop = FALSE])) ||
       any(param.i[, "sd"] <= 0)){
      stop("RMixtComp SEM returned invalid Gaussian parameters.")
    }
    means[, i] <- param.i[, "mean"]
    for(l in seq_len(L)){
      covs[i, i, l] <- param.i[l, "sd"]^2
    }
  }

  covariance.repair <- .stabilize.LPA.covariances(
    covs, "V0", fallback = stats::cov(response)
  )
  covs <- covariance.repair$covs
  expectation <- lpa_expectation_cpp(
    response, means, covs, P.Z, repair = FALSE
  )
  if(!isTRUE(expectation$valid)){
    stop("RMixtComp returned invalid stabilized LPA covariance matrices.")
  }
  P.Z.Xn <- expectation$posterior

  list(
    fit = fit.cur,
    means = means,
    covs = covs,
    P.Z = P.Z,
    P.Z.Xn = P.Z.Xn,
    Z = max.col(P.Z.Xn, ties.method = "first"),
    Log.Lik = expectation$Log.Lik,
    covariance.repaired = covariance.repair$repaired,
    covariance.repair = covariance.repair$diagnostics
  )
}

RMixtComp.LPA <- function(response, L, constraint = "V0", vis = TRUE,
                          control.RMixtComp){
  if(!requireNamespace("RMixtComp", quietly = TRUE) ||
     !requireNamespace("RMixtCompUtilities", quietly = TRUE)){
    stop(
      "method = 'RMixtComp' requires the optional packages 'RMixtComp' and ",
      "'RMixtCompUtilities'.",
      call. = FALSE
    )
  }

  response <- as.matrix(response)
  N <- nrow(response)
  I <- ncol(response)
  if(anyNA(response) || any(!is.finite(response))){
    stop("RMixtComp SEM requires complete finite LPA responses.")
  }
  if(!is.character(constraint) || length(constraint) != 1L || constraint != "V0"){
    stop("RMixtComp SEM for LPA supports only constraint = 'V0'.")
  }
  if(N < L){
    stop("RMixtComp requires N >= L so every latent profile can be non-empty.")
  }

  variable.names <- if(is.null(colnames(response))){
    paste0("V", seq_len(I))
  }else{
    make.unique(colnames(response))
  }
  data.RMixtComp <- as.data.frame(response, check.names = FALSE)
  names(data.RMixtComp) <- variable.names
  model <- as.list(setNames(rep("Gaussian", I), variable.names))

  fit <- .RMixtComp.fit(
    data.RMixtComp, model, L, control.RMixtComp
  )
  result <- .extract.RMixtComp.LPA(fit, response, L, variable.names)

  means <- result$means
  covs <- result$covs
  P.Z <- result$P.Z
  P.Z.Xn <- result$P.Z.Xn
  Z <- result$Z
  Log.Lik <- result$Log.Lik
  npar <- get.npar.LPA(I, L, constraint)
  AIC <- -2 * Log.Lik + 2 * npar
  BIC <- -2 * Log.Lik + npar * log(N)

  if(vis){
    .print.estimation.summary("RMixtComp SEM", Log.Lik, BIC)
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
    Log.Lik.nrep = Log.Lik,
    covariance.repaired = result$covariance.repaired,
    covariance.repair = result$covariance.repair
  )
}
