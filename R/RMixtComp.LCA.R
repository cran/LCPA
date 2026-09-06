.extract.RMixtComp.LCA <- function(fit.cur, response, L, poly.value,
                                   variable.names){
  N <- nrow(response)
  I <- ncol(response)
  P.Z <- as.numeric(RMixtCompUtilities::getProportion(fit.cur))
  if(length(P.Z) != L || any(!is.finite(P.Z)) || any(P.Z <= 0)){
    stop("RMixtComp SEM returned invalid latent class proportions.")
  }
  P.Z <- P.Z / sum(P.Z)

  P.Z.Xn <- as.matrix(RMixtCompUtilities::getTik(fit.cur, log = FALSE))
  if(!identical(dim(P.Z.Xn), c(N, as.integer(L))) ||
     any(!is.finite(P.Z.Xn)) || any(rowSums(P.Z.Xn) <= 0)){
    stop("RMixtComp SEM returned invalid posterior class probabilities.")
  }
  P.Z.Xn <- P.Z.Xn / rowSums(P.Z.Xn)

  poly.max <- max(poly.value)
  par <- array(NA_real_, dim = c(L, I, poly.max))
  for(i in seq_len(I)){
    par.i <- as.matrix(RMixtCompUtilities::getParam(fit.cur, variable.names[i]))
    modality.names <- paste0("modality ", seq.int(0, poly.value[i] - 1L))
    if(all(modality.names %in% colnames(par.i))){
      par.i <- par.i[, modality.names, drop = FALSE]
    }else{
      par.i <- par.i[, seq_len(poly.value[i]), drop = FALSE]
    }
    if(nrow(par.i) != L || any(!is.finite(par.i)) || any(par.i < 0) ||
       any(rowSums(par.i) <= 0)){
      stop("RMixtComp SEM returned invalid multinomial parameters.")
    }
    par[, i, seq_len(poly.value[i])] <- par.i / rowSums(par.i)
  }

  list(
    fit = fit.cur,
    par = par,
    P.Z = P.Z,
    P.Z.Xn = P.Z.Xn,
    Z = max.col(P.Z.Xn, ties.method = "first"),
    Log.Lik = get.Log.Lik.LCA(response, par, P.Z)
  )
}

RMixtComp.LCA <- function(response, L, poly.value, vis = TRUE,
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
    stop("RMixtComp SEM requires complete finite LCA responses.")
  }
  if(N < L){
    stop("RMixtComp requires N >= L so every latent class can be non-empty.")
  }

  variable.names <- if(is.null(colnames(response))){
    paste0("I", seq_len(I))
  }else{
    make.unique(colnames(response))
  }
  data.RMixtComp <- as.data.frame(lapply(seq_len(I), function(i){
    factor(response[, i], levels = seq.int(0, poly.value[i] - 1L))
  }), check.names = FALSE)
  names(data.RMixtComp) <- variable.names
  model <- as.list(setNames(rep("Multinomial", I), variable.names))

  fit <- .RMixtComp.fit(
    data.RMixtComp, model, L, control.RMixtComp
  )
  result <- .extract.RMixtComp.LCA(
    fit, response, L, poly.value, variable.names
  )

  par <- result$par
  P.Z <- result$P.Z
  P.Z.Xn <- result$P.Z.Xn
  Z <- result$Z
  Log.Lik <- result$Log.Lik
  npar <- get.npar.LCA(poly.value, L)
  AIC <- -2 * Log.Lik + 2 * npar
  BIC <- -2 * Log.Lik + npar * log(N)

  if(vis){
    .print.estimation.summary("RMixtComp SEM", Log.Lik, BIC)
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
    Log.Lik.nrep = Log.Lik
  )
}
