Rmixmod.LPA <- function(response, L, constraint = "VV", nrep = 20,
                        starts = 100, maxiter.warmup = 20, vis = TRUE,
                        control.Rmixmod){

  if(!requireNamespace("Rmixmod", quietly = TRUE)){
    stop("method = 'Rmixmod' requires the optional package 'Rmixmod'.", call. = FALSE)
  }
  response <- as.matrix(response)
  N <- nrow(response)
  I <- ncol(response)
  if(any(!is.finite(response))){
    stop("Rmixmod requires complete finite LPA responses.")
  }

  model.names <- c(
    V0 = "Gaussian_pk_Lk_Bk",
    EE = "Gaussian_pk_L_C",
    VV = "Gaussian_pk_Lk_Ck"
  )
  constraint <- as.character(constraint)[1]
  if(!(constraint %in% names(model.names))){
    stop("Rmixmod SEM for LPA supports constraints 'V0', 'EE', and 'VV'.")
  }

  nrep <- as.integer(nrep)
  starts <- as.integer(starts)
  maxiter.warmup <- as.integer(maxiter.warmup)
  path <- match.arg(as.character(control.Rmixmod$path)[1L], c("LCPA", "Rmixmod"))
  is.lcpa <- path == "LCPA"
  algorithm <- "SEM"
  if(is.lcpa && (nrep < 1L || starts < nrep || maxiter.warmup < 1L)){
    stop("Rmixmod requires starts >= nrep >= 1 and maxiter.warmup >= 1.")
  }
  if(N < L){
    stop("Rmixmod requires N >= L so every latent profile can be non-empty.")
  }
  model <- Rmixmod::mixmodGaussianModel(listModels = unname(model.names[constraint]))
  data.Rmixmod <- as.data.frame(response, check.names = FALSE)
  names(data.Rmixmod) <- if(is.null(colnames(response))) paste0("V", seq_len(I)) else colnames(response)

  seed <- sample.int(.Machine$integer.max, 1L)

  random.parameters <- function(){
    q10 <- apply(response, 2, stats::quantile, probs = 0.10,
                 names = FALSE, type = 7)
    q90 <- apply(response, 2, stats::quantile, probs = 0.90,
                 names = FALSE, type = 7)
    means <- matrix(NA_real_, L, I)
    for(i in seq_len(I)){
      means[, i] <- runif(L, q10[i], q90[i])
    }
    list(
      means = means,
      covs = replicate(L, diag(I), simplify = FALSE),
      P.Z = as.numeric(rdirichlet(1L, rep(3, L)))
    )
  }

  covariance.fallback <- matrix(stats::cov(response), I, I)
  if(any(!is.finite(covariance.fallback))) covariance.fallback <- diag(I)
  regularize.covariance <- function(covariance, diagonal=FALSE){
    covariance <- matrix(covariance, I, I)
    if(diagonal) covariance <- diag(diag(covariance), I)
    .stabilize.LPA.covariances(
      array(covariance, dim = c(I, I, 1L)), "VV",
      fallback = covariance.fallback
    )$covs[, , 1L]
  }

  estimate.parameters <- function(labels){
    frequency <- tabulate(labels, nbins = L)
    means <- matrix(NA_real_, L, I)
    covs.ini <- vector("list", L)
    for(l in seq_len(L)){
      response.l <- response[labels == l, , drop=FALSE]
      means[l, ] <- colMeans(response.l)
      dev <- sweep(response.l, 2, means[l, ], "-")
      covs.ini[[l]] <- regularize.covariance(
        crossprod(dev) / frequency[l], diagonal = constraint == "V0"
      )
    }
    if(constraint == "EE"){
      residual <- response - means[labels, , drop=FALSE]
      covariance <- regularize.covariance(crossprod(residual) / N)
      covs.ini <- replicate(L, covariance, simplify = FALSE)
    }
    list(means = means, covs = covs.ini, P.Z = frequency / N)
  }

  as.Rmixmod.parameter <- function(parameters){
    methods::new(
      "GaussianParameter",
      proportions = parameters$P.Z,
      mean = parameters$means,
      variance = parameters$covs
    )
  }

  expectation.step <- function(parameters){
    covs <- array(unlist(parameters$covs, use.names = FALSE), c(I, I, L))
    result <- lpa_expectation_cpp(
      response, parameters$means, covs, parameters$P.Z, repair = FALSE
    )
    if(!isTRUE(result$valid)) stop("Invalid LPA expectation step")
    list(P.Z.Xn = result$posterior, Log.Lik = result$Log.Lik)
  }

  stochastic.step <- function(P.Z.Xn){
    cumulative <- P.Z.Xn
    if(L > 1L){
      for(l in 2:L){
        cumulative[, l] <- cumulative[, l - 1L] + cumulative[, l]
      }
    }
    labels <- rowSums(cumulative < runif(N)) + 1L
    labels <- pmin(labels, L)

    frequency <- tabulate(labels, nbins = L)
    while(any(frequency == 0L)){
      target <- which(frequency == 0L)[1]
      candidates <- which(frequency[labels] > 1L)
      if(!length(candidates)){
        stop("SEM warm-up could not maintain non-empty latent profiles.")
      }
      selected <- candidates[which.max(P.Z.Xn[candidates, target])]
      donor <- labels[selected]
      labels[selected] <- target
      frequency[donor] <- frequency[donor] - 1L
      frequency[target] <- frequency[target] + 1L
    }
    labels
  }

  run.sem.warmup <- function(seed.cur){
    set.seed(seed.cur)
    parameters <- random.parameters()
    for(iter in seq_len(maxiter.warmup)){
      estep <- expectation.step(parameters)
      labels <- stochastic.step(estep$P.Z.Xn)
      parameters <- estimate.parameters(labels)
    }
    estep <- expectation.step(parameters)
    list(
      Log.Lik = estep$Log.Lik,
      parameter = as.Rmixmod.parameter(parameters)
    )
  }

  run.mixmod <- function(strategy, seed.cur, fail.silently = TRUE){
    run.cur <- function(){
      Rmixmod::mixmodCluster(
        data = data.Rmixmod,
        nbCluster = as.integer(L),
        models = model,
        strategy = strategy,
        criterion = "BIC",
        seed = seed.cur
      )
    }
    if(!fail.silently){
      return(run.cur())
    }
    tryCatch({
      withCallingHandlers(
        run.cur(),
        warning = function(w){
          if(grepl("All models got errors!", conditionMessage(w), fixed = TRUE)){
            invokeRestart("muffleWarning")
          }
        }
      )
    }, error = function(e){
      NULL
    })
  }

  is.valid.fit <- function(fit.cur){
    tryCatch({
      result.cur <- fit.cur["bestResult"]
      parameters.cur <- fit.cur["parameters"]
      means.cur <- as.matrix(parameters.cur["mean"])
      variance.cur <- parameters.cur["variance"]
      P.Z.cur <- as.numeric(parameters.cur["proportions"])
      post.cur <- as.matrix(fit.cur["proba"])

      identical(result.cur@error, "No error") &&
        is.finite(as.numeric(result.cur@likelihood)) &&
        identical(dim(means.cur), c(as.integer(L), I)) &&
        all(is.finite(means.cur)) &&
        length(variance.cur) == L &&
        all(vapply(variance.cur, function(x){ all(is.finite(as.matrix(x))) }, logical(1))) &&
        length(P.Z.cur) == L && all(is.finite(P.Z.cur)) && all(P.Z.cur > 0) &&
        identical(dim(post.cur), c(N, as.integer(L))) &&
        all(is.finite(post.cur)) && all(rowSums(post.cur) > 0)
    }, error = function(e){
      FALSE
    })
  }

  corrected.Log.Lik <- function(fit.cur){
    parameters.cur <- fit.cur["parameters"]
    means.cur <- as.matrix(parameters.cur["mean"])
    variance.cur <- parameters.cur["variance"]
    covs.cur <- array(NA_real_, dim = c(I, I, L))
    for(l in seq_len(L)) covs.cur[, , l] <- as.matrix(variance.cur[[l]])
    covariance.cur <- .stabilize.LPA.covariances(
      covs.cur, constraint, fallback = stats::cov(response)
    )$covs
    P.Z.cur <- as.numeric(parameters.cur["proportions"])
    P.Z.cur <- P.Z.cur / sum(P.Z.cur)
    expectation.cur <- lpa_expectation_cpp(
      response, means.cur, covariance.cur, P.Z.cur, repair = FALSE
    )
    if(!isTRUE(expectation.cur$valid)) -Inf else expectation.cur$Log.Lik
  }

  fit <- NULL
  best.Log.Lik <- -Inf
  if(!is.lcpa){
    strategy <- .Rmixmod.strategy(control.Rmixmod)
    if(!methods::is(strategy, "Strategy")){
      stop("control.Rmixmod$strategy must be an Rmixmod Strategy object.")
    }
    algorithm <- as.character(methods::slot(strategy, "algo"))
    fit.cur <- run.mixmod(strategy, as.integer(seed), fail.silently = FALSE)
    if(!is.null(fit.cur) && is.valid.fit(fit.cur)){
      Log.Lik.cur <- corrected.Log.Lik(fit.cur)
      if(is.finite(Log.Lik.cur)){
        fit <- fit.cur
        best.Log.Lik <- Log.Lik.cur
      }
    }
    Log.Lik.nrep <- best.Log.Lik
  }else{
    warmup.Log.Lik <- rep(-Inf, starts)
    warmup.parameters <- vector("list", starts)
    warmup.progress.state <- .new.progress.state()
    for(s in seq_len(starts)){
      seed.cur <- as.integer((as.double(seed) + s - 2) %% .Machine$integer.max + 1)
      warmup.cur <- tryCatch(
        run.sem.warmup(seed.cur),
        error = function(e){ NULL }
      )

      if(!is.null(warmup.cur) && is.finite(warmup.cur$Log.Lik)){
        warmup.Log.Lik[s] <- warmup.cur$Log.Lik
        warmup.parameters[[s]] <- warmup.cur$parameter
      }
      if(vis){
        .print.estimation.progress(
          "Warm", s, starts, warmup.Log.Lik[s],
          algorithm = "Rmixmod SEM", iterations = maxiter.warmup,
          progress.state = warmup.progress.state
        )
      }
    }
    if(vis){
      .end.estimation.progress()
    }

    warmup.position <- order(warmup.Log.Lik, decreasing = TRUE)
    warmup.position <- warmup.position[is.finite(warmup.Log.Lik[warmup.position])]
    if(length(warmup.position) < nrep){
      stop(
        "Rmixmod warm-up produced only ", length(warmup.position),
        " finite LPA solutions from ", starts, " starts; nrep = ", nrep, "."
      )
    }
    warmup.position <- warmup.position[seq_len(nrep)]

    Log.Lik.nrep <- rep(-Inf, nrep)
    replication.progress.state <- .new.progress.state()
    for(r in seq_len(nrep)){
      seed.cur <- as.integer((as.double(seed) + starts + r - 2) %% .Machine$integer.max + 1)
      maxiter <- control.Rmixmod$maxiter
      if(is.null(maxiter)) maxiter <- 1000L
      strategy <- Rmixmod::mixmodStrategy(
        algo = "SEM",
        nbTry = 1L,
        initMethod = "parameter",
        parameter = warmup.parameters[[warmup.position[r]]],
        nbIterationInAlgo = as.integer(maxiter)
      )
      fit.cur <- run.mixmod(strategy, seed.cur)

      if(!is.null(fit.cur) && is.valid.fit(fit.cur)){
        Log.Lik.nrep[r] <- corrected.Log.Lik(fit.cur)
        if(Log.Lik.nrep[r] > best.Log.Lik){
          fit <- fit.cur
          best.Log.Lik <- Log.Lik.nrep[r]
        }
      }
      if(vis){
        .print.estimation.progress(
          "Rep", r, nrep, Log.Lik.nrep[r], best.Log.Lik,
          progress.state = replication.progress.state
        )
      }
    }
    if(vis){
      .end.estimation.progress()
    }
  }
  if(is.null(fit)){
    stop("Rmixmod ", paste(algorithm, collapse = " -> "),
         " failed to produce a finite LPA solution. Sample size may be too small.")
  }

  parameters <- fit["parameters"]
  means <- as.matrix(parameters["mean"])
  variance <- parameters["variance"]
  covs <- array(NA_real_, dim = c(I, I, L))
  for(l in seq_len(L)){
    covs[, , l] <- as.matrix(variance[[l]])
  }

  P.Z <- as.numeric(parameters["proportions"])
  P.Z <- P.Z / sum(P.Z)
  covariance.repair <- .stabilize.LPA.covariances(
    covs, constraint, fallback = stats::cov(response)
  )
  covs <- covariance.repair$covs
  expectation <- lpa_expectation_cpp(
    response, means, covs, P.Z, repair = FALSE
  )
  if(!isTRUE(expectation$valid)){
    stop("Rmixmod returned invalid stabilized LPA covariance matrices.")
  }
  P.Z.Xn <- expectation$posterior
  Z <- max.col(P.Z.Xn, ties.method = "first")
  Log.Lik <- expectation$Log.Lik
  npar <- get.npar.LPA(I, L, constraint)
  AIC <- -2 * Log.Lik + 2 * npar
  BIC <- -2 * Log.Lik + npar * log(N)

  if(vis){
    .print.estimation.summary(
      paste("Rmixmod", paste(algorithm, collapse = " -> ")),
      Log.Lik, BIC
    )
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
    Log.Lik.nrep = Log.Lik.nrep,
    covariance.repaired = covariance.repair$repaired,
    covariance.repair = covariance.repair$diagnostics
  )
}
