Rmixmod.LCA <- function(response, L, poly.value, nrep = 20,
                        starts = 100, maxiter.warmup = 20, vis = TRUE,
                        control.Rmixmod){
  if(!requireNamespace("Rmixmod", quietly = TRUE)){
    stop("method = 'Rmixmod' requires the optional package 'Rmixmod'.", call. = FALSE)
  }
  response <- as.matrix(response)
  N <- nrow(response)
  I <- ncol(response)
  poly.value <- as.integer(poly.value)
  response.int <- matrix(as.integer(response), N, I)

  invalid <- vapply(seq_len(I), function(i){
    any(!is.finite(response[, i])) ||
      any(!(response[, i] %in% 0:(poly.value[i] - 1L)))
  }, logical(1))
  if(any(invalid)){
    stop("Rmixmod requires complete LCA responses coded from 0 to C_i - 1.")
  }

  data.Rmixmod <- as.data.frame(lapply(seq_len(I), function(i){
    factor(response[, i], levels = 0:(poly.value[i] - 1L))
  }), check.names = FALSE)
  names(data.Rmixmod) <- if(is.null(colnames(response))) paste0("I", seq_len(I)) else colnames(response)

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
    stop("Rmixmod requires N >= L so every latent class can be non-empty.")
  }
  model <- Rmixmod::mixmodMultinomialModel(listModels = "Binary_pk_Ekjh")

  seed <- sample.int(.Machine$integer.max, 1L)

  random.parameters <- function(){
    poly.max <- max(poly.value)
    par <- array(0, dim = c(L, I, poly.max))
    for(i in seq_len(I)){
      par[, i, seq_len(poly.value[i])] <-
        rdirichlet(L, rep(3, poly.value[i]))
    }
    list(
      par = par,
      P.Z = as.numeric(rdirichlet(1L, rep(3, L)))
    )
  }

  estimate.parameters <- function(labels){
    maximization <- lca_maximization_cpp(
      response.int, .class.indicator(labels, L), poly.value, 1e-8
    )
    list(par = maximization$par, P.Z = maximization$P.Z)
  }

  expectation.step <- function(parameters){
    result <- lca_expectation_cpp(response.int, parameters$par, parameters$P.Z)
    if(!isTRUE(result$valid)) stop("Invalid LCA expectation step")
    list(P.Z.Xn = result$posterior, Log.Lik = result$Log.Lik)
  }

  stochastic.step <- function(P.Z.Xn){
    cumulative <- P.Z.Xn
    if(L > 1L){
      for(l in 2:L){
        cumulative[, l] <- cumulative[, l - 1L] + cumulative[, l]
      }
    }
    labels <- pmin(rowSums(cumulative < runif(N)) + 1L, L)
    frequency <- tabulate(labels, nbins = L)
    while(any(frequency == 0L)){
      target <- which(frequency == 0L)[1]
      candidates <- which(frequency[labels] > 1L)
      if(!length(candidates)){
        stop("SEM warm-up could not maintain non-empty latent classes.")
      }
      selected <- candidates[which.max(P.Z.Xn[candidates, target])]
      donor <- labels[selected]
      labels[selected] <- target
      frequency[donor] <- frequency[donor] - 1L
      frequency[target] <- 1L
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
    list(Log.Lik = estep$Log.Lik, labels = labels)
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
      center.cur <- as.matrix(parameters.cur["center"])
      scatter.cur <- parameters.cur["scatter"]
      P.Z.cur <- as.numeric(parameters.cur["proportions"])
      post.cur <- as.matrix(fit.cur["proba"])

      if(!identical(result.cur@error, "No error") ||
         !is.finite(as.numeric(result.cur@likelihood)) ||
         !identical(dim(center.cur), c(as.integer(L), I)) ||
         length(scatter.cur) != L ||
         length(P.Z.cur) != L || any(!is.finite(P.Z.cur)) || any(P.Z.cur <= 0) ||
         !identical(dim(post.cur), c(N, as.integer(L))) ||
         any(!is.finite(post.cur)) || any(rowSums(post.cur) <= 0)){
        return(FALSE)
      }

      for(l in seq_len(L)){
        scatter.l <- as.matrix(scatter.cur[[l]])
        for(i in seq_len(I)){
          prob <- scatter.l[i, seq_len(poly.value[i])]
          center.li <- as.integer(center.cur[l, i])
          if(any(!is.finite(prob)) || center.li < 1L || center.li > poly.value[i]){
            return(FALSE)
          }
          prob[center.li] <- 1 - prob[center.li]
          if(!is.finite(sum(prob)) || sum(prob) <= 0){
            return(FALSE)
          }
        }
      }
      TRUE
    }, error = function(e){
      FALSE
    })
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
      fit <- fit.cur
      best.Log.Lik <- as.numeric(fit.cur["likelihood"])
    }
    Log.Lik.nrep <- best.Log.Lik
  }else{
    warmup.Log.Lik <- rep(-Inf, starts)
    warmup.labels <- vector("list", starts)
    warmup.progress.state <- .new.progress.state()
    for(s in seq_len(starts)){
      seed.cur <- as.integer((as.double(seed) + s - 2) %% .Machine$integer.max + 1)
      warmup.cur <- tryCatch(
        run.sem.warmup(seed.cur),
        error = function(e){ NULL }
      )

      if(!is.null(warmup.cur) && is.finite(warmup.cur$Log.Lik)){
        warmup.Log.Lik[s] <- warmup.cur$Log.Lik
        warmup.labels[[s]] <- warmup.cur$labels
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
        " finite LCA solutions from ", starts, " starts; nrep = ", nrep, "."
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
        initMethod = "label",
        labels = warmup.labels[[warmup.position[r]]],
        nbIterationInAlgo = as.integer(maxiter)
      )
      methods::slot(strategy, "initMethod") <- "partition"
      fit.cur <- run.mixmod(strategy, seed.cur)

      if(!is.null(fit.cur) && is.valid.fit(fit.cur)){
        Log.Lik.nrep[r] <- as.numeric(fit.cur["likelihood"])
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
         " failed to produce a finite LCA solution. Sample size may be too small.")
  }

  parameters <- fit["parameters"]
  center <- as.matrix(parameters["center"])
  scatter <- parameters["scatter"]
  poly.max <- max(poly.value)
  par <- array(NA_real_, dim = c(L, I, poly.max))
  for(l in seq_len(L)){
    scatter.l <- as.matrix(scatter[[l]])
    for(i in seq_len(I)){
      prob <- scatter.l[i, seq_len(poly.value[i])]
      center.li <- as.integer(center[l, i])
      prob[center.li] <- 1 - prob[center.li]
      par[l, i, seq_len(poly.value[i])] <- prob / sum(prob)
    }
  }

  P.Z <- as.numeric(parameters["proportions"])
  P.Z <- P.Z / sum(P.Z)
  P.Z.Xn <- as.matrix(fit["proba"])
  P.Z.Xn <- P.Z.Xn / rowSums(P.Z.Xn)
  Z <- max.col(P.Z.Xn, ties.method = "first")
  Log.Lik <- get.Log.Lik.LCA(response, par, P.Z)
  npar <- get.npar.LCA(poly.value, L)
  AIC <- -2 * Log.Lik + 2 * npar
  BIC <- -2 * Log.Lik + npar * log(N)

  if(vis){
    .print.estimation.summary(
      paste("Rmixmod", paste(algorithm, collapse = " -> ")),
      Log.Lik, BIC
    )
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
    Log.Lik.nrep = Log.Lik.nrep
  )
}
