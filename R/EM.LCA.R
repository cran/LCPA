EM.LCA <- function(response, L=5, par.ini="random", starts=50, maxiter.warmup=100,
                   nrep=1, vis=FALSE, maxiter=2000, tol = 1e-4){

  .validate.training.stages(starts, maxiter.warmup, nrep)
  adjust.response.obj <- adjust.response(response)
  response <- adjust.response.obj$response
  poly.max <- adjust.response.obj$poly.max
  poly.value <- adjust.response.obj$poly.value
  poly.orig <- adjust.response.obj$poly.orig

  npar <- get.npar.LCA(poly.value, L)

  Y <- as.matrix(response)
  N <- nrow(Y)
  I <- ncol(Y)
  Y.int <- matrix(as.integer(Y), N, I)
  stabilize.posterior <- function(posterior){
    empty <- !is.finite(colSums(posterior)) | colSums(posterior) <= 1e-10
    if(any(empty)){
      posterior[, empty] <- pmax(posterior[, empty, drop = FALSE], min(1e-8, 0.01 / L))
      posterior <- posterior / rowSums(posterior)
    }
    posterior
  }

  run.EM <- function(par.ini.current, r=0, best_BIC=Inf, warmup=FALSE) {
    par.cur <- NULL
    attempts <- 0
    maxattempts <- 100
    iteration.progress.state <- .new.progress.state()

    while(is.null(par.cur) && attempts < maxattempts){
      attempts <- attempts + 1
      tryCatch({
        if(is.character(par.ini.current)){
          if(par.ini.current == "kmeans"){
            Kmeans.LCA.obj <- Kmeans.LCA(response, L, starts=1)
            par.cur <- par.pre <- Kmeans.LCA.obj$params$par
            P.Z.pre <- Kmeans.LCA.obj$params$P.Z
          } else if(par.ini.current == "random"){
            par.pre <- array(NA, dim=c(L, I, poly.max))
            for(i in 1:I){
              for(l in 1:L){
                par.pre[l , i, 1:poly.value[i]] <- rdirichlet(n = 1, alpha = rep(3, poly.value[i]))
              }
            }
            par.cur <- par.pre
            P.Z.pre <- rdirichlet(n = 1, alpha = rep(3, L))
          }
        } else {
          par.cur <- par.pre <- par.ini.current$par
          P.Z.pre <- par.ini.current$P.Z
        }

        maxiter.current <- if(warmup) maxiter.warmup else maxiter
        expectation <- lca_expectation_cpp(Y.int, par.pre, as.numeric(P.Z.pre))
        if(!isTRUE(expectation$valid)) stop("Invalid initial LCA solution")
        P.Z.Xn <- stabilize.posterior(expectation$posterior)
        L.X.pre <- expectation$Log.Lik
        Log.Lik.history <- numeric(maxiter.current + 1L)
        Log.Lik.history[1L] <- L.X.pre
        iter <- 0

        while(iter < maxiter.current){
          iter <- iter + 1
          maximization <- lca_maximization_cpp(Y.int, P.Z.Xn, poly.value, 1e-10)
          par.cur <- maximization$par
          P.Z.cur <- maximization$P.Z

          expectation <- lca_expectation_cpp(Y.int, par.cur, P.Z.cur)
          if(!isTRUE(expectation$valid)) stop("Invalid LCA expectation step")
          P.Z.Xn.cur <- stabilize.posterior(expectation$posterior)
          L.X.cur <- expectation$Log.Lik

          Log.Lik.history[iter + 1L] <- L.X.cur

          maxchg <- abs(L.X.cur - L.X.pre)
          AIC <- -2*L.X.cur + 2*npar
          BIC <- -2*L.X.cur + npar*log(N)
          Deviance <- -2*L.X.cur
          if(BIC < best_BIC){
            best_BIC <- BIC
          }

          if(vis && iter > 1 && r == 0){
            .print.iteration.progress(
              iter, maxchg, BIC,
              progress.prefix = .estimation.output.prefix(),
              progress.state = iteration.progress.state
            )
          }

          P.Z.Xn <- P.Z.Xn.cur
          if(maxchg < tol)
            break

          par.pre <- par.cur
          P.Z.pre <- P.Z.cur
          L.X.pre <- L.X.cur

          if (iter == maxiter.current && vis && r == 0){
            message(
              '\n', .estimation.output.prefix(),
              'Maximum number of iterations reached; convergence may not have been achieved\n'
            )
          }
        }

        Log.Lik <- L.X.cur

        P.Z.cur <- as.table(P.Z.cur)
        names(P.Z.cur) <- .latent.group.names(L, "LCA")

        res = list(
          params = list(par = par.cur, P.Z = P.Z.cur),
          npar = npar,
          Log.Lik = Log.Lik,
          AIC=AIC,
          BIC=BIC,
          best_BIC=best_BIC,
          P.Z.Xn = P.Z.Xn,
          P.Z = P.Z.cur,
          Z = max.col(P.Z.Xn, ties.method = "first"),
          probability = NULL,
          Log.Lik.history = Log.Lik.history[seq_len(iter + 1L)]
        )

        return(res)
      }, error = function(e){
        par.cur <- NULL
      })
    }

    if(is.null(par.cur)) {
      stop("Failed to initialize after ", maxattempts, " attempts. Try different initial values.")
    }
  }

  if (is.null(par.ini)){
    par.ini <- "random"
  }
  if(is.character(par.ini)){
    if(length(par.ini) != 1L || !par.ini %in% c("random", "kmeans")){
      stop("par.ini must be 'random', 'kmeans', or a parameter list.")
    }
  }else if(!is.list(par.ini) || !all(c("par", "P.Z") %in% names(par.ini))){
    stop("par.ini must contain par and P.Z.")
  }

  warmup.results <- list(params = vector("list", starts), Log.Lik = rep(-Inf, starts))
  best.BIC.warmup <- Inf
  warmup.progress.state <- .new.progress.state()
  for (s in seq_len(starts)) {
    res.warmup <- run.EM(
      par.ini, r = s, best_BIC = best.BIC.warmup, warmup = TRUE
    )
    best.BIC.warmup <- res.warmup$best_BIC
    warmup.results$params[[s]] <- res.warmup$params
    warmup.results$Log.Lik[s] <- res.warmup$Log.Lik
    if(vis){
      .print.estimation.progress(
        "Warm", s, starts, res.warmup$Log.Lik,
        algorithm = "EM",
        iterations = length(res.warmup$Log.Lik.history) - 1L,
        progress.state = warmup.progress.state
      )
    }
  }
  if(vis) .end.estimation.progress()

  warmup.position <- order(warmup.results$Log.Lik, decreasing = TRUE)
  warmup.position <- warmup.position[
    is.finite(warmup.results$Log.Lik[warmup.position])
  ]
  if(length(warmup.position) < nrep){
    stop(
      "EM warm-up produced only ", length(warmup.position),
      " finite LCA solutions from ", starts, " starts; nrep = ", nrep, "."
    )
  }
  warmup.position <- warmup.position[seq_len(nrep)]

  results <- vector("list", nrep)
  Log.Lik.nrep <- rep(-Inf, nrep)
  best.BIC.main <- Inf
  best.Log.Lik <- -Inf
  replication.progress.state <- .new.progress.state()
  for (r in seq_len(nrep)) {
    par.ini.current <- warmup.results$params[[warmup.position[r]]]
    res.cur <- run.EM(
      par.ini.current, r = r, best_BIC = best.BIC.main, warmup = FALSE
    )
    results[[r]] <- res.cur
    Log.Lik.nrep[r] <- res.cur$Log.Lik
    best.BIC.main <- min(best.BIC.main, res.cur$BIC)
    best.Log.Lik <- max(best.Log.Lik, Log.Lik.nrep[r])
    if(vis){
      .print.estimation.progress(
        "Rep", r, nrep, Log.Lik.nrep[r], best.Log.Lik,
        progress.state = replication.progress.state
      )
    }
  }
  if(vis) .end.estimation.progress()

  best.idx <- which.max(Log.Lik.nrep)
  res <- results[[best.idx]]
  res$Log.Lik.nrep <- Log.Lik.nrep

  res$probability <- NULL

  names(res$params$P.Z) <- 1:L

  if (vis) .print.estimation.summary("EM", res$Log.Lik, res$BIC)
  return(res)
}
