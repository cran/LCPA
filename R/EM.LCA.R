EM.LCA <- function(response, L=5, par.ini=NULL, nrep=1, starts=50, maxiter.wa=100, vis=FALSE, maxiter=2000, tol = 1e-4){

  adjust.response.obj <- adjust.response(response)
  response <- adjust.response.obj$response
  poly.max <- adjust.response.obj$poly.max
  poly.value <- adjust.response.obj$poly.value
  poly.orig <- adjust.response.obj$poly.orig

  npar <- sum(poly.value * L - 1) + L-1

  Y <- as.matrix(response)
  N <- nrow(Y)
  I <- ncol(Y)
  Y_hot <- array(0, dim = c(N, I, poly.max))
  idx <- cbind(c(row(Y)), c(col(Y)), c(Y+1))
  Y_hot[idx] <- 1

  int_width <- ceiling(log10(N * I * L))
  total_width <- int_width + 5
  fmt_string_maxchg <- sprintf("%%%d.%df", total_width, 5)

  int_width <- ceiling(log10(N * I * L)) + 1L
  total_width <- int_width + 3
  fmt_string_BIC <- sprintf("%%%d.%df", total_width, 2)

  run_em_once <- function(init_method, r=0, nrep=0, best_BIC=Inf, wa=FALSE) {
    par.cur <- NULL
    retry_count <- 0
    max_retry <- 100

    while(is.null(par.cur) && retry_count < max_retry){
      retry_count <- retry_count + 1
      tryCatch({
        if(is.character(init_method)){
          if(init_method == "kmeans"){
            Kmeans.LCA.obj <- Kmeans.LCA(response, L, nrep=100)
            par.cur <- par.pre <- Kmeans.LCA.obj$params$par
            P.Z.pre <- Kmeans.LCA.obj$params$P.Z
          } else if(init_method == "random"){
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
          par.cur <- par.pre <- init_method$par
          P.Z.pre <- init_method$P.Z
        }

        P.Z.Xn <- matrix(1/L, N, L)
        Rlj <- Ilj <- array(NA, dim=c(L, I, poly.max))

        iter <- 0
        L.X.pre <- -Inf
        Log.Lik.history <- numeric(if(wa) maxiter.wa else maxiter)

        maxiter.once <- if(wa) maxiter.wa else maxiter

        while(iter < maxiter.once){
          iter <- iter + 1
          par_pre_vec <- as.vector(par.pre)
          P_Z_pre_vec <- as.vector(P.Z.pre)
          Y_int <- matrix(as.integer(Y), nrow = N, ncol = I)

          e_step_res <- em_e_step(Y_int, par_pre_vec, P_Z_pre_vec, N, I, L, poly.max)
          L.Xi.Z <- e_step_res$L_Xi_Z
          P.Z.Xn <- e_step_res$P_Z_Xn
          L.Xi_vec <- e_step_res$L_Xi_vec

          Y_hot_vec <- as.vector(Y_hot)

          m_step_res <- em_m_step(P.Z.Xn, Y_hot_vec, poly.value, N, I, L, poly.max)
          Rlj <- array(m_step_res$Rlj, dim = c(L, I, poly.max))
          Ilj <- array(m_step_res$Ilj, dim = c(L, I, poly.max))

          par.cur <- (Rlj + 1e-10) / (Ilj + 2e-10)

          L.X.cur <- sum(log(L.Xi_vec))
          Log.Lik.history[iter] <- L.X.cur

          maxchg <- if(iter > 1) abs(L.X.cur - L.X.pre) else Inf
          AIC <- -2*L.X.cur + 2*npar
          BIC <- -2*L.X.cur + npar*log(N)
          Deviance <- -2*L.X.cur
          if(BIC < best_BIC){
            best_BIC <- BIC
          }

          if(vis && iter > 1 && r){
            if(wa){
              cat('\rWarm', paste0(sprintf("%2d", r), "/", sprintf("%2d", if(wa) starts else nrep)), '| Iter =', sprintf("%4d", iter),
                  '  \u0394Log.Lik =', sprintf(fmt_string_maxchg, maxchg),
                  '  BIC =', sprintf(fmt_string_BIC, BIC), '  Best_BIC =', sprintf(fmt_string_BIC, best_BIC))
            } else {
              cat('\rRep ', paste0(sprintf("%2d", r), "/", sprintf("%2d", nrep)), '| Iter =', sprintf("%4d", iter),
                  '  \u0394Log.Lik =', sprintf(fmt_string_maxchg, maxchg),
                  '  BIC =', sprintf(fmt_string_BIC, BIC), '  Best_BIC =', sprintf(fmt_string_BIC, best_BIC))
            }
          } else if(vis && iter > 1){
            cat('\rIter =', sprintf("%4d", iter), '  \u0394Log.Lik =', sprintf(fmt_string_maxchg, maxchg),
                '  BIC =', sprintf(fmt_string_BIC, BIC))
          }

          if(maxchg < tol)
            break

          par.pre <- par.cur
          P.Z.pre <- colSums(P.Z.Xn) / sum(P.Z.Xn)
          L.X.pre <- L.X.cur

          if (iter == maxiter.once && vis && r == 0){
            message('\nMaximum number of iterations reached; convergence may not have been achieved\n')
          }
        }

        Log.Lik <- L.X.cur

        P.Z.pre <- as.table(P.Z.pre)
        names(P.Z.pre) <- paste0("Class ", 1:L)

        res = list(
          params = list(par = par.cur, P.Z = P.Z.pre),
          npar = npar,
          Log.Lik = Log.Lik,
          AIC=AIC,
          BIC=BIC,
          best_BIC=best_BIC,
          P.Z.Xn = P.Z.Xn,
          P.Z = P.Z.pre,
          Z = apply(P.Z.Xn, 1, which.max),
          probability = NULL,
          Log.Lik.history = Log.Lik.history[1:iter]
        )

        return(res)
      }, error = function(e){
        par.cur <- NULL
      })
    }

    if(is.null(par.cur)) {
      stop("Failed to initialize after ", max_retry, " attempts. Try different initial values.")
    }
  }

  if (is.null(par.ini)){
    par.ini <- "random"
  }

  results.wa <- list(params = list(), BIC = numeric(0))
  best_BIC_wa <- Inf

  if (starts >= 1 && any(par.ini == "random")) {

    for (s in 1:starts) {
      res_wa <- run_em_once("random", r = s, nrep = starts, best_BIC = best_BIC_wa, wa = TRUE)
      best_BIC_wa <- res_wa$best_BIC

      if (length(results.wa$BIC) < nrep) {
        results.wa$params[[length(results.wa$BIC) + 1]] <- res_wa$params
        results.wa$BIC <- c(results.wa$BIC, res_wa$BIC)
      } else {
        worst_idx <- which.max(results.wa$BIC)
        if (res_wa$BIC < results.wa$BIC[worst_idx]) {
          results.wa$params[[worst_idx]] <- res_wa$params
          results.wa$BIC[worst_idx] <- res_wa$BIC
        }
      }
    }
    if(vis){
      cat("\n")
    }
  }

  if (!is.null(results.wa$BIC) && length(results.wa$BIC) > 0 && any(par.ini == "random")) {
    results <- vector("list", nrep)
    Log.Lik.nrep <- numeric(nrep)
    best_BIC_main <- Inf

    for (r in 1:nrep) {
      init_vals <- list(par = results.wa$params[[r]]$par, P.Z = results.wa$params[[r]]$P.Z)
      res <- run_em_once(init_vals, r = r, nrep = nrep, best_BIC = best_BIC_main, wa = FALSE)
      results[[r]] <- res
      Log.Lik.nrep[r] <- res$Log.Lik
      if (res$BIC < best_BIC_main) best_BIC_main <- res$BIC
    }

    best_idx <- which.max(Log.Lik.nrep)
    res <- results[[best_idx]]
    res$Log.Lik.nrep <- Log.Lik.nrep

  } else if (any(par.ini == "random")) {
    results <- vector("list", nrep)
    Log.Lik.nrep <- numeric(nrep)
    best_BIC <- Inf

    for (r in 1:nrep) {
      res <- run_em_once("random", r = r, nrep = nrep, best_BIC = best_BIC)
      results[[r]] <- res
      Log.Lik.nrep[r] <- res$Log.Lik
      best_BIC <- res$best_BIC
    }

    best_idx <- which.max(Log.Lik.nrep)
    res <- results[[best_idx]]
    res$Log.Lik.nrep <- Log.Lik.nrep

  } else if (any(par.ini == "kmeans")) {
    Kmeans.LCA.obj <- Kmeans.LCA(response, L, nrep = 100)
    par.ini.KM <- list(par = Kmeans.LCA.obj$params$par, P.Z = Kmeans.LCA.obj$params$P.Z)
    res <- run_em_once(par.ini.KM)
    res$Log.Lik.nrep <- res$Log.Lik

  } else {
    res <- run_em_once(par.ini)
    res$Log.Lik.nrep <- res$Log.Lik
  }

  res$probability <- NULL

  names(res$params$P.Z) <- 1:L

  if (vis) cat("\n\n")
  return(res)
}
