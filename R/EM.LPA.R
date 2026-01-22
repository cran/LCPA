
#' @importFrom stats cov sd kmeans
EM.LPA <- function(response, L = 2, par.ini = NULL, constraint = "VV",
                   nrep = 1, starts=50, maxiter.wa=100, vis = TRUE,
                   maxiter = 2000, tol = 1e-4) {

  if (is.vector(response) && is.numeric(response)) response <- matrix(response, ncol = 1)

  if (is.data.frame(response)) response <- as.matrix(response)
  if (!is.matrix(response)) stop("response must be matrix/data.frame/vector")
  N <- nrow(response)
  I <- ncol(response)
  if (N == 0) stop("empty response")
  jitter = 1e-10

  int_width <- ceiling(log10(N * I * L))
  total_width <- int_width + 5
  fmt_string_maxchg <- sprintf("%%%d.%df", total_width, 5)

  int_width <- ceiling(log10(N * I * L)) + 1L
  total_width <- int_width + 3
  fmt_string_BIC <- sprintf("%%%d.%df", total_width, 2)

  if (is.null(par.ini)) {
    par.ini <- "random"
  }

  covs.global <- tryCatch({
    cov(response)
  }, error = function(e) diag(pmax(apply(response, 2, var, na.rm = TRUE), jitter)))

  if (any(!is.finite(covs.global)) || !is.matrix(covs.global)) {
    vars <- apply(response, 2, var, na.rm = TRUE)
    if (any(!is.finite(vars))) vars <- rep(1, I)
    covs.global <- diag(pmax(vars, jitter))
  }
  covs.global <- (covs.global + t(covs.global)) / 2
  eig0 <- eigen(covs.global, symmetric = TRUE)
  if (min(eig0$values) < jitter) {
    covs.global <- eig0$vectors %*% diag(pmax(eig0$values, jitter)) %*% t(eig0$vectors)
    covs.global <- (covs.global + t(covs.global)) / 2
  }
  covs.global <- as.matrix(covs.global)

  npar <- get.npar.LPA(I, L, constraint)

  results.wa <- list(params=NULL, BIC=NULL)
  best_BIC <- Inf

  run_em_once <- function(init_vals, r = 0, best_BIC = Inf, wa=FALSE) {
    means.cur <- init_vals$means
    covs.cur <- init_vals$covs
    P.Z.cur <- init_vals$P.Z

    if(wa){
      maxiter.once <- maxiter.wa
    }else{
      maxiter.once <- maxiter
    }

    if(wa){
      nrep.once <- starts
    }else{
      nrep.once <- nrep
    }

    Log.Lik <- -Inf
    Log.Lik.history <- numeric(maxiter)
    P.Z.Xn <- matrix(1/L, nrow = N, ncol = L)

    prev_means <- means.cur
    prev_covs <- covs.cur
    prev_P.Z <- P.Z.cur

    tresponse <- t(response)

    for (iter in 1:maxiter.once) {
      logres <- matrix(NA_real_, nrow = N, ncol = L)

      for (l in 1:L) {
        covs.l <- covs.cur[,,l]
        mean.l <- means.cur[l,]
        lp <- logpdf_component(mean.l, covs.l, tresponse, jitter)
        logres[,l] <- lp + log(P.Z.cur[l] + 1e-12)
      }

      rowmax <- apply(logres, 1, max)
      nonfinite_rows <- !is.finite(rowmax)
      if (any(nonfinite_rows)) {
        finite_vals <- logres[is.finite(logres)]
        if (length(finite_vals) > 0) {
          rowmax[nonfinite_rows] <- max(finite_vals)
        } else {
          rowmax[nonfinite_rows] <- 0
        }
      }
      exp_rel <- exp(logres - matrix(rowmax, N, L))
      row_sums <- rowSums(exp_rel)
      bad <- !is.finite(row_sums) | row_sums < 1e-20
      if (any(bad)) {
        exp_rel[bad,] <- 1/L
        row_sums[bad] <- 1
      }
      P.Z.Xn <- exp_rel / row_sums

      nk <- colSums(P.Z.Xn)
      empty_clusters <- which(nk < 1e-5)
      if (length(empty_clusters) > 0) {
        non_empty <- setdiff(1:L, empty_clusters)
        P.Z.Xn[, empty_clusters] <- 0
        row_sums <- rowSums(P.Z.Xn[, non_empty, drop = FALSE])
        row_sums[row_sums < 1e-20] <- 1
        P.Z.Xn[, non_empty] <- P.Z.Xn[, non_empty, drop = FALSE] / row_sums
        nk <- colSums(P.Z.Xn)
      }

      P.Z.cur <- pmax(nk / N, 1e-12)
      P.Z.cur <- P.Z.cur / sum(P.Z.cur)

      means_new <- (t(P.Z.Xn) %*% response) / matrix(nk + 1e-12, nrow = L, ncol = I)
      nonfinite_means <- apply(means_new, 1, function(x) any(!is.finite(x)))
      if (any(nonfinite_means)) {
        means_new[nonfinite_means, ] <- prev_means[nonfinite_means, ]
      }
      means.cur <- means_new

      covs.new <- array(0, dim = c(I, I, L))
      for (l in 1:L) {
        wc <- P.Z.Xn[, l]
        totw <- sum(wc)
        mean.l <- means.cur[l,]
        dev <- sweep(response, 2, mean.l, "-")
        WDev <- dev * wc
        S.l <- crossprod(dev, WDev) / (totw + 1e-12)
        diag(S.l) <- pmax(diag(S.l), jitter)
        covs.new[,,l] <- matrix(S.l, nrow = I, ncol = I)
      }

      if(I == 1 && any(constraint %in% c("UE", "UV"))) {
        if(constraint == "UE") {
          var.total <- sum(nk * sapply(1:L, function(l) covs.new[,,l])) / sum(nk)
          var.total <- pmax(var.total, jitter)
          for(l in 1:L) {
            covs.new[,,l] <- matrix(var.total, nrow = I, ncol = I)
          }
        } else if(constraint == "UV") {
          for(l in 1:L) {
            covs.new[,,l] <- matrix(pmax(covs.new[,,l], jitter), nrow = I, ncol = I)
          }
        }
      } else if(any(constraint %in% c("E0", "V0", "EE", "VV", "VE", "EV"))) {
        if (constraint == "E0") {
          S.l <- array(0, dim = c(I, I))
          covs.new <- array(0, dim = c(I, I, L))
          totw = 0
          for (l in 1:L) {
            wc <- P.Z.Xn[, l]
            totw <- totw + sum(wc)

            mean.l <- means.cur[l, ]
            dev <- sweep(response, 2, mean.l, "-")
            WDev <- dev * wc
            S.l <- S.l + crossprod(dev, WDev)
          }
          S.l <- S.l / totw
          for (l in 1:L) diag(covs.new[, , l]) <- diag(S.l)
        } else if (constraint == "V0") {
          covs.new <- array(0, dim = c(I, I, L))
          for (l in 1:L) {
            wc <- P.Z.Xn[, l]
            totw <- sum(wc)
            mean.l <- means.cur[l,]
            dev <- sweep(response, 2, mean.l, "-")
            WDev <- dev * wc
            S.l <- crossprod(dev, WDev) / (totw + 1e-12)
            diag(S.l) <- pmax(diag(S.l), jitter)
            diag(covs.new[, , l]) <- diag(S.l)
          }
        } else if (constraint == "EE") {
          S.l <- array(0, dim = c(I, I))
          covs.new <- array(0, dim = c(I, I, L))
          totw = 0
          for (l in 1:L) {
            wc <- P.Z.Xn[, l]
            totw <- totw + sum(wc)

            mean.l <- means.cur[l, ]
            dev <- sweep(response, 2, mean.l, "-")
            WDev <- dev * wc
            S.l <- S.l + crossprod(dev, WDev)
          }
          S.l <- S.l / totw
          for (l in 1:L) covs.new[, , l] <- S.l
        } else if (constraint == "VE" |constraint == "EV") {

          S.l <- array(0, dim = c(I, I))
          covs.EE <- covs.VV <- array(0, dim = c(I, I, L))
          totw = 0
          for (l in 1:L) {
            wc <- P.Z.Xn[, l]
            totw <- totw + sum(wc)

            mean.l <- means.cur[l, ]
            dev <- sweep(response, 2, mean.l, "-")
            WDev <- dev * wc
            S.l <- S.l + crossprod(dev, WDev)
          }
          S.l <- S.l / totw
          for (l in 1:L) covs.EE[, , l] <- S.l

          for (l in 1:L) {
            covc <- covs.new[,,l]
            diag(covc) <- pmax(diag(covc), jitter)
            eigc <- eigen(covc, symmetric = TRUE)
            if (min(eigc$values) < jitter) {
              covc <- eigc$vectors %*% diag(pmax(eigc$values, jitter)) %*% t(eigc$vectors)
            }
            covs.VV[,,l] <- covc
          }

          if (constraint == "VE"){
            for (l in 1:L) {
              covs.new[, , l] <- covs.EE[, , l]
              diag(covs.new[, , l]) <- diag(covs.VV[, , l])
            }
          }else{
            for (l in 1:L) {
              covs.new[, , l] <- covs.VV[, , l]
              diag(covs.new[, , l]) <- diag(covs.EE[, , l])
            }
          }

        } else{
          for (l in 1:L) {
            covc <- covs.new[,,l]
            diag(covc) <- pmax(diag(covc), jitter)
            eigc <- eigen(covc, symmetric = TRUE)
            if (min(eigc$values) < jitter) {
              covc <- eigc$vectors %*% diag(pmax(eigc$values, jitter)) %*% t(eigc$vectors)
            }
            covs.new[,,l] <- covc
          }
        }
      }else{
        S.l <- array(0, dim = c(I, I))
        covs.EE <- array(0, dim = c(I, I, L))
        totw = 0
        for (l in 1:L) {
          wc <- P.Z.Xn[, l]
          totw <- totw + sum(wc)

          mean.l <- means.cur[l, ]
          dev <- sweep(response, 2, mean.l, "-")
          WDev <- dev * wc
          S.l <- S.l + crossprod(dev, WDev)
        }
        S.l <- S.l / totw
        for (l in 1:L) covs.EE[, , l] <- S.l

        for (l in 1:L) {
          covc <- covs.new[,,l]
          diag(covc) <- pmax(diag(covc), jitter)
          eigc <- eigen(covc, symmetric = TRUE)
          if (min(eigc$values) < jitter) {
            covc <- eigc$vectors %*% diag(pmax(eigc$values, jitter)) %*% t(eigc$vectors)
          }
          covs.new[,,l] <- covc
        }

        for(v in constraint){
          covs.new[v[1], v[2], ] <- covs.new[v[2], v[1], ] <- covs.EE[v[1], v[2], ]
        }
      }

      valid_update <- TRUE
      for (l in 1:L) {
        try_chol <- tryCatch(chol(covs.new[,,l]), error = function(e) NULL)
        if (is.null(try_chol)) {
          valid_update <- FALSE
          break
        }
      }

      if (!valid_update) {
        shared_mask <- matrix(FALSE, I, I)

        if (is.list(constraint)) {
          for (v in constraint) {
            i <- v[1]; j <- v[2]
            shared_mask[i, j] <- TRUE
            shared_mask[j, i] <- TRUE
          }
        } else if (constraint %in% c("E0", "EE")) {
          shared_mask[,] <- TRUE
        } else if (constraint == "VE") {
          shared_mask[lower.tri(shared_mask, diag = FALSE)] <- TRUE
          shared_mask[upper.tri(shared_mask, diag = FALSE)] <- TRUE
        } else if (constraint == "EV") {
          diag(shared_mask) <- TRUE
        } else if (constraint == "V0") {
          shared_mask[!diag(I)] <- TRUE
        }
        diag_shared <- diag(shared_mask)
        fixed_classes <- logical(L)
        max_fix_attempts <- 100

        for (attempt in 1:max_fix_attempts) {
          for (l in 1:L) {
            if (fixed_classes[l]) next

            cov_mat <- covs.new[, , l]
            cov_mat <- (cov_mat + t(cov_mat)) / 2
            eig <- eigen(cov_mat, symmetric = TRUE, only.values = TRUE)
            min_eig <- min(eig$values)

            if (min_eig > jitter) {
              fixed_classes[l] <- TRUE
              covs.new[, , l] <- cov_mat
              next
            }

            perturbation <- matrix(0, I, I)
            nonshared_vars <- which(!diag_shared)
            if (length(nonshared_vars) > 0) {
              required_increase <- abs(min_eig) + jitter
              var_increase <- rep(required_increase / length(nonshared_vars), length(nonshared_vars))
              diag(perturbation)[nonshared_vars] <- var_increase
            }

            if (min_eig <= jitter && length(nonshared_vars) == 0) {
              nonshared_offdiag <- which(!shared_mask & !diag(I), arr.ind = TRUE)
              if (nrow(nonshared_offdiag) > 0) {
                upper_idx <- nonshared_offdiag[nonshared_offdiag[,1] < nonshared_offdiag[,2], , drop = FALSE]
                if (nrow(upper_idx) > 0) {
                  perturbation[upper_idx] <- runif(nrow(upper_idx), -1e-4, 1e-4)
                  perturbation[upper_idx[,2], upper_idx[,1]] <- perturbation[upper_idx]
                }
              }
            }

            cov_mat <- cov_mat + perturbation
            cov_mat <- (cov_mat + t(cov_mat)) / 2

            eig <- eigen(cov_mat, symmetric = TRUE, only.values = TRUE)
            if (min(eig$values) > jitter) {
              fixed_classes[l] <- TRUE
              covs.new[, , l] <- cov_mat
            }
          }

          if (all(fixed_classes)) break

          if (attempt == max_fix_attempts) {
            for (l in which(!fixed_classes)) {
              cov_mat <- covs.new[, , l]
              cov_mat <- (cov_mat + t(cov_mat)) / 2
              eig <- eigen(cov_mat, symmetric = TRUE)
              new_vals <- pmax(eig$values, jitter)
              cov_new <- eig$vectors %*% diag(new_vals) %*% t(eig$vectors)
              cov_new <- (cov_new + t(cov_new)) / 2

              for (i in 1:I) {
                for (j in 1:I) {
                  if (shared_mask[i, j]) {
                    cov_new[i, j] <- cov_mat[i, j]
                  }
                }
              }
              covs.new[, , l] <- cov_new
            }
          }
        }

        valid_update <- TRUE
        for (l in 1:L) {
          try_chol <- tryCatch(chol(covs.new[, , l]), error = function(e) NULL)
          if (is.null(try_chol)) {
            valid_update <- FALSE
            break
          }
        }

        if (!valid_update) {
          covs.new <- prev_covs
          means.cur <- prev_means
          P.Z.cur <- prev_P.Z
        } else {
          covs.cur <- covs.new
          prev_means <- means.cur
          prev_covs <- covs.new
          prev_P.Z <- P.Z.cur
        }
      } else {
        covs.cur <- covs.new
        prev_means <- means.cur
        prev_covs <- covs.new
        prev_P.Z <- P.Z.cur
      }

      mx <- rowmax
      ll <- sum(mx + log(rowSums(exp(logres - mx))))
      if (!is.finite(ll)) {
        ll <- Log.Lik - 1e5
      }

      Log.Lik.history[iter] <- ll
      maxchg <- abs(ll - Log.Lik)
      AIC <- -2 * ll + 2 * npar
      BIC <- -2 * ll + npar * log(N)

      if (BIC < best_BIC){
        best_BIC <- BIC
      }

      if (vis) {
        if(vis && iter > 1 && r){
          if(wa){
            cat('\rWarm', paste0(sprintf("%2d", r), "/", sprintf("%2d", nrep.once)), '| Iter =', sprintf("%4d", iter),
                '  \u0394Log.Lik =', sprintf(fmt_string_maxchg, maxchg),
                '  BIC =', sprintf(fmt_string_BIC, BIC), '  Best_BIC =', sprintf(fmt_string_BIC, best_BIC))
          }else{
            cat('\rRep ', paste0(sprintf("%2d", r), "/", sprintf("%2d", nrep.once)), '| Iter =', sprintf("%4d", iter),
                '  \u0394Log.Lik =', sprintf(fmt_string_maxchg, maxchg),
                '  BIC =', sprintf(fmt_string_BIC, BIC), '  Best_BIC =', sprintf(fmt_string_BIC, best_BIC))
          }

        }else if(vis && iter > 1){
          cat('\rIter =', sprintf("%4d", iter), '  \u0394Log.Lik =', sprintf(fmt_string_maxchg, maxchg),
              '  BIC =', sprintf(fmt_string_BIC, BIC))
        }
      }

      if (maxchg < tol) {
        Log.Lik.history <- Log.Lik.history[1:iter]
        break
      }

      Log.Lik <- ll

      if (iter == maxiter.once && vis && r == 0) {
        message('\nMaximum number of iterations reached; convergence may not have been achieved\n')
      }
    }

    colnames(P.Z.Xn) <- paste0("Class ", 1:L)
    rownames(means.cur) <- paste0("Class ", 1:L)

    P.Z.Xn[is.na(P.Z.Xn)] <- 1/L
    P.Z.Xn[P.Z.Xn < 1e-20] <- 1e-20
    row_sums <- rowSums(P.Z.Xn)
    row_sums[row_sums < 1e-20] <- 1
    P.Z.Xn <- P.Z.Xn / row_sums

    Z <- apply(P.Z.Xn, 1, which.max)

    P.Z.cur <- as.table(P.Z.cur)
    names(P.Z.cur) <- paste0("Class ", 1:L)

    res <- list(
      params = list(means = means.cur, covs = covs.cur, P.Z = P.Z.cur),
      npar = npar,
      Log.Lik = Log.Lik,
      AIC = AIC,
      BIC = BIC,
      best_BIC=best_BIC,
      P.Z.Xn = P.Z.Xn,
      P.Z = P.Z.cur,
      Z = Z,
      Log.Lik.history = Log.Lik.history[1:iter]
    )

    return(res)
  }

  if(starts >= 1 && any(par.ini %in% c("random", "kmeans"))){
    for (s in 1:starts) {
      if (par.ini == "random") {
        means_init <- matrix(rnorm(L * I, mean = mean(response), sd = sd(as.vector(response))), nrow = L, ncol = I)
        P.Z_init <- rep(1/L, L)
        covs_init <- array(0, dim = c(I, I, L))
        for (l in 1:L) {
          # covs_init[,,l] <- covs.global
          covs_init[,,l] <- diag(I)
        }
        init_vals <- list(means = means_init, covs = covs_init, P.Z = P.Z_init)
      } else {
        kmeans_attempts <- 0
        while (kmeans_attempts < 5) {
          try_kmeans <- tryCatch({
            kmeans(response, centers = L, nstart = 1, iter.max = 500)
          }, error = function(e) NULL)
          if (!is.null(try_kmeans)) break
          kmeans_attempts <- kmeans_attempts + 1
        }

        if (is.null(try_kmeans)) {
          cluster_assignments <- sample(1:L, N, replace = TRUE)
          means_init <- matrix(apply(response, 2, mean), nrow = L, ncol = I, byrow = TRUE)
          P.Z_init <- rep(1/L, L)
        } else {
          kmeans.obj <- try_kmeans
          means_init <- matrix(kmeans.obj$centers, nrow = L, ncol = I)
          P.Z_init <- pmax(kmeans.obj$size / N, 1e-12)
          P.Z_init <- P.Z_init / sum(P.Z_init)
          cluster_assignments <- kmeans.obj$cluster
        }

        covs_init <- array(0, dim = c(I, I, L))
        for (l in 1:L) {
          idx <- which(cluster_assignments == l)
          if (length(idx) > 1) {
            covc <- cov(response[idx,,drop=FALSE])
            covc <- (covc + t(covc))/2
            diag(covc) <- pmax(diag(covc), jitter)
            covs_init[,,l] <- covc
          } else {
            covs_init[,,l] <- covs.global
          }
        }
        init_vals <- list(means = means_init, covs = covs_init, P.Z = P.Z_init)
      }

      res <- run_em_once(init_vals, s, best_BIC, wa=TRUE)
      best_BIC <- res$best_BIC

      if(length(results.wa$BIC) < nrep){
        results.wa$params[[length(results.wa$BIC)+1]] <- res$params
        results.wa$BIC <- c(results.wa$BIC, res$BIC)
      }else{
        BIC.wa.posi <- which.max(results.wa$BIC)
        BIC.wa <- max(results.wa$BIC)
        if(BIC.wa > res$BIC){
          results.wa$params[[BIC.wa.posi]] <- res$params
          results.wa$BIC[[BIC.wa.posi]] <- res$BIC
        }
      }
    }
    if(vis){
      cat("\n")
    }
  }

  if (is.null(results.wa$BIC) && nrep >= 1 && any(par.ini %in% c("random", "kmeans"))) {
    results <- NULL
    Log.Lik.nrep <- numeric(nrep)
    best_BIC <- Inf

    for (r in 1:nrep) {
      if (par.ini == "random") {
        means_init <- matrix(rnorm(L * I, mean = mean(response), sd = sd(as.vector(response))), nrow = L, ncol = I)
        P.Z_init <- rep(1/L, L)
        covs_init <- array(0, dim = c(I, I, L))
        for (l in 1:L) {
          covs_init[,,l] <- diag(I)
        }
        init_vals <- list(means = means_init, covs = covs_init, P.Z = P.Z_init)
      } else {
        kmeans_attempts <- 0
        while (kmeans_attempts < 5) {
          try_kmeans <- tryCatch({
            kmeans(response, centers = L, nstart = 1, iter.max = 500)
          }, error = function(e) NULL)
          if (!is.null(try_kmeans)) break
          kmeans_attempts <- kmeans_attempts + 1
        }

        if (is.null(try_kmeans)) {
          cluster_assignments <- sample(1:L, N, replace = TRUE)
          means_init <- matrix(apply(response, 2, mean), nrow = L, ncol = I, byrow = TRUE)
          P.Z_init <- rep(1/L, L)
        } else {
          kmeans.obj <- try_kmeans
          means_init <- matrix(kmeans.obj$centers, nrow = L, ncol = I)
          P.Z_init <- pmax(kmeans.obj$size / N, 1e-12)
          P.Z_init <- P.Z_init / sum(P.Z_init)
          cluster_assignments <- kmeans.obj$cluster
        }

        covs_init <- array(0, dim = c(I, I, L))
        for (l in 1:L) {
          idx <- which(cluster_assignments == l)
          if (length(idx) > 1) {
            covc <- cov(response[idx,,drop=FALSE])
            covc <- (covc + t(covc))/2
            diag(covc) <- pmax(diag(covc), jitter)
            covs_init[,,l] <- covc
          } else {
            covs_init[,,l] <- covs.global
          }
        }
        init_vals <- list(means = means_init, covs = covs_init, P.Z = P.Z_init)
      }

      res <- run_em_once(init_vals, r, best_BIC, wa=FALSE)
      results[[r]] <- res
      Log.Lik.nrep[r] <- res$Log.Lik

      best_BIC <- res$best_BIC
    }

    best_idx <- which.max(Log.Lik.nrep)
    res <- results[[best_idx]]
    res$Log.Lik.nrep <- Log.Lik.nrep

  }else if(!is.null(results.wa$BIC) && nrep >= 1 && any(par.ini %in% c("random", "kmeans"))){
    results <- NULL
    Log.Lik.nrep <- numeric(nrep)
    best_BIC <- Inf

    for(r in 1:nrep){
      res <- run_em_once(results.wa$params[[r]], r, best_BIC, wa=FALSE)
      results[[r]] <- res
      Log.Lik.nrep[r] <- res$Log.Lik

      best_BIC <- res$best_BIC
    }
    best_idx <- which.max(Log.Lik.nrep)
    res <- results[[best_idx]]
    res$Log.Lik.nrep <- Log.Lik.nrep

  } else {
    init_vals <- list(means = par.ini$means, covs = par.ini$covs, P.Z = par.ini$P.Z)
    res <- run_em_once(init_vals, 0, Inf, wa=FALSE)
    res$Log.Lik.nrep <- res$Log.Lik
  }

  if (vis) cat("\n\n")
  return(res)
}
