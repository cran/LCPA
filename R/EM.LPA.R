
#' @importFrom stats cov sd
EM.LPA <- function(response, L = 2, par.ini = "random", constraint = "VV",
                   starts=50, maxiter.warmup=100, nrep = 1, vis = TRUE,
                   maxiter = 2000, tol = 1e-4) {

  if (is.vector(response) && is.numeric(response)) response <- matrix(response, ncol = 1)
  if (is.data.frame(response)) response <- as.matrix(response)
  if (!is.matrix(response) || !is.numeric(response)) stop("response must be a numeric matrix/data.frame/vector")
  if (any(!is.finite(response))) stop("response must not contain missing or non-finite values")

  N <- nrow(response)
  I <- ncol(response)
  if (N == 0) stop("empty response")
  if (length(L) != 1 || L < 1 || L != as.integer(L) || L > N) stop("L must be an integer between 1 and nrow(response)")
  .validate.training.stages(starts, maxiter.warmup, nrep)
  if (maxiter < 1 || tol <= 0) stop("maxiter and tol must be positive")
  jitter = 1e-10

  constraint <- .validate.LPA.constraint(constraint, I)

  if (is.null(par.ini)) par.ini <- "random"
  if (is.character(par.ini) && (length(par.ini) != 1 || !par.ini %in% c("random", "kmeans"))) {
    stop("par.ini must be 'random', 'kmeans', or a parameter list")
  }
  if (!is.character(par.ini) &&
      (!is.list(par.ini) || !all(c("means", "covs", "P.Z") %in% names(par.ini)))) {
    stop("par.ini must contain means, covs, and P.Z")
  }

  covs.global <- matrix(cov(response), nrow = I, ncol = I)
  if (any(!is.finite(covs.global))) {
    vars <- apply(response, 2, var)
    vars[!is.finite(vars) | vars <= 0] <- 1
    covs.global <- diag(vars, I)
  }
  covs.global <- (covs.global + t(covs.global)) / 2
  repair_covariance <- function(covariance, fallback = covs.global){
    .stabilize.LPA.covariances(
      array(covariance, dim = c(I, I, 1L)), "VV", fallback = fallback
    )$covs[, , 1L]
  }
  covs.global <- repair_covariance(covs.global, diag(I))

  npar <- get.npar.LPA(I, L, constraint)
  expectation.step <- function(means, covs, P.Z) {
    if (any(!is.finite(P.Z)) || any(P.Z <= 0)) return(NULL)
    P.Z <- P.Z / sum(P.Z)
    result <- lpa_expectation_cpp(response, means, covs, P.Z, repair = FALSE)
    if(!isTRUE(result$valid) || !is.finite(result$Log.Lik)) return(NULL)
    P.Z.Xn <- result$posterior
    empty <- !is.finite(colSums(P.Z.Xn)) | colSums(P.Z.Xn) <= jitter
    if (any(empty)) {
      P.Z.Xn[, empty] <- pmax(P.Z.Xn[, empty, drop = FALSE], min(1e-8, 0.01 / L))
      P.Z.Xn <- P.Z.Xn / rowSums(P.Z.Xn)
    }
    list(P.Z.Xn = P.Z.Xn, Log.Lik = result$Log.Lik)
  }

  covariance_is_pd <- function(covs) {
    all(vapply(1:L, function(l) {
      covs.l <- covs[, , l]
      all(is.finite(covs.l)) &&
        max(abs(covs.l - t(covs.l))) < 1e-8 &&
        !is.null(tryCatch(chol(covs.l), error = function(e) NULL))
    }, logical(1)))
  }

  update_constrained_covariances <- function(scatter, nk, covs.start) {
    pairs <- which(lower.tri(matrix(TRUE, I, I), diag = TRUE), arr.ind = TRUE)
    pair_keys <- apply(pairs, 1, function(x) paste(sort(x), collapse = ":"))

    if (is.character(constraint)) {
      shared <- if (constraint == "VE") pairs[, 1] != pairs[, 2] else pairs[, 1] == pairs[, 2]
    } else {
      constraint_keys <- unique(vapply(constraint, function(x) {
        paste(sort(as.integer(x)), collapse = ":")
      }, character(1)))
      shared <- pair_keys %in% constraint_keys
    }

    param_index <- matrix(0L, nrow = nrow(pairs), ncol = L)
    ntheta <- 0L
    for (p in 1:nrow(pairs)) {
      if (shared[p]) {
        ntheta <- ntheta + 1L
        param_index[p, ] <- ntheta
      } else {
        param_index[p, ] <- ntheta + seq_len(L)
        ntheta <- ntheta + L
      }
    }

    make_covs <- function(theta) {
      covs <- array(0, dim = c(I, I, L))
      for (p in 1:nrow(pairs)) {
        i <- pairs[p, 1]
        j <- pairs[p, 2]
        for (l in 1:L) covs[i, j, l] <- covs[j, i, l] <- theta[param_index[p, l]]
      }
      covs
    }

    make_theta <- function(covs) {
      theta <- numeric(ntheta)
      for (p in 1:nrow(pairs)) {
        i <- pairs[p, 1]
        j <- pairs[p, 2]
        if (shared[p]) {
          theta[param_index[p, 1]] <- sum(nk * covs[i, j, ]) / sum(nk)
        } else {
          theta[param_index[p, ]] <- covs[i, j, ]
        }
      }
      theta
    }

    objective <- function(theta) {
      covs <- make_covs(theta)
      value <- 0
      for (l in 1:L) {
        R <- tryCatch(chol(covs[, , l]), error = function(e) NULL)
        if (is.null(R)) return(.Machine$double.xmax^0.25)
        inv <- chol2inv(R)
        value <- value + nk[l] * 2 * sum(log(diag(R))) + sum(scatter[, , l] * inv)
      }
      value / 2
    }

    gradient <- function(theta) {
      covs <- make_covs(theta)
      grad <- numeric(ntheta)
      for (l in 1:L) {
        R <- tryCatch(chol(covs[, , l]), error = function(e) NULL)
        if (is.null(R)) return(rep(0, ntheta))
        inv <- chol2inv(R)
        G <- (nk[l] * inv - inv %*% scatter[, , l] %*% inv) / 2
        for (p in 1:nrow(pairs)) {
          i <- pairs[p, 1]
          j <- pairs[p, 2]
          value <- if (i == j) G[i, j] else 2 * G[i, j]
          idx <- param_index[p, l]
          grad[idx] <- grad[idx] + value
        }
      }
      grad
    }

    unconstrained <- array(0, dim = c(I, I, L))
    safe <- array(0, dim = c(I, I, L))
    for (l in 1:L) unconstrained[, , l] <- scatter[, , l] / nk[l]

    for (i in 1:I) {
      p <- which(pairs[, 1] == i & pairs[, 2] == i)
      values <- unconstrained[i, i, ]
      if (shared[p]) values[] <- sum(scatter[i, i, ]) / sum(nk)
      if (any(!is.finite(values))) values <- covs.start[i, i, ]
      values[values <= 0] <- covs.start[i, i, ][values <= 0]
      if (shared[p]) values[] <- sum(nk * values) / sum(nk)
      safe[i, i, ] <- values
    }

    theta_candidates <- list(make_theta(covs.start), make_theta(unconstrained), make_theta(safe))
    valid_candidates <- vapply(theta_candidates, function(theta) covariance_is_pd(make_covs(theta)), logical(1))
    if (!any(valid_candidates)) return(NULL)
    theta_candidates <- theta_candidates[valid_candidates]
    objective_values <- vapply(theta_candidates, objective, numeric(1))
    theta.start <- theta_candidates[[which.min(objective_values)]]
    objective.start <- min(objective_values)

    fit <- tryCatch(
      stats::optim(
        theta.start,
        objective,
        gradient,
        method = "BFGS",
        control = list(maxit = 100, reltol = min(tol / 10, 1e-8))
      ),
      error = function(e) NULL
    )

    if (is.null(fit)) return(make_covs(theta.start))
    covs.fit <- make_covs(fit$par)
    if (!covariance_is_pd(covs.fit) ||
        objective(fit$par) > objective.start + 1e-8 * (1 + abs(objective.start))) {
      return(make_covs(theta.start))
    }
    covs.fit
  }

  update_covariances <- function(scatter, nk, covs.start) {
    pooled <- apply(scatter, c(1, 2), sum) / sum(nk)

    if (I == 1) {
      shared <- if (is.character(constraint)) {
        constraint %in% c("UE", "E0", "EE", "EV")
      } else {
        any(vapply(constraint, function(x) all(as.integer(x) == 1L), logical(1)))
      }
      values <- if (shared) rep(pooled[1, 1], L) else scatter[1, 1, ] / nk
      if (any(!is.finite(values))) values <- rep(covs.global[1, 1], L)
      covs <- array(values, dim = c(1, 1, L))
      for (l in 1:L) covs[, , l] <- repair_covariance(covs[, , l])
      return(covs)
    }

    covs <- array(0, dim = c(I, I, L))
    if (is.character(constraint) && constraint == "E0") {
      for (l in 1:L) covs[, , l] <- repair_covariance(diag(diag(pooled), I))
    } else if (is.character(constraint) && constraint == "V0") {
      for (l in 1:L){
        covariance <- diag(diag(scatter[, , l]) / nk[l], I)
        covs[, , l] <- repair_covariance(covariance)
      }
    } else if (is.character(constraint) && constraint == "EE") {
      pooled <- repair_covariance(pooled)
      for (l in 1:L) covs[, , l] <- pooled
    } else if (is.character(constraint) && constraint == "VV") {
      for (l in 1:L) covs[, , l] <- repair_covariance(scatter[, , l] / nk[l])
    } else {
      return(update_constrained_covariances(scatter, nk, covs.start))
    }

    if (!covariance_is_pd(covs)) return(NULL)
    covs
  }

  run.EM <- function(par.ini.current, r = 0, best_BIC = Inf, warmup=FALSE) {
    iteration.progress.state <- .new.progress.state()
    means.cur <- as.matrix(par.ini.current$means)
    covs.cur <- par.ini.current$covs
    P.Z.cur <- as.numeric(par.ini.current$P.Z)

    if (length(dim(means.cur)) != 2 || any(dim(means.cur) != c(L, I))) stop("initial means must be an L x I matrix")
    if (any(!is.finite(means.cur))) stop("initial means must be finite")
    if (length(dim(covs.cur)) != 3 || any(dim(covs.cur) != c(I, I, L))) stop("initial covs must be an I x I x L array")
    if (length(P.Z.cur) != L || any(!is.finite(P.Z.cur)) || any(P.Z.cur <= 0)) {
      stop("initial P.Z must contain L positive probabilities")
    }
    P.Z.cur <- P.Z.cur / sum(P.Z.cur)
    for (l in 1:L) covs.cur[, , l] <- repair_covariance(covs.cur[, , l])

    maxiter.current <- if (warmup) maxiter.warmup else maxiter
    Log.Lik.history <- numeric(maxiter.current)

    estep.cur <- expectation.step(means.cur, covs.cur, P.Z.cur)
    if (is.null(estep.cur)) {
      return(list(
        params = list(means = means.cur, covs = covs.cur, P.Z = P.Z.cur),
        npar = npar, Log.Lik = -Inf, AIC = Inf, BIC = Inf, best_BIC = best_BIC,
        P.Z.Xn = matrix(1 / L, N, L), P.Z = P.Z.cur,
        Z = rep(1L, N), Log.Lik.history = -Inf
      ))
    }

    Log.Lik <- estep.cur$Log.Lik
    P.Z.Xn <- estep.cur$P.Z.Xn
    AIC <- -2 * Log.Lik + 2 * npar
    BIC <- -2 * Log.Lik + npar * log(N)
    iter <- 0L

    for (iter in 1:maxiter.current) {
      nk <- colSums(P.Z.Xn)
      if (any(!is.finite(nk)) || any(nk <= jitter)) {
        break
      }

      P.Z.new <- nk / N
      means.new <- sweep(t(P.Z.Xn) %*% response, 1, nk, "/")

      scatter <- array(0, dim = c(I, I, L))
      for (l in 1:L) {
        dev <- sweep(response, 2, means.new[l, ], "-")
        scatter[, , l] <- crossprod(dev, dev * P.Z.Xn[, l])
      }

      covs.new <- update_covariances(scatter, nk, covs.cur)
      if (is.null(covs.new)) {
        break
      }

      estep.new <- expectation.step(means.new, covs.new, P.Z.new)
      if (is.null(estep.new)) {
        break
      }

      ll <- estep.new$Log.Lik
      maxchg <- abs(ll - Log.Lik)
      means.cur <- means.new
      covs.cur <- covs.new
      P.Z.cur <- P.Z.new
      P.Z.Xn <- estep.new$P.Z.Xn
      Log.Lik <- ll
      Log.Lik.history[iter] <- ll
      AIC <- -2 * ll + 2 * npar
      BIC <- -2 * ll + npar * log(N)

      if (BIC < best_BIC) best_BIC <- BIC

      if (vis && iter > 1 && r == 0) {
        .print.iteration.progress(
          iter, maxchg, BIC,
          progress.prefix = .estimation.output.prefix(),
          progress.state = iteration.progress.state
        )
      }

      if (maxchg < tol) break
      if (iter == maxiter.current && vis && r == 0) {
        message(
          '\n', .estimation.output.prefix(),
          'Maximum number of iterations reached; convergence may not have been achieved\n'
        )
      }
    }

    colnames(P.Z.Xn) <- .latent.group.names(L, "LPA")
    rownames(means.cur) <- .latent.group.names(L, "LPA")
    Z <- max.col(P.Z.Xn, ties.method = "first")

    P.Z.cur <- as.table(P.Z.cur)
    names(P.Z.cur) <- .latent.group.names(L, "LPA")

    list(
      params = list(means = means.cur, covs = covs.cur, P.Z = P.Z.cur),
      npar = npar,
      Log.Lik = Log.Lik,
      AIC = AIC,
      BIC = BIC,
      best_BIC=best_BIC,
      P.Z.Xn = P.Z.Xn,
      P.Z = P.Z.cur,
      Z = Z,
      Log.Lik.history = Log.Lik.history[seq_len(iter)]
    )
  }

  make_initial_values <- function() {
    if (!is.character(par.ini)) {
      return(par.ini)
    }
    if (par.ini == "random") {
      sd.response <- sd(as.vector(response))
      if (!is.finite(sd.response) || sd.response <= 0) sd.response <- 1
      means_init <- matrix(rnorm(L * I, mean = mean(response), sd = sd.response), nrow = L, ncol = I)
      P.Z_init <- rep(1 / L, L)
      covs_init <- array(0, dim = c(I, I, L))
      for (l in 1:L) covs_init[, , l] <- covs.global
      return(list(means = means_init, covs = covs_init, P.Z = P.Z_init))
    }

    Kmeans.LPA(response, L, constraint = constraint, starts = 1)$params
  }

  best_BIC <- Inf
  warmup.params <- vector("list", starts)
  warmup.Log.Lik <- rep(-Inf, starts)
  warmup.progress.state <- .new.progress.state()
  for (s in seq_len(starts)) {
    res.warmup <- run.EM(make_initial_values(), s, best_BIC, warmup=TRUE)
    best_BIC <- res.warmup$best_BIC
    warmup.params[[s]] <- res.warmup$params
    warmup.Log.Lik[s] <- res.warmup$Log.Lik
    if(vis){
      .print.estimation.progress(
        "Warm", s, starts, res.warmup$Log.Lik,
        algorithm = "EM",
        iterations = length(res.warmup$Log.Lik.history),
        progress.state = warmup.progress.state
      )
    }
  }
  if(vis) .end.estimation.progress()

  warmup.position <- order(warmup.Log.Lik, decreasing = TRUE)
  warmup.position <- warmup.position[is.finite(warmup.Log.Lik[warmup.position])]
  if(length(warmup.position) < nrep){
    stop(
      "EM warm-up produced only ", length(warmup.position),
      " finite LPA solutions from ", starts, " starts; nrep = ", nrep, "."
    )
  }
  warmup.position <- warmup.position[seq_len(nrep)]

  results <- vector("list", nrep)
  Log.Lik.nrep <- rep(-Inf, nrep)
  best.Log.Lik <- -Inf
  replication.progress.state <- .new.progress.state()
  for (r in seq_len(nrep)) {
    par.ini.current <- warmup.params[[warmup.position[r]]]
    results[[r]] <- run.EM(par.ini.current, r, best_BIC, warmup=FALSE)
    Log.Lik.nrep[r] <- results[[r]]$Log.Lik
    best_BIC <- results[[r]]$best_BIC
    best.Log.Lik <- max(best.Log.Lik, Log.Lik.nrep[r])
    if(vis){
      .print.estimation.progress(
        "Rep", r, nrep, Log.Lik.nrep[r], best.Log.Lik,
        progress.state = replication.progress.state
      )
    }
  }
  if(vis) .end.estimation.progress()

  if (!any(is.finite(Log.Lik.nrep))) {
    stop("all EM refinements failed because a class became empty or a covariance matrix became singular")
  }
  best_idx <- which.max(Log.Lik.nrep)
  res <- results[[best_idx]]
  res$Log.Lik.nrep <- Log.Lik.nrep

  if (vis) .print.estimation.summary("EM", res$Log.Lik, res$BIC)
  return(res)
}
