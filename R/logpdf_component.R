logpdf_component <- function(mean.l, covs.l, tresponse_local, jitter_local = 1e-10) {
  N <- ncol(tresponse_local)
  I <- nrow(tresponse_local)
  min_jitter <- max(jitter_local, sqrt(.Machine$double.eps))

  if (I == 1) {
    var_val <- if (is.matrix(covs.l)) covs.l[1, 1] else as.numeric(covs.l)

    if (!is.finite(var_val) || var_val <= 0) {
      var_val <- min_jitter
    } else {
      var_val <- max(var_val, min_jitter)
    }

    dev <- as.vector(tresponse_local - mean.l)
    log_var <- log(var_val)
    quad <- (dev^2) / var_val

    const_term <- -0.5 * log(2 * pi)
    lp <- const_term - 0.5 * (log_var + quad)

    lp[!is.finite(lp)] <- -1e300
    return(lp)
  }

  covs.l <- (covs.l + t(covs.l)) * 0.5

  R <- tryCatch(chol(covs.l), error = function(e) NULL)

  if (is.null(R)) {
    diag_vals <- diag(covs.l)
    off_sum <- rowSums(abs(covs.l)) - abs(diag_vals)
    gersh_bound <- min(diag_vals - off_sum)

    jitter_needed <- max(min_jitter, -gersh_bound + min_jitter)
    cov_perturb <- covs.l
    diag(cov_perturb) <- diag_vals + jitter_needed

    R <- tryCatch(chol(cov_perturb), error = function(e) NULL)

    if (is.null(R)) {
      eig <- eigen(covs.l, symmetric = TRUE)
      eig_vals <- pmax(eig$values, min_jitter)
      covs.l <- eig$vectors %*% (eig_vals * t(eig$vectors))
      covs.l <- (covs.l + t(covs.l)) * 0.5
      R <- chol(covs.l)
    } else {
      covs.l <- cov_perturb
    }
  }

  dev_mat <- tresponse_local - matrix(mean.l, nrow = I, ncol = N, byrow = FALSE)

  tryCatch({
    sol <- forwardsolve(t(R), dev_mat)
    quad <- colSums(sol^2)
  }, error = function(e) {
    quad <- rep(Inf, N)
  })

  diag_R <- diag(R)
  if (any(diag_R <= 0)) {
    diag_R[diag_R <= 0] <- min_jitter
    logdet <- 2 * sum(log(diag_R))
  } else {
    logdet <- 2 * sum(log(diag_R))
  }

  const_term <- -0.5 * I * log(2 * pi) - 0.5 * logdet
  lp <- const_term - 0.5 * quad

  if (any(!is.finite(lp))) {
    bad_idx <- !is.finite(lp)
    lp[bad_idx] <- -1e300
    warning(sum(bad_idx), " non-finite log-density values replaced with -1e300")
  }

  return(lp)
}
