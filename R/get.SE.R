#' Compute Standard Errors
#'
#' Computes approximate standard errors (SEs) for estimated parameters in Latent Class Analysis (LCA)
#' or Latent Profile Analysis (LPA) models using two methods:
#' \itemize{
#'   \item \code{"Bootstrap"}: Non-parametric bootstrap with label-switching correction.
#'                             McLachlan & Peel (2004) suggest that 50–100 replicates often provide adequate accuracy
#'                             for practical purposes, though more (e.g., 500–1000) may be preferred for publication-quality inference.
#'   \item \code{"Obs"}: Numerical evaluation of the observed information matrix (Hessian of negative log-likelihood)
#' }
#' Users should note that \code{\link[LCPA]{get.SE}} computes standard errors based on the observed
#' information matrix via numerical differentiation, which may lack precision and often yields
#' ill-conditioned matrices. Therefore, we recommend using \code{method = "Bootstrap"}.
#'
#' @param object An object of class \code{"LCA"} or \code{"LPA"} returned by \code{\link{LCA}} or \code{\link{LPA}}.
#' @param method Character specifying SE calculation method: \code{"Obs"} or \code{"Bootstrap"} (default).
#' @param n.Bootstrap Integer. Number of bootstrap replicates when \code{method="Bootstrap"} (default=100).
#' @param vis Logical. If \code{TRUE}, displays progress information during estimation (default: \code{TRUE}).
#'
#' @return A list of class \code{"SE"} containing:
#'   \describe{
#'     \item{\code{se}}{Named list of SEs matching parameter structure of input model:
#'       \itemize{
#'         \item LPA: \code{means} (matrix: classes x variables),
#'                            \code{covs} (array: vars x vars x classes),
#'                            \code{P.Z} (vector: class prob SEs)
#'         \item LCA: \code{par} (array: classes x indicators x categories),
#'                            \code{P.Z} (vector: class prob SEs)
#'         \item Critical Note for \code{"Obs"} method:
#'               Only \emph{free parameters} have non-zero SEs. Non-free parameters (e.g., last class probability in \code{P.Z} due
#'               to sum-to-1 constraint; last category probability in LCA indicators) have SE=0.
#'               Bootstrap provides SEs for all parameters.
#'       }
#'     }
#'     \item{\code{vcov}}{\code{NULL} for bootstrap. For \code{"Obs"}: variance-covariance matrix (may be regularized). Diagonal
#'                       contains squared SEs of free parameters.}
#'     \item{\code{hessian}}{\code{NULL} for bootstrap. For \code{"Obs"}: observed information matrix (pre-regularization).
#'                           Dimension = number of free parameters.}
#'     \item{\code{diagnostics}}{Method-specific diagnostics:
#'       \itemize{
#'         \item Bootstrap: \code{n.Bootstrap.requested}, \code{n.Bootstrap.completed}
#'         \item Obs: Hessian computation details, condition number, regularization status, step sizes
#'       }
#'     }
#'     \item{\code{call}}{Function call that generated the object}
#'     \item{\code{arguments}}{List of input arguments}
#'   }
#'
#'
#' @references
#' McLachlan, G. J., & Peel, D. (2004). Finite Mixture Models. Wiley.
#' https://books.google.com.sg/books?id=c2_fAox0DQoC
#'
#' @examples
#' \donttest{
#' library(LCPA)
#' set.seed(123)
#'
#' # LPA with Bootstrap (minimal replicates for example)
#' lpa_data <- sim.LPA(N = 500, I = 4, L = 3)
#' lpa_fit <- LPA(lpa_data$response, L = 3)
#' se_boot <- get.SE(lpa_fit, method = "Bootstrap", n.Bootstrap = 10)
#'
#' print(se_boot)
#' extract(se_boot, "covs")
#'
#'
#' # LCA with Observed Information (note zeros for constrained parameters)
#' lca_data <- sim.LCA(N = 500, I = 4, L = 3, poly.value = 5)
#' lca_fit <- LCA(lca_data$response, L = 3)
#' se_obs <- get.SE(lca_fit, method = "Obs")
#'
#' print(se_obs)
#' extract(se_obs, "par")
#' }
#'
#' @importFrom numDeriv hessian
#' @importFrom Matrix nearPD
#' @importFrom clue solve_LSAP
#' @importFrom stats sd plogis qlogis
#' @export
get.SE <- function(object, method = "Bootstrap", n.Bootstrap = 100, vis = TRUE) {
  if (!inherits(object, c("LCA", "LPA"))) {
    stop("object must be of class 'LCA' or 'LPA'")
  }

  call <- match.call()

  model.type <- class(object)[1]
  response <- object$arguments$response
  params <- object$params
  N <- nrow(response)

  if (model.type == "LCA") {
    adjust.response.obj <- adjust.response(response)
    poly.value <- adjust.response.obj$poly.value
    poly.max <- adjust.response.obj$poly.max
  } else {
    poly.value <- NULL
    poly.max <- NULL
  }

  I <- ncol(response)
  L <- length(params$P.Z)

  if (method == "Bootstrap") {
    if(model.type == "LPA"){
      P.Z.Bootstrap <- array(0, dim = c(L, n.Bootstrap))
      means.Bootstrap <- array(0, dim = c(L, I, n.Bootstrap))
      covs.Bootstrap <- array(0, dim = c(I, I, L, n.Bootstrap))
      par.Bootstrap <- NULL
    }else{
      P.Z.Bootstrap <- array(0, dim = c(L, n.Bootstrap))
      par.Bootstrap <- array(0, dim = c(L, I, poly.max, n.Bootstrap))
      means.Bootstrap <- NULL
      covs.Bootstrap <- NULL

      par.t <- NULL
      for(po in 1:poly.max){
        par.t <- cbind(par.t, params$par[ , , po])
      }
      posi <- apply(par.t, 2, function(x) !any(is.na(x)))
      par.t <- t(par.t)[posi, ]
    }
    bs <- 0
    iter <- 0
    while (iter < n.Bootstrap){
      iter <- iter + 1
      if(vis){
        cat("\rRunning bootstrap (", paste0(iter, "/", n.Bootstrap), "replicates) ...")
      }
      indices <- sample(seq_len(N), size = N, replace = TRUE)
      response.cur <- response[indices, , drop = FALSE]
      object.cur <- tryCatch(
        update(object, response = response.cur, par.ini = "random", vis = FALSE),
        error = function(e) NULL
      )

      if (is.null(object.cur)) {
        next
      } else {
        params.cur <- object.cur$params
        if(model.type == "LPA"){
          means.cur <- params.cur$means
          covs.cur <- params.cur$covs
          P.Z.cur <- params.cur$P.Z

          dist.mat <- distance.matrix(t(params$means), t(means.cur))
          assignment <- clue::solve_LSAP(dist.mat)
          means.Bootstrap[, , iter]  <- means.cur[as.numeric(assignment), ]
          covs.Bootstrap[, , , iter]  <- covs.cur[, , as.numeric(assignment)]
          P.Z.Bootstrap[, iter]  <- P.Z.cur[as.numeric(assignment)]
        }else{
          par.cur <- params.cur$par
          P.Z.cur <- params.cur$P.Z

          par.e <- NULL
          for(po in 1:poly.max){
            par.e <- cbind(par.e, params.cur$par[ , , po])
          }
          par.e <- t(par.e)[posi, ]
          dist.mat <- distance.matrix(par.t, par.e)
          assignment <- clue::solve_LSAP(dist.mat)

          par.Bootstrap[, , , iter]  <- par.cur[as.numeric(assignment), , ]
          P.Z.Bootstrap[, iter]  <- P.Z.cur[as.numeric(assignment)]
        }
        bs <- bs + 1
      }
    }

    if(model.type == "LPA"){
      P.Z.se <- apply(P.Z.Bootstrap, 1, sd, na.rm = TRUE)
      means.se <- apply(means.Bootstrap, c(1, 2), sd, na.rm = TRUE)
      covs.se <- apply(covs.Bootstrap, c(1, 2, 3), sd, na.rm = TRUE)
    }else if(model.type == "LCA"){
      P.Z.se <- apply(P.Z.Bootstrap, 1, sd, na.rm = TRUE)
      par.se <- apply(par.Bootstrap, c(1, 2, 3), sd, na.rm = TRUE)
    }

    if (model.type == "LPA") {
      se.list <- list(means=means.se, covs=covs.se, P.Z=P.Z.se)
    } else if (model.type == "LCA") {
      se.list <- list(par=par.se, P.Z=P.Z.se)
    }

    diagnostics <- list(
      method = "Bootstrap",
      n.Bootstrap.requested = n.Bootstrap,
      n.Bootstrap.completed = bs
    )

  } else {
    if (model.type == "LPA") {
      logit_PZ <- if (L > 1) log(params$P.Z[-L] / params$P.Z[L]) else numeric(0)
      means_vec <- as.vector(t(params$means))
      covs_vec <- unlist(lapply(1:L, function(l) {
        cov_mat <- params$covs[, , l]
        cov_mat[upper.tri(cov_mat, diag = TRUE)]
      }))
      params.vec <- c(logit_PZ, means_vec, covs_vec)

      compute_step_sizes1 <- function(params.vec, L, I) {
        d_vec <- rep(1e-4, length(params.vec))
        total_params <- length(params.vec)

        if (L > 1) {
          pz_idx <- 1:(L-1)
          d_vec[pz_idx] <- pmax(abs(params.vec[pz_idx]) * 0.01, 1e-5)
        }

        means_start <- if (L > 1) L else 1
        means_end <- means_start + L * I - 1
        if (means_start <= total_params) {
          mean_vals <- params.vec[means_start:means_end]
          d_vec[means_start:means_end] <- pmax(abs(mean_vals) * 0.1, 0.01)  # 10% relative step
          d_vec[means_start:means_end] <- pmin(d_vec[means_start:means_end], 0.1)  # Cap at 0.1
        }

        covs_start <- means_end + 1
        if (covs_start <= total_params) {
          cov_vals <- params.vec[covs_start:total_params]
          d_vec[covs_start:total_params] <- pmax(abs(cov_vals) * 0.1, 0.001)  # 10% relative step
          d_vec[covs_start:total_params] <- pmin(d_vec[covs_start:total_params], 0.1)  # Cap at 0.1
        }
        d_vec
      }

      d_vec <- compute_step_sizes1(params.vec, L, I)

      func.nll.logit <- function(par, model.type, response, I, L, poly.value = NULL) {
        if (L > 1) {
          logit_PZ <- par[1:(L-1)]
          P.Z <- c(exp(logit_PZ), 1)
          P.Z <- P.Z / sum(P.Z)
        } else {
          P.Z <- rep(1, L)
        }

        means_mat <- matrix(par[(L):(L + L*I - 1)], nrow = L, byrow = TRUE)
        covs_start <- L + L*I
        covs_list <- vector("list", L)
        idx <- covs_start
        for (l in 1:L) {
          n_elems <- I * (I + 1) / 2
          vec <- par[idx:(idx + n_elems - 1)]
          mat <- matrix(0, I, I)
          mat[upper.tri(mat, diag = TRUE)] <- vec
          mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
          covs_list[[l]] <- mat
          idx <- idx + n_elems
        }
        covs_array <- array(unlist(covs_list), dim = c(I, I, L))

        -get.Log.Lik.LPA(response, P.Z, means_mat, covs_array)
      }

      hess_result <- tryCatch({
        numDeriv::hessian(
          func = function(x) func.nll.logit(x, model.type, response, I, L, poly.value),
          x = params.vec,
          method = "Richardson",
          method.args = list(eps = 1, d = d_vec, zero.tol = 1e-8, r = 8)
        )
      }, error = function(e) {
        matrix(0, ncol = length(params.vec), nrow = length(params.vec))
      })

      if (!is.matrix(hess_result) || nrow(hess_result) != length(params.vec)) {
        hess.matrix <- matrix(0, ncol = length(params.vec), nrow = length(params.vec))
      } else {
        hess.matrix <- hess_result
      }

      suppressWarnings(cov.mat <- tryCatch({
        hess.matrix <- (hess.matrix + t(hess.matrix)) / 2
        eig <- eigen(hess.matrix, symmetric = TRUE)
        vals <- eig$values
        vecs <- eig$vectors
        max_eig <- max(abs(vals))
        abs_threshold <- 1e-6
        rel_threshold <- 1e-4 * max_eig
        threshold <- max(abs_threshold, rel_threshold)
        n_neg <- sum(vals < 0)
        n_small <- sum(vals > 0 & vals < threshold)
        frac_trunc <- (n_neg + n_small) / length(vals)
        vals[vals < threshold] <- threshold
        hess.pos <- vecs %*% diag(vals) %*% t(vecs)

        if (min(eigen(hess.pos, symmetric = TRUE, only.values = TRUE)$values) < 1e-8) {
          hess.pos <- Matrix::nearPD(hess.pos,
                                     corr = FALSE,
                                     eig.tol = 1e-6,
                                     conv.tol = 1e-8,
                                     keepDiag = TRUE)$mat
          regularization_note <- "Eigenvalue truncation + nearPD fallback"
        } else {
          regularization_note <- "Eigenvalue truncation only"
        }
        chol_test <- try(chol(hess.pos), silent = TRUE)
        if (inherits(chol_test, "try-error")) {
          hess.pos <- Matrix::nearPD(hess.pos,
                                     corr = FALSE,
                                     eig.tol = 1e-6,
                                     conv.tol = 1e-8,
                                     keepDiag = TRUE)$mat
          regularization_note <- "Cholesky failed, nearPD used"
        }

        solve(hess.pos)
      }, error = function(e) {
        matrix(0, ncol = length(params.vec), nrow = length(params.vec))
      }))
      se.vec.logit <- sqrt(pmax(1e-6, diag(cov.mat)))

      se.vec <- numeric(length(se.vec.logit))
      if (L > 1) {
        PZ_free <- pmin(pmax(params$P.Z[-L], 1e-6), 1-1e-6)
        se.vec[1:(L-1)] <- PZ_free * (1 - PZ_free) * se.vec.logit[1:(L-1)]
      }

      if (length(se.vec.logit) > (L - 1)) {
        se.vec[L:length(se.vec)] <- se.vec.logit[L:length(se.vec.logit)]
      }

      num_pz <- L - 1
      num_means <- L * I
      num_covs <- L * I * (I + 1) / 2

      se.means <- matrix(se.vec[(num_pz + 1):(num_pz + num_means)],
                         nrow = L, byrow = TRUE)

      se.covs <- array(0, dim = c(I, I, L))
      start <- num_pz + num_means + 1
      for (l in 1:L) {
        n_elems <- I * (I + 1) / 2
        vec <- se.vec[start:(start + n_elems - 1)]
        mat <- matrix(0, I, I)
        mat[upper.tri(mat, diag = TRUE)] <- vec
        mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
        se.covs[, , l] <- mat
        start <- start + n_elems
      }

      se.pz <- numeric(L)
      if (L > 1) {
        se.pz[-L] <- se.vec[1:(L-1)]
        se.pz[L] <- 0
      } else {
        se.pz <- c(0)
      }

      se.list <- list(means = se.means, covs = se.covs, P.Z = se.pz)
      final_cond_num <- if (!exists("hess.pos")) NA else kappa(hess.pos)

      diagnostics <- list(
        method = "Obs",
        hessian_cond_number_initial = kappa(hess.matrix),
        hessian_cond_number_final = final_cond_num,
        eigenvalue_truncation_fraction = frac_trunc,
        regularization_method = regularization_note,
        step_sizes = d_vec,
        note = "Use Bootstrap method for more reliable SEs when condition number > 1e6"
      )

    } else {  # LCA Obs branch
      free_par_prob <- c()
      param_indices <- list()
      idx_counter <- 1
      for (i in 1:I) {
        K_i <- poly.value[i]
        if (K_i > 1) {
          for (l in 1:L) {
            free_probs <- params$par[l, i, 1:(K_i-1)]
            free_par_prob <- c(free_par_prob, free_probs)
            for (k in 1:(K_i-1)) {
              param_indices[[idx_counter]] <- c(l, i, k)
              idx_counter <- idx_counter + 1
            }
          }
        }
      }

      free_par_logit <- qlogis(pmax(pmin(free_par_prob, 1 - 1e-8), 1e-8))

      logit_PZ <- if (L > 1) log(params$P.Z[-L] / params$P.Z[L]) else numeric(0)
      params.vec <- c(logit_PZ, free_par_logit)

      compute_step_sizes2 <- function(params.vec, model.type, L, I, poly.value = NULL) {
        d_vec <- rep(1e-4, length(params.vec))
        if (L > 1) {
          d_vec[1:(L-1)] <- pmax(abs(params.vec[1:(L-1)]) * 0.01, 1e-5)
        }
        if (length(params.vec) > (L-1)) {
          d_vec[-(1:(L-1))] <- 0.1
        }
        d_vec
      }

      d_vec <- compute_step_sizes2(params.vec, model.type, L, I, poly.value)

      func.nll.logit <- function(par, model.type, response, I, L, poly.value = NULL) {
        if (L > 1) {
          logit_PZ <- par[1:(L-1)]
          P.Z <- c(exp(logit_PZ), 1)
          P.Z <- P.Z / sum(P.Z)
        } else {
          P.Z <- rep(1, L)
        }

        par_start <- if (L > 1) L else 1
        par_free_logit <- par[par_start:length(par)]

        par.array <- array(0, dim = c(L, I, poly.max))
        idx <- 1
        for (i in 1:I) {
          K_i <- poly.value[i]
          n_free <- K_i - 1
          for (l in 1:L) {
            if (n_free > 0) {
              logits <- par_free_logit[idx:(idx + n_free - 1)]
              probs_free <- plogis(logits)
              probs_free <- pmin(pmax(probs_free, 1e-6), 1 - 1e-6)
              total_free <- sum(probs_free)
              if (total_free >= 1) {
                probs_free <- probs_free / (total_free + 1e-6) * 0.999
              }
              par.array[l, i, 1:n_free] <- probs_free
              par.array[l, i, K_i] <- 1 - sum(par.array[l, i, 1:n_free])
              idx <- idx + n_free
            } else {
              par.array[l, i, 1] <- 1
            }
          }
        }

        -get.Log.Lik.LCA(response, P.Z, par.array)
      }

      hess_result <- tryCatch({
        numDeriv::hessian(
          func = func.nll.logit,
          x = params.vec,
          method = "Richardson",
          method.args = list(eps = 1, d = d_vec, zero.tol = 1e-8, r = 8),
          model.type = model.type,
          response = response,
          I = I,
          L = L,
          poly.value = poly.value
        )
      }, error = function(e) {
        matrix(0, ncol = length(params.vec), nrow = length(params.vec))
      })

      if (!is.matrix(hess_result) || nrow(hess_result) != length(params.vec)) {
        hess.matrix <- matrix(0, ncol = length(params.vec), nrow = length(params.vec))
      } else {
        hess.matrix <- hess_result
      }

      suppressWarnings(cov.mat <- tryCatch({
        hess.matrix <- (hess.matrix + t(hess.matrix)) / 2
        eig <- eigen(hess.matrix, symmetric = TRUE, only.values = TRUE)
        min_eig <- min(eig$values)

        if (min_eig < 1e-5) {
          max_diag <- max(diag(hess.matrix))
          perturb <- max(1e-5 * max_diag, 1e-7)
          diag(hess.matrix) <- diag(hess.matrix) + perturb
        }

        chol_test <- try(chol(hess.matrix), silent = TRUE)
        if (inherits(chol_test, "try-error")) {
          hess.matrix <- Matrix::nearPD(hess.matrix,
                                        corr = FALSE,
                                        eig.tol = 1e-6,
                                        conv.tol = 1e-8)$mat
        }

        solve(hess.matrix)
      }, error = function(e) {
        matrix(0, ncol = length(params.vec), nrow = length(params.vec))
      }))

      se.vec.logit <- sqrt(pmax(1e-6, diag(cov.mat)))

      se.vec.prob <- numeric(length(se.vec.logit))
      if (L > 1) {
        PZ_free <- params$P.Z[-L]
        se.vec.prob[1:(L-1)] <- PZ_free * (1 - PZ_free) * se.vec.logit[1:(L-1)]
      }

      if (length(param_indices) > 0) {
        start_idx <- if (L > 1) L else 1
        for (i in seq_along(param_indices)) {
          l <- param_indices[[i]][1]
          item <- param_indices[[i]][2]
          cat_idx <- param_indices[[i]][3]
          prob_val <- params$par[l, item, cat_idx]
          deriv <- prob_val * (1 - prob_val)
          se.vec.prob[start_idx + i - 1] <- deriv * se.vec.logit[start_idx + i - 1]
        }
      }

      se.pz <- numeric(L)
      if (L > 1) {
        se.pz[-L] <- se.vec.prob[1:(L-1)]
        se.pz[L] <- 0
      } else {
        se.pz <- c(0)
      }

      se.array <- array(0, dim = c(L, I, poly.max))
      if (length(param_indices) > 0) {
        start_idx <- if (L > 1) L else 1
        for (i in seq_along(param_indices)) {
          l <- param_indices[[i]][1]
          item <- param_indices[[i]][2]
          cat_idx <- param_indices[[i]][3]
          se.array[l, item, cat_idx] <- se.vec.prob[start_idx + i - 1]
        }
      }

      se.list <- list(par = se.array, P.Z = se.pz)

      diagnostics <- list(
        method = "Obs",
        hessian_cond_number = if (!is.null(hess.matrix)) kappa(hess.matrix) else NA,
        step_sizes = d_vec,
        note = "SEs computed in logit space + delta method for improved accuracy in LCA models."
      )
    }
  }

  if (model.type == "LPA") {
    if (!is.null(se.list$means)) {
      rownames(se.list$means) <- rownames(params$means)
      colnames(se.list$means) <- colnames(params$means)
    }
    if (!is.null(se.list$covs)) {
      dimnames(se.list$covs) <- dimnames(params$covs)
    }
    names(se.list$P.Z) <- names(params$P.Z)
  } else if (model.type == "LCA") {
    dimnames(se.list$par) <- dimnames(params$par)
    names(se.list$P.Z) <- names(params$P.Z)
  }

  res <- list(
    se = se.list,
    vcov = if (method == "Obs") cov.mat else NULL,
    hessian = if (method == "Obs") hess.matrix else NULL,
    diagnostics = diagnostics
  )

  res$call <- call
  res$arguments = list(
    object = object, method = method, n.Bootstrap = n.Bootstrap
  )

  class(res) <- "SE"
  return(res)
}
