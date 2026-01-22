#' Compute Posterior Latent Profile Probabilities Based on Fixed Parameters
#'
#' Computes posterior probabilities of latent profile membership while simultaneously
#' re-estimating profile prevalences via an EM algorithm. Unlike standard posterior
#' computation, this function iteratively updates profile prevalences (\eqn{\pi_l})
#' using fixed profile characteristics (means and covariances).
#'
#' @param response Numeric matrix (\eqn{N \times I}) of continuous responses.
#'   Missing values are not allowed. Data should typically be standardized prior to analysis.
#' @param means Numeric matrix (\eqn{L \times I}) of fixed profile means where:
#'   \itemize{
#'     \item \eqn{L} = number of latent profiles
#'     \item \eqn{I} = number of observed variables
#'   }
#'   Row \eqn{l} contains profile-specific means \eqn{\boldsymbol{\mu}_l}.
#' @param covs 3D array (\eqn{I \times I \times L}) of fixed profile covariance matrices where:
#'   \itemize{
#'     \item \code{covs[, , l]} = profile-specific covariance matrix \eqn{\boldsymbol{\Sigma}_l}
#'   }
#'   Each slice must be symmetric and positive definite (after jittering).
#' @param tol Convergence tolerance for absolute change in log-likelihood. Default: 1e-10.
#' @param maxiter Maximum EM iterations. Default: 2000.
#' @param vis Logical: show iteration progress? Default: TRUE.
#'
#' @return Numeric matrix (\eqn{N \times L}) of posterior probabilities.
#'   Rows sum to 1. Columns named "Class.1", "Class.2", etc.
#'
#' @section Algorithm:
#' The function implements an EM algorithm with:
#' \describe{
#'   \item{E-step}{Compute posterior probabilities for observation \eqn{n} and profile \eqn{l}:
#'     \deqn{
#'       P(Z_n = l \mid \mathbf{x}_n) =
#'       \frac{\pi_l^{(t)} \cdot \mathcal{N}(\mathbf{x}_n \mid \boldsymbol{\mu}_l, \boldsymbol{\Sigma}_l)}
#'            {\sum_{k=1}^L \pi_k^{(t)} \cdot \mathcal{N}(\mathbf{x}_n \mid \boldsymbol{\mu}_k, \boldsymbol{\Sigma}_k)}
#'     }}
#'   \item{M-step}{Update profile prevalences:
#'     \deqn{
#'       \pi_l^{(t+1)} = \frac{1}{N} \sum_{n=1}^N P(Z_n = l \mid \mathbf{x}_n)
#'     }}
#' }
#' Convergence is determined by the absolute change in log-likelihood between iterations.
#'
#' @details
#' \enumerate{
#'   \item Numerical stability:
#'        \itemize{
#'          \item Covariance matrices are jittered with \code{tol} for positive definiteness
#'          \item Log-space computation with log-sum-exp trick
#'          \item Uniform probabilities used as fallback for non-finite densities
#'        }
#'   \item Profile prevalences are initialized uniformly (\eqn{\pi_l^{(0)} = 1/L}).
#'   \item Termination occurs when:
#'        \itemize{
#'          \item \eqn{|\log L^{(t)} - \log L^{(t-1)}| < \code{tol}} (log-likelihood change)
#'          \item Maximum iterations reached
#'        }
#' }
#'
#' @examples
#' \donttest{
#' library(LCPA)
#' set.seed(123)
#' data.obj <- sim.LPA(N = 300, I = 2, L = 2, constraint = "VV")  # From LCPA
#' fit <- LPA(data.obj$response, L = 2, method = "EM", nrep = 5)  # From LCPA
#'
#' P.Z.Xn <- get.P.Z.Xn.LPA(
#'   response = data.obj$response,
#'   means = fit$params$means,    # Fixed profile means
#'   covs = fit$params$covs       # Fixed profile covariances
#' )
#' head(P.Z.Xn)
#' }
#'
#' @export
get.P.Z.Xn.LPA <- function (response, means, covs, tol = 1e-10, maxiter=2000, vis=TRUE){
  if (!is.matrix(response))
    stop("response must be a matrix")

  I <- ncol(response)
  L <- nrow(means)

  N <- nrow(response)
  tresponse <- t(response)
  logres <- matrix(NA_real_, nrow = N, ncol = L)

  nt_width <- ceiling(log10(N * I * L))
  total_width <- nt_width + 5
  fmt_string_maxchg <- sprintf("%%%d.%df", total_width, 5)

  P.Z <- rep(1/L, L)
  iter <- 1
  L.X.pre <- -Inf
  while(iter < maxiter){

    iter <- iter + 1

    for (l in 1:L) {
      covs.l <- covs[, , l]
      mean.l <- means[l, ]
      lp <- logpdf_component(mean.l, covs.l, tresponse, tol)
      logres[, l] <- lp + log(P.Z[l] + 1e-12)
    }

    rowmax <- apply(logres, 1, max)
    nonfinite_rows <- !is.finite(rowmax)
    if (any(nonfinite_rows)) {
      finite_vals <- logres[is.finite(logres)]
      if (length(finite_vals) > 0) {
        rowmax[nonfinite_rows] <- max(finite_vals)
      }
      else {
        rowmax[nonfinite_rows] <- 0
      }
    }

    exp_rel <- exp(logres - matrix(rowmax, N, L))
    row_sums <- rowSums(exp_rel)
    bad <- !is.finite(row_sums) | row_sums < 1e-20
    if (any(bad)) {
      exp_rel[bad, ] <- 1/L
      row_sums[bad] <- 1
    }
    P.Z.Xn <- exp_rel/row_sums

    P.Z <- apply(P.Z.Xn, 2, sum) / N

    L.X.cur <- sum(rowmax + log(rowSums(exp(logres - rowmax))))
    maxchg <- if(iter > 1) abs(L.X.cur - L.X.pre) else Inf
    if(vis){
      cat('\rIter =', sprintf("%4d", iter), '  \u0394Log.Lik =', sprintf(fmt_string_maxchg, maxchg))
    }
    if(maxchg < tol){
      break
    }
    L.X.pre <- L.X.cur
  }

  colnames(P.Z.Xn) <- paste0("Class.", 1:L)

  if(vis){
    cat("\n\n")
  }
  return(P.Z.Xn)
}
