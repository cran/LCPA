#' Compute Posterior Latent Profile Probabilities Based on Fixed Parameters
#'
#' Computes posterior probabilities of latent profile membership for each observation
#' using fixed profile parameters (means, covariances) and fixed prior probabilities.
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
#'   Each slice must be symmetric and positive definite.
#' @param P.Z Vector of length \eqn{L} with fixed profile prior probabilities (\eqn{\pi_l}).
#'   These values are used directly without re-estimation.
#'
#' @return Numeric matrix (\eqn{N \times L}) of posterior probabilities.
#'   Rows sum to 1. Columns are named "Class.1", "Class.2", etc.
#'
#' @details
#' Unlike an EM algorithm, this function does NOT iteratively update profile prevalences.
#' It performs a single E-step calculation:
#' \deqn{
#'   P(Z_n = l \mid \mathbf{x}_n) =
#'   \frac{\pi_l \cdot \mathcal{N}(\mathbf{x}_n \mid \boldsymbol{\mu}_l, \boldsymbol{\Sigma}_l)}
#'        {\sum_{k=1}^L \pi_k \cdot \mathcal{N}(\mathbf{x}_n \mid \boldsymbol{\mu}_k, \boldsymbol{\Sigma}_k)}
#' }
#' where \eqn{\pi_l} are the fixed priors provided in the \code{P.Z} argument.
#'
#'
#' @examples
#' \donttest{
#' library(LCPA)
#' set.seed(123)
#' data.obj <- sim.LPA(N = 300, I = 2, L = 2, constraint = "VV")
#' fit <- LPA(data.obj$response, L = 2, method = "EM", nrep = 5)
#'
#' # Calculate posteriors using fixed parameters from a fitted model
#' P.Z.Xn <- get.P.Z.Xn.LPA(
#'   response = data.obj$response,
#'   means = fit$params$means,
#'   covs = fit$params$covs,
#'   P.Z = fit$params$P.Z
#' )
#' head(P.Z.Xn)
#' }
#'
#' @export
get.P.Z.Xn.LPA <- function (response, means, covs, P.Z){
  if (!is.matrix(response))
    stop("response must be a matrix")

  I <- ncol(response)
  L <- nrow(means)
  N <- nrow(response)

  tresponse <- t(response)
  logres <- matrix(NA_real_, nrow = N, ncol = L)

  # Compute log-likelihood contributions for each class
  # Formula: log(P(Z=l)) + log(N(x | mu_l, Sigma_l))
  for (l in 1:L) {
    covs.l <- covs[, , l]
    mean.l <- means[l, ]
    # Assumes logpdf_component handles the multivariate normal log-density calculation
    lp <- logpdf_component(mean.l, covs.l, tresponse)
    logres[, l] <- lp + log(P.Z[l] + 1e-12)
  }

  # Log-Sum-Exp trick for numerical stability
  rowmax <- apply(logres, 1, max)

  # Handle cases where max is non-finite (e.g., all -Inf)
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

  # Convert to probabilities: exp(log_val - max_log_val) / sum(...)
  exp_rel <- exp(logres - matrix(rowmax, N, L))
  row_sums <- rowSums(exp_rel)

  # Handle numerical errors where sum is 0 or NaN
  bad <- !is.finite(row_sums) | row_sums < 1e-20
  if (any(bad)) {
    exp_rel[bad, ] <- 1/L
    row_sums[bad] <- 1
  }

  P.Z.Xn <- exp_rel / row_sums

  colnames(P.Z.Xn) <- paste0("Class.", 1:L)

  return(P.Z.Xn)
}
