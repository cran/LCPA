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
#'   Rows sum to 1. Columns are named `"Profile 1"`, `"Profile 2"`, and so on.
#'
#' @details
#' Unlike an EM algorithm, this function does NOT iteratively update profile prevalences.
#' It performs a single E-step calculation:
#' \deqn{
#'   \tau_{nl}=P(Z_n=l\mid\mathbf{X}_n)=
#'   \frac{\pi_l\mathcal{N}(\mathbf{X}_n\mid
#'   \boldsymbol{\mu}_l,\boldsymbol{\Sigma}_l)}
#'   {\sum_{h=1}^L\pi_h\mathcal{N}(\mathbf{X}_n\mid
#'   \boldsymbol{\mu}_h,\boldsymbol{\Sigma}_h)}
#' }
#' where the profile probabilities, means, and covariance matrices are fixed by
#' the \code{P.Z}, \code{means}, and \code{covs} arguments.
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
  if(!is.numeric(response) || any(!is.finite(response)) ||
     !is.matrix(means) || ncol(means) != I || any(!is.finite(means)) ||
     length(dim(covs)) != 3L || !identical(dim(covs), c(I, I, L)) ||
     any(!is.finite(covs))){
    stop("response, means, and covs must contain compatible finite numeric values")
  }
  if(length(P.Z) != L || any(!is.finite(P.Z)) || any(P.Z <= 0) || sum(P.Z) <= 0){
    stop("P.Z must contain L finite positive values")
  }
  stabilized <- .stabilize.LPA.covariances(
    covs, "VV", fallback = stats::cov(response)
  )
  if(stabilized$repaired) covs <- stabilized$covs
  expectation <- .lpa.expectation.reference(
    response, means, covs, as.numeric(P.Z) / sum(P.Z)
  )
  if(!isTRUE(expectation$valid)){
    stop("Unable to obtain positive-definite profile covariance matrices")
  }
  P.Z.Xn <- expectation$posterior
  colnames(P.Z.Xn) <- .latent.group.names(L, "LPA")
  P.Z.Xn
}
