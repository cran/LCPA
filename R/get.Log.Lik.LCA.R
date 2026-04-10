#' Calculate Log-Likelihood for Latent Class Analysis
#'
#' Computes the log-likelihood of observed categorical data under a Latent Class Analysis (LCA) model
#' given class probabilities and conditional response probabilities. The calculation assumes local independence
#' of responses conditional on latent class membership.
#'
#' @param response A numeric matrix of dimension \eqn{N \times I} containing discrete responses.
#'   Values can be any categorical encoding (e.g., 1/2/3, A/B/C, or 0/1). The function automatically:
#'   \itemize{
#'     \item Converts all responses to 0-based integer encoding internally
#'     \item Determines the maximum number of categories (\eqn{K_{\max}}) across indicators
#'   }
#' @param P.Z A numeric vector of length \eqn{L} containing prior probabilities for latent classes.
#'   Must satisfy:
#'   \itemize{
#'     \item \eqn{\sum_{l=1}^L \pi_l = 1}
#'     \item \eqn{\pi_l > 0} for all \eqn{l = 1, \dots, L}
#'   }
#' @param par A 3-dimensional array of dimension \eqn{L \times I \times K_{\max}} containing conditional probabilities,
#'   where \eqn{par[l, i, k]} represents \eqn{P(X_i = k-1 \mid Z=l)} (after internal 0-based re-encoding).
#'   Must satisfy:
#'   \itemize{
#'     \item For each class \eqn{l} and indicator \eqn{i}: \eqn{\sum_{k=1}^{K_i} par[l,i,k] = 1}
#'     \item Probabilities for non-existent categories (where \eqn{k > K_i}) are ignored but must be present in the array
#'   }
#'
#' @return A single numeric value representing the total log-likelihood:
#'   \deqn{\log \mathcal{L} = \sum_{n=1}^N \log \left[ \sum_{l=1}^L \pi_l \prod_{i=1}^I P(X_{ni} = x_{ni} \mid Z=l) \right]}
#'   where \eqn{x_{ni}} is the standardized (0-based) response for person \eqn{n} on indicator \eqn{i}.
#'
#' @details The log-likelihood calculation follows these steps:
#'
#' \itemize{
#'   \item Response Standardization:
#'   Original responses are converted to 0-based integers
#'   using \code{\link[LCPA]{adjust.response}}.
#'   For example, original values \{1,2,5\} become \{0,1,2\}
#'   (ordered and relabeled sequentially).
#'
#'   \item Class-Specific Likelihood:
#'   For each observation \eqn{n} and class \eqn{l}, compute:
#'   \deqn{P(\mathbf{X}_n \mid Z_n=l) = \prod_{i=1}^I P(X_{ni} = x_{ni} \mid Z_n=l)}
#'   where \eqn{x_{ni}} is the standardized response value, and probabilities are taken from \code{par[l, i, x_{ni}+1]}.
#'
#'   \item Marginal Likelihood:
#'   For each observation \eqn{n}, combine class-specific likelihoods weighted by class probabilities:
#'   \deqn{P(\mathbf{X}_n) = \sum_{l=1}^L \pi_l \cdot P(\mathbf{X}_n \mid Z_n=l)}
#'
#'   \item Log Transformation:
#'   Sum log-transformed marginal likelihoods across all observations:
#'   \deqn{\log \mathcal{L} = \sum_{n=1}^N \log P(\mathbf{X}_n)}
#' }
#'
#'
#' @export
get.Log.Lik.LCA <- function(response, P.Z, par){
  adjust.response.obj <- adjust.response(response)
  response <- adjust.response.obj$response
  poly.max <- adjust.response.obj$poly.max
  poly.value <- adjust.response.obj$poly.value
  poly.orig <- adjust.response.obj$poly.orig

  Y <- as.matrix(response)
  N <- nrow(Y)
  I <- ncol(Y)
  L <- length(P.Z)

  par_vec <- as.vector(par)
  P_Z_vec <- as.vector(P.Z)
  Y_int <- matrix(as.integer(Y), nrow = N, ncol = I)
  e_step_res <- em_e_step(Y_int, par_vec, P_Z_vec, N, I, L, poly.max)
  L.Xi.Z <- e_step_res$L_Xi_Z
  P.Z.Xn <- e_step_res$P_Z_Xn
  L.Xi_vec <- e_step_res$L_Xi_vec
  total_log_lik <- sum(log(L.Xi_vec))

  return(total_log_lik)
}
