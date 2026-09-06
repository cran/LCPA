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
#' @param par A 3-dimensional array of dimension \eqn{L \times I \times K_{\max}} containing conditional probabilities,
#'   where \eqn{par[l, i, q]} represents \eqn{P(X_i = q-1 \mid Z=l)}
#'   (after internal 0-based re-encoding).
#'   Must satisfy:
#'   \itemize{
#'     \item For each class \eqn{l} and indicator \eqn{i}:
#'       \eqn{\sum_{q=1}^{K_i} par[l,i,q] = 1}
#'     \item Probabilities for non-existent categories (where \eqn{q > K_i})
#'       are ignored but must be present in the array
#'   }
#' @param P.Z A numeric vector of length \eqn{L} containing prior probabilities for latent classes.
#'   Must satisfy:
#'   \itemize{
#'     \item \eqn{\sum_{l=1}^L \pi_l = 1}
#'     \item \eqn{\pi_l > 0} for all \eqn{l = 1, \dots, L}
#'   }
#'
#' @return A single numeric value equal to the total observed-data
#'   log-likelihood \eqn{\log\mathcal{L}_{\mathrm{LCA}}} defined below.
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
#'   \item Class-Conditional Likelihood Contribution:
#'   For each participant \eqn{n} and class \eqn{l}, local independence gives
#'   \deqn{P(\mathbf{X}_n=\mathbf{x}_n\mid Z_n=l)
#'   =\prod_{i=1}^I P(X_{ni}=x_{ni}\mid Z_n=l),}
#'   where \eqn{x_{ni}} is the standardized response value, and probabilities
#'   are taken from \code{par[l, i, x_{ni}+1]}.
#'
#'   \item Participant-Level Marginal Likelihood Contribution:
#'   For each participant \eqn{n}, combine class-specific likelihoods weighted
#'   by class probabilities:
#'   \deqn{P(\mathbf{X}_n=\mathbf{x}_n)=
#'   \sum_{l=1}^L\pi_l\prod_{i=1}^I
#'   P(X_{ni}=x_{ni}\mid Z_n=l).}
#'
#'   \item Total Observed-Data Log-Likelihood:
#'   Sum the log marginal contributions across participants:
#'   \deqn{\log\mathcal{L}_{\mathrm{LCA}}=
#'   \sum_{n=1}^N\log\left\{\sum_{l=1}^L\pi_l
#'   \prod_{i=1}^I P(X_{ni}=x_{ni}\mid Z_n=l)\right\}.}
#' }
#'
#'
#' @export
get.Log.Lik.LCA <- function(response, par, P.Z){
  adjust.response.obj <- adjust.response(response)
  response <- adjust.response.obj$response
  poly.max <- adjust.response.obj$poly.max
  poly.value <- adjust.response.obj$poly.value
  poly.orig <- adjust.response.obj$poly.orig

  Y <- as.matrix(response)
  Y.int <- matrix(as.integer(Y), nrow(Y), ncol(Y))
  expectation <- lca_expectation_cpp(Y.int, par, as.numeric(P.Z))
  if(!isTRUE(expectation$valid)){
    stop("At least one response pattern has zero probability in every class")
  }
  expectation$Log.Lik
}
