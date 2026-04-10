#' Compute Posterior Latent Class Probabilities Based on Fixed Parameters
#'
#' Computes posterior probabilities of latent class membership for each observation
#' using fixed conditional response probabilities (\code{par}) and fixed class prior
#' probabilities (\code{P.Z}).
#'
#' @param response Numeric matrix (\eqn{N \times I}) of categorical responses.
#'   Categories are automatically remapped to 0-based integers internally via
#'   \code{\link[LCPA]{adjust.response}}.
#' @param par 3D array (\eqn{L \times I \times K_{\max}}) of fixed conditional response
#'   probabilities where:
#'   \itemize{
#'     \item \eqn{L} = number of latent classes
#'     \item \eqn{I} = number of indicators
#'     \item \eqn{K_{\max}} = maximum categories across indicators
#'   }
#'   \code{par[l, i, k]} = \eqn{P(X_i = k-1 \mid Z = l)} (using 1-based indexing for the array dimension corresponding to category k-1).
#'
#' @param P.Z Vector of length \eqn{L} with fixed class prior probabilities (\eqn{\pi_l}).
#'   These values are used directly without re-estimation.
#'
#' @return Numeric matrix (\eqn{N \times L}) of posterior probabilities.
#'   Rows sum to 1. Columns are named "Class.1", "Class.2", etc.
#'
#' @details
#' Unlike an EM algorithm, this function does NOT iteratively update class prevalences.
#' It performs a single calculation step based on Bayes' theorem:
#' \deqn{
#'   P(Z_n = l \mid \mathbf{X}_n) =
#'   \frac{\pi_l \prod_{i=1}^I P(X_{ni} = x_{ni} \mid Z_n = l)}
#'        {\sum_{k=1}^L \pi_k \prod_{i=1}^I P(X_{ni} = x_{ni} \mid Z_n = k)}
#' }
#' where \eqn{\pi_l} are the fixed priors provided in the \code{P.Z} argument.
#'
#'
#' @examples
#' \donttest{
#' library(LCPA)
#' set.seed(123)
#' # Simulate data
#' data.obj <- sim.LCA(N = 200, I = 3, L = 2, IQ = 0.85)
#'
#' # Fit a model to get parameters
#' fit <- LCA(data.obj$response, L = 2, method = "EM", nrep = 5)
#'
#' # Calculate posteriors using fixed parameters from the fitted model
#' P.Z.Xn <- get.P.Z.Xn.LCA(
#'   response = data.obj$response,
#'   par = fit$params$par,
#'   P.Z = fit$params$P.Z
#' )
#' head(P.Z.Xn)
#' }
#'
#' @export
get.P.Z.Xn.LCA <- function(response, par, P.Z){
  adjust.response.obj <- adjust.response(response)
  response <- adjust.response.obj$response
  poly.max <- adjust.response.obj$poly.max

  Y <- as.matrix(response)
  N <- nrow(Y)
  I <- ncol(Y)

  Y_hot <- array(0, dim = c(N, I, poly.max))
  idx <- cbind(c(row(Y)), c(col(Y)), c(Y+1))
  Y_hot[idx] <- 1

  L <- dim(par)[1]

  P.Z.Xn <- matrix(1/L, N, L)
  L.Xi.Z <- matrix(1, L, N)

  for(p in 1:N){
    for(i in 1:I){
      L.Xi.Z[ , p] <- L.Xi.Z[ , p] * par[ , i, Y[p, i]+1]
    }

    L.Xi <- sum(L.Xi.Z[, p] * P.Z)
    P.Z.Xn[p, ] <- (L.Xi.Z[, p] * P.Z + 1e-50) / (L.Xi + 2e-50)
  }

  colnames(P.Z.Xn) <- paste0("Class.", 1:L)

  return(P.Z.Xn)
}
