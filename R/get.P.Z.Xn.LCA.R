#' Compute Posterior Latent Class Probabilities Based on Fixed Parameters
#'
#' Computes posterior probabilities of latent class membership while simultaneously
#' re-estimating class prevalences via an EM algorithm. Unlike standard posterior
#' computation, this function iteratively updates class prevalences (\eqn{\pi_l})
#' using fixed conditional response probabilities (\code{par}).
#'
#' @param response Numeric matrix (\eqn{N \times I}) of categorical responses.
#'   Categories are automatically remapped to 0-based integers via
#'   \code{\link[LCPA]{adjust.response}}.
#' @param par 3D array (\eqn{L \times I \times K_{\max}}) of fixed conditional response
#'   probabilities where:
#'   \itemize{
#'     \item \eqn{L} = number of latent classes
#'     \item \eqn{I} = number of items
#'     \item \eqn{K_{\max}} = maximum categories across items
#'   }
#'   \code{par[l, i, k]} = \eqn{P(X_i = k-1 \mid Z = l)} (0-based indexing).
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
#'   \item{E-step}{Compute posterior probabilities for observation \eqn{n} and class \eqn{l}:
#'     \deqn{
#'       P(Z_n = l \mid \mathbf{X}_n) =
#'       \frac{\pi_l^{(t)} \prod_{i=1}^I P(X_{ni} = x_{ni} \mid Z_n = l)}
#'            {\sum_{k=1}^L \pi_k^{(t)} \prod_{i=1}^I P(X_{ni} = x_{ni} \mid Z_n = k)}
#'     }}
#'   \item{M-step}{Update class prevalences:
#'     \deqn{
#'       \pi_l^{(t+1)} = \frac{1}{N} \sum_{n=1}^N P(Z_n = l \mid \mathbf{X}_n)
#'     }}
#' }
#' Convergence is determined by the absolute change in log-likelihood between iterations.
#'
#' @details
#' \enumerate{
#'   \item Response categories are standardized to 0-based integers using
#'         \code{\link[LCPA]{adjust.response}}.
#'   \item Class prevalences are initialized uniformly (\eqn{\pi_l^{(0)} = 1/L}).
#'   \item Numerical stability: Small constants (1e-50) prevent division by zero.
#'   \item Termination occurs when:
#'        \itemize{
#'          \item \eqn{|\log L^{(t)} - \log L^{(t-1)}| < \code{tol}} (log-likelihood change)
#'          \item Maximum iterations reached
#'        }
#' }
#'
#'
#' @examples
#' \donttest{
#' library(LCPA)
#' set.seed(123)
#' data.obj <- sim.LCA(N = 200, I = 3, L = 2, IQ = 0.85)  # From LCPA
#' fit <- LCA(data.obj$response, L = 2, method = "EM", nrep = 5)  # From LCPA
#'
#' P.Z.Xn <- get.P.Z.Xn.LCA(
#'   response = data.obj$response,
#'   par = fit$params$par  # Fixed conditional probabilities
#' )
#' head(P.Z.Xn)
#' }
#'
#' @export
#'
get.P.Z.Xn.LCA <- function(response, par, tol=1e-10, maxiter=2000, vis=TRUE){
  adjust.response.obj <- adjust.response(response)
  response <- adjust.response.obj$response
  poly.max <- adjust.response.obj$poly.max
  poly.value <- adjust.response.obj$poly.value
  poly.orig <- adjust.response.obj$poly.orig

  Y <- as.matrix(response)
  N <- nrow(Y)
  I <- ncol(Y)
  Y_hot <- array(0, dim = c(N, I, poly.max))
  idx <- cbind(c(row(Y)), c(col(Y)), c(Y+1))
  Y_hot[idx] <- 1
  N <- nrow(Y)
  I <- ncol(Y)
  L <- dim(par)[1]


  int_width <- ceiling(log10(N * I * L))
  total_width <- int_width + 5
  fmt_string_maxchg <- sprintf("%%%d.%df", total_width, 5)

  P.Z <- rep(1/L, L)
  iter <- 0
  L.X.pre <- -Inf
  while (iter < maxiter) {
    P.Z.Xn <- matrix(1/L, N, L)
    L.Xi.Z <- matrix(1, L, N)
    iter <- iter + 1

    for(p in 1:N){
      for(i in 1:I){
        L.Xi.Z[ , p] <- L.Xi.Z[ , p] * par[ , i, Y[p, i]+1]
      }
      L.Xi <- sum(L.Xi.Z[, p] * P.Z)
      P.Z.Xn[p, ] <- (L.Xi.Z[, p] * P.Z+1e-50) / (L.Xi+2e-50)
    }
    P.Z <- apply(P.Z.Xn, 2, sum) / N

    L.X.cur <- sum(log(colSums(L.Xi.Z * matrix(P.Z, nrow = L, ncol = N, byrow = FALSE) + 1e-300)))
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
