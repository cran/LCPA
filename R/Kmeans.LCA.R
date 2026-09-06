#' Initialize LCA Parameters via K-means Clustering
#'
#' Performs hard clustering of observations using K-means algorithm to generate
#' initial parameter estimates for Latent Class Analysis (LCA) models. This
#' provides a data-driven initialization strategy that often outperforms random
#' starts when the number of observed categorical variables \eqn{I} is large
#' (i.e., \eqn{I > 50}).
#'
#' @param response A numeric matrix of dimension \eqn{N \times I}, where \eqn{N} is the number of observations
#'   and \eqn{I} is the number of observed categorical variables. Each column must contain nominal-scale
#'   discrete responses (e.g., integers representing categories). Non-sequential category values are
#'   automatically re-encoded to sequential integers starting from 1.
#' @param L Integer specifying the number of latent classes. Must be \eqn{2 \leq L < N}.
#' @param starts Integer specifying the number of random starts for K-means algorithm
#'   (default: 1). The solution with the lowest within-cluster sum of squares is retained.
#'
#' @details
#' The function executes the following steps:
#' \itemize{
#'   \item Data preprocessing: Automatically adjusts non-sequential category values
#'         to sequential integers (e.g., categories \{1,3,5\} become \{1,2,3\}) using internal adjustment routines.
#'   \item K-means clustering: Scales variables to mean=0 and SD=1 before clustering.
#'         Uses Lloyd's algorithm with Euclidean distance.
#'   \item Parameter estimation:
#'     \itemize{
#'       \item For each cluster \eqn{l}, computes empirical response probabilities
#'             \eqn{P(X_i=q\mid Z=l)} for all indicators \eqn{i} and response
#'             categories \eqn{q}.
#'     }
#'   \item Posterior probabilities: Constructs hard-classification matrix where
#'         \eqn{P(Z_n=l\mid\mathbf{X}_n)=1} for the assigned cluster and 0 otherwise.
#' }
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{params}}{List of initialized parameters:
#'     \describe{
#'       \item{\code{par}}{An \eqn{L \times I \times K_{\max}} array of initial conditional probabilities,
#'                         where \eqn{K_{\max}} is the maximum number of categories across indicators.
#'                         Dimension order: latent classes (1:L), indicators (1:I), response categories (1:K_max).}
#'       \item{\code{P.Z}}{Numeric vector of length \eqn{L} containing initial class prior probabilities
#'                         derived from cluster proportions.}
#'     }
#'   }
#'   \item{\code{P.Z.Xn}}{An \eqn{N \times L} matrix of posterior class probabilities. Contains
#'                         hard assignments (0/1 values) based on K-means cluster memberships.}
#' }
#' @note This function is primarily designed as an initialization method for
#'   \code{\link[LCPA]{LCA}()} and not for final model estimation.
#'
#' @examples
#' # Simulate response data
#' set.seed(123)
#' response <- matrix(sample(1:4, 200, replace = TRUE), ncol = 5)
#'
#' # Generate K-means initialization for 3-class LCA
#' init_params <- Kmeans.LCA(response, L = 3, starts = 5)
#'
#' # Inspect initial class probabilities
#' print(init_params$params$P.Z)
#' @export
#' @importFrom stats kmeans
Kmeans.LCA <- function(response, L, starts=1){
  adjust.response.obj <- adjust.response(response)
  response <- adjust.response.obj$response
  poly.value <- adjust.response.obj$poly.value

  cluster.res <- kmeans(
    x = scale(response),
    centers = L,
    iter.max = 1000,
    nstart = starts,
    algorithm = "Lloyd"
  )

  Z <- as.integer(cluster.res$cluster)
  P.Z.Xn <- .class.indicator(Z, L)
  maximization <- lca_maximization_cpp(
    matrix(as.integer(response), nrow(response), ncol(response)),
    P.Z.Xn, as.integer(poly.value), 1e-10
  )
  par <- maximization$par
  P.Z <- as.table(maximization$P.Z)
  names(P.Z) <- 1:L

  res <- list(params=list(
    par=par, P.Z=P.Z
  ),
  P.Z.Xn=P.Z.Xn)

  return(res)
}
