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
#' @param nrep Integer specifying the number of random starts for K-means algorithm
#'   (default: 10). The solution with the lowest within-cluster sum of squares is retained.
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
#'             \eqn{P(X_i=k|Z=l)} for all indicators \eqn{i} and categories \eqn{k}.
#'       \item Handles singleton clusters by assigning near-deterministic probabilities
#'             (e.g., \eqn{1-10^{-10}} for observed category, \eqn{10^{-10}} for others).
#'     }
#'   \item Posterior probabilities: Constructs hard-classification matrix where
#'         \eqn{P(Z=l|\mathbf{X}_n)=1} for the assigned cluster and 0 otherwise.
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
#' @note
#' \itemize{
#'   \item Requires at least one observation per cluster. If a cluster has only one observation,
#'         probabilities are set to avoid zero values (using \eqn{10^{-10}}) for numerical stability.
#'   \item Data scaling is applied internally. Variables with zero variance are automatically
#'         excluded from clustering.
#'   \item This function is primarily designed as an initialization method for \code{\link{LCA}} and not
#'         intended for final model estimation.
#' }
#'
#' @examples
#' # Simulate response data
#' set.seed(123)
#' response <- matrix(sample(1:4, 200, replace = TRUE), ncol = 5)
#'
#' # Generate K-means initialization for 3-class LCA
#' init_params <- Kmeans.LCA(response, L = 3, nrep = 5)
#'
#' # Inspect initial class probabilities
#' print(init_params$params$P.Z)
#' @export
#' @importFrom stats kmeans
Kmeans.LCA <- function(response, L, nrep=10){
  adjust.response.obj <- adjust.response(response)
  response <- adjust.response.obj$response
  poly.max <- adjust.response.obj$poly.max
  poly.value <- adjust.response.obj$poly.value
  poly.orig <- adjust.response.obj$poly.orig

  N <- nrow(response)
  I <- ncol(response)
  par <- array(NA, dim=c(L, I, poly.max))

  cluster.res <- kmeans(
    x = scale(response),
    centers = L,
    iter.max = 1000,
    nstart = nrep,
    algorithm = "Lloyd"
  )

  Z <- cluster.res$cluster
  for(l in 1:L){

    participants.cur <- response[which(Z == l), ]

    if(length(which(Z == l)) > 1){
      for(i in 1:I){
        prop.cur <- table(participants.cur[ , i])
        prop.cur <- prop.cur / sum(prop.cur)

        prop.i <- rep(0, poly.value[i])
        prop.i[as.numeric(names(prop.cur))+1] <- prop.cur
        par[l, i, 1:poly.value[i]] <- prop.i
      }
    }else{
      for(i in 1:I){
        par[l, i, participants.cur[i]] <- 1-1e-10
        par[l, i, c(1:poly.value[i])[participants.cur[i]]] <- 1e-10
      }
    }
  }

  P.Z.Xn <- t(sapply(Z, function(x) {
    temp <- rep(0, L)
    temp[x] <- 1
    return(temp)
  }))

  P.Z <- as.table(colSums(P.Z.Xn) / sum(P.Z.Xn))
  names(P.Z) <- 1:L

  res <- list(params=list(
    par=par, P.Z=P.Z
  ),
  P.Z.Xn=P.Z.Xn)

  return(res)
}
