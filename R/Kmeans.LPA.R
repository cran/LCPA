#' Initialize LPA Parameters via K-means Clustering
#'
#' Performs hard clustering with K-means and converts the resulting partition into
#' initial means, covariance matrices, profile proportions, and hard posterior
#' assignments for Latent Profile Analysis (LPA).
#'
#' @param response A numeric matrix or data frame of dimension \eqn{N \times I} containing
#'   continuous indicators. Missing and non-finite values are not allowed. The indicators
#'   must be standardized before input using \code{\link[base]{scale}} or
#'   \code{\link[LCPA]{normalize}}; no additional scaling is performed internally.
#' @param L Positive integer specifying the number of latent profiles. It must be smaller
#'   than the number of observations.
#' @param constraint Covariance structure used to construct the initial profile covariance
#'   matrices. Named options are \code{"UE"}, \code{"UV"}, \code{"E0"}, \code{"V0"},
#'   \code{"EE"}, \code{"VV"}, \code{"VE"}, and \code{"EV"}. A custom list of index
#'   pairs may be supplied to constrain selected variance or covariance elements equal
#'   across profiles, using the same semantics as \code{\link[LCPA]{LPA}}.
#' @param starts Positive integer specifying the number of internal K-means random starts
#'   (default: 1). The solution with the lowest within-cluster sum of squares is retained.
#' @details
#' The function performs four operations:
#' \itemize{
#'   \item Runs K-means on the supplied standardized indicators using Lloyd's algorithm.
#'   \item Uses the cluster centers as initial profile means and cluster proportions as
#'     initial profile probabilities.
#'   \item Computes within-cluster maximum-likelihood covariance matrices, then applies
#'     the requested named or custom equality constraints.
#' }
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{params}}{A list with:
#'     \describe{
#'       \item{\code{means}}{An \eqn{L \times I} matrix of initial profile means.}
#'       \item{\code{covs}}{An \eqn{I \times I \times L} array of constrained,
#'         positive-definite initial covariance matrices.}
#'       \item{\code{P.Z}}{A numeric vector of length \eqn{L} containing initial profile
#'         proportions.}
#'     }
#'   }
#'   \item{\code{P.Z.Xn}}{An \eqn{N \times L} hard-assignment matrix whose rows contain
#'     one 1 and \eqn{L-1} zeros.}
#' }
#'
#' @note This function constructs initialization values only; it does not fit a final LPA
#'   model. When called directly with \code{starts > 1}, K-means selects one solution from
#'   those internal starts. The package-level \code{LPA(..., par.ini = "kmeans")} workflow
#'   instead calls \code{Kmeans.LPA(..., starts = 1)} separately for every outer start.
#'
#' @examples
#' set.seed(123)
#' response <- scale(matrix(rnorm(300), ncol = 3))
#' initial <- Kmeans.LPA(response, L = 2, constraint = "V0")
#' initial$params$means
#' initial$params$P.Z
#'
#' @export
#' @importFrom stats cov kmeans
Kmeans.LPA <- function(response, L, constraint = "VV", starts = 1){
  if(is.vector(response) && is.numeric(response)) response <- matrix(response, ncol = 1L)
  if(is.data.frame(response)) response <- as.matrix(response)
  if(!is.matrix(response) || !is.numeric(response)){
    stop("response must be a numeric matrix, data frame, or vector")
  }
  if(any(!is.finite(response))){
    stop("response must not contain missing or non-finite values")
  }

  N <- nrow(response)
  I <- ncol(response)
  if(N < 2L || I < 1L){
    stop("response must contain at least two observations and one indicator")
  }
  if(length(L) != 1L || !is.finite(L) || L < 1L || L != as.integer(L) || L >= N){
    stop("L must be a positive integer smaller than nrow(response)")
  }
  if(length(starts) != 1L || !is.finite(starts) || starts < 1L ||
     starts != as.integer(starts)){
    stop("starts must be a positive integer")
  }
  constraint <- .validate.LPA.constraint(constraint, I)

  cluster.res <- tryCatch(
    stats::kmeans(
      response, centers = as.integer(L), iter.max = 1000L,
      nstart = as.integer(starts), algorithm = "Lloyd"
    ),
    error = function(error){
      stop("K-means LPA initialization failed: ", conditionMessage(error), call. = FALSE)
    }
  )

  Z <- as.integer(cluster.res$cluster)
  P.Z.Xn <- .class.indicator(Z, L)
  means <- matrix(cluster.res$centers, nrow = L, ncol = I)
  P.Z <- as.table(tabulate(Z, nbins = L) / N)
  names(P.Z) <- seq_len(L)

  covariance.global <- matrix(stats::cov(response), I, I)
  if(any(!is.finite(covariance.global))){
    variances <- apply(response, 2L, stats::var)
    variances[!is.finite(variances)] <- 1
    covariance.global <- diag(pmax(variances, .Machine$double.eps), I)
  }
  covariance.global <- (covariance.global + t(covariance.global)) / 2

  covs <- array(0, dim = c(I, I, L))
  scatter <- array(0, dim = c(I, I, L))
  sizes <- tabulate(Z, nbins = L)
  for(l in seq_len(L)){
    centered <- sweep(response[Z == l, , drop = FALSE], 2L, means[l, ], "-")
    scatter[, , l] <- crossprod(centered)
    covs[, , l] <- if(sizes[l] > 1L){
      scatter[, , l] / sizes[l]
    }else{
      covariance.global
    }
  }

  pooled <- apply(scatter, c(1L, 2L), sum) / sum(sizes)
  if(is.character(constraint)){
    if(constraint == "E0"){
      for(l in seq_len(L)) covs[, , l] <- diag(diag(pooled), I)
    }else if(constraint == "V0"){
      for(l in seq_len(L)) covs[, , l] <- diag(diag(covs[, , l]), I)
    }else if(constraint %in% c("EE", "UE")){
      for(l in seq_len(L)) covs[, , l] <- pooled
    }else if(constraint == "VE" && I > 1L){
      off.diagonal <- row(pooled) != col(pooled)
      for(l in seq_len(L)){
        covariance.l <- covs[, , l]
        covariance.l[off.diagonal] <- pooled[off.diagonal]
        covs[, , l] <- covariance.l
      }
    }else if(constraint == "EV"){
      for(l in seq_len(L)) diag(covs[, , l]) <- diag(pooled)
    }
  }else{
    for(indices in constraint){
      i <- as.integer(indices[1L])
      j <- as.integer(indices[2L])
      covs[i, j, ] <- covs[j, i, ] <- pooled[i, j]
    }
  }

  covs <- .stabilize.LPA.covariances(
    covs, constraint, fallback = covariance.global
  )$covs

  list(
    params = list(means = means, covs = covs, P.Z = P.Z),
    P.Z.Xn = P.Z.Xn
  )
}
