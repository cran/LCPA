#' Compute Classification Error Probability (CEP) Matrices
#'
#' Computes the Classification Error Probability (CEP) matrices (Liang et al., 2023) used in the bias-corrected
#' three-step estimation of Latent Class/Profile Analysis with Covariates.
#'
#' @param P.Z.Xns A list of length \eqn{T} (number of time points). Each element is an
#'   \eqn{N \times L} matrix of posterior probabilities
#'   \eqn{P(Z_{it} = l \mid X_i)} from the first-step model.
#'   \itemize{
#'     \item Rows correspond to individuals (\eqn{i = 1, \dots, N});
#'     \item Columns correspond to latent classes (\eqn{l = 1, \dots, L});
#'     \item Each row must sum to 1.
#'   }
#'   The list must be ordered chronologically (e.g., time 1 to \eqn{T}).
#' @param time.cross Logical. If \code{TRUE} (default), returns a list where every element
#'   is the same pooled CEP matrix (averaged across all time points). If \code{FALSE},
#'   returns time-specific CEP matrices.
#'
#' @return A named list of length \eqn{T}. Each element is an \eqn{L \times L} matrix:
#'   \itemize{
#'     \item Row \eqn{l}: true latent class;
#'     \item Column \eqn{l'}: individuals assigned to class \eqn{l'};
#'     \item Entry \eqn{(l, l')}: estimated
#'           \eqn{P(\text{assigned class} = l' \mid \text{true class} = l)}.
#'   }
#'
#'   When \code{time.cross = TRUE}, all matrices in the list are identical.
#'   Names are \code{"t1"}, \code{"t2"}, \dots, \code{"tT"}.
#'
#' @details
#' The CEP matrix at time \eqn{t} gives the probability that an individual truly belongs
#' to latent class \eqn{l'} given that they were assigned (via modal assignment)
#' to class \eqn{l} at time \eqn{t}.
#'
#' Formally, for time point \eqn{t}:
#' \deqn{
#' \mathrm{CEP}_t(l, l') =
#' P(Z_t = l \mid \hat{Z}_t = l')
#' =
#' \frac{
#'   \sum_{i:\,\hat{z}_{it} = l'}
#'   P(Z_{it} = l \mid X_i)
#' }{
#'   N \, \hat{\pi}_{tl}
#' }
#' }
#'
#' where:
#' \itemize{
#'   \item \eqn{Z_{it}} is the true latent class of individual \eqn{i} at time \eqn{t};
#'   \item \eqn{P(Z_{it} = l \mid X_i)} is the posterior probability from the first-step model;
#'   \item \eqn{\hat{z}_{it} = \arg\max_l P(Z_{it} = l' \mid X_i)}
#'         is the modal (most likely) assigned class;
#'   \item \eqn{\hat{\pi}_{tl} = \frac{1}{N} \sum_{i=1}^N I(\hat{z}_{it} = l)}
#'         is the observed proportion assigned to class \eqn{l} at time \eqn{t};
#'   \item \eqn{N} is the total sample size.
#' }
#'
#' If \code{time.cross = TRUE} (default), a single pooled CEP matrix is computed by
#' aggregating counts across all time points. This assumes the classification error
#' structure is invariant over time (i.e., measurement invariance), as in
#' Liang et al. (2023). The same pooled matrix is then returned for every time point.
#'
#' @note
#' \itemize{
#'   \item Assumes complete data (no missing values in posterior matrices).
#'   \item All matrices in \code{P.Z.Xns} must have identical dimensions
#'         (same \eqn{N} and \eqn{L}).
#'   \item Assignment is based on modal class (\code{which.max}).
#'   \item If no individual is assigned to a class at a time point,
#'         division by zero may occur.
#' }

#'
#' @references
#' Liang, Q., la Torre, J. d., & Law, N. (2023). Latent Transition Cognitive Diagnosis Model With Covariates: A Three-Step Approach. Journal of Educational and Behavioral Statistics, 48(6), 690-718. https://doi.org/10.3102/10769986231163320
#'
#' @examples
#' # Simulate posterior probabilities for 2 time points, 3 classes, 100 individuals
#' set.seed(123)
#' N <- 100; L <- 3; times <- 2
#' P.Z.Xns <- replicate(times,
#'   t(apply(matrix(runif(N * L), N, L), 1, function(x) x / sum(x))),
#'   simplify = FALSE)
#'
#' # Compute time-specific CEP matrices
#' cep_time_specific <- get.CEP(P.Z.Xns, time.cross = FALSE)
#'
#' # Compute time-invariant (pooled) CEP matrix
#' cep_pooled <- get.CEP(P.Z.Xns, time.cross = TRUE)
#'
#' @export
#'
get.CEP <- function(P.Z.Xns, time.cross=TRUE){

  if(inherits(P.Z.Xns, "list")){
    times <- length(P.Z.Xns)
  }else{
    stop("P.Z.Xns must be a list of length 'times', each element being an N (number of observations) by L (number of latent classes) matrix!")
  }

  N <- nrow(P.Z.Xns[[1]])
  L <- ncol(P.Z.Xns[[1]])

  Zs <- P.Zs <- list()
  for(t in 1:times){
    Zs[[t]] <- apply(P.Z.Xns[[t]], 1, which.max)
    P.Zs[[t]] <- colSums(P.Z.Xns[[t]]) / N
  }

  CEP <- list()
  CEP.sum <- matrix(0, L, L)
  P.Zs.sum <- rep(0, L)
  for(t in 1:times){
    CEP.cur <- matrix(0, L, L)
    for(l in 1:L){
      for(ll in 1:L){
        CEP.cur[l, ll] <- sum(P.Z.Xns[[t]][Zs[[t]] == l, ll])
      }
    }
    CEP[[t]] <- CEP.cur / (N * matrix(P.Zs[[t]], L, L, byrow = FALSE))
    CEP.sum <- CEP.sum + CEP.cur
    P.Zs.sum <- P.Zs.sum + P.Zs[[t]]
  }

  if(time.cross){
    CEP.ave <- CEP.sum / (N*P.Zs.sum)
    for(t in 1:times){
      CEP[[t]] <- CEP.ave
    }
  }

  for(t in 1:times){
    if(!is.null(colnames(P.Z.Xns[[t]]))){
      rownames(CEP[[t]]) <- paste0("True.", colnames(P.Z.Xns[[t]]))
      colnames(CEP[[t]]) <- paste0("Pred.", colnames(P.Z.Xns[[t]]))
    }else{
      rownames(CEP[[t]]) <- paste0("True.", 1:L)
      colnames(CEP[[t]]) <- paste0("Pred.", 1:L)
    }
  }

  names(CEP) <- paste0("t", 1:times)

  return(CEP)
}
