#' Calculate Average Posterior Probability (AvePP)
#'
#' Computes the average posterior probability for the most likely class assignment
#' in latent class/profile analysis. This metric quantifies classification precision.
#' The total average posterior probability \eqn{\geq 0.70} (Nylund-Gibson & Choi, 2018) indicate adequate classification quality.
#'
#' @param object An object of class \code{"LCA"} or \code{"LPA"} returned by
#'   \code{\link[LCPA]{LCA}} or \code{\link[LCPA]{LPA}}, or any object containing:
#'   \itemize{
#'     \item \code{P.Z.Xn}: \eqn{N \times L} matrix of posterior class probabilities, where:
#'       \itemize{
#'         \item \eqn{N} = Total number of observations (\eqn{n = 1, 2, \dots, N})
#'         \item \eqn{L} = Number of latent classes (\eqn{l = 1, 2, \dots, L})
#'         \item Element \eqn{p_{nl} = P(Z_n = l \mid \mathbf{X}_n)} denotes the posterior probability
#'           that observation \eqn{n} belongs to class \eqn{l} given observed data \eqn{\mathbf{X}_n}
#'       }
#'   }
#'
#' @return A \eqn{(L+1) \times (L+1)} matrix with the following structure:
#'   \itemize{
#'     \item Rows: Represent each latent class (1 to L) and a final "Total" row.
#'     \item Columns: Represent each latent class (1 to L) and a final "Total" column.
#'     \item Diagonal elements \eqn{\text{ave}[l,l]}: Average posterior probability for observations assigned to class \eqn{l}.
#'       That is, \deqn{\overline{P}_{ll} = \frac{1}{N_l} \sum_{n: \hat{z}_n = l} p_{nl},}
#'       where \eqn{N_l} is the number of observations assigned to class \eqn{l}, and \eqn{\hat{z}_n = \arg\max_{l'} p_{nl'}}.
#'     \item Off-diagonal elements \eqn{\text{ave}[l,k]} (\eqn{l \ne k}): Average posterior probability of class \eqn{k}
#'       among observations assigned to class \eqn{l}. Useful for assessing classification confusion.
#'       \deqn{\overline{P}_{lk} = \frac{1}{N_l} \sum_{n: \hat{z}_n = l} p_{nk}.}
#'     \item Bottom-right corner \eqn{\text{ave}[L+1,L+1]}: Overall average posterior probability across all observations,
#'       \deqn{\overline{P}_{\text{total}} = \frac{1}{N} \sum_{n=1}^N \max_{l} p_{nl}.}
#'   }
#'
#' @note Classification quality is considered acceptable if \eqn{\overline{P}_{\text{total}} \geq 0.70} (Nylund-Gibson & Choi, 2018).
#'
#' @references
#' Nylund-Gibson, K., & Choi, A. Y. (2018). Ten frequently asked questions about latent class analysis. Translational Issues in Psychological Science, 4(4), 440-461. https://doi.org/10.1037/tps0000176
#'
#' @examples
#' # Example with simulated data
#' set.seed(123)
#' data.obj <- sim.LCA(N = 500, I = 4, L = 2, IQ=0.9)
#' response <- data.obj$response
#'
#' # Fit 2-class model with EM algorithm
#' \donttest{
#' fit.em <- LCA(response, L = 2, method = "EM", nrep = 10)
#'
#' AvePP_value <- get.AvePP(fit.em)
#' print(AvePP_value)
#'
#' }
#'
#' @export
get.AvePP <- function(object) {
  if (!("P.Z.Xn" %in% names(object))) {
    stop("Input object missing required element 'P.Z.Xn' (posterior probabilities)")
  }

  P.Z.Xn <- object$P.Z.Xn

  if (!is.matrix(P.Z.Xn)) {
    stop("'P.Z.Xn' must be a matrix")
  }

  max.pp <- apply(P.Z.Xn, 1, max)
  ave.total <- mean(max.pp)

  L <-ncol(P.Z.Xn)
  ave <- matrix(NA, L+1, L+1)
  Z <- apply(P.Z.Xn, 1, which.max)
  for(l in 1:L){
    ave[l, 1:L] <- colMeans(P.Z.Xn[which(Z == l), , drop=FALSE], na.rm = TRUE)
  }
  ave[L+1, L+1] <- ave.total

  if(is.null(colnames(P.Z.Xn))){
    rownames(ave) <- colnames(ave) <- c(paste0("Class.", 1:L), "Total")
  }else{
    rownames(ave) <- colnames(ave) <- c(colnames(P.Z.Xn), "Total")
  }

  return(ave)
}
