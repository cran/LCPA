#' Calculate Number of Free Parameters in Latent Transition Analysis
#'
#' Computes the total number of free parameters in a Latent Transition Analysis (LTA) model
#' estimated via the three-step approach. The count depends on the number of latent classes,
#' the number of time points, the number of covariates at each time point, and whether
#' transition coefficients are constrained to be equal across time.
#'
#' @param covariates.ncol An integer vector of length \eqn{T} (number of time points).
#'   Each element equals \eqn{U_t+1}: \eqn{U_t} observed covariates plus the
#'   all-ones intercept column at time \eqn{t}.
#' @param L Integer scalar. Number of latent classes (\eqn{L \geq 2}).
#' @param covariates.time.cross Logical. If \code{TRUE}, transition coefficients are
#'   constrained to be identical across all transitions (time-invariant effects). This requires
#'   that the number of covariates is the same for all time points after the
#'   first (i.e., \eqn{U_2=U_3=\cdots=U_T}).
#'   If \code{FALSE} (default), each transition has its own set of coefficients.
#'
#' @return Integer representing the total number of free parameters:
#'   \deqn{npar = (U_1+1)(L-1) + \begin{cases}
#'     L(L-1)(U_2+1) & \text{if } T>1 \text{ and time-invariant effects} \\
#'     \sum_{t=2}^T L(L-1)(U_t+1) & \text{if } T>1 \text{ and time-varying effects} \\
#'     0 & \text{if } T=1
#'   \end{cases}}
#'   where:
#'   \itemize{
#'     \item \emph{time-invariant effects} corresponds to \code{covariates.time.cross = TRUE}
#'     \item \emph{time-varying effects} corresponds to \code{covariates.time.cross = FALSE}
#'   }
#'
#' @details Parameterization:
#'   \describe{
#'     \item{Initial status model (time 1):}{
#'       Multinomial logit model with \eqn{L} classes (one class is the reference).
#'       Number of free parameters: \eqn{(U_1+1)(L-1)}.
#'     }
#'     \item{Transition models (time \eqn{t \to t+1}):}{
#'       For each transition, a multinomial logit model conditioned on previous class.
#'       For each origin class \eqn{k} and destination class \eqn{l} (\eqn{l \neq L}),
#'       there is a coefficient vector of length \eqn{U_{t+1}+1}.
#'       Total per transition: \eqn{L(L-1)(U_{t+1}+1)} parameters.
#'       The constraint \code{covariates.time.cross} determines whether these parameters
#'       are shared across transitions.
#'     }
#'   }
#'
#' @note Critical assumptions:
#'   \itemize{
#'     \item One latent class is selected as the reference category for all logits; this choice does not change the count.
#'     \item When \code{covariates.time.cross = TRUE}, it is assumed that all time points after the first
#'           have identical covariate structures
#'           (\eqn{U_2=U_3=\cdots=U_T}). If violated, the function
#'           requires the transition design matrices at times 2 through \eqn{T} to have the same number of columns.
#'     \item For \eqn{T=1}, no transition parameters are estimated (pure latent class/profile analysis).
#'   }
#'
#' @examples
#' # Example 1: 2 time points, 2 classes, time-invariant transition coefficients
#' #   Time1: 2 covariates (intercept + 1 predictor)
#' #   Time2: 3 covariates; this determines the shared transition block
#' covariates.ncol <- c(2, 3)
#' L <- 2
#' get.npar.LTA(covariates.ncol, L, covariates.time.cross = TRUE)
#'
#' # Example 2: Same as above but time-varying coefficients
#' get.npar.LTA(covariates.ncol, L, covariates.time.cross = FALSE)
#'
#' # Example 3: 3 time points, 3 classes, time-invariant coefficients
#' covariates.ncol <- c(2, 2, 2)  # All time points have identical covariates
#' L <- 3
#' get.npar.LTA(covariates.ncol, L, covariates.time.cross = TRUE)
#'
#' # Example 4: 3 time points, 3 classes, time-varying coefficients
#' covariates.ncol <- c(2, 3, 4)
#' L <- 3
#' get.npar.LTA(covariates.ncol, L, covariates.time.cross = FALSE)
#'
#' # Example 5: Single time point (equivalent to LCA)
#' covariates.ncol <- c(3)
#' L <- 4
#' get.npar.LTA(covariates.ncol, L)
#'
#' @export
get.npar.LTA <- function(covariates.ncol, L, covariates.time.cross = FALSE) {
  times <- length(covariates.ncol)
  n_beta <- covariates.ncol[1] * (L - 1)
  n_gama <- 0

  if (times > 1) {
    if (covariates.time.cross) {
      M_transfer <- covariates.ncol[2]
      n_gama <- L * (L - 1) * M_transfer
    } else {
      for (t in 2:times) {
        n_gama <- n_gama + L * (L - 1) * covariates.ncol[t]
      }
    }
  }

  total_npar <- n_beta + n_gama
  return(total_npar)
}
