#' Calculate Number of Free Parameters in Latent Transition Analysis
#'
#' Computes the total number of free parameters in a Latent Transition Analysis (LTA) model
#' estimated via the three-step approach. The count depends on the number of latent classes,
#' the number of time points, the number of covariates at each time point, and whether
#' transition coefficients are constrained to be equal across time.
#'
#' @param covariates.ncol An integer vector of length \eqn{T} (number of time points).
#'   Each element \eqn{M_t} represents the number of covariates (columns) for time point \eqn{t}.
#'   Must include an intercept column (all \code{1}s) as the first covariate.
#' @param L Integer scalar. Number of latent classes (\eqn{L \geq 2}).
#' @param covariates.timeCross Logical. If \code{TRUE}, transition coefficients are
#'   constrained to be identical across all transitions (time-invariant effects). This requires
#'   that the number of covariates is the same for all time points after the first (i.e., \eqn{M_2 = M_3 = \dots = M_T}).
#'   If \code{FALSE} (default), each transition has its own set of coefficients.
#'
#' @return Integer representing the total number of free parameters:
#'   \deqn{npar = M_1 \times (L-1) + \begin{cases}
#'     L \times (L-1) \times M_2 & \text{if } T>1 \text{ and time-invariant effects} \\
#'     \sum_{t=2}^T L \times (L-1) \times M_t & \text{if } T>1 \text{ and time-varying effects} \\
#'     0 & \text{if } T=1
#'   \end{cases}}
#'   where:
#'   \itemize{
#'     \item \emph{time-invariant effects} corresponds to \code{covariates.timeCross = TRUE}
#'     \item \emph{time-varying effects} corresponds to \code{covariates.timeCross = FALSE}
#'   }
#'
#' @details Parameterization:
#'   \describe{
#'     \item{Initial status model (time 1):}{
#'       Multinomial logit model with \eqn{L} classes (last class is reference).
#'       Number of free parameters: \eqn{M_1 \times (L-1)}.
#'     }
#'     \item{Transition models (time \eqn{t \to t+1}):}{
#'       For each transition, a multinomial logit model conditioned on previous class.
#'       For each origin class \eqn{k} and destination class \eqn{l} (\eqn{l \neq L}),
#'       there is a coefficient vector of length \eqn{M_{t+1}}.
#'       Total per transition: \eqn{L \times (L-1) \times M_{t+1}} parameters.
#'       The constraint \code{covariates.timeCross} determines whether these parameters
#'       are shared across transitions.
#'     }
#'   }
#'
#' @note Critical assumptions:
#'   \itemize{
#'     \item The last latent class (\eqn{L}) is always the reference category for all logits.
#'     \item When \code{covariates.timeCross = TRUE}, it is assumed that all time points after the first
#'           have identical covariate structures (\eqn{M_2 = M_3 = \dots = M_T}). If violated, the function
#'           uses \eqn{M_1} for all transitions to match \code{\link[LCPA]{LTA}}'s internal behavior.
#'     \item For \eqn{T=1}, no transition parameters are estimated (pure latent class/profile analysis).
#'   }
#'
#' @examples
#' # Example 1: 2 time points, 2 classes, time-invariant transition coefficients
#' #   Time1: 2 covariates (intercept + 1 predictor)
#' #   Time2: 3 covariates (but constrained to match Time1 due to timeCross=TRUE)
#' covariates.ncol <- c(2, 3)
#' L <- 2
#' get.npar.LTA(covariates.ncol, L, covariates.timeCross = TRUE)
#'
#' # Example 2: Same as above but time-varying coefficients
#' get.npar.LTA(covariates.ncol, L, covariates.timeCross = FALSE)
#'
#' # Example 3: 3 time points, 3 classes, time-invariant coefficients
#' covariates.ncol <- c(2, 2, 2)  # All time points have identical covariates
#' L <- 3
#' get.npar.LTA(covariates.ncol, L, covariates.timeCross = TRUE)
#'
#' # Example 4: 3 time points, 3 classes, time-varying coefficients
#' covariates.ncol <- c(2, 3, 4)
#' L <- 3
#' get.npar.LTA(covariates.ncol, L, covariates.timeCross = FALSE)
#'
#' # Example 5: Single time point (equivalent to LCA)
#' covariates.ncol <- c(3)
#' L <- 4
#' get.npar.LTA(covariates.ncol, L)
#'
#' @export
get.npar.LTA <- function(covariates.ncol, L, covariates.timeCross = FALSE) {
  times <- length(covariates.ncol)
  n_beta <- covariates.ncol[1] * (L - 1)
  n_gama <- 0

  if (times > 1) {
    if (covariates.timeCross) {
      M_transfer <- covariates.ncol[1]
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
