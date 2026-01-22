#' Calculate Number of Free Parameters in Latent Profile Analysis
#'
#' Computes the total number of free parameters in an LPA model based on the number of observed
#' variables (\eqn{I}), number of latent profiles (\eqn{L}), and covariance structure constraints.
#'
#' @param I Integer specifying the number of continuous observed variables.
#' @param L Integer specifying the number of latent profiles.
#' @param constraint Character string specifying covariance structure constraints. Supported options:
#'   \describe{
#'     \item{Univariate case (\eqn{I = 1}):}{
#'       \describe{
#'         \item{\code{"UE"}}{Equal variance across all profiles (1 shared variance parameter).}
#'         \item{\code{"UV"}}{Varying variances across profiles (\eqn{L} profile-specific variance parameters).}
#'       }
#'     }
#'     \item{Multivariate case (\eqn{I > 1}):}{
#'       \describe{
#'         \item{\code{"E0"}}{Equal variances across profiles, zero covariances. Requires \eqn{I} variance parameters.}
#'         \item{\code{"V0"}}{Varying variances across profiles, zero covariances. Requires \eqn{L \times I} variance parameters.}
#'         \item{\code{"EE"}}{Equal variances and equal covariances across profiles (homogeneous covariance matrix). Requires \eqn{\frac{I(I+1)}{2}} parameters.}
#'         \item{\code{"VE"}}{Varying variances per profile, but \emph{equal covariances} across profiles. Requires \eqn{L \times I + \frac{I(I-1)}{2}} parameters.}
#'         \item{\code{"EV"}}{Equal variances across profiles, but \emph{varying covariances} per profile. Requires \eqn{I + L \times \frac{I(I-1)}{2}} parameters.}
#'         \item{\code{"VV"}}{Varying variances and varying covariances across profiles (heterogeneous covariance matrices). Requires \eqn{L \times \frac{I(I+1)}{2}} parameters.}
#'         \item{\code{list}}{Custom constraints. Each element is a 2-element integer vector specifying variables whose covariance
#'                            parameters are constrained equal across all classes. The constraint applies to:
#'                            \itemize{
#'                              \item Variances: When both indices are identical (e.g., \code{c(3,3)} forces variance of variable 3 to be equal across classes)
#'                              \item Covariances: When indices differ (e.g., \code{c(1,2)} forces covariance between variables 1 and 2 to be equal across classes)
#'                            }
#'                            Constraints are symmetric (e.g., \code{c(1,2)} automatically constrains \code{c(2,1)}). All unconstrained parameters
#'                            vary freely across classes while maintaining positive definiteness.
#'                            }
#'       }
#'     }
#'   }
#'   Default: \code{"VV"}.
#'
#' @return Integer representing the total number of free parameters in the model:
#'   \deqn{\text{npar} = \underbrace{L \times I}_{\text{means}} + \underbrace{(L-1)}_{\text{class proportions}} + \underbrace{\text{covariance parameters}}_{\text{depends on constraint}}}
#'
#' @details Parameter count breakdown:
#'   \enumerate{
#'     \item Fixed components (always present):
#'       \itemize{
#'         \item Profile-specific means: \eqn{L \times I} parameters
#'         \item Independent class proportions: \eqn{L-1} parameters (since \eqn{\sum_{l=1}^L \pi_l = 1})
#'       }
#'     \item Covariance parameters (varies by constraint):
#'       \describe{
#'         \item{\eqn{I = 1}:}{
#'           \itemize{
#'             \item \code{"UE"}: 1 shared variance parameter
#'             \item \code{"UV"}: \eqn{L} profile-specific variance parameters
#'           }
#'         }
#'         \item{\eqn{I > 1}:}{
#'           \itemize{
#'             \item \code{"E0"}: \eqn{I} shared variance parameters (no covariances)
#'             \item \code{"V0"}: \eqn{L \times I} profile-specific variance parameters (no covariances)
#'             \item \code{"EE"}: \eqn{\frac{I(I+1)}{2}} parameters for one shared full covariance matrix
#'             \item \code{"VE"}: \eqn{L \times I} diagonal variances (free per profile) + \eqn{\frac{I(I-1)}{2}} off-diagonal covariances (shared across profiles)
#'             \item \code{"EV"}: \eqn{I} diagonal variances (shared across profiles) + \eqn{L \times \frac{I(I-1)}{2}} off-diagonal covariances (free per profile)
#'             \item \code{"VV"}: \eqn{L \times \frac{I(I+1)}{2}} parameters for \eqn{L} distinct full covariance matrices
#'           }
#'         }
#'       }
#'   }
#'
#' @note Important considerations:
#'   \itemize{
#'     \item For \eqn{I = 1}, only \code{"UE"} and \code{"UV"} are meaningful; \code{"EE"}, \code{"E0"}, \code{"VV"}, \code{"V0"}, etc., are treated as \code{"UE"} or \code{"UV"} respectively.
#'     \item Covariance parameters count only free elements in symmetric matrices (diagonal + upper triangle).
#'     \item If an user-defined \code{constraint} is provided, the function defaults to \code{"VV"}
#'           behavior but subtracts \eqn{(L-1) \times \text{length(constraint)}}.
#'     \item \code{"VE"} and \code{"EV"} constraints require \eqn{I > 1} to be meaningful (otherwise no covariances exist).
#'   }
#'
#' @examples
#' # Univariate examples (I=1)
#' get.npar.LPA(I = 1, L = 2, constraint = "UE")
#' get.npar.LPA(I = 1, L = 3, constraint = "UV")
#'
#' # Multivariate examples (I=3)
#' get.npar.LPA(I = 3, L = 2, constraint = "E0")
#' get.npar.LPA(I = 3, L = 2, constraint = "V0")
#' get.npar.LPA(I = 3, L = 2, constraint = "EE")
#' get.npar.LPA(I = 3, L = 2, constraint = "VV")
#' get.npar.LPA(I = 3, L = 2, constraint = "VE")
#' get.npar.LPA(I = 3, L = 2, constraint = "EV")
#'
#' # User defined example
#' get.npar.LPA(I = 3, L = 2, constraint = list(c(1, 2), c(3, 3)))
#'
#' @export
get.npar.LPA <- function(I, L, constraint = "VV") {
  npar <- L * I + (L - 1)

  if (I == 1 && any(constraint %in% c("UE", "UV"))) {
    if (constraint == "UE") {
      npar <- npar + 1
    } else if (constraint == "UV") {
      npar <- npar + L
    }
  }
  else if (any(constraint %in% c("E0", "V0", "EE", "VV", "VE", "EV"))) {
    if (constraint == "E0") {
      npar <- npar + I
    } else if (constraint == "V0") {
      npar <- npar + L * I
    } else if (constraint == "EE") {
      npar <- npar + I * (I + 1) / 2
    } else if (constraint == "VV") {
      npar <- npar + L * I * (I + 1) / 2
    }else if (constraint == "VE") {
      npar <- npar + L * I + I * (I - 1) / 2
    }else if (constraint == "EV") {
      npar <- npar + I + L * I * (I - 1) / 2
    }
  }else{
    npar <- npar + L * I * (I + 1) / 2
    npar = npar - (L-1) * length(constraint)
  }

  return(npar)
}
