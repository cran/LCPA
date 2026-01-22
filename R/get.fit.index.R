#' Calculate Fit Indices
#'
#' Computes a comprehensive set of model fit indices for objects returned by
#' \code{\link[LCPA]{LCA}} or \code{\link[LCPA]{LPA}}. These indices balance model
#' fit (log-likelihood) with model complexity (number of parameters) to facilitate
#' model selection. All indices are derived from the observed-data log-likelihood
#' and parameter count.
#'
#' @param object An object of class \code{"LCA"} or \code{"LPA"} returned by
#'   \code{\link[LCPA]{LCA}}, \code{\link[LCPA]{LPA}} or any object containing:
#'   \itemize{
#'     \item \code{Log.Lik}: Log-likelihood value
#'     \item \code{npar}: Number of free parameters
#'     \item \eqn{N} = Total number of observations (\eqn{n = 1, 2, \dots, N})
#'   }
#'
#' @return An object of class \code{"fit.index"} containing:
#'   \describe{
#'     \item{npar}{Number of free parameters in the model}
#'     \item{Log.Lik}{Log-likelihood of the model: \eqn{\log \mathcal{L}}}
#'     \item{-2LL}{Deviance statistic: \eqn{-2 \log \mathcal{L}}
#'       \deqn{-2 \sum_{n=1}^N \log \left[ \sum_{l=1}^L \pi_l \cdot f(\mathbf{x}_n \mid \boldsymbol{\theta}_l) \right]},
#'       where \eqn{\pi_l} is the prior probability of class \eqn{l},
#'       \eqn{f(\cdot)} is the probability density/mass function (multivariate normal for LPA,
#'       multinomial for LCA), and \eqn{\boldsymbol{\theta}_l} are class-specific parameters.}
#'     \item{AIC}{Akaike Information Criterion:
#'       \eqn{\mathrm{AIC} = -2 \log \mathcal{L} + 2k},
#'       where \eqn{npar} = number of free parameters. Lower values indicate better fit.}
#'     \item{BIC}{Bayesian Information Criterion:
#'       \eqn{\mathrm{BIC} = -2 \log \mathcal{L} + npar \log(N)},
#'       where \eqn{N} = sample size. Incorporates stronger penalty for complexity than AIC.}
#'     \item{SIC}{Sample-Size Adjusted BIC:
#'       \eqn{\mathrm{SIC} = -\frac{1}{2} \mathrm{BIC}}.
#'       Equivalent to \eqn{\log \mathcal{L} - \frac{npar}{2} \log(N)}. Often used in latent class modeling.}
#'     \item{CAIC}{Consistent AIC:
#'       \eqn{\mathrm{CAIC} = -2 \log \mathcal{L} + npar \left[ \log(N) + 1 \right]}.
#'       Consistent version of AIC that converges to true model as \eqn{N \to \infty}.}
#'     \item{AWE}{Approximate Weight of Evidence:
#'       \eqn{\mathrm{AWE} = -2 \log \mathcal{L} + 1.5k \left[ \log(N) + 1 \right]}.
#'       Penalizes complexity more heavily than CAIC.}
#'     \item{SABIC}{Sample-Size Adjusted BIC:
#'       \eqn{\mathrm{SABIC} = -2 \log \mathcal{L} + npar \log \left( \frac{N + 2}{24} \right)}.
#'       Recommended for latent class/profile analysis with moderate sample sizes.}
#'   }
#'
#' @examples
#' \donttest{
#' # Fit LPA model
#' set.seed(123)
#' data.obj <- sim.LPA(N = 100, I = 3, L = 2, constraint = "E0")
#' fit <- LPA(data.obj$response, L = 2, constraint = "VV", method = "EM")
#'
#' # Compute fit indices
#' fit_indices <- get.fit.index(fit)
#'
#' fit_indices
#'
#' extract(fit_indices, "SABIC")
#' }
#' @export
get.fit.index <- function(object){

  call <- match.call()

  Log.Lik <- object$Log.Lik
  npar <- object$npar

  if(!is.null(object$arguments$response)){
    N <- nrow(object$arguments$response)
  }

  if(!is.null(object$N)){
    N <- object$N
  }

  `-2LL` <- -2 * Log.Lik

  AIC <- -2 * Log.Lik + 2 * npar

  BIC <- -2 * Log.Lik + log(N) * npar

  CAIC <- -2 * Log.Lik + npar * (log(N)+1)

  SABIC <- -2 * Log.Lik + npar * log((N+2)/24)

  AWE <- -2 * Log.Lik + npar * (log(N)+1.5)

  SIC = -0.5 * BIC

  res <- list(npar = npar, Log.Lik = Log.Lik,
              `-2LL`=`-2LL`, AIC=AIC, BIC=BIC, SIC=SIC, CAIC=CAIC,
              AWE=AWE, SABIC=SABIC,
              arguments = list(
                object=object
              ))
  res$call <- call

  class(res) <- "fit.index"

  return(res)
}
