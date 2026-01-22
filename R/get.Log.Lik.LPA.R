#' Calculate Log-Likelihood for Latent Profile Analysis
#'
#' Computes the log-likelihood of observed continuous data under a Latent Profile Analysis (LPA) model
#' with multivariate normal distributions within each latent profile. Implements robust numerical
#' techniques to handle near-singular covariance matrices.
#'
#' @param response A numeric matrix of dimension \eqn{N \times I} containing continuous observations.
#'   Rows represent observations, columns represent variables. Missing values are not permitted.
#' @param P.Z A numeric vector of length \eqn{L} containing prior probabilities for latent profiles.
#'   Must satisfy:
#'   \itemize{
#'     \item \eqn{\sum_{l=1}^L \pi_l = 1}
#'     \item \eqn{\pi_l > 0} for all \eqn{l = 1, \dots, L}
#'   }
#' @param means A matrix of dimension \eqn{L \times I} where row \eqn{l} contains the mean vector
#'   \eqn{\boldsymbol{\mu}_l} for profile \eqn{l}.
#' @param covs An array of dimension \eqn{I \times I \times L} where slice \eqn{l} contains the
#'   covariance matrix \eqn{\boldsymbol{\Sigma}_l} for profile \eqn{l}. Must be symmetric positive semi-definite.
#' @param jitter A small positive constant (default: 1e-10) added to diagonal elements of
#'   covariance matrices to ensure numerical stability during Cholesky decomposition.
#'
#' @return A single numeric value representing the total log-likelihood:
#'   \deqn{\log \mathcal{L} = \sum_{n=1}^N \log \left[ \sum_{l=1}^L \pi_l \cdot \mathcal{N}(\mathbf{x}_n \mid \boldsymbol{\mu}_l, \boldsymbol{\Sigma}_l) \right]}
#'   where \eqn{\mathcal{N}(\cdot)} denotes the multivariate normal density function.
#'
#' @details The log-likelihood calculation follows these steps:
#'
#' \itemize{
#'   \item Covariance Stabilization:
#'   Each covariance matrix \eqn{\boldsymbol{\Sigma}_l} is symmetrized as \eqn{(\boldsymbol{\Sigma}_l + \boldsymbol{\Sigma}_l^\top)/2}.
#'   If Cholesky decomposition fails:
#'   \itemize{
#'     \item Add \code{jitter} to diagonal elements iteratively (up to 10 attempts, scaling jitter by 10x each attempt)
#'     \item Fall back to diagonal covariance matrix if decomposition still fails
#'   }
#'
#'   \item Profile-Specific Density for observation \eqn{n} in profile \eqn{l}:
#'   \deqn{\log f(\mathbf{x}_n \mid Z_n=l) = -\frac{I}{2}\log(2\pi) - \frac{1}{2}\log|\boldsymbol{\Sigma}_l| - \frac{1}{2}(\mathbf{x}_n - \boldsymbol{\mu}_l)^\top \boldsymbol{\Sigma}_l^{-1} (\mathbf{x}_n - \boldsymbol{\mu}_l)}
#'   Computed efficiently using Cholesky decomposition \eqn{\boldsymbol{\Sigma}_l = \mathbf{R}^\top\mathbf{R}} where applicable.
#'
#'   \item Joint Probability for observation \eqn{n} and profile \eqn{l}:
#'   \deqn{\log[\pi_l \cdot f(\mathbf{x}_n \mid Z_n=l)] = \log(\pi_l) + \log f(\mathbf{x}_n \mid Z_n=l)}
#'   \eqn{\log(\pi_l)} uses \eqn{\log(\pi_l + 10^{-12})} to avoid undefined values.
#'
#'   \item Marginal Likelihood per observation using log-sum-exp trick for numerical stability:
#'   \deqn{\log f(\mathbf{x}_n) = a_{\max} + \log\left( \sum_{l=1}^L \exp\left\{ \log[\pi_l \cdot f(\mathbf{x}_n \mid Z_n=l)] - a_{\max} \right\} \right)}
#'   where \eqn{a_{\max} = \max_l \log[\pi_l \cdot f(\mathbf{x}_n \mid Z_n=l)]}.
#'
#'   \item Total Log-Likelihood: Sum of \eqn{\log f(\mathbf{x}_n)} across all observations \eqn{n=1,\dots,N}.
#' }
#'
#' @note Critical implementation details:
#'   \itemize{
#'     \item Cholesky Decomposition: For non-degenerate cases (\eqn{I>1}), used to compute:
#'       \eqn{\log|\boldsymbol{\Sigma}_l| = 2\sum_{i=1}^I \log(r_{ii})} and quadratic form \eqn{\|\mathbf{R}^{-\top}(\mathbf{x}_n - \boldsymbol{\mu}_l)\|^2}
#'     \item Univariate Handling: When \eqn{I=1}, computes density directly without decomposition
#'     \item Numerical Safeguards:
#'       \itemize{
#'         \item Densities clamped to \eqn{-10^{10}} when non-finite
#'         \item Marginal likelihoods clamped to \eqn{-10^{10}} when non-finite
#'         \item Explicit dimension checks for \code{means} and \code{covs}
#'       }
#'     \item Assumptions:
#'       \itemize{
#'         \item Multivariate normality within profiles
#'         \item No missing data in \code{response}
#'         \item Positive-definite covariance matrices (after stabilization)
#'       }
#'   }
#'
#' @export
get.Log.Lik.LPA <- function(response, P.Z, means, covs, jitter = 1e-10) {
  if (!is.matrix(response)) stop("response must be a matrix")
  if (!is.numeric(P.Z) || any(P.Z <= 0) || abs(sum(P.Z) - 1) > 1e-10) {
    stop("P.Z must be a valid probability vector (positive values summing to 1)")
  }
  if (!is.matrix(means) || nrow(means) != length(P.Z)) {
    stop("means must be a matrix with rows corresponding to classes")
  }
  I <- ncol(response)
  L <- length(P.Z)

  if (!is.array(covs) || length(dim(covs)) != 3 ||
      dim(covs)[1] != I || dim(covs)[2] != I || dim(covs)[3] != L) {
    stop("covs must be an I x I x L array matching the dimensions of response and classes")
  }

  N <- nrow(response)
  tresponse <- t(response)
  const.term <- -0.5 * I * log(2 * pi)
  log_densities <- matrix(NA_real_, nrow = N, ncol = L)

  for (l in 1:L) {
    covs.l <- covs[,,l]
    mean.l <- means[l,]
    covs.l <- (covs.l + t(covs.l)) / 2
    if (I == 1) {
      var_val <- covs.l[1, 1]
      var_val <- max(var_val, jitter)
      dev <- tresponse - mean.l
      quad <- (dev^2) / var_val
      logdet <- log(var_val)
      lp <- const.term - 0.5 * (logdet + quad)

    } else {
      R <- tryCatch(chol(covs.l), error = function(e) NULL)

      if (is.null(R)) {
        diag(covs.l) <- pmax(diag(covs.l), jitter)
        add <- jitter
        tries <- 0
        while (is.null(R) && tries < 10) {
          diag(covs.l) <- diag(covs.l) + add
          R <- tryCatch(chol(covs.l), error = function(e) NULL)
          add <- add * 10
          tries <- tries + 1
        }
        if (is.null(R)) {
          covs.l <- diag(pmax(diag(covs.l), jitter), I)
          R <- chol(covs.l)
        }
      }

      dev <- tresponse - matrix(mean.l, nrow = I, ncol = N, byrow = FALSE)

      sol <- tryCatch({
        backsolve(R, dev, transpose = TRUE)
      }, error = function(e) {
        matrix(NA_real_, nrow = I, ncol = N)
      })

      quad <- colSums(sol^2, na.rm = FALSE)
      logdet <- 2 * sum(log(diag(R)))
      lp <- const.term - 0.5 * (logdet + quad)
    }

    invalid <- !is.finite(lp) | is.na(lp)
    if (any(invalid)) lp[invalid] <- -1e10

    log_densities[, l] <- lp
  }

  log_prior <- log(P.Z + 1e-10)  # 避免log(0)
  log_joint <- sweep(log_densities, 2, log_prior, "+")  # log(P.Z[l] * density)

  row_max <- apply(log_joint, 1, max)
  log_sum_exp <- row_max + log(rowSums(exp(log_joint - matrix(row_max, N, L))))

  invalid_rows <- !is.finite(log_sum_exp) | is.na(log_sum_exp)
  if (any(invalid_rows)) {
    log_sum_exp[invalid_rows] <- -1e10
  }

  total_log_lik <- sum(log_sum_exp)

  return(total_log_lik)
}
