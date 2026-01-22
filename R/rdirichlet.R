#' Generate Random Samples from the Dirichlet Distribution
#'
#' @description
#' \code{rdirichlet} generates \code{n} random observations from a Dirichlet distribution
#' with a specified concentration parameter vector \code{alpha}.
#'
#' @details
#' The Dirichlet distribution is a family of continuous multivariate probability distributions
#' parameterized by a vector \eqn{\alpha} of positive reals. It is the multivariate
#' generalization of the beta distribution and is commonly used as a conjugate prior
#' to the multinomial distribution in Bayesian statistics.
#'
#' \bold{Probability Density Function:}
#'
#' For a vector \eqn{x = (x_1, \dots, x_K)} on the unit simplex (where \eqn{\sum x_i = 1}
#' and \eqn{x_i \ge 0}), the density is given by:
#'
#' \deqn{f(x_1, \dots, x_K; \alpha_1, \dots, \alpha_K) = \frac{1}{B(\alpha)} \prod_{i=1}^{K} x_i^{\alpha_i - 1}}
#'
#' where the normalizing constant \eqn{B(\alpha)} is the multivariate beta function:
#'
#' \deqn{B(\alpha) = \frac{\prod_{i=1}^{K} \Gamma(\alpha_i)}{\Gamma(\sum_{i=1}^{K} \alpha_i)}}
#'
#' \bold{Simulation Method:}
#'
#' The function utilizes the property that if \eqn{Y_1, \dots, Y_K} are independent
#' Gamma random variables such that \eqn{Y_i \sim Gamma(shape = \alpha_i, rate = 1)}, then:
#'
#' \deqn{X_i = \frac{Y_i}{\sum_{j=1}^{K} Y_j}}
#'
#' The resulting vector \eqn{(X_1, \dots, X_K)} follows a Dirichlet distribution with parameters \eqn{\alpha}.
#'
#' @param n Integer. The number of random vectors to generate.
#' @param alpha Numeric vector. The concentration parameters (must be positive).
#' The length of this vector determines the number of dimensions \eqn{K}.
#'
#' @return A matrix with \code{n} rows and \code{length(alpha)} columns.
#' Each row sums to 1, representing a single sample from the Dirichlet distribution.
#'
#' @examples
#' # Generate 5 samples from a 3-dimensional Dirichlet distribution
#' set.seed(123)
#' alpha_params <- c(1, 2, 5)
#' result <- rdirichlet(n = 5, alpha = alpha_params)
#' print(result)
#'
#' # Check that rows sum to 1
#' rowSums(result)
#'
#' @importFrom stats rgamma
#' @export
rdirichlet <- function(n, alpha) {
  l <- length(alpha)
  # Generate Gamma distributed variables
  x <- matrix(rgamma(l * n, shape = alpha, rate = 1), ncol = l, byrow = TRUE)
  # Normalize each row to sum to 1 (Project onto the simplex)
  return(sweep(x, 1, rowSums(x), "/"))
}
