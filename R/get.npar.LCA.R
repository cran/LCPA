#' Calculate Number of Free Parameters in Latent Class Analysis
#'
#' Computes the total number of free parameters in an LCA model based on the number of categories
#' per observed variable and the number of latent classes. This follows standard LCA parameterization
#' with local independence assumption.
#'
#' @param poly.value A numeric vector of length \eqn{I} where each element \eqn{K_i} represents
#'   the number of response categories for observed variable \eqn{i}.
#' @param L Integer specifying the number of latent classes.
#'
#' @return Integer representing the total number of free parameters in the model:
#'   \deqn{\text{npar} = \sum_{i=1}^I \underbrace{(L \times K_i - 1)}_{\text{free parameters}} + \underbrace{(L-1)}_{\text{class proportions}}}
#'
#' @details Parameter count derivation:
#'   \describe{
#'     \item{Fixed components (always present):}{
#'       \itemize{
#'         \item Conditional response probabilities: \eqn{\sum_{i=1}^I (L \times K_i - 1)} parameters
#'         \item Independent class proportions: \eqn{L-1} parameters (since \eqn{\sum_{l=1}^L \pi_l = 1})
#'       }
#'     }
#'     \item{Per-variable parameterization:}{
#'       For each observed variable \eqn{i} with \eqn{K_i} categories:
#'       \itemize{
#'         \item Each latent class requires \eqn{K_i} conditional probabilities \eqn{P(X_i=k|Z=l)}
#'         \item With constraints \eqn{\sum_{k=1}^{K_i} P(X_i=k|Z=l) = 1} for each class \eqn{l}
#'         \item Global constraints reduce total parameters to \eqn{L \times K_i - 1} per variable
#'       }
#'     }
#'   }
#'
#' @examples
#' # Example 1: 3 binary variables (K_i=2), 2 latent classes
#' poly.value <- c(2, 2, 2)  # Three binary variables
#' L <- 2
#' npar <- sum(poly.value * L - 1) + (L - 1)  # = (4-1)+(4-1)+(4-1) + 1 = 3+3+3+1 = 10
#' get.npar.LCA(poly.value, L)  # Returns 10
#'
#' # Example 2: Mixed variable types (binary, ternary, quaternary)
#' poly.value <- c(2, 3, 4)  # Variables with 2, 3, and 4 categories
#' L <- 3
#' npar <- sum(poly.value * L - 1) + (L - 1)  # = (6-1)+(9-1)+(12-1) + 2 = 5+8+11+2 = 26
#' get.npar.LCA(poly.value, L)  # Returns 26
#'
#' # Example 3: Single polytomous variable with 5 categories, 4 latent classes
#' poly.value <- 5
#' L <- 4
#' npar <- sum(poly.value * L - 1) + (L - 1)  # = (20-1) + 3 = 19+3 = 22
#' get.npar.LCA(poly.value, L)  # Returns 22
#'
#' @export
get.npar.LCA <- function(poly.value, L) {
  npar <- sum(poly.value * L - 1) + (L - 1)
  return(npar)
}
