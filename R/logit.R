#' Compute the Logistic (Sigmoid) Function
#'
#' This function computes the logistic (also known as sigmoid) transformation of the input.
#' The logistic function maps real-valued numbers to the open interval (0, 1), and is widely
#' used in machine learning, statistical modeling (e.g., logistic regression), and neural networks
#' as an activation function or link function.
#'
#' The logistic function is defined as:
#' \deqn{ \mathrm{logit}^{-1}(x) = \frac{1}{1 + e^{-x}} }
#'
#' Note: Despite the name "logit", this function actually computes the *inverse logit* (i.e., the
#' logistic function). The true logit function is the inverse: \eqn{\log(p / (1 - p))}.
#' However, in many applied contexts—especially in software—the term "logit" is sometimes
#' informally used to refer to the sigmoid. For clarity, this implementation follows the
#' conventional definition of the logistic/sigmoid function.
#'
#' @param x A numeric vector, matrix, or array. Accepts any real number, including \code{Inf}
#'   and \code{-Inf}. Missing values (\code{NA}) are preserved.
#'
#' @return A numeric object of the same dimension as \code{x}, where each element is the
#'   logistic transformation of the corresponding input:
#'   \itemize{
#'     \item If \code{x = 0}, returns \code{0.5}
#'     \item As \code{x -> Inf}, output approaches \code{1}
#'     \item As \code{x -> -Inf}, output approaches \code{0}
#'     \item \code{NA} values remain \code{NA}
#'   }
#'
#' @examples
#' logit(0)        # 0.5
#' logit(c(-Inf, 0, Inf))  # c(0, 0.5, 1)
#' logit(c(-2, -1, 0, 1, 2))
#'
#' @export
logit <- function(x) {
  # Input validation: ensure x is numeric
  if (!is.numeric(x)) {
    stop("Input 'x' must be numeric.")
  }
  out <- ifelse(
    x >= 0,
    1 / (1 + exp(-x)),
    exp(x) / (1 + exp(x))
  )

  # Preserve NA values
  out[is.na(x)] <- NA

  return(out)
}
