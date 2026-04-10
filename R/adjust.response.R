#' Adjust Categorical Response Data for Polytomous Indicators
#'
#' Standardizes polytomous response data by converting raw category values to
#' consecutive integers starting from 0. Records original category values for
#' potential reverse transformation. Handles varying numbers of response
#' categories across indicators.
#'
#' @param response A matrix or data frame containing response data where:
#'   \itemize{
#'     \item Rows represent respondents (\eqn{N} observations)
#'     \item Columns represent indicators/items/questions (\eqn{I} indicators)
#'     \item Cells contain raw response values (numeric)
#'   }
#'   Non-numeric columns will be coerced to numeric with warning.
#'
#' @return A named list containing:
#'   \describe{
#'     \item{\code{poly.orig}}{\eqn{I \times K_{max}} matrix. Original sorted category values for each indicator.
#'       Rows correspond to indicators, columns to category positions. Empty cells filled with \code{NA}.}
#'     \item{\code{poly.value}}{Integer vector of length \eqn{I}. Number of unique response categories per indicator.}
#'     \item{\code{poly.max}}{Scalar integer. Maximum number of categories across all indicators, i.e., \eqn{K_{max}}.}
#'     \item{\code{response}}{\eqn{N \times I} matrix. Adjusted response data where original values are replaced by
#'                            zero-based category indices (0 to \eqn{k-1} for \eqn{k} categories).}
#'   }
#'
#' @details The function processes each indicator column independently:
#'   \enumerate{
#'     \item Extracts unique response values and sorts them in ascending order
#'     \item Maps smallest value to 0, second smallest to 1, etc.
#'     \item Records original values in \code{poly.orig} for possible reverse transformation
#'     \item Handles indicators with different numbers of categories through NA-padding
#'   }
#'   Missing values (\code{NA}) in input are preserved as \code{NA} in output.
#'
#' @examples
#' # Simulate response data with 3 indicators and varying categories
#' set.seed(123)
#' resp <- data.frame(
#'   indicator1 = sample(1:3, 10, replace = TRUE),
#'   indicator2 = sample(c(0, 5, 10), 10, replace = TRUE),
#'   indicator3 = sample(1:2, 10, replace = TRUE)
#' )
#'
#' # Apply adjustment
#' adjusted <- adjust.response(resp)
#'
#' # Inspect results
#' str(adjusted)
#' print(adjusted$poly.orig)  # Original category values
#' print(adjusted$response)   # Standardized responses
#'
#' @export
adjust.response <- function(response) {
  if (!is.matrix(response)) {
    response <- as.matrix(response)
  }
  result <- adjust_response_cpp(response)
  return(result)
}
