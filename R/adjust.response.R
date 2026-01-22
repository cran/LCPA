#' Adjust Categorical Response Data for Polytomous Items
#'
#' Standardizes polytomous response data by converting raw category values to
#' consecutive integers starting from 0. Records original category values for
#' potential reverse transformation. Handles varying numbers of response
#' categories across items.
#'
#' @param response A matrix or data frame containing response data where:
#'   \itemize{
#'     \item Rows represent respondents (\eqn{N} observations)
#'     \item Columns represent items/questions (\eqn{I} items)
#'     \item Cells contain raw response values (numeric)
#'   }
#'   Non-numeric columns will be coerced to numeric with warning.
#'
#' @return A named list containing:
#'   \describe{
#'     \item{\code{poly.orig}}{\eqn{I \times K_{max}} matrix. Original sorted category values for each item.
#'       Rows correspond to items, columns to category positions. Empty cells filled with \code{NA}.}
#'     \item{\code{poly.value}}{Integer vector of length \eqn{I}. Number of unique response categories per item.}
#'     \item{\code{poly.max}}{Scalar integer. Maximum number of categories across all items, i.e., \eqn{K_{max}}.}
#'     \item{\code{response}}{\eqn{N \times I} matrix. Adjusted response data where original values are replaced by
#'                            zero-based category indices (0 to \eqn{k-1} for \eqn{k} categories).}
#'   }
#'
#' @details The function processes each item column independently:
#'   \enumerate{
#'     \item Extracts unique response values and sorts them in ascending order
#'     \item Maps smallest value to 0, second smallest to 1, etc.
#'     \item Records original values in \code{poly.orig} for possible reverse transformation
#'     \item Handles items with different numbers of categories through NA-padding
#'   }
#'   Missing values (\code{NA}) in input are preserved as \code{NA} in output.
#'
#' @examples
#' # Simulate response data with 3 items and varying categories
#' set.seed(123)
#' resp <- data.frame(
#'   item1 = sample(1:3, 10, replace = TRUE),
#'   item2 = sample(c(0, 5, 10), 10, replace = TRUE),
#'   item3 = sample(1:2, 10, replace = TRUE)
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
