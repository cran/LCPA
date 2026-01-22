#' Validate response matrix against expected polytomous category counts
#'
#' Checks whether each column in the response matrix contains exactly the number
#' of unique response categories specified in \code{poly.value}. Handles edge cases
#' where all items have identical category counts efficiently.
#'
#' @param response A numeric matrix of dimension \eqn{N \times I}, where:
#'   \itemize{
#'     \item \eqn{N}: Number of subjects/observations (rows)
#'     \item \eqn{I}: Number of items/variables (columns)
#'   }
#'   Each cell contains the observed response value for a subject on an item.
#' @param poly.value An integer vector of length \eqn{I} specifying the expected number
#'   of unique response categories (levels) for each corresponding item in
#'   \code{response}. Values must be positive integers.
#'
#' @return Logical value indicating validation status:
#'   \itemize{
#'     \item \code{TRUE} if either:
#'       \itemize{
#'         \item All columns have identical numbers of unique values (regardless of \code{poly.value} specification)
#'         \item Each column's unique value count matches its corresponding \code{poly.value} entry
#'       }
#'     \item \code{FALSE} if any column's unique value count mismatches its specified \code{poly.value}
#'       (when columns have varying category counts)
#'   }
#'
#' @note This function contains a specific behavior: When all items have identical numbers of
#'   unique response categories, it returns \code{TRUE} immediately without validating against
#'   \code{poly.value}. This may lead to unexpected results if \code{poly.value} contains
#'   inconsistent expectations. Users should ensure \code{poly.value} accurately reflects
#'   their measurement model.
#'
#' @examples
#' # Valid case: Matching category counts
#' resp_matrix <- matrix(c(1,1,2,2, 1,2,3,1), ncol = 2)
#' check.response(resp_matrix, poly.value = c(2, 3))  # Returns TRUE
#'
#' # Invalid case: Mismatched category counts
#' check.response(resp_matrix, poly.value = c(2, 2))  # Returns FALSE
#'
#' # Special case: Uniform category counts bypass poly.value check
#' uniform_resp <- matrix(rep(1:2, each = 4), ncol = 2)
#' check.response(uniform_resp, poly.value = c(2, 5))  # Returns TRUE (bypass behavior)
#'
#' @export
check.response <- function(response, poly.value) {
  poly.value.response <- apply(response, 2, unique)
  if(any(class(poly.value.response) == "matrix")){
    return(TRUE)
  } else {
    for(i in 1:ncol(response)){
      if(length(poly.value.response[[i]]) != poly.value[i]){
        return(FALSE)
      }
    }
  }
  return(TRUE)
}
