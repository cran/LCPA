#' Compute Posterior Latent Class Probabilities Based on Fixed Parameters
#'
#' Computes posterior probabilities of latent class membership for each observation
#' using fixed conditional response probabilities (\code{par}) and fixed class prior
#' probabilities (\code{P.Z}).
#'
#' @param response Numeric matrix (\eqn{N \times I}) of categorical responses.
#' @param par 3D array (\eqn{L \times I \times K_{\max}}) of fixed conditional response
#'   probabilities where:
#'   \itemize{
#'     \item \eqn{L} = number of latent classes
#'     \item \eqn{I} = number of indicators
#'     \item \eqn{K_{\max}} = maximum categories across indicators
#'   }
#'   \code{par[l, i, q]} = \eqn{P(X_i=q-1\mid Z=l)} (using 1-based array
#'   indexing for the dimension corresponding to category \eqn{q-1}).
#'
#' @param P.Z Vector of length \eqn{L} with fixed class prior probabilities (\eqn{\pi_l}).
#'   These values are used directly without re-estimation.
#' @param category.levels List of length \eqn{I}. Element \eqn{i} contains the
#'   ordered response categories for indicator \eqn{i} as fixed in Step 1.
#'   Categories absent from the current sample remain in this mapping.
#'
#' @return Numeric matrix (\eqn{N \times L}) of posterior probabilities.
#'   Rows sum to 1. Columns are named `"Class 1"`, `"Class 2"`, and so on.
#'
#' @details
#' Unlike an EM algorithm, this function does NOT iteratively update class prevalences.
#' It performs a single calculation step based on Bayes' theorem:
#' \deqn{
#'   \tau_{nl}=P(Z_n=l\mid\mathbf{X}_n)=
#'   \frac{\pi_l\prod_{i=1}^I
#'   P(X_{ni}=x_{ni}\mid Z_n=l)}
#'   {\sum_{h=1}^L\pi_h\prod_{i=1}^I
#'   P(X_{ni}=x_{ni}\mid Z_n=h)}
#' }
#' where the class probabilities and conditional response probabilities are
#' fixed by the \code{P.Z} and \code{par} arguments.
#'
#'
#' @examples
#' \donttest{
#' library(LCPA)
#' set.seed(123)
#' # Simulate data
#' data.obj <- sim.LCA(N = 200, I = 3, L = 2, IQ = 0.85)
#'
#' # Fit a model to get parameters
#' fit <- LCA(data.obj$response, L = 2, method = "EM", nrep = 5)
#'
#' # Calculate posteriors using fixed parameters from the fitted model
#' P.Z.Xn <- get.P.Z.Xn.LCA(
#'   response = data.obj$response,
#'   par = fit$params$par,
#'   P.Z = fit$params$P.Z,
#'   category.levels = fit$params$category.levels
#' )
#' head(P.Z.Xn)
#' }
#'
#' @export
get.P.Z.Xn.LCA <- function(response, par, P.Z, category.levels){
  Y <- .LCA.encode.response(response, category.levels)
  N <- nrow(Y)
  I <- ncol(Y)
  if (length(dim(par)) != 3L || dim(par)[2L] != I) {
    stop("par must be an L by I by K array matching response")
  }
  L <- dim(par)[1]
  if (length(P.Z) != L || any(!is.finite(P.Z)) || any(P.Z < 0) || sum(P.Z) <= 0) {
    stop("P.Z must contain L finite non-negative values with a positive sum")
  }
  if (any(Y + 1L > dim(par)[3L])) {
    stop("The Step 1 category mapping exceeds the category dimension of par")
  }

  expectation <- lca_expectation_cpp(Y, par, as.numeric(P.Z) / sum(P.Z))
  if(!isTRUE(expectation$valid)){
    stop("At least one response pattern has zero probability in every class")
  }
  P.Z.Xn <- expectation$posterior
  colnames(P.Z.Xn) <- .latent.group.names(L, "LCA")
  P.Z.Xn
}

.LCA.category.levels <- function(response) {
  response <- as.matrix(response)
  if (!is.numeric(response) || any(!is.finite(response))) {
    stop("LCA responses must be a finite numeric matrix")
  }
  lapply(seq_len(ncol(response)), function(i) {
    sort(unique(response[, i]))
  })
}

.LCA.encode.response <- function(response, category.levels) {
  response <- as.matrix(response)
  if (!is.numeric(response) || any(!is.finite(response))) {
    stop("LCA responses must be a finite numeric matrix")
  }
  I <- ncol(response)
  if (!is.list(category.levels) || length(category.levels) != I) {
    stop("category.levels must be a list with one element per indicator")
  }

  encoded <- matrix(NA_integer_, nrow(response), I, dimnames = dimnames(response))
  for (i in seq_len(I)) {
    levels.cur <- category.levels[[i]]
    if (!is.numeric(levels.cur) || length(levels.cur) < 1L ||
        any(!is.finite(levels.cur)) || anyDuplicated(levels.cur)) {
      stop("Every element of category.levels must contain unique finite numeric values")
    }
    index <- match(response[, i], levels.cur)
    if (anyNA(index)) {
      unknown <- unique(response[is.na(index), i])
      stop(
        "Indicator ", i, " contains categories not present in the Step 1 mapping: ",
        paste(unknown, collapse = ", ")
      )
    }
    encoded[, i] <- index - 1L
  }
  encoded
}
