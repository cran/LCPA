#' Likelihood Ratio Test
#'
#' Conducts a likelihood ratio test to compare the fit of two
#' models. The test evaluates whether a model with more parameters
#' provides a significantly better fit than a model with fewer parameters.
#'
#' @param object1 Fitted model object with fewer parameters (i.e., fewer \code{npar}, small model).
#' @param object2 Fitted model object with more parameters (i.e., more \code{npar}, large model).
#'
#' @return An object of class \code{"htest"} containing:
#' \itemize{
#'   \item \code{statistic}: VLMR adjusted test statistic
#'   \item \code{parameter}: Degrees of freedom (\eqn{df = npar_2 - npar_1})
#'   \item \code{p.value}: P-value from \eqn{\chi^2_df} distribution
#'   \item \code{method}: Name of the test
#'   \item \code{data.name}: Model comparison description
#' }
#'
#' @details
#' Note that since the small model may be nested within the large model, the result
#' of \code{\link[LCPA]{LRT.test}} may not be accurate and is provided for reference only.
#' More reliable conclusions should be based on a combination of fit indices (i.e., \code{\link[LCPA]{get.fit.index}}),
#' classification accuracy measures (i.e., \code{\link[LCPA]{get.entropy}}, \code{\link[LCPA]{get.AvePP}}), and a bootstrapped
#' likelihood-ratio test (i.e., BLRT, \code{\link[LCPA]{LRT.test.Bootstrap}}, which is very time-consuming).
#' Above all and the most important criterion, is that the better model is the one that aligns with theoretical
#' expectations and offers clear interpretability.
#'
#' The \code{\link[LCPA]{LRT.test}} test statistic is defined as:
#' \itemize{
#'   \item The models must be \emph{nested} (i.e., the model with fewer parameters is a constrained version of the more one).
#'   \item Both models must be fit on the identical dataset with the same response variables.
#'   \item The test statistic asymptotically follows a chi-square distribution.
#' }
#'
#' \deqn{LRT = -2 \times (\text{LogLik}_{1} - \text{LogLik}_{2})}
#' where:
#' \itemize{
#'   \item \eqn{\text{LogLik}_{1}}: Log-likelihood of the smaller model (fewer parameters).
#'   \item \eqn{\text{LogLik}_{2}}: Log-likelihood of the larger model (more parameters).
#' }
#' Under the null hypothesis (\code{H_0}: small model is true), LRT asymptotically follows
#' a chi-square distribution with \eqn{df} degrees of freedom.
#'
#' @importFrom stats pchisq
#' @export
#'
LRT.test <- function(object1, object2) {
  if (!identical(class(object1), class(object2))) {
    stop("Model classes must be identical. Both objects must be either 'LCA' or 'LPA' type.")
  }

  valid_classes <- c("LCA", "LPA")
  if (!(class(object1) %in% valid_classes)) {
    stop("Invalid model class. Only supported for 'LCA' or 'LPA' class.")
  }

  N <- nrow(object1$arguments$response)

  npar1 <- object1$npar
  npar2 <- object2$npar

  if (npar2 < npar1) {
    model1 <- object2
    model2 <- object1
    L1 <- object2$arguments$L
    L2 <- object1$arguments$L
    npar1 <- object2$npar
    npar2 <- object1$npar
  } else {
    model1 <- object1
    model2 <- object2
    L1 <- object1$arguments$L
    L2 <- object2$arguments$L
  }

  fit.index1 <- get.fit.index(model1)
  fit.index2 <- get.fit.index(model2)

  LRT.statistic <- -2 * (fit.index1$Log.Lik - fit.index2$Log.Lik)
  df <- abs(npar2 - npar1)
  p.value <- pchisq(q = LRT.statistic, df = df, lower.tail = FALSE)

  res <- list(
    statistic = c(`LRT` = LRT.statistic),
    parameter = c(df = df),
    p.value = c(p.value = p.value),
    method = "Likelihood Ratio Test",
    data.name = paste("\nModel with", L1, "classes ( npar=", npar1, ") \n              vs\n Model with", L2, "classes ( npar=", npar2, ")")
  )
  class(res) <- "htest"

  return(res)
}
