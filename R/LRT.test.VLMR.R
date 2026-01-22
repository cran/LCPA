#' Lo-Mendell-Rubin likelihood ratio test
#'
#' Implements the ad-hoc adjusted likelihood ratio test (LRT) described in Formula 15 of Lo, Mendell,
#' & Rubin (2001), also known as VLMR LRT (Vuong-Lo-Mendell-Rubin Adjusted LRT). This method is typically
#' used to compare models with \code{L-1} and \code{L} classes. If the difference in the number of
#' classes exceeds 1, conclusions should be interpreted with extreme caution.
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
#' of \code{\link[LCPA]{LRT.test.VLMR}} may not be accurate and is provided for reference only.
#' More reliable conclusions should be based on a combination of fit indices (i.e., \code{\link[LCPA]{get.fit.index}}),
#' classification accuracy measures (i.e., \code{\link[LCPA]{get.entropy}}, \code{\link[LCPA]{get.AvePP}}), and a bootstrapped
#' likelihood-ratio test (i.e., BLRT, \code{\link[LCPA]{LRT.test.Bootstrap}}, which is very time-consuming).
#' Above all and the most important criterion, is that the better model is the one that aligns with theoretical
#' expectations and offers clear interpretability.
#'
#' The \code{\link[LCPA]{LRT.test.VLMR}} statistic is defined as:
#' \deqn{VLMR = \frac{LRT}{c} \quad \text{where} \quad c = 1 + \frac{1}{df \cdot \log(N)}}
#' where:
#' \itemize{
#'   \item \eqn{LRT} is the standard likelihood ratio statistic. see \code{\link[LCPA]{LRT.test}}
#'   \item \eqn{df = npar_2 - npar_1} is the difference in number of free parameters between models.
#'   \item \eqn{N} is the sample size.
#' }
#' Under the null hypothesis (\code{H_0}: small model is true), VLMR asymptotically follows
#' a chi-square distribution with \eqn{df} degrees of freedom.
#'
#' @references
#' Lo, Y., Mendell, N. R., & Rubin, D. B. (2001). Testing the number of components in a normal mixture. Biometrika, 88(3), 767-778. https://doi.org/10.1093/biomet/88.3.767
#'
#' Nylund-Gibson, K., & Choi, A. Y. (2018). Ten frequently asked questions about latent class analysis. Translational Issues in Psychological Science, 4(4), 440-461. https://doi.org/10.1037/tps0000176
#'
#' @importFrom stats pchisq
#' @export
LRT.test.VLMR <- function(object1, object2) {
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

  c.factor <- 1 + 1 / (df * log(N))

  modlr.statistic <- LRT.statistic / c.factor

  p.value <- pchisq(q = modlr.statistic, df = df, lower.tail = FALSE)

  res <- list(
    statistic = c(`Adjusted LRT` = modlr.statistic),
    parameter = c(df = df),
    p.value = c(p.value = p.value),
    method = "Lo-Mendell-Rubin Likelihood Ratio Test",
    data.name = paste("\nModel with", L1, "classes ( npar=", npar1, ") \n              vs\n Model with", L2, "classes ( npar=", npar2, ")")
  )
  class(res) <- "htest"

  return(res)
}
