#' Likelihood Ratio Test
#'
#' Conducts a likelihood ratio test to compare the fit of two LCA or LPA
#' models with any numbers of latent classes/profiles, including equal class
#' counts. The test evaluates whether a model with more parameters provides a
#' significantly better fit than a model with fewer parameters.
#'
#' @param object1 Fitted LCA or LPA model. When both models have the same
#'   number of free parameters, this is treated as the null model.
#' @param object2 Fitted LCA or LPA model of the same type. When both models
#'   have the same number of free parameters, this is treated as the
#'   alternative model.
#'
#' @return An object of class \code{"htest"} containing:
#' \itemize{
#'   \item \code{statistic}: Standard likelihood ratio test statistic
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
#' Models may have any class counts; they do not need to differ by exactly one
#' class. If both models have the same number of free parameters, the
#' likelihood-ratio statistic is returned but the chi-square p-value is
#' \code{NA} because its reference distribution has zero degrees of freedom.
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

  models <- .order.LRT.models(object1, object2)
  model1 <- models$null
  model2 <- models$alternative
  L1 <- model1$arguments$L
  L2 <- model2$arguments$L
  npar1 <- model1$npar
  npar2 <- model2$npar

  fit.index1 <- get.fit.index(model1)
  fit.index2 <- get.fit.index(model2)

  LRT.statistic <- -2 * (fit.index1$Log.Lik - fit.index2$Log.Lik)
  df <- npar2 - npar1
  p.value <- if(df > 0){
    pchisq(q = LRT.statistic, df = df, lower.tail = FALSE)
  }else{
    NA_real_
  }

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

.order.LRT.models <- function(object1, object2) {
  npar <- c(object1$npar, object2$npar)
  if(length(npar) != 2L || any(!is.finite(npar))){
    stop("Both models must have a finite number of free parameters")
  }
  if(npar[1L] <= npar[2L]){
    list(null = object1, alternative = object2)
  }else{
    list(null = object2, alternative = object1)
  }
}
