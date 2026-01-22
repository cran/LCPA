#' Bootstrap Likelihood Ratio Test
#'
#' Conducts a bootstrap likelihood ratio test to compare the fit of two nested models.
#' This test evaluates whether a model with more parameters provides a significantly
#' better fit than a model with fewer parameters by approximating the null distribution
#' through parametric bootstrapping.
#'
#' @param object1 Fitted model object with fewer parameters (i.e., fewer \code{npar}, small model).
#' @param object2 Fitted model object with more parameters (i.e., more \code{npar}, large model).
#' @param n.Bootstrap Number of bootstrap replicates (default = 100). Higher values increase
#'                    accuracy but computation time (we suggest that \code{n.Bootstrap = 1000}).
#' @param vis Logical. If \code{TRUE}, displays progress information during bootstrap (default: \code{TRUE}).
#'
#' @return An object of class \code{"htest"} containing:
#' \itemize{
#'   \item \code{statistic}: Observed likelihood ratio test statistic
#'   \item \code{parameter}: Degrees of freedom (reported as \code{NA} since p-value is bootstrap-derived)
#'   \item \code{p.value}: Bootstrap p-value
#'   \item \code{method}: Name of the test ("Bootstrap Likelihood Ratio Test")
#'   \item \code{data.name}: Model comparison description
#' }
#'
#' @details
#' Note that even the result of \code{\link[LCPA]{LRT.test.Bootstrap}} should not be taken as the sole
#' criterion; fit indices (e.g., \code{\link[LCPA]{get.fit.index}}) and classification accuracy
#' measures (e.g., \code{\link[LCPA]{get.entropy}}, \code{\link[LCPA]{get.AvePP}}) must be considered together.
#' Above all and the most important criterion, is that the better model is the one that aligns with theoretical
#' expectations and offers clear interpretability.
#'
#' The \code{\link[LCPA]{LRT.test.Bootstrap}} statistic is calculated as:
#' \deqn{LRT = -2 \times (\text{LogLik}_{1} - \text{LogLik}_{2})}
#' where:
#' \itemize{
#'   \item \eqn{\text{LogLik}_{1}}: Log-likelihood of the smaller model (fewer parameters).
#'   \item \eqn{\text{LogLik}_{2}}: Log-likelihood of the larger model (more parameters).
#' }
#'
#' The \code{\link[LCPA]{LRT.test.Bootstrap}} function employs a parametric bootstrap procedure to
#' empirically estimate the distribution of the LRT statistic under the null hypothesis (that the
#' smaller model is sufficient). The specific steps are as follows:
#' \enumerate{
#'   \item Parameter Extraction: The estimated parameters (\code{params}) from the smaller model
#'   (\code{object1}) are treated as the true population values (ground truth).
#'   \item Data Simulation: The function invokes \code{\link[LCPA]{sim.LCA}} or \code{\link[LCPA]{sim.LPA}}
#'   to generate \code{n.Bootstrap} independent datasets. Each dataset maintains the same sample size (\code{N})
#'   and number of indicators (\code{I}) as the original empirical data.
#'   \item Model Re-fitting: For each simulated dataset, both the small model and the large model are
#'   re-fitted. To ensure consistency, the estimation settings (e.g., convergence criteria, iterations)
#'   are identical to those used for the original models, with the exception that
#'   \code{\link[LCPA]{LRT.test.Bootstrap}} forces \code{par.ini = "random"} to avoid local maxima.
#'   \item Distribution Generation: This process generates \code{n.Bootstrap} pairs of
#'   \eqn{\text{LogLik}_{1, boot}} and \eqn{\text{LogLik}_{2, boot}}, which are used to compute
#'   a collection of bootstrap LRT statistics:
#'   \eqn{LRT_{boot} = -2 \times (\text{LogLik}_{1, boot} - \text{LogLik}_{2, boot})}.
#'   \item P-value Calculation: The bootstrap \eqn{p}-value is calculated as the proportion of
#'   simulated \eqn{LRT_{boot}} values that are greater than or equal to the observed \eqn{LRT}
#'   statistic from the original data.
#' }
#'
#' This method is particularly recommended for Latent Class and Latent Profile Analysis because
#' the traditional Chi-square distribution for LRT often fails to hold due to parameters
#' being on the boundary of the parameter space (e.g., probabilities near 0 or 1).
#'
#' @importFrom stats pchisq
#' @export
#'
LRT.test.Bootstrap <- function(object1, object2, n.Bootstrap=100, vis=TRUE) {
  if (!identical(class(object1), class(object2))) {
    stop("Model classes must be identical. Both objects must be either 'LCA' or 'LPA' type.")
  }

  valid_classes <- c("LCA", "LPA")
  if (!(class(object1) %in% valid_classes)) {
    stop("Invalid model class. Only supported for 'LCA' or 'LPA' class.")
  }

  N <- nrow(object1$arguments$response)

  if (object1$npar > object2$npar) {
    model1 <- object2
    model2 <- object1
  } else {
    model1 <- object1
    model2 <- object2
  }
  L1 <- model1$arguments$L
  L2 <- model2$arguments$L
  npar1 <- model1$npar
  npar2 <- model2$npar

  LRT.statistic <- -2 * (model1$Log.Lik - model2$Log.Lik)

  params1 <- model1$params
  params2 <- model2$params

  response <- model1$arguments$response
  N <- nrow(response)
  I <- ncol(response)
  nrep <- 1
  starts <- 1
  maxiter.wa <- 1

  if(inherits(model1, "LCA")){
    poly.value <- rep(0, I)
    poly.temp <- vector("list", I)
    for (i in 1:I) {
      valid_responses <- response[, i][!is.na(response[, i])]
      poly.temp[[i]] <- sort(unique(valid_responses))
      poly.value[i] <- length(poly.temp[[i]])
    }
  }

  LogLik <- matrix(0, n.Bootstrap, 2)
  for(r in 1:n.Bootstrap){
    if(vis){
      cat('\r', paste0(L1, " classes (npar=", npar1, ") vs ", L2, " classes (npar=", npar2, "):"), ' bootstrapping r =', sprintf("%4d", r), '/', n.Bootstrap)
    }
    if(inherits(model1, "LCA")){
      data.obj.cur <- sim.LCA(N=N, I=I, L=L1, poly.value=poly.value, params=params1)
    }else if(inherits(model1, "LPA")){
      data.obj.cur <- sim.LPA(N=N, I=I, L=L1, constraint=model1$arguments$constraint, params=params1)
    }
    response.cur <- data.obj.cur$response

    model1.cur <- update(model1, response=response.cur, par.ini="random", vis=FALSE)
    model2.cur <- update(model2, response=response.cur, par.ini="random", vis=FALSE)

    LogLik[r, 1] <- model1.cur$Log.Lik
    LogLik[r, 2] <- model2.cur$Log.Lik
  }
  LRT.Bootstrap <- -2 * (LogLik[, 1] - LogLik[, 2])

  p.value <- (sum(LRT.Bootstrap >= LRT.statistic) + 1e-8) / (n.Bootstrap +  1e-8)

  res <- list(
    statistic = c(`LRT` = LRT.statistic),
    parameter = c(df = NA),
    p.value = c(p.value = p.value),
    method = "Bootstrap Likelihood Ratio Test",
    data.name = paste0("\nModel with ", L1, " classes (npar=", npar1, ") \n              vs\n Model with ", L2, " classes (npar=", npar2, ")")
  )
  class(res) <- "htest"

  if(vis){
    cat("\n")
  }
  return(res)
}
