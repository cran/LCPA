#' @title Model Comparison Tool
#'
#' @description
#' Compares two nested latent class/profile models using multiple fit indices, likelihood ratio tests, and classification metrics.
#'
#' @param object1 An object of class \code{\link[LCPA]{LCA}} or \code{\link[LCPA]{LPA}}, representing
#'                the first latent class/profile model.
#' @param object2 An object of class \code{\link[LCPA]{LCA}} or \code{\link[LCPA]{LPA}}, representing
#'                the second latent class/profile model. Must be of the same type as \code{object1}.
#' @param n.Bootstrap Integer specifying the number of bootstrap replications for the parametric
#'                    bootstrap likelihood ratio test (BLRT). Default is \code{0} (no bootstrap test performed).
#'
#' @return An object of class \code{compare.model} containing:
#' \describe{
#'   \item{\code{npar}}{Named vector with number of free parameters for each model}
#'   \item{\code{entropy}}{Named vector with entropy values (classification accuracy measure) for each model}
#'   \item{\code{AvePP}}{List containing average posterior probabilities per latent class/profile}
#'   \item{\code{fit.index}}{List of \code{\link[LCPA]{get.fit.index}} objects for both models}
#'   \item{\code{BF}}{Bayes Factor for model comparison (based on SIC)}
#'   \item{\code{LRT.obj}}{Likelihood ratio test (LRT) results}
#'   \item{\code{LRT.VLMR.obj}}{Vuong-Lo-Mendell-Rubin (VLMR) adjusted LRT results}
#'   \item{\code{LRT.Bootstrap.obj}}{Bootstrap LRT results (if \code{n.Bootstrap > 0})}
#'   \item{\code{call}}{The matched function call}
#'   \item{\code{arguments}}{List containing the original arguments passed to the function}
#' }
#'
#' @details
#' This function performs comprehensive model comparison between two nested LCA/LPA models. Key features include:
#' \itemize{
#'   \item Automatically orders models by parameter count (smaller model first)
#'   \item Computes multiple fit indices via \code{\link[LCPA]{get.fit.index}}
#'   \item Calculates classification quality metrics (entropy, average posterior probabilities)
#'   \item Performs three types of likelihood ratio tests:
#'     \itemize{
#'       \item Standard LRT, see \code{\link[LCPA]{LRT.test}}
#'       \item VLMR adjusted LRT, see \code{\link[LCPA]{LRT.test.VLMR}}
#'       \item Parametric bootstrap LRT (computationally intensive but robust), see \code{\link[LCPA]{LRT.test.Bootstrap}}
#'     }
#'   \item Computes Bayes Factor using Sample-Size Adjusted BIC (SIC)
#' }
#'
#' \strong{Important Requirements}:
#' \itemize{
#'   \item Both models must be of the same type (\code{LCA} or \code{LPA})
#'   \item Models must be nested (one model is a constrained version of the other)
#'   \item \code{n.Bootstrap > 0} requires significant computational resources
#' }
#'
#' @examples
#' library(LCPA)
#' set.seed(123)
#'
#' data.obj <- sim.LPA(N = 500, I = 5, L = 3, constraint = "V0")
#' response <- data.obj$response
#'
#' # need Mplus
#' \dontrun{
#' # Compare 3-class vs 4-class LPA models
#' object1 <- LPA(response, L = 3, method = "Mplus", constraint = "V0")
#' object2 <- LPA(response, L = 4, method = "Mplus", constraint = "V0")
#'
#' compare.model.obj <- compare.model(object1, object2)
#'
#' print(compare.model.obj)
#' }
#'
#' @seealso
#' \code{\link[LCPA]{LCA}}, \code{\link[LCPA]{LPA}}, \code{\link[LCPA]{get.fit.index}},
#' \code{\link[LCPA]{extract}}, \code{\link[LCPA]{LRT.test}}, \code{\link[LCPA]{LRT.test.VLMR}}
#'
#' @export
compare.model <- function(object1, object2, n.Bootstrap=0){

  if (!identical(class(object1), class(object2))) {
    stop("Model classes must be identical. Both objects must be either 'LCA' or 'LPA' type.")
  }

  valid_classes <- c("LCA", "LPA")
  if (!(class(object1) %in% valid_classes)) {
    stop("Invalid model class. Only supported for 'LCA' or 'LPA' class.")
  }

  call <- match.call()

  model.type <- class(object1)
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

  entropy1 <- get.entropy(model1)
  entropy2 <- get.entropy(model2)
  AvePP1 <- get.AvePP(model1)
  AvePP2 <- get.AvePP(model2)
  fit.index1 <- get.fit.index(model1)
  fit.index2 <- get.fit.index(model2)

  BF <- exp(fit.index2$SIC - fit.index1$SIC)

  LRT.obj <- LRT.test(object1, object2)
  LRT.VLMR.obj <- LRT.test.VLMR(object1, object2)
  if(n.Bootstrap){
    LRT.Bootstrap.obj <- LRT.test.Bootstrap(object1, object2, n.Bootstrap)
  }else{
    LRT.Bootstrap.obj <- NULL
  }

  N <- c(model1=nrow(object1$arguments$response), model2=nrow(object2$arguments$response))
  I <- c(model1=ncol(object1$arguments$response), model2=ncol(object2$arguments$response))
  L <- c(model1=object1$arguments$L, model2=object2$arguments$L)
  res <- list(N=N, I=I, L=L, npar=c(model1=npar1, model2=npar2),
              entropy=c(model1=entropy1, model2=entropy2),
              AvePP=list(model1=AvePP1, model2=AvePP2),
              fit.index=list(model1=fit.index1, model2=fit.index2),
              BF=BF,
              LRT.obj=LRT.obj, LRT.VLMR.obj=LRT.VLMR.obj,
              LRT.Bootstrap.obj=LRT.Bootstrap.obj)

  res$call <- call
  res$arguments = list( object1=object1, object2=object2, n.Bootstrap=n.Bootstrap)

  class(res) <- "compare.model"

  return(res)
}
