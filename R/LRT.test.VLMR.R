#' Vuong-Lo-Mendell-Rubin likelihood ratio test
#'
#' Computes the Mplus TECH11 implementation of the Vuong-Lo-Mendell-Rubin
#' likelihood ratio test (VLMR; Vermunt, 2024) and the Lo-Mendell-Rubin adjusted test (aLMR; Lo et al., 2001)
#' for LCA or LPA models with any numbers of latent classes/profiles, including
#' equal class counts.
#'
#' @param object1 Fitted \code{LCA} or \code{LPA} model. When both models have
#'   the same number of free parameters, this is treated as the null model.
#' @param object2 Fitted model of the same type and fitted to the same data.
#'   When both models have the same number of free parameters, this is treated
#'   as the alternative model.
#'
#' @return An object of classes \code{"VLMR"} and \code{"htest"} containing:
#' \itemize{
#'   \item \code{statistic}: Unadjusted VLMR likelihood ratio statistic.
#'   \item \code{p.value}: VLMR p-value from the Mplus weighted chi-square distribution.
#'   \item \code{adjusted.statistic}: LMR statistic adjusted using Equation
#'     (15) of Lo et al. (2001).
#'   \item \code{adjusted.p.value}: Adjusted LMR p-value from the same weighted chi-square distribution.
#'   \item \code{correction.factor}: Correction factor from Equation (15) of
#'     Lo et al. (2001).
#'   \item \code{distribution}: Mean, standard deviation, and eigenvalue weights
#'     of the estimated VLMR reference distribution.
#'   \item \code{diagnostics}: Score, information-matrix, eigenvalue, and Imhof
#'     integration diagnostics.
#' }
#'
#' @details
#' This function reports two related tests. The unadjusted VLMR is based on the
#' likelihood-ratio test for nested or overlapping models developed by Vuong
#' (1989), which Lo, Mendell, and Rubin (2001) applied to the comparison of
#' \eqn{K}- and \eqn{K+1}-component mixture models. The adjusted LMR test
#' (aLMR) is the modification proposed by Lo et al. (2001). The robust
#' reference distribution follows the Mplus TECH11 reconstruction described by
#' Vermunt (2024) and is calculated from observation-level scores and
#' observed-information matrices as follows.
#'
#' Let \eqn{H_0} denote the null model and \eqn{H_1} the alternative model,
#' ordered by parameter count as described above. For model
#' \eqn{r \in \{0,1\}}, let
#' \eqn{\ell_{rn}} be participant \eqn{n}'s log-likelihood contribution,
#' \eqn{s_{rn}} its score vector, and
#' \eqn{I_r=-\sum_{n=1}^N\partial^2\ell_{rn}/
#' \partial\theta_r\partial\theta_r^{\mathsf T}} the observed information
#' matrix. Define
#' \deqn{
#'   B_r = \sum_{n=1}^{N} s_{rn}s_{rn}^{\mathsf T}, \qquad
#'   B_{10} = \sum_{n=1}^{N} s_{1n}s_{0n}^{\mathsf T},
#' }
#' and the robust sandwich covariance matrix
#' \deqn{
#'   R_r = I_r^{-1}B_rI_r^{-1}.
#' }
#' Following the Mplus implementation identified by Vermunt (2024), the
#' reference-distribution weights are the eigenvalues of
#' \deqn{
#'   W_{\mathrm{Mplus}} =
#'   \left[\begin{array}{cc}
#'     B_1R_1 & B_{10}R_0 \\
#'     -B_{10}^{\mathsf T}R_1 & -B_0R_0
#'   \end{array}\right].
#' }
#' If these eigenvalues are \eqn{\lambda_1,\ldots,\lambda_m}, the unadjusted
#' statistic and its estimated null distribution are
#' \deqn{
#'   \mathrm{LR} = 2\{\ell(\widehat\theta_1)-\ell(\widehat\theta_0)\},
#'   \qquad
#'   Q = \sum_{j=1}^{m}\lambda_j\chi^2_{1,j}.
#' }
#' Thus, the VLMR p-value is \eqn{\Pr(Q \geq \mathrm{LR})}. It is evaluated
#' using Imhof's (1961) method, retaining negative eigenvalue weights. The
#' reported mean and standard deviation of the reference distribution are
#' \eqn{\sum_j\lambda_j} and
#' \eqn{\{2\sum_j\lambda_j^2\}^{1/2}}, respectively (Vermunt, 2024).
#'
#' The aLMR adjustment is specifically Equation (15) of Lo et al. (2001):
#' \deqn{
#'   \mathrm{aLMR} = \frac{\mathrm{LR}}{c}, \qquad
#'   c = 1 + \frac{1}{(p-q)\log(N)},
#' }
#' where \eqn{p-q} is the difference in the numbers of free parameters. The
#' adjusted statistic uses the same weighted chi-square reference distribution
#' as the unadjusted VLMR statistic (Vermunt, 2024). Consequently,
#' \code{p.value} is the unadjusted VLMR result derived from Vuong (1989) and
#' Lo et al. (2001), whereas \code{adjusted.p.value} is the aLMR result based
#' on Equation (15) of Lo et al. (2001).
#' Models may have any class counts; they do not need to differ by exactly one
#' class. When their numbers of free parameters are equal, the unadjusted VLMR
#' is computed, but the aLMR correction factor, adjusted statistic, and adjusted
#' p-value are \code{NA} because Equation (15) contains \eqn{p-q} in the
#' denominator.
#'
#' The result matches the TECH11 calculation conditional on the two supplied
#' maximum-likelihood solutions. Local maxima, singular information matrices,
#' and boundary solutions can invalidate the comparison.
#'
#' @references
#' Imhof, J. P. (1961). Computing the distribution of quadratic forms in normal
#' variables. *Biometrika, 48*(3--4), 419--426.
#' \doi{10.1093/biomet/48.3-4.419}
#'
#' Lo, Y., Mendell, N. R., & Rubin, D. B. (2001). Testing the number of
#' components in a normal mixture. *Biometrika, 88*(3), 767--778.
#' \doi{10.1093/biomet/88.3.767}
#'
#' Vermunt, J. K. (2024). The Vuong-Lo-Mendell-Rubin test for latent class and
#' latent profile analysis: A note on the different implementations in Mplus
#' and LatentGOLD. *Methodology, 20*(1), 72--83.
#' \doi{10.5964/meth.12467}
#'
#' Vuong, Q. H. (1989). Likelihood ratio tests for model selection and
#' non-nested hypotheses. *Econometrica, 57*(2), 307--333.
#' \doi{10.2307/1912557}
#'
#' @importFrom CompQuadForm imhof
#' @export
LRT.test.VLMR <- function(object1, object2) {
  model.class <- class(object1)[1L]
  if(!identical(model.class, class(object2)[1L])){
    stop("Model classes must be identical")
  }
  if(!model.class %in% c("LCA", "LPA")){
    stop("Only LCA and LPA models are supported")
  }

  response1 <- as.matrix(object1$arguments$response)
  response2 <- as.matrix(object2$arguments$response)
  if(!isTRUE(all.equal(response1, response2, tolerance = 0,
                       check.attributes = FALSE))){
    stop("Models must be fitted to the same observations in the same order")
  }
  if(model.class == "LPA" &&
     !identical(object1$arguments$constraint, object2$arguments$constraint)){
    stop("LPA models must use the same covariance constraint")
  }

  models <- .order.LRT.models(object1, object2)
  model0 <- models$null
  model1 <- models$alternative

  N <- nrow(response1)
  df <- model1$npar - model0$npar

  LRT.statistic <- 2 * (model1$Log.Lik - model0$Log.Lik)
  if(LRT.statistic < 0){
    warning("The alternative model has a lower log-likelihood; inspect model ordering and local maxima")
  }

  alternative <- .VLMR.model.components(model1)
  null <- .VLMR.model.components(model0)
  distribution <- .VLMR.weights(
    alternative$score.observation, alternative$information,
    null$score.observation, null$information
  )
  VLMR.p <- .VLMR.p.value(LRT.statistic, distribution$eigenvalues)

  if(df > 0){
    correction.factor <- 1 + 1 / (df * log(N))
    adjusted.statistic <- LRT.statistic / correction.factor
    adjusted.p <- .VLMR.p.value(adjusted.statistic, distribution$eigenvalues)
  }else{
    correction.factor <- NA_real_
    adjusted.statistic <- NA_real_
    adjusted.p <- list(p.value = NA_real_, error = NA_real_)
  }

  res <- list(
    statistic = c(`VLMR LRT` = LRT.statistic),
    parameter = c(df = df),
    p.value = c(p.value = VLMR.p$p.value),
    adjusted.statistic = c(`Adjusted LMR LRT` = adjusted.statistic),
    adjusted.p.value = c(p.value = adjusted.p$p.value),
    correction.factor = correction.factor,
    distribution = list(
      mean = sum(distribution$eigenvalues),
      sd = sqrt(2 * sum(distribution$eigenvalues^2)),
      eigenvalues = distribution$eigenvalues
    ),
    diagnostics = list(
      score.max.abs = c(
        null = max(abs(null$score)),
        alternative = max(abs(alternative$score))
      ),
      information.condition = distribution$information.condition,
      information.reciprocal.condition = distribution$information.reciprocal.condition,
      information.generalized.inverse = distribution$information.generalized.inverse,
      covariance.repaired = c(
        null = isTRUE(model0$covariance.repaired),
        alternative = isTRUE(model1$covariance.repaired)
      ),
      eigenvalue.max.imaginary = distribution$eigenvalue.max.imaginary,
      imhof.error = c(VLMR = VLMR.p$error, adjusted.LMR = adjusted.p$error)
    ),
    method = "Vuong-Lo-Mendell-Rubin Likelihood Ratio Test (Mplus TECH11)",
    data.name = paste0(
      "Model with ", model0$arguments$L, " classes (npar=", model0$npar,
      ") vs model with ", model1$arguments$L, " classes (npar=", model1$npar, ")"
    )
  )
  class(res) <- c("VLMR", "htest")
  res
}

.VLMR.model.components <- function(object) {
  if(inherits(object, "LCA")){
    specification <- lca.parameterization(object)
    if(any(!is.finite(specification$theta))){
      stop("VLMR requires interior LCA probability estimates")
    }
    components <- lca.score.components(
      specification$theta, specification, compute.information = TRUE
    )
    information <- components$information
  }else{
    specification <- lpa.parameterization(object)
    if(any(!is.finite(specification$theta))){
      stop("VLMR requires finite LPA parameter estimates")
    }
    components <- lpa.score.components(
      specification$theta, specification
    )
    if(is.null(components)){
      stop("Fitted LPA covariance matrices must be positive definite")
    }
    information <- louis.SE.LPA(object)$hessian
  }

  if(length(specification$theta) != object$npar){
    stop("The fitted parameter count does not match the VLMR parameterization")
  }
  if(any(!is.finite(components$score.observation)) ||
     any(!is.finite(information))){
    stop("VLMR score and information matrices must be finite")
  }
  list(
    score.observation = components$score.observation,
    score = components$score,
    information = information
  )
}

.VLMR.robust.components <- function(score.observation, information) {
  information <- (information + t(information)) / 2
  reciprocal.condition <- rcond(information)
  use.ginv <- !is.finite(reciprocal.condition) ||
    reciprocal.condition < sqrt(.Machine$double.eps)
  information.inverse <- NULL
  if(!use.ginv){
    information.inverse <- tryCatch(
      solve(information),
      error = function(e) NULL
    )
  }
  if(is.null(information.inverse) || any(!is.finite(information.inverse))){
    information.inverse <- MASS::ginv(information)
    use.ginv <- TRUE
  }
  if(any(!is.finite(information.inverse))){
    stop("VLMR generalized inverse of the observed information matrix is not finite")
  }
  B <- crossprod(score.observation)
  list(
    B = B,
    covariance = information.inverse %*% B %*% information.inverse,
    condition = kappa(information, exact = TRUE),
    reciprocal.condition = reciprocal.condition,
    generalized.inverse = use.ginv
  )
}

.VLMR.weights <- function(score.alternative, information.alternative,
                          score.null, information.null) {
  if(nrow(score.alternative) != nrow(score.null)){
    stop("VLMR score matrices must contain the same observations")
  }
  alternative <- .VLMR.robust.components(
    score.alternative, information.alternative
  )
  null <- .VLMR.robust.components(score.null, information.null)
  B.cross <- crossprod(score.alternative, score.null)

  W <- rbind(
    cbind(
      alternative$B %*% alternative$covariance,
      B.cross %*% null$covariance
    ),
    cbind(
      -t(B.cross) %*% alternative$covariance,
      -null$B %*% null$covariance
    )
  )
  eigenvalues <- eigen(W, only.values = TRUE)$values
  eigenvalue.scale <- max(1, max(abs(Re(eigenvalues))))
  imaginary <- max(abs(Im(eigenvalues)))
  if(imaginary > 1e-6 * eigenvalue.scale){
    warning("VLMR matrix has non-negligible complex eigenvalues")
  }

  list(
    eigenvalues = Re(eigenvalues),
    information.condition = c(
      null = null$condition,
      alternative = alternative$condition
    ),
    information.reciprocal.condition = c(
      null = null$reciprocal.condition,
      alternative = alternative$reciprocal.condition
    ),
    information.generalized.inverse = c(
      null = null$generalized.inverse,
      alternative = alternative$generalized.inverse
    ),
    eigenvalue.max.imaginary = imaginary
  )
}

.VLMR.p.value <- function(statistic, eigenvalues) {
  threshold <- sqrt(.Machine$double.eps) * max(1, max(abs(eigenvalues)))
  weights <- eigenvalues[abs(eigenvalues) > threshold]
  if(length(weights) == 0L){
    return(list(p.value = NA_real_, error = NA_real_))
  }
  result <- imhof(statistic, weights)
  list(
    p.value = min(max(result$Qq, 0), 1),
    error = result$abserr
  )
}

#' @export
print.VLMR <- function(x, digits = getOption("digits"), ...) {
  statistic.digits <- max(1L, digits - 2L)
  p.digits <- max(1L, digits - 3L)
  cat("\n\t", x$method, "\n\n", sep = "")
  cat("data:  ", x$data.name, "\n", sep = "")
  cat(
    "VLMR LRT = ", format(x$statistic, digits = statistic.digits),
    ", p-value = ", format.pval(x$p.value, digits = p.digits), "\n",
    sep = ""
  )
  cat(
    "Adjusted LMR LRT = ",
    format(x$adjusted.statistic, digits = statistic.digits),
    ", p-value = ", format.pval(x$adjusted.p.value, digits = p.digits), "\n",
    sep = ""
  )
  cat(
    "Reference distribution: mean = ",
    format(x$distribution$mean, digits = statistic.digits),
    ", SD = ", format(x$distribution$sd, digits = statistic.digits),
    "\n\n", sep = ""
  )
  invisible(x)
}
