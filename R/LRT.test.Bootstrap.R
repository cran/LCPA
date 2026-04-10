#' Bootstrap Likelihood Ratio Test for Latent Class/Profile Models
#'
#' Conducts a bootstrap likelihood ratio test (BLRT) to compare two nested latent class analysis (LCA)
#' or latent profile analysis (LPA) models. The test evaluates whether a model with more classes
#' (the alternative model) provides significantly better fit than a model with fewer classes
#' (the null model). Implements both fixed-replicate and sequential stopping procedures for computational efficiency.
#'
#' @param object1 Fitted model object representing the null model (e.g., \code{L-1} classes).
#'                Must be of class \code{"LCA"} or \code{"LPA"}.
#' @param object2 Fitted model object representing the alternative model (e.g., \code{L} classes).
#'                Must be of the same class as \code{object1}.
#' @param n.Bootstrap Maximum number of bootstrap replicates (default = 100).
#'                    If \code{use.sequential = FALSE}, exactly this many replicates are performed.
#'                    If \code{use.sequential = TRUE} (default), this is an upper bound; the algorithm may stop early.
#'                    McLachlan & Peel (2000) suggest 100 replicates as typically sufficient for stable p-values,
#'                    especially near the conventional \eqn{\alpha = 0.05} threshold.
#'                    This default also aligns with the implementation in Nylund et al. (2007).
#' @param vis Logical. If \code{TRUE} (default), displays real-time progress during bootstrapping.
#' @param use.sequential Logical. If \code{TRUE} (default), applies the sequential stopping rule from
#'                       Nylund et al. (2007) to terminate early when the decision is statistically clear.
#'                       If \code{FALSE}, performs exactly \code{n.Bootstrap} replicates (traditional fixed bootstrap).
#'
#' @return An object of class \code{"htest"} containing:
#' \itemize{
#'   \item \code{statistic}: Observed likelihood ratio test statistic,
#'         \eqn{-2 \times (\text{logLik}_{\text{null}} - \text{logLik}_{\text{alt}})}.
#'   \item \code{parameter}: Degrees of freedom (set to \code{NA}; the bootstrap p-value does not rely on asymptotic chi-square distribution).
#'   \item \code{p.value}: Bootstrap p-value, computed as the proportion of bootstrap LRT statistics
#'         greater than or equal to the observed LRT statistic.
#'   \item \code{method}: Character string indicating the procedure used:
#'         either "Bootstrap LRT with Sequential Stopping" or "Bootstrap LRT (Fixed Replicates)".
#'   \item \code{data.name}: Descriptive string summarizing the model comparison,
#'         including the actual number of bootstrap replicates performed and, if applicable,
#'         the reason for early stopping under the sequential rule.
#'   \item \code{LRT.Bootstrap}: Numeric vector of length \eqn{R} (where \eqn{R} is the actual number of
#'         bootstrap replicates performed) containing the bootstrap LRT statistics:
#'         \eqn{LRT_{\text{boot}} = -2 \times (\text{logLik}_{\text{null,boot}} - \text{logLik}_{\text{alt,boot}})}.
#'         These are generated under the null hypothesis (i.e., data simulated from the smaller model).
#'         This vector can be used for diagnostic purposes, such as inspecting the empirical null distribution
#'         or computing alternative p-value estimates.
#' }
#'
#' @details
#' Core Workflow (Parametric Bootstrap):
#' \itemize{
#'   \item Parameter Extraction: Parameters from the null model (\code{object1}) are treated as the population truth.
#'   \item Data Simulation: Generates datasets using \code{\link[LCPA]{sim.LCA}} or \code{\link[LCPA]{sim.LPA}},
#'         preserving the original sample size and indicator structure.
#'   \item Model Re-fitting: For each simulated dataset, both the null and alternative models are re-estimated.
#'         Initial parameters are set to random values (\code{par.ini = "random"}) to avoid convergence to local maxima;
#'         the EM algorithm includes fallback strategies upon non-convergence.
#'   \item LRT Distribution: Computes bootstrap LRT statistics:
#'         \eqn{LRT_{\text{boot}} = -2 \times (\text{logLik}_{\text{null,boot}} - \text{logLik}_{\text{alt,boot}})}.
#'   \item P-value Calculation: The p-value is the proportion of bootstrap LRTs that are
#'         greater than or equal to the observed LRT statistic.
#' }
#'
#' Sequential Stopping Rule (Nylund et al., 2007):
#' When \code{use.sequential = TRUE}, the algorithm checks three stopping criteria after each replication:
#' \itemize{
#'   \item Upper Stopping Points: Stop if the current estimated p-value (\eqn{\hat{p}}) is
#'         greater than or equal to a pre-specified threshold at replication \eqn{n_i},
#'         indicating insufficient evidence against the null model.
#'   \item Lower Stopping Points: Stop if \eqn{\hat{p} \leq p_i} at replication \eqn{n_i},
#'         indicating strong evidence in favor of the alternative model.
#'   \item Conditional Lower Stopping Points: Stop if \eqn{\hat{p} = 0} \emph{and}
#'         the observed LRT exceeds the mean of the current bootstrap LRTs by more than \eqn{s} standard deviations.
#'         This leverages the approximate normality of the LRT distribution under the null (Nylund et al., 2007).
#' }
#' The exact stopping thresholds follow Appendix A of Nylund et al. (2007):
#' \itemize{
#'   \item Upper: (n, 2/n) for n = 2–3; (n, 3/n) for n = 4–9; (n, 4/n) for n = 10–17;
#'         (n, 5/n) for n = 18–26; (n, 6/n) for n = 27–99; (100, 0).
#'   \item Lower: (49, 0), (78, 1/78).
#'   \item Conditional: (5, 0, 20), (10, 0, 10), (20, 0, 5),
#'         where the third value denotes the required z-score threshold.
#' }
#' This rule ensures >95% agreement with the decision that would be made using an infinite number of replicates
#' when the true p-value is below 0.10. Near \eqn{p = 0.05}, the algorithm continues until the maximum
#' number of replicates to avoid premature conclusions.
#'
#' Critical Interpretation Notes:
#' \itemize{
#'   \item The BLRT should \strong{not} be used in isolation. Always integrate results with:
#'         \itemize{
#'           \item Information criteria: e.g., BIC, aBIC, CAIC via \code{\link[LCPA]{get.fit.index}}.
#'           \item Classification quality metrics: e.g., entropy (\code{\link[LCPA]{get.entropy}}),
#'                 average posterior probability (\code{\link[LCPA]{get.AvePP}}).
#'           \item Substantive interpretability and theoretical coherence of the extracted classes.
#'         }
#'   \item The BLRT is especially valuable in LCA/LPA contexts where the standard asymptotic LRT is invalid
#'         due to boundary parameter issues (e.g., class probabilities approaching 0 or 1; McLachlan & Peel, 2000).
#'   \item Early termination under sequential stopping (e.g., at \eqn{n \ll 100}) yields reliable p-values
#'         when the evidence is strong (\eqn{p \ll 0.05} or \eqn{p \gg 0.05}).
#'         However, if stopping occurs near \code{n.Bootstrap} with \eqn{p \approx 0.05},
#'         the result should be considered inconclusive.
#'   \item In small samples (\eqn{N < 300}) or with poorly separated classes,
#'         the BLRT may exhibit inflated Type I error rates; corroborate findings with BIC and other indices.
#' }
#'
#' @references
#' McLachlan, G. J., & Peel, D. (2000). Finite Mixture Models. John Wiley & Sons.
#'
#' Nylund, K. L., Asparouhov, T., & Muthén, B. O. (2007). Deciding on the Number of Classes in Latent Class Analysis and Growth Mixture Modeling: A Monte Carlo Simulation Study. Structural Equation Modeling: A Multidisciplinary Journal, 14(4), 535–569. https://doi.org/10.1080/10705510701575396
#'
#' @importFrom stats pchisq
#' @export
#'
LRT.test.Bootstrap <- function(object1, object2, n.Bootstrap = 100, vis = TRUE, use.sequential = TRUE) {
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

  if (inherits(model1, "LCA")) {
    poly.value <- rep(0, I)
    poly.temp <- vector("list", I)
    for (i in 1:I) {
      valid_responses <- response[, i][!is.na(response[, i])]
      poly.temp[[i]] <- sort(unique(valid_responses))
      poly.value[i] <- length(poly.temp[[i]])
    }
  }

  if (use.sequential) {
    upper.n <- c(2:3, 4:9, 10:17, 18:26, 27:99, 100)
    upper.p <- c(
      2 / (2:3),
      3 / (4:9),
      4 / (10:17),
      5 / (18:26),
      6 / (27:99),
      0
    )
    lower.n <- c(49, 78)
    lower.p <- c(0, 1/78)

    cond.n <- c(5, 10, 20)
    cond.p <- c(0, 0, 0)
    cond.s <- c(20, 10, 5)

    stop.reason <- ""

  } else {
    stop.reason <- "Fixed bootstrap (no sequential stopping)"
  }

  LRT.vec <- numeric(n.Bootstrap)
  actual.n <- n.Bootstrap

  current.p.hat <- 0
  for (r in 1:n.Bootstrap) {
    if (vis) {
      cat('\r', paste0(L1, " classes (npar=", npar1, ") vs ",
                       L2, " classes (npar=", npar2, "): ",
                       "bootstrapping r=",
                       sprintf(sprintf("%%%dd", ceiling(log10(r+1))), r), "/",
                       sprintf(sprintf("%%%dd", ceiling(log10(n.Bootstrap+1))), n.Bootstrap),
                       sprintf(", current p.hat=%.3f", current.p.hat)))
    }

    if (inherits(model1, "LCA")) {
      data.obj.cur <- sim.LCA(N = N, I = I, L = L1, poly.value = poly.value, params = params1)
    } else if (inherits(model1, "LPA")) {
      data.obj.cur <- sim.LPA(N = N, I = I, L = L1, constraint = model1$arguments$constraint, params = params1)
    }
    response.cur <- data.obj.cur$response

    model1.cur <- tryCatch({
      update(model1, response = response.cur, par.ini = "random", vis = FALSE)
    }, error = function(e) {
      update(model1, response = response.cur, par.ini = "random", method = "EM", vis = FALSE)
    })

    model2.cur <- tryCatch({
      update(model2, response = response.cur, par.ini = "random", vis = FALSE)
    }, error = function(e) {
      update(model2, response = response.cur, par.ini = "random", method = "EM", vis = FALSE)
    })

    LRT.vec[r] <- -2 * (model1.cur$Log.Lik - model2.cur$Log.Lik)

    current.p.hat <- mean(LRT.vec[1:r] >= LRT.statistic)
    if (use.sequential) {
      triggered <- FALSE

      if (r %in% upper.n) {
        p.thresh <- upper.p[which(upper.n == r)]
        if (current.p.hat >= p.thresh) {
          stop.reason <- paste0("Upper stop (n=", r, ", p.hat=", round(current.p.hat, 4), " >= ", round(p.thresh, 4), ")")
          actual.n <- r
          triggered <- TRUE
        }
      }

      if (!triggered && r %in% lower.n) {
        p.thresh <- lower.p[which(lower.n == r)]
        if (current.p.hat <= p.thresh) {
          stop.reason <- paste0("Lower stop (n=", r, ", p.hat=", round(current.p.hat, 4), " <= ", round(p.thresh, 4), ")")
          actual.n <- r
          triggered <- TRUE
        }
      }

      if (!triggered && r %in% cond.n) {
        idx <- which(cond.n == r)
        if (current.p.hat <= cond.p[idx]) {
          boot.sd <- sd(LRT.vec[1:r])
          if (boot.sd > 0) {
            z.score <- (LRT.statistic - mean(LRT.vec[1:r])) / boot.sd
            if (z.score > cond.s[idx]) {
              stop.reason <- paste0("Conditional lower stop (n=", r, ", z=", round(z.score, 2), " > ", cond.s[idx], ")")
              actual.n <- r
              triggered <- TRUE
            }
          }
        }
      }

      if (triggered) break
      if (r == n.Bootstrap) stop.reason <- "Max bootstrap reached"
    }
  }

  if (use.sequential && actual.n < n.Bootstrap) {
    LRT.Bootstrap <- LRT.vec[1:actual.n]
  } else {
    LRT.Bootstrap <- LRT.vec
    actual.n <- n.Bootstrap
  }
  p.value <- mean(LRT.Bootstrap >= LRT.statistic)

  base.name <- paste0(
    "\nModel with ", L1, " classes (npar=", npar1, ") \n              vs\n Model with ",
    L2, " classes (npar=", npar2, ")"
  )

  if (use.sequential) {
    data.name.str <- paste0(base.name, "\n[Sequential stopping: n=", actual.n, ", ", stop.reason, "]")
    method.str <- "Bootstrap Likelihood Ratio Test (Sequential Stopping)"
  } else {
    data.name.str <- paste0(base.name, "\n[Fixed bootstrap: n=", n.Bootstrap, "]")
    method.str <- "Bootstrap Likelihood Ratio Test (Fixed Bootstrap)"
  }

  res <- list(
    statistic = c(`LRT` = LRT.statistic),
    parameter = c(df = NA),
    p.value = c(p.value = p.value),
    method = method.str,
    data.name = data.name.str,
    LRT.Bootstrap = LRT.Bootstrap
  )
  class(res) <- "htest"

  if (vis) {
    cat("\n")
    if (use.sequential && stop.reason != "Max bootstrap reached") {
      cat("Stopping reason:", stop.reason, "\n")
    }
  }

  return(res)
}
