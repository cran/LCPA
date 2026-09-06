#' @title S3 Methods: summary
#'
#' @description
#' Generates structured, comprehensive summaries of objects produced by the \code{LCPA} package.
#' This generic function dispatches to class-specific methods that extract and organize key information
#' including model configurations, fit statistics, parameter estimates, simulation truths, and diagnostics.
#' Designed for programmatic access and downstream reporting.
#'
#' @param object An object of one of the following classes:
#'   \itemize{
#'     \item Model objects: \code{\link[LCPA]{LCA}}, \code{\link[LCPA]{LPA}},
#'       \code{\link[LCPA]{LCPA}}, \code{\link[LCPA]{LTA}}
#'     \item Simulation objects: \code{\link[LCPA]{sim.LCA}}, \code{\link[LCPA]{sim.LPA}},
#'       \code{\link[LCPA]{sim.LTA}}
#'     \item Fit/comparison objects: \code{\link[LCPA]{get.fit.index}}, \code{\link[LCPA]{compare.model}}
#'     \item Standard error objects: \code{\link[LCPA]{get.SE}}
#'   }
#' @param digits Number of decimal places for numeric output (default: 4). Applied universally across all methods.
#' @param I.max Maximum number of variables/items to display for \code{LCA},
#'   \code{LPA}, \code{sim.LCA}, \code{sim.LPA}, and \code{sim.LTA} summaries
#'   (default: 5).
#' @param L.max Maximum number of latent classes/profiles to display before truncation (\code{sim.LTA} only; default: 5).
#'   Useful when models have many latent groups. Ignored for other classes.
#' @param ... Reserved for S3 method compatibility; no additional arguments are used.
#'
#' @return A structured list whose S3 class identifies the corresponding
#'   summary method, such as \code{\link[LCPA]{summary.LCA}} or
#'   \code{\link[LCPA]{summary.LPA}}. Every summary method returns its object
#'   visibly, so an interactive call to \code{\link[base]{summary}()} dispatches
#'   automatically to the corresponding \code{\link[LCPA]{print}} method.
#'
#' @details
#' Each method returns a class-specific list designed both for its corresponding
#' print method and for programmatic access:
#' \describe{
#'   \item{\code{LCA}}{A \code{summary.LCA} object with:
#'     \describe{
#'       \item{\code{call}}{Original fitting call.}
#'       \item{\code{model.config}}{Number of latent classes and estimation method.}
#'       \item{\code{data.info}}{Sample size, item count, number of categories per item, and whether category counts are uniform.}
#'       \item{\code{fit.stats}}{Log-likelihood, AIC, BIC, entropy, and number of free parameters.}
#'       \item{\code{class.probs}}{Data frame containing class labels, modal-assignment counts, and estimated class proportions.}
#'       \item{\code{item.probs}}{Conditional response-probability matrices for the first \code{I.max} items.}
#'       \item{\code{convergence}}{Backend-specific algorithm, iteration, tolerance, initialization, and diagnostic information.}
#'       \item{\code{replication}}{Replication count and best BIC when replication selection applies; otherwise \code{NULL}.}
#'       \item{\code{digits}, \code{I.max.shown}, \code{total.items}}{Formatting and truncation metadata.}
#'     }}
#'
#'   \item{\code{LPA}}{A \code{summary.LPA} object with:
#'     \describe{
#'       \item{\code{call}}{Original fitting call.}
#'       \item{\code{model.config}}{Number of profiles, requested covariance constraint, its expanded description, and estimation method.}
#'       \item{\code{data.info}}{Sample size, variable count, and multivariate-normal distribution label.}
#'       \item{\code{fit.stats}}{Log-likelihood, AIC, BIC, entropy, and number of free parameters.}
#'       \item{\code{class.probs}}{Data frame containing profile labels, modal-assignment counts, and estimated profile proportions.}
#'       \item{\code{class.means}}{Profile-specific means for the first \code{I.max} variables.}
#'       \item{\code{convergence}}{Backend-specific algorithm, iteration, tolerance, initialization, and diagnostic information.}
#'       \item{\code{replication}}{Replication count and best BIC when replication selection applies; otherwise \code{NULL}.}
#'       \item{\code{digits}, \code{I.max.shown}, \code{total.vars}}{Formatting and truncation metadata.}
#'     }}
#'
#'   \item{\code{LCPA}}{A \code{summary.LCPA} object with:
#'     \describe{
#'       \item{\code{call}}{Original fitting call.}
#'       \item{\code{model.config}}{Analysis path, number of classes/profiles,
#'         model and three-step methods, Step 1 source, dependent-variable
#'         structure where applicable, and classification-error handling.}
#'       \item{\code{data.info}}{Sample size and number of response variables.}
#'       \item{\code{fit.stats}}{For XZ, log-likelihood, AIC, BIC, and number of free parameters.}
#'       \item{\code{class.probs}}{Data frame containing class probabilities, proportions, and modal-assignment frequencies.}
#'       \item{\code{coefficients}}{For XZ, the non-reference-class coefficient table with estimates, standard errors, 95 percent confidence limits, z statistics, and two-sided p-values.}
#'       \item{\code{dependent.variables}}{For ZY, fitted conditional
#'         distributions nested by model and dependent variable. Gaussian
#'         entries contain class/profile-specific means and variances, their
#'         standard errors and covariance matrices, and separate omnibus Wald
#'         tests. Categorical entries contain class/profile-specific category
#'         probabilities, standard errors, covariance matrices, and an omnibus
#'         test of equality of the conditional distributions.}
#'       \item{\code{covariates.names}, \code{ref.class}}{Displayed covariate names and the multinomial-logit reference class.}
#'       \item{\code{convergence}}{Overall and model-specific Step 3 convergence and iteration information.}
#'       \item{\code{digits}, \code{vars.to.show}, \code{total.vars}, \code{has.covariates}}{Formatting and covariate metadata.}
#'     }}
#'
#'   \item{\code{LTA}}{A \code{summary.LTA} object with:
#'     \describe{
#'       \item{\code{call}}{Original fitting call.}
#'       \item{\code{model.config}}{Number of time points and classes, model type, Step 1 source, reference class, covariate mode, classification-error handling, and transition mode.}
#'       \item{\code{data.info}}{Sample size, response-variable count, and number of time points.}
#'       \item{\code{fit.stats}}{For XZ, log-likelihood, AIC, BIC, and number of free parameters.}
#'       \item{\code{class.probs}}{Time-indexed data frames containing class probabilities, proportions, and modal-assignment frequencies.}
#'       \item{\code{initial.model}}{For XZ, the initial-status coefficient table, covariate names, and reference class.}
#'       \item{\code{transition.models}}{For XZ, time-invariant or time-indexed transition coefficient tables with origin class, destination class, covariate, estimate, standard error, confidence limits, z statistic, and p-value.}
#'       \item{\code{dependent.variables}}{For ZY, state- or path-specific
#'         Gaussian means and variances or categorical probabilities, together
#'         with their standard errors, covariance matrices, confidence-interval
#'         inputs, group masses, and omnibus Wald tests.}
#'       \item{\code{convergence}}{Overall and model-specific Step 3 convergence and iteration information.}
#'       \item{\code{digits}, \code{total.vars}, \code{covariates.time.cross}, \code{ref.class}}{Formatting, covariate, and reference-class metadata.}
#'     }}
#'
#'   \item{\code{sim.LCA}}{A \code{summary.sim.LCA} object with:
#'     \describe{
#'       \item{\code{call}}{Original simulation call.}
#'       \item{\code{config}}{Sample size, item count, class count, category counts, category-count uniformity, item quality, and generating distribution.}
#'       \item{\code{class.probs}}{True class probabilities and realized frequencies.}
#'       \item{\code{item.probs}}{True conditional response probabilities for the first \code{I.max} items.}
#'       \item{\code{digits}, \code{I.max.shown}, \code{total.vars}}{Formatting and truncation metadata.}
#'     }}
#'
#'   \item{\code{sim.LPA}}{A \code{summary.sim.LPA} object with:
#'     \describe{
#'       \item{\code{call}}{Original simulation call.}
#'       \item{\code{config}}{Sample size, variable count, profile count, constraint specification and description, and generating distribution.}
#'       \item{\code{class.probs}}{True profile probabilities and realized frequencies.}
#'       \item{\code{class.means}}{True profile means for the first \code{I.max} variables.}
#'       \item{\code{constraint}}{Expanded description of the covariance constraint.}
#'       \item{\code{digits}, \code{I.max.shown}, \code{total.vars}}{Formatting and truncation metadata.}
#'     }}
#'
#'   \item{\code{sim.LTA}}{A \code{summary.sim.LTA} object with:
#'     \describe{
#'       \item{\code{call}}{Original simulation call.}
#'       \item{\code{config}}{Sample size, variable count, class count, time points, model type, generating distribution, coefficient reference class, and LPA constraint when applicable.}
#'       \item{\code{class.probs}}{Time-indexed true class probabilities and realized frequencies.}
#'       \item{\code{item.probs}, \code{class.means}}{Truncated true measurement parameters for LCA or LPA simulations, respectively.}
#'       \item{\code{transition}}{Fixed-rate or covariate-dependent transition specification, including beta/gamma parameters and time indices when present.}
#'       \item{\code{covariates}}{Time-indexed covariate summaries containing minima, maxima, and means, or \code{NULL}.}
#'       \item{\code{digits}, \code{I.max.shown}, \code{L.max.shown}, \code{total.vars}, \code{total.classes}}{Formatting and truncation metadata.}
#'     }}
#'
#'   \item{\code{\link[LCPA:get.fit.index]{fit.index}}}{A
#'     \code{\link[LCPA]{summary.fit.index}} object with:
#'     \describe{
#'       \item{\code{call}}{Call that produced the fit-index object.}
#'       \item{\code{data.info}}{List containing the sample size \code{N}.}
#'       \item{\code{fit.table}}{Data frame with \code{Statistic}, \code{Value}, and \code{Description} columns for
#'         \code{npar}, \code{Log.Lik}, \code{-2LL}, AIC, BIC, SIC, CAIC, AWE, and SABIC.}
#'       \item{\code{digits}}{Requested numeric precision.}
#'     }}
#'
#'   \item{\code{compare.model}}{A \code{summary.compare.model} object with:
#'     \describe{
#'       \item{\code{call}}{Call that produced the model comparison.}
#'       \item{\code{data.info}}{Lists the named model-specific sample sizes and indicator counts and the two class counts.}
#'       \item{\code{fit.table}}{Side-by-side table of class count, parameter count, log-likelihood, -2LL, AIC, BIC, SIC, CAIC, AWE, and SABIC.}
#'       \item{\code{model.comparison}}{Data frame comparing class counts, parameter counts, diagonal average posterior probabilities, and entropy.}
#'       \item{\code{BF}, \code{BF.interpretation}}{Bayes factor computed from SIC and its evidence label.}
#'       \item{\code{LRT.table}}{Separate rows for the standard LRT, VLMR, adjusted LMR, and BLRT when available, with statistics, degrees of freedom, p-values, and significance symbols.}
#'       \item{\code{LRT.objects}}{Named list containing the unmodified hypothesis-test objects used to build \code{LRT.table}.}
#'       \item{\code{digits}}{Requested numeric precision.}
#'     }}
#'
#'   \item{\code{\link[LCPA:get.SE]{SE}}}{A \code{\link[LCPA]{summary.SE}} object with:
#'     \describe{
#'       \item{\code{call}}{Call that produced the standard-error object.}
#'       \item{\code{method}}{Selected \code{"Bootstrap"}, \code{"Obs"}, or \code{"Louis"} method.}
#'       \item{\code{diagnostics}}{Complete method-specific diagnostic list from \code{\link[LCPA]{get.SE}()}.}
#'       \item{\code{type}}{\code{"LCA"}, \code{"LPA"}, or \code{"Unknown"}, inferred from the standard-error components.}
#'       \item{\code{L}, \code{I}}{Number of classes/profiles and variables/items.}
#'       \item{\code{nonzero.counts}}{Counts of nonzero standard errors for \code{P.Z} and, as applicable, \code{par}, \code{means}, and \code{covs}.}
#'       \item{\code{total.P.Z}}{Total number of class-proportion standard errors.}
#'     }}
#' }
#'
#' @name summary
NULL

.summarize.Rmixmod <- function(object, arguments){
  control <- arguments$control.Rmixmod
  path <- if(is.null(control$path)) "LCPA" else control$path

  if(identical(path, "Rmixmod")){
    strategy <- methods::slot(object$model, "strategy")
    algorithm <- as.character(methods::slot(strategy, "algo"))
    iterations <- as.integer(methods::slot(strategy, "nbIterationInAlgo"))
    tol <- as.numeric(methods::slot(strategy, "epsilonInAlgo"))
    deterministic <- which(algorithm != "SEM")

    return(list(
      algorithm = paste0(paste(algorithm, collapse = " -> "), " (Rmixmod)"),
      iterations = paste(paste0(algorithm, "=", iterations), collapse = " -> "),
      tol = if(length(deterministic)) tol[deterministic] else NULL,
      initialization = paste0(
        methods::slot(strategy, "initMethod"), " (",
        as.integer(methods::slot(strategy, "nbTryInInit")), " tries; ",
        as.integer(methods::slot(strategy, "nbIterationInInit")), " iterations)"
      ),
      criterion = "log-likelihood",
      note = "Native Rmixmod strategy; starts, maxiter.warmup, and nrep are not used"
    ))
  }

  final.iterations <- control$maxiter
  if(is.null(final.iterations)) final.iterations <- 1000L
  list(
    algorithm = "Stochastic EM (SEM; Rmixmod)",
    iterations = paste0(arguments$maxiter.warmup, " warm-up + ", final.iterations, " final SEM"),
    initialization = paste0(arguments$starts, " starts -> ", arguments$nrep, " SEM runs"),
    criterion = "log-likelihood",
    note = "SEM uses a fixed iteration count; epsilon convergence is not defined"
  )
}

#' @describeIn summary Summary method for \code{LCA} objects
#' @importFrom utils tail
#' @exportS3Method summary LCA
summary.LCA <- function(object, digits = 4, I.max = 5, ...) {
  call_info <- object$call
  arguments <- object$arguments
  params <- object$params
  I <- ncol(arguments$response)
  N <- nrow(arguments$response)

  P.Z.Xn <- object$P.Z.Xn
  L <- length(params$P.Z)
  entropy.i <- -rowSums(P.Z.Xn * log(P.Z.Xn + 1e-10))
  entropy <- 1 - sum(entropy.i) / (N * log(L))

  poly.value <- sapply(object$probability, ncol)
  poly.value.uniform <- length(unique(poly.value)) == 1

  P.Z <- params$P.Z
  class.probs <- data.frame(
    Class = names(P.Z),
    Count = as.numeric(table(factor(object$Z, levels = seq_len(L)))),
    Proportion = sprintf("%.1f%%", P.Z * 100),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  items_to_show <- min(I.max, I)
  item.probs <- if (I > 0) {
    lapply(seq_len(items_to_show), function(i) {
      prob_mat <- round(object$probability[[i]], digits)
      dimnames(prob_mat) <- list(
        .latent.group.names(nrow(prob_mat), "LCA"),
        paste("Cat", 1:ncol(prob_mat))
      )
      prob_mat
    })
  } else {
    list()
  }

  convergence <- list()
  if (arguments$method == "EM") {
    iter.n <- length(object$Log.Lik.history) - 1
    tol <- if (!is.null(arguments$control.EM) && !is.null(arguments$control.EM$tol)) {
      arguments$control.EM$tol
    } else {
      1e-4
    }
    ll.initial <- object$Log.Lik.history[1]
    ll.final <- tail(object$Log.Lik.history, 1)
    ll.delta <- abs(ll.final - ll.initial)

    convergence <- list(
      algorithm = "Expectation-Maximization (EM)",
      iterations = iter.n,
      tol = tol,
      Log.Lik.change = ll.delta,
      Log.Lik.initial = ll.initial,
      Log.Lik.final = ll.final
    )
  } else if (arguments$method == "NNE") {
    iter.n <- length(object$Log.Lik.history) - 1
    patience.early <- if (!is.null(arguments$control.NNE) && !is.null(arguments$control.NNE$patience.early)) {
      arguments$control.NNE$patience.early
    } else {
      5
    }
    ll.initial <- object$Log.Lik.history[1]
    ll.final <- tail(object$Log.Lik.history, 1)
    ll.delta <- abs(ll.final - ll.initial)

    convergence <- list(
      algorithm = "Neural Network Estimation (NNE)",
      iterations = iter.n,
      patience = patience.early,
      Log.Lik.change = ll.delta,
      Log.Lik.initial = ll.initial,
      Log.Lik.final = ll.final,
      hardware = arguments$control.NNE$device
    )
  } else if (arguments$method == "Rmixmod") {
    convergence <- .summarize.Rmixmod(object, arguments)
  } else if (arguments$method == "flexmix") {
    convergence <- list(
      algorithm = "Stochastic EM (SEM; flexmix)",
      iterations = paste0(
        arguments$maxiter.warmup, " warm-up + ",
        arguments$control.flexmix$maxiter, " final SEM"
      ),
      tol = if (arguments$control.flexmix$tol > 0) {
        arguments$control.flexmix$tol
      } else {
        NULL
      },
      initialization = paste0(arguments$starts, " starts -> ", arguments$nrep, " SEM runs"),
      criterion = "log-likelihood",
      note = paste0(
        "Fixed iteration count; likelihood tolerance disabled; ",
        "one stochastic classification step per SEM iteration"
      )
    )
  } else if (arguments$method == "RMixtComp") {
    convergence <- list(
      algorithm = "Stochastic EM (SEM; RMixtComp)",
      iterations = paste0(
        arguments$control.RMixtComp$maxiter.burnin, " burn-in + ",
        arguments$control.RMixtComp$maxiter, " recorded SEM + ",
        arguments$control.RMixtComp$maxiter.gibbs.burnin, " Gibbs burn-in + ",
        arguments$control.RMixtComp$maxiter.gibbs, " recorded Gibbs"
      ),
      initialization = paste0(
        arguments$control.RMixtComp$n.init.per.class,
        " observations per class; ", arguments$control.RMixtComp$nrep,
        " native run(s)"
      ),
      criterion = arguments$control.RMixtComp$criterion,
      note = "Direct native RMixtComp flow; no LCPA warm-up or nrep promotion"
    )
  } else if (arguments$method == "Mplus") {
    convergence <- list(
      algorithm = "Mplus (External Estimation)",
      note = "Convergence diagnostics unavailable for external estimators"
    )
  }

  replication <- if (arguments$nrep > 1 &&
                     !arguments$method %in% c("Mplus", "RMixtComp") &&
                     !(arguments$method == "Rmixmod" &&
                       identical(arguments$control.Rmixmod$path, "Rmixmod"))) {
    list(
      nrep = arguments$nrep,
      best_BIC = object$best_BIC
    )
  } else {
    NULL
  }

  summary_obj <- list(
    call = call_info,
    model.config = list(
      L = arguments$L,
      method = arguments$method
    ),
    data.info = list(
      N = N,
      I = I,
      poly.value = poly.value,
      categories.uniform = poly.value.uniform
    ),
    fit.stats = list(
      LogLik = object$Log.Lik,
      AIC = object$AIC,
      BIC = object$BIC,
      entropy = entropy,
      npar = object$npar
    ),
    class.probs = class.probs,
    item.probs = item.probs,
    convergence = convergence,
    replication = replication,
    digits = digits,
    I.max.shown = items_to_show,
    total.items = I
  )

  class(summary_obj) <- "summary.LCA"
  summary_obj
}

#' @describeIn summary Summary method for \code{LPA} objects
#' @importFrom utils tail
#' @exportS3Method summary LPA
summary.LPA <- function(object, digits = 4, I.max = 5, ...) {
  call_info <- object$call
  arguments <- object$arguments
  params <- object$params
  N <- nrow(arguments$response)
  I <- ncol(arguments$response)
  L <- arguments$L

  item.names <- colnames(arguments$response)
  if(is.null(item.names)){
    item.names <- paste0("V", 1:I)
  }

  P.Z.Xn <- object$P.Z.Xn
  entropy.i <- -rowSums(P.Z.Xn * log(P.Z.Xn + 1e-10))
  entropy <- 1 - sum(entropy.i) / (N * log(L))

  P.Z <- params$P.Z
  class.probs <- data.frame(
    Profile = names(P.Z),
    Count = as.numeric(table(factor(object$Z, levels = seq_len(L)))),
    Proportion = sprintf("%.1f%%", P.Z * 100),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  vars_to_show <- min(I.max, I)
  class.means <- round(params$means[, seq_len(vars_to_show), drop = FALSE], digits)

  convergence <- list()
  if (arguments$method %in% c("EM", "NNE") && !is.null(object$Log.Lik.history)) {
    iter.n <- length(object$Log.Lik.history) - 1
    ll.initial <- object$Log.Lik.history[1]
    ll.final <- tail(object$Log.Lik.history, 1)
    ll.delta <- abs(ll.final - ll.initial)

    if (arguments$method == "EM") {
      tol <- if (!is.null(arguments$control.EM) && !is.null(arguments$control.EM$tol)) {
        arguments$control.EM$tol
      } else {
        1e-4
      }

      convergence <- list(
        algorithm = "Expectation-Maximization (EM)",
        iterations = iter.n,
        tol = tol,
        Log.Lik.change = ll.delta,
        Log.Lik.initial = ll.initial,
        Log.Lik.final = ll.final
      )
    } else if (arguments$method == "NNE") {
      patience.early <- if (!is.null(arguments$control.NNE) && !is.null(arguments$control.NNE$patience.early)) {
        arguments$control.NNE$patience.early
      } else {
        5
      }

      convergence <- list(
        algorithm = "Neural Network Estimation (NNE)",
        iterations = iter.n,
        patience = patience.early,
        Log.Lik.change = ll.delta,
        Log.Lik.initial = ll.initial,
        Log.Lik.final = ll.final,
        hardware = arguments$control.NNE$device
      )
    }
  } else if (arguments$method == "Rmixmod") {
    convergence <- .summarize.Rmixmod(object, arguments)
  } else if (arguments$method == "flexmix") {
    convergence <- list(
      algorithm = "Stochastic EM (SEM; flexmix)",
      iterations = paste0(
        arguments$maxiter.warmup, " warm-up + ",
        arguments$control.flexmix$maxiter, " final SEM"
      ),
      tol = if (arguments$control.flexmix$tol > 0) {
        arguments$control.flexmix$tol
      } else {
        NULL
      },
      initialization = paste0(arguments$starts, " starts -> ", arguments$nrep, " SEM runs"),
      criterion = "log-likelihood",
      note = paste0(
        "Fixed iteration count; likelihood tolerance disabled; ",
        "one stochastic classification step per SEM iteration"
      )
    )
  } else if (arguments$method == "RMixtComp") {
    convergence <- list(
      algorithm = "Stochastic EM (SEM; RMixtComp)",
      iterations = paste0(
        arguments$control.RMixtComp$maxiter.burnin, " burn-in + ",
        arguments$control.RMixtComp$maxiter, " recorded SEM + ",
        arguments$control.RMixtComp$maxiter.gibbs.burnin, " Gibbs burn-in + ",
        arguments$control.RMixtComp$maxiter.gibbs, " recorded Gibbs"
      ),
      initialization = paste0(
        arguments$control.RMixtComp$n.init.per.class,
        " observations per class; ", arguments$control.RMixtComp$nrep,
        " native run(s)"
      ),
      criterion = arguments$control.RMixtComp$criterion,
      note = "Direct native RMixtComp flow; no LCPA warm-up or nrep promotion"
    )
  } else if (arguments$method == "Mplus") {
    convergence <- list(
      algorithm = "Mplus (External Estimation)",
      note = "Convergence diagnostics unavailable for external estimators"
    )
  }

  # Replication info
  replication <- if (arguments$nrep > 1 &&
                     !arguments$method %in% c("Mplus", "RMixtComp") &&
                     !(arguments$method == "Rmixmod" &&
                       identical(arguments$control.Rmixmod$path, "Rmixmod"))) {
    finite.Log.Lik <- object$Log.Lik.nrep[is.finite(object$Log.Lik.nrep)]
    list(
      nrep = arguments$nrep,
      best_BIC = if(length(finite.Log.Lik) > 0) -2 * max(finite.Log.Lik) + log(N) * object$npar else NA_real_
    )
  } else {
    NULL
  }

  # Covariance structure description
  if(any(arguments$constraint %in% c("VV", "VE", "EV", "EE", "V0", "E0", "UV", "UE"))){
    constraint <- switch(arguments$constraint,
                            VV = "   Free variance, free covariance",
                            VE = "   Free variance, shared covariance",
                            EV = "   Shared variance, free covariance",
                            EE = "   Shared variance, shared covariance",
                            V0 = "   Free variance, covariance = 0",
                            E0 = "   Shared variance, covariance = 0",
                            UV  = "   Univariate free variance",
                            UE  = "   Univariate shared variance")
  }else {
    variance_constraints <- c()
    covariance_constraints <- c()
    for (const in arguments$constraint) {
      if (length(const) == 2 && is.numeric(const)) {
        var1 <- const[1]
        var2 <- const[2]
        if (var1 != var2) {
          cov_desc <- paste0("    ", item.names[var1], " & ", item.names[var2], " share covariance")
          covariance_constraints <- c(covariance_constraints, cov_desc)
        }
        else {
          var_desc <- paste0("    ", item.names[var1], " shares variance")
          variance_constraints <- c(variance_constraints, var_desc)
        }
      } else {
        warning(paste("Invalid constraint format:", paste(const, collapse = ",")))
      }
    }
    all_constraints <- c(
      if (length(variance_constraints) > 0) variance_constraints else "No shared variance constraints",
      if (length(covariance_constraints) > 0) covariance_constraints else "No shared covariance constraints"
    )
    constraint <- paste(all_constraints, collapse = "\n")
  }

  # Build summary object
  summary_obj <- list(
    call = call_info,
    model.config = list(
      L = L,
      constraint = arguments$constraint,
      constraint.description = constraint,
      method = arguments$method
    ),
    data.info = list(
      N = N,
      I = I,
      distribution = "Multivariate Normal"
    ),
    fit.stats = list(
      LogLik = object$Log.Lik,
      AIC = object$AIC,
      BIC = object$BIC,
      entropy = entropy,
      npar = object$npar
    ),
    class.probs = class.probs,
    class.means = class.means,
    convergence = convergence,
    replication = replication,
    digits = digits,
    I.max.shown = vars_to_show,
    total.vars = I
  )

  class(summary_obj) <- "summary.LPA"
  summary_obj
}

.summarize.ZY <- function(object, digits, longitudinal = FALSE) {
  arguments <- object$arguments
  times <- length(object$P.Z.Xns)
  L <- ncol(object$P.Z.Xns[[1L]])
  class.probs <- lapply(seq_len(times), function(t) {
    .latent.group.columns(data.frame(
      Class = .latent.group.names(L, arguments$type.model),
      Probability = round(object$P.Zs[[t]], digits),
      Proportion = sprintf("%.1f%%", object$P.Zs[[t]] * 100),
      Frequency = as.vector(table(factor(object$Zs[[t]], levels = seq_len(L))))
    ), arguments$type.model)
  })
  names(class.probs) <- paste0("Time ", seq_len(times))
  dependent.model.info <- do.call(rbind, lapply(
    names(object$dependent.variables), function(group.name) {
      group <- object$dependent.variables[[group.name]]
      do.call(rbind, lapply(names(group), function(dependent.variable.name) {
        model <- group[[dependent.variable.name]]
        data.frame(
          Model = group.name,
          Dependent.Variable = dependent.variable.name,
          Family = model$family,
          Observations = model$observations,
          Omitted = model$omitted,
          Converged = model$converged,
          Iterations = model$iterations,
          stringsAsFactors = FALSE
        )
      }))
    }
  ))
  list(
    call = object$call,
    type.analysis = "ZY",
    model.config = list(
      times = times,
      L = L,
      type = ifelse(arguments$type.model == "LCA",
                    "Latent Class Analysis (categorical)",
                    "Latent Profile Analysis (continuous)"),
      type.model = arguments$type.model,
      method.3step = arguments$method.3step,
      method.model = arguments$method.model,
      method.regression = arguments$method.regression,
      method.SE = arguments$method.SE,
      dependent.variable.structure = if(longitudinal) {
        arguments$dependent.variable.structure
      } else {
        "State"
      },
      dependent.variable.time = if(longitudinal) arguments$dependent.variable.time else 1L,
      dependent.variable.time.cross = if(longitudinal) arguments$dependent.variable.time.cross else FALSE,
      step1.source = if (!is.null(arguments$control.model$params)) {
        "User-supplied fixed parameters"
      } else if (longitudinal && isTRUE(arguments$step1.pool)) {
        "Pooled responses across all time points"
      } else if (longitudinal) {
        "Time point 1 (default)"
      } else {
        "Cross-sectional response data"
      },
      npar = object$npar,
      CEP.handling = ifelse(arguments$CEP.error,
                            ifelse(longitudinal && arguments$CEP.time.cross,
                                   "Classification error correction (time-invariant CEP)",
                                   "Classification error correction"),
                            "No classification error correction (naive modal assignment)")
    ),
    data.info = list(
      N = nrow(object$P.Z.Xns[[1L]]),
      variables = if(longitudinal) ncol(arguments$responses[[1L]]) else ncol(arguments$response),
      times = times,
      dependent.variable.models = nrow(dependent.model.info)
    ),
    class.probs = class.probs,
    dependent.variables = object$dependent.variables,
    latent.paths = object$latent.paths,
    SE.diagnostics = object$SE.diagnostics,
    diagnostics = object$diagnostics,
    convergence = list(
      converged = object$converged,
      models = dependent.model.info
    ),
    digits = digits
  )
}

#' @describeIn summary Summary method for \code{LTA} objects
#' @exportS3Method summary LTA
summary.LTA <- function(object, digits = 4, ...) {
  if(identical(object$type.analysis, "ZY")){
    summary_obj <- .summarize.ZY(object, digits, longitudinal = TRUE)
    class(summary_obj) <- "summary.LTA"
    return(summary_obj)
  }
  arguments <- object$arguments
  times <- length(object$P.Zs)
  L <- ncol(object$beta)
  ref.class <- arguments$ref.class
  type <- arguments$type.model
  N <- nrow(arguments$responses[[1]])
  I <- ncol(arguments$responses[[1]])
  covariates.time.cross <- arguments$covariates.time.cross

  # Model configuration
  model_config <- list(
    times = times,
    L = L,
    type.model = type,
    type = ifelse(type == "LCA", "Latent Class Analysis (categorical)",
                        "Latent Profile Analysis (continuous)"),
    method.3step = arguments$method.3step,
    method.model = arguments$method.model,
    method.regression = arguments$method.regression,
    method.SE = arguments$method.SE,
    step1.source = if (!is.null(arguments$control.model$params)) {
      "User-supplied fixed parameters"
    } else if (isTRUE(arguments$step1.pool)) {
      "Pooled responses across all time points"
    } else {
      "Time point 1 (default)"
    },
    ref.class = ref.class,
    covariates.mode = if (!is.null(arguments$covariates)) {
      if (covariates.time.cross) "Time-invariant covariates" else "Time-varying covariates"
    } else "No covariates",
    CEP.handling = ifelse(arguments$CEP.error,
                          ifelse(arguments$CEP.time.cross,
                                 "Classification error correction (time-invariant CEP)",
                                 "Classification error correction (time-varying CEP)"),
                          "No classification error correction (naive modal assignment)"),
    transition.mode = ifelse(covariates.time.cross,
                             "Time-invariant transition effects (coefficients constant across time)",
                             "Time-varying transition effects (coefficients differ by time point)")
  )

  # Class probabilities over time
  class_probs <- lapply(1:times, function(t) {
    .latent.group.columns(data.frame(
      Class = .latent.group.names(L, type),
      Probability = round(object$P.Zs[[t]], digits),
      Proportion = sprintf("%.1f%%", object$P.Zs[[t]] * 100),
      Frequency = as.vector(table(factor(object$Zs[[t]], levels = 1:L)))
    ), type)
  })
  names(class_probs) <- paste0("Time ", 1:times)

  # Extract covariate names (unified if time-cross)
  cov_names_list <- vector("list", times)
  for (t in 1:times) {
    cov_mat <- arguments$covariates[[t]]
    cov_names <- colnames(cov_mat)
    if (is.null(cov_names) || any(cov_names == "")) {
      cov_names <- character(ncol(cov_mat))
      cov_names[1] <- "Intercept"
      if (ncol(cov_mat) > 1) {
        cov_names[2:ncol(cov_mat)] <- paste0("Cov", 1:(ncol(cov_mat)-1))
      }
    }
    cov_names_list[[t]] <- cov_names
  }

  if (covariates.time.cross) {
    # Use first time point's covariate names for all
    cov_names_list <- replicate(times, cov_names_list[[1]], simplify = FALSE)
  }

  # Initial state coefficients table (always needed)
  cov_names_initial <- cov_names_list[[1]]
  vars_to_show_initial <- length(cov_names_initial)

  initial_coef <- data.frame(
    Class = character(),
    Covariate = character(),
    Estimate = numeric(),
    Std.Error = numeric(),
    lower.95 = numeric(),
    upper.95 = numeric(),
    z.value = numeric(),
    p.value = numeric(),
    stringsAsFactors = FALSE
  )

  non_ref_classes <- setdiff(1:L, ref.class)
  for (cls in non_ref_classes) {
    for (j in 1:vars_to_show_initial) {
      coef_val <- object$beta[j, cls]
      se_val <- if (!is.null(object$beta.se)) object$beta.se[j, cls] else NA
      z_val <- if (!is.null(object$beta.Z.sta)) object$beta.Z.sta[j, cls] else NA
      p_val <- if (!is.null(object$beta.p.value.tail2)) object$beta.p.value.tail2[j, cls] else NA

      # Calculate 95% confidence intervals safely
      lower.95 <- if (!is.na(se_val) && is.finite(se_val)) coef_val - 1.96 * se_val else NA
      upper.95 <- if (!is.na(se_val) && is.finite(se_val)) coef_val + 1.96 * se_val else NA

      initial_coef <- rbind(initial_coef, data.frame(
        Class = .latent.group.names(L, type)[cls],
        Covariate = cov_names_initial[j],
        Estimate = round(coef_val, digits),
        Std.Error = ifelse(is.na(se_val) || !is.finite(se_val), NA, round(se_val, digits)),
        lower.95 = ifelse(is.na(se_val) || !is.finite(se_val), NA, round(lower.95, digits)),
        upper.95 = ifelse(is.na(se_val) || !is.finite(se_val), NA, round(upper.95, digits)),
        z.value = ifelse(is.na(z_val) || !is.finite(z_val), NA, round(z_val, digits)),
        p.value = ifelse(is.na(p_val) || !is.finite(p_val), NA, sprintf("%.4f", p_val))
      ))
    }
  }
  initial_coef <- .latent.group.columns(initial_coef, type)

  # Transition coefficients: time-invariant vs time-varying
  if (covariates.time.cross) {
    # Only one transition model needed (all time points share same coefficients)
    trans_df <- data.frame(
      From.Class = character(),
      To.Class = character(),
      Covariate = character(),
      Estimate = numeric(),
      Std.Error = numeric(),
      lower.95 = numeric(),
      upper.95 = numeric(),
      z.value = numeric(),
      p.value = numeric(),
      stringsAsFactors = FALSE
    )

    cov_names_trans <- cov_names_list[[2]]  # Same for all time points
    vars_to_show_trans <- length(cov_names_trans)

    # CRITICAL FIX: Correct reference class handling for transitions
    # For each FROM class, all transitions are relative to the SAME reference destination class
    for (from_cls in 1:L) {
      for (to_cls in non_ref_classes) {
        for (j in 1:vars_to_show_trans) {
          # Use first transition (t1->t2) as representative
          coef_val <- object$gamma[[1]][[from_cls]][[to_cls]][j]
          se_val <- if (!is.null(object$gamma.se)) object$gamma.se[[1]][[from_cls]][[to_cls]][j] else NA
          z_val <- if (!is.null(object$gamma.Z.sta)) object$gamma.Z.sta[[1]][[from_cls]][[to_cls]][j] else NA
          p_val <- if (!is.null(object$gamma.p.value.tail2)) object$gamma.p.value.tail2[[1]][[from_cls]][[to_cls]][j] else NA

          # Calculate 95% confidence intervals safely
          lower.95 <- if (!is.na(se_val) && is.finite(se_val)) coef_val - 1.96 * se_val else NA
          upper.95 <- if (!is.na(se_val) && is.finite(se_val)) coef_val + 1.96 * se_val else NA

          trans_df <- rbind(trans_df, data.frame(
            From.Class = .latent.group.names(L, type)[from_cls],
            To.Class = .latent.group.names(L, type)[to_cls],
            Covariate = cov_names_trans[j],
            Estimate = round(coef_val, digits),
            Std.Error = ifelse(is.na(se_val) || !is.finite(se_val), NA, round(se_val, digits)),
            lower.95 = ifelse(is.na(se_val) || !is.finite(se_val), NA, round(lower.95, digits)),
            upper.95 = ifelse(is.na(se_val) || !is.finite(se_val), NA, round(upper.95, digits)),
            z.value = ifelse(is.na(z_val) || !is.finite(z_val), NA, round(z_val, digits)),
            p.value = ifelse(is.na(p_val) || !is.finite(p_val), NA, sprintf("%.4f", p_val))
          ))
        }
      }
    }

    trans_df <- .latent.group.columns(trans_df, type)

    # Add reference class explanation to dataframe
    attr(trans_df, "ref_class_note") <- paste0(
      "All transitions are relative to reference destination ",
      .latent.group.names(L, type)[ref.class]
    )

    # Single transition model with special label
    transition_coefs <- list(`Time-invariant effects` = trans_df)
  } else {
    # Time-varying effects: show each transition phase
    transition_coefs <- list()
    for (t in 1:(times-1)) {
      cov_names_trans <- cov_names_list[[t+1]]
      vars_to_show_trans <- length(cov_names_trans)

      trans_df <- data.frame(
        From.Class = character(),
        To.Class = character(),
        Covariate = character(),
        Estimate = numeric(),
        Std.Error = numeric(),
        lower.95 = numeric(),
        upper.95 = numeric(),
        z.value = numeric(),
        p.value = numeric(),
        stringsAsFactors = FALSE
      )

      # CRITICAL FIX: Correct reference class handling for transitions
      for (from_cls in 1:L) {
        for (to_cls in non_ref_classes) {
          for (j in 1:vars_to_show_trans) {
            coef_val <- object$gamma[[t]][[from_cls]][[to_cls]][j]
            se_val <- if (!is.null(object$gamma.se)) object$gamma.se[[t]][[from_cls]][[to_cls]][j] else NA
            z_val <- if (!is.null(object$gamma.Z.sta)) object$gamma.Z.sta[[t]][[from_cls]][[to_cls]][j] else NA
            p_val <- if (!is.null(object$gamma.p.value.tail2)) object$gamma.p.value.tail2[[t]][[from_cls]][[to_cls]][j] else NA

            # Calculate 95% confidence intervals safely
            lower.95 <- if (!is.na(se_val) && is.finite(se_val)) coef_val - 1.96 * se_val else NA
            upper.95 <- if (!is.na(se_val) && is.finite(se_val)) coef_val + 1.96 * se_val else NA

            trans_df <- rbind(trans_df, data.frame(
              From.Class = .latent.group.names(L, type)[from_cls],
              To.Class = .latent.group.names(L, type)[to_cls],
              Covariate = cov_names_trans[j],
              Estimate = round(coef_val, digits),
              Std.Error = ifelse(is.na(se_val) || !is.finite(se_val), NA, round(se_val, digits)),
              lower.95 = ifelse(is.na(se_val) || !is.finite(se_val), NA, round(lower.95, digits)),
              upper.95 = ifelse(is.na(se_val) || !is.finite(se_val), NA, round(upper.95, digits)),
              z.value = ifelse(is.na(z_val) || !is.finite(z_val), NA, round(z_val, digits)),
              p.value = ifelse(is.na(p_val) || !is.finite(p_val), NA, sprintf("%.4f", p_val))
            ))
          }
        }
      }
      trans_df <- .latent.group.columns(trans_df, type)

      # Add reference class explanation
      attr(trans_df, "ref_class_note") <- paste0(
        "All transitions are relative to reference destination ",
        .latent.group.names(L, type)[ref.class]
      )
      transition_coefs[[t]] <- trans_df
    }
    names(transition_coefs) <- paste0("Time ", 1:(times-1), " -> Time ", 2:times)
  }

  # Convergence information
  log_lik_at_start <- if (length(object$Log.Lik.history) > 0) object$Log.Lik.history[1L] else NA
  log_lik_at_end <- object$Log.Lik
  Log.Lik.change <- if (!is.na(log_lik_at_start)) round(log_lik_at_end - log_lik_at_start, digits) else NA

  convergence_info <- list(
    iterations = object$iterations,
    converged = object$converged,
    note = sprintf("  Log-likelihood change: |%.2f - %.2f| = %.4f\n",
                           log_lik_at_start, log_lik_at_end, Log.Lik.change)
  )

  # Build summary object
  summary_obj <- list(
    call = object$call,
    model.config = model_config,
    data.info = list(
      N = N,
      variables = I,
      times = times
    ),
    fit.stats = list(
      LogLik = object$Log.Lik,
      AIC = object$AIC,
      BIC = object$BIC,
      npar = object$npar
    ),
    class.probs = class_probs,
    initial.model = list(
      coefficients = initial_coef,
      covariates.names = cov_names_initial,
      ref.class = ref.class
    ),
    transition.models = transition_coefs,
    convergence = convergence_info,
    digits = digits,
    total.vars = if (!is.null(arguments$covariates)) max(sapply(arguments$covariates, ncol)) else 1,
    covariates.time.cross = covariates.time.cross,
    ref.class = ref.class  # Store reference class at top level for easy access
  )

  class(summary_obj) <- "summary.LTA"
  summary_obj
}

#' @describeIn summary Summary method for \code{LCPA} objects
#' @exportS3Method summary LCPA
summary.LCPA <- function(object, digits = 4, ...) {
  if(identical(object$type.analysis, "ZY")){
    summary_obj <- .summarize.ZY(object, digits, longitudinal = FALSE)
    class(summary_obj) <- "summary.LCPA"
    return(summary_obj)
  }
  arguments <- object$arguments
  L <- ncol(object$beta)
  ref.class <- arguments$ref.class
  type <- arguments$type.model
  N <- nrow(object$P.Z.Xn)  # Single time point
  I <- ncol(arguments$response)

  # Model configuration
  model_config <- list(
    L = L,
    type.model = type,
    type = ifelse(type == "LCA", "Latent Class Analysis (categorical)",
                        "Latent Profile Analysis (continuous)"),
    method.3step = arguments$method.3step,
    method.model = arguments$method.model,
    method.regression = arguments$method.regression,
    method.SE = arguments$method.SE,
    ref.class = ref.class,
    covariates.mode = if (!is.null(arguments$covariates)) {
      "Covariates included in class membership model"
    } else "No covariates (intercept-only model)",
    CEP.handling = ifelse(arguments$CEP.error,
                          "Classification error correction applied",
                          "No classification error correction (naive modal assignment)")
  )

  # Class probabilities (single time point)
  class_probs_df <- .latent.group.columns(data.frame(
    Class = .latent.group.names(L, type),
    Probability = round(object$P.Z, digits),
    Proportion = sprintf("%.1f%%", object$P.Z * 100),
    Frequency = as.vector(table(factor(object$Z, levels = 1:L)))
  ), type)

  # Extract covariate names
  cov_mat <- arguments$covariates
  if (is.null(cov_mat)) {
    cov_names <- "Intercept"
  } else {
    cov_names <- colnames(cov_mat)
    if (is.null(cov_names) || any(cov_names == "")) {
      cov_names <- character(ncol(cov_mat))
      cov_names[1] <- "Intercept"
      if (ncol(cov_mat) > 1) {
        cov_names[2:ncol(cov_mat)] <- paste0("Cov", 1:(ncol(cov_mat)-1))
      }
    }
  }

  vars_to_show <- length(cov_names)

  # Coefficients table (beta)
  coef_df <- data.frame(
    Class = character(),
    Covariate = character(),
    Estimate = numeric(),
    Std.Error = numeric(),
    lower.95 = numeric(),
    upper.95 = numeric(),
    z.value = numeric(),
    p.value = numeric(),
    stringsAsFactors = FALSE
  )

  non_ref_classes <- setdiff(1:L, ref.class)
  for (cls in non_ref_classes) {
    for (j in 1:vars_to_show) {
      coef_val <- object$beta[j, cls]
      se_val <- if (!is.null(object$beta.se)) object$beta.se[j, cls] else NA
      z_val <- if (!is.null(object$beta.Z.sta)) object$beta.Z.sta[j, cls] else NA
      p_val <- if (!is.null(object$beta.p.value.tail2)) object$beta.p.value.tail2[j, cls] else NA

      # Calculate 95% confidence intervals safely
      lower.95 <- if (!is.na(se_val) && is.finite(se_val)) coef_val - 1.96 * se_val else NA
      upper.95 <- if (!is.na(se_val) && is.finite(se_val)) coef_val + 1.96 * se_val else NA

      coef_df <- rbind(coef_df, data.frame(
        Class = .latent.group.names(L, type)[cls],
        Covariate = cov_names[j],
        Estimate = round(coef_val, digits),
        Std.Error = ifelse(is.na(se_val) || !is.finite(se_val), NA, round(se_val, digits)),
        lower.95 = ifelse(is.na(se_val) || !is.finite(se_val), NA, round(lower.95, digits)),
        upper.95 = ifelse(is.na(se_val) || !is.finite(se_val), NA, round(upper.95, digits)),
        z.value = ifelse(is.na(z_val) || !is.finite(z_val), NA, round(z_val, digits)),
        p.value = ifelse(is.na(p_val) || !is.finite(p_val), NA, sprintf("%.4f", p_val))
      ))
    }
  }
  coef_df <- .latent.group.columns(coef_df, type)

  # Convergence information
  log_lik_at_start <- if (length(object$Log.Lik.history) > 0) object$Log.Lik.history[1L] else NA
  log_lik_at_end <- object$Log.Lik
  Log.Lik.change <- if (!is.na(log_lik_at_start)) round(log_lik_at_end - log_lik_at_start, digits) else NA

  convergence_info <- list(
    iterations = object$iterations,
    converged = object$converged,
    note = if (!is.na(Log.Lik.change)) {
      sprintf("  Log-likelihood change: |%.2f - %.2f| = %.4f\n",
              log_lik_at_start, log_lik_at_end, Log.Lik.change)
    } else {
      "  Log-likelihood history unavailable\n"
    }
  )

  # Build summary object
  summary_obj <- list(
    call = object$call,
    model.config = model_config,
    data.info = list(
      N = N,
      variables = I
    ),
    fit.stats = list(
      LogLik = object$Log.Lik,
      AIC = object$AIC,
      BIC = object$BIC,
      npar = object$npar
    ),
    class.probs = class_probs_df,
    coefficients = coef_df,
    covariates.names = cov_names,
    ref.class = ref.class,
    convergence = convergence_info,
    digits = digits,
    vars.to.show = vars_to_show,
    total.vars = length(cov_names),
    has.covariates = !is.null(arguments$covariates)
  )

  class(summary_obj) <- "summary.LCPA"
  summary_obj
}

#' @describeIn summary Summary method for \code{sim.LCA} objects
#' @exportS3Method summary sim.LCA
summary.sim.LCA <- function(object, digits = 4, I.max = 5, ...) {
  N <- nrow(object$response)
  I <- ncol(object$response)
  L <- length(object$P.Z)

  vars_to_show <- min(I.max, I)
  item.probs <- vector("list", vars_to_show)

  for (i in 1:vars_to_show) {
    prob_mat <- matrix(NA, nrow = L, ncol = object$poly.value[i])
    for (l in 1:L) {
      prob_mat[l, ] <- object$par[l, i, 1:object$poly.value[i]]
    }
    dimnames(prob_mat) <- list(
      .latent.group.names(L, "LCA"),
      paste0("Cat", 0:(object$poly.value[i]-1))
    )
    item.probs[[i]] <- round(prob_mat, digits)
  }

  summary_obj <- list(
    call = object$call,
    config = list(
      N = N,
      I = I,
      L = L,
      poly.value = object$poly.value,
      categories.uniform = length(unique(object$poly.value)) == 1,
      IQ = object$arguments$IQ,
      distribution = object$arguments$distribution
    ),
    class.probs = data.frame(
      Class = .latent.group.names(L, "LCA"),
      Probability = round(object$P.Z, digits),
      Frequency = as.vector(table(factor(object$Z, levels = seq_len(L))))
    ),
    item.probs = item.probs,
    I.max.shown = vars_to_show,
    total.vars = I,
    digits = digits
  )

  class(summary_obj) <- "summary.sim.LCA"
  return(summary_obj)
}

#' @describeIn summary Summary method for \code{sim.LPA} objects
#' @exportS3Method summary sim.LPA
summary.sim.LPA <- function(object, digits = 4, I.max = 5, ...) {
  N <- nrow(object$response)
  I <- ncol(object$response)
  L <- length(object$P.Z)
  var_names <- colnames(object$response)

  vars_to_show <- min(I.max, I)
  selected_vars <- var_names[1:vars_to_show]

  class.means <- object$means[, selected_vars, drop = FALSE]

  constraint_desc <- describe_constraint(object$constraint, I)

  summary_obj <- list(
    call = object$call,
    config = list(
      N = N,
      I = I,
      L = L,
      constraint = object$constraint,
      constraint_desc = constraint_desc,
      distribution = object$arguments$distribution
    ),
    class.probs = data.frame(
      Profile = .latent.group.names(L, "LPA"),
      Probability = object$P.Z,
      Frequency = as.vector(table(factor(object$Z, levels = seq_len(L))))
    ),
    class.means = class.means,
    constraint = constraint_desc$details,
    I.max.shown = vars_to_show,
    total.vars = I,
    digits = digits
  )

  class(summary_obj) <- "summary.sim.LPA"
  return(summary_obj)
}

describe_constraint <- function(constraint, I) {
  if (is.character(constraint)) {
    desc_map <- c(
      VV = "Varying full covariance matrices (heterogeneous)",
      VE = "Varying variances, equal correlations",
      EV = "Equal variances, varying covariances",
      EE = "Equal full covariance matrices (homogeneous)",
      V0 = "Varying variances, zero covariances (diagonal)",
      E0 = "Equal variances, zero covariances (diagonal)",
      UV  = "Varying variances (univariate)",
      UE  = "Equal variances (univariate)"
    )

    description <- desc_map[constraint]
    if (is.null(description)) {
      description <- sprintf("Custom constraint: %s", constraint)
      details <- "See object$constraint for specification"
      example <- NULL
    } else {
      if (constraint %in% c("VV", "VE", "EV", "EE", "V0", "E0") && I > 1) {
        details <- switch(constraint,
                          VV = "All covariance parameters vary freely across profiles",
                          VE = "Diagonal elements vary across profiles; off-diagonal elements share correlation structure",
                          EV = "Diagonal elements equal across profiles; off-diagonal elements vary freely",
                          EE = "All covariance parameters equal across profiles",
                          V0 = "All off-diagonal elements constrained to zero; diagonal elements vary freely",
                          E0 = "All off-diagonal elements constrained to zero; diagonal elements equal across profiles"
        )
      } else if (constraint %in% c("UV", "UE") && I == 1) {
        details <- switch(constraint,
                          UV = "Variance parameters vary freely across profiles",
                          UE = "Variance parameter equal across all profiles"
        )
      } else {
        details <- "Constraint applied appropriately for data dimensionality"
      }
      example <- constraint
    }
  } else if (is.list(constraint)) {
    description <- "Custom constraints"

    # 生成详细描述
    var_constraints <- character(0)
    cov_constraints <- character(0)

    for (con in constraint) {
      i <- con[1]
      j <- con[2]

      if (i == j) {
        var_constraints <- c(var_constraints,
                             sprintf("Variance of %s equal across profiles",
                                     if(I==1) "variable" else paste0("UV", i)))
      } else {
        cov_constraints <- c(cov_constraints,
                             sprintf("Covariance between V%d and V%d equal across profiles", i, j))
      }
    }

    details <- character(0)
    if (length(var_constraints) > 0) {
      details <- c(details, "Shared variance constraints:")
      details <- c(details, paste0("  - ", var_constraints))
    } else {
      details <- c(details, "No shared variance constraints")
    }

    if (length(cov_constraints) > 0) {
      details <- c(details, "Shared covariance constraints:")
      details <- c(details, paste0("  - ", cov_constraints))
    } else {
      details <- c(details, "No shared covariance constraints")
    }

    details <- paste(details, collapse = "\n")
    example <- "custom"
  } else {
    description <- "Unknown constraint type"
    details <- "Check object$constraint for details"
    example <- NULL
  }

  list(description = description, details = details, example = example)
}

#' @describeIn summary Summary method for \code{sim.LTA} objects
#' @exportS3Method summary sim.LTA
summary.sim.LTA <- function(object, digits = 4, I.max = 5, L.max = 5, ...) {
  # Extract basic information
  N <- nrow(object$responses[[1]])
  I <- ncol(object$responses[[1]])
  L <- length(object$P.Zs[[1]])
  times <- length(object$responses)
  type <- if (!is.null(object$par)) "LCA" else "LPA"

  # Determine how many classes and variables to show
  classes_to_show <- min(L.max, L)
  vars_to_show <- min(I.max, I)

  # Class probabilities at each time point
  class_probs <- lapply(1:times, function(t) {
    .latent.group.columns(data.frame(
      Class = .latent.group.names(L, type),
      Probability = round(object$P.Zs[[t]], digits),
      Frequency = as.vector(table(factor(object$Zs[[t]], levels = seq_len(L))))
    ), type)
  })
  names(class_probs) <- paste0("Time ", 1:times)

  # Item probabilities or class means depending on type
  if (type == "LCA") {
    item_probs <- vector("list", times)
    for (t in 1:times) {
      item_probs[[t]] <- vector("list", vars_to_show)
      for (i in 1:vars_to_show) {
        prob_mat <- matrix(NA, nrow = classes_to_show, ncol = object$poly.value[i])
        for (l in 1:classes_to_show) {
          prob_mat[l, ] <- object$par[l, i, 1:object$poly.value[i]]
        }
        dimnames(prob_mat) <- list(
          .latent.group.names(classes_to_show, type),
          paste0("Cat", 0:(object$poly.value[i]-1))
        )
        item_probs[[t]][[i]] <- round(prob_mat, digits)
      }
      names(item_probs[[t]]) <- paste0("Item", 1:vars_to_show)
    }
    names(item_probs) <- paste0("Time ", 1:times)

    class_means <- NULL
  } else { # LPA
    class_means <- lapply(1:times, function(t) {
      # Use consistent class ordering across time points
      means_mat <- object$means[1:classes_to_show, 1:vars_to_show, drop = FALSE]
      round(means_mat, digits)
    })
    names(class_means) <- paste0("Time ", 1:times)

    item_probs <- NULL
  }

  # Transition information
  transition_info <- list()
  if (!is.null(object$rate)) {
    # Non-covariate mode: fixed transition probabilities
    transition_info$mode <- "fixed"
    transition_info$rate <- lapply(1:(times-1), function(t) {
      rate_mat <- round(object$rate[[t]], digits)
      dimnames(rate_mat) <- list(
        paste("From", .latent.group.names(L, type)),
        paste("To", .latent.group.names(L, type))
      )
      rate_mat
    })
    # Store time point indices directly for printing
    transition_info$times <- data.frame(
      from = 1:(times-1),
      to = 2:times
    )
  } else if (!is.null(object$beta) && !is.null(object$gamma)) {
    # Covariate mode: regression coefficients
    transition_info$mode <- "covariate"
    transition_info$beta <- round(object$beta, digits)
    dimnames(transition_info$beta) <- list(
      colnames(object$covariates[[1]]),
      .latent.group.names(ncol(object$beta), type)
    )

    gamma_display <- vector("list", times-1)
    for (t in 1:(times-1)) {
      # Get number of covariates from the first gamma coefficient
      num_covariates <- length(object$gamma[[t]][[1]][[1]])
      gamma_t <- vector("list", L)
      for (l in 1:L) {
        # Create matrix with proper dimensions
        gamma_mat <- matrix(NA, nrow = num_covariates, ncol = L)
        for (k in 1:L) {
          gamma_mat[, k] <- object$gamma[[t]][[l]][[k]]
        }
        gamma_t[[l]] <- round(gamma_mat, digits)
        dimnames(gamma_t[[l]]) <- list(
          colnames(object$covariates[[t + 1L]]),
          paste("To", .latent.group.names(L, type))
        )
      }
      names(gamma_t) <- paste("From", .latent.group.names(L, type))
      gamma_display[[t]] <- gamma_t
    }
    names(gamma_display) <- paste0("T", 1:(times-1), "->T", 2:times)
    transition_info$gamma <- gamma_display

    # Store time point indices directly for printing
    transition_info$times <- data.frame(
      from = 1:(times-1),
      to = 2:times
    )
  }

  # Create summary object
  summary_obj <- list(
    call = object$call,
    config = list(
      N = N,
      I = I,
      L = L,
      times = times,
      type = type,
      distribution = object$arguments$distribution,
      ref.class = object$ref.class,
      constraint = if (type == "LPA") object$arguments$constraint else NULL
    ),
    class.probs = class_probs,
    item.probs = if (type == "LCA") item_probs else NULL,
    class.means = if (type == "LPA") class_means else NULL,
    transition = transition_info,
    covariates = if (!is.null(object$covariates)) {
      lapply(1:length(object$covariates), function(t) {
        cov_names <- if (!is.null(colnames(object$covariates[[t]]))) {
          colnames(object$covariates[[t]])
        } else {
          paste0("Cov", 1:ncol(object$covariates[[t]]))
        }
        cov_summary <- data.frame(
          Variable = cov_names,
          Min = apply(object$covariates[[t]], 2, min),
          Max = apply(object$covariates[[t]], 2, max),
          Mean = apply(object$covariates[[t]], 2, mean)
        )
        rownames(cov_summary) <- NULL
        cov_summary
      })
    } else {
      NULL
    },
    I.max.shown = vars_to_show,
    L.max.shown = classes_to_show,
    total.vars = I,
    total.classes = L,
    digits = digits
  )

  class(summary_obj) <- "summary.sim.LTA"
  return(summary_obj)
}

#' @describeIn summary Summary method for \code{\link[LCPA:get.fit.index]{fit.index}} objects
#' @exportS3Method summary fit.index
summary.fit.index <- function(object, digits = 4, ...) {

  N <- object$N
  npar <- object$npar
  fmt_num <- function(x) sprintf(paste0("%.", digits, "f"), x)
  fit_stats <- data.frame(
    Statistic = c("npar", "Log.Lik", "-2LL", "AIC", "BIC", "SIC", "CAIC", "AWE", "SABIC"),
    Value = c(
      object$npar,
      fmt_num(object$Log.Lik),
      fmt_num(object[["-2LL"]]),
      fmt_num(object$AIC),
      fmt_num(object$BIC),
      fmt_num(object$SIC),
      fmt_num(object$CAIC),
      fmt_num(object$AWE),
      fmt_num(object$SABIC)
    ),
    Description = c(
      "-",
      "Higher is better",
      "Lower is better",
      "Lower is better",
      "Lower is better",
      "-0.5 * BIC",
      "Lower is better",
      "Lower is better",
      "Lower is better"
    ),
    stringsAsFactors = FALSE
  )

  res <- list(
    call = object$call,
    data.info = list(N = N),
    fit.table = fit_stats,
    digits = digits
  )

  class(res) <- "summary.fit.index"
  res
}

#' @describeIn summary Summary method for \code{compare.model} objects
#' @exportS3Method summary compare.model
summary.compare.model <- function(object, digits = 4, ...) {

  fit.index1 <- object$fit.index$model1
  fit.index2 <- object$fit.index$model2
  N <- object$N
  I <- object$I
  L1 <- object$L[1]
  L2 <- object$L[2]
  fmt_num <- function(x) sprintf(paste0("%.", digits, "f"), x)

  fit_stats <- data.frame(
    Statistic = c("class", "npar", "Log.Lik", "-2LL", "AIC", "BIC",
                  "SIC", "CAIC", "AWE", "SABIC"),
    mode1 = c(
      as.character(L1),
      as.character(fit.index1$npar),
      fmt_num(fit.index1$Log.Lik),
      fmt_num(fit.index1[["-2LL"]]),
      fmt_num(fit.index1$AIC),
      fmt_num(fit.index1$BIC),
      fmt_num(fit.index1$SIC),
      fmt_num(fit.index1$CAIC),
      fmt_num(fit.index1$AWE),
      fmt_num(fit.index1$SABIC)
    ),
    mode2 = c(
      as.character(L2),
      as.character(fit.index2$npar),
      fmt_num(fit.index2$Log.Lik),
      fmt_num(fit.index2[["-2LL"]]),
      fmt_num(fit.index2$AIC),
      fmt_num(fit.index2$BIC),
      fmt_num(fit.index2$SIC),
      fmt_num(fit.index2$CAIC),
      fmt_num(fit.index2$AWE),
      fmt_num(fit.index2$SABIC)
    ),
    Description = c("-", "-", "Higher is better", "Lower is better",
                    "Lower is better", "Lower is better",
                    "-0.5 * BIC",
                    "Lower is better", "Lower is better", "Lower is better"),
    stringsAsFactors = FALSE
  )

  AvePP <- object$AvePP
  entropy <- object$entropy

  model.comparison <- data.frame(
    Classes = c(L1, L2),
    npar = c(fit.index1$npar, fit.index2$npar),
    AvePP = fmt_num(c(AvePP$model1[L1+1, L1+1], AvePP$model2[L2+1, L2+1])),
    Entropy = fmt_num(c(entropy[1], entropy[2])),
    stringsAsFactors = FALSE
  )

  BF <- object$BF
  bf_interpretation <- if (!is.na(BF)) {
    if (BF >= 10) {
      "Strong evidence for the larger model (Model 2)"
    } else if (BF >= 5) {
      "Moderate evidence for the larger model (Model 2)"
    } else if (BF >= 3) {
      "Weak evidence for the larger model (Model 2)"
    } else {
      "Strong evidence for the smaller model (Model 1)"
    }
  } else {
    "Bayes factor not available"
  }

  LRT.list <- list()
  test_names <- character(0)
  stats <- numeric(0)
  dfs <- numeric(0)
  pvals <- numeric(0)

  if (!is.null(object$LRT.obj)) {
    LRT.list[["LRT"]] <- object$LRT.obj
    test_names <- c(test_names, "Standard LRT")
    stats <- c(stats, object$LRT.obj$statistic)
    dfs <- c(dfs, object$LRT.obj$parameter)
    pvals <- c(pvals, object$LRT.obj$p.value)
  }

  if (!is.null(object$LRT.VLMR.obj)) {
    LRT.list[["VLMR"]] <- object$LRT.VLMR.obj
    test_names <- c(test_names, "VLMR LRT", "Adjusted LMR LRT")
    stats <- c(
      stats, object$LRT.VLMR.obj$statistic,
      object$LRT.VLMR.obj$adjusted.statistic
    )
    dfs <- c(
      dfs, object$LRT.VLMR.obj$parameter,
      object$LRT.VLMR.obj$parameter
    )
    pvals <- c(
      pvals, object$LRT.VLMR.obj$p.value,
      object$LRT.VLMR.obj$adjusted.p.value
    )
  }

  if (!is.null(object$LRT.Bootstrap.obj)) {
    LRT.list[["Bootstrap"]] <- object$LRT.Bootstrap.obj
    test_names <- c(test_names, "Bootstrap LRT")
    stats <- c(stats, object$LRT.Bootstrap.obj$statistic)
    dfs <- c(dfs, object$LRT.Bootstrap.obj$parameter)
    pvals <- c(pvals, object$LRT.Bootstrap.obj$p.value)
  }

  LRT.table <- if (length(test_names) > 0) {
    sig_labels <- ifelse(pvals < 0.001, "***",
                         ifelse(pvals < 0.01, "**",
                                ifelse(pvals < 0.05, "*", "")))

    data.frame(
      Test = test_names,
      Statistic = fmt_num(stats),
      DF = dfs,
      `p-value` = fmt_num(pvals),
      Sig = sig_labels,
      stringsAsFactors = FALSE
    )
  } else {
    NULL
  }

  res <- list(
    call = object$call,
    data.info = list(N = N, I = I, L = c(L1, L2)),
    fit.table = fit_stats,
    model.comparison = model.comparison,
    BF = BF,
    BF.interpretation = bf_interpretation,
    LRT.table = LRT.table,
    LRT.objects = LRT.list,
    digits = digits
  )

  class(res) <- "summary.compare.model"
  res
}

#' @describeIn summary Summary method for \code{\link[LCPA]{summary.SE}} objects
#' @exportS3Method summary SE
summary.SE <- function(object, ...) {
  # Determine model type
  type <- if (!is.null(object$se$means)) "LPA" else if (!is.null(object$se$par)) "LCA" else "Unknown"

  # Extract dimensions
  L <- length(object$se$P.Z)
  I <- if (type == "LPA") {
    ncol(object$se$means)
  } else if (type == "LCA") {
    dim(object$se$par)[2]
  } else {
    NA
  }

  # Count non-zero SEs
  nonzero.counts <- list(
    P.Z = sum(object$se$P.Z != 0, na.rm = TRUE)
  )
  if (type == "LPA") {
    nonzero.counts$means <- sum(object$se$means != 0, na.rm = TRUE)
    nonzero.counts$covs <- sum(object$se$covs != 0, na.rm = TRUE)
  } else if (type == "LCA") {
    nonzero.counts$par <- sum(object$se$par != 0, na.rm = TRUE)
  }

  res <- list(
    call = object$call,
    method = object$diagnostics$method,
    diagnostics = object$diagnostics,
    type = type,
    L = L,
    I = I,
    nonzero.counts = nonzero.counts,
    total.P.Z = length(object$se$P.Z)
  )
  class(res) <- "summary.SE"
  res
}
