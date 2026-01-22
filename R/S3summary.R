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
#' @param I.max Maximum number of variables/items to display (\code{LCA}, \code{LPA}, \code{sim.LCA}, \code{sim.LPA}, \code{sim.LTA},
#'   \code{LCPA}, \code{LTA}, and \code{compare.model} only; default: 5). Controls verbosity for high-dimensional outputs.
#' @param L.max Maximum number of latent classes/profiles to display before truncation (\code{sim.LTA} only; default: 5).
#'   Useful when models have many latent groups. Ignored for other classes.
#' @param ... Additional arguments passed to or from other methods (currently ignored).
#'
#' @return Invisibly returns a structured list containing summary components.
#'   The exact structure depends on the class of \code{object}. All returned objects carry an appropriate
#'   S3 class (e.g., \code{summary.LCA}, \code{summary.LPA}) for use with corresponding \code{print} methods.
#'
#' @details
#' Each method returns a named list with class-specific components optimized for structured access:
#'
#' \describe{
#'   \item{\strong{\code{LCA}}}{Returns a \code{summary.LCA} object with components:
#'     \describe{
#'       \item{\code{call}}{Original function call.}
#'       \item{\code{model.config}}{List: \code{latent_classes}, \code{method}.}
#'       \item{\code{data.info}}{List: \code{N}, \code{I}, \code{poly.value}, \code{uniform_categories}.}
#'       \item{\code{fit.stats}}{List: \code{LogLik}, \code{AIC}, \code{BIC}, \code{entropy}, \code{npar}.}
#'       \item{\code{class.probs}}{Data frame: \code{Class}, \code{Count}, \code{Proportion}.}
#'       \item{\code{item.probs}}{List of matrices (first \code{I.max} items) with conditional probabilities per class/category.}
#'       \item{\code{convergence}}{List: algorithm, iterations, tolerance, loglik change, hardware (if applicable).}
#'       \item{\code{replication}}{List: \code{nrep}, \code{best_BIC} (if multiple replications performed).}
#'       \item{\code{digits}, \code{I.max.shown}, \code{total.items}}{Metadata for printing/formatting.}
#'     }}
#'
#'   \item{\strong{\code{LPA}}}{Returns a \code{summary.LPA} object with components:
#'     \describe{
#'       \item{\code{call}}{Original function call.}
#'       \item{\code{model.config}}{List: \code{latent_profiles}, \code{constraint}, \code{cov_structure}, \code{method}.}
#'       \item{\code{data.info}}{List: \code{N}, \code{I}, \code{distribution}.}
#'       \item{\code{fit.stats}}{List: \code{LogLik}, \code{AIC}, \code{BIC}, \code{entropy}, \code{npar}.}
#'       \item{\code{class.probs}}{Data frame: \code{Profile}, \code{Count}, \code{Proportion}.}
#'       \item{\code{class.means}}{Matrix (first \code{I.max} variables) of profile-specific means.}
#'       \item{\code{convergence}}{List: algorithm, iterations, tolerance, loglik change, hardware (if applicable).}
#'       \item{\code{replication}}{List: \code{nrep}, \code{best_BIC} (if multiple replications performed).}
#'       \item{\code{digits}, \code{I.max.shown}, \code{total.vars}}{Metadata for printing/formatting.}
#'     }}
#'
#'   \item{\strong{\code{LCPA}}}{Returns a \code{summary.LCPA} object with components:
#'     \describe{
#'       \item{\code{call}}{Original function call.}
#'       \item{\code{model.config}}{List: \code{latent_classes}, \code{model_type}, \code{reference_class},
#'         \code{covariates_mode}, \code{CEP_handling}.}
#'       \item{\code{data.info}}{List: \code{sample_size}, \code{variables}.}
#'       \item{\code{fit.stats}}{List: \code{LogLik}, \code{AIC}, \code{BIC}, \code{npar}.}
#'       \item{\code{class.probs}}{Data frame: \code{Class}, \code{Probability}, \code{Proportion}, \code{Frequency}.}
#'       \item{\code{coefficients}}{Data frame: regression coefficients for non-reference classes (Estimate, Std_Error, z_value, p_value, 95% CI).}
#'       \item{\code{reference_class}}{Integer: reference class for multinomial logit.}
#'       \item{\code{convergence}}{List: \code{iterations}, \code{coveraged}, \code{converg_note}.}
#'       \item{\code{digits}, \code{I.max.shown}, \code{total.vars}, \code{has.covariates}}{Metadata for printing/formatting.}
#'     }}
#'
#'   \item{\strong{\code{LTA}}}{Returns a \code{summary.LTA} object with components:
#'     \describe{
#'       \item{\code{call}}{Original function call.}
#'       \item{\code{model.config}}{List: \code{time_points}, \code{latent_classes}, \code{model_type},
#'         \code{reference_class}, \code{covariates_mode}, \code{CEP_handling}, \code{transition_mode}.}
#'       \item{\code{data.info}}{List: \code{sample_size}, \code{variables}, \code{time_points}.}
#'       \item{\code{fit.stats}}{List: \code{LogLik}, \code{AIC}, \code{BIC}, \code{npar}.}
#'       \item{\code{class.probs}}{List of data frames (per time point): \code{Class}, \code{Probability}, \code{Proportion}, \code{Frequency}.}
#'       \item{\code{initial_model}}{List: \code{coefficients} (data frame), \code{covariate_names}, \code{reference_class}.}
#'       \item{\code{transition_models}}{Named list of data frames: transition coefficients per time interval (From_Class, To_Class, Estimate, Std_Error, etc.).}
#'       \item{\code{reference_class}}{Integer: reference destination class for transitions.}
#'       \item{\code{convergence}}{List: \code{iterations}, \code{coveraged}, \code{converg_note}.}
#'       \item{\code{digits}, \code{I.max.shown}, \code{total.vars}, \code{covariates.timeCross}}{Metadata for printing/formatting.}
#'     }}
#'
#'   \item{\strong{\code{sim.LCA}}}{Returns a \code{summary.sim.LCA} object with components:
#'     \describe{
#'       \item{\code{call}}{Original simulation call.}
#'       \item{\code{config}}{List: \code{N}, \code{I}, \code{L}, \code{poly.value}, \code{uniform_categories}, \code{IQ}, \code{distribution}.}
#'       \item{\code{class.probs}}{Data frame: \code{Class}, \code{Probability}, \code{Frequency}.}
#'       \item{\code{item.probs}}{List of matrices (first \code{I.max} items) with true conditional probabilities per class/category.}
#'       \item{\code{digits}, \code{I.max.shown}, \code{total.vars}}{Metadata for printing/formatting.}
#'     }}
#'
#'   \item{\strong{\code{sim.LPA}}}{Returns a \code{summary.sim.LPA} object with components:
#'     \describe{
#'       \item{\code{call}}{Original simulation call.}
#'       \item{\code{config}}{List: \code{N}, \code{I}, \code{L}, \code{constraint}, \code{constraint_desc}, \code{distribution}.}
#'       \item{\code{class.probs}}{Data frame: \code{Profile}, \code{Probability}, \code{Frequency}.}
#'       \item{\code{class.means}}{Matrix (first \code{I.max} variables) of true profile-specific means.}
#'       \item{\code{cov_structure}}{Character: detailed description of covariance constraints.}
#'       \item{\code{digits}, \code{I.max.shown}, \code{total.vars}}{Metadata for printing/formatting.}
#'     }}
#'
#'   \item{\strong{\code{sim.LTA}}}{Returns a \code{summary.sim.LTA} object with components:
#'     \describe{
#'       \item{\code{call}}{Original simulation call.}
#'       \item{\code{config}}{List: \code{N}, \code{I}, \code{L}, \code{times}, \code{type}, \code{distribution}, \code{constraint} (if LPA).}
#'       \item{\code{class.probs}}{List of data frames (per time point): \code{Class}, \code{Probability}, \code{Frequency}.}
#'       \item{\code{item.probs}}{Nested list (by time/item) of true conditional probabilities (if \code{type="LCA"}).}
#'       \item{\code{class.means}}{List of matrices (by time) of true profile means (if \code{type="LPA"}).}
#'       \item{\code{transition}}{List: \code{mode} ("fixed" or "covariate"), \code{rate} or \code{beta}/\code{gamma} coefficients, \code{time_points}.}
#'       \item{\code{covariates}}{List of data frames (per time point) with covariate summaries (Min, Max, Mean), if present.}
#'       \item{\code{digits}, \code{I.max.shown}, \code{L.max.shown}, \code{total.vars}, \code{total.classes}}{Metadata for printing/formatting.}
#'     }}
#'
#'   \item{\strong{\code{fit.index}}}{Returns a \code{summary.fit.index} object with components:
#'     \describe{
#'       \item{\code{call}}{Function call that generated the fit indices.}
#'       \item{\code{data.info}}{List: \code{N}.}
#'       \item{\code{fit.table}}{Data frame: \code{Statistic}, \code{Value}, \code{Description} for -2LL, AIC, BIC, SIC, CAIC, AWE, SABIC.}
#'       \item{\code{digits}}{Numeric: precision used for formatting.}
#'     }}
#'
#'   \item{\strong{\code{compare.model}}}{Returns a \code{summary.compare.model} object with components:
#'     \describe{
#'       \item{\code{call}}{Function call that generated the comparison.}
#'       \item{\code{data.info}}{List: \code{N}, \code{I}, \code{L} (named vector for two models).}
#'       \item{\code{fit.table}}{Data frame comparing fit indices for both models.}
#'       \item{\code{model_comparison}}{Data frame: \code{Classes}, \code{npar}, \code{AvePP}, \code{Entropy}.}
#'       \item{\code{BF}}{Numeric: Bayes Factor value (if computed).}
#'       \item{\code{BF_interpretation}}{Character: interpretive guidance for Bayes Factor.}
#'       \item{\code{lrt_table}}{Data frame: \code{Test}, \code{Statistic}, \code{DF}, \code{p-value}, \code{Sig} (significance markers).}
#'       \item{\code{lrt_objects}}{List: raw hypothesis test objects for further inspection.}
#'       \item{\code{digits}}{Numeric: precision used for formatting.}
#'     }}
#'
#'   \item{\strong{\code{SE}}}{Returns a \code{summary.SE} object with components:
#'     \describe{
#'       \item{\code{call}}{Original function call.}
#'       \item{\code{method}}{Character: "Obs" or "Bootstrap".}
#'       \item{\code{diagnostics}}{List: method-specific diagnostic info (e.g., n.Bootstrap, hessian_cond_number).}
#'       \item{\code{model_type}}{Character: "LCA" or "LPA".}
#'       \item{\code{L}}{Integer: number of latent classes/profiles.}
#'       \item{\code{I}}{Integer: number of variables/items (NA if unknown).}
#'       \item{\code{nonzero_counts}}{List: counts of non-zero SEs by parameter type (P.Z, means/par, covs).}
#'       \item{\code{total_PZ}}{Integer: total number of class probability parameters.}
#'     }}
#' }
#'
#' @name summary
NULL

#' @describeIn summary Summary method for \code{LCA} objects
#' @importFrom utils tail
#' @export
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
    Count = as.numeric(table(object$Z)),
    Proportion = sprintf("%.1f%%", P.Z * 100),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  items_to_show <- min(I.max, I)
  item.probs <- if (I > 0) {
    lapply(seq_len(items_to_show), function(i) {
      prob_mat <- round(object$probability[[i]], digits)
      dimnames(prob_mat) <- list(
        paste("Class", 1:nrow(prob_mat)),
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
      tolerance = tol,
      loglik_change = ll.delta,
      loglik_initial = ll.initial,
      loglik_final = ll.final
    )
  } else if (arguments$method == "NNE") {
    iter.n <- length(object$Log.Lik.history) - 1
    maxiter_early <- if (!is.null(arguments$control.NNE) && !is.null(arguments$control.NNE$maxiter.early)) {
      arguments$control.NNE$maxiter.early
    } else {
      5
    }
    ll.initial <- object$Log.Lik.history[1]
    ll.final <- tail(object$Log.Lik.history, 1)
    ll.delta <- abs(ll.final - ll.initial)

    convergence <- list(
      algorithm = "Neural Network Estimation (NNE)",
      iterations = iter.n,
      early_stop_threshold = maxiter_early,
      loglik_change = ll.delta,
      loglik_initial = ll.initial,
      loglik_final = ll.final,
      hardware = arguments$control.NNE$device
    )
  } else if (arguments$method == "Mplus") {
    convergence <- list(
      algorithm = "Mplus (External Estimation)",
      note = "Convergence diagnostics unavailable for external estimators"
    )
  }

  replication <- if (arguments$nrep > 1 && arguments$method != "Mplus") {
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
      latent_classes = arguments$L,
      method = arguments$method
    ),
    data.info = list(
      N = N,
      I = I,
      poly.value = poly.value,
      uniform_categories = poly.value.uniform
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
  invisible(summary_obj)
}

#' @describeIn summary Summary method for \code{LPA} objects
#' @importFrom utils tail
#' @export
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
    Count = as.numeric(table(object$Z)),
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
        tolerance = tol,
        loglik_change = ll.delta,
        loglik_initial = ll.initial,
        loglik_final = ll.final
      )
    } else if (arguments$method == "NNE") {
      maxiter_early <- if (!is.null(arguments$control.NNE) && !is.null(arguments$control.NNE$maxiter.early)) {
        arguments$control.NNE$maxiter.early
      } else {
        5
      }

      convergence <- list(
        algorithm = "Neural Network Estimation (NNE)",
        iterations = iter.n,
        early_stop_threshold = maxiter_early,
        loglik_change = ll.delta,
        loglik_initial = ll.initial,
        loglik_final = ll.final,
        hardware = arguments$control.NNE$device
      )
    }
  } else if (arguments$method == "Mplus") {
    convergence <- list(
      algorithm = "Mplus (External Estimation)",
      note = "Convergence diagnostics unavailable for external estimators"
    )
  }

  # Replication info
  replication <- if (arguments$nrep > 1 && arguments$method != "Mplus") {
    list(
      nrep = arguments$nrep,
      best_BIC = min(object$Log.Lik.nrep) * -2 + log(N) * object$npar  # Recalculate BIC from stored log-likelihoods
    )
  } else {
    NULL
  }

  # Covariance structure description
  if(any(arguments$constraint %in% c("VV", "VE", "EV", "EE", "V0", "E0", "UV", "UE"))){
    cov_structure <- switch(arguments$constraint,
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
    cov_structure <- paste(all_constraints, collapse = "\n")
  }

  # Build summary object
  summary_obj <- list(
    call = call_info,
    model.config = list(
      latent_profiles = L,
      constraint = arguments$constraint,
      cov_structure = cov_structure,
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
  invisible(summary_obj)
}

#' @describeIn summary Summary method for \code{LTA} objects
#' @export
summary.LTA <- function(object, digits = 4, I.max = 5, ...) {
  arguments <- object$arguments
  times <- length(object$P.Zs)
  L <- ncol(object$beta)
  ref.class <- arguments$ref.class
  type <- arguments$type
  N <- nrow(arguments$responses[[1]])
  I <- ncol(arguments$responses[[1]])
  covariates.timeCross <- arguments$covariates.timeCross

  # Model configuration
  model_config <- list(
    time_points = times,
    latent_classes = L,
    model_type = ifelse(type == "LCA", "Latent Class Analysis (categorical)",
                        "Latent Profile Analysis (continuous)"),
    reference_class = ref.class,
    covariates_mode = if (!is.null(arguments$covariates)) {
      if (covariates.timeCross) "Time-invariant covariates" else "Time-varying covariates"
    } else "No covariates",
    CEP_handling = ifelse(arguments$CEP.error,
                          ifelse(arguments$CEP.timeCross,
                                 "Classification error correction (time-invariant CEP)",
                                 "Classification error correction (time-varying CEP)"),
                          "No classification error correction (naive modal assignment)"),
    transition_mode = ifelse(covariates.timeCross,
                             "Time-invariant transition effects (coefficients constant across time)",
                             "Time-varying transition effects (coefficients differ by time point)")
  )

  # Class probabilities over time
  class_probs <- lapply(1:times, function(t) {
    data.frame(
      Class = as.character(1:L),
      Probability = round(object$P.Zs[[t]], digits),
      Proportion = sprintf("%.1f%%", object$P.Zs[[t]] * 100),
      Frequency = as.vector(table(factor(object$Zs[[t]], levels = 1:L)))
    )
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

  if (covariates.timeCross) {
    # Use first time point's covariate names for all
    cov_names_list <- replicate(times, cov_names_list[[1]], simplify = FALSE)
  }

  # Initial state coefficients table (always needed)
  cov_names_initial <- cov_names_list[[1]]
  vars_to_show_initial <- min(I.max, length(cov_names_initial))

  initial_coef <- data.frame(
    Class = character(),
    Covariate = character(),
    Estimate = numeric(),
    Std_Error = numeric(),
    lower_95 = numeric(),
    upper_95 = numeric(),
    z_value = numeric(),
    p_value = numeric(),
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
      lower_95 <- if (!is.na(se_val) && is.finite(se_val)) coef_val - 1.96 * se_val else NA
      upper_95 <- if (!is.na(se_val) && is.finite(se_val)) coef_val + 1.96 * se_val else NA

      initial_coef <- rbind(initial_coef, data.frame(
        Class = paste0("Class ", cls),
        Covariate = cov_names_initial[j],
        Estimate = round(coef_val, digits),
        Std_Error = ifelse(is.na(se_val) || !is.finite(se_val), NA, round(se_val, digits)),
        lower_95 = ifelse(is.na(se_val) || !is.finite(se_val), NA, round(lower_95, digits)),
        upper_95 = ifelse(is.na(se_val) || !is.finite(se_val), NA, round(upper_95, digits)),
        z_value = ifelse(is.na(z_val) || !is.finite(z_val), NA, round(z_val, digits)),
        p_value = ifelse(is.na(p_val) || !is.finite(p_val), NA, sprintf("%.4f", p_val))
      ))
    }
  }

  # Transition coefficients: time-invariant vs time-varying
  if (covariates.timeCross) {
    # Only one transition model needed (all time points share same coefficients)
    trans_df <- data.frame(
      From_Class = character(),
      To_Class = character(),
      Covariate = character(),
      Estimate = numeric(),
      Std_Error = numeric(),
      lower_95 = numeric(),
      upper_95 = numeric(),
      z_value = numeric(),
      p_value = numeric(),
      stringsAsFactors = FALSE
    )

    cov_names_trans <- cov_names_list[[2]]  # Same for all time points
    vars_to_show_trans <- min(I.max, length(cov_names_trans))

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
          lower_95 <- if (!is.na(se_val) && is.finite(se_val)) coef_val - 1.96 * se_val else NA
          upper_95 <- if (!is.na(se_val) && is.finite(se_val)) coef_val + 1.96 * se_val else NA

          trans_df <- rbind(trans_df, data.frame(
            From_Class = paste0("Class ", from_cls),
            To_Class = paste0("Class ", to_cls),
            Covariate = cov_names_trans[j],
            Estimate = round(coef_val, digits),
            Std_Error = ifelse(is.na(se_val) || !is.finite(se_val), NA, round(se_val, digits)),
            lower_95 = ifelse(is.na(se_val) || !is.finite(se_val), NA, round(lower_95, digits)),
            upper_95 = ifelse(is.na(se_val) || !is.finite(se_val), NA, round(upper_95, digits)),
            z_value = ifelse(is.na(z_val) || !is.finite(z_val), NA, round(z_val, digits)),
            p_value = ifelse(is.na(p_val) || !is.finite(p_val), NA, sprintf("%.4f", p_val))
          ))
        }
      }
    }

    # Add reference class explanation to dataframe
    attr(trans_df, "ref_class_note") <- paste0("All transitions are relative to reference destination Class ", ref.class)

    # Single transition model with special label
    transition_coefs <- list(`Time-invariant effects` = trans_df)
  } else {
    # Time-varying effects: show each transition phase
    transition_coefs <- list()
    for (t in 1:(times-1)) {
      cov_names_trans <- cov_names_list[[t+1]]
      vars_to_show_trans <- min(I.max, length(cov_names_trans))

      trans_df <- data.frame(
        From_Class = character(),
        To_Class = character(),
        Covariate = character(),
        Estimate = numeric(),
        Std_Error = numeric(),
        lower_95 = numeric(),
        upper_95 = numeric(),
        z_value = numeric(),
        p_value = numeric(),
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
            lower_95 <- if (!is.na(se_val) && is.finite(se_val)) coef_val - 1.96 * se_val else NA
            upper_95 <- if (!is.na(se_val) && is.finite(se_val)) coef_val + 1.96 * se_val else NA

            trans_df <- rbind(trans_df, data.frame(
              From_Class = paste0("Class ", from_cls),
              To_Class = paste0("Class ", to_cls),
              Covariate = cov_names_trans[j],
              Estimate = round(coef_val, digits),
              Std_Error = ifelse(is.na(se_val) || !is.finite(se_val), NA, round(se_val, digits)),
              lower_95 = ifelse(is.na(se_val) || !is.finite(se_val), NA, round(lower_95, digits)),
              upper_95 = ifelse(is.na(se_val) || !is.finite(se_val), NA, round(upper_95, digits)),
              z_value = ifelse(is.na(z_val) || !is.finite(z_val), NA, round(z_val, digits)),
              p_value = ifelse(is.na(p_val) || !is.finite(p_val), NA, sprintf("%.4f", p_val))
            ))
          }
        }
      }
      # Add reference class explanation
      attr(trans_df, "ref_class_note") <- paste0("All transitions are relative to reference destination Class ", ref.class)
      transition_coefs[[t]] <- trans_df
    }
    names(transition_coefs) <- paste0("Time ", 1:(times-1), " -> Time ", 2:times)
  }

  # Convergence information
  log_lik_at_start <- if (length(object$Log.Lik.history) > 0) object$Log.Lik.history[2] else NA
  log_lik_at_end <- object$Log.Lik
  loglik_change <- if (!is.na(log_lik_at_start)) round(log_lik_at_end - log_lik_at_start, digits) else NA

  convergence_info <- list(
    iterations = object$iterations,
    coveraged = object$coveraged,
    converg_note = sprintf("  Log-likelihood change: |%.2f - %.2f| = %.4f\n",
                           log_lik_at_start, log_lik_at_end, loglik_change)
  )

  # Build summary object
  summary_obj <- list(
    call = object$call,
    model.config = model_config,
    data.info = list(
      sample_size = N,
      variables = I,
      time_points = times
    ),
    fit.stats = list(
      LogLik = object$Log.Lik,
      AIC = object$AIC,
      BIC = object$BIC,
      npar = object$npar
    ),
    class.probs = class_probs,
    initial_model = list(
      coefficients = initial_coef,
      covariate_names = cov_names_initial,
      reference_class = ref.class
    ),
    transition_models = transition_coefs,
    convergence = convergence_info,
    digits = digits,
    I.max.shown = I.max,
    total.vars = if (!is.null(arguments$covariates)) max(sapply(arguments$covariates, ncol)) else 1,
    covariates.timeCross = covariates.timeCross,
    reference_class = ref.class  # Store reference class at top level for easy access
  )

  class(summary_obj) <- "summary.LTA"
  invisible(summary_obj)
}

#' @describeIn summary Summary method for \code{LCPA} objects
#' @export
summary.LCPA <- function(object, digits = 4, I.max = 5, ...) {
  arguments <- object$arguments
  L <- ncol(object$beta)
  ref.class <- arguments$ref.class
  type <- arguments$type
  N <- nrow(object$P.Z.Xn)  # Single time point
  I <- ncol(arguments$response)

  # Model configuration
  model_config <- list(
    latent_classes = L,
    model_type = ifelse(type == "LCA", "Latent Class Analysis (categorical)",
                        "Latent Profile Analysis (continuous)"),
    reference_class = ref.class,
    covariates_mode = if (!is.null(arguments$covariate)) {
      "Covariates included in class membership model"
    } else "No covariates (intercept-only model)",
    CEP_handling = ifelse(arguments$CEP.error,
                          "Classification error correction applied",
                          "No classification error correction (naive modal assignment)")
  )

  # Class probabilities (single time point)
  class_probs_df <- data.frame(
    Class = as.character(1:L),
    Probability = round(object$P.Z, digits),
    Proportion = sprintf("%.1f%%", object$P.Z * 100),
    Frequency = as.vector(table(factor(object$Z, levels = 1:L)))
  )

  # Extract covariate names
  cov_mat <- arguments$covariate
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

  vars_to_show <- min(I.max, length(cov_names))

  # Coefficients table (beta)
  coef_df <- data.frame(
    Class = character(),
    Covariate = character(),
    Estimate = numeric(),
    Std_Error = numeric(),
    lower_95 = numeric(),
    upper_95 = numeric(),
    z_value = numeric(),
    p_value = numeric(),
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
      lower_95 <- if (!is.na(se_val) && is.finite(se_val)) coef_val - 1.96 * se_val else NA
      upper_95 <- if (!is.na(se_val) && is.finite(se_val)) coef_val + 1.96 * se_val else NA

      coef_df <- rbind(coef_df, data.frame(
        Class = paste0("Class ", cls),
        Covariate = cov_names[j],
        Estimate = round(coef_val, digits),
        Std_Error = ifelse(is.na(se_val) || !is.finite(se_val), NA, round(se_val, digits)),
        lower_95 = ifelse(is.na(se_val) || !is.finite(se_val), NA, round(lower_95, digits)),
        upper_95 = ifelse(is.na(se_val) || !is.finite(se_val), NA, round(upper_95, digits)),
        z_value = ifelse(is.na(z_val) || !is.finite(z_val), NA, round(z_val, digits)),
        p_value = ifelse(is.na(p_val) || !is.finite(p_val), NA, sprintf("%.4f", p_val))
      ))
    }
  }

  # Convergence information
  log_lik_at_start <- if (length(object$Log.Lik.history) > 0) object$Log.Lik.history[2] else NA
  log_lik_at_end <- object$Log.Lik
  loglik_change <- if (!is.na(log_lik_at_start)) round(log_lik_at_end - log_lik_at_start, digits) else NA

  convergence_info <- list(
    iterations = object$iterations,
    coveraged = object$coveraged,
    converg_note = if (!is.na(loglik_change)) {
      sprintf("  Log-likelihood change: |%.2f - %.2f| = %.4f\n",
              log_lik_at_start, log_lik_at_end, loglik_change)
    } else {
      "  Log-likelihood history unavailable\n"
    }
  )

  # Build summary object
  summary_obj <- list(
    call = object$call,
    model.config = model_config,
    data.info = list(
      sample_size = N,
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
    covariate_names = cov_names,
    reference_class = ref.class,
    convergence = convergence_info,
    digits = digits,
    I.max.shown = vars_to_show,
    total.vars = length(cov_names),
    has.covariates = !is.null(arguments$covariate)
  )

  class(summary_obj) <- "summary.LCPA"
  invisible(summary_obj)
}

#' @describeIn summary Summary method for \code{sim.LCA} objects
#' @export
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
      paste0("Class ", 1:L),
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
      uniform_categories = length(unique(object$poly.value)) == 1,
      IQ = object$arguments$IQ,
      distribution = object$arguments$distribution
    ),
    class.probs = data.frame(
      Class = paste0("L", 1:L),
      Probability = round(object$P.Z, digits),
      Frequency = as.vector(table(object$Z))
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
#' @export
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
      Profile = names(object$P.Z),
      Probability = object$P.Z,
      Frequency = as.vector(table(object$Z))
    ),
    class.means = class.means,
    cov_structure = constraint_desc$details,
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
#' @export
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
    data.frame(
      Class = paste0("L", 1:L),
      Probability = round(object$P.Zs[[t]], digits),
      Frequency = as.vector(table(object$Zs[[t]]))
    )
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
          paste0("Class ", 1:classes_to_show),
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
        paste0("From L", 1:L),
        paste0("To L", 1:L)
      )
      rate_mat
    })
    # Store time point indices directly for printing
    transition_info$time_points <- data.frame(
      from = 1:(times-1),
      to = 2:times
    )
  } else if (!is.null(object$beta) && !is.null(object$gamma)) {
    # Covariate mode: regression coefficients
    transition_info$mode <- "covariate"
    transition_info$beta <- round(object$beta, digits)
    dimnames(transition_info$beta) <- list(
      colnames(object$covariates[[1]]),
      paste0("Class ", 1:ncol(object$beta))
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
          colnames(object$covariates[[t]]),
          paste0("To L", 1:L)
        )
      }
      names(gamma_t) <- paste0("From L", 1:L)
      gamma_display[[t]] <- gamma_t
    }
    names(gamma_display) <- paste0("T", 1:(times-1), "->T", 2:times)
    transition_info$gamma <- gamma_display

    # Store time point indices directly for printing
    transition_info$time_points <- data.frame(
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

#' @describeIn summary Summary method for \code{fit.index} objects
#' @export
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
  invisible(res)
}

#' @describeIn summary Summary method for \code{compare.model} objects
#' @export
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
      fmt_num(fit.index1$npar),
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
      fmt_num(fit.index2$npar),
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

  model_comparison <- data.frame(
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

  lrt_list <- list()
  test_names <- character(0)
  stats <- numeric(0)
  dfs <- numeric(0)
  pvals <- numeric(0)

  if (!is.null(object$LRT.obj)) {
    lrt_list[["LRT"]] <- object$LRT.obj
    test_names <- c(test_names, "Standard LRT")
    stats <- c(stats, object$LRT.obj$statistic)
    dfs <- c(dfs, object$LRT.obj$parameter)
    pvals <- c(pvals, object$LRT.obj$p.value)
  }

  if (!is.null(object$LRT.VLMR.obj)) {
    lrt_list[["VLMR"]] <- object$LRT.VLMR.obj
    test_names <- c(test_names, "VLMR-adjusted LRT")
    stats <- c(stats, object$LRT.VLMR.obj$statistic)
    dfs <- c(dfs, object$LRT.VLMR.obj$parameter)
    pvals <- c(pvals, object$LRT.VLMR.obj$p.value)
  }

  if (!is.null(object$LRT.Bootstrap.obj)) {
    lrt_list[["Bootstrap"]] <- object$LRT.Bootstrap.obj
    test_names <- c(test_names, "Bootstrap LRT")
    stats <- c(stats, object$LRT.Bootstrap.obj$statistic)
    dfs <- c(dfs, object$LRT.Bootstrap.obj$parameter)
    pvals <- c(pvals, object$LRT.Bootstrap.obj$p.value)
  }

  lrt_table <- if (length(test_names) > 0) {
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
    model_comparison = model_comparison,
    BF = BF,
    BF_interpretation = bf_interpretation,
    lrt_table = lrt_table,
    lrt_objects = lrt_list,
    digits = digits
  )

  class(res) <- "summary.compare.model"
  invisible(res)
}

#' @describeIn summary Summary method for \code{summary.SE} objects
#' @export
summary.SE <- function(object, ...) {
  # Determine model type
  model_type <- if (!is.null(object$se$means)) "LPA" else if (!is.null(object$se$par)) "LCA" else "Unknown"

  # Extract dimensions
  L <- length(object$se$P.Z)
  I <- if (model_type == "LPA") {
    ncol(object$se$means)
  } else if (model_type == "LCA") {
    dim(object$se$par)[2]
  } else {
    NA
  }

  # Count non-zero SEs
  nonzero_counts <- list(
    P.Z = sum(object$se$P.Z != 0, na.rm = TRUE)
  )
  if (model_type == "LPA") {
    nonzero_counts$means <- sum(object$se$means != 0, na.rm = TRUE)
    nonzero_counts$covs <- sum(object$se$covs != 0, na.rm = TRUE)
  } else if (model_type == "LCA") {
    nonzero_counts$par <- sum(object$se$par != 0, na.rm = TRUE)
  }

  res <- list(
    call = object$call,
    method = object$diagnostics$method,
    diagnostics = object$diagnostics,
    model_type = model_type,
    L = L,
    I = I,
    nonzero_counts = nonzero_counts,
    total_PZ = length(object$se$P.Z)
  )
  class(res) <- "summary.SE"
  res
}
