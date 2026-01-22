#' @title  S3 Methods: extract
#'
#' @description
#' A generic S3 extractor function designed to retrieve internal components from various model and simulation objects
#' produced by the \code{LCPA} package. This function provides a consistent interface across different classes,
#' allowing users to access estimated parameters, fit statistics, simulation truths, standard errors, and more.
#'
#' @param object An object of one of the following classes:
#'   \itemize{
#'     \item \code{\link[LCPA]{LCA}} — Latent Class Analysis model results.
#'     \item \code{\link[LCPA]{LPA}} — Latent Profile Analysis model results.
#'     \item \code{\link[LCPA]{LCPA}} — Latent Class/Profile Analysis with covariates.
#'     \item \code{\link[LCPA]{LTA}} — Latent Transition Analysis model results.
#'     \item \code{\link[LCPA]{sim.LCA}} — Simulated LCA data with known truth.
#'     \item \code{\link[LCPA]{sim.LPA}} — Simulated LPA data with known truth.
#'     \item \code{\link[LCPA]{sim.LTA}} — Simulated LTA data with known truth.
#'     \item \code{\link[LCPA]{get.fit.index}} — Model fit indices object.
#'     \item \code{\link[LCPA]{compare.model}} — Model comparison results.
#'     \item \code{\link[LCPA]{get.SE}} — Standard error estimation results.
#'   }
#' @param what A character string specifying the name of the component to extract.
#'   Valid choices depend on the class of \code{object}. See Details section for full listings.
#' @param ... Additional arguments passed to methods (currently ignored).
#'
#' @return The requested component. Return type varies depending on \code{what} and the class of \code{object}.
#'   If an invalid \code{what} is provided, an informative error is thrown listing valid options.
#'
#' @details
#' This function supports extraction from ten primary object classes. Below are available components for each:
#'
#' \describe{
#'   \item{\strong{\code{LCA}}}{Latent Class Analysis model results. Available components:
#'     \describe{
#'       \item{\code{params}}{List containing all estimated model parameters.}
#'       \item{\code{par}}{3D array (\eqn{L \times I \times K_{\max}}) of conditional response probabilities.}
#'       \item{\code{P.Z}}{Vector of length \eqn{L} with latent class prior probabilities.}
#'       \item{\code{npar}}{Number of free parameters in the model.}
#'       \item{\code{Log.Lik}}{Log-likelihood of the final model.}
#'       \item{\code{AIC}}{Akaike Information Criterion.}
#'       \item{\code{BIC}}{Bayesian Information Criterion.}
#'       \item{\code{best_BIC}}{Best BIC value across replication runs (if \code{nrep > 1}).}
#'       \item{\code{P.Z.Xn}}{\eqn{N \times L} matrix of posterior class probabilities.}
#'       \item{\code{Z}}{Vector of length \eqn{N} with MAP-classified latent class memberships.}
#'       \item{\code{probability}}{List of formatted conditional probability matrices per item.}
#'       \item{\code{Log.Lik.history}}{Vector tracking log-likelihood at each EM iteration.}
#'       \item{\code{Log.Lik.nrep}}{Vector of log-likelihoods from each replication run.}
#'       \item{\code{model}}{Trained neural network model object (only when \code{method="NNE"}).}
#'       \item{\code{call}}{The original function call used for model estimation.}
#'       \item{\code{arguments}}{List containing all input arguments passed to the \code{LCA} function.}
#'     }}
#'
#'   \item{\strong{\code{LPA}}}{Latent Profile Analysis model results. Available components:
#'     \describe{
#'       \item{\code{params}}{List containing all estimated model parameters.}
#'       \item{\code{means}}{\eqn{L \times I} matrix of estimated mean vectors for each profile.}
#'       \item{\code{covs}}{\eqn{I \times I \times L} array of estimated covariance matrices.}
#'       \item{\code{P.Z}}{Vector of length \eqn{L} with profile prior probabilities.}
#'       \item{\code{npar}}{Number of free parameters (depends on \code{constraint}).}
#'       \item{\code{Log.Lik}}{Log-likelihood of the final model.}
#'       \item{\code{AIC}}{Akaike Information Criterion.}
#'       \item{\code{BIC}}{Bayesian Information Criterion.}
#'       \item{\code{best_BIC}}{Best BIC value across replication runs (if \code{nrep > 1}).}
#'       \item{\code{P.Z.Xn}}{\eqn{N \times L} matrix of posterior profile probabilities.}
#'       \item{\code{Z}}{Vector of length \eqn{N} with MAP-classified profile memberships.}
#'       \item{\code{Log.Lik.history}}{Vector tracking log-likelihood at each EM iteration.}
#'       \item{\code{Log.Lik.nrep}}{Vector of log-likelihoods from each replication run.}
#'       \item{\code{model}}{Trained model object (neural network or Mplus).}
#'       \item{\code{call}}{The original function call used for model estimation.}
#'       \item{\code{arguments}}{List containing all input arguments passed to the \code{LPA} function.}
#'       \item{\code{constraint}}{Covariance structure constraints applied during estimation (from original arguments).}
#'     }}
#'
#'   \item{\strong{\code{LCPA}}}{Latent Class/Profile Analysis (with covariates). Available components:
#'     \describe{
#'       \item{\code{beta}}{Initial class coefficients (p1 x L matrix).}
#'       \item{\code{beta.se}}{Standard errors for beta.}
#'       \item{\code{beta.Z.sta}}{Z-statistics for beta.}
#'       \item{\code{beta.p.value.tail1}}{One-tailed p-values for beta.}
#'       \item{\code{beta.p.value.tail2}}{Two-tailed p-values for beta.}
#'       \item{\code{P.Z.Xn}}{Posterior probabilities (N x L).}
#'       \item{\code{P.Z}}{Prior proportions (length L).}
#'       \item{\code{Z}}{Modal class assignments (length N).}
#'       \item{\code{npar}}{Number of free parameters.}
#'       \item{\code{Log.Lik}}{Log-likelihood.}
#'       \item{\code{AIC}}{AIC.}
#'       \item{\code{BIC}}{BIC.}
#'       \item{\code{iterations}}{Optimization iterations in Step 3.}
#'       \item{\code{coveraged}}{Logical: did optimization converge early?}
#'       \item{\code{params}}{Step 1 model parameters (LCA/LPA output).}
#'       \item{\code{call}}{Function call.}
#'       \item{\code{arguments}}{Input arguments list.}
#'     }}
#'
#'   \item{\strong{\code{LTA}}}{Latent Transition Analysis model results. Available components:
#'     \describe{
#'       \item{\code{beta}}{Initial class coefficients (p1 x L matrix).}
#'       \item{\code{gamma}}{Transition coefficients (nested list).}
#'       \item{\code{beta.se}}{Standard errors for beta.}
#'       \item{\code{gamma.se}}{Standard errors for gamma.}
#'       \item{\code{beta.Z.sta}}{Z-statistics for beta.}
#'       \item{\code{gamma.Z.sta}}{Z-statistics for gamma.}
#'       \item{\code{beta.p.value.tail1}}{One-tailed p-values for beta.}
#'       \item{\code{gamma.p.value.tail1}}{One-tailed p-values for gamma.}
#'       \item{\code{beta.p.value.tail2}}{Two-tailed p-values for beta.}
#'       \item{\code{gamma.p.value.tail2}}{Two-tailed p-values for gamma.}
#'       \item{\code{P.Z.Xns}}{List of posterior probabilities per time (each N x L).}
#'       \item{\code{P.Zs}}{List of prior proportions per time (each length L).}
#'       \item{\code{Zs}}{List of modal class assignments per time (each length N).}
#'       \item{\code{npar}}{Number of free parameters.}
#'       \item{\code{Log.Lik}}{Log-likelihood.}
#'       \item{\code{AIC}}{AIC.}
#'       \item{\code{BIC}}{BIC.}
#'       \item{\code{iterations}}{Optimization iterations in Step 3.}
#'       \item{\code{coveraged}}{Logical: did optimization converge early?}
#'       \item{\code{params}}{Step 1 model parameters (LCA/LPA output).}
#'       \item{\code{call}}{Function call.}
#'       \item{\code{arguments}}{Input arguments list.}
#'     }}
#'
#'   \item{\strong{\code{sim.LCA}}}{Simulated Latent Class Analysis data. Available components:
#'     \describe{
#'       \item{\code{response}}{Integer matrix (\eqn{N \times I}) of simulated categorical observations.}
#'       \item{\code{par}}{Array (\eqn{L \times I \times P_{\max}}) of true class-specific category probabilities.}
#'       \item{\code{Z}}{Integer vector (length \eqn{N}) of true latent class assignments.}
#'       \item{\code{P.Z}}{Numeric vector (length \eqn{L}) of true class proportions.}
#'       \item{\code{poly.value}}{Integer vector (length \eqn{I}) specifying categories per variable.}
#'       \item{\code{P.Z.Xn}}{Binary matrix (\eqn{N \times L}) of true class membership indicators.}
#'       \item{\code{call}}{The original function call used for simulation.}
#'       \item{\code{arguments}}{List containing all input arguments passed to \code{\link[LCPA]{sim.LCA}}.}
#'     }}
#'
#'   \item{\strong{\code{sim.LPA}}}{Simulated Latent Profile Analysis data. Available components:
#'     \describe{
#'       \item{\code{response}}{Numeric matrix (\eqn{N \times I}) of simulated continuous observations.}
#'       \item{\code{means}}{\eqn{L \times I} matrix of true class-specific means.}
#'       \item{\code{covs}}{\eqn{I \times I \times L} array of true class-specific covariance matrices.}
#'       \item{\code{P.Z.Xn}}{\eqn{N \times L} matrix of true class membership indicators.}
#'       \item{\code{P.Z}}{Numeric vector (length \eqn{L}) of true class proportions.}
#'       \item{\code{Z}}{Integer vector (length \eqn{N}) of true profile assignments.}
#'       \item{\code{constraint}}{Original constraint specification passed to \code{\link[LCPA]{sim.LPA}}.}
#'       \item{\code{call}}{The original function call used for simulation.}
#'       \item{\code{arguments}}{List containing all input arguments passed to \code{\link[LCPA]{sim.LPA}}.}
#'     }}
#'
#'   \item{\strong{\code{sim.LTA}}}{Simulated Latent Transition Analysis data. Available components:
#'     \describe{
#'       \item{\code{responses}}{List of response matrices per time point.}
#'       \item{\code{Zs}}{List of true latent class assignments per time.}
#'       \item{\code{P.Zs}}{List of true class proportions per time.}
#'       \item{\code{par}}{True conditional probabilities (for categorical items).}
#'       \item{\code{means}}{True profile means (for continuous variables).}
#'       \item{\code{covs}}{True covariance matrices per class and time.}
#'       \item{\code{poly.value}}{Categories per variable (for categorical items).}
#'       \item{\code{rate}}{Transition rate matrix or structure.}
#'       \item{\code{covariates}}{Simulated covariate matrix.}
#'       \item{\code{beta}}{True initial class coefficients.}
#'       \item{\code{gamma}}{True transition coefficients.}
#'       \item{\code{call}}{Original simulation function call.}
#'       \item{\code{arguments}}{Input arguments used in simulation.}
#'     }}
#'
#'   \item{\strong{\code{fit.index}}}{Model fit indices object. Available components:
#'     \describe{
#'       \item{\code{npar}}{Number of free parameters in the model.}
#'       \item{\code{Log.Lik}}{Log-likelihood of the model.}
#'       \item{\code{-2LL}}{Deviance statistic (-2 times log-likelihood).}
#'       \item{\code{AIC}}{Akaike Information Criterion.}
#'       \item{\code{BIC}}{Bayesian Information Criterion.}
#'       \item{\code{SIC}}{Sample-Size Adjusted BIC (-0.5 * BIC).}
#'       \item{\code{CAIC}}{Consistent AIC.}
#'       \item{\code{AWE}}{Approximate Weight of Evidence.}
#'       \item{\code{SABIC}}{Sample-Size Adjusted BIC (alternative formulation).}
#'       \item{\code{call}}{Original function call that generated the fit indices.}
#'       \item{\code{arguments}}{List containing input arguments (includes original model object).}
#'     }}
#'
#'   \item{\strong{\code{compare.model}}}{Model comparison results. Available components:
#'     \describe{
#'       \item{\code{npar}}{Named numeric vector with free parameters for each model (\code{model1}, \code{model2}).}
#'       \item{\code{entropy}}{Named numeric vector with entropy values for each model.}
#'       \item{\code{AvePP}}{List of average posterior probabilities per class/profile for each model.}
#'       \item{\code{fit.index}}{List of \code{\link[LCPA]{get.fit.index}} objects for both models.}
#'       \item{\code{BF}}{Bayes Factor comparing models (based on SIC differences).}
#'       \item{\code{LRT.obj}}{Standard likelihood ratio test results (requires nested models).}
#'       \item{\code{LRT.VLMR.obj}}{Vuong-Lo-Mendell-Rubin adjusted likelihood ratio test results.}
#'       \item{\code{LRT.Bootstrap.obj}}{Parametric bootstrap likelihood ratio test results (if \code{n.Bootstrap > 0}).}
#'       \item{\code{call}}{The matched function call used for comparison.}
#'       \item{\code{arguments}}{List containing original input arguments (\code{object1}, \code{object2}, \code{n.Bootstrap}).}
#'     }}
#'
#'   \item{\strong{\code{SE}}}{Standard error estimation results. Available components:
#'     \describe{
#'       \item{\code{se}}{List containing standard errors for parameters (components depend on model type).}
#'       \item{\code{vcov}}{Variance-covariance matrix (only for \code{method="Obs"}).}
#'       \item{\code{hessian}}{Observed information matrix (only for \code{method="Obs"}).}
#'       \item{\code{diagnostics}}{Method-specific diagnostic information (e.g., estimation method).}
#'       \item{\code{call}}{Function call that generated the object.}
#'       \item{\code{arguments}}{Input arguments used in estimation.}
#'       \item{\code{means}}{Standard errors for profile means (LPA models only — accessed via \code{se} list).}
#'       \item{\code{covs}}{Standard errors for covariance parameters (LPA models only — accessed via \code{se} list).}
#'       \item{\code{P.Z}}{Standard errors for class proportions (both LCA/LPA — accessed via \code{se} list).}
#'       \item{\code{par}}{Standard errors for conditional probabilities (LCA models only — accessed via \code{se} list).}
#'     }}
#' }
#'
#' @section Usage Notes:
#' \itemize{
#'   \item For \code{LCA}, \code{LPA}, \code{LCPA}, and \code{LTA} objects, components reflect \emph{estimated} parameters.
#'   \item For \code{sim.LCA}, \code{sim.LPA}, and \code{sim.LTA} objects, components reflect \emph{true} data-generating parameters.
#'   \item In \code{SE} objects:
#'     \itemize{
#'       \item Top-level components like \code{vcov} and \code{hessian} are only available when \code{method = "Obs"}.
#'         Requesting them under \code{Bootstrap} triggers a warning and returns \code{NULL}.
#'       \item Parameter-specific SEs (e.g., \code{means}, \code{par}) are stored within the \code{se} list.
#'         You can extract them directly by name (e.g., \code{extract(se_obj, "means")}).
#'       \item Attempting to extract unavailable parameter SEs (e.g., \code{par} from an LPA model) triggers an error with available options.
#'     }
#'   \item For \code{fit.index} and \code{compare.model} objects, valid components are dynamically determined from the object’s names.
#'   \item All methods ignore additional arguments (\code{...}).
#' }
#'
#' @examples
#' set.seed(123)
#'
#' # Simulate LPA data: 500 observations, 3 continuous variables, 2 latent profiles
#' # Constraint "E0": Equal variances across classes, zero covariances
#' data.obj <- sim.LPA(N = 500, I = 3, L = 2, constraint = "E0")
#'
#' # Extract the simulated response matrix (N x I) for model fitting
#' response <- extract(data.obj, "response")
#'
#' # Extract the TRUE covariance matrices (I x I x L array)
#' extract(data.obj, "covs")
#'
#' # Fit an LPA model to the simulated data using the SAME constraint ("E0")
#' fit_E0 <- LPA(response, L = 2, constraint = "E0")
#'
#' # Extract the ESTIMATED covariance matrices from the fitted model
#' extract(fit_E0, "covs")
#'
#' # Simulate LCA data: 30 observations, 5 categorical items, 3 latent classes
#' sim_data <- sim.LCA(N = 30, I = 5, L = 3)
#'
#' # Extract the TRUE conditional probability array
#' extract(sim_data, "par")
#'
#' @name extract
#' @export
extract <- function(object, what, ...) {
  UseMethod("extract")
}

#' @describeIn extract Extract fields from a \code{LCA} object
#' @export
extract.LCA <- function(object, what, ...) {
  choices <- c("params", "par", "P.Z", "npar", "Log.Lik", "AIC", "BIC",
               "best_BIC", "P.Z.Xn", "Z", "probability", "Log.Lik.history",
               "Log.Lik.nrep", "model", "call", "arguments")

  if (!what %in% choices) {
    stop(sprintf("'%s' is not a valid field for LCA objects. Choose from: %s",
                 what, paste(choices, collapse = ", ")), call. = FALSE)
  }

  switch(what,
         params          = object$params,
         par             = object$params$par,
         P.Z             = object$params$P.Z,
         npar            = object$npar,
         Log.Lik         = object$Log.Lik,
         AIC             = object$AIC,
         BIC             = object$BIC,
         best_BIC        = object$best_BIC,
         P.Z.Xn          = object$P.Z.Xn,
         Z               = object$Z,
         probability     = object$probability,
         Log.Lik.history = object$Log.Lik.history,
         Log.Lik.nrep    = object$Log.Lik.nrep,
         model           = object$model,
         call            = object$call,
         arguments       = object$arguments)
}

#' @describeIn extract Extract fields from a \code{LPA} object
#' @export
extract.LPA <- function(object, what, ...) {
  choices <- c("params", "means", "covs", "P.Z", "npar", "Log.Lik", "AIC", "BIC",
               "best_BIC", "P.Z.Xn", "Z", "Log.Lik.history", "Log.Lik.nrep",
               "model", "call", "arguments", "constraint")

  if (!what %in% choices) {
    stop(sprintf("'%s' is not a valid field for LPA objects. Choose from: %s",
                 what, paste(choices, collapse = ", ")), call. = FALSE)
  }

  switch(what,
         params          = object$params,
         means           = object$params$means,
         covs            = object$params$covs,
         P.Z             = object$params$P.Z,
         npar            = object$npar,
         Log.Lik         = object$Log.Lik,
         AIC             = object$AIC,
         BIC             = object$BIC,
         best_BIC        = object$best_BIC,
         P.Z.Xn          = object$P.Z.Xn,
         Z               = object$Z,
         Log.Lik.history = object$Log.Lik.history,
         Log.Lik.nrep    = object$Log.Lik.nrep,
         model           = object$model,
         call            = object$call,
         arguments       = object$arguments,
         constraint      = object$arguments$constraint)
}

#' @describeIn extract Extract fields from a \code{LCPA} object
#' @export
extract.LCPA <- function(object, what, ...) {
  choices <- c(
    "beta",                   # Initial class coefficients (p1 x L matrix)
    "beta.se",                # Standard errors for beta
    "beta.Z.sta",             # Z-statistics for beta
    "beta.p.value.tail1",     # One-tailed p-values for beta
    "beta.p.value.tail2",     # Two-tailed p-values for beta
    "P.Z.Xn",                 # posterior probabilities (N x L)
    "P.Z",                    # prior proportions (length L)
    "Z",                      # modal class assignments (length N)
    "npar",                   # Number of free parameters
    "Log.Lik",                # Log-likelihood
    "AIC",                    # AIC
    "BIC",                    # BIC
    "iterations",             # Optimization iterations in Step 3
    "coveraged",              # Logical: did optimization converge early?
    "params",                 # Step 1 model parameters (LCA/LPA output)
    "call",                   # Function call
    "arguments"               # Input arguments list
  )

  if (!what %in% choices) {
    stop(sprintf("'%s' is not a valid field for LTA objects. Choose from: %s",
                 what, paste(choices, collapse = ", ")), call. = FALSE)
  }

  switch(what,
         beta                   = object$beta,
         beta.se                = object$beta.se,
         beta.Z.sta             = object$beta.Z.sta,
         beta.p.value.tail1     = object$beta.p.value.tail1,
         beta.p.value.tail2     = object$beta.p.value.tail2,
         P.Z.Xn                = object$P.Z.Xn,
         P.Z                   = object$P.Z,
         Z                     = object$Z,
         npar                   = object$npar,
         Log.Lik                = object$Log.Lik,
         AIC                    = object$AIC,
         BIC                    = object$BIC,
         iterations             = object$iterations,
         coveraged              = object$coveraged,
         params                 = object$params,
         call                   = object$call,
         arguments              = object$arguments)
}

#' @describeIn extract Extract fields from a \code{LTA} object
#' @export
extract.LTA <- function(object, what, ...) {
  choices <- c(
    "beta",                   # Initial class coefficients (p1 x L matrix)
    "gamma",                   # Transition coefficients (nested list)
    "beta.se",                # Standard errors for beta
    "gamma.se",                # Standard errors for gamma
    "beta.Z.sta",             # Z-statistics for beta
    "gamma.Z.sta",             # Z-statistics for gamma
    "beta.p.value.tail1",     # One-tailed p-values for beta
    "gamma.p.value.tail1",     # One-tailed p-values for gamma
    "beta.p.value.tail2",     # Two-tailed p-values for beta
    "gamma.p.value.tail2",     # Two-tailed p-values for gamma
    "P.Z.Xns",                # List of posterior probabilities per time (N x L)
    "P.Zs",                   # List of prior proportions per time (length L)
    "Zs",                     # List of modal class assignments per time (length N)
    "npar",                   # Number of free parameters
    "Log.Lik",                # Log-likelihood
    "AIC",                    # AIC
    "BIC",                    # BIC
    "iterations",             # Optimization iterations in Step 3
    "coveraged",              # Logical: did optimization converge early?
    "params",                 # Step 1 model parameters (LCA/LPA output)
    "call",                   # Function call
    "arguments"               # Input arguments list
  )

  if (!what %in% choices) {
    stop(sprintf("'%s' is not a valid field for LTA objects. Choose from: %s",
                 what, paste(choices, collapse = ", ")), call. = FALSE)
  }

  switch(what,
         beta                   = object$beta,
         gamma                   = object$gamma,
         beta.se                = object$beta.se,
         gamma.se                = object$gamma.se,
         beta.Z.sta             = object$beta.Z.sta,
         gamma.Z.sta             = object$gamma.Z.sta,
         beta.p.value.tail1     = object$beta.p.value.tail1,
         gamma.p.value.tail1     = object$gamma.p.value.tail1,
         beta.p.value.tail2     = object$beta.p.value.tail2,
         gamma.p.value.tail2     = object$gamma.p.value.tail2,
         P.Z.Xns                = object$P.Z.Xns,
         P.Zs                   = object$P.Zs,
         Zs                     = object$Zs,
         npar                   = object$npar,
         Log.Lik                = object$Log.Lik,
         AIC                    = object$AIC,
         BIC                    = object$BIC,
         iterations             = object$iterations,
         coveraged              = object$coveraged,
         params                 = object$params,
         call                   = object$call,
         arguments              = object$arguments)
}

#' @describeIn extract Extract fields from a \code{sim.LCA} object
#' @export
extract.sim.LCA <- function(object, what, ...) {
  choices <- c("response", "par", "Z", "P.Z", "poly.value", "P.Z.Xn", "call", "arguments")

  if (!what %in% choices) {
    stop(sprintf("'%s' is not a valid field for sim.LCA objects. Choose from: %s",
                 what, paste(choices, collapse = ", ")), call. = FALSE)
  }

  switch(what,
         response   = object$response,
         par        = object$par,
         Z          = object$Z,
         P.Z        = object$P.Z,
         poly.value = object$poly.value,
         P.Z.Xn     = object$P.Z.Xn,
         call       = object$call,
         arguments  = object$arguments)
}

#' @describeIn extract Extract fields from a \code{sim.LPA} object
#' @export
extract.sim.LPA <- function(object, what, ...) {
  choices <- c("response", "means", "covs", "P.Z.Xn", "P.Z", "Z",
               "constraint", "call", "arguments")

  if (!what %in% choices) {
    stop(sprintf("'%s' is not a valid field for sim.LPA objects. Choose from: %s",
                 what, paste(choices, collapse = ", ")), call. = FALSE)
  }

  switch(what,
         response   = object$response,
         means      = object$means,
         covs       = object$covs,
         P.Z.Xn     = object$P.Z.Xn,
         P.Z        = object$P.Z,
         Z          = object$Z,
         constraint = object$constraint,
         call       = object$call,
         arguments  = object$arguments)
}

#' @describeIn extract Extract fields from a \code{sim.LTA} object
#' @export
extract.sim.LTA <- function(object, what, ...) {
  choices <- c(
    "responses", "Zs", "P.Zs",  "par", "means", "covs", "poly.value", "rate", "covariates",
    "beta",  "gamma", "call", "arguments")

  if (!what %in% choices) {
    stop(sprintf("'%s' is not a valid field for sim.LTA objects. Choose from: %s",
                 what, paste(choices, collapse = ", ")), call. = FALSE)
  }

  switch(what,
         responses   = object$responses,
         Zs          = object$Zs,
         P.Zs        = object$P.Zs,
         par         = object$par,
         means       = object$means,
         covs        = object$covs,
         poly.value  = object$poly.value,
         rate        = object$rate,
         covariates  = object$covariates,
         beta        = object$beta,
         gamma       = object$gamma,
         call        = object$call,
         arguments   = object$arguments)
}

#' @describeIn extract Extractor method for \code{fit.index} objects
#' @export
extract.fit.index <- function(object, what, ...) {
  valid_components <- names(object)

  if (!what %in% valid_components) {
    stop(
      sprintf(
        "Invalid component '%s' requested for 'fit.index' object.\nValid components are: %s",
        what,
        paste(valid_components, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  return(object[[what]])
}

#' @describeIn extract Extract fields from a \code{compare.model} object
#' @export
extract.compare.model <- function(object, what, ...) {
  valid_components <- names(object)

  if (!what %in% valid_components) {
    stop(
      sprintf(
        "Invalid component '%s' requested for 'compare.model' object.\nValid components are: %s",
        what,
        paste(valid_components, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  return(object[[what]])
}

#' @describeIn extract Extract fields from a \code{SE} object
#' @export
extract.SE <- function(object, what, ...) {
  # Top-level components (directly in SE object)
  top_components <- c("se", "vcov", "hessian", "diagnostics", "call", "arguments")

  # Parameter-level components (inside se list)
  param_components <- c("means", "covs", "P.Z", "par")

  # Check if requested component exists at top level
  if (what %in% top_components) {
    # Handle method-specific warnings
    if (what %in% c("vcov", "hessian") && object$diagnostics$method == "Bootstrap") {
      warning(sprintf(
        "Component '%s' is NULL for Bootstrap method. Only available for 'Obs' method.",
        what
      ), call. = FALSE)
    }
    return(object[[what]])
  }

  # Check if requested component exists in se list
  if (what %in% param_components) {
    if (!exists("se", where = object) || is.null(object$se)) {
      stop("Internal error: 'se' component missing from SE object", call. = FALSE)
    }

    if (!what %in% names(object$se)) {
      model_type <- if (!is.null(object$se$means)) "LPA" else if (!is.null(object$se$par)) "LCA" else "unknown"
      available <- names(object$se)
      stop(sprintf(
        "Component '%s' not available for this %s model.\nAvailable parameter components: %s",
        what,
        model_type,
        paste(available, collapse = ", ")
      ), call. = FALSE)
    }

    return(object$se[[what]])
  }

  # Invalid component
  valid_all <- c(top_components, param_components)
  stop(
    sprintf(
      "Invalid component '%s' requested for 'SE' object.\nValid components are: %s",
      what,
      paste(valid_all, collapse = ", ")
    ),
    call. = FALSE
  )
}
