#' Latent Class/Profile Analysis with Covariates
#'
#' Implements the three-step estimation method (Vermunt, 2010; Liang et al., 2023) for latent class/profile analysis
#' with covariates, treating latent class membership as an observed variable with measurement error.
#' This is mathematically equivalent to a latent transition analysis (\code{\link[LCPA]{LTA}}) with \code{times=1}.
#'
#' @param response A matrix or data frame of observed responses.
#'                  Rows of the matrix represent individuals/participants/observations (\eqn{N}), columns of the
#'                  matrix represent observed indicators/items/variables (\eqn{I}).
#'                  For \code{type = "LCA"}: indicators must be binary or categorical (coded as integers starting from 0).
#'                  For \code{type = "LPA"}: indicators must be continuous (numeric), and the response matrix must be
#'                  standardized using \code{\link[base]{scale}} or \code{\link[LCPA]{normalize}} prior to input.
#' @param L Integer scalar. Number of latent classes/profiles. Must satisfy \eqn{L \geq 2}.
#' @param ref.class Integer \eqn{L \geq ref.class \geq 1}. Specifies which latent class to use as the reference category.
#'                  Default is \code{L} (last class). Coefficients for the reference class are fixed to zero.
#'                  When \code{is.sort=TRUE}, classes are first ordered by decreasing \code{P.Z} (class 1 has highest probability),
#'                  then \code{ref.class} refers to the position in this sorted order.
#' @param type Character string. Specifies the type of latent variable model for Step 1:
#'             \itemize{
#'               \item \code{"LCA"} — Latent Class Analysis for categorical indicators.
#'               \item \code{"LPA"} — Latent Profile Analysis for continuous indicators.
#'             }
#'             See \code{\link[LCPA]{LCA}} and \code{\link[LCPA]{LPA}} for details.
#' @param covariate Optional. A matrix or data frame of covariates for modeling latent class membership.
#'                  Must include an intercept column (all 1s) as the first column.
#'                  If \code{NULL} (default), only intercept terms are used (i.e., no covariates).
#'                  Dimension is \eqn{N \times p} where \eqn{p} is the number of covariates including intercept.
#' @param CEP.error Logical. If \code{TRUE} (recommended), incorporates classification uncertainty via
#'                  estimated Classification Error Probability (\code{\link[LCPA]{get.CEP}}) matrices from Step 1. If \code{FALSE},
#'                  uses identity CEP matrices (equivalent to naive modal assignment; introduces bias).
#' @param par.ini Specification for parameter initialization. Options include:
#'   \itemize{
#'     \item \code{"random"}: Completely random initialization (default).
#'     \item \code{"kmeans"}: Initializes parameters via K-means clustering on observed data (McLachlan & Peel, 2004).
#'     \item A \code{list} for LCA containing:
#'       \describe{
#'         \item{\code{par}}{An \eqn{L \times I \times K_{\max}} array of initial conditional probabilities for
#'                           each latent class, indicator, and response category (where \eqn{K_{\max}} is the maximum
#'                           number of categories across indicators).}
#'         \item{\code{P.Z}}{A numeric vector of length \eqn{L} specifying initial prior probabilities for latent classes.}
#'       }
#'     \item A \code{list} for LPA containing:
#'       \describe{
#'         \item{\code{means}}{An \eqn{L \times I} matrix of initial mean vectors for each profile.}
#'         \item{\code{covs}}{An \eqn{I \times I \times L} array of initial covariance matrices for each profile.}
#'         \item{\code{P.Z}}{A numeric vector of length \eqn{L} specifying initial prior probabilities for profiles.}
#'       }
#'   }
#' @param params Optional \code{list} of pre-estimated Step 1 parameters. If \code{NULL} (default),
#'               Step 1 models are estimated internally. If provided, no LCA or LPA parameter estimation
#'               will be performed; instead, the parameters provided in \code{params} will be used as
#'               fixed values. Additionally, \code{params} must contain:
#'               \itemize{
#'                 \item A \code{list} for LCA containing:
#'                    \describe{
#'                      \item{\code{par}}{An \eqn{L \times I \times K_{\max}} array of initial conditional probabilities for
#'                                        each latent class, indicator, and response category (where \eqn{K_{\max}} is the maximum
#'                                        number of categories across indicators).}
#'                      \item{\code{P.Z}}{A numeric vector of length \eqn{L} specifying initial prior probabilities for latent classes.}
#'                    }
#'                 \item A \code{list} for LPA containing:
#'                    \describe{
#'                      \item{\code{means}}{An \eqn{L \times I} matrix of initial mean vectors for each profile.}
#'                      \item{\code{covs}}{An \eqn{I \times I \times L} array of initial covariance matrices for each profile.}
#'                      \item{\code{P.Z}}{A numeric vector of length \eqn{L} specifying initial prior probabilities for profiles.}
#'                    }
#'               }
#' @param is.sort A logical value. If \code{TRUE} (Default), the latent classes will be ordered in descending
#'                order according to \code{P.Z}. All other parameters will be adjusted accordingly
#'                based on the reordered latent classes.
#' @param constraint Character (LPA only). Specifies structure of within-class covariance matrices:
#'                  \itemize{
#'                    \item \code{"VV"} — Class-varying variances and covariances (unconstrained; default).
#'                    \item \code{"EE"} — Equal variances and covariances across all classes (homoscedastic).
#'                  }
#' @param method Character. Estimation algorithm for Step 1 models:
#'               \itemize{
#'                 \item \code{"EM"} — Expectation-Maximization (default; robust and widely used).
#'                 \item \code{"Mplus"} — Interfaces with Mplus software (requires external installation).
#'                 \item \code{"NNE"} — Neural Network Estimator (experimental; uses transformer + simulated
#'                                      annealing, more reliable than both \code{"EM"} and \code{"Mplus"}).
#'               }
#' @param tol Convergence tolerance for log-likelihood difference (default: 1e-4).
#' @param lower The upper bound for the estimation of regression coefficients, default is -10
#' @param upper The lower bound for the estimation of regression coefficients, default is 10
#' @param method.SE Character. Method for estimating standard errors of parameter estimates:
#'               \itemize{
#'                 \item \code{"Obs"} — Approximates the observed information matrix via numerical differentiation (Richardson's method).
#'                       Standard errors are obtained from the inverse Hessian. May fail or be unreliable in small samples or with
#'                       complex likelihood surfaces.
#'                 \item \code{"Bootstrap"} — Uses nonparametric bootstrap resampling to estimate empirical sampling variability.
#'                       More robust to model misspecification and small-sample bias. Computationally intensive but recommended when
#'                       asymptotic assumptions are questionable.
#'               }
#'               Default is \code{"Bootstrap"}.
#' @param n.Bootstrap Integer. Number of bootstrap replicates used when \code{method.SE = "Bootstrap"}.
#'                    Default is 100. McLachlan & Peel (2004) suggest that 50–100 replicates often provide adequate accuracy
#'                    for practical purposes, though more (e.g., 500–1000) may be preferred for publication-quality inference.
#'                    Each replicate involves re-estimating the full three-step LTA model on a resampled dataset.
#' @param maxiter Maximum number of iterations for optimizing the regression coefficients. Default: 5000.
#' @param nrep Integer controlling replication behavior:
#'   \itemize{
#'     \item If \code{par.ini = "random"}, number of random initializations.
#'     \item If \code{par.ini = "kmeans"}, number of K-means runs for initialization.
#'     \item For \code{method="Mplus"}, controls number of random starts in Mplus via \code{STARTS} option.
#'     \item Best solution is selected by log-likelihood/BIC across replications.
#'     \item Ignored for user-provided initial parameters.
#'   }
#' @param starts Number of random initializations to explore during warm-up phase (default: 100).
#' @param maxiter.wa Maximum number of training iterations allowed per warm-up run. After completion,
#'                   the top \code{nrep} solutions (by log-likelihood) are promoted to final training phase (default: 20).
#' @param vis Logical. If \code{TRUE}, displays progress information during estimation (default: \code{TRUE}).
#' @param control.EM List of control parameters for EM algorithm:
#'   \describe{
#'     \item{\code{maxiter}}{Maximum iterations (default: 2000).}
#'     \item{\code{tol}}{Convergence tolerance for log-likelihood difference (default: 1e-4).}
#'   }
#' @param control.NNE List of control parameters for NNE algorithm:
#'   \describe{
#'     \item{\code{hidden.layers}}{Integer vector specifying layer sizes in fully-connected network (default: \code{c(16,16)}).}
#'     \item{\code{activation.function}}{Activation function (e.g., \code{"tanh"}, default: \code{"tanh"}).}
#'     \item{\code{use.attention}}{Whether to enable the self-attention mechanism (i.e., transformer encoder) (default: \code{TRUE}).}
#'     \item{\code{d.model}}{Dimensionality of transformer encoder embeddings (default: 8).}
#'     \item{\code{nhead}}{Number of attention heads in transformer (default: 2).}
#'     \item{\code{dim.feedforward}}{Dimensionality of transformer feedforward network (default: 16).}
#'     \item{\code{eps}}{Small constant for numerical stability (default: 1e-8).}
#'     \item{\code{lambda}}{A factor for slight regularization of all parameters (default: 1e-5).}
#'     \item{\code{initial.temperature}}{Initial temperature for simulated annealing (default: 1000).}
#'     \item{\code{cooling.rate}}{Cooling rate per iteration in simulated annealing (default: 0.5).}
#'     \item{\code{maxiter.sa}}{Maximum iterations for simulated annealing (default: 1000).}
#'     \item{\code{threshold.sa}}{Minimum temperature threshold for annealing (default: 1e-10).}
#'     \item{\code{maxiter}}{Maximum training epochs (default: 1000).}
#'     \item{\code{maxiter.early}}{Patience parameter for early stopping (default: 50).}
#'     \item{\code{maxcycle}}{Maximum cycles for optimization (default: 10).}
#'     \item{\code{lr}}{Learning rate, controlling the step size of neural network parameter updates (default: 0.025).}
#'     \item{\code{scheduler.patience}}{Patience for learning rate decay (if the loss function does not improve for more than `patience` consecutive epochs, the learning rate will be reduced) (default: 10).}
#'     \item{\code{scheduler.factor}}{Learning rate decay factor; the new learning rate equals the original learning rate multiplied by `scheduler.factor` (default: 0.70).}
#'     \item{\code{plot.interval}}{Interval (in epochs) for plotting training diagnostics (default: 100).}
#'     \item{\code{device}}{Specifies the hardware device; can be \code{"CPU"} (default) or \code{"GPU"}. If the GPU is not available, it automatically falls back to CPU.}
#'   }
#' @param control.Mplus List of control parameters for Mplus estimation:
#'   \describe{
#'     \item{\code{maxiter}}{Maximum iterations for Mplus optimization (default: 2000).}
#'     \item{\code{tol}}{Convergence tolerance for log-likelihood difference (default: 1e-4).}
#'     \item{\code{files.path}}{Character string specifying the directory path where Mplus will write its intermediate files
#'                              (e.g., \code{.inp} model input, \code{.dat} data file, \code{.out} output, and saved posterior probabilities).
#'                              This argument is \strong{required} — if \code{NULL} (default), the function throws an error.
#'                              The specified directory must exist and be writable; if it does not exist, the function attempts to create it recursively.
#'                              A unique timestamped subdirectory (e.g., \code{"Mplus_LPA_YYYY-MM-DD_HH-MM-SS"} or
#'                              \code{"Mplus_LCA_YYYY-MM-DD_HH-MM-SS"}) will be created within this path
#'                              to store all run-specific files and avoid naming conflicts. See in \code{\link[LCPA]{LCA}} and \code{\link[LCPA]{LPA}}.}
#'     \item{\code{files.clean}}{Logical. If \code{TRUE} (default), all intermediate files and the temporary working directory
#'                               created for this run are deleted upon successful completion or error exit (via \code{on.exit()}).
#'                               If \code{FALSE}, all generated files are retained in \code{files.path} (or the auto-generated temp dir)
#'                               for inspection or debugging. Note: when \code{files.path = NULL}, even if \code{files.clean = FALSE},
#'                               the temporary directory may still be cleaned up by the system later — for guaranteed persistence,
#'                               specify a custom \code{files.path}.}
#'   }
#'
#' @return An object of class \code{LCPA}, a named list containing:
#' \describe{
#'   \item{\code{beta}}{Matrix of size \eqn{p \times L}. Coefficients for class membership multinomial logit model.
#'     Columns 1 to \eqn{L-1} are free parameters; column \eqn{L} (reference class) is constrained to \eqn{\boldsymbol{\beta}_L = \mathbf{0}}.}
#'   \item{\code{beta.se}}{Standard errors for \code{beta} (if Hessian is invertible). Same dimensions as \code{beta}.
#'     May contain \code{NA} if variance-covariance matrix is not positive definite.}
#'   \item{\code{beta.Z.sta}}{Z-statistics for testing null hypothesis that each beta coefficient equals zero.
#'     Computed as \code{beta / beta.se}. Same structure as \code{beta}.}
#'   \item{\code{beta.p.value.tail1}}{One-tailed p-values based on standard normal distribution: \eqn{P(Z < -|z|)}.
#'     Useful for directional hypotheses. Same structure as \code{beta}.}
#'   \item{\code{beta.p.value.tail2}}{Two-tailed p-values: \eqn{2 \times P(Z < -|z|)}.
#'     Standard test for non-zero effect. Same structure as \code{beta}.}
#'   \item{\code{P.Z.Xn}}{Matrix of size \eqn{N \times L} of posterior class probabilities
#'     \eqn{P(Z_n=l \mid \mathbf{X}_n)} for each individual \eqn{n} and class \eqn{l}.}
#'   \item{\code{P.Z}}{Vector of length \eqn{L} containing prior class proportions
#'     \eqn{P(Z = l)} estimated at Step 1.}
#'   \item{\code{Z}}{Vector of length \eqn{N} containing modal class assignments
#'     (MAP classifications) \eqn{\hat{z}_n} for each individual.}
#'   \item{\code{npar}}{Number of free parameters in the model (depends on covariates).}
#'   \item{\code{Log.Lik}}{Observed-data log-likelihood value at convergence.}
#'   \item{\code{Log.Lik.history}}{Vector tracking log-likelihood at each iteration.}
#'   \item{\code{AIC}}{Akaike Information Criterion value.}
#'   \item{\code{BIC}}{Bayesian Information Criterion value.}
#'   \item{\code{iterations}}{Integer. Number of optimization iterations in Step 3.}
#'   \item{\code{coveraged}}{Logical. \code{TRUE} if optimization terminated before reaching \code{maxiter} (suggesting convergence).
#'                           Note: This is a heuristic indicator; formal convergence diagnostics should check Hessian properties.}
#'   \item{\code{params}}{List. Step 1 model parameters (output from \code{LCA()} or \code{LPA()}).}
#'   \item{\code{call}}{The matched function call.}
#'   \item{\code{arguments}}{List of all input arguments passed to the function (useful for reproducibility).}
#' }
#'
#' @section Methodology Overview:
#' The three-step procedure follows the same principles as LTA but for a single time point:
#'
#' Step 1 — Unconditional Latent Class/Profile Model:
#' Fit an unconditional LCA or LPA model (ignoring covariates). Obtain posterior class membership probabilities
#' \eqn{P(Z_n=l \mid \mathbf{X}_n)} for each individual \eqn{n} and class \eqn{l} using Bayes' theorem.
#'
#' Step 2 — Classification Error Probabilities (equal to \code{\link[LCPA]{get.CEP}}):
#' Compute the \eqn{L \times L} CEP matrix where element \eqn{(k,l)} estimates:
#' \deqn{
#'   \text{CEP}(k,l) = P(\hat{Z}_n = l \mid Z_n = k)
#' }
#' using a non-parametric approximation:
#' \deqn{
#'   \widehat{\text{CEP}}(k,l) = \frac{ \sum_{n=1}^N \mathbb{I}(\hat{z}_n = l) \cdot P(Z_n=k \mid \mathbf{X}_n) }{ \sum_{n=1}^N P(Z_n=k \mid \mathbf{X}_n) }
#' }
#' where \eqn{\hat{z}_n} is the modal class assignment.
#'
#' Step 3 — Class Membership Model with Measurement Error Correction:
#' Estimate the multinomial logit model for class membership:
#' \deqn{
#' P(Z_n = l \mid \mathbf{X}_n) = \frac{\exp(\boldsymbol{\beta}_l^\top \mathbf{X}_n)}{\sum_{k=1}^L \exp(\boldsymbol{\beta}_k^\top \mathbf{X}_n)}
#' }
#' where \eqn{\mathbf{X}_n = (1, W_{n1}, \dots, W_{nM})^\top} is the covariate vector for individual \eqn{n}
#' (with intercept as first column), and \eqn{\boldsymbol{\beta}_l = (\beta_{l0}, \beta_{l1}, \dots, \beta_{lM})^\top}
#' contains intercept and regression coefficients. Class \eqn{L} is the reference category
#' (\eqn{\boldsymbol{\beta}_L = \mathbf{0}}).
#'
#' The observed-data likelihood integrates over latent classes:
#' \deqn{
#' \log \mathcal{L}(\boldsymbol{\beta}) =
#' \sum_{n=1}^N \log \left[
#'   \sum_{l=1}^L \text{CEP}(l, \hat{z}_n) \cdot P(Z_n=l \mid \mathbf{X}_n)
#' \right]
#' }
#' Parameters \eqn{\boldsymbol{\beta}} are estimated via maximum likelihood using the BOBYQA algorithm.
#'
#' @section Important Implementation Details:
#' \itemize{
#'   \item Reference Class: Coefficients for the reference class (\code{ref.class}) are ALWAYS fixed to
#'                          zero (\eqn{\boldsymbol{\beta}_{ref.class} = \mathbf{0}}) in the multinomial
#'                          logit model.
#'   \item CEP Matrices: When \code{CEP.error = TRUE}, misclassification probabilities are estimated
#'                       non-parametrically using Step 1 posterior probabilities. This corrects for
#'                       classification uncertainty. See in \code{\link[LCPA]{get.CEP}}.
#'   \item Covariate Requirements: Covariate matrix MUST include an intercept column (all 1s) as the first
#'                                 column. Dimensions must be \eqn{N \times (M+1)}, where \eqn{M} represents
#'                                 the number of covariate and \eqn{1} is the Intercept.
#'   \item \strong{Optimization & Standard Errors}:
#'         \itemize{
#'           \item Step 3 uses BOBYQA algorithm (\code{nloptr::nloptr}) for stable optimization with box constraints.
#'           \item For \code{method.SE = "Obs"}: Standard errors derived from inverse Hessian (\code{\link[numDeriv]{hessian}}). If Hessian is singular:
#'                 \itemize{
#'                   \item Uses Moore-Penrose pseudoinverse (\code{\link[MASS]{ginv}})
#'                   \item Sets negative variances to \code{NA}
#'                 }
#'           \item For \code{method.SE = "Bootstrap"}: Each replicate independently re-estimates Steps 1-3.
#'                 Failed bootstrap runs yield \code{NA} in SEs and derived statistics. Progress messages include
#'                 replicate index and optimization diagnostics.
#'         }
#'   \item \strong{Computational Notes}:
#'         \itemize{
#'           \item Step 1 complexity increases with \eqn{L} and \eqn{I}.
#'           \item Bootstrap is computationally intensive: 100 replicates = 100 full re-estimations of Steps 1-3.
#'         }
#'   \item \strong{Bootstrap Reproducibility}: Always set a seed (e.g., \code{set.seed(123)}) before
#'                                             calling \code{LCPA()} when using \code{method.SE = "Bootstrap"}.
#'                                             Monitor convergence in bootstrap runs via progress messages.
#' }
#'
#' @references
#' Vermunt, J. K. (2010). Latent class modeling with covariates: Two improved three-step approaches. Political Analysis, 18(4), 450–469. https://doi.org/10.1093/pan/mpq025
#'
#' Liang, Q., la Torre, J. d., & Law, N. (2023). Latent Transition Cognitive Diagnosis Model With Covariates: A Three-Step Approach. Journal of Educational and Behavioral Statistics, 48(6), 690-718. https://doi.org/10.3102/10769986231163320
#'
#'
#' @examples
#' library(LCPA)
#'
#' set.seed(123)
#' N <- 2000  # Sample size
#' L <- 3    # Number of latent classes
#' I <- 6    # Number of indicators
#'
#' # Create covariates (intercept + 2 covariates + 1 interaction)
#'  Intercept = rep(1, N)
#'  X1 <- rnorm(N)
#'  X2 <- rbinom(N, 1, 0.5)
#'  X1.X2 <- X1 * X2
#' covariate <- cbind(Intercept, X1, X2, X1.X2)
#'
#' # Simulate data for LPA
#' sim_data <- sim.LTA(
#'   N = N, I = I, L = L, times = 1, type = "LPA",
#'   covariates = list(covariate), is.sort=TRUE,
#'   beta = matrix(c(
#'    -0.2, 0.0, -0.1,  ## fix reference class to class 2
#'     0.2, 0.0, -0.3,
#'     0.8, 0.0, -0.6,
#'    -0.1, 0.0,  0.3
#'   ), ncol = L, byrow = TRUE)
#' )
#' response <- sim_data$responses[[1]]
#'
#' ## It is strongly recommended to perform the following
#' ## standardization to obtain more stable results when LPA.
#' ## Standardization is not performed here in order to
#' ## compare estimated values with true values.
#' # response <- normalize(response)
#'
#' # Fit cross-sectional LPA with covariates
#' ## fix reference class to class 2
#' # need Mplus
#' \dontrun{
#' fit <- LCPA(
#'   response = response,
#'   L = L, ref.class = 2,
#'   type = "LPA", is.sort=TRUE,
#'   covariate = covariate,
#'   method.SE = "Obs",
#'   CEP.error = TRUE,
#'   method = "Mplus",
#'   control.Mplus = list(files.path = ""),
#'   vis = TRUE
#' )
#' print(fit)
#' }
#'
#' @importFrom nloptr nloptr
#' @importFrom Matrix nearPD
#' @importFrom numDeriv hessian
#' @importFrom stats pnorm
#' @importFrom MASS ginv
#'
#' @export
#'
LCPA <- function(response, L = 2,
                 ref.class = L, type = "LCA",
                 covariate = NULL,
                 CEP.error = TRUE,
                 par.ini = "random",
                 params = NULL, is.sort = TRUE,
                 constraint = "VV",
                 method = "EM", tol = 1e-4,
                 lower=-10, upper=10,
                 method.SE = "Bootstrap", n.Bootstrap=100,
                 maxiter = 5000, nrep = 20,
                 starts = 100, maxiter.wa = 20,
                 vis = TRUE,
                 control.EM = NULL,
                 control.Mplus = NULL,
                 control.NNE = NULL) {

  call <- match.call()

  if (ref.class < 1 || ref.class > L) {
    stop("ref.class must be between 1 and L")
  }

  # Default control parameters
  default_control.EM <- list(maxiter = 2000, tol = 1e-4)
  default_control.Mplus <- list(maxiter = 2000, tol = 1e-4, files.path = NULL, files.clean = TRUE)
  default_control.NNE <- list(
    hidden.layers = c(16, 16),
    activation.function = "tanh",
    use.attention=TRUE,
    d.model = 8,
    nhead = 2,
    dim.feedforward = 16,
    eps = 1e-8,
    lambda = 1e-5,
    initial.temperature = 1000,
    cooling.rate = 0.5,
    maxiter.sa = 1000,
    threshold.sa = 1e-10,
    maxiter = 1000,
    maxiter.early = 50,
    maxcycle = 10,
    lr = 0.025,
    scheduler.patience = 10,
    scheduler.factor = 0.80,
    plot.interval = 100,
    device = "CPU"
  )

  merge_and_clean_control <- function(user_control, default_control) {
    if (is.null(user_control)) return(default_control)
    merged <- modifyList(default_control, user_control)
    merged[names(default_control)]
  }

  control.EM <- merge_and_clean_control(control.EM, default_control.EM)
  control.Mplus <- merge_and_clean_control(control.Mplus, default_control.Mplus)
  control.NNE <- merge_and_clean_control(control.NNE, default_control.NNE)

  # Convert to matrix if needed
  response <- as.matrix(response)
  N <- nrow(response)
  I <- ncol(response)

  if(vis){
    cat(paste0("Starting the first step for ", type, " ...\n\n"))
  }

  # Step 1: Unconditional LCA/LPA
  if (is.null(params)) {
    if (type == "LCA") {
      adj <- adjust.response(response)
      response <- adj$response
      LCPA.obj <- LCA(response, L = L, par.ini = par.ini,
                      method = method, is.sort = is.sort, nrep = nrep,
                      starts = starts, maxiter.wa = maxiter.wa,
                      vis = vis,
                      control.EM = control.EM,
                      control.Mplus = control.Mplus,
                      control.NNE = control.NNE)
    } else {
      LCPA.obj <- LPA(response, L = L, par.ini = par.ini,
                      constraint = constraint,
                      method = method, is.sort = is.sort, nrep = nrep,
                      starts = starts, maxiter.wa = maxiter.wa,
                      vis = vis,
                      control.EM = control.EM,
                      control.Mplus = control.Mplus,
                      control.NNE = control.NNE)
    }
    params <- LCPA.obj$params
  }

  if(vis){
    cat(paste0("Starting the second step for ", type, " ...\n\n"))
  }
  if (type == "LCA") {
    P.Z.Xns <- list(get.P.Z.Xn.LCA(response = response, par = params$par, P.Z=params$P.Z))
  } else {
    P.Z.Xns <- list(get.P.Z.Xn.LPA(response = response, means = params$means, covs = params$covs, P.Z=params$P.Z))
  }
  P.Zs <- list(colSums(P.Z.Xns[[1]]) / sum(P.Z.Xns[[1]]))
  Zs <- list(apply(P.Z.Xns[[1]], 1, which.max))

  if (CEP.error) {
    CEP <- get.CEP(P.Z.Xns, time.cross=FALSE)
  } else {
    CEP <- list(diag(L))
  }

  if (is.null(covariate)) {
    covariate <- matrix(1, nrow = N, ncol = 1)
    colnames(covariate) <- "Intercept"
  } else {
    covariate <- as.matrix(covariate)
    if (ncol(covariate) == 0) {
      covariate <- matrix(1, nrow = N, ncol = 1)
      colnames(covariate) <- "Intercept"
    }
  }
  covariates <- list(covariate)
  p <- ncol(covariate)
  num_beta <- p * (L - 1)
  num_gamma <- 0
  init.par <- rep(0, num_beta + num_gamma)

  if(vis){
    cat(paste0("Starting the third step for ", type, " ...\n\n"))
  }

  npar <- length(init.par)
  lb <- rep(lower, npar)
  ub <- rep(upper, npar)

  int_width <- ceiling(log10(N * I * L))
  total_width <- int_width + 5
  fmt_string_maxchg <- sprintf("%%%d.%df", total_width, 5)

  int_width <- ceiling(log10(N * I * L)) + 1L
  total_width <- int_width + 3
  fmt_string_BIC <- sprintf("%%%d.%df", total_width, 2)

  Log.Lik.history <- c(0)
  make_loglik_with_print <- function(vis, ref.class) {
    iter <- 0
    function(par, CEP, P.Z.Xns, Zs, covariates, covariates.timeCross) {
      iter <<- iter + 1
      Log.Lik <- get.Log.Lik.LTA.optim(par, CEP, P.Z.Xns, Zs, covariates, FALSE, ref.class)

      Log.Lik.history <<- c(Log.Lik.history, -Log.Lik)
      maxchg <- abs(Log.Lik.history[iter+1] - Log.Lik.history[iter])
      BIC <- -2 * Log.Lik.history[iter+1] + npar * log(N)
      if(vis){
        cat('\r  Iter =', sprintf("%4d", iter), '  \u0394Log.Lik =', sprintf(fmt_string_maxchg, maxchg),
            '  BIC =', sprintf(fmt_string_BIC, BIC))
      }
      return(Log.Lik)
    }
  }
  get.Log.Lik.LTA.print <- make_loglik_with_print(vis, ref.class)

  nlopt_res <- nloptr::nloptr(
    x0 = init.par,
    eval_f = function(x) get.Log.Lik.LTA.print(x, CEP, P.Z.Xns, Zs, covariates, FALSE),
    lb = lb,
    ub = ub,
    opts = list(
      algorithm = "NLOPT_LN_BOBYQA",
      initial_step = rep(0.5, length(init.par)),
      xtol_rel = 1e-8,
      ftol_rel = 1e-8,
      maxeval = maxiter,
      print_level = 0
    )
  )
  if(vis){
    cat("\n\n")
  }

  refined_init <- nlopt_res$solution

  if(method.SE == "Obs"){
    if(vis){
      cat("  Calculating Observed Information Matrix for Standard Errors ...\n")
    }
    h <- tryCatch({
      base_eps <- .Machine$double.eps^(1/3)
      adaptive_eps <- base_eps * pmax(1, abs(refined_init))
      numDeriv::hessian(
        func = function(x) {
          get.Log.Lik.LTA.optim(
            x,
            CEP,
            P.Z.Xns,
            Zs,
            covariates,
            FALSE,
            ref.class
          )
        },
        x = refined_init,
        method = "Richardson",
        method.args = list(
          eps = adaptive_eps,
          r   = 6,
          v   = 2,
          zero.tol = .Machine$double.eps * 100
        )
      )

    }, error = function(e1) {
      tryCatch({
        warning("Hessian (r = 6) failed: ", conditionMessage(e1),
                ". Retrying with r = 4.")
        numDeriv::hessian(
          func = function(x) {
            get.Log.Lik.LTA.optim(
              x,
              CEP,
              P.Z.Xns,
              Zs,
              covariates,
              FALSE,
              ref.class
            )
          },
          x = refined_init,
          method = "Richardson",
          method.args = list(
            eps = adaptive_eps,
            r   = 4,
            v   = 2,
            zero.tol = .Machine$double.eps * 100
          )
        )

      }, error = function(e2) {
        tryCatch({
          warning("Hessian (r = 4) failed: ", conditionMessage(e2),
                  ". Retrying with r = 2.")
          numDeriv::hessian(
            func = function(x) {
              get.Log.Lik.LTA.optim(
                x,
                CEP,
                P.Z.Xns,
                Zs,
                covariates,
                FALSE,
                ref.class
              )
            },
            x = refined_init,
            method = "Richardson",
            method.args = list(
              eps = adaptive_eps,
              r   = 2,
              v   = 2,
              zero.tol = .Machine$double.eps * 100
            )
          )

        }, error = function(e3) {
          warning("All Hessian attempts failed: ", conditionMessage(e3))
          message(
            "Suggestions:\n",
            "1. Check local identifiability / label switching\n",
            "2. Inspect parameter scaling (very large |theta|)\n",
            "3. Try profile likelihood or reparameterization\n",
            "4. Consider observed information via EM if applicable\n",
            "Current parameters:\n",
            paste(round(refined_init, 4), collapse = ", ")
          )
          NULL
        })
      })
    })

    if (!is.null(h)) {
      near_h <- Matrix::nearPD(h, corr = FALSE,
                               eig.tol = 1e-10,
                               conv.tol = 1e-8)$mat
      vcov <- tryCatch({
        solve(near_h)
      }, error = function(e) {
        MASS::ginv(as.matrix(near_h))
      })
      diag_vcov <- diag(vcov)
      if (any(diag_vcov < 0)) {
        diag_vcov[diag_vcov < 0] <- NA
      }

      se.vec <- sqrt(diag_vcov)
    }else{
      se.vec <- NULL
    }
  }else{
    if(vis){
      cat("  Bootstrapping for Standard Errors ...\n")
    }
    par.Bootstrap <- matrix(0, n.Bootstrap, npar)
    for(bs in 1:n.Bootstrap){
      covariates.cur <- vector("list", 1)
      samples.cur <- sample(1:N, N, replace = TRUE)
      covariates.cur[[1]] <- covariates[[1]][samples.cur, , drop=FALSE]

      P.Z.Xns.cur <- Zs.cur <- P.Zs.cur <- vector("list", 1)
      if(type == "LCA"){
        P.Z.Xns.cur[[1]] <- get.P.Z.Xn.LCA(response=response[samples.cur, ], par=params$par, P.Z=params$P.Z)
      } else {
        P.Z.Xns.cur[[1]] <- get.P.Z.Xn.LPA(response=response[samples.cur, ], means=params$means, covs=params$covs, P.Z=params$P.Z)
      }

      P.Zs.cur[[1]] <- colSums(P.Z.Xns.cur[[1]]) / sum(P.Z.Xns.cur[[1]])
      Zs.cur[[1]] <- apply(P.Z.Xns.cur[[1]], 1, which.max)

      if(CEP.error){
        CEP.cur <- get.CEP(P.Z.Xns.cur, time.cross=FALSE)
      } else {
        CEP.cur <- replicate(1, diag(L), simplify=FALSE)
      }

      init.par <- rnorm(npar, mean=0, sd=abs(refined_init) * 0.1)
      Log.Lik.history.cur <- c(0)
      make_loglik_with_print_Bootstrap <- function(vis, ref.class, bs, n.Bootstrap) {
        iter <- 0
        function(par, CEP.cur, P.Z.Xns.cur, Zs.cur, covariates.cur, covariates.timeCross, bs, n.Bootstrap) {
          iter <<- iter + 1
          Log.Lik <- get.Log.Lik.LTA.optim(par, CEP.cur, P.Z.Xns.cur, Zs.cur, covariates.cur, FALSE, ref.class)

          Log.Lik.history.cur <<- c(Log.Lik.history.cur, -Log.Lik)
          maxchg <- abs(Log.Lik.history.cur[iter+1] - Log.Lik.history.cur[iter])
          BIC <- -2 * Log.Lik.history.cur[iter+1] + npar * log(N)
          if(vis){
            cat('\r  Bootstrap =', sprintf("%4d/%4d", bs, n.Bootstrap), ' | Iter =', sprintf("%4d", iter),
                '  \u0394Log.Lik =', sprintf(fmt_string_maxchg, maxchg),
                '  BIC =', sprintf(fmt_string_BIC, BIC))
          }
          return(Log.Lik)
        }
      }
      get.Log.Lik.LTA.print <- make_loglik_with_print_Bootstrap(vis, ref.class, bs, n.Bootstrap)
      nlopt_res <- nloptr::nloptr(
        x0 = init.par,
        eval_f = function(x) get.Log.Lik.LTA.print(x, CEP.cur, P.Z.Xns.cur, Zs.cur, covariates.cur, FALSE, bs, n.Bootstrap),
        lb = lb,
        ub = ub,
        opts = list(
          algorithm = "NLOPT_LN_BOBYQA",
          initial_step = rep(0.5, length(init.par)),
          xtol_rel = 1e-8,
          ftol_rel = 1e-8,
          maxeval = maxiter,
          print_level = 0
        )
      )
      par.Bootstrap[bs, ] <- nlopt_res$solution
    }

    se.vec <- apply(par.Bootstrap, 2, sd, na.rm = TRUE)
    if(vis){
      cat("\n\n")
    }
  }

  params.LTA.obj <- LTA.vector.to.parameters(refined_init, covariates, L, ref.class)

  if(!is.null(se.vec)){
    se.obj <- LTA.vector.to.parameters(se.vec, covariates, L, ref.class)
    Z.sta.vec <- refined_init/se.vec
    p.value.tail1 <- pnorm(-abs(Z.sta.vec))
    p.value.tail2 <- pnorm(-abs(Z.sta.vec)) * 2
    Z.sta.obj <- LTA.vector.to.parameters(Z.sta.vec, covariates, L, ref.class)
    p.value.tail1.obj <- LTA.vector.to.parameters(p.value.tail1, covariates, L, ref.class)
    p.value.tail2.obj <- LTA.vector.to.parameters(p.value.tail1, covariates, L, ref.class)
  }

  beta <- params.LTA.obj$beta
  beta.se <- se.obj$beta
  beta.Z.sta <- Z.sta.obj$beta
  beta.p.value.tail1 <- p.value.tail1.obj$beta
  beta.p.value.tail2 <- p.value.tail2.obj$beta

  covariates.ncol <- unlist(lapply(covariates, ncol))
  npar <- get.npar.LTA(covariates.ncol, L, FALSE)

  Log.Lik = -get.Log.Lik.LTA.optim(refined_init, CEP, P.Z.Xns, Zs, covariates, FALSE, ref.class)
  AIC <- -2 * Log.Lik + 2 * npar
  BIC <- -2 * Log.Lik + npar * log(N)

  # Prepare result object
  res <- list(
    beta = beta,
    beta.se = beta.se,
    beta.Z.sta = beta.Z.sta,
    beta.p.value.tail1 = beta.p.value.tail1,
    beta.p.value.tail2 = beta.p.value.tail2,
    npar = npar,
    Log.Lik = Log.Lik,
    AIC = AIC,
    BIC = BIC,
    P.Z.Xn = P.Z.Xns[[1]],
    P.Z = P.Zs[[1]],
    Z = Zs[[1]],
    Log.Lik.history = Log.Lik.history[1:nlopt_res$iterations],
    iterations = nlopt_res$iterations,
    coveraged = nlopt_res$iterations < maxiter,
    params = params,
    call = call,
    arguments = list(
      response = response, L = L, type = type,
      covariate = covariate,
      CEP.error = CEP.error,
      par.ini = par.ini,
      params = params,
      constraint = constraint,
      method = method, tol = tol,
      maxiter = maxiter, is.sort = is.sort,
      nrep = nrep, starts = starts, maxiter.wa = maxiter.wa,
      vis = vis,
      control.EM = control.EM,
      control.Mplus = control.Mplus,
      control.NNE = control.NNE,
      ref.class = ref.class
    )
  )

  class(res) <- "LCPA"
  return(res)
}
