#' Latent Class/Profile Analysis with Covariates
#'
#' Implements the three-step estimation method (Vermunt, 2010; Liang et al., 2023) for latent class/profile analysis
#' with covariates, treating latent class membership as an observed variable with measurement error.
#' This is mathematically equivalent to a latent transition analysis (\code{\link[LCPA]{LTA}}) with \code{times=1}.
#'
#' @param response A matrix or data frame of observed responses.
#'                  Rows of the matrix represent individuals/participants/observations (\eqn{N}), columns of the
#'                  matrix represent observed indicators/items/variables (\eqn{I}).
#'                  For \code{type.model = "LCA"}: indicators must be binary or categorical (coded as integers starting from 0).
#'                  For \code{type.model = "LPA"}: indicators must be continuous (numeric), and the response matrix must be
#'                  standardized using \code{\link[base]{scale}} or \code{\link[LCPA]{normalize}} prior to input.
#' @param L Integer scalar. Number of latent classes/profiles. Must satisfy \eqn{L \geq 2}.
#' @param ref.class Integer \eqn{L \geq ref.class \geq 1}. Specifies which latent class to use as the reference category.
#'                  Default is \code{L} (last class). Coefficients for the reference class are fixed to zero.
#'                  When \code{is.sort=TRUE}, classes are first ordered by decreasing \code{P.Z} (class 1 has highest probability),
#'                  then \code{ref.class} refers to the position in this sorted order.
#' @param type.model Character string. Specifies the type of latent variable model for Step 1:
#'             \itemize{
#'               \item \code{"LCA"} — Latent Class Analysis for categorical indicators.
#'               \item \code{"LPA"} — Latent Profile Analysis for continuous indicators.
#'             }
#'             See \code{\link[LCPA]{LCA}} and \code{\link[LCPA]{LPA}} for details.
#' @param covariates Optional. A matrix or data frame of covariates for modeling latent class membership.
#'                  Must include an intercept column (all 1s) as the first column.
#'                  All non-intercept covariates must be standardized before analysis, for example
#'                  with \code{scale()} or \code{\link[LCPA]{normalize}}. This places the covariate regression
#'                  coefficients on comparable scales and prevents their magnitudes from being driven by the
#'                  covariates' original measurement units. Do not standardize the intercept column; construct
#'                  interaction terms from the standardized covariates.
#'                  If \code{NULL} (default), only intercept terms are used (i.e., no covariates).
#'                  Dimension is \eqn{N\times(U+1)}, where \eqn{U} is the number
#'                  of observed covariates and the additional column is the intercept.
#' @param CEP.error Logical. If \code{TRUE} (recommended), incorporates classification uncertainty via
#'                  estimated Classification Error Probability (\code{\link[LCPA]{get.CEP}}) matrices from Step 1. If \code{FALSE},
#'                  uses identity CEP matrices (equivalent to naive modal assignment; introduces bias).
#' @param par.ini Specification for parameter initialization. Options include:
#'   \itemize{
#'     \item \code{"random"}: Completely random initialization (default).
#'     \item \code{"kmeans"}: Initializes parameters via K-means clustering on observed data (McLachlan & Peel, 2000).
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
#'   For methods that expose this initialization interface, \code{par.ini} is used only
#'   to construct the \code{starts} warm-up initializations for the internally fitted
#'   Step 1 model. It does not initialize the \code{nrep} refinement runs, which continue
#'   directly from the selected warm-up states. With \code{par.ini = "kmeans"}, each
#'   outer warm-up start performs exactly one K-means run, so \code{starts} is the number
#'   of separately initialized K-means outputs passed to warm-up training.
#'   Backends such as Mplus and Rmixmod may instead use their native initialization
#'   mechanisms and are not required to implement K-means initialization.
#'   If \code{"kmeans"} is requested for a method that does not support it, \code{par.ini}
#'   is automatically changed to \code{"random"}.
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
#'                      \item{\code{category.levels}}{Fixed ordered response categories for each indicator.}
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
#'                  The complete named set is \code{"UE"}, \code{"UV"}, \code{"E0"},
#'                  \code{"V0"}, \code{"EE"}, \code{"VV"}, \code{"VE"}, and \code{"EV"};
#'                  custom equality lists follow \code{\link[LCPA]{LPA}()} and its backend restrictions.
#' @param method.model Character. Estimation algorithm for Step 1 models:
#'               \itemize{
#'                 \item \code{"EM"} — Expectation-Maximization (default; robust and widely used).
#'                 \item \code{"Mplus"} — Interfaces with Mplus software (requires external installation).
#'                 \item \code{"NNE"} — Neural Network Estimator (experimental), using feed-forward layers,
#'                                      optional transformer attention, gradient optimization, and simulated annealing.
#'                 \item \code{"Rmixmod"}: Package-staged SEM or a native EM/CEM/SEM strategy
#'                   through the optional \code{Rmixmod} backend.
#'                 \item \code{"RMixtComp"}: Stochastic EM (SEM) through the optional \code{RMixtComp} backend.
#'                 \item \code{"flexmix"}: Stochastic EM (SEM) through the optional \code{flexmix} backend.
#'               }
#' @param tol Convergence tolerance for log-likelihood difference (default: 1e-4).
#' @param method.regression Step 3 derivative method. \code{"Analytic"}
#'   uses the compiled exact score; \code{"Numeric"} uses Richardson numerical
#'   differentiation (default: \code{"Analytic"}).
#' @param lower Lower bound for regression coefficients (default: -10).
#' @param upper Upper bound for regression coefficients (default: 10).
#' @param method.SE Character. Method for estimating standard errors of parameter estimates:
#'               \itemize{
#'                 \item \code{"Analytic"}: Analytic observed information based on Louis' identity.
#'                       The analytic information is used as the sandwich bread; the empirical
#'                       individual score and the influence of re-estimating CEP form the meat.
#'                 \item \code{"Numeric"} — Approximates the observed information matrix via numerical differentiation (Richardson's method).
#'                       The numerical information is used with the same empirical sandwich meat as \code{"Analytic"}.
#'                 \item \code{"Bootstrap"} — Uses nonparametric bootstrap resampling to estimate empirical sampling variability.
#'                       More robust to model misspecification and small-sample bias. Computationally intensive but recommended when
#'                       asymptotic assumptions are questionable.
#'               }
#'               Default is \code{"Bootstrap"}.
#' @param nrep.bootstrap Integer. Number of bootstrap replicates used when \code{method.SE = "Bootstrap"}.
#'                    Default is 100. McLachlan & Peel (2000) suggest that 50–100 replicates often provide adequate accuracy
#'                    for practical purposes, though more (e.g., 500–1000) may be preferred for publication-quality inference.
#'                    Each replicate samples \eqn{N} individuals with replacement, keeps the Step 1 measurement parameters fixed,
#'                    recomputes posterior assignments and CEP, and re-estimates Step 3. Only successful
#'                    optimization replicates contribute to the empirical covariance matrix.
#' @param maxiter Maximum number of iterations for optimizing the regression coefficients. Default: 5000.
#' @param starts Positive integer. Number of Step 1 warm-up analyses to run (default: 100).
#'   Each analysis is initialized by the selected method and trained for at most
#'   \code{maxiter.warmup} iterations, producing exactly \code{starts} warm-up solutions.
#'   For methods supporting \code{par.ini = "kmeans"}, each outer warm-up start performs
#'   exactly one K-means run, so \code{starts} is also the number of K-means outputs.
#' @param maxiter.warmup Positive integer. Maximum number of training iterations for
#'   each of the \code{starts} Step 1 warm-up analyses (default: 20). This limit applies
#'   only to warm-up and does not limit the subsequent refinement phase.
#' @param nrep Positive integer not exceeding \code{starts}. Number of Step 1 refinement
#'   analyses (default: 20). The \code{nrep} warm-up solutions with the largest
#'   log-likelihoods are continued from their saved states until the full-training
#'   stopping rule is met; no new initialization occurs in this phase. The refined
#'   solution with the largest log-likelihood is used as the final Step 1 result.
#'   These three staged-training arguments are not used by \code{method = "RMixtComp"}
#'   or by \code{method = "Rmixmod"} with \code{control.Rmixmod$path = "Rmixmod"};
#'   those paths use their documented native controls instead.
#' @param vis Logical. If \code{TRUE}, displays the three-step stage headings and
#'            indented process information for Step 1 model estimation, Step 2
#'            posterior classification and CEP-matrix preparation, and Step 3
#'            optimization and standard-error estimation (default: \code{TRUE}).
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
#'     \item{\code{patience.early}}{Patience parameter for early stopping (default: 100).}
#'     \item{\code{maxcycle}}{Maximum cycles for optimization (default: 20).}
#'     \item{\code{lr}}{Learning rate, controlling the step size of neural network parameter updates (default: 0.025).}
#'     \item{\code{scheduler.patience}}{Patience for learning rate decay (if the loss function does not improve for more than `patience` consecutive epochs, the learning rate will be reduced) (default: 10).}
#'     \item{\code{scheduler.factor}}{Learning rate decay factor; the new learning rate equals the original learning rate multiplied by \code{scheduler.factor} (default: 0.80).}
#'     \item{\code{plot.interval}}{Interval (in epochs) for plotting training diagnostics (default: 200).}
#'     \item{\code{device}}{Specifies the hardware device; can be \code{"CPU"} (default) or \code{"GPU"}. If the GPU is not available, it automatically falls back to CPU.}
#'   }
#' @param control.Mplus List of control parameters for Mplus estimation:
#'   \describe{
#'     \item{\code{maxiter}}{Maximum iterations for Mplus optimization (default: 2000).}
#'     \item{\code{tol}}{Convergence tolerance for log-likelihood difference (default: 1e-4).}
#'     \item{\code{files.path}}{A character string specifying the directory under which Mplus writes
#'       intermediate files, including model input, data, output, and saved posterior probabilities.
#'       When \code{method.model = "Mplus"}, \code{NULL} (the default) is invalid and produces an error.
#'       A non-empty path is created recursively when necessary and must be writable. Within it, the
#'       Step 1 function creates a unique timestamped \code{"Mplus_LCA_YYYY-MM-DD_HH-MM-SS"} or
#'       \code{"Mplus_LPA_YYYY-MM-DD_HH-MM-SS"} subdirectory to isolate all files from the current run.
#'       If \code{files.path = ""}, the timestamped subdirectory is created directly under R's current
#'       working directory, \code{\link[base]{getwd}()}.}
#'     \item{\code{files.clean}}{Logical. If \code{TRUE} (default), all intermediate files and the temporary working directory
#'       created for the Step 1 run are deleted on successful completion or error exit via \code{on.exit()}.
#'       If \code{FALSE}, the complete timestamped working directory is retained under \code{files.path},
#'       or under \code{\link[base]{getwd}()} when \code{files.path = ""}, for inspection and debugging.}
#'   }
#'
#' @param control.Rmixmod Optional Rmixmod controls passed to \code{\link[LCPA]{LCA}()} or \code{\link[LCPA]{LPA}()}.
#'   Supports \code{path}, \code{algorithm}, \code{nrep}, \code{method.init}, \code{starts},
#'   \code{maxiter.init}, \code{maxiter}, \code{tol.init}, \code{tol}, \code{par.ini},
#'   \code{labels.ini}, and \code{strategy}; see \code{\link[LCPA]{LCA}()} and \code{\link[LCPA]{LPA}()}
#'   for their complete semantics. The default \code{path = "LCPA"} uses SEM, whereas native
#'   \code{path = "Rmixmod"} supports EM, CEM, SEM, and ordered combinations. LPA Step 1
#'   models support \code{"V0"}, \code{"EE"}, and \code{"VV"}.
#' @param control.flexmix Optional flexmix SEM controls passed to \code{\link[LCPA]{LCA}()} or \code{\link[LCPA]{LPA}()}.
#'   Supports \code{maxiter}, \code{minprior}, and \code{tol};
#'   \code{classify = "SEM"} is fixed by LCPA.
#' @param control.RMixtComp Optional RMixtComp SEM controls passed to \code{\link[LCPA]{LCA}()} or \code{\link[LCPA]{LPA}()}.
#'   See those functions for the full native control list. RMixtComp is used only for SEM,
#'   without LCPA's external warm-up or promoted-replication stages; for LPA, it requires
#'   \code{constraint = "V0"}.
#'
#' @section Random-number reproducibility:
#' Except for \code{method = "NNE"}, which intentionally uses its fixed backend seed,
#' all stochastic stages and supported backends are driven from R's current random-number
#' generator. The user only needs to call \code{set.seed()} immediately before \code{\link[LCPA]{LCPA}()};
#' no flexmix, Rmixmod, RMixtComp, or Mplus-specific seed setting is required.
#'
#' @return An object of class \code{LCPA}, a named list containing:
#' \describe{
#'   \item{\code{beta}}{Matrix of size \eqn{(U+1)\times L}. Coefficients for class membership multinomial logit model.
#'     All non-reference columns are free; column \code{ref.class} is constrained to zero.}
#'   \item{\code{beta.se}}{Standard errors for \code{beta}. The fixed reference-class column is \code{NA};
#'     free coefficients may also be \code{NA} when the information matrix is unavailable or ill-conditioned.}
#'   \item{\code{beta.Z.sta}}{Z-statistics for testing null hypothesis that each beta coefficient equals zero.
#'     Computed as \code{beta / beta.se}. Same structure as \code{beta}.}
#'   \item{\code{beta.p.value.tail1}}{One-tailed p-values based on standard normal distribution: \eqn{P(Z < -|z|)}.
#'     Useful for directional hypotheses. Same structure as \code{beta}.}
#'   \item{\code{beta.p.value.tail2}}{Two-tailed p-values: \eqn{2 \times P(Z < -|z|)}.
#'     Standard test for non-zero effect. Same structure as \code{beta}.}
#'   \item{\code{vcov}}{Variance-covariance matrix of the free Step 3 coefficients. For bootstrap,
#'     this is the empirical covariance matrix across successful replicates; for
#'     \code{"Numeric"} and \code{"Analytic"}, it is the sandwich covariance matrix.}
#'   \item{\code{information}}{Observed information matrix for
#'     \code{method.SE = "Numeric"} or \code{"Analytic"}; \code{NULL} for
#'     bootstrap.}
#'   \item{\code{SE.diagnostics}}{Conditioning diagnostics for observed-information methods,
#'     or the bootstrap replication count.}
#'   \item{\code{bound.diagnostics}}{Step 3 coefficient indices within \code{1e-4} of an optimization bound.}
#'   \item{\code{P.Z.Xn}}{Matrix of size \eqn{N \times L} of posterior class probabilities
#'     \eqn{P(Z_n=l \mid \mathbf{X}_n)} for each participant \eqn{n} and class \eqn{l}.}
#'   \item{\code{P.Z}}{Vector of length \eqn{L} containing prior class proportions
#'     \eqn{P(Z = l)} estimated at Step 1.}
#'   \item{\code{Z}}{Vector of length \eqn{N} containing modal class assignments
#'     (MAP classifications) \eqn{\widehat{Z}_n} for each participant.}
#'   \item{\code{npar}}{Number of free parameters in the model (depends on covariates).}
#'   \item{\code{Log.Lik}}{Observed-data log-likelihood value at convergence.}
#'   \item{\code{Log.Lik.history}}{Vector tracking log-likelihood at each iteration.}
#'   \item{\code{AIC}}{Akaike Information Criterion value.}
#'   \item{\code{BIC}}{Bayesian Information Criterion value.}
#'   \item{\code{iterations}}{Integer. Number of optimization iterations in Step 3.}
#'   \item{\code{converged}}{Logical indicator based on a successful NLopt termination code.}
#'   \item{\code{params}}{List. Step 1 model parameters (output from \code{\link[LCPA]{LCA}()} or \code{\link[LCPA]{LPA}()}).}
#'   \item{\code{call}}{The matched function call.}
#'   \item{\code{arguments}}{List of all input arguments passed to the function (useful for reproducibility).}
#' }
#'
#' @section Methodology Overview:
#' The three-step procedure follows the same principles as LTA but for a single time point:
#'
#' Step 1 — Unconditional Latent Class/Profile Model:
#' Fit an unconditional LCA or LPA model (ignoring covariates). Obtain posterior class membership probabilities
#' \eqn{P(Z_n=l \mid \mathbf{X}_n)} for each participant \eqn{n} and class
#' \eqn{l} using Bayes' theorem.
#'
#' Step 2 — Classification Error Probabilities (equal to \code{\link[LCPA]{get.CEP}}):
#' Compute the \eqn{L\times L} matrix with
#' \eqn{\mathrm{CEP}(l,k)=P(\widehat{Z}_n=k\mid Z_n=l)}. Its modal assignment,
#' posterior-weight estimator, and matrix orientation are defined in
#' \code{\link[LCPA]{get.CEP}()}.
#'
#' Step 3 — Class Membership Model with Measurement Error Correction:
#' Estimate the multinomial logit model for class membership:
#' \deqn{
#' P(Z_n=l\mid\boldsymbol{\zeta}_n) =
#' \frac{\exp(\boldsymbol{\beta}_l^\top\boldsymbol{\zeta}_n)}
#' {\sum_{h=1}^L\exp(\boldsymbol{\beta}_h^\top\boldsymbol{\zeta}_n)}
#' }
#' where
#' \eqn{\boldsymbol{\zeta}_n=(1,\zeta_{n1},\ldots,\zeta_{nU})^\top} is the
#' covariate vector for participant \eqn{n}, with
#' \eqn{u=1,\ldots,U}, and
#' \eqn{\boldsymbol{\beta}_l=(\beta_{l0},\beta_{l1},\ldots,\beta_{lU})^\top}
#' contains the corresponding intercept and slopes. The class selected by \code{ref.class} is the
#' reference category and its coefficient vector is fixed to zero.
#'
#' The observed-data likelihood integrates over latent classes:
#' \deqn{
#' \log \mathcal{L}(\boldsymbol{\beta}) =
#' \sum_{n=1}^N \log \left[
#'   \sum_{l=1}^L \mathrm{CEP}(l,\widehat{Z}_n)
#'   P(Z_n=l\mid\boldsymbol{\zeta}_n)
#' \right]
#' }
#' Parameters \eqn{\boldsymbol{\beta}} are estimated via maximum likelihood using the L-BFGS algorithm.
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
#'                                 column. Dimensions must be \eqn{N\times(U+1)},
#'                                 where \eqn{U} is the number of covariates.
#'   \item Covariate standardization: All non-intercept covariates must be standardized before
#'         analysis (e.g., with \code{scale()} or \code{\link[LCPA]{normalize}}). Otherwise, the estimated
#'         covariate regression coefficients remain tied to incompatible original measurement units, creating
#'         a scale problem that prevents meaningful comparison of coefficient magnitudes. The intercept must
#'         remain equal to 1, and interaction terms should be formed after standardizing their component covariates.
#'   \item Optimization and standard errors:
#'         \itemize{
#'           \item Step 3 uses the L-BFGS algorithm (\code{nloptr::nloptr}) for stable optimization with box constraints.
#'           \item For \code{method.SE = "Analytic"} or \code{"Numeric"}: the analytic or numerical observed
#'                 information supplies the sandwich bread. The meat uses individual empirical scores and,
#'                 when \code{CEP.error = TRUE}, the CEP ratio-estimator influence function. Thus both target
#'                 the same fixed-Step-1 nonparametric sampling scheme as bootstrap.
#'           \item For \code{method.SE = "Bootstrap"}: each replicate keeps Step 1 fixed, recomputes Step 2
#'                 assignments and CEP, and re-estimates Step 3. Only successful replicates are retained.
#'         }
#'   \item Computational notes:
#'         \itemize{
#'           \item Step 1 complexity increases with \eqn{L} and \eqn{I}.
#'           \item Bootstrap is computationally intensive because Steps 2 and 3 are repeated for every replicate.
#'         }
#'   \item Bootstrap reproducibility: Always set a seed (e.g., \code{set.seed(123)}) before
#'                                             calling \code{\link[LCPA]{LCPA}()} when using \code{method.SE = "Bootstrap"}.
#'                                             Monitor convergence in bootstrap runs via progress messages.
#' }
#'
#' @references
#' Liang, Q., de la Torre, J., & Law, N. (2023). Latent transition cognitive
#' diagnosis model with covariates: A three-step approach. *Journal of
#' Educational and Behavioral Statistics, 48*(6), 690--718.
#' \doi{10.3102/10769986231163320}
#'
#' Louis, T. A. (1982). Finding the observed information matrix when using the EM algorithm.
#' *Journal of the Royal Statistical Society: Series B (Methodological),
#' 44*(2), 226--233.
#' \doi{10.1111/j.2517-6161.1982.tb01203.x}
#'
#' Vermunt, J. K. (2010). Latent class modeling with covariates: Two improved
#' three-step approaches. *Political Analysis, 18*(4), 450--469.
#' \doi{10.1093/pan/mpq025}
#'
#' @examples
#'
#' \donttest{
#' ## long time
#' library(LCPA)
#'
#' set.seed(1245)
#' N <- 2000  # Sample size
#' L <- 3    # Number of latent classes
#' I <- 6    # Number of indicators
#'
#' # Create standardized covariates (intercept + 2 covariates + 1 interaction)
#'  Intercept = rep(1, N)
#'  X1 <- as.numeric(scale(rnorm(N)))
#'  X2 <- as.numeric(scale(rbinom(N, 1, 0.5)))
#'  X1.X2 <- X1 * X2
#' covariates <- cbind(Intercept, X1, X2, X1.X2)
#'
#' # Simulate data for LPA
#' sim_data <- sim.LTA(
#'   N = N, I = I, L = L, times = 1, type = "LPA",
#'   ref.class = 2,
#'   covariates = list(covariates), is.sort=TRUE,
#'   beta = matrix(c(
#'    -0.2, 0.0, -0.1,  ## fix reference class to class 2
#'     0.2, 0.0, -0.3,
#'     0.8, 0.0, -0.6,
#'    -0.3, 0.0,  0.3
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
#' fit <- LCPA(
#'   response = response,
#'   L = L, ref.class = 2,
#'   type.model = "LPA", is.sort=TRUE,
#'   covariates = covariates,
#'   method.SE = "Analytic",
#'   CEP.error = TRUE,
#'   method = "EM",
#'   vis = TRUE
#' )
#' print(fit)
#' }
#'
#' # Optional external Step 1 backends retain the same Step 2 and Step 3 workflow.
#' \dontrun{
#' # need Mplus
#' fit.mplus <- LCPA(
#'   response, L = L, ref.class = 2, type.model = "LPA", covariates = covariates,
#'   method.SE = "Analytic", method.model = "Mplus",
#'   control.Mplus = list(files.path = "", files.clean = TRUE), vis = TRUE
#' )
#'
#' # need Rmixmod
#' fit.rmixmod <- LCPA(
#'   response, L = L, ref.class = 2, type.model = "LPA", covariates = covariates,
#'   method = "Rmixmod", constraint = "VV", nrep = 3,
#'   control.Rmixmod = list(maxiter = 200), vis = TRUE
#' )
#' }
#'
#' @importFrom nloptr nloptr
#' @importFrom Matrix nearPD
#' @importFrom numDeriv hessian
#' @importFrom stats pnorm complete.cases
#' @importFrom MASS ginv
#'
#' @noRd
ML.XZ.LCPA <- function(response, L = 2,
                 ref.class = L, type.model = "LCA",
                 covariates = NULL,
                 CEP.error = TRUE,
                 par.ini = "random",
                 params = NULL, is.sort = TRUE,
                 constraint = "VV",
                 method.model = "EM", tol = 1e-4,
                 method.regression = "Analytic",
                 lower=-10, upper=10,
                 method.SE = "Bootstrap", nrep.bootstrap=100,
                 maxiter = 5000, starts = 100,
                 maxiter.warmup = 20, nrep = 20,
                 vis = TRUE,
                 control.EM = NULL,
                 control.Mplus = NULL,
                 control.NNE = NULL,
                 control.flexmix = NULL,
                 control.Rmixmod = NULL,
                 control.RMixtComp = NULL,
                 .progress.path = "X -> Z",
                 .progress.step3 = TRUE) {

  call <- match.call()
  method.model <- match.arg(method.model, c("EM", "NNE", "Mplus", "flexmix", "Rmixmod", "RMixtComp"))
  if(is.null(params)) par.ini <- .normalize.par.ini(par.ini, method.model)
  type.model <- match.arg(type.model, c("LCA", "LPA"))
  method.regression <- match.arg(method.regression, c("Analytic", "Numeric"))
  method.SE <- match.arg(method.SE, c("Analytic", "Numeric", "Bootstrap"))
  method.SE.internal <- c(Analytic = "Louis", Numeric = "Obs", Bootstrap = "Bootstrap")[[method.SE]]
  gradient.cur <- tolower(method.regression)
  vis.step3 <- isTRUE(vis) && isTRUE(.progress.step3)

  if (ref.class < 1 || ref.class > L) {
    stop("ref.class must be between 1 and L")
  }
  if (method.SE == "Bootstrap" && nrep.bootstrap < 2L) {
    stop("nrep.bootstrap must be at least 2 when method.SE = 'Bootstrap'")
  }

  # Default control parameters
  default_control.EM <- list(maxiter = 2000, tol = 1e-4)
  default_control.Mplus <- list(maxiter = 2000, tol = 1e-4, files.path = NULL, files.clean = TRUE)
  default_control.Rmixmod <- .default.Rmixmod.control(maxiter = 1000L)
  default_control.RMixtComp <- .default.RMixtComp.control()
  default_control.flexmix <- .default.flexmix.control()
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
    patience.early = 100,
    maxcycle = 20,
    lr = 0.025,
    scheduler.patience = 10,
    scheduler.factor = 0.80,
    plot.interval = 200,
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
  control.flexmix <- merge_and_clean_control(control.flexmix, default_control.flexmix)
  control.Rmixmod <- merge_and_clean_control(control.Rmixmod, default_control.Rmixmod)
  control.RMixtComp <- merge_and_clean_control(control.RMixtComp, default_control.RMixtComp)

  # Convert to matrix if needed
  response <- as.matrix(response)
  N <- nrow(response)
  I <- ncol(response)

  .three.step.progress.start(
    1L, "LCPA", type.model = type.model, vis = vis
  )

  # Step 1: Unconditional LCA/LPA
  params.supplied <- !is.null(params)
  if (is.null(params)) {
    LCPA.obj <- .with.estimation.output.prefix({
      if (type.model == "LCA") {
        LCA(response, L = L, par.ini = par.ini,
            method = method.model, is.sort = is.sort, nrep = nrep,
            starts = starts, maxiter.warmup = maxiter.warmup,
            vis = vis,
            control.EM = control.EM,
            control.Mplus = control.Mplus,
            control.NNE = control.NNE,
            control.flexmix = control.flexmix,
            control.Rmixmod = control.Rmixmod,
            control.RMixtComp = control.RMixtComp)
      } else {
        LPA(response, L = L, par.ini = par.ini,
            constraint = constraint,
            method = method.model, is.sort = is.sort, nrep = nrep,
            starts = starts, maxiter.warmup = maxiter.warmup,
            vis = vis,
            control.EM = control.EM,
            control.Mplus = control.Mplus,
            control.NNE = control.NNE,
            control.flexmix = control.flexmix,
            control.Rmixmod = control.Rmixmod,
            control.RMixtComp = control.RMixtComp)
      }
    })
    params <- LCPA.obj$params
  } else if(vis) {
    cat("  Supplied Step 1 parameters were used.\n")
  }

  if (type.model == "LCA") {
    if (is.null(params$category.levels)) {
      stop("LCA Step 1 parameters must include category.levels")
    }
  }
  if(params.supplied && is.sort){
    position <- order(params$P.Z, decreasing = TRUE)
    params <- .three.step.reorder.params(params, position, type.model)
  }

  .three.step.progress.start(
    2L, "LCPA", path = .progress.path, vis = vis
  )
  P.Z.Xns <- list(.three.step.posterior(response, params, type.model))
  P.Zs <- list(colSums(P.Z.Xns[[1]]) / sum(P.Z.Xns[[1]]))
  Zs <- list(max.col(P.Z.Xns[[1]], ties.method = "first"))

  if (CEP.error) {
    CEP <- get.CEP(P.Z.Xns, CEP.time.cross =FALSE)
  } else {
    CEP <- list(diag(L))
  }
  if(vis){
    cat("  Posterior probabilities and modal classifications were computed.\n")
    cat("  The CEP matrix was computed.\n")
  }

  if (is.null(covariates)) {
    covariates <- matrix(1, nrow = N, ncol = 1)
    colnames(covariates) <- "Intercept"
  } else {
    covariates <- as.matrix(covariates)
    if (ncol(covariates) == 0) {
      covariates <- matrix(1, nrow = N, ncol = 1)
      colnames(covariates) <- "Intercept"
    }
  }
  if(nrow(covariates) != N || any(!is.finite(covariates))){
    stop("covariates must contain N rows and only finite values")
  }
  covariates.list <- list(covariates)
  p <- ncol(covariates)
  npar.beta <- p * (L - 1)
  npar.gamma <- 0
  par.ini.step3 <- rep(0, npar.beta + npar.gamma)

  .three.step.progress.start(
    3L, "LCPA", path = .progress.path, vis = vis.step3
  )

  npar <- length(par.ini.step3)
  lb <- rep(lower, npar)
  ub <- rep(upper, npar)

  optimization <- .three.step.optimize(
    par.ini.step3, CEP, P.Z.Xns, Zs, covariates.list,
    covariates.time.cross = FALSE, ref.class = ref.class,
    lower = lower, upper = upper, tol = tol, maxiter = maxiter,
    gradient = gradient.cur, vis = vis.step3, progress.prefix = "  "
  )
  optimization.result <- optimization$result
  Log.Lik.history <- optimization$Log.Lik.history
  optimization.result.main <- optimization.result
  if(vis.step3){
    cat("\n")
  }

  params.step3 <- optimization.result$solution
  bound.diagnostics <- .three.step.bound.diagnostics(params.step3, lb, ub)

  vcov.step3 <- NULL
  information.step3 <- NULL
  SE.diagnostics <- list(method = method.SE)
  if (method.SE.internal != "Bootstrap") {
    derivatives <- .three.step.derivatives(
      params.step3, CEP, Zs, covariates.list, FALSE, ref.class, FALSE
    )
  }
  if(method.SE.internal == "Louis"){
    if(vis.step3){
      cat("  Calculating Analytic Observed Information via Louis' Identity ...\n")
    }
    se.result <- .three.step.bootstrap.target.se(
      derivatives$information, derivatives,
      P.Z.Xns, Zs, CEP,
      CEP.time.cross = FALSE,
      CEP.reestimated = CEP.error,
      label = "Louis"
    )
    se.vec <- se.result$se
    vcov.step3 <- se.result$vcov
    information.step3 <- se.result$information
    SE.diagnostics <- se.result$diagnostics
  }else if(method.SE.internal == "Obs"){
    if(vis.step3){
      cat("  Calculating Observed Information Matrix for Standard Errors ...\n")
    }
    h <- .three.step.numeric.information(
      params.step3,
      function(x) {
        get.Log.Lik.LTA.optim(
          x, CEP, P.Z.Xns, Zs, covariates.list, FALSE, ref.class
        )
      }
    )

    if (!is.null(h)) {
      se.result <- .three.step.bootstrap.target.se(
        h, derivatives,
        P.Z.Xns, Zs, CEP,
        CEP.time.cross = FALSE,
        CEP.reestimated = CEP.error,
        label = "Obs"
      )
      se.vec <- se.result$se
      vcov.step3 <- se.result$vcov
      information.step3 <- se.result$information
      SE.diagnostics <- se.result$diagnostics
    }else{
      se.vec <- rep(NA_real_, npar)
      SE.diagnostics <- list(method = "Obs", invertible = FALSE)
    }
  }else{
    if(vis.step3){
      cat("  Bootstrapping for Standard Errors ...\n")
    }
    params.bootstrap <- matrix(NA_real_, nrep.bootstrap, npar)
    progress.state <- .new.progress.state()
    for(bs in 1:nrep.bootstrap){
      covariates.cur <- vector("list", 1)
      samples.cur <- .three.step.bootstrap.indices(N)
      covariates.cur[[1]] <- covariates.list[[1]][samples.cur, , drop=FALSE]

      P.Z.Xns.cur <- Zs.cur <- P.Zs.cur <- vector("list", 1)
      P.Z.Xns.cur[[1]] <- .three.step.posterior(
        response[samples.cur, , drop = FALSE], params, type.model
      )

      P.Zs.cur[[1]] <- colSums(P.Z.Xns.cur[[1]]) / sum(P.Z.Xns.cur[[1]])
      Zs.cur[[1]] <- max.col(P.Z.Xns.cur[[1]], ties.method = "first")

      if(CEP.error){
        CEP.cur <- get.CEP(P.Z.Xns.cur, CEP.time.cross =FALSE)
      } else {
        CEP.cur <- replicate(1, diag(L), simplify=FALSE)
      }

      optimization <- .three.step.optimize(
        params.step3, CEP.cur, P.Z.Xns.cur, Zs.cur, covariates.cur,
        covariates.time.cross = FALSE, ref.class = ref.class,
        lower = lower, upper = upper, tol = tol, maxiter = maxiter,
        gradient = gradient.cur, vis = vis.step3,
        progress.prefix = paste0("  Bootstrap = ", bs, "/", nrep.bootstrap, " | "),
        progress.state = progress.state
      )
      optimization.result <- optimization$result
      if (optimization.result$status %in% 1:4 &&
          all(is.finite(optimization.result$solution)) && is.finite(optimization.result$objective)) {
        params.bootstrap[bs, ] <- optimization.result$solution
      }
    }

    successful <- complete.cases(params.bootstrap)
    if (sum(successful) < 2L) {
      warning("Fewer than two bootstrap replicates converged; bootstrap standard errors are unavailable")
      se.vec <- rep(NA_real_, npar)
    } else {
      se.vec <- apply(params.bootstrap[successful, , drop = FALSE], 2, sd)
      vcov.step3 <- stats::cov(params.bootstrap[successful, , drop = FALSE])
    }
    SE.diagnostics <- list(
      method = "Bootstrap", successful = sum(successful),
      attempted = nrep.bootstrap
    )
    if(vis.step3){
      cat("\n\n")
    }
  }

  SE.diagnostics$method <- method.SE
  params.LTA.obj <- LTA.vector.to.parameters(params.step3, covariates.list, L, ref.class)

  if(!is.null(se.vec)){
    se.obj <- .three.step.reference.na(
      LTA.vector.to.parameters(se.vec, covariates.list, L, ref.class), ref.class
    )
    Z.sta.vec <- params.step3/se.vec
    p.value.tail1 <- pnorm(-abs(Z.sta.vec))
    p.value.tail2 <- pnorm(-abs(Z.sta.vec)) * 2
    Z.sta.obj <- .three.step.reference.na(
      LTA.vector.to.parameters(Z.sta.vec, covariates.list, L, ref.class), ref.class
    )
    p.value.tail1.obj <- .three.step.reference.na(
      LTA.vector.to.parameters(p.value.tail1, covariates.list, L, ref.class), ref.class
    )
    p.value.tail2.obj <- .three.step.reference.na(
      LTA.vector.to.parameters(p.value.tail2, covariates.list, L, ref.class), ref.class
    )
  }

  beta <- params.LTA.obj$beta
  beta.se <- se.obj$beta
  beta.Z.sta <- Z.sta.obj$beta
  beta.p.value.tail1 <- p.value.tail1.obj$beta
  beta.p.value.tail2 <- p.value.tail2.obj$beta

  covariates.ncol <- unlist(lapply(covariates.list, ncol))
  npar <- get.npar.LTA(covariates.ncol, L, FALSE)

  Log.Lik = -get.Log.Lik.LTA.optim(params.step3, CEP, P.Z.Xns, Zs, covariates.list, FALSE, ref.class)
  AIC <- -2 * Log.Lik + 2 * npar
  BIC <- -2 * Log.Lik + npar * log(N)

  # Prepare result object
  res <- list(
    beta = beta,
    beta.se = beta.se,
    beta.Z.sta = beta.Z.sta,
    beta.p.value.tail1 = beta.p.value.tail1,
    beta.p.value.tail2 = beta.p.value.tail2,
    vcov = vcov.step3,
    information = information.step3,
    SE.diagnostics = SE.diagnostics,
    bound.diagnostics = bound.diagnostics,
    npar = npar,
    Log.Lik = Log.Lik,
    AIC = AIC,
    BIC = BIC,
    P.Z.Xn = P.Z.Xns[[1]],
    P.Z = P.Zs[[1]],
    Z = Zs[[1]],
    Log.Lik.history = Log.Lik.history,
    iterations = optimization.result.main$iterations,
    converged = optimization.result.main$status %in% 1:4,
    params = params,
    call = call,
    arguments = list(
      response = response, L = L, type.model = type.model,
      covariates = covariates,
      CEP.error = CEP.error,
      par.ini = par.ini,
      params = params,
      constraint = constraint,
      method.model = method.model, tol = tol,
      method.regression = method.regression,
      lower = lower, upper = upper,
      method.SE = method.SE, nrep.bootstrap = nrep.bootstrap,
      maxiter = maxiter, is.sort = is.sort,
      nrep = nrep, starts = starts, maxiter.warmup = maxiter.warmup,
      vis = vis,
      control.EM = control.EM,
      control.Mplus = control.Mplus,
      control.NNE = control.NNE,
      control.flexmix = control.flexmix,
      control.Rmixmod = control.Rmixmod,
      control.RMixtComp = control.RMixtComp,
      ref.class = ref.class
    )
  )

  class(res) <- "LCPA"
  return(res)
}
