#' Latent Transition Analysis (LTA)
#'
#' Implements the three-step estimation method (Vermunt, 2010; Liang et al., 2023) for Latent Transition Analysis (LTA),
#' treating latent class memberships at each time point as observed variables with measurement error.
#' Classification uncertainty from Step 1 (latent class/profile analysis) is explicitly incorporated
#' into the transition model estimation in Step 3, ensuring asymptotically unbiased estimates of
#' transition probabilities and covariate effects. This avoids the bias introduced by "hard" modal-class assignment.
#'
#' @param responses A \code{list} of response matrices or data frames. Each matrix corresponds to one time point.
#'                  Rows of each matrix represent individuals/participants/observations (\eqn{N}), columns of each
#'                  matrix represent observed items/variables (\eqn{I}).
#'                  For \code{type.model = "LCA"}: items must be binary or categorical (coded as integers starting from 0).
#'                  For \code{type.model = "LPA"}: items must be continuous (numeric), and each response matrix must be
#'                  standardized using \code{\link[base]{scale}} or \code{\link[LCPA]{normalize}} prior to input.
#' @param L Integer scalar. Number of latent classes/profiles at each time point. Must satisfy \eqn{L \geq 2}.
#' @param ref.class Integer \eqn{L \geq ref.class \geq 1}. Specifies which first-time-point latent class
#'                  to use as the reference category. Default is \code{L} (last class). Coefficients for
#'                  the reference class are fixed to zero. When \code{is.sort=TRUE}, the first-time-point
#'                  classes are first ordered by decreasing \code{P.Z}; \code{ref.class} refers to the
#'                  position in this established order, which is retained at every later time point.
#' @param type.model Character string. Specifies the type of latent variable model for Step 1:
#'             \itemize{
#'               \item \code{"LCA"} — Latent Class Analysis for categorical items.
#'               \item \code{"LPA"} — Latent Profile Analysis for continuous items.
#'             }
#'             See \code{\link[LCPA]{LCA}} and \code{\link[LCPA]{LPA}} for details.
#' @param covariates Optional. A \code{list} of matrices/data frames (length = number of time points).
#'                   Each matrix contains covariates for modeling transitions or initial status.
#'                   Must include an intercept column (all 1s) as the first column.
#'                   All non-intercept covariates must be standardized before analysis, for example
#'                   with \code{scale()} or \code{\link[LCPA]{normalize}}. This places the covariate regression
#'                   coefficients on comparable scales and prevents their magnitudes from being driven by the
#'                   covariates' original measurement units. Do not standardize the intercept column; construct
#'                   interaction terms from the standardized covariates.
#'                   If \code{NULL} (default), only intercept terms are used (i.e., no covariates).
#'                   For time \eqn{t}, dimension is \eqn{N\times(U_t+1)}.
#'                   Covariates can vary across time.
#' @param CEP.time.cross Logical. If \code{TRUE}, assumes measurement invariance and uses the same
#'                      Classification Error Probability (\code{\link[LCPA]{get.CEP}}) matrix across all time points.
#'                      Requires that item parameters are invariant over time (not checked internally).
#'                      Default is \code{FALSE}.
#' @param CEP.error Logical. If \code{TRUE} (recommended), incorporates classification uncertainty via
#'                  estimated CEP matrices from Step 1. If \code{FALSE}, uses identity CEP matrices
#'                  (equivalent to naive modal assignment; introduces bias and not recommended).
#' @param covariates.time.cross Logical. If \code{TRUE}, forces the use of identical \eqn{\gamma} parameters
#'                             across all time points (i.e., a time-invariant probability transition matrix).
#'                             In this case, users should ensure that the covariate matrices at different time points
#'                             have the same dimensions (values may differ) to
#'                             match the fixed form of the
#'                             \eqn{\boldsymbol{\gamma}_{klt}} coefficients.
#'                             Default is \code{FALSE}, allowing for potentially different probability transition matrices across time points.
#' @param par.ini Specification for parameter initialization. Options include:
#'   \itemize{
#'     \item \code{"random"}: Completely random initialization (default).
#'     \item \code{"kmeans"}: Initializes parameters via K-means clustering on observed data (McLachlan & Peel, 2000).
#'     \item A \code{list} for LCA containing:
#'       \describe{
#'         \item{\code{par}}{An \eqn{L \times I \times K_{\max}} array of initial conditional probabilities for
#'                           each latent class, item, and response category (where \eqn{K_{\max}} is the maximum
#'                           number of categories across items).}
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
#'                                        each latent class, item, and response category (where \eqn{K_{\max}} is the maximum
#'                                        number of categories across items).}
#'                      \item{\code{P.Z}}{A numeric vector of length \eqn{L} specifying initial prior probabilities for latent classes.}
#'                      \item{\code{category.levels}}{Fixed ordered response categories for each item.}
#'                    }
#'                 \item A \code{list} for LPA containing:
#'                    \describe{
#'                      \item{\code{means}}{An \eqn{L \times I} matrix of initial mean vectors for each profile.}
#'                      \item{\code{covs}}{An \eqn{I \times I \times L} array of initial covariance matrices for each profile.}
#'                      \item{\code{P.Z}}{A numeric vector of length \eqn{L} specifying initial prior probabilities for profiles.}
#'                    }
#'               }
#' @param step1.pool Logical. If \code{FALSE} (default), estimates the common Step 1 measurement
#'                   model from the first time point and applies it to every time point. If \code{TRUE},
#'                   estimates one common measurement model from the row-bound responses across all
#'                   time points and then splits the posterior probabilities back by time point.
#'                   This option is ignored when \code{params} is supplied because those Step 1
#'                   parameters are used as fixed values.
#' @param is.sort A logical value. If \code{TRUE} (default), the first-time-point latent classes are
#'                ordered by decreasing \code{P.Z} and that order is retained at every later time point.
#'                If \code{FALSE}, the original Step 1 order is retained.
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
#'                 \item \code{"Rmixmod"} — Package-staged SEM or a native EM/CEM/SEM strategy
#'                   through \code{Rmixmod}.
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
#'                       More robust to model misspecification and small-sample bias. Computationally intensive but recommended when asymptotic assumptions are questionable.
#'               }
#'               Default is \code{"Bootstrap"}.
#' @param nrep.bootstrap Integer. Number of bootstrap replicates used when \code{method.SE = "Bootstrap"}.
#'                    Default is 100. McLachlan & Peel (2000) suggest that 50–100 replicates often provide adequate accuracy
#'                    for practical purposes, though more (e.g., 500–1000) may be preferred for publication-quality inference.
#'                    Each replicate resamples individuals jointly over time, recomputes posterior assignments
#'                    and CEP from the fixed Step 1 measurement parameters, and re-estimates Step 3.
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
#'   These three staged-training arguments are not used when
#'   \code{control.Rmixmod$path = "Rmixmod"}; that path uses its documented native controls instead.
#' @param vis Logical. If \code{TRUE}, displays the three-step stage headings and
#'            indented process information for Step 1 model estimation, Step 2
#'            time-specific posterior classification and CEP-matrix preparation,
#'            and Step 3 optimization and standard-error estimation (default: \code{TRUE}).
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
#' @param control.Rmixmod List of control parameters for Rmixmod estimation:
#'   \describe{
#'     \item{\code{path}}{Execution path: \code{"LCPA"} (default) uses the package's
#'       warm-up and promoted SEM runs; \code{"Rmixmod"} delegates one native EM, CEM,
#'       SEM, or combined strategy and ignores the public staged-training arguments.}
#'     \item{\code{algorithm}, \code{nrep}, \code{method.init}, \code{starts},
#'       \code{maxiter.init}, \code{maxiter}, \code{tol.init}, \code{tol},
#'       \code{par.ini}, \code{labels.ini}}{Package-standard controls translated to
#'       \code{Rmixmod::mixmodStrategy()} when \code{path = "Rmixmod"}. Under
#'       \code{path = "LCPA"}, only \code{maxiter} is used, with 1000 iterations when omitted.
#'       \code{algorithm} accepts \code{"EM"}, \code{"CEM"}, and \code{"SEM"}, including
#'       ordered combinations supported by Rmixmod.}
#'     \item{\code{strategy}}{Optional pre-built Rmixmod \code{Strategy} object for
#'       \code{path = "Rmixmod"}; it takes precedence over individual strategy controls.}
#'   }
#'   For LPA Step 1 models, Rmixmod supports \code{"V0"}, \code{"EE"}, and \code{"VV"} covariance constraints.
#'
#' @section Random-number reproducibility:
#' Except for \code{method = "NNE"}, which intentionally uses its fixed backend seed,
#' all stochastic stages and supported Step 1 backends are driven from R's current
#' random-number generator. The user only needs to call \code{set.seed()} immediately
#' before \code{\link[LCPA]{LTA}()}; no Rmixmod or Mplus-specific seed setting is required.
#'
#' @return An object of class \code{LTA}, a named list containing:
#' \describe{
#'   \item{\code{beta}}{Matrix of size \eqn{(U_1+1)\times L}. Coefficients for initial class membership multinomial logit model.
#'     All non-reference columns are free; column \code{ref.class} is constrained to zero.}
#'   \item{\code{gamma}}{List of length \eqn{T-1}. Each element \code{gamma[[t]]} (for transition from time \eqn{t} to \eqn{t+1})
#'     is a nested list: \code{gamma[[t]][[from_class]][[to_class]]} returns a
#'     coefficient vector of length \eqn{U_{t+1}+1}.
#'     Coefficients for transitions to \code{ref.class} are fixed to zero for every origin class.}
#'   \item{\code{beta.se}}{Standard errors for \code{beta}. The fixed reference-class column is \code{NA};
#'     free coefficients may also be \code{NA} if the information matrix is ill-conditioned.}
#'   \item{\code{gamma.se}}{Standard errors for \code{gamma}; coefficients for transitions to the fixed
#'     reference class are \code{NA}.}
#'   \item{\code{beta.Z.sta}}{Z-statistics for testing null hypothesis that each beta coefficient equals zero.
#'     Computed as \code{beta / beta.se}. Same structure as \code{beta}.}
#'   \item{\code{gamma.Z.sta}}{Z-statistics for gamma coefficients. Same nested structure as \code{gamma}.
#'     Used for testing significance of transition effects.}
#'   \item{\code{beta.p.value.tail1}}{One-tailed p-values based on standard normal distribution: \eqn{P(Z < -|z|)}.
#'     Useful for directional hypotheses. Same structure as \code{beta}.}
#'   \item{\code{gamma.p.value.tail1}}{One-tailed p-values for gamma coefficients. Same nested structure as \code{gamma}.}
#'   \item{\code{beta.p.value.tail2}}{Two-tailed p-values: \eqn{2 \times P(Z < -|z|)}.
#'     Standard test for non-zero effect. Same structure as \code{beta}.}
#'   \item{\code{gamma.p.value.tail2}}{Two-tailed p-values for gamma coefficients. Same nested structure as \code{gamma}.}
#'   \item{\code{vcov}}{Variance-covariance matrix of the free Step 3 coefficients. For bootstrap,
#'     this is the empirical covariance matrix across successful replicates.}
#'   \item{\code{information}}{Observed information matrix for \code{method.SE = "Numeric"} or \code{"Analytic"};
#'     \code{NULL} for bootstrap.}
#'   \item{\code{SE.diagnostics}}{Conditioning or bootstrap-convergence diagnostics for the selected SE method.}
#'   \item{\code{bound.diagnostics}}{Step 3 coefficient indices within \code{1e-4} of an optimization bound.}
#'   \item{\code{P.Z.Xns}}{List of length \eqn{T}. Each element is an \eqn{N \times L} matrix of posterior class probabilities
#'     \eqn{P(Z_{nt}=l \mid \mathbf{X}_{nt})} for each participant \eqn{n} at time \eqn{t}.}
#'   \item{\code{P.Zs}}{List of length \eqn{T}. Each element is a vector of length \eqn{L} containing prior class proportions
#'     \eqn{P(Z_t = l)} estimated at Step 1 for time \eqn{t}.}
#'   \item{\code{Zs}}{List of length \eqn{T}. Each element is a vector of length \eqn{N} containing modal class assignments
#'     (MAP classifications) \eqn{\widehat{Z}_{nt}} for each participant at time \eqn{t}.}
#'   \item{\code{npar}}{Number of free parameters in the model (depends on \code{covariates}).}
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
#' The three-step LTA proceeds as follows:
#'
#' Step 1 — Unconditional Latent Class/Profile Model:
#' By default, fit an unconditional LCA or LPA model at the first time point and apply its common
#' measurement parameters at every time point. With \code{step1.pool = TRUE}, fit the common measurement
#' model to the row-bound responses from all time points instead.
#' Obtain posterior class membership probabilities
#' \eqn{P(Z_{nt}=l \mid \mathbf{X}_{nt})} for each participant \eqn{n}
#' and class \eqn{l} using Bayes' theorem.
#'
#' Step 2 — Classification Error Probabilities (equal to \code{\link[LCPA]{get.CEP}}):
#' Compute the time-specific \eqn{L\times L} matrices with
#' \eqn{\mathrm{CEP}_t(l,k)=P(\widehat{Z}_{nt}=k\mid Z_{nt}=l)}. Their modal
#' assignments, posterior-weight estimators, orientation, and optional
#' cross-time pooling are defined in \code{\link[LCPA]{get.CEP}()}.
#'
#' Step 3 — Transition Model with Measurement Error Correction:
#' Estimate the multinomial logit models for:
#' \itemize{
#'   \item Initial class membership (time 1):
#'     \eqn{P(Z_{n1}=l\mid\boldsymbol{\zeta}_{n1})=
#'     \frac{\exp(\boldsymbol{\beta}_l^\top\boldsymbol{\zeta}_{n1})}
#'     {\sum_{h=1}^L\exp(\boldsymbol{\beta}_h^\top\boldsymbol{\zeta}_{n1})}}
#'   \item Transitions (time \eqn{t > 1}):
#'     \eqn{P(Z_{nt}=l\mid Z_{n,t-1}=k,\boldsymbol{\zeta}_{nt})=
#'     \frac{\exp(\boldsymbol{\gamma}_{klt}^\top\boldsymbol{\zeta}_{nt})}
#'     {\sum_{h=1}^L\exp(\boldsymbol{\gamma}_{kht}^\top
#'     \boldsymbol{\zeta}_{nt})}}
#' }
#' where
#' \eqn{\boldsymbol{\zeta}_{n1}=
#' (1,\zeta_{n11},\ldots,\zeta_{n1U_1})^\top} is the covariate vector for
#' participant \eqn{n} at time 1. The coefficient vector
#' \eqn{\boldsymbol{\beta}_l=
#' (\beta_{l0},\beta_{l1},\ldots,\beta_{lU_1})^\top} contains the
#' corresponding intercept and slopes.
#' The class selected by \code{ref.class} is the reference class and its beta vector is fixed to zero.
#' At time \eqn{t},
#' \eqn{\boldsymbol{\zeta}_{nt}=
#' (1,\zeta_{nt1},\ldots,\zeta_{ntU_t})^\top}, and
#' \eqn{\boldsymbol{\gamma}_{klt}=
#' (\gamma_{klt0},\gamma_{klt1},\ldots,\gamma_{kltU_t})^\top} contains the
#' corresponding intercept and slopes. Coefficients for transitions to
#' \code{ref.class} are fixed to zero.
#'
#' The full observed-data likelihood integrates over all possible latent class paths:
#' \deqn{
#' \begin{aligned}
#' \log \mathcal{L}(\boldsymbol{\beta},\boldsymbol{\gamma}) &=
#' \sum_{n=1}^N \log \Biggl[
#'   \sum_{z_{n1}=1}^L\cdots\sum_{z_{nT}=1}^L
#'   \Bigl(\prod_{t=1}^T
#'   \mathrm{CEP}_t(z_{nt},\widehat{Z}_{nt})\Bigr) \cdot \\
#'   &\quad P(Z_{n1}=z_{n1}\mid\boldsymbol{\zeta}_{n1}) \cdot \\
#'   &\quad \prod_{t=2}^T P(Z_{nt}=z_{nt}\mid
#'   Z_{n,t-1}=z_{n,t-1},\boldsymbol{\zeta}_{nt})
#' \Biggr]
#' \end{aligned}
#' }
#' The coefficients \eqn{\boldsymbol{\beta}} and \eqn{\boldsymbol{\gamma}}
#' are estimated via maximum likelihood using the L-BFGS algorithm
#' (box-constrained gradient optimization). The destination selected by
#' \code{ref.class} has zero beta and gamma coefficient vectors.
#'
#' @section Bootstrap Standard Error Estimation:
#' When \code{method.SE = "Bootstrap"}, standard errors are estimated using a nonparametric bootstrap procedure:
#' \enumerate{
#'   \item Draw \eqn{B} (=\code{nrep.bootstrap}) independent samples of size \eqn{N} with replacement from the original data.
#'   \item For each bootstrap sample \eqn{b=1,\dots,B}, retain the fitted Step 1 measurement parameters,
#'         recompute posterior assignments and CEP, and re-estimate the Step 3 parameter vector
#'         \eqn{\hat{\boldsymbol{\theta}}^{(b)}}.
#'   \item Compute the bootstrap standard error for each parameter as the sample standard deviation across replicates:
#'         \deqn{
#'           \widehat{\mathrm{SE}}_{\mathrm{boot}}(\hat{\theta}_j) = \sqrt{ \frac{1}{B-1} \sum_{b=1}^B \left( \hat{\theta}_j^{(b)} - \bar{\theta}_j \right)^2 },
#'         }
#'         where \eqn{\bar{\theta}_j = \frac{1}{B}\sum_{b=1}^B \hat{\theta}_j^{(b)}}.
#' }
#' This approach does not rely on large-sample normality or correct specification of the information matrix,
#' making it particularly suitable for complex models like LTA where analytic derivatives are difficult or unstable.
#' However, it increases computational cost linearly with \eqn{B}.
#'
#' @section Important Implementation Details:
#' \itemize{
#'   \item Reference Class: The class selected by \code{ref.class} is treated as the reference category.
#'         Its corresponding destination coefficients in \code{beta} and \code{gamma} are fixed to zero.
#'   \item CEP Matrices: When \code{CEP.error = TRUE}, misclassification probabilities are estimated
#'                       non-parametrically using Step 1 posterior probabilities. This corrects for
#'                       classification uncertainty. Setting \code{CEP.time.cross = TRUE} assumes these
#'                       error structures are identical across time (measurement invariance).
#'                       See in \code{\link[LCPA]{get.CEP}}.
#'   \item Covariate Handling: Covariates for initial status (time 1) and transitions (time \eqn{t \geq 2}) can differ.
#'         For transitions to time \eqn{t}, the covariate matrix must have
#'         dimensions \eqn{N\times(U_t+1)}, i.e., an intercept column of all
#'         \eqn{1} plus \eqn{U_t} covariates.
#'   \item Covariate standardization: All non-intercept covariates must be standardized before
#'         analysis (e.g., with \code{scale()} or \code{\link[LCPA]{normalize}}). Otherwise, the estimated
#'         covariate regression coefficients remain tied to incompatible original measurement units, creating
#'         a scale problem that prevents meaningful comparison of coefficient magnitudes. The intercept must
#'         remain equal to 1, and interaction terms should be formed after standardizing their component covariates.
#'         When \code{covariates.time.cross = TRUE}, use the same centering and scaling constants at every time point
#'         so that the equality constraint applies to coefficients defined on the same covariate scale.
#'   \item Optimization: Step 3 uses box-constrained BOBYQA via \code{\link[nloptr]{nloptr}}.
#'         \code{method.SE = "Analytic"} computes the sandwich bread analytically and \code{method.SE = "Numeric"}
#'         approximates it numerically. Both use the same empirical score and CEP influence-function meat,
#'         targeting the fixed-Step-1 nonparametric bootstrap sampling scheme.
#'   \item Computational Complexity: Step 3 likelihood and gradient evaluation use scaled
#'         forward--backward recursions with \eqn{O(NTL^2)} time rather than explicit
#'         enumeration of \eqn{L^T} latent paths.
#'   \item Bootstrap Computation: Each bootstrap iteration keeps the Step 1 measurement parameters fixed,
#'         recomputes \code{P.Z.Xns} and \code{CEP}, and re-estimates the transition parameters. Failed
#'         optimization runs are excluded rather than entered as zero vectors.
#'         To ensure reproducibility, set a seed before calling \code{\link[LCPA]{LTA}()} when using \code{method.SE = "Bootstrap"}.
#'         Progress messages during bootstrapping include current replicate index and optimization diagnostics.
#'         Users should monitor convergence in each bootstrap run; failed runs will result in \code{NA} entries in SEs and derived statistics.
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
#' McLachlan, G. J., & Peel, D. (2000). *Finite mixture models*. John Wiley & Sons.
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
#' N <- 2000 ## sample size
#' L <- 3   ## number of latent class
#' I <- 6   ## number of variables/items
#'
#' ## Covariates at time point T1
#' covariates.inter <- rep(1, N)                 # Intercept term is always 1 for each individual
#' covariates.X11 <- as.numeric(scale(rnorm(N))) # Standardized continuous covariate
#' # Combine into covariates at T1
#' covariates.T1 <- cbind(Intercept=covariates.inter, X1=covariates.X11)
#' ## Covariates at time point T2
#' covariates.inter <- rep(1, N)                 # Intercept term is always 1 for each individual
#' covariates.X21 <- as.numeric(scale(rnorm(N))) # Standardized continuous covariate
#' # Combine into covariates at T1
#' covariates.T2 <- cbind(Intercept=covariates.inter, X1=covariates.X21)
#'
#' # Combine into final covariates list
#' covariates <- list(t1=covariates.T1, t2=covariates.T2)
#'
#' ## Simulate beta coefficients
#' ## fix reference class to class 2
#' beta <- matrix(c( 0.7, 0.0, -0.1,
#'                   0.2, 0.0,  0.5), ncol=L, byrow=TRUE)
#'
#' ## Simulate gamma coefficients
#' gamma <- list(
#'   lapply(1:L, function(l) {
#'     lapply(1:L, function(k) if(k != 2)
#'            runif(2, -2.0, 2.0) else c(0, 0)) # Class 2 as reference
#'   })
#' )
#'
#' ## Simulate the data
#' sim_custom <- sim.LTA(
#'   N=N, I=I, L=L, times=2, type="LCA", IQ=0.9,
#'   ref.class=2,
#'   covariates=covariates,
#'   beta=beta,
#'   gamma=gamma
#' )
#' summary(sim_custom)
#' responses <- sim_custom$responses
#' covariates <- sim_custom$covariates
#'
#' ## fix reference class to class 2
#' LTA.obj <- LTA(responses, L=L, ref.class=2, type.model="LCA",
#'                covariates=covariates,
#'                method.SE="Bootstrap", nrep.bootstrap=100,
#'                CEP.time.cross =FALSE, CEP.error=TRUE, covariates.time.cross =FALSE,
#'                par.ini = "random", method="EM", vis = TRUE)
#'
#' print(LTA.obj)
#' }
#'
#' @importFrom nloptr nloptr
#' @importFrom Matrix nearPD
#' @importFrom numDeriv hessian
#' @importFrom stats pnorm complete.cases
#' @importFrom MASS ginv
#'
#' @noRd
ML.XZ.LTA <- function(responses, L=2,
                ref.class = L, type.model="LCA",
                covariates=NULL,
                CEP.time.cross =FALSE,
                CEP.error=TRUE,
                covariates.time.cross =FALSE,
                par.ini = "random",
                params=NULL, step1.pool=FALSE, is.sort=TRUE,
                constraint = "VV",
                method.model="EM", tol=1e-4,
                method.regression="Analytic",
                lower=-10, upper=10,
                method.SE="Bootstrap", nrep.bootstrap=100,
                maxiter=5000, starts=100,
                maxiter.warmup=20, nrep = 20,
                vis = TRUE,
                control.EM=NULL,
                control.Mplus=NULL,
                control.NNE=NULL,
                control.Rmixmod=NULL,
                .progress.path = "X -> Z",
                .progress.step3 = TRUE){

  call <- match.call()
  method.model <- match.arg(method.model, c("EM", "NNE", "Mplus", "Rmixmod"))
  if(is.null(params)) par.ini <- .normalize.par.ini(par.ini, method.model)
  type.model <- match.arg(type.model, c("LCA", "LPA"))
  method.regression <- match.arg(method.regression, c("Analytic", "Numeric"))
  method.SE <- match.arg(method.SE, c("Analytic", "Numeric", "Bootstrap"))
  method.SE.internal <- c(Analytic = "Louis", Numeric = "Obs", Bootstrap = "Bootstrap")[[method.SE]]
  gradient.cur <- tolower(method.regression)
  step1.pool <- isTRUE(step1.pool)
  params.input <- params
  vis.step3 <- isTRUE(vis) && isTRUE(.progress.step3)

  if(ref.class < 1 || ref.class > L) {
    stop("ref.class must be between 1 and L")
  }
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("tol must be a positive finite number")
  }
  if (method.SE == "Bootstrap" && nrep.bootstrap < 2L) {
    stop("nrep.bootstrap must be at least 2 when method.SE = 'Bootstrap'")
  }

  default_control.EM <- list(maxiter=2000, tol=1e-4)
  default_control.Mplus <- list(maxiter=2000, tol=1e-4, files.path = NULL, files.clean = TRUE)
  default_control.Rmixmod <- .default.Rmixmod.control(maxiter = 1000L)
  default_control.NNE <- list(
    hidden.layers = c(16, 16),
    activation.function = "tanh",
    use.attention = TRUE,
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
    device="CPU"
  )

  merge_and_clean_control <- function(user_control, default_control) {
    if (is.null(user_control)) return(default_control)
    merged <- modifyList(default_control, user_control)
    merged[names(default_control)]
  }

  control.EM     <- merge_and_clean_control(control.EM, default_control.EM)
  control.Mplus  <- merge_and_clean_control(control.Mplus, default_control.Mplus)
  control.NNE    <- merge_and_clean_control(control.NNE, default_control.NNE)
  control.Rmixmod <- merge_and_clean_control(control.Rmixmod, default_control.Rmixmod)

  if (!is.list(responses) || length(responses) < 1L) {
    stop("responses must be a non-empty list of response matrices")
  }
  responses <- lapply(responses, as.matrix)
  times <- length(responses)
  N <- nrow(responses[[1]])
  I <- ncol(responses[[1]])
  if (any(vapply(responses, nrow, integer(1)) != N) ||
      any(vapply(responses, ncol, integer(1)) != I)) {
    stop("All response matrices must have the same numbers of individuals and indicators")
  }
  response <- if(step1.pool) do.call(rbind, responses) else responses[[1]]

  .three.step.progress.start(
    1L, "LTA", type.model = type.model, vis = vis
  )

  if(is.null(params)){
    LCPA.obj <- .with.estimation.output.prefix({
      if(type.model == "LCA"){
        LCA(response, L=L, par.ini=par.ini,
            method=method.model, is.sort=FALSE, nrep=nrep,
            starts=starts, maxiter.warmup=maxiter.warmup,
            vis=vis,
            control.EM=control.EM,
            control.Mplus=control.Mplus,
            control.NNE=control.NNE,
            control.Rmixmod=control.Rmixmod)
      } else {
        LPA(response, L=L, par.ini=par.ini,
            constraint=constraint,
            method=method.model, is.sort=FALSE, nrep=nrep,
            starts=starts, maxiter.warmup=maxiter.warmup,
            vis=vis,
            control.EM=control.EM,
            control.Mplus=control.Mplus,
            control.NNE=control.NNE,
            control.Rmixmod=control.Rmixmod)
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
  if(is.sort){
    posterior.t1 <- .three.step.posterior(responses[[1]], params, type.model)
    P.Z.t1 <- colSums(posterior.t1) / sum(posterior.t1)
    position <- order(P.Z.t1, decreasing = TRUE)
    params <- .three.step.reorder.params(params, position, type.model)
  }

  .three.step.progress.start(
    2L, "LTA", path = .progress.path, vis = vis
  )

  P.Z.Xns <- Zs <- P.Zs <- vector("list", times)
  if(step1.pool){
    response.pooled <- do.call(rbind, responses)
    posterior.pooled <- .three.step.posterior(response.pooled, params, type.model)
    for(t in seq_len(times)){
      index.t <- ((t - 1L) * N + 1L):(t * N)
      P.Z.Xns[[t]] <- posterior.pooled[index.t, , drop=FALSE]
    }
  } else {
    for(t in seq_len(times)){
      P.Z.Xns[[t]] <- .three.step.posterior(responses[[t]], params, type.model)
    }
  }
  for(t in seq_len(times)){
    P.Zs[[t]] <- colSums(P.Z.Xns[[t]]) / sum(P.Z.Xns[[t]])
    Zs[[t]] <- max.col(P.Z.Xns[[t]], ties.method = "first")
  }

  if(CEP.error){
    CEP <- get.CEP(P.Z.Xns, CEP.time.cross =CEP.time.cross)
  } else {
    CEP <- replicate(times, diag(L), simplify=FALSE)
  }
  if(vis){
    cat("  Posterior probabilities and modal classifications were computed.\n")
    cat(sprintf("  CEP matrices were computed for %d time points.\n", times))
  }

  if(is.null(covariates)){
    covariates <- lapply(responses, function(x){
      matrix(1, nrow=nrow(x), ncol=1)
    })
  }
  if (!is.list(covariates) || length(covariates) != times) {
    stop("covariates must be a list with one design matrix per time point")
  }
  covariates <- lapply(covariates, as.matrix)
  if (any(vapply(covariates, nrow, integer(1)) != N) ||
      any(!vapply(covariates, function(x) all(is.finite(x)), logical(1)))) {
    stop("Every covariate design matrix must contain N rows and only finite values")
  }
  if (covariates.time.cross && times > 2L &&
      length(unique(vapply(covariates[-1L], ncol, integer(1)))) != 1L) {
    stop("covariates.time.cross = TRUE requires identical transition designs at times 2, ..., T")
  }

  p1 <- ncol(covariates[[1]])
  npar.beta <- p1 * (L - 1)
  npar.gamma <- 0L
  if(times > 1L){
    if(covariates.time.cross){
      npar.gamma <- L * (L - 1) * ncol(covariates[[2]])
    } else {
      npar.gamma <- sum(vapply(2:times, function(t){
        L * (L - 1) * ncol(covariates[[t]])
      }, numeric(1)))
    }
  }
  par.ini.step3 <- rep(0, npar.beta + npar.gamma)

  .three.step.progress.start(
    3L, "LTA", path = .progress.path, vis = vis.step3
  )

  npar <- length(par.ini.step3)
  lb <- rep(lower, npar)
  ub <- rep(upper, npar)

  optimization <- .three.step.optimize(
    par.ini.step3, CEP, P.Z.Xns, Zs, covariates,
    covariates.time.cross = covariates.time.cross, ref.class = ref.class,
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
      params.step3, CEP, Zs, covariates,
      covariates.time.cross, ref.class, CEP.time.cross
    )
  }
  if(method.SE.internal == "Louis"){
    if(vis.step3){
      cat("  Calculating Analytic Observed Information via Louis' Identity ...\n")
    }
    se.result <- .three.step.bootstrap.target.se(
      derivatives$information, derivatives,
      P.Z.Xns, Zs, CEP,
      CEP.time.cross = CEP.time.cross,
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
          x, CEP, P.Z.Xns, Zs, covariates,
          covariates.time.cross, ref.class
        )
      }
    )

    if (!is.null(h)) {
      se.result <- .three.step.bootstrap.target.se(
        h, derivatives,
        P.Z.Xns, Zs, CEP,
        CEP.time.cross = CEP.time.cross,
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

      covariates.cur <- vector("list", times)
      samples.cur <- .three.step.bootstrap.indices(N)
      for(t in 1:times){
        covariates.cur[[t]] <- covariates[[t]][samples.cur, , drop=FALSE]
      }

      P.Z.Xns.cur <- Zs.cur <- P.Zs.cur <- vector("list", times)
      if(step1.pool){
        response.cur <- do.call(rbind, lapply(responses, function(x){
          x[samples.cur, , drop=FALSE]
        }))
        posterior.cur <- .three.step.posterior(response.cur, params, type.model)
        for(t in seq_len(times)){
          index.t <- ((t - 1L) * N + 1L):(t * N)
          P.Z.Xns.cur[[t]] <- posterior.cur[index.t, , drop=FALSE]
        }
      } else {
        for(t in seq_len(times)){
          P.Z.Xns.cur[[t]] <- .three.step.posterior(
            responses[[t]][samples.cur, , drop = FALSE], params, type.model
          )
        }
      }
      for(t in seq_len(times)){
        P.Zs.cur[[t]] <- colSums(P.Z.Xns.cur[[t]]) / sum(P.Z.Xns.cur[[t]])
        Zs.cur[[t]] <- max.col(P.Z.Xns.cur[[t]], ties.method = "first")
      }

      if(CEP.error){
        CEP.cur <- get.CEP(P.Z.Xns.cur, CEP.time.cross =CEP.time.cross)
      } else {
        CEP.cur <- replicate(times, diag(L), simplify=FALSE)
      }

      optimization <- .three.step.optimize(
        params.step3, CEP.cur, P.Z.Xns.cur, Zs.cur, covariates.cur,
        covariates.time.cross = covariates.time.cross, ref.class = ref.class,
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
      cat("\n")
    }
  }

  SE.diagnostics$method <- method.SE
  params.LTA.obj <- LTA.vector.to.parameters(
    params.step3, covariates, L, ref.class, covariates.time.cross
  )

  if(!is.null(se.vec)){
    se.obj <- .three.step.reference.na(
      LTA.vector.to.parameters(
        se.vec, covariates, L, ref.class, covariates.time.cross
      ), ref.class
    )
    Z.sta.vec <- params.step3/se.vec
    p.value.tail1 <- pnorm(-abs(Z.sta.vec))
    p.value.tail2 <- pnorm(-abs(Z.sta.vec)) * 2
    Z.sta.obj <- .three.step.reference.na(
      LTA.vector.to.parameters(
        Z.sta.vec, covariates, L, ref.class, covariates.time.cross
      ), ref.class
    )
    p.value.tail1.obj <- .three.step.reference.na(
      LTA.vector.to.parameters(
        p.value.tail1, covariates, L, ref.class, covariates.time.cross
      ), ref.class
    )
    p.value.tail2.obj <- .three.step.reference.na(
      LTA.vector.to.parameters(
        p.value.tail2, covariates, L, ref.class, covariates.time.cross
      ), ref.class
    )
  }else{
    Z.sta.vec <- NULL
    p.value.tail1 <- NULL
    p.value.tail2 <- NULL
    Z.sta.obj <- NULL
    p.value.tail1.obj <- NULL
    p.value.tail2.obj <- NULL
  }

  beta <- params.LTA.obj$beta
  gamma <- params.LTA.obj$gamma
  beta.se <- se.obj$beta
  gamma.se <- se.obj$gamma
  beta.Z.sta <- Z.sta.obj$beta
  gamma.Z.sta <- Z.sta.obj$gamma
  beta.p.value.tail1 <- p.value.tail1.obj$beta
  gamma.p.value.tail1 <- p.value.tail1.obj$gamma
  beta.p.value.tail2 <- p.value.tail2.obj$beta
  gamma.p.value.tail2 <- p.value.tail2.obj$gamma

  covariates.ncol <- unlist(lapply(covariates, ncol))
  npar <- get.npar.LTA(covariates.ncol, L, covariates.time.cross)

  Log.Lik = -get.Log.Lik.LTA.optim(params.step3, CEP, P.Z.Xns, Zs, covariates, covariates.time.cross, ref.class)
  AIC <- -2 * Log.Lik + 2 * npar
  BIC <- -2 * Log.Lik + npar * log(N)

  res <- list(
    beta = beta, gamma = gamma,
    beta.se = beta.se, gamma.se = gamma.se,
    beta.Z.sta=beta.Z.sta, gamma.Z.sta=gamma.Z.sta,
    beta.p.value.tail1=beta.p.value.tail1, gamma.p.value.tail1=gamma.p.value.tail1,
    beta.p.value.tail2=beta.p.value.tail2, gamma.p.value.tail2=gamma.p.value.tail2,
    vcov = vcov.step3,
    information = information.step3,
    SE.diagnostics = SE.diagnostics,
    bound.diagnostics = bound.diagnostics,
    npar = npar,
    Log.Lik = Log.Lik,
    AIC = AIC,
    BIC = BIC,
    P.Z.Xns=P.Z.Xns,
    P.Zs=P.Zs,
    Zs=Zs,
    Log.Lik.history = Log.Lik.history,
    iterations = optimization.result.main$iterations,
    converged = optimization.result.main$status %in% 1:4,
    params = params,
    call = call,
    arguments = list(
      responses=responses, L=L, type.model=type.model,
      covariates=covariates,
      CEP.time.cross =CEP.time.cross,
      CEP.error=CEP.error,
      covariates.time.cross =covariates.time.cross,
      par.ini = par.ini,
      params=params.input,
      step1.pool=step1.pool,
      constraint = constraint,
      method.model=method.model, tol=tol,
      method.regression=method.regression,
      lower=lower, upper=upper,
      method.SE=method.SE, nrep.bootstrap=nrep.bootstrap,
      maxiter=maxiter, is.sort=is.sort,
      nrep = nrep, starts=starts, maxiter.warmup=maxiter.warmup,
      vis = vis,
      control.EM=control.EM,
      control.Mplus=control.Mplus,
      control.NNE=control.NNE,
      control.Rmixmod=control.Rmixmod,
      ref.class = ref.class
    )
  )

  class(res) <- "LTA"
  return(res)
}

.make.latent.paths <- function(L, times) {
  make_latent_paths_cpp(L, times)
}

LTA.vector.to.parameters <- function(params, covariates, L, ref.class,
                                     covariates.time.cross = FALSE) {
  if (ref.class < 1 || ref.class > L) {
    stop("ref.class must be between 1 and L")
  }

  params.expanded <- .three.step.expand.parameters(
    params, covariates, L, ref.class, covariates.time.cross
  )

  lta_vector_to_parameters_cpp(
    params = params.expanded,
    covariates_list = covariates,
    L = L,
    ref_class = ref.class
  )
}

get.Log.Lik.LTA.optim <- function(params, CEP, P.Z.Xns, Zs, covariates,
                                  covariates.time.cross = FALSE, ref.class,
                                  compute.gradient = FALSE) {
  if (length(CEP) != length(covariates) || length(Zs) != length(covariates)) {
    stop("CEP, Zs, and covariates must have the same length")
  }

  L <- ncol(P.Z.Xns[[1]])
  parameter.index <- .three.step.parameter.index(
    covariates, L, ref.class, covariates.time.cross
  )
  params.expanded <- unname(params[parameter.index])

  result <- get_log_lik_lta_optim_cpp(
    params = params.expanded,
    CEP_list = CEP,
    Zs_list = Zs,
    covariates_list = covariates,
    covariates_time_cross = covariates.time.cross,
    ref_class = ref.class,
    compute_gradient = compute.gradient
  )

  if (!is.finite(result$objective)) {
    if(compute.gradient){
      return(list(
        objective = 1e10,
        gradient = rep(0, length(params))
      ))
    }
    return(1e10)
  }

  if(!compute.gradient){
    return(result$objective)
  }
  gradient <- rowsum(
    matrix(result$gradient, ncol = 1L), parameter.index,
    reorder = FALSE
  )[, 1L]
  list(objective = result$objective, gradient = unname(gradient))
}
