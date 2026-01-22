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
#'                  For \code{type = "LCA"}: items must be binary or categorical (coded as integers starting from 0).
#'                  For \code{type = "LPA"}: items must be continuous (numeric), and each response matrix must be
#'                  standardized using \code{\link[base]{scale}} or \code{\link[LCPA]{normalize}} prior to input.
#' @param L Integer scalar. Number of latent classes/profiles at each time point. Must satisfy \eqn{L \geq 2}.
#' @param ref.class Integer \eqn{L \geq ref.class \geq 1}. Specifies which latent class to use as the reference category.
#'                  Default is \code{L} (last class). Coefficients for the reference class are fixed to zero.
#'                  When \code{is.sort=TRUE}, classes are first ordered by decreasing \code{P.Z} (class 1 has highest probability),
#'                  then \code{ref.class} refers to the position in this sorted order.
#' @param type Character string. Specifies the type of latent variable model for Step 1:
#'             \itemize{
#'               \item \code{"LCA"} — Latent Class Analysis for categorical items.
#'               \item \code{"LPA"} — Latent Profile Analysis for continuous items.
#'             }
#'             See \code{\link[LCPA]{LCA}} and \code{\link[LCPA]{LPA}} for details.
#' @param covariates Optional. A \code{list} of matrices/data frames (length = number of time points).
#'                   Each matrix contains covariates for modeling transitions or initial status.
#'                   Must include an intercept column (all 1s) as the first column.
#'                   If \code{NULL} (default), only intercept terms are used (i.e., no covariates).
#'                   For time \eqn{t}, dimension is \eqn{N \times p_t}. Covariates can vary across time.
#' @param CEP.timeCross Logical. If \code{TRUE}, assumes measurement invariance and uses the same
#'                      Classification Error Probability (\code{\link[LCPA]{get.CEP}}) matrix across all time points.
#'                      Requires that item parameters are invariant over time (not checked internally).
#'                      Default is \code{FALSE}.
#' @param CEP.error Logical. If \code{TRUE} (recommended), incorporates classification uncertainty via
#'                  estimated CEP matrices from Step 1. If \code{FALSE}, uses identity CEP matrices
#'                  (equivalent to naive modal assignment; introduces bias and not recommended).
#' @param covariates.timeCross Logical. If \code{TRUE}, forces the use of identical \eqn{\gamma} parameters
#'                             across all time points (i.e., a time-invariant probability transition matrix).
#'                             In this case, users should ensure that the covariate matrices at different time points
#'                             have the same dimensions (values may differ) to match the fixed form of the \eqn{\gamma_{lkt}} coefficients.
#'                             Default is \code{FALSE}, allowing for potentially different probability transition matrices across time points.
#' @param par.ini Specification for parameter initialization. Options include:
#'   \itemize{
#'     \item \code{"random"}: Completely random initialization (default).
#'     \item \code{"kmeans"}: Initializes parameters via K-means clustering on observed data (McLachlan & Peel, 2004).
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
#' @param method.SE Character. Method for estimating standard errors of parameter estimates:
#'               \itemize{
#'                 \item \code{"Obs"} — Approximates the observed information matrix via numerical differentiation (Richardson's method).
#'                       Standard errors are obtained from the inverse Hessian. May fail or be unreliable in small samples or with complex likelihood surfaces.
#'                 \item \code{"Bootstrap"} — Uses nonparametric bootstrap resampling to estimate empirical sampling variability.
#'                       More robust to model misspecification and small-sample bias. Computationally intensive but recommended when asymptotic assumptions are questionable.
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
#'     \item{\code{hidden.layers}}{Integer vector specifying layer sizes in fully-connected network (default: \code{c(12,12)}).}
#'     \item{\code{activation.function}}{Activation function (e.g., \code{"tanh"}, default: \code{"tanh"}).}
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
#' @return An object of class \code{LTA}, a named list containing:
#' \describe{
#'   \item{\code{beta}}{Matrix of size \eqn{p_1 \times L}. Coefficients for initial class membership multinomial logit model.
#'     Columns 1 to \eqn{L-1} are free parameters; column \eqn{L} (reference class) is constrained to \eqn{\boldsymbol{\beta}_L = \mathbf{0}}.}
#'   \item{\code{gamma}}{List of length \eqn{T-1}. Each element \code{gamma[[t]]} (for transition from time \eqn{t} to \eqn{t+1})
#'     is a nested list: \code{gamma[[t]][[from_class]][[to_class]]} returns coefficient vector of length \eqn{p_{t+1}}.
#'     The last class (\eqn{L}) is reference → coefficients fixed to \eqn{\boldsymbol{\gamma}_{kl,t+1} = \mathbf{0}} for all \eqn{k}.}
#'   \item{\code{beta.se}}{Standard errors for \code{beta} (if Hessian is invertible). Same dimensions as \code{beta}.
#'     May contain \code{NA} if variance-covariance matrix is not positive definite.}
#'   \item{\code{gamma.se}}{Standard errors for \code{gamma}, same nested structure. May contain \code{NA}s.}
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
#'   \item{\code{P.Z.Xns}}{List of length \eqn{T}. Each element is an \eqn{N \times L} matrix of posterior class probabilities
#'     \eqn{P(Z_{nt}=l \mid \mathbf{X}_{nt})} for each individual \eqn{n} at time \eqn{t}.}
#'   \item{\code{P.Zs}}{List of length \eqn{T}. Each element is a vector of length \eqn{L} containing prior class proportions
#'     \eqn{P(Z_t = l)} estimated at Step 1 for time \eqn{t}.}
#'   \item{\code{Zs}}{List of length \eqn{T}. Each element is a vector of length \eqn{N} containing modal class assignments
#'     (MAP classifications) \eqn{\hat{z}_{nt}} for each individual at time \eqn{t}.}
#'   \item{\code{npar}}{Number of free parameters in the model (depends on \code{covariates}).}
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
#' The three-step LTA proceeds as follows:
#'
#' Step 1 — Unconditional Latent Class/Profile Model:
#' At each time point \eqn{t}, fit an unconditional LCA or LPA model (ignoring transitions and covariates).
#' Obtain posterior class membership probabilities \eqn{P(Z_{nt}=l \mid \mathbf{X}_{nt})} for each individual \eqn{n}
#' and class \eqn{l} using Bayes' theorem.
#'
#' Step 2 — Classification Error Probabilities (equal to \code{\link[LCPA]{get.CEP}}):
#' Compute the \eqn{L \times L} CEP matrix for each time point \eqn{t}, where element \eqn{(k,l)} estimates:
#' \deqn{
#'   \text{CEP}_t(k,l) = P(\hat{Z}_{nt} = l \mid Z_{nt} = k)
#' }
#' using a non-parametric approximation based on posterior weights:
#' \deqn{
#'   \widehat{\text{CEP}}_t(k,l) = \frac{ \sum_{n=1}^N \mathbb{I}(\hat{z}_{nt} = l) \cdot P(Z_{nt}=k \mid \mathbf{X}_{nt}) }{ \sum_{n=1}^N P(Z_{nt}=k \mid \mathbf{X}_{nt}) }
#' }
#' where \eqn{\hat{z}_{nt}} is the modal (most likely) class assignment for individual \eqn{n} at time \eqn{t}.
#'
#' Step 3 — Transition Model with Measurement Error Correction:
#' Estimate the multinomial logit models for:
#' \itemize{
#'   \item Initial class membership (time 1): \eqn{P(Z_{n1} = l \mid \mathbf{X}_{n1}) = \frac{\exp(\boldsymbol{\beta}_l^\top \mathbf{X}_{n1})}{\sum_{k=1}^L \exp(\boldsymbol{\beta}_k^\top \mathbf{X}_{n1})}}
#'   \item Transitions (time \eqn{t > 1}): \eqn{P(Z_{nt} = l \mid Z_{n,t-1} = k, \mathbf{X}_{nt}) = \frac{\exp(\boldsymbol{\gamma}_{kl t}^\top \mathbf{X}_{nt})}{\sum_{j=1}^L \exp(\boldsymbol{\gamma}_{kj t}^\top \mathbf{X}_{nt})}}
#' }
#' where \eqn{\mathbf{X}_{n1} = (X_{n10}, X_{n11}, \dots, X_{n1M})^\top} is the covariate vector for observation/participant \eqn{n} at time 1,
#' with \eqn{X_{n10} = 1} (intercept term) and \eqn{X_{n1m}} (\eqn{m=1,\dots,M}) representing the value of the \eqn{m}-th covariate.
#' The coefficient vector \eqn{\boldsymbol{\beta}_l = (\beta_{l0}, \beta_{l1}, \dots, \beta_{lM})^\top} corresponds element-wise to \eqn{\mathbf{X}_{n1}},
#' where \eqn{\beta_{l0}} is the intercept and \eqn{\beta_{lm}} (\eqn{m \geq 1}) are regression coefficients for covariates.
#' Class \eqn{L} is the reference class (\eqn{\boldsymbol{\beta}_L = \mathbf{0}}).
#' \eqn{\mathbf{X}_{nt} = (X_{nt0}, X_{nt1}, \dots, X_{ntM})^\top} is the covariate vector at time \eqn{t},
#' with \eqn{X_{nt0} = 1} (intercept) and \eqn{X_{ntm}} (\eqn{m=1,\dots,M}) as the \eqn{m}-th covariate value.
#' The coefficient vector \eqn{\boldsymbol{\gamma}_{lkt} = (\gamma_{lkt0}, \gamma_{lkt1}, \dots, \gamma_{lktM})^\top}
#' corresponds element-wise to \eqn{\mathbf{X}_{nt}}, where \eqn{\gamma_{lkt0}} is the intercept and \eqn{\gamma_{lktm}} (\eqn{m \geq 1})
#' are regression coefficients. Class \eqn{L} is the reference class (\eqn{\boldsymbol{\gamma}_{lLt} = \mathbf{0}} for all \eqn{l}).
#'
#' The full observed-data likelihood integrates over all possible latent class paths \eqn{\mathbf{z}_n = (z_{n1},\dots,z_{nT})}:
#' \deqn{
#' \begin{aligned}
#' \log \mathcal{L}(\boldsymbol{\theta}) &=
#' \sum_{n=1}^N \log \Biggl[
#'   \sum_{\mathbf{z}_n \in \{1,\dots,L\}^T}
#'   \Bigl( \prod_{t=1}^T \text{CEP}_t(z_{nt}, \hat{z}_{nt}) \Bigr) \cdot \\
#'   &\quad P(Z_{n1}=z_{n1} \mid \mathbf{X}_{n1}) \cdot \\
#'   &\quad \prod_{t=2}^T P(Z_{nt}=z_{nt} \mid Z_{n,t-1}=z_{n,t-1}, \mathbf{X}_{nt})
#' \Biggr]
#' \end{aligned}
#' }
#' Parameters \eqn{\boldsymbol{\theta} = \{\boldsymbol{\beta}, \boldsymbol{\gamma}\}} are estimated via maximum likelihood using
#' the BOBYQA algorithm (box-constrained derivative-free optimization). Reference class \eqn{L} satisfies \eqn{\boldsymbol{\beta}_L = \mathbf{0}} and \eqn{\boldsymbol{\gamma}_{kl t} = \mathbf{0}} for all \eqn{k} when \eqn{l = L}.
#'
#' @section Bootstrap Standard Error Estimation:
#' When \code{method.SE = "Bootstrap"}, standard errors are estimated using a nonparametric bootstrap procedure:
#' \enumerate{
#'   \item Draw \eqn{B} (=\code{n.Bootstrap}) independent samples of size \eqn{N} with replacement from the original data.
#'   \item For each bootstrap sample \eqn{b=1,\dots,B}, re-estimate the full three-step LTA model (Steps 1–3), yielding parameter vector \eqn{\hat{\boldsymbol{\theta}}^{(b)}}.
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
#'   \item Reference Class: The last latent class (\eqn{L}) is always treated as the reference category.
#'         All corresponding coefficients in \code{beta} and \code{gamma} are fixed to zero (\eqn{\boldsymbol{\beta}_L = \mathbf{0}}, \eqn{\boldsymbol{\gamma}_{kl t} = \mathbf{0}} for \eqn{l=L}).
#'   \item CEP Matrices: When \code{CEP.error = TRUE}, misclassification probabilities are estimated
#'                       non-parametrically using Step 1 posterior probabilities. This corrects for
#'                       classification uncertainty. Setting \code{CEP.timeCross = TRUE} assumes these
#'                       error structures are identical across time (measurement invariance).
#'                       See in \code{\link[LCPA]{get.CEP}}.
#'   \item Covariate Handling: Covariates for initial status (time 1) and transitions (time \eqn{t \geq 2}) can differ.
#'         For transitions to time \eqn{t}, the covariate matrix must have dimensions \eqn{N \times (M_{t}+1)},
#'         i.e., an intercept column of all \eqn{1} plus \eqn{M_{t}} columns of covariates in time \eqn{t}.
#'   \item Optimization: Step 3 uses L-BFGS-B via \code{\link[nloptr]{nloptr}} to ensure numerical stability.
#'         Standard errors are derived from the inverse Hessian (via \code{\link[numDeriv]{hessian}}). If singular,
#'         Moore-Penrose pseudoinverse (\code{\link[MASS]{ginv}}) is used, and negative variances are set to \code{NA}.
#'   \item Computational Complexity: Likelihood evaluation requires enumerating \eqn{L^T} possible latent paths.
#'   \item Bootstrap Computation: Each bootstrap iteration re-fits Steps 1–3 independently, including re-estimation of
#'         \code{P.Z.Xns}, \code{CEP} matrices and transition parameters.
#'         To ensure reproducibility, set a seed before calling \code{LTA()} when using \code{method.SE = "Bootstrap"}.
#'         Progress messages during bootstrapping include current replicate index and optimization diagnostics.
#'         Users should monitor convergence in each bootstrap run; failed runs will result in \code{NA} entries in SEs and derived statistics.
#' }
#'
#' @references
#' Liang, Q., la Torre, J. d., & Law, N. (2023). Latent Transition Cognitive Diagnosis Model With Covariates: A Three-Step Approach. Journal of Educational and Behavioral Statistics, 48(6), 690-718. https://doi.org/10.3102/10769986231163320
#'
#' McLachlan, G. J., & Peel, D. (2004). Finite Mixture Models. Wiley. https://books.google.com.sg/books?id=c2_fAox0DQoC
#'
#' Vermunt, J. K. (2010). Latent class modeling with covariates: Two improved three-step approaches. Political Analysis, 18(4), 450–469. https://doi.org/10.1093/pan/mpq025
#'
#' @examples
#'
#' \donttest{
#'
#' library(LCPA)
#'
#' set.seed(123)
#' N <- 2000 ## sample size
#' L <- 3   ## number of latent class
#' I <- 6   ## number of variables/items
#'
#' ## Covariates at time point T1
#' covariates.inter <- rep(1, N)                 # Intercept term is always 1 for each individual
#' covariates.X11 <- rnorm(N)                    # Covariate X1 is a continuous variable
#' # Combine into covariates at T1
#' covariates.T1 <- cbind(Intercept=covariates.inter, X1=covariates.X11)
#' ## Covariates at time point T2
#' covariates.inter <- rep(1, N)                 # Intercept term is always 1 for each individual
#' covariates.X21 <- rnorm(N)                    # Covariate X1 is a continuous variable
#' # Combine into covariates at T1
#' covariates.T2 <- cbind(Intercept=covariates.inter, X1=covariates.X21)
#'
#' # Combine into final covariates list
#' covariates <- list(t1=covariates.T1, t2=covariates.T2)
#'
#' ## Simulate beta coefficients
#' # last column is zero because the last category is used as reference
#' ## fix reference class to class 2
#' beta <- matrix(c( 0.8, 0.0,  0.1,
#'                  -0.3, 0.0, -0.5), ncol=L, byrow=TRUE)
#'
#' ## Simulate gamma coefficients
#' gamma <- list(
#'   lapply(1:L, function(l) {
#'     lapply(1:L, function(k) if(k < L)
#'            runif(2, -2.0, 2.0) else c(0, 0)) # Last class as reference
#'   })
#' )
#'
#' ## Simulate the data
#' sim_custom <- sim.LTA(
#'   N=N, I=I, L=L, times=2, type="LCA", IQ=0.9,
#'   covariates=covariates,
#'   beta=beta,
#'   gamma=gamma
#' )
#' summary(sim_custom)
#' responses <- sim_custom$responses
#' covariates <- sim_custom$covariates
#'
#' ## fix reference class to class 2
#' LTA.obj <- LTA(responses, L=L, ref.class=2, type="LCA",
#'                covariates=covariates,
#'                method.SE="Bootstrap", n.Bootstrap=10,
#'                CEP.timeCross=FALSE, CEP.error=TRUE, covariates.timeCross=FALSE,
#'                par.ini = "random", method="EM", vis = TRUE)
#'
#' print(LTA.obj)
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
LTA <- function(responses, L=2,
                ref.class = L, type="LCA",
                covariates=NULL,
                CEP.timeCross=FALSE,
                CEP.error=TRUE,
                covariates.timeCross=FALSE,
                par.ini = "random",
                params=NULL, is.sort=TRUE,
                constraint = "VV",
                method="EM", tol=1e-4,
                method.SE="Bootstrap", n.Bootstrap=100,
                maxiter=5000, nrep = 20,
                starts=100, maxiter.wa=20,
                vis = TRUE,
                control.EM=NULL,
                control.Mplus=NULL,
                control.NNE=NULL){

  call <- match.call()

  if(ref.class < 1 || ref.class > L) {
    stop("ref.class must be between 1 and L")
  }

  default_control.EM <- list(maxiter=2000, tol=1e-4)
  default_control.Mplus <- list(maxiter=2000, tol=1e-4, files.path = NULL, files.clean = TRUE)
  default_control.NNE <- list(
    hidden.layers = c(12, 12),
    activation.function = "tanh",
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

  response <- as.matrix(responses[[1]])
  N <- nrow(response)
  I <- ncol(response)
  times <- length(responses)

  if(vis){
    cat("Starting the first step for LTA ...\n\n")
  }

  if(is.null(params)){
    if(type == "LCA"){
      adj <- adjust.response(response)
      response <- adj$response
      LCPA.obj <- LCA(response, L=L, par.ini=par.ini,
                      method=method, is.sort=is.sort, nrep=nrep,
                      starts=starts, maxiter.wa=maxiter.wa,
                      vis=vis,
                      control.EM=control.EM,
                      control.Mplus=control.Mplus,
                      control.NNE=control.NNE)
    } else {
      LCPA.obj <- LPA(response, L=L, par.ini=par.ini,
                      constraint=constraint,
                      method=method, is.sort=is.sort, nrep=nrep,
                      starts=starts, maxiter.wa=maxiter.wa,
                      vis=vis,
                      control.EM=control.EM,
                      control.Mplus=control.Mplus,
                      control.NNE=control.NNE)
    }
    params <- LCPA.obj$params
  }

  if(vis){
    cat("Starting the second step for LTA ...\n\n")
  }

  P.Z.Xns <- Zs <- P.Zs <- vector("list", times)
  for(t in 1:times){
    if(type == "LCA"){
      P.Z.Xns[[t]] <- get.P.Z.Xn.LCA(response=responses[[t]], par=params$par, vis=vis)
    } else {
      P.Z.Xns[[t]] <- get.P.Z.Xn.LPA(response=responses[[t]], means=params$means, covs=params$covs, vis=vis)
    }
    P.Zs[[t]] <- colSums(P.Z.Xns[[t]]) / sum(P.Z.Xns[[t]])
    Zs[[t]] <- apply(P.Z.Xns[[t]], 1, which.max)
  }

  if(CEP.error){
    CEP <- get.CEP(P.Z.Xns, time.cross=CEP.timeCross)
  } else {
    CEP <- replicate(times, diag(L), simplify=FALSE)
  }

  if(is.null(covariates)){
    covariates <- lapply(responses, function(x){
      matrix(1, nrow=nrow(x), ncol=1)
    })
  }

  p1 <- ncol(covariates[[1]])
  num_beta <- p1 * (L - 1)
  num_gamma <- sum(sapply(2:times, function(t){
    L * (L - 1) * ncol(covariates[[t]])
  }))
  init.par <- rep(0, num_beta + num_gamma)

  if(vis){
    cat("Starting the third step for LTA ...\n\n")
  }

  npar <- length(init.par)
  lb <- rep(-5, npar)
  ub <- rep( 5, npar)

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
      Log.Lik <- get.Log.Lik.LTA.optim(par, CEP, P.Z.Xns, Zs, covariates, covariates.timeCross, ref.class)

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
    eval_f = function(x) get.Log.Lik.LTA.print(x, CEP, P.Z.Xns, Zs, covariates, covariates.timeCross),
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
            covariates.timeCross,
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
              covariates.timeCross,
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
                covariates.timeCross,
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

      samples.cur <- sample(1:N, N, replace = TRUE)
      P.Z.Xns.cur <- Zs.cur <- P.Zs.cur <- covariates.cur <- vector("list", times)
      for(t in 1:times){
        if(type == "LCA"){
          P.Z.Xns.cur[[t]] <- get.P.Z.Xn.LCA(response=responses[[t]][samples.cur, ], par=params$par, vis=FALSE)
        } else {
          P.Z.Xns.cur[[t]] <- get.P.Z.Xn.LPA(response=responses[[t]][samples.cur, ], means=params$means, covs=params$covs, vis=FALSE)
        }

        covariates.cur[[t]] <- covariates[[t]][samples.cur, , drop=FALSE]
        P.Zs.cur[[t]] <- colSums(P.Z.Xns.cur[[t]]) / sum(P.Z.Xns.cur[[t]])
        Zs.cur[[t]] <- apply(P.Z.Xns.cur[[t]], 1, which.max)
      }

      if(CEP.error){
        CEP.cur <- get.CEP(P.Z.Xns.cur, time.cross=CEP.timeCross)
      } else {
        CEP.cur <- replicate(times, diag(L), simplify=FALSE)
      }

      init.par <- rnorm(npar, mean=refined_init, sd=abs(refined_init) * 0.1)
      Log.Lik.history.cur <- c(0)
      make_loglik_with_print_Bootstrap <- function(vis, ref.class, bs, n.Bootstrap) {
        iter <- 0
        function(par, CEP.cur, P.Z.Xns.cur, Zs.cur, covariates.cur, covariates.timeCross, bs, n.Bootstrap) {
          iter <<- iter + 1
          Log.Lik <- get.Log.Lik.LTA.optim(par, CEP.cur, P.Z.Xns.cur, Zs.cur, covariates.cur, covariates.timeCross, ref.class)

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
        eval_f = function(x) get.Log.Lik.LTA.print(x, CEP.cur, P.Z.Xns.cur, Zs.cur, covariates.cur, covariates.timeCross, bs, n.Bootstrap),
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
  npar <- get.npar.LTA(covariates.ncol, L, covariates.timeCross)

  Log.Lik = -get.Log.Lik.LTA.optim(refined_init, CEP, P.Z.Xns, Zs, covariates, covariates.timeCross, ref.class)
  AIC <- -2 * Log.Lik + 2 * npar
  BIC <- -2 * Log.Lik + npar * log(N)

  res <- list(
    beta = beta, gamma = gamma,
    beta.se = beta.se, gamma.se = gamma.se,
    beta.Z.sta=beta.Z.sta, gamma.Z.sta=gamma.Z.sta,
    beta.p.value.tail1=beta.p.value.tail1, gamma.p.value.tail1=gamma.p.value.tail1,
    beta.p.value.tail2=beta.p.value.tail2, gamma.p.value.tail2=gamma.p.value.tail2,
    npar = npar,
    Log.Lik = Log.Lik,
    AIC = AIC,
    BIC = BIC,
    P.Z.Xns=P.Z.Xns,
    P.Zs=P.Zs,
    Zs=Zs,
    Log.Lik.history = Log.Lik.history[1:nlopt_res$iterations],
    iterations = nlopt_res$iterations,
    coveraged = nlopt_res$iterations < maxiter,
    params = params,
    call = call,
    arguments = list(
      responses=responses, L=L, type=type,
      covariates=covariates,
      CEP.timeCross=CEP.timeCross,
      CEP.error=CEP.error,
      covariates.timeCross=covariates.timeCross,
      par.ini = par.ini,
      params=params,
      constraint = constraint,
      method=method, tol=tol,
      method.SE=method.SE, n.Bootstrap=n.Bootstrap,
      maxiter=maxiter, is.sort=is.sort,
      nrep = nrep, starts=starts, maxiter.wa=maxiter.wa,
      vis = vis,
      control.EM=control.EM,
      control.Mplus=control.Mplus,
      control.NNE=control.NNE,
      ref.class = ref.class
    )
  )

  class(res) <- "LTA"
  return(res)
}

.make.latent.paths <- function(L, times) {
  make_latent_paths_cpp(L, times)
}

LTA.vector.to.parameters <- function(params, covariates, L, ref.class) {
  if (ref.class < 1 || ref.class > L) {
    stop("ref.class must be between 1 and L")
  }

  lta_vector_to_parameters_cpp(
    params = params,
    covariates_list = covariates,
    L = L,
    ref_class = ref.class
  )
}

get.Log.Lik.LTA.optim <- function(init.par, CEP, P.Z.Xns, Zs, covariates,
                                  covariates.timeCross = TRUE, ref.class) {
  if (length(CEP) != length(covariates) || length(Zs) != length(covariates)) {
    stop("CEP, Zs, and covariates must have the same length")
  }

  latent.paths <- .make.latent.paths(L = ncol(P.Z.Xns[[1]]), times = length(covariates))

  result <- get_log_lik_lta_optim_cpp(
    init_par = init.par,
    CEP_list = CEP,
    Zs_list = Zs,
    covariates_list = covariates,
    covariates_timeCross = covariates.timeCross,
    ref_class = ref.class,
    latent_paths = latent.paths
  )

  if (!is.finite(result)) {
    warning("Non-finite likelihood returned from C++")
    return(1e10)
  }

  result
}
