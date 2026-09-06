#' Fit Latent Profile Analysis
#'
#' This function estimates parameters of a Latent Profile Analysis (LPA) model for continuous observed variables
#' using the Expectation-Maximization (EM) algorithm, stochastic EM (SEM) through \code{flexmix}
#' or \code{RMixtComp}, native EM/CEM/SEM strategies through \code{Rmixmod},
#' Neural Network Estimation (NNE), or external Mplus software.
#'
#' @param response A numeric matrix of dimension \eqn{N \times I}, where \eqn{N} is the number of participants
#'                 and \eqn{I} is the number of continuous observed indicators. Missing values are not allowed.
#'                 Note that \code{response} must be standardized using \code{\link[base]{scale}} or
#'                 \code{\link[LCPA]{normalize}} before input.
#' @param L Integer specifying the number of latent profiles (default: 2).
#' @param constraint Character string specifying covariance structure constraints:
#'   \describe{
#'     \item{\code{"VV"}}{Varying variances and varying covariances across profiles (heterogeneous full covariance; Default).}
#'     \item{\code{"VE"}}{Varying variances but equal covariances across profiles.}
#'     \item{\code{"EV"}}{Equal variances but varying covariances across profiles.}
#'     \item{\code{"EE"}}{Equal variances and equal covariances across profiles (homogeneous full covariance).}
#'     \item{\code{"E0"}}{Equal variances across profiles, zero covariances (diagonal with shared variances).}
#'     \item{\code{"V0"}}{Varying variances across profiles, zero covariances (diagonal with free variances).}
#'     \item{\code{"UE"}}{Univariate response only: equal variance across profiles.}
#'     \item{\code{"UV"}}{Univariate response only: varying variances across profiles.}
#'     \item{\code{list}}{Custom constraints. Each element is a 2-element integer vector specifying variables whose covariance
#'                        parameters are constrained equal across all classes. The constraint applies to:
#'                        \itemize{
#'                          \item Variances: When both indices are identical (e.g., \code{c(3,3)} forces variance of variable 3 to be equal across classes).
#'                          \item Covariances: When indices differ (e.g., \code{c(1,2)} forces covariance between variables 1 and 2 to be equal across classes).
#'                        }
#'                        Constraints are symmetric (e.g., \code{c(1,2)} automatically constrains \code{c(2,1)}). All unconstrained parameters
#'                        vary freely across classes while maintaining positive definiteness.
#'                        }
#'   }
#'   For \code{method = "flexmix"}, all eight named structures and custom constraint lists are supported
#'   through LCPA's joint Gaussian M-step driver inside the flexmix SEM loop.
#'   For \code{method = "Rmixmod"}, only \code{"V0"}, \code{"EE"}, and \code{"VV"} are supported;
#'   custom constraint lists are not supported.
#'   For \code{method = "RMixtComp"}, only the locally independent \code{"V0"} structure is supported;
#'   custom constraint lists are not supported.
#'   For \code{method = "Mplus"}, at least two indicators are required; the supported named structures are
#'   \code{"E0"}, \code{"V0"}, \code{"EE"}, \code{"VE"}, \code{"EV"}, and \code{"VV"}, together with
#'   custom constraint lists.
#' @param method Character string specifying estimation algorithm:
#'   \itemize{
#'     \item \code{"EM"}: Expectation-Maximization algorithm (Default).
#'     \item \code{"NNE"}: Neural Network Estimation (experimental), using feed-forward layers,
#'                         optional transformer attention, gradient optimization, and simulated annealing.
#'                         See \code{\link[LCPA]{install_python_dependencies}}.
#'     \item \code{"Mplus"}: Calls external Mplus software for estimation.
#'                           Uses Mplus defaults for optimization unless overridden by \code{control.Mplus}.
#'     \item \code{"flexmix"}: Stochastic EM (SEM) through \code{flexmix}, using LCPA's
#'                              warm-up and promoted-replication scheme. Requires \code{flexmix}.
#'     \item \code{"Rmixmod"}: LCPA warm-up plus SEM, or a native Rmixmod EM, CEM, or SEM strategy.
#'                              Requires the \code{Rmixmod} package.
#'     \item \code{"RMixtComp"}: Stochastic EM (SEM) estimation through \code{RMixtComp}.
#'                                Requires \code{RMixtComp} and \code{RMixtCompUtilities}.
#'                                No non-SEM RMixtComp algorithm is exposed.
#'   }
#' @param par.ini Specification for parameter initialization. Options include:
#'   \itemize{
#'     \item \code{"random"}: Random initialization of means and covariances (default).
#'     \item \code{"kmeans"}: Initializes parameters via K-means clustering on observed data (McLachlan & Peel, 2000).
#'     \item A \code{list} containing exactly three elements:
#'       \describe{
#'         \item{\code{means}}{An \eqn{L \times I} matrix of initial mean vectors for each profile.}
#'         \item{\code{covs}}{An \eqn{I \times I \times L} array of initial covariance matrices for each profile.}
#'         \item{\code{P.Z}}{A numeric vector of length \eqn{L} specifying initial prior probabilities for profiles.}
#'       }
#'   }
#'   For methods that expose this initialization interface, \code{par.ini} is used only
#'   to construct the \code{starts} warm-up initializations. It does not initialize the
#'   \code{nrep} refinement runs, which continue directly from the selected warm-up states.
#'   With \code{par.ini = "kmeans"}, each EM warm-up start calls
#'   \code{\link[LCPA]{Kmeans.LPA}} exactly once with one internal K-means start; NNE applies
#'   the same single-start K-means rule in its Python backend. Thus, \code{starts} is the number
#'   of separately initialized K-means outputs passed to warm-up training. Backends such as
#'   Mplus and Rmixmod may instead use their
#'   native initialization mechanisms and are not required to implement K-means initialization.
#'   If \code{"kmeans"} is requested for a method that does not support it, \code{par.ini}
#'   is automatically changed to \code{"random"}.
#' @param is.sort A logical value. If \code{TRUE} (Default), the latent classes will be ordered in descending
#'                order according to \code{P.Z}. All other parameters will be adjusted accordingly
#'                based on the reordered latent classes.
#' @param starts Positive integer. Number of warm-up analyses to run (default: 100).
#'   Each analysis is initialized by the selected method and trained for at most
#'   \code{maxiter.warmup} iterations, producing exactly \code{starts} warm-up solutions.
#'   With \code{par.ini = "kmeans"}, each EM warm-up start calls
#'   \code{\link[LCPA]{Kmeans.LPA}} exactly once with \code{starts = 1}; NNE applies the
#'   equivalent single-start rule in Python. The public \code{starts} argument is therefore
#'   also the number of K-means outputs.
#' @param maxiter.warmup Positive integer. Maximum number of training iterations for
#'   each of the \code{starts} warm-up analyses (default: 20). This limit applies only
#'   to warm-up and does not limit the subsequent refinement phase.
#' @param nrep Positive integer not exceeding \code{starts}. Number of refinement
#'   analyses (default: 20). The \code{nrep} warm-up solutions with the largest
#'   log-likelihoods are continued from their saved states until the full-training
#'   stopping rule is met; no new initialization occurs in this phase. The refined
#'   solution with the largest log-likelihood is returned as the final result.
#'   These three staged-training arguments are not used by \code{method = "RMixtComp"}
#'   or by \code{method = "Rmixmod"} with \code{control.Rmixmod$path = "Rmixmod"};
#'   those paths use their documented native controls instead.
#' @param vis Logical. If \code{TRUE}, displays carriage-return-updated \code{Warm} and
#'   \code{Rep} lines when the backend exposes those stages, followed by one final fit-summary
#'   line (default: \code{TRUE}). Each stage occupies one console line and ends with one newline.
#' @param control.EM List of control parameters for EM algorithm:
#'   \describe{
#'     \item{\code{maxiter}}{Maximum iterations (default: 2000).}
#'     \item{\code{tol}}{Convergence tolerance for log-likelihood difference (default: 1e-4).}
#'   }
#' @param control.Mplus List of control parameters for Mplus estimation:
#'   \describe{
#'     \item{\code{maxiter}}{Maximum iterations for Mplus optimization (default: 2000).}
#'     \item{\code{tol}}{Convergence tolerance for log-likelihood difference (default: 1e-4).}
#'     \item{\code{files.path}}{A character string specifying the directory under which Mplus writes
#'       intermediate files, including model input, data, output, and saved posterior probabilities.
#'       The effective default is \code{""}. A non-empty path is created recursively when necessary
#'       and must be writable. Within it, the function creates a unique timestamped subdirectory named
#'       \code{"Mplus_LPA_YYYY-MM-DD_HH-MM-SS"} to isolate all files from the current run.
#'       If \code{files.path = ""}, that timestamped subdirectory is created directly under R's current
#'       working directory, \code{\link[base]{getwd}()}. Explicit \code{NULL} is invalid.}
#'     \item{\code{files.clean}}{Logical. If \code{TRUE} (default), all intermediate files and the temporary working directory
#'       created for the run are deleted on successful completion or error exit via \code{on.exit()}.
#'       If \code{FALSE}, the complete timestamped working directory is retained under \code{files.path},
#'       or under \code{\link[base]{getwd}()} when \code{files.path = ""}, for inspection and debugging.}
#'   }
#' @param control.NNE List of control parameters for NNE algorithm:
#'   \describe{
#'     \item{\code{hidden.layers}}{Integer vector specifying layer sizes in fully-connected network (default: \code{c(16,16)}).}
#'     \item{\code{activation.function}}{Activation function (e.g., \code{"tanh"}, default: \code{"tanh"}).}
#'     \item{\code{use.attention}}{Whether to enable the self-attention mechanism (i.e., transformer encoder) (default: \code{TRUE}).}
#'     \item{\code{d.model}}{Dimensionality of transformer encoder embeddings (default: 8).}
#'     \item{\code{nhead}}{Number of attention heads in transformer (default: 2).}
#'     \item{\code{dim.feedforward}}{Dimensionality of transformer feedforward network (default: 16).}
#'     \item{\code{eps}}{Positive offset used in the NNE objective (default: 1e-8).}
#'     \item{\code{lambda}}{Coefficient of the NNE parameter penalty (default: 1e-5).}
#'     \item{\code{initial.temperature}}{Initial temperature for simulated annealing (default: 1000).}
#'     \item{\code{cooling.rate}}{Cooling rate per iteration in simulated annealing (default: 0.5).}
#'     \item{\code{maxiter.sa}}{Maximum iterations for simulated annealing (default: 1000).}
#'     \item{\code{threshold.sa}}{Minimum temperature threshold for annealing (default: 1e-10).}
#'     \item{\code{maxiter}}{Maximum training epochs (default: 1000).}
#'     \item{\code{patience.early}}{Maximum consecutive iterations without improvement before early stopping (default: 100).}
#'     \item{\code{maxcycle}}{Maximum cycles for optimization (default: 10).}
#'     \item{\code{lr}}{Learning rate, controlling the step size of neural network parameter updates (default: 0.025).}
#'     \item{\code{scheduler.patience}}{Patience for learning rate decay (if the loss function does not improve for more than `patience` consecutive epochs, the learning rate will be reduced) (default: 10).}
#'     \item{\code{scheduler.factor}}{Learning rate decay factor; the new learning rate equals the original learning rate multiplied by \code{scheduler.factor} (default: 0.80).}
#'     \item{\code{plot.interval}}{Interval (in epochs) for plotting training diagnostics (default: 100).}
#'     \item{\code{device}}{Specifies the hardware device; can be \code{"CPU"} (default) or \code{"GPU"}. If the GPU is not available, it automatically falls back to CPU.}
#'   }
#' @param control.flexmix List of control parameters for flexmix SEM estimation:
#'   \describe{
#'     \item{\code{maxiter}}{Number of SEM iterations in every promoted full run (default: 1000).}
#'     \item{\code{minprior}}{Minimum component prior accepted by flexmix (default: 0, so LCPA does not intentionally remove requested profiles).}
#'     \item{\code{tol}}{Relative likelihood-change threshold used by flexmix (default: 0, which enforces the fixed SEM iteration count).}
#'   }
#' @param control.Rmixmod List of control parameters for Rmixmod estimation:
#'   \describe{
#'     \item{\code{path}}{Execution path: \code{"LCPA"} (default) preserves the package's
#'       \code{starts}/\code{maxiter.warmup}/\code{nrep} warm-up plus pure SEM procedure; \code{"Rmixmod"}
#'       delegates one native EM, CEM, SEM, or combined strategy and ignores those three arguments.}
#'     \item{\code{algorithm}, \code{nrep}, \code{method.init}, \code{starts},
#'       \code{maxiter.init}, \code{maxiter}, \code{tol.init}, \code{tol},
#'       \code{par.ini}, \code{labels.ini}}{Package-standard controls translated to the corresponding
#'       arguments of \code{Rmixmod::mixmodStrategy()} when \code{path="Rmixmod"}. Unspecified arguments
#'       retain the installed Rmixmod version's defaults. Under \code{path="LCPA"}, only
#'       \code{maxiter} is used, with 1000 iterations when omitted. \code{algorithm} accepts
#'       \code{"EM"}, \code{"CEM"}, and \code{"SEM"}, including ordered combinations supported by Rmixmod.}
#'     \item{\code{strategy}}{Optional pre-built Rmixmod \code{Strategy} object for
#'       \code{path="Rmixmod"}. When supplied, it takes precedence over the individual strategy controls.}
#'   }
#' @param control.RMixtComp List of control parameters for RMixtComp SEM estimation:
#'   \describe{
#'     \item{\code{maxiter.burnin}}{Number of native SEM burn-in iterations (default: 50).}
#'     \item{\code{maxiter}}{Number of recorded post-burn-in SEM iterations (default: 50).}
#'     \item{\code{maxiter.gibbs.burnin}, \code{maxiter.gibbs}}{Numbers of burn-in and
#'       recorded iterations in RMixtComp's subsequent fixed-parameter Gibbs stage
#'       (defaults: 50 and 50).}
#'     \item{\code{n.init.per.class}}{Number of observations per class used by
#'       RMixtComp's native parameter initialization (default: 50).}
#'     \item{\code{maxattempts.sem}}{Maximum number of SEM attempts (default: 20).}
#'     \item{\code{confidence.level}, \code{stable.ratio}, \code{n.stable}}{
#'       Native RMixtComp SEM controls (defaults: 0.95, 0.99, and 20).}
#'     \item{\code{criterion}}{RMixtComp model-selection criterion: \code{"BIC"}
#'       (default) or \code{"ICL"}.}
#'     \item{\code{nrep}}{Number of native RMixtComp SEM runs for the requested number
#'       of profiles; RMixtComp retains the run with the largest observed likelihood (default: 1).}
#'     \item{\code{ncores}}{Number of cores used by RMixtComp to parallelize
#'       \code{nrep}; must not exceed \code{nrep} (default: 1).}
#'   }
#'
#' @section Random-number reproducibility:
#' Except for \code{method = "NNE"}, which intentionally uses its fixed backend seed,
#' every stochastic estimator is driven from R's current random-number generator.
#' The user only needs to call \code{set.seed()} immediately before \code{\link[LCPA]{LPA}()} to
#' reproduce EM, K-means, flexmix, Rmixmod, RMixtComp, and Mplus estimation. LCPA
#' automatically passes an R-derived seed to backends with independent random streams;
#' no backend-specific seed setting is required.
#'
#' @return An object of class \code{"LPA"} containing:
#'   \describe{
#'     \item{\code{params}}{List with estimated profile parameters:
#'       \describe{
#'         \item{\code{means}}{\eqn{L \times I} matrix of estimated mean vectors for each profile.}
#'         \item{\code{covs}}{\eqn{I \times I \times L} array of estimated covariance matrices for each profile.}
#'         \item{\code{P.Z}}{Vector of length \eqn{L} with profile prior probabilities.}
#'       }
#'     }
#'     \item{\code{npar}}{Number of free parameters in the model (depends on \code{constraint}).}
#'     \item{\code{Log.Lik}}{Log-likelihood of the final model.}
#'     \item{\code{AIC}}{Akaike Information Criterion value.}
#'     \item{\code{BIC}}{Bayesian Information Criterion value.}
#'     \item{\code{best_BIC}}{Best BIC value across \code{nrep} runs when applicable;
#'       for native Rmixmod and RMixtComp paths, the selected native fit's BIC.}
#'     \item{\code{P.Z.Xn}}{\eqn{N \times L} matrix of posterior profile probabilities for each observation.}
#'     \item{\code{P.Z}}{Vector of length \eqn{L} containing the prior probabilities/structural parameters/proportions for each latent class.}
#'     \item{\code{Z}}{Vector of length \eqn{N} with MAP-classified profile memberships.}
#'     \item{\code{Log.Lik.history}}{Vector tracking log-likelihood at each EM iteration (only for \code{method="EM"}).}
#'     \item{\code{Log.Lik.nrep}}{Vector of log-likelihoods from each replication run.
#'       For native Rmixmod and RMixtComp paths, this is the selected native fit's scalar log-likelihood.}
#'     \item{\code{model}}{The optimal model object:
#'       \itemize{
#'         \item For \code{method="NNE"}: Trained neural network model.
#'         \item For \code{method="Mplus"}: Estimated Mplus model.
#'         \item For \code{method="flexmix"}: Selected \code{flexmix} SEM object.
#'         \item For \code{method="Rmixmod"}: Selected \code{MixmodCluster} object.
#'         \item For \code{method="RMixtComp"}: Selected \code{MixtCompLearn} SEM object.
#'       }
#'     }
#'     \item{\code{call}}{Matched function call.}
#'     \item{\code{arguments}}{A list containing all effective input arguments.}
#'   }
#'
#' @section Notation:
#' Write the response matrix as
#' \eqn{\mathbf{X}=(X_{ni})_{N\times I}}, where
#' \eqn{n=1,2,\ldots,N} indexes participants and
#' \eqn{i=1,2,\ldots,I} indexes observed indicators. The response vector for
#' participant \eqn{n} is
#' \eqn{\mathbf{X}_n=(X_{n1},\ldots,X_{nI})^\top}. The latent profile variable
#' is \eqn{Z_n\in\{1,2,\ldots,L\}}, and
#' \eqn{l=1,2,\ldots,L} indexes a particular latent profile.
#'
#' With profile mean \eqn{\boldsymbol{\mu}_l} and covariance
#' \eqn{\boldsymbol{\Sigma}_l}, the observed-data log-likelihood is
#' \deqn{\log\mathcal{L}_{\mathrm{LPA}}=
#' \sum_{n=1}^N\log\left\{\sum_{l=1}^L\pi_l
#' \mathcal{N}(\mathbf{X}_n\mid
#' \boldsymbol{\mu}_l,\boldsymbol{\Sigma}_l)\right\}.}
#'
#' @section EM Algorithm:
#' When \code{method = "EM"}, parameter estimation uses the Expectation-Maximization (EM) algorithm to maximize the observed-data log-likelihood:
#'
#' \deqn{\log\mathcal{L}_{\mathrm{LPA}}=
#' \sum_{n=1}^N\log\left\{\sum_{l=1}^L\pi_l
#' \mathcal{N}(\mathbf{X}_n\mid
#' \boldsymbol{\mu}_l,\boldsymbol{\Sigma}_l)\right\}.}
#'
#' The algorithm iterates between two steps until convergence (change in log-likelihood < \code{tol} or max iterations reached):
#'
#' \describe{
#'   \item{E-step:}{
#'     Compute posterior class probabilities (responsibilities) for participant
#'     \eqn{n} and class \eqn{l}:
#'     \deqn{\tau_{nl}=P(Z_n=l\mid\mathbf{X}_n)=
#'     \frac{\pi_l\mathcal{N}(\mathbf{X}_n\mid
#'     \boldsymbol{\mu}_l,\boldsymbol{\Sigma}_l)}
#'     {\sum_{h=1}^L\pi_h\mathcal{N}(\mathbf{X}_n\mid
#'     \boldsymbol{\mu}_h,\boldsymbol{\Sigma}_h)}.}
#'     where \eqn{\mathcal{N}(\cdot)} is the multivariate normal density, \eqn{\pi_l} is the prior class probability, and \eqn{\boldsymbol{\mu}_l},
#'     \eqn{\boldsymbol{\Sigma}_l} are current parameters.
#'   }
#'
#'   \item{M-step:}{
#'     Update parameters using responsibilities \eqn{\tau_{nl}}:
#'     \itemize{
#'       \item Class probabilities: \eqn{\pi_l^{\text{new}} = \frac{1}{N}\sum_{n=1}^N \tau_{nl}}
#'       \item Class means: \eqn{\boldsymbol{\mu}_l^{\text{new}} =
#'         \frac{\sum_{n=1}^N \tau_{nl} \mathbf{X}_n}
#'         {\sum_{n=1}^N \tau_{nl}}}
#'       \item Class covariances: Updated under constraints:
#'         \describe{
#'           \item{\code{"VV"}}{\eqn{\boldsymbol{\Sigma}_l^{\text{new}} =
#'             \frac{\sum_{n=1}^N \tau_{nl}(\mathbf{X}_n-\boldsymbol{\mu}_l^{\text{new}})
#'             (\mathbf{X}_n-\boldsymbol{\mu}_l^{\text{new}})^\top}
#'             {\sum_{n=1}^N \tau_{nl}}}}
#'           \item{\code{"EE"}}{Shared covariance:
#'             \eqn{\boldsymbol{\Sigma}^{\text{new}} =
#'             \frac{\sum_{l=1}^L\sum_{n=1}^N\tau_{nl}
#'             (\mathbf{X}_n-\boldsymbol{\mu}_l^{\text{new}})
#'             (\mathbf{X}_n-\boldsymbol{\mu}_l^{\text{new}})^\top}
#'             {\sum_{l=1}^L\sum_{n=1}^N\tau_{nl}}}}
#'           \item{\code{"VE"} / \code{"EV"}}{Hybrid constraints (e.g., \code{"VE"}: varying variances, equal covariances).
#'                                            The covariance part of the standard EM Q-function is maximized numerically under the exact equality constraints.}
#'           \item{Custom constraints}{User-specified variances/covariances (e.g., \code{list(c(1,2), c(2, 2))}, meaning the
#'                                     covariates of observed variable 1 and observed variable 2 are equal across latent classes,
#'                                     and the variance of observed variable 2 is equal across classes) are estimated by maximizing
#'                                     the covariance part of the standard EM Q-function under the specified equalities.}
#'         }
#'     }
#'   }
#' }
#'
#' @section Neural Network Estimation (NNE):
#' When \code{method = "NNE"}, parameters are estimated using a hybrid neural network architecture
#' combining fully-connected layers with transformer-based attention mechanisms. This approach jointly
#' optimizes profile parameters and posterior probabilities through stochastic optimization with
#' simulated annealing. See \code{\link[LCPA]{install_python_dependencies}}. Key components include:
#'
#' Architecture:
#' \describe{
#'   \item{Input Representation:}{
#'     Continuous observed indicators \eqn{\mathbf{X}_n \in \mathbb{R}^I} are standardized
#'     (mean-centered and scaled to unit variance) during training. No encoding is required.
#'   }
#'   \item{Feature Encoder (Feedforward Network):}{
#'     A multi-layer perceptron with architecture defined by \code{hidden.layers} and \code{activation.function}
#'     maps the continuous input vector into a latent space of dimension \code{d.model}. This layer learns non-linear
#'     feature combinations predictive of latent profile membership.
#'   }
#'   \item{Attention Refiner (Transformer Encoder)}{
#'     A transformer encoder with \code{nhead} attention heads that learns latent class prior probabilities
#'     \eqn{\boldsymbol{\pi} = (\pi_1, \pi_2, \dots, \pi_L)} directly from observed responses.
#'   }
#'   \item{Parameter Head (Means & Covariances):}{
#'     Two separate projection heads branch from the transformer output:
#'     \itemize{
#'       \item Means Head: Linear projection to \eqn{L \times I} matrix \eqn{\boldsymbol{\mu}_l}.
#'       \item Covariance Head: Outputs one covariance matrix
#'         \eqn{\boldsymbol{\Sigma}_l} for each profile.
#'     }
#'   }
#' }
#'
#' Constraint handling:
#' \itemize{
#'   \item Covariance constraints (\code{constraint}) are enforced after activation;
#'     variances or covariances marked for equality are shared across profiles.
#'   \item Custom constraints: e.g., \code{list(c(1,2), c(3,3))}, force equality of specific covariance elements
#'           across profiles, with symmetry (\eqn{\sigma_{12} = \sigma_{21}}) automatically enforced.
#' }
#'
#' @section Mplus:
#' When \code{method = "Mplus"}, estimation is delegated to external Mplus software.
#' The function automates the entire workflow:
#'
#' Workflow:
#' \describe{
#'   \item{Working Directory Setup}{Creates a timestamped \code{"Mplus_LPA_YYYY-MM-DD_HH-MM-SS"}
#'     directory under \code{control.Mplus$files.path}, or under the current working directory when that path is empty, to store:
#'     \itemize{
#'       \item Mplus input syntax (\code{.inp})
#'       \item Data file in Mplus format (\code{.dat})
#'       \item Posterior probabilities output (\code{.dat})
#'     }
#'     Files are automatically deleted after estimation unless \code{control.Mplus$files.clean = FALSE}.
#'   }
#'
#'   \item{Syntax Generation}{Constructs Mplus syntax with:
#'     \itemize{
#'       \item \code{CLASSES = c1(L)} specification for \eqn{L} latent classes
#'       \item \code{ANALYSIS} block with optimization controls:
#'         \describe{
#'           \item{\code{TYPE = mixture}}{Standard mixture modeling setup}
#'           \item{\code{STARTS = starts nrep}}{Random \code{starts} and final stage optimizations}
#'           \item{\code{STSEED}}{Random-start seed drawn from R's current random-number generator}
#'           \item{\code{STITERATIONS = maxiter.warmup}}{max itertions during \code{starts}.}
#'           \item{\code{MITERATIONS = maxiter}}{Maximum EM iterations}
#'           \item{\code{CONVERGENCE = tol}}{Log-likelihood convergence tolerance}
#'         }
#'       \item \code{MODEL} block reflecting the specified \code{constraint} structure
#'     }
#'   }
#'
#'   \item{Execution}{Calls Mplus via \code{MplusAutomation::mplusModeler()}}, which:
#'     \itemize{
#'       \item Writes data to disk in Mplus-compatible format
#'       \item Invokes the Mplus executable (requires valid license)
#'       \item Captures convergence status and output
#'     }
#' }
#'
#' Constraint handling:
#' \itemize{
#'   \item Covariance restrictions are encoded directly in the generated Mplus model syntax before estimation;
#'     they are not imposed by post-estimation averaging.
#'   \item The named structures are translated as follows:
#'     \itemize{
#'       \item \code{"E0"}: variances carry common equality labels across profiles and all covariances are fixed to zero.
#'       \item \code{"V0"}: variances are freely estimated within profile and all covariances are fixed to zero.
#'       \item \code{"EE"}: variances and covariances carry common equality labels across profiles.
#'       \item \code{"VE"}: variances vary across profiles while covariances carry common equality labels.
#'       \item \code{"EV"}: variances carry common equality labels while covariances vary across profiles.
#'       \item \code{"VV"}: all variances and covariances vary across profiles.
#'     }
#'   \item For a custom list, \code{c(i, i)} assigns a common cross-profile label to the variance of
#'     variable \eqn{i}, whereas \code{c(i, j)} assigns one to the covariance between variables
#'     \eqn{i} and \eqn{j}. Unlisted parameters remain profile-specific, and covariance symmetry is
#'     represented by a single Mplus \code{WITH} parameter for each variable pair.
#'   \item \code{"UE"} and \code{"UV"} are not available through this Mplus backend because it requires
#'     at least two indicators. They remain available for the univariate EM and flexmix paths.
#' }
#'
#' @section flexmix Stochastic EM:
#' With \code{method = "flexmix"}, each SEM iteration performs one stochastic
#' classification draw from the current posterior probabilities before the M-step. LCPA runs
#' exactly \code{starts} short trajectories of \code{maxiter.warmup} iterations, promotes the best
#' \code{nrep} trajectories by observed log-likelihood, continues each for
#' \code{control.flexmix$maxiter} SEM iterations, and retains the largest-likelihood final state
#' across the promoted runs.
#' This selection compares the final state returned by each flexmix SEM run; flexmix's
#' \code{classify = "SEM"} does not retain the largest-likelihood state visited within a run.
#' The stochastic classification and iteration control remain those of flexmix. LCPA supplies a
#' joint Gaussian M-step driver so that \code{"UE"}, \code{"UV"}, \code{"E0"}, \code{"V0"},
#' \code{"EE"}, \code{"VV"}, \code{"VE"}, \code{"EV"}, and custom equality lists obey the same
#' covariance definitions used by LCPA's EM estimator. Setting \code{control.flexmix$tol = 0}
#' prevents likelihood-based early termination. Only SEM is exposed; \code{par.ini} is not used.
#' @section Rmixmod stochastic strategies:
#' When \code{method = "Rmixmod"}, estimation uses \code{Rmixmod::mixmodCluster()}. The covariance
#' constraints map to Rmixmod Gaussian models as follows: \code{"V0"} to
#' \code{"Gaussian_pk_Lk_Bk"}, \code{"EE"} to \code{"Gaussian_pk_L_C"}, and \code{"VV"} to
#' \code{"Gaussian_pk_Lk_Ck"}.
#'
#' With \code{path="LCPA"}, the function generates exactly \code{starts} random balanced
#' partitions and runs exactly \code{maxiter.warmup} consecutive stochastic E-S-M iterations for every
#' warm-up start. The best \code{nrep} warm-up parameter sets are then passed directly to independent
#' Rmixmod SEM runs, and the finite solution with the largest final log-likelihood is retained.
#' Rmixmod's internal \code{"smallEM"} and \code{"SEMMax"} initialization searches are not used;
#' in particular, \code{"SEMMax"} is an initialization search rather than one consecutive SEM
#' trajectory for each user-level start. \code{par.ini} is not used for this method.
#' SEM stops after \code{control.Rmixmod$maxiter} iterations; an epsilon convergence
#' criterion is not defined for SEM in Rmixmod. Unlike flexmix's \code{classify = "SEM"}, Rmixmod
#' retains the largest-likelihood parameter state visited within each SEM run, so equal iteration
#' counts do not imply identical final-state selection.
#' Custom covariance-constraint lists are not accepted by the Rmixmod backend.
#'
#' With \code{path="Rmixmod"}, LCPA translates all non-\code{NULL} package-standard strategy controls
#' to \code{Rmixmod::mixmodStrategy()} and runs one native \code{mixmodCluster()} call. Its
#' \code{algorithm} may contain \code{"EM"}, \code{"CEM"}, \code{"SEM"}, or an ordered combination
#' of these algorithms. LCPA does
#' not add its outer \code{starts}, \code{maxiter.warmup}, or \code{nrep};
#' \code{control.Rmixmod$nrep} controls complete
#' strategy repetitions. The published strategy of Mulder et al.
#' (2015) uses 200 SEM iterations followed by EM with a relative likelihood-change tolerance of
#' \code{1e-5}. For reproducibility across Rmixmod versions, the example and simulation scripts
#' explicitly pin the contemporaneous documented defaults: \code{smallEM}, 50 initialization tries,
#' 5 initialization iterations, \code{tol.init=0.001}, and a 200-iteration EM limit.
#'
#' @section RMixtComp Stochastic EM:
#' When \code{method = "RMixtComp"}, LCPA calls \code{RMixtComp::mixtCompLearn()} in
#' classic, non-hierarchical learning mode with a univariate Gaussian model for every indicator.
#' This integration exposes only RMixtComp's stochastic EM (SEM) algorithm; it does not
#' introduce any other RMixtComp estimation algorithm. Conditional independence implies
#' diagonal, class-varying covariance matrices, so only \code{constraint = "V0"} is accepted.
#' Custom covariance-constraint lists are not accepted by the RMixtComp backend.
#' LCPA does not add its own \code{starts}, \code{maxiter.warmup}, or \code{nrep} stages to this
#' backend. Each native run performs RMixtComp initialization, SEM burn-in, recorded SEM
#' iterations, and then the fixed-parameter Gibbs burn-in and recorded Gibbs iterations.
#' One stochastic S-step is performed per SEM iteration, between the E-step and M-step.
#' Native repetition and parallelization are controlled only by
#' \code{control.RMixtComp$nrep} and \code{control.RMixtComp$ncores}.
#' RMixtComp's native random stream is independent of R, so LCPA passes it one integer drawn
#' from R's current random stream. The user only needs an external \code{set.seed()} for reproducibility.
#' The default algorithm controls reproduce \code{RMixtCompUtilities::createAlgo()}
#' defaults (version 4.1.4 or later).
#' \code{par.ini} is not used by this backend.
#'
#' @references
#' Biernacki, C. (2015). MixtComp software: Model-based clustering/imputation
#' with mixed data, missing data and uncertain data. *MISSDATA 2015*.
#' \url{https://inria.hal.science/hal-01253393}
#'
#' Leisch, F. (2004). FlexMix: A general framework for finite mixture models and latent
#' class regression in R. *Journal of Statistical Software, 11*(8), 1--18.
#' \doi{10.18637/jss.v011.i08}
#'
#' McLachlan, G. J., & Peel, D. (2000). *Finite mixture models*. John Wiley & Sons.
#'
#' Mulder, V. L., Lacoste, M., Martin, M. P., Richer-de-Forges, A., & Arrouays, D. (2015).
#' Understanding large-extent controls of soil organic carbon storage in relation to soil depth
#' and soil-landscape systems. *Global Biogeochemical Cycles, 29*(8), 1210--1229.
#' \doi{10.1002/2015GB005178}
#'
#' @examples
#'
#' library(LCPA)
#'
#' # Simulate bivariate continuous data for 2 profiles
#' set.seed(123)
#' data.obj <- sim.LPA(N = 500, I = 3, L = 2, constraint = "VV")
#' response <- data.obj$response
#'
#' ## It is strongly recommended to perform the following
#' ## standardization to obtain more stable results.
#' ## Standardization is not performed here in order to
#' ## compare estimated values with true values.
#' # response <- normalize(response)
#'
#' # Fit 2-profile model with VV constraint (default)
#' fit_vv <- LPA(response, L = 2, constraint = "VV")
#'
#' # Fit 2-profile model with E0 constraint using neural network estimation
#' # need Python
#' \dontrun{
#' fit_e0_nne <- LPA(response, L = 2, constraint = "E0", method = "NNE", nrep = 2)
#' }
#'
#' # Fit 2-profile model using Mplus
#' # Requires Mplus to be installed and available in system PATH.
#' # An empty 'files.path' instead uses the current working directory.
#' # This example creates a timestamped subdirectory
#' # (e.g., "Mplus_LPA_YYYY-MM-DD_HH-MM-SS") under './inst'
#' # to store all temporary Mplus files (.inp, .dat, .out, etc.).
#' # The 'inst' directory will be created if it does not exist.
#' # Setting files.clean=FALSE means temporary files will be preserved after execution.
#' \dontrun{
#' fit_mplus <- LPA(response, L = 2, method = "Mplus", constraint = list(c(1, 2), c(3, 3)),
#'                  control.Mplus = list(files.path = "inst", files.clean=FALSE))
#' }
#'
#' # Fit an EE model with flexmix SEM and the joint constrained M-step
#' # need flexmix
#' \dontrun{
#'   fit_flexmix <- LPA(response, L = 2, constraint = "EE", method = "flexmix",
#'                      nrep = 2, starts = 5, maxiter.warmup = 5,
#'                      control.flexmix = list(maxiter = 50))
#' }
#'
#' # Fit 2-profile model with the published Rmixmod SEM-to-EM strategy
#' # need Rmixmod
#' \dontrun{
#'   fit_rmixmod <- LPA(response, L = 2, constraint = "VV", method = "Rmixmod",
#'                      control.Rmixmod = list(path = "Rmixmod",
#'                                               algorithm = c("SEM", "EM"),
#'                                               nrep = 1,
#'                                               method.init = "smallEM",
#'                                               starts = 50,
#'                                               maxiter.init = 5,
#'                                               tol.init = 0.001,
#'                                               maxiter = c(200, 200),
#'                                               tol = c(NA, 1e-5)))
#' }
#'
#' @import reticulate
#' @export
#'
LPA <- function(response,
                L = 2, constraint = "VV",
                method="EM", par.ini = "random", is.sort=TRUE,
                starts=100, maxiter.warmup=20, nrep = 20,
                vis = TRUE,
                control.EM=NULL,
                control.Mplus=NULL,
                control.NNE=NULL,
                control.flexmix=NULL,
                control.Rmixmod=NULL,
                control.RMixtComp=NULL){

  call <- match.call()
  method <- match.arg(method, c("EM", "NNE", "Mplus", "flexmix", "Rmixmod", "RMixtComp"))
  par.ini <- .normalize.par.ini(par.ini, method)

  response <- as.matrix(response)
  N <- nrow(response)
  I <- ncol(response)

  default_control.EM <- list(maxiter=2000, tol=1e-4)
  default_control.Mplus <- list(
    maxiter=2000, tol=1e-4,
    files.path = "", files.clean = TRUE
  )
  default_control.Rmixmod <- .default.Rmixmod.control()
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
    device="CPU"
  )

  merge_and_clean_control <- function(user_control, default_control) {
    if (is.null(user_control)) return(default_control)
    valid_names <- intersect(names(user_control), names(default_control))
    cleaned_user_control <- user_control[valid_names]
    merged_control <- default_control
    for (param_name in names(cleaned_user_control)) {
      user_value <- cleaned_user_control[[param_name]]
      if (!is.null(user_value)) {
        merged_control[[param_name]] <- user_value
      }
    }
    return(merged_control)
  }

  control.EM <- merge_and_clean_control(control.EM, default_control.EM)
  control.Mplus <- merge_and_clean_control(control.Mplus, default_control.Mplus)
  control.NNE <- merge_and_clean_control(control.NNE, default_control.NNE)
  control.flexmix <- merge_and_clean_control(control.flexmix, default_control.flexmix)
  control.Rmixmod <- merge_and_clean_control(control.Rmixmod, default_control.Rmixmod)
  control.RMixtComp <- merge_and_clean_control(control.RMixtComp, default_control.RMixtComp)

  if(method != "RMixtComp" &&
     !(method == "Rmixmod" && identical(control.Rmixmod$path, "Rmixmod"))){
    .validate.training.stages(starts, maxiter.warmup, nrep)
  }

  method.requested <- method
  fit.requested <- function(){
  if(method == "NNE"){

    py_file <- system.file(
      "python", "Net_LPA.py",
      package = "LCPA",
      mustWork = TRUE
    )
    reticulate::source_python(py_file)

    res <- NN_LPA(response,
                  L=L,
                  par_ini=par.ini,
                  constraint=constraint,
                  nrep=nrep,
                  starts=starts,
                  maxiter_warmup=maxiter.warmup,
                  vis=vis,
                  hidden_layers=control.NNE$hidden.layers,
                  activation_function=control.NNE$activation.function,
                  use_attention=control.NNE$use.attention,
                  d_model=control.NNE$d.model,
                  nhead=control.NNE$nhead,
                  dim_feedforward=control.NNE$dim.feedforward,
                  eps=control.NNE$eps,
                  lambda_=control.NNE$lambda,
                  initial_temperature=control.NNE$initial.temperature,
                  cooling_rate=control.NNE$cooling.rate,
                  maxiter_sa=control.NNE$maxiter.sa,
                  threshold_sa=control.NNE$threshold.sa,
                  maxiter=control.NNE$maxiter,
                  patience_early=control.NNE$patience.early,
                  maxcycle=control.NNE$maxcycle,
                  lr = control.NNE$lr,
                  scheduler_patience = control.NNE$scheduler.patience,
                  scheduler_factor = control.NNE$scheduler.factor,
                  plot_interval = control.NNE$plot.interval,
                  device=control.NNE$device,
                  output_prefix=.estimation.output.prefix())

    res$Log.Lik.history <- unlist(res$Log.Lik.history)[2*(1:(length(unlist(res$Log.Lik.history))/2))]
    res$Log.Lik.nrep <- unlist(res$Log.Lik.nrep)

  }else if(method == "EM"){
    res <- EM.LPA(response, L = L, par.ini=par.ini, constraint = constraint, nrep=nrep,
                  starts=starts, maxiter.warmup=maxiter.warmup,
                  vis = vis,
                  maxiter = control.EM$maxiter,
                  tol = control.EM$tol)
  }else if(method == "Mplus"){
    res <- Mplus.LPA(response, L = L, constraint = constraint, nrep = nrep,
                     starts=starts, maxiter.warmup=maxiter.warmup,
                     vis = vis,
                     maxiter = control.Mplus$maxiter,
                     tol = control.Mplus$tol,
                     files.path = control.Mplus$files.path,
                     files.clean = control.Mplus$files.clean)
  }else if(method == "flexmix"){
    res <- flexmix.LPA(response, L = L, constraint = constraint,
                       nrep = nrep, starts = starts, maxiter.warmup = maxiter.warmup,
                       vis = vis,
                       control.flexmix = control.flexmix)
  }else if(method == "Rmixmod"){
    res <- Rmixmod.LPA(response, L = L, constraint = constraint,
                       nrep = nrep, starts = starts, maxiter.warmup = maxiter.warmup,
                       vis = vis,
                       control.Rmixmod = control.Rmixmod)
  }else if(method == "RMixtComp"){
    res <- RMixtComp.LPA(response, L = L, constraint = constraint,
                         vis = vis,
                         control.RMixtComp = control.RMixtComp)
  }
    res
  }

  estimation.error <- NULL
  res <- if(method == "EM"){
    fit.requested()
  }else{
    suppressWarnings(tryCatch(
      fit.requested(),
      error = function(e){
        estimation.error <<- conditionMessage(e)
        NULL
      }
    ))
  }

  finalize.result <- function(result){
    .stabilize.LPA.result(result, response, constraint)
  }
  is.valid.result <- function(result){
    valid <- is.list(result) && is.list(result$params) &&
      is.matrix(result$params$means) && all(dim(result$params$means) == c(L, I)) &&
      length(dim(result$params$covs)) == 3L && all(dim(result$params$covs) == c(I, I, L)) &&
      length(result$params$P.Z) == L && all(is.finite(result$params$P.Z)) &&
      all(result$params$P.Z > 0) &&
      is.matrix(result$P.Z.Xn) && all(dim(result$P.Z.Xn) == c(N, L)) &&
      all(is.finite(result$P.Z.Xn)) && all(rowSums(result$P.Z.Xn) > 0) &&
      length(result$Log.Lik) == 1L && is.finite(result$Log.Lik) &&
      all(is.finite(result$params$means)) && all(is.finite(result$params$covs))
    if(!valid) return(FALSE)
    all(vapply(seq_len(L), function(l){
      covariance <- result$params$covs[, , l]
      scale <- max(1, max(abs(covariance)))
      max(abs(covariance - t(covariance))) <=
        sqrt(.Machine$double.eps) * scale &&
        !is.null(tryCatch(chol(covariance), error = function(e){ NULL }))
    }, logical(1)))
  }

  res <- finalize.result(res)
  valid.result <- is.valid.result(res)

  if(method != "EM" && !valid.result){
    warning(method.requested, " failed; EM was used instead.", call. = FALSE)
    method <- "EM"
    res <- EM.LPA(response, L = L, par.ini = par.ini, constraint = constraint,
                  nrep = nrep, starts = starts, maxiter.warmup = maxiter.warmup,
                  vis = vis, maxiter = control.EM$maxiter, tol = control.EM$tol)
    res <- finalize.result(res)
    valid.result <- is.valid.result(res)
  }

  if(!valid.result) stop("EM failed to return a valid LPA solution", call. = FALSE)

  res$requested.method <- method.requested
  res$estimation.method <- method
  if(!is.null(estimation.error)) res$fallback.reason <- estimation.error

  if (is.sort) {
    posi <- order(res$params$P.Z, decreasing = TRUE)
    res$P.Z <- res$params$P.Z[posi]
    res$params$P.Z <- res$params$P.Z[posi]
    res$params$means <- res$params$means[posi, , drop = FALSE]
    res$params$covs  <- res$params$covs[, , posi, drop = FALSE]
    res$P.Z.Xn <- res$P.Z.Xn[, posi, drop = FALSE]
    res$Z <- match(res$Z, posi)
  }


  res$call <- call
  res$arguments <- list(
    response = response,
    L = L,
    constraint = constraint,
    method = method,
    par.ini = par.ini,
    is.sort = is.sort,
    starts = starts,
    maxiter.warmup = maxiter.warmup,
    nrep = nrep,
    vis = vis,
    control.EM = control.EM,
    control.Mplus = control.Mplus,
    control.NNE = control.NNE,
    control.flexmix = control.flexmix,
    control.Rmixmod = control.Rmixmod,
    control.RMixtComp = control.RMixtComp
  )

  rownames(res$params$means) <- .latent.group.names(L, "LPA")
  if(!is.null(colnames(response))){
    colnames(res$params$means) <- colnames(response)
  }else{
    colnames(res$params$means) <- paste0("V", 1:I)
  }

  covs.dimnames <- list(
    colnames(res$params$means), colnames(res$params$means),
    .latent.group.names(L, "LPA")
  )
  dimnames(res$params$covs) <- covs.dimnames

  colnames(res$P.Z.Xn) <- names(res$P.Z) <- names(res$params$P.Z) <-
    .latent.group.names(L, "LPA")

  class(res) <- "LPA"

  return(res)
}
