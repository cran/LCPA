#' Fit Latent Profile Analysis
#'
#' This function estimates parameters of a Latent Profile Analysis (LPA) model for continuous observed variables
#' using one of three methods: Expectation-Maximization (EM) algorithm, Neural Network Estimation (NNE), or
#' external Mplus software. It supports flexible covariance structures and initialization strategies.
#'
#' @param response A numeric matrix of dimension \eqn{N \times I}, where \eqn{N} is the number of individuals/participants/observations
#'                 and \eqn{I} is the number of continuous observed items/variables. Missing values are not allowed.
#'                 Note that \code{response} must be standardized using \code{\link[base]{scale}} or
#'                 \code{\link[LCPA]{normalize}} before input.
#' @param L Integer specifying the number of latent profiles (default: 2).
#' @param par.ini Specification for parameter initialization. Options include:
#'   \itemize{
#'     \item \code{"random"}: Random initialization of means and covariances (default).
#'     \item \code{"kmeans"}: Initializes parameters via K-means clustering on observed data (McLachlan & Peel, 2004).
#'     \item A \code{list} containing exactly three elements:
#'       \describe{
#'         \item{\code{means}}{An \eqn{L \times I} matrix of initial mean vectors for each profile.}
#'         \item{\code{covs}}{An \eqn{I \times I \times L} array of initial covariance matrices for each profile.}
#'         \item{\code{P.Z}}{A numeric vector of length \eqn{L} specifying initial prior probabilities for profiles.}
#'       }
#'   }
#' @param constraint Character string specifying covariance structure constraints:
#'   \describe{
#'     \item{\code{"VV"}}{Varying variances and varying covariances across profiles (heterogeneous full covariance; Default).}
#'     \item{\code{"VE"}}{Varying variances but equal correlations across profiles.}
#'     \item{\code{"EV"}}{Equal variances but varying covariances across profiles.}
#'     \item{\code{"EE"}}{Equal variances and equal covariances across profiles (homogeneous full covariance).}
#'     \item{\code{"E0"}}{Equal variances across profiles, zero covariances (diagonal with shared variances).}
#'     \item{\code{"V0"}}{Varying variances across profiles, zero covariances (diagonal with free variances).}
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
#' @param method Character string specifying estimation algorithm:
#'   \itemize{
#'     \item \code{"EM"}: Expectation-Maximization algorithm (Default).
#'     \item \code{"NNE"}: Neural Network Estimation with transformer architecture (experimental; uses transformer + simulated
#'                         annealing, more reliable than both \code{"EM"} and \code{"Mplus"}). See \code{\link[LCPA]{install_python_dependencies}}.
#'     \item \code{"Mplus"}: Calls external Mplus software for estimation.
#'                           Uses Mplus defaults for optimization unless overridden by \code{control.Mplus}.
#'   }
#' @param is.sort A logical value. If \code{TRUE} (Default), the latent classes will be ordered in descending
#'                order according to \code{P.Z}. All other parameters will be adjusted accordingly
#'                based on the reordered latent classes.
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
#'                              A unique timestamped subdirectory (e.g., \code{"Mplus_LPA_YYYY-MM-DD_HH-MM-SS"}) will be created within this path
#'                              to store all run-specific files and avoid naming conflicts.}
#'     \item{\code{files.clean}}{Logical. If \code{TRUE} (default), all intermediate files and the temporary working directory
#'                               created for this run are deleted upon successful completion or error exit (via \code{on.exit()}).
#'                               If \code{FALSE}, all generated files are retained in \code{files.path} (or the auto-generated temp dir)
#'                               for inspection or debugging. Note: when \code{files.path = NULL}, even if \code{files.clean = FALSE},
#'                               the temporary directory may still be cleaned up by the system later — for guaranteed persistence,
#'                               specify a custom \code{files.path}.}
#'   }
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
#'     \item{\code{best_BIC}}{Best BIC value across \code{nrep} runs (if applicable).}
#'     \item{\code{P.Z.Xn}}{\eqn{N \times L} matrix of posterior profile probabilities for each observation.}
#'     \item{\code{P.Z}}{Vector of length \eqn{L} containing the prior probabilities/structural parameters/proportions for each latent class.}
#'     \item{\code{Z}}{Vector of length \eqn{N} with MAP-classified profile memberships.}
#'     \item{\code{Log.Lik.history}}{Vector tracking log-likelihood at each EM iteration (only for \code{method="EM"}).}
#'     \item{\code{Log.Lik.nrep}}{Vector of log-likelihoods from each replication run.}
#'     \item{\code{model}}{The optimal model object:
#'       \itemize{
#'         \item For \code{method="NNE"}: Trained neural network model.
#'         \item For \code{method="Mplus"}: Estimated Mplus model.
#'       }
#'     }
#'   }
#'
#' @section EM Algorithm:
#' When \code{method = "EM"}, parameter estimation uses the Expectation-Maximization (EM) algorithm to maximize the observed-data log-likelihood:
#'
#' \deqn{
#'  \log \mathcal{L} = \sum_{n=1}^N \log \left[ \sum_{l=1}^L \pi_l \cdot \mathcal{N}(\mathbf{x}_n \mid \boldsymbol{\mu}_l, \boldsymbol{\Sigma}_l) \right]
#' }
#'
#' The algorithm iterates between two steps until convergence (change in log-likelihood < \code{tol} or max iterations reached):
#'
#' \describe{
#'   \item{E-step:}{
#'     Compute posterior class probabilities (responsibilities) for observation \eqn{n} and class \eqn{l}:
#'     \deqn{\tau_{nl} = \frac{\pi_l \cdot \mathcal{N}(\mathbf{x}_n \mid \boldsymbol{\mu}_l, \boldsymbol{\Sigma}_l)}{\sum_{k=1}^L \pi_k \cdot \mathcal{N}(\mathbf{x}_n \mid \boldsymbol{\mu}_k, \boldsymbol{\Sigma}_k)}}
#'     where \eqn{\mathcal{N}(\cdot)} is the multivariate normal density, \eqn{\pi_l} is the prior class probability, and \eqn{\boldsymbol{\mu}_l},
#'     \eqn{\boldsymbol{\Sigma}_l} are current parameters. Numerical stability is ensured via the log-sum-exp trick.
#'   }
#'
#'   \item{M-step:}{
#'     Update parameters using responsibilities \eqn{\tau_{nl}}:
#'     \itemize{
#'       \item Class probabilities: \eqn{\pi_l^{\text{new}} = \frac{1}{N}\sum_{n=1}^N \tau_{nl}}
#'       \item Class means: \eqn{\boldsymbol{\mu}_l^{\text{new}} = \frac{\sum_{n=1}^N \tau_{nl} \mathbf{x}_n}{\sum_{n=1}^N \tau_{nl}}}
#'       \item Class covariances: Updated under constraints:
#'         \describe{
#'           \item{\code{"VV"}}{\eqn{\boldsymbol{\Sigma}_l^{\text{new}} = \frac{\sum_{n=1}^N \tau_{nl} (\mathbf{x}_n - \boldsymbol{\mu}_l^{\text{new}}) (\mathbf{x}_n - \boldsymbol{\mu}_l^{\text{new}})^\top}{\sum_{n=1}^N \tau_{nl}}}}
#'           \item{\code{"EE"}}{Shared covariance: \eqn{\boldsymbol{\Sigma}^{\text{new}} = \frac{\sum_{l=1}^L \sum_{n=1}^N \tau_{nl} (\mathbf{x}_n - \boldsymbol{\mu}_l^{\text{new}}) (\mathbf{x}_n - \boldsymbol{\mu}_l^{\text{new}})^\top}{\sum_{l=1}^L \sum_{n=1}^N \tau_{nl}}}}
#'           \item{\code{"VE"} (default) / \code{"EV"}}{Hybrid constraints (e.g., \code{"VE"}: varying variances, equal correlations).
#'                                            Off-diagonal elements use weighted averages across classes; diagonals retain class-specific values.}
#'           \item{Custom constraints}{User-specified variances/covariances (e.g., \code{list(c(1,2), c(2, 2))}, meaning the
#'                                     covariates of observed variable 1 and observed variable 2 are equal across latent classes,
#'                                     and the variance of observed variable 2 is equal across classes) are forced equal across
#'                                     classes via weighted averaging.}
#'         }
#'     }
#'     Covariance matrices are regularized to ensure positive definiteness:
#'     \itemize{
#'       \item Eigenvalues < \code{jitter} (1e-10) are replaced with \code{jitter}
#'       \item Failed Cholesky decompositions trigger diagonal jittering or perturbation of non-constrained elements
#'     }
#'   }
#' }
#'
#' \strong{Edge Handling}:
#' \itemize{
#'   \item Empty classes (\eqn{\sum_n \tau_{nl} < 10^{-5}}) are reinitialized by redistributing responsibilities.
#'   \item Non-finite likelihoods trigger fallback to previous valid parameters or covariance perturbation.
#'   \item Univariate cases (\eqn{I=1}) bypass Cholesky decomposition for direct variance updates.
#' }
#'
#' @section Neural Network Estimation (NNE):
#' When \code{method = "NNE"}, parameters are estimated using a hybrid neural network architecture
#' combining fully-connected layers with transformer-based attention mechanisms. This approach jointly
#' optimizes profile parameters and posterior probabilities through stochastic optimization with
#' simulated annealing. See \code{\link[LCPA]{install_python_dependencies}}. Key components include:
#'
#' \strong{Architecture}:
#' \describe{
#'   \item{Input Representation:}{
#'     Continuous observed variables \eqn{\mathbf{x}_n \in \mathbb{R}^I} are standardized (mean-centered, unit-variance)
#'     internally during training to improve numerical stability. No encoding is needed — raw values are fed directly.
#'   }
#'   \item{Feature Encoder (Feedforward Network):}{
#'     A multi-layer perceptron with architecture defined by \code{hidden.layers} and \code{activation.function}
#'     maps the continuous input vector into a latent space of dimension \code{d.model}. This layer learns non-linear
#'     feature combinations predictive of latent profile membership.
#'   }
#'   \item{Attention Refiner (Transformer Encoder)}{
#'     A transformer encoder with \code{nhead} attention heads that learns latent class prior probabilities
#'     \eqn{\boldsymbol{\pi} = (\pi_1, \pi_2, \dots, \pi_L)} directly from observed responses.
#'     \itemize{
#'       \item Input: \code{response} matrix (\eqn{N \times I}), where \eqn{N} = observations, \eqn{I} = continuous variables.
#'       \item Mechanism: Self-attention dynamically weighs variable importance during profile assignment, capturing complex
#'                        multivariate interactions.
#'       \item Output: Class prior vector \eqn{\boldsymbol{\pi}} computed as the mean of posteriors:
#'                     \deqn{\pi_l = \frac{1}{N}\sum_{n=1}^N attention(\mathbf{X}_n)}
#'                     This ensures probabilistic consistency with the mixture model framework.
#'     }
#'   }
#'   \item{Parameter Head (Means & Covariances):}{
#'     Two separate projection heads branch from the transformer output:
#'     \itemize{
#'       \item Means Head: Linear projection to \eqn{L \times I} matrix \eqn{\boldsymbol{\mu}_l}.
#'       \item Covariance Head: Outputs lower-triangular elements of Cholesky factors \eqn{\mathbf{L}_l}
#'               for each profile. Diagonal elements are passed through \code{softplus} to ensure positivity;
#'               off-diagonals use \code{tanh} scaled by 1.2 to bound magnitude and promote stability.
#'               The full covariance is reconstructed via \eqn{\boldsymbol{\Sigma}_l = \mathbf{L}_l \mathbf{L}_l^\top}.
#'     }
#'     After reconstruction, covariance constraints (e.g., \code{"EE"}, \code{"V0"}, or custom lists) are applied by
#'     averaging constrained elements across profiles and re-symmetrizing.
#'   }
#' }
#'
#' \strong{Optimization Strategy}:
#' \itemize{
#'   \item \strong{Hybrid Training Protocol}: Alternates between:
#'     \itemize{
#'       \item \emph{Gradient-based phase}: AdamW optimizer minimizes negative log-likelihood with weight decay regularization:
#'         \deqn{-\log \mathcal{L} + \lambda \|\boldsymbol{\theta}\|_2^2}
#'         where \eqn{\lambda} is controlled by \code{lambda} (default: 1e-5). Learning rate decays adaptively when
#'         loss plateaus (controlled by \code{scheduler.patience} and \code{scheduler.factor}).
#'       \item \emph{Simulated annealing phase}: After gradient-based early stopping (\code{maxiter.early}), parameters
#'         are perturbed with noise scaled by temperature:
#'         \deqn{\theta_{\text{new}} = \theta_{\text{current}} + \mathcal{N}(0, \theta_{\text{current}} \times \frac{T}{T_0})}
#'         Temperature \eqn{T} decays geometrically (\eqn{T \leftarrow T \times \text{cooling.rate}}) from
#'         \code{initial.temperature} until \code{threshold.sa} is reached. This escapes poor local minima.
#'     }
#'     Each full cycle (gradient descent + annealing) repeats up to \code{maxcycle} times.
#'   \item \strong{Model Selection}: Across \code{nrep} random restarts (using Dirichlet-distributed initializations
#'     or K-means), the solution with lowest BIC is retained.
#'   \item \strong{Diagnostics}: Training loss, annealing points, and global best solution are plotted when \code{vis=TRUE}.
#' }
#'
#' \strong{Constraint Handling}:
#' \itemize{
#'   \item Covariance constraints (\code{constraint}) are enforced \emph{after} activation via:
#'     \itemize{
#'       \item Shared Parameters: Variances/covariances marked for equality are replaced by
#'                                their average across profiles.
#'       \item Positive Definiteness: Non-positive definite matrices are corrected via eigenvalue clamping,
#'                                           diagonal jittering, or adaptive Cholesky decomposition.
#'     }
#'   \item Custom constraints: e.g., \code{list(c(1,2), c(3,3))}, force equality of specific covariance elements
#'           across profiles, with symmetry (\eqn{\sigma_{12} = \sigma_{21}}) automatically enforced.
#' }
#'
#' @section Mplus:
#' When \code{method = "Mplus"}, estimation is delegated to external Mplus software.
#' The function automates the entire workflow:
#'
#' \strong{Workflow}:
#' \describe{
#'   \item{Temporary Directory Setup}{Creates \code{inst/Mplus} to store:
#'     \itemize{
#'       \item Mplus input syntax (\code{.inp})
#'       \item Data file in Mplus format (\code{.dat})
#'       \item Posterior probabilities output (\code{.dat})
#'     }
#'     Files are automatically deleted after estimation unless \code{control.Mplus$clean.files = FALSE}.
#'   }
#'
#'   \item{Syntax Generation}{Constructs Mplus syntax with:
#'     \itemize{
#'       \item \code{CLASSES = c1(L)} specification for \eqn{L} latent classes
#'       \item \code{ANALYSIS} block with optimization controls:
#'         \describe{
#'           \item{\code{TYPE = mixture}}{Standard mixture modeling setup}
#'           \item{\code{STARTS = starts nrep}}{Random \code{starts} and final stage optimizations}
#'           \item{\code{STITERATIONS = maxiter.wa}}{max itertions during \code{starts}.}
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
#' \strong{Constraint Handling}:
#' \itemize{
#'   \item Covariance constraints (\code{constraint}) are enforced \emph{after} activation via:
#'     \itemize{
#'       \item Shared Parameters: Variances/covariances marked for equality are replaced by
#'                                their average across profiles.
#'       \item Positive Definiteness: Non-positive definite matrices are corrected via eigenvalue clamping,
#'                                           diagonal jittering, or adaptive Cholesky decomposition.
#'     }
#'   \item Custom constraints: e.g., \code{list(c(1,2), c(3,3))}, force equality of specific covariance elements
#'           across profiles, with symmetry (\eqn{\sigma_{12} = \sigma_{21}}) automatically enforced.
#' }
#'
#' @references
#' McLachlan, G. J., & Peel, D. (2004). Finite Mixture Models. Wiley. https://books.google.com.sg/books?id=c2_fAox0DQoC
#'
#' @examples
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
#' # NOTE: 'files.path' in control.Mplus is REQUIRED — the function will
#' # throw an error if not provided.
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
#' @import reticulate
#' @export
#'
LPA <- function(response,
                L = 2, par.ini = "random",
                constraint = "VV",
                method="EM", is.sort=TRUE,
                nrep = 20, starts=100, maxiter.wa=20,
                vis = TRUE,
                control.EM=NULL,
                control.Mplus=NULL,
                control.NNE=NULL){

  call <- match.call()

  response <- as.matrix(response)
  N <- nrow(response)
  I <- ncol(response)

  default_control.EM <- list(maxiter=2000, tol=1e-4)
  default_control.Mplus <- list(maxiter=2000, tol=1e-4, files.path = NULL, files.clean = TRUE)
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
                  maxiter_wa=maxiter.wa,
                  vis=vis,
                  hidden_layers=control.NNE$hidden.layers,
                  activation_function=control.NNE$activation.function,
                  use_attention=control.NNE$use.attention,
                  d_model=control.NNE$d.model,
                  nhead=control.NNE$nhead,
                  dim_feedforward=control.NNE$dim.feedforward,
                  eps=control.NNE$eps,
                  Lambda=control.NNE$lambda,
                  initial_temperature=control.NNE$initial.temperature,
                  cooling_rate=control.NNE$cooling.rate,
                  maxiter_sa=control.NNE$maxiter.sa,
                  threshold_sa=control.NNE$threshold.sa,
                  maxiter=control.NNE$maxiter,
                  maxiter_early=control.NNE$maxiter.early,
                  maxcycle=control.NNE$maxcycle,
                  lr = control.NNE$lr,
                  scheduler_patience = control.NNE$scheduler.patience,
                  scheduler_factor = control.NNE$scheduler.factor,
                  plot_interval = control.NNE$plot.interval,
                  device=control.NNE$device)

    res$Log.Lik.history <- unlist(res$Log.Lik.history)[2*(1:(length(unlist(res$Log.Lik.history))/2))]
    res$Log.Lik.nrep <- unlist(res$Log.Lik.nrep)

  }else if(method == "EM"){
    res <- EM.LPA(response, L = L, par.ini=par.ini, constraint = constraint, nrep=nrep,
                  starts=starts, maxiter.wa=maxiter.wa,
                  vis = vis,
                  maxiter = control.EM$maxiter,
                  tol = control.EM$tol)
  }else if(method == "Mplus"){
    res <- Mplus.LPA(response, L = L, constraint = constraint, nrep = nrep,
                     starts=starts, maxiter.wa=maxiter.wa,
                     vis = vis,
                     maxiter = control.Mplus$maxiter,
                     tol = control.Mplus$tol,
                     files.path = control.Mplus$files.path,
                     files.clean = control.Mplus$files.clean)
  }

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
  res$arguments = list(
    response=response,
    L=L, par.ini=par.ini, constraint=constraint,
    method=method, is.sort=is.sort,
    nrep=nrep, starts=starts, maxiter.wa=maxiter.wa,
    vis=vis,
    control.EM=control.EM,
    control.NNE=control.NNE,
    control.Mplus=control.Mplus
  )

  rownames(res$params$means) <- paste0("Class.", 1:L)
  if(!is.null(colnames(response))){
    colnames(res$params$means) <- colnames(response)
  }else{
    colnames(res$params$means) <- paste0("V", 1:I)
  }

  covs.dimnames <- list(colnames(res$params$means), colnames(res$params$means), paste0("Class.", 1:L))
  dimnames(res$params$covs) <- covs.dimnames

  colnames(res$P.Z.Xn) <- names(res$P.Z) <- names(res$params$P.Z) <- paste0("Class.", 1:L)

  class(res) <- "LPA"

  return(res)
}
