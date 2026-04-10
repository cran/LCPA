#' Fit Latent Class Analysis Models
#'
#' This function estimates parameters of a Latent Class Analysis (LCA; Hagenaars & McCutcheon, 2002) model using
#' either the Expectation-Maximization (EM) algorithm or Neural Network Estimation (NNE).
#' It supports flexible initialization strategies and provides comprehensive model diagnostics.
#'
#' @param response A numeric matrix of dimension \eqn{N \times I}, where \eqn{N} is the number of individuals/participants/observations
#'                 and \eqn{I} is the number of observed categorical items/variables/indicators. Each column must contain nominal-scale
#'                 discrete responses (e.g., integers representing categories).
#' @param L Integer specifying the number of latent classes (default: 2).
#' @param par.ini Specification for parameter initialization. Options include:
#'   \itemize{
#'     \item \code{"random"}: Completely random initialization (default).
#'     \item \code{"kmeans"}: Initializes parameters via K-means clustering on observed data (McLachlan & Peel, 2004).
#'     \item A \code{list} containing:
#'       \describe{
#'         \item{\code{par}}{An \eqn{L \times I \times K_{\max}} array of initial conditional probabilities for
#'                           each latent class, indicator, and response category (where \eqn{K_{\max}} is the maximum
#'                           number of categories across indicators).}
#'         \item{\code{P.Z}}{A numeric vector of length \eqn{L} specifying initial prior probabilities for latent classes.}
#'       }
#'   }
#' @param method Character string specifying estimation algorithm:
#'   \itemize{
#'     \item \code{"EM"}: Expectation-Maximization algorithm (default).
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
##'     \item{\code{lr}}{Learning rate, controlling the step size of neural network parameter updates (default: 0.025).}
#'     \item{\code{scheduler.patience}}{Patience for learning rate decay (if the loss function does not improve for more than `patience` consecutive epochs, the learning rate will be reduced) (default: 10).}
#'     \item{\code{scheduler.factor}}{Learning rate decay factor; the new learning rate equals the original learning rate multiplied by `scheduler.factor` (default: 0.70).}
#'     \item{\code{plot.interval}}{Interval (in epochs) for plotting training diagnostics (default: 100).}
#'     \item{\code{device}}{Specifies the hardware device; can be \code{"CPU"} (default) or \code{"GPU"}. If the GPU is not available, it automatically falls back to CPU.}
#'   }
#' @param control.Mplus List of control parameters for Mplus estimation:
#'   \describe{
#'     \item{\code{maxiter}}{Maximum iterations for Mplus optimization (default: 2000).}
#'     \item{\code{tol}}{Convergence tolerance for log-likelihood difference (default: 1e-4).}
#'     \item{\code{files.path}}{A character string specifying the directory path where Mplus will write its intermediate files
#'                              (e.g., \code{.inp} model input, \code{.dat} data file, \code{.out} output, and saved posterior probabilities).
#'                              This argument is \strong{required} — if \code{NULL} (the default), the function throws an error.
#'                              The specified directory must exist and be writable; if it does not exist, the function attempts to create it recursively.
#'                              A unique timestamped subdirectory (e.g., \code{"Mplus_LCA_YYYY-MM-DD_HH-MM-SS"}) will be created within this path
#'                              to store all run-specific files and avoid naming conflicts.
#'                              If it is an empty string (\code{""}), the timestamped subdirectory \code{"Mplus_LCA_YYYY-MM-DD_HH-MM-SS"}
#'                              will be created directly under R's current working directory (\code{\link[base]{getwd}()}).}
#'     \item{\code{files.clean}}{Logical. If \code{TRUE} (default), all intermediate files and the temporary working directory
#'                               created for this run are deleted upon successful completion or error exit (via \code{on.exit()}).
#'                               If \code{FALSE}, all generated files are retained in \code{files.path} (or the auto-generated temp dir)
#'                               for inspection or debugging. Note: when \code{files.path = NULL}, even if \code{files.clean = FALSE},
#'                               the temporary directory may still be cleaned up by the system later — for guaranteed persistence,
#'                               specify a custom \code{files.path}.}
#'   }
#'
#' @return An object of class \code{"LCA"} containing:
#'   \describe{
#'     \item{\code{params}}{List with estimated parameters:
#'       \describe{
#'         \item{\code{par}}{\eqn{L \times I \times K_{\max}} array of conditional response probabilities per latent class.}
#'         \item{\code{P.Z}}{Vector of length \eqn{L} with latent class prior probabilities.}
#'       }
#'     }
#'     \item{\code{npar}}{Number of free parameters in the model. see \code{\link[LCPA]{get.npar.LCA}}}
#'     \item{\code{Log.Lik}}{Log-likelihood of the final model. see \code{\link[LCPA]{get.Log.Lik.LCA}}}
#'     \item{\code{AIC}}{Akaike Information Criterion value.}
#'     \item{\code{BIC}}{Bayesian Information Criterion value.}
#'     \item{\code{best_BIC}}{Best BIC value across \code{nrep} runs (if applicable).}
#'     \item{\code{P.Z.Xn}}{\eqn{N \times L} matrix of posterior class probabilities for each observation.}
#'     \item{\code{P.Z}}{Vector of length \eqn{L} containing the prior probabilities/structural parameters/proportions for each latent class.}
#'     \item{\code{Z}}{Vector of length \eqn{N} with MAP-classified latent class memberships.}
#'     \item{\code{probability}}{Same as \code{params$par} (redundant storage for convenience).}
#'     \item{\code{Log.Lik.history}}{Vector tracking log-likelihood at each EM iteration.}
#'     \item{\code{Log.Lik.nrep}}{Vector of log-likelihoods from each replication run.}
#'     \item{\code{model}}{The optimal neural network model object (only for \code{method="NNE"}).
#'                 Contains the trained transformer architecture corresponding to \code{best_loss}.
#'                 This object can be used for further predictions or model inspection.}
#'     \item{\code{arguments}}{A list containing all input arguments}
#'   }
#'
#' @section EM Algorithm:
#' When \code{method = "EM"}, parameters are estimated via the Expectation-Maximization algorithm, which iterates between:
#' \itemize{
#'   \item \strong{E-step:} Compute posterior class probabilities given current parameters:
#'     \deqn{P(Z_n = l \mid \mathbf{X}_n) = \frac{\pi_l \prod_{i=1}^I P(X_{ni} = x_{ni} \mid Z_n=l)}{\sum_{k=1}^L \pi_k \prod_{i=1}^I P(X_{ni} = x_{ni} \mid Z_n=k)}}
#'     where \eqn{x_{ni}} is the standardized (0-based) response for person \eqn{n} on indicator \eqn{i} (see \code{\link[LCPA]{adjust.response}}).
#'   \item \strong{M-step:} Update parameters by maximizing expected complete-data log-likelihood:
#'     \itemize{
#'       \item Class probabilities: \eqn{\pi_l^{\text{new}} = \frac{1}{N} \sum_{n=1}^N P(Z_n = l \mid \mathbf{X}_n)}
#'       \item Conditional probabilities: \eqn{P(X_i = k \mid Z=l)^{\text{new}} = \frac{\sum_{n:x_{ni}=k} P(Z_n = l \mid \mathbf{X}_n)}{\sum_{n=1}^N P(Z_n = l \mid \mathbf{X}_n)}}
#'     }
#'   \item \strong{Convergence}: Stops when \eqn{|\log\mathcal{L}^{(t)} - \log\mathcal{L}^{(t-1)}| < \texttt{tol}} or maximum iterations reached.
#' }
#'
#' @section Neural Network Estimation (NNE):
#' When \code{method = "NNE"}, parameters are estimated using a hybrid neural network architecture
#' that combines feedforward layers with transformer-based attention mechanisms. This approach jointly
#' optimizes profile parameters and posterior probabilities through stochastic optimization enhanced
#' with simulated annealing. See \code{\link[LCPA]{install_python_dependencies}}. Key components include:
#'
#' \strong{Architecture}:
#' \describe{
#'   \item{Input Representation}{
#'     Observed categorical responses are converted to 0-based integer indices per indicator (not one-hot encoded).
#'     For example, original responses \eqn{[1, 2, 4]} become \eqn{[0, 1, 2]}.
#'   }
#'   \item{Feature Estimator (Feedforward Network)}{
#'     A fully-connected neural network with layer sizes specified by \code{hidden.layers} and activation
#'     function \code{activation.function} processes the integer-indexed responses. This network outputs
#'     unnormalized logits for posterior class membership (\eqn{N \times L} matrix).
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
#'   \item{Profile Parameter Estimation}{
#'     Global conditional probability parameters (\eqn{P(X_i = k \mid Z = l)}) are stored as learnable
#'     parameters \code{par} (an \eqn{L \times I \times K_{\max}} tensor). A \emph{masked softmax} is applied
#'     along categories to enforce:
#'     \itemize{
#'       \item Probabilities sum to 1 within each indicator-class pair
#'       \item Non-existent categories (beyond indicator's actual max response) are masked to zero probability
#'     }
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
#'       \item \code{CATEGORICAL} declaration for all indicator variables
#'       \item \code{ANALYSIS} block with optimization controls:
#'         \describe{
#'           \item{\code{TYPE = mixture}}{Standard mixture modeling setup}
#'           \item{\code{STARTS = starts nrep}}{Random \code{starts} and final stage optimizations}
#'           \item{\code{STITERATIONS = maxiter.wa}}{max itertions during \code{starts}.}
#'           \item{\code{MITERATIONS = maxiter}}{Maximum EM iterations}
#'           \item{\code{CONVERGENCE = tol}}{Log-likelihood convergence tolerance}
#'         }
#'       \item \code{MODEL} block with \code{\%OVERALL\%}
#'     }
#'   }
#'
#'   \item{Execution}{Calls Mplus via \code{MplusAutomation::mplusModeler()}, which:
#'     \itemize{
#'       \item Converts R data to Mplus-compatible format with automatic recoding
#'       \item Invokes Mplus executable (requires valid license and system PATH configuration)
#'     }
#'   }
#' }
#'
#' @references
#' Hagenaars, J. A. , & McCutcheon, A. L. (2002). Applied Latent Class Analysis. United Kingdom: Cambridge University Press.
#'
#' McLachlan, G. J., & Peel, D. (2004). Finite Mixture Models. Wiley.
#' https://books.google.com.sg/books?id=c2_fAox0DQoC
#'
#'
#' @examples
#' library(LCPA)
#'
#' # Example with simulated data
#' set.seed(123)
#' data.obj <- sim.LCA(N = 500, I = 4, L = 2, IQ=0.9)
#' response <- data.obj$response
#'
#' # Fit 2-class model with EM algorithm
#' \donttest{
#' fit.em <- LCA(response, L = 2, method = "EM", nrep = 10)
#' }
#'
#' # Fit 2-profile model using Mplus
#' # need Mplus
#' # NOTE: 'files.path' in control.Mplus is REQUIRED — function will error if not provided.
#' # Example creates a timestamped subfolder (e.g., "Mplus_LCA_YYYY-MM-DD_HH-MM-SS") under './'
#' # to store all temporary Mplus files (.inp, .dat, .out, etc.).
#' \dontrun{
#' fit.mplus <- LCA(response, L = 2, method = "Mplus", nrep = 3,
#'                  control.Mplus = list(files.path = ""))
#' }
#'
#' # Fit 2-class model with neural network estimation
#' # need Python
#' \dontrun{
#' fit.nne <- LCA(response, L = 2, method = "NNE", nrep = 3)
#' }
#'
#' @import reticulate
#' @importFrom utils modifyList
#' @export
#'
LCA <- function(response,
                L=2, par.ini="random",
                method="EM", is.sort=TRUE,
                nrep=20, starts=100, maxiter.wa=20,
                vis=TRUE,
                control.EM=NULL,
                control.Mplus=NULL,
                control.NNE=NULL){

  call <- match.call()

  response <- as.matrix(response)

  adjust.response.obj <- adjust.response(response)
  response <- adjust.response.obj$response
  poly.max <- adjust.response.obj$poly.max
  poly.value <- adjust.response.obj$poly.value
  poly.orig <- adjust.response.obj$poly.orig
  I <- ncol(response)
  N <- nrow(response)

  default_control.EM <- list(maxiter=2000, tol=1e-4)
  default_control.Mplus <- list(maxiter=2000, tol=1e-4, files.path = "", files.clean = TRUE)
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
      "python", "Net_LCA.py",
      package = "LCPA",
      mustWork = TRUE
    )
    reticulate::source_python(py_file)

    res <- NN_LCA(response,
                  L=L,
                  par_ini=par.ini,
                  nrep=nrep,
                  starts=starts,
                  maxiter_wa=maxiter.wa,
                  vis=vis,
                  hidden_layers=control.NNE$hidden.layers,
                  d_model=control.NNE$d.model,
                  nhead=control.NNE$nhead,
                  dim_feedforward=control.NNE$dim.feedforward,
                  eps=control.NNE$eps,
                  Lambda=control.NNE$lambda,
                  activation_function=control.NNE$activation.function,
                  use_attention=control.NNE$use.attention,
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

    if (!is.null(res$Log.Lik.history)) {
      res$Log.Lik.history <- unlist(res$Log.Lik.history)[2*(1:(length(unlist(res$Log.Lik.history))/2))]
    }
    if (!is.null(res$Log.Lik.nrep)) {
      res$Log.Lik.nrep <- unlist(res$Log.Lik.nrep)
    }

  } else if(method == "EM"){
    res <- EM.LCA(response, L=L, par.ini=par.ini, nrep=nrep,
                  starts=starts, maxiter.wa=maxiter.wa,
                  vis=vis,
                  maxiter = control.EM$maxiter,
                  tol = control.EM$tol)
  }else if(method == "Mplus"){
    res <- Mplus.LCA(response, L=L, nrep=nrep,
                     starts=starts, maxiter.wa=maxiter.wa,
                     vis=vis,
                     maxiter = control.Mplus$maxiter,
                     tol = control.Mplus$tol,
                     files.path = control.Mplus$files.path,
                     files.clean = control.Mplus$files.clean)
  }

  if (is.sort) {
    posi <- order(res$params$P.Z, decreasing = TRUE)
    res$P.Z <- res$params$P.Z[posi]
    res$params$P.Z <- res$params$P.Z[posi]
    res$params$par <- res$params$par[posi, , , drop = FALSE]
    res$P.Z.Xn <- res$P.Z.Xn[, posi, drop = FALSE]
    res$Z <- match(res$Z, posi)
  }


  res$call <- call
  res$arguments = list(
    response=response,
    L=L, par.ini=par.ini,
    method=method, is.sort=is.sort,
    nrep=nrep, starts=starts, maxiter.wa=maxiter.wa,
    vis=vis,
    control.EM=control.EM,
    control.Mplus=control.Mplus,
    control.NNE=control.NNE
  )

  probability <- vector("list", I)
  for(i in 1:I){
    probability[[i]] <- res$params$par[ , i, 1:poly.value[i]]
    if(L == 1){
      probability[[i]] <- matrix(probability[[i]], nrow=1, ncol=poly.value[i])
    }
    rownames(probability[[i]]) <- paste0("class.", 1:L)
    colnames(probability[[i]]) <- paste0("category.", poly.orig[i, 1:poly.value[i]])

  }
  if(!is.null(colnames(response))){
    names(probability) <- colnames(response)
  } else {
    names(probability) <- paste0("indicator.", 1:I)
  }
  res$probability <- probability

  dimnames(res$params$par) <- list(rownames(probability[[1]]), names(probability),
                                   paste0("category.", 1:poly.max))

  colnames(res$P.Z.Xn) <- names(res$P.Z) <- names(res$params$P.Z) <- paste0("class.", 1:L)

  if(vis){
    cat("\n")
  }

  class(res) <- "LCA"
  return(res)
}
