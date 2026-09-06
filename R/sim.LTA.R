#' Simulate Data for Latent Transition Analysis (LTA)
#'
#' Simulates longitudinal latent class/profile data where initial class membership and transition
#' probabilities may be influenced by time-varying covariates. Supports both Latent Class Analysis (LCA)
#' for categorical outcomes and Latent Profile Analysis (LPA) for continuous outcomes. Measurement
#' invariance is assumed by default (identical indicator parameters across time).
#'
#' @param N Integer; sample size.
#' @param I Integer; number of observed indicators/items/indicators per time point.
#' @param L Integer; number of latent classes/profiles.
#' @param times Integer; number of time points (must be \eqn{\geq 1}).
#' @param type Character; type of latent model. \code{"LCA"} for categorical indicators (default),
#'   \code{"LPA"} for continuous indicators.
#' @param distribution Character; distribution of initial class probabilities when not using covariates or \code{params}.
#'   Options: \code{"uniform"} (equal probabilities) or \code{"random"} (Dirichlet-distributed, default).
#' @param constraint Character; covariance structure for LPA (\code{type="LPA"} only). Options:
#'   \code{"UE"} and \code{"UV"} for univariate responses; \code{"E0"},
#'   \code{"V0"}, \code{"EE"}, \code{"EV"}, \code{"VE"}, and \code{"VV"}
#'   for multivariate responses. See \code{\link[LCPA]{sim.LPA}}. The default
#'   is \code{"VV"}.
#' @param poly.value Integer; number of categories for polytomous LCA indicators (default: 5).
#' @param IQ Character; method for generating indicator discrimination in LCA. \code{"random"} (default) or fixed values.
#' @param mean.range Numeric vector; range for randomly generated class means in LPA (default: \code{c(-2, 2)}).
#' @param covs.range Numeric vector; range for covariance matrix diagonals in LPA (default: \code{c(0.01, 4)}).
#' @param params List or NULL; pre-specified measurement and initial-class parameters in the final
#'   first-time-point class order (see Details).
#' @param is.sort A logical value. If \code{TRUE} (default), internally generated first-time-point
#'   classes are established in decreasing model-implied probability order. Supplied parameters
#'   already use this final order and are never reordered; a supplied \code{beta}, \code{params$P.Z},
#'   or \code{params$Z} must therefore imply decreasing first-time-point class proportions. If
#'   \code{FALSE}, the supplied or generated order is retained. No later time point is reordered.
#' @param rate List of matrices or NULL; transition probability matrices in the final class order for non-covariate mode.
#'   Each matrix is \eqn{L \times L} with rows summing to 1. If \code{NULL} (default), matrices are
#'   generated with 0.7 diagonal probability and uniform off-diagonals. Ignored when \code{times=1}.
#' @param covariates List of matrices or NULL; covariate matrices for each time point. Each matrix must have
#'   dimensions \eqn{N\times(U_t+1)} and include an intercept column (first
#'   column must be all 1s). If \code{NULL},
#'   covariate mode is disabled. See Details for automatic coefficient generation.
#' @param ref.class Integer between 1 and \code{L}; reference class in the class order established
#'   at the first time point. When \code{is.sort=TRUE}, this is the position after ordering the
#'   first-time-point classes by decreasing \code{P.Z}. The same class order is retained at every
#'   later time point.
#' @param beta Matrix or NULL; initial state regression coefficients of dimension
#'   \eqn{(U_1+1)\times L}
#'   in the final first-time-point class order. Column \code{ref.class} must be zero. Supplied
#'   coefficients determine the logits and are returned unchanged.
#'   If \code{NULL} and covariates are used, coefficients are randomly generated from \eqn{\text{Uniform}(-1, 1)}.
#' @param gamma List or NULL; transition regression coefficients. Must be a list of length \code{times-1}.
#'   Each element \eqn{t} is a list of length \eqn{L} (previous state). Each sub-list contains \eqn{L} vectors
#'   (next state). The supplied coefficients refer to the class order established at the first time
#'   point and are returned unchanged. Ignored when \code{times=1}.
#'   If \code{NULL} and covariates are used with \code{times>=2}, coefficients are randomly generated from
#'   \eqn{\text{Uniform}(-1, 1)} for non-reference classes.
#'
#' @return A list of class \code{"sim.LTA"} containing:
#' \describe{
#'   \item{\code{responses}}{List of length \code{times}; observed data matrices (\eqn{N \times I}).}
#'   \item{\code{Zs}}{List of length \code{times}; true latent class memberships (\eqn{N \times 1} vectors).}
#'   \item{\code{P.Zs}}{List of length \code{times}; marginal class probabilities at each time.}
#'   \item{\code{par}}{Indicator parameters for LCA (if \code{type="LCA"}).}
#'   \item{\code{means}}{Class means for LPA (if \code{type="LPA"}).}
#'   \item{\code{covs}}{Class covariance matrices for LPA (if \code{type="LPA"}).}
#'   \item{\code{poly.value}}{Category counts for LCA (if \code{type="LCA"}).}
#'   \item{\code{rate}}{True transition matrices (non-covariate mode only; \code{NULL} when \code{times=1}).}
#'   \item{\code{covariates}}{List of covariate matrices used (covariate mode only).}
#'   \item{\code{beta}}{True initial state coefficients (covariate mode only).}
#'   \item{\code{gamma}}{True transition coefficients (covariate mode only; \code{NULL} when \code{times=1}).}
#'   \item{\code{ref.class}}{Reference class of the returned \code{beta} and \code{gamma} coefficients.}
#'   \item{\code{call}}{Function call.}
#'   \item{\code{arguments}}{Input arguments.}
#' }
#'
#' @section Model Specification:
#' \describe{
#'   \item{Initial Class Probabilities (with covariates):}{
#'     For participant \eqn{n} at time 1, the probability of belonging to
#'     latent class \eqn{l} is
#'     \deqn{P(Z_{n1}=l\mid\boldsymbol{\zeta}_{n1}) =
#'       \frac{\exp(\boldsymbol{\beta}_l^\top\boldsymbol{\zeta}_{n1})}
#'            {\sum_{h=1}^L\exp(\boldsymbol{\beta}_h^\top
#'            \boldsymbol{\zeta}_{n1})}.}
#'     Here
#'     \eqn{\boldsymbol{\zeta}_{n1}=
#'     (1,\zeta_{n11},\ldots,\zeta_{n1U_1})^\top}; the leading 1 is the
#'     intercept, and \eqn{u=1,\ldots,U_1} indexes observed covariates. The
#'     coefficient vector
#'     \eqn{\boldsymbol{\beta}_l=(\beta_{l0},\beta_{l1},\ldots,
#'     \beta_{lU_1})^\top} has the corresponding intercept and slopes.
#'     The class selected by \code{ref.class} is the reference class and has a zero coefficient vector.
#'   }
#'   \item{Transition Probabilities (with covariates and times>=2):}{
#'     For participant \eqn{n} transitioning from class \eqn{k} at time
#'     \eqn{t-1} to class \eqn{l} at time \eqn{t} (\eqn{t\geq2}),
#'     \deqn{P(Z_{nt}=l\mid Z_{n,t-1}=k,\boldsymbol{\zeta}_{nt}) =
#'       \frac{\exp(\boldsymbol{\gamma}_{klt}^\top\boldsymbol{\zeta}_{nt})}
#'            {\sum_{h=1}^L\exp(\boldsymbol{\gamma}_{kht}^\top
#'            \boldsymbol{\zeta}_{nt})}.}
#'     Here
#'     \eqn{\boldsymbol{\zeta}_{nt}=
#'     (1,\zeta_{nt1},\ldots,\zeta_{ntU_t})^\top}, and
#'     \eqn{\boldsymbol{\gamma}_{klt}=
#'     (\gamma_{klt0},\gamma_{klt1},\ldots,\gamma_{kltU_t})^\top} contains the
#'     corresponding intercept and slopes. The destination class selected by
#'     \code{ref.class} has a zero
#'     coefficient vector for every origin class.
#'   }
#'   \item{Without Covariates or When times=1:}{
#'     Initial probabilities follow a multinomial distribution with probabilities \eqn{\boldsymbol{\pi} = (\pi_1, \dots, \pi_L)}.
#'     When \eqn{times \geq 2}, transitions follow a Markov process with fixed
#'     probabilities
#'     \eqn{\tau_{kl}^{(t)}=P(Z_{nt}=l\mid Z_{n,t-1}=k)}, where
#'     \eqn{\sum_{l=1}^L\tau_{kl}^{(t)}=1} for each origin class \eqn{k} and
#'     time \eqn{t}.
#'   }
#' }
#'
#'
#' @details
#' Covariate Requirements:
#' \itemize{
#'   \item Covariate matrices must include an intercept (first column = 1). If omitted, the function adds an intercept
#'     and issues a warning.
#'   \item When \code{covariates} is provided but \code{beta} or \code{gamma} is \code{NULL}, coefficients are
#'     randomly generated from \eqn{\text{Uniform}(-1, 1)} (non-reference classes only).
#'   \item All supplied simulation parameters use the final first-time-point class order. They are
#'     never reordered or reparameterized. Internally generated initial coefficients are ordered and
#'     reparameterized before they are returned. The same order is used for every transition.
#' }
#'
#' Parameter Compatibility:
#' \itemize{
#'   \item Use \code{params} to fix indicator parameters (LCA) or class means/covariances (LPA) across simulations.
#'   \item In non-covariate mode, \code{rate} must be a list of \eqn{(times-1)} valid transition matrices (ignored when \code{times=1}).
#'   \item In covariate mode with \code{times>=2}, all three (\code{covariates}, \code{beta}, \code{gamma}) must be consistent in dimensions.
#' }
#'
#' @examples
#' ####################### Example 1: Single time point (times=1) ######################
#' library(LCPA)
#' set.seed(123)
#' sim_single <- sim.LTA(N = 200, I = 4, L = 3, times = 1, type = "LCA")
#' print(sim_single)
#'
#' ####################### Example 2: LPA without covariates ######################
#' set.seed(123)
#' sim_lta <- sim.LTA(N = 200, I = 3, L = 3, times = 3, type = "LPA", constraint = "VE")
#' print(sim_lta)
#'
#' ################## Example 3: With custom covariates (times>=2) ######################
#' set.seed(123)
#' N <- 200 ## sample size
#'
#' ## Covariates at time point T1
#' covariates.inter <- rep(1, N) # Intercept term is always 1 for each n
#' covariates.X1 <- rnorm(N)     # Covariate X1 is a continuous variable
#' covariates.X2 <- rbinom(N, 1, 0.5) # Covariate X2 is a binary variable
#' covariates.X1.X2 <- covariates.X1 * covariates.X2 # Interaction between covariates X1 and X2
#' covariates.T1 <- cbind(inter=covariates.inter, X1=covariates.X1,
#'                        X2=covariates.X2, X1.X2=covariates.X1.X2) # Combine into covariates at T1
#'
#' ## Covariates at time point T2
#' covariates.inter <- rep(1, N) # Intercept term is always 1 for each n
#' covariates.X1 <- rnorm(N)     # Covariate X1 is a continuous variable
#' covariates.X2 <- rbinom(N, 1, 0.5) # Covariate X2 is a binary variable
#' covariates.X1.X2 <- covariates.X1 * covariates.X2 # Interaction between covariates X1 and X2
#' covariates.T2 <- cbind(inter=covariates.inter, X1=covariates.X1,
#'                        X2=covariates.X2, X1.X2=covariates.X1.X2) # Combine into covariates at T2
#'
#' covariates <- list(t1=covariates.T1, t2=covariates.T2) # Combine into final covariates list
#'
#' ## Simulate beta coefficients
#' # 3x3 matrix (last column is zero because the last category is used as reference)
#' beta <- matrix(c( 1.3,  0.5, 0.0,
#'                   0.1,  0.4, 0.0,
#'                  -0.6, -0.8, 0.0,
#'                  -0.3, -0.2, 0.0), ncol=3, byrow=TRUE)
#'
#' ## Simulate gamma coefficients (only needed when times>=2)
#' gamma <- list(
#'   lapply(1:3, function(l) {
#'     lapply(1:3, function(k) if(k < 3)
#'            runif(4, -1.0, 1.0) else c(0, 0, 0, 0)) # Last class as reference
#'   })
#' )
#'
#' ## Simulate the data
#' sim_custom <- sim.LTA(
#'   N=N, I=4, L=3, times=2, type="LPA",
#'   ref.class=3,
#'   covariates=covariates,
#'   beta=beta,
#'   gamma=gamma
#' )
#'
#' summary(sim_custom)
#'
#' @export
sim.LTA <- function(N=500, I=5, L=3, times=2, type="LCA",
                    distribution="random",
                    constraint = "VV", poly.value=5, IQ="random",
                    mean.range = c(-2, 2), covs.range = c(0.01, 4),
                    params=NULL, is.sort=TRUE,
                    rate=NULL, covariates = NULL, ref.class=L,
                    beta = NULL, gamma = NULL) {

  call <- match.call()

  if (times < 1) stop("times must be at least 1")
  if (L < 2) stop("L must be at least 2")
  if (length(ref.class) != 1L || !is.numeric(ref.class) ||
      !is.finite(ref.class) || ref.class != as.integer(ref.class) ||
      ref.class < 1L || ref.class > L) {
    stop("ref.class must be an integer between 1 and L")
  }
  ref.class <- as.integer(ref.class)

  use_covariates <- !is.null(covariates)
  beta.supplied <- !is.null(beta)

  if (use_covariates) {
    if (!is.list(covariates) || length(covariates) != times) {
      stop("covariates must be a list of length 'times'")
    }
    for (t in 1:times) {
      if (!is.matrix(covariates[[t]]) || nrow(covariates[[t]]) != N) {
        stop(sprintf("covariates[[%d]] must be a %d x p matrix", t, N))
      }
      if (any(abs(covariates[[t]][, 1] - 1) > 1e-6)) {
        warning(sprintf("First column of covariates[[%d]] not all 1s; adding intercept", t))
        covariates[[t]] <- cbind(1, covariates[[t]])
      }
    }

    p1 <- ncol(covariates[[1]])
    if (is.null(beta)) {
      beta_mat <- matrix(0, p1, L)
      free.classes <- setdiff(seq_len(L), ref.class)
      beta_mat[, free.classes] <- matrix(
        runif(p1 * (L - 1L), -1, 1), p1, L - 1L
      )
    } else {
      if (!is.matrix(beta) || nrow(beta) != p1 || ncol(beta) != L) {
        stop(sprintf("beta must be %d x %d matrix (p1 x L)", p1, L))
      }
      if(any(abs(beta[, ref.class]) > sqrt(.Machine$double.eps))){
        stop("The ref.class column of beta must be zero")
      }
      beta_mat <- beta
    }

    gamma_list <- NULL
    if (times >= 2) {
      gamma_list <- vector("list", times-1)
      for (t in 1:(times-1)) {
        pt <- ncol(covariates[[t+1]])
        if (is.null(gamma)) {
          gamma_t <- vector("list", L)
          for (l in 1:L) {
            gamma_t[[l]] <- vector("list", L)
            for (k in 1:L) {
              gamma_t[[l]][[k]] <- if (k != ref.class) runif(pt, -1, 1) else rep(0, pt)
            }
          }
          gamma_list[[t]] <- gamma_t
        } else {
          if (length(gamma) != times-1) {
            stop("gamma must be a list of length times-1")
          }
          gamma_t <- gamma[[t]]
          if (length(gamma_t) != L) {
            stop(sprintf("gamma[[%d]] must contain L=%d elements (one per current state)", t, L))
          }

          gamma_valid <- vector("list", L)
          for (l in 1:L) {
            if (length(gamma_t[[l]]) != L) {
              stop(sprintf("gamma[[%d]][[%d]] must contain L=%d elements (one per next state)", t, l, L))
            }
            gamma_valid[[l]] <- vector("list", L)
            for (k in 1:L) {
              coef_vec <- gamma_t[[l]][[k]]
              if (!is.numeric(coef_vec) || length(coef_vec) != pt) {
                stop(sprintf("gamma[[%d]][[%d]][[%d]] must be numeric vector of length %d",
                             t, l, k, pt))
              }
              if(k == ref.class && any(abs(coef_vec) > sqrt(.Machine$double.eps))){
                stop(sprintf(
                  "gamma[[%d]][[%d]][[%d]] must be zero because ref.class=%d",
                  t, l, k, ref.class
                ))
              }
              gamma_valid[[l]][[k]] <- coef_vec
            }
          }
          gamma_list[[t]] <- gamma_valid
        }
      }
    } else {
      if (!is.null(gamma)) {
        warning("gamma parameter is ignored when times=1 in covariate mode")
      }
    }

  } else {
    if (!is.null(beta) || !is.null(gamma)) {
      warning("beta/gamma ignored in non-covariate mode; provide 'covariates' to enable covariate mode")
    }
    if(times > 1){
      if (is.null(rate)) {
        rate <- vector("list", times - 1)
        for (t in 1:(times-1)) {
          rate[[t]] <- matrix((1-0.7)/(L-1), L, L)
          diag(rate[[t]]) <- 0.7
        }
      } else if (length(rate) != times-1) {
        stop("Length of 'rate' must equal times-1")
      }

      for (t in 1:(times-1)) {
        mat <- rate[[t]]
        if (!is.matrix(mat) || nrow(mat) != L || ncol(mat) != L) {
          stop(sprintf("rate[[%d]] must be L x L matrix (L=%d)", t, L))
        }
        if (any(mat < -1e-8)) stop(sprintf("All entries in rate[[%d]] must be non-negative", t))

        row_sums <- rowSums(mat)
        if (any(abs(row_sums - 1) > 1e-6)) {
          warning(sprintf("Standardizing rows of rate[[%d]] to sum to 1", t))
          mat <- mat / row_sums
        }
        rate[[t]] <- mat
      }
    }
  }

  responses <- vector("list", times)
  Zs <- vector("list", times)
  P.Zs <- vector("list", times)

  position <- seq_len(L)
  if(use_covariates){
    P.Z.model <- mean_multinomial_logit_probability_cpp(
      covariates[[1]], beta_mat
    )
    if(is.sort){
      if(beta.supplied){
        .simulation.require.sorted(P.Z.model, TRUE, "beta")
      }else{
        position <- order(P.Z.model, decreasing = TRUE)
      }
    }
    Zs[[1]] <- sample_multinomial_logit_cpp(covariates[[1]], beta_mat)
    if(!identical(position, seq_len(L))){
      Zs[[1]] <- match(Zs[[1]], position)
      beta_mat <- .simulation.relabel.coefficients(
        beta_mat, position, ref.class
      )
    }
  }else if(!is.null(params) && !is.null(params$Z)){
    Zs[[1]] <- as.integer(params$Z)
    if(length(Zs[[1]]) != N || !all(Zs[[1]] %in% seq_len(L))){
      stop("params$Z must be length N with values in 1:L")
    }
    .simulation.require.sorted(
      .class.proportions(Zs[[1]], L), is.sort, "params$Z"
    )
  }else{
    if(!is.null(params) && !is.null(params$P.Z)){
      P.Z1 <- .simulation.probability(params$P.Z, L, "params$P.Z")
      .simulation.require.sorted(P.Z1, is.sort, "params$P.Z")
    }else if(distribution == "uniform"){
      P.Z1 <- rep(1 / L, L)
    }else if(distribution == "random"){
      P.Z1 <- as.numeric(rdirichlet(1L, rep(3, L)))
      if(is.sort) P.Z1 <- sort(P.Z1, decreasing = TRUE)
    }else{
      stop("distribution must be 'uniform' or 'random'")
    }
    Zs[[1]] <- .simulation.sample.classes(P.Z1, N)
  }

  P.Zs[[1]] <- .class.proportions(Zs[[1]], L)

  if(times > 1){
    if(use_covariates){
      for (t in 2:times) {
        Zs[[t]] <- sample_transition_logit_cpp(
          Zs[[t-1]], covariates[[t]], gamma_list[[t-1]]
        )
      }
    } else {
      for (t in 2:times) {
        Zs[[t]] <- sample_markov_cpp(Zs[[t-1]], rate[[t-1]])
      }
    }

    for (t in 2:times) {
      P.Zs[[t]] <- .class.proportions(Zs[[t]], L)
    }
  }

  latent.group.names <- .latent.group.names(L, type)
  P.Zs <- lapply(P.Zs, function(P.Z){
    names(P.Z) <- latent.group.names
    P.Z
  })
  if(!use_covariates && length(rate)){
    rate <- lapply(rate, function(rate.cur){
      dimnames(rate.cur) <- list(latent.group.names, latent.group.names)
      rate.cur
    })
  }

  params_t1 <- params
  if (is.null(params_t1)) params_t1 <- list()
  params_t1$Z <- Zs[[1]]

  if (type == "LCA") {
    data.obj.cur <- sim.LCA(N=N, I=I, L=L, poly.value=poly.value, IQ=IQ,
                            distribution=distribution, params=params_t1, is.sort=FALSE)
    par <- data.obj.cur$par
    poly.value <- data.obj.cur$poly.value
    means <- covs <- NULL
  } else if (type == "LPA") {
    data.obj.cur <- sim.LPA(N=N, I=I, L=L, constraint=constraint, distribution=distribution,
                            mean.range=mean.range, covs.range=covs.range, params=params_t1, is.sort=FALSE)
    means <- data.obj.cur$means
    covs <- data.obj.cur$covs
    par <- NULL
    poly.value <- NULL
  } else {
    stop("type must be 'LCA' or 'LPA'")
  }
  responses[[1]] <- data.obj.cur$response

  if(times > 1){
    for (t in 2:times) {
      Z.cur <- Zs[[t]]

      if (type == "LCA") {
        params.new <- list(par = par, Z = Z.cur, P.Z = NULL)
      } else {
        params.new <- list(means = means, covs = covs, Z = Z.cur, P.Z = NULL)
      }
      updated.obj <- update(data.obj.cur, params = params.new, is.sort=FALSE)
      responses[[t]] <- updated.obj$response
    }
  }

  names(Zs) <- names(P.Zs) <- names(responses) <- paste0("t", 1:times)

  res <- list(
    responses = responses,
    Zs = Zs,
    P.Zs = P.Zs,
    par = if (type == "LCA") par else NULL,
    means = if (type == "LPA") means else NULL,
    covs = if (type == "LPA") covs else NULL,
    poly.value = if (type == "LCA") poly.value else NULL,
    rate = if (!use_covariates) rate else NULL,
    covariates = if (use_covariates) covariates else NULL,
    beta = if (use_covariates) beta_mat else NULL,
    gamma = if (use_covariates && times >= 2) gamma_list else NULL,
    ref.class = ref.class,
    call = call,
    arguments = list(
      N = N,
      I = I,
      L = L,
      times = times,
      type = type,
      distribution = distribution,
      constraint = constraint,
      poly.value = poly.value,
      IQ = IQ,
      mean.range = mean.range,
      covs.range = covs.range,
      params = params,
      is.sort = is.sort,
      rate = rate,
      covariates = covariates,
      ref.class = ref.class,
      beta = beta,
      gamma = gamma
    )
  )

  class(res) <- "sim.LTA"
  return(res)
}
