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
#' @param distribution Character; distribution of initial class probabilities when not using covariates or \code{params}.
#'   Options: \code{"uniform"} (equal probabilities) or \code{"random"} (Dirichlet-distributed, default).
#' @param times Integer; number of time points (must be \eqn{\geq 1}).
#' @param type Character; type of latent model. \code{"LCA"} for categorical indicators (default),
#'   \code{"LPA"} for continuous indicators.
#' @param rate List of matrices or NULL; transition probability matrices for non-covariate mode.
#'   Each matrix is \eqn{L \times L} with rows summing to 1. If \code{NULL} (default), matrices are
#'   generated with 0.7 diagonal probability and uniform off-diagonals. **Ignored when \code{times=1}**.
#' @param constraint Character; covariance structure for LPA (\code{type="LPA"} only). Options:
#'   \code{"VV"} (unstructured, default), \code{"VE"} (diagonal variance), \code{"EE"} (equal variance).
#' @param mean.range Numeric vector; range for randomly generated class means in LPA (default: \code{c(-2, 2)}).
#' @param covs.range Numeric vector; range for covariance matrix diagonals in LPA (default: \code{c(0.01, 4)}).
#' @param poly.value Integer; number of categories for polytomous LCA indicators (default: 5).
#' @param IQ Character; method for generating indicator discrimination in LCA. \code{"random"} (default) or fixed values.
#' @param params List or NULL; pre-specified parameters for reproducibility (see Details).
#' @param is.sort A logical value. If \code{TRUE} (Default), the latent classes will be ordered in descending
#'                order according to \code{P.Z}. All other parameters will be adjusted accordingly
#'                based on the reordered latent classes.
#' @param covariates List of matrices or NULL; covariate matrices for each time point. Each matrix must have
#'   dimensions \eqn{N \times p_t} and include an intercept column (first column must be all 1s). If \code{NULL},
#'   covariate mode is disabled. See Details for automatic coefficient generation.
#' @param beta Matrix or NULL; initial state regression coefficients of dimension \eqn{p_1 \times L}.
#'   Columns correspond to classes 1 to \eqn{L} (last class \eqn{L} is reference and must be zero).
#'   If \code{NULL} and covariates are used, coefficients are randomly generated from \eqn{\text{Uniform}(-1, 1)}.
#' @param gamma List or NULL; transition regression coefficients. Must be a list of length \code{times-1}.
#'   Each element \eqn{t} is a list of length \eqn{L} (previous state). Each sub-list contains \eqn{L} vectors
#'   (next state), where the last vector (reference class) is always \eqn{\mathbf{0}}. **Ignored when \code{times=1}**.
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
#'   \item{\code{rate}}{True transition matrices (non-covariate mode only; \code{NULL} when \code{times=1}).}
#'   \item{\code{covariates}}{List of covariate matrices used (covariate mode only).}
#'   \item{\code{beta}}{True initial state coefficients (covariate mode only).}
#'   \item{\code{gamma}}{True transition coefficients (covariate mode only; \code{NULL} when \code{times=1}).}
#'   \item{\code{call}}{Function call.}
#'   \item{\code{arguments}}{Input arguments.}
#' }
#'
#' @section Model Specification:
#' \describe{
#'   \item{Initial Class Probabilities (with covariates):}{
#'     For observation/participant \eqn{n} at time 1, the probability of belonging to latent class \eqn{l} is:
#'     \deqn{P(Z_{n1} = l \mid \mathbf{X}_{n1}) =
#'       \frac{\exp(\boldsymbol{\beta}_l^\top \mathbf{X}_{n1})}
#'            {\sum_{k=1}^L \exp(\boldsymbol{\beta}_k^\top \mathbf{X}_{n1})}}
#'     where \eqn{\mathbf{X}_{n1} = (X_{n10}, X_{n11}, \dots, X_{n1M})^\top} is the covariate vector for observation/participant \eqn{n} at time 1,
#'     with \eqn{X_{n10} = 1} (intercept term) and \eqn{X_{n1m}} (\eqn{m=1,\dots,M}) representing the value of the \eqn{m}-th covariate.
#'     The coefficient vector \eqn{\boldsymbol{\beta}_l = (\beta_{l0}, \beta_{l1}, \dots, \beta_{lM})^\top} corresponds element-wise to \eqn{\mathbf{X}_{n1}},
#'     where \eqn{\beta_{l0}} is the intercept and \eqn{\beta_{lm}} (\eqn{m \geq 1}) are regression coefficients for covariates.
#'     Class \eqn{L} is the reference class (\eqn{\boldsymbol{\beta}_L = \mathbf{0}}).
#'   }
#'   \item{Transition Probabilities (with covariates and times>=2):}{
#'     For observation/participant \eqn{n} transitioning from class \eqn{l} at time \eqn{t-1} to class \eqn{k} at time \eqn{t} (\eqn{t \geq 2}):
#'     \deqn{P(Z_{nt} = k \mid Z_{n,t-1} = l, \mathbf{X}_{nt}) =
#'       \frac{\exp(\boldsymbol{\gamma}_{lkt}^\top \mathbf{X}_{nt})}
#'            {\sum_{j=1}^L \exp(\boldsymbol{\gamma}_{ljt}^\top \mathbf{X}_{nt})}}
#'     where \eqn{\mathbf{X}_{nt} = (X_{nt0}, X_{nt1}, \dots, X_{ntM})^\top} is the covariate vector at time \eqn{t},
#'     with \eqn{X_{nt0} = 1} (intercept) and \eqn{X_{ntm}} (\eqn{m=1,\dots,M}) as the \eqn{m}-th covariate value.
#'     The coefficient vector \eqn{\boldsymbol{\gamma}_{lkt} = (\gamma_{lkt0}, \gamma_{lkt1}, \dots, \gamma_{lktM})^\top}
#'     corresponds element-wise to \eqn{\mathbf{X}_{nt}}, where \eqn{\gamma_{lkt0}} is the intercept and \eqn{\gamma_{lktm}} (\eqn{m \geq 1})
#'     are regression coefficients. Class \eqn{L} is the reference class (\eqn{\boldsymbol{\gamma}_{lLt} = \mathbf{0}} for all \eqn{l}).
#'   }
#'   \item{Without Covariates or When times=1:}{
#'     Initial probabilities follow a multinomial distribution with probabilities \eqn{\boldsymbol{\pi} = (\pi_1, \dots, \pi_L)}.
#'     When \eqn{times \geq 2}, transitions follow a Markov process with fixed probabilities \eqn{\tau_{lk}^{(t)} = P(Z_t = k \mid Z_{t-1} = l)},
#'     where \eqn{\sum_{k=1}^L \tau_{lk}^{(t)} = 1} for each \eqn{l} and \eqn{t}.
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
#'   \item The reference class (\eqn{L}) always has zero coefficients (\eqn{\boldsymbol{\beta}_L = \mathbf{0}},
#'     \eqn{\boldsymbol{\gamma}_{l,L} = \mathbf{0}}).
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
#' beta <- matrix(c( 0.8, -0.5, 0.0,
#'                  -0.3, -0.4, 0.0,
#'                   0.2,  0.8, 0.0,
#'                  -0.1,  0.2, 0.0), ncol=3, byrow=TRUE)
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
#'   covariates=covariates,
#'   beta=beta,
#'   gamma=gamma
#' )
#'
#' summary(sim_custom)
#'
#' @export
sim.LTA <- function(N=500, I=5, L=3, distribution="random",
                    times=2, type="LCA", rate=NULL,
                    constraint = "VV", mean.range = c(-2, 2), covs.range = c(0.01, 4),
                    poly.value=5, IQ="random",
                    params=NULL, is.sort=TRUE,
                    covariates = NULL,
                    beta = NULL,
                    gamma = NULL) {

  call <- match.call()

  if (times < 1) stop("times must be at least 1")
  if (L < 2) stop("L must be at least 2")

  use_covariates <- !is.null(covariates)

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
      beta_mat[, 1:(L-1)] <- matrix(runif(p1 * (L-1), -1, 1), p1, L-1)
    } else {
      if (!is.matrix(beta) || nrow(beta) != p1 || ncol(beta) != L) {
        stop(sprintf("beta must be %d x %d matrix (p1 x L)", p1, L))
      }
      beta_mat <- beta
    }

    # Modified block: Handle gamma only when times >= 2
    gamma_list <- NULL  # Default to NULL
    if (times >= 2) {
      gamma_list <- vector("list", times-1)
      for (t in 1:(times-1)) {
        pt <- ncol(covariates[[t+1]])
        if (is.null(gamma)) {
          gamma_t <- vector("list", L)
          for (l in 1:L) {
            gamma_t[[l]] <- vector("list", L)
            for (k in 1:L) {
              gamma_t[[l]][[k]] <- if (k < L) runif(pt, -1, 1) else rep(0, pt)
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
              gamma_valid[[l]][[k]] <- coef_vec
            }
          }
          gamma_list[[t]] <- gamma_valid
        }
      }
    } else {
      # times == 1: Explicitly ignore gamma
      if (!is.null(gamma)) {
        warning("gamma parameter is ignored when times=1 in covariate mode")
      }
    }
    # End modified block

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

      rate.accumulate <- vector("list", times-1)
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
        rate.accumulate[[t]] <- t(apply(mat, 1, cumsum))
      }
    }
  }

  responses <- vector("list", times)
  Zs <- vector("list", times)
  P.Zs <- vector("list", times)

  if (use_covariates) {
    X1 <- covariates[[1]]
    prob1 <- matrix(0, N, L)
    for (l in 1:L) {
      eta <- X1 %*% beta_mat[, l]
      prob1[, l] <- exp(eta)
    }
    prob1 <- prob1 / rowSums(prob1)
    Zs[[1]] <- apply(prob1, 1, function(p) sample(1:L, 1, prob = p))

    # When times=1, this loop is skipped automatically
    if(times > 1){
      for (t in 2:times) {
        Z_prev <- Zs[[t-1]]
        Z_cur <- integer(N)
        Xt <- covariates[[t]]
        gamma_t <- gamma_list[[t-1]]

        for (n in 1:N) {
          l_prev <- Z_prev[n]
          eta_vec <- numeric(L)
          for (k in 1:L) {
            coef_vec <- gamma_t[[l_prev]][[k]]
            eta_vec[k] <- sum(coef_vec * Xt[n, ])
          }
          prob_vec <- exp(eta_vec) / sum(exp(eta_vec))
          Z_cur[n] <- sample(1:L, 1, prob = prob_vec)
        }
        Zs[[t]] <- Z_cur
      }
    }
  } else {
    if (!is.null(params) && !is.null(params$Z)) {
      Zs[[1]] <- params$Z
      if (length(Zs[[1]]) != N || !all(Zs[[1]] %in% 1:L)) {
        stop("params$Z must be length N with values in 1:L")
      }
      if (use_covariates) {
        warning("params$Z ignored in covariate mode; states generated from covariates")
      }
    } else {
      if (distribution == "uniform") {
        P.Z1 <- rep(1/L, L)
      } else if (distribution == "random") {
        P.Z1 <- rgamma(L, 3); P.Z1 <- P.Z1 / sum(P.Z1)
      } else {
        stop("distribution must be 'uniform' or 'random'")
      }
      Zs[[1]] <- sample(1:L, N, replace = TRUE, prob = P.Z1)
    }

    if(times > 1){
      for (t in 2:times) {
        Z_prev <- Zs[[t-1]]
        Z_cur <- integer(N)

        for (l in 1:L) {
          idx <- which(Z_prev == l)
          if (length(idx) > 0) {
            u <- runif(length(idx))
            breaks <- c(0, rate.accumulate[[t-1]][l, 1:(L-1)], 1)
            Z_cur[idx] <- findInterval(u, breaks, rightmost.closed = TRUE)
          }
        }
        Zs[[t]] <- Z_cur
      }
    }
  }

  for (t in 1:times) {
    P.Zs[[t]] <- as.numeric(table(factor(Zs[[t]], levels = 1:L))) / N
  }

  if(is.sort){
    posi <- order(P.Zs[[1]], decreasing = TRUE)
    P.Zs[[1]]     <- P.Zs[[1]][posi]
    Zs[[1]] <- match(Zs[[1]], posi)
    beta <- beta[, posi]
    if(times > 1){
      for(t in 2:times){
        P.Zs[[t]]     <- P.Zs[[t]][posi]
        Zs[[t]] <- match(Zs[[t]], posi)
        gamma.temp <- gamma
        for(l in 1:L){
          for(ll in 1:L){
            gamma[[t-1]][[l]][[ll]] <- gamma.temp[[t-1]][[ posi[l] ]][[ posi[ll] ]]
          }
        }
      }
    }
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
    call = call,
    arguments = list(
      N = N, I = I, L = L, distribution = distribution,
      times = times, type = type, rate = rate,
      constraint = constraint, mean.range = mean.range, covs.range = covs.range,
      poly.value = poly.value, IQ = IQ, params = params, is.sort=is.sort,
      covariates = covariates, beta = beta, gamma = gamma
    )
  )

  class(res) <- "sim.LTA"
  return(res)
}
