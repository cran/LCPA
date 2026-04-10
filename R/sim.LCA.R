#' Simulate Data for Latent Class Analysis
#'
#' Generates synthetic multivariate categorical data from a latent class model with \code{L} latent classes.
#' Each observed variable follows a multinomial distribution within classes, with flexible control over
#' class separation via the \code{IQ} parameter and class size distributions.
#'
#' @param N Integer; total number of observations to simulate. Must be \eqn{> L}. Default: \code{1000}.
#' @param I Integer; number of categorical observed variables. Must be \eqn{\geq 1}. Default: \code{10}.
#' @param L Integer; number of latent classes. Must be \eqn{\geq 2} when \code{IQ} is numeric. Default: \code{3}.
#' @param poly.value Integer or integer vector; number of categories (levels) for each observed variable.
#'   If scalar, all variables share the same number of categories. If vector, must have length \code{I}.
#'   Minimum valid value is \code{2} when \code{IQ} is numeric. Default: \code{5}.
#' @param IQ Character or numeric; controls category probability distributions:
#'   \describe{
#'     \item{\code{"random"}}{(default) Dirichlet-distributed probabilities (\eqn{\alpha=3}).}
#'     \item{\code{Numeric}}{in \eqn{(0.5, 1)}. Forces high discriminative power (see details in section below).}
#'   }
#' @param distribution Character; distribution of class sizes. Options: \code{"random"} (default) or \code{"uniform"}.
#' @param params List with fixed parameters for simulation:
#'   \describe{
#'     \item{\code{par}}{\eqn{L \times I \times K_{\max}} array of conditional response probabilities per latent class.}
#'     \item{\code{P.Z}}{Vector of length \eqn{L} with latent class prior probabilities.}
#'     \item{\code{Z}}{Vector of length \eqn{N} containing the latent classes of observations. A fixed
#'                     observation classes \code{Z} is applied directly to simulate data only when \code{P.Z}
#'                     is \code{NULL} and \code{Z} is a \code{N} length vector.}
#'   }
#' @param is.sort A logical value. If \code{TRUE} (Default), the latent classes will be ordered in descending
#'                order according to \code{P.Z}. All other parameters will be adjusted accordingly
#'                based on the reordered latent classes.
#'
#' @return A list containing:
#'   \describe{
#'     \item{response}{Integer matrix (\eqn{N \times I}) of simulated observations. Rows are observations (named \code{"O1"}, \code{"O2"}, ...),
#'       columns are variables named \code{"I1"}, \code{"I2"}, ... Values range from \code{0} to \code{poly.value[i]-1}.}
#'     \item{par}{Array (\eqn{L \times I \times K}) of true class-specific category probabilities,
#'       where \eqn{K = \text{max}(poly.value)} (i.e., the maximum number of categories across variables).
#'       Dimensions: classes x variables x categories.
#'       Note: For variables with \code{poly.value[i] < K}, unused category dimensions contain \code{NA}.
#'       Dimension names: \code{"L1"}, \code{"L2"}, ... (classes); \code{"I1"}, \code{"I2"}, ... (variables);
#'       \code{"poly0"}, \code{"poly1"}, ... (categories).}
#'     \item{Z}{Integer vector (length \eqn{N}) of true class assignments (1 to L). Named with observation IDs (e.g., \code{"O1"}).}
#'     \item{P.Z}{Numeric vector (length \eqn{L}) of true class proportions. Named with class labels (e.g., \code{"L1"}).}
#'     \item{poly.value}{Integer vector (length \eqn{I}) specifying number of categories per variable.}
#'     \item{P.Z.Xn}{Binary matrix (\eqn{N \times L}) of true class membership indicators (one-hot encoded).
#'       Row \code{i}, column \code{l} = 1 if observation \code{i} belongs to class \code{l}, else 0.
#'       Row/column names match \code{Z} and class labels.}
#'     \item{arguments}{A list containing all input arguments.}
#'   }
#'
#' @section Indicator Quality (IQ) Parameter:
#' Controls the discriminative power of observed variables:
#' \describe{
#'   \item{\code{IQ = "random"}}{(Default)
#'     Category probabilities for each variable-class combination are drawn from a symmetric Dirichlet distribution
#'     (\eqn{\alpha = 3}), resulting in moderate class separation.
#'   }
#'   \item{\code{IQ = numeric}}{ (0.5 < IQ < 1)
#'     Forces high discriminative power for each variable:
#'     \enumerate{
#'       \item For each variable, two categories per class are assigned extreme probabilities:
#'         one category gets probability \eqn{IQ}, another gets \eqn{1-IQ}.
#'       \item Remaining categories share the residual probability \eqn{1 - IQ - (1-IQ) = 0}.
#'         \emph{Note: This requires \code{poly.value} >= 2 for all variables.}
#'       \item Category assignments are randomized within classes to avoid structural patterns.
#'     }
#'     Higher \code{IQ} values (closer to 1) yield stronger class separation but increase simulation failure risk.
#'   }
#' }
#'
#' @section Class Size Distribution:
#' \describe{
#'   \item{\code{"random"}}{(Default) Class proportions drawn from Dirichlet distribution (\eqn{\alpha = 3} for all classes),
#'     ensuring no empty classes. Sizes are rounded to integers with adjustment for exact \code{N}.}
#'   \item{\code{"uniform"}}{Equal probability of class membership (\eqn{1/L} per class), sampled with replacement.
#'     May produce empty classes if \code{N} is small relative to \code{L}.}
#' }
#'
#' @section Response Validation:
#' The simulation enforces a critical constraint: \emph{every category of every observed variable must appear
#' at least once in the dataset}. If initial generation violates this (e.g., a rare category is missing),
#' parameters and responses are regenerated until satisfied. This ensures compatibility with standard LCA estimation.
#'
#' @details
#' \strong{Probability Generation:}
#' \itemize{
#'   \item Dirichlet Sampling (\code{IQ="random"}):
#'     For each variable-class combination, probabilities are drawn from
#'     \eqn{\text{Dirichlet}(\alpha_1=3, \dots, \alpha_k=3)} where \eqn{k = \text{poly.value}[i]}.
#'   \item High-Discrimination Mode (\code{IQ=numeric}):
#'     For each variable \code{i}:
#'     \enumerate{
#'       \item Generate special probabilities \code{par.special} of length \code{L} containing:
#'         \eqn{IQ}, \eqn{1-IQ}, and \eqn{L-2} values uniformly sampled from \eqn{[1-IQ, IQ]}.
#'       \item For each class \code{l}, assign \code{par.special[l]} to one category, distribute
#'         remaining probability \eqn{1 - \text{par.special}[l]} uniformly (via Dirichlet) across other categories.
#'       \item Shuffle category assignments to avoid position bias.
#'     }
#' }
#'
#' \strong{Data Generation:}
#' \itemize{
#'   \item Class assignments \code{Z} are generated first according to \code{distribution}.
#'   \item For each observation \code{p} and variable \code{i}:
#'     \enumerate{
#'       \item Retrieve cumulative probabilities for class \code{Z[p]}
#'       \item Draw uniform random number \eqn{u \sim \text{Uniform}(0,1)}
#'       \item Assign category \code{k} where \eqn{P(\text{category} \leq k-1) < u \leq P(\text{category} \leq k)}
#'     }
#'   \item Entire dataset is regenerated if any category of any variable has zero observations.
#' }
#'
#' \strong{Critical Constraints:}
#' \itemize{
#'   \item When \code{IQ} is numeric: \eqn{0.5 < IQ < 1} and \code{min(poly.value) >= 2}
#'   \item \code{N} must be sufficiently large to observe all categories, especially when \code{IQ} is high
#'     or \code{poly.value} is large. Simulation may fail for small \code{N}.
#'   \item For \code{distribution="uniform"}, empty classes possible when \eqn{N < L}.
#' }
#'
#' @importFrom stats runif
#'
#' @examples
#' # Example 1: Default settings (moderate separation, random class sizes)
#' sim_data <- sim.LCA(N = 30, I = 5, L = 3)
#'
#' # Example 2: High-discrimination indicators (IQ=0.85), uniform class sizes
#' sim_high_disc <- sim.LCA(
#'   N = 30,
#'   I = 4,
#'   L = 2,
#'   poly.value = c(3,4,3,5),  # Variable category counts
#'   IQ = 0.85,
#'   distribution = "uniform"
#' )
#'
#' # Example 3: Binary indicators (poly.value=2) with high separation
#' sim_binary <- sim.LCA(N = 300, I = 10, L = 2, poly.value = 2, IQ = 0.9)
#'
#' @export
sim.LCA <- function(N=1000, I=10, L=3, poly.value=5, IQ="random", distribution="random", params=NULL, is.sort=TRUE){

  call <- match.call()

  if(length(poly.value) != I){
    poly.value <- rep(poly.value[1], I)
  }
  poly.max <- max(poly.value)

  is.availiable <- FALSE

  while(!is.availiable){
    if(!is.null(params$Z)){
      Z <- params$Z
      P.Z <- table(params$Z) / sum(Z)
    }else if(!is.null(params$P.Z)) {
      P.Z <- params$P.Z
      Z.frequency <- ceiling(P.Z*N)
      l.max <- which.max(Z.frequency)
      Z.frequency[l.max] <- N - sum(Z.frequency[c(1:L)[-l.max]])
      Z <- NULL
      for(l in 1:L){ Z <- c(Z, rep(l, Z.frequency[l])) }
      Z <- sample(Z, N, replace = FALSE)
    }else{
      if(distribution == "random"){
        Z.frequency <- ceiling(rdirichlet(n = 1, alpha = rep(3, L))*N)

        l.max <- which.max(Z.frequency)
        Z.frequency[l.max] <- N - sum(Z.frequency[c(1:L)[-l.max]])
        Z <- NULL
        for(l in 1:L){ Z <- c(Z, rep(l, Z.frequency[l])) }
        Z <- sample(Z, N, replace = FALSE)
      }else if(distribution == "uniform"){
        Z <- sample(1:L, N, replace = TRUE)
      }
      P.Z <- rep(0, L)
      P.Z[as.numeric(names(table(Z)))] <- table(Z) / N
    }

    if(is.null(params$par)){
      par <- array(NA, dim=c(L, I, poly.max))
      if(is.character(IQ)){
        if(IQ == "random"){
          for(i in 1:I){
            for(l in 1:L){
              par[l , i, 1:poly.value[i]] <- rdirichlet(n=1, alpha = rep(3, poly.value[i]))
            }
          }
        }
      }else if(is.numeric(IQ)){
        for(i in 1:I){
          par.special <- sample(c(runif(L-2, 1-IQ, IQ), IQ, 1-IQ), L, replace = FALSE)
          for(l in 1:L){
            par[l , i, 1:(poly.value[i]-1)] <- rdirichlet(n=1, alpha=rep(3, (poly.value[i]-1))) * (1 - par.special[l])
            par[l , i, poly.value[i]] <- par.special[l]
            par[l , i, 1:poly.value[i]] <- sample(par[l , i, 1:poly.value[i]], poly.value[i], replace = FALSE)
          }
        }
      }
    }else{
      par <- params$par
    }

    LCA.R <- function(Z, par){
      P.LCA <- par[Z, , ]
      return(P.LCA)
    }

    P.LCA <- LCA.R(Z, par)
    P.random <- P.LCA.cumulative <- array(0, dim=c(N, I, poly.max))
    for(p in 1:N){
      for(i in 1:I){
        for(po in 1:poly.value[i]){
          P.LCA.cumulative[p , i, po] <- sum(P.LCA[p , i, 1:po])
        }
      }
    }

    P.random[, , 1] <- matrix(runif(N*I, 0, 1), N, I)
    for(po in 2:poly.max){
      P.random[ , , po] <- P.random[ , , 1]
    }
    response <- matrix(0, N, I)
    for(i in 1:I){
      temp <- P.random[ , i, 1:poly.value[i]] >= P.LCA.cumulative[ , i, 1:poly.value[i]]
      response[ , i] <- apply(temp, 1, sum)
    }

    is.availiable <- check.response(response, poly.value)
  }

  P.Z.Xn <- t(sapply(1:N, function(i){
    temp <- rep(0, L)
    temp[Z[i]] <-  1
    temp
  }))

  if (is.sort) {
    posi <- order(P.Z, decreasing = TRUE)
    P.Z     <- P.Z[posi]
    par     <- par[posi, , , drop = FALSE]
    P.Z.Xn  <- P.Z.Xn[, posi, drop = FALSE]
    Z <- match(Z, posi)
  }

  dimnames(par) <- list(paste0("Class.", 1:L),
                        paste0("I", 1:I),
                        paste0("poly", 0:(poly.max-1)))
  colnames(P.Z.Xn) <- paste0("Class.", 1:L)
  names(P.Z) <- paste0("Class.", 1:L)
  names(Z) <- paste0("O", 1:N)

  res <- list(response=response, par=par, Z=Z, P.Z=P.Z,
              poly.value=poly.value, P.Z.Xn=P.Z.Xn)
  res$call <- call
  res$arguments = list(
    N=N, I=I, L=L, poly.value=poly.value, IQ=IQ,
    distribution=distribution, params=params, is.sort=is.sort
  )

  class(res) <- "sim.LCA"

  return(res)
}
