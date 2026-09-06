#' Simulate Data for Latent Profile Analysis
#'
#' Generates synthetic multivariate continuous data from a latent profile model with \code{L} latent classes.
#' Supports flexible covariance structure constraints (including custom equality constraints) and
#' class size distributions. All covariance matrices are ensured to be positive definite.
#'
#' @param N Integer; total number of observations to simulate. Must be \eqn{\geq} \code{L} (Default = 1000).
#' @param I Integer; number of continuous observed variables. Must be \eqn{\geq 1} (Default = 5).
#' @param L Integer; number of latent profiles (classes). Must be \eqn{\geq 1} (Default = 2).
#' @param constraint Character string or list specifying covariance constraints. See detailed description below.
#'   Default is \code{"VV"} (fully heterogeneous covariances).
#' @param distribution Character; distribution of class sizes. Options: \code{"random"} (default) or \code{"uniform"}.
#' @param mean.range Numeric vector of length 2; range for sampling class-specific means.
#'   Each variable's means are sampled uniformly from \code{mean.range[1]} to \code{mean.range[2]}.
#'   Default: \code{c(-2, 2)}.
#' @param covs.range Numeric vector of length 2; range for sampling variance parameters (diagonal elements).
#'   Must satisfy \code{covs.range[1] > 0} and \code{covs.range[2] > covs.range[1]}. Off-diagonal covariances
#'   are derived from correlations scaled by these variances. Default: \code{c(0.01, 4)}.
#' @param params List with fixed parameters for simulation:
#'   \describe{
#'     \item{\code{means}}{\eqn{L \times I} matrix of class-specific means in the final class order.}
#'     \item{\code{covs}}{\eqn{I \times I \times L} array of class-specific
#'       covariance matrices in the final class order.}
#'     \item{\code{P.Z}}{Vector of length \eqn{L} with latent class prior probabilities in the final class order.}
#'     \item{\code{Z}}{Vector of length \eqn{N} containing fixed class
#'       assignments in the final class order. When supplied, \code{Z} takes precedence over
#'       \code{P.Z}.}
#'   }
#' @param is.sort A logical value. If \code{TRUE} (default), internally generated class probabilities
#'   are ordered decreasingly before class-specific parameters and observations are generated.
#'   Supplied parameters are already defined in this final order and are never reordered; supplied
#'   \code{P.Z} or \code{Z} must therefore follow decreasing class proportions.
#'
#' @return A list containing:
#'   \describe{
#'     \item{response}{Numeric matrix (\eqn{N \times I}) of simulated observations. Rows are observations,
#'       columns are variables named \code{"V1"}, \code{"V2"}, ..., or \code{"V"} for univariate data.}
#'     \item{means}{Numeric matrix (\eqn{L \times I}) of true class-specific means.
#'       Row names: `"Profile 1"`, `"Profile 2"`, and so on; column names
#'       match `response`.}
#'     \item{covs}{Array (\eqn{I \times I \times L}) of true class-specific covariance matrices.
#'       Dimensions: variables x variables x classes. Constrained parameters have identical values across class slices.
#'       Dimension names match \code{response} and class labels.}
#'     \item{P.Z.Xn}{Numeric matrix (\eqn{N \times L}) of true class membership probabilities (one-hot encoded).
#'       Row \code{i}, column \code{l} = 1 if observation \code{i} belongs to class \code{l}, else 0.
#'       Row names: \code{"O1"}, \code{"O2"}, ...; column names:
#'       `"Profile 1"`, `"Profile 2"`, and so on.}
#'     \item{P.Z}{Numeric vector (length \eqn{L}) of true class proportions.
#'       Named `"Profile 1"`, `"Profile 2"`, and so on.}
#'     \item{Z}{Integer vector (length \eqn{N}) of true class assignments (1 to L).
#'       Named with observation IDs (e.g., \code{"O1"}).}
#'     \item{constraint}{Original constraint specification (character string or list) passed to the function.}
#'     \item{call}{Matched simulation call.}
#'     \item{arguments}{List of the effective simulation arguments.}
#'   }
#'
#' @section Covariance Constraints:
#' The \code{constraint} parameter controls equality constraints on covariance parameters across classes:
#' \describe{
#'   \item{Predefined Constraints (Character Strings):}{
#'     \describe{
#'       \item{\code{"UE"} (Univariate only)}{Equal variance across all classes.}
#'       \item{\code{"UV"} (Univariate only)}{Varying variances across classes.}
#'       \item{\code{"E0"}}{Equal variances across classes, zero covariances (diagonal matrix with shared variances).}
#'       \item{\code{"V0"}}{Varying variances across classes, zero covariances (diagonal matrix with free variances).}
#'       \item{\code{"EE"}}{Equal full covariance matrix across all classes (homogeneous).}
#'       \item{\code{"EV"}}{Equal variances but varying covariances (equal diagonal, free off-diagonal).}
#'       \item{\code{"VE"}}{Varying variances but equal correlations (free diagonal, equal correlation structure).}
#'       \item{\code{"VV"}}{Varying full covariance matrices across classes (heterogeneous; default).}
#'     }
#'   }
#'   \item{Custom Constraints (List of integer vectors):}{
#'     Each element specifies a pair of variables whose covariance parameters are constrained equal across classes:
#'     \describe{
#'       \item{\code{c(i,i)}}{Constrains variance of variable \code{i} to be equal across all classes.}
#'       \item{\code{c(i,j)}}{Constrains covariance between variables \code{i} and \code{j} to be equal across all classes
#'         (symmetric: automatically includes \code{c(j,i)}).}
#'     }
#'     Unconstrained parameters vary freely, and the returned covariance matrices
#'     are positive definite.
#'     Critical requirements for custom constraints:
#'     \describe{
#'       \item{At least one variance must be unconstrained if any off-diagonal covariance is unconstrained.}{}
#'       \item{All indices must be between 1 and \code{I}.}{}
#'       \item{For univariate data (\code{I=1}), only \code{list(c(1,1))} is valid.}{}
#'     }
#'   }
#' }
#'
#' @section Class Size Distribution:
#' \describe{
#'   \item{\code{"random"}}{(Default) Class proportions drawn from Dirichlet distribution (\eqn{\alpha = 3} for all classes),
#'     ensuring no empty classes. Sizes are rounded to integers with adjustment for exact \code{N}.}
#'   \item{\code{"uniform"}}{Equal probability of class membership (\eqn{1/L} per class), sampled with replacement.}
#' }
#'
#' @details
#' Mean Generation: For each variable, \eqn{3L} candidate means are sampled uniformly from \code{mean.range}.
#' \eqn{L} distinct means are selected without replacement to ensure separation between classes.
#'
#' Covariance Generation:
#' \itemize{
#'   \item Univariate case (\code{I=1}): Constraints \code{"UE"} and \code{"UV"} are enforced automatically.
#'     Predefined constraints like \code{"E0"} map to \code{"UE"}.
#' }
#'
#' Class Assignment:
#' \itemize{
#'   \item \code{"random"}: Uses Dirichlet distribution (\eqn{\alpha = 3}) to avoid extremely small classes.
#'     Sizes are rounded and adjusted to sum exactly to \code{N}.
#'   \item \code{"uniform"}: Simple random sampling with equal probability. May produce empty classes if \code{N} is small.
#' }
#'
#' Data Generation: Observations are simulated using \code{mvtnorm::rmvnorm} per class.
#' Final data and class labels are shuffled to remove ordering artifacts.
#'
#' @importFrom Matrix nearPD
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats runif rnorm
#'
#' @examples
#' # Example 1: Bivariate data, 3 classes, heterogeneous covariances (default)
#' sim_data <- sim.LPA(N = 500, I = 2, L = 3, constraint = "VV")
#'
#' # Example 2: Univariate data, equal variances
#' # 'E0' automatically maps to 'UE' for I=1
#' sim_uni <- sim.LPA(N = 200, I = 1, L = 2, constraint = "E0")
#'
#' # Example 3: Custom constraints
#' # - Equal covariance between V1 and V2 across classes
#' # - Equal variance for V3 across classes
#' sim_custom <- sim.LPA(
#'   N = 300,
#'   I = 3,
#'   L = 4,
#'   constraint = list(c(1, 2), c(3, 3))
#' )
#'
#' # Example 4: VE constraint (varying variances, equal correlations)
#' sim_ve <- sim.LPA(N = 400, I = 3, L = 3, constraint = "VE")
#'
#' # Example 5: Uniform class sizes
#' sim_uniform <- sim.LPA(N = 300, I = 4, L = 5, distribution = "uniform")
#'
#' @export
sim.LPA <- function(N = 1000, I = 5, L = 2, constraint = "VV", distribution = "random",
                    mean.range = c(-2, 2), covs.range = c(0.01, 4), params = NULL, is.sort=TRUE) {

  call <- match.call()

  if (covs.range[1] <= 0) stop("covs.range[1] must be > 0 for positive definiteness")
  if (covs.range[2] <= covs.range[1]) stop("covs.range[2] must be > covs.range[1]")

  requested_constraint <- constraint
  custom_constraint <- FALSE

  if (is.list(constraint)) {
    custom_constraint <- TRUE
    if (!all(sapply(constraint, function(x) is.numeric(x) && length(x) == 2 && all(x == as.integer(x))))) {
      stop("When constraint is a list, each element must be a two-integer vector (e.g., list(c(1,2), c(3,3))).")
    }
    all_indices <- unlist(constraint)
    if (any(all_indices < 1) || any(all_indices > I)) {
      stop("Index in constraint list out of bounds (must be between 1 and I).")
    }
    if (I == 1) {
      if (!all(all_indices == 1)) {
        stop("For I=1, constraint list must only contain pairs of 1 (e.g., list(c(1,1))).")
      }
    }
  } else {
    if (I > 1 && any(constraint %in% c("UE", "UV"))) {
      stop("Single-variable constraints ('UE','UV') require I = 1")
    }
    if (I == 1 && !any(constraint %in% c("UE", "UV"))) {
      if (any(constraint %in% c("E0", "EE", "EV"))) {
        constraint <- "UE"
        message(sprintf("I=1: Mapping to univariate constraint 'UE' from '%s'", requested_constraint))
      } else {
        constraint <- "UV"
        message(sprintf("I=1: Mapping to univariate constraint 'UV' from '%s'", requested_constraint))
      }
    }
  }

  if (is.null(params$means)) {
    means <- matrix(0, L, I)
    for (i in 1:I) {
      means.temp <- runif(L * 3, min = mean.range[1], max = mean.range[2])
      posi <- sample(c(1:L) * 3, L, replace = FALSE)
      means[, i] <- means.temp[posi]
    }
  } else {
    means <- params$means
  }

  is.univariate <- (I == 1)

  generate_positive_definite_cov <- function(I, covs.range) {
    R <- sim.correlation(I)
    sds <- runif(I, sqrt(covs.range[1]), sqrt(covs.range[2]))
    D <- diag(sds)
    Sigma <- D %*% R %*% D
    return(Sigma)
  }

  if (is.null(params$covs)) {
    generate_covs <- function() {
      covs_attempt <- array(0, dim = c(I, I, L))
      dimnames(covs_attempt) <- list(NULL, NULL, .latent.group.names(L, "LPA"))
      is.univariate_local <- (I == 1)

      if (custom_constraint) {
        S0 <- generate_positive_definite_cov(I, covs.range)
        mask <- matrix(FALSE, I, I)
        for (pair in constraint) {
          i <- pair[1]
          j <- pair[2]
          if (i == j) {
            mask[i, i] <- TRUE
          } else {
            mask[i, j] <- TRUE
            mask[j, i] <- TRUE
          }
        }
        mask[lower.tri(mask)] <- t(mask)[lower.tri(mask)]
        diag_logical <- matrix(FALSE, I, I)
        diag(diag_logical) <- TRUE
        nonshared_offdiag_exists <- any(!mask & !diag_logical)
        nonshared_diag_indices_global <- which(!mask[cbind(1:I, 1:I)])

        if (nonshared_offdiag_exists && length(nonshared_diag_indices_global) == 0) {
          stop("Cannot have non-shared covariances without at least one non-shared variance for adjustment.")
        }

        for (l in 1:L) {
          M <- S0

          if (I > 1) {
            upper_tri_nonshared <- which(!mask & upper.tri(mask, diag = FALSE), arr.ind = TRUE)
            if (nrow(upper_tri_nonshared) > 0) {
              factors <- runif(nrow(upper_tri_nonshared), 0.9, 1.1)
              for (k in 1:nrow(upper_tri_nonshared)) {
                i <- upper_tri_nonshared[k, 1]
                j <- upper_tri_nonshared[k, 2]
                factor <- factors[k]
                M[i, j] <- M[i, j] * factor
                M[j, i] <- M[j, i] * factor
              }
            }
          }

          diag_indices <- 1:I
          nonshared_diag_indices <- diag_indices[!mask[cbind(diag_indices, diag_indices)]]

          maxattempts.adjust <- 1000
          attempt <- 0
          while (!is.positive.definite(M) && attempt < maxattempts.adjust &&
                 length(nonshared_diag_indices) > 0) {
            projected <- as.matrix(Matrix::nearPD(
              M, corr = FALSE, keepDiag = FALSE,
              eig.tol = .Machine$double.eps^(2 / 3),
              posd.tol = sqrt(.Machine$double.eps)
            )$mat)
            adjustment <- pmax(diag(projected) - diag(M), 0)
            diag(M)[nonshared_diag_indices] <-
              diag(M)[nonshared_diag_indices] + adjustment[nonshared_diag_indices]
            attempt <- attempt + 1
          }

          if (!is.positive.definite(M)) {
            stop("Failed to make covariance matrix positive definite for class ", l,
                 " after ", maxattempts.adjust, " attempts. Consider wider covs.range or reducing I.")
          }
          covs_attempt[, , l] <- M
        }
      } else {
        if (is.univariate_local) {
          if (constraint == "UE") {
            var_shared <- runif(1, covs.range[1], covs.range[2])
            covs_attempt[1, 1, ] <- var_shared
          } else {
            covs_attempt[1, 1, ] <- runif(L, covs.range[1], covs.range[2])
          }
        } else {
          if (constraint == "EE") {
            shared_cov <- generate_positive_definite_cov(I, covs.range)
            for (l in 1:L) covs_attempt[, , l] <- shared_cov
          } else if (constraint == "E0") {
            shared_vars <- runif(I, covs.range[1], covs.range[2])
            shared_cov <- diag(shared_vars)
            for (l in 1:L) covs_attempt[, , l] <- shared_cov
          } else if (constraint == "VE") {
            maxattempts.ve <- 1000
            base_cov <- NULL
            non_diag_matrix <- NULL
            min_diag_val <- NULL

            for (attempt in 1:maxattempts.ve) {
              base_cov <- generate_positive_definite_cov(I, covs.range)
              non_diag_matrix <- base_cov
              diag(non_diag_matrix) <- 0
              e_vals <- eigen(non_diag_matrix, symmetric = TRUE, only.values = TRUE)$values
              lambda_min <- min(e_vals)
              numerical.margin <- sqrt(.Machine$double.eps) *
                max(1, max(abs(e_vals)))
              min_diag_val <- max(covs.range[1], -lambda_min + numerical.margin)

              if (min_diag_val <= covs.range[2]) break
              if (attempt == maxattempts.ve) {
                stop("Failed to generate VE constraint matrix within covs.range. ",
                     "Try increasing covs.range[2] or reducing dimension I.")
              }
            }

            for (l in 1:L) {
              diag_vals <- runif(I, min_diag_val, covs.range[2])
              cov_mat <- non_diag_matrix
              diag(cov_mat) <- diag_vals
              covs_attempt[, , l] <- cov_mat
            }
          } else if (constraint == "EV") {
            sds_shared <- sqrt(runif(I, covs.range[1], covs.range[2]))
            for (l in 1:L) {
              R_l <- sim.correlation(I)
              covs_attempt[, , l] <- diag(sds_shared) %*% R_l %*% diag(sds_shared)
            }
          } else if (constraint == "VV") {
            for (l in 1:L) {
              covs_attempt[, , l] <- generate_positive_definite_cov(I, covs.range)
            }
          } else if (constraint == "V0") {
            for (l in 1:L) {
              vars_l <- runif(I, covs.range[1], covs.range[2])
              covs_attempt[, , l] <- diag(vars_l)
            }
          }
        }
      }
      return(covs_attempt)
    }

    maxattempts.total <- 1000
    covs <- NULL
    for (attempt in 1:maxattempts.total) {
      covs_temp <- try(generate_covs(), silent = TRUE)
      if (!inherits(covs_temp, "try-error")) {
        covs <- covs_temp
        break
      } else {
        if (attempt < maxattempts.total) {
          message(sprintf("Attempt %d/%d: Failed to generate positive definite covariance matrices. Retrying...",
                          attempt, maxattempts.total))
        }
      }
    }
    if (is.null(covs)) {
      stop("Failed to generate positive definite covariance matrices after ", maxattempts.total, " attempts.")
    }
  } else {
    covs <- params$covs
  }

  if(!is.null(params$Z)){
    Z <- as.integer(params$Z)
    if(length(Z) != N || any(!Z %in% seq_len(L))){
      stop("params$Z must contain N class labels between 1 and L")
    }
    P.Z <- .class.proportions(Z, L)
    .simulation.require.sorted(P.Z, is.sort, "params$Z")
  }else if(!is.null(params$P.Z)){
    P.Z <- .simulation.probability(params$P.Z, L, "params$P.Z")
    .simulation.require.sorted(P.Z, is.sort, "params$P.Z")
    Z <- .simulation.sample.classes.from.pool(P.Z, N)
  }else{
    if (distribution == "random") {
      P.Z <- as.numeric(rdirichlet(1, rep(3, L)))
    } else if (distribution == "uniform") {
      P.Z <- rep(1 / L, L)
    } else {
      stop("Invalid 'distribution'. Choose 'random' or 'uniform'.")
    }
    if(is.sort) P.Z <- sort(P.Z, decreasing = TRUE)
    Z <- .simulation.sample.classes.from.pool(P.Z, N)
  }

  var_names <- if (is.univariate) "V" else paste0("V", 1:I)
  response <- matrix(0, nrow = N, ncol = I)
  colnames(response) <- var_names

  for (l in 1:L) {
    n_l <- sum(Z == l)
    if (n_l == 0) next

    covs.cur <- as.matrix(covs[, , l], I, I)
    mean_l <- means[l, , drop = FALSE]

    if (I > 1 && !is.positive.definite(covs.cur)) {
        nearPD_result <- Matrix::nearPD(
          covs.cur, corr = FALSE, keepDiag = FALSE,
          eig.tol = .Machine$double.eps^(2 / 3),
          posd.tol = sqrt(.Machine$double.eps)
        )
        covs.cur <- as.matrix(nearPD_result$mat)
        if(!is.positive.definite(covs.cur)){
          stop("Matrix::nearPD failed to return a positive-definite covariance matrix")
        }
        covs[, , l] <- covs.cur
    }

    class_data <- mvtnorm::rmvnorm(n = n_l, mean = as.numeric(mean_l), sigma = (covs.cur + t(covs.cur)) / 2)
    response[Z == l, ] <- class_data
  }
  P.Z.Xn <- .class.indicator(Z, L)
  P.Z <- .class.proportions(Z, L)


  colnames(P.Z.Xn) <- .latent.group.names(L, "LPA")
  rownames(P.Z.Xn) <- paste0("O", 1:N)
  names(P.Z) <- colnames(P.Z.Xn)
  rownames(means) <- .latent.group.names(L, "LPA")
  colnames(means) <- var_names
  dimnames(covs) <- list(var_names, var_names, .latent.group.names(L, "LPA"))
  names(Z) <- paste0("O", 1:N)

  res <- list(
    response = response,
    means = means,
    covs = covs,
    P.Z.Xn = P.Z.Xn,
    P.Z = P.Z,
    Z = Z,
    constraint = requested_constraint,
    call = call,
    arguments = list(
      N = N, I = I, L = L, constraint = constraint, distribution = distribution,
      mean.range = mean.range, covs.range = covs.range, params=params, is.sort=is.sort
    )
  )

  class(res) <- "sim.LPA"
  return(res)
}
