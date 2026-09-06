#' Calculate Log-Likelihood for Latent Profile Analysis
#'
#' Computes the log-likelihood of observed continuous data under a Latent Profile Analysis (LPA) model
#' with multivariate normal distributions within each latent profile.
#'
#' @param response A numeric matrix of dimension \eqn{N \times I} containing continuous observations.
#'   Rows represent observations, columns represent variables. Missing values are not permitted.
#' @param means A matrix of dimension \eqn{L \times I} where row \eqn{l} contains the mean vector
#'   \eqn{\boldsymbol{\mu}_l} for profile \eqn{l}.
#' @param covs An array of dimension \eqn{I \times I \times L} where slice \eqn{l} contains the
#'   covariance matrix \eqn{\boldsymbol{\Sigma}_l} for profile \eqn{l}. Must be symmetric positive semi-definite.
#' @param P.Z A numeric vector of length \eqn{L} containing prior probabilities for latent profiles.
#'   Must satisfy:
#'   \itemize{
#'     \item \eqn{\sum_{l=1}^L \pi_l = 1}
#'     \item \eqn{\pi_l > 0} for all \eqn{l = 1, \dots, L}
#'   }
#'
#' @return A single numeric value equal to the total observed-data
#'   log-likelihood \eqn{\log\mathcal{L}_{\mathrm{LPA}}} defined below.
#'
#' @details The log-likelihood calculation follows these steps:
#'
#' \itemize{
#'   \item Class-Conditional Likelihood Contribution:
#'   For participant \eqn{n} in profile \eqn{l}, the multivariate normal log
#'   density is
#'   \deqn{\log\mathcal{N}(\mathbf{X}_n\mid
#'   \boldsymbol{\mu}_l,\boldsymbol{\Sigma}_l) =
#'   -\frac{I}{2}\log(2\pi)-\frac{1}{2}\log|\boldsymbol{\Sigma}_l|
#'   -\frac{1}{2}(\mathbf{X}_n-\boldsymbol{\mu}_l)^\top
#'   \boldsymbol{\Sigma}_l^{-1}(\mathbf{X}_n-\boldsymbol{\mu}_l).}
#'   \item Participant-Level Marginal Likelihood Contribution:
#'   The marginal density for participant \eqn{n} is
#'   \deqn{P(\mathbf{X}_n)=\sum_{l=1}^L\pi_l
#'   \mathcal{N}(\mathbf{X}_n\mid
#'   \boldsymbol{\mu}_l,\boldsymbol{\Sigma}_l).}
#'   \item Total Observed-Data Log-Likelihood:
#'   \deqn{\log\mathcal{L}_{\mathrm{LPA}}=
#'   \sum_{n=1}^N\log\left\{\sum_{l=1}^L\pi_l
#'   \mathcal{N}(\mathbf{X}_n\mid
#'   \boldsymbol{\mu}_l,\boldsymbol{\Sigma}_l)\right\}.}
#' }
#'
#' @export
get.Log.Lik.LPA <- function(response, means, covs, P.Z) {
  if (!is.matrix(response) || !is.numeric(response) || any(!is.finite(response))) {
    stop("response must be a finite numeric matrix")
  }
  if (!is.numeric(P.Z) || any(!is.finite(P.Z)) || any(P.Z <= 0) ||
      abs(sum(P.Z) - 1) > 1e-5) {
    stop("P.Z must be a valid probability vector (positive values summing to 1)")
  }
  if (!is.matrix(means) || nrow(means) != length(P.Z) ||
      ncol(means) != ncol(response) || any(!is.finite(means))) {
    stop("means must be a matrix with rows corresponding to classes")
  }
  I <- ncol(response)
  L <- length(P.Z)

  if (!is.array(covs) || length(dim(covs)) != 3 ||
      dim(covs)[1] != I || dim(covs)[2] != I || dim(covs)[3] != L) {
    stop("covs must be an I x I x L array matching the dimensions of response and classes")
  }

  if(any(!is.finite(covs))){
    stop("covs must contain finite values")
  }
  stabilized <- .stabilize.LPA.covariances(
    covs, "VV", fallback = stats::cov(response)
  )
  if(stabilized$repaired) covs <- stabilized$covs
  expectation <- .lpa.expectation.reference(response, means, covs, P.Z)
  if(!isTRUE(expectation$valid)){
    stop("Unable to obtain positive-definite profile covariance matrices")
  }
  expectation$Log.Lik
}

.lpa.expectation.reference <- function(response, means, covs, P.Z){
  python.configured <- nzchar(Sys.getenv("RETICULATE_PYTHON")) ||
    reticulate::py_available(initialize = FALSE)
  if(!python.configured){
    return(lpa_expectation_cpp(
      response, means, covs, P.Z, repair = FALSE
    ))
  }
  torch <- tryCatch(
    reticulate::import("torch", convert = FALSE),
    error = function(e){ NULL }
  )
  if(is.null(torch)){
    return(lpa_expectation_cpp(
      response, means, covs, P.Z, repair = FALSE
    ))
  }

  result <- tryCatch({
    response.torch <- torch$tensor(
      reticulate::r_to_py(response), dtype = torch$float32
    )
    means.torch <- torch$tensor(
      reticulate::r_to_py(means), dtype = torch$float32
    )
    covs.torch <- torch$tensor(
      reticulate::r_to_py(covs), dtype = torch$float32
    )
    P.Z.torch <- torch$tensor(
      reticulate::r_to_py(as.numeric(P.Z)), dtype = torch$float32
    )
    N <- nrow(response)
    I <- ncol(response)
    L <- length(P.Z)
    eps <- 1e-6
    covariance.eps <- 1e-5
    dev <- response.torch$unsqueeze(1L) - means.torch$unsqueeze(0L)
    covs.batch <- covs.torch$permute(c(2L, 0L, 1L))
    off.diagonal <- covs
    for(l in seq_len(L)) diag(off.diagonal[, , l]) <- 0
    if(I == 1L || all(off.diagonal == 0)){
      variances <- torch$diagonal(
        covs.batch, dim1 = 1L, dim2 = 2L
      ) + covariance.eps * 10
      diag.chol <- torch$sqrt(variances)
      y <- dev$permute(c(1L, 2L, 0L)) / diag.chol$unsqueeze(2L)
    }else{
      eye <- torch$eye(I, dtype = torch$float32)$expand(c(L, I, I))
      covs.batch <- covs.batch + covariance.eps * 10 * eye
      chol.batch <- tryCatch(
        torch$linalg$cholesky(covs.batch),
        error = function(e){
          torch$linalg$cholesky(covs.batch + covariance.eps * 100 * eye)
        }
      )
      diag.chol <- torch$diagonal(chol.batch, dim1 = 1L, dim2 = 2L)
      y <- torch$linalg$solve_triangular(
        chol.batch, dev$permute(c(1L, 2L, 0L)), upper = FALSE
      )
    }
    logdet <- 2 * torch$sum(torch$log(diag.chol), dim = 1L)
    quad <- torch$sum(y^2, dim = 1L)
    constant <- -0.5 * I * log(2 * pi)
    log.pdf <- constant - 0.5 * (logdet$unsqueeze(1L) + quad)
    log.pdf <- torch$clamp(log.pdf, min = -1e8, max = 1e8)$t()
    log.joint <- log.pdf + torch$log(P.Z.torch + eps)
    row.max <- torch$max(log.joint, dim = 1L, keepdim = TRUE)$values
    row.max <- torch$where(
      torch$isfinite(row.max), row.max, torch$zeros_like(row.max)
    )
    relative <- torch$exp(log.joint - row.max)
    denominator <- torch$sum(relative, dim = 1L, keepdim = TRUE)
    all.zero <- denominator < eps
    denominator <- torch$where(
      all.zero, torch$ones_like(denominator) * eps, denominator
    )
    log.sum.exp <- row.max + torch$log(denominator)
    log.sum.exp <- torch$where(
      all.zero, torch$full_like(log.sum.exp, -Inf), log.sum.exp
    )
    log.sum.exp <- torch$where(
      torch$isfinite(log.sum.exp), log.sum.exp,
      torch$full_like(log.sum.exp, -1e8)
    )
    posterior <- relative / denominator
    list(
      posterior = reticulate::py_to_r(posterior$detach()$cpu()$numpy()),
      Log.Lik = reticulate::py_to_r(torch$sum(log.sum.exp)$item()),
      valid = TRUE
    )
  }, error = function(e){ NULL })
  if(!is.null(result)) return(result)
  lpa_expectation_cpp(response, means, covs, P.Z, repair = FALSE)
}

.stabilize.LPA.covariances <- function(covs, constraint, fallback = NULL){
  dimensions <- dim(covs)
  if(length(dimensions) != 3L || dimensions[1L] != dimensions[2L]){
    stop("covs must be an I x I x L array")
  }
  I <- dimensions[1L]
  L <- dimensions[3L]
  if(is.null(fallback)) fallback <- diag(I)
  fallback <- matrix(fallback, I, I)
  if(any(!is.finite(fallback))) fallback <- diag(I)
  fallback <- (fallback + t(fallback)) / 2
  minimum.eigenvalue <- function(covariance){
    if(any(!is.finite(covariance))) return(-Inf)
    values <- tryCatch(
      eigen((covariance + t(covariance)) / 2,
            symmetric = TRUE, only.values = TRUE)$values,
      error = function(e){ numeric(0) }
    )
    if(length(values) != I) -Inf else min(values)
  }

  is.positive.definite <- function(covariance){
    if(any(!is.finite(covariance))) return(FALSE)
    scale <- max(1, max(abs(covariance)))
    if(max(abs(covariance - t(covariance))) >
       sqrt(.Machine$double.eps) * scale) return(FALSE)
    !is.null(tryCatch(chol(covariance), error = function(e){ NULL }))
  }

  repair.one <- function(covariance){
    covariance <- matrix(covariance, I, I)
    if(any(!is.finite(covariance))) covariance <- fallback
    if(is.positive.definite(covariance)) return(covariance)
    covariance <- (covariance + t(covariance)) / 2
    covariance <- tryCatch(
      as.matrix(suppressWarnings(Matrix::nearPD(
        covariance, base.matrix = TRUE,
        eig.tol = .Machine$double.eps^(2 / 3),
        conv.tol = sqrt(.Machine$double.eps),
        posd.tol = sqrt(.Machine$double.eps)
      ))$mat),
      error = function(e){ fallback }
    )
    covariance <- (covariance + t(covariance)) / 2
    if(!is.positive.definite(covariance)){
      decomposition <- eigen(covariance, symmetric = TRUE)
      tolerance <- sqrt(.Machine$double.eps) *
        max(1, max(abs(decomposition$values)))
      covariance <- decomposition$vectors %*%
        diag(pmax(decomposition$values, tolerance), I) %*%
        t(decomposition$vectors)
      covariance <- (covariance + t(covariance)) / 2
    }
    covariance
  }

  original <- covs
  minimum.before <- vapply(seq_len(L), function(l){
    minimum.eigenvalue(covs[, , l])
  }, numeric(1))
  needs.repair <- !vapply(seq_len(L), function(l){
    is.positive.definite(covs[, , l])
  }, logical(1))
  if(!any(needs.repair)){
    return(list(
      covs = original,
      repaired = FALSE,
      diagnostics = list(
        minimum.eigenvalue.before = minimum.before,
        minimum.eigenvalue.after = minimum.before,
        repair.method = NULL
      )
    ))
  }
  for(l in seq_len(L)) covs[, , l] <- repair.one(covs[, , l])

  if(is.character(constraint)){
    if(constraint == "E0"){
      shared.variance <- vapply(seq_len(I), function(i){
        mean(covs[i, i, ])
      }, numeric(1))
      for(l in seq_len(L)) covs[, , l] <- diag(shared.variance, I)
    }else if(constraint == "V0"){
      for(l in seq_len(L)) covs[, , l] <- diag(diag(covs[, , l]), I)
    }else if(constraint %in% c("EE", "UE")){
      pooled <- apply(covs, c(1L, 2L), mean)
      for(l in seq_len(L)) covs[, , l] <- pooled
    }else if(constraint == "VE"){
      for(i in seq_len(I - 1L)){
        for(j in (i + 1L):I){
          value <- mean(covs[i, j, ])
          covs[i, j, ] <- covs[j, i, ] <- value
        }
      }
    }else if(constraint == "EV"){
      for(i in seq_len(I)) covs[i, i, ] <- mean(covs[i, i, ])
    }
  }else if(is.list(constraint)){
    for(pair in constraint){
      i <- as.integer(pair[1L])
      j <- as.integer(pair[2L])
      value <- mean(covs[i, j, ])
      covs[i, j, ] <- covs[j, i, ] <- value
    }
  }

  for(l in seq_len(L)) covs[, , l] <- repair.one(covs[, , l])
  list(
    covs = covs,
    repaired = TRUE,
    diagnostics = list(
      minimum.eigenvalue.before = minimum.before,
      minimum.eigenvalue.after = vapply(seq_len(L), function(l){
        minimum.eigenvalue(covs[, , l])
      }, numeric(1)),
      repair.method = "Matrix::nearPD"
    )
  )
}

.stabilize.LPA.result <- function(object, response, constraint){
  if(!is.list(object) || !is.list(object$params) ||
     is.null(object$params$covs) || is.null(object$params$means) ||
     is.null(object$params$P.Z)){
    return(object)
  }
  response <- as.matrix(response)
  stabilized <- .stabilize.LPA.covariances(
    object$params$covs, constraint,
    fallback = stats::cov(response)
  )
  if(!stabilized$repaired) return(object)

  P.Z <- pmax(as.numeric(object$params$P.Z), 1e-12)
  if(any(!is.finite(P.Z)) || sum(P.Z) <= 0) return(object)
  P.Z <- P.Z / sum(P.Z)
  expectation <- .lpa.expectation.reference(
    response, as.matrix(object$params$means), stabilized$covs, P.Z
  )
  if(!isTRUE(expectation$valid)) return(object)

  object$params$covs <- stabilized$covs
  object$params$P.Z <- object$P.Z <- P.Z
  object$P.Z.Xn <- expectation$posterior
  object$Z <- max.col(object$P.Z.Xn, ties.method = "first")
  object$Log.Lik <- expectation$Log.Lik
  if(!is.null(object$npar)){
    object$AIC <- -2 * object$Log.Lik + 2 * object$npar
    object$BIC <- -2 * object$Log.Lik + object$npar * log(nrow(response))
    object$best_BIC <- object$BIC
  }
  if(!is.null(object$Log.Lik.history) && length(object$Log.Lik.history) > 0L){
    if(stabilized$repaired){
      object$Log.Lik.history <- object$Log.Lik
    }else{
      object$Log.Lik.history[length(object$Log.Lik.history)] <- object$Log.Lik
    }
  }
  if(!is.null(object$Log.Lik.nrep) && length(object$Log.Lik.nrep) > 0L){
    if(stabilized$repaired){
      object$Log.Lik.nrep <- object$Log.Lik
    }else if(length(object$Log.Lik.nrep) == 1L){
      object$Log.Lik.nrep <- object$Log.Lik
    }else{
      object$Log.Lik.nrep[which.max(object$Log.Lik.nrep)] <- object$Log.Lik
    }
  }
  object$covariance.repaired <- isTRUE(object$covariance.repaired) ||
    stabilized$repaired
  if(stabilized$repaired || is.null(object$covariance.repair)){
    object$covariance.repair <- stabilized$diagnostics
  }
  object
}
