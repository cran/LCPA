#' Compute Standard Errors
#'
#' Computes standard errors (SEs) for parameters estimated by \code{\link[LCPA]{LCA}()} or \code{\link[LCPA]{LPA}()}.
#' Three methods are available:
#' \itemize{
#'   \item \code{"Bootstrap"}: Nonparametric bootstrap with label-switching correction.
#'   \item \code{"Obs"}: Numerical observed information from adaptive central differences
#'     of the analytic observed-data score.
#'   \item \code{"Louis"}: Analytic observed information based on Louis' identity.
#' }
#'
#' @param object An object of class \code{"LCA"} or \code{"LPA"} returned by \code{\link[LCPA]{LCA}()} or \code{\link[LCPA]{LPA}()}.
#' @param method Character specifying \code{"Bootstrap"}, \code{"Obs"}, or \code{"Louis"}.
#' @param nrep.bootstrap Integer. Number of successful bootstrap replicates when
#'   \code{method = "Bootstrap"} (default: 100).
#' @param vis Logical. If \code{TRUE}, displays bootstrap progress (default: \code{TRUE}).
#'
#' @return A list of class \code{"SE"} containing:
#'   \describe{
#'     \item{\code{se}}{SEs with the same parameter structure and dimnames as the fitted model.}
#'     \item{\code{vcov}}{Covariance matrix on the unconstrained estimation scale for
#'       \code{"Obs"} and \code{"Louis"}; \code{NULL} for bootstrap.}
#'     \item{\code{hessian}}{Observed information matrix on the unconstrained estimation scale for
#'       \code{"Obs"} and \code{"Louis"}; \code{NULL} for bootstrap.}
#'     \item{\code{diagnostics}}{Method-specific diagnostics.}
#'     \item{\code{call}}{Function call that generated the object.}
#'     \item{\code{arguments}}{Input arguments.}
#'   }
#'
#' @details
#' Class proportions and LCA conditional probabilities are represented by additive log-ratios
#' with the final category as reference. Standard errors for every probability, including the
#' reference probability, are obtained using the full softmax Jacobian.
#'
#' For LPA, the observed-information parameterization includes only covariance parameters that
#' are free under the fitted \code{constraint}; variances are represented on the log scale.
#' This avoids the nonidentified Hessian directions produced by differentiating all class-specific
#' covariance elements under equality or zero constraints.
#'
#' The analytic \code{"Louis"} method evaluates the observed information from the conditional
#' complete-data score and Hessian. For both LCA and LPA it accounts for posterior class uncertainty
#' and therefore is not the naive complete-data information matrix. The numerical \code{"Obs"}
#' method provides an independent check by differentiating the analytic observed-data score.
#'
#' @references
#' Louis, T. A. (1982). Finding the observed information matrix when using the EM algorithm.
#' *Journal of the Royal Statistical Society: Series B (Methodological),
#' 44*(2), 226--233.
#' \doi{10.1111/j.2517-6161.1982.tb01203.x}
#'
#' McLachlan, G. J., & Peel, D. (2000). *Finite mixture models*. John Wiley & Sons.
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' lca.data <- sim.LCA(N = 300, I = 4, L = 2, poly.value = 3)
#' lca.fit <- LCA(lca.data$response, L = 2, nrep = 3, vis = FALSE)
#' se.louis <- get.SE(lca.fit, method = "Louis")
#'
#' lpa.data <- sim.LPA(N = 300, I = 3, L = 2, constraint = "V0")
#' lpa.fit <- LPA(lpa.data$response, L = 2, constraint = "V0", nrep = 3, vis = FALSE)
#' se.louis.lpa <- get.SE(lpa.fit, method = "Louis")
#' }
#'
#' @importFrom clue solve_LSAP
#' @importFrom stats sd
#' @export
get.SE <- function(object, method = "Bootstrap", nrep.bootstrap = 100, vis = TRUE) {
  if(!inherits(object, c("LCA", "LPA"))){
    stop("object must be of class 'LCA' or 'LPA'")
  }

  call <- match.call()
  method <- match.arg(method, c("Bootstrap", "Obs", "Louis"))
  type <- class(object)[1]

  if(method == "Bootstrap"){
    result <- bootstrap.SE(object, nrep.bootstrap = nrep.bootstrap, vis = vis)
  }else if(method == "Louis"){
    result <- if(type == "LCA") louis.SE.LCA(object) else louis.SE.LPA(object)
  }else if(type == "LCA"){
    result <- observed.SE.LCA(object)
  }else{
    result <- observed.SE.LPA(object)
  }

  res <- list(
    se = result$se,
    vcov = result$vcov,
    hessian = result$hessian,
    diagnostics = result$diagnostics,
    call = call,
    arguments = list(object = object, method = method, nrep.bootstrap = nrep.bootstrap)
  )
  class(res) <- "SE"
  res
}

softmax.reference <- function(eta){
  if(length(eta) == 0){
    return(1)
  }
  eta <- c(eta, 0)
  eta <- eta - max(eta)
  prob <- exp(eta)
  prob / sum(prob)
}

softmax.jacobian <- function(prob){
  K <- length(prob)
  if(K == 1){
    return(matrix(numeric(0), nrow = 1, ncol = 0))
  }
  diag(prob)[, seq_len(K - 1), drop = FALSE] - outer(prob, prob[seq_len(K - 1)])
}

invert.information <- function(information){
  information <- (information + t(information)) / 2
  eigenvalues.initial <- eigen(
    information, symmetric = TRUE, only.values = TRUE
  )$values
  repaired <- is.null(tryCatch(
    chol(information), error = function(e){ NULL }
  ))
  information.used <- if(repaired){
    as.matrix(Matrix::nearPD(
      information, base.matrix = TRUE,
      eig.tol = .Machine$double.eps^(2 / 3),
      conv.tol = sqrt(.Machine$double.eps),
      posd.tol = sqrt(.Machine$double.eps)
    )$mat)
  }else{
    information
  }
  decomposition <- chol(information.used)
  vcov <- chol2inv(decomposition)
  eigenvalues.final <- eigen(
    information.used, symmetric = TRUE, only.values = TRUE
  )$values
  scale <- max(1, max(abs(eigenvalues.initial)))
  adjusted <- abs(eigenvalues.final - eigenvalues.initial) >
    sqrt(.Machine$double.eps) * scale

  list(
    vcov = vcov,
    information = information,
    information.used = information.used,
    condition.initial = if(any(abs(eigenvalues.initial) < .Machine$double.eps)) Inf else
      max(abs(eigenvalues.initial)) / min(abs(eigenvalues.initial)),
    condition.final = max(eigenvalues.final) / min(eigenvalues.final),
    adjusted.fraction = mean(adjusted),
    threshold = if(repaired) min(eigenvalues.final) else NA_real_
  )
}

numeric.information.from.score <- function(theta, score.function,
                                           rel.step = 1e-4, max.halvings = 12L){
  q <- length(theta)
  score.estimate <- tryCatch(score.function(theta), error = function(e) NULL)
  if(is.null(score.estimate) || length(score.estimate) != q ||
     any(!is.finite(score.estimate))){
    stop("analytic score is not finite at the fitted parameter estimates")
  }

  score.jacobian <- matrix(NA_real_, q, q)
  step.sizes <- numeric(q)
  for(j in seq_len(q)){
    h <- rel.step * (abs(theta[j]) + 1)
    success <- FALSE
    for(attempt in 0:max.halvings){
      theta.plus <- theta.minus <- theta
      theta.plus[j] <- theta.plus[j] + h
      theta.minus[j] <- theta.minus[j] - h
      score.plus <- tryCatch(score.function(theta.plus), error = function(e) NULL)
      score.minus <- tryCatch(score.function(theta.minus), error = function(e) NULL)

      theta.plus[j] <- theta[j] + h / 2
      theta.minus[j] <- theta[j] - h / 2
      score.plus.half <- tryCatch(score.function(theta.plus), error = function(e) NULL)
      score.minus.half <- tryCatch(score.function(theta.minus), error = function(e) NULL)
      scores <- c(score.plus, score.minus, score.plus.half, score.minus.half)
      if(length(scores) == 4L * q && all(is.finite(scores))){
        derivative.full <- (score.plus - score.minus) / (2 * h)
        derivative.half <- (score.plus.half - score.minus.half) / h
        score.jacobian[, j] <- (4 * derivative.half - derivative.full) / 3
        step.sizes[j] <- h / 2
        success <- TRUE
        break
      }
      h <- h / 2
    }
    if(!success){
      stop("could not obtain a finite central difference for parameter ", j)
    }
  }

  information.raw <- -score.jacobian
  asymmetry <- max(abs(information.raw - t(information.raw))) /
    max(1, max(abs(information.raw)))
  list(
    information = (information.raw + t(information.raw)) / 2,
    score = score.estimate,
    step.sizes = step.sizes,
    relative.asymmetry = asymmetry
  )
}

lca.parameterization <- function(object){
  params <- object$params
  if (is.null(params$category.levels)) {
    stop("LCA parameters must include category.levels")
  }
  response <- .LCA.encode.response(
    object$arguments$response, params$category.levels
  )
  poly.value <- lengths(params$category.levels)
  poly.max <- max(poly.value)
  L <- length(params$P.Z)
  I <- ncol(response)

  theta <- if(L > 1) log(params$P.Z[seq_len(L - 1)] / params$P.Z[L]) else numeric(0)
  blocks <- vector("list", I * L)
  block.cur <- 0L
  position <- length(theta) + 1L
  for(i in seq_len(I)){
    K <- poly.value[i]
    for(l in seq_len(L)){
      block.cur <- block.cur + 1L
      index <- if(K > 1) position:(position + K - 2L) else integer(0)
      blocks[[block.cur]] <- list(i = i, l = l, K = K, index = index)
      if(K > 1){
        prob <- params$par[l, i, seq_len(K)]
        theta <- c(theta, log(prob[seq_len(K - 1)] / prob[K]))
        position <- position + K - 1L
      }
    }
  }

  list(theta = theta, blocks = blocks, poly.value = poly.value,
       poly.max = poly.max, L = L, I = I, response = response)
}

unpack.LCA <- function(theta, specification){
  L <- specification$L
  I <- specification$I
  P.Z <- softmax.reference(if(L > 1) theta[seq_len(L - 1)] else numeric(0))
  par <- array(NA_real_, dim = c(L, I, specification$poly.max))
  for(block in specification$blocks){
    prob <- softmax.reference(theta[block$index])
    par[block$l, block$i, seq_len(block$K)] <- prob
  }
  list(P.Z = P.Z, par = par)
}

lca.SE.from.vcov <- function(vcov, object, specification){
  params <- object$params
  L <- specification$L
  I <- specification$I
  se.P.Z <- numeric(L)
  if(L > 1){
    J <- softmax.jacobian(params$P.Z)
    index <- seq_len(L - 1)
    se.P.Z <- sqrt(pmax(diag(J %*% vcov[index, index, drop = FALSE] %*% t(J)), 0))
  }

  se.par <- array(NA_real_, dim = dim(params$par), dimnames = dimnames(params$par))
  for(block in specification$blocks){
    if(block$K == 1){
      se.par[block$l, block$i, 1] <- 0
    }else{
      prob <- params$par[block$l, block$i, seq_len(block$K)]
      J <- softmax.jacobian(prob)
      V <- vcov[block$index, block$index, drop = FALSE]
      se.par[block$l, block$i, seq_len(block$K)] <- sqrt(pmax(diag(J %*% V %*% t(J)), 0))
    }
  }
  names(se.P.Z) <- names(params$P.Z)
  list(par = se.par, P.Z = se.P.Z)
}

lca.score.components <- function(theta, specification,
                                 compute.information = FALSE){
  params <- unpack.LCA(theta, specification)
  result <- lca_score_information_cpp(
    specification$response, as.vector(params$par), params$P.Z,
    as.integer(specification$poly.value), compute.information
  )
  if(any(!is.finite(result$score.observation))){
    stop("non-finite LCA posterior probabilities")
  }
  result$params <- params
  result
}

lca.observed.score <- function(theta, specification){
  lca.score.components(theta, specification)$score
}

observed.SE.LCA <- function(object){
  specification <- lca.parameterization(object)
  numerical <- numeric.information.from.score(
    specification$theta,
    function(theta) lca.observed.score(theta, specification)
  )
  inverse <- invert.information(numerical$information)

  list(
    se = lca.SE.from.vcov(inverse$vcov, object, specification),
    vcov = inverse$vcov,
    hessian = inverse$information,
    diagnostics = list(
      method = "Obs",
      parameterization = "additive log-ratio; analytic-score central differences",
      score_max_abs = max(abs(numerical$score)),
      information_relative_asymmetry = numerical$relative.asymmetry,
      hessian_cond_number_initial = inverse$condition.initial,
      hessian_cond_number_final = inverse$condition.final,
      eigenvalue_adjustment_fraction = inverse$adjusted.fraction,
      eigenvalue_threshold = inverse$threshold,
      step_sizes = numerical$step.sizes
    )
  )
}

louis.SE.LCA <- function(object){
  specification <- lca.parameterization(object)
  components <- lca.score.components(
    specification$theta, specification, compute.information = TRUE
  )
  inverse <- invert.information(components$information)
  list(
    se = lca.SE.from.vcov(inverse$vcov, object, specification),
    vcov = inverse$vcov,
    hessian = inverse$information,
    diagnostics = list(
      method = "Louis",
      parameterization = "additive log-ratio",
      score_max_abs = max(abs(components$score)),
      hessian_cond_number_initial = inverse$condition.initial,
      hessian_cond_number_final = inverse$condition.final,
      eigenvalue_adjustment_fraction = inverse$adjusted.fraction,
      eigenvalue_threshold = inverse$threshold
    )
  )
}

lpa.parameterization <- function(object){
  params <- object$params
  response <- object$arguments$response
  constraint <- object$arguments$constraint
  L <- length(params$P.Z)
  I <- ncol(response)

  theta <- if(L > 1) log(params$P.Z[seq_len(L - 1)] / params$P.Z[L]) else numeric(0)
  mean.index <- length(theta) + seq_len(L * I)
  theta <- c(theta, as.vector(t(params$means)))
  covariance.specification <- list()

  is.zero <- function(a, b){
    is.character(constraint) && a != b && constraint %in% c("E0", "V0", "UE", "UV")
  }
  is.shared <- function(a, b){
    if(is.character(constraint)){
      if(constraint %in% c("EE", "UE")) return(TRUE)
      if(constraint == "E0") return(a == b)
      if(constraint == "VE") return(a != b)
      if(constraint == "EV") return(a == b)
      return(FALSE)
    }
    key <- sort(c(a, b))
    any(vapply(constraint, function(x){ identical(sort(as.integer(x)), key) }, logical(1)))
  }

  for(a in seq_len(I)){
    for(b in a:I){
      if(is.zero(a, b)) next
      classes <- if(is.shared(a, b)) list(seq_len(L)) else lapply(seq_len(L), function(l) l)
      for(class.index in classes){
        value <- mean(params$covs[a, b, class.index])
        theta <- c(theta, if(a == b) log(value) else value)
        covariance.specification[[length(covariance.specification) + 1L]] <- list(
          a = a, b = b, classes = class.index, diagonal = a == b,
          index = length(theta)
        )
      }
    }
  }

  list(theta = theta, mean.index = mean.index,
       covariance.specification = covariance.specification,
       L = L, I = I, response = response, constraint = constraint)
}

unpack.LPA <- function(theta, specification){
  L <- specification$L
  I <- specification$I
  P.Z <- softmax.reference(if(L > 1) theta[seq_len(L - 1)] else numeric(0))
  means <- matrix(theta[specification$mean.index], nrow = L, byrow = TRUE)
  covs <- array(0, dim = c(I, I, L))
  for(entry in specification$covariance.specification){
    value <- if(entry$diagonal) exp(theta[entry$index]) else theta[entry$index]
    covs[entry$a, entry$b, entry$classes] <- value
    covs[entry$b, entry$a, entry$classes] <- value
  }
  list(P.Z = P.Z, means = means, covs = covs)
}

lpa.score.components <- function(theta, specification){
  params <- unpack.LPA(theta, specification)
  response <- specification$response
  N <- nrow(response)
  L <- specification$L
  I <- specification$I
  q <- length(theta)
  log.joint <- matrix(NA_real_, N, L)
  class.info <- vector("list", L)
  complete.score <- vector("list", L)

  for(l in seq_len(L)){
    R <- tryCatch(chol(params$covs[, , l]), error = function(e) NULL)
    if(is.null(R)) return(NULL)
    A <- chol2inv(R)
    dev <- sweep(response, 2, params$means[l, ], "-")
    u <- dev %*% A
    log.joint[, l] <- log(params$P.Z[l]) -
      0.5 * (I * log(2 * pi) + 2 * sum(log(diag(R))) + rowSums(u * dev))

    score.l <- matrix(0, N, q)
    if(L > 1){
      mixing.index <- seq_len(L - 1)
      score.l[, mixing.index] <- matrix(
        as.numeric(mixing.index == l) - params$P.Z[mixing.index],
        N, length(mixing.index), byrow = TRUE
      )
    }
    mean.index <- specification$mean.index[(l - 1L) * I + seq_len(I)]
    score.l[, mean.index] <- u

    entry.ids <- which(vapply(
      specification$covariance.specification,
      function(entry) l %in% entry$classes,
      logical(1)
    ))
    D <- vector("list", length(entry.ids))
    for(j in seq_along(entry.ids)){
      entry <- specification$covariance.specification[[entry.ids[j]]]
      D[[j]] <- matrix(0, I, I)
      if(entry$diagonal){
        D[[j]][entry$a, entry$a] <- params$covs[entry$a, entry$a, l]
      }else{
        D[[j]][entry$a, entry$b] <- 1
        D[[j]][entry$b, entry$a] <- 1
      }
      score.l[, entry$index] <- 0.5 * (
        rowSums((u %*% D[[j]]) * u) - sum(diag(A %*% D[[j]]))
      )
    }
    complete.score[[l]] <- score.l
    class.info[[l]] <- list(
      A = A, u = u, mean.index = mean.index,
      entry.ids = entry.ids, D = D
    )
  }

  row.max <- apply(log.joint, 1, max)
  if(any(!is.finite(row.max))) return(NULL)
  posterior <- exp(log.joint - row.max)
  row.total <- rowSums(posterior)
  if(any(!is.finite(row.total)) || any(row.total <= 0)) return(NULL)
  posterior <- posterior / row.total

  score.observation <- matrix(0, N, q)
  for(l in seq_len(L)){
    score.observation <- score.observation +
      complete.score[[l]] * posterior[, l]
  }

  list(
    params = params,
    posterior = posterior,
    complete.score = complete.score,
    score.observation = score.observation,
    score = colSums(score.observation),
    class.info = class.info
  )
}

lpa.observed.score <- function(theta, specification){
  components <- lpa.score.components(theta, specification)
  if(is.null(components)){
    return(rep(NA_real_, length(theta)))
  }
  components$score
}

lpa.SE.from.vcov <- function(vcov, object, specification){
  params <- object$params
  L <- specification$L
  I <- specification$I
  se.P.Z <- numeric(L)
  if(L > 1){
    J <- softmax.jacobian(params$P.Z)
    index <- seq_len(L - 1)
    se.P.Z <- sqrt(pmax(diag(J %*% vcov[index, index, drop = FALSE] %*% t(J)), 0))
  }
  names(se.P.Z) <- names(params$P.Z)

  se.means <- matrix(sqrt(pmax(diag(vcov)[specification$mean.index], 0)),
                     nrow = L, byrow = TRUE, dimnames = dimnames(params$means))
  se.covs <- array(0, dim = dim(params$covs), dimnames = dimnames(params$covs))
  for(entry in specification$covariance.specification){
    derivative <- if(entry$diagonal) params$covs[entry$a, entry$b, entry$classes] else rep(1, length(entry$classes))
    value <- abs(derivative) * sqrt(max(vcov[entry$index, entry$index], 0))
    se.covs[entry$a, entry$b, entry$classes] <- value
    se.covs[entry$b, entry$a, entry$classes] <- value
  }
  list(means = se.means, covs = se.covs, P.Z = se.P.Z)
}

observed.SE.LPA <- function(object){
  specification <- lpa.parameterization(object)
  numerical <- numeric.information.from.score(
    specification$theta,
    function(theta) lpa.observed.score(theta, specification)
  )
  inverse <- invert.information(numerical$information)

  list(
    se = lpa.SE.from.vcov(inverse$vcov, object, specification),
    vcov = inverse$vcov,
    hessian = inverse$information,
    diagnostics = list(
      method = "Obs",
      parameterization = paste(
        "constraint-aware log-variance/covariance;",
        "analytic-score central differences"
      ),
      free_parameters = length(specification$theta),
      expected_free_parameters = object$npar,
      score_max_abs = max(abs(numerical$score)),
      information_relative_asymmetry = numerical$relative.asymmetry,
      hessian_cond_number_initial = inverse$condition.initial,
      hessian_cond_number_final = inverse$condition.final,
      eigenvalue_adjustment_fraction = inverse$adjusted.fraction,
      eigenvalue_threshold = inverse$threshold,
      step_sizes = numerical$step.sizes
    )
  )
}

louis.SE.LPA <- function(object){
  specification <- lpa.parameterization(object)
  components <- lpa.score.components(specification$theta, specification)
  if(is.null(components)){
    stop("fitted LPA covariance matrices must be positive definite")
  }

  N <- nrow(specification$response)
  L <- specification$L
  I <- specification$I
  q <- length(specification$theta)
  expected.hessian <- matrix(0, q, q)
  if(L > 1){
    mixing.index <- seq_len(L - 1)
    prob <- components$params$P.Z[mixing.index]
    expected.hessian[mixing.index, mixing.index] <- -N *
      (diag(prob, nrow = length(prob)) - outer(prob, prob))
  }

  for(l in seq_len(L)){
    weight <- components$posterior[, l]
    nk <- sum(weight)
    info <- components$class.info[[l]]
    A <- info$A
    u <- info$u
    mean.index <- info$mean.index
    expected.hessian[mean.index, mean.index] <-
      expected.hessian[mean.index, mean.index] - nk * A

    weighted.u <- colSums(u * weight)
    for(j in seq_along(info$entry.ids)){
      entry.j <- specification$covariance.specification[[info$entry.ids[j]]]
      index.j <- entry.j$index
      mean.cov <- -A %*% info$D[[j]] %*% weighted.u
      expected.hessian[mean.index, index.j] <-
        expected.hessian[mean.index, index.j] + mean.cov
      expected.hessian[index.j, mean.index] <-
        expected.hessian[index.j, mean.index] + as.vector(mean.cov)

      for(k in j:length(info$entry.ids)){
        entry.k <- specification$covariance.specification[[info$entry.ids[k]]]
        index.k <- entry.k$index
        D.j <- info$D[[j]]
        D.k <- info$D[[k]]
        value.n <- -2 * rowSums((u %*% (D.j %*% A %*% D.k)) * u) +
          sum(diag(A %*% D.k %*% A %*% D.j))
        if(j == k && entry.j$diagonal){
          value.n <- value.n + rowSums((u %*% D.j) * u) -
            sum(diag(A %*% D.j))
        }
        value <- 0.5 * sum(weight * value.n)
        expected.hessian[index.j, index.k] <-
          expected.hessian[index.j, index.k] + value
        if(index.k != index.j){
          expected.hessian[index.k, index.j] <-
            expected.hessian[index.k, index.j] + value
        }
      }
    }
  }

  second.moment <- matrix(0, q, q)
  for(l in seq_len(L)){
    weighted.score <- components$complete.score[[l]] *
      sqrt(components$posterior[, l])
    second.moment <- second.moment + crossprod(weighted.score)
  }
  missing.information <- second.moment -
    crossprod(components$score.observation)
  observed.information <- -expected.hessian - missing.information
  inverse <- invert.information(observed.information)

  list(
    se = lpa.SE.from.vcov(inverse$vcov, object, specification),
    vcov = inverse$vcov,
    hessian = inverse$information,
    diagnostics = list(
      method = "Louis",
      parameterization = "constraint-aware log-variance/covariance",
      free_parameters = length(specification$theta),
      expected_free_parameters = object$npar,
      score_max_abs = max(abs(components$score)),
      hessian_cond_number_initial = inverse$condition.initial,
      hessian_cond_number_final = inverse$condition.final,
      eigenvalue_adjustment_fraction = inverse$adjusted.fraction,
      eigenvalue_threshold = inverse$threshold
    )
  )
}

bootstrap.SE <- function(object, nrep.bootstrap, vis){
  nrep.bootstrap <- as.integer(nrep.bootstrap)[1]
  if(is.na(nrep.bootstrap) || nrep.bootstrap < 2L){
    stop("nrep.bootstrap must be an integer greater than or equal to 2")
  }

  type <- class(object)[1]
  response <- object$arguments$response
  params <- object$params
  N <- nrow(response)
  I <- ncol(response)
  L <- length(params$P.Z)
  maxattempts <- max(5L * nrep.bootstrap, nrep.bootstrap + 10L)

  if(type == "LPA"){
    P.Z.Bootstrap <- matrix(NA_real_, L, nrep.bootstrap)
    means.Bootstrap <- array(NA_real_, c(L, I, nrep.bootstrap))
    covs.Bootstrap <- array(NA_real_, c(I, I, L, nrep.bootstrap))
  }else{
    if (is.null(params$category.levels)) {
      stop("LCA parameters must include category.levels")
    }
    poly.value <- lengths(params$category.levels)
    poly.max <- max(poly.value)
    P.Z.Bootstrap <- matrix(NA_real_, L, nrep.bootstrap)
    par.Bootstrap <- array(NA_real_, c(L, I, poly.max, nrep.bootstrap))
    par.reference <- do.call(rbind, lapply(seq_len(I), function(i){
      t(matrix(params$par[, i, seq_len(poly.value[i]), drop = FALSE], nrow = L))
    }))
  }

  completed <- 0L
  attempts <- 0L
  while(completed < nrep.bootstrap && attempts < maxattempts){
    attempts <- attempts + 1L
    if(vis){
      cat("\rRunning bootstrap (", paste0(completed + 1L, "/", nrep.bootstrap),
          " successful replicates; attempt ", attempts, ") ...", sep = "")
    }
    indices <- sample.int(N, size = N, replace = TRUE)
    response.cur <- response[indices, , drop = FALSE]

    if(type == "LCA"){
      categories.present <- vapply(seq_len(I), function(i){
        length(unique(response.cur[, i])) == poly.value[i]
      }, logical(1))
      if(!all(categories.present)) next
    }

    object.cur <- tryCatch(
      update(object, response = response.cur, par.ini = "random", vis = FALSE),
      error = function(e) NULL
    )
    if(is.null(object.cur)) next

    completed <- completed + 1L
    params.cur <- object.cur$params
    if(type == "LPA"){
      dist.mat <- distance.matrix(t(params$means), t(params.cur$means))
      assignment <- as.numeric(clue::solve_LSAP(dist.mat))
      means.Bootstrap[, , completed] <- params.cur$means[assignment, , drop = FALSE]
      covs.Bootstrap[, , , completed] <- params.cur$covs[, , assignment, drop = FALSE]
      P.Z.Bootstrap[, completed] <- params.cur$P.Z[assignment]
    }else{
      par.current <- do.call(rbind, lapply(seq_len(I), function(i){
        t(matrix(params.cur$par[, i, seq_len(poly.value[i]), drop = FALSE], nrow = L))
      }))
      dist.mat <- distance.matrix(par.reference, par.current)
      assignment <- as.numeric(clue::solve_LSAP(dist.mat))
      par.Bootstrap[, , , completed] <- params.cur$par[assignment, , , drop = FALSE]
      P.Z.Bootstrap[, completed] <- params.cur$P.Z[assignment]
    }
  }
  if(vis) cat("\n")
  if(completed < 2L){
    stop("fewer than two bootstrap fits completed successfully")
  }
  if(completed < nrep.bootstrap){
    warning("Only ", completed, " of ", nrep.bootstrap,
            " bootstrap fits completed after ", attempts, " attempts.")
  }

  keep <- seq_len(completed)
  se.P.Z <- apply(P.Z.Bootstrap[, keep, drop = FALSE], 1, sd)
  names(se.P.Z) <- names(params$P.Z)
  if(type == "LPA"){
    se.means <- apply(means.Bootstrap[, , keep, drop = FALSE], c(1, 2), sd)
    se.covs <- apply(covs.Bootstrap[, , , keep, drop = FALSE], c(1, 2, 3), sd)
    dimnames(se.means) <- dimnames(params$means)
    dimnames(se.covs) <- dimnames(params$covs)
    se <- list(means = se.means, covs = se.covs, P.Z = se.P.Z)
  }else{
    se.par <- apply(par.Bootstrap[, , , keep, drop = FALSE], c(1, 2, 3), sd, na.rm = TRUE)
    dimnames(se.par) <- dimnames(params$par)
    se <- list(par = se.par, P.Z = se.P.Z)
  }

  list(
    se = se,
    vcov = NULL,
    hessian = NULL,
    diagnostics = list(
      method = "Bootstrap",
      nrep.bootstrap.requested = nrep.bootstrap,
      nrep.bootstrap.completed = completed,
      attempts = attempts,
      maxattempts = maxattempts
    )
  )
}
