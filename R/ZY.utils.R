.ZY.CEP.inverse <- function(CEP, tolerance = sqrt(.Machine$double.eps)) {
  CEP <- as.matrix(CEP)
  singular.values <- svd(CEP, nu = 0L, nv = 0L)$d
  largest <- max(singular.values)
  threshold <- tolerance * largest
  rank <- if(is.finite(largest) && largest > 0){
    sum(singular.values > threshold)
  }else{
    0L
  }
  reciprocal.condition <- if(rank == nrow(CEP)){
    min(singular.values) / largest
  }else{
    0
  }
  generalized <- rank < nrow(CEP) ||
    !is.finite(reciprocal.condition) || reciprocal.condition <= tolerance
  inverse <- if(generalized){
    MASS::ginv(CEP, tol = tolerance)
  }else{
    solve(CEP)
  }
  list(
    inverse = inverse,
    method = if(generalized) "MASS::ginv" else "solve",
    rank = rank,
    condition = if(reciprocal.condition > 0){
      1 / reciprocal.condition
    }else{
      Inf
    },
    reciprocal.condition = reciprocal.condition,
    threshold = threshold
  )
}

.ZY.stabilize.weight.mass <- function(weight, fallback) {
  weight <- as.matrix(weight)
  fallback <- as.matrix(fallback)
  if(!identical(dim(weight), dim(fallback))){
    stop("weight and fallback must have identical dimensions")
  }
  fallback.mass <- colSums(fallback)
  class.mass <- colSums(weight)
  minimum.mass <- max(1e-8, 1e-6 * sum(fallback.mass))
  regularization <- 0
  unstable <- !is.finite(class.mass) | class.mass <= minimum.mass
  if(any(unstable)){
    required <- (minimum.mass - class.mass[unstable]) /
      (fallback.mass[unstable] - class.mass[unstable])
    required[!is.finite(required)] <- 1
    regularization <- min(1, max(0, required) + sqrt(.Machine$double.eps))
    weight <- (1 - regularization) * weight + regularization * fallback
    class.mass <- colSums(weight)
    if(any(!is.finite(class.mass)) || any(class.mass <= minimum.mass)){
      regularization <- 1
      weight <- fallback
      class.mass <- fallback.mass
    }
  }
  list(
    weight = weight,
    class.mass = class.mass,
    regularization = regularization,
    minimum.mass = minimum.mass
  )
}

.ZY.weights <- function(posterior, CEP = NULL, CEP.error = TRUE) {
  posterior <- as.matrix(posterior)
  modal <- max.col(posterior, ties.method = "first")
  if(is.null(CEP)){
    CEP <- if(CEP.error){
      get.CEP(list(posterior), CEP.time.cross = FALSE)[[1L]]
    }else{
      diag(ncol(posterior))
    }
  }
  inverse <- .ZY.CEP.inverse(CEP)
  case.weights <- inverse$inverse[modal, , drop = FALSE]
  row.mass <- rowSums(case.weights)
  stable.row <- is.finite(row.mass) & abs(row.mass) > sqrt(.Machine$double.eps)
  if(any(stable.row)){
    case.weights[stable.row, ] <- case.weights[stable.row, , drop = FALSE] /
      row.mass[stable.row]
  }
  if(any(!stable.row)){
    case.weights[!stable.row, ] <- posterior[!stable.row, , drop = FALSE]
  }

  stabilized <- .ZY.stabilize.weight.mass(case.weights, posterior)
  case.weights <- stabilized$weight
  colnames(case.weights) <- colnames(posterior)
  list(
    CEP = CEP,
    case.weights = case.weights,
    inverse.method = inverse$method,
    rank = inverse$rank,
    condition = inverse$condition,
    reciprocal.condition = inverse$reciprocal.condition,
    inverse.threshold = inverse$threshold,
    weight.regularization = stabilized$regularization,
    reconstruction.error = max(abs(
      colSums(case.weights) - colSums(posterior)
    ))
  )
}

.ZY.project.simplex <- function(value) {
  value <- as.numeric(value)
  if(any(!is.finite(value))) return(rep(1 / length(value), length(value)))
  ordered <- sort(value, decreasing = TRUE)
  cumulative <- cumsum(ordered)
  active <- which(ordered - (cumulative - 1) / seq_along(ordered) > 0)
  if(!length(active)) return(rep(1 / length(value), length(value)))
  threshold <- (cumulative[max(active)] - 1) / max(active)
  pmax(value - threshold, 0)
}

.ZY.dependent.variables <- function(dependent.variables, N) {
  if(is.data.frame(dependent.variables)){
    res <- dependent.variables
  }else if(is.matrix(dependent.variables)){
    res <- as.data.frame(dependent.variables, stringsAsFactors = FALSE)
  }else{
    res <- data.frame(Dependent.Variable = dependent.variables)
  }
  if(nrow(res) != N){
    stop("dependent.variables must contain one row for every fitted observation")
  }
  if(ncol(res) < 1L){
    stop("dependent.variables must contain at least one observed dependent variable")
  }
  if(is.null(names(res)) || any(names(res) == "")){
    names(res) <- paste0("Dependent.Variable.", seq_len(ncol(res)))
  }
  res
}

.ZY.family <- function(family, dependent.variables) {
  if(length(family) != 1L && length(family) != ncol(dependent.variables)){
    stop("family must have length one or one value per dependent variable")
  }
  family <- rep(family, length.out = ncol(dependent.variables))
  family <- vapply(family, match.arg, character(1),
                   choices = c("gaussian", "categorical"))
  names(family) <- names(dependent.variables)
  family
}

.ZY.path.matrix <- function(matrices) {
  if(!is.list(matrices) || !length(matrices)){
    stop("matrices must be a non-empty list")
  }
  matrices <- lapply(matrices, as.matrix)
  N <- nrow(matrices[[1L]])
  L <- ncol(matrices[[1L]])
  if(L < 2L || any(vapply(matrices, nrow, integer(1)) != N) ||
     any(vapply(matrices, ncol, integer(1)) != L)){
    stop("All matrices must have the same dimensions and at least two columns")
  }
  value <- matrices[[1L]]
  if(length(matrices) > 1L){
    for(t in 2:length(matrices)){
      value <- do.call(cbind, lapply(seq_len(L), function(l){
        value * matrices[[t]][, l]
      }))
    }
  }
  paths <- make_latent_paths_cpp(L, length(matrices)) + 1L
  path.names <- .latent.path.names(paths)
  colnames(value) <- path.names
  colnames(paths) <- paste0("t", seq_len(length(matrices)))
  rownames(paths) <- path.names
  list(value = value, paths = paths)
}

.ZY.path.modal <- function(modal, L) {
  if(!is.list(modal) || !length(modal)){
    stop("modal must be a non-empty list")
  }
  N <- length(modal[[1L]])
  if(any(lengths(modal) != N)){
    stop("All modal-state vectors must have the same length")
  }
  index <- rep.int(1L, N)
  multiplier <- 1L
  for(t in seq_along(modal)){
    value <- as.integer(modal[[t]])
    if(anyNA(value) || any(value < 1L | value > L)){
      stop("Modal states must be integers between 1 and L")
    }
    index <- index + (value - 1L) * multiplier
    multiplier <- multiplier * L
  }
  index
}

.ZY.LTA.dependent.variables <- function(dependent.variables, family, dependent.variable.time, times, N) {
  dependent.variables.list <- vector("list", times)
  if(is.list(dependent.variables) && !is.data.frame(dependent.variables)){
    if(length(dependent.variables) != times){
      stop("dependent.variables must be a list with one element per time point")
    }
    selected <- if(is.null(dependent.variable.time)){
      which(!vapply(dependent.variables, is.null, logical(1)))
    }else{
      as.integer(dependent.variable.time)
    }
    if(!length(selected) || anyNA(selected) ||
       any(selected < 1L | selected > times) || anyDuplicated(selected)){
      stop("dependent.variable.time must contain distinct valid time-point indices")
    }
    for(t in selected){
      if(is.null(dependent.variables[[t]])){
        stop("Selected dependent-variable time points cannot be NULL")
      }
      dependent.variables.list[[t]] <- .ZY.dependent.variables(dependent.variables[[t]], N)
    }
  }else{
    selected <- if(is.null(dependent.variable.time)) times else as.integer(dependent.variable.time)
    if(length(selected) != 1L){
      stop("dependent.variable.time must select one time point when dependent.variables is not a list")
    }
    if(is.na(selected) || selected < 1L || selected > times){
      stop("dependent.variable.time must contain a valid time-point index")
    }
    dependent.variables.list[[selected]] <- .ZY.dependent.variables(dependent.variables, N)
  }
  if(!length(selected) || any(selected < 1L | selected > times)){
    stop("dependent.variable.time must contain valid time-point indices")
  }

  family.list <- vector("list", times)
  if(is.list(family)){
    if(length(family) != times){
      stop("family must be a list with one element per time point")
    }
    for(t in selected){
      family.list[[t]] <- .ZY.family(family[[t]], dependent.variables.list[[t]])
    }
  }else{
    for(t in selected){
      family.list[[t]] <- .ZY.family(family, dependent.variables.list[[t]])
    }
  }
  list(dependent.variables = dependent.variables.list, family = family.list, selected = selected)
}

.ZY.cluster.score <- function(score, cluster) {
  if(is.null(cluster)) return(score)
  rowsum(score, cluster, reorder = FALSE)
}

.ZY.omnibus <- function(estimate, vcov, contrast) {
  difference <- as.vector(contrast %*% estimate)
  covariance <- contrast %*% vcov %*% t(contrast)
  df <- if(all(is.finite(covariance))) qr(covariance)$rank else NA_integer_
  statistic <- if(!is.na(df) && df > 0L){
    as.numeric(t(difference) %*% MASS::ginv(covariance) %*% difference)
  }else{
    NA_real_
  }
  res <- list(
    statistic = c(Wald = statistic),
    parameter = c(df = df),
    p.value = c(p.value = if(!is.na(df) && df > 0L) stats::pchisq(
      statistic, df, lower.tail = FALSE
    ) else NA_real_),
    method = "Wald test of equality across latent groups",
    data.name = "observed dependent variable"
  )
  class(res) <- "htest"
  res
}

.ZY.fit.gaussian <- function(dependent.variable, weight, cluster = NULL,
                             method.regression = "Analytic",
                             method.SE = "Analytic") {
  complete <- !is.na(dependent.variable)
  if(!is.numeric(dependent.variable)){
    stop("gaussian dependent variables must be numeric")
  }
  dependent.variable.cur <- as.numeric(dependent.variable[complete])
  if(!length(dependent.variable.cur) || any(!is.finite(dependent.variable.cur))){
    stop("gaussian dependent variables must contain finite observed values")
  }
  weight.cur <- weight[complete, , drop = FALSE]
  cluster.cur <- if(is.null(cluster)) NULL else cluster[complete]
  stabilized <- .ZY.stabilize.weight.mass(
    weight.cur,
    matrix(1 / ncol(weight.cur), nrow(weight.cur), ncol(weight.cur))
  )
  weight.cur <- stabilized$weight
  class.mass <- stabilized$class.mass
  estimate <- colSums(weight.cur * dependent.variable.cur) / class.mass
  G <- length(estimate)
  index.mean <- seq_len(G)
  index.variance <- G + seq_len(G)
  residual <- sweep(
    matrix(dependent.variable.cur, nrow(weight.cur), ncol(weight.cur)),
    2, estimate, "-"
  )
  variance.fallback <- mean((dependent.variable.cur - mean(dependent.variable.cur))^2)
  if(!is.finite(variance.fallback) || variance.fallback <= 0) variance.fallback <- 1
  variance.floor <- max(.Machine$double.eps, variance.fallback * 1e-8)
  variance.raw <- colSums(weight.cur * residual^2) / class.mass
  variance <- pmax(variance.raw, variance.floor)
  iterations <- 0L
  estimating.equation <- function(value) {
    mean.cur <- value[index.mean]
    variance.cur <- value[index.variance]
    residual.cur <- sweep(
      matrix(dependent.variable.cur, nrow(weight.cur), ncol(weight.cur)),
      2, mean.cur, "-"
    )
    c(
      colSums(weight.cur * residual.cur),
      colSums(weight.cur * sweep(residual.cur^2, 2, variance.cur, "-"))
    )
  }
  if(method.regression == "Numeric"){
    optimization <- stats::optim(
      c(estimate, variance),
      function(value) sum((estimating.equation(value) /
                             rep(class.mass, 2L))^2),
      method = "L-BFGS-B",
      lower = c(rep(-Inf, G), rep(variance.floor, G))
    )
    estimate <- optimization$par[index.mean]
    variance <- optimization$par[index.variance]
    iterations <- unname(optimization$counts["function"])
  }
  residual <- sweep(
    matrix(dependent.variable.cur, nrow(weight.cur), ncol(weight.cur)), 2, estimate, "-"
  )
  score <- .ZY.cluster.score(cbind(
    weight.cur * residual,
    weight.cur * sweep(residual^2, 2, variance, "-")
  ), cluster.cur)
  bread <- if(method.SE == "Numeric"){
    -numDeriv::jacobian(estimating.equation, c(estimate, variance))
  }else{
    value <- matrix(0, 2L * G, 2L * G)
    value[index.mean, index.mean] <- diag(class.mass, nrow = G)
    value[index.variance, index.mean] <- diag(
      2 * colSums(weight.cur * residual), nrow = G
    )
    value[index.variance, index.variance] <- diag(class.mass, nrow = G)
    value
  }
  inverse.bread <- MASS::ginv(bread)
  vcov.parameters <- inverse.bread %*% crossprod(score) %*% t(inverse.bread)
  vcov <- vcov.parameters[index.mean, index.mean, drop = FALSE]
  variance.vcov <- vcov.parameters[index.variance, index.variance, drop = FALSE]
  mean.variance.vcov <- vcov.parameters[index.mean, index.variance, drop = FALSE]
  names(estimate) <- colnames(weight)
  names(variance) <- names(estimate)
  dimnames(vcov) <- list(names(estimate), names(estimate))
  dimnames(variance.vcov) <- list(names(estimate), names(estimate))
  dimnames(mean.variance.vcov) <- list(names(estimate), names(estimate))
  se <- sqrt(pmax(diag(vcov), 0))
  names(se) <- names(estimate)
  variance.se <- sqrt(pmax(diag(variance.vcov), 0))
  names(variance.se) <- names(estimate)
  contrast <- cbind(diag(ncol(weight) - 1L), -1)
  list(
    estimate = estimate,
    se = se,
    vcov = vcov,
    vcov.parameters = vcov.parameters,
    information = bread,
    omnibus = .ZY.omnibus(estimate, vcov, contrast),
    class.mass = class.mass,
    weight.regularization = stabilized$regularization,
    prior = class.mass / sum(class.mass),
    variance = variance,
    variance.se = variance.se,
    variance.vcov = variance.vcov,
    mean.variance.vcov = mean.variance.vcov,
    omnibus.variance = .ZY.omnibus(variance, variance.vcov, contrast),
    variance.regularization = max(abs(variance - variance.raw)),
    admissible = TRUE,
    class.shift = 0,
    observations = sum(complete),
    omitted = sum(!complete),
    levels = NULL,
    converged = TRUE,
    iterations = iterations
  )
}

.ZY.fit.categorical <- function(dependent.variable, weight, cluster = NULL,
                                levels.dependent.variable = NULL,
                                method.regression = "Analytic",
                                method.SE = "Analytic") {
  complete <- !is.na(dependent.variable)
  if(is.null(levels.dependent.variable)){
    levels.dependent.variable <- levels(factor(dependent.variable[complete]))
  }
  dependent.variable.factor <- factor(dependent.variable[complete], levels = levels.dependent.variable)
  K <- nlevels(dependent.variable.factor)
  if(K < 2L){
    stop("categorical dependent variables must contain at least two observed categories")
  }
  weight.cur <- weight[complete, , drop = FALSE]
  cluster.cur <- if(is.null(cluster)) NULL else cluster[complete]
  stabilized <- .ZY.stabilize.weight.mass(
    weight.cur,
    matrix(1 / ncol(weight.cur), nrow(weight.cur), ncol(weight.cur))
  )
  weight.cur <- stabilized$weight
  class.mass <- stabilized$class.mass
  indicator <- diag(K)[as.integer(dependent.variable.factor), , drop = FALSE]
  estimate.raw <- sweep(t(weight.cur) %*% indicator, 1, class.mass, "/")
  dimnames(estimate.raw) <- list(colnames(weight), levels.dependent.variable)
  estimate <- t(apply(estimate.raw, 1L, .ZY.project.simplex))
  dimnames(estimate) <- dimnames(estimate.raw)
  regularized <- max(abs(estimate - estimate.raw)) > sqrt(.Machine$double.eps)
  iterations <- 0L

  estimating.equation <- function(value){
    value <- matrix(value, nrow = ncol(weight.cur), byrow = TRUE)
    unlist(lapply(seq_len(ncol(weight.cur)), function(l){
      colSums(weight.cur[, l] * sweep(indicator, 2, value[l, ], "-"))
    }), use.names = FALSE)
  }
  if(method.regression == "Numeric"){
    optimization <- stats::optim(
      as.vector(t(estimate)),
      function(value) sum((estimating.equation(value) /
                             rep(class.mass, each = K))^2),
      method = "BFGS"
    )
    estimate.numeric <- matrix(
      optimization$par, nrow = ncol(weight.cur), byrow = TRUE,
      dimnames = dimnames(estimate)
    )
    estimate <- t(apply(estimate.numeric, 1L, .ZY.project.simplex))
    dimnames(estimate) <- dimnames(estimate.raw)
    regularized <- regularized ||
      max(abs(estimate - estimate.numeric)) > sqrt(.Machine$double.eps)
    iterations <- unname(optimization$counts["function"])
  }

  score <- matrix(0, nrow(weight.cur), ncol(weight.cur) * K)
  for(l in seq_len(ncol(weight.cur))){
    index <- (l - 1L) * K + seq_len(K)
    score[, index] <- weight.cur[, l] * sweep(
      indicator, 2, estimate[l, ], "-"
    )
  }
  score <- .ZY.cluster.score(score, cluster.cur)
  bread <- if(method.SE == "Numeric"){
    -numDeriv::jacobian(estimating.equation, as.vector(t(estimate)))
  }else{
    diag(rep(class.mass, each = K), nrow = ncol(score))
  }
  inverse.bread <- MASS::ginv(bread)
  vcov <- inverse.bread %*% crossprod(score) %*% inverse.bread
  se <- matrix(
    sqrt(pmax(diag(vcov), 0)), nrow = ncol(weight.cur), byrow = TRUE,
    dimnames = dimnames(estimate)
  )
  class.contrast <- cbind(diag(ncol(weight.cur) - 1L), -1)
  contrast <- kronecker(class.contrast, diag(K))
  list(
    estimate = estimate,
    se = se,
    vcov = vcov,
    information = bread,
    omnibus = .ZY.omnibus(as.vector(t(estimate)), vcov, contrast),
    class.mass = class.mass,
    weight.regularization = stabilized$regularization,
    prior = class.mass / sum(class.mass),
    variance = NULL,
    admissible = TRUE,
    regularized = regularized,
    class.shift = 0,
    observations = sum(complete),
    omitted = sum(!complete),
    levels = levels.dependent.variable,
    converged = TRUE,
    iterations = iterations
  )
}

.ZY.fit.dependent.variables <- function(dependent.variables, family, weight, cluster = NULL,
                             levels.dependent.variables = NULL,
                             method.regression = "Analytic",
                             method.SE = "Analytic") {
  res <- vector("list", ncol(dependent.variables))
  names(res) <- names(dependent.variables)
  for(j in seq_len(ncol(dependent.variables))){
    if(family[j] == "gaussian"){
      res[[j]] <- .ZY.fit.gaussian(
        dependent.variables[[j]], weight, cluster, method.regression, method.SE
      )
    }else{
      levels.cur <- if(is.null(levels.dependent.variables)) NULL else levels.dependent.variables[[j]]
      res[[j]] <- .ZY.fit.categorical(
        dependent.variables[[j]], weight, cluster, levels.cur,
        method.regression, method.SE
      )
    }
    res[[j]]$family <- family[j]
  }
  res
}

.ZY.vector <- function(models) {
  unlist(lapply(models, function(x) {
    if(identical(unname(x$family), "gaussian")){
      c(Mean = x$estimate, Variance = x$variance)
    }else{
      as.vector(t(x$estimate))
    }
  }), use.names = TRUE)
}

.ZY.flatten.models <- function(models) {
  if(!length(models)) return(list())
  is.model <- vapply(models, function(x){
    is.list(x) && !is.null(x$family)
  }, logical(1))
  if(all(is.model)) return(models)
  unlist(lapply(models, .ZY.flatten.models), recursive = FALSE)
}

.ZY.progress.models <- function(models, method.3step, method.SE,
                                vis = TRUE) {
  if(!vis) return(invisible(NULL))
  models <- .ZY.flatten.models(models)
  n.models <- length(models)
  if(!n.models) return(invisible(NULL))
  family <- vapply(models, function(x) as.character(x$family)[1L], character(1))
  family.count <- table(factor(family, levels = c("gaussian", "categorical")))
  family.output <- paste(
    paste0(
      c("Gaussian", "categorical")[family.count > 0L], " = ",
      as.integer(family.count[family.count > 0L])
    ),
    collapse = ", "
  )
  converged <- sum(vapply(models, function(x) isTRUE(x$converged), logical(1)))
  iterations <- sum(vapply(models, function(x){
    if(is.null(x$iterations) || !is.finite(x$iterations)) 0 else x$iterations
  }, numeric(1)))
  output <- sprintf(
    "  %s dependent-variable models = %d (%s) | Converged = %d/%d | Iterations = %d",
    method.3step, n.models, family.output, converged, n.models, iterations
  )
  cat(output, "\n", sep = "")
  if(method.SE != "Bootstrap"){
    cat(sprintf("  %s Standard Errors were computed.\n", method.SE))
  }
  invisible(NULL)
}

.ZY.progress.bootstrap.complete <- function(SE.diagnostics, vis = TRUE) {
  if(!vis) return(invisible(NULL))
  if(!is.null(SE.diagnostics$groups)){
    diagnostics <- SE.diagnostics$groups
    successful <- vapply(diagnostics, `[[`, numeric(1), "successful")
    attempted <- vapply(diagnostics, `[[`, numeric(1), "attempted")
    cat(sprintf(
      "  Bootstrap Standard Errors were computed from %d--%d/%d successful replications across dependent-variable groups.\n",
      min(successful), max(successful), max(attempted)
    ))
  }else{
    cat(sprintf(
      "  Bootstrap Standard Errors were computed from %d/%d successful replications.\n",
      SE.diagnostics$successful, SE.diagnostics$attempted
    ))
  }
  invisible(NULL)
}

.ZY.bootstrap <- function(models, estimates.bootstrap) {
  successful <- complete.cases(estimates.bootstrap)
  if(sum(successful) < 2L){
    vcov <- matrix(NA_real_, ncol(estimates.bootstrap), ncol(estimates.bootstrap))
  }else{
    vcov <- stats::cov(estimates.bootstrap[successful, , drop = FALSE])
  }
  index <- 0L
  for(j in seq_along(models)){
    if(identical(unname(models[[j]]$family), "gaussian")){
      G <- length(models[[j]]$estimate)
      index.mean <- index + seq_len(G)
      index.variance <- index + G + seq_len(G)
      index.cur <- c(index.mean, index.variance)
      vcov.cur <- vcov[index.cur, index.cur, drop = FALSE]
      models[[j]]$reported.vcov <- vcov.cur
      models[[j]]$vcov <- vcov[index.mean, index.mean, drop = FALSE]
      models[[j]]$variance.vcov <- vcov[index.variance, index.variance, drop = FALSE]
      models[[j]]$mean.variance.vcov <- vcov[index.mean, index.variance, drop = FALSE]
      dimnames(models[[j]]$vcov) <- list(
        names(models[[j]]$estimate), names(models[[j]]$estimate)
      )
      dimnames(models[[j]]$variance.vcov) <- list(
        names(models[[j]]$variance), names(models[[j]]$variance)
      )
      dimnames(models[[j]]$mean.variance.vcov) <- list(
        names(models[[j]]$estimate), names(models[[j]]$variance)
      )
      models[[j]]$se <- setNames(
        sqrt(pmax(diag(models[[j]]$vcov), 0)), names(models[[j]]$estimate)
      )
      models[[j]]$variance.se <- setNames(
        sqrt(pmax(diag(models[[j]]$variance.vcov), 0)),
        names(models[[j]]$variance)
      )
      contrast <- cbind(diag(G - 1L), -1)
      models[[j]]$omnibus <- .ZY.omnibus(
        models[[j]]$estimate, models[[j]]$vcov, contrast
      )
      models[[j]]$omnibus.variance <- .ZY.omnibus(
        models[[j]]$variance, models[[j]]$variance.vcov, contrast
      )
      index <- index + 2L * G
    }else{
      npar <- length(models[[j]]$estimate)
      index.cur <- index + seq_len(npar)
      vcov.cur <- vcov[index.cur, index.cur, drop = FALSE]
      se.cur <- sqrt(pmax(diag(vcov.cur), 0))
      models[[j]]$se <- matrix(
        se.cur, nrow = nrow(models[[j]]$estimate), byrow = TRUE,
        dimnames = dimnames(models[[j]]$estimate)
      )
      class.contrast <- cbind(diag(nrow(models[[j]]$estimate) - 1L), -1)
      contrast <- kronecker(class.contrast, diag(ncol(models[[j]]$estimate)))
      estimate.cur <- as.vector(t(models[[j]]$estimate))
      models[[j]]$vcov <- vcov.cur
      models[[j]]$omnibus <- .ZY.omnibus(estimate.cur, vcov.cur, contrast)
      index <- index + npar
    }
  }
  list(
    models = models,
    vcov = vcov,
    diagnostics = list(
      method = "Bootstrap",
      successful = sum(successful),
      attempted = nrow(estimates.bootstrap)
    )
  )
}
