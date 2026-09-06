.ML.ZY.classification <- function(modal, CEP) {
  t(CEP[, modal, drop = FALSE])
}

.ML.ZY.softmax <- function(eta) {
  eta <- c(eta, 0)
  eta <- eta - max(eta)
  probability <- exp(eta)
  probability / sum(probability)
}

.ML.ZY.gaussian.inference <- function(mean, sd, vcov.parameters,
                                       index.mean, index.log.sd,
                                       group.names) {
  estimate <- setNames(mean, group.names)
  variance <- setNames(sd^2, group.names)
  vcov <- vcov.parameters[index.mean, index.mean, drop = FALSE]
  log.sd.vcov <- vcov.parameters[index.log.sd, index.log.sd, drop = FALSE]
  derivative.variance <- diag(2 * variance, nrow = length(variance))
  variance.vcov <- derivative.variance %*% log.sd.vcov %*%
    derivative.variance
  mean.variance.vcov <- vcov.parameters[
    index.mean, index.log.sd, drop = FALSE
  ] %*% derivative.variance
  dimnames(vcov) <- list(group.names, group.names)
  dimnames(variance.vcov) <- list(group.names, group.names)
  dimnames(mean.variance.vcov) <- list(group.names, group.names)
  contrast <- cbind(diag(length(group.names) - 1L), -1)
  list(
    estimate = estimate,
    se = setNames(sqrt(pmax(diag(vcov), 0)), group.names),
    vcov = vcov,
    variance = variance,
    variance.se = setNames(
      sqrt(pmax(diag(variance.vcov), 0)), group.names
    ),
    variance.vcov = variance.vcov,
    mean.variance.vcov = mean.variance.vcov,
    omnibus = .ZY.omnibus(estimate, vcov, contrast),
    omnibus.variance = .ZY.omnibus(variance, variance.vcov, contrast)
  )
}

.ML.ZY.path.specification <- function(L, times) {
  position <- L
  transition.index <- vector("list", max(times - 1L, 0L))
  if(times > 1L){
    for(t in 2:times){
      transition.index[[t - 1L]] <- vector("list", L)
      for(k in seq_len(L)){
        transition.index[[t - 1L]][[k]] <-
          position + seq_len(L - 1L) - 1L
        position <- position + L - 1L
      }
    }
  }
  list(
    initial.index = seq_len(L - 1L),
    transition.index = transition.index,
    npar = position - 1L
  )
}

.ML.ZY.path.initial <- function(modal, paths, L, specification) {
  assigned <- paths[modal, , drop = FALSE]
  initial.count <- tabulate(assigned[, 1L], nbins = L) + 0.5
  value <- log(initial.count[-L] / initial.count[L])
  if(ncol(paths) > 1L){
    for(t in 2:ncol(paths)){
      for(k in seq_len(L)){
        count <- tabulate(
          assigned[assigned[, t - 1L] == k, t], nbins = L
        ) + 0.5
        value[specification$transition.index[[t - 1L]][[k]]] <-
          log(count[-L] / count[L])
      }
    }
  }
  value
}

.ML.ZY.path.structure <- function(params, paths, L, specification) {
  times <- ncol(paths)
  probability.initial <- .ML.ZY.softmax(
    params[specification$initial.index]
  )
  probability.transition <- vector("list", max(times - 1L, 0L))
  if(times > 1L){
    for(t in 2:times){
      probability.transition[[t - 1L]] <- t(vapply(
        seq_len(L), function(k){
          .ML.ZY.softmax(
            params[specification$transition.index[[t - 1L]][[k]]]
          )
        }, numeric(L)
      ))
    }
  }

  log.prior <- log(pmax(probability.initial[paths[, 1L]], 1e-300))
  if(times > 1L){
    for(t in 2:times){
      log.prior <- log.prior + log(pmax(
        probability.transition[[t - 1L]][
          cbind(paths[, t - 1L], paths[, t])
        ],
        1e-300
      ))
    }
  }
  prior <- exp(log.prior - max(log.prior))
  prior <- prior / sum(prior)

  score <- matrix(0, nrow(paths), specification$npar)
  score[, specification$initial.index] <- sweep(
    outer(paths[, 1L], seq_len(L - 1L), `==`) + 0,
    2, probability.initial[-L], "-"
  )
  if(times > 1L){
    for(t in 2:times){
      for(k in seq_len(L)){
        selected <- paths[, t - 1L] == k
        score[selected, specification$transition.index[[t - 1L]][[k]]] <-
          sweep(
            outer(paths[selected, t], seq_len(L - 1L), `==`) + 0,
            2, probability.transition[[t - 1L]][k, -L], "-"
          )
      }
    }
  }
  list(
    prior = prior,
    probability.initial = probability.initial,
    probability.transition = probability.transition,
    score = score
  )
}

.ML.ZY.path.structural.information <- function(responsibility, paths,
                                                structure, specification, L) {
  information <- matrix(0, specification$npar, specification$npar)
  probability <- structure$probability.initial[-L]
  information[specification$initial.index, specification$initial.index] <-
    nrow(responsibility) * (
      diag(probability, nrow = L - 1L) - tcrossprod(probability)
    )
  if(ncol(paths) > 1L){
    for(t in 2:ncol(paths)){
      for(k in seq_len(L)){
        index <- specification$transition.index[[t - 1L]][[k]]
        probability <- structure$probability.transition[[t - 1L]][k, -L]
        mass <- sum(responsibility[, paths[, t - 1L] == k, drop = FALSE])
        information[index, index] <- mass * (
          diag(probability, nrow = L - 1L) - tcrossprod(probability)
        )
      }
    }
  }
  information
}

.ML.ZY.group.specification <- function(dependent.group, complete, G.path,
                                        group.names = NULL) {
  N <- length(complete)
  if(is.null(dependent.group)) dependent.group <- seq_len(G.path)
  if(is.matrix(dependent.group)){
    if(!identical(dim(dependent.group), c(N, G.path))){
      stop("dependent.group has incompatible dimensions")
    }
    dependent.group <- dependent.group[complete, , drop = FALSE]
  }else{
    if(length(dependent.group) != G.path){
      stop("dependent.group must contain one value per latent path")
    }
    dependent.group <- matrix(
      dependent.group, sum(complete), G.path, byrow = TRUE
    )
  }
  if(anyNA(dependent.group) || any(!is.finite(dependent.group)) ||
     any(dependent.group != as.integer(dependent.group)) ||
     any(dependent.group < 1L)){
    stop("dependent.group must contain positive integer group indices")
  }
  dependent.group <- matrix(
    as.integer(dependent.group), nrow(dependent.group), ncol(dependent.group)
  )
  G <- max(dependent.group)
  if(!identical(sort(unique(as.vector(dependent.group))), seq_len(G))){
    stop("dependent.group indices must be consecutive from 1")
  }
  if(is.null(group.names)) group.names <- paste0("Group.", seq_len(G))
  if(length(group.names) != G || anyNA(group.names) || any(!nzchar(group.names))){
    stop("group.names must contain one non-empty name per dependent-variable group")
  }
  list(value = dependent.group, G = G, names = as.character(group.names))
}

.ML.ZY.group.responsibility <- function(responsibility, dependent.group, G) {
  value <- matrix(0, nrow(responsibility), G)
  for(g in seq_len(G)){
    value[, g] <- rowSums(responsibility * (dependent.group == g))
  }
  value
}

.ML.ZY.group.prior <- function(prior, dependent.group, G) {
  vapply(seq_len(G), function(g){
    mean(rowSums(sweep(dependent.group == g, 2, prior, "*")))
  }, numeric(1))
}

.ML.ZY.path.fit.gaussian <- function(dependent.variable, classification,
                                      paths, L,
                                      method.regression = "Analytic",
                                      method.SE = "Analytic",
                                      modal.original = NULL,
                                      dependent.group = NULL,
                                      group.names = NULL) {
  complete <- !is.na(dependent.variable)
  y <- as.numeric(dependent.variable[complete])
  classification <- classification[complete, , drop = FALSE]
  modal.original <- modal.original[complete]
  if(!length(y) || any(!is.finite(y))){
    stop("gaussian dependent variables must contain finite observed values")
  }
  G.path <- ncol(classification)
  if(nrow(paths) != G.path){
    stop("classification columns must correspond to latent paths")
  }
  if(is.null(dependent.group) && is.null(group.names)){
    group.names <- colnames(classification)
  }
  group <- .ML.ZY.group.specification(
    dependent.group, complete, G.path, group.names
  )
  dependent.group <- group$value
  G <- group$G
  modal.group <- dependent.group[cbind(seq_along(y), modal.original)]
  specification <- .ML.ZY.path.specification(L, ncol(paths))
  mean.overall <- mean(y)
  sd.overall <- stats::sd(y)
  if(!is.finite(sd.overall) || sd.overall <= 0) sd.overall <- 1
  means <- vapply(seq_len(G), function(g){
    value <- mean(y[modal.group == g])
    if(is.finite(value)) value else mean.overall
  }, numeric(1))
  sds <- vapply(seq_len(G), function(g){
    value <- stats::sd(y[modal.group == g])
    if(is.finite(value) && value > 0) value else sd.overall
  }, numeric(1))
  par.ini <- c(
    .ML.ZY.path.initial(modal.original, paths, L, specification),
    means, log(sds)
  )
  index.structure <- seq_len(specification$npar)
  index.mean <- specification$npar + seq_len(G)
  index.log.sd <- specification$npar + G + seq_len(G)

  evaluate <- function(params, gradient = FALSE, components = FALSE) {
    structure <- .ML.ZY.path.structure(
      params[index.structure], paths, L, specification
    )
    mean.cur <- params[index.mean]
    sd.cur <- exp(params[index.log.sd])
    log.component <- sapply(seq_len(G.path), function(path){
      group.cur <- dependent.group[, path]
      log(pmax(classification[, path], 1e-300)) +
        log(pmax(structure$prior[path], 1e-300)) +
        stats::dnorm(
          y, mean.cur[group.cur], sd.cur[group.cur], log = TRUE
        )
    })
    maximum <- apply(log.component, 1L, max)
    scaled <- exp(log.component - maximum)
    denominator <- rowSums(scaled)
    responsibility <- scaled / denominator
    responsibility.group <- .ML.ZY.group.responsibility(
      responsibility, dependent.group, G
    )
    objective <- -sum(maximum + log(denominator))
    if(!gradient && !components) return(objective)

    score <- matrix(0, length(y), length(params))
    score[, index.structure] <- responsibility %*% structure$score
    for(g in seq_len(G)){
      residual <- y - mean.cur[g]
      score[, index.mean[g]] <-
        responsibility.group[, g] * residual / sd.cur[g]^2
      score[, index.log.sd[g]] <- responsibility.group[, g] *
        (-1 + residual^2 / sd.cur[g]^2)
    }
    list(
      objective = objective,
      gradient = -colSums(score),
      responsibility = responsibility,
      responsibility.group = responsibility.group,
      structure = structure,
      mean = mean.cur,
      sd = sd.cur
    )
  }

  lower <- c(
    rep(-10, specification$npar),
    rep(min(y) - 10 * sd.overall, G), rep(-10, G)
  )
  upper <- c(
    rep(10, specification$npar),
    rep(max(y) + 10 * sd.overall, G), rep(10, G)
  )
  optimization <- stats::optim(
    par.ini, fn = function(x) evaluate(x),
    gr = if(method.regression == "Analytic"){
      function(x) evaluate(x, gradient = TRUE)$gradient
    }else{
      function(x) numDeriv::grad(function(value) evaluate(value), x)
    },
    method = "L-BFGS-B", lower = lower, upper = upper,
    control = list(maxit = 5000, factr = 1e7)
  )
  fit <- evaluate(optimization$par, gradient = TRUE, components = TRUE)

  if(method.SE == "Analytic"){
    information <- matrix(0, length(optimization$par), length(optimization$par))
    information[index.structure, index.structure] <-
      .ML.ZY.path.structural.information(
        fit$responsibility, paths, fit$structure, specification, L
      )
    for(g in seq_len(G)){
      residual <- y - fit$mean[g]
      mass <- fit$responsibility.group[, g]
      information[index.mean[g], index.mean[g]] <-
        sum(mass) / fit$sd[g]^2
      cross.value <- sum(mass * 2 * residual / fit$sd[g]^2)
      information[index.mean[g], index.log.sd[g]] <- cross.value
      information[index.log.sd[g], index.mean[g]] <- cross.value
      information[index.log.sd[g], index.log.sd[g]] <-
        sum(mass * 2 * residual^2 / fit$sd[g]^2)
    }
    for(i in seq_along(y)){
      complete.score <- matrix(0, G.path, length(optimization$par))
      complete.score[, index.structure] <- fit$structure$score
      for(path in seq_len(G.path)){
        g <- dependent.group[i, path]
        residual <- y[i] - fit$mean[g]
        complete.score[path, index.mean[g]] <- residual / fit$sd[g]^2
        complete.score[path, index.log.sd[g]] <-
          -1 + residual^2 / fit$sd[g]^2
      }
      mean.score <- colSums(complete.score * fit$responsibility[i, ])
      information <- information - crossprod(
        sweep(complete.score, 2, mean.score, "-") *
          sqrt(fit$responsibility[i, ])
      )
    }
  }else{
    information <- numDeriv::hessian(
      function(x) evaluate(x), optimization$par
    )
  }
  vcov.parameters <- MASS::ginv(information)
  inference <- .ML.ZY.gaussian.inference(
    fit$mean, fit$sd, vcov.parameters, index.mean, index.log.sd,
    group$names
  )
  prior <- setNames(
    .ML.ZY.group.prior(fit$structure$prior, dependent.group, G),
    group$names
  )
  list(
    estimate = inference$estimate,
    se = inference$se,
    vcov = inference$vcov,
    vcov.parameters = vcov.parameters,
    information = information,
    parameters = optimization$par,
    structural.parameters = optimization$par[index.structure],
    omnibus = inference$omnibus,
    class.mass = prior * sum(complete),
    prior = prior,
    variance = inference$variance,
    variance.se = inference$variance.se,
    variance.vcov = inference$variance.vcov,
    mean.variance.vcov = inference$mean.variance.vcov,
    omnibus.variance = inference$omnibus.variance,
    admissible = TRUE,
    class.shift = mean(
      max.col(fit$responsibility.group, ties.method = "first") != modal.group
    ),
    observations = sum(complete),
    omitted = sum(!complete),
    levels = NULL,
    family = "gaussian",
    converged = optimization$convergence == 0L,
    iterations = unname(optimization$counts["function"])
  )
}

.ML.ZY.path.fit.categorical <- function(dependent.variable, classification,
                                         paths, L,
                                         method.regression = "Analytic",
                                         method.SE = "Analytic",
                                         levels.dependent.variable = NULL,
                                         modal.original = NULL,
                                         dependent.group = NULL,
                                         group.names = NULL) {
  complete <- !is.na(dependent.variable)
  if(is.null(levels.dependent.variable)){
    levels.dependent.variable <- levels(factor(dependent.variable[complete]))
  }
  y <- factor(
    dependent.variable[complete], levels = levels.dependent.variable
  )
  K <- nlevels(y)
  if(K < 2L){
    stop("categorical dependent variables must contain at least two observed categories")
  }
  indicator <- diag(K)[as.integer(y), , drop = FALSE]
  classification <- classification[complete, , drop = FALSE]
  modal.original <- modal.original[complete]
  G.path <- ncol(classification)
  if(nrow(paths) != G.path){
    stop("classification columns must correspond to latent paths")
  }
  if(is.null(dependent.group) && is.null(group.names)){
    group.names <- colnames(classification)
  }
  group <- .ML.ZY.group.specification(
    dependent.group, complete, G.path, group.names
  )
  dependent.group <- group$value
  G <- group$G
  modal.group <- dependent.group[cbind(seq_along(y), modal.original)]
  specification <- .ML.ZY.path.specification(L, ncol(paths))
  probability.dependent.variable <- matrix(0, G, K)
  for(g in seq_len(G)){
    count <- colSums(indicator[modal.group == g, , drop = FALSE]) + 0.5
    probability.dependent.variable[g, ] <- count / sum(count)
  }
  par.ini <- c(
    .ML.ZY.path.initial(modal.original, paths, L, specification),
    unlist(lapply(seq_len(G), function(g){
      log(probability.dependent.variable[g, -K] /
            probability.dependent.variable[g, K])
    }))
  )
  index.structure <- seq_len(specification$npar)
  index.dependent.variable <- lapply(seq_len(G), function(g){
    specification$npar + (g - 1L) * (K - 1L) + seq_len(K - 1L)
  })

  evaluate <- function(params, gradient = FALSE, components = FALSE) {
    structure <- .ML.ZY.path.structure(
      params[index.structure], paths, L, specification
    )
    probability.dependent.variable <- t(vapply(
      index.dependent.variable, function(index){
        .ML.ZY.softmax(params[index])
      }, numeric(K)
    ))
    log.component <- sapply(seq_len(G.path), function(path){
      group.cur <- dependent.group[, path]
      log(pmax(classification[, path], 1e-300)) +
        log(pmax(structure$prior[path], 1e-300)) +
        log(pmax(
          probability.dependent.variable[
            cbind(group.cur, as.integer(y))
          ],
          1e-300
        ))
    })
    maximum <- apply(log.component, 1L, max)
    scaled <- exp(log.component - maximum)
    denominator <- rowSums(scaled)
    responsibility <- scaled / denominator
    responsibility.group <- .ML.ZY.group.responsibility(
      responsibility, dependent.group, G
    )
    objective <- -sum(maximum + log(denominator))
    if(!gradient && !components) return(objective)

    score <- matrix(0, length(y), length(params))
    score[, index.structure] <- responsibility %*% structure$score
    for(g in seq_len(G)){
      score[, index.dependent.variable[[g]]] <-
        responsibility.group[, g] * sweep(
          indicator[, -K, drop = FALSE], 2,
          probability.dependent.variable[g, -K], "-"
        )
    }
    list(
      objective = objective,
      gradient = -colSums(score),
      responsibility = responsibility,
      responsibility.group = responsibility.group,
      structure = structure,
      probability.dependent.variable = probability.dependent.variable
    )
  }

  optimization <- stats::optim(
    par.ini, fn = function(x) evaluate(x),
    gr = if(method.regression == "Analytic"){
      function(x) evaluate(x, gradient = TRUE)$gradient
    }else{
      function(x) numDeriv::grad(function(value) evaluate(value), x)
    },
    method = "L-BFGS-B",
    lower = rep(-10, length(par.ini)), upper = rep(10, length(par.ini)),
    control = list(maxit = 5000, factr = 1e7)
  )
  fit <- evaluate(optimization$par, gradient = TRUE, components = TRUE)
  if(method.SE == "Analytic"){
    information <- matrix(0, length(optimization$par), length(optimization$par))
    information[index.structure, index.structure] <-
      .ML.ZY.path.structural.information(
        fit$responsibility, paths, fit$structure, specification, L
      )
    for(g in seq_len(G)){
      probability <- fit$probability.dependent.variable[g, -K]
      information[index.dependent.variable[[g]],
                  index.dependent.variable[[g]]] <-
        sum(fit$responsibility.group[, g]) * (
          diag(probability, nrow = K - 1L) - tcrossprod(probability)
        )
    }
    for(i in seq_along(y)){
      complete.score <- matrix(0, G.path, length(optimization$par))
      complete.score[, index.structure] <- fit$structure$score
      for(path in seq_len(G.path)){
        g <- dependent.group[i, path]
        complete.score[path, index.dependent.variable[[g]]] <-
          indicator[i, -K] - fit$probability.dependent.variable[g, -K]
      }
      mean.score <- colSums(complete.score * fit$responsibility[i, ])
      information <- information - crossprod(
        sweep(complete.score, 2, mean.score, "-") *
          sqrt(fit$responsibility[i, ])
      )
    }
  }else{
    information <- numDeriv::hessian(
      function(x) evaluate(x), optimization$par
    )
  }
  vcov.parameters <- MASS::ginv(information)
  transform <- matrix(0, G * K, length(optimization$par))
  for(g in seq_len(G)){
    probability <- fit$probability.dependent.variable[g, ]
    derivative <- rbind(
      diag(probability[-K], nrow = K - 1L) - tcrossprod(probability[-K]),
      -probability[K] * probability[-K]
    )
    rows <- (g - 1L) * K + seq_len(K)
    transform[rows, index.dependent.variable[[g]]] <- derivative
  }
  vcov <- transform %*% vcov.parameters %*% t(transform)
  estimate <- fit$probability.dependent.variable
  dimnames(estimate) <- list(group$names, levels.dependent.variable)
  se <- matrix(
    sqrt(pmax(diag(vcov), 0)), nrow = G, byrow = TRUE,
    dimnames = dimnames(estimate)
  )
  class.contrast <- cbind(diag(G - 1L), -1)
  contrast <- kronecker(class.contrast, diag(K))
  prior <- setNames(
    .ML.ZY.group.prior(fit$structure$prior, dependent.group, G),
    group$names
  )
  list(
    estimate = estimate,
    se = se,
    vcov = vcov,
    vcov.parameters = vcov.parameters,
    information = information,
    parameters = optimization$par,
    structural.parameters = optimization$par[index.structure],
    omnibus = .ZY.omnibus(as.vector(t(estimate)), vcov, contrast),
    class.mass = prior * sum(complete),
    variance = NULL,
    prior = prior,
    admissible = TRUE,
    class.shift = mean(
      max.col(fit$responsibility.group, ties.method = "first") != modal.group
    ),
    observations = sum(complete),
    omitted = sum(!complete),
    levels = levels.dependent.variable,
    family = "categorical",
    converged = optimization$convergence == 0L,
    iterations = unname(optimization$counts["function"])
  )
}

.ML.ZY.path.fit.dependent.variables <- function(
    dependent.variables, family, classification, paths, L,
    method.regression = "Analytic", method.SE = "Analytic",
    levels.dependent.variables = NULL, modal.original = NULL,
    dependent.group = NULL, group.names = NULL) {
  res <- vector("list", ncol(dependent.variables))
  names(res) <- names(dependent.variables)
  for(j in seq_len(ncol(dependent.variables))){
    if(family[j] == "gaussian"){
      res[[j]] <- .ML.ZY.path.fit.gaussian(
        dependent.variables[[j]], classification, paths, L,
        method.regression, method.SE, modal.original,
        dependent.group, group.names
      )
    }else{
      levels.cur <- if(is.null(levels.dependent.variables)){
        NULL
      }else levels.dependent.variables[[j]]
      res[[j]] <- .ML.ZY.path.fit.categorical(
        dependent.variables[[j]], classification, paths, L,
        method.regression, method.SE, levels.cur, modal.original,
        dependent.group, group.names
      )
    }
  }
  res
}

.ML.ZY.fit.gaussian <- function(dependent.variable, classification, prior,
                                method.regression = "Analytic",
                                method.SE = "Analytic",
                                modal.original = NULL) {
  complete <- !is.na(dependent.variable)
  y <- as.numeric(dependent.variable[complete])
  classification <- classification[complete, , drop = FALSE]
  if(!length(y) || any(!is.finite(y))){
    stop("gaussian dependent variables must contain finite observed values")
  }
  L <- ncol(classification)
  if(is.null(modal.original)){
    modal.original <- max.col(classification, ties.method = "first")
  }else{
    modal.original <- modal.original[complete]
  }
  mean.overall <- mean(y)
  sd.overall <- stats::sd(y)
  if(!is.finite(sd.overall) || sd.overall <= 0) sd.overall <- 1
  means <- vapply(seq_len(L), function(l){
    value <- mean(y[modal.original == l])
    if(is.finite(value)) value else mean.overall
  }, numeric(1))
  sds <- vapply(seq_len(L), function(l){
    value <- stats::sd(y[modal.original == l])
    if(is.finite(value) && value > 0) value else sd.overall
  }, numeric(1))
  prior <- pmax(prior / sum(prior), 1e-8)
  par.ini <- c(log(prior[-L] / prior[L]), means, log(sds))
  index.prior <- seq_len(L - 1L)
  index.mean <- (L - 1L) + seq_len(L)
  index.log.sd <- (2L * L - 1L) + seq_len(L)

  evaluate <- function(params, gradient = FALSE, components = FALSE) {
    probability <- .ML.ZY.softmax(params[index.prior])
    mean.cur <- params[index.mean]
    sd.cur <- exp(params[index.log.sd])
    log.component <- sapply(seq_len(L), function(l){
      log(pmax(classification[, l], 1e-300)) + log(probability[l]) +
        stats::dnorm(y, mean.cur[l], sd.cur[l], log = TRUE)
    })
    maximum <- apply(log.component, 1, max)
    scaled <- exp(log.component - maximum)
    denominator <- rowSums(scaled)
    responsibility <- scaled / denominator
    objective <- -sum(maximum + log(denominator))
    if(!gradient && !components) return(objective)

    score <- matrix(0, length(y), length(params))
    score[, index.prior] <- sweep(
      responsibility[, -L, drop = FALSE], 2, probability[-L], "-"
    )
    for(l in seq_len(L)){
      residual <- y - mean.cur[l]
      score[, index.mean[l]] <- responsibility[, l] * residual / sd.cur[l]^2
      score[, index.log.sd[l]] <- responsibility[, l] *
        (-1 + residual^2 / sd.cur[l]^2)
    }
    res <- list(
      objective = objective,
      gradient = -colSums(score),
      score = score,
      responsibility = responsibility,
      probability = probability,
      mean = mean.cur,
      sd = sd.cur
    )
    res
  }

  lower <- c(rep(-10, L - 1L), rep(min(y) - 10 * sd.overall, L), rep(-10, L))
  upper <- c(rep(10, L - 1L), rep(max(y) + 10 * sd.overall, L), rep(10, L))
  optimization <- stats::optim(
    par.ini,
    fn = function(x) evaluate(x),
    gr = if(method.regression == "Analytic"){
      function(x) evaluate(x, gradient = TRUE)$gradient
    }else NULL,
    method = "L-BFGS-B", lower = lower, upper = upper,
    control = list(maxit = 5000, factr = 1e7)
  )
  fit <- evaluate(optimization$par, gradient = TRUE, components = TRUE)

  if(method.SE == "Analytic"){
    probability <- fit$probability
    information <- matrix(0, length(optimization$par), length(optimization$par))
    V.prior <- diag(probability[-L], nrow = L - 1L) -
      tcrossprod(probability[-L])
    for(i in seq_along(y)){
      expected.negative.hessian <- matrix(
        0, length(optimization$par), length(optimization$par)
      )
      expected.negative.hessian[index.prior, index.prior] <- V.prior
      complete.score <- matrix(0, L, length(optimization$par))
      for(l in seq_len(L)){
        complete.score[l, index.prior] <- -probability[-L]
        if(l < L) complete.score[l, index.prior[l]] <-
          complete.score[l, index.prior[l]] + 1
        residual <- y[i] - fit$mean[l]
        complete.score[l, index.mean[l]] <- residual / fit$sd[l]^2
        complete.score[l, index.log.sd[l]] <- -1 + residual^2 / fit$sd[l]^2
        expected.negative.hessian[index.mean[l], index.mean[l]] <-
          expected.negative.hessian[index.mean[l], index.mean[l]] +
          fit$responsibility[i, l] / fit$sd[l]^2
        cross.value <- 2 * residual / fit$sd[l]^2
        expected.negative.hessian[index.mean[l], index.log.sd[l]] <-
          expected.negative.hessian[index.mean[l], index.log.sd[l]] +
          fit$responsibility[i, l] * cross.value
        expected.negative.hessian[index.log.sd[l], index.mean[l]] <-
          expected.negative.hessian[index.mean[l], index.log.sd[l]]
        expected.negative.hessian[index.log.sd[l], index.log.sd[l]] <-
          expected.negative.hessian[index.log.sd[l], index.log.sd[l]] +
          fit$responsibility[i, l] * 2 * residual^2 / fit$sd[l]^2
      }
      mean.score <- colSums(complete.score * fit$responsibility[i, ])
      variance.score <- crossprod(
        sweep(complete.score, 2, mean.score, "-") * sqrt(fit$responsibility[i, ])
      )
      information <- information + expected.negative.hessian - variance.score
    }
  }else{
    information <- numDeriv::hessian(function(x) evaluate(x), optimization$par)
  }
  vcov.parameters <- MASS::ginv(information)
  inference <- .ML.ZY.gaussian.inference(
    fit$mean, fit$sd, vcov.parameters, index.mean, index.log.sd,
    colnames(classification)
  )
  list(
    estimate = inference$estimate,
    se = inference$se,
    vcov = inference$vcov,
    vcov.parameters = vcov.parameters,
    information = information,
    omnibus = inference$omnibus,
    class.mass = setNames(
      fit$probability * sum(complete), names(inference$estimate)
    ),
    prior = setNames(fit$probability, names(inference$estimate)),
    variance = inference$variance,
    variance.se = inference$variance.se,
    variance.vcov = inference$variance.vcov,
    mean.variance.vcov = inference$mean.variance.vcov,
    omnibus.variance = inference$omnibus.variance,
    admissible = TRUE,
    class.shift = mean(
      max.col(fit$responsibility, ties.method = "first") != modal.original
    ),
    observations = sum(complete),
    omitted = sum(!complete),
    levels = NULL,
    family = "gaussian",
    converged = optimization$convergence == 0L,
    iterations = unname(optimization$counts["function"])
  )
}

.ML.ZY.fit.categorical <- function(dependent.variable, classification, prior,
                                   method.regression = "Analytic",
                                   method.SE = "Analytic",
                                   levels.dependent.variable = NULL,
                                   modal.original = NULL) {
  complete <- !is.na(dependent.variable)
  if(is.null(levels.dependent.variable)) levels.dependent.variable <- levels(factor(dependent.variable[complete]))
  y <- factor(dependent.variable[complete], levels = levels.dependent.variable)
  K <- nlevels(y)
  if(K < 2L) stop("categorical dependent variables must contain at least two observed categories")
  indicator <- diag(K)[as.integer(y), , drop = FALSE]
  classification <- classification[complete, , drop = FALSE]
  L <- ncol(classification)
  if(is.null(modal.original)){
    modal.original <- max.col(classification, ties.method = "first")
  }else{
    modal.original <- modal.original[complete]
  }
  prior <- pmax(prior / sum(prior), 1e-8)
  probability.dependent.variable <- matrix(0, L, K)
  for(l in seq_len(L)){
    count <- colSums(indicator[modal.original == l, , drop = FALSE]) + 0.5
    probability.dependent.variable[l, ] <- count / sum(count)
  }
  par.ini <- c(
    log(prior[-L] / prior[L]),
    unlist(lapply(seq_len(L), function(l){
      log(probability.dependent.variable[l, -K] / probability.dependent.variable[l, K])
    }))
  )
  index.prior <- seq_len(L - 1L)
  index.dependent.variable <- lapply(seq_len(L), function(l){
    (L - 1L) + (l - 1L) * (K - 1L) + seq_len(K - 1L)
  })

  evaluate <- function(params, gradient = FALSE, components = FALSE) {
    probability <- .ML.ZY.softmax(params[index.prior])
    probability.dependent.variable <- t(vapply(index.dependent.variable, function(index){
      .ML.ZY.softmax(params[index])
    }, numeric(K)))
    log.component <- sapply(seq_len(L), function(l){
      log(pmax(classification[, l], 1e-300)) + log(probability[l]) +
        log(pmax(probability.dependent.variable[l, as.integer(y)], 1e-300))
    })
    maximum <- apply(log.component, 1, max)
    scaled <- exp(log.component - maximum)
    denominator <- rowSums(scaled)
    responsibility <- scaled / denominator
    objective <- -sum(maximum + log(denominator))
    if(!gradient && !components) return(objective)

    score <- matrix(0, length(y), length(params))
    score[, index.prior] <- sweep(
      responsibility[, -L, drop = FALSE], 2, probability[-L], "-"
    )
    for(l in seq_len(L)){
      score[, index.dependent.variable[[l]]] <- responsibility[, l] * sweep(
        indicator[, -K, drop = FALSE], 2, probability.dependent.variable[l, -K], "-"
      )
    }
    list(
      objective = objective,
      gradient = -colSums(score),
      score = score,
      responsibility = responsibility,
      probability = probability,
      probability.dependent.variable = probability.dependent.variable
    )
  }

  optimization <- stats::optim(
    par.ini,
    fn = function(x) evaluate(x),
    gr = if(method.regression == "Analytic"){
      function(x) evaluate(x, gradient = TRUE)$gradient
    }else NULL,
    method = "L-BFGS-B",
    lower = rep(-10, length(par.ini)), upper = rep(10, length(par.ini)),
    control = list(maxit = 5000, factr = 1e7)
  )
  fit <- evaluate(optimization$par, gradient = TRUE, components = TRUE)
  if(method.SE == "Analytic"){
    information <- matrix(0, length(optimization$par), length(optimization$par))
    V.prior <- diag(fit$probability[-L], nrow = L - 1L) -
      tcrossprod(fit$probability[-L])
    for(i in seq_along(y)){
      expected.negative.hessian <- matrix(
        0, length(optimization$par), length(optimization$par)
      )
      expected.negative.hessian[index.prior, index.prior] <- V.prior
      complete.score <- matrix(0, L, length(optimization$par))
      for(l in seq_len(L)){
        complete.score[l, index.prior] <- -fit$probability[-L]
        if(l < L) complete.score[l, index.prior[l]] <-
          complete.score[l, index.prior[l]] + 1
        complete.score[l, index.dependent.variable[[l]]] <-
          indicator[i, -K] - fit$probability.dependent.variable[l, -K]
        V.dependent.variable <- diag(
          fit$probability.dependent.variable[l, -K], nrow = K - 1L
        ) -
          tcrossprod(fit$probability.dependent.variable[l, -K])
        expected.negative.hessian[index.dependent.variable[[l]], index.dependent.variable[[l]]] <-
          expected.negative.hessian[index.dependent.variable[[l]], index.dependent.variable[[l]]] +
          fit$responsibility[i, l] * V.dependent.variable
      }
      mean.score <- colSums(complete.score * fit$responsibility[i, ])
      variance.score <- crossprod(
        sweep(complete.score, 2, mean.score, "-") * sqrt(fit$responsibility[i, ])
      )
      information <- information + expected.negative.hessian - variance.score
    }
  }else{
    information <- numDeriv::hessian(function(x) evaluate(x), optimization$par)
  }
  vcov.parameters <- MASS::ginv(information)
  transform <- matrix(0, L * K, length(optimization$par))
  for(l in seq_len(L)){
    q <- fit$probability.dependent.variable[l, ]
    derivative <- rbind(
      diag(q[-K], nrow = K - 1L) - tcrossprod(q[-K]),
      -q[K] * q[-K]
    )
    rows <- (l - 1L) * K + seq_len(K)
    transform[rows, index.dependent.variable[[l]]] <- derivative
  }
  vcov <- transform %*% vcov.parameters %*% t(transform)
  estimate <- fit$probability.dependent.variable
  dimnames(estimate) <- list(colnames(classification), levels.dependent.variable)
  se <- matrix(
    sqrt(pmax(diag(vcov), 0)), nrow = L, byrow = TRUE,
    dimnames = dimnames(estimate)
  )
  class.contrast <- cbind(diag(L - 1L), -1)
  contrast <- kronecker(class.contrast, diag(K))
  list(
    estimate = estimate,
    se = se,
    vcov = vcov,
    vcov.parameters = vcov.parameters,
    information = information,
    omnibus = .ZY.omnibus(as.vector(t(estimate)), vcov, contrast),
    class.mass = setNames(
      fit$probability * sum(complete), rownames(estimate)
    ),
    variance = NULL,
    prior = setNames(fit$probability, rownames(estimate)),
    admissible = TRUE,
    class.shift = mean(
      max.col(fit$responsibility, ties.method = "first") != modal.original
    ),
    observations = sum(complete),
    omitted = sum(!complete),
    levels = levels.dependent.variable,
    family = "categorical",
    converged = optimization$convergence == 0L,
    iterations = unname(optimization$counts["function"])
  )
}

.ML.ZY.fit.dependent.variables <- function(dependent.variables, family, classification, prior,
                                method.regression = "Analytic",
                                method.SE = "Analytic",
                                levels.dependent.variables = NULL,
                                modal.original = NULL) {
  res <- vector("list", ncol(dependent.variables))
  names(res) <- names(dependent.variables)
  for(j in seq_len(ncol(dependent.variables))){
    if(family[j] == "gaussian"){
      res[[j]] <- .ML.ZY.fit.gaussian(
        dependent.variables[[j]], classification, prior,
        method.regression, method.SE, modal.original
      )
    }else{
      levels.cur <- if(is.null(levels.dependent.variables)) NULL else levels.dependent.variables[[j]]
      res[[j]] <- .ML.ZY.fit.categorical(
        dependent.variables[[j]], classification, prior,
        method.regression, method.SE, levels.cur, modal.original
      )
    }
  }
  res
}
