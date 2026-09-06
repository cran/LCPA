.default.flexmix.control <- function(){
  list(
    maxiter = 1000L,
    minprior = 0,
    tol = 0
  )
}

.fit.flexmix.SEM <- function(formula, data, model, L, nrep, starts,
                             maxiter.warmup, vis, control.flexmix, type){
  nrep <- as.integer(nrep)
  starts <- as.integer(starts)
  maxiter.warmup <- as.integer(maxiter.warmup)
  maxiter <- as.integer(control.flexmix$maxiter)
  minprior <- as.numeric(control.flexmix$minprior)
  tol <- as.numeric(control.flexmix$tol)

  .validate.training.stages(starts, maxiter.warmup, nrep)
  if(length(maxiter) != 1L || is.na(maxiter) || maxiter < 1L){
    stop("control.flexmix$maxiter must be a positive integer.")
  }
  if(length(minprior) != 1L || !is.finite(minprior) || minprior < 0 || minprior >= 1){
    stop("control.flexmix$minprior must be in [0, 1).")
  }
  if(length(tol) != 1L || !is.finite(tol) || tol < 0){
    stop("control.flexmix$tol must be a non-negative number.")
  }

  seed <- sample.int(.Machine$integer.max, 1L)

  run.flexmix <- function(maxiter, seed.cur, cluster = NULL){
    set.seed(seed.cur)
    tryCatch({
      model.cur <- if(is.function(model)) model() else model
      fit.cur <- flexmix::flexmix(
        formula = formula,
        data = data,
        k = as.integer(L),
        cluster = cluster,
        model = model.cur,
        control = list(
          classify = "SEM",
          iter.max = as.integer(maxiter),
          minprior = minprior,
          tolerance = tol,
          verbose = 0
        )
      )
      if(is.function(model)) fit.cur <- .finalize.flexmix.LPA(fit.cur)
      fit.cur
      },
      error = function(e){ NULL }
    )
  }

  is.valid.fit <- function(fit.cur){
    if(is.null(fit.cur)) return(FALSE)
    tryCatch({
      posterior.cur <- as.matrix(flexmix::posterior(fit.cur))
      prior.cur <- as.numeric(methods::slot(fit.cur, "prior"))
      identical(as.integer(methods::slot(fit.cur, "k")), as.integer(L)) &&
        identical(dim(posterior.cur), c(nrow(data), as.integer(L))) &&
        all(is.finite(posterior.cur)) && all(rowSums(posterior.cur) > 0) &&
        length(prior.cur) == L && all(is.finite(prior.cur)) && all(prior.cur > 0) &&
        is.finite(as.numeric(methods::slot(fit.cur, "logLik")))
    }, error = function(e){ FALSE })
  }

  warmup.Log.Lik <- rep(-Inf, starts)
  warmup.fits <- vector("list", starts)
  warmup.progress.state <- .new.progress.state()
  for(s in seq_len(starts)){
    seed.cur <- as.integer((as.double(seed) + s - 2) %% .Machine$integer.max + 1)
    set.seed(seed.cur)
    cluster <- sample(rep(seq_len(L), length.out = nrow(data)))
    fit.cur <- run.flexmix(maxiter.warmup, seed.cur, cluster = cluster)
    if(is.valid.fit(fit.cur)){
      warmup.fits[[s]] <- fit.cur
      warmup.Log.Lik[s] <- as.numeric(methods::slot(fit.cur, "logLik"))
    }
    if(vis){
      .print.estimation.progress(
        "Warm", s, starts, warmup.Log.Lik[s],
        algorithm = "flexmix SEM", iterations = maxiter.warmup,
        progress.state = warmup.progress.state
      )
    }
  }
  if(vis) .end.estimation.progress()

  warmup.position <- order(warmup.Log.Lik, decreasing = TRUE)
  warmup.position <- warmup.position[is.finite(warmup.Log.Lik[warmup.position])]
  if(length(warmup.position) < nrep){
    stop(
      "flexmix warm-up produced only ", length(warmup.position),
      " finite ", type, " solutions from ", starts, " starts; nrep = ", nrep, "."
    )
  }
  warmup.position <- warmup.position[seq_len(nrep)]

  fit <- NULL
  best.Log.Lik <- -Inf
  Log.Lik.nrep <- rep(-Inf, nrep)
  replication.progress.state <- .new.progress.state()
  for(r in seq_len(nrep)){
    seed.cur <- as.integer((as.double(seed) + starts + r - 2) %% .Machine$integer.max + 1)
    warmup.fit <- warmup.fits[[warmup.position[r]]]
    cluster <- flexmix::posterior(warmup.fit)
    fit.cur <- run.flexmix(maxiter, seed.cur, cluster = cluster)
    if(!is.valid.fit(fit.cur)) fit.cur <- warmup.fit
    if(is.valid.fit(fit.cur)){
      Log.Lik.nrep[r] <- as.numeric(methods::slot(fit.cur, "logLik"))
      if(Log.Lik.nrep[r] > best.Log.Lik){
        fit <- fit.cur
        best.Log.Lik <- Log.Lik.nrep[r]
      }
    }
    if(vis){
      .print.estimation.progress(
        "Rep", r, nrep, Log.Lik.nrep[r], best.Log.Lik,
        progress.state = replication.progress.state
      )
    }
  }
  if(vis) .end.estimation.progress()

  if(is.null(fit)){
    stop("flexmix SEM failed to produce a finite ", type, " solution.")
  }
  list(fit = fit, Log.Lik.nrep = Log.Lik.nrep)
}

.repair.flexmix.covariance <- function(covariance, fallback){
  I <- nrow(covariance)
  .stabilize.LPA.covariances(
    array(covariance, dim = c(I, I, 1L)), "VV", fallback = fallback
  )$covs[, , 1L]
}

.flexmix.covariance.is.pd <- function(covs){
  L <- dim(covs)[3L]
  all(vapply(seq_len(L), function(l){
    covariance <- covs[, , l]
    all(is.finite(covariance)) &&
      max(abs(covariance - t(covariance))) < 1e-8 &&
      !is.null(tryCatch(chol(covariance), error = function(e){ NULL }))
  }, logical(1)))
}

.update.flexmix.covariances <- function(scatter, nk, covs.start,
                                        constraint){
  I <- dim(scatter)[1L]
  L <- dim(scatter)[3L]
  pooled <- apply(scatter, c(1, 2), sum) / sum(nk)
  fallback <- apply(covs.start, c(1, 2), mean)
  fallback <- .repair.flexmix.covariance(fallback, diag(I))

  if(I == 1L){
    shared <- if(is.character(constraint)){
      constraint %in% c("UE", "E0", "EE", "EV")
    }else{
      any(vapply(constraint, function(x){ all(as.integer(x) == 1L) }, logical(1)))
    }
    values <- if(shared) rep(pooled[1, 1], L) else scatter[1, 1, ] / nk
    if(any(!is.finite(values))) values <- rep(fallback[1, 1], L)
    covs <- array(values, dim = c(1, 1, L))
    for(l in seq_len(L)){
      covs[, , l] <- .repair.flexmix.covariance(covs[, , l], fallback)
    }
    return(covs)
  }

  covs <- array(0, dim = c(I, I, L))
  if(is.character(constraint) && constraint == "E0"){
    for(l in seq_len(L)){
      covs[, , l] <- .repair.flexmix.covariance(diag(diag(pooled), I), fallback)
    }
  }else if(is.character(constraint) && constraint == "V0"){
    for(l in seq_len(L)){
      covariance <- diag(diag(scatter[, , l]) / nk[l], I)
      covs[, , l] <- .repair.flexmix.covariance(covariance, fallback)
    }
  }else if(is.character(constraint) && constraint == "EE"){
    pooled <- .repair.flexmix.covariance(pooled, fallback)
    for(l in seq_len(L)) covs[, , l] <- pooled
  }else if(is.character(constraint) && constraint == "VV"){
    for(l in seq_len(L)){
      covs[, , l] <- .repair.flexmix.covariance(
        scatter[, , l] / nk[l], fallback
      )
    }
  }else{
    pairs <- which(lower.tri(matrix(TRUE, I, I), diag = TRUE), arr.ind = TRUE)
    pair.keys <- apply(pairs, 1, function(x){ paste(sort(x), collapse = ":") })
    if(is.character(constraint)){
      shared <- if(constraint == "VE") pairs[, 1] != pairs[, 2] else pairs[, 1] == pairs[, 2]
    }else{
      constraint.keys <- unique(vapply(constraint, function(x){
        paste(sort(as.integer(x)), collapse = ":")
      }, character(1)))
      shared <- pair.keys %in% constraint.keys
    }

    param.index <- matrix(0L, nrow = nrow(pairs), ncol = L)
    ntheta <- 0L
    for(p in seq_len(nrow(pairs))){
      if(shared[p]){
        ntheta <- ntheta + 1L
        param.index[p, ] <- ntheta
      }else{
        param.index[p, ] <- ntheta + seq_len(L)
        ntheta <- ntheta + L
      }
    }

    make.covs <- function(theta){
      result <- array(0, dim = c(I, I, L))
      for(p in seq_len(nrow(pairs))){
        i <- pairs[p, 1]
        j <- pairs[p, 2]
        for(l in seq_len(L)){
          result[i, j, l] <- result[j, i, l] <- theta[param.index[p, l]]
        }
      }
      result
    }

    make.theta <- function(covs){
      theta <- numeric(ntheta)
      for(p in seq_len(nrow(pairs))){
        i <- pairs[p, 1]
        j <- pairs[p, 2]
        if(shared[p]){
          theta[param.index[p, 1L]] <- sum(nk * covs[i, j, ]) / sum(nk)
        }else{
          theta[param.index[p, ]] <- covs[i, j, ]
        }
      }
      theta
    }

    objective <- function(theta){
      covs <- make.covs(theta)
      value <- 0
      for(l in seq_len(L)){
        R <- tryCatch(chol(covs[, , l]), error = function(e){ NULL })
        if(is.null(R)) return(.Machine$double.xmax^0.25)
        inverse <- chol2inv(R)
        value <- value + nk[l] * 2 * sum(log(diag(R))) +
          sum(scatter[, , l] * inverse)
      }
      value / 2
    }

    gradient <- function(theta){
      covs <- make.covs(theta)
      grad <- numeric(ntheta)
      for(l in seq_len(L)){
        R <- tryCatch(chol(covs[, , l]), error = function(e){ NULL })
        if(is.null(R)) return(rep(0, ntheta))
        inverse <- chol2inv(R)
        G <- (nk[l] * inverse - inverse %*% scatter[, , l] %*% inverse) / 2
        for(p in seq_len(nrow(pairs))){
          i <- pairs[p, 1]
          j <- pairs[p, 2]
          value <- if(i == j) G[i, j] else 2 * G[i, j]
          index <- param.index[p, l]
          grad[index] <- grad[index] + value
        }
      }
      grad
    }

    unconstrained <- array(0, dim = c(I, I, L))
    safe <- array(0, dim = c(I, I, L))
    for(l in seq_len(L)) unconstrained[, , l] <- scatter[, , l] / nk[l]
    for(i in seq_len(I)){
      p <- which(pairs[, 1] == i & pairs[, 2] == i)
      values <- unconstrained[i, i, ]
      if(shared[p]) values[] <- sum(scatter[i, i, ]) / sum(nk)
      if(any(!is.finite(values))) values <- covs.start[i, i, ]
      values[values <= 0] <- covs.start[i, i, ][values <= 0]
      if(shared[p]) values[] <- sum(nk * values) / sum(nk)
      safe[i, i, ] <- values
    }

    candidates <- list(make.theta(covs.start), make.theta(unconstrained), make.theta(safe))
    valid <- vapply(candidates, function(theta){
      .flexmix.covariance.is.pd(make.covs(theta))
    }, logical(1))
    if(!any(valid)) return(NULL)
    candidates <- candidates[valid]
    values <- vapply(candidates, objective, numeric(1))
    theta.start <- candidates[[which.min(values)]]
    objective.start <- min(values)
    fit <- tryCatch(
      stats::optim(
        theta.start, objective, gradient, method = "BFGS",
        control = list(maxit = 100, reltol = 1e-8)
      ),
      error = function(e){ NULL }
    )
    if(is.null(fit)) return(make.covs(theta.start))
    covs <- make.covs(fit$par)
    if(!.flexmix.covariance.is.pd(covs) ||
       objective(fit$par) > objective.start + 1e-8 * (1 + abs(objective.start))){
      covs <- make.covs(theta.start)
    }
  }

  if(!.flexmix.covariance.is.pd(covs)) return(NULL)
  covs
}

.flexmix.LCA.model <- function(response, poly.value){
  response <- as.matrix(response)
  storage.mode(response) <- "integer"
  N <- nrow(response)
  I <- ncol(response)
  poly.value <- as.integer(poly.value)
  poly.max <- max(poly.value)
  category.offset <- c(0L, cumsum(poly.value))[seq_len(I)]
  category.index <- sweep(response + 1L, 2L, category.offset, "+")
  category.indicator <- matrix(0, N, sum(poly.value))
  category.indicator[cbind(rep(seq_len(N), I), as.vector(category.index))] <- 1
  item.index <- rep(seq_len(I), poly.value)

  model <- methods::new(
    "FLXMC", weighted = TRUE, formula = . ~ ., dist = "multinomial",
    name = "LCPA categorical latent class model"
  )
  methods::slot(model, "defineComponent") <- function(para){
    probability <- para$probability
    probability.vector <- probability[cbind(item.index, sequence(poly.value))]
    log.probability <- log(probability.vector)
    logLik <- function(x, y){
      rowSums(matrix(log.probability[category.index], N, I))
    }
    expected <- vapply(seq_len(I), function(i){
      sum((0:(poly.value[i] - 1L)) * probability[i, seq_len(poly.value[i])])
    }, numeric(1))
    predict <- function(x, ...){
      matrix(expected, nrow(x), I, byrow = TRUE)
    }
    methods::new(
      "FLXcomponent", parameters = list(probability = probability),
      df = sum(poly.value - 1L), logLik = logLik, predict = predict
    )
  }
  methods::slot(model, "fit") <- function(x, y, w, ...){
    nk <- sum(w)
    if(!is.finite(nk) || nk <= 0){
      stop("empty class in flexmix SEM S-step")
    }
    counts <- as.numeric(crossprod(w, category.indicator))
    probability <- matrix(NA_real_, I, poly.max)
    position <- 0L
    for(i in seq_len(I)){
      index <- position + seq_len(poly.value[i])
      probability.i <- pmax(counts[index] / nk, .Machine$double.xmin)
      probability[i, seq_len(poly.value[i])] <- probability.i / sum(probability.i)
      position <- position + poly.value[i]
    }
    methods::slot(model, "defineComponent")(list(probability = probability))
  }
  model
}

.flexmix.LPA.model <- function(response, L, constraint){
  response <- as.matrix(response)
  N <- nrow(response)
  I <- ncol(response)
  state <- new.env(parent = emptyenv())
  state$index <- 0L
  state$weights <- matrix(0, N, L)
  covariance.global <- matrix(stats::cov(response), I, I)
  if(any(!is.finite(covariance.global))){
    variances <- apply(response, 2, stats::var)
    variances[!is.finite(variances) | variances <= 0] <- 1
    covariance.global <- diag(variances, I)
  }
  covariance.global <- (covariance.global + t(covariance.global)) / 2
  covariance.global <- .repair.flexmix.covariance(covariance.global, diag(I))
  state$params <- list(
    means = matrix(colMeans(response), L, I, byrow = TRUE),
    covs = array(rep(covariance.global, L), dim = c(I, I, L))
  )
  covariance.df <- get.npar.LPA(I, L, constraint) - L * I - (L - 1L)

  model <- methods::new(
    "FLXMC", weighted = TRUE, formula = . ~ ., dist = "mvnorm",
    name = "LCPA constrained Gaussian clustering"
  )
  methods::slot(model, "defineComponent") <- function(para){
    methods::new("FLXcomponent", parameters = list(), df = 0)
  }
  methods::slot(model, "fit") <- function(x, y, w, ...){
    state$index <- state$index + 1L
    component <- state$index
    state$weights[, component] <- w
    if(component == L){
      nk <- colSums(state$weights)
      empty <- which(!is.finite(nk) | nk <= 0)
      if(length(empty) > 0L){
        available <- order(apply(state$weights, 1L, max))
        for(j in seq_along(empty)){
          state$weights[available[(j - 1L) %% N + 1L], empty[j]] <- 1
        }
        nk <- colSums(state$weights)
      }
      means <- sweep(t(state$weights) %*% response, 1, nk, "/")
      scatter <- array(0, dim = c(I, I, L))
      for(l in seq_len(L)){
        dev <- sweep(response, 2, means[l, ], "-")
        scatter[, , l] <- crossprod(dev, dev * state$weights[, l])
      }
      covs <- .update.flexmix.covariances(
        scatter, nk, state$params$covs, constraint
      )
      if(is.null(covs)) stop("invalid constrained covariance in flexmix SEM M-step")
      state$params <- list(means = means, covs = covs)
      state$index <- 0L
    }

    component.index <- component
    logLik <- function(x, y){
      mvn_log_density_cpp(
        as.matrix(y), state$params$means[component.index, ],
        state$params$covs[, , component.index], FALSE
      )
    }
    predict <- function(x, ...){
      matrix(state$params$means[component.index, ], nrow(x), I, byrow = TRUE)
    }
    methods::new(
      "FLXcomponent",
      parameters = list(component = component.index),
      df = I + if(component.index == 1L) covariance.df else 0,
      logLik = logLik,
      predict = predict
    )
  }
  model
}

.finalize.flexmix.LPA <- function(fit){
  model <- methods::slot(fit, "model")[[1L]]
  state <- get("state", envir = environment(methods::slot(model, "fit")), inherits = TRUE)
  components <- methods::slot(fit, "components")
  for(l in seq_len(methods::slot(fit, "k"))){
    component <- components[[l]][[1L]]
    methods::slot(component, "parameters") <- list(
      center = state$params$means[l, ],
      cov = state$params$covs[, , l]
    )
    components[[l]][[1L]] <- component
  }
  methods::slot(fit, "components") <- components
  fit
}
