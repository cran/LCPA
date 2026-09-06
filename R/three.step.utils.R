.three.step.control.model.names <- c(
  "params", "constraint", "is.sort", "par.ini", "starts",
  "maxiter.warmup", "nrep"
)

.three.step.control.model <- function(control.model = NULL){
  defaults <- list(
    params = NULL,
    constraint = "VV",
    is.sort = TRUE,
    par.ini = "random",
    starts = 100,
    maxiter.warmup = 20,
    nrep = 20
  )
  if(is.null(control.model)) return(defaults)
  if(!is.list(control.model)){
    stop("control.model must be a named list", call. = FALSE)
  }
  control.names <- names(control.model)
  if(is.null(control.names) || any(!nzchar(control.names))){
    stop("Every control.model element must be named", call. = FALSE)
  }
  if(anyDuplicated(control.names)){
    stop("control.model elements must have unique names", call. = FALSE)
  }
  unknown <- setdiff(control.names, .three.step.control.model.names)
  if(length(unknown)){
    stop(
      "Unknown control.model element(s): ", paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }
  defaults[control.names] <- control.model
  defaults
}

.three.step.softmax <- function(eta) {
  eta <- eta - max(eta)
  prob <- exp(eta)
  prob / sum(prob)
}

.three.step.bootstrap.indices <- function(N) {
  sample.int(N, N, replace = TRUE)
}

.three.step.progress.start <- function(step, analysis, path = NULL,
                                       type.model = NULL, vis = TRUE) {
  if(!vis) return(invisible(NULL))
  ordinal <- c("first", "second", "third")
  if(!step %in% seq_along(ordinal)) stop("step must be 1, 2, or 3")
  prefix <- if(step == 1L) "" else "\n"
  subject <- if(step == 1L){
    paste0(analysis, ": ", type.model, " measurement model")
  }else{
    paste0(analysis, ": ", path)
  }
  cat(sprintf(
    "%sStarting the %s step for %s ...\n",
    prefix, ordinal[step], subject
  ))
  invisible(NULL)
}

.three.step.progress.bootstrap <- function(current, total, successful,
                                            progress.state, vis = TRUE) {
  if(!vis) return(invisible(NULL))
  .print.progress.line(
    sprintf(
      "  Bootstrap = %d/%d | Successful = %d",
      current, total, successful
    ),
    progress.state
  )
}

.three.step.progress.regression <- function(fit, method.3step, method.SE,
                                            model.label, vis = TRUE) {
  if(!vis) return(invisible(NULL))
  fits <- if(is.null(fit$fits)) list(fit) else fit$fits
  n.models <- length(fits)
  converged <- sum(vapply(fits, function(x) isTRUE(x$converged), logical(1)))
  iterations <- sum(vapply(fits, function(x){
    if(is.null(x$iterations) || !is.finite(x$iterations)) 0 else x$iterations
  }, numeric(1)))
  output <- sprintf(
    "  %s %s = %d | Converged = %d/%d | Iterations = %d",
    method.3step, model.label, n.models, converged, n.models, iterations
  )
  cat(output, "\n", sep = "")
  if(method.SE != "Bootstrap"){
    cat(sprintf("  %s Standard Errors were computed.\n", method.SE))
  }
  invisible(NULL)
}

.three.step.posterior <- function(response, params, type) {
  if(type == "LCA"){
    get.P.Z.Xn.LCA(
      response = response, par = params$par, P.Z = params$P.Z,
      category.levels = params$category.levels
    )
  }else{
    get.P.Z.Xn.LPA(
      response = response, means = params$means, covs = params$covs,
      P.Z = params$P.Z
    )
  }
}

.three.step.reorder.params <- function(params, position, type) {
  params$P.Z <- params$P.Z[position]
  if(type == "LCA"){
    params$par <- params$par[position, , , drop = FALSE]
  }else{
    params$means <- params$means[position, , drop = FALSE]
    params$covs <- params$covs[, , position, drop = FALSE]
  }
  params
}

.three.step.specification <- function(covariates, L, ref.class,
                                      covariates.time.cross = FALSE) {
  times <- length(covariates)
  if (times < 1L) {
    stop("covariates must contain at least one design matrix")
  }
  if (ref.class < 1L || ref.class > L) {
    stop("ref.class must be between 1 and L")
  }

  free.classes <- setdiff(seq_len(L), ref.class)
  position <- 1L
  p1 <- ncol(covariates[[1]])
  beta.index <- matrix(NA_integer_, p1, L)
  for (class.cur in free.classes) {
    beta.index[, class.cur] <- position:(position + p1 - 1L)
    position <- position + p1
  }

  gamma.index <- vector("list", max(times - 1L, 0L))
  if (times > 1L) {
    transition.p <- vapply(covariates[-1L], ncol, integer(1))
    if (covariates.time.cross && length(unique(transition.p)) != 1L) {
      stop("covariates.time.cross = TRUE requires identical numbers of transition covariates at times 2, ..., T")
    }

    for (t in 2:times) {
      if (covariates.time.cross && t > 2L) {
        gamma.index[[t - 1L]] <- gamma.index[[1L]]
        next
      }

      pt <- ncol(covariates[[t]])
      gamma.t <- vector("list", L)
      for (from.class in seq_len(L)) {
        index.cur <- matrix(NA_integer_, pt, L)
        for (to.class in free.classes) {
          index.cur[, to.class] <- position:(position + pt - 1L)
          position <- position + pt
        }
        gamma.t[[from.class]] <- index.cur
      }
      gamma.index[[t - 1L]] <- gamma.t
    }
  }

  list(
    beta.index = beta.index,
    gamma.index = gamma.index,
    free.classes = free.classes,
    npar = position - 1L,
    times = times
  )
}

.three.step.parameter.index <- function(covariates, L, ref.class,
                                        covariates.time.cross = FALSE) {
  specification <- .three.step.specification(
    covariates, L, ref.class, covariates.time.cross
  )
  index <- as.vector(
    specification$beta.index[, specification$free.classes, drop = FALSE]
  )
  if (specification$times > 1L) {
    for (t in 2:specification$times) {
      for (from.class in seq_len(L)) {
        index <- c(
          index,
          as.vector(specification$gamma.index[[t - 1L]][[from.class]][,
            specification$free.classes, drop = FALSE])
        )
      }
    }
  }
  index
}

.three.step.expand.parameters <- function(params, covariates, L, ref.class,
                                          covariates.time.cross = FALSE) {
  specification <- .three.step.specification(
    covariates, L, ref.class, covariates.time.cross
  )
  if (length(params) != specification$npar) {
    stop("params has length ", length(params), "; expected ", specification$npar)
  }
  index <- .three.step.parameter.index(
    covariates, L, ref.class, covariates.time.cross
  )
  unname(params[index])
}

.three.step.reference.na <- function(parameters, ref.class) {
  parameters$beta[, ref.class] <- NA_real_
  if (length(parameters$gamma) > 0L) {
    for (t in seq_along(parameters$gamma)) {
      for (from.class in seq_along(parameters$gamma[[t]])) {
        parameters$gamma[[t]][[from.class]][[ref.class]][] <- NA_real_
      }
    }
  }
  parameters
}

.three.step.bound.diagnostics <- function(params, lower, upper,
                                          bound.tol = 1e-4) {
  lower.hit <- is.finite(lower) & abs(params - lower) <= bound.tol
  upper.hit <- is.finite(upper) & abs(params - upper) <= bound.tol
  index <- which(lower.hit | upper.hit)
  side <- character(length(index))
  side[lower.hit[index]] <- "lower"
  side[upper.hit[index]] <- "upper"
  side[lower.hit[index] & upper.hit[index]] <- "both"

  list(
    hit = length(index) > 0L,
    n.hit = length(index),
    index = index,
    side = side,
    bound.tol = bound.tol
  )
}

.three.step.optimize <- function(par.ini, CEP, P.Z.Xns, Zs, covariates,
                                 covariates.time.cross, ref.class,
                                 lower, upper, tol, maxiter,
                                 gradient = c("analytic", "numeric"),
                                 vis = FALSE, progress.prefix = "",
                                 progress.state = NULL) {
  gradient <- match.arg(gradient)
  Log.Lik.history <- numeric(0)
  N <- nrow(covariates[[1L]])
  npar <- length(par.ini)
  if(is.null(progress.state)){
    progress.state <- .new.progress.state()
  }
  objective <- function(par) {
    get.Log.Lik.LTA.optim(
      par, CEP, P.Z.Xns, Zs, covariates,
      covariates.time.cross, ref.class
    )
  }
  evaluate <- function(par) {
    result <- if(gradient == "analytic") {
      get.Log.Lik.LTA.optim(
        par, CEP, P.Z.Xns, Zs, covariates,
        covariates.time.cross, ref.class, compute.gradient = TRUE
      )
    } else {
      list(
        objective = objective(par),
        gradient = numDeriv::grad(objective, par, method = "Richardson")
      )
    }

    Log.Lik <- -result$objective
    Log.Lik.history <<- c(Log.Lik.history, Log.Lik)
    if(vis){
      iteration <- length(Log.Lik.history)
      maxchg <- if(iteration > 1L) {
        abs(Log.Lik.history[iteration] - Log.Lik.history[iteration - 1L])
      } else {
        Inf
      }
      BIC <- -2 * Log.Lik + npar * log(N)
      .print.iteration.progress(
        iteration, maxchg, BIC,
        progress.prefix = progress.prefix,
        progress.state = progress.state
      )
    }
    result
  }

  result <- nloptr::nloptr(
    x0 = par.ini,
    eval_f = evaluate,
    lb = rep(lower, length(par.ini)),
    ub = rep(upper, length(par.ini)),
    opts = list(
      algorithm = "NLOPT_LD_LBFGS",
      xtol_rel = tol,
      ftol_rel = tol,
      maxeval = maxiter,
      print_level = 0
    )
  )
  list(result = result, Log.Lik.history = Log.Lik.history)
}

.three.step.derivatives <- function(params, CEP, Zs, covariates,
                                    covariates.time.cross = FALSE,
                                    ref.class,
                                    CEP.time.cross = FALSE) {
  L <- nrow(CEP[[1]])
  N <- nrow(covariates[[1]])
  times <- length(covariates)
  specification <- .three.step.specification(
    covariates, L, ref.class, covariates.time.cross
  )
  if (length(params) != specification$npar) {
    stop("params has length ", length(params), "; expected ", specification$npar)
  }
  if (length(CEP) != times || length(Zs) != times) {
    stop("CEP, Zs, and covariates must have the same length")
  }
  if (any(vapply(covariates, nrow, integer(1)) != N) ||
      any(vapply(Zs, length, integer(1)) != N)) {
    stop("All covariate matrices and modal-class vectors must contain the same N individuals")
  }
  if (any(vapply(CEP, nrow, integer(1)) != L) ||
      any(vapply(CEP, ncol, integer(1)) != L)) {
    stop("Every CEP matrix must be L by L")
  }

  latent.paths <- .make.latent.paths(L, times)
  n.paths <- nrow(latent.paths)
  q <- specification$npar
  free.classes <- specification$free.classes
  log.lik <- 0
  score <- numeric(q)
  score.individual <- matrix(0, N, q)
  information <- matrix(0, q, q)
  cep.npar <- L * L * if (CEP.time.cross) 1L else times
  cep.cross <- matrix(0, q, cep.npar)

  for (n in seq_len(N)) {
    x1 <- covariates[[1]][n, ]
    eta1 <- numeric(L)
    for (class.cur in free.classes) {
      eta1[class.cur] <- sum(
        x1 * params[specification$beta.index[, class.cur]]
      )
    }
    initial.prob <- .three.step.softmax(eta1)

    transition.prob <- vector("list", max(times - 1L, 0L))
    if (times > 1L) {
      for (t in 2:times) {
        xt <- covariates[[t]][n, ]
        probability.cur <- matrix(0, L, L)
        for (from.class in seq_len(L)) {
          eta.cur <- numeric(L)
          for (to.class in free.classes) {
            eta.cur[to.class] <- sum(
              xt * params[specification$gamma.index[[t - 1L]][[from.class]][, to.class]]
            )
          }
          probability.cur[from.class, ] <- .three.step.softmax(eta.cur)
        }
        transition.prob[[t - 1L]] <- probability.cur
      }
    }

    log.weight <- numeric(n.paths)
    complete.score <- matrix(0, n.paths, q)
    beta.index <- as.vector(
      specification$beta.index[, free.classes, drop = FALSE]
    )
    for (class.cur in free.classes) {
      index.cur <- specification$beta.index[, class.cur]
      complete.score[, index.cur] <- outer(
        as.numeric(latent.paths[, 1L] + 1L == class.cur) - initial.prob[class.cur],
        x1
      )
    }

    for (path.cur in seq_len(n.paths)) {
      class.cur <- latent.paths[path.cur, 1L] + 1L
      emission <- CEP[[1]][class.cur, Zs[[1]][n]]
      if (!is.finite(emission) || emission <= 0) {
        log.weight[path.cur] <- -Inf
      } else {
        log.weight[path.cur] <- log(initial.prob[class.cur]) + log(emission)
      }
    }

    if (times > 1L) {
      for (t in 2:times) {
        xt <- covariates[[t]][n, ]
        probability.cur <- transition.prob[[t - 1L]]
        for (from.class in seq_len(L)) {
          active <- as.numeric(latent.paths[, t - 1L] + 1L == from.class)
          for (to.class in free.classes) {
            index.cur <- specification$gamma.index[[t - 1L]][[from.class]][, to.class]
            complete.score[, index.cur] <- complete.score[, index.cur, drop = FALSE] +
              outer(
                active * (as.numeric(latent.paths[, t] + 1L == to.class) -
                            probability.cur[from.class, to.class]),
                xt
              )
          }
        }

        for (path.cur in seq_len(n.paths)) {
          if (!is.finite(log.weight[path.cur])) {
            next
          }
          from.class <- latent.paths[path.cur, t - 1L] + 1L
          to.class <- latent.paths[path.cur, t] + 1L
          emission <- CEP[[t]][to.class, Zs[[t]][n]]
          if (!is.finite(emission) || emission <= 0) {
            log.weight[path.cur] <- -Inf
          } else {
            log.weight[path.cur] <- log.weight[path.cur] +
              log(probability.cur[from.class, to.class]) + log(emission)
          }
        }
      }
    }

    max.log.weight <- max(log.weight)
    if (!is.finite(max.log.weight)) {
      stop("All latent paths have zero probability for individual ", n,
           "; check CEP orientation and zero cells")
    }
    scaled.weight <- exp(log.weight - max.log.weight)
    posterior.path <- scaled.weight / sum(scaled.weight)
    log.lik <- log.lik + max.log.weight + log(sum(scaled.weight))

    complete.information <- matrix(0, q, q)
    initial.jacobian <- diag(
      initial.prob[free.classes], nrow = length(free.classes)
    ) -
      tcrossprod(initial.prob[free.classes])
    complete.information[beta.index, beta.index] <-
      kronecker(initial.jacobian, tcrossprod(x1))

    if (times > 1L) {
      for (t in 2:times) {
        xt <- covariates[[t]][n, ]
        probability.cur <- transition.prob[[t - 1L]]
        for (from.class in seq_len(L)) {
          origin.prob <- sum(
            posterior.path[latent.paths[, t - 1L] + 1L == from.class]
          )
          index.cur <- as.vector(
            specification$gamma.index[[t - 1L]][[from.class]][,
              free.classes, drop = FALSE]
          )
          transition.jacobian <-
            diag(probability.cur[from.class, free.classes],
                 nrow = length(free.classes)) -
            tcrossprod(probability.cur[from.class, free.classes])
          complete.information[index.cur, index.cur] <-
            complete.information[index.cur, index.cur] +
            origin.prob * kronecker(transition.jacobian, tcrossprod(xt))
        }
      }
    }

    score.cur <- colSums(complete.score * posterior.path)
    score.individual[n, ] <- score.cur
    for (t in seq_len(times)) {
      observed.class <- Zs[[t]][n]
      offset <- if (CEP.time.cross) 0L else (t - 1L) * L * L
      for (class.cur in seq_len(L)) {
        cep.cur <- CEP[[t]][class.cur, observed.class]
        if (!is.finite(cep.cur) || cep.cur <= 0) {
          next
        }
        state.indicator <- as.numeric(
          latent.paths[, t] + 1L == class.cur
        )
        state.prob <- sum(posterior.path * state.indicator)
        score.state <- colSums(
          complete.score * (posterior.path * state.indicator)
        )
        index.cur <- offset + class.cur + (observed.class - 1L) * L
        cep.cross[, index.cur] <- cep.cross[, index.cur] +
          (score.state - score.cur * state.prob) / cep.cur
      }
    }
    weighted.score <- complete.score * sqrt(posterior.path)
    missing.information <- crossprod(weighted.score) - tcrossprod(score.cur)
    score <- score + score.cur
    information <- information + complete.information - missing.information
  }

  list(
    logLik = log.lik,
    score = score,
    score.individual = score.individual,
    cep.cross = cep.cross,
    information = (information + t(information)) / 2
  )
}

.three.step.cep.influence <- function(P.Z.Xns, Zs, CEP,
                                      CEP.time.cross = FALSE) {
  times <- length(P.Z.Xns)
  N <- nrow(P.Z.Xns[[1]])
  L <- ncol(P.Z.Xns[[1]])
  if (length(Zs) != times || length(CEP) != times ||
      any(vapply(P.Z.Xns, nrow, integer(1)) != N) ||
      any(vapply(P.Z.Xns, ncol, integer(1)) != L) ||
      any(vapply(Zs, length, integer(1)) != N)) {
    stop("P.Z.Xns, Zs, and CEP must describe the same N, L, and time points")
  }

  influence <- matrix(
    0, N, L * L * if (CEP.time.cross) 1L else times
  )
  if (CEP.time.cross) {
    for (class.cur in seq_len(L)) {
      denominator <- sum(vapply(
        P.Z.Xns, function(x) mean(x[, class.cur]), numeric(1)
      ))
      if (!is.finite(denominator) || denominator <= 0) {
        stop("CEP influence is undefined for a class with zero posterior prevalence")
      }
      for (observed.class in seq_len(L)) {
        contribution <- numeric(N)
        for (t in seq_len(times)) {
          contribution <- contribution +
            P.Z.Xns[[t]][, class.cur] *
            (as.numeric(Zs[[t]] == observed.class) -
               CEP[[t]][class.cur, observed.class])
        }
        index.cur <- class.cur + (observed.class - 1L) * L
        influence[, index.cur] <- contribution / denominator
      }
    }
  } else {
    for (t in seq_len(times)) {
      offset <- (t - 1L) * L * L
      for (class.cur in seq_len(L)) {
        posterior.cur <- P.Z.Xns[[t]][, class.cur]
        denominator <- mean(posterior.cur)
        if (!is.finite(denominator) || denominator <= 0) {
          stop("CEP influence is undefined for a class with zero posterior prevalence")
        }
        for (observed.class in seq_len(L)) {
          index.cur <- offset + class.cur + (observed.class - 1L) * L
          influence[, index.cur] <- posterior.cur *
            (as.numeric(Zs[[t]] == observed.class) -
               CEP[[t]][class.cur, observed.class]) / denominator
        }
      }
    }
  }
  influence
}

.three.step.bootstrap.target.se <- function(information, derivatives,
                                            P.Z.Xns, Zs, CEP,
                                            CEP.time.cross = FALSE,
                                            CEP.reestimated = TRUE,
                                            label) {
  conditional <- .three.step.information.se(information, label)
  if (is.null(conditional$vcov)) {
    return(conditional)
  }

  N <- nrow(derivatives$score.individual)
  score.adjusted <- derivatives$score.individual
  if (CEP.reestimated) {
    cep.influence <- .three.step.cep.influence(
      P.Z.Xns, Zs, CEP, CEP.time.cross
    )
    if (ncol(derivatives$cep.cross) != ncol(cep.influence)) {
      stop("CEP derivative and influence-function dimensions do not match")
    }
    cep.gradient <- derivatives$cep.cross / N
    score.adjusted <- score.adjusted +
      cep.influence %*% t(cep.gradient)
  }
  max.score.mean <- max(abs(colMeans(score.adjusted)))
  score.adjusted <- sweep(
    score.adjusted, 2L, colMeans(score.adjusted), "-"
  )
  meat <- crossprod(score.adjusted)
  vcov <- conditional$vcov %*% meat %*% conditional$vcov
  vcov <- (vcov + t(vcov)) / 2
  variance <- diag(vcov)
  variance[variance < 0] <- NA_real_

  list(
    se = sqrt(variance),
    vcov = vcov,
    information = conditional$information,
    conditional.vcov = conditional$vcov,
    meat = meat,
    diagnostics = c(
      conditional$diagnostics,
      list(
        target = if (CEP.reestimated) {
          "Nonparametric bootstrap with fixed Step 1 and re-estimated CEP"
        } else {
          "Nonparametric bootstrap with fixed Step 1 and fixed identity CEP"
        },
        sandwich = TRUE,
        cep.reestimated = CEP.reestimated,
        max.adjusted.score.mean = max.score.mean
      )
    )
  )
}

.three.step.information.se <- function(information, label) {
  if (is.null(information) || any(!is.finite(information))) {
    return(list(
      se = rep(NA_real_, nrow(information)), vcov = NULL,
      information = information,
      diagnostics = list(method = label, invertible = FALSE)
    ))
  }

  inverse <- invert.information(information)
  variance <- diag(inverse$vcov)
  variance[variance < 0] <- NA_real_
  list(
    se = sqrt(variance),
    vcov = inverse$vcov,
    information = inverse$information,
    diagnostics = c(list(method = label, invertible = TRUE), inverse[
      c("condition.initial", "condition.final", "adjusted.fraction", "threshold")
    ])
  )
}

.three.step.numeric.information <- function(params, objective) {
  tryCatch(
    numDeriv::hessian(
      func = objective,
      x = params,
      method = "Richardson",
      method.args = list(r = 4, v = 2)
    ),
    error = function(e1) {
      tryCatch(
        numDeriv::hessian(func = objective, x = params, method = "simple"),
        error = function(e2) NULL
      )
    }
  )
}
