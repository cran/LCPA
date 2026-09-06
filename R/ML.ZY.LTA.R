#' Longitudinal classification-error-corrected ML distal model
#'
#' @details
#' The observed likelihood sums over every complete latent-state path. The
#' classification term is the product of the time-specific CEP probabilities,
#' and the latent-path probability is determined by the initial-state and
#' transition probabilities. With `dependent.variable.structure = "State"`,
#' paths sharing the selected-time latent state share the same dependent-
#' variable parameters. With `dependent.variable.structure = "Path"`, each
#' complete path has its own dependent-variable parameters.
#' Time-specific models retain separate parameters. With
#' `dependent.variable.time.cross = TRUE`, selected time records share
#' dependent-variable parameters;
#' participant bootstrap is recommended and preserves their dependence by
#' resampling subjects across all times. Class-shift diagnostics compare the
#' third-step responsibilities with the Step 2 modal assignments.
#'
#' @noRd
ML.ZY.LTA <- function(responses, dependent.variables, L = 2,
                      type.model = "LCA",
                      family = "gaussian",
                      dependent.variable.structure = c("State", "Path"),
                      dependent.variable.time = NULL,
                      dependent.variable.time.cross = FALSE,
                      CEP.time.cross = FALSE,
                      CEP.error = TRUE,
                      par.ini = "random",
                      params = NULL, step1.pool = FALSE, is.sort = TRUE,
                      constraint = "VV",
                      method.model = "EM",
                      method.regression = "Analytic",
                      method.SE = "Bootstrap", nrep.bootstrap = 100,
                      starts = 100, maxiter.warmup = 20, nrep = 20,
                      vis = TRUE,
                      control.EM = NULL,
                      control.Mplus = NULL,
                      control.NNE = NULL,
                      control.Rmixmod = NULL) {

  call <- match.call()
  dependent.variable.structure <- match.arg(dependent.variable.structure)
  method.regression <- match.arg(method.regression, c("Analytic", "Numeric"))
  method.SE <- match.arg(method.SE, c("Analytic", "Numeric", "Bootstrap"))
  if(method.SE == "Bootstrap" && nrep.bootstrap < 2L){
    stop("nrep.bootstrap must be at least 2 when method.SE = 'Bootstrap'")
  }

  res <- ML.XZ.LTA(
    responses = responses, L = L,
    ref.class = L, type.model = type.model,
    covariates = NULL,
    CEP.time.cross = CEP.time.cross, CEP.error = CEP.error,
    covariates.time.cross = FALSE,
    par.ini = par.ini, params = params, step1.pool = step1.pool,
    is.sort = is.sort, constraint = constraint,
    method.model = method.model,
    method.regression = "Analytic", method.SE = "Analytic",
    starts = starts, maxiter.warmup = maxiter.warmup, nrep = nrep,
    vis = vis,
    control.EM = control.EM, control.Mplus = control.Mplus,
    control.NNE = control.NNE, control.Rmixmod = control.Rmixmod,
    .progress.path = "Z -> Y", .progress.step3 = FALSE
  )

  times <- length(res$P.Z.Xns)
  N <- nrow(res$P.Z.Xns[[1L]])
  dependent.variable.specification <- .ZY.LTA.dependent.variables(
    dependent.variables, family, dependent.variable.time, times, N
  )
  selected <- dependent.variable.specification$selected
  dependent.variables.list <- dependent.variable.specification$dependent.variables
  family.list <- dependent.variable.specification$family
  res$CEP <- if(CEP.error){
    get.CEP(res$P.Z.Xns, CEP.time.cross = CEP.time.cross)
  }else replicate(times, diag(L), simplify = FALSE)
  grouping <- function(posterior.cur, CEP.cur){
    modal.time <- lapply(
      posterior.cur, max.col, ties.method = "first"
    )
    classification.time <- lapply(seq_len(times), function(t){
      value <- .ML.ZY.classification(modal.time[[t]], CEP.cur[[t]])
      colnames(value) <- colnames(posterior.cur[[t]])
      value
    })
    classification.path <- .ZY.path.matrix(classification.time)
    posterior.path <- .ZY.path.matrix(posterior.cur)
    modal.path <- .ZY.path.modal(modal.time, L)
    list(
      classification = if(dependent.variable.structure == "State"){
        classification.time
      }else{
        replicate(times, classification.path$value, simplify = FALSE)
      },
      posterior = if(dependent.variable.structure == "State"){
        posterior.cur
      }else{
        replicate(times, posterior.path$value, simplify = FALSE)
      },
      modal = if(dependent.variable.structure == "State"){
        modal.time
      }else{
        replicate(times, modal.path, simplify = FALSE)
      },
      classification.path = classification.path$value,
      posterior.path = posterior.path$value,
      modal.path = modal.path,
      paths = classification.path$paths
    )
  }
  grouping.initial <- grouping(res$P.Z.Xns, res$CEP)
  classifications <- grouping.initial$classification
  if(vis){
    cat(sprintf(
      "  Classification-error-corrected ML weights were prepared for the %s dependent-variable model.\n",
      tolower(dependent.variable.structure)
    ))
  }
  levels.list <- lapply(selected, function(t){
    lapply(seq_along(dependent.variables.list[[t]]), function(j){
      if(family.list[[t]][j] == "categorical"){
        levels(factor(dependent.variables.list[[t]][[j]]))
      }else NULL
    })
  })

  fit.dependent.variables <- function(dependent.variables.cur, posterior.cur, CEP.cur,
                           selected.cur, levels.cur = NULL,
                           method.SE.cur = "Analytic"){
    grouping.cur <- grouping(posterior.cur, CEP.cur)
    fit.group <- function(dependent.variables.fit, family.fit,
                          classification.fit, modal.fit, levels.fit,
                          dependent.group.fit, group.names.fit){
      .ML.ZY.path.fit.dependent.variables(
        dependent.variables.fit, family.fit, classification.fit,
        grouping.cur$paths, L, method.regression, method.SE.cur,
        levels.fit, modal.fit, dependent.group.fit, group.names.fit
      )
    }
    if(dependent.variable.time.cross){
      names.cur <- lapply(dependent.variables.cur[selected.cur], names)
      if(length(unique(vapply(names.cur, paste, character(1), collapse = "\r"))) != 1L){
        stop("dependent.variable.time.cross = TRUE requires identical dependent-variable columns across selected time points")
      }
      family.key <- vapply(family.list[selected.cur], paste, character(1), collapse = "\r")
      if(length(unique(family.key)) != 1L){
        stop("dependent.variable.time.cross = TRUE requires identical dependent-variable families across selected time points")
      }
      dependent.variables.pooled <- do.call(rbind, dependent.variables.cur[selected.cur])
      classification.pooled <- do.call(rbind, replicate(
        length(selected.cur), grouping.cur$classification.path,
        simplify = FALSE
      ))
      modal.pooled <- rep(grouping.cur$modal.path, length(selected.cur))
      if(dependent.variable.structure == "Path"){
        dependent.group.pooled <- seq_len(nrow(grouping.cur$paths))
        group.names.pooled <- colnames(grouping.cur$classification.path)
      }else{
        dependent.group.pooled <- do.call(rbind, lapply(selected.cur, function(t){
          matrix(
            grouping.cur$paths[, t], N, nrow(grouping.cur$paths),
            byrow = TRUE
          )
        }))
        group.names.pooled <- .latent.group.names(L, type.model)
      }
      list(pooled = fit.group(
        dependent.variables.pooled, family.list[[selected.cur[1L]]],
        classification.pooled, modal.pooled,
        if(is.null(levels.cur)) NULL else levels.cur[[1L]],
        dependent.group.pooled, group.names.pooled
      ))
    }else{
      value <- lapply(seq_along(selected.cur), function(j){
        t <- selected.cur[j]
        if(dependent.variable.structure == "Path"){
          dependent.group.cur <- seq_len(nrow(grouping.cur$paths))
          group.names.cur <- colnames(grouping.cur$classification.path)
        }else{
          dependent.group.cur <- grouping.cur$paths[, t]
          group.names.cur <- .latent.group.names(L, type.model)
        }
        fit.group(
          dependent.variables.cur[[t]], family.list[[t]],
          grouping.cur$classification.path, grouping.cur$modal.path,
          if(is.null(levels.cur)) NULL else levels.cur[[j]],
          dependent.group.cur, group.names.cur
        )
      })
      names(value) <- paste0("t", selected.cur)
      value
    }
  }

  method.SE.fit <- if(method.SE == "Bootstrap") "Analytic" else method.SE
  .three.step.progress.start(
    3L, "LTA", path = "Z -> Y", vis = vis
  )
  models <- fit.dependent.variables(
    dependent.variables.list, res$P.Z.Xns, res$CEP, selected,
    levels.list, method.SE.fit
  )
  .ZY.progress.models(models, "ML", method.SE, vis)
  vcov <- lapply(models, function(x) NULL)
  SE.diagnostics <- list(
    method = method.SE,
    class.shift = lapply(models, function(x){
      vapply(x, `[[`, numeric(1), "class.shift")
    })
  )

  if(method.SE == "Bootstrap"){
    if(vis) cat("  Bootstrapping for Standard Errors ...\n")
    estimates.bootstrap <- lapply(models, function(x){
      matrix(NA_real_, nrep.bootstrap, length(.ZY.vector(x)))
    })
    progress.state <- .new.progress.state()
    successful.bootstrap <- 0L
    for(bs in seq_len(nrep.bootstrap)){
      samples.cur <- .three.step.bootstrap.indices(N)
      posterior.cur <- lapply(res$P.Z.Xns, function(x){
        x[samples.cur, , drop = FALSE]
      })
      CEP.cur <- if(CEP.error){
        get.CEP(posterior.cur, CEP.time.cross = CEP.time.cross)
      }else{
        replicate(times, diag(L), simplify = FALSE)
      }
      dependent.variables.cur <- dependent.variables.list
      for(t in selected){
        dependent.variables.cur[[t]] <- dependent.variables.list[[t]][samples.cur, , drop = FALSE]
      }
      models.cur <- tryCatch(
        fit.dependent.variables(
          dependent.variables.cur, posterior.cur, CEP.cur, selected,
          levels.list, "Analytic"
        ),
        error = function(e) NULL
      )
      if(!is.null(models.cur)){
        for(g in seq_along(models)){
          estimates.bootstrap[[g]][bs, ] <- .ZY.vector(models.cur[[g]])
        }
        successful.bootstrap <- successful.bootstrap + as.integer(all(
          vapply(models.cur, function(x){
            all(is.finite(.ZY.vector(x)))
          }, logical(1))
        ))
      }
      .three.step.progress.bootstrap(
        bs, nrep.bootstrap, successful.bootstrap, progress.state, vis
      )
    }
    if(vis) .end.estimation.progress()
    diagnostics.bootstrap <- vector("list", length(models))
    names(diagnostics.bootstrap) <- names(models)
    for(g in seq_along(models)){
      bootstrap <- .ZY.bootstrap(models[[g]], estimates.bootstrap[[g]])
      models[[g]] <- bootstrap$models
      vcov[[g]] <- bootstrap$vcov
      diagnostics.bootstrap[[g]] <- bootstrap$diagnostics
    }
    SE.diagnostics <- list(
      method = "Bootstrap", groups = diagnostics.bootstrap,
      class.shift = lapply(models, function(x){
        vapply(x, `[[`, numeric(1), "class.shift")
      })
    )
    .ZY.progress.bootstrap.complete(SE.diagnostics, vis)
  }

  class.shift <- unlist(lapply(models, function(x){
    vapply(x, `[[`, numeric(1), "class.shift")
  }), use.names = FALSE)
  if(any(class.shift > 0.20)){
    warning("The ML distal model changed more than 20% of modal assignments; inspect class-shift diagnostics and consider method.3step = 'BCH'")
  }
  if(dependent.variable.time.cross && method.SE != "Bootstrap" && length(selected) > 1L){
    warning("Analytic and numeric ML standard errors treat pooled time records as independent; use method.SE = 'Bootstrap' for participant-level uncertainty")
  }
  res$dependent.variables <- models
  res$beta <- NULL
  res$gamma <- NULL
  res$beta.se <- NULL
  res$gamma.se <- NULL
  res$beta.Z.sta <- NULL
  res$gamma.Z.sta <- NULL
  res$beta.p.value.tail1 <- NULL
  res$gamma.p.value.tail1 <- NULL
  res$beta.p.value.tail2 <- NULL
  res$gamma.p.value.tail2 <- NULL
  res$vcov <- vcov
  res$SE.diagnostics <- SE.diagnostics
  res$classification <- classifications
  res$latent.paths <- if(dependent.variable.structure == "Path"){
    grouping.initial$paths
  }else NULL
  res$diagnostics$class.shift <- class.shift
  res$npar <- sum(vapply(models, function(group){
    sum(vapply(group, function(x) length(x$parameters), numeric(1)))
  }, numeric(1)))
  res$converged <- all(vapply(models, function(group){
    all(vapply(group, `[[`, logical(1), "converged"))
  }, logical(1)))
  res$call <- call
  res$arguments$dependent.variables <- dependent.variables.list
  res$arguments$family <- family.list
  res$arguments$dependent.variable.structure <- dependent.variable.structure
  res$arguments$dependent.variable.time <- selected
  res$arguments$dependent.variable.time.cross <- dependent.variable.time.cross
  res$arguments$CEP.time.cross <- CEP.time.cross
  res$arguments$CEP.error <- CEP.error
  res$arguments$method.regression <- method.regression
  res$arguments$method.SE <- method.SE
  res$arguments$nrep.bootstrap <- nrep.bootstrap
  res
}
