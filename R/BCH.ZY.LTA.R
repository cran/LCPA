#' Longitudinal BCH distal dependent-variable model
#'
#' @details
#' A separate BCH correction is constructed at each time point unless
#' `CEP.time.cross = TRUE`, in which case a posterior-mass-weighted common CEP
#' is used. State-specific external dependent-variable parameters are estimated at the
#' selected `dependent.variable.time` values. For
#' `dependent.variable.structure = "Path"`, the product of the time-specific
#' BCH weights defines one weight for every complete latent-state path. With
#' `dependent.variable.time.cross = TRUE`, records are
#' stacked and common parameters are estimated with participant-cluster
#' sandwich covariance. Participant bootstrap resamples the same individuals at
#' every time and recomputes all CEP matrices and BCH weights.
#'
#' @noRd
BCH.ZY.LTA <- function(responses, dependent.variables, L = 2,
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
  method.model <- match.arg(method.model, c("EM", "NNE", "Mplus", "Rmixmod"))
  if(is.null(params)) par.ini <- .normalize.par.ini(par.ini, method.model)
  type.model <- match.arg(type.model, c("LCA", "LPA"))
  dependent.variable.structure <- match.arg(dependent.variable.structure)
  method.regression <- match.arg(method.regression, c("Analytic", "Numeric"))
  method.SE <- match.arg(method.SE, c("Analytic", "Numeric", "Bootstrap"))
  step1.pool <- isTRUE(step1.pool)
  if(method.SE == "Bootstrap" && nrep.bootstrap < 2L){
    stop("nrep.bootstrap must be at least 2 when method.SE = 'Bootstrap'")
  }

  default_control.EM <- list(maxiter = 2000, tol = 1e-4)
  default_control.Mplus <- list(maxiter = 2000, tol = 1e-4, files.path = NULL, files.clean = TRUE)
  default_control.Rmixmod <- .default.Rmixmod.control(maxiter = 1000L)
  default_control.NNE <- list(
    hidden.layers = c(16, 16),
    activation.function = "tanh",
    use.attention = TRUE,
    d.model = 8,
    nhead = 2,
    dim.feedforward = 16,
    eps = 1e-8,
    lambda = 1e-5,
    initial.temperature = 1000,
    cooling.rate = 0.5,
    maxiter.sa = 1000,
    threshold.sa = 1e-10,
    maxiter = 1000,
    patience.early = 100,
    maxcycle = 20,
    lr = 0.025,
    scheduler.patience = 10,
    scheduler.factor = 0.80,
    plot.interval = 200,
    device = "CPU"
  )

  merge_and_clean_control <- function(user_control, default_control) {
    if (is.null(user_control)) return(default_control)
    merged <- modifyList(default_control, user_control)
    merged[names(default_control)]
  }

  control.EM <- merge_and_clean_control(control.EM, default_control.EM)
  control.Mplus <- merge_and_clean_control(control.Mplus, default_control.Mplus)
  control.NNE <- merge_and_clean_control(control.NNE, default_control.NNE)
  control.Rmixmod <- merge_and_clean_control(control.Rmixmod, default_control.Rmixmod)

  if(!is.list(responses) || length(responses) < 1L){
    stop("responses must be a non-empty list of response matrices")
  }
  responses <- lapply(responses, as.matrix)
  times <- length(responses)
  N <- nrow(responses[[1]])
  I <- ncol(responses[[1]])
  if(any(vapply(responses, nrow, integer(1)) != N) ||
     any(vapply(responses, ncol, integer(1)) != I)){
    stop("All response matrices must have the same numbers of individuals and indicators")
  }

  dependent.variable.specification <- .ZY.LTA.dependent.variables(
    dependent.variables, family, dependent.variable.time, times, N
  )
  dependent.variables.list <- dependent.variable.specification$dependent.variables
  family.list <- dependent.variable.specification$family
  selected <- dependent.variable.specification$selected

  response <- if(step1.pool) do.call(rbind, responses) else responses[[1]]
  .three.step.progress.start(
    1L, "LTA", type.model = type.model, vis = vis
  )
  if(is.null(params)){
    LCPA.obj <- .with.estimation.output.prefix({
      if(type.model == "LCA"){
        LCA(response, L = L, par.ini = par.ini,
            method = method.model, is.sort = FALSE, nrep = nrep,
            starts = starts, maxiter.warmup = maxiter.warmup,
            vis = vis,
            control.EM = control.EM,
            control.Mplus = control.Mplus,
            control.NNE = control.NNE,
            control.Rmixmod = control.Rmixmod)
      }else{
        LPA(response, L = L, par.ini = par.ini,
            constraint = constraint,
            method = method.model, is.sort = FALSE, nrep = nrep,
            starts = starts, maxiter.warmup = maxiter.warmup,
            vis = vis,
            control.EM = control.EM,
            control.Mplus = control.Mplus,
            control.NNE = control.NNE,
            control.Rmixmod = control.Rmixmod)
      }
    })
    params <- LCPA.obj$params
  }else if(vis){
    cat("  Supplied Step 1 parameters were used.\n")
  }
  if(type.model == "LCA" && is.null(params$category.levels)){
    stop("LCA Step 1 parameters must include category.levels")
  }
  if(is.sort){
    posterior.t1 <- .three.step.posterior(responses[[1]], params, type.model)
    P.Z.t1 <- colSums(posterior.t1) / sum(posterior.t1)
    position <- order(P.Z.t1, decreasing = TRUE)
    params <- .three.step.reorder.params(params, position, type.model)
  }

  .three.step.progress.start(
    2L, "LTA", path = "Z -> Y", vis = vis
  )
  P.Z.Xns <- vector("list", times)
  if(step1.pool){
    posterior.pooled <- .three.step.posterior(
      do.call(rbind, responses), params, type.model
    )
    for(t in seq_len(times)){
      index.t <- ((t - 1L) * N + 1L):(t * N)
      P.Z.Xns[[t]] <- posterior.pooled[index.t, , drop = FALSE]
    }
  }else{
    for(t in seq_len(times)){
      P.Z.Xns[[t]] <- .three.step.posterior(responses[[t]], params, type.model)
    }
  }
  for(t in seq_len(times)){
    if(is.null(colnames(P.Z.Xns[[t]]))){
      colnames(P.Z.Xns[[t]]) <- .latent.group.names(L, type.model)
    }
  }
  P.Zs <- lapply(P.Z.Xns, function(x) colSums(x) / sum(x))
  Zs <- lapply(P.Z.Xns, max.col, ties.method = "first")
  CEP <- if(CEP.error){
    get.CEP(P.Z.Xns, CEP.time.cross = CEP.time.cross)
  }else{
    replicate(times, diag(L), simplify = FALSE)
  }
  bch <- lapply(seq_len(times), function(t){
    .ZY.weights(P.Z.Xns[[t]], CEP[[t]], CEP.error)
  })
  weights.time <- lapply(bch, `[[`, "case.weights")
  latent.paths <- NULL
  path.regularization <- 0
  if(dependent.variable.structure == "Path"){
    path.weight <- .ZY.path.matrix(weights.time)
    path.posterior <- .ZY.path.matrix(P.Z.Xns)
    stabilized.path <- .ZY.stabilize.weight.mass(
      path.weight$value, path.posterior$value
    )
    weights.path <- stabilized.path$weight
    colnames(weights.path) <- colnames(path.weight$value)
    case.weights <- replicate(times, weights.path, simplify = FALSE)
    latent.paths <- path.weight$paths
    path.regularization <- stabilized.path$regularization
  }else{
    case.weights <- weights.time
  }
  if(vis){
    cat("  Posterior probabilities and modal classifications were computed.\n")
    cat(sprintf(
      "  CEP matrices and BCH weights were computed for %d time points.\n",
      times
    ))
    if(dependent.variable.structure == "Path"){
      cat("  Complete latent-state-path weights were prepared.\n")
    }else{
      cat(sprintf(
        "  State-specific dependent-variable weights were prepared for %d selected time point%s.\n",
        length(selected), if(length(selected) == 1L) "" else "s"
      ))
    }
  }

  fit.dependent.variables <- function(dependent.variables.cur, family.cur, weights.cur,
                           selected.cur, levels.cur = NULL) {
    if(dependent.variable.time.cross){
      names.cur <- lapply(dependent.variables.cur[selected.cur], names)
      if(length(unique(vapply(names.cur, paste, character(1), collapse = "\r"))) != 1L){
        stop("dependent.variable.time.cross = TRUE requires identical dependent-variable columns across selected time points")
      }
      family.key <- vapply(family.cur[selected.cur], paste, character(1), collapse = "\r")
      if(length(unique(family.key)) != 1L){
        stop("dependent.variable.time.cross = TRUE requires identical dependent-variable families across selected time points")
      }
      dependent.variables.pooled <- do.call(rbind, dependent.variables.cur[selected.cur])
      weights.pooled <- do.call(rbind, weights.cur[selected.cur])
      cluster <- rep(seq_len(nrow(dependent.variables.cur[[selected.cur[1L]]])), length(selected.cur))
      list(pooled = .ZY.fit.dependent.variables(
        dependent.variables.pooled, family.cur[[selected.cur[1L]]], weights.pooled,
        cluster = cluster,
        levels.dependent.variables = if(is.null(levels.cur)) NULL else levels.cur[[1L]],
        method.regression = method.regression,
        method.SE = if(method.SE == "Bootstrap") "Analytic" else method.SE
      ))
    }else{
      res <- lapply(seq_along(selected.cur), function(j){
        t <- selected.cur[j]
        .ZY.fit.dependent.variables(
          dependent.variables.cur[[t]], family.cur[[t]], weights.cur[[t]],
          levels.dependent.variables = if(is.null(levels.cur)) NULL else levels.cur[[j]],
          method.regression = method.regression,
          method.SE = if(method.SE == "Bootstrap") "Analytic" else method.SE
        )
      })
      names(res) <- paste0("t", selected.cur)
      res
    }
  }

  .three.step.progress.start(
    3L, "LTA", path = "Z -> Y", vis = vis
  )
  models <- fit.dependent.variables(
    dependent.variables.list, family.list, case.weights, selected
  )
  .ZY.progress.models(models, "BCH", method.SE, vis)
  levels.models <- lapply(models, function(x) lapply(x, `[[`, "levels"))
  vcov <- lapply(models, function(x) NULL)
  SE.diagnostics <- list(method = method.SE)

  if(method.SE == "Bootstrap"){
    if(vis) cat("  Bootstrapping for Standard Errors ...\n")
    estimates.bootstrap <- lapply(models, function(x){
      matrix(NA_real_, nrep.bootstrap, length(.ZY.vector(x)))
    })
    progress.state <- .new.progress.state()
    successful.bootstrap <- 0L
    for(bs in seq_len(nrep.bootstrap)){
      samples.cur <- .three.step.bootstrap.indices(N)
      P.Z.Xns.cur <- lapply(P.Z.Xns, function(x){
        x[samples.cur, , drop = FALSE]
      })
      CEP.cur <- if(CEP.error){
        get.CEP(P.Z.Xns.cur, CEP.time.cross = CEP.time.cross)
      }else{
        replicate(times, diag(L), simplify = FALSE)
      }
      weights.time.cur <- lapply(seq_len(times), function(t){
        .ZY.weights(P.Z.Xns.cur[[t]], CEP.cur[[t]], CEP.error)$case.weights
      })
      weights.cur <- if(dependent.variable.structure == "Path"){
        path.weight.cur <- .ZY.path.matrix(weights.time.cur)$value
        path.posterior.cur <- .ZY.path.matrix(P.Z.Xns.cur)$value
        stabilized.cur <- .ZY.stabilize.weight.mass(
          path.weight.cur, path.posterior.cur
        )$weight
        colnames(stabilized.cur) <- colnames(path.weight.cur)
        replicate(times, stabilized.cur, simplify = FALSE)
      }else{
        weights.time.cur
      }
      dependent.variables.cur <- dependent.variables.list
      for(t in selected){
        dependent.variables.cur[[t]] <- dependent.variables.list[[t]][samples.cur, , drop = FALSE]
      }
      models.cur <- tryCatch(
        fit.dependent.variables(
          dependent.variables.cur, family.list, weights.cur, selected, levels.models
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
    SE.diagnostics <- list(method = "Bootstrap", groups = diagnostics.bootstrap)
    .ZY.progress.bootstrap.complete(SE.diagnostics, vis)
  }

  group.count <- ncol(case.weights[[selected[1L]]])
  npar <- sum(vapply(models, function(group){
    sum(vapply(group, function(x){
      if(x$family == "gaussian"){
        2L * group.count
      }else{
        group.count * (ncol(x$estimate) - 1L)
      }
    }, numeric(1)))
  }, numeric(1)))
  res <- list(
    beta = NULL,
    gamma = NULL,
    dependent.variables = models,
    vcov = vcov,
    SE.diagnostics = SE.diagnostics,
    npar = npar,
    P.Z.Xns = P.Z.Xns,
    P.Zs = P.Zs,
    Zs = Zs,
    CEP = CEP,
    case.weights = case.weights,
    latent.paths = latent.paths,
    diagnostics = list(
      CEP.inverse.method = vapply(bch, `[[`, character(1), "inverse.method"),
      CEP.rank = vapply(bch, `[[`, integer(1), "rank"),
      CEP.condition = vapply(bch, `[[`, numeric(1), "condition"),
      CEP.reciprocal.condition = vapply(
        bch, `[[`, numeric(1), "reciprocal.condition"
      ),
      weight.regularization = vapply(
        bch, `[[`, numeric(1), "weight.regularization"
      ),
      negative.weight.proportion = vapply(
        case.weights, function(x) mean(x < 0), numeric(1)
      ),
      path.weight.regularization = path.regularization,
      class.mass.reconstruction.error = vapply(
        bch, `[[`, numeric(1), "reconstruction.error"
      )
    ),
    params = params,
    step1.pool = step1.pool,
    converged = TRUE,
    call = call,
    arguments = list(
      responses = responses, dependent.variables = dependent.variables.list,
      L = L, type.model = type.model, family = family.list,
      dependent.variable.structure = dependent.variable.structure,
      dependent.variable.time = selected,
      dependent.variable.time.cross = dependent.variable.time.cross,
      CEP.time.cross = CEP.time.cross,
      CEP.error = CEP.error,
      par.ini = par.ini, params = params,
      step1.pool = step1.pool, is.sort = is.sort,
      constraint = constraint,
      method.model = method.model,
      method.regression = method.regression, method.SE = method.SE,
      nrep.bootstrap = nrep.bootstrap,
      starts = starts, maxiter.warmup = maxiter.warmup, nrep = nrep,
      vis = vis,
      control.EM = control.EM,
      control.Mplus = control.Mplus,
      control.NNE = control.NNE,
      control.Rmixmod = control.Rmixmod
    )
  )
  class(res) <- "LTA"
  res
}
