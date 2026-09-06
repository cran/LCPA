#' BCH distal dependent-variable model for latent classes or profiles
#'
#' @details
#' The dependent variable is external to the measurement model. For
#' \eqn{\mathrm{CEP}(l,k)=P(\widehat{Z}_n=k\mid Z_n=l)}, continuous
#' class means solve
#' \deqn{\sum_{n=1}^N
#' (\mathrm{CEP}^{-1})_{\widehat{Z}_n,l}
#' (Y_{nv}-\mu_{lv})=0,}
#' class variances solve the corresponding weighted centered-square equation,
#' and categorical probabilities replace `Y` by category indicators. Analytic
#' and numeric joint sandwich covariance is available. Bootstrap resamples observations,
#' recomputes CEP/BCH weights, and therefore propagates Step 2 sampling error.
#' Singular CEP matrices use a generalized inverse. If inverse weighting gives
#' inadmissible categorical probabilities, they are projected to the probability
#' simplex and the fitted variable records `regularized = TRUE`.
#'
#' @noRd
BCH.ZY.LCPA <- function(response, dependent.variables, L = 2,
                    type.model = "LCA",
                    family = "gaussian",
                    CEP.error = TRUE,
                    par.ini = "random",
                    params = NULL, is.sort = TRUE,
                    constraint = "VV",
                    method.model = "EM",
                    method.regression = "Analytic",
                    method.SE = "Bootstrap", nrep.bootstrap = 100,
                    starts = 100, maxiter.warmup = 20, nrep = 20,
                    vis = TRUE,
                    control.EM = NULL,
                    control.Mplus = NULL,
                    control.NNE = NULL,
                    control.flexmix = NULL,
                    control.Rmixmod = NULL,
                    control.RMixtComp = NULL) {

  call <- match.call()
  method.model <- match.arg(method.model, c("EM", "NNE", "Mplus", "flexmix", "Rmixmod", "RMixtComp"))
  if(is.null(params)) par.ini <- .normalize.par.ini(par.ini, method.model)
  type.model <- match.arg(type.model, c("LCA", "LPA"))
  method.regression <- match.arg(method.regression, c("Analytic", "Numeric"))
  method.SE <- match.arg(method.SE, c("Analytic", "Numeric", "Bootstrap"))
  if(method.SE == "Bootstrap" && nrep.bootstrap < 2L){
    stop("nrep.bootstrap must be at least 2 when method.SE = 'Bootstrap'")
  }

  default_control.EM <- list(maxiter = 2000, tol = 1e-4)
  default_control.Mplus <- list(maxiter = 2000, tol = 1e-4, files.path = NULL, files.clean = TRUE)
  default_control.Rmixmod <- .default.Rmixmod.control(maxiter = 1000L)
  default_control.RMixtComp <- .default.RMixtComp.control()
  default_control.flexmix <- .default.flexmix.control()
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
  control.flexmix <- merge_and_clean_control(control.flexmix, default_control.flexmix)
  control.Rmixmod <- merge_and_clean_control(control.Rmixmod, default_control.Rmixmod)
  control.RMixtComp <- merge_and_clean_control(control.RMixtComp, default_control.RMixtComp)

  response <- as.matrix(response)
  N <- nrow(response)
  dependent.variables <- .ZY.dependent.variables(dependent.variables, N)
  family <- .ZY.family(family, dependent.variables)

  .three.step.progress.start(
    1L, "LCPA", type.model = type.model, vis = vis
  )
  params.supplied <- !is.null(params)
  if(is.null(params)){
    LCPA.obj <- .with.estimation.output.prefix({
      if(type.model == "LCA"){
        LCA(response, L = L, par.ini = par.ini,
            method = method.model, is.sort = is.sort, nrep = nrep,
            starts = starts, maxiter.warmup = maxiter.warmup,
            vis = vis,
            control.EM = control.EM,
            control.Mplus = control.Mplus,
            control.NNE = control.NNE,
            control.flexmix = control.flexmix,
            control.Rmixmod = control.Rmixmod,
            control.RMixtComp = control.RMixtComp)
      }else{
        LPA(response, L = L, par.ini = par.ini,
            constraint = constraint,
            method = method.model, is.sort = is.sort, nrep = nrep,
            starts = starts, maxiter.warmup = maxiter.warmup,
            vis = vis,
            control.EM = control.EM,
            control.Mplus = control.Mplus,
            control.NNE = control.NNE,
            control.flexmix = control.flexmix,
            control.Rmixmod = control.Rmixmod,
            control.RMixtComp = control.RMixtComp)
      }
    })
    params <- LCPA.obj$params
  }else if(vis){
    cat("  Supplied Step 1 parameters were used.\n")
  }
  if(type.model == "LCA" && is.null(params$category.levels)){
    stop("LCA Step 1 parameters must include category.levels")
  }
  if(params.supplied && is.sort){
    position <- order(params$P.Z, decreasing = TRUE)
    params <- .three.step.reorder.params(params, position, type.model)
  }

  .three.step.progress.start(
    2L, "LCPA", path = "Z -> Y", vis = vis
  )
  posterior <- .three.step.posterior(response, params, type.model)
  if(is.null(colnames(posterior))){
    colnames(posterior) <- .latent.group.names(ncol(posterior), type.model)
  }
  P.Z <- colSums(posterior) / sum(posterior)
  Z <- max.col(posterior, ties.method = "first")
  bch <- .ZY.weights(posterior, CEP.error = CEP.error)
  if(vis){
    cat("  Posterior probabilities and modal classifications were computed.\n")
    cat("  The CEP matrix and BCH weights were computed.\n")
  }

  .three.step.progress.start(
    3L, "LCPA", path = "Z -> Y", vis = vis
  )
  levels.dependent.variables <- lapply(seq_along(dependent.variables), function(j){
    if(family[j] == "categorical") levels(factor(dependent.variables[[j]])) else NULL
  })
  models <- .ZY.fit.dependent.variables(
    dependent.variables, family, bch$case.weights,
    levels.dependent.variables = levels.dependent.variables,
    method.regression = method.regression,
    method.SE = if(method.SE == "Bootstrap") "Analytic" else method.SE
  )
  .ZY.progress.models(models, "BCH", method.SE, vis)
  vcov <- NULL
  SE.diagnostics <- list(method = method.SE)

  if(method.SE == "Bootstrap"){
    if(vis) cat("  Bootstrapping for Standard Errors ...\n")
    estimates.bootstrap <- matrix(
      NA_real_, nrep.bootstrap, length(.ZY.vector(models))
    )
    progress.state <- .new.progress.state()
    successful.bootstrap <- 0L
    for(bs in seq_len(nrep.bootstrap)){
      samples.cur <- .three.step.bootstrap.indices(N)
      posterior.cur <- posterior[samples.cur, , drop = FALSE]
      dependent.variables.cur <- dependent.variables[samples.cur, , drop = FALSE]
      models.cur <- tryCatch({
        bch.cur <- .ZY.weights(posterior.cur, CEP.error = CEP.error)
        .ZY.fit.dependent.variables(
          dependent.variables.cur, family, bch.cur$case.weights,
          levels.dependent.variables = levels.dependent.variables,
          method.regression = method.regression,
          method.SE = "Analytic"
        )
      }, error = function(e) NULL)
      if(!is.null(models.cur)){
        estimates.bootstrap[bs, ] <- .ZY.vector(models.cur)
        successful.bootstrap <- successful.bootstrap +
          as.integer(all(is.finite(estimates.bootstrap[bs, ])))
      }
      .three.step.progress.bootstrap(
        bs, nrep.bootstrap, successful.bootstrap, progress.state, vis
      )
    }
    if(vis) .end.estimation.progress()
    bootstrap <- .ZY.bootstrap(models, estimates.bootstrap)
    models <- bootstrap$models
    vcov <- bootstrap$vcov
    SE.diagnostics <- bootstrap$diagnostics
    .ZY.progress.bootstrap.complete(SE.diagnostics, vis)
  }

  npar <- sum(vapply(seq_along(models), function(j){
    if(family[j] == "gaussian") 2L * L else L * (ncol(models[[j]]$estimate) - 1L)
  }, numeric(1)))
  res <- list(
    beta = NULL,
    gamma = NULL,
    dependent.variables = list(t1 = models),
    vcov = list(t1 = vcov),
    SE.diagnostics = SE.diagnostics,
    npar = npar,
    P.Z.Xn = posterior,
    P.Z = P.Z,
    Z = Z,
    P.Z.Xns = list(posterior),
    P.Zs = list(P.Z),
    Zs = list(Z),
    CEP = list(bch$CEP),
    case.weights = list(bch$case.weights),
    diagnostics = list(
      CEP.inverse.method = bch$inverse.method,
      CEP.rank = bch$rank,
      CEP.condition = bch$condition,
      CEP.reciprocal.condition = bch$reciprocal.condition,
      weight.regularization = bch$weight.regularization,
      negative.weight.proportion = mean(bch$case.weights < 0),
      class.mass.reconstruction.error = bch$reconstruction.error
    ),
    params = params,
    converged = TRUE,
    call = call,
    arguments = list(
      response = response, dependent.variables = dependent.variables,
      L = L, type.model = type.model, family = family,
      CEP.error = CEP.error,
      par.ini = par.ini, params = params,
      is.sort = is.sort, constraint = constraint,
      method.model = method.model,
      method.regression = method.regression, method.SE = method.SE,
      nrep.bootstrap = nrep.bootstrap,
      starts = starts, maxiter.warmup = maxiter.warmup, nrep = nrep,
      vis = vis,
      control.EM = control.EM,
      control.Mplus = control.Mplus,
      control.NNE = control.NNE,
      control.flexmix = control.flexmix,
      control.Rmixmod = control.Rmixmod,
      control.RMixtComp = control.RMixtComp
    )
  )
  class(res) <- "LCPA"
  res
}
