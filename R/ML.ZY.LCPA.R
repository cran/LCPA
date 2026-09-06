#' Classification-error-corrected ML distal dependent-variable model
#'
#' @details
#' For a Gaussian dependent variable, the observed likelihood contribution is
#' \deqn{\sum_{l=1}^L\mathrm{CEP}(l,\widehat{Z}_n)\pi_l
#' \mathcal{N}(Y_{nv}\mid\mu_{lv},\sigma_{lv}^2).}
#' For a categorical dependent variable, replace the normal density by
#' \eqn{\prod_q p_{lvq}^{\mathbb{1}(Y_{nv}=q)}}.
#' Gaussian dependent variables estimate class-specific means and unequal variances;
#' categorical dependent variables estimate class-specific multinomial probabilities.
#' Analytic gradients and Louis observed information are implemented directly;
#' numeric alternatives use finite-difference optimization/Hessians. Because Y
#' updates third-step class responsibilities, the proportion of changed modal
#' assignments is returned as a class-shift diagnostic. Bootstrap recomputes CEP
#' for every resample while holding Step 1 measurement parameters fixed.
#'
#' @noRd
ML.ZY.LCPA <- function(response, dependent.variables, L = 2,
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
  method.regression <- match.arg(method.regression, c("Analytic", "Numeric"))
  method.SE <- match.arg(method.SE, c("Analytic", "Numeric", "Bootstrap"))
  if(method.SE == "Bootstrap" && nrep.bootstrap < 2L){
    stop("nrep.bootstrap must be at least 2 when method.SE = 'Bootstrap'")
  }

  res <- ML.XZ.LCPA(
    response = response, L = L,
    ref.class = L, type.model = type.model,
    covariates = NULL, CEP.error = CEP.error,
    par.ini = par.ini, params = params, is.sort = is.sort,
    constraint = constraint, method.model = method.model,
    method.regression = "Analytic", method.SE = "Analytic",
    starts = starts, maxiter.warmup = maxiter.warmup, nrep = nrep,
    vis = vis,
    control.EM = control.EM, control.Mplus = control.Mplus,
    control.NNE = control.NNE, control.flexmix = control.flexmix,
    control.Rmixmod = control.Rmixmod,
    control.RMixtComp = control.RMixtComp,
    .progress.path = "Z -> Y", .progress.step3 = FALSE
  )

  posterior <- res$P.Z.Xn
  modal <- res$Z
  dependent.variables <- .ZY.dependent.variables(dependent.variables, nrow(posterior))
  family <- .ZY.family(family, dependent.variables)
  res$P.Z.Xns <- list(posterior)
  res$P.Zs <- list(res$P.Z)
  res$Zs <- list(modal)
  res$CEP <- if(CEP.error){
    get.CEP(res$P.Z.Xns, CEP.time.cross = FALSE)
  }else list(diag(L))
  classification <- .ML.ZY.classification(modal, res$CEP[[1L]])
  colnames(classification) <- colnames(posterior)
  if(vis){
    cat("  Classification-error-corrected ML weights were prepared.\n")
  }
  levels.dependent.variables <- lapply(seq_along(dependent.variables), function(j){
    if(family[j] == "categorical") levels(factor(dependent.variables[[j]])) else NULL
  })
  method.SE.fit <- if(method.SE == "Bootstrap") "Analytic" else method.SE
  .three.step.progress.start(
    3L, "LCPA", path = "Z -> Y", vis = vis
  )
  models <- .ML.ZY.fit.dependent.variables(
    dependent.variables, family, classification, res$P.Z,
    method.regression, method.SE.fit, levels.dependent.variables, modal
  )
  .ZY.progress.models(models, "ML", method.SE, vis)
  vcov <- NULL
  SE.diagnostics <- list(
    method = method.SE,
    class.shift = vapply(models, `[[`, numeric(1), "class.shift")
  )

  if(method.SE == "Bootstrap"){
    if(vis) cat("  Bootstrapping for Standard Errors ...\n")
    estimates.bootstrap <- matrix(
      NA_real_, nrep.bootstrap, length(.ZY.vector(models))
    )
    N <- nrow(posterior)
    progress.state <- .new.progress.state()
    successful.bootstrap <- 0L
    for(bs in seq_len(nrep.bootstrap)){
      samples.cur <- .three.step.bootstrap.indices(N)
      posterior.cur <- posterior[samples.cur, , drop = FALSE]
      dependent.variables.cur <- dependent.variables[samples.cur, , drop = FALSE]
      modal.cur <- max.col(posterior.cur, ties.method = "first")
      CEP.cur <- if(CEP.error){
        get.CEP(list(posterior.cur), CEP.time.cross = FALSE)[[1L]]
      }else{
        diag(L)
      }
      classification.cur <- .ML.ZY.classification(modal.cur, CEP.cur)
      colnames(classification.cur) <- colnames(posterior.cur)
      models.cur <- tryCatch(
        .ML.ZY.fit.dependent.variables(
          dependent.variables.cur, family, classification.cur,
          colSums(posterior.cur) / sum(posterior.cur),
          method.regression, "Analytic", levels.dependent.variables, modal.cur
        ),
        error = function(e) NULL
      )
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
    SE.diagnostics$class.shift <- vapply(
      models, `[[`, numeric(1), "class.shift"
    )
    .ZY.progress.bootstrap.complete(SE.diagnostics, vis)
  }

  class.shift <- vapply(models, `[[`, numeric(1), "class.shift")
  if(any(class.shift > 0.20)){
    warning("The ML distal model changed more than 20% of modal assignments; inspect class-shift diagnostics and consider method.3step = 'BCH'")
  }
  res$dependent.variables <- list(t1 = models)
  res$beta <- NULL
  res$gamma <- NULL
  res$beta.se <- NULL
  res$beta.Z.sta <- NULL
  res$beta.p.value.tail1 <- NULL
  res$beta.p.value.tail2 <- NULL
  res$vcov <- list(t1 = vcov)
  res$SE.diagnostics <- SE.diagnostics
  res$classification <- list(classification)
  res$diagnostics$class.shift <- class.shift
  res$npar <- sum(vapply(models, function(x) length(x$parameters), numeric(1)))
  res$converged <- all(vapply(models, `[[`, logical(1), "converged"))
  res$call <- call
  res$arguments$dependent.variables <- dependent.variables
  res$arguments$family <- family
  res$arguments$CEP.error <- CEP.error
  res$arguments$method.regression <- method.regression
  res$arguments$method.SE <- method.SE
  res$arguments$nrep.bootstrap <- nrep.bootstrap
  res
}
