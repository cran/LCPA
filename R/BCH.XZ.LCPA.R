#' BCH three-step model for covariates predicting latent membership
#'
#' @details
#' For the orientation
#' \eqn{\mathrm{CEP}(l,k)=P(\widehat{Z}_n=k\mid Z_n=l)}, the coefficient for
#' class \eqn{l\ne l_0} solves
#' \deqn{\mathbf{U}_l(\boldsymbol{\beta}_l)=
#' \sum_{n=1}^N\boldsymbol{\zeta}_n
#' \left\{(\mathrm{CEP}^{-1})_{\widehat{Z}_n,l}
#' -P(Z_n=l\mid\boldsymbol{\zeta}_n)
#' \sum_{h=1}^L(\mathrm{CEP}^{-1})_{\widehat{Z}_n,h}\right\}=0.}
#' The analytic SE uses a robust sandwich with multinomial bread and empirical
#' score meat. The numeric SE replaces the bread by the numerical Hessian.
#' Bootstrap resamples individuals and recomputes the CEP inverse before every
#' refit. Negative BCH weights are retained rather than truncated.
#'
#' @noRd
BCH.XZ.LCPA <- function(response, L = 2,
                        ref.class = L, type.model = "LCA",
                        covariates = NULL,
                        CEP.error = TRUE,
                        par.ini = "random",
                        params = NULL, is.sort = TRUE,
                        constraint = "VV",
                        method.model = "EM", tol = 1e-4,
                        method.regression = "Analytic",
                        lower = -10, upper = 10,
                        method.SE = "Bootstrap", nrep.bootstrap = 100,
                        maxiter = 5000, starts = 100,
                        maxiter.warmup = 20, nrep = 20,
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
    response = response, L = L, ref.class = ref.class,
    type.model = type.model, covariates = covariates,
    CEP.error = CEP.error,
    par.ini = par.ini, params = params, is.sort = is.sort,
    constraint = constraint, method.model = method.model,
    tol = tol, method.regression = "Analytic",
    lower = lower, upper = upper,
    method.SE = "Analytic", maxiter = maxiter,
    starts = starts, maxiter.warmup = maxiter.warmup, nrep = nrep,
    vis = vis,
    control.EM = control.EM, control.Mplus = control.Mplus,
    control.NNE = control.NNE, control.flexmix = control.flexmix,
    control.Rmixmod = control.Rmixmod,
    control.RMixtComp = control.RMixtComp,
    .progress.path = "X -> Z", .progress.step3 = FALSE
  )
  covariates <- res$arguments$covariates
  posterior <- res$P.Z.Xn
  CEP <- if(CEP.error){
    get.CEP(list(posterior), CEP.time.cross = FALSE)
  }else list(diag(L))
  bch <- .ZY.weights(posterior, CEP[[1L]], CEP.error)
  if(vis) cat("  BCH weights were computed.\n")
  .three.step.progress.start(
    3L, "LCPA", path = "X -> Z", vis = vis
  )
  method.SE.fit <- if(method.SE == "Bootstrap") "Analytic" else method.SE
  fit <- .BCH.XZ.multinomial(
    covariates, bch$case.weights, ref.class,
    method.regression, method.SE.fit,
    lower, upper, tol, maxiter
  )
  .three.step.progress.regression(
    fit, "BCH", method.SE, "latent-class regression models", vis
  )
  vcov <- fit$vcov
  SE.diagnostics <- list(method = method.SE)

  if(method.SE == "Bootstrap"){
    if(vis) cat("  Bootstrapping for Standard Errors ...\n")
    params.bootstrap <- matrix(NA_real_, nrep.bootstrap, length(fit$params))
    N <- nrow(posterior)
    progress.state <- .new.progress.state()
    successful.bootstrap <- 0L
    for(bs in seq_len(nrep.bootstrap)){
      samples.cur <- .three.step.bootstrap.indices(N)
      posterior.cur <- posterior[samples.cur, , drop = FALSE]
      CEP.cur <- if(CEP.error){
        get.CEP(list(posterior.cur), CEP.time.cross = FALSE)[[1L]]
      }else diag(L)
      fit.cur <- tryCatch(
        .BCH.XZ.multinomial(
          covariates[samples.cur, , drop = FALSE],
          .ZY.weights(posterior.cur, CEP.cur, CEP.error)$case.weights,
          ref.class, method.regression, "Analytic",
          lower, upper, tol, maxiter
        ),
        error = function(e) NULL
      )
      if(!is.null(fit.cur) && fit.cur$converged){
        params.bootstrap[bs, ] <- fit.cur$params
        successful.bootstrap <- successful.bootstrap +
          as.integer(all(is.finite(params.bootstrap[bs, ])))
      }
      .three.step.progress.bootstrap(
        bs, nrep.bootstrap, successful.bootstrap, progress.state, vis
      )
    }
    if(vis) .end.estimation.progress()
    successful <- complete.cases(params.bootstrap)
    vcov <- if(sum(successful) >= 2L){
      stats::cov(params.bootstrap[successful, , drop = FALSE])
    }else{
      matrix(NA_real_, length(fit$params), length(fit$params))
    }
    SE.diagnostics <- list(
      method = "Bootstrap", successful = sum(successful),
      attempted = nrep.bootstrap
    )
    if(vis){
      cat(sprintf(
        "  Bootstrap Standard Errors were computed from %d/%d successful replications.\n",
        SE.diagnostics$successful, SE.diagnostics$attempted
      ))
    }
  }

  statistics <- .BCH.XZ.statistics(
    fit$params, sqrt(pmax(diag(vcov), 0)), list(covariates),
    L, ref.class, FALSE
  )
  res$beta <- statistics$estimate$beta
  res$beta.se <- statistics$se$beta
  res$beta.Z.sta <- statistics$z$beta
  res$beta.p.value.tail1 <- statistics$tail1$beta
  res$beta.p.value.tail2 <- statistics$tail2$beta
  res$vcov <- vcov
  res$information <- fit$information
  res$SE.diagnostics <- SE.diagnostics
  res$bound.diagnostics <- fit$bound.diagnostics
  res$npar <- length(fit$params)
  res$Log.Lik <- -fit$objective
  res$AIC <- 2 * fit$objective + 2 * res$npar
  res$BIC <- 2 * fit$objective + log(nrow(posterior)) * res$npar
  res$iterations <- fit$iterations
  res$converged <- fit$converged
  res$CEP <- CEP
  res$case.weights <- list(bch$case.weights)
  res$diagnostics <- list(
    CEP.inverse.method = bch$inverse.method,
    CEP.rank = bch$rank,
    CEP.condition = bch$condition,
    CEP.reciprocal.condition = bch$reciprocal.condition,
    weight.regularization = bch$weight.regularization,
    negative.weight.proportion = mean(bch$case.weights < 0),
    class.mass.reconstruction.error = bch$reconstruction.error
  )
  res$call <- call
  res$arguments$method.regression <- method.regression
  res$arguments$method.SE <- method.SE
  res$arguments$nrep.bootstrap <- nrep.bootstrap
  res
}
