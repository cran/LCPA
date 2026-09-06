#' Longitudinal BCH model for initial-state and transition regressions
#'
#' @details
#' At time `t`, case weights use the corresponding row of the stable inverse of
#' `CEP[[t]]`; singular matrices use `MASS::ginv()`. Initial-state
#' regression uses `w[i1, c]`. If classification errors are conditionally
#' independent over time given the latent-state path, the inverse correction for
#' the joint pair `(Z[t-1], Z[t])` is the Kronecker product of the two inverse
#' CEP matrices, hence the transition pseudo-weight is
#' `w[i, t-1, from] * w[i, t, to]`. All initial and transition scores are joined
#' before forming the participant-cluster sandwich covariance, preserving
#' cross-equation and repeated-measure score covariance. With
#' `covariates.time.cross = TRUE`, transition records are pooled while sharing
#' one coefficient block for each origin state. Participant bootstrap recomputes
#' every CEP and every joint pseudo-weight.
#'
#' @noRd
BCH.XZ.LTA <- function(responses, L = 2,
                       ref.class = L, type.model = "LCA",
                       covariates = NULL,
                       CEP.time.cross = FALSE,
                       CEP.error = TRUE,
                       covariates.time.cross = FALSE,
                       par.ini = "random",
                       params = NULL, step1.pool = FALSE, is.sort = TRUE,
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
                       control.Rmixmod = NULL) {

  call <- match.call()
  method.regression <- match.arg(method.regression, c("Analytic", "Numeric"))
  method.SE <- match.arg(method.SE, c("Analytic", "Numeric", "Bootstrap"))
  if(method.SE == "Bootstrap" && nrep.bootstrap < 2L){
    stop("nrep.bootstrap must be at least 2 when method.SE = 'Bootstrap'")
  }
  res <- ML.XZ.LTA(
    responses = responses, L = L, ref.class = ref.class,
    type.model = type.model, covariates = covariates,
    CEP.time.cross = CEP.time.cross, CEP.error = CEP.error,
    covariates.time.cross = covariates.time.cross,
    par.ini = par.ini, params = params, step1.pool = step1.pool,
    is.sort = is.sort, constraint = constraint,
    method.model = method.model, tol = tol,
    method.regression = "Analytic", lower = lower, upper = upper,
    method.SE = "Analytic", maxiter = maxiter,
    starts = starts, maxiter.warmup = maxiter.warmup, nrep = nrep,
    vis = vis,
    control.EM = control.EM, control.Mplus = control.Mplus,
    control.NNE = control.NNE, control.Rmixmod = control.Rmixmod,
    .progress.path = "X -> Z", .progress.step3 = FALSE
  )
  covariates <- res$arguments$covariates
  posterior <- res$P.Z.Xns
  times <- length(posterior)
  N <- nrow(posterior[[1L]])
  CEP <- if(CEP.error){
    get.CEP(posterior, CEP.time.cross = CEP.time.cross)
  }else replicate(times, diag(L), simplify = FALSE)
  bch <- lapply(seq_len(times), function(t){
    .ZY.weights(posterior[[t]], CEP[[t]], CEP.error)
  })
  weights <- lapply(bch, `[[`, "case.weights")
  if(vis){
    cat(sprintf("  BCH weights were computed for %d time points.\n", times))
  }

  fit.all <- function(covariates.cur, weights.cur, method.SE.cur){
    fits <- list(.BCH.XZ.multinomial(
      covariates.cur[[1L]], weights.cur[[1L]], ref.class,
      method.regression, method.SE.cur,
      lower, upper, tol, maxiter,
      cluster = seq_len(nrow(covariates.cur[[1L]]))
    ))
    if(times > 1L){
      if(covariates.time.cross){
        for(from.class in seq_len(L)){
          design.cur <- do.call(rbind, covariates.cur[2:times])
          weight.cur <- do.call(rbind, lapply(2:times, function(t){
            weights.cur[[t - 1L]][, from.class] * weights.cur[[t]]
          }))
          fits[[length(fits) + 1L]] <- .BCH.XZ.multinomial(
            design.cur, weight.cur, ref.class,
            method.regression, method.SE.cur,
            lower, upper, tol, maxiter,
            cluster = rep(seq_len(nrow(covariates.cur[[1L]])), times - 1L)
          )
        }
      }else{
        for(t in 2:times){
          for(from.class in seq_len(L)){
            weight.cur <- weights.cur[[t - 1L]][, from.class] *
              weights.cur[[t]]
            fits[[length(fits) + 1L]] <- .BCH.XZ.multinomial(
              covariates.cur[[t]], weight.cur, ref.class,
              method.regression, method.SE.cur,
              lower, upper, tol, maxiter,
              cluster = seq_len(nrow(covariates.cur[[t]]))
            )
          }
        }
      }
    }
    params.cur <- unlist(lapply(fits, `[[`, "params"), use.names = FALSE)
    information.cur <- .BCH.XZ.block.diagonal(
      lapply(fits, `[[`, "information")
    )
    score.cur <- do.call(cbind, lapply(fits, `[[`, "score"))
    inverse.information <- MASS::ginv(information.cur)
    vcov.cur <- inverse.information %*% crossprod(score.cur) %*%
      inverse.information
    list(
      fits = fits, params = params.cur,
      information = information.cur, vcov = vcov.cur,
      objective = sum(vapply(fits, `[[`, numeric(1), "objective")),
      iterations = sum(vapply(fits, `[[`, numeric(1), "iterations")),
      converged = all(vapply(fits, `[[`, logical(1), "converged"))
    )
  }

  .three.step.progress.start(
    3L, "LTA", path = "X -> Z", vis = vis
  )
  method.SE.fit <- if(method.SE == "Bootstrap") "Analytic" else method.SE
  fit <- fit.all(covariates, weights, method.SE.fit)
  .three.step.progress.regression(
    fit, "BCH", method.SE, "latent-state regression models", vis
  )
  vcov <- fit$vcov
  SE.diagnostics <- list(
    method = method.SE,
    longitudinal.BCH.assumption =
      "Classification errors are conditionally independent across time given the latent-state path"
  )

  if(method.SE == "Bootstrap"){
    if(vis) cat("  Bootstrapping for Standard Errors ...\n")
    params.bootstrap <- matrix(NA_real_, nrep.bootstrap, length(fit$params))
    progress.state <- .new.progress.state()
    successful.bootstrap <- 0L
    for(bs in seq_len(nrep.bootstrap)){
      samples.cur <- .three.step.bootstrap.indices(N)
      posterior.cur <- lapply(posterior, function(x){
        x[samples.cur, , drop = FALSE]
      })
      CEP.cur <- if(CEP.error){
        get.CEP(posterior.cur, CEP.time.cross = CEP.time.cross)
      }else replicate(times, diag(L), simplify = FALSE)
      weights.cur <- lapply(seq_len(times), function(t){
        .ZY.weights(posterior.cur[[t]], CEP.cur[[t]], CEP.error)$case.weights
      })
      covariates.cur <- lapply(covariates, function(x){
        x[samples.cur, , drop = FALSE]
      })
      fit.cur <- tryCatch(
        fit.all(covariates.cur, weights.cur, "Analytic"),
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
      attempted = nrep.bootstrap,
      longitudinal.BCH.assumption =
        "Classification errors are conditionally independent across time given the latent-state path"
    )
    if(vis){
      cat(sprintf(
        "  Bootstrap Standard Errors were computed from %d/%d successful replications.\n",
        SE.diagnostics$successful, SE.diagnostics$attempted
      ))
    }
  }

  statistics <- .BCH.XZ.statistics(
    fit$params, sqrt(pmax(diag(vcov), 0)), covariates,
    L, ref.class, covariates.time.cross
  )
  res$beta <- statistics$estimate$beta
  res$gamma <- statistics$estimate$gamma
  res$beta.se <- statistics$se$beta
  res$gamma.se <- statistics$se$gamma
  res$beta.Z.sta <- statistics$z$beta
  res$gamma.Z.sta <- statistics$z$gamma
  res$beta.p.value.tail1 <- statistics$tail1$beta
  res$gamma.p.value.tail1 <- statistics$tail1$gamma
  res$beta.p.value.tail2 <- statistics$tail2$beta
  res$gamma.p.value.tail2 <- statistics$tail2$gamma
  res$vcov <- vcov
  res$information <- fit$information
  res$SE.diagnostics <- SE.diagnostics
  res$bound.diagnostics <- lapply(fit$fits, `[[`, "bound.diagnostics")
  res$npar <- length(fit$params)
  res$Log.Lik <- -fit$objective
  res$AIC <- 2 * fit$objective + 2 * res$npar
  res$BIC <- 2 * fit$objective + log(N) * res$npar
  res$iterations <- fit$iterations
  res$converged <- fit$converged
  res$CEP <- CEP
  res$case.weights <- weights
  res$diagnostics <- list(
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
      weights, function(x) mean(x < 0), numeric(1)
    ),
    class.mass.reconstruction.error = vapply(
      bch, `[[`, numeric(1), "reconstruction.error"
    ),
    longitudinal.BCH.assumption =
      "Classification errors are conditionally independent across time given the latent-state path"
  )
  res$call <- call
  res$arguments$method.regression <- method.regression
  res$arguments$method.SE <- method.SE
  res$arguments$nrep.bootstrap <- nrep.bootstrap
  res
}
