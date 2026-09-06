.default.RMixtComp.control <- function(){
  list(
    maxiter.burnin = 50L,
    maxiter = 50L,
    maxiter.gibbs.burnin = 50L,
    maxiter.gibbs = 50L,
    n.init.per.class = 50L,
    maxattempts.sem = 20L,
    confidence.level = 0.95,
    stable.ratio = 0.99,
    n.stable = 20L,
    criterion = "BIC",
    nrep = 1L,
    ncores = 1L
  )
}

check.RMixtComp.control <- function(control.RMixtComp){
  integer.names <- c(
    "maxiter.burnin", "maxiter", "maxiter.gibbs.burnin", "maxiter.gibbs",
    "n.init.per.class", "maxattempts.sem", "n.stable", "nrep", "ncores"
  )
  minimum <- setNames(rep(1L, length(integer.names)), integer.names)
  minimum["maxiter.burnin"] <- 0L
  invalid.integer <- vapply(integer.names, function(name){
    value <- control.RMixtComp[[name]]
    length(value) != 1L || !is.finite(value) || value < minimum[name] ||
      value > .Machine$integer.max || value != floor(value)
  }, logical(1))
  if(any(invalid.integer)){
    stop(
      "RMixtComp controls must be integers with maxiter.burnin >= 0 and all other iteration/count controls >= 1: ",
      paste(integer.names[invalid.integer], collapse = ", "), "."
    )
  }
  if(control.RMixtComp$ncores > control.RMixtComp$nrep){
    stop("control.RMixtComp$ncores cannot exceed control.RMixtComp$nrep.")
  }
  if(length(control.RMixtComp$confidence.level) != 1L ||
     !is.finite(control.RMixtComp$confidence.level) ||
     control.RMixtComp$confidence.level <= 0 || control.RMixtComp$confidence.level >= 1){
    stop("control.RMixtComp$confidence.level must be strictly between 0 and 1.")
  }
  if(length(control.RMixtComp$stable.ratio) != 1L ||
     !is.finite(control.RMixtComp$stable.ratio) ||
     control.RMixtComp$stable.ratio <= 0 || control.RMixtComp$stable.ratio > 1){
    stop("control.RMixtComp$stable.ratio must be in (0, 1].")
  }
  if(length(control.RMixtComp$criterion) != 1L ||
     !control.RMixtComp$criterion %in% c("BIC", "ICL")){
    stop("control.RMixtComp$criterion must be 'BIC' or 'ICL'.")
  }
  invisible(control.RMixtComp)
}

.RMixtComp.algorithm.values <- function(control.RMixtComp){
  list(
    nbBurnInIter = control.RMixtComp$maxiter.burnin,
    nbIter = control.RMixtComp$maxiter,
    nbGibbsBurnInIter = control.RMixtComp$maxiter.gibbs.burnin,
    nbGibbsIter = control.RMixtComp$maxiter.gibbs,
    nInitPerClass = control.RMixtComp$n.init.per.class,
    nSemTry = control.RMixtComp$maxattempts.sem,
    confidenceLevel = control.RMixtComp$confidence.level,
    ratioStableCriterion = control.RMixtComp$stable.ratio,
    nStableCriterion = control.RMixtComp$n.stable
  )
}

.RMixtComp.fit <- function(data.RMixtComp, model, L, control.RMixtComp){
  check.RMixtComp.control(control.RMixtComp)

  seed <- sample.int(.Machine$integer.max, 1L)
  seed.old <- Sys.getenv("MC_DETERMINISTIC", unset = NA_character_)
  on.exit({
    if(is.na(seed.old)){
      Sys.unsetenv("MC_DETERMINISTIC")
    }else{
      Sys.setenv(MC_DETERMINISTIC = seed.old)
    }
  }, add = TRUE)
  Sys.setenv(MC_DETERMINISTIC = as.character(as.integer(seed)))

  algorithm.RMixtComp <- do.call(
    RMixtCompUtilities::createAlgo,
    .RMixtComp.algorithm.values(control.RMixtComp)
  )
  fit <- RMixtComp::mixtCompLearn(
    data = data.RMixtComp,
    model = model,
    algo = algorithm.RMixtComp,
    nClass = as.integer(L),
    criterion = control.RMixtComp$criterion,
    hierarchicalMode = "no",
    nRun = as.integer(control.RMixtComp$nrep),
    nCore = as.integer(control.RMixtComp$ncores),
    verbose = FALSE
  )

  Log.Lik <- suppressWarnings(as.numeric(fit$mixture$lnObservedLikelihood))
  if(length(Log.Lik) != 1L || !is.finite(Log.Lik)){
    error <- if(!is.null(fit$warnLog)) fit$warnLog else
      "RMixtComp did not return a finite observed-data log-likelihood."
    stop(error, call. = FALSE)
  }
  fit
}
