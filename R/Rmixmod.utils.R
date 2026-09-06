.default.Rmixmod.control <- function(maxiter = NULL){
  list(
    path = "LCPA",
    algorithm = NULL,
    nrep = NULL,
    method.init = NULL,
    starts = NULL,
    maxiter.init = NULL,
    tol.init = NULL,
    maxiter = maxiter,
    tol = NULL,
    par.ini = NULL,
    labels.ini = NULL,
    strategy = NULL
  )
}

.Rmixmod.strategy <- function(control.Rmixmod){
  if(!is.null(control.Rmixmod$strategy)){
    return(control.Rmixmod$strategy)
  }
  strategy.args <- list(
    algo = control.Rmixmod$algorithm,
    nbTry = control.Rmixmod$nrep,
    initMethod = control.Rmixmod$method.init,
    nbTryInInit = control.Rmixmod$starts,
    nbIterationInInit = control.Rmixmod$maxiter.init,
    nbIterationInAlgo = control.Rmixmod$maxiter,
    epsilonInInit = control.Rmixmod$tol.init,
    epsilonInAlgo = control.Rmixmod$tol,
    parameter = control.Rmixmod$par.ini,
    labels = control.Rmixmod$labels.ini
  )
  strategy.args <- strategy.args[!vapply(strategy.args, is.null, logical(1))]
  do.call(Rmixmod::mixmodStrategy, strategy.args)
}
