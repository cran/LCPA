#' @importFrom methods is new slot
#' @importFrom reticulate source_python
#' @importFrom scales breaks_width
#' @importFrom stats dnorm optim quantile
#' @importFrom utils globalVariables packageDescription
#' @import Rcpp
#' @useDynLib LCPA, .registration = TRUE
#'
NULL

utils::globalVariables(c("NN_LCA", "NN_LPA"))

.validate.training.stages <- function(starts, maxiter.warmup, nrep){
  values <- list(
    starts = starts,
    maxiter.warmup = maxiter.warmup,
    nrep = nrep
  )
  invalid <- vapply(values, function(value){
    length(value) != 1L || !is.finite(value) || value < 1 ||
      value > .Machine$integer.max || value != floor(value)
  }, logical(1))
  if(any(invalid)){
    stop(
      paste(names(values)[invalid], collapse = ", "),
      " must be positive integers."
    )
  }
  if(starts < nrep){
    stop("starts must be greater than or equal to nrep.")
  }
  invisible(NULL)
}

.normalize.par.ini <- function(par.ini, method){
  if(is.character(par.ini) && length(par.ini) == 1L &&
     identical(par.ini, "kmeans") && !method %in% c("EM", "NNE")){
    warning(
      "par.ini = 'kmeans' is unavailable for method = '", method,
      "'; using par.ini = 'random'.",
      call. = FALSE
    )
    return("random")
  }
  par.ini
}

.validate.LPA.constraint <- function(constraint, I){
  named.constraints <- c("UE", "UV", "E0", "V0", "EE", "VV", "VE", "EV")
  if(is.character(constraint)){
    if(length(constraint) != 1L || !constraint %in% named.constraints){
      stop("unsupported covariance constraint")
    }
    if(I > 1L && constraint %in% c("UE", "UV")){
      stop("UE and UV are only available for univariate response")
    }
    return(constraint)
  }
  if(!is.list(constraint)){
    stop("constraint must be a character string or a list")
  }
  valid <- vapply(constraint, function(indices){
    length(indices) == 2L && is.numeric(indices) && all(is.finite(indices)) &&
      all(indices == as.integer(indices)) && all(indices >= 1L) && all(indices <= I)
  }, logical(1))
  if(!all(valid)){
    stop("each custom constraint must contain two valid variable indices")
  }
  constraint.keys <- vapply(constraint, function(indices){
    paste(sort(as.integer(indices)), collapse = ":")
  }, character(1))
  constraint[!duplicated(constraint.keys)]
}

.class.indicator <- function(classes, L) {
  indicator <- matrix(0, length(classes), L)
  indicator[cbind(seq_along(classes), classes)] <- 1
  indicator
}

.class.proportions <- function(classes, L) {
  tabulate(classes, nbins = L) / length(classes)
}

.latent.group.label <- function(type.model, plural = FALSE) {
  type.model <- match.arg(type.model, c("LCA", "LPA"))
  label <- if(type.model == "LCA") "Class" else "Profile"
  if(plural) paste0(label, "s") else label
}

.latent.group.names <- function(L, type.model) {
  paste(.latent.group.label(type.model), seq_len(L))
}

.latent.group.columns <- function(x, type.model) {
  label <- .latent.group.label(type.model)
  replacement <- c(
    Class = label,
    From.Class = paste0("From.", label),
    To.Class = paste0("To.", label)
  )
  position <- match(names(replacement), names(x), nomatch = 0L)
  names(x)[position[position > 0L]] <- replacement[position > 0L]
  x
}

.latent.path.names <- function(paths) {
  paste("Path", apply(paths, 1L, paste, collapse = "-"))
}

.estimation.output.prefix <- function(){
  getOption("LCPA.estimation.output.prefix", "")
}

.with.estimation.output.prefix <- function(expr, prefix = "  "){
  old.options <- options(
    LCPA.estimation.output.prefix = paste0(.estimation.output.prefix(), prefix)
  )
  on.exit(options(old.options), add = TRUE)
  force(expr)
}

.new.progress.state <- function(){
  state <- new.env(parent = emptyenv())
  state$width <- 0L
  state
}

.print.progress.line <- function(output, progress.state){
  output.width <- nchar(output, type = "width")
  cat(
    "\r", output,
    strrep(" ", max(0L, progress.state$width - output.width)), sep = ""
  )
  progress.state$width <- output.width
  invisible(NULL)
}

.print.estimation.progress <- function(stage, current, total, Log.Lik,
                                       best.Log.Lik = NULL,
                                       algorithm = NULL, iterations = NULL,
                                       progress.state){
  stage <- match.arg(stage, c("Warm", "Rep"))
  output <- sprintf(
    "%s%s %d/%d", .estimation.output.prefix(), stage, current, total
  )
  if(!is.null(algorithm) && !is.null(iterations)){
    output <- paste0(output, " | ", algorithm, " iterations = ", iterations)
  }
  output <- paste0(output, " | Log-likelihood = ", sprintf("%.5f", Log.Lik))
  if(!is.null(best.Log.Lik)){
    output <- paste0(output, " | Best = ", sprintf("%.5f", best.Log.Lik))
  }
  .print.progress.line(output, progress.state)
}

.print.iteration.progress <- function(iteration, maxchg, BIC,
                                      progress.prefix = "",
                                      progress.state){
  output <- paste0(
    progress.prefix, "Iter = ", iteration,
    " | \u0394Log.Lik = ", sprintf("%.5f", maxchg),
    " | BIC = ", sprintf("%.2f", BIC)
  )
  .print.progress.line(output, progress.state)
  invisible(NULL)
}

.end.estimation.progress <- function(){
  cat("\n")
  invisible(NULL)
}

.print.estimation.summary <- function(algorithm, Log.Lik, BIC){
  cat(sprintf("%s%s: Log-likelihood = %.5f | BIC = %.2f\n",
              .estimation.output.prefix(), algorithm, Log.Lik, BIC))
  invisible(NULL)
}
