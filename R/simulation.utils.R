.simulation.probability <- function(probability, L, name){
  probability <- as.numeric(probability)
  if(length(probability) != L || any(!is.finite(probability)) ||
     any(probability < 0) || sum(probability) <= 0){
    stop(name, " must contain L finite non-negative values with a positive sum")
  }
  probability / sum(probability)
}

.simulation.require.sorted <- function(probability, is.sort, name){
  if(is.sort && any(diff(probability) > sqrt(.Machine$double.eps))){
    position <- order(probability, decreasing = TRUE)
    stop(
      name, " must already follow decreasing class proportions when is.sort=TRUE; ",
      "supplied or implied proportions are c(",
      paste(sprintf("%.6f", probability), collapse = ", "),
      ") and their decreasing order is c(", paste(position, collapse = ", "),
      "). Supply parameters in that final order or use is.sort=FALSE"
    )
  }
  invisible(NULL)
}

.simulation.sample.classes <- function(probability, N){
  sample.int(length(probability), N, replace = TRUE, prob = probability)
}

.simulation.sample.classes.from.pool <- function(probability, N){
  sizes.class <- as.integer(round(probability * N))
  difference <- N - sum(sizes.class)
  if(difference != 0L){
    largest.class <- which.max(sizes.class)
    sizes.class[largest.class] <- sizes.class[largest.class] + difference
  }
  class.pool <- rep.int(seq_along(probability), sizes.class)
  sample(class.pool, N, replace = TRUE)
}

.simulation.relabel.coefficients <- function(coefficients, position, ref.class){
  coefficients <- coefficients[, position, drop = FALSE]
  sweep(coefficients, 1L, coefficients[, ref.class], "-")
}
