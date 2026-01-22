
#' @importFrom clue solve_LSAP
get.index.LCA <- function(res.t, res.e){

  N <- length(res.t$Z)
  L <- dim(res.t$par)[1]
  I <- dim(res.t$par)[2]
  poly.max <- dim(res.t$par)[3]

  if(!is.null(res.e$par)){
    par.t <- par.e <- NULL
    for(po in 1:poly.max){
      par.t <- cbind(par.t, res.t$par[ , , po])
      par.e <- cbind(par.e, res.e$par[ , , po])
    }
    posi <- apply(par.t, 2, function(x) !any(is.na(x)))
    par.t <- t(par.t)[posi, ]
    par.e[which(is.infinite(par.e))] <- NA
    par.e <- t(par.e)[posi, ]
    posi <- apply(par.e, 1, function(x) !any(is.na(x)))
    par.e <- par.e[posi, ]
    par.t <- par.t[posi, ]

    # dist.mat <- distance.matrix(par.t, par.e)
    # assignment <- clue::solve_LSAP(dist.mat)
    MAP.e <- matrix(apply(res.e$P.Z.Xn, 1, max), N, L, byrow=FALSE) == res.e$P.Z.Xn
    dist.mat <- distance.matrix(res.t$P.Z.Xn, MAP.e)
    assignment <- clue::solve_LSAP(dist.mat)

    par.e.aligned  <- par.e[, as.numeric(assignment)]

    par.dif <- par.t - par.e.aligned
    par.MSE <- mean(par.dif^2, na.rm = TRUE)
  }else{
    par.MSE <- NA
  }

  if(!is.null(res.e$P.Z.Xn)){
    Z.t <- res.t$Z
    P.Z.Xn <- res.e$P.Z.Xn[, as.numeric(assignment)]
    Z.e <- apply(P.Z.Xn, 1, which.max)
    Z.dif <- Z.t != Z.e
    acc <- 1 - mean(Z.dif)
    res <- c(par.MSE, acc)
  }else{
    res <- c(par.MSE, NA)
  }

  names(res) <- c("MSE", "acc")

  return(res)
}

#' @importFrom clue solve_LSAP
 get.index.LPA <- function(res.t, res.e){

  N <- nrow(res.t$response)
  I <- ncol(res.t$means)
  L <- nrow(res.t$means)

  means.t <- t(res.t$means)
  means.e <- t(res.e$params$means)
  covs.t <- res.t$covs
  covs.e <- res.e$params$covs

  MAP.e <- matrix(apply(res.e$P.Z.Xn, 1, max), N, L, byrow=FALSE) == res.e$P.Z.Xn

  covs.t[which(covs.t == 0)] <- NA

  if(!is.null(means.e)){
    dist.mat <- distance.matrix(res.t$P.Z.Xn, MAP.e)
    assignment <- clue::solve_LSAP(dist.mat)
    means.e.aligned  <- means.e[, as.numeric(assignment)]

    means.dif <- means.t - means.e.aligned
    MSE.m <- mean(means.dif^2, na.rm = TRUE)

    covs.e.aligned <- covs.e[, , as.numeric(assignment)]

    covs.dif <- covs.t - covs.e.aligned
    MSE.c <- mean(covs.dif^2, na.rm = TRUE)

  }else{
    MSE.m <- NA
    MSE.c <- NA
  }

  if(!is.null(res.e$P.Z.Xn) & !is.null(means.e)){
    Z.t <- res.t$Z
    P.Z.Xn <- res.e$P.Z.Xn[, as.numeric(assignment)]
    Z.e <- apply(P.Z.Xn, 1, which.max)
    Z.dif <- Z.t != Z.e
    acc <- 1 - mean(Z.dif)
    res <- c(MSE.m, MSE.c, acc)
  }else{
    res <- c(MSE.m, MSE.c, NA)
  }

  names(res) <- c("MSE.m", "MSE.c", "acc")

  return(res)
}

 distance.matrix <- function(A, B) {
  L <- ncol(A)
  dist.mat <- matrix(0, L, L)
  for (i in 1:L) {
    for (j in 1:L) {
      dist.mat[i, j] <- sum((A[, i] - B[, j])^2)
    }
  }
  return(dist.mat)
}


get.runs <- function(posi, trial.list){
  number.trial <- length(trial.list)
  trials.length <- length(trial.list[[number.trial]])
  for(i in (number.trial-1):1)
    trials.length <- c(trials.length[1] * length(trial.list[[i]]), trials.length)
  for(i in 1:number.trial)
    trials.length[i] <- trials.length[i] / length(trial.list[[i]])

  runs <- rep(1, number.trial)
  for(i in 1:number.trial){
    if(posi > trials.length[i]){
      runs[i] <- ceiling(posi / trials.length[i])
      posi <- posi - (runs[i] - 1) * trials.length[i]
      if(posi == 0 & i < number.trial)
        for(j in (i+1):number.trial)
          runs[j] <- length(trial.list[[j]])
    }
  }
  return(runs)
}

format_mplus_vars_auto <- function(var_names, indent = "  ", max_line_length = 80) {
  current_line <- indent
  lines <- character(0)

  for (v in var_names) {
    proposed <- paste(current_line, v)
    if (nchar(proposed) > max_line_length) {
      lines <- c(lines, current_line)
      current_line <- paste(indent, v)
    } else {
      current_line <- proposed
    }
  }
  lines <- c(lines, current_line)  # add last line

  return(paste(lines, collapse = "\n"))
}

is.valid.params <- function(params, model.type) {
  if (any(params$P.Z <= 0) || any(params$P.Z >= 1) || abs(sum(params$P.Z) - 1) > 1e-6) {
    return(FALSE)
  }
  if (model.type == "LPA") {
    L <- length(params$P.Z)
    I <- ncol(params$means)
    if (!is.numeric(params$means) || any(is.na(params$means))) {
      return(FALSE)
    }
    for (l in 1:L) {
      cov_mat <- params$covs[, , l]
      if (!isSymmetric(cov_mat, tol = 1e-8) || !is.positive.definite(cov_mat)) {
        return(FALSE)
      }
    }
  } else if (model.type == "LCA") {
    L <- length(params$P.Z)
    I <- dim(params$par)[2]
    poly.max <- dim(params$par)[3]

    for (l in 1:L) {
      for (i in 1:I) {
        probs <- params$par[l, i, ]
        probs <- probs[!is.na(probs)]  # Handle unused categories
        if (length(probs) > 0 && (any(probs <= 0) || any(probs >= 1) || abs(sum(probs) - 1) > 1e-6)) {
          return(FALSE)
        }
      }
    }
  }

  return(TRUE)
}

params.to.row <- function(params, model.type, poly.value = NULL) {
  if (model.type == "LPA") {
    means <- params$means
    covs <- params$covs
    pz <- params$P.Z

    L <- nrow(means)
    I <- ncol(means)
    means.vec <- as.vector(t(means))
    covs.vec <- c()
    for (class.idx in 1:L) {
      mat <- covs[, , class.idx]
      covs.vec <- c(covs.vec, mat[upper.tri(mat, diag = TRUE)])
    }
    pz.vec <- if (L > 1) as.vector(pz[-L]) else numeric(0)

    return(c(means.vec, covs.vec, pz.vec))

  } else if (model.type == "LCA") {
    if (is.null(poly.value)) {
      stop("poly.value must be provided for LCA model")
    }
    P.Z <- params$P.Z
    par.array <- params$par

    L_class <- length(P.Z)
    I_var <- dim(par.array)[2]
    stopifnot(length(poly.value) == I_var)
    pz.vec <- if (L_class > 1) P.Z[-L_class] else numeric(0)
    par.vec <- c()
    for (l in 1:L_class) {
      for (i in 1:I_var) {
        K_i <- poly.value[i]
        if (K_i > 1) {
          probs.vec <- par.array[l, i, 1:(K_i - 1)]
          par.vec <- c(par.vec, probs.vec)
        }
      }
    }

    return(c(pz.vec, par.vec))
  } else {
    stop("Unsupported model type: ", model.type)
  }
}

row.to.params <- function(vec, model.type, I, L, poly.value = NULL, validate = FALSE) {
  # +++ ONLY CHANGE: Handle complex numbers safely +++
  vec <- Re(vec)

  if (model.type == "LPA") {
    num_means <- L * I
    num_cov_elements <- I * (I + 1) / 2
    num_covs <- L * num_cov_elements
    num_pz <- max(L - 1, 0)  # Handle L=1 case
    expected.length <- num_means + num_covs + num_pz

    if (length(vec) != expected.length) {
      stop(sprintf("Parameter vector length mismatch. Expected %d, got %d",
                   expected.length, length(vec)))
    }

    idx <- 1
    means.vec <- vec[idx:(idx + num_means - 1)]
    idx <- idx + num_means
    means.mat <- matrix(means.vec, nrow = L, ncol = I, byrow = TRUE)

    covs.vec <- vec[idx:(idx + num_covs - 1)]
    idx <- idx + num_covs
    covs.arr <- array(0, dim = c(I, I, L))
    pos <- 1
    for (class.idx in 1:L) {
      mat <- matrix(0, nrow = I, ncol = I)
      upper_idx <- upper.tri(mat, diag = TRUE)
      num_elements <- sum(upper_idx)

      if (pos + num_elements - 1 <= length(covs.vec)) {
        mat[upper_idx] <- covs.vec[pos:(pos + num_elements - 1)]
        mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
        pos <- pos + num_elements
      } else {
        stop(sprintf("Insufficient values for class %d covariance matrix", class.idx))
      }
      if (!is.positive.definite(mat)) {
        if (validate) {
          return(NULL)  # Trigger invalid parameter handling
        }
        warning(sprintf("Covariance matrix for class %d is not positive definite. Correcting...", class.idx))
        suppressWarnings(pd <- Matrix::nearPD(mat, corr = FALSE, maxit = 1000))
        mat <- as.matrix(pd$mat)
      }

      covs.arr[, , class.idx] <- mat
    }
    if (num_pz > 0) {
      pz.vec <- vec[idx:(idx + num_pz - 1)]
      pz.vec <- c(pz.vec, 1 - sum(pz.vec))
    } else {
      pz.vec <- 1  # Only one class
    }
    if (any(pz.vec <= 0) || abs(sum(pz.vec) - 1) > 1e-6) {
      if (validate) {
        return(NULL)  # Trigger invalid parameter handling
      }
      pz.vec[pz.vec <= 0] <- .Machine$double.eps
      pz.vec <- pz.vec / sum(pz.vec)
    }
    pz.named <- setNames(pz.vec, paste0("Class.", 1:L))

    return(list(
      means = means.mat,
      covs = covs.arr,
      P.Z = pz.named
    ))

  } else if (model.type == "LCA") {
    if (is.null(poly.value)) {
      stop("poly.value must be provided for LCA model")
    }
    stopifnot(length(poly.value) == I)
    num_pz <- max(L - 1, 0)
    if (num_pz > 0) {
      pz.vec <- vec[1:num_pz]
      par.vec <- vec[-(1:num_pz)]
    } else {
      pz.vec <- numeric(0)
      par.vec <- vec
    }

    if (length(pz.vec) > 0) {
      pz.vec <- c(pz.vec, 1 - sum(pz.vec))
      if (any(pz.vec <= 0) || abs(sum(pz.vec) - 1) > 1e-6) {
        if (validate) {
          return(NULL)  # Trigger invalid parameter handling
        }
        pz.vec[pz.vec <= 0] <- .Machine$double.eps
        pz.vec <- pz.vec / sum(pz.vec)
      }
    } else {
      pz.vec <- 1  # Only one class
    }
    pz.named <- setNames(pz.vec, paste0("Class.", 1:L))
    poly.max <- max(poly.value)
    par.array <- array(NA, dim = c(L, I, poly.max))

    pos <- 1
    for (l in 1:L) {
      for (i in 1:I) {
        K_i <- poly.value[i]
        if (K_i > 1) {
          num.vec <- K_i - 1
          if (pos + num.vec - 1 > length(par.vec)) {
            stop(sprintf("Parameter vector too short for item %d in class %d", i, l))
          }
          probs.vec <- par.vec[pos:(pos + num.vec - 1)]
          pos <- pos + num.vec
          probs.full <- c(probs.vec, 1 - sum(probs.vec))
          if (any(probs.full <= 0) || abs(sum(probs.full) - 1) > 1e-6) {
            if (validate) {
              return(NULL)  # Trigger invalid parameter handling
            }
            probs.full[probs.full <= 0] <- .Machine$double.eps
            probs.full <- probs.full / sum(probs.full)
          }

          par.array[l, i, 1:K_i] <- probs.full
        } else if (K_i == 1) {
          par.array[l, i, 1] <- 1
        }
      }
    }

    return(list(
      par = par.array,
      P.Z = pz.named
    ))
  } else {
    stop("Unsupported model type: ", model.type)
  }
}

is.positive.definite <- function(mat) {
  if (!is.matrix(mat) || nrow(mat) != ncol(mat)) return(FALSE)
  tryCatch({
    chol(mat)
    TRUE
  }, error = function(e) FALSE)
}

