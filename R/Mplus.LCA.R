#' @importFrom MplusAutomation mplusObject mplusModeler
#' @importFrom dplyr filter
#' @importFrom stats na.omit

Mplus.LCA <- function(response, L = 2,
                      nrep = 10, starts = 200, maxiter.warmup=20,
                      vis = TRUE,
                      maxiter = 2000, tol = 1e-4,
                      files.path = NULL,
                      files.clean = TRUE) {

  if (!is.matrix(response) && !is.data.frame(response)) {
    stop("response must be a matrix or data frame")
  }
  response <- as.matrix(response)
  if (!is.numeric(response)) stop("response must be numeric")
  if (ncol(response) < 1) stop("At least one indicator variable required")
  L <- as.integer(L)

  if (is.null(files.path)) {
    stop("No valid 'files.path' provided!", call. = FALSE)
  }
  if(files.path != ""){
    if (!dir.exists(files.path)) {
      dir.create(
        files.path,
        recursive    = TRUE,
        showWarnings = FALSE
      )
      if (!dir.exists(files.path)) {
        stop("Failed to create: ", paste0(getwd(), "/", files.path), call. = FALSE)
      }
    }
    timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
    temp_dir  <- file.path(files.path, paste0("Mplus_LCA_", timestamp))
  }else{
    timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
    temp_dir <- paste0("Mplus_LCA_", timestamp)
  }
  if (!dir.exists(temp_dir)) {
    dir.create(
      temp_dir,
      recursive    = TRUE,
      showWarnings = FALSE
    )
    if (!dir.exists(temp_dir)) {
      stop("Failed to create: ", paste0(getwd(), "/", temp_dir), call. = FALSE)
    }
  }
  if(vis){
    cat(.estimation.output.prefix(), "Temporary working: ", paste0(getwd(), "/", temp_dir), "\n", sep = "")
  }
  if (isTRUE(files.clean)) {
    on.exit({
      if (dir.exists(temp_dir)) {
        for (i in 1:5) {
          unlink(temp_dir, recursive = TRUE, force = TRUE)
          if (!dir.exists(temp_dir)) break
          Sys.sleep(0.2)
        }
        if (dir.exists(temp_dir)) {
          warning(
            "Failed to clean up: ", paste0(getwd(), "/", temp_dir),
            "\nPlease remove it manually if safe to do so.",
            call. = FALSE
          )
        }else{
          if(vis){
            cat(.estimation.output.prefix(), "Successed to clean up: ", paste0(getwd(), "/", temp_dir), "\n", sep = "")
          }
        }
      }
    }, add = TRUE)
  }

  variable.names <- colnames(response)
  if (is.null(variable.names)) {
    variable.names <- paste0("V", seq_len(ncol(response)))
  }

  standardize_varnames <- function(names) {
    names <- gsub("^[^a-zA-Z]", "V", names)
    names <- gsub("[^a-zA-Z0-9_]", "_", names)
    names <- make.unique(names, sep = "_")
    return(names)
  }

  variable.names.standardized <- standardize_varnames(variable.names)
  colnames(response) <- variable.names.standardized
  variable.names <- variable.names.standardized

  df <- as.data.frame(response)

  poly.value <- sapply(df, function(x) length(unique(na.omit(x))))
  if (any(poly.value == 1)) stop("Some variables have only 1 level; invalid for LCA.")

  variable.names.formatted <- format_mplus_vars_auto(variable.names)
  seed <- sample.int(.Machine$integer.max, 1L)

  variable_str <- paste0("CLASSES = c1(", L, ");\n",
                         "CATEGORICAL = ", variable.names.formatted, ";\n",
                         "ANALYSIS:\n",
                         "  TYPE = mixture;\n",
                         "  STARTS = ", starts, " ", nrep, ";\n",
                         "  STSEED = ", seed, ";\n",
                         "  STITERATIONS = ", maxiter.warmup, ";\n",
                         "  MITERATIONS = ", maxiter, ";\n",
                         "  CONVERGENCE = ", tol, ";")

  model_str <- "%OVERALL%"

  output_str <- "  TECH8;"

  post_file <- file.path(temp_dir, "posterior.dat")
  savedata_str <- paste0('  FILE = "', post_file, '";\n',
                         "  SAVE = CPROBABILITIES;")

  title_str <- sprintf("LCA with %d classes", L)

  mobj <- MplusAutomation::mplusObject(
    TITLE    = title_str,
    VARIABLE = variable_str,
    MODEL    = model_str,
    OUTPUT   = output_str,
    SAVEDATA = savedata_str,
    rdata    = df
  )

  modelout_path <- file.path(temp_dir, "lca_model.inp")
  dataout_path <- file.path(temp_dir, "lca_data.dat")

  if (vis) {
    cat(.estimation.output.prefix(), "Running Mplus ...\n", sep = "")
  }

  Mplus.obj <- suppressMessages(suppressWarnings(
    MplusAutomation::mplusModeler(
      mobj,
      dataout   = dataout_path,
      modelout  = modelout_path,
      run       = TRUE,
      writeData = "always",
      check     = FALSE,
      quiet     = TRUE
    )
  ))

  if(is.null(Mplus.obj) ||
     is.null(Mplus.obj$results) ||
     is.null(Mplus.obj$results$parameters)){
    stop("Mplus reported an error in parameter estimation, please switch to using method = 'EM' or method = 'NNE'")
  }

  N <- nrow(response)
  I <- ncol(response)
  P.Z <- Mplus.obj$results$class_counts$modelEstimated$proportion
  params_df <- Mplus.obj$results$parameters$probability.scale
  poly.max <- max(response) + 1
  poly.value <- apply(response, 2, function(x){length(unique(x))})

  prob_data <- dplyr::filter(
    params_df,
    .data[["param"]] %in% variable.names,
    .data[["LatentClass"]] %in% 1:L
  )

  par <- array(
    NA_real_,
    dim = c(L, I, poly.max),
    dimnames = list(
      .latent.group.names(L, "LCA"),
      variable.names,
      paste0("Cat.", 1:poly.max)
    )
  )

  for (i in 1:nrow(prob_data)) {
    row <- prob_data[i, ]
    cls_idx <- as.integer(row$LatentClass)
    variable.index <- match(row$param, variable.names)
    cat_idx <- as.integer(row$category)

    if (is.na(variable.index)) next

    if (cls_idx < 1 || cls_idx > L) next
    if (variable.index < 1 || variable.index > I) next
    if (cat_idx < 1 || cat_idx > poly.max) next

    est_val <- suppressWarnings(as.numeric(row$est))
    if (is.na(est_val) || est_val < 0 || est_val > 1) {
      warning(paste("Invalid probability at row", i,
                    ": Class", cls_idx, "Var", variable.index, "Cat", cat_idx,
                    "Value =", row$est))
      next
    }
    par[cls_idx, variable.index, cat_idx] <- est_val
  }

  category.levels <- lapply(seq_len(I), function(i) seq_len(poly.value[i]) - 1L)
  P.Z.Xn <- get.P.Z.Xn.LCA(response, par, P.Z, category.levels)

  Log.Lik <- get.Log.Lik.LCA(response, par, P.Z)
  npar <- get.npar.LCA(poly.value, L)
  AIC <- -2 * Log.Lik + 2 * npar
  BIC <- -2 * Log.Lik + npar * log(N)

  if (vis) {
    cat(sprintf("%sMplus Model: %s\n%sLog-likelihood = %.5f | BIC = %.2f\n",
                .estimation.output.prefix(), title_str,
                .estimation.output.prefix(), Log.Lik, BIC))
  }

  res = list(
    params = list(par = par, P.Z = P.Z),
    npar = npar,
    Log.Lik = Log.Lik,
    AIC=AIC,
    BIC=BIC,
    P.Z.Xn = P.Z.Xn,
    P.Z = P.Z,
    Z = max.col(P.Z.Xn, ties.method = "first"),
    probability = NULL
  )

  return(res)
}

format_mplus_vars_auto <- function(variable.names, indent = "  ", max.line.length = 70) {
  result_lines <- c()
  current_line <- indent

  for (var in variable.names) {
    sep <- if (nchar(current_line) == nchar(indent)) "" else " "
    candidate <- paste0(current_line, sep, var)

    if (nchar(candidate) > max.line.length) {
      result_lines <- c(result_lines, current_line)
      current_line <- paste0(indent, var)
    } else {
      current_line <- candidate
    }
  }

  result_lines <- c(result_lines, current_line)
  return(paste(result_lines, collapse = "\n"))
}
