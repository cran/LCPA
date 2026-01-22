#' @importFrom MplusAutomation mplusObject mplusModeler
#' @importFrom dplyr filter
#' @importFrom stats na.omit
Mplus.LCA <- function(response, L = 2,
                      nrep = 10, starts = 200, maxiter.wa=20,
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
    cat("Temporary working: ", paste0(getwd(), "/", temp_dir), "\n")
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
            cat("Successed to clean up: ", paste0(getwd(), "/", temp_dir), "\n")
          }
        }
      }
    }, add = TRUE)
  }

  var_names <- colnames(response)
  if (is.null(var_names)) {
    var_names <- paste0("V", seq_len(ncol(response)))
  }
  colnames(response) <- var_names
  var_names <- format_mplus_vars_auto(var_names)
  df <- as.data.frame(response)

  poly.value <- sapply(df, function(x) length(unique(na.omit(x))))
  if (any(poly.value == 1)) stop("Some variables have only 1 level; invalid for LCA.")

  variable_str <- paste0("CLASSES = c1(", L, ");\n",
                         "CATEGORICAL = ", paste(var_names, collapse = " "), ";\n",
                         "ANALYSIS:\n",
                         "  TYPE = mixture;\n",
                         "  STARTS = ", starts, " ", nrep, ";\n",
                         "  STITERATIONS = ", maxiter.wa, ";\n",
                         "  MITERATIONS = ", maxiter, ";\n",
                         "  CONVERGENCE = ", tol, ";")

  model_str <- "%OVERALL%"

  output_str <- "  TECH11 TECH14;"

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
    cat("Running Mplus ...\n")
  }

  Mplus.obj <- suppressWarnings(
    MplusAutomation::mplusModeler(
      mobj,
      dataout   = dataout_path,
      modelout  = modelout_path,
      run       = TRUE,
      writeData = "always",
      check     = FALSE,
      quiet     = TRUE
    )
  )

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
    grepl("^V[0-9]+$", .data[["param"]]),
    .data[["LatentClass"]] %in% 1:L
  )

  par <- array(
    NA_real_,
    dim = c(L, I, poly.max),
    dimnames = list(
      paste0("Class.", 1:L),
      paste0("V", 1:I),
      paste0("Cat.", 1:poly.max)
    )
  )
  for (i in 1:nrow(prob_data)) {
    row <- prob_data[i, ]
    cls_idx <- as.integer(row$LatentClass)
    var_idx <- as.integer(gsub("V", "", row$param))
    cat_idx <- as.integer(row$category)
    if (cls_idx < 1 || cls_idx > L) next
    if (var_idx < 1 || var_idx > I) next
    if (cat_idx < 1 || cat_idx > poly.max) next
    est_val <- suppressWarnings(as.numeric(row$est))
    if (is.na(est_val) || est_val < 0 || est_val > 1) {
      warning(paste("Invalid probability at row", i,
                    ": Class", cls_idx, "Var", var_idx, "Cat", cat_idx,
                    "Value =", row$est))
      next
    }
    par[cls_idx, var_idx, cat_idx] <- est_val
  }

  P.Z.Xn <- matrix(1/L, N, L)
  L.Xi.Z <- matrix(1, L, N)
  for(p in 1:N){
    for(i in 1:I){
      L.Xi.Z[ , p] <- L.Xi.Z[ , p] * par[ , i, response[p, i]+1]
    }
    L.Xi <- sum(L.Xi.Z[, p] * P.Z)
    P.Z.Xn[p, ] <- (L.Xi.Z[, p] * P.Z) / (L.Xi + 1e-300)
  }

  Log.Lik <- get.Log.Lik.LCA(response, P.Z, par)
  npar <- get.npar.LCA(poly.value, L)
  AIC <- -2 * Log.Lik + 2 * npar
  BIC <- -2 * Log.Lik + npar * log(N)

  if (vis) {
    cat(sprintf("Mplus Model: %s\nLog-likelihood = %.5f | BIC = %.2f\n",
                title_str, Log.Lik, BIC))
  }

  res = list(
    params = list(par = par, P.Z = P.Z),
    npar = npar,
    Log.Lik = Log.Lik,
    AIC=AIC,
    BIC=BIC,
    P.Z.Xn = P.Z.Xn,
    P.Z = P.Z,
    Z = apply(P.Z.Xn, 1, which.max),
    probability = NULL
  )

  return(res)
}

