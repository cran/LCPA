#' @importFrom MplusAutomation mplusObject mplusModeler
#' @importFrom dplyr filter
#' @importFrom stats var

Mplus.LPA <- function(response, L = 2, constraint = "VV",
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
  if (ncol(response) < 2) stop("At least two indicator variables required")

  L <- as.integer(L)
  I <- ncol(response)
  N <- nrow(response)

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
    temp_dir  <- file.path(files.path, paste0("Mplus_LPA_", timestamp))
  }else{
    timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
    temp_dir <- paste0("Mplus_LPA_", timestamp)
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

  orig_variable.names <- colnames(response)
  if (is.null(orig_variable.names)) {
    orig_variable.names <- paste0("V", seq_len(I))
  }
  standardize_varnames <- function(names) {
    names_std <- tolower(names)
    names_std <- gsub("^[^a-z]", "v", names_std)
    names_std <- gsub("[^a-z0-9_]", "_", names_std)
    names_std <- make.unique(names_std, sep = "_")
    return(names_std)
  }

  variable.names.standardized <- standardize_varnames(orig_variable.names)
  colnames(response) <- variable.names.standardized
  variable.names <- variable.names.standardized

  df <- as.data.frame(response)
  seed <- sample.int(.Machine$integer.max, 1L)

  variable_str <- paste0("CLASSES = c1(", L, ");\n",
                         "ANALYSIS:\n",
                         "  TYPE = mixture;\n",
                         "  STARTS = ", starts, " ", nrep, ";\n",
                         "  STSEED = ", seed, ";\n",
                         "  STITERATIONS = ", maxiter.warmup, ";\n",
                         "  MITERATIONS = ", maxiter, ";\n",
                         "  CONVERGENCE = ", tol, ";")

  model_lines <- generate_mplus_model(constraint, variable.names, L, class.variable = "c1")
  model_str <- paste(model_lines, collapse = "\n")
  output_str <- "  TECH8;"

  post_file <- file.path(temp_dir, "posterior.dat")
  savedata_str <- paste0('  FILE = "', post_file, '";\n',
                         "  SAVE = CPROBABILITIES;")

  title_str <- sprintf("LPA with %d classes (%s constraint)", L,
                       ifelse(is.character(constraint), constraint, "custom"))

  mobj <- MplusAutomation::mplusObject(
    TITLE    = title_str,
    VARIABLE = variable_str,
    MODEL    = model_str,
    OUTPUT   = output_str,
    SAVEDATA = savedata_str,
    rdata    = df
  )

  modelout_path <- file.path(temp_dir, "lpa_model.inp")
  dataout_path <- file.path(temp_dir, "lpa_data.dat")

  if (vis) {
    cat(.estimation.output.prefix(), "Running Mplus ...\n", sep = "")
  }

  Mplus.obj <- MplusAutomation::mplusModeler(
    mobj,
    dataout  = dataout_path,
    modelout = modelout_path,
    run      = TRUE,
    writeData = "always",
    check     = FALSE,
    quiet     = TRUE
  )

  out_file <- sub("\\.inp$", ".out", modelout_path)
  for (i in 1:20) {
    if (file.exists(out_file)) break
    Sys.sleep(0.1)
  }

  if (!file.exists(out_file)) {
    stop("Mplus output file not found: ", out_file)
  }

  if(is.null(Mplus.obj) ||
     is.null(Mplus.obj$results) ||
     is.null(Mplus.obj$results$parameters)){
    stop("Mplus reported an error in parameter estimation, please switch to using method = 'EM' or method = 'NNE'")
  }

  analysis_result <- Mplus.obj$results
  params_df <- analysis_result$parameters$unstandardized

  variable.names.upper <- toupper(variable.names)

  mean_data <- dplyr::filter(params_df, .data[["paramHeader"]] == "Means" & .data[["LatentClass"]] %in% 1:L)
  var_data  <- dplyr::filter(params_df, .data[["paramHeader"]] == "Variances" & .data[["param"]] %in% variable.names.upper & .data[["LatentClass"]] %in% 1:L)
  cov_data <- dplyr::filter(
    params_df,
    grepl("\\.WITH$", .data[["paramHeader"]]) &
      .data[["LatentClass"]] %in% 1:L
  )

  if (nrow(mean_data) > 0) {
    mean_data$param <- tolower(mean_data$param)
  }
  if (nrow(var_data) > 0) {
    var_data$param <- tolower(var_data$param)
  }
  if (nrow(cov_data) > 0) {
    cov_data$param <- tolower(cov_data$param)
    cov_data$paramHeader <- tolower(cov_data$paramHeader)
  }

  means <- matrix(0, nrow = L, ncol = I, dimnames = list(1:L, variable.names))
  covs <- array(0, dim = c(I, I, L), dimnames = list(variable.names, variable.names, 1:L))

  for (i in 1:nrow(mean_data)) {
    cls <- as.integer(sub("C1#", "", mean_data$LatentClass[i]))
    var <- mean_data$param[i]
    if (var %in% colnames(means)) {
      est_val <- as.numeric(mean_data$est[i])
      if (!is.na(est_val) && is.finite(est_val)) {
        means[cls, var] <- est_val
      }
    }
  }

  for (i in 1:nrow(var_data)) {
    cls <- as.integer(sub("C1#", "", var_data$LatentClass[i]))
    var <- var_data$param[i]
    if (var %in% rownames(covs) && var %in% colnames(covs)) {
      est_val <- as.numeric(var_data$est[i])
      if (!is.na(est_val) && is.finite(est_val)) {
        covs[var, var, cls] <- est_val
      }
    }
  }

  for (i in 1:nrow(cov_data)) {
    cls_raw <- as.character(cov_data$LatentClass[i])
    cls <- as.integer(sub("C1#", "", cls_raw))

    raw_var1 <- sub("\\.with$", "", as.character(cov_data$paramHeader[i]), ignore.case = TRUE)
    raw_var1 <- tolower(raw_var1)

    raw_var2 <- tolower(as.character(cov_data$param[i]))
    variable1.index <- match(raw_var1, variable.names)
    variable2.index <- match(raw_var2, variable.names)

    if (!is.na(variable1.index) && !is.na(variable2.index)) {
      est_val <- as.numeric(cov_data$est[i])

      if (!is.na(est_val) && is.finite(est_val)) {
        covs[variable1.index, variable2.index, cls] <- est_val
        covs[variable2.index, variable1.index, cls] <- est_val
      }
    }
  }

  P.Z <- analysis_result$class_counts$posteriorProb[, 3]
  if(is.null(P.Z)){
    P.Z <- rep(1/L, L)
  }
  P.Z <- pmax(as.numeric(P.Z), 1e-12)
  P.Z <- P.Z / sum(P.Z)

  covariance.repair <- .stabilize.LPA.covariances(
    covs, constraint, fallback = stats::cov(response)
  )
  covs <- covariance.repair$covs

  P.Z.Xn <- get.P.Z.Xn.LPA(response, means, covs, P.Z)
  nk <- colSums(P.Z.Xn)
  empty_clusters <- which(nk < 1e-5)
  if (length(empty_clusters) > 0) {
    non_empty <- setdiff(1:L, empty_clusters)
    P.Z.Xn[, empty_clusters] <- 0
    row_sums <- rowSums(P.Z.Xn[, non_empty, drop = FALSE])
    row_sums[row_sums < 1e-20] <- 1
    P.Z.Xn[, non_empty] <- P.Z.Xn[, non_empty, drop = FALSE] / row_sums
    nk <- colSums(P.Z.Xn)
  }

  Z <- max.col(P.Z.Xn, ties.method = "first")
  P.Z <- pmax(nk / N, 1e-12)
  P.Z <- P.Z / sum(P.Z)

  Log.Lik <- get.Log.Lik.LPA(response, means, covs, P.Z)
  npar <- get.npar.LPA(I, L, constraint)

  AIC <- -2 * Log.Lik + 2 * npar
  BIC <- -2 * Log.Lik + npar * log(N)

  if (vis) {
    cat(sprintf("%sMplus Model: %s\n%sLog-likelihood = %.5f | BIC = %.2f\n",
                .estimation.output.prefix(), title_str,
                .estimation.output.prefix(), Log.Lik, BIC))
  }

  dimnames(means) <- list(1:L, orig_variable.names)
  dimnames(covs) <- list(orig_variable.names, orig_variable.names, 1:L)

  res <- list(
    params = list(means = means, covs = covs, P.Z = P.Z),
    model = Mplus.obj,
    Log.Lik = Log.Lik,
    npar = npar,
    AIC = AIC,
    BIC = BIC,
    P.Z.Xn = P.Z.Xn,
    P.Z = P.Z,
    Z = Z,
    covariance.repaired = covariance.repair$repaired,
    covariance.repair = covariance.repair$diagnostics
  )

  return(res)
}

#' @importFrom utils combn
#' @importFrom stats setNames
generate_mplus_model <- function(constraint, variable.names, L, class.variable = "c1") {
  I <- length(variable.names)
  variable.pairs <- t(combn(variable.names, 2))

  if (is.character(constraint) && constraint == "VE") {
    overall_lines <- c(
      "%OVERALL%",
      apply(variable.pairs, 1, function(pair) {
        pair_idx <- which(apply(variable.pairs, 1, function(x) all(sort(x) == sort(pair))))
        sprintf("%s WITH %s (c%d);", pair[1], pair[2], pair_idx)
      })
    )

    class_lines <- lapply(1:L, function(k) {
      c(
        paste0("%", class.variable, "#", k, "%"),
        paste0("[", paste(variable.names, collapse = " "), "];"),
        paste0(variable.names, ";", collapse = "\n")
      )
    })

    return(c(
      overall_lines,
      unlist(class_lines)
    ))
  }

  if (is.character(constraint)) {
    switch(constraint,
           "E0" = {
             overall_lines <- c(
               sapply(seq_along(variable.names), function(i) {
                 sprintf("%s (%d);", variable.names[i], i)
               }),
               apply(variable.pairs, 1, function(pair) {
                 sprintf("%s WITH %s@0;", pair[1], pair[2])
               })
             )

             class_lines <- lapply(1:L, function(k) {
               c(
                 paste0("%", class.variable, "#", k, "%"),
                 paste0("[", variable.names, "];")
               )
             })

             return(c("%OVERALL%", overall_lines, unlist(class_lines)))
           },
           "V0" = {
             class_lines <- lapply(1:L, function(k) {
               c(
                 paste0("%", class.variable, "#", k, "%"),
                 paste0("[", variable.names, "];"),
                 paste0(variable.names, ";"),
                 apply(variable.pairs, 1, function(pair) {
                   sprintf("%s WITH %s@0;", pair[1], pair[2])
                 })
               )
             })

             return(c(unlist(class_lines)))
           },
           "EE" = {
             overall_lines <- character(0)
             overall_lines <- c(overall_lines, sapply(seq_along(variable.names), function(i) {
               sprintf("%s (%d);", variable.names[i], i)
             }))

             cov_labels <- (I + 1):(I + nrow(variable.pairs))
             overall_lines <- c(overall_lines, apply(variable.pairs, 1, function(pair) {
               idx <- which(apply(variable.pairs, 1, function(x) all(sort(x) == sort(pair))))
               sprintf("%s WITH %s (%d);", pair[1], pair[2], cov_labels[idx])
             }))

             class_lines <- lapply(1:L, function(k) {
               c(
                 paste0("%", class.variable, "#", k, "%"),
                 paste0("[", variable.names, "];")
               )
             })

             return(c("%OVERALL%", overall_lines, unlist(class_lines)))
           },
           "EV" = {
             overall_lines <- sapply(seq_along(variable.names), function(i) {
               sprintf("%s (%d);", variable.names[i], i)
             })

             class_lines <- lapply(1:L, function(k) {
               cov_lines <- apply(variable.pairs, 1, function(pair) {
                 paste0(pair[1], " WITH ", pair[2], ";")
               })

               c(
                 paste0("%", class.variable, "#", k, "%"),
                 paste0("[", variable.names, "];"),
                 cov_lines
               )
             })

             return(c("%OVERALL%", overall_lines, unlist(class_lines)))
           },
           "VV" = {
             class_lines <- lapply(1:L, function(k) {
               cov_lines <- apply(variable.pairs, 1, function(pair) {
                 paste0(pair[1], " WITH ", pair[2], ";")
               })

               c(
                 paste0("%", class.variable, "#", k, "%"),
                 paste0("[", variable.names, "];"),
                 paste0(variable.names, ";"),
                 cov_lines
               )
             })

             return(c(unlist(class_lines)))
           },
           {
             stop("Unsupported constraint type: ", constraint)
           }
    )
  } else if (is.list(constraint)) {
    constraints <- lapply(constraint, function(x) {
      if (length(x) != 2 || !is.numeric(x)) {
        stop("Each constraint must be a numeric vector of length 2")
      }
      x <- as.integer(x)
      if (any(x < 1) || any(x > I)) {
        stop("Constraint indices must be between 1 and ", I)
      }
      sort(x)
    })

    constraint_keys <- sapply(constraints, function(x) paste(sort(x), collapse = ":"))
    unique_idx <- !duplicated(constraint_keys)
    if (sum(!unique_idx) > 0) {
      warning(sum(!unique_idx), " duplicate constraints removed")
    }
    constraints <- constraints[unique_idx]
    n_constraints <- length(constraints)
    labels <- 1000 + seq_len(n_constraints) - 1
    constraint_to_label <- setNames(labels, sapply(constraints, paste, collapse = ":"))

    var_constraints <- rep(NA_integer_, I)
    cov_constraints <- matrix(NA_integer_, nrow = I, ncol = I)

    for (con in constraints) {
      key <- paste(con, collapse = ":")
      label <- constraint_to_label[key]
      i <- con[1]
      j <- con[2]
      if (i == j) {
        if (!is.na(var_constraints[i])) {
          warning("Variable ", i, " already has a variance constraint; using first constraint")
        } else {
          var_constraints[i] <- label
        }
      } else {
        if (!is.na(cov_constraints[i, j])) {
          warning("Covariance between variables ", i, " and ", j,
                  " already constrained; using first constraint")
        } else {
          cov_constraints[i, j] <- label
          cov_constraints[j, i] <- label
        }
      }
    }

    class_lines <- lapply(1:L, function(k) {
      lines <- c(
        paste0("%", class.variable, "#", k, "%"),
        paste0("[", paste(variable.names, collapse = " "), "];")
      )

      for (i in seq_len(I)) {
        variable.name <- variable.names[i]
        label <- var_constraints[i]
        if (!is.na(label)) {
          lines <- c(lines, sprintf("%s (%d);", variable.name, label))
        } else {
          lines <- c(lines, sprintf("%s;", variable.name))
        }
      }

      for (pair_idx in seq_len(nrow(variable.pairs))) {
        v1 <- variable.pairs[pair_idx, 1]
        v2 <- variable.pairs[pair_idx, 2]
        idx1 <- match(v1, variable.names)
        idx2 <- match(v2, variable.names)
        label <- cov_constraints[idx1, idx2]
        if (!is.na(label)) {
          lines <- c(lines, sprintf("%s WITH %s (%d);", v1, v2, label))
        } else {
          lines <- c(lines, sprintf("%s WITH %s;", v1, v2))
        }
      }

      lines
    })

    return(c(unlist(class_lines)))
  } else {
    stop("Invalid constraint specification. Must be character string or list of integer vectors.")
  }
}
