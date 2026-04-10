#' @title Install Required Python Dependencies for Neural Latent Variable Models
#'
#' @description
#' Checks whether essential Python packages required to run neural latent variable models
#' (e.g., LCAnet, LPAnet) are installed in the current Python environment. If any are missing,
#' the user is interactively prompted to install them via \code{reticulate::py_install()}.
#' The targeted packages are:
#' \itemize{
#'   \item \code{numpy} — Fundamental package for numerical computing in Python.
#'   \item \code{torch} — PyTorch deep learning framework (supports CPU/GPU computation).
#'   \item \code{matplotlib} — 2D plotting and visualization library.
#'   \item \code{scikit-learn} — Machine learning utilities.
#'   \item \code{scipy} — Scientific computing and advanced linear algebra routines.
#'   \item \code{six} — Python 3 compatibility library.
#' }
#'
#' For \code{torch}, users can choose between CPU-only or GPU-enabled versions (with CUDA support).
#' Available CUDA versions are filtered by OS compatibility.
#'
#' This function is especially useful when deploying models that bridge R and Python via \pkg{reticulate},
#' ensuring all backend dependencies are met before model execution.
#'
#' @details
#' The function performs the following steps for each dependency:
#' \enumerate{
#'   \item Uses \code{reticulate::py_module_available()} to test if the module is importable.
#'   \item If not available, prints a message describing the package's purpose.
#'   \item Prompts the user interactively (via \code{readline}) whether to proceed with installation.
#'   \item For \code{torch}, offers CPU/GPU choice and CUDA version selection if GPU is chosen.
#'   \item Installs the package using \code{reticulate::py_install()} with appropriate index URL if needed.
#'   \item Returns a logical list indicating initial installation status of each package.
#' }
#'
#' \strong{Note:} This function requires \pkg{reticulate} to be loaded and a valid Python environment configured.
#' It does NOT automatically install \pkg{reticulate} or configure Python — that must be done separately.
#'
#' @return
#' A named list of logical values indicating whether each package was already installed before running this function:
#' \item{numpy_installed}{Logical. Was \code{numpy} already available?}
#' \item{torch_installed}{Logical. Was \code{torch} already available?}
#' \item{matplotlib_installed}{Logical. Was \code{matplotlib} already available?}
#' \item{sklearn_installed}{Logical. Was \code{scikit-learn} already available?}
#' \item{scipy_installed}{Logical. Was \code{scipy} already available?}
#' \item{six_installed}{Logical. Was \code{six} already available?}
#'
#' @examples
#' library(reticulate)
#'
#' # Ensure reticulate is loaded and Python is configured
#' # need python
#' \dontrun{
#' # Run dependency installer
#' deps <- install_python_dependencies()
#'
#' # Check which were missing
#' print(deps)
#' }
#'
#' @importFrom reticulate py_module_available py_install
#' @importFrom utils packageVersion
#' @export
install_python_dependencies <- function() {
  py_exe <- py_config()$python

  if (is.null(py_exe) || !file.exists(py_exe)) {
    stop("Unable to obtain a valid Python executable path, please call use_python() first to specify.")
  }

  message("Using Python: ", py_exe)

  numpy_installed   <- py_module_available("numpy")
  torch_installed   <- py_module_available("torch")
  matplotlib_installed <- py_module_available("matplotlib")
  sklearn_installed <- py_module_available("sklearn")
  scipy_installed   <- py_module_available("scipy")
  six_installed     <- py_module_available("six")

  os_type    <- Sys.info()["sysname"]
  is_windows <- grepl("Windows", os_type, ignore.case = TRUE)
  is_mac     <- grepl("Darwin", os_type, ignore.case = TRUE)
  is_linux   <- grepl("Linux", os_type, ignore.case = TRUE)

  install_with_pip <- function(packages, extra_args = character(0)) {
    args <- c("-m", "pip", "install", "--user")
    args <- c(args, extra_args, packages)
    result <- system2(py_exe, args = args, stdout = TRUE, stderr = TRUE)

    cat(paste(result, collapse = "\n"), "\n")

    if (any(grepl("Successfully installed", tail(result, 10), ignore.case = TRUE))) {
      message("Installation successful.")
    } else {
      message("Installation may have failed, please check the output above.")
    }
  }

  if (!numpy_installed) {
    message("'numpy' not detected.")
    install_numpy <- readline(prompt = "Install 'numpy'? (y/n): ")
    if (tolower(install_numpy) == "y") {
      install_with_pip("numpy")
      message("numpy installation completed (or already exists).")
    } else {
      message("numpy installation skipped.")
    }
  }

  if (!torch_installed) {
    message("'torch' not detected.")
    install_torch <- readline(prompt = "Install 'torch'? (y/n): ")
    if (tolower(install_torch) == "y") {
      if (is_mac) {
        message("macOS environment, installing CPU/MPS version by default.")
        install_with_pip(c("torch", "torchvision"))
      } else {
        gpu_choice <- readline(prompt = "Install GPU version? (y=GPU / n=CPU): ")
        if (tolower(gpu_choice) == "y" && (is_windows || is_linux)) {
          cuda_options <- list(
            cu128 = "https://download.pytorch.org/whl/cu128 ",
            cu126 = "https://download.pytorch.org/whl/cu126 ",
            cu130 = "https://download.pytorch.org/whl/cu130 "
          )
          message("Available CUDA versions:")
          for (name in names(cuda_options)) message("  ", name)

          selected_cuda <- readline(prompt = "Enter CUDA version (e.g., cu128, cu126, cu130): ")
          selected_cuda <- trimws(selected_cuda)

          if (selected_cuda %in% names(cuda_options)) {
            index_url <- cuda_options[[selected_cuda]]
            message(sprintf("Installing torch with CUDA %s support...", gsub("cu", "", selected_cuda)))
            install_with_pip(
              c("torch", "torchvision"),
              extra_args = c("--index-url", index_url)
            )
          } else {
            message("Invalid CUDA version, installing CPU version instead.")
            install_with_pip(c("torch", "torchvision"))
          }
        } else {
          message("Installing CPU version of torch...")
          install_with_pip(c("torch", "torchvision"))
        }
      }
    } else {
      message("torch installation skipped.")
    }
  }

  if (!matplotlib_installed) {
    message("'matplotlib' not detected.")
    install_matplotlib <- readline(prompt = "Install 'matplotlib'? (y/n): ")
    if (tolower(install_matplotlib) == "y") {
      install_with_pip("matplotlib")
      message("matplotlib installation completed.")
    } else {
      message("matplotlib installation skipped.")
    }
  }

  # scikit-learn
  if (!sklearn_installed) {
    message("'scikit-learn' not detected.")
    install_sklearn <- readline(prompt = "Install 'scikit-learn'? (y/n): ")
    if (tolower(install_sklearn) == "y") {
      install_with_pip("scikit-learn")
      message("scikit-learn installation completed.")
    } else {
      message("scikit-learn installation skipped.")
    }
  }

  if (!scipy_installed) {
    message("'scipy' not detected.")
    install_scipy <- readline(prompt = "Install 'scipy'? (y/n): ")
    if (tolower(install_scipy) == "y") {
      install_with_pip("scipy")
      message("scipy installation completed.")
    } else {
      message("scipy installation skipped.")
    }
  }

  if (!six_installed) {
    message("'six' not detected.")
    install_six <- readline(prompt = "Install 'six'? (y/n): ")
    if (tolower(install_six) == "y") {
      install_with_pip("six")
      message("six installation completed.")
    } else {
      message("six installation skipped.")
    }
  }

  return(list(
    numpy_installed   = py_module_available("numpy"),
    torch_installed   = py_module_available("torch"),
    matplotlib_installed = py_module_available("matplotlib"),
    sklearn_installed = py_module_available("sklearn"),
    scipy_installed   = py_module_available("scipy"),
    six_installed     = py_module_available("six")
  ))
}
