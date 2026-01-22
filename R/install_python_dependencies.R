#' @title Install Required Python Dependencies for Neural Latent Variable Models
#'
#' @description
#' Checks whether five essential Python packages required to run neural latent variable models
#' (e.g., LCAnet, LPAnet) are installed in the current Python environment. If any are missing,
#' the user is interactively prompted to install them via \code{reticulate::py_install()}.
#' The targeted packages are:
#' \itemize{
#'   \item \code{numpy} — Fundamental package for numerical computing in Python.
#'   \item \code{torch} — PyTorch deep learning framework (supports CPU/GPU computation).
#'   \item \code{matplotlib} — 2D plotting and visualization library.
#'   \item \code{scikit-learn} — Machine learning utilities (used here primarily for KMeans initialization).
#'   \item \code{scipy} — Scientific computing and advanced linear algebra routines.
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
#'   \item If not available, prints a message describing the package’s purpose.
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

  numpy_installed <- py_module_available("numpy")
  torch_installed <- py_module_available("torch")
  matplotlib_installed <- py_module_available("matplotlib")
  sklearn_installed <- py_module_available("sklearn")
  scipy_installed <- py_module_available("scipy")

  # Helper: get OS type
  os_type <- Sys.info()["sysname"]
  is_windows <- grepl("Windows", os_type, ignore.case = TRUE)
  is_mac <- grepl("Darwin", os_type, ignore.case = TRUE)
  is_linux <- grepl("Linux", os_type, ignore.case = TRUE)

  if (!numpy_installed) {
    message("The 'numpy' library is not installed. (Purpose: Numerical computations)")
    install_numpy <- readline(prompt = "Do you want to install 'numpy'? (y/n): ")
    if (tolower(install_numpy) == "y") {
      py_install("numpy")
      message("numpy has been installed.")
    } else {
      message("numpy installation skipped.")
    }
  }

  if (!torch_installed) {
    message("The 'torch' library is not installed. (Purpose: Neural network estimation - CPU or GPU version)")

    install_torch <- readline(prompt = "Do you want to install 'torch'? (y/n): ")
    if (tolower(install_torch) == "y") {

      # Ask CPU vs GPU
      if (is_mac) {
        message("You are on macOS. GPU support requires Apple Silicon (MPS) or external GPU via ROCm (not officially supported).")
        message("PyTorch on macOS defaults to CPU/MPS backend. Installing standard version.")
        py_install(c("torch", "torchvision"))
        message("torch (CPU/MPS) has been installed.")
      } else {
        gpu_choice <- readline(prompt = "Install GPU version? (y for GPU / n for CPU): ")

        if (tolower(gpu_choice) == "y" && !is_mac) {
          # Define available CUDA versions with URLs
          cuda_options <- list(
            cu128 = "https://download.pytorch.org/whl/cu128",
            cu126 = "https://download.pytorch.org/whl/cu126",
            cu130 = "https://download.pytorch.org/whl/cu130"
          )

          # Filter by OS (Windows/Linux only; macOS skipped above)
          if (is_windows || is_linux) {
            message("Available CUDA versions:")
            for (name in names(cuda_options)) {
              message(sprintf("  %s", name))
            }

            selected_cuda <- readline(prompt = "Enter CUDA version (e.g., cu128, cu126, cu130): ")
            selected_cuda <- trimws(selected_cuda)

            if (selected_cuda %in% names(cuda_options)) {
              index_url <- cuda_options[[selected_cuda]]
              message(sprintf("Installing torch with CUDA %s support...", gsub("cu", "", selected_cuda)))
              py_install(
                packages = c("torch", "torchvision"),
                pip = TRUE,
                extra_args = c("--index-url", index_url)
              )
              message(sprintf("torch (CUDA %s) has been installed.", gsub("cu", "", selected_cuda)))
            } else {
              message("Invalid CUDA version. Installing CPU version instead.")
              py_install(c("torch", "torchvision"))
            }
          } else {
            message("Unsupported OS for CUDA. Installing CPU version.")
            py_install(c("torch", "torchvision"))
          }
        } else {
          message("Installing CPU version of torch...")
          py_install(c("torch", "torchvision"))
          message("torch (CPU) has been installed.")
        }
      }
    } else {
      message("torch installation skipped.")
    }
  }

  if (!matplotlib_installed) {
    message("The 'matplotlib' library is not installed. (Purpose: Visualization)")
    install_matplotlib <- readline(prompt = "Do you want to install 'matplotlib'? (y/n): ")
    if (tolower(install_matplotlib) == "y") {
      py_install("matplotlib")
      message("matplotlib has been installed.")
    } else {
      message("matplotlib installation skipped.")
    }
  }

  if (!sklearn_installed) {
    message("The 'scikit-learn' library is not installed. (Purpose: K-means initialization)")
    install_sklearn <- readline(prompt = "Do you want to install 'scikit-learn'? (y/n): ")
    if (tolower(install_sklearn) == "y") {
      py_install("scikit-learn")
      message("scikit-learn has been installed.")
    } else {
      message("scikit-learn installation skipped.")
    }
  }

  if (!scipy_installed) {
    message("The 'scipy' library is not installed. (Purpose: Linear algebra operations)")
    install_scipy <- readline(prompt = "Do you want to install 'scipy'? (y/n): ")
    if (tolower(install_scipy) == "y") {
      py_install("scipy")
      message("scipy has been installed.")
    } else {
      message("scipy installation skipped.")
    }
  }

  return(list(
    numpy_installed = numpy_installed,
    torch_installed = torch_installed,
    matplotlib_installed = matplotlib_installed,
    sklearn_installed = sklearn_installed,
    scipy_installed = scipy_installed
  ))
}
