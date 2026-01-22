#' Generate a Random Correlation Matrix via C-Vine Partial Correlations
#'
#' This function generates a random \eqn{I \times I} correlation matrix using the C-vine partial correlation
#' parameterization described in Joe & Kurowicka (2026). The method constructs the matrix recursively using
#' partial correlations organized in a C-vine structure, with distributional properties controlled by LKJ
#' concentration and skewness parameters.
#'
#' @param I Dimension of the correlation matrix (must be \eqn{I \geq 1}).
#' @param eta LKJ concentration parameter (\eqn{\eta > 0}). When \eqn{\eta = 1} and \eqn{\text{skew} = 0}, the
#'            distribution is uniform over correlation matrices. Larger \eqn{\eta} values concentrate
#'            mass near the identity matrix. Critical for positive definiteness: Requires \eqn{\eta > (I-2)/2}
#'            to theoretically guarantee positive definiteness (Theorem 1, Joe & Kurowicka 2026). Default is 1.
#' @param skew Skewness parameter (\eqn{-1 < \text{skew} < 1}). Controls asymmetry in the partial correlation distribution:
#'   - \eqn{\text{skew} > 0}: Biased toward positive partial correlations
#'   - \eqn{\text{skew} < 0}: Biased toward negative partial correlations
#'   - \eqn{\text{skew} = 0}: Symmetric distribution (default)
#' @param positive Logical. If \code{TRUE}, restricts partial correlations to \eqn{(0,1)} and enforces positive definiteness. Default is \code{FALSE}.
#' @param permute Logical. If \code{TRUE}, applies a random permutation to rows/columns to ensure exchangeability (invariance to variable ordering). Default is \code{TRUE}.
#' @param maxiter Integer. Maximum number of generation attempts before numerical adjustment when \code{positive = TRUE}. Default is 10.
#'
#' @return An \eqn{I \times I} positive definite correlation matrix with unit diagonal.
#'
#' @details
#' The algorithm follows four key steps:
#'
#' 1. Partial correlation sampling:
#' For tree level \eqn{k = 1, \dots, I-1} and node \eqn{j = k+1, \dots, I}, partial correlations \eqn{\rho_{k,j \mid 1:(k-1)}} are sampled as:
#' \deqn{
#'   \alpha_k = \eta + \frac{I - k - 1}{2}, \quad
#'   a_k = \alpha_k (1 + \text{skew}), \quad
#'   b_k = \alpha_k (1 - \text{skew})
#' }
#' \itemize{
#'   \item If \code{positive = FALSE}:
#'   \deqn{\rho_{k,j} \sim 2 \cdot \mathrm{Beta}(a_k, b_k) - 1}
#'   \item If \code{positive = TRUE}:
#'   \deqn{\rho_{k,j} \sim \mathrm{Beta}(a_k, b_k)}
#' }
#'
#' 2. Recursive matrix construction (C-vine):
#' The correlation matrix \eqn{\mathbf{R}} is built without matrix inversion using backward recursion:
#' \itemize{
#'   \item Tree 1 (raw correlations): \eqn{R_{1j} = \rho_{1,j}} for \eqn{j = 2,\dots,I}
#'   \item Trees \eqn{l \geq 2}: For pairs \eqn{(l,j)} where \eqn{l = 2,\dots,I-1} and \eqn{j = l+1,\dots,I}:
#'   \deqn{
#'     c \gets \rho_{l,j \mid 1:(l-1)} \\
#'     \text{for } k = l-1 \text{ down to } 1: \\
#'     \quad c \gets c \cdot \sqrt{(1 - \rho_{k,l}^2)(1 - \rho_{k,j}^2)} + \rho_{k,l} \cdot \rho_{k,j} \\
#'     R_{lj} \gets c
#'   }
#' }
#' This implements the dynamic programming approach from Joe & Kurowicka (2026, Section 2.1).
#'
#' 3. Positive definiteness enforcement (when \code{positive = TRUE}):
#' \itemize{
#'   \item Attempt up to \code{maxiter} generations
#'   \item On failure, project to nearest positive definite correlation matrix using \code{\link[Matrix]{nearPD}} with \code{corr = TRUE}
#'   \item Final matrix has minimum eigenvalue > 1e-8
#' }
#'
#' 4. Exchangeability (optional):
#' If \code{permute = TRUE}, rows/columns are randomly permuted before returning the matrix.
#'
#' @note When \code{positive = TRUE}, the function guarantees positive definiteness either through direct generation
#' (with retries) or numerical projection. The theoretical guarantee \eqn{\eta > (I-2)/2} is recommended for high dimensions.
#'
#' @references
#' Joe, H., & Kurowicka, D. (2026). Random correlation matrices generated via partial correlation C-vines. Journal of Multivariate Analysis, 211, 105519. https://doi.org/10.1016/j.jmva.2025.105519
#'
#' @examples
#' # Default 3x3 correlation matrix
#' sim.correlation(3)
#'
#' # 5x5 matrix concentrated near identity (eta=3)
#' sim.correlation(5, eta = 3)
#'
#' # Skewed toward positive correlations (no permutation)
#' sim.correlation(4, skew = 0.7, permute = FALSE)
#'
#' # Positive partial correlations (enforced positive definiteness)
#' R <- sim.correlation(6, positive = TRUE)
#' min(eigen(R, symmetric = TRUE, only.values = TRUE)$values)  # > 1e-8
#'
#' # High-dimensional case (I=20) with theoretical guarantee
#' R <- sim.correlation(20, eta = 10)  # eta=10 > (20-2)/2=9
#' min(eigen(R, symmetric = TRUE, only.values = TRUE)$values)
#'
#' @importFrom stats rbeta
#' @export
sim.correlation <- function(I,
                            eta = 1,
                            skew = 0,
                            positive = FALSE,
                            permute = TRUE,
                            maxiter = 10) {

  # Handle trivial case immediately
  if (I < 2) return(matrix(1, 1, 1))

  # Validate input parameters
  if (eta <= 0) stop("eta must be > 0")
  if (abs(skew) >= 1) stop("skew must satisfy -1 < skew < 1")
  if (maxiter < 1) stop("maxiter must be at least 1")

  # Determine if positive definiteness enforcement is needed
  enforce_pd <- positive
  attempt <- 1
  success <- FALSE
  R_final <- NULL

  # Attempt matrix generation with retries for positive definite case
  while (attempt <= maxiter && !success) {
    # Initialize matrices
    R_current <- diag(I)  # Working correlation matrix
    partial_cor <- matrix(0, nrow = I, ncol = I)  # Stores partial correlations rho_{k,j}

    # Step 1: Generate partial correlations for each tree level
    for (tree_level in 1:(I - 1)) {
      # Tree-specific concentration parameter
      alpha_k <- eta + (I - tree_level - 1) / 2
      # Skewed Beta distribution parameters
      a_k <- alpha_k * (1 + skew)
      b_k <- alpha_k * (1 - skew)

      # Validate distribution parameters
      if (a_k <= 0 || b_k <= 0) {
        stop(sprintf("Invalid Beta parameters at tree %d: a_k=%.3f, b_k=%.3f. ",
                     tree_level, a_k, b_k),
             "Increase eta or reduce |skew| for high dimensions.")
      }

      # Generate partial correlations for current tree
      for (node in (tree_level + 1):I) {
        if (positive) {
          # Sample from Beta(a_k, b_k) on (0,1)
          rho_val <- rbeta(1, a_k, b_k)
        } else {
          # Transform Beta to (-1,1) interval
          rho_val <- 2 * rbeta(1, a_k, b_k) - 1
        }
        partial_cor[tree_level, node] <- rho_val
      }
    }

    # Step 2: Construct correlation matrix via C-vine recursion
    # Tree 1: Direct correlations with root node (variable 1)
    for (j in 2:I) {
      R_current[1, j] <- R_current[j, 1] <- partial_cor[1, j]
    }

    # Higher trees (l >= 2): Recursive update through conditioning sets
    if (I > 2) {
      for (l in 2:(I - 1)) {  # Current row index
        for (j in (l + 1):I) {  # Column index > l
          current_val <- partial_cor[l, j]  # Start with tree-l partial correlation

          # Backward recursion through previous trees (k = l-1 down to 1)
          for (k in (l - 1):1) {
            rho_kl <- partial_cor[k, l]  # Partial correlation (k,l)
            rho_kj <- partial_cor[k, j]  # Partial correlation (k,j)

            # Recursive update (Joe & Kurowicka 2026, Section 2.1)
            current_val <- current_val * sqrt(1 - rho_kl^2) * sqrt(1 - rho_kj^2) +
              rho_kl * rho_kj
          }
          # Assign symmetric entries
          R_current[l, j] <- R_current[j, l] <- current_val
        }
      }
    }

    # Numerical safeguards: enforce symmetry and unit diagonal
    R_current <- (R_current + t(R_current)) / 2
    diag(R_current) <- 1

    # Step 3: Check positive definiteness if required
    if (enforce_pd) {
      min_eigenval <- min(eigen(R_current, symmetric = TRUE, only.values = TRUE)$values)
      if (min_eigenval > 1e-8) {
        R_final <- R_current
        success <- TRUE
      }
    } else {
      R_final <- R_current
      success <- TRUE
    }
    attempt <- attempt + 1
  }

  # Step 4: Numerical adjustment if PD enforcement failed
  if (enforce_pd && !success) {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
      stop("Package 'Matrix' required for positive definiteness adjustment. ",
           "Install with: install.packages('Matrix')")
    }

    warning(sprintf("Failed to generate PD matrix after %d attempts. ", maxiter),
            "Applying numerical projection via Matrix::nearPD(). ",
            "Result may have small numerical deviations.")

    # Project to nearest PD correlation matrix
    pd_result <- Matrix::nearPD(R_current,
                                corr = TRUE,
                                keepDiag = TRUE,
                                doDykstra = TRUE,
                                conv.tol = 1e-7,
                                maxit = 1000)
    R_final <- as.matrix(pd_result$mat)

    # Re-enforce exact unit diagonal after projection
    diag(R_final) <- 1

    # Final eigenvalue check
    min_eig <- min(eigen(R_final, symmetric = TRUE, only.values = TRUE)$values)
    if (min_eig <= 1e-8) {
      warning(sprintf("Numerical adjustment yielded min eigenvalue = %.3e. ",
                      min_eig),
              "Consider increasing eta or maxiter for strict PD guarantee.")
    }
  }

  # Step 5: Ensure exchangeability via random permutation
  if (permute) {
    perm_order <- sample(I)
    R_final <- R_final[perm_order, perm_order]
  }

  return(R_final)
}
