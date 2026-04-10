#' Align Latent Class/Profile Models via Optimal Permutation
#'
#' This function reorders the latent classes/profiles of \code{object2} to best match those in \code{object1}
#' by minimizing the total assignment cost based on posterior class membership (MAP classification).
#' It uses the Linear Sum Assignment Problem (LSAP) solver to find the optimal one-to-one mapping between latent classes.
#' Useful for comparing or averaging models across replications, initializations, or algorithms where class labels may be permuted.
#'
#' @param object1 An object of class \code{"LCA"} or \code{"LPA"}, typically the reference model.
#' @param object2 An object of class \code{"LCA"} or \code{"LPA"}, whose latent classes will be reordered to align with \code{object1}.
#'
#' @return A modified version of \code{object2}, with all parameters and posterior probabilities reordered
#'         to best match the latent class structure of \code{object1}. The returned object retains its original class
#'         (\code{"LCA"} or \code{"LPA"}) and includes aligned:
#'         \itemize{
#'           \item Prior probabilities (\code{P.Z})
#'           \item Posterior probabilities (\code{P.Z.Xn})
#'           \item MAP classifications (\code{Z})
#'           \item Profile/Class-specific parameters:
#'             \itemize{
#'               \item For \code{"LPA"}: \code{means}, \code{covs}
#'               \item For \code{"LCA"}: \code{par}, \code{probability}
#'             }
#'           \item All relevant \code{dimnames} and \code{names} are synchronized with \code{object1}
#'         }
#'
#' @details
#' The alignment is performed by:
#' \enumerate{
#'   \item Computing Maximum A Posteriori (MAP) classification matrices for both models.
#'   \item Calculating a distance matrix between classes (typically Euclidean distance between binary MAP vectors).
#'   \item Solving the Linear Sum Assignment Problem (LSAP) via \code{\link[clue]{solve_LSAP}} to find the permutation minimizing total mismatch cost.
#'   \item Reordering all class-specific components of \code{object2} according to this optimal assignment.
#' }
#'
#' @note
#' - Both models must have identical numbers of observations (\eqn{N}), latent classes (\eqn{L}), and indicators (\eqn{I}).
#' - Designed for use after fitting multiple models (e.g., different random starts) to ensure consistent class labeling.
#' - Does not modify \code{object1}; only returns a reordered \code{object2}.
#'
#' @examples
#' \dontrun{
#' # need Mplus and Pyrthon
#'
#' library(LCPA)
#' set.seed(123)
#'
#' data.obj <- sim.LCA(N = 500, I = 4, L = 3)
#'
#' # Fit two models with different random seeds
#' fit1 <- LCA(data.obj$response, L = 3, method = "Mplus", nrep = 1)
#' fit2 <- LCA(data.obj$response, L = 3, method = "NNE", nrep = 1)
#'
#' # Align fit2 to fit1's class ordering
#' fit2_aligned <- adjust.model(fit1, fit2)
#'
#' # Compare prior probabilities before and after
#' print("Before alignment:")
#' print(fit2$params$P.Z)
#' print("After alignment:")
#' print(fit2_aligned$params$P.Z)
#'
#' }
#'
#' @importFrom clue solve_LSAP
#' @export
#'
adjust.model <- function(object1, object2) {

  # Validate input classes
  model.type1 <- class(object1)[1]
  model.type2 <- class(object2)[1]

  if (model.type1 != model.type2) {
    stop("'object1' and 'object2' must be of the same model type (both 'LCA' or both 'LPA').",
         call. = FALSE)
  }

  model.type <- model.type1
  if (!(model.type %in% c("LCA", "LPA"))) {
    stop("Only 'LCA' and 'LPA' model types are supported.", call. = FALSE)
  }

  # Validate dimensions
  N1 <- nrow(object1$P.Z.Xn)
  N2 <- nrow(object2$P.Z.Xn)
  if (N1 != N2) {
    stop(sprintf("Number of observations must match: object1 has %d, object2 has %d.", N1, N2),
         call. = FALSE)
  }
  N <- N1

  L1 <- ncol(object1$P.Z.Xn)
  L2 <- ncol(object2$P.Z.Xn)
  if (L1 != L2) {
    stop(sprintf("Number of latent classes must match: object1 has %d, object2 has %d.", L1, L2),
         call. = FALSE)
  }
  L <- L1

  I1 <- ncol(object1$arguments$response)
  I2 <- ncol(object2$arguments$response)
  if (I1 != I2) {
    stop(sprintf("Number of indicators must match: object1 has %d, object2 has %d.", I1, I2),
         call. = FALSE)
  }
  I <- I1

  # Compute MAP classification matrices (binary indicator: 1 if assigned to that class)
  MAP1 <- matrix(apply(object1$P.Z.Xn, 1, max), N, L, byrow = FALSE) == object1$P.Z.Xn
  MAP2 <- matrix(apply(object2$P.Z.Xn, 1, max), N, L, byrow = FALSE) == object2$P.Z.Xn

  # Compute distance matrix between classes (rows = object1 classes, cols = object2 classes)
  dist.mat <- distance.matrix(MAP1, MAP2)

  # Solve LSAP to get optimal assignment (which class in object2 maps to which in object1)
  assignment <- as.numeric(clue::solve_LSAP(dist.mat))  # assignment[j] = i means class j in object2 → class i in object1

  # Create inverse mapping: from object2's original index to new aligned index
  # If assignment = c(2, 3, 1), then reorder object2 so that its class 1→2, 2→3, 3→1
  # We need: new_object2_class[i] = old_object2_class[assignment[i]]
  # So we reorder using assignment as the new order
  # Example: if assignment = c(2,3,1), then we want object2[ ,c(2,3,1)] — so just use assignment directly

  # Reorder object2 to match object1's class ordering
  if (model.type == "LPA") {
    object2$params$means  <- object2$params$means[assignment, , drop = FALSE]
    object2$params$covs   <- object2$params$covs[, , assignment, drop = FALSE]
    object2$params$P.Z    <- object2$params$P.Z[assignment]

    dimnames(object2$params$means) <- dimnames(object1$params$means)
    dimnames(object2$params$covs)  <- dimnames(object1$params$covs)
    names(object2$params$P.Z)      <- names(object1$params$P.Z)
  } else if (model.type == "LCA") {
    object2$params$par    <- object2$params$par[assignment, , , drop = FALSE]
    object2$params$P.Z    <- object2$params$P.Z[assignment]

    dimnames(object2$params$par) <- dimnames(object1$params$par)
    names(object2$params$P.Z)    <- names(object1$params$P.Z)

    for (i in 1:I) {
      object2$probability[[i]] <- object2$probability[[i]][assignment, , drop = FALSE]
      dimnames(object2$probability[[i]]) <- dimnames(object1$probability[[i]])
    }
  }

  # Reorder posterior probabilities and MAP assignments
  object2$P.Z.Xn <- object2$P.Z.Xn[, assignment, drop = FALSE]
  dimnames(object2$P.Z.Xn) <- dimnames(object1$P.Z.Xn)

  object2$P.Z <- object2$P.Z[assignment]
  names(object2$P.Z) <- names(object1$P.Z)

  # Update Z (MAP class memberships): map old class labels to new ones
  # Since we reordered classes: old class 'k' is now at position 'which(assignment == k)'
  # But easier: create reverse lookup table
  reverse_map <- integer(L)
  reverse_map[assignment] <- 1:L
  object2$Z <- reverse_map[object2$Z]

  return(object2)
}
