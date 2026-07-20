#' Load or install a package
#'
#' Attempts to load a package; installs it first if not available.
#' Supports installation from CRAN or GitHub.
#'
#' @param pkg Character string. The name of the package to load.
#' @param github Character string or NULL. If non-NULL, the GitHub repo string
#'   (e.g., `"KlugerLab/GeneTrajectory"`) used to install the package via
#'   `devtools::install_github()`.
#'
#' @return Called for its side effect of attaching the package.
#' @keywords internal
load_or_install <- function(pkg, github = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!is.null(github)) {
      if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
      devtools::install_github(github)
    } else {
      install.packages(pkg)
    }
  }
  library(pkg, character.only = TRUE)
}


#' Inverse square root of a symmetric matrix
#'
#' Computes the inverse square root of a symmetric positive (semi-)definite
#' matrix using its eigen-decomposition.
#'
#' @param L_sym A named, symmetric numeric matrix.
#'
#' @return A symmetric numeric matrix of the same dimensions as `L_sym`,
#'   representing \eqn{L_{\text{sym}}^{-1/2}}.
#' @keywords internal
inv_sqrt <- function(L_sym){
  E.list <- rARPACK::eigs_sym(L_sym, k = nrow(L_sym))
  L_sym_inv_sqrt <- E.list$vectors %*% diag(1/sqrt(E.list$values)) %*% t(E.list$vectors)
  rownames(L_sym_inv_sqrt) = rownames(L_sym)
  colnames(L_sym_inv_sqrt) = colnames(L_sym)
  L_sym_inv_sqrt
}

#' Run a truncated symmetric eigendecomposition
#'
#' Wrapper around [rARPACK::eigs_sym()] that retains row names and returns
#' eigenvalues in descending order.
#'
#' @param K A named, symmetric positive semi-definite numeric matrix.
#' @param eig_keep Integer. Number of eigenpairs to compute. Defaults to
#'   `nrow(K)` (all eigenpairs).
#'
#' @return A list with components:
#'   \describe{
#'     \item{values}{Numeric vector of eigenvalues (descending).}
#'     \item{vectors}{Matrix of eigenvectors (columns), with row names from `K`.}
#'   }
#'
#' @export
#'
#' @importFrom rARPACK eigs_sym
#'
#' @examples
#' \dontrun{
#' evd <- RunSVD(P_diff)
#' }
RunSVD = function(K, eig_keep = nrow(K), seed = 233){
  set.seed(seed)
  E.list <- rARPACK::eigs_sym(K, k = eig_keep, which = "LM")
  rownames(E.list$vectors) = rownames(K)
  return(E.list)
}
