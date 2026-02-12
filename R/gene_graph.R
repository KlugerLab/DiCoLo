#' Build a symmetrized Gaussian affinity matrix
#'
#' Converts a pairwise distance matrix into a symmetrised Gaussian kernel
#' affinity matrix, then applies a symmetric normalisation (like a normalised
#' graph Laplacian step).
#'
#' @param dist.mat A square numeric matrix of pairwise distances. Row and
#'   column names should be gene (or sample) identifiers.
#' @param K Integer. The bandwidth is set to the distance to the K-th nearest
#'   neighbour. Ignored when `sigma` is provided. Default `10`.
#' @param sigma Numeric or NULL. Fixed kernel bandwidth. If NULL (default),
#'   the bandwidth is determined adaptively per row using the K-th nearest
#'   neighbour distance.
#' @param nEV Integer. Kept for API compatibility; not currently used internally.
#'   Default `30`.
#' @param t Integer. Kept for API compatibility; not currently used internally.
#'   Default `1`.
#'
#' @return A symmetric, normalised affinity matrix with the same row/column
#'   names as `dist.mat`.
#'
#' @keywords internal
affinity.mtx <- function(dist.mat, K = 10, sigma = NULL, nEV = 30, t = 1) {
  K <- min(K, nrow(dist.mat))
  dists <- as.matrix(dist.mat)
  sigma.list <- c()
  if (is.null(sigma)) {
    for (i in 1:nrow(dists)) {
      sigma.list[i] <- sort(dists[i, ])[K+1]
    }
  }
  if (!is.null(sigma)) {
    sigma.list <- rep(sigma, nrow(dists))
  }
  affinity.matrix <- matrix(0, nrow = nrow(dists), ncol = ncol(dists))
  for (i in 1:nrow(affinity.matrix)) {
    dist.vec <- dists[i, ]
    dist.vec[is.na(dist.vec)] <- 10^6
    affinity.matrix[i, ] <- exp(-(dist.vec/sigma.list[i])^2 )
  }
  affinity.matrix.2 <- (affinity.matrix + t(affinity.matrix))/2
  
  normalized.vec <- 1/sqrt(rowSums(affinity.matrix.2))
  D <- diag(normalized.vec)
  affinity.matrix.3 <- D %*% affinity.matrix.2 %*% D
  rownames(affinity.matrix.3) = colnames(affinity.matrix.3) = rownames(dist.mat)
  return(affinity.matrix.3)
}


#' Compute a positive semi-definite gene graph operator
#'
#' Builds a symmetrised Gaussian affinity matrix from a gene distance matrix
#' and projects any negative eigenvalues up to a small positive threshold to
#' ensure positive semi-definiteness.
#'
#' @param gene.dist.mat A square numeric distance matrix (e.g., from
#'   [LoadGeneEMD()]).
#' @param knn Integer. Number of nearest neighbours used to set the adaptive
#'   kernel bandwidth. Default `10`.
#' @param shrink_thred Numeric. Eigenvalues below this threshold are replaced
#'   by this value. Default `1e-10`.
#'
#' @return A symmetric positive semi-definite matrix (the graph operator) with
#'   the same row/column names as `gene.dist.mat`.
#'
#' @export
#'
#' @importFrom rARPACK eigs_sym
#'
#' @examples
#' \dontrun{
#' G <- ComputeGraphOperator(gene.dist.mat, knn = 10)
#' }
ComputeGraphOperator <- function(gene.dist.mat, knn = 10, shrink_thred = 1e-10) {
  gene_graph <- affinity.mtx(gene.dist.mat, K = knn)

  E.list <- rARPACK::eigs_sym(gene_graph, k = nrow(gene_graph))
  if (min(E.list$values) > 0) {
    return(gene_graph)
  }

  gene_graph_shrink <- E.list$vectors %*%
    diag(pmax(E.list$values, shrink_thred)) %*%
    t(E.list$vectors)

  rownames(gene_graph_shrink) <- colnames(gene_graph_shrink) <- rownames(gene_graph)
  return(gene_graph_shrink)
}


#' Compute the differential co-localization operator
#'
#' Computes the differential operator \eqn{P_{\text{diff}} = P_{\text{proj}}^{-1/2}
#' P_{\text{comp}} P_{\text{proj}}^{-1/2}} using regularisation.
#' The regularisation parameter is chosen automatically at the knee point
#' of the eigenvalue spectrum of `gene_graph_proj`.
#'
#' @param gene_graph_comp A symmetric positive semi-definite matrix representing
#'   the comparison graph operator (e.g., from [ComputeGraphOperator()]).
#' @param gene_graph_proj A symmetric positive semi-definite matrix representing
#'   the projection/reference graph operator.
#' @param regulator_min Numeric. Minimum allowed regularisation parameter to
#'   avoid numerical instability. Default `1e-2`.
#'
#' @return A symmetric numeric matrix — the differential operator — with the
#'   same dimensions as the input matrices.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' P_diff <- ComputeDifferentialOperator(gene_graph_comp, gene_graph_proj)
#' }
ComputeDifferentialOperator <- function(gene_graph_comp, gene_graph_proj,
                                        regulator_min = 1e-2) {
  E.list    <- RunSVD(gene_graph_proj)
  regulator <- E.list$values[FindKneePoint(E.list$values, plot.fig = FALSE)]
  regulator <- max(regulator, regulator_min)

  P_inv  <- inv_sqrt(gene_graph_proj + regulator * diag(nrow(gene_graph_proj)))
  P_diff <- P_inv %*% gene_graph_comp %*% P_inv

  return(P_diff)
}


