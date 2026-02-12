#' Find the knee point of an eigenvalue curve
#'
#' Identifies the index of maximum perpendicular distance from the line
#' connecting the first and last points of the eigenvalue sequence.
#'
#' @param eigval Numeric vector of eigenvalues (typically in descending order).
#' @param plot.fig Logical. If `TRUE`, plots the eigenvalue curve and marks
#'   the knee point. Default `FALSE`.
#'
#' @return An integer giving the index of the knee point.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' k <- FindKneePoint(evd$values, plot.fig = TRUE)
#' }
FindKneePoint <- function(eigval, plot.fig = FALSE) {
  Vecs <- eigval
  curve_data <- as.matrix(data.frame(x = seq_along(Vecs), y = Vecs))
  line_start  <- curve_data[1, ]
  line_end    <- curve_data[nrow(curve_data), ]
  dx <- line_end["x"] - line_start["x"]
  dy <- line_end["y"] - line_start["y"]
  norm_factor <- sqrt(dx^2 + dy^2)
  A <- dy
  B <- -dx
  C <- dx * line_start["y"] - dy * line_start["x"]
  distances <- abs(A * curve_data[, "x"] + B * curve_data[, "y"] + C) / norm_factor
  knee_point <- which.max(distances)

  if (isTRUE(plot.fig)) {
    plot(Vecs, pch = 20, xlab = "Index", ylab = "EigVal",
         cex.lab = 1.4, cex.axis = 1.2, cex.main = 1.5, cex = 1)
    abline(v = knee_point, col = "grey", lwd = 2, lty = "dashed")
    points(seq_len(knee_point - 1), Vecs[seq_len(knee_point - 1)],
           col = "red", pch = 20, cex = 1)
  }

  return(knee_point)
}


#' Select significant genes per eigenvector using local FDR
#'
#' For each column (eigenvector) of a gene-by-eigenvector matrix, z-transforms
#' the loadings (using median/MAD for robustness) and applies `locfdr` to
#' identify significantly loaded genes. Returns a binary indicator matrix.
#'
#' @param V A numeric matrix with genes as rows and eigenvectors as columns.
#'   Row names should be gene names.
#' @param lfdr_thresh Numeric. Local FDR threshold for significance. Default
#'   is `0.2`. Automatically capped at the empirical FDR of the first bin.
#' @param min_genes Integer. Minimum number of significant genes required in
#'   a direction (positive or negative) to retain them. Default is `3`.
#' @param max_genes Integer. Maximum number of significant genes to keep per
#'   direction per eigenvector. Default is `50`.
#'
#' @return A binary matrix of the same dimensions as `V` (genes × eigenvectors),
#'   where `1` indicates the gene is significantly loaded on that eigenvector.
#'
#' @export
#'
#' @importFrom locfdr locfdr
#'
#' @examples
#' \dontrun{
#' gene_indicator <- SelectSignificantGenes(V, lfdr_thresh = 0.2)
#' }
SelectSignificantGenes <- function(V, lfdr_thresh = 0.2,
                                   min_genes = 3,
                                   max_genes = 50) {
  stopifnot(is.matrix(V))
  if (is.null(colnames(V))) colnames(V) <- paste0("EV", seq_len(ncol(V)))
  
  eigenvec_indicator <- matrix(0, nrow(V), ncol(V))
  rownames(eigenvec_indicator) <- rownames(V)
  colnames(eigenvec_indicator) <- colnames(V)
  
  for (j in seq_len(ncol(V))) {
    v <- V[, j]
    names(v) <- rownames(V)
    
    # z-transform (MAD robust)
    med <- median(v, na.rm = TRUE)
    s   <- mad(v, center = med, constant = 1.4826, na.rm = TRUE)
    if (s < 1e-12) s <- sd(v)
    z <- (v - med) / s
    
    # Local FDR
    lf   <- suppressWarnings(locfdr(z, nulltype = 0, plot = 0))
    lfdr <- lf$fdr
    names(lfdr) <- names(v)
    
    lfdr_thresh <- min(lfdr_thresh, lf$Efdr[1])
    sel <- names(lfdr)[lfdr < lfdr_thresh]
    if (length(sel) == 0) next
    
    pos <- sel[z[sel] > 0]
    neg <- sel[z[sel] < 0]
    
    if (length(pos) > max_genes) {
      pos <- pos[order(z[pos], decreasing = TRUE)[1:max_genes]]
    }
    if (length(neg) > max_genes) {
      neg <- neg[order(-z[neg], decreasing = TRUE)[1:max_genes]]
    }
    
    if (length(pos) > min_genes) {
      eigenvec_indicator[pos, j] <- 1
    }
    if (length(neg) > min_genes) {
      eigenvec_indicator[neg, j] <- 1
    }
  }
  
  return(eigenvec_indicator)
}

#' Remove outlier genes from clusters
#'
#' For each cluster, computes each member's average distance to other cluster
#' members and flags as `NA` those that are more than `z_thresh` MADs above
#' the median average distance.
#'
#' @param gene_partition A named factor of cluster assignments (e.g., output
#'   of [ClusterGenes()]).
#' @param distM A square numeric distance matrix with names matching
#'   `names(gene_partition)`.
#' @param z_thresh Numeric. Number of MADs above the median to use as the
#'   outlier threshold. Default `3`.
#'
#' @return The input `gene_partition` factor with outlier genes set to `NA`.
#' @keywords internal
tighten_modules <- function(gene_partition, distM, z_thresh = 3) {
  unique_clusters <- levels(gene_partition)
  new_partition   <- gene_partition
  
  for (cl in unique_clusters) {
    members <- names(gene_partition[gene_partition == cl])
    if (length(members) <= 2) next
    
    sub_dist <- distM[members, members]
    avg_d    <- apply(sub_dist, 1, function(x) mean(x[-which.min(x == 0)]))
    mu       <- median(avg_d)
    sigma    <- mad(avg_d)
    outliers <- names(avg_d[avg_d > mu + z_thresh * sigma])
    
    new_partition[outliers] <- NA
  }
  
  return(new_partition)
}


#' Cluster genes using hierarchical clustering and dynamic tree cutting
#'
#' Performs agglomerative hierarchical clustering on a gene distance matrix,
#' cuts the tree dynamically with `dynamicTreeCut::cutreeDynamic`, tightens
#' clusters by removing outliers with [tighten_modules()], and drops
#' singleton/small clusters.
#'
#' @param gene_dist A `dist` object (or symmetric matrix coercible to one)
#'   of pairwise gene distances.
#' @param clustering_method Character. Linkage method passed to [hclust()].
#'   Default `"average"`.
#' @param min_gene Integer. Minimum cluster size to retain. Default `2`.
#' @param deepSplit Integer (0–4). Controls the sensitivity of dynamic tree
#'   cutting; higher values produce more, smaller clusters. Default `1`.
#' @param z_thresh Numeric. MAD outlier threshold passed to [tighten_modules()].
#'   Default `3`.
#'
#' @return A named factor of cluster assignments (levels re-indexed as
#'   consecutive integers). Genes removed as singletons or outliers are
#'   dropped from the result.
#'
#' @export
#'
#' @importFrom dynamicTreeCut cutreeDynamic
#'
#' @examples
#' \dontrun{
#' gene_partition <- ClusterGenes(as.dist(gene.dist.mat))
#' }
ClusterGenes <- function(gene_dist, clustering_method = "average",
                         min_gene = 2, deepSplit = 1, z_thresh = 3) {
  gene_hree      <- hclust(gene_dist, method = clustering_method)
  gene_partition <- dynamicTreeCut::cutreeDynamic(
    dendro              = gene_hree,
    distM               = as.matrix(gene_dist),
    deepSplit           = deepSplit,
    pamStage            = TRUE,
    pamRespectsDendro   = TRUE,
    minClusterSize      = min_gene
  )
  names(gene_partition) <- labels(gene_dist)
  gene_partition <- as.factor(gene_partition)
  gene_partition <- tighten_modules(gene_partition, as.matrix(gene_dist), z_thresh)
  gene_partition <- factor(gene_partition, levels = seq(nlevels(gene_partition)))
  
  # Drop singletons / small clusters
  small <- names(table(gene_partition))[table(gene_partition) < min_gene]
  gene_partition[gene_partition %in% small] <- NA
  gene_partition <- gene_partition[!is.na(gene_partition)]
  gene_partition <- droplevels(gene_partition)
  levels(gene_partition) <- seq_len(nlevels(gene_partition))
  
  return(gene_partition)
}


#' Rank gene modules by eigenvector importance
#'
#' Assigns a weighted importance score to each gene module based on which
#' eigenvectors (and their eigenvalue weights) its member genes are loaded on.
#' The module order is then re-indexed from most to least important.
#'
#' @param gene_partition A named factor of cluster assignments (output of
#'   [ClusterGenes()]).
#' @param gene_indicator A binary matrix (genes × eigenvectors) indicating
#'   significant gene loadings (output of [SelectSignificantGenes()]).
#' @param eigen_values Numeric vector of eigenvalues corresponding to the
#'   columns of `gene_indicator` (e.g., `RunSVD(P_diff)$values`).
#'
#' @return A named integer factor of module assignments, re-ranked so that
#'   module `1` is the most important.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' gene_partition_ranked <- RankGeneModules(gene_partition, gene_indicator, evd$values)
#' }
RankGeneModules <- function(gene_partition, gene_indicator, eigen_values) {
  eigen_values <- eigen_values[seq_len(ncol(gene_indicator))]
  eigen_vec_weight <- eigen_values / sum(eigen_values)
  names(eigen_vec_weight) <- colnames(gene_indicator)
  
  common_genes   <- intersect(names(gene_partition), rownames(gene_indicator))
  gene_partition <- gene_partition[common_genes]
  gene_indicator <- gene_indicator[common_genes, ]
  gene_indicator <- apply(gene_indicator, 1, function(x) {
    colnames(gene_indicator)[min(which(x == 1))]
  })
  
  modules        <- split(names(gene_partition), gene_partition)
  module_weights <- sapply(modules, function(genes) {
    eig_prop <- table(gene_indicator[genes]) / length(genes)
    w <- sum(eig_prop * eigen_vec_weight[names(eig_prop)])
    return(w)
  })
  
  module_order <- names(sort(module_weights, decreasing = TRUE))
  name_map     <- setNames(seq_along(module_order), module_order)
  gene_partition <- setNames(name_map[as.character(gene_partition)], names(gene_partition))
  gene_partition <- as.factor(gene_partition)
  
  return(gene_partition)
}

