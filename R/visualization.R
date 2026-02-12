#' Compute a 2-D t-SNE embedding of genes from a distance matrix
#'
#' Wraps [Rtsne::Rtsne()] for a distance-based t-SNE. The random seed is
#' fixed to `233` for reproducibility.
#'
#' @param dist A `dist` object (or symmetric matrix coercible to one) of
#'   pairwise gene distances.
#' @param perplexity Numeric. t-SNE perplexity parameter. Default `30`.
#'
#' @return A data frame with columns `tSNE_1` and `tSNE_2` and row names
#'   matching `labels(dist)`.
#'
#' @export
#'
#' @importFrom Rtsne Rtsne
#'
#' @examples
#' \dontrun{
#' embedding <- ObtainGeneTSNE(as.dist(gene.dist.mat))
#' }
ObtainGeneTSNE <- function(dist, perplexity = 30) {
  set.seed(233)
  df <- data.frame(
    Rtsne::Rtsne(dist, is_distance = TRUE, dim = 2,
                 perplexity = perplexity, max.iter = 500)$Y
  )
  colnames(df) <- c("tSNE_1", "tSNE_2")
  rownames(df) <- labels(dist)
  return(df)
}


#' Visualise a gene t-SNE coloured by module or continuous score
#'
#' Plots a 2-D t-SNE of genes coloured by a discrete gene partition (modules)
#' or a continuous numeric score. If no pre-computed embedding or partition is
#' provided, they are computed internally.
#'
#' @param dist A `dist` object of pairwise gene distances. Required if
#'   `gene_embedding` or `gene_partition` is NULL.
#' @param gene_embedding Data frame with columns `tSNE_1` and `tSNE_2` (output
#'   of [ObtainGeneTSNE()]). If NULL the embedding is computed from `dist`.
#' @param gene_partition Named factor or character/logical vector of gene module
#'   assignments. If NULL, [ClusterGenes()] is called on `dist`.
#' @param perplexity Numeric. Passed to [ObtainGeneTSNE()] when
#'   `gene_embedding` is NULL. Default `30`.
#' @param text Logical. If `TRUE` (default) and `gene_partition` is a factor,
#'   module centroid labels are overlaid.
#' @param module_color Named character vector mapping module levels to colours.
#'   If NULL, colours are generated automatically from `RColorBrewer::brewer.pal`.
#'
#' @return A `ggplot2` object.
#'
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual scale_color_gradient2
#'   theme_minimal ggtitle labs geom_text element_rect
#' @importFrom dplyr group_by summarise %>%
#' @importFrom RColorBrewer brewer.pal
#'
#' @examples
#' \dontrun{
#' p <- VisualizeGeneTSNE(dist = gene_dist, gene_partition = gene_partition)
#' print(p)
#' }
VisualizeGeneTSNE <- function(dist = NULL, gene_embedding = NULL,
                              gene_partition = NULL, perplexity = 30, text = TRUE,
                              module_color = NULL) {
  if (is.null(gene_partition)) {
    gene_partition <- ClusterGenes(dist)
  }
  if (is.null(gene_embedding)) {
    df <- ObtainGeneTSNE(dist, perplexity = perplexity)
  } else {
    df <- gene_embedding
  }
  
  if (is.character(gene_partition) | is.logical(gene_partition)) {
    gene_partition <- as.factor(gene_partition)
  }
  df[names(gene_partition), "group"] <- gene_partition
  
  if (is.factor(gene_partition)) {
    if (is.null(module_color)) {
      cols_all <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(nlevels(gene_partition) + 5)
      cols_all <- cols_all[cols_all != "#999999"]
      module_color <- setNames(cols_all[seq_len(nlevels(gene_partition))],
                               levels(gene_partition))
    }
    centroids <- df %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(tSNE_1 = mean(tSNE_1), tSNE_2 = mean(tSNE_2))
    
    p <- ggplot2::ggplot(
      df[order(df$group, na.last = FALSE), ],
      ggplot2::aes(x = tSNE_1, y = tSNE_2, color = group)
    ) +
      ggplot2::geom_point(size = 1) +
      ggplot2::scale_color_manual(values = module_color, na.value = "lightgrey") +
      ggplot2::theme_minimal() +
      ggplot2::ggtitle("TSNE Visualization of Genes") +
      ggplot2::labs(color = "Modules")
    
    if (text) {
      p <- p + ggplot2::geom_text(
        data = centroids,
        ggplot2::aes(label = group), vjust = -1, color = "black"
      )
    }
  } else {
    p <- ggplot2::ggplot(
      df[order(abs(df$"group"), na.last = FALSE), ],
      ggplot2::aes(x = tSNE_1, y = tSNE_2, color = group)
    ) +
      ggplot2::geom_point(size = 1) +
      ggplot2::scale_color_gradient2(low = "blue", mid = "white", high = "red",
                                     na.value = "lightgrey") +
      ggplot2::theme_minimal() +
      ggplot2::ggtitle("TSNE Visualization of Genes") +
      ggplot2::theme(plot.background = ggplot2::element_rect(fill = "lightgrey", color = NA))
  }
  
  return(p)
}