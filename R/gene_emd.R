#' Select common variable genes across two Seurat objects
#'
#' Identifies the union of top variable genes in two Seurat objects and
#' then restricts to genes expressed in both samples. Genes are filtered
#' to those with an expression percentage between 0.5% and 50% before
#' calling `FindVariableFeatures`.
#'
#' @param srat1 A `Seurat` object (sample 1).
#' @param srat2 A `Seurat` object (sample 2).
#' @param ngenes Integer. Number of top variable features to select per sample.
#'   Default is `500`.
#'
#' @return A character vector of gene names present and expressed in both
#'   samples.
#'
#' @export
#'
#' @importFrom Seurat DefaultAssay GetAssayData VariableFeatures NormalizeData FindVariableFeatures
#' @importFrom dplyr %>%
#'
#' @examples
#' \dontrun{
#' common_genes <- SelectCommonGenes(srat1, srat2, ngenes = 500)
#' }
SelectCommonGenes <- function(srat1, srat2, ngenes = 500){
  data_S_ls = list(srat1,srat2)
  assay_name = "RNA"; layer = "data"
  common_genes = Reduce(union,lapply(data_S_ls,function(data_S){
    DefaultAssay(data_S) <- assay_name
    expr_percent <- apply(as.matrix(GetAssayData(data_S,assay = assay_name,layer = layer)) > 0, 1, sum)/ncol(data_S)
    genes <- names(expr_percent)[which(expr_percent > 0.005 & expr_percent < 0.5)]
    VariableFeatures(data_S) = VariableFeatures(subset(data_S,features = genes) %>% NormalizeData() %>% FindVariableFeatures())
    genes <- data_S[[assay_name]]@var.features[1:ngenes]
    genes
  }))
  # Make sure genes express in all samples
  common_genes = Reduce(intersect,lapply(data_S_ls,function(data_S){
    selected_genes = common_genes[common_genes%in%rownames(data_S)]
    expr_percent <- apply(as.matrix(data_S[[assay_name]]@data[selected_genes, ]) > 0, 1, sum)/ncol(data_S)
    selected_genes[expr_percent!=0]
  }))
  return(common_genes)
}


#' Compute gene Earth Mover's Distances via a Python backend
#'
#' Runs diffusion maps and coarse-graining on a Seurat object, writes
#' intermediate files to disk, and launches a Python script
#' (`gene_distance_cal_parallel.py`) in the background to compute pairwise
#' gene EMDs via optimal transport. Results are saved to `dir.path` and can
#' be loaded with [LoadGeneEMD()].
#'
#' @param srat A `Seurat` object. Must have a normalised `RNA` assay.
#' @param common_genes Character vector of gene names to include in the
#'   computation (e.g., the output of [SelectCommonGenes()]).
#' @param dir.path Character string. Path to the directory where intermediate
#'   and output files will be written.
#' @param script_dir Character string or NULL. Path to the directory containing
#'   `gene_distance_cal_parallel.py`. If NULL (default), uses the bundled
#'   Python script included with the package.
#' @param npc Integer. Number of diffusion map components to use. Default `10`.
#' @param K Integer. Number of nearest neighbours for the cell graph. Default `10`.
#' @param reduction Character string. Name of the dimensionality reduction to
#'   use when computing graph distances. Default `"dm"`.
#'
#' @return Invisibly returns `NULL`. Results are written to `dir.path`.
#'
#' @export
#'
#' @importFrom GeneTrajectory RunDM GetGraphDistance CoarseGrain
#' @importFrom Matrix writeMM Matrix
#'
#' @examples
#' \dontrun{
#' ComputeGeneEMD(srat, common_genes, dir.path = "results/emd/")
#' }
ComputeGeneEMD <- function(srat, common_genes, dir.path, script_dir = NULL,
                           npc = 10, K = 10, reduction = "dm", seed = 233) {
  set.seed(seed)
  genes <- common_genes
  srat  <- GeneTrajectory::RunDM(srat, dims = 1:npc, K = K)

  cell.graph.dist <- tryCatch({
    res <- GeneTrajectory::GetGraphDistance(srat, K = K,
                                            reduction = reduction,
                                            dims = 1:npc)
    message("Graph is connected. Proceeding with the rest of the script.")
    res
  }, error = function(e) {
    message("Error detected: ", e$message)
    NULL
  })

  if (is.null(cell.graph.dist)) {
    return(invisible(NULL))
  }

  cg_output <- GeneTrajectory::CoarseGrain(srat, cell.graph.dist, genes,
                                            N = min(ncol(srat) - 1, 1000),
                                            reduction = reduction,
                                            dims = 1:npc)

  dir.path <- file.path(dir.path, "")
  if (!dir.exists(dir.path)) {
    dir.create(dir.path, recursive = TRUE)
  }

  write.table(cg_output[["features"]], file.path(dir.path,"gene_names.csv"), row.names = F, col.names = F, sep = ",")
  write.table(cg_output[["graph.dist"]], file.path(dir.path,"ot_cost.csv"), row.names = F, col.names = F, sep = ",")
  Matrix::writeMM(Matrix::Matrix(cg_output[["gene.expression"]], sparse = T), file.path(dir.path,"gene_expression.mtx"))

  if (is.null(script_dir)) {
    py_script <- system.file("python", "gene_distance_cal_parallel.py", 
                             package = "DiCoLo")
    if (py_script == "") {
      stop("Python script not found. Please check package installation.")
    }
  } else{
    py_script <- file.path(script_dir, "gene_distance_cal_parallel.py")
  }
  
  py <- Sys.which("python3")
  if (py == "") py <- Sys.which("python")
  if (py == "") stop("No Python found. Please install Python and ensure it is on your PATH.")

  system(sprintf("nohup %s %s %s &", py, py_script, dir.path))
  message(sprintf("OT running in background. Load results from:\n%s", dir.path))

  return(invisible(NULL))
}


#' Load a gene EMD distance matrix from disk
#'
#' Reads the pairwise gene Earth Mover's Distance matrix produced by the
#' Python optimal-transport backend (see [ComputeGeneEMD()]) and deduplicates
#' row/column names.
#'
#' @param folder.path Character string. Path to the directory containing the
#'   `emd.csv` output file written by the Python backend.
#'
#' @return A square numeric matrix of gene–gene EMDs with gene names as
#'   row and column names, or `NULL` if the file cannot be loaded.
#'
#' @export
#'
#' @importFrom GeneTrajectory LoadGeneDistMat
#'
#' @examples
#' \dontrun{
#' gene.dist.mat <- LoadGeneEMD("results/emd/")
#' }
LoadGeneEMD <- function(folder.path) {
  gene.dist.mat <- tryCatch(
    {
      GeneTrajectory::LoadGeneDistMat(folder.path, file_name = "emd.csv")
    },
    error = function(e) {
      message(e$message)
      return(NULL)
    }
  )

  if (is.null(gene.dist.mat)) {
    return(NULL)
  }

  gene.dist.mat <- gene.dist.mat[
    unique(rownames(gene.dist.mat)),
    unique(colnames(gene.dist.mat))
  ]

  return(gene.dist.mat)
}


