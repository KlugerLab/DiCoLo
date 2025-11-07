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
load_or_install("GeneTrajectory", "KlugerLab/GeneTrajectory")
load_or_install("dplyr")
load_or_install("locfdr")
load_or_install("RColorBrewer")

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

ComputeGeneEMD <- function(srat, common_genes, dir.path){
  genes = common_genes
  npc = 10
  srat <- GeneTrajectory::RunDM(srat, dims = 1:npc)
  # check if the cell graph is connected
  cell.graph.dist <- tryCatch({
    res = GeneTrajectory::GetGraphDistance(srat, K = 10, dims = 1:npc)
    message("Graph is connected. Proceeding with the rest of the script.")
    res
  }, error = function(e) {
    message("Error detected: ", e$message)
    NULL
  })
  if(is.null(cell.graph.dist)){return(NULL)}
  cg_output <- GeneTrajectory::CoarseGrain(srat, cell.graph.dist, genes, N = min(ncol(srat)-1,1000), dims = 1:npc)


  dir.path = file.path(dir.path,"")
  if(!dir.exists(dir.path)){
    dir.create(dir.path, recursive = TRUE)
  }
  write.table(cg_output[["features"]], file.path(dir.path,"gene_names.csv"), row.names = F, col.names = F, sep = ",")
  write.table(cg_output[["graph.dist"]], file.path(dir.path,"ot_cost.csv"), row.names = F, col.names = F, sep = ",")
  Matrix::writeMM(Matrix::Matrix(cg_output[["gene.expression"]], sparse = T), file.path(dir.path,"gene_expression.mtx"))

  script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
  py_script = file.path(script_dir,"gene_distance_cal_parallel.py")
  
  py <- Sys.which("python3")
  if (py == "") py <- Sys.which("python")
  if (py == "") stop("No Python found. Please install Python and ensure it is on your PATH.")
  
  system(sprintf("nohup %s %s %s &", py, py_script, dir.path))
  message(sprintf("OT running backend, please load result from \n%s",dir.path))
  return(NULL)
}

LoadGeneEMD <- function(folder.path){
  gene.dist.mat <- tryCatch(
    {
      LoadGeneDistMat(folder.path, file_name = "emd.csv")
    },
    error = function(e) {
      message(e$message)
      return(NULL)
    }
  )
  if (is.null(gene.dist.mat)) {
    return(NULL)
  }
  gene.dist.mat <- gene.dist.mat[unique(rownames(gene.dist.mat)),
                                 unique(colnames(gene.dist.mat))]
  return(gene.dist.mat)
}

inv_sqrt <- function(L_sym){
  E.list <- rARPACK::eigs_sym(L_sym, k = nrow(L_sym))
  L_sym_inv_sqrt <- E.list$vectors %*% diag(1/sqrt(E.list$values)) %*% t(E.list$vectors)
  rownames(L_sym_inv_sqrt) = rownames(L_sym)
  colnames(L_sym_inv_sqrt) = colnames(L_sym)
  L_sym_inv_sqrt
}

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

ComputeGraphOperator = function(gene.dist.mat, knn = 10, shrink_thred = 1e-10){
  # Normalize graph laplacian
  gene_graph = affinity.mtx(gene.dist.mat, K = knn)
  
  E.list <- rARPACK::eigs_sym(gene_graph, k = nrow(gene_graph))
  # cat("before_shrink: ",min(E.list$values),"\n")
  if(min(E.list$values) > 0){return(gene_graph)}
  
  gene_graph_shrink <- E.list$vectors %*% diag(pmax(E.list$values,shrink_thred)) %*% t(E.list$vectors)
  rownames(gene_graph_shrink) = colnames(gene_graph_shrink) = rownames(gene_graph)
  return(gene_graph_shrink)
}

ComputeDifferentialOperator = function(gene_graph_comp, gene_graph_proj, 
                    regulator_min = 1e-2){
  # estimate regulator
  E.list <- RunSVD(gene_graph_proj)
  regulator <- E.list$values[FindKneePoint(E.list$values, plot.fig = FALSE)]
  regulator = max(regulator,regulator_min)
  # message(sprintf("Regulator: %s",regulator))
  
  P_inv <- inv_sqrt(gene_graph_proj + regulator * diag(nrow(gene_graph_proj)))
  P_diff <- P_inv %*% gene_graph_comp %*% P_inv
  
  return(P_diff)
}

RunSVD = function(K, eig_keep = nrow(K)){
  E.list <- rARPACK::eigs_sym(K, k = eig_keep, which = "LM")
  rownames(E.list$vectors) = rownames(K)
  return(E.list)
}

FindKneePoint = function(eigval, plot.fig = FALSE){
  Vecs = eigval
  curve_data = as.matrix(data.frame(x = 1:length(Vecs), y = Vecs))
  line_start <- curve_data[1, ]
  line_end <- curve_data[nrow(curve_data), ]
  dx <- line_end["x"] - line_start["x"]
  dy <- line_end["y"] - line_start["y"]
  norm_factor <- sqrt(dx^2 + dy^2)
  A <- dy
  B <- -dx
  C <- dx * line_start["y"] - dy * line_start["x"]
  distances <- abs(A * curve_data[, "x"] + B * curve_data[, "y"] + C)/norm_factor
  knee_point = which.max(distances)
  if (plot.fig == TRUE) {
    plot(Vecs, pch = 20, xlab = "Index", ylab = "EigVal",
         cex.lab = 1.4,   # axis label size
         cex.axis = 1.2,  # tick label size
         cex.main = 1.5,
         cex = 1)
    abline(v = knee_point, col = "grey", lwd = 2, lty = 'dashed')
    points(1:(knee_point-1), Vecs[1:(knee_point-1)], col = "red", pch = 20,cex = 1)
  }
  return(knee_point)
}

ObtainGeneTSNE = function(dist, perplexity = 30){
  set.seed(233)
  df = data.frame(Rtsne::Rtsne(dist, is_distance = TRUE, dim = 2, 
                               perplexity = perplexity, max.iter = 500)$Y)
  colnames(df) = c("tSNE_1", "tSNE_2")
  rownames(df) = labels(dist)
  return(df)
}

VisualizeGeneTSNE = function(dist = NULL, gene_embedding = NULL, gene_partition = NULL, 
                             filtered = TRUE, perplexity = 30, text = TRUE, module_color = NULL) 
{
  if (is.null(gene_partition)) {
    gene_partition = ClusterGenes(dist, filtered = filtered)
  }
  if(is.null(gene_embedding)){
    df = ObtainGeneTSNE(dist, perplexity = perplexity)
  }else{
    df = gene_embedding
  }
  if(is.character(gene_partition) | is.logical(gene_partition)){
    gene_partition = as.factor(gene_partition)
  }
  df[names(gene_partition), "group"] = gene_partition
  if(is.factor(gene_partition)){
    if(is.null(module_color)){
      cols_all = colorRampPalette(brewer.pal(9, "Set1"))(nlevels(gene_partition)+5)
      cols_all <- cols_all[cols_all != "#999999"]
      module_color = setNames(cols_all[1:nlevels(gene_partition)], 
                              levels(gene_partition))
    }
    
    centroids <- df %>% group_by(group) %>% summarise(tSNE_1 = mean(tSNE_1), 
                                                      tSNE_2 = mean(tSNE_2))
    p = ggplot(df[order(df$group,na.last = FALSE),], aes(x = tSNE_1, y = tSNE_2, color = group)) + 
      geom_point(size = 1) + scale_color_manual(values = module_color,
                                                na.value = "lightgrey") + 
      theme_minimal() + ggtitle("TSNE Visualization of Genes") + 
      labs(color = "Modules")
    if(text){
      p = p + geom_text(data = centroids, aes(label = group), vjust = -1, 
                        color = "black")
    }
  }else{
    p = ggplot(df[order(abs(df$'group'), na.last = FALSE),], aes(x = tSNE_1, y = tSNE_2, color = group)) + 
      geom_point(size = 1) + scale_color_gradient2(low = "blue", mid = "white", high = "red", na.value = "lightgrey") + theme_minimal() + 
      ggtitle("TSNE Visualization of Genes") + 
      theme(plot.background = element_rect(fill = "lightgrey", color = NA))
  }
  return(p)
}

SelectSignificantGenes <- function(V, lfdr_thresh = 0.2,
                                  min_genes = 3,
                                  max_genes = 50) {
  stopifnot(is.matrix(V))
  if (is.null(colnames(V))) colnames(V) <- paste0("EV", seq_len(ncol(V)))
  # modules <- list()
  eigenvec_indicator = matrix(0,nrow(V),ncol(V))
  rownames(eigenvec_indicator) = rownames(V); colnames(eigenvec_indicator) = colnames(V)
  
  for (j in seq_len(ncol(V))) {
    v <- V[, j]
    names(v) <- rownames(V)
    
    ## z-transform (MAD robust)
    med <- median(v, na.rm = TRUE)
    s   <- mad(v, center = med, constant = 1.4826, na.rm = TRUE)
    if (s < 1e-12) s <- sd(v)  # fallback if too flat
    z <- (v - med) / s
    
    ## run local FDR
    lf <- suppressWarnings(locfdr(z, nulltype = 0, plot = 0))
    lfdr <- lf$fdr
    names(lfdr) <- names(v)
    
    ## select by lfdr
    lfdr_thresh = min(lfdr_thresh,lf$Efdr[1])
    # cat(lfdr_thresh,"\n")
    sel <- names(lfdr)[lfdr < lfdr_thresh]
    if (length(sel) == 0) next
    
    pos <- sel[z[sel] > 0]
    neg <- sel[z[sel] < 0]
    
    ## keep most extreme if too many
    if (length(pos) > max_genes) {
      pos <- pos[order(z[pos], decreasing = TRUE)[1:max_genes]]
    }
    if (length(neg) > max_genes) {
      neg <- neg[order(-z[neg], decreasing = TRUE)[1:max_genes]]
    }
    
    if (length(pos) > min_genes){
      # modules[[paste0(colnames(V)[j], "+")]] <- pos
      eigenvec_indicator[pos,j] = 1
    }
    if (length(neg) > min_genes){
      # modules[[paste0(colnames(V)[j], "-")]] <- neg
      eigenvec_indicator[neg,j] = 1
    }
  }
  return(eigenvec_indicator)
  # return(list(modules = modules,
  #             gene_indicator = eigenvec_indicator))
}


tighten_modules <- function(gene_partition, distM, z_thresh = 3) {
  unique_clusters = levels(gene_partition)
  new_partition = gene_partition
  for (cl in unique_clusters) {
    members <- names(gene_partition[gene_partition == cl])
    if (length(members) <= 2) next
    
    sub_dist <- distM[members, members]
    avg_d <- apply(sub_dist, 1, function(x) mean(x[-which.min(x == 0)]))
    mu <- median(avg_d)
    sigma <- mad(avg_d)  # robust spread
    outliers <- names(avg_d[avg_d > mu + z_thresh * sigma])
    
    new_partition[outliers] <- NA  # mark as unassigned
  }
  return(new_partition)
}


ClusterGenes <- function (gene_dist, clustering_method = "average", 
                                   min_gene = 2, deepSplit = 1, z_thresh = 3) 
{
  gene_hree = hclust(gene_dist, method = clustering_method)
  gene_partition = dynamicTreeCut::cutreeDynamic(dendro = gene_hree, 
                                                 distM = as.matrix(gene_dist), 
                                                 deepSplit = deepSplit, 
                                                 pamStage = TRUE, 
                                                 pamRespectsDendro = TRUE, 
                                                 minClusterSize = min_gene)
  names(gene_partition) = labels(gene_dist)
  gene_partition = as.factor(gene_partition)
  gene_partition = tighten_modules(gene_partition, as.matrix(gene_dist), z_thresh)
  
  gene_partition = factor(gene_partition, levels = seq(nlevels(gene_partition)))
  # filter-out singleton
  gene_partition[gene_partition %in% names(table(gene_partition))[table(gene_partition) < min_gene]] = NA
  gene_partition = gene_partition[!is.na(gene_partition)]
  gene_partition <- droplevels(gene_partition)
  levels(gene_partition) = 1:nlevels(gene_partition)
  return(gene_partition)
}

RankGeneModules <- function(gene_partition,gene_indicator,eigen_values){
  eigen_values = eigen_values[1:ncol(gene_indicator)]
  eigen_vec_weight = eigen_values/sum(eigen_values); names(eigen_vec_weight) = colnames(gene_indicator)
  # ensure all are named consistently
  common_genes <- intersect(names(gene_partition), rownames(gene_indicator))
  gene_partition <- gene_partition[common_genes]
  gene_indicator <- gene_indicator[common_genes,]
  gene_indicator <- apply(gene_indicator,1,function(x) colnames(gene_indicator)[min(which(x == 1))])
  
  # for each module m, compute c_{m,k} = fraction of its genes from eigenvector k
  modules <- split(names(gene_partition), gene_partition)
  module_weights <- sapply(modules, function(genes) {
    eig_prop <- table(gene_indicator[genes]) / length(genes)
    # multiply proportions by eigenvector weights (lambda_k)
    
    w <- sum(eig_prop * eigen_vec_weight[names(eig_prop)])
    return(w)
  })
  module_order = names(sort(module_weights,decreasing = T))
  name_map = setNames(1:length(module_order),module_order)
  gene_partition = setNames(name_map[as.character(gene_partition)],names(gene_partition))
  gene_partition = as.factor(gene_partition)
  return(gene_partition)
}
