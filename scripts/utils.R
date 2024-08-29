# required packages ####
library(Seurat)
library(ggplot2)
library(ggpubr)
library(ggbreak)
library(ggfortify)
library(mosaic)
library(RColorBrewer)
library(dplyr)
library(EnhancedVolcano)
library(fgsea)
library(qusage)
library(data.table)
library(pheatmap)
library(tidyverse)
library(clustree)
library(celldex)
library(SingleR)
library(DoubletFinder)
library(harmony)
library(anndata)
library(chisq.posthoc.test)
library(presto)
library(CellChat)
library(patchwork)
library(NMF)
library(ggalluvial)
library(reticulate)
library(ComplexHeatmap)
library(monocle3)
library(SeuratWrappers)
library(magrittr)
library(infercnv)
library(sctransform)
library(glmGamPoi)
library(scCustomize)
library(parallel)
library(BiocParallel)
library(fpc)
library(genomicInstability)
library(ineq)
library(plyr)
library(cowplot)
library(sceasy)
library(dendsort)
library(lmerTest)
library(emmeans)
library(rstatix)
library(expm)
# database of color vectors ####
treatment_colors = c("dodgerblue2", "orange", "grey80")
patient_colors = c("#E41A1C", "#2840ad", "#4DAF4A")
group_colors = c("#f2797b", "#E41A1C", "#610b0c", "#4869f0", "#2847c9", "#2840ad","#122680", "#031252","#a1f09e", "#73d66f", "#4DAF4A", "#378c34", "#246b21", "#134511")

broadcelltype_colors_alpha = c("#f032e6", "#000075", "#f58231",  "#aaffc3", "#e6194B", "#bfef45",
                               "#469990", "#42d4f4", "#911eb4","#3cb44b", "#ffe119", "#dcbeff")
broadcelltype_order_alpha = c("B", "Basal_Tumor", "CAF", "Cycling_Epithelial", "Endothelial", "LEC", 
                              "Luminal_Progenitor", "Luminal_Tumor", "Macrophage", "MEC", "PVL", "T")

broadcelltype_colors_rainbow <- c("#e6194B", "#f58231", "#ffe119", "#bfef45", "#3cb44b", "#aaffc3", "#469990", "#42d4f4", "#000075", "#911eb4", "#f032e6",  "#dcbeff")
broadcelltype_order_rainbow = c("Endothelial", "CAF", "PVL", "LEC", "MEC", "Cycling_Epithelial", "Luminal_Progenitor", "Luminal_Tumor", "Basal_Tumor", "Macrophage", "B", "T")

# all fine tune together
all_fine_anno_colors_rainbow = c("#e6194B", "#ffac63", "#f58231", "#fff27b", "#ccad00", "#bfef45", "#3cb44b", 
                                 "#aaffc3", "#469990", "#2C8DC9", "#42d4f4", "#afebfa", "#86b1f7", "#4363d8", 
                                 "#000075", "#7979d9", "#7740c9", "#3c1361", "#dcbeff", "#c68ff7", 
                                 "#b380ff", "#911eb4", "#f032e6", "#f5a4f1", "#f0e8ff")
all_fine_anno_order_rainbow = c("Endothelial", "iCAF", "myCAF", "imPVL", "dPVL", "LEC", "MEC",
                                "Cycling_Epithelial", "Luminal_Progenitor", "Luminal_Tumor2", "Luminal_Tumor3", "Luminal_Tumor1", "Luminal_Tumor5", "Luminal_Tumor4", 
                                "Basal_Tumor1", "Basal_Tumor2", "Teff", "ChopT", "Tcm", "Tem",
                                "MAIT", "Macrophage", "Plasma", "B", "Th1")
all_fine_anno_short_order_rainbow = c("Endothelial", "iCAF", "myCAF", "imPVL", "dPVL", "LEC", "MEC",
                                      "CYC", "LP", "LT2", "LT3", "LT1", "LT5", "LT4", 
                                      "BT1", "BT2", "Teff", "ChopT", "Tcm", "Tem",
                                      "MAIT", "Macrophage", "Plasma", "B", "Th1")

all_fine_anno_colors_alpha = c("#f032e6", "#000075","#7979d9", "#3c1361", "#aaffc3","#ccad00",
                               "#e6194B", "#ffac63", "#fff27b", "#bfef45", "#469990", "#afebfa", 
                               "#2C8DC9", "#42d4f4", "#4363d8", "#86b1f7", "#911eb4", "#b380ff", 
                               "#3cb44b", "#f58231", "#f5a4f1", "#dcbeff", "#7740c9", "#c68ff7", "#f0e8ff")

all_fine_anno_order_alpha = c("B", "Basal_Tumor1", "Basal_Tumor2", "ChopT", "Cycling_Epithelial", "dPVL",
                              "Endothelial", "iCAF", "imPVL", "LEC", "Luminal_Progenitor", "Luminal_Tumor1", 
                              "Luminal_Tumor2", "Luminal_Tumor3", "Luminal_Tumor4", "Luminal_Tumor5", "Macrophage", "MAIT", 
                              "MEC", "myCAF", "Plasma", "Tcm", "Teff", "Tem", "Th1")

# epithelial specific colors
epi_fine_anno_colors_rainbow = c("#bfef45", "#3cb44b", "#469990", "#aaffc3", "#2C8DC9", "#42d4f4", "#afebfa", "#86b1f7", "#4363d8", "#000075", "#7979d9")
epi_fine_anno_order_rainbow = c("LEC", "MEC", "Luminal_Progenitor", "Cycling_Epithelial", "Luminal_Tumor2", "Luminal_Tumor3", "Luminal_Tumor1", "Luminal_Tumor5", "Luminal_Tumor4", "Basal_Tumor1", "Basal_Tumor2")

epi_fine_anno_colors_alpha = c("#000075","#7979d9", "#aaffc3", "#bfef45", "#469990", "#afebfa", "#2C8DC9", "#42d4f4", "#4363d8", "#86b1f7", "#3cb44b")
epi_fine_anno_order_alpha = c("Basal_Tumor1", "Basal_Tumor2", "Cycling_Epithelial", "LEC", "Luminal_Progenitor", "Luminal_Tumor1", "Luminal_Tumor2", "Luminal_Tumor3", "Luminal_Tumor4", "Luminal_Tumor5", "MEC")

epi_fine_anno_short_colors_rainbow = c("#bfef45", "#3cb44b", "#469990", "#aaffc3", "#2C8DC9", "#42d4f4", "#afebfa", "#86b1f7", "#4363d8", "#000075", "#7979d9")
epi_fine_anno_short_order_rainbow = c("LEC", "MEC", "LP", "CYC", "LT2", "LT3", "LT1", "LT5", "LT4", "BT1", "BT2")

epi_fine_anno_short_colors_alpha = c("#000075","#7979d9", "#aaffc3", "#bfef45", "#469990", "#afebfa", "#2C8DC9", "#42d4f4", "#4363d8", "#86b1f7", "#3cb44b")
epi_fine_anno_short_order_alpha = c("BT1", "BT2", "CYC", "LEC", "LP", "LT1", "LT2", "LT3", "LT4", "LT5", "MEC")

# stromal specific colors
stro_fine_anno_colors_rainbow = c("#e6194B", "#ffac63", "#f58231", "#fff27b", "#ccad00")
stro_fine_anno_order_rainbow = c("Endothelial", "iCAF", "myCAF", "imPVL", "dPVL")

stro_fine_anno_colors_alpha = c("#ccad00", "#e6194B", "#ffac63", "#fff27b", "#f58231")
stro_fine_anno_order_alpha = c("dPVL", "Endothelial",  "iCAF", "imPVL", "myCAF")

# immune specific colors
imm_fine_anno_colors_rainbow = c("#dcbeff", "#c68ff7", "#7740c9", "#f0e8ff", "#3c1361", "#b380ff", "#f032e6", "#f5a4f1", "#911eb4")
imm_fine_anno_order_rainbow = c("Tcm", "Tem", "Teff", "Th1", "ChopT", "MAIT", "B", "Plasma", "Macrophage")

imm_fine_anno_colors_alpha = c("#f032e6", "#3c1361", "#911eb4", "#b380ff", "#f5a4f1", "#dcbeff", "#7740c9", "#c68ff7", "#f0e8ff")
imm_fine_anno_order_alpha = c("B", "ChopT", "Macrophage", "MAIT", "Plasma", "Tcm", "Teff", "Tem", "Th1")

# preprocessing functions ####
create_seurat_object <- function(data_dir, project_name) {
  data <- Seurat::Read10X(data.dir = data_dir)
  seurat_obj <- CreateSeuratObject(counts = data, project = project_name, min.cells = 3, min.features = 100)
  return(seurat_obj)
}
process_and_doublet_finder <- function(obj, nExp_multiplier) {
  # objects need to be processed through standard workflow before calculating doublets
  obj_dd <- subset(obj, subset = nFeature_RNA > 50 & nFeature_RNA < 10000)
  obj_dd <- NormalizeData(obj_dd)
  obj_dd <- FindVariableFeatures(obj_dd, selection.method = "vst", nfeatures = 2000)
  obj_dd <- ScaleData(obj_dd)
  obj_dd <- RunPCA(obj_dd, npcs = 25)
  obj_dd <- FindNeighbors(obj_dd, dims = 1:15)
  obj_dd <- FindClusters(obj_dd, resolution = 0.5)
  obj_dd <- RunUMAP(obj_dd, dims = 1:15)
  
  # calculate pK
  sweep.res.list <- paramSweep(obj_dd, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  # determine the optimal pK value
  optimal_pK <- bcmvn[which.max(bcmvn$BCmetric), "pK"]
  print(paste("Optimal pK value for current object:", optimal_pK))
  
  # calculate homotypic proportion and expected number of doublets
  annotations <- obj_dd@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(nExp_multiplier * nrow(obj_dd@meta.data))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  
  # apply DoubletFinder
  obj_dd <- doubletFinder(obj_dd, 
                          PCs = 1:15, 
                          pN = 0.25, 
                          pK = as.numeric(optimal_pK),
                          nExp = nExp_poi.adj, 
                          reuse.pANN = FALSE, 
                          sct = FALSE)
  
  return(obj_dd)
}
doublet_dimplot_check <- function(doublet_detected_objs) {
  plot_list <- list()
  
  for (obj_name in names(doublet_detected_objs)) {
    obj <- doublet_detected_objs[[obj_name]]
    
    # find the doublet column
    group_by_column <- grep("^DF.classifications_", colnames(obj@meta.data), value = TRUE)
    
    # generate DimPlot if the column is found
    if (length(group_by_column) > 0) {
      plot <- scCustomize::DimPlot_scCustom(obj, reduction = "umap", group.by = group_by_column[1], colors_use = c("red", "grey80"))
      plot_list[[obj_name]] <- plot
    } else {
      stop("No column starting with 'DF.classifications_' found in object: ", obj_name)
    }
  }
  
  # combine all plots using cowplot
  combined_plot <- plot_grid(plotlist = plot_list, ncol = 3)
  
  return(combined_plot)
}
print_metadata_table <- function(seurat_obj, column_prefix) {
  # find the column name that starts with the given prefix
  column_name <- grep(column_prefix, colnames(seurat_obj@meta.data), value = TRUE)
  
  if (length(column_name) > 0) {
    for (col in column_name) {
      cat("Table for column:", col, "\n")
      print(table(seurat_obj@meta.data[[col]]))
      cat("\n")
    }
  } else {
    print(paste("No columns starting with", column_prefix, "found in metadata"))
  }
}
update_doublets <- function(raw, doublet) {
  for (obj_name in names(raw)) {
    raw_obj <- raw[[obj_name]]$object
    doublet_obj <- doublet[[obj_name]]
    
    if (!is.null(raw_obj) && !is.null(doublet_obj)) {
      # find the correct DF.classifications column in doublet_obj
      classification_col <- grep("^DF.classifications", colnames(doublet_obj@meta.data), value = TRUE)
      
      if (length(classification_col) > 0) {
        classification_col <- classification_col[1]
        
        # subset the doublet detected metadata
        doublet_cells <- rownames(subset(doublet_obj@meta.data, doublet_obj@meta.data[[classification_col]] == "Doublet"))
        
        # update the doublet column in the raw object metadata
        raw_obj@meta.data$doublet = paste("single")
        raw_obj@meta.data$doublet[rownames(raw_obj@meta.data) %in% doublet_cells] <- "doublet"
        
        # assign the updated raw object back to the list
        raw_objs[[obj_name]]$object <- raw_obj
      }
    }
  }
  
  return(raw_objs)
}
qc_metrics <- function(obj){
  obj <- scCustomize::Add_Mito_Ribo(object = obj, species = "Human")
  obj <- scCustomize::Add_Cell_Complexity(object = obj)
  
  return(obj)
}
qc_plots <- function(obj, pdo_name) {
  # create plots
  p1 <- scCustomize::QC_Plots_Genes(seurat_object = obj, low_cutoff = 100, high_cutoff = 8500)
  p2 <- scCustomize::QC_Plots_UMIs(seurat_object = obj, low_cutoff = 500, high_cutoff = 50000)
  p3 <- scCustomize::QC_Plots_Mito(seurat_object = obj, high_cutoff = 20)
  p4 <- scCustomize::QC_Plots_Complexity(seurat_object = obj, high_cutoff = 0.75)
  
  qc_plots <- wrap_plots(p1, p2, p3, p4, ncol = 4)
  
  p1 <- scCustomize::QC_Plot_UMIvsGene(seurat_object = obj, 
                                       low_cutoff_gene = 100, high_cutoff_gene = 8500, 
                                       low_cutoff_UMI = 500, high_cutoff_UMI = 50000)
  p2 <- scCustomize::QC_Plot_GenevsFeature(seurat_object = obj, feature1 = "percent_mito", 
                                           low_cutoff_gene = 100, high_cutoff_gene = 8500, 
                                           high_cutoff_feature = 20)
  
  umi_gene_plots <- wrap_plots(p1, p2)
  
  # combine all plots 
  combined_plots <- wrap_plots(qc_plots, umi_gene_plots)
  
  # display plots in the RStudio plots window
  print(combined_plots + plot_annotation(title = paste("QC Plots:", pdo_name)), ncol = 1)
}
qc_filtering <- function(obj, obj_name) {
  # subset by filters determined
  obj <- subset(obj, subset = nFeature_RNA > 100 & nFeature_RNA < 9000 & percent_mito < 20 & doublet=="single")
  
  # print the number of rows in the meta.data slot
  cat("Number of cells in", obj_name, ":", nrow(obj@meta.data), "\n")
  
  # return the subsetted object
  return(obj)
}

# annotation functions ####
add_singleR_metadata <- function(seurat_obj, predicted_celltypes) {
  # add barcode prefix to the metadata
  seurat_obj@meta.data$barcode <- rownames(seurat_obj@meta.data)
  
  # prepare the labels data frame
  to_merge_swarbrick_labels <- dplyr::select(as.data.frame(predicted_celltypes), X, labels, pruned.labels)
  colnames(to_merge_swarbrick_labels)=c("barcode", "singleR_labels", "singleR_pruned.labels")
  
  # merge the metadata with the labels
  meta_data <- dplyr::left_join(seurat_obj@meta.data, to_merge_swarbrick_labels, by = "barcode")
  rownames(meta_data) <- meta_data$barcode
  seurat_obj@meta.data <- meta_data
  
  return(seurat_obj)
}
get_cells_for_highlight <- function(group_name, seurat_obj) {
  cells <- WhichCells(seurat_obj, idents = c(group_name))
  return(cells)
}

# integration functions ####
do_clustree_integrated <- function(obj, res_vector = NULL, res_seq = NULL){
  # set resolutions to use
  if(is.null(res_vector)){
    if(!is.null(res_seq)){
      res_vector = seq(res_seq[1], res_seq[2], res_seq[3])
    } else stop("Please input res_vector or res_seq parameters")
  }
  
  # perform multiple clusterings at resolutions specified
  for (i in res_vector){
    obj <- FindClusters(obj, resolution = i)
  }
  
  # plot clustree graph
  clustree::clustree(obj@meta.data, prefix = "SCT_snn_res.")
}
integration_simpleplots <- function(obj, reduction = NULL, cluster_labels) {
  # generate the individual plots
  plot1 <- DimPlot_scCustom(obj, reduction = reduction, group.by = cluster_labels)+NoLegend()
  plot2 <- DimPlot_scCustom(obj, reduction = reduction, group.by = "patientID", colors_use = patient_colors)
  plot3 <- DimPlot_scCustom(obj, reduction = reduction, group.by = "group")
  plot4 <- DimPlot_scCustom(obj, reduction = reduction, group.by = "treatment", colors_use = treatment_colors)
  plot5 <- DimPlot_scCustom(obj, reduction = reduction, group.by = "singleR_majorlabels")
  # combine the plots into a single cowplot
  combined_plot <- plot_grid(plot1, plot2, plot3, plot4, plot5, ncol = 3)
  
  return(combined_plot)
}
quick_process <- function(seurat_obj){
  seurat_obj <- SCTransform(seurat_obj, vars.to.regress = "percent_mito", verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
  
  return(seurat_obj)
}
count_filtered_cells <- function(markers, log2FC_condition) {
  if (log2FC_condition == "greater") {
    return(nrow(subset(markers, avg_log2FC > 1 & p_val_adj < 0.05)))
  } else if (log2FC_condition == "less") {
    return(nrow(subset(markers, avg_log2FC < 1 & p_val_adj < 0.05)))
  }
}

# CNV functions ####
process_inferCNV <- function(seurat_obj, obj_name, base_path) {
  # select the necessary columns
  cnv_anno <- dplyr::select(seurat_obj[[]], barcode, coarse_anno)
  cnv_counts <- GetAssayData(seurat_obj, layer="counts")
  
  # adjust barcode
  cnv_anno$barcode <- paste0("X", cnv_anno$barcode) # needed because R has an in-built behavior to add an X to the beginning of colnames
  cnv_anno$barcode <- gsub("-1", ".1", cnv_anno$barcode) # needed because R has an in-built behavior to change this when saving colnames
  
  # define file paths
  counts_file <- paste0(base_path, "/", obj_name, "_counts.txt")
  anno_file <- paste0(base_path, "/", obj_name, "_anno.txt")
  
  # write files
  write.table(cnv_counts, file = counts_file,
              row.names = TRUE, col.names = TRUE, sep = "\t")
  write.table(cnv_anno, file = anno_file,
              row.names = FALSE, col.names = FALSE, sep = "\t")
  
  # create infercnvObject
  infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = counts_file,
    annotations_file = anno_file,
    delim = "\t",
    gene_order_file = paste0(base_path, "/hg38_gencode_v27.txt"),
    ref_group_names = c("Immune")
  )
  
  # run inferCNV
  infercnv_obj_hmm <- infercnv::run(
    infercnv_obj,
    cutoff = 0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
    out_dir = paste0(base_path, "/run_", obj_name, "_hmm/"),
    cluster_by_groups = FALSE,
    denoise = TRUE,
    HMM = TRUE
  )
  
  return(infercnv_obj_hmm)
}
scale_between_minus1_and_1 <- function(x) {
  min_x <- min(x)
  max_x <- max(x)
  scaled_x <- 2 * (x - min_x) / (max_x - min_x) - 1
  return(scaled_x)
}
calculate_gi_score <- function(cnv_data) {
  scaled_cnv <- scale_between_minus1_and_1(cnv_data)
  gi_scores <- colMeans(scaled_cnv^2)
  names(gi_scores) <- colnames(cnv_data) 
  return(gi_scores)
}
evaluate_thresholds_f1 <- function(score_df, x_multiplier, y_multiplier) {
  # group into each cluster and calculate mean gi and correlation score
  cluster_means <- score_df %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(
      mean_gi_score = mean(gi_score),
      mean_correlation_score = mean(correlation_score)
    ) %>%
    dplyr::arrange(mean_gi_score + mean_correlation_score)
  
  first_cancer_cluster <- cluster_means$cluster[2]
  first_cancer_df <- score_df %>% dplyr::filter(cluster == first_cancer_cluster)
  
  x_int <- mean(first_cancer_df$gi_score) - (x_multiplier * sd(first_cancer_df$gi_score))
  y_int <- mean(first_cancer_df$correlation_score) - (y_multiplier * sd(first_cancer_df$correlation_score))
  
  score_df <- score_df %>%
    dplyr::mutate(
      predicted_label = case_when(
        gi_score < x_int & correlation_score < y_int ~ "normal",
        TRUE ~ "cancer"
      )
    )
  
  # constructing confusion matrix
  confusion <- table(score_df$true_label, score_df$predicted_label)
  
  # check if labels exist in confusion matrix
  # Ensure all categories are represented in the confusion matrix
  if (!"cancer" %in% rownames(confusion)) {confusion <- rbind(confusion, cancer = rep(0, ncol(confusion)))}
  if (!"normal" %in% rownames(confusion)) {confusion <- rbind(confusion, normal = rep(0, ncol(confusion)))}
  if (!"cancer" %in% colnames(confusion)) {confusion <- cbind(confusion, cancer = rep(0, nrow(confusion)))}
  if (!"normal" %in% colnames(confusion)) {confusion <- cbind(confusion, normal = rep(0, nrow(confusion)))}
  
  # calculating f1 score
  precision <- ifelse(sum(confusion[2, ]) > 0, confusion[2, 2] / sum(confusion[2, ]), 0)
  recall <- ifelse(sum(confusion[, 2]) > 0, confusion[2, 2] / sum(confusion[, 2]), 0)
  f1_score <- ifelse(precision + recall > 0, 2 * (precision * recall) / (precision + recall), 0)
  
  return(list(
    x_multiplier = x_multiplier,
    y_multiplier = y_multiplier,
    precision = precision,
    recall = recall,
    f1_score = f1_score,
    x_int = x_int,
    y_int = y_int
  ))
}
evaluate_thresholds_silwidth <- function(score_df, x_multiplier, y_multiplier) {
  cluster_means <- score_df %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(
      mean_gi_score = mean(gi_score),
      mean_correlation_score = mean(correlation_score)
    ) %>%
    dplyr::arrange(mean_gi_score + mean_correlation_score)
  
  first_cancer_cluster <- cluster_means$cluster[2]
  first_cancer_df <- score_df %>% dplyr::filter(cluster == first_cancer_cluster)
  
  x_int <- mean(first_cancer_df$gi_score) - (x_multiplier * sd(first_cancer_df$gi_score))
  y_int <- mean(first_cancer_df$correlation_score) - (y_multiplier * sd(first_cancer_df$correlation_score))
  
  score_df <- score_df %>%
    dplyr::mutate(
      predicted_label = case_when(
        gi_score < x_int & correlation_score < y_int ~ "normal",
        TRUE ~ "cancer"
      )
    )
  
  # Evaluate separation metrics: silhouette score for clusters
  pamk_result <- pamk(score_df[, c("gi_score", "correlation_score")], krange = 2:4)
  silhouette_result <- pam(score_df[, c("gi_score", "correlation_score")], pamk_result$nc)
  
  avg_silhouette_width <- mean(silhouette_result$silinfo$widths[, "sil_width"])
  
  return(list(
    x_multiplier = x_multiplier,
    y_multiplier = y_multiplier,
    avg_silhouette_width = avg_silhouette_width,
    x_int = x_int,
    y_int = y_int
  ))
}
process_cnv_results_f1 <- function(infercnv_obj, base_path, obj_name, x_multipliers = seq(0.5, 2, by = 0.25), y_multipliers = seq(0.5, 2, by = 0.25)) {
  cnv_data <- infercnv_obj@expr.data
  # calculate genomic instability score
  gi_scores <- calculate_gi_score(cnv_data)
  
  # determine top 5% cells for CNV profile
  top_5_percent_threshold <- quantile(gi_scores, 0.95)
  top_5_percent_cells <- names(subset(gi_scores, gi_scores >= top_5_percent_threshold))
  top_5_percent_cells_no_na <- top_5_percent_cells[!is.na(top_5_percent_cells)]
  
  # create average CNV profile across the top 5% most unstable cells
  avg_cnv_profile <- rowMeans(cnv_data[,top_5_percent_cells_no_na], na.rm = TRUE)
  
  # calculate correlation scores
  correlation_scores <- sapply(1:ncol(cnv_data), function(i) {
    cell_cnv_profile <- cnv_data[, i]
    cor(cell_cnv_profile, avg_cnv_profile, use = "complete.obs")
  })
  names(correlation_scores) <- colnames(cnv_data)
  
  # combine scores into a data frame
  score_df <- data.frame(
    cell = colnames(cnv_data),
    gi_score = gi_scores,
    correlation_score = correlation_scores
  )
  
  # initialize true labels
  true_labels <- c()
  
  # extract the subclusters under the 'Immune' category
  immune_subclusters <- infercnv_obj@tumor_subclusters$subclusters$Immune
  
  # loop through each subcluster to collect cell IDs
  for (sc_name in names(immune_subclusters)) {
    cell_ids <- immune_subclusters[[sc_name]]
    true_labels <- c(true_labels, cell_ids)
  }
  
  # Create a named vector for true labels
  score_df$true_label <- ifelse(score_df$cell %in% names(true_labels), "normal", "cancer")
  
  # run silhouette cluster analysis
  pamk_result <- fpc::pamk(score_df[, c("gi_score", "correlation_score")], krange = 2:4)
  num_clusters <- pamk_result$nc
  silhouette_result <- pam(score_df, num_clusters)
  saveRDS(silhouette_result, file = paste0(base_path, "/", obj_name, "_silhouette_result.rds"))
  
  # initialize classification
  if (num_clusters > 1) {
    cluster_info <- data.frame(
      #extract cluster info
      cell = names(silhouette_result$clustering),
      cluster = silhouette_result$clustering
    )
    score_df <- merge(score_df, cluster_info, by = "cell")
    
    results <- list()
    
    # run multipliers to determine best fit
    for (x_mult in x_multipliers) {
      for (y_mult in y_multipliers) {
        result <- evaluate_thresholds_f1(score_df, x_mult, y_mult)
        results <- append(results, list(result))
      }
    }
    
    # x_int and y_int with highest f1 score should be best fit
    results_df <- do.call(rbind, lapply(results, as.data.frame))
    best_result <- results_df %>%
      arrange(desc(f1_score)) %>%
      slice_head(n=1)
    
    x_int <- best_result$x_int
    y_int <- best_result$y_int
    
    # call normal vs cancer based on best fit xy int
    score_df <- score_df %>%
      mutate(
        normal_cell_call = case_when(
          gi_score < x_int & correlation_score < y_int ~ "normal",
          TRUE ~ "cancer"
        )
      )
    
    # to count how many normal/cancer were called
    labeled_table <- table(score_df$true_label, score_df$normal_cell_call)
    labeled_df <- as.data.frame(labeled_table)
    colnames(labeled_df) <- c("TrueLabel", "PredictedLabel", "Count")
    
    # calculate total and misdiagnosed counts
    total_normals <- sum(labeled_table["normal", ])
    misdiagnosed_normals <- labeled_table["normal", "cancer"]
    total_cancers <- sum(labeled_table["cancer", ])
    misdiagnosed_cancers <- labeled_table["cancer", "normal"]
    
    # visualizations
    plot_path <- paste0(base_path, "/plots/")
    dir.create(plot_path, showWarnings = FALSE)
    
    p1 <- ggplot(score_df, aes(x = gi_score, y = correlation_score, color = true_label)) +
      geom_point() +
      scale_color_manual(values = c("normal" = "blue", "cancer" = "red")) +
      geom_vline(xintercept = x_int, linetype = "dashed") +
      geom_hline(yintercept = y_int, linetype = "dashed") +
      ggtitle(paste("True Label for", obj_name)) +
      theme_minimal()
    
    ggsave(paste0(plot_path, obj_name, "_true_label_plot.png"), p1,
           width = 10, height = 8, dpi = 150, units = "in")
    
    p2 <- ggplot(score_df, aes(x = gi_score, y = correlation_score, color = normal_cell_call)) +
      geom_point() +
      scale_color_manual(values = c("normal" = "blue", "cancer" = "red")) +
      geom_vline(xintercept = x_int, linetype = "dashed") +
      geom_hline(yintercept = y_int, linetype = "dashed") +
      ggtitle(paste("Classification for", obj_name)) +
      theme_minimal()
    
    ggsave(paste0(plot_path, obj_name, "_classification_plot.png"), p2,
           width = 10, height = 8, dpi = 150, units = "in")
    
    p3 <- ggplot(score_df, aes(x = gi_score, y = correlation_score, color = as.factor(cluster))) +
      geom_point() +
      scale_color_manual(values = c("magenta", "cyan", "green")) +
      geom_vline(xintercept = x_int, linetype = "dashed") +
      geom_hline(yintercept = y_int, linetype = "dashed") +
      ggtitle(paste("Clusters for", obj_name)) +
      theme_minimal()
    
    ggsave(paste0(plot_path, obj_name, "_clusters_plot.png"), p3,
           width = 10, height = 8, dpi = 150, units = "in")
    
    
    p4 <- ggplot(labeled_df, aes(x = TrueLabel, y = Count, fill = PredictedLabel)) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_manual(values = c("normal" = "blue", "cancer" = "red")) +
      ggtitle(paste("Comparison of Predicted and True Labels for", obj_name)) +
      theme_minimal() +
      geom_text(aes(label = Count), position = position_dodge(width = 0.9), vjust = -0.5) +
      labs(x = "True Label", y = "Count", fill = "Predicted Label") +
      annotate("text", x = 1, y = total_cancers, label = paste("Misdiagnosed Normals: ", misdiagnosed_normals, "/", total_normals, sep=""), hjust = 0, vjust = -1.5, color = "black") +
      annotate("text", x = 1, y = (total_cancers+total_normals), label = paste("Misdiagnosed Cancers: ", misdiagnosed_cancers, "/", total_cancers, sep=""), hjust = 0, vjust = -1.5, color = "black")
    
    ggsave(paste0(plot_path, obj_name, "_comparison_plot.png"), p4,
           width = 10, height = 8, dpi = 150, units = "in")
    
  } else {
    message(paste0("Sample ", obj_name, " needs to be recorded as mono-cluster and script rerun!"))
  }
  
  write.csv(score_df, file = paste0(base_path, "/", obj_name, "_cell_classification.csv"))
  
  return(list(
    classification_results = score_df,
    best_multipliers = list(
      x_multiplier = best_result$x_multiplier,
      y_multiplier = best_result$y_multiplier
    )
  ))
}
process_cnv_results_silwidth <- function(infercnv_obj, base_path, obj_name, x_multipliers = seq(0.5, 2, by = 0.25), y_multipliers = seq(0.5, 2, by = 0.25)) {
  cnv_data <- infercnv_obj@expr.data
  
  # calculate genomic instability score
  gi_scores <- calculate_gi_score(cnv_data)
  
  # determine top 5% cells for CNV profile
  top_5_percent_threshold <- quantile(gi_scores, 0.95)
  top_5_percent_cells <- names(subset(gi_scores, gi_scores >= top_5_percent_threshold))
  top_5_percent_cells_no_na <- top_5_percent_cells[!is.na(top_5_percent_cells)]
  
  # create average CNV profile across the top 5% most unstable cells
  avg_cnv_profile <- rowMeans(cnv_data[,top_5_percent_cells_no_na], na.rm = TRUE)
  
  # calculate correlation scores
  correlation_scores <- sapply(1:ncol(cnv_data), function(i) {
    cell_cnv_profile <- cnv_data[, i]
    cor(cell_cnv_profile, avg_cnv_profile, use = "complete.obs")
  })
  names(correlation_scores) <- colnames(cnv_data)
  
  # combine scores into a data frame
  score_df <- data.frame(
    cell = colnames(cnv_data),
    gi_score = gi_scores,
    correlation_score = correlation_scores
  )
  
  # run silhouette cluster analysis
  pamk_result <- fpc::pamk(score_df[, c("gi_score", "correlation_score")], krange = 2:4)
  num_clusters <- pamk_result$nc
  silhouette_result <- pam(score_df, num_clusters)
  saveRDS(silhouette_result, file = paste0(base_path, "/", obj_name, "_silhouette_result.rds"))
  
  # initialize classification
  score_df$normal_cell_call <- "cancer"
  
  # initialize classification
  if (num_clusters > 1) {
    cluster_info <- data.frame(
      #extract cluster info
      cell = names(silhouette_result$clustering),
      cluster = silhouette_result$clustering
    )
    score_df <- merge(score_df, cluster_info, by = "cell")
    
    results <- list()
    
    # run multipliers to determine best fit
    for (x_mult in x_multipliers) {
      for (y_mult in y_multipliers) {
        result <- evaluate_thresholds_silwidth(score_df, x_mult, y_mult)
        results <- append(results, list(result))
      }
    }
    
    results_df <- do.call(rbind, lapply(results, as.data.frame))
    best_result <- results_df %>%
      arrange(desc(avg_silhouette_width)) %>%
      slice_head(n=1)
    
    # define x and y intercepts for normal vs cancer classification
    x_int <- best_result$x_int
    y_int <- best_result$y_int
    
    # classify cells
    score_df <- score_df %>%
      mutate(
        normal_cell_call = case_when(
          gi_score < x_int & correlation_score < y_int ~ "normal",
          TRUE ~ "cancer"
        )
      )
    
    #visualizations
    plot_path <- paste0(base_path, "/plots/")
    dir.create(plot_path, showWarnings = FALSE)
    
    p1 <- ggplot(score_df, aes(x = gi_score, y = correlation_score, color = normal_cell_call)) +
      geom_point() +
      scale_color_manual(values = c("normal" = "blue", "cancer" = "red")) +
      geom_vline(xintercept = x_int, linetype = "dashed") +
      geom_hline(yintercept = y_int, linetype = "dashed") +
      ggtitle(paste("Classification for", obj_name)) +
      theme_minimal()
    
    ggsave(paste0(plot_path, obj_name, "_classification_plot.png"), p1,
           width = 10, height = 8, dpi = 150, units = "in")
    
    p2 <- ggplot(score_df, aes(x = gi_score, y = correlation_score, color = as.factor(cluster))) +
      geom_point() +
      scale_color_manual(values = c("magenta", "cyan", "green")) +
      geom_vline(xintercept = x_int, linetype = "dashed") +
      geom_hline(yintercept = y_int, linetype = "dashed") +
      ggtitle(paste("Clusters for", obj_name)) +
      theme_minimal()
    
    ggsave(paste0(plot_path, obj_name, "_clusters_plot.png"), p2,
           width = 10, height = 8, dpi = 150, units = "in")
  } else {
    message(paste0("Sample ", obj_name, " needs to be recorded as mono-cluster and script rerun!"))
  }
  
  write.csv(score_df, file = paste0(base_path, "/", obj_name, "_cell_classification.csv"))
  
  return(list(
    classification_results = score_df,
    best_multipliers = list(
      x_multiplier = best_result$x_multiplier,
      y_multiplier = best_result$y_multiplier
    )
  ))
}
reformatting_calls_for_merge_f1 <- function(df){
  df<-df[,2:5]
  colnames(df) = c("barcode", "gi_score", "correlation_score", "normal_cell_call_f1")
  df$barcode = sub("^X", "", df$barcode)
  df$barcode = sub("\\.(\\d+)$", "-\\1", df$barcode)
  return(df)
}
reformatting_calls_for_merge_silwidth <- function(df){
  df<-df[,2:5]
  colnames(df) = c("barcode", "gi_score", "correlation_score", "normal_cell_call_silwidth")
  df$barcode = sub("^X", "", df$barcode)
  df$barcode = sub("\\.(\\d+)$", "-\\1", df$barcode)
  return(df)
}

# analysis functions ####
combined_volcano_plot <- function(data_list, x = 'avg_log2FC', y = 'p_val_adj', pCutoff = 0.05, FCcutoff = 1, col = c('black', 'black', 'black', 'red3'), ncol = 1) {
  plots <- lapply(names(data_list), function(name) {
    data <- data_list[[name]]
    EnhancedVolcano(data, lab = data$X, x = x, y = y, pCutoff = pCutoff, FCcutoff = FCcutoff, col = col, title = name, legendPosition = 'none', 
                    labSize = 4, axisLabSize = 10, titleLabSize = 12, subtitleLabSize = 10, captionLabSize = 6)
  })
  
  combined_plot <- plot_grid(plotlist = plots, ncol = ncol)
  
  return(combined_plot)
}
create_gsea_vectors <- function(data_list) {
  gsea_data_list <- lapply(data_list, function(data) {
    gsea_data <- setNames(data$avg_log2FC, data$X)
    gsea_data <- sort(gsea_data, decreasing = TRUE)
    return(gsea_data)
  })
  
  return(gsea_data_list)
}
run_through_fgsea <- function(msigdb_list, named_degs_list, minSize = 10, maxSize = 500, nPermSimple = 50000, gseaParam = 0.5) {
  # initialize empty results list
  fgsea_results <- list()
  
  # iterate over msigdb and degs
  for (i in seq_along(msigdb_list)) {
    msigdb_name <- names(msigdb_list)[i]
    for (j in seq_along(named_degs_list)) {
      degs_name <- names(named_degs_list)[j]
      
      # print progress message since this takes a little bit of time
      message(paste("Processing:", degs_name, "through", msigdb_name))
      
      # running fgsea
      fgsea_res <- fgsea(pathways = msigdb_list[[i]], 
                         stats    = named_degs_list[[j]],
                         minSize  = minSize,
                         maxSize  = maxSize, 
                         nPermSimple = nPermSimple)
      fgsea_res$Enrichment = ifelse(fgsea_res$NES > 0, "Up-regulated", "Down-regulated")
      fgsea_results[[paste0(msigdb_name, "_", degs_name)]] <- fgsea_res
    }
  }
  return(fgsea_results)
}
generate_gsea_table_plots <- function(treatment_comparison, msigdb_list, named_degs_list, fgsea_res, search_string) {
  get_top_pathways <- function(category, comparison) {
    topPathwaysUp <- fgsea_res[[paste0(category, "_", comparison)]][ES > 0][head(order(pval), n=10), pathway]
    topPathwaysDown <- fgsea_res[[paste0(category, "_", comparison)]][ES < 0][head(order(pval), n=10), pathway]
    topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
    return(topPathways)
  }
  
  plot_gsea <- function(category, topPathways, comparison) {
    plotGseaTable(msigdb_list[[category]][topPathways], named_degs_list[[comparison]], fgsea_res[[paste0(category, "_", comparison)]], gseaParam=0.5)
  }
  
  # generate full comparisons
  suffixes <- c("allsamples", "epi", "tumor", "stro", "imm", "fullreps")
  comparisons <- paste0(treatment_comparison, ".", suffixes)
  
  # filter the comparisons to include only those that are present in the data
  available_comparisons <- comparisons[comparisons %in% names(named_degs_list)]
  
  categories <- c("all", "hallmarks", "c2", "c5", "c6", "c7")
  filtered_categories <- grep(search_string, categories, value = TRUE)
  
  plots <- list()
  
  for (comparison in available_comparisons) {
    for (cat in filtered_categories) {
      topPathways <- get_top_pathways(cat, comparison)
      plot_name <- paste0(cat, "_", comparison)
      plots[[plot_name]] <- plot_gsea(cat, topPathways, comparison)
    }
  }
  
  combined_plot <- plot_grid(plotlist = plots, labels = "AUTO", ncol = 1)
  return(combined_plot)
}
find_gene_modules_corrected <- function(cds, 
                                        reduction_method = c("UMAP"), max_components = 2, umap.metric = "cosine",
                                        umap.min_dist = 0.1, umap.n_neighbors = 15L, umap.fast_sgd = FALSE, 
                                        umap.nn_method = "annoy", k = 20, leiden_iter = 1,
                                        partition_qval = 0.05, weight = FALSE, resolution = NULL, 
                                        random_seed = 0L, cores=1, verbose = FALSE,
                                        preprocess_method = c('PCA', 'LSI'), 
                                        nn_control = list(),
                                        ...) {
  method = 'leiden'
  
  nn_control_default <- monocle3:::get_global_variable('nn_control_annoy_euclidean')
  nn_control <- monocle3:::set_nn_control(mode=3,
                                          nn_control=nn_control,
                                          nn_control_default=nn_control_default,
                                          nn_index=NULL,
                                          k=k,
                                          verbose=verbose)
  
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(preprocess_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "preprocess_method must be one of 'PCA' or 'LSI'")
  preprocess_method <- match.arg(preprocess_method)
  
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'UMAP', 'PCA' or 'tSNE'")
  reduction_method <- match.arg(reduction_method)
  
  assertthat::assert_that(methods::is(cds, "cell_data_set"))
  assertthat::assert_that(is.character(reduction_method))
  assertthat::assert_that(assertthat::is.count(k))
  assertthat::assert_that(is.logical(weight))
  assertthat::assert_that(assertthat::is.count(leiden_iter))
  ## TO DO what is resolution?
  assertthat::assert_that(is.numeric(partition_qval))
  assertthat::assert_that(is.logical(verbose))
  assertthat::assert_that(!is.null(SingleCellExperiment::reducedDims(cds)[[reduction_method]]),
                          msg = paste("No dimensionality reduction for",
                                      reduction_method, "calculated.",
                                      "Please run reduce_dimension with",
                                      "reduction_method =", reduction_method,
                                      "before running cluster_cells"))
  
  # preprocess_mat is gene_loading matrix. The gene_loadings were calculated only
  # for preprocess_method='PCA' in preprocess_cds() but I extend this to 'LSI' and
  # calculate gene_loadings here.
  preprocess_mat <- cds@reduce_dim_aux$gene_loadings # changing this so that it is not dependent on preprocessing which doesnt work with Seurat's irlba requirements
  preprocess_mat <- preprocess_mat[intersect(rownames(cds), row.names(preprocess_mat)),]
  
  # uwot::umap uses a random number generator
  if( random_seed != 0L )
    set.seed( random_seed )
  
  umap_res = uwot::umap(as.matrix(preprocess_mat),
                        n_components = max_components,
                        metric = umap.metric,
                        min_dist = umap.min_dist,
                        n_neighbors = umap.n_neighbors,
                        fast_sgd = umap.fast_sgd,
                        n_threads=cores,
                        verbose=verbose,
                        nn_method= umap.nn_method,
                        ...)
  
  row.names(umap_res) <- row.names(preprocess_mat)
  if(ncol(umap_res) < 1) warning('bad loop: ncol(umap_res) < 1')
  colnames(umap_res) <- paste0('dim_', 1:ncol(umap_res))
  reduced_dim_res <- umap_res
  
  if(verbose)
    message("Running leiden clustering algorithm ...")
  
  cluster_result <- monocle3:::leiden_clustering(data=reduced_dim_res,
                                                 pd=rowData(cds)[row.names(reduced_dim_res),,drop=FALSE],
                                                 weight=weight, nn_index=NULL, k=k, nn_control=nn_control,
                                                 num_iter=leiden_iter, resolution_parameter=resolution,
                                                 random_seed=random_seed, verbose=verbose, ...)
  
  cluster_graph_res <- monocle3:::compute_partitions(cluster_result$g, # this doesn't work!!
                                                     cluster_result$optim_res,
                                                     partition_qval, verbose)
  
  partitions <-
    igraph::components(cluster_graph_res$cluster_g)$membership[
      cluster_result$optim_res$membership]
  names(partitions) <- row.names(reduced_dim_res)
  partitions <- as.factor(partitions)
  
  gene_module_df <- tibble::tibble(id = row.names(preprocess_mat),
                                   module = factor(
                                     igraph::membership(cluster_result$optim_res)),
                                   supermodule = partitions)
  gene_module_df <- tibble::as_tibble(cbind(gene_module_df, umap_res))
  
  return(gene_module_df)
}
run_through_fgsea_modules <- function(msigdb_list, named_degs_list, minSize = 10, maxSize = 500, nPermSimple = 50000, gseaParam = 0.5) {
  # initialize empty results list
  fgsea_results <- list()
  
  # iterate over msigdb and degs
  for (i in seq_along(msigdb_list)) {
    msigdb_name <- names(msigdb_list)[i]
    for (j in seq_along(named_degs_list)) {
      degs_name <- names(named_degs_list)[j]
      
      # print progress message since this takes a little bit of time
      message(paste("Processing:", degs_name, "through", msigdb_name))
      
      # running fgsea
      fgsea_res <- fgsea(pathways = msigdb_list[[i]], 
                         stats    = named_degs_list[[j]],
                         minSize  = minSize,
                         maxSize  = maxSize, 
                         nPermSimple = nPermSimple, 
                         scoreType = "pos")
      fgsea_res$Enrichment = ifelse(fgsea_res$NES > 0, "Up-regulated", "Down-regulated")
      fgsea_results[[paste0(msigdb_name, "_", degs_name)]] <- fgsea_res
    }
  }
  return(fgsea_results)
}
calculate_lmer_modules <- function(data) {
  # initialize lists to store results
  anova_results <- list()
  summary_results <- list()
  emmeans_results <- list()
  contrasts_results <- list()
  
  # loop over the modules and fit the models
  for (i in 1:8) {
    # construct the module name
    module_name <- paste0("Module", i)
    
    # fit the linear mixed-effects model
    lm_model <- lmerTest::lmer(data = data, 
                               formula = as.formula(paste(module_name, "~ treatment + correlation_score + correlation_score*patientID + (1|patientID)")))
    
    # perform ANOVA and save the result
    anova_result <- anova(lm_model)
    anova_results[[module_name]] <- anova_result
    
    # save the summary
    summary_result <- summary(lm_model)
    summary_results[[module_name]] <- summary_result
    
    # perform emmeans analysis and save the result
    emmeans_result <- emmeans(lm_model, pairwise ~ treatment | correlation_score)
    emmeans_results[[module_name]] <- emmeans_result
    
    # format contrasts for plot
    contrasts_emm <- emmeans_result$contrasts %>%
      summary(infer = TRUE, type = 'response') %>%
      rbind() %>%
      as.data.frame()
    
    # significance stars
    contrasts_emm$p.val_corrected <- contrasts_emm$p.value
    contrasts_emm$p.val_corrected[contrasts_emm$p.value > 0.05] <- "ns"
    contrasts_emm$p.val_corrected[contrasts_emm$p.value <= 0.05] <- "*"
    contrasts_emm$p.val_corrected[contrasts_emm$p.value <= 0.01] <- "**"
    contrasts_emm$p.val_corrected[contrasts_emm$p.value < 0.001] <- "***"
    
    contrasts_results[[module_name]] <- contrasts_emm
  }
  
  # combine all results into a single list
  results <- list(anova = anova_results, summary = summary_results, emmeans = emmeans_results, contrasts = contrasts_results)
  
  return(results)
}
plotEnrichment_customcolor <- function (pathway, stats, gseaParam = 1, ticksSize = 0.2, custom_color) {
  pd <- plotEnrichmentData(pathway = pathway, stats = stats, 
                           gseaParam = gseaParam)
  with(pd, ggplot(data = curve) + geom_line(aes(x = rank, y = ES), linewidth = 2,
                                            color = custom_color) + geom_segment(data = ticks, mapping = aes(x = rank, 
                                                                                                             y = -spreadES/16, xend = rank, yend = spreadES/16), linewidth = ticksSize) + 
         geom_hline(yintercept = posES, colour = "red", linetype = "dashed") + 
         geom_hline(yintercept = negES, colour = "red", linetype = "dashed") + 
         geom_hline(yintercept = 0, colour = "black") + theme(panel.background = element_blank(), 
                                                              panel.grid.major = element_line(color = "grey92")) + 
         labs(x = "rank", y = "enrichment score"))
}