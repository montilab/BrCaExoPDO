# 3.0 load packages and data ====
source("../scripts/utils.R")
save_path <- "path/to/save/dir/" # please update file path to save in your desired location
readRDS(file=file.path(save_path, "anno_obj_list.rds")) 

# 3.1 merge into one object ====
# merge together
pdo <- merge(anno_obj_list$p318C, y = c(anno_obj_list$p318N, anno_obj_list$p318T, 
                                        anno_obj_list$p377C, anno_obj_list$p377N, anno_obj_list$p377T,
                                        anno_obj_list$p409C, anno_obj_list$p409N, anno_obj_list$p409T,
                                        anno_obj_list$p409C2, anno_obj_list$p409N2, anno_obj_list$p409T2),
             add.cell.ids = c("318C", "318N", "318T", 
                              "377C", "377N", "377T",
                              "409C", "409N", "409T",
                              "409C2", "409N2", "409T2"),
             project = "pdo")
pdo #24,599 cells (already undergone filtering and doublet removal)

# 3.2 add group labels ====
pdo@meta.data$barcode = rownames(pdo@meta.data)

pdo@meta.data$group = rownames(pdo@meta.data)
pdo@meta.data$group = gsub("_.*", "", pdo@meta.data$group)

pdo@meta.data$treatment[pdo@meta.data$group %in% c("318C", "377C", "409C", "409C2")] <- "UT"
pdo@meta.data$treatment[pdo@meta.data$group %in% c("318N", "377N", "409N", "409N2")] <- "NDexo"
pdo@meta.data$treatment[pdo@meta.data$group %in% c("318T", "377T", "409T", "409T2")] <- "T2Dexo"

pdo@meta.data$patientID[pdo@meta.data$group %in% c("318C", "318N", "318T")] <- "318"
pdo@meta.data$patientID[pdo@meta.data$group %in% c("377C", "377N", "377T")] <- "377"
pdo@meta.data$patientID[pdo@meta.data$group %in% c("409C", "409N", "409T", "409C2", "409N2", "409T2")] <- "409"

pdo@meta.data$batch[pdo@meta.data$group %in% c("318C", "318N", "318T")] <- "rep1_318"
pdo@meta.data$batch[pdo@meta.data$group %in% c("377C", "377N", "377T")] <- "rep1_377"
pdo@meta.data$batch[pdo@meta.data$group %in% c("409C", "409N", "409T")] <- "rep1_409"
pdo@meta.data$batch[pdo@meta.data$group %in% c("409C2", "409N2", "409T2")] <- "rep2_409"

# adding in singleR major labels
major_labels = dplyr::select(swarbrick_metadata, celltype_minor, celltype_major)
major_labels = unique(major_labels)

for (category in names(table(major_labels$celltype_major))) {
  pdo@meta.data$singleR_majorlabels[pdo@meta.data$singleR_labels %in% subset(major_labels, celltype_major == category)$celltype_minor] <- category
}

# 3.3 check QC now that everything is together ====
p1 <- scCustomize::QC_Plots_Genes(seurat_object = pdo, low_cutoff = 100, high_cutoff = 8500)
p2 <- scCustomize::QC_Plots_UMIs(seurat_object = pdo, low_cutoff = 500, high_cutoff = 50000)
p3 <- scCustomize::QC_Plots_Mito(seurat_object = pdo, high_cutoff = 20)
p4 <- scCustomize::QC_Plots_Complexity(seurat_object = pdo, high_cutoff = 0.75)
wrap_plots(p1, p2, p3, p4, ncol = 4) #looks great

pdo@meta.data %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(
    min_value = min(log10GenesPerUMI),
    max_value = max(log10GenesPerUMI),
    median_value = median(log10GenesPerUMI)
  )

# 3.4 isolating compartments #### 
pdo@meta.data$coarse_anno[pdo@meta.data$singleR_majorlabels %in% 
                            c("Cancer Epithelial", "Normal Epithelial")] <- "Epithelial"
pdo@meta.data$coarse_anno[pdo@meta.data$singleR_majorlabels %in% 
                            c("CAFs", "Endothelial", "PVL")] <- "Stromal"
pdo@meta.data$coarse_anno[pdo@meta.data$singleR_majorlabels %in% 
                            c("B-cells", "Myeloid", "Plasmablasts", "T-cells")] <- "Immune"
pdo_epi = subset(x = pdo, subset = coarse_anno == "Epithelial")
pdo_stro = subset(x = pdo, subset = coarse_anno == "Stromal")
pdo_imm = subset(x = pdo, subset = coarse_anno == "Immune")
# 3.5 save objects ====
saveRDS(pdo, file = save_path, "pdo_unintegrated.rds")
saveRDS(pdo_epi, file = save_path, "pdo_epi_unintegrated.rds")
saveRDS(pdo_stro, file = save_path, "pdo_stro_unintegrated.rds")
saveRDS(pdo_imm, file = save_path, "pdo_imm_unintegrated.rds")
