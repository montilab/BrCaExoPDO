# 6.0 load packages and data ====
source("../scripts/utils.R")
save_path <- "path/to/save/dir/" # please update file path to save in your desired location
pdo_unintegrated <- readRDS(file.path(save_path, "pdo_unintegrated.rds"))
pdo <- readRDS(file.path(save_path, "pdo_integrated.rds"))
pdo_epi_postintegration <- readRDS(file.path(save_path, "pdo_epi_split_post_integrated.rds"))
pdo_epi <- readRDS(file.path(save_path, "pdo_epi_integrated.rds"))
pdo_stro_postintegration <- readRDS(file.path(save_path, "pdo_stro_split_post_integrated.rds"))
pdo_stro <- readRDS(file.path(save_path, "pdo_stro_integrated.rds"))
pdo_imm_postintegration <- readRDS(file.path(save_path, "pdo_imm_split_post_integrated.rds"))
pdo_imm <- readRDS(file.path(save_path, "pdo_imm_integrated.rds"))

# 6.1 cell cycle scoring ====
pdo <- CellCycleScoring(pdo,
                        s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, 
                        set.ident = TRUE, assay = 'SCT')

scCustomize::DimPlot_scCustom(pdo, reduction = "umap.scvi") 
FeaturePlot(pdo, features = "S.Score", reduction = "umap.scvi")
FeaturePlot(pdo, features = "G2M.Score", reduction = "umap.scvi")

FeaturePlot(pdo, features = "percent_mito", reduction = "umap.scvi")
FeaturePlot(pdo, features = "percent_ribo", reduction = "umap.scvi")
# 6.2 stromal fine tune annotations ####
plot1 <- scCustomize::DimPlot_scCustom(pdo_stro_postintegration, reduction = "umap", figure_plot = T, group.by = "patientID", colors_use = patient_colors)
plot2 <- scCustomize::DimPlot_scCustom(pdo_stro, reduction = "umap.scvi", figure_plot = T, group.by = "patientID", colors_use = patient_colors)
plot_grid(plot1, plot2) # separately integrated looks better

pdo_stro <- CellCycleScoring(pdo_stro,
                             s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, 
                             set.ident = TRUE, assay = 'SCT')

plot1 <- scCustomize::DimPlot_scCustom(pdo_stro, reduction = "umap.scvi", figure_plot = T, group.by = "scvi_clusters")
plot2 <- scCustomize::DimPlot_scCustom(pdo_stro, reduction = "umap.scvi", figure_plot = T, group.by = "singleR_labels")
plot3 <- scCustomize::DimPlot_scCustom(pdo_stro, reduction = "umap.scvi") 
plot4 <- scCustomize::DimPlot_scCustom(pdo_stro, reduction = "umap.scvi", group.by = "patientID", colors_use = patient_colors) 
plot5 <- scCustomize::DimPlot_scCustom(pdo_stro, reduction = "umap.scvi", group.by = "treatment", colors_use = treatment_colors) 

plot_grid(plot1, plot2, plot3, plot4, plot5, ncol = 3)

# finding marker genes for each cluster
Idents(object = pdo_stro) <- "scvi_clusters"
pdo_stro <- PrepSCTFindMarkers(pdo_stro, assay = "SCT", verbose = TRUE) # needs to be done on SCT merged data
# stro.markers <- FindAllMarkers(pdo_stro, only.pos = TRUE,
#                                min.pct = 0.25, logfc.threshold = 0.25,
#                                test.use = "MAST", latent.vars = "batch")
stro.markers = read.csv("../results/stro_FindAllMarkers.csv")

stro.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.25) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5
DoHeatmap(pdo_stro, features = top5$gene) + NoLegend()

#to explore marker genes
stro.markers %>% 
  dplyr::filter(cluster == 0) %>% # Change this to explore markers of other clusters
  dplyr::arrange(p_val_adj, desc(avg_log2FC)) %>%
  dplyr::slice_head(n=20) %>%
  dplyr::pull(gene)

# STROMAL ANNOTATIONS
# basics:
# CAF - PDGFRA and COL1A1
# iCAF - ALDH1A1, KLF4 and LEPR, CXCL12 and C3
# myCAF - ACTA2, TAGLN, FAP and COL1A1
# PVL - MCAM, ACTA2 and PDGFRB
# imPVL - PDGFRB, ALDH1A1, CD44, CSPG4, RGS5 and CD36, ICAM1, VCAM1 and ITGB1
# dPVL - MYH11 and ACTA2
# endothelial - PECAM1 and CD34

# 0 - imPVL (RGS5, CD36)
# 1 - myCAF (FN1, COL1A2)
# 2 - dPVL
# 3 - dPVL
# 4 - myCAF (FN1, COL6A1)
# 5 - iCAF (BASP1, ANXA1, MMP2)
# 6 - iCAF
# 7 - imPVL
# 8 - endothelial (PECAM1, CD34)

pdo_stro[[]]$fine_anno[pdo_stro[[]]$scvi_clusters %in% c("0", "7")] <- "imPVL"
pdo_stro[[]]$fine_anno[pdo_stro[[]]$scvi_clusters %in% c("1", "4")] <- "myCAF"
pdo_stro[[]]$fine_anno[pdo_stro[[]]$scvi_clusters %in% c("2", "3")] <- "dPVL"
pdo_stro[[]]$fine_anno[pdo_stro[[]]$scvi_clusters %in% c("5", "6")] <- "iCAF"
pdo_stro[[]]$fine_anno[pdo_stro[[]]$scvi_clusters %in% c("8")] <- "Endothelial"

# adding labels to main object
for (category in names(table(pdo_stro@meta.data$fine_anno))) {
  pdo@meta.data$fine_anno[pdo@meta.data$barcode %in% subset(pdo_stro@meta.data, fine_anno == category)$barcode] <- category
}

# 6.3 immune fine tune annotations ####
plot1 <- scCustomize::DimPlot_scCustom(pdo_imm_postintegration, reduction = "umap", figure_plot = T, group.by = "patientID", colors_use = patient_colors)
plot2 <- scCustomize::DimPlot_scCustom(pdo_imm, reduction = "umap.scvi", figure_plot = T, group.by = "patientID", colors_use = patient_colors)
plot_grid(plot1, plot2) # look similar

plot1 <- scCustomize::DimPlot_scCustom(pdo_imm_postintegration, reduction = "umap", figure_plot = T, group.by = "SCT_snn_res.1")
plot2 <- scCustomize::DimPlot_scCustom(pdo_imm_postintegration, reduction = "umap", figure_plot = T, group.by = "singleR_labels")
plot3 <- scCustomize::DimPlot_scCustom(pdo_imm_postintegration, reduction = "umap", group.by = "patientID", colors_use = patient_colors) 
plot4 <- scCustomize::DimPlot_scCustom(pdo_imm_postintegration, reduction = "umap", group.by = "treatment", colors_use = treatment_colors) 
plot_grid(plot1, plot2, plot3, plot4, ncol = 2)

# finding marker genes for each cluster
Idents(object = pdo_imm_postintegration) <- "SCT_snn_res.1"
# imm.markers <- FindAllMarkers(pdo_imm_postintegration, only.pos = TRUE,
#                               min.pct = 0.25, logfc.threshold = 0.25,
#                               test.use = "MAST", latent.vars = "batch")
imm.markers = read.csv("../results/imm_FindAllMarkers.csv")

imm.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.25) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5
DoHeatmap(pdo_imm_postintegration, features = top5$gene) + NoLegend()

# to explore marker genes
imm.markers %>% 
  dplyr::filter(cluster == 0) %>% # Change this to explore markers of other clusters
  dplyr::arrange(p_val_adj, desc(avg_log2FC)) %>%
  dplyr::slice_head(n=20) %>%
  dplyr::pull(gene)

# IMMUNE ANNOTATIONS
# 0 - CTL (CCL5, GZMK, CD8A, KLRK1, CD8B, PIK3R1)
# 1 - Tcm (IL7R, SELL+, CCR7+, CXCR3+)
# 2 - MAIT (KLRB1 aka CD161, RORA) 
# 3 - Tem (IL7R, SELL low, CCR7 low)
# 4 - B cells (CD79A, CD74, MS4A1, BANK1, CD37)
# 5 - stressed T (DDIT3, XBP1, lncRNAs)
# 6 - Th1 (IL32, GZMB, GZMA, IFNG, STAT1)
# 7 - plasma (IGLC2, IGLC3, IGHA1, IGLC1, MS4A1, BANK1)
# 8 - Macrophage (CD68, LYZ)
# 9 - plasma (JCHAIN, MZB1)

pdo_imm_postintegration[[]]$fine_anno[pdo_imm_postintegration[[]]$SCT_snn_res.1 %in% c("0")] <- "Teff"
pdo_imm_postintegration[[]]$fine_anno[pdo_imm_postintegration[[]]$SCT_snn_res.1 %in% c("1")] <- "Tcm"
pdo_imm_postintegration[[]]$fine_anno[pdo_imm_postintegration[[]]$SCT_snn_res.1 %in% c("2")] <- "MAIT"
pdo_imm_postintegration[[]]$fine_anno[pdo_imm_postintegration[[]]$SCT_snn_res.1 %in% c("3")] <- "Tem"
pdo_imm_postintegration[[]]$fine_anno[pdo_imm_postintegration[[]]$SCT_snn_res.1 %in% c("4")] <- "B"
pdo_imm_postintegration[[]]$fine_anno[pdo_imm_postintegration[[]]$SCT_snn_res.1 %in% c("5")] <- "ChopT"
pdo_imm_postintegration[[]]$fine_anno[pdo_imm_postintegration[[]]$SCT_snn_res.1 %in% c("6")] <- "Th1"
pdo_imm_postintegration[[]]$fine_anno[pdo_imm_postintegration[[]]$SCT_snn_res.1 %in% c("7", "9")] <- "Plasma"
pdo_imm_postintegration[[]]$fine_anno[pdo_imm_postintegration[[]]$SCT_snn_res.1 %in% c("8")] <- "Macrophage"

# adding labels to main object
for (category in names(table(pdo_imm_postintegration@meta.data$fine_anno))) {
  pdo@meta.data$fine_anno[pdo@meta.data$barcode %in% subset(pdo_imm_postintegration@meta.data, fine_anno == category)$barcode] <- category
}

DimPlot_scCustom(pdo_imm_postintegration, group.by = "fine_anno", reduction = "umap", colors_use = imm_fine_anno_colors_alpha, figure_plot = T)

FeaturePlot_scCustom(pdo_imm_postintegration, reduction = "umap",
                     features = c("CD3D", "CD8A", "CD4", 
                                  "IL7R", "SELL", "CCR7",
                                  "IFNG", "LTA", "KLRB1",
                                  "DDIT3", "PPP1R15A"), figure_plot = T, num_columns = 3)

# 6.4 epithelial fine tune annotations ####
plot1 <- scCustomize::DimPlot_scCustom(pdo_epi_postintegration, reduction = "umap", figure_plot = T, group.by = "patientID", colors_use = patient_colors)
plot2 <- scCustomize::DimPlot_scCustom(pdo_epi, reduction = "umap.scvi", figure_plot = T, group.by = "patientID", colors_use = patient_colors)
plot_grid(plot1, plot2) # separately integrated looks better

# updating obj with normal/tumor calls
for (barcode in unique(pdo@meta.data$barcode)) {
  # find gi_score, correlation_score, f1 call, silwidth call
  gi_score <- pdo@meta.data[pdo@meta.data$barcode == barcode, 'gi_score']
  correlation_score <- pdo@meta.data[pdo@meta.data$barcode == barcode, 'correlation_score']
  f1_call <- pdo@meta.data[pdo@meta.data$barcode == barcode, 'f1_call']
  silwidth_call <- pdo@meta.data[pdo@meta.data$barcode == barcode, 'silwidth_call']
  manual_call <- pdo@meta.data[pdo@meta.data$barcode == barcode, 'manual_call']
  
  # assign to the corresponding rows in destination 
  pdo_epi@meta.data$gi_score[pdo_epi@meta.data$barcode == barcode] <- gi_score
  pdo_epi@meta.data$correlation_score[pdo_epi@meta.data$barcode == barcode] <- correlation_score
  pdo_epi@meta.data$f1_call[pdo_epi@meta.data$barcode == barcode] <- f1_call
  pdo_epi@meta.data$silwidth_call[pdo_epi@meta.data$barcode == barcode] <- silwidth_call
  pdo_epi@meta.data$manual_call[pdo_epi@meta.data$barcode == barcode] <- manual_call
}

# adding in checks for normal/tumor
# compare to cancer gene census breast cancer genes (in gene module score) and msigDB to see how well it did calling tumor
# please export CSV from COSMIC cancer gene census: https://cancer.sanger.ac.uk/census
cgc <- read.csv("path_to/cosmic_cgc.csv") # please update file path
cgc = dplyr::select(cgc, Gene.Symbol, Genome.Location, Tier, Hallmark, Chr.Band, Tumour.Types.Somatic.)
cgc = subset(cgc, Tumour.Types.Somatic. %like% "breast")
cgc_genes = c(cgc$Gene.Symbol)

# check against msigDB signatures
# please export gmt file of GOBP_POSITIVE_REGULATION_OF_MAMMARY_GLAND_EPITHELIAL_CELL_PROLIFERATION: https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/GOBP_POSITIVE_REGULATION_OF_MAMMARY_GLAND_EPITHELIAL_CELL_PROLIFERATION.html
gobp_mamm_epi_prolif = read.gmt("path_to/gobp_gmt.gmt") # please update file path

bc_tumor_genesets = list(cgc_genes, 
                         gobp_mamm_epi_prolif$GOBP_POSITIVE_REGULATION_OF_MAMMARY_GLAND_EPITHELIAL_CELL_PROLIFERATIONP)

pdo_epi <- Seurat::AddModuleScore(pdo_epi, features = bc_tumor_genesets)

pdo_epi@meta.data <- pdo_epi@meta.data %>% dplyr::rename("cgc_exp" = `Cluster1`,
                                                         "gobp_mamm_epi_prolif" = `Cluster2`) 

scCustomize::FeaturePlot_scCustom(pdo_epi, features = c("correlation_score", "cgc_exp", 
                                                        "gobp_mamm_epi_prolif"), 
                                  figure_plot = T, order = T, num_columns = 2, reduction = "umap.scvi")

pdo_epi <- CellCycleScoring(pdo_epi,
                            s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, 
                            set.ident = TRUE, assay = 'SCT')

plot1 <- scCustomize::DimPlot_scCustom(pdo_epi, reduction = "umap.scvi", figure_plot = T, group.by = "scvi_clusters")
plot2 <- scCustomize::DimPlot_scCustom(pdo_epi, reduction = "umap.scvi", figure_plot = T, group.by = "singleR_labels")
plot3 <- scCustomize::DimPlot_scCustom(pdo_epi, reduction = "umap.scvi") 
plot4 <- scCustomize::DimPlot_scCustom(pdo_epi, reduction = "umap.scvi", group.by = "patientID", colors_use = patient_colors) 
plot5 <- scCustomize::DimPlot_scCustom(pdo_epi, reduction = "umap.scvi", group.by = "treatment", colors_use = treatment_colors) 

plot_grid(plot1, plot2, plot3, plot4, plot5, ncol = 3)

plot1 <- scCustomize::DimPlot_scCustom(pdo_epi, reduction = "umap.scvi", figure_plot = T, group.by = "f1_call", colors_use = c("red", "blue"))
plot2 <- scCustomize::DimPlot_scCustom(pdo_epi, reduction = "umap.scvi", figure_plot = T, group.by = "silwidth_call", colors_use = c("red", "blue"))
plot3 <- scCustomize::DimPlot_scCustom(pdo_epi, reduction = "umap.scvi", figure_plot = T, group.by = "manual_call", colors_use = c("red", "blue"))

plot_grid(plot1, plot2, plot3, ncol = 3)

# finding marker genes for each cluster
Idents(object = pdo_epi) <- "scvi_clusters"
pdo_epi <- PrepSCTFindMarkers(pdo_epi, assay = "SCT", verbose = TRUE) # needs to be done on SCT merged data
# epi.markers <- FindAllMarkers(pdo_epi, only.pos = TRUE,
#                                  min.pct = 0.25, logfc.threshold = 0.25,
#                                  test.use = "MAST", latent.vars = "batch", recorrect_umi = F)
epi.markers = read.csv("../results/epi_FindAllMarkers.csv")

epi.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.25) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5
DoHeatmap(pdo_epi, features = top5$gene) + NoLegend()

#to explore marker genes
epi.markers %>% 
  dplyr::filter(cluster == 0) %>% # Change this to explore markers of other clusters
  dplyr::arrange(p_val_adj, desc(avg_log2FC)) %>%
  dplyr::slice_head(n=20) %>%
  dplyr::pull(gene)

# EPITHELIAL ANNOTATIONS

# 0 - myoepithelial or basal (KRT14, ITGA6 (basal epi marker), )
# 1 - luminal progenitor (ELF5, KIT, LTF, SCGB1D2, ALDH1A3)
# 2 - luminal tumor (TRPS1, MALAT1, NEAT1, metabolic activity)
# 3 - lumtum (based on ESR PGR ERBB2 in main umap)
# 4 - myoepithelial (ACTA2, MYLK, TAGLN, SYNM, JUNB, FOS, SPARCL1)
# 5 - cycling epithelial (MKI67, TOP2A, cell cycle score)
# 6 - lumtum (cgc_exp, TFF3, AGR2, GALNT6, XBP1, DHCR24, ASPH, TFAP2B and TBX3)
# 7 - LEC (silwidth and manual call, cgc_exp and gobp_mamm_prolif)
# 8 - LEC (silwidth and manual call, cgc_exp and gobp_mamm_prolif)
# 9 - immune evasive tumor? (HLA-A, HLA-B, HLA-C and B2M)
# 10 - normal... LEC or MEC? (silwidth call, cgc_exp and gobp_mamm_prolif - MEC based off clustree)
# 11 - MEC or basal tumor
# 12 - LEC (silwidth and manual call, cgc_exp and gobp_mamm_prolif)
# 13 - immune evasive tumor? (MUC1, AGR2, gobp_mamm_epi_prolif, cgc_exp, SERPINA3, CFB, and CD46)

pdo_epi[[]]$fine_anno[pdo_epi[[]]$scvi_clusters %in% c("0")] <- "Basal_Tumor1"
pdo_epi[[]]$fine_anno[pdo_epi[[]]$scvi_clusters %in% c("11")] <- "Basal_Tumor2"
pdo_epi[[]]$fine_anno[pdo_epi[[]]$scvi_clusters %in% c("1")] <- "Luminal_Progenitor"
pdo_epi[[]]$fine_anno[pdo_epi[[]]$scvi_clusters %in% c("2")] <- "Luminal_Tumor1"
pdo_epi[[]]$fine_anno[pdo_epi[[]]$scvi_clusters %in% c("3")] <- "Luminal_Tumor2"
pdo_epi[[]]$fine_anno[pdo_epi[[]]$scvi_clusters %in% c("6")] <- "Luminal_Tumor3"
pdo_epi[[]]$fine_anno[pdo_epi[[]]$scvi_clusters %in% c("9")] <- "Luminal_Tumor4"
pdo_epi[[]]$fine_anno[pdo_epi[[]]$scvi_clusters %in% c("13")] <- "Luminal_Tumor5"
pdo_epi[[]]$fine_anno[pdo_epi[[]]$scvi_clusters %in% c("4", "10")] <- "MEC"
pdo_epi[[]]$fine_anno[pdo_epi[[]]$scvi_clusters %in% c("5")] <- "Cycling_Epithelial"
pdo_epi[[]]$fine_anno[pdo_epi[[]]$scvi_clusters %in% c("7", "8", "12")] <- "LEC"

# adding labels to main object
for (category in names(table(pdo_epi@meta.data$fine_anno))) {
  pdo@meta.data$fine_anno[pdo@meta.data$barcode %in% subset(pdo_epi@meta.data, fine_anno == category)$barcode] <- category
}

DimPlot_scCustom(pdo, group.by = "fine_anno", reduction = "umap.scvi")

DimPlot_scCustom(pdo_epi, group.by = "fine_anno", reduction = "umap.scvi",
                 figure_plot = T, colors_use = epi_fine_anno_colors_alpha)
DimPlot_scCustom(pdo_epi, group.by = "treatment", reduction = "umap.scvi",
                 figure_plot = T, colors_use = treatment_colors)

# checking csc
# https://geneglobe.qiagen.com/us/product-groups/rt2-profiler-pcr-arrays/PAHS-176Z
csc_genes <- c("ABCB5","ALCAM", "ATXN1", "BMI1", "CD24", "CD38", "CD44", "ENG", "ETFA", 
               "FLOT2", "GATA3", "ITGA2", "ITGA4", "ITGA6", "ITGB1", "KIT", "MS4A1", "MUC1", "PROM1", "PTPRC", "THY1", # general csc markers
               "EGF",  "ERBB2", "KITLG", "LIN28B", "NOS2", # cell proliferation
               "BMP7", "DNMT1", "FGFR2", # self renewal
               "KLF4", "LIN28A", "MYC", "NANOG", "POU5F1", "SOX2", # pluripotency 
               "FOXP1", "HDAC1", "MYCN", "SIRT1", "WNT1", # asymmetric division
               "AXL", "ID1", "CXCL8", "KLF17", "PLAT", "PLAUR", "SNAI1", "TWIST1", "TWIST2", "ZEB1", "ZEB2", # migration/metastasis
               "DLL1", "DLL4",  "GATA3", "JAG1", "JAK2",  "LATS1", "MAML1", "MERTK", 
               "NOTCH1", "NOTCH2", "SAV1",  "TAZ", "WWC1", "YAP1", # signal transduction
               "ABCG2","ATM", "CHEK1", "DDR1", "DKK1", "EPCAM", "FZD7",  "GSK3B", "IKBKB", 
               "NFKB1",  "SMO", "STAT3","TGFBR1","WEE1"  # cancer therapeutic targets
)

pdo <- Seurat::AddModuleScore(pdo, features = list(csc_genes))

pdo@meta.data <- pdo@meta.data %>% dplyr::rename("csc_exp" = `Cluster1`) 

scCustomize::FeaturePlot_scCustom(pdo, features = c("csc_exp"), 
                                  figure_plot = T, reduction = "umap.scvi")

# 6.5 cleaning up metadata and simple overview figures ####
# adding broad anno
pdo@meta.data$finalized_fine_anno = pdo@meta.data$fine_anno
pdo@meta.data$finalized_broad_anno[pdo@meta.data$finalized_fine_anno %in% c("B", "Plasma")] <- "B"
pdo@meta.data$finalized_broad_anno[pdo@meta.data$finalized_fine_anno %in% c("Basal_Tumor1", "Basal_Tumor2")] <- "Basal_Tumor"
pdo@meta.data$finalized_broad_anno[pdo@meta.data$finalized_fine_anno %in% c("iCAF", "myCAF")] <- "CAF"
pdo@meta.data$finalized_broad_anno[pdo@meta.data$finalized_fine_anno %in% c("Cycling_Epithelial")] <- "Cycling_Epithelial"
pdo@meta.data$finalized_broad_anno[pdo@meta.data$finalized_fine_anno %in% c("Endothelial")] <- "Endothelial"
pdo@meta.data$finalized_broad_anno[pdo@meta.data$finalized_fine_anno %in% c("LEC")] <- "LEC"
pdo@meta.data$finalized_broad_anno[pdo@meta.data$finalized_fine_anno %in% c("Luminal_Progenitor")] <- "Luminal_Progenitor"
pdo@meta.data$finalized_broad_anno[pdo@meta.data$finalized_fine_anno %in% c("Luminal_Tumor1", "Luminal_Tumor2", "Luminal_Tumor3", "Luminal_Tumor4", "Luminal_Tumor5")] <- "Luminal_Tumor"
pdo@meta.data$finalized_broad_anno[pdo@meta.data$finalized_fine_anno %in% c("Macrophage")] <- "Macrophage"
pdo@meta.data$finalized_broad_anno[pdo@meta.data$finalized_fine_anno %in% c("MEC")] <- "MEC"
pdo@meta.data$finalized_broad_anno[pdo@meta.data$finalized_fine_anno %in% c("dPVL", "imPVL")] <- "PVL"
pdo@meta.data$finalized_broad_anno[pdo@meta.data$finalized_fine_anno %in% c("ChopT", "MAIT", "Tcm", "Teff", "Tem", "Th1")] <- "T"

# removing a few unnecessary columns
meta_data <- dplyr::select(pdo@meta.data, -orig.ident, -singleR_majorlabels, -seurat_clusters,
                           -old.ident, -fine_anno)

rownames(meta_data) <- meta_data$barcode
pdo@meta.data <- meta_data

# nice figures
DimPlot_scCustom(pdo, reduction = "umap.scvi", group.by = "finalized_broad_anno", 
                 colors_use = broadcelltype_colors_alpha, figure_plot = T)

DimPlot_scCustom(pdo, reduction = "umap.scvi", group.by = "treatment", 
                 colors_use = treatment_colors, figure_plot = T)

# feature plot simple figures
extra_yelred <- colorRampPalette(brewer.pal(9,"YlOrRd"))
FeaturePlot_scCustom(pdo, reduction = "umap.scvi", features = "correlation_score", figure_plot = T,
                     colors_use = extra_yelred(100), na_color = "#FFEDA0")

oncotype = c("MKI67", "AURKA", "BIRC5", "CCND1", "MYBL2",
             "STMY3", "CTSV", "GRB7", "ERBB2", "ESR1", "PGR",
             "BCL1", "SCUBE2", "GSTM1", "BAG1", "CD68")
pdo <- Seurat::AddModuleScore(pdo, features = list(oncotype))
pdo@meta.data <- pdo@meta.data %>% dplyr::rename("oncotype_exp" = `Cluster1`)

plot1 <- ggplot(subset(pdo@meta.data, coarse_anno =="Epithelial"), aes(y=oncotype_exp, x = treatment))+
  geom_violin(aes(fill = treatment))+
  geom_boxplot(fill = "white", outlier.shape = NA, width = 0.2)+
  theme_minimal()+
  xlab("") + ylab("Oncotype Expression")+
  theme(legend.position = "none")+
  scale_fill_manual(values = treatment_colors)+
  stat_compare_means(comparisons = list(c("UT", "NET"), c("UT", "DET"), c("NET", "DET")),
                     method="wilcox.test",label = "p.signif", 
                     symnum.args = list(cutpoints = c(0, 0.01, 0.05)))
plot2 <- ggplot(subset(pdo@meta.data, coarse_anno =="Epithelial"), aes(y=csc_exp, x = treatment))+
  geom_violin(aes(fill = treatment))+
  geom_boxplot(fill = "white", outlier.shape = NA, width = 0.2)+
  theme_minimal()+
  xlab("") + ylab("CSC Marker Expression")+
  theme(legend.position = "none")+
  scale_fill_manual(values = treatment_colors)+
  stat_compare_means(comparisons = list(c("UT", "NET"), c("UT", "DET"), c("NET", "DET")),
                     method="wilcox.test",label = "p.signif", 
                     symnum.args = list(cutpoints = c(0, 0.01, 0.05)))
plot3 <- ggplot(subset(pdo@meta.data, coarse_anno =="Epithelial"), aes(y=cgc_exp, x = treatment))+
  geom_violin(aes(fill = treatment))+
  geom_boxplot(fill = "white", outlier.shape = NA, width = 0.2)+
  theme_minimal()+
  xlab("") + ylab("COSMIC CGC Expression")+
  theme(legend.position = "none")+
  scale_fill_manual(values = treatment_colors)+
  stat_compare_means(comparisons = list(c("UT", "NET"), c("UT", "DET"), c("NET", "DET")),
                     method="wilcox.test",label = "p.signif", 
                     symnum.args = list(cutpoints = c(0, 0.01, 0.05)))

plot_grid(plot1, plot2, plot3, ncol = 3)

# 6.6 save progress ====
saveRDS(pdo, file = file.path(save_path, "pdo_integrated.rds"))
saveRDS(pdo_epi, file = file.path(save_path, "pdo_epi_integrated.rds"))
saveRDS(pdo_stro, file = file.path(save_path, "pdo_stro_integrated.rds"))
saveRDS(pdo_imm, file = file.path(save_path, "pdo_imm_integrated.rds"))
saveRDS(pdo_epi_postintegration, file = file.path(save_path, "pdo_epi_split_post_integrated.rds"))
saveRDS(pdo_stro_postintegration, file = file.path(save_path, "pdo_stro_split_post_integrated.rds"))
saveRDS(pdo_imm_postintegration, file = file.path(save_path, "pdo_imm_split_post_integrated.rds"))