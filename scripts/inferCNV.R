# 5.0 load packages and data ====
source("../scripts/utils.R")
save_path <- "path/to/save/dir/" # please update file path to save in your desired location
pdo_unintegrated <- readRDS(file.path(save_path, "pdo_unintegrated.rds"))
pdo <- readRDS(file.path(save_path, "pdo_integrated.rds"))

# 5.1 normal/tumor inferCNV calling ####
# creating patient subsets
cnv318 = subset(x = pdo_unintegrated, subset = patientID == "318")
cnv377 = subset(x = pdo_unintegrated, subset = patientID == "377")
cnv409 = subset(x = pdo_unintegrated, subset = patientID == "409")

cnv318 = subset(x = cnv318, subset = coarse_anno %in% c("Immune", "Epithelial"))
cnv377 = subset(x = cnv377, subset = coarse_anno %in% c("Immune", "Epithelial"))
cnv409 = subset(x = cnv409, subset = coarse_anno %in% c("Immune", "Epithelial"))

# process each object
results_409 <- process_inferCNV(cnv409, "p409", file.path(save_path, "inferCNV/"))
results_377 <- process_inferCNV(cnv377, "p377", file.path(save_path, "inferCNV/"))
results_318 <- process_inferCNV(cnv318, "p318", file.path(save_path, "inferCNV/"))

# 5.2 benchmarking normal v tumor ====
# via: https://www.nature.com/articles/s41588-021-00911-1
# developed by: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6703186/
# calculate genomic instability score: use the inferred CNV changes to compute a score for each cell
# determine top Cells for CNV Profile: identify the top 5% of cells with the highest genomic instability scores to create an average CNV profile
# calculate correlation scores: correlate each cell’s CNV profile with the average profile to obtain correlation scores
# plot and cluster cells: use these scores to plot cells and apply clustering for classification
# set thresholds for classification: define thresholds for identifying normal and tumor cells based on clustering results

results_409 <- readRDS(file.path(save_dir, "inferCNV/run_p409_hmm/run.final.infercnv_obj"))
results_377 <- readRDS(file.path(save_dir, "inferCNV/run_p377_hmm/run.final.infercnv_obj"))
results_318 <- readRDS(file.path(save_dir, "inferCNV/run_p318_hmm/run.final.infercnv_obj"))

# process each inferCNV result with each method
results_409_called <- process_cnv_results_f1(results_409, file.path(save_path, "tumor_calling_f1/"), "p409")
results_318_called <- process_cnv_results_f1(results_318, file.path(save_path, "tumor_calling_f1/"), "p418")
results_377_called <- process_cnv_results_f1(results_377, file.path(save_path, "tumor_calling_f1/"), "p377")

results_409_called <- process_cnv_results_silwidth(results_409, file.path(save_path, "tumor_calling_silwidth/"), "p409")
results_318_called <- process_cnv_results_silwidth(results_318, file.path(save_path, "tumor_calling_silwidth/"), "p318")
results_377_called <- process_cnv_results_silwidth(results_377, file.path(save_path, "tumor_calling_silwidth/"), "p377")

# 5.3 adding calls to main object ====
# load in called normal/tumor 
p318_cell_classification_f1 <- read.csv(file.path(save_path, "inferCNV/tumor_calling_f1/p318_cell_classification.csv"))
p377_cell_classification_f1 <- read.csv(file.path(save_path, "inferCNV/tumor_calling_f1/p377_cell_classification.csv"))
p409_cell_classification_f1 <- read.csv(file.path(save_path, "inferCNV/tumor_calling_f1/p409_cell_classification.csv"))
p318_cell_classification_silwidth <- read.csv(file.path(save_path, "inferCNV/tumor_calling_silwidth/p318_cell_classification.csv"))
p377_cell_classification_silwidth <- read.csv(file.path(save_path, "inferCNV/tumor_calling_silwidth/p377_cell_classification.csv"))
p409_cell_classification_silwidth <- read.csv(file.path(save_path, "inferCNV/tumor_calling_silwidth/p409_cell_classification.csv"))

p318_cell_classification_f1 = reformatting_calls_for_merge_f1(p318_cell_classification_f1)
p377_cell_classification_f1 = reformatting_calls_for_merge_f1(p377_cell_classification_f1)
p409_cell_classification_f1 = reformatting_calls_for_merge_f1(p409_cell_classification_f1)
p318_cell_classification_silwidth = reformatting_calls_for_merge_silwidth(p318_cell_classification_silwidth)
p377_cell_classification_silwidth = reformatting_calls_for_merge_silwidth(p377_cell_classification_silwidth)
p409_cell_classification_silwidth = reformatting_calls_for_merge_silwidth(p409_cell_classification_silwidth)

cell_classification_f1 = do.call("rbind", list(p318_cell_classification_f1,
                                               p377_cell_classification_f1,
                                               p409_cell_classification_f1))
cell_classification_silwidth = do.call("rbind", list(p318_cell_classification_silwidth,
                                                     p377_cell_classification_silwidth,
                                                     p409_cell_classification_silwidth))

cell_classification = merge(cell_classification_f1, cell_classification_silwidth, by = c("barcode", "gi_score", "correlation_score"))

cell_classification = subset(cell_classification, barcode %in% pdo@meta.data$barcode)

cell_classification$coarse_anno[cell_classification$barcode %in% subset(pdo@meta.data, coarse_anno=="Epithelial")$barcode] <- "Epithelial"
cell_classification$coarse_anno[cell_classification$barcode %in% subset(pdo@meta.data, coarse_anno=="Immune")$barcode] <- "Immune"
cell_classification$coarse_anno[cell_classification$barcode %in% subset(pdo@meta.data, coarse_anno=="Stromal")$barcode] <- "Stromal"

# adding a column to see how well f1 vs silwidth calls match
cell_classification$match = rep("mismatch")
cell_classification$match[cell_classification$normal_cell_call_f1 == cell_classification$normal_cell_call_silwidth] <- "match"

# adding scores and calls to meta data
for (barcode in unique(cell_classification$barcode)) {
  # find gi_score, correlation_score, f1 call, silwidth call
  gi_score <- cell_classification[cell_classification$barcode == barcode, 'gi_score']
  correlation_score <- cell_classification[cell_classification$barcode == barcode, 'correlation_score']
  f1_call <- cell_classification[cell_classification$barcode == barcode, 'normal_cell_call_f1']
  silwidth_call <- cell_classification[cell_classification$barcode == barcode, 'normal_cell_call_silwidth']
  
  # assign to the corresponding rows in destination 
  pdo@meta.data$gi_score[pdo@meta.data$barcode == barcode] <- gi_score
  pdo@meta.data$correlation_score[pdo@meta.data$barcode == barcode] <- correlation_score
  pdo@meta.data$f1_call[pdo@meta.data$barcode == barcode] <- f1_call
  pdo@meta.data$silwidth_call[pdo@meta.data$barcode == barcode] <- silwidth_call
}

plot1 <- scCustomize::DimPlot_scCustom(pdo, group.by = "f1_call", colors_use = c("red", "blue"), reduction = "umap.scvi")
plot2 <- scCustomize::DimPlot_scCustom(pdo, group.by = "silwidth_call", colors_use = c("red", "blue"), reduction = "umap.scvi")
plot_grid(plot1, plot2)

# 5.4 comparing calls to cosmic and canonical marker sets ====
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
                         gobp_mamm_epi_prolif$GOBP_POSITIVE_REGULATION_OF_MAMMARY_GLAND_EPITHELIAL_CELL_PROLIFERATION)

pdo <- Seurat::AddModuleScore(pdo, features = bc_tumor_genesets)

pdo@meta.data <- pdo@meta.data %>% dplyr::rename("cgc_exp" = `Cluster1`,
                                                 "gobp_mamm_epi_prolif" = `Cluster2`) 

scCustomize::FeaturePlot_scCustom(pdo, features = c("correlation_score", "cgc_exp", 
                                                    "gobp_mamm_epi_prolif"), 
                                  figure_plot = T, order = T, num_columns = 2, reduction = "umap.scvi")

ggplot(subset(pdo@meta.data, !is.na(silwidth_call) & patientID =="318"), aes(gi_score)) + 
  geom_vline(xintercept = 0.1916, col = "black", lwd = 1) +
  geom_histogram(bins = 2000, col = patient_colors[1]) +
  theme_pubr()

plot1 <- ggplot(subset(pdo@meta.data, !is.na(silwidth_call) & patientID =="318"), aes(gi_score)) + 
  geom_vline(xintercept = 0.1916, col = "black", lwd = 1) +
  geom_density(fill = patient_colors[1], col = patient_colors[1], alpha = 0.5) +
  theme_pubr()
plot2 <- ggplot(subset(pdo@meta.data, !is.na(silwidth_call) & patientID =="318"), aes(correlation_score)) + 
  geom_vline(xintercept = 0.695, col = "black", lwd = 1) +
  geom_density(fill = patient_colors[1], col = patient_colors[1], alpha = 0.5) +
  theme_pubr()
plot3 <- ggplot(subset(pdo@meta.data, !is.na(silwidth_call) & patientID =="377"), aes(gi_score)) + 
  geom_vline(xintercept = 0.4873, col = "black", lwd = 1) +
  geom_density(fill = patient_colors[2], col = patient_colors[2], alpha = 0.5) +
  theme_pubr()
plot4 <- ggplot(subset(pdo@meta.data, !is.na(silwidth_call) & patientID =="377"), aes(correlation_score)) + 
  geom_vline(xintercept = 0.12, col = "black", lwd = 1) +
  geom_density(fill = patient_colors[2], col = patient_colors[2], alpha = 0.5) +
  theme_pubr()
plot5 <- ggplot(subset(pdo@meta.data, !is.na(silwidth_call) & patientID =="409"), aes(gi_score)) + 
  geom_vline(xintercept = 0.062, col = "black", lwd = 1) +
  geom_density(fill = patient_colors[3], col = patient_colors[3], alpha = 0.5) +
  theme_pubr()
plot6 <- ggplot(subset(pdo@meta.data, !is.na(silwidth_call) & patientID =="409"), aes(correlation_score)) + 
  geom_vline(xintercept = 0.41, col = "black", lwd = 1) +
  geom_density(fill = patient_colors[3], col = patient_colors[3], alpha = 0.5) +
  theme_pubr()

plot_grid(plot1, plot2, plot3, plot4, plot5, plot6, ncol = 2)

pdo@meta.data$manual_call[pdo@meta.data$coarse_anno=="Epithelial"] <- "cancer"
pdo@meta.data$manual_call[pdo@meta.data$patientID=="318" & pdo@meta.data$gi_score < 0.1916 & pdo@meta.data$correlation_score < 0.695] <- "normal"
pdo@meta.data$manual_call[pdo@meta.data$patientID=="377" & pdo@meta.data$gi_score < 0.4873 & pdo@meta.data$correlation_score < 0.12] <- "normal"
pdo@meta.data$manual_call[pdo@meta.data$patientID=="409" & pdo@meta.data$gi_score < 0.062 & pdo@meta.data$correlation_score < 0.41] <- "normal"
pdo@meta.data$manual_call[pdo@meta.data$coarse_anno == "Stromal"] <- NA

DimPlot_scCustom(pdo, group.by = "manual_call", reduction = "umap.scvi", colors_use = c("red", "blue"))

# 5.5 save progress ====
saveRDS(pdo, file = file.path(save_path, "pdo_integrated.rds"))