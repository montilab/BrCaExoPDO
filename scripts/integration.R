# 4.0 load packages and data ====
source("../scripts/utils.R")
save_path <- "path/to/save/dir/" # please update file path to save in your desired location
pdo <- readRDS(file.path(save_path, "pdo_unintegrated.rds"))
pdo_epi <- readRDS(file.path(save_path, "pdo_epi_unintegrated.rds"))
pdo_stro <- readRDS(file.path(save_path, "pdo_stro_unintegrated.rds"))
pdo_imm <- readRDS(file.path(save_path, "pdo_imm_unintegrated.rds"))

# 4.1 one line integration ####
pdo <- SCTransform(pdo, vars.to.regress = "percent_mito")
pdo <- RunPCA(pdo)
options(future.globals.maxSize = 8000 * 1024^2)

# cca integration
pdo <- IntegrateLayers(object = pdo, method = CCAIntegration, normalization.method = "SCT", k.weight = 10,
                       orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)
write.csv(pdo@reductions$integrated.cca@cell.embeddings,
          file.path(save_path, "/pdo_cca_embedding_mitoregress.csv"))

# rpca integration
pdo <- IntegrateLayers(object = pdo, method = RPCAIntegration, normalization.method = "SCT", k.weight = 10,
                       orig.reduction = "pca", new.reduction = "integrated.rpca", verbose = FALSE)
write.csv(pdo@reductions$integrated.rpca@cell.embeddings,
          file.path(save_path, "pdo_rpca_embedding_mitoregress.csv"))

# harmony
pdo <- IntegrateLayers(object = pdo, method = HarmonyIntegration, normalization.method = "SCT", k.weight = 10,
                       orig.reduction = "pca", new.reduction = "harmony", verbose = FALSE)
write.csv(pdo@reductions$harmony@cell.embeddings,
          file.path(save_path, "pdo_harmony_embedding_mitoregress.csv"))

# scVI integration
reticulate::use_condaenv("path_to/scvi_env/", required = TRUE) # please update directory path to your scVI conda environment
reticulate::py_config()
# remotes::install_github(repo="satijalab/seurat-wrappers", ref = remotes::github_pull(184)) # had trouble with code not working, updating seuratwrappers to this pull request fixed it

pdo <- IntegrateLayers(object = pdo, method = scVIIntegration, k.weight = 10, normalization.method = "SCT",
                       orig.reduction = "pca", new.reduction = "integrated.scvi", 
                       conda_env = "path_to/scvi_env/", verbose = FALSE) # please update directory path to your scVI conda environment
write.csv(pdo@reductions$integrated.scvi@cell.embeddings,
          file.path(save_path, "pdo_scvi_embedding_mitoregress.csv"))

# 4.2 benchmarking via scIB ####
# https://scib-metrics.readthedocs.io/en/latest/index.html
# this was built for python - run in r via reticulate

# checking/specifying environment
reticulate::use_condaenv("path_to/scib_env/", required = TRUE) # please update directory path to your scVI conda environment
reticulate::py_config()

# saving seurat object as h5ad (needed for python work)
loompy <- reticulate::import('loompy')
pdo_adata <- pdo
pdo_adata[["RNA3"]] <- as(object = pdo_adata[["RNA"]], Class = "Assay")
DefaultAssay(pdo_adata) <- "RNA3"
pdo_adata[["RNA"]] <- NULL
pdo_adata <- RenameAssays(object = pdo_adata, RNA3 = 'RNA')
sceasy::convertFormat(pdo_adata, from="seurat", to="anndata",
                      outFile=file.path(save_path, 'pdo_mitoregress.h5ad'))

# import necessary python modules
sc <- reticulate::import("scanpy")
np <- reticulate::import("numpy")
pd <- reticulate::import("pandas")
rich <- reticulate::import("rich")
plt <- reticulate::import("matplotlib.pyplot")
Benchmarker <- reticulate::import("scib_metrics.benchmark")$Benchmarker
BioConservation <- reticulate::import("scib_metrics.benchmark")$BioConservation
BatchCorrection <- reticulate::import("scib_metrics.benchmark")$BatchCorrection

# load the data 
adata <- sc$read_h5ad(file.path(save_path, "pdo_mitoregress.h5ad"))

# update embeddings in anndata object
adata$obsm['harmony'] <- adata$obsm["X_harmony"]
adata$obsm['rpca'] <- adata$obsm["X_integrated.rpca"]
adata$obsm["scVI"] <- adata$obsm["X_integrated.scvi"]
adata$obsm['cca'] <- adata$obsm["X_integrated.cca"]
adata$obsm["Unintegrated"] <- adata$obsm["X_pca"]

# define metrics
biocons <- BioConservation(
  isolated_labels=TRUE, 
  nmi_ari_cluster_labels_leiden=FALSE, 
  nmi_ari_cluster_labels_kmeans=TRUE, 
  silhouette_label=TRUE, 
  clisi_knn=TRUE
)
batchcons <- BatchCorrection(
  silhouette_batch=TRUE, 
  ilisi_knn=TRUE, 
  kbet_per_label=TRUE, 
  graph_connectivity=TRUE, 
  pcr_comparison=TRUE
)

# benchmark
bm <- Benchmarker(adata,
                  batch_key="group", label_key="singleR_labels",
                  bio_conservation_metrics=biocons, batch_correction_metrics=batchcons,
                  embedding_obsm_keys=c("Unintegrated", "scVI", "harmony", "rpca", "cca"),
                  pre_integrated_embedding_obsm_key="X_pca", n_jobs=-1)

adata$obs[['group']] <- as.character(adata$obs[['group']])
adata$obs[['singleR_labels']] <- as.character(adata$obs[['singleR_labels']])
bm$benchmark()

# plot and save results
bm$plot_results_table(min_max_scale=FALSE, save_dir=file.path(save_path, "benchmarking/"))
df <- bm$get_results(min_max_scale=FALSE)
df <- df %>% dplyr::mutate(across(where(is.list), ~sapply(., function(x) paste(x, collapse = ", "))))
write.csv(df, file = file.path(save_path, "bm_singleR_mitoregress.csv"), row.names = FALSE)

# 4.3 calculating neighbors for each method ####
pdo <- FindNeighbors(pdo, reduction = "pca", dims = 1:30)
do_clustree_integrated(obj = pdo, res_seq = c(0.1, 1, 0.1))
pdo <- FindClusters(pdo, resolution = 0.5, cluster.name = "unintegrated_clusters")
pdo <- RunUMAP(pdo, dims = 1:30, reduction = "pca", reduction.name = "umap.pca")

pdo <- FindNeighbors(pdo, reduction = "integrated.cca", dims = 1:30)
do_clustree_integrated(obj = pdo, res_seq = c(0.1, 1, 0.1))
pdo <- FindClusters(pdo, resolution = 0.4, cluster.name = "cca_clusters")
pdo <- RunUMAP(pdo, dims = 1:30, reduction = "integrated.cca", reduction.name = "umap.cca")

pdo <- FindNeighbors(pdo, reduction = "integrated.rpca", dims = 1:30)
do_clustree_integrated(obj = pdo, res_seq = c(0.1, 1, 0.1))
pdo <- FindClusters(pdo, resolution = 0.5, cluster.name = "rpca_clusters")
pdo <- RunUMAP(pdo, dims = 1:30, reduction = "integrated.rpca", reduction.name = "umap.rpca")

pdo <- FindNeighbors(pdo, reduction = "integrated.scvi", dims = 1:30)
do_clustree_integrated(obj = pdo, res_seq = c(0.1, 1, 0.1))
pdo <- FindClusters(pdo, resolution = 0.5, cluster.name = "scvi_clusters")
pdo <- RunUMAP(pdo, dims = 1:30, reduction = "integrated.scvi", reduction.name = "umap.scvi")

pdo <- FindNeighbors(pdo, reduction = "harmony", dims = 1:30)
do_clustree_integrated(obj = pdo, res_seq = c(0.1, 1, 0.1))
pdo <- FindClusters(pdo, resolution = 0.4, cluster.name = "harmony_clusters")
pdo <- RunUMAP(pdo, dims = 1:30, reduction = "harmony", reduction.name = "umap.harmony")

pdo <- JoinLayers(pdo, assay = "RNA")

# 4.4 visualizations ####
# visualizations to compare integration methods
integration_simpleplots(pdo, reduction = "umap.pca", cluster_labels = "unintegrated_clusters")
integration_simpleplots(pdo, reduction = "umap.cca", cluster_labels = "cca_clusters")
integration_simpleplots(pdo, reduction = "umap.rpca", cluster_labels = "rpca_clusters")
integration_simpleplots(pdo, reduction = "umap.harmony", cluster_labels = "harmony_clusters")
integration_simpleplots(pdo, reduction = "umap.scvi", cluster_labels = "scvi_clusters")

# 4.5 coarse annotation split pre integration ====
options(future.globals.maxSize = 8000 * 1024^2)
reticulate::use_condaenv("path_to/scvi_env/", required = TRUE) # please update directory path to your scVI conda environment
reticulate::py_config()

# sct and integration processing
pdo_epi <- SCTransform(pdo_epi, verbose = FALSE, vars.to.regress = "percent_mito")
pdo_epi <- RunPCA(pdo_epi)
pdo_epi <- IntegrateLayers(object = pdo_epi, method = scVIIntegration, k.weight = 10, normalization.method = "SCT",
                           orig.reduction = "pca", new.reduction = "integrated.scvi", 
                           conda_env = "path_to/scvi_env/", verbose = FALSE) # please update directory path to your scVI conda environment
write.csv(pdo_epi@reductions$integrated.scvi@cell.embeddings,
          file.path(save_path, "epi_scvi_embedding.csv"))

pdo_stro <- SCTransform(pdo_stro, verbose = FALSE, min_cells = 2, vars.to.regress = "percent_mito")
pdo_stro <- RunPCA(pdo_stro)
pdo_stro <- IntegrateLayers(object = pdo_stro, method = scVIIntegration, k.weight = 10, normalization.method = "SCT",
                            orig.reduction = "pca", new.reduction = "integrated.scvi", 
                            conda_env = "path_to/scvi_env/", verbose = FALSE) # please update directory path to your scVI conda environment
write.csv(pdo_stro@reductions$integrated.scvi@cell.embeddings,
          file.path(save_path, "stro_scvi_embedding.csv"))

# counts.9 doesn't have any immune cells
pdo_imm@assays$RNA@layers$counts.9 <- NULL
pdo_imm@assays$RNA$counts.9 <- NULL
pdo_imm <- SCTransform(pdo_imm, verbose = FALSE, min_cells = 2, vars.to.regress = "percent_mito")
pdo_imm <- RunPCA(pdo_imm)
pdo_imm <- IntegrateLayers(object = pdo_imm, method = scVIIntegration, k.weight = 10, normalization.method = "SCT",
                           orig.reduction = "pca", new.reduction = "integrated.scvi", 
                           conda_env = "path_to/scvi_env/", verbose = FALSE) # please update directory path to your scVI conda environment
write.csv(pdo_imm@reductions$integrated.scvi@cell.embeddings,
          file.path(save_path, "imm_scvi_embedding.csv"))

# finding right resolution for each object
pdo_epi <- FindNeighbors(pdo_epi, reduction = "pca", dims = 1:30)
do_clustree_integrated(obj = pdo_epi, res_seq = c(0.1, 1, 0.1))
pdo_epi <- FindClusters(pdo_epi, resolution = 0.6, cluster.name = "unintegrated_clusters")
pdo_epi <- RunUMAP(pdo_epi, dims = 1:30, reduction = "pca", reduction.name = "umap.pca")
pdo_epi <- FindNeighbors(pdo_epi, reduction = "integrated.scvi", dims = 1:30)
do_clustree_integrated(obj = pdo_epi, res_seq = c(0.1, 1, 0.1))
pdo_epi <- FindClusters(pdo_epi, resolution = 0.5, cluster.name = "scvi_clusters")
pdo_epi <- RunUMAP(pdo_epi, dims = 1:30, reduction = "integrated.scvi", reduction.name = "umap.scvi")

pdo_stro <- FindNeighbors(pdo_stro, reduction = "pca", dims = 1:30)
do_clustree_integrated(obj = pdo_stro, res_seq = c(0.1, 1, 0.1))
pdo_stro <- FindClusters(pdo_stro, resolution = 0.4, cluster.name = "unintegrated_clusters")
pdo_stro <- RunUMAP(pdo_stro, dims = 1:30, reduction = "pca", reduction.name = "umap.pca")
pdo_stro <- FindNeighbors(pdo_stro, reduction = "integrated.scvi", dims = 1:30)
do_clustree_integrated(obj = pdo_stro, res_seq = c(0.1, 1, 0.1))
pdo_stro <- FindClusters(pdo_stro, resolution = 0.7, cluster.name = "scvi_clusters")
pdo_stro <- RunUMAP(pdo_stro, dims = 1:30, reduction = "integrated.scvi", reduction.name = "umap.scvi")

pdo_imm <- FindNeighbors(pdo_imm, reduction = "pca", dims = 1:30)
do_clustree_integrated(obj = pdo_imm, res_seq = c(0.1, 1, 0.1))
pdo_imm <- FindClusters(pdo_imm, resolution = 1.2, cluster.name = "unintegrated_clusters")
pdo_imm <- RunUMAP(pdo_imm, dims = 1:30, reduction = "pca", reduction.name = "umap.pca")
pdo_imm <- FindNeighbors(pdo_imm, reduction = "integrated.scvi", dims = 1:30)
do_clustree_integrated(obj = pdo_imm, res_seq = c(0.1, 1, 0.1))
pdo_imm <- FindClusters(pdo_imm, resolution = 1.2, cluster.name = "scvi_clusters")
pdo_imm <- RunUMAP(pdo_imm, dims = 1:30, reduction = "integrated.scvi", reduction.name = "umap.scvi")

# joining layers
pdo_epi <- JoinLayers(pdo_epi, assay = "RNA")
pdo_stro <- JoinLayers(pdo_stro, assay = "RNA")
pdo_imm <- JoinLayers(pdo_imm, assay = "RNA")

# visualizing 
plot1 <- scCustomize::DimPlot_scCustom(pdo_epi, reduction = "umap.pca", figure_plot = T, group.by = "unintegrated_clusters")
plot2 <- scCustomize::DimPlot_scCustom(pdo_epi, reduction = "umap.scvi", figure_plot = T, group.by = "scvi_clusters")
plot_grid(plot1, plot2)

plot1 <- scCustomize::DimPlot_scCustom(pdo_stro, reduction = "umap.pca", figure_plot = T, group.by = "unintegrated_clusters")
plot2 <- scCustomize::DimPlot_scCustom(pdo_stro, reduction = "umap.scvi", figure_plot = T, group.by = "scvi_clusters")
plot_grid(plot1, plot2)

plot1 <- scCustomize::DimPlot_scCustom(pdo_imm, reduction = "umap.pca", figure_plot = T, group.by = "unintegrated_clusters")
plot2 <- scCustomize::DimPlot_scCustom(pdo_imm, reduction = "umap.scvi", figure_plot = T, group.by = "scvi_clusters")
plot_grid(plot1, plot2)

# checking my subsets
FeaturePlot(pdo_epi, features = c("EPCAM", "KRT8", "KRT18", # epithelial
                                  "PECAM1", "VWF", "CDH5", # endothelial 
                                  "PTPRC", "CD33", "CSF1R", # general immune/myeloid
                                  "FAP", "COL1A1", "PDPN"), # fibroblast
            ncol = 3, reduction = "umap.scvi")

FeaturePlot(pdo_stro, features = c("EPCAM", "KRT8", "KRT18", # epithelial
                                   "PECAM1", "VWF", "CDH5", # endothelial 
                                   "PTPRC", "CD33", "CSF1R", # general immune/myeloid
                                   "FAP", "COL1A1", "PDPN"), # fibroblast
            ncol = 3, reduction = "umap.scvi")

FeaturePlot(pdo_imm, features = c("EPCAM", "KRT8", "KRT18", # epithelial
                                  "PECAM1", "VWF", "CDH5", # endothelial 
                                  "PTPRC", "CD33", "CSF1R", # general immune/myeloid
                                  "FAP", "COL1A1", "PDPN"), # fibroblast
            ncol = 3, reduction = "umap.scvi")

# 4.6 coarse annotation split post integration ====
pdo_epi_postintegration = subset(x = pdo, subset = coarse_anno == "Epithelial")
pdo_stro_postintegration = subset(x = pdo, subset = coarse_anno == "Stromal")
pdo_imm_postintegration = subset(x = pdo, subset = coarse_anno == "Immune")

pdo_epi_postintegration = quick_process(pdo_epi_postintegration)
pdo_stro_postintegration = quick_process(pdo_stro_postintegration)
pdo_imm_postintegration = quick_process(pdo_imm_postintegration)

# finding right resolution for each object
do_clustree_integrated(obj = pdo_epi_postintegration, res_seq = c(0.1, 1, 0.1))
do_clustree_integrated(obj = pdo_stro_postintegration, res_seq = c(0.1, 1, 0.1))
do_clustree_integrated(obj = pdo_imm_postintegration, res_seq = c(0.1, 1, 0.1))

pdo_epi_postintegration <- FindClusters(pdo_epi_postintegration, resolution = 0.5) 
pdo_epi_postintegration <- RunUMAP(pdo_epi_postintegration, dims = 1:20)

pdo_stro_postintegration <- FindClusters(pdo_stro_postintegration, resolution = 0.6)
pdo_stro_postintegration <- RunUMAP(pdo_stro_postintegration, dims = 1:20)

pdo_imm_postintegration <- FindClusters(pdo_imm_postintegration, resolution = 1)
pdo_imm_postintegration <- RunUMAP(pdo_imm_postintegration, dims = 1:20)

scCustomize::DimPlot_scCustom(pdo_epi_postintegration, reduction = "umap", figure_plot = T)
scCustomize::DimPlot_scCustom(pdo_stro_postintegration, reduction = "umap", figure_plot = T)
scCustomize::DimPlot_scCustom(pdo_imm_postintegration, reduction = "umap", figure_plot = T)

# 4.7 saving objects ####
saveRDS(pdo, file = file.path(save_path, "pdo_integrated.rds"))
saveRDS(pdo_epi, file = file.path(save_path, "pdo_epi_integrated.rds"))
saveRDS(pdo_stro, file = file.path(save_path, "pdo_stro_integrated.rds"))
saveRDS(pdo_imm, file = file.path(save_path, "pdo_imm_integrated.rds"))
saveRDS(pdo_epi_postintegration, file = file.path(save_path, "pdo_epi_split_post_integrated.rds"))
saveRDS(pdo_stro_postintegration, file = file.path(save_path, "pdo_stro_split_post_integrated.rds"))
saveRDS(pdo_imm_postintegration, file = file.path(save_path, "pdo_imm_split_post_integrated.rds"))