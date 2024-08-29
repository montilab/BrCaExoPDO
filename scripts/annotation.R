# 2.0 load packages and data ====
source("../scripts/utils.R")
save_path <- "path/to/save/dir/" # please update file path to save in your desired location
readRDS(file=file.path(save_path, "subset_objs.rds")) 

# 2.1 creating reference ====
# paper: https://www.nature.com/articles/s41588-021-00911-1
# please download the data from Gene Expression Omnibus under accession number GSE176078 and update file paths as needed
swarbrick_metadata <- read.csv("/path_to/metadata.csv")
swarbrick_barcodes <- read_tsv("/path_to/count_matrix_barcodes.tsv", col_names = F)
swarbrick_genes <- read_tsv("/path_to/count_matrix_genes.tsv", col_names = F)
swarbrick__sparse <- Matrix::readMM("/path_to/count_matrix_sparse.mtx")

rownames(swarbrick__sparse) = swarbrick_genes$X1
colnames(swarbrick__sparse) = swarbrick_barcodes$X1

swarbrick <- SingleCellExperiment::SingleCellExperiment(assays=list(logcounts=swarbrick__sparse),
                                                        rowData = swarbrick_genes, colData=swarbrick_barcodes,
                                                        metadata=swarbrick_metadata)
# 2.2 running parallel workflow ====
for (obj_name in names(subset_objs)) {
  obj <- sct_obj_list[[obj_name]]
  
  pred <- SingleR(
    test = obj@assays$RNA$counts, # can be on raw or log normalized - but github of singleR did specify not sct normalized
    ref = swarbrick@assays@data$logcounts,
    labels = swarbrick@metadata$celltype_minor,
    de.method = "wilcox",
    fine.tune = TRUE
  )
  
  output_file <- paste0(save_path, "predicted_celltypes_swarbrick_", obj_name, ".csv")
  
  write.csv(pred, output_file, row.names = TRUE)
}
# 2.3 adding labels to seurat objects ====
#load in matrices - saved in github
pred_p318C <- read.csv("../results/annotations/predicted_celltypes_p318C.csv")
pred_p318N <- read.csv("../results/annotations/predicted_celltypes_p318N.csv")
pred_p318T <- read.csv("../results/annotations/predicted_celltypes_p318T.csv")
pred_p377C <- read.csv("../results/annotations/predicted_celltypes_p377C.csv")
pred_p377N <- read.csv("../results/annotations/predicted_celltypes_p377N.csv")
pred_p377T <- read.csv("../results/annotations/predicted_celltypes_p377T.csv")
pred_p409C <- read.csv("../results/annotations/predicted_celltypes_p409C.csv")
pred_p409N <- read.csv("../results/annotations/predicted_celltypes_p409N.csv")
pred_p409T <- read.csv("../results/annotations/predicted_celltypes_p409T.csv")
pred_p409C2 <- read.csv("../results/annotations/predicted_celltypes_p409C2.csv")
pred_p409N2 <- read.csv("../results/annotations/predicted_celltypes_p409N2.csv")
pred_p409T2 <- read.csv("../results/annotations/predicted_celltypes_p409T2.csv")

# apply the function to each seurat object in the list
anno_obj_list <- setNames(
  lapply(names(subset_objs), function(name) {
    seurat_obj <- subset_objs[[name]]
    predicted_celltypes <- get(paste0("pred_", name))
    add_singleR_metadata(seurat_obj, predicted_celltypes)
  }),
  names(subset_objs)
)

# 2.4 save progress ####
saveRDS(anno_obj_list, file=file.path(save_path, "anno_obj_list.rds"))