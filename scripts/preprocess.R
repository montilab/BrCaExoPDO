# 1.0 load packages and initial data ====
source("../scripts/utils.R")
save_path <- "path/to/save/dir/" # please update file path to save in your desired location
options(Seurat.object.assay.version = "v5")
basefile_path <- "/path/to/files" # please download files from GEO and update this line to point to your directory

data_dirs <- c(
  "p318C" = file.path(basefile_path, "318C"),
  "p318N" = file.path(basefile_path, "318N"),
  "p318T" = file.path(basefile_path, "318T"),
  "p377C" = file.path(basefile_path, "377C"),
  "p377N" = file.path(basefile_path, "377N"),
  "p377T" = file.path(basefile_path, "377T"),
  "p409C" = file.path(basefile_path, "409C"),
  "p409N" = file.path(basefile_path, "409N"),
  "p409T" = file.path(basefile_path, "409T"),
  "p409C2" = file.path(basefile_path, "409C2"),
  "p409N2" = file.path(basefile_path, "409N2"),
  "p409T2" = file.path(basefile_path, "409T2")
  )

raw_objs <- lapply(names(data_dirs), function(name) {
  create_seurat_object(data_dirs[[name]], "pdo")
})

names(raw_objs) <- names(data_dirs)
sapply(raw_objs, ncol) # print number of cells per obj

# 1.1 doublet detection ====
expected_rates <- c(
  p318C = 0.015, # from 10X; expected doublet rate when recovering ~2000 cells is ~1.5%
  p318N = 0.016, # expected doublet rate when recovering ~2400 cells is ~1.6%
  p318T = 0.023, # expected doublet rate when recovering ~2600 cells is ~2.3%
  p377C = 0.046, # expected doublet rate when recovering ~6152 cells is ~4.6%
  p377N = 0.061, # expected doublet rate when recovering ~7891 cells is ~6.1%
  p377T = 0.046, # expected doublet rate when recovering ~5658 cells is ~4.6% 
  p409C = 0.004, # expected doublet rate when recovering ~71 cells is ~0.4%
  p409N = 0.004, # expected doublet rate when recovering ~71 cells is ~0.4%
  p409T = 0.004, # expected doublet rate when recovering ~71 cells is ~0.4%
  p409C2 = 0.004, # expected doublet rate when recovering ~71 cells is ~0.4%
  p409N2 = 0.004, # expected doublet rate when recovering ~71 cells is ~0.4%
  p409T2 = 0.004 # expected doublet rate when recovering ~71 cells is ~0.4%
)

raw_objs <- lapply(names(raw_objs), function(name) {
  list(object = raw_objs[[name]], expected_rate = expected_rates[[name]])
})

doublet_detected_objs <- lapply(seq_along(raw_objs), function(x) {
  obj <- raw_objs[[x]]$object
  expected_rate <- raw_objs[[x]]$expected_rate
  cat(paste("Processing object", x, "of", length(raw_objs), "\n"))
  processed_obj <- process_and_doublet_finder(obj, expected_rate)
  cat(paste("Finished processing object", x, "with optimal pK value\n\n"))
  return(processed_obj)
})

names(doublet_detected_objs) = c("p318C", "p318N", "p318T", "p377C", "p377N", "p377T",
                                 "p409C", "p409N", "p409T", "p409C2", "p409N2", "p409T2") # double checked the order and this is correct

# checking results
doublet_dimplot_check(doublet_detected_objs)

for (obj in doublet_detected_objs) {
  print_metadata_table(obj, "DF.classifications_")
}

# adding doublets as column in real data
raw_doublet_detected_objs <- update_doublets(raw = raw_objs, doublet = doublet_detected_objs)

# 1.2 QC visualizations ====
for (i in seq_along(raw_doublet_detected_objs)) {
  raw_doublet_detected_objs[[i]] <- qc_metrics(raw_doublet_detected_objs[[i]]$object)
}

for (i in seq_along(raw_doublet_detected_objs)) {
  qc_plots(raw_doublet_detected_objs[[i]], names(raw_doublet_detected_objs)[i])
}

# 1.3 low quality cell filtering ====
subset_objs <- mapply(qc_filtering, raw_doublet_detected_objs, names(raw_doublet_detected_objs), SIMPLIFY = FALSE)

# number of cells in p318C : 2007 
# number of cells in p318N : 2384 
# number of cells in p318T : 2562 
# number of cells in p377C2 : 5396 
# number of cells in p377N2 : 6915 
# number of cells in p377T2 : 4953 
# number of cells in p409C : 67 
# number of cells in p409N : 84 
# number of cells in p409T : 54 
# number of cells in p409C2 : 97 
# number of cells in p409N2 : 40 
# number of cells in p409T2 : 40 
# 1.4 save progress ####
saveRDS(subset_objs, file=file.path(save_path, "subset_objs.rds"))