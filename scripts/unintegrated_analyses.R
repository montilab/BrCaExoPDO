# 7.0 load packages and data ====
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

for (category in names(table(pdo@meta.data$finalized_fine_anno))) {
  pdo_unintegrated@meta.data$finalized_fine_anno[pdo_unintegrated@meta.data$barcode %in% subset(pdo@meta.data, finalized_fine_anno == category)$barcode] <- category
}
for (category in names(table(pdo@meta.data$finalized_broad_anno))) {
  pdo_unintegrated@meta.data$finalized_broad_anno[pdo_unintegrated@meta.data$barcode %in% subset(pdo@meta.data, finalized_broad_anno == category)$barcode] <- category
}

#loading msigdb files
# please export gmt files from msigdb: https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
all_msigdb = read.gmt("path_to/msigdb.v2022.1.Hs.symbols.gmt")
hallmarks = read.gmt("path_to/h.all.v2022.1.Hs.symbols.gmt")
c2_cp = read.gmt("path_to/c2.cp.v2022.1.Hs.symbols.gmt")
c5_go = read.gmt("path_to/c5.go.v2022.1.Hs.symbols.gmt")
c6_onco = read.gmt("path_to/c6.all.v2022.1.Hs.symbols.gmt")
c7_imm = read.gmt("path_to/c7.all.v2022.1.Hs.symbols.gmt")

msigdb_of_interest = list(all = all_msigdb, hallmarks = hallmarks, 
                          c2 = c2_cp, c5 = c5_go, c6 = c6_onco, c7 = c7_imm)

# 7.1 all cell differential expression ####
Idents(object = pdo_unintegrated) <- "treatment"

# setting up marker genes
# n.c.allcells.markers_fullreps <- FindMarkers(pdo_unintegrated_fullreps, ident.1 = "NET", ident.2 = "UT", test.use="MAST",  latent.vars="patient_replicate")
# t.c.allcells.markers_fullreps <- FindMarkers(pdo_unintegrated_fullreps, ident.1 = "DET", ident.2 = "UT", test.use="MAST", latent.vars="patient_replicate")
# t.n.allcells.markers_fullreps <- FindMarkers(pdo_unintegrated_fullreps, ident.1 = "DET", ident.2 = "NET", test.use="MAST", latent.vars="patient_replicate")
n.c.allcells = read.csv("../results/analyses/n.c.allcells.markers.csv")
t.c.allcells = read.csv("../results/analyses/t.c.allcells.markers.csv")
t.n.allcells = read.csv("../results/analyses/t.n.allcells.markers.csv")

markers_list = list(n.c.allsamples = n.c.allcells, 
                    t.c.allsamples = t.c.allcells, 
                    t.n.allsamples = t.n.allcells)

# running setup functions
combined_volcano_plot(markers_list, ncol = 3)
allcells_for_gsea = create_gsea_vectors(markers_list)
fgsea_results_allcells = run_through_fgsea(msigdb_list = msigdb_of_interest, 
                                           named_degs_list = allcells_for_gsea) # this takes a few minutes

saveRDS(fgsea_results_allcells, file.path(save_path, "fgsea_results_allcells.rds"))

generate_gsea_table_plots(treatment_comparison = "t.n", search_string = "hallmarks",
                          msigdb_list = msigdb_of_interest, 
                          named_degs_list = allcells_for_gsea, 
                          fgsea_res = fgsea_results_allcells)

# 7.2 all cell proportion shift ####
prop_calcs = pdo@meta.data %>%
  dplyr::group_by(treatment, finalized_broad_anno) %>%
  dplyr::summarise(count = n())
prop_calcs$finalized_broad_anno = factor(prop_calcs$finalized_broad_anno, levels = broadcelltype_order_rainbow)
ggplot(prop_calcs, aes(x=treatment, y=count))+
  geom_bar(aes(fill=finalized_broad_anno), stat = "identity", position = "fill")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")+
  labs(x="", y = "Proportion", fill = "")+
  scale_fill_manual(values = broadcelltype_colors_rainbow)

# chi square
prop_stat = pdo@meta.data %>%
  dplyr::group_by(treatment, finalized_broad_anno) %>%
  dplyr::summarise(count = n()) %>%
  spread(finalized_broad_anno, count, fill = 0)
prop_stat = as.data.frame(prop_stat)
rownames(prop_stat) = prop_stat$treatment
prop_stat$treatment <- NULL
chisq.test(prop_stat) #p=2.2e-16
chisq.posthoc.test(prop_stat, method = "bonferroni")

# 7.3 epi differential expression ####
Idents(object = pdo_epi) <- "treatment"

# setting up marker genes
# n.c.epi.markers <- FindMarkers(pdo_epi, ident.1 = "NET", ident.2 = "UT", test.use="MAST", latent.vars="patient_replicate")
# t.c.epi.markers <- FindMarkers(pdo_epi, ident.1 = "DET", ident.2 = "UT", test.use="MAST", latent.vars="patient_replicate")
# t.n.epi.markers <- FindMarkers(pdo_epi, ident.1 = "DET", ident.2 = "NET", test.use="MAST", latent.vars="patient_replicate")
n.c.epi = read.csv("../results/analyses/n.c.epi.markers.csv")
t.c.epi = read.csv("../results/analyses/t.c.epi.markers.csv")
t.n.epi = read.csv("../results/analyses/t.n.epi.markers.csv")

markers_list = list(n.c.allsamples = n.c.epi, 
                    t.c.allsamples = t.c.epi, 
                    t.n.allsamples = t.n.epi)

combined_volcano_plot(markers_list, ncol = 3)
epi_for_gsea = create_gsea_vectors(markers_list)
fgsea_results_epi = run_through_fgsea(msigdb_list = msigdb_of_interest, 
                                      named_degs_list = epi_for_gsea) # this takes a few minutes

saveRDS(fgsea_results_epi, file.path(save_path, "fgsea_results_epi.rds"))

generate_gsea_table_plots(treatment_comparison = "t.n", search_string = "c5",
                          msigdb_list = msigdb_of_interest, 
                          named_degs_list = epi_for_gsea, 
                          fgsea_res = fgsea_results_epi)

# 7.4 epi proportion shift ####
prop_calcs = pdo@meta.data %>%
  dplyr::filter(coarse_anno =="Epithelial") %>%
  dplyr::group_by(treatment, finalized_fine_anno) %>%
  dplyr::summarise(count = n())
prop_calcs$finalized_fine_anno = factor(prop_calcs$finalized_fine_anno, levels = epi_fine_anno_order_rainbow)
ggplot(prop_calcs, aes(x=treatment, y=count))+
  geom_bar(aes(fill=finalized_fine_anno), stat = "identity", position = "fill")+theme_pubr(legend="right")+
  theme_pubr(legend = "right")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x="", y = "Proportion", fill = "")+
  scale_fill_manual(values = epi_fine_anno_colors_rainbow)

# chi square
prop_stat = pdo@meta.data %>%
  dplyr::filter(coarse_anno =="Epithelial") %>%
  dplyr::group_by(treatment, finalized_fine_anno) %>%
  dplyr::summarise(count = n()) %>%
  spread(finalized_fine_anno, count, fill = 0)
prop_stat = as.data.frame(prop_stat)
rownames(prop_stat) = prop_stat$treatment
prop_stat$treatment <- NULL
chisq.test(prop_stat) #p=2.2e-16
chisq.posthoc.test(prop_stat, method = "bonferroni")

# 7.5 tumor differential expression - finding gsea of markers ####
pdo_lumt <- subset(x=pdo_unintegrated, subset = finalized_broad_anno=="Luminal_Tumor")
Idents(object = pdo_lumt) <- "finalized_fine_anno"

# setting up marker genes
# lumtum.markers.1 <- FindMarkers(pdo_lumt, ident.1 = "Luminal_Tumor1", test.use="MAST", latent.vars="patient_replicate", recorrect_umi = F)
# lumtum.markers.3 <- FindMarkers(pdo_lumt, ident.1 = "Luminal_Tumor3", test.use="MAST", latent.vars="patient_replicate", recorrect_umi = F)
lumtum1 = read.csv("../results/analyses/lumtum1.markers.csv")
lumtum3 = read.csv("../results/analyses/lumtum3.markers.csv")

markers_list = list(lumtum_1 = lumtum1,
                    lumtum_3 = lumtum3)

combined_volcano_plot(markers_list, ncol = 2)
lumtum_for_gsea = create_gsea_vectors(markers_list)
fgsea_results_lumtum = run_through_fgsea(msigdb_list = msigdb_of_interest, 
                                         named_degs_list = lumtum_for_gsea) # this takes a few minutes

# saving fgsea results because they take a while
saveRDS(fgsea_results_lumtum, file.path(save_path, "fgsea_results_lumtum.rds"))

# 7.6 lumtumor proportion shift ####
prop_calcs = pdo@meta.data %>%
  dplyr::filter(finalized_broad_anno %in% c("Luminal_Tumor")) %>%
  dplyr::group_by(treatment, finalized_fine_anno) %>%
  dplyr::summarise(count = n())
prop_calcs$finalized_fine_anno = factor(prop_calcs$finalized_fine_anno, levels = epi_fine_anno_order_rainbow)
ggplot(prop_calcs, aes(x=treatment, y=count))+
  geom_bar(aes(fill=finalized_fine_anno), stat = "identity", position = "fill", color = "white")+theme_pubr(legend="right")+
  theme_minimal()+
  theme(legend.position = "none")+
  labs(x="", y = "Proportion", fill = "")+
  scale_fill_manual(values = epi_fine_anno_colors_rainbow[5:10])

# chi square
prop_stat = pdo@meta.data %>%
  dplyr::filter(finalized_broad_anno %in% c("Luminal_Tumor")) %>%
  # dplyr::filter(treatment %in% c("NET", "UT")) %>% to check specific comparisons
  dplyr::group_by(treatment, finalized_fine_anno) %>%
  dplyr::summarise(count = n()) %>%
  spread(finalized_fine_anno, count, fill = 0)
prop_stat = as.data.frame(prop_stat)
rownames(prop_stat) = prop_stat$treatment
prop_stat$treatment <- NULL
chisq.test(prop_stat) #p=2.2e-16
chisq.posthoc.test(prop_stat, method = "bonferroni")

# 7.7 stro differential expression ####
Idents(object = pdo_stro) <- "treatment"

# setting up marker genes
# n.c.stro.markers <- FindMarkers(pdo_stro, ident.1 = "NET", ident.2 = "UT", test.use="MAST", latent.vars="patient_replicate")
# t.c.stro.markers <- FindMarkers(pdo_stro, ident.1 = "DET", ident.2 = "UT", test.use="MAST", latent.vars="patient_replicate")
# t.n.stro.markers <- FindMarkers(pdo_stro, ident.1 = "DET", ident.2 = "NET", test.use="MAST", latent.vars="patient_replicate")
n.c.stro = read.csv("../results/analyses/n.c.stro.markers.csv")
t.c.stro = read.csv("../results/analyses/t.c.stro.markers.csv")
t.n.stro = read.csv("../results/analyses/t.n.stro.markers.csv")

markers_list = list(n.c.allsamples = n.c.stro, 
                    t.c.allsamples = t.c.stro, 
                    t.n.allsamples = t.n.stro)

# running setup functions - found in 6.2
combined_volcano_plot(markers_list, ncol = 3)
stro_for_gsea = create_gsea_vectors(markers_list)
fgsea_results_stro = run_through_fgsea(msigdb_list = msigdb_of_interest, 
                                       named_degs_list = stro_for_gsea) # this takes a few minutes

# saving fgsea results because they take a while
saveRDS(fgsea_results_stro, file.path(save_path, "fgsea_results_stro.rds"))

# 7.8 imm differential expression ####
Idents(object = pdo_imm) <- "treatment"
# n.c.imm.markers_fullreps <- FindMarkers(pdo_imm_fullreps, ident.1 = "NET", ident.2 = "UT", test.use="MAST",  latent.vars="patient_replicate")
# t.c.imm.markers_fullreps <- FindMarkers(pdo_imm_fullreps, ident.1 = "DET", ident.2 = "UT", test.use="MAST", latent.vars="patient_replicate")
# t.n.imm.markers_fullreps <- FindMarkers(pdo_imm_fullreps, ident.1 = "DET", ident.2 = "NET", test.use="MAST", latent.vars="patient_replicate")
n.c.imm = read.csv("../results/analyses/n.c.imm.markers.csv")
t.c.imm = read.csv("../results/analyses/t.c.imm.markers.csv")
t.n.imm = read.csv("../results/analyses/t.n.imm.markers.csv")

markers_list = list(n.c.allsamples = n.c.imm, 
                    t.c.allsamples = t.c.imm, 
                    t.n.allsamples = t.n.imm)

# running setup functions - found in 6.2
combined_volcano_plot(markers_list, ncol = 3)
imm_for_gsea = create_gsea_vectors(markers_list)
fgsea_results_imm = run_through_fgsea(msigdb_list = msigdb_of_interest, 
                                      named_degs_list = imm_for_gsea) # this takes a few minutes

# saving fgsea results because they take a while
saveRDS(fgsea_results_imm, file.path(save_path, "HOME/ennisc/pdo_sc/results/2024_05_15/objs/fgsea_results_imm.rds"))

# 7.9 imm proportion shift ####
prop_calcs = pdo@meta.data %>%
  dplyr::filter(coarse_anno =="Immune") %>%
  dplyr::group_by(treatment, v2_fine_anno) %>%
  dplyr::summarise(count = n())
prop_calcs$v2_fine_anno = factor(prop_calcs$v2_fine_anno, levels = imm_fine_anno_order_rainbow)
ggplot(prop_calcs, aes(x=treatment, y=count))+
  geom_bar(aes(fill=v2_fine_anno), stat = "identity", position = "fill")+theme_pubr(legend="right")+
  theme_pubr(legend = "right")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x="", y = "Proportion", fill = "")+
  scale_fill_manual(values = imm_fine_anno_colors_rainbow)

# chi square
prop_stat = pdo@meta.data %>%
  dplyr::filter(coarse_anno =="Immune") %>%
  dplyr::group_by(treatment, v2_fine_anno) %>%
  dplyr::summarise(count = n()) %>%
  spread(fine_anno, count, fill = 0)
prop_stat = as.data.frame(prop_stat)
rownames(prop_stat) = prop_stat$treatment
prop_stat$treatment <- NULL
chisq.test(prop_stat) #p=2.2e-16
chisq.posthoc.test(prop_stat, method = "bonferroni")

