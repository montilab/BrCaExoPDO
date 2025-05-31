# 7.0 load packages and data ====
source("../scripts/utils.R")
save_path <- "path/to/save/dir/" # please update file path to save in your desired location
pdo_unintegrated <- readRDS(file.path(save_path, "pdo_unintegrated.rds"))
pdo <- readRDS(file.path(save_path, "pdo_integrated.rds"))
brca_atlas <- read.csv(file.path(save_path, "cellxgene.csv")) # please download atlas metadata, preprint here: https://www.biorxiv.org/content/10.1101/2025.03.13.643025v2.full
imm_atlas <- read.RDS(file.path(save_path, "cellxgene.rds")) # please download atlas immune obj, preprint here: https://www.biorxiv.org/content/10.1101/2025.03.13.643025v2.full
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

apop_er_stress = read.gmt("path_to/GOBP_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_IN_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS.v2023.2.Hs.gmt")
hallmark_hypoxia = read.gmt("path_to/HALLMARK_HYPOXIA.v2023.2.Hs.gmt")
hallmark_upr = read.gmt("path_to/HALLMARK_UNFOLDED_PROTEIN_RESPONSE.v2023.2.Hs.gmt")
t_cytotoxic = read.gmt("path_to/BIOCARTA_TCYTOTOXIC_PATHWAY.v2024.1.Hs.gmt")

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

# 7.10 proportion comparisons to atlas ====
pdo_summary <- pdo@meta.data %>%
  dplyr::group_by(patientID, coarse_anno) %>% 
  dplyr::summarize(count = n(), .groups = 'drop') %>%
  dplyr::group_by(patientID) %>%
  dplyr::mutate(proportion = count / sum(count)) 

atlas_summary <- brca_atlas %>%
  # create a new column with broad categories
  dplyr::mutate(broad_category = case_when(
    # epithelial: usually the malignant (tumor) cells
    author_cell_type == "Malignant" ~ "Epithelial",
    
    # stromal: includes CAFs, endothelial cells, myofibroblasts, vascular smooth muscle, etc.
    author_cell_type %in% c("CAFs (CA12+)", "CAFs (COL11A1+)", "CAFs (DPT+)", 
                            "CAFs (LAMP5)", "CAFs (PI16+)", "Myofibroblast", 
                            "Prolif Stromal", "Endo Arter.", "Endo Capil.", 
                            "Endo Imm.", "Endo Lymph.", "Endo Ven.", "VSMC", 
                            "PCs (ECM)") ~ "Stromal",
    
    # immune: includes B cells, T cells, dendritic cells, macrophages, NK cells, etc.
    author_cell_type %in% c("B_Cell", "CD4_Naive", "CD4_Tem", "CD4_Tfh", "CD4_Th", 
                            "CD4_Treg", "CD8_Isg", "CD8_Tem", "CD8_Tex", "cDC1", 
                            "cDC2", "Mac_Col.", "Mac_Lipo", "Mac_Prolif.", "Mast", 
                            "mDC", "NK", "pDC", "Plasma_IgG", "T_Prolif.") ~ "Immune",
    
    # default case if none of the above match
    TRUE ~ "Other"
  )) %>%
  # group by donor and the new broad category
  dplyr::group_by(donor_id, broad_category) %>% 
  dplyr::summarise(count = n(), .groups = 'drop') %>%
  dplyr::group_by(donor_id) %>%
  dplyr::mutate(proportion = count / sum(count))

# merge together
pdo_summary$patientID <- paste0("pdo_", pdo_summary$patientID)
colnames(atlas_summary) = c("patientID", "coarse_anno", "count", "proportion")

merge_summary = rbind(pdo_summary, atlas_summary)

# reorder for descending perc epi
patient_order <- merge_summary %>%
  dplyr::filter(coarse_anno == "Epithelial") %>%
  dplyr::arrange(desc(proportion)) %>%
  dplyr::pull(patientID)

merge_summary <- merge_summary %>%
  dplyr::mutate(patientID = factor(patientID, levels = patient_order))

merge_summary = na.omit(merge_summary)

highlight_ids <- c("pdo_318", "pdo_377", "pdo_409") 

# convert patientID to numeric based on the factor order.
highlight_positions <- data.frame(patientID = highlight_ids) %>%
  dplyr::mutate(patientID = factor(patientID, levels = levels(merge_summary$patientID))) %>%
  dplyr::mutate(x = as.numeric(patientID))

# plot
ggplot(merge_summary, aes(x = patientID, y = proportion, fill = coarse_anno)) +
  geom_bar(stat = "identity", width = 0.9) +
  # overlay a rectangle only for highlighted patients
  geom_rect(data = highlight_positions,
            aes(xmin = x - 0.45, xmax = x + 0.45, ymin = 0, ymax = 1),
            fill = NA, color = "black", size = 1, inherit.aes = FALSE) +
  labs(x = "", y = "Proportion", fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "top")+
  scale_fill_brewer(palette = "Set1")
# 7.11 projection of signatures in TCGA ====
tn_df <- read.csv("../results/analyses/t.n.allcells.markers.csv")
tn_sig <- subset(tn_df, p_val_adj < 0.05)
tn_gs <- list(pos = head(tn_sig[order(-tn_sig$avg_log2FC), ]$X, 100),
              neg = head(tn_sig[order(tn_sig$avg_log2FC), ]$X, 100))

gsva_tn <- gsva_data(sigs_list = tn_gs, brca_data = "TCGA", adjust_prolif = TRUE, adjust_inflam = TRUE)
gsva_combined_tcga <- add_combined_signature(gsva_tn)
gsva_combined_tcga <- add_combined_score(gsva_combined_tcga)
gsva_cox_tn <- gsva_cox_fit(gsva_tn, adjust_age = TRUE, adjust_prolif = TRUE, adjust_inflam = TRUE)
tn_cox_corrected_tcga_pam50 <- gsva_cox_pam50(gsva_combined_tcga, brca_data = "TCGA", adjust_age = TRUE, adjust_pam50 = TRUE, five_year = TRUE, adjust_prolif = TRUE, adjust_inflam = TRUE)

gsva_tn_metabric <- gsva_data(sigs_list = tn_gs, brca_data = "METABRIC", adjust_prolif = TRUE, adjust_inflam = TRUE)
gsva_tn_metabric$age_at_index <- gsva_tn_metabric$AGE_AT_DIAGNOSIS
gsva_combined_mb <- add_combined_signature(gsva_tn_metabric)
gsva_combined_mb <- add_combined_score(gsva_combined_mb)
gsva_cox_tn_metabric <- gsva_cox_fit(gsva_tn_metabric, adjust_age = TRUE, adjust_prolif = TRUE, adjust_inflam = TRUE)
tn_cox_corrected_mb_pam50 <- gsva_cox_pam50(gsva_combined_mb, brca_data = "METABRIC", adjust_age = TRUE, adjust_pam50 = TRUE, five_year = TRUE, adjust_prolif = TRUE, adjust_inflam = TRUE)

# plot theming constants
time_breaks <- seq(0, 1825, by = 365)
time_labels <- c("0", "1 year", "2 year", "3 year", "4 year", "5 year")
common_theme <- theme_minimal() + 
  theme(axis.ticks = element_blank(),
        legend.position = "bottom")

color_labels <- c("Low: NDexo-like", "High: T2Dexo-like")
color_values <- rev(treatment_colors[1:2])

# TCGA Cox model and plot
cox_model_tcga <- coxph(Surv(as.numeric(time), vital_status_5) ~ age_at_index + strata(combined_group) + prolif + inflam + subtype_m_rna, 
                        data = Biobase::pData(gsva_combined_tcga))

cox_fit_tcga <- summary(tn_cox_corrected_tcga_pam50$combined)
HR_tcga <- round(cox_fit_tcga$conf.int[, "exp(coef)"], 2)
lower_CI_tcga <- round(cox_fit_tcga$conf.int[, "lower .95"], 2)
upper_CI_tcga <- round(cox_fit_tcga$conf.int[, "upper .95"], 2)
pval_tcga <- signif(cox_fit_tcga$coefficients["combined", 5], 3)

plot_tcga <- ggadjustedcurves(fit = cox_model_tcga, 
                              data = Biobase::pData(gsva_combined_tcga),
                              method = "conditional",
                              variable = "combined_group",
                              xlab = "Days", ylab = "Overall Survival Probability",
                              legend.title = "Combined Score", 
                              ggtheme = theme_minimal())

plot_tcga <- plot_tcga +
  annotate("text", x = 50, y = 0.7, label = format_annotation(HR_tcga, lower_CI_tcga, upper_CI_tcga, pval_tcga),
           size = 4, hjust = 0) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0.6, 1)) +
  scale_x_continuous(breaks = time_breaks, labels = time_labels) +
  coord_cartesian(xlim = c(0, 1825)) +
  labs(x = NULL) +
  common_theme +
  scale_color_manual(values = color_values, labels = color_labels)

# METABRIC Cox model and plot
cox_model_mb <- coxph(Surv(as.numeric(time), vital_status_5) ~ AGE_AT_DIAGNOSIS + strata(combined_group) + prolif + inflam + Pam50_SUBTYPE, 
                      data = Biobase::pData(gsva_combined_mb))

cox_fit_mb <- summary(tn_cox_corrected_mb_pam50$combined)
HR_mb <- round(cox_fit_mb$conf.int[, "exp(coef)"], 2)
lower_CI_mb <- round(cox_fit_mb$conf.int[, "lower .95"], 2)
upper_CI_mb <- round(cox_fit_mb$conf.int[, "upper .95"], 2)
pval_mb <- signif(cox_fit_mb$coefficients["combined", 5], 2)

plot_mb <- ggadjustedcurves(fit = cox_model_mb, 
                            data = Biobase::pData(gsva_combined_mb),
                            method = "conditional",
                            variable = "combined_group",
                            xlab = "Days", ylab = "Overall Survival Probability",
                            legend.title = "Combined Score", 
                            ggtheme = theme_minimal())

plot_mb <- plot_mb +
  annotate("text", x = 50, y = 0.7, label = format_annotation(HR_mb, lower_CI_mb, upper_CI_mb, pval_mb),
           size = 4, hjust = 0) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0.6, 1)) +
  scale_x_continuous(breaks = time_breaks, labels = time_labels) +
  coord_cartesian(xlim = c(0, 1825)) +
  labs(x = NULL) +
  common_theme +
  scale_color_manual(values = color_values, labels = color_labels)

# TCGA  CI plot
surv_tcga_fit <- survival::survfit(cox_model_tcga, data = Biobase::pData(gsva_combined_tcga))

ci_plot_tcga <- ggsurvplot(surv_tcga_fit, data = Biobase::pData(gsva_combined_tcga), conf.int = TRUE,
                           palette = color_values, ggtheme = theme_minimal(),
                           xlab = "", ylab = "Overall Survival Probability",
                           legend.labs = color_labels, legend.title = "Combined Score")

ci_plot_tcga$plot <- ci_plot_tcga$plot +
  annotate("text", x = 0.5 * 50, y = 0.5,
           label = format_annotation(HR_tcga, lower_CI_tcga, upper_CI_tcga, pval_tcga),
           size = 4, hjust = 0) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0.45, 1)) +
  scale_x_continuous(breaks = time_breaks, labels = time_labels) +
  coord_cartesian(xlim = c(0, 1825)) +
  common_theme

# METABRIC  CI plot
surv_mb_fit <- survival::survfit(cox_model_mb, data = Biobase::pData(gsva_combined_mb))

ci_plot_mb <- ggsurvplot(surv_mb_fit, data = Biobase::pData(gsva_combined_mb), conf.int = TRUE,
                         palette = color_values, ggtheme = theme_minimal(),
                         xlab = "", ylab = "Overall Survival Probability",
                         legend.labs = color_labels, legend.title = "Combined Score")

ci_plot_mb$plot <- ci_plot_mb$plot +
  annotate("text", x = 0.5 * 50, y = 0.5,
           label = format_annotation(HR_mb, lower_CI_mb, upper_CI_mb, pval_mb),
           size = 4, hjust = 0) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0.45, 1)) +
  scale_x_continuous(breaks = time_breaks, labels = time_labels) +
  coord_cartesian(xlim = c(0, 1825)) +
  common_theme

# final combined plots
plot_tcga_clean <- plot_tcga +
  labs(x = NULL, y = "Overall Survival Probability", title = "TCGA") +
  theme(plot.title = element_text(hjust = 0.5))

plot_mb_clean <- plot_mb +
  labs(x = NULL, y = NULL, title = "METABRIC") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5))

ci_plot_tcga_clean <- ci_plot_tcga$plot +
  labs(x = NULL, y = "Overall Survival Probability", title = "TCGA") +
  theme(plot.title = element_text(hjust = 0.5))

ci_plot_mb_clean <- ci_plot_mb$plot +
  labs(x = NULL, y = NULL, title = "METABRIC") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5))

# Combine: title above each plot + plots side-by-side
final_simple_plot <- (plot_tcga_clean | plot_mb_clean) +
  plot_layout(guides = "collect", heights = c(0.08, 1)) &
  theme(legend.position = "bottom")

final_ci_plot <- (ci_plot_tcga_clean | ci_plot_mb_clean) +
  plot_layout(guides = "collect", heights = c(0.08, 1)) &
  theme(legend.position = "bottom")

final_simple_plot
final_ci_plot
# 7.12 projection of chopT signature in atlas ====
t_cell_ids = rownames(subset(imm_atlas@meta.data, cluster_annot %in% c("CD4_Naive", "CD4_Tem", "CD4_Tfh", "CD4_Th", "CD4_Treg", "CD8_Isg", "CD8_Tem", "CD8_Tex", "T_Prolif.")))
t_cells <- subset(imm_atlas, cells = t_cell_ids)

t.n.imm = read.csv("../results/analyses/t.n.imm.markers.csv")
chopt_genes <- subset(t.n.imm, cluster == "ChopT")$gene

my_genesets <- list(
  "chopt" = chopt_genes, 
  "t_cytotoxic" = t_cytotoxic$BIOCARTA_TCYTOTOXIC_PATHWAY,
  "apop_er_stress" = apop_er_stress$GOBP_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_IN_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS,
  "hypoxia" = hallmark_hypoxia$HALLMARK_HYPOXIA,
  "upr" = hallmark_upr$HALLMARK_UNFOLDED_PROTEIN_RESPONSE,
)

t_cells <- AddModuleScore(object = t_cells, features = my_genesets, name = c(names(my_genesets)))

cor_results_seurat <- get_cor_tests(c("hypoxia", "upr", "t_cytotoxic", "t_effector"), 
                                    c("chopt"), t_cells@meta.data)

ggplot(t_cells@meta.data, aes(y = t_effector, x = chopt)) +
  geom_point(alpha = 0.2, color = "grey80") +
  geom_smooth(color = "red", method = "lm", se = TRUE) +
  ylab("T Effector Score") + xlab("ChopT Score") +
  stat_cor(method = "spearman", label.x.npc = "center", label.y.npc = "top") +
  theme_minimal()
ggplot(t_cells@meta.data, aes(y = t_cytotoxic, x = chopt)) +
  geom_point(alpha = 0.2, color = "grey80") +
  geom_smooth(color = "red", method = "lm", se = TRUE) +
  ylab("T Cytotoxic Score") + xlab("ChopT Score") +
  stat_cor(method = "spearman", label.x.npc = "center", label.y.npc = "top") +
  theme_minimal()
ggplot(t_cells@meta.data, aes(y = upr, x = chopt)) +
  geom_point(alpha = 0.2, color = "grey80") +
  geom_smooth(color = "red", method = "lm", se = TRUE) +
  ylab("Unfolded Protein Response Score") + xlab("ChopT Score") +
  stat_cor(method = "spearman", label.x.npc = "center", label.y.npc = "top") +
  theme_minimal()
ggplot(t_cells@meta.data, aes(y = hypoxia, x = chopt)) +
  geom_point(alpha = 0.2, color = "grey80") +
  geom_smooth(color = "red", method = "lm", se = TRUE) +
  ylab("Hypoxia Score") + xlab("ChopT Score") +
  stat_cor(method = "spearman", label.x.npc = "center", label.y.npc = "top") +
  theme_minimal()
