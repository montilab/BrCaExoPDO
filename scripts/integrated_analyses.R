# 8.0 load packages and data ====
source("../scripts/utils.R")
save_path <- "path/to/save/dir/" # please update file path to save in your desired location
pdo <- readRDS(file.path(save_path, "pdo_integrated.rds"))
pdo_epi <- readRDS(file.path(save_path, "pdo_epi_integrated.rds"))

pdo@meta.data$finalized_fine_anno_short = pdo@meta.data$finalized_fine_anno
pdo@meta.data$finalized_fine_anno_short[pdo@meta.data$finalized_fine_anno =="Luminal_Tumor1"] <- "LT1"
pdo@meta.data$finalized_fine_anno_short[pdo@meta.data$finalized_fine_anno =="Luminal_Tumor2"] <- "LT2"
pdo@meta.data$finalized_fine_anno_short[pdo@meta.data$finalized_fine_anno =="Luminal_Tumor3"] <- "LT3"
pdo@meta.data$finalized_fine_anno_short[pdo@meta.data$finalized_fine_anno =="Luminal_Tumor4"] <- "LT4"
pdo@meta.data$finalized_fine_anno_short[pdo@meta.data$finalized_fine_anno =="Luminal_Tumor5"] <- "LT5"
pdo@meta.data$finalized_fine_anno_short[pdo@meta.data$finalized_fine_anno =="Basal_Tumor1"] <- "BT1"
pdo@meta.data$finalized_fine_anno_short[pdo@meta.data$finalized_fine_anno =="Basal_Tumor2"] <- "BT2"
pdo@meta.data$finalized_fine_anno_short[pdo@meta.data$finalized_fine_anno =="Cycling_Epithelial"] <- "CYC"
pdo@meta.data$finalized_fine_anno_short[pdo@meta.data$finalized_fine_anno =="Luminal_Progenitor"] <- "LP"

pdo_circ <- read.csv("../results/analyses/organoid_circularity.csv")

# 8.1 ligand receptor analysis  ####
#vignette for basic: https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html
#vignette for comparing conditions: https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets.html

#set cellchat database
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
CellChat::showDatabaseCategory(CellChatDB)

# use a subset of CellChatDB for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB - possible switch to liana::select_resource("Consensus"))

#make cellchat object 
#need to make separate for each condition then can merge together
pdo_net = subset(x=pdo, subset = treatment =="NET")
pdo_det = subset(x=pdo, subset = treatment =="DET")

#preprocess net
dataNET = pdo_net@assays$SCT$data # normalized data matrix
metaNET = pdo_net[[]] # a dataframe with rownames containing cell mata data
chatNET <- createCellChat(object = dataNET, meta = metaNET, group.by = "finalized_broad_anno")
chatNET@DB <- CellChatDB.use # set the used database in the object
chatNET <- subsetData(chatNET) # This step is necessary even if using the whole database
chatNET <- identifyOverExpressedGenes(chatNET)
chatNET <- identifyOverExpressedInteractions(chatNET) 
#compute communication probability and infer cellular communication network
chatNET <- computeCommunProb(chatNET, type = "triMean") #this takes a few minutes
chatNET <- computeCommunProbPathway(chatNET) #infer the cell-cell communication at a signaling pathway level
chatNET <- aggregateNet(chatNET) #calculate the aggregated cell-cell communication network
#centrality and communication patterns
chatNET <- netAnalysis_computeCentrality(chatNET, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
selectK(chatNET, pattern = "outgoing") #identify and visualize outgoing communication pattern of secreting cells
chatNET <- identifyCommunicationPatterns(chatNET, pattern = "outgoing", k = 7)
selectK(chatNET, pattern = "incoming")
chatNET <- identifyCommunicationPatterns(chatNET, pattern = "incoming", k = 7)
chatNET <- computeNetSimilarity(chatNET, type = "functional")
chatNET <- computeNetSimilarity(chatNET, type = "structural")

#preprocess net
dataDET = pdo_det@assays$SCT$data # normalized data matrix
metaDET = pdo_det[[]] # a dataframe with rownames containing cell mata data
chatDET <- createCellChat(object = dataDET, meta = metaDET, group.by = "finalized_broad_anno")
chatDET@DB <- CellChatDB.use # set the used database in the object
chatDET <- subsetData(chatDET) # This step is necessary even if using the whole database
chatDET <- identifyOverExpressedGenes(chatDET)
chatDET <- identifyOverExpressedInteractions(chatDET) 
#compute communication probability and infer cellular communication network
chatDET <- computeCommunProb(chatDET, type = "triMean") #this takes a few minutes
chatDET <- computeCommunProbPathway(chatDET) #infer the cell-cell communication at a signaling pathway level
chatDET <- aggregateNet(chatDET) #calculate the aggregated cell-cell communication network
#centrality and communication patterns
chatDET <- netAnalysis_computeCentrality(chatDET, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
selectK(chatDET, pattern = "outgoing") #identify and visualize outgoing communication pattern of secreting cells
chatDET <- identifyCommunicationPatterns(chatDET, pattern = "outgoing", k = 5)
selectK(chatDET, pattern = "incoming")
chatDET <- identifyCommunicationPatterns(chatDET, pattern = "incoming", k = 5)
chatDET <- computeNetSimilarity(chatDET, type = "functional")
chatDET <- computeNetSimilarity(chatDET, type = "structural")

# save individual chat objects
saveRDS(chatNET, file = file.path(save_path, "cellchat_nd.rds"))
saveRDS(chatDET, file = file.path(save_path, "cellchat_t2d.rds"))

#merge - only doing ND vs T2D for comparison plots
object.list <- list(NET = chatNET, DET = chatDET)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

#identify altered interactions and cell populations
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), color.use = rev(treatment_colors[1:2]))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight", color.use = rev(treatment_colors[1:2]))
gg1 + gg2 #compare total number of interactions and interaction strength

#circle plot showing differential number of interactions/interaction strength among different cell pops
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, color.use = broadcelltype_colors_alpha, arrow.size = 0.01)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", color.use = broadcelltype_colors_alpha, arrow.size = 0.01)
#heatmap
gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

#major sources and targets in 2d space
#identify cell populations with significant changes in sending/receiving signals
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)

#signaling changes of specific cell populations
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Basal_Tumor")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "T", signaling.exclude = c("COLLAGEN", "LAMININ"))
patchwork::wrap_plots(plots = list(gg1,gg2))

#Identify altered signaling with distinct network architecture and interaction strength
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
rankNet(cellchat, mode = "comparison", measure = "weight", 
        sources.use = NULL, targets.use = NULL, 
        stacked = T, do.stat = TRUE, color.use = c("orange", "dodgerblue2"), 
        cutoff.pvalue = 0.01, font.size = 6)

# figure out counts of pathways and genes
table(cellchat@netP$DET$pathways %in% cellchat@netP$NET$pathways)
table(cellchat@netP$NET$pathways %in% cellchat@netP$DET$pathways)
length(unique(c(cellchat@netP$NET$pathways,cellchat@netP$DET$pathways)))
length(unique(c(cellchat@net$NET$LRs,cellchat@net$DET$LRs)))

#compare outgoing (or incoming) signaling patterns associated with each cell population
i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 10, height = 15)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 10, height = 15)
ComplexHeatmap::draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 10, height = 15, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 10, height = 15, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 10, height = 15, color.heatmap = "OrRd", font.size = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 10, height = 15, color.heatmap = "OrRd", font.size = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

#Identify dysfunctional signaling by comparing the communication probabities
netVisual_bubble(cellchat, sources.use = c(7:9),  targets.use = c(12), 
                 comparison = c(1, 2), angle.x = 45, 
                 color.text = c("#E41A1C", "#377EB8"), font.size = 6) #sources.use is to determine the sending cell type
netVisual_bubble(cellchat, sources.use = c(12),  targets.use = c(7:9), comparison = c(1, 2), angle.x = 45, color.text = c("#E41A1C", "#377EB8"), font.size = 4) #sources.use is to determine the sending cell type

#Identify dysfunctional signaling by comparing the communication probabities
netVisual_bubble(cellchat, sources.use = 12,  targets.use = c(1,2,3,4,5,6,7,8,9,10,11,12),
                 comparison = c(1, 2), angle.x = 45, thresh = 0.05,
                 color.text = c("orange", "dodgerblue2"), font.size = 6, remove.isolate = T) #sources.use is to determine the sending cell type

#identify dysfunctional signaling
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "DET"
# define a char name used for storing the results of differential expression analysis
features.name = paste0(pos.dataset, ".merged")

# perform differential expression analysis 
# Of note, compared to CellChat version < v2, CellChat v2 now performs an ultra-fast Wilcoxon test using the presto package, which gives smaller values of logFC. Thus we here set a smaller value of thresh.fc compared to the original one (thresh.fc = 0.1). Users can also provide a vector and dataframe of customized DEGs by modifying the cellchat@var.features$LS.merged and cellchat@var.features$LS.merged.info. 

cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,thresh.p = 0.05, group.DE.combined = FALSE) 
#> Use the joint cell labels from the merged CellChat object

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "DET",ligand.logFC = 0.05, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated receptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "NET",ligand.logFC = -0.05, receptor.logFC = NULL)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = c(7:9),  targets.use = c(12), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = 7,  targets.use = c(2, 6:11), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> You may try the function `netVisual_chord_cell` for visualizing individual signaling pathway

#Visually compare cell-cell communication using Hierarchy plot, Circle plot or Chord diagram
pathways.show <- c("AREG") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}

#> Do heatmap based on a single object
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

# 8.2 pseudotime analysis ####
#convert into monocle object using seurat wrappers and run basic preprocessing
pdo_epi@reductions$umap <- pdo_epi@reductions$umap.scvi # monocle requires this format specifically
pdo_lum <- subset(x = pdo_epi, subset = finalized_fine_anno %in% c("LT1", "LT2", "LT3", "LT4", "LT5", "LP", "LEC"))
mon_pdo <- SeuratWrappers::as.cell_data_set(pdo_lum)
mon_pdo <- cluster_cells(mon_pdo)
p1 <- plot_cells(mon_pdo, show_trajectory_graph = FALSE)
p2 <- plot_cells(mon_pdo, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)

#generate trajectory path
mon_pdo <- learn_graph(mon_pdo)
plot_cells(mon_pdo, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)

plot_cells(mon_pdo,
           color_cells_by = "finalized_fine_anno",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)+
  scale_color_manual(values = c(epi_fine_anno_short_colors_alpha[4:10]))

normal_epithelial_cells = c(subset(pdo_lum@meta.data, finalized_fine_anno %in% c("LEC"))$barcode)

mon_pdo <- order_cells(mon_pdo, root_cells = normal_epithelial_cells)
plot_cells(mon_pdo, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE)

#can i determine if one cell type is more on the beginning/end of pseudotime?
pseudotime_labels = as.data.frame(mon_pdo@principal_graph_aux@listData$UMAP$pseudotime)
colnames(pseudotime_labels) = c("monocle3_pseudotime")
pseudotime_labels$barcode = rownames(pseudotime_labels)

meta_data <- dplyr::left_join(pdo_lum@meta.data, pseudotime_labels, by="barcode")
rownames(meta_data) = meta_data$barcode
pdo_lum@meta.data <- meta_data

ggplot(pdo_lum@meta.data, aes(y=reorder(finalized_fine_anno, -monocle3_pseudotime), x=monocle3_pseudotime))+
  geom_violin(aes(fill=finalized_fine_anno), scale = "width", bw = 0.2)+
  geom_boxplot(fill = "white", width = 0.1, outlier.shape = NA)+
  theme_minimal()+
  theme(legend.position = "none")+
  scale_fill_manual(values = c(epi_fine_anno_short_colors_alpha[4:10]))+
  xlab("Pseudotime")+ylab("")

ggplot(pdo_lum@meta.data, aes(x=treatment, y=monocle3_pseudotime))+
  geom_violin(aes(fill=treatment), bw = 0.2)+
  geom_boxplot(fill = "white", width = 0.1, outlier.shape = NA)+
  theme_minimal()+
  theme(legend.position = "none")+
  xlab("")+ylab("Pseudotime")+
  scale_fill_manual(values = treatment_colors)+
  stat_compare_means(comparisons = list(c("UT", "NET"), c("UT", "DET"), c("NET", "DET")), # switch for rank thing monti suggested
                     method="wilcox.test",label = "p.signif", method.args = list(alternative = "less"), hide.ns = T)

saveRDS(mon_pdo, file = file.path(save_path, "mon_pdo_lum.rds"))

#differential expression - what genes define pseudotime?
graph_test_results <- graph_test(mon_pdo, neighbor_graph = "principal_graph", cores = 4) #this takes some time, can make faster with cores = 4 (for example)
head(graph_test_results, error=FALSE, message=FALSE, warning=FALSE)

# one or the other - depending on if read in or calculated
deg_ids <- subset(graph_test_results[order(graph_test_results$morans_I, decreasing = TRUE),], q_value < 0.05)$X
deg_ids <- rownames(subset(graph_test_results[order(graph_test_results$morans_I, decreasing = TRUE),], q_value < 0.05))

monocle_degs = pdo_lum@assays$SCT$data[rownames(pdo_lum@assays$SCT$data) %in% deg_ids,]
monocle_degs =as.data.frame(t(as.matrix(monocle_degs)))
monocle_degs$barcode = rownames(monocle_degs)
monocle_degs = merge(monocle_degs, dplyr::select(pdo_lum[[]], monocle3_pseudotime, barcode, finalized_fine_anno), by="barcode")

#modules
# gene_modules <- find_gene_modules_corrected(mon_pdo[deg_ids,],
#                                             resolution=10^(-2))
# table(gene_modules$module) # 8 gene modules

gene_modules <- read.csv("../results/analyses/monocle_gene_modules.csv")
rownames(gene_modules) <- gene_modules$X
gene_modules$X <- NULL

# adding gene module score to main pdo to see if treatment is different
pdo <- Seurat::AddModuleScore(pdo, features = list(mod1 = subset(gene_modules, module == "1")$id,
                                                   mod2 = subset(gene_modules, module == "2")$id,
                                                   mod3 = subset(gene_modules, module == "3")$id,
                                                   mod4 = subset(gene_modules, module == "4")$id,
                                                   mod5 = subset(gene_modules, module == "5")$id,
                                                   mod6 = subset(gene_modules, module == "6")$id,
                                                   mod7 = subset(gene_modules, module == "7")$id,
                                                   mod8 = subset(gene_modules, module == "8")$id))

pdo@meta.data <- pdo@meta.data %>% dplyr::rename("Module1" = `Cluster1`,
                                                 "Module2" = `Cluster2`,
                                                 "Module3" = `Cluster3`,
                                                 "Module4" = `Cluster4`,
                                                 "Module5" = `Cluster5`,
                                                 "Module6" = `Cluster6`,
                                                 "Module7" = `Cluster7`,
                                                 "Module8" = `Cluster8`) 

ggplot(subset(pdo@meta.data, finalized_fine_anno_short %in% c("LT1", "LT2", "LT3", "LT4", "LT5", "LP", "LEC")), 
       aes(x=treatment, y=Module5))+
  geom_violin(aes(fill=treatment), bw = 0.2)+
  geom_boxplot(fill = "white", width = 0.1, outlier.shape = NA)+
  theme_minimal()+
  theme(legend.position = "none")+
  xlab("")+ylab("Module 5 Score")+
  scale_fill_manual(values = treatment_colors)

# loop to figure out which are significant
for (i in 1:8) {
  # construct the module name
  module_name <- paste0("Module", i)
  
  # fit the linear mixed-effects model
  lm_model <- lmerTest::lmer(data = pdo@meta.data, 
                             formula = as.formula(paste(module_name, "~ treatment + correlation_score + correlation_score*patientID + (1|patientID)")))
  
  # perform ANOVA
  anova_result <- anova(lm_model)
  print(paste("ANOVA for", module_name))
  print(anova_result)
  
  # print the summary
  summary_result <- summary(lm_model)
  print(paste("Summary for", module_name))
  print(summary_result)
  
  # perform emmeans analysis
  emmeans_result <- emmeans(lm_model, pairwise ~ treatment | correlation_score)
  print(paste("emmeans for", module_name))
  print(emmeans_result)
} # all cells
# mod1: ANOVA(< 2e-16), emmeans(DETvNET=<.0001, DETvUT=0.0009, NETvUT=<.0001)
# mod2: ANOVA(< 2e-16), emmeans(DETvNET=<.0001, DETvUT=<.0001, NETvUT=<.0001)
# mod3: ANOVA(< 2e-16), emmeans(DETvNET=<.0001, DETvUT=<.0001, NETvUT=<.0001)
# mod4: ANOVA(< 2e-16), emmeans(DETvNET=<.0001, DETvUT=<.0001, NETvUT=<.0001)
# mod5: ANOVA(< 2e-16), emmeans(DETvNET=<.0001, DETvUT=<.0001, NETvUT=<.0001)
# mod6: ANOVA(< 2e-16), emmeans(DETvNET=<.0001, DETvUT=0.0008, NETvUT=<.0001)
# mod7: ANOVA(< 2e-16), emmeans(DETvNET=<.0001, DETvUT=<.0001, NETvUT=<.0001)
# mod8: ANOVA(< 2e-16), emmeans(DETvNET=<.0001, DETvUT=0.7308, NETvUT=<.0001)

for (i in 1:8) {
  # construct the module name
  module_name <- paste0("Module", i)
  
  # fit the linear mixed-effects model
  lm_model <- lmerTest::lmer(data = subset(pdo@meta.data, finalized_fine_anno_short %in% c("LT1", "LT2", "LT3", "LT4", "LT5", "LP", "LEC")), 
                             formula = as.formula(paste(module_name, "~ treatment + correlation_score + correlation_score*patientID + (1|patientID)")))
  
  # perform ANOVA
  anova_result <- anova(lm_model)
  print(paste("ANOVA for", module_name))
  print(anova_result)
  
  # print the summary
  summary_result <- summary(lm_model)
  print(paste("Summary for", module_name))
  print(summary_result)
  
  # perform emmeans analysis
  emmeans_result <- emmeans(lm_model, pairwise ~ treatment | correlation_score)
  print(paste("emmeans for", module_name))
  print(emmeans_result)
} # only luminal like epi
# mod1: ANOVA(< 2e-16), emmeans(DETvNET=<.0001, DETvUT=<.0001, NETvUT=<.0001)
# mod2: ANOVA(< 2e-16), emmeans(DETvNET=0.9541, DETvUT=<.0001, NETvUT=<.0001)
# mod3: ANOVA(< 2e-16), emmeans(DETvNET=<.0001, DETvUT=<.0001, NETvUT=<.0001)
# mod4: ANOVA(< 2e-16), emmeans(DETvNET=0.4763, DETvUT=<.0001, NETvUT=<.0001)
# mod5: ANOVA(< 2e-16), emmeans(DETvNET=<.0001, DETvUT=<.0001, NETvUT=<.0001) ** UP in t2d
# mod6: ANOVA(< 2e-16), emmeans(DETvNET=<.0001, DETvUT=0.0135, NETvUT=<.0001)
# mod7: ANOVA(< 2e-16), emmeans(DETvNET=<.0001, DETvUT=0.6920, NETvUT=<.0001)
# mod8: ANOVA(< 2e-16), emmeans(DETvNET=<.0001, DETvUT=<.0001, NETvUT=<.0001)

# can i do gsea on the modules to figure out what they are?
# ranked by morans i from graph_test_results
mod_1_moran <- subset(graph_test_results, X %in% subset(gene_modules, module == "1")$id)
mod_2_moran <- subset(graph_test_results, X %in% subset(gene_modules, module == "2")$id)
mod_3_moran <- subset(graph_test_results, X %in% subset(gene_modules, module == "3")$id)
mod_4_moran <- subset(graph_test_results, X %in% subset(gene_modules, module == "4")$id)
mod_5_moran <- subset(graph_test_results, X %in% subset(gene_modules, module == "5")$id)
mod_6_moran <- subset(graph_test_results, X %in% subset(gene_modules, module == "6")$id)
mod_7_moran <- subset(graph_test_results, X %in% subset(gene_modules, module == "7")$id)
mod_8_moran <- subset(graph_test_results, X %in% subset(gene_modules, module == "8")$id)

gsea_mod_1 <- setNames(mod_1_moran$morans_test_statistic, mod_1_moran$X)
gsea_mod_1 <- sort(gsea_mod_1, decreasing = TRUE)
gsea_mod_2 <- setNames(mod_2_moran$morans_test_statistic, mod_2_moran$X)
gsea_mod_2 <- sort(gsea_mod_2, decreasing = TRUE)
gsea_mod_3 <- setNames(mod_3_moran$morans_test_statistic, mod_3_moran$X)
gsea_mod_3 <- sort(gsea_mod_3, decreasing = TRUE)
gsea_mod_4 <- setNames(mod_4_moran$morans_test_statistic, mod_4_moran$X)
gsea_mod_4 <- sort(gsea_mod_4, decreasing = TRUE)
gsea_mod_5 <- setNames(mod_5_moran$morans_test_statistic, mod_5_moran$X)
gsea_mod_5 <- sort(gsea_mod_5, decreasing = TRUE)
gsea_mod_6 <- setNames(mod_6_moran$morans_test_statistic, mod_6_moran$X)
gsea_mod_6 <- sort(gsea_mod_6, decreasing = TRUE)
gsea_mod_7 <- setNames(mod_7_moran$morans_test_statistic, mod_7_moran$X)
gsea_mod_7 <- sort(gsea_mod_7, decreasing = TRUE)
gsea_mod_8 <- setNames(mod_8_moran$morans_test_statistic, mod_8_moran$X)
gsea_mod_8 <- sort(gsea_mod_8, decreasing = TRUE)

modules_for_gsea = list(module1 = gsea_mod_1, module2 = gsea_mod_2, module3 = gsea_mod_3,
                        module4 = gsea_mod_4, module5 = gsea_mod_5, module6 = gsea_mod_6,
                        module7 = gsea_mod_7, module8 = gsea_mod_8)


fgsea_results_modules = run_through_fgsea_modules(msigdb_list = msigdb_of_interest,
                                                  named_degs_list = modules_for_gsea) 

genes_in_modules <- list(
  module1 = names(gsea_mod_1), module2 = names(gsea_mod_2),
  module3 = names(gsea_mod_3), module4 = names(gsea_mod_4),
  module5 = names(gsea_mod_5), module6 = names(gsea_mod_6),
  module7 = names(gsea_mod_7), module8 = names(gsea_mod_8)
)

max_length <- max(sapply(genes_in_modules, length))
genes_in_modules <- lapply(genes_in_modules, function(x) c(x, rep(NA, max_length - length(x))))
genes_in_modules <- as.data.frame(genes_in_modules) # supplied this to david to define function

# module 1 involved in cell cycle
# module 2 involved in apoptosis, protein folding
# module 3 involved in fatty acid/lipid metabolism, mitochondrial respiration
# module 4 involved in splicing, DNA damage, stress response, t cell mediated toxicity
# module 5 involved in splicing, mitochondrial respiration, EMT
# module 6 involved in translation regulation
# module 7 involved in protein biosynthesis
# module 8 involved in cholesterol homeostasis/metabolism

# adding pseudotime to main object
for (pseudo in pdo_lum@meta.data$monocle3_pseudotime) {
  pdo@meta.data$monocle3_pseudotime[pdo@meta.data$barcode %in% subset(pdo_lum@meta.data, monocle3_pseudotime == pseudo)$barcode] <- pseudo
}

ggplot(subset(pdo@meta.data, finalized_fine_anno_short %in% c("LT1", "LT2", "LT3", "LT4", "LT5", "LP", "LEC")), 
       aes(x=monocle3_pseudotime, y=Module5))+
  geom_point(aes(color = treatment), alpha = 0.1)+
  geom_smooth(aes(color = treatment), method = "lm", linewidth = 2)+
  theme_minimal()+
  theme(legend.position = "none")+
  scale_color_manual(values = treatment_colors)

# module dynamics through pseudotime
module_dynamics = dplyr::select(subset(pdo@meta.data, finalized_fine_anno_short %in% c("LT1", "LT2", "LT3", "LT4", "LT5", "LP", "LEC")),
                                barcode, treatment, finalized_fine_anno_short, monocle3_pseudotime, Module1, Module2, Module3, Module4, Module5, Module6, Module7, Module8)
module_dynamics = reshape2::melt(module_dynamics,id.vars = c("barcode", "treatment", "finalized_fine_anno_short", "monocle3_pseudotime"))

ggplot(module_dynamics, aes(x=monocle3_pseudotime, y=value))+
  geom_point(aes(color = treatment))+
  geom_smooth(color = "black", method = "lm")+
  facet_wrap(~variable, scales = "free_y", ncol = 4)+
  theme_minimal()+
  theme(legend.position = "none")+
  scale_color_manual(values = treatment_colors)

ggplot(module_dynamics, aes(x=monocle3_pseudotime, y=value))+
  geom_point(aes(color = finalized_fine_anno_short))+
  geom_smooth(color = "black", method = "lm")+
  facet_wrap(~variable, scales = "free_y", ncol = 4)+
  theme_minimal()+
  theme(legend.position = "none")+
  scale_color_manual(values = epi_fine_anno_colors_alpha[4:10])+
  xlab("Pseudotime")+ylab("Module Score")

# 7.6 circularity analysis ####
pdo_circ$corrected_treatment[pdo_circ$treatment =="C"] <- "UT"
pdo_circ$corrected_treatment[pdo_circ$treatment =="N"] <- "NET"
pdo_circ$corrected_treatment[pdo_circ$treatment =="T"] <- "DET"

ggplot(pdo_circ, aes(x=corrected_treatment, y=circularity))+
  geom_boxplot()+
  theme_pubr()+
  stat_compare_means(label.y = (-0.03), label.x=1.5)+
  stat_compare_means(comparisons = list(c("C", "N"), 
                                        c("C", "T"), 
                                        c("N", "T")))+
  xlab("Treatment")+
  ylab("Circularity Score")

ggplot(pdo_circ, aes(x=corrected_treatment, y=circularity))+
  geom_violin(aes(fill = corrected_treatment))+
  geom_boxplot(fill = "white", width = 0.1, outlier.shape = NA)+
  theme_minimal()+
  stat_compare_means(comparisons = list(c("UT", "NET"), 
                                        c("UT", "DET"), 
                                        c("NET", "DET")),
                     method="wilcox.test",
                     label = "p.signif",
                     symnum.args = list(cutpoints = c(0, 0.005, 0.01)),
                     hide.ns = T)+
  theme(legend.position = "none")+
  xlab("")+ ylab("Circularity Score")+
  scale_fill_manual(values = treatment_colors)
