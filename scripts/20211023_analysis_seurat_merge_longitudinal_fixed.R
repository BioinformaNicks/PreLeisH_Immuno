#################################################################
##                    Loading in the packages                   #
#################################################################

library(here)
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(clusterProfiler))
suppressMessages(library(enrichplot))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(tidyr))
suppressMessages(library(scRepertoire)) #dev branch
suppressMessages(library(Seurat))
suppressMessages(library(reticulate))
suppressMessages(library(SingleR))
suppressMessages(library(scGate))
#suppressMessages(library(harmony))
suppressMessages(library(pheatmap))
suppressMessages(library(ggsci))
suppressMessages(library(ggpubr))
suppressMessages(library(gridExtra))
suppressMessages(library(patchwork))
suppressMessages(library(geomtextpath))
suppressMessages(library(gplots))
suppressMessages(library(ggrepel))
suppressMessages(library(ggraph))
suppressMessages(library(tidydr))
suppressMessages(library(hrbrthemes))
suppressMessages(library(extrafont))
suppressMessages(library(Cairo))

options(ggrepel.max.overlaps = Inf)
# hrbrthemes::import_roboto_condensed()
# extrafont::loadfonts(device="win")

##################################################################
##                    Setting global variables                   #
##################################################################

text_font <- 'Roboto Condensed'

colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", 
                                        "#0348A6", "#7AC5FF", "#C6FDEC"))

batch1_dir_path <- 'E:/PreLeisH_scRNAseq/First_Batch/Multi/'
batch1_ids <- c("0114UV", "0104UV", "0114W4", "0104W4")
batch2_dir_path <- 'E:/PreLeisH_scRNAseq/Second_Batch/Multi/'
batch2_ids <- c("0117M6", "0170M6", "0117W4", "0170W4")

flowcyto_markers <- c('FOXP3', 'CXCR3', 'IFNG', 'IL10', 'LAMP1', 'FAS', 
                      'B3GAT1', 'KLRG1', 
                      'LAG3', 'HAVCR2',  'TIGIT', 'PDCD1', 'BTLA')

##################################################################
##                      Load the GE datasets                     #
##################################################################

#Create an empty list filled with two elements (one for each batch)
#each of these two elements contain a list with the length of sample per batch
batch_list <- list(vector(mode = 'list', length = length(batch1_ids)),
                   vector(mode = 'list', length = length(batch2_ids)))

GE_data_loader <- function(batch_dir, batch_ids, batch_num) {
  
  disease_groups <- c("ncVL_HIV", "cVL_HIV", "ncVL_HIV", "cVL_HIV")
  treatment_timepoints <- c('D0', 'D0', 'EOT', 'EOT')
  
  for (i in seq_along(batch_ids)) {
    
    d10x <- Read10X(paste0(batch_dir,batch_ids[i], "/outs/per_sample_outs/", batch_ids[i], "/count/sample_feature_bc_matrix"))
    colnames(d10x) <- paste(batch_ids[i], colnames(d10x), sep = '_')
    colnames(d10x) <- stringr::str_remove(colnames(d10x), "-1")
    
    batch_list[[batch_num]][[i]] <- CreateSeuratObject(counts = d10x, min.cells = 3, min.features = 200)
    
    batch_num_to_add <- data.frame(batch_id = rep(batch_num,
                                                  nrow(batch_list[[batch_num]][[i]]@meta.data)),
                                   disease_group = rep(disease_groups[i],
                                                       nrow(batch_list[[batch_num]][[i]]@meta.data)),
                                   treatment_timepoint = rep(treatment_timepoints[i],
                                                             nrow(batch_list[[batch_num]][[i]]@meta.data)),
                                   stringsAsFactors = F)
    
    row.names(batch_num_to_add) <- row.names(batch_list[[batch_num]][[i]]@meta.data)
    
    batch_list[[batch_num]][[i]] <- AddMetaData(object = batch_list[[batch_num]][[i]],
                                                metadata = batch_num_to_add, col.name = NULL)
    
    batch_list[[batch_num]][[i]] <- NormalizeData(batch_list[[batch_num]][[i]],normalization.method = "LogNormalize", verbose = FALSE)
    batch_list[[batch_num]][[i]] <- FindVariableFeatures(batch_list[[batch_num]][[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE) #according to paper COVID19
    batch_list[[batch_num]][[i]] <- ScaleData(object = batch_list[[batch_num]][[i]], verbose = FALSE)
    
  }
  return(batch_list)
}

batch_list <- GE_data_loader(batch1_dir_path, batch1_ids, 1)
batch_list <- GE_data_loader(batch2_dir_path, batch2_ids, 2)

batch_list <- unlist(batch_list, recursive = TRUE, use.names = TRUE)
names(batch_list) <- c(batch1_ids, batch2_ids)

#################################################################
##                    Integrate all samples                     #
#################################################################

# sample_list <- lapply(X = sample_list, FUN = function(x) {
#   x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
#   x <- subset(x, subset = nFeature_RNA >= 400 & nFeature_RNA <= 3000 & percent.mt < 10)
#   x <- NormalizeData(x,normalization.method = "LogNormalize", scale.factor = 10000)
#   x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000) #according to paper COVID19
# })


long_combined <- merge(batch_list[[1]], batch_list[c(2:length(batch_list))])
rm(batch_list)
gc()

#################################################################
##                        Quality Control                       #
#################################################################

long_combined[["percent.mt"]] <- PercentageFeatureSet(long_combined, pattern = "^MT-")
long_combined[["percent.rb"]] <- PercentageFeatureSet(long_combined, pattern = "^RP[SL]")
long_combined <- subset(long_combined,
                        subset = nFeature_RNA >= 800 & nFeature_RNA <= 3000 & percent.mt > 2.5 & percent.mt < 12.5 & percent.rb > 5)

#################################################################
##                    Load the VDJ datasets                     #
#################################################################

contig_list <- list(vector(mode = 'list', length = length(batch1_ids)), vector(mode = 'list', length = length(batch2_ids)))

vdj_loading <- function(batch_dir, batch_ids, batch_num) {
  
  for (i in seq_along(batch_ids)) {
    
    vdj_contigs <- read.csv(paste0(batch_dir,batch_ids[i], "/outs/per_sample_outs/", batch_ids[i], "/vdj_t/filtered_contig_annotations.csv"))
    contig_list[[batch_num]][[i]] <- vdj_contigs
  }
  return(contig_list)
}

contig_list <- vdj_loading(batch1_dir_path, batch1_ids, 1)
contig_list <- vdj_loading(batch2_dir_path, batch2_ids, 2)
contig_list <- unlist(contig_list, recursive = F, use.names = TRUE)

for (i in seq_along(contig_list)) {
  contig_list[[i]][,"barcode"] <- stringr::str_remove(contig_list[[i]][,"barcode"], "-1")
}

combined <- combineTCR(contig_list, samples = c(c(batch1_ids, batch2_ids)), cells = "T-AB")

##################################################################
##                Add VDJ data to Seurat object                  #
##################################################################

long_combined <- combineExpression(combined, long_combined, proportion = FALSE,
                                     cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))

#################################################################
##            Calculating UMAP and finding clusters             #
#################################################################

long_combined <- NormalizeData(long_combined,normalization.method = "LogNormalize")
long_combined <- FindVariableFeatures(long_combined, selection.method = "vst", nfeatures = 2000) #according to paper COVID19
long_combined <- ScaleData(object = long_combined, verbose = FALSE)
long_combined <- RunPCA(object = long_combined, npcs = 40, verbose = FALSE)
#ElbowPlot(long_combined)
long_combined <- RunUMAP(object = long_combined, reduction = "pca", dims = 1:15)
long_combined <- FindNeighbors(object = long_combined, dims = 1:15, force.recalc = T)
long_combined <- FindClusters(object = long_combined, resolution = 0.55)

#################################################################
##              Save object for Azimuth predictions             #
#################################################################

###Saving object to run Azimuth pipeline (https://satijalab.org/azimuth/)
#saveRDS(long_combined, file = "immune_combined_long_antonio.rds")

##################################################################
##  Add Azimuth celltype predictions to Seurat object metadata   #
##################################################################

azimuth_predictions <- read.delim(here::here('data', 'processed_data', 'azimuth_pred_long.tsv'), row.names = 1)
long_combined <- AddMetaData(object = long_combined, metadata = azimuth_predictions)

#################################################################
##              SingleR for celltype identification             #
#################################################################

# hpca_se <- celldex::HumanPrimaryCellAtlasData()
# singler_pred_immune <- SingleR(test = GetAssayData(long_combined), ref = hpca_se, labels = hpca_se$label.fine)
# long_combined[["SingleR.labels"]] <- singler_pred_immune$labels
# rm(hpca_se, singler_pred_immune)

##################################################################
##                    Identify Cluster Markers                   #
##################################################################
#Could be used for manual cluster celltype identification

# find markers for every cluster compared to all remaining cells, report only the positive ones
# pbmc.markers <- FindAllMarkers(long_combined, only.pos = TRUE, min.pct = 0.25, logfc.threshtop30_markers <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)old = 0.25) %>% mutate(diff = pct.1 - pct.2)
# 
#clus_marker_heatmap <- DoHeatmap(long_combined, features = top30_markers$gene) + NoLegend()

#################################################################
##                      Celltype labeling                       #
#################################################################

#CD4 and CD8 naive could not be clearly separated here through clustering, but it is easy to make the distinction
cd4_naive <- subset(long_combined, idents = 6)
cd4_naive <- WhichCells(cd4_naive, slot = 'counts', expression = CD4 > 0 & CD8A <= 0 & CD8B <= 0)
Idents(long_combined, cells = cd4_naive) <- '17'

###Half of cluster 4 has VDJ data, other half doesn't. Very sizable portion. Makes labeling hard.
split_4 <- subset(long_combined, idents = 4)
split_4 <- WhichCells(subset(split_4, cloneType %in% c(NA, 'NA')))
Idents(long_combined, cells = split_4) <- '18'

#There are some cells in this cluster that have no VDJ data and that are fully CD3 negative, so these are probably NK cells
#Also lot of differential NK markers as compared to other clusters except clus 5 (NK), here almost no diff markers
cd3neg_cells_in_18 <- subset(long_combined, idents = 18)
cd3neg_cells_in_18 <- WhichCells(cd3neg_cells_in_18, slot = 'counts',
                                 expression = CD3E <= 0 & CD3D <= 0 & CD3G <= 0)
Idents(long_combined, cells = cd3neg_cells_in_18) <- '5'
#There are also some cells that are mostly CD3 negative with maybe one exception i.e. CD3D & E negative but CD3G positive
#Probably also NK cells based on differential markers
cd3neg_cells_in_18 <- subset(long_combined, idents = 18)
cd3neg_cells_in_18 <- WhichCells(cd3neg_cells_in_18, slot = 'counts',
                                 expression = ((CD3D <= 0 & CD3E <= 0) | (CD3D <= 0 & CD3G <= 0) | (CD3E <= 0 & CD3G <= 0)))
Idents(long_combined, cells = cd3neg_cells_in_18) <- '5'
#So there are some fully CD3 positive and fully CD8 positive cells, these should be grouped with the CD8 TEM's
cd8pos_cd3pos_cells_in_18 <- subset(long_combined, idents = 18)
cd8pos_cd3pos_cells_in_18 <- WhichCells(cd8pos_cd3pos_cells_in_18, slot = 'counts',
                                        expression = (CD3D > 0 & CD3E > 0 & CD3G > 0) & (CD8A > 0 & CD8B > 0))
Idents(long_combined, cells = cd8pos_cd3pos_cells_in_18) <- '4'
#So there are some mostly CD3 positive and fully CD8 positive cells, these should be grouped with the CD8 TEM's
cd8pos_cd3pos_cells_in_18 <- subset(long_combined, idents = 18)
cd8pos_cd3pos_cells_in_18 <- WhichCells(cd8pos_cd3pos_cells_in_18, slot = 'counts',
                                        expression = ((CD3D > 0 & CD3E > 0) | (CD3D > 0 & CD3G > 0) | (CD3E > 0 & CD3G > 0)) & (CD8A > 0 & CD8B > 0))
Idents(long_combined, cells = cd8pos_cd3pos_cells_in_18) <- '4'
#There are some mostly CD3 positive but fully CD8 negative, these are probably gdT's
cd8neg_cd3pos_cells_in_18 <- subset(long_combined, idents = 18)
cd8neg_cd3pos_cells_in_18 <- WhichCells(cd8neg_cd3pos_cells_in_18, slot = 'counts',
                                        expression = ((CD3D > 0 & CD3E > 0) | (CD3D > 0 & CD3G > 0) | (CD3E > 0 & CD3G > 0)) & ((CD8A <= 0 & CD8B <= 0) | (CD8A > 0 & CD8B <= 0)))
Idents(long_combined, cells = cd8neg_cd3pos_cells_in_18) <- '18' #This ended up being 92% of the cells, and the remaining 8% was just CD8 positive but no real diff markers
#Ended up just including the cluster 18 as gdT's

#gdt_only <- WhichCells(gdt, slot = 'counts', expression = (TRGV2 > 0 | TRGV3 > 0 | TRGV4 > 0 | TRGV8 > 0 | TRGV9 > 0 | TRGV10 > 0 | TRGC1 > 0 | TRGC2 > 0 | TRDV1 > 0 | TRDV3 > 0 | TRDC > 0))

#Setting up new celltype labels
new.cluster.ids <- c("0" = "CD8+ TEM",
                     "1" = "CD8+ TCM",
                     "2" = "CD4+ TCM",
                     "3" = "CD4+ TEM",
                     "4" = "CD8+ TEM",
                     "5" = "NK Cells",
                     "6" = "CD8+ T Naive",
                     "7" = "CD4+ TEM",
                     "8" = "CD14+ Monocytes",
                     "9" = "CD16+ Monocytes",
                     "10" = "B Intermediate",
                     "11" = "B Naive",
                     "12" = "CD14+ Monocytes",
                     "13" = "Proliferating CD8+ T & NK",
                     "14" = "Platelets",
                     "15" = "Erythrocytes",
                     "16" = "Plasmablasts",
                     "17" = "CD4+ T Naive",
                     "18" = 'gdT')

long_combined <- RenameIdents(long_combined, new.cluster.ids)
long_combined[["stashed.ident"]] <- Idents(long_combined)


scGate_modelsDB <- get_scGateDB(here("data", "meta_data", "scGateDB"), force_update = F)

hspc_model <- gating_model(name = 'HSPC', signature = c('SPINK2', 'SOX4', 'CD34'))
platelet_model <- gating_model(name = 'Platelets', signature = c('PPBP', 'PF4', 'TUBB1'))
erythrocytes_model <- gating_model(name = 'Erythrocytes', signature = c('HBB', 'HBA1', 'HBA2', 'HBD', 'ALAS2'))
cd4_prolif_model <- gating_model(model = scGate_modelsDB$human$generic$CD4T, name = 'CD4Prolif', level = 5, signature = c('MKI67', 'STMN1'))
cd8_prolif_model <- gating_model(model = scGate_modelsDB$human$generic$CD8T, name = 'CD8Prolif', level = 5, signature = c('MKI67', 'STMN1'))
NK_prolif_model <- gating_model(model = scGate_modelsDB$human$generic$NK, name = 'Prolif', level = 4, signature = c('MKI67', 'STMN1'))
CD56dim_NK_model <- gating_model(model = scGate_modelsDB$human$generic$NK, name = 'CD56dim', level = 4, signature = c('FCGR3A', 'NCAM1-'))
CD56dim_NK_model <- gating_model(model = CD56dim_NK_model, name = 'noProlif', level = 4, negative = T, signature = c('MKI67', 'STMN1'))
CD56bright_NK_model <- gating_model(model = scGate_modelsDB$human$generic$NK, name = 'CD56bright', level = 4, signature = c('NCAM1', 'XCL1', 'XCL2', 'FCER1G', 'FCGR3A-'))
CD56bright_NK_model <- gating_model(model = CD56bright_NK_model, name = 'noProlif', level = 4, negative = T, signature = c('MKI67', 'STMN1'))

long_combined <- scGate(long_combined, CD56dim_NK_model)
cd56dim_nk_cells <- WhichCells(subset(long_combined, subset = is.pure == 'Pure'))
long_combined <- scGate(long_combined, CD56bright_NK_model)
cd56bright_nk_cells <- WhichCells(subset(long_combined, subset = is.pure == 'Pure'))
long_combined <- scGate(long_combined, NK_prolif_model)
nk_prolif_cells <- WhichCells(subset(long_combined, subset = is.pure == 'Pure'))
long_combined <- scGate(long_combined, cd8_prolif_model)
cd8_prolif_cells <- WhichCells(subset(long_combined, subset = is.pure == 'Pure'))
long_combined <- scGate(long_combined, cd4_prolif_model)
cd4_prolif_cells <- WhichCells(subset(long_combined, subset = is.pure == 'Pure'))

long_combined <- scGate(long_combined, hspc_model)
hspcs <- WhichCells(subset(long_combined, subset = is.pure == 'Pure'))
long_combined <- scGate(long_combined, platelet_model)
platelets <- WhichCells(subset(long_combined, subset = is.pure == 'Pure'))
long_combined <- scGate(long_combined, erythrocytes_model)
erythrocytes <- WhichCells(subset(long_combined, subset = is.pure == 'Pure'))

Idents(long_combined) <- "stashed.ident"
Idents(long_combined, cells = hspcs) <- "HSPCs"
Idents(long_combined, cells = platelets) <- "Platelets"
Idents(long_combined, cells = erythrocytes) <- "Erythrocytes"
Idents(long_combined, cells = cd4_prolif_cells) <- "CD4+ T Proliferating"
Idents(long_combined, cells = cd8_prolif_cells) <- "CD8+ T Proliferating"
Idents(long_combined, cells = nk_prolif_cells) <- "NK Proliferating"
Idents(long_combined, cells = WhichCells(subset(long_combined, idents = 'Proliferating CD8+ T & NK', subset = cloneType %in% c(NA, 'NA')))) <- 'NK Proliferating'
Idents(long_combined, cells = WhichCells(subset(long_combined, idents = 'Proliferating CD8+ T & NK'))) <- 'CD8+ T Proliferating'
Idents(long_combined, cells = cd56dim_nk_cells) <- "CD56dim NK"
Idents(long_combined, cells = cd56bright_nk_cells) <- "CD56bright NK"
Idents(long_combined, cells = WhichCells(subset(long_combined, idents = 'NK Cells'))) <- "CD56dim NK"

Idents(long_combined) <- factor(Idents(long_combined), levels = c("CD4+ T Proliferating", "CD4+ T Naive", "CD4+ TCM", "CD4+ TEM",
                                                                  "CD8+ T Proliferating", "CD8+ T Naive", "CD8+ TCM", "CD8+ TEM",
                                                                  'gdT', "NK Proliferating", "CD56dim NK", "CD56bright NK",
                                                                  "CD14+ Monocytes", "CD16+ Monocytes",
                                                                  "B Naive", "B Intermediate", "Plasmablasts",
                                                                  "Platelets", "Erythrocytes", "HSPCs"))

long_combined$final.ident <- Idents(long_combined)

##################################################################
##                      Cellular composition                     #
##################################################################

long_combined$disease_group <- factor(long_combined$disease_group, levels = c('ncVL_HIV', 'cVL_HIV'))
long_combined$treatment_timepoint <- as.factor(long_combined$treatment_timepoint)
long_combined$treatment_timepoint <- factor(long_combined$treatment_timepoint, levels = c('D0', 'EOT'))
long_combined$cloneType <- as.factor(long_combined$cloneType)
long_combined$cloneType <- factor(long_combined$cloneType, levels = c("Hyperexpanded (100 < X <= 500)", "Large (20 < X <= 100)", "Medium (5 < X <= 20)", "Small (1 < X <= 5)", "Single (0 < X <= 1)", NA))
long_combined$dg_tt <- factor(paste0(long_combined$disease_group, '_', long_combined$treatment_timepoint), levels = c('ncVL_HIV_D0', 'ncVL_HIV_EOT', 'cVL_HIV_D0', 'cVL_HIV_EOT'))

##----------------------------------------------------------------
##                  UMAPs and frequency tables                   -
##----------------------------------------------------------------

celltype_colors <- c("CD4+ T Proliferating" = "#FF9D9A", "CD4+ T Naive" = "#D37295", "CD4+ TCM" ="#FABFD2", "CD4+ TEM" = "#E15759",
                     "CD8+ T Proliferating" = "#499894", "CD8+ T Naive" = "#86BCB6", "CD8+ TCM" = "#A0CBE8", "CD8+ TEM" = "#4E79A7",
                     'gdT' = "#e377c2", "NK Proliferating" = "#F1CE63", "CD56dim NK" = "#FFBE7D", "CD56bright NK" = "#F28E2B",
                     "CD14+ Monocytes" = "#98df8a", "CD16+ Monocytes" = "#2ca02c",
                     "B Naive" = "#D7B5A6", "B Intermediate" = "#B6992D", "Plasmablasts" = "#8c564b",
                     "Platelets" = "#D4A6C8", "Erythrocytes" ="#17becf", "HSPCs" ="#B07AA1")

longitudinal_umap_axes <- DimPlot(long_combined, label = T, repel = T, label.size = 8, pt.size = 1.2,
                                  cols = celltype_colors) +
  ggtitle("Longitudinal comparison VL-HIV patients (VL treatment)") +
  theme_void() + theme(legend.position = 'none',
                     panel.grid = element_blank(),
                     plot.title = element_text(hjust = 0.5, family = text_font, size = 36, face = 'bold')) + 
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 7.5)))

# freq_table_dg <- tableGrob(rbind(table(long_combined@active.ident, long_combined$disease_group),`Total Cells` = colSums(table(long_combined@active.ident, long_combined$disease_group))), theme = ttheme('classic', base_size = 15))
# freq_table_tt <- tableGrob(rbind(table(long_combined@active.ident, long_combined$treatment_timepoint),`Total Cells` = colSums(table(long_combined@active.ident, long_combined$treatment_timepoint))), theme = ttheme('classic', base_size = 15))

##---------------------------------------------------------------
##                  Cellular composition graphs                 -
##---------------------------------------------------------------

patient_table <- as_tibble(rbind(table(long_combined@active.ident,
                                       long_combined$orig.ident)),
                           rownames = 'Celltype') %>%
  pivot_longer(!Celltype, names_to = 'PatientID', values_to = 'Count_per_patient') %>%
  inner_join(., long_combined[[]], by = c('PatientID' = 'orig.ident')) %>%
  select(c(Celltype, PatientID, disease_group, treatment_timepoint, Count_per_patient)) %>% 
  distinct()

total_per_patient_table <- patient_table %>% 
  group_by(PatientID) %>% 
  summarise(Total = sum(Count_per_patient))

patient_table <- inner_join(patient_table, total_per_patient_table, by = c("PatientID")) %>%
  mutate(Percentage = (Count_per_patient / Total) * 100, Condition = paste0(disease_group, '_', treatment_timepoint)) %>%
  group_by(Celltype, Condition) %>% mutate(Mean_Percentage = mean(Percentage))

mean_table <- patient_table %>% select(., c(Celltype, Condition, Mean_Percentage)) %>% distinct()

temp_long_combined <- long_combined
temp_long_combined$condition <- paste0(temp_long_combined$disease_group, '_', temp_long_combined$treatment_timepoint)
Idents(temp_long_combined) <- 'condition'

total_cellcounts_per_group <- as_tibble(rbind(table(Idents(temp_long_combined)))) %>% select(c(1,3,2,4)) %>% 
  pivot_longer(everything(), names_to = 'condition', values_to = 'Cells')
rm(temp_long_combined)
gc()

stacked_barchart_longitudinal <- ggplot(mean_table,
                                        aes(factor(Condition, levels = c("ncVL_HIV_D0", "ncVL_HIV_EOT", "cVL_HIV_D0", "cVL_HIV_EOT")), Mean_Percentage,
                                            fill = factor(Celltype, levels = c("CD4+ T Proliferating", "CD4+ T Naive", "CD4+ TCM", "CD4+ TEM",
                                                                                          "CD8+ T Proliferating", "CD8+ T Naive", "CD8+ TCM", "CD8+ TEM",
                                                                                          'gdT', "NK Proliferating", "CD56dim NK", "CD56bright NK",
                                                                                          "CD14+ Monocytes", "CD16+ Monocytes",
                                                                                          "B Naive", "B Intermediate", "Plasmablasts",
                                                                                          "Platelets", "Erythrocytes", "HSPCs")))) +
  geom_bar(stat = "identity") +
  labs(y = 'Percentage') +
  scale_x_discrete(expand = c(0.08, 0.23),
                   labels = paste0(levels(factor(mean_table$Condition, levels = c("ncVL_HIV_D0", "ncVL_HIV_EOT", "cVL_HIV_D0", "cVL_HIV_EOT"))), ' (n = ', total_cellcounts_per_group$Cells, ')'))+
  scale_y_continuous(expand = c(0.00, 0.01)) +
  ggtitle("Cellular composition of VL-HIV patients") +
  theme_classic() + theme(axis.title = element_text(size = 32, family = text_font),
                          axis.title.x = element_blank(),
                          axis.text = element_text(size = 26, family = text_font, color = 'black'),
                          axis.ticks.length = unit(.2, "cm"),
                          legend.position = 'right',
                          legend.title = element_blank(),
                          legend.text = element_text(size = 26, family = text_font),
                          plot.title = element_text(hjust = 0.5, family = text_font, size = 36, face = 'bold')) + 
  guides(fill = guide_legend(ncol = 1, override.aes = list(size = 7.5))) +
  ggpubr::rotate_x_text(angle = 40) +
  scale_fill_manual(values = unname(celltype_colors), limits = levels(factor(mean_table$Celltype, levels = c("CD4+ T Proliferating", "CD4+ T Naive", "CD4+ TCM", "CD4+ TEM",
                                                                                                             "CD8+ T Proliferating", "CD8+ T Naive", "CD8+ TCM", "CD8+ TEM",
                                                                                                             'gdT', "NK Proliferating", "CD56dim NK", "CD56bright NK",
                                                                                                             "CD14+ Monocytes", "CD16+ Monocytes",
                                                                                                             "B Naive", "B Intermediate", "Plasmablasts",
                                                                                                             "Platelets", "Erythrocytes", "HSPCs"))))

#-----------Saving cross_sectional + longitudinal at the same time
layout <- "
AAAABB
CCCCDD
"

arranged_plot <- cross_sect_umap_axes + stacked_barchart_cross_sectional +
  longitudinal_umap_axes + stacked_barchart_longitudinal +
  #celltype_barchart_cross_sectional + celltype_barchart_longitudinal +
  plot_layout(design = layout) & plot_annotation(tag_levels = c('a')) & theme(plot.tag = element_text(size = 36, face = 'bold'))

ggsave(here("analyses", "final_output", 'Figure3_UMAP_Composition.pdf'), arranged_plot, height = 17.5, width = 30, device = cairo_pdf)
ggsave(here("analyses", "final_output", 'Figure3_UMAP_Composition.png'), arranged_plot, height = 17.5, width = 30, type = 'cairo-png')

#################################################################
##                    Differential Expression                   #
#################################################################

##----------------------------------------------------------------
##                Differential Expression CD4+ T                 -
##----------------------------------------------------------------

allcd4 <- names(long_combined@active.ident[grepl('CD4', long_combined@active.ident)])
cd4_subset <- subset(long_combined, cells = allcd4)

Idents(cd4_subset) <- 'treatment_timepoint'
CD4_D0_cVL_vs_ncVL <- FindMarkers(cd4_subset, ident.1 = 'ncVL_HIV', ident.2 = 'cVL_HIV', group.by = 'disease_group', subset.ident = 'D0', min.pct = 0.1, logfc.threshold = 0.25, only.pos = F) %>% mutate(gene = rownames(.), diff = pct.1 - pct.2)
CD4_EOT_cVL_vs_ncVL <- FindMarkers(cd4_subset, ident.1 = 'ncVL_HIV', ident.2 = 'cVL_HIV', group.by = 'disease_group', subset.ident = 'EOT', min.pct = 0.1, logfc.threshold = 0.25, only.pos = F) %>% mutate(gene = rownames(.), diff = pct.1 - pct.2)

Idents(cd4_subset) <- 'disease_group'
CD4_cVL_EOT_vs_D0 <- FindMarkers(cd4_subset, ident.1 = 'EOT', ident.2 = 'D0', group.by = 'treatment_timepoint', subset.ident = 'ncVL_HIV', min.pct = 0.1, logfc.threshold = 0.25, only.pos = F) %>% mutate(gene = rownames(.), diff = pct.1 - pct.2)
CD4_ncVL_EOT_vs_D0 <- FindMarkers(cd4_subset, ident.1 = 'EOT', ident.2 = 'D0', group.by = 'treatment_timepoint', subset.ident = 'cVL_HIV', min.pct = 0.1, logfc.threshold = 0.25, only.pos = F) %>% mutate(gene = rownames(.), diff = pct.1 - pct.2)


##----------------------------------------------------------------
##                Differential Expression CD8+ T                 -
##----------------------------------------------------------------

allcd8 <- names(long_combined@active.ident[grepl('CD8', long_combined@active.ident)])
cd8_subset <- subset(long_combined, cells = allcd8)

Idents(cd8_subset) <- 'treatment_timepoint'
CD8_D0_cVL_vs_ncVL <- FindMarkers(cd8_subset, ident.1 = 'ncVL_HIV', ident.2 = 'cVL_HIV', group.by = 'disease_group', subset.ident = 'D0', min.pct = 0.1, logfc.threshold = 0.25, only.pos = F) %>% mutate(gene = rownames(.), diff = pct.1 - pct.2)
CD8_EOT_cVL_vs_ncVL <- FindMarkers(cd8_subset, ident.1 = 'ncVL_HIV', ident.2 = 'cVL_HIV', group.by = 'disease_group', subset.ident = 'EOT', min.pct = 0.1, logfc.threshold = 0.25, only.pos = F) %>% mutate(gene = rownames(.), diff = pct.1 - pct.2)

Idents(cd8_subset) <- 'disease_group'
CD8_cVL_EOT_vs_D0 <- FindMarkers(cd8_subset, ident.1 = 'EOT', ident.2 = 'D0', group.by = 'treatment_timepoint', subset.ident = 'ncVL_HIV', min.pct = 0.1, logfc.threshold = 0.25, only.pos = F) %>% mutate(gene = rownames(.), diff = pct.1 - pct.2)
CD8_ncVL_EOT_vs_D0 <- FindMarkers(cd8_subset, ident.1 = 'EOT', ident.2 = 'D0', group.by = 'treatment_timepoint', subset.ident = 'cVL_HIV', min.pct = 0.1, logfc.threshold = 0.25, only.pos = F) %>% mutate(gene = rownames(.), diff = pct.1 - pct.2)


##################################################################
##        Generating DE VolcanoPlots & Enrichment Networks       #
##################################################################

enrichment_plotter <- function(marker_list, x_axis_text_size) {
  
  m_t2g <- msigdbr::msigdbr(species = "Homo sapiens", category = c('H')) %>% 
    dplyr::select(gs_name, entrez_gene) %>% dplyr::distinct(gs_name, entrez_gene)
  
  go_m_t2g <- msigdbr::msigdbr(species = "Homo sapiens", category = c('C5'), subcategory = 'BP') %>% filter(., gs_subcat != 'HPO') %>%
    dplyr::select(gs_name, entrez_gene) %>% dplyr::distinct(gs_name, entrez_gene)
  
  m_t2g <- bind_rows(m_t2g, go_m_t2g)
  
  go_jointable <- msigdbr::msigdbr(species = "Homo sapiens", category = c('C5'), subcategory = 'BP') %>% filter(., gs_subcat != 'HPO')
  
  enrich_result <- lapply(seq_along(marker_list), FUN = function(x) {
    tmp <- filter(marker_list[[x]], p_val_adj < 0.05) %>% filter(., !str_starts(gene, '^RP')) %>% filter(., !str_starts(gene, '^MT'))
    
    gene_list <- row.names(tmp)
    marker_expression <- tmp$avg_log2FC
    names(marker_expression) <- gene_list
    if (length(marker_expression) > 1) {
      gene_aliases <- AnnotationDbi::select(org.Hs.eg.db, gene_list, c("ENTREZID"), "ALIAS")
      gene_aliases <- inner_join(gene_aliases, tibble::enframe(marker_expression, name = 'gene', value = 'fold_change'),
                                 by = c('ALIAS' = 'gene'))
      
      enrich_ready <- gene_aliases$fold_change
      names(enrich_ready) <- gene_aliases$ENTREZID
      
      enriched <- enricher(names(enrich_ready), TERM2GENE=m_t2g)
      
      enricher_result <- enriched@result
      enricher_result <- enricher_result %>% mutate(go_jointable[match(enricher_result$ID, go_jointable$gs_name), c('gs_name', 'gs_subcat', 'gs_exact_source')])
      hallmark_result <- filter(enricher_result, is.na(gs_name)) %>% select(!c('gs_name', 'gs_subcat', 'gs_exact_source'))
      go_result <- filter(enricher_result, !is.na(gs_name))
      go_result <- mutate(go_result, ONTOLOGY = str_split(go_result$gs_subcat, ':', simplify = T)[,2])
      go_result$ID <- go_result$gs_exact_source
      go_result <- go_result %>% select(c("ONTOLOGY","ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count" ))
      rownames(go_result) <- go_result$ID
      go_result_bp <- filter(go_result, ONTOLOGY == 'BP')

      enriched@ontology <- 'BP'
      enriched@keytype <- "ENTREZID"
      enriched@organism <- "Homo sapiens"

      enriched@result <- go_result_bp
      enriched <- clusterProfiler::simplify(enriched, cutoff = 0.6)
      simplified_bp <- enriched@result

      total_simplified <- simplified_bp %>% select(!ONTOLOGY) %>% mutate(ID = Description)
      rownames(total_simplified) <- total_simplified$ID
      total_simplified <- bind_rows(total_simplified, hallmark_result)

      enriched@result <- total_simplified
      
      enriched@result <- enriched@result[order(enriched@result$p.adjust),]
      
      enriched <- setReadable(enriched, org.Hs.eg.db, keyType = "ENTREZID")

      hp <- heatplot(enriched, showCategory = 15, foldChange = enrich_ready) +
        scale_y_discrete(label = function(y) gsub("_", " ", stringr::str_trunc(y, 40))) +
        ggtitle(names(marker_list)[x]) +
        theme_bw() + 
        #scale_fill_gradient2(high = "#9a133d", low = "#1a318b", mid = "#ffffff") +
        scale_fill_gradient2(name = "FoldChange", low = "#1F77B4FF", high = "#B2182B", mid = 'white') +
        #scale_fill_distiller(palette = 'RdBu') +
        theme(plot.title = element_text(hjust = 0.5, family = text_font, size = 36, face = 'bold'),
              text = element_text(family = text_font),
              axis.title = element_blank(),
              axis.text.x = element_text(family = text_font, size = x_axis_text_size[x], face = 'bold'),
              axis.text.y = element_text(family = text_font, size = 17, face = 'bold'),
              legend.title = element_text(family = text_font, size = 21, face = 'bold'),
              legend.text = element_text(family = text_font, size = 18)) + ggpubr::rotate_x_text(angle = 90)
      
      dp <- mutate(enriched, qscore = -log(p.adjust, base = 10)) %>%
        barplot(x = 'qscore', showCategory = 15) +
        #dotplot(enriched, showCategory = nrow(filter(enriched@result, p.adjust < 0.05))) +
        scale_y_discrete(label = function(y) gsub("_", " ", stringr::str_trunc(y, 40))) +
        ggtitle(names(marker_list)[x]) +
        theme_bw() + 
        theme(plot.title = element_text(hjust = 0.5, family = text_font, size = 36, face = 'bold'),
              text = element_text(family = text_font),
              axis.title.x = element_text(family = text_font, size = 20, face = 'bold'),
              axis.text.x = element_text(family = text_font, size = x_axis_text_size[x], face = 'bold'),
              axis.title.y = element_blank(),
              axis.text.y = element_text(family = text_font, size = 14, face = 'bold'),
              legend.title = element_text(family = text_font, size = 21, face = 'bold'),
              legend.text = element_text(family = text_font, size = 18))
      
      #hp$layers[[1]]$aes_params$colour <- 'black'
      
      return(list(hp, dp))
    } else {
      return(NULL)
    }
    
  })
  names(enrich_result) <- names(marker_list)
  return(enrich_result)
  
}

cd4_cVL_vs_ncVL_list <- list('CD4+ T Non-Chronic VL-HIV vs Chronic VL-HIV at D0' = CD4_D0_cVL_vs_ncVL,
                        'CD4+ T Non-Chronic VL-HIV vs Chronic VL-HIV at EOT' = CD4_EOT_cVL_vs_ncVL)

cd4_EOT_vs_D0_list <- list('CD4+ T Non-Chronic VL-HIV EOT vs D0' = CD4_cVL_EOT_vs_D0,
                           'CD4+ T Chronic VL-HIV EOT vs D0' = CD4_ncVL_EOT_vs_D0)

cd8_cVL_vs_ncVL_list <- list('CD8+ T Non-Chronic VL-HIV vs Chronic VL-HIV at D0' = CD8_D0_cVL_vs_ncVL,
                            'CD8+ T Non-Chronic VL-HIV vs Chronic VL-HIV at EOT' = CD8_EOT_cVL_vs_ncVL)

cd8_EOT_vs_D0_list <- list('CD8+ T Non-Chronic VL-HIV EOT vs D0' = CD8_cVL_EOT_vs_D0,
                           'CD8+ T Chronic VL-HIV EOT vs D0' = CD8_ncVL_EOT_vs_D0)




cd4_cVL_vs_ncVL_enrichplots <- enrichment_plotter(cd4_cVL_vs_ncVL_list, c(16, 14))
cd8_cVL_vs_ncVL_enrichplots <- enrichment_plotter(cd8_cVL_vs_ncVL_list, c(13, 13))

cd4_EOT_vs_D0_enrichplots <- enrichment_plotter(cd4_EOT_vs_D0_list, c(17, 16))
cd8_EOT_vs_D0_enrichplots <- enrichment_plotter(cd8_EOT_vs_D0_list, c(16, 17))


layout <- "
AAAAAA
BBBBBB
CCCCCC
DDDDDD
"

cd4_enrichment_arranged <- cd4_cVL_vs_ncVL_enrichplots$`CD4+ T Non-Chronic VL-HIV vs Chronic VL-HIV at D0`[[1]] + cd4_cVL_vs_ncVL_enrichplots$`CD4+ T Non-Chronic VL-HIV vs Chronic VL-HIV at EOT`[[1]] +
  cd4_EOT_vs_D0_enrichplots$`CD4+ T Non-Chronic VL-HIV EOT vs D0`[[1]] + cd4_EOT_vs_D0_enrichplots$`CD4+ T Chronic VL-HIV EOT vs D0`[[1]] +
  plot_layout(design = layout) + plot_annotation(tag_levels = c('a')) & theme(plot.tag = element_text(size = 40, face = 'bold'))

cd8_enrichment_arranged <- cd8_cVL_vs_ncVL_enrichplots$`CD8+ T Non-Chronic VL-HIV vs Chronic VL-HIV at D0`[[1]] + cd8_cVL_vs_ncVL_enrichplots$`CD8+ T Non-Chronic VL-HIV vs Chronic VL-HIV at EOT`[[1]] +
  cd8_EOT_vs_D0_enrichplots$`CD8+ T Non-Chronic VL-HIV EOT vs D0`[[1]] + cd8_EOT_vs_D0_enrichplots$`CD8+ T Chronic VL-HIV EOT vs D0`[[1]] +
  plot_layout(design = layout) + plot_annotation(tag_levels = c('a')) & theme(plot.tag = element_text(size = 40, face = 'bold'))

cd4_cd8_cVL_vs_ncVL_arranged <- cd4_cVL_vs_ncVL_enrichplots$`CD4+ T Non-Chronic VL-HIV vs Chronic VL-HIV at D0`[[1]] + cd4_cVL_vs_ncVL_enrichplots$`CD4+ T Non-Chronic VL-HIV vs Chronic VL-HIV at EOT`[[1]] +
  cd8_cVL_vs_ncVL_enrichplots$`CD8+ T Non-Chronic VL-HIV vs Chronic VL-HIV at D0`[[1]] + cd8_cVL_vs_ncVL_enrichplots$`CD8+ T Non-Chronic VL-HIV vs Chronic VL-HIV at EOT`[[1]] +
  plot_layout(design = layout) + plot_annotation(tag_levels = c('a')) & theme(plot.tag = element_text(size = 40, face = 'bold'))

cd4_cd8_EOT_vs_D0_arranged <- cd4_EOT_vs_D0_enrichplots$`CD4+ T Non-Chronic VL-HIV EOT vs D0`[[1]] + cd4_EOT_vs_D0_enrichplots$`CD4+ T Chronic VL-HIV EOT vs D0`[[1]] +
  cd8_EOT_vs_D0_enrichplots$`CD8+ T Non-Chronic VL-HIV EOT vs D0`[[1]] + cd8_EOT_vs_D0_enrichplots$`CD8+ T Chronic VL-HIV EOT vs D0`[[1]] +
  plot_layout(design = layout) + plot_annotation(tag_levels = c('a')) & theme(plot.tag = element_text(size = 40, face = 'bold'))

cd4_cd8_EOT_vs_D0_score_arranged <- cd4_EOT_vs_D0_enrichplots$`CD4+ T Non-Chronic VL-HIV EOT vs D0`[[2]] + cd4_EOT_vs_D0_enrichplots$`CD4+ T Chronic VL-HIV EOT vs D0`[[2]] +
  cd8_EOT_vs_D0_enrichplots$`CD8+ T Non-Chronic VL-HIV EOT vs D0`[[2]] + cd8_EOT_vs_D0_enrichplots$`CD8+ T Chronic VL-HIV EOT vs D0`[[2]] +
  plot_layout(design = layout) + plot_annotation(tag_levels = c('a')) & theme(plot.tag = element_text(size = 40, face = 'bold'))


ggsave(here("analyses", "final_output", "Figure8_SingleCell_DE.pdf"), cd4_cd8_EOT_vs_D0_arranged, dev = Cairo, type = 'pdf', height = 22, width = 30)
ggsave(here("analyses", "final_output", "Figure8_SingleCell_DE.png"), cd4_cd8_EOT_vs_D0_arranged, type = 'cairo-png', height = 22, width = 30)

ggsave(here("analyses", "final_output", "supplemental", "SuppFigure12_SingleCell_DEScore.pdf"), cd4_cd8_EOT_vs_D0_score_arranged, dev = Cairo, type = 'pdf', height = 20, width = 27.5)
ggsave(here("analyses", "final_output", "supplemental", "SuppFigure12_SingleCell_DEScore.png"), cd4_cd8_EOT_vs_D0_score_arranged, type = 'cairo-png', height = 20, width = 27.5)

ggsave(here("analyses", "final_output", "supplemental", "SuppFigure13_SingleCell_DE.pdf"), cd4_cd8_cVL_vs_ncVL_arranged, dev = Cairo, type = 'pdf', height = 22, width = 30)
ggsave(here("analyses", "final_output", "supplemental", "SuppFigure13_SingleCell_DE.png"), cd4_cd8_cVL_vs_ncVL_arranged, type = 'cairo-png', height = 22, width = 30)

##################################################################
##                          TCR Analysis                         #
##################################################################

longitudinal_split_axes_umap <- DimPlot(long_combined, label = F, repel = F, label.size = 7.5, pt.size = 1.2,
                                        cols = celltype_colors, split.by = 'dg_tt') +
  theme_void() + 
  theme(plot.title = element_blank(), legend.text = element_text(size=27), text = element_text(family = text_font),
        strip.text = element_text(size = 27, family = text_font, face = 'bold')) +
  guides(color = guide_legend(override.aes = list(size = 7.5)))

longitudinal_clonotype_axes <- DimPlot(long_combined, group.by = 'cloneType', split.by = 'dg_tt', pt.size = 1.2) + 
  scale_colour_manual(values = colorblind_vector(5), na.value = 'grey') + 
  theme_void() + 
  theme(plot.title = element_blank(), legend.text = element_text(size=27), text = element_text(family = text_font),
                       strip.text = element_text(size = 27, family = text_font, face = 'bold')) +
  guides(color = guide_legend(override.aes = list(size = 7.5)))


layout <- "
AAA
BBB
"

combined_UMAPs <- longitudinal_split_axes_umap + longitudinal_clonotype_axes + plot_layout(design = layout)

# ggsave(here("analyses", "final_output", 'SingleCell_Figure7_ClonoTypeOverlay.png'), combined_UMAPs, height = 20, width = 30, type = 'cairo-png')

long_combined$cluster <- Idents(long_combined)
long_tcr_by_cluster <- expression2List(long_combined, split.by = NULL)
long_tcr_by_cluster <- long_tcr_by_cluster[(grepl("T", names(long_tcr_by_cluster)) & names(long_tcr_by_cluster) != c('gdT'))]
long_tcr_by_individual <- expression2List(long_combined, split.by = 'orig.ident')
long_tcr_by_disease_group <- expression2List(long_combined, split.by = 'disease_group')
long_tcr_by_treatment_timepoint <- expression2List(long_combined, split.by = 'treatment_timepoint')

##----------------------------------------------------------------
##                  Clonotype expansion alluvial                 -
##----------------------------------------------------------------

NonChronic_subset <- subset(long_combined, subset = dg_tt %in% c("ncVL_HIV_D0", "ncVL_HIV_EOT"))
Chronic_subset <- subset(long_combined, subset = dg_tt %in% c("cVL_HIV_D0", "cVL_HIV_EOT"))

alluvial_colors <- c("#660d20", "#a00e00", "#eb1e2c", "#bD0a36", "#F83B4A", 
                      "#329941", "#459471", "#51b364", "#90DB2C", "#A4D56E", 
                      "#1170aa", "#3C8BCB", "#5fa2ce", "#98d9e4", "#a3cce9",
                      "#E9530F", "#e9a00e", "#fbe183", "#fc7D0b", "#eebe04",
                      "#8E79CB", "#C04FFA",  "#ED21CD", "#C15DCC", "#E842A0",
                      "#472c0b", "#7e5522", "#a97f2f", "#8c564b", "#D7B5A6")

NonChronic_alluvial <- compareClonotypes(NonChronic_subset, split.by = 'dg_tt', numbers = 17, chain = 'both', cloneCall = 'aa') + #26 for top 25, 13 for new update. New update: 17 for top 30
  coord_cartesian(ylim = c(0, 0.45)) +
 scale_fill_manual(values = alluvial_colors) +
  theme(axis.text.y = element_text(size = 18, family = text_font, color = 'black'),
        axis.title.y = element_text(size = 20, family = text_font),
        axis.text.x = element_text(color='black', size=18, family=text_font, face = 'bold'),
        axis.title.x = element_blank(),
        legend.title = element_text(size = 18, family = text_font, face = 'bold', hjust = 0.5),
        legend.text = element_text(size = 11, family = text_font),
        legend.position = 'bottom') +
  guides(fill=guide_legend(ncol=3, title.position = 'top'))

Chronic_alluvial <- scRepertoire::compareClonotypes(Chronic_subset, split.by = 'dg_tt', numbers = 21, chain = 'both', cloneCall = 'aa') + #33 for top 25, 17 for new update. New update: 21 for top 30
  coord_cartesian(ylim = c(0, 0.45)) +
  scale_fill_manual(values = alluvial_colors) +
  theme(axis.text.y = element_text(size = 18, family = text_font, color = 'black'),
        axis.title.y = element_text(size = 20, family = text_font),
        axis.text.x = element_text(color='black', size=18, family=text_font, face = 'bold'),
        axis.title.x = element_blank(),
        legend.title = element_text(size = 18, family = text_font, face = 'bold', hjust = 0.5),
        legend.text = element_text(size = 11, family = text_font),
        legend.position = 'bottom') +
  guides(fill=guide_legend(ncol=3, title.position = 'top'))
  
#arranged_alluvial <- NonChronic_alluvial + Chronic_alluvial + plot_annotation(title = "T cell clonotype expansion in VL-HIV patients after treatment", theme = theme(plot.title = element_text(size = 30, family = text_font, face = 'bold', hjust = 0.5)))

NonChronic_alluvial_data <- NonChronic_alluvial$data %>% pivot_wider(. ,names_from = Sample, values_from = Proportion) %>%
  replace(is.na(.), 0)

Chronic_alluvial_data <- Chronic_alluvial$data %>% pivot_wider(. ,names_from = Sample, values_from = Proportion) %>%
  replace(is.na(.), 0) #%>% pivot_longer(., cols = c("cVL_HIV_D0", "cVL_HIV_EOT"), names_to = 'Sample', values_to = "Proportion")

long_combined <- highlightClonotypes(long_combined, cloneCall = 'aa', sequence = NonChronic_alluvial$data$Clonotypes)

long_combined$aahigh <- ifelse(is.na(long_combined$highlight), NA, long_combined$CTaa)

NonChronic_highlighted_subset <- subset(long_combined, subset = aahigh %in% c(NA, 'NA'), invert = T)

NonChronic_highlighted <- occupiedscRepertoire(NonChronic_highlighted_subset, proportion = F) +
  ggtitle("Non-Chronic VL-HIV: Phenotype of top 30 clonotypes") +
  scale_fill_manual(values = colorblind_vector(5)) +
  scale_x_discrete(limits = c("CD4+ TEM", "CD8+ T Proliferating", "CD8+ TCM", "CD8+ TEM")) +
  theme(axis.text.y = element_text(size = 18, family = text_font, color = 'black'),
        axis.title.y = element_text(size = 20, family = text_font),
        axis.text.x = element_text(color='black', size=18, angle = 60, vjust = 0.5, family=text_font),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 18, family = text_font),
        plot.title = element_text(hjust = 0.5, family = text_font, size = 30, face = 'bold'))
  
# NonChronic_highlighted <- DimPlot(long_combined, group.by = 'aahigh', pt.size = 1.2,
#                              cols = alluvial_colors) +
#   stat_ellipse(data = cd4_subset@meta.data, aes(x = Embeddings(cd4_subset[['umap']])[, "UMAP_1"], y = Embeddings(cd4_subset[['umap']])[, "UMAP_2"], label = 'CD4+ T cells', linewidth = 1, size = 20), geom = 'labelpath', hjust = 0.25, level = 0.93) +
#   stat_ellipse(data = cd8_subset@meta.data, aes(x = Embeddings(cd8_subset[['umap']])[, "UMAP_1"], y = Embeddings(cd8_subset[['umap']])[, "UMAP_2"], label = 'CD8+ T cells', linewidth = 1, size = 20), geom = 'labelpath', hjust = 0.75, level = 0.85) +
#   ggtitle("Non-Chronic VL-HIV") +
#   theme_void() + 
#   theme(legend.position = 'none',
#         plot.title = element_text(hjust = 0.5, family = text_font, size = 28, face = 'bold'))# + guides(color = guide_legend(nrow = 5, override.aes = list(size = 7.5)))

long_combined <- highlightClonotypes(long_combined, cloneCall = 'aa', sequence = Chronic_alluvial$data$Clonotypes)
long_combined$aahigh <- ifelse(is.na(long_combined$highlight), NA, long_combined$CTaa)

Chronic_highlighted_subset <- subset(long_combined, subset = aahigh %in% c(NA, 'NA'), invert = T)

Chronic_highlighted <- occupiedscRepertoire(Chronic_highlighted_subset, proportion = F) +
  ggtitle(" Chronic VL-HIV: Phenotype of top 30 clonotypes") +
  scale_fill_manual(values = colorblind_vector(5)[2:5]) +
  scale_x_discrete(limits = c("CD4+ TEM", "CD8+ T Proliferating", "CD8+ TCM", "CD8+ TEM")) +
  theme(axis.text.y = element_text(size = 18, family = text_font, color = 'black'),
        axis.title.y = element_text(size = 20, family = text_font),
        axis.text.x = element_text(color='black', size=18, angle = 60, vjust = 0.5, family=text_font),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 18, family = text_font),
        plot.title = element_text(hjust = 0.5, family = text_font, size = 30, face = 'bold'))


# Chronic_highlighted <- DimPlot(long_combined, group.by = 'aahigh', pt.size = 1.2,
#                                  cols = alluvial_colors) +
#   stat_ellipse(data = cd4_subset@meta.data, aes(x = Embeddings(cd4_subset[['umap']])[, "UMAP_1"], y = Embeddings(cd4_subset[['umap']])[, "UMAP_2"], label = 'CD4+ T cells', linewidth = 1, size = 20), geom = 'labelpath', hjust = 0.25, level = 0.93) +
#   stat_ellipse(data = cd8_subset@meta.data, aes(x = Embeddings(cd8_subset[['umap']])[, "UMAP_1"], y = Embeddings(cd8_subset[['umap']])[, "UMAP_2"], label = 'CD8+ T cells', linewidth = 1, size = 20), geom = 'labelpath', hjust = 0.75, level = 0.85) +
#   ggtitle("Chronic VL-HIV") +
#   theme_void() +
#   theme(legend.position = 'none',
#         plot.title = element_text(hjust = 0.5, family = text_font, size = 28, face = 'bold')) #+ guides(color = guide_legend(nrow = 5, override.aes = list(size = 7.5)))

layout <- "
AAABBB
CCCDDD
"

arranged_clonotype <- NonChronic_alluvial + Chronic_alluvial + NonChronic_highlighted + Chronic_highlighted +
  plot_annotation(title = "T cell clonotype expansion in VL-HIV patients after treatment", tag_levels = c('a'), theme = theme(plot.title = element_text(size = 30, family = text_font, face = 'bold', hjust = 0.5))) &
  theme(plot.tag = element_text(size = 36, face = 'bold'))

ggsave(here("analyses", "final_output", 'Figure9_ClonoTypes.pdf'), arranged_clonotype, height = 15, width = 27.5, dev = Cairo, type = 'pdf')
ggsave(here("analyses", "final_output", 'Figure9_ClonoTypes.png'), arranged_clonotype, height = 15, width = 27.5, type = 'cairo-png')

#arranged_clonotype <- arranged_alluvial + NonChronic_highlighted + Chronic_highlighted

# NonChronic_alluvial_foldchanges <- NonChronic_alluvial$data %>% pivot_wider(. ,names_from = Sample, values_from = Proportion) %>% 
#   group_by(Clonotypes) %>% 
#   mutate(Log2FC = log2(ncVL_HIV_EOT / ncVL_HIV_D0)) %>% 
#   pivot_longer(., cols = c("ncVL_HIV_D0", "ncVL_HIV_EOT"), names_to = 'Sample', values_to = "Proportion")
# 
# Chronic_alluvial_foldchanges <- Chronic_alluvial$data %>% pivot_wider(. ,names_from = Sample, values_from = Proportion) %>%
#   replace(is.na(.), 0) %>%
#   group_by(Clonotypes) %>%
#   mutate(Log2FC = log2((cVL_HIV_EOT + 1) / (cVL_HIV_D0 + 1))) %>%
#   pivot_longer(., cols = c("cVL_HIV_D0", "cVL_HIV_EOT"), names_to = 'Sample', values_to = "Proportion")


###TEST

# occupiedscRepertoire(subset(long_combined, CTaa %in% Chronic_alluvial_data$Clonotypes), proportion = T) +
#   labs(y = 'Proportion') +
#   ggtitle("T cell clonotype expansion by phenotype") +
#   scale_fill_manual(values = colorblind_vector(5)) +
#   theme(axis.text.y = element_text(size = 18, family = text_font, color = 'black'),
#         axis.title.y = element_text(size = 20, family = text_font),
#         axis.text.x = element_text(color='black', size=18, angle = 60, vjust = 0.5, family=text_font),
#         axis.title.x = element_blank(),
#         legend.title = element_blank(),
#         legend.text = element_text(size = 18, family = text_font),
#         plot.title = element_text(hjust = 0.5, family = text_font, size = 30, face = 'bold')) + geom_text(aes(label = total, y = 1.025), fontface = 'bold')

# test_subset <- subset(long_combined, CTaa %in% c(NonChronic_alluvial_data$Clonotypes, Chronic_alluvial_data$Clonotypes))
# Idents(test_subset) <- "disease_group"
# test <- FindMarkers(test_subset, ident.1 = 'ncVL_HIV', ident.2 = 'cVL_HIV', group.by = 'disease_group', min.pct = 0.1, logfc.threshold = 0.25, only.pos = F) %>% mutate(gene = rownames(.), diff = pct.1 - pct.2)
# test_enrichplot <- enrichment_plotter(list(test), c(14))
