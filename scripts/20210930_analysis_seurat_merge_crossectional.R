#################################################################
##                    Loading in the packages                   #
#################################################################

suppressMessages(library(here))
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

##################################################################
##                    Setting global variables                   #
##################################################################

text_font <- 'Roboto Condensed'

colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", 
                                        "#7AC5FF", "#C6FDEC", "#0348A6"))

batch1_dir_path <- 'E:/PreLeisH_scRNAseq/First_Batch/Multi/'
batch1_ids <- c("HEC1", "0001M3", "0123M3", "0114UV", "0104UV")
batch2_dir_path <- 'E:/PreLeisH_scRNAseq/Second_Batch/Multi/'
batch2_ids <- c("HEC2", "0090M3", "0120M3", "0117M6", "0170M6")

flowcyto_markers <- c('FOXP3', 'CXCR3', 'IFNG', 'IL10', 'LAMP1', 'FAS', 
                      'B3GAT1', 'KLRG1', 
                      'LAG3', 'HAVCR2',  'TIGIT', 'PDCD1', 'BTLA', 
                      'CD28', 'CD83', 'PRDX1', 'CCL5')

##################################################################
##                      Load the GE datasets                     #
##################################################################

#Create an empty list filled with two elements (one for each batch)
#each of these two elements contain a list with the length of sample per batch
batch_list <- list(vector(mode = 'list', length = length(batch1_ids)),
                   vector(mode = 'list', length = length(batch2_ids)))

GE_data_loader <- function(batch_dir, batch_ids, batch_num) {
  
  disease_groups <- c("Healthy", "HIV", "aL_HIV", "VL_HIV", "VL_HIV")
  
  for (i in seq_along(batch_ids)) {
    
    d10x <- Read10X(paste0(batch_dir,batch_ids[i], "/outs/per_sample_outs/", batch_ids[i], "/count/sample_feature_bc_matrix"))
    colnames(d10x) <- paste(batch_ids[i], colnames(d10x), sep = '_')
    colnames(d10x) <- stringr::str_remove(colnames(d10x), "-1")
    
    batch_list[[batch_num]][[i]] <- CreateSeuratObject(counts = d10x, min.cells = 3, min.features = 200)
    
    batch_num_to_add <- data.frame(batch_id = rep(batch_num, nrow(batch_list[[batch_num]][[i]]@meta.data)),
                                   disease_group = rep(disease_groups[i], nrow(batch_list[[batch_num]][[i]]@meta.data)),
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

# batch_list <- lapply(X = batch_list, FUN = function(x) {
#   x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
#   x <- subset(x, subset = nFeature_RNA >= 400 & nFeature_RNA <= 3000 & percent.mt < 10)
#   x <- NormalizeData(x,normalization.method = "LogNormalize", scale.factor = 10000)
#   x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000) #according to paper COVID19
# })


csect_combined <- merge(batch_list[[1]], batch_list[c(2:10)])
rm(batch_list)

csect_combined[["percent.mt"]] <- PercentageFeatureSet(csect_combined, pattern = "^MT-")
csect_combined[["percent.rb"]] <- PercentageFeatureSet(csect_combined, pattern = "^RP[SL]")
csect_combined <- subset(csect_combined,
                          subset = nFeature_RNA >= 800 & nFeature_RNA <= 3000 & percent.mt > 2.5 & percent.mt < 10)

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

combined <- combineTCR(contig_list, samples = c(c(batch1_ids, batch2_ids)), ID = c(c(batch1_ids, batch2_ids)), cells = "T-AB")

for (i in seq_along(combined)) {
  combined[[i]][,"barcode"] <- paste0(stringr::str_split(combined[[i]][,"barcode"], "_", simplify = T)[,1],
                                      '_',
                                      stringr::str_split(combined[[i]][,"barcode"], "_", simplify = T)[,3])
}

names(combined) <- stringr::str_split(names(combined), "_", simplify = T)[,1]

##################################################################
##                Add VDJ data to Seurat object                  #
##################################################################

csect_combined <- combineExpression(combined, csect_combined, proportion = FALSE,
                                     cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))

#################################################################
##            Calculating UMAP and finding clusters             #
#################################################################

csect_combined <- NormalizeData(csect_combined,normalization.method = "LogNormalize")
csect_combined <- FindVariableFeatures(csect_combined, selection.method = "vst", nfeatures = 2000) #according to paper COVID19
csect_combined <- ScaleData(object = csect_combined, verbose = FALSE)
csect_combined <- RunPCA(object = csect_combined, npcs = 40, verbose = FALSE)
#ElbowPlot(csect_combined)
csect_combined <- RunUMAP(object = csect_combined, reduction = "pca", dims = 1:15)
csect_combined <- FindNeighbors(object = csect_combined, dims = 1:15, force.recalc = T)
csect_combined <- FindClusters(object = csect_combined, resolution = 0.8)

#################################################################
##              Save object for Azimuth predictions             #
#################################################################

###Saving object to run Azimuth pipeline (https://satijalab.org/azimuth/)
#saveRDS(csect_combined, file = "immune_combined_crosssectional.rds")

##################################################################
##  Add Azimuth celltype predictions to Seurat object metadata   #
##################################################################

azimuth_predictions <- read.delim(here::here('data', 'processed_data', 'azimuth_pred_crosssectional.tsv'), row.names = 1)
csect_combined <- AddMetaData(object = csect_combined, metadata = azimuth_predictions)

#################################################################
##              SingleR for celltype identification             #
#################################################################

# hpca_se <- celldex::HumanPrimaryCellAtlasData()
# singler_pred_immune <- SingleR(test = GetAssayData(csect_combined), ref = hpca_se, labels = hpca_se$label.fine)
# csect_combined[["SingleR.labels"]] <- singler_pred_immune$labels
# rm(hpca_se, singler_pred_immune)

##################################################################
##                    Identify Cluster Markers                   #
##################################################################
#Could be used for manual cluster celltype identification

# find markers for every cluster compared to all remaining cells, report only the positive ones
# pbmc.markers <- FindAllMarkers(csect_combined, only.pos = TRUE, logfc.threshold = 0.25) %>% mutate(diff = pct.1 - pct.2)
# top30_FC_markers <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
# top30_diff_markers <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 30, wt = diff)
#clus_marker_heatmap <- DoHeatmap(csect_combined, features = top30_markers$gene) + NoLegend()

#################################################################
##                      Celltype labeling                       #
#################################################################

#Some plasma cells have been clustered with proliferating CD8 + NK cells
#Even though on the UMAP they form their own mini-cluster next to B cells
#So manually relabel them as a new cluster 22 and then rename them all
plasma <- c("HEC1_AAAGCAAAGATACACA", "HEC1_ACGAGCCCATCGGAAG", "HEC1_AGGCCACTCCCTGACT",
            "0001M3_CGCGTTTCAAGTACCT", "0123M3_GTGCATAGTCGACTAT", "0123M3_TTGGAACCATGGTCAT",
            "0104UV_AAATGCCCACGTCTCT", "0104UV_AACTCAGGTGGACGAT", "0104UV_AAGGAGCGTTGTGGAG",
            "0104UV_ACGCCGAAGGCCCGTT", "0104UV_ACTTTCACAAGCGTAG", "0104UV_AGAGCGAAGCAGATCG",
            "0104UV_AGAGTGGAGCCAGAAC", "0104UV_CAACCAAAGCTAGGCA", "0104UV_CACACAAAGACACGAC",
            "0104UV_CATGACAGTCGACTGC", "0104UV_CCTTCCCAGGCTAGCA", "0104UV_CGTTGGGGTTGGTTTG",
            "0104UV_CTCAGAAGTCAGAGGT", "0104UV_CTGAAACTCCCACTTG", "0104UV_CTTAACTCATCCGGGT",
            "0104UV_TCTCATACATGACGGA", "HEC2_AGTCTTTCAGTATCTG", "0090M3_AACCGCGGTCCGACGT",
            "0090M3_GATCTAGAGTCCTCCT", "0117M6_TCCACACGTCCAGTGC", "0170M6_CGAATGTGTCGTTGTA",
            "0170M6_CGCCAAGCACAGACTT", "0170M6_CGTAGGCGTCGACTGC", "0170M6_CTTACCGAGAAGGACA")

Idents(csect_combined, cells = plasma) <- '22'

cd8_in_cd4_naive <- WhichCells(csect_combined, slot = 'counts', idents = 6, expression = (CD8A > 0 | CD8B > 0) & CD4 <= 0)
Idents(csect_combined, cells = cd8_in_cd4_naive) <- '8'

cd8_in_cd4_prolif <- WhichCells(csect_combined, slot = 'counts', idents = 7, expression = (CD8A > 0 | CD8B > 0) & CD4 <= 0)
Idents(csect_combined, cells = cd8_in_cd4_prolif) <- '0'

cd8_in_cd4_TEM <- WhichCells(csect_combined, slot = 'counts', idents = c(3), expression = (CD8A > 0 | CD8B > 0) & CD4 <= 0)
cd8_in_cd4_TCM <- WhichCells(csect_combined, slot = 'counts', idents = c(4), expression = (CD8A > 0 | CD8B > 0) & CD4 <= 0)
Idents(csect_combined, cells = cd8_in_cd4_TEM) <- '0'
Idents(csect_combined, cells = cd8_in_cd4_TCM) <- '0'

#Systematic gating to split gdTs & possible NKs from CD8's
non_VDJ_carrying <- subset(csect_combined, idents = c(9, 11))
non_VDJ_carrying <- subset(non_VDJ_carrying, cloneType %in% c(NA, 'NA'))
Idents(csect_combined, cells = WhichCells(non_VDJ_carrying)) <- '23'
#Fully CD3 negative without VDJ data should be NK, very little so shouldnt matter much
NK_in_CD8 <- subset(csect_combined, idents = 23)
NK_in_CD8 <- WhichCells(NK_in_CD8, slot = 'counts',
                        expression = CD3D <= 0 & CD3E <= 0 & CD3G <= 0)
Idents(csect_combined, cells = NK_in_CD8) <- '2'
# #Mostly CD3 negative without VDJ data should probably also be NK, not big problem if missclassified
NK_in_CD8 <- subset(csect_combined, idents = '23')
NK_in_CD8 <- WhichCells(NK_in_CD8, slot = 'counts',
                        expression = ((CD3D <= 0 & CD3E <= 0) | (CD3D <= 0 & CD3G <= 0) | (CD3E <= 0 & CD3G <= 0)))
Idents(csect_combined, cells = NK_in_CD8) <- '2'
#Fully CD3 positive and fully CD8 positive should be grouped in with either cluster 9 or 11, does not matter: both will be CD8+ TEM
cd8pos_cd3pos_cells_in_23 <- subset(csect_combined, idents = '23')
cd8pos_cd3pos_cells_in_23 <- WhichCells(cd8pos_cd3pos_cells_in_23, slot = 'counts',
                                        expression = (CD3D > 0 & CD3E > 0 & CD3G > 0) & (CD8A > 0 & CD8B > 0))
Idents(csect_combined, cells = cd8pos_cd3pos_cells_in_23) <- '11'
#So there are some mostly CD3 positive and fully CD8 positive cells, these should be grouped with the CD8 TEM's
cd8pos_cd3pos_cells_in_23 <- subset(csect_combined, idents = 23)
cd8pos_cd3pos_cells_in_23 <- WhichCells(cd8pos_cd3pos_cells_in_23, slot = 'counts',
                                        expression = ((CD3D > 0 & CD3E > 0) | (CD3D > 0 & CD3G > 0) | (CD3E > 0 & CD3G > 0)) & (CD8A > 0 & CD8B > 0))
Idents(csect_combined, cells = cd8pos_cd3pos_cells_in_23) <- '11'
#Fully CD3 positive without VDJ data and fully CD8 negative should be gdT
#The following steps are just to check whether they are truly gdT based on markers, they won't be split off from cluster 23
gdt_in_CD8 <- subset(csect_combined, idents = 23)
gdt_in_CD8 <- WhichCells(gdt_in_CD8, slot = 'counts',
                         expression = (CD3E > 0 & CD3D > 0 & CD3G > 0) & (CD8A <= 0 & CD8B <= 0))
Idents(csect_combined, cells = gdt_in_CD8) <- '23'
#Mostly CD3 positive without VDJ data and CD8A positive but CD8B negative could possibly be gdT
gdt_in_CD8 <- subset(csect_combined, idents = 23)
gdt_in_CD8 <- WhichCells(gdt_in_CD8, slot = 'counts',
                         expression = ((CD3D > 0 & CD3E > 0) | (CD3D > 0 & CD3G > 0) | (CD3E > 0 & CD3G > 0)) & (CD8A >= 0 & CD8B <= 0))
Idents(csect_combined, cells = gdt_in_CD8) <- '23'


#Gating to isolate CD8's from NK cells. gdT's will be too hard / not possible
cd8_in_NK <- subset(csect_combined, idents = 2)
cd8_in_NK <- subset(cd8_in_NK, cloneType %in% c(NA, 'NA'), invert = T)
Idents(csect_combined, cells = WhichCells(cd8_in_NK)) <- '24'

#Setting up new celltype labels
new.cluster.ids <- c("0" = "CD8+ TCM",
                     "1" = "CD8+ TEM",
                     "2" = "CD56dim NK",
                     "3" = "CD4+ TCM",
                     "4" = "CD4+ TCM",
                     "5" = "CD14+ Monocytes",
                     "6" = "CD4+ T Naive",
                     "7" = "CD4+ TEM",
                     "8" = "CD8+ T Naive",
                     "9" = "CD8+ TEM",
                     "10" = "B Intermediate",
                     "11" = "CD8+ TEM",
                     "12" = "CD16+ Monocytes",
                     "13" = "B Naive",
                     "14" = "CD14+ Monocytes",
                     "15" = "CD56bright NK",
                     "16" = "CD14+ Monocytes",
                     "17" = "Platelets",
                     "18" = "Proliferating CD8+ T & NK",
                     "19" = "pDC",
                     "20" = "Erythrocytes",
                     "21" = "HSPCs",
                     "22" = "Plasmablasts",
                     "23" = "gdT",
                     "24" = "CD8+ TEM")

csect_combined <- RenameIdents(csect_combined, new.cluster.ids)
csect_combined[["stashed.ident"]] <- Idents(csect_combined)


scGate_modelsDB <- get_scGateDB(here("data", "meta_data", "scGateDB"), force_update = F)

hspc_model <- gating_model(name = 'HSPC', signature = c('SPINK2', 'SOX4', 'CD34'))
platelet_model <- gating_model(name = 'Platelets', signature = c('PPBP', 'PF4', 'TUBB1'))
cd4_prolif_model <- gating_model(model = scGate_modelsDB$human$generic$CD4T, name = 'CD4Prolif', level = 5, signature = c('MKI67', 'STMN1'))
cd8_prolif_model <- gating_model(model = scGate_modelsDB$human$generic$CD8T, name = 'CD8Prolif', level = 5, signature = c('MKI67', 'STMN1'))
NK_prolif_model <- gating_model(model = scGate_modelsDB$human$generic$NK, name = 'Prolif', level = 4, signature = c('MKI67', 'STMN1'))

csect_combined <- scGate(csect_combined, NK_prolif_model)
nk_prolif_cells <- WhichCells(subset(csect_combined, subset = is.pure == 'Pure'))
csect_combined <- scGate(csect_combined, cd8_prolif_model)
cd8_prolif_cells <- WhichCells(subset(csect_combined, subset = is.pure == 'Pure'))
csect_combined <- scGate(csect_combined, cd4_prolif_model)
cd4_prolif_cells <- WhichCells(subset(csect_combined, subset = is.pure == 'Pure'))

csect_combined <- scGate(csect_combined, hspc_model)
hspcs <- WhichCells(subset(csect_combined, subset = is.pure == 'Pure'))
csect_combined <- scGate(csect_combined, platelet_model)
platelets <- WhichCells(subset(csect_combined, subset = is.pure == 'Pure'))

Idents(csect_combined) <- "stashed.ident"
Idents(csect_combined, cells = hspcs) <- "HSPCs"
Idents(csect_combined, cells = platelets) <- "Platelets"
Idents(csect_combined, cells = cd4_prolif_cells) <- "CD4+ T Proliferating"
Idents(csect_combined, cells = cd8_prolif_cells) <- "CD8+ T Proliferating"
Idents(csect_combined, cells = nk_prolif_cells) <- "NK Proliferating"
Idents(csect_combined, cells = WhichCells(subset(csect_combined, idents = 'Proliferating CD8+ T & NK', subset = cloneType %in% c(NA, 'NA')))) <- 'NK Proliferating'
Idents(csect_combined, cells = WhichCells(subset(csect_combined, idents = 'Proliferating CD8+ T & NK'))) <- 'CD8+ T Proliferating'

Idents(csect_combined) <- factor(Idents(csect_combined), levels = c("CD4+ T Proliferating", "CD4+ T Naive", "CD4+ TCM", "CD4+ TEM",
                                                                  "CD8+ T Proliferating", "CD8+ T Naive", "CD8+ TCM", "CD8+ TEM",
                                                                  'gdT', "NK Proliferating", "CD56dim NK", "CD56bright NK",
                                                                  "CD14+ Monocytes", "CD16+ Monocytes",
                                                                  "B Naive", "B Intermediate", "Plasmablasts",
                                                                  "Platelets", "Erythrocytes", "HSPCs", "pDC"))

#################################################################
##              Save UMAP figures for publication               #
#################################################################

csect_combined$disease_group <- as.factor(csect_combined$disease_group)
csect_combined$disease_group <- factor(csect_combined$disease_group, levels = c('Healthy', 'HIV', 'aL_HIV', 'VL_HIV'))
csect_combined$cloneType <- as.factor(csect_combined$cloneType)
csect_combined$cloneType <- factor(csect_combined$cloneType, levels = c("Large (20 < X <= 100)", "Medium (5 < X <= 20)", "Small (1 < X <= 5)", "Single (0 < X <= 1)", NA))

##---------------------------------------------------------------
##                  Cellular composition graphs                 -
##---------------------------------------------------------------

##----------------------------------------------------------------
##                  UMAPs and frequency tables                   -
##----------------------------------------------------------------

celltype_colors <- c("CD4+ T Proliferating" = "#FF9D9A", "CD4+ T Naive" = "#D37295", "CD4+ TCM" ="#FABFD2", "CD4+ TEM" = "#E15759",
                     "CD8+ T Proliferating" = "#499894", "CD8+ T Naive" = "#86BCB6", "CD8+ TCM" = "#A0CBE8", "CD8+ TEM" = "#4E79A7",
                     'gdT' = "#e377c2", "NK Proliferating" = "#F1CE63", "CD56dim NK" = "#FFBE7D", "CD56bright NK" = "#F28E2B",
                     "CD14+ Monocytes" = "#98df8a", "CD16+ Monocytes" = "#2ca02c",
                     "B Naive" = "#D7B5A6", "B Intermediate" = "#B6992D", "Plasmablasts" = "#8c564b",
                     "Platelets" = "#D4A6C8", "Erythrocytes" ="#17becf", "HSPCs" = "#B07AA1", "pDC" = "#bcbd22")

cross_sect_umap_axes <- DimPlot(csect_combined, label = T, repel = T, label.size = 8, pt.size = 1.2,
                                cols = celltype_colors) +
  ggtitle("UMAP representation of immune cells of all patient groups") +
  theme_void() + theme(legend.position = 'none',
                     panel.grid = element_blank(),
                     plot.title = element_text(hjust = 0.5, family = text_font, size = 36, face = 'bold')) + 
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 7.5)))

# freq_table <- tableGrob(rbind(table(csect_combined@active.ident, csect_combined$disease_group),`Total Cells` = colSums(table(csect_combined@active.ident, csect_combined$disease_group))),
#                         theme = ttheme('classic', base_size = 15))

patient_table <- as_tibble(rbind(table(csect_combined@active.ident,
                                       csect_combined$orig.ident)),
                           rownames = 'Celltype') %>%
  pivot_longer(!Celltype, names_to = 'PatientID', values_to = 'Count_per_patient') %>%
  inner_join(., csect_combined[[]], by = c('PatientID' = 'orig.ident')) %>%
  select(c(Celltype, PatientID, disease_group, Count_per_patient)) %>% 
  distinct()

total_per_patient_table <- patient_table %>% 
  group_by(PatientID) %>% 
  summarise(Total = sum(Count_per_patient))

patient_table <- inner_join(patient_table, total_per_patient_table, by = c("PatientID")) %>%
  mutate(Percentage = (Count_per_patient / Total) * 100) %>%
  group_by(Celltype, disease_group) %>% mutate(Mean_Percentage = mean(Percentage)) %>% ungroup()

mean_table <- patient_table %>% select(., c(Celltype, disease_group, Mean_Percentage)) %>% distinct()

temp_csect_combined <- csect_combined
Idents(temp_csect_combined) <- 'disease_group'

total_cellcounts_per_group <- as_tibble(rbind(table(Idents(temp_csect_combined)))) %>% pivot_longer(everything(), names_to = 'disease_group', values_to = 'Cells')
rm(temp_csect_combined)

stacked_barchart_cross_sectional <- ggplot(mean_table,
  aes(factor(disease_group, levels = c("Healthy", "HIV", "aL_HIV", "VL_HIV")),
  Mean_Percentage,
  fill = factor(Celltype, levels = c("CD4+ T Proliferating", "CD4+ T Naive", "CD4+ TCM", "CD4+ TEM",
                                                "CD8+ T Proliferating", "CD8+ T Naive", "CD8+ TCM", "CD8+ TEM",
                                                'gdT', "NK Proliferating", "CD56dim NK", "CD56bright NK",
                                                "CD14+ Monocytes", "CD16+ Monocytes",
                                                "B Naive", "B Intermediate", "Plasmablasts",
                                                "Platelets", "Erythrocytes", "HSPCs", "pDC")))) +
  geom_bar(stat = "identity") +
  labs(y = 'Percentage') +
  scale_x_discrete(expand = c(0.08, 0.23),
                   labels = paste0(levels(mean_table$disease_group), ' (n = ', total_cellcounts_per_group$Cells, ')')) +
  scale_y_continuous(expand = c(0.00, 0.01)) +
  ggtitle("Cellular composition of all patient groups") +
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
                                                                                                    "Platelets", "Erythrocytes", "HSPCs", "pDC"))))

#-----------Saving cross_sectional + longitudinal at the same time
layout <- "
AAAABB
CCCCDD
"

arranged_plot <- cross_sect_umap_axes + stacked_barchart_cross_sectional +
  longitudinal_umap_axes + stacked_barchart_longitudinal +
  #celltype_barchart_cross_sectional + celltype_barchart_longitudinal +
  plot_layout(design = layout) & plot_annotation(tag_levels = c('a')) & theme(plot.tag = element_text(size = 36, face = 'bold'))

# ggsave(here("analyses", "final_output", 'SingleCell_Figure3_UMAP_Composition.pdf'), arranged_plot, height = 17.5, width = 30, device = cairo_pdf)
# ggsave(here("analyses", "final_output", 'SingleCell_Figure3_UMAP_Composition.png'), arranged_plot, height = 17.5, width = 30, type = 'cairo-png')

##################################################################
##                Gene Set Enrichments All Cells                 #
##################################################################

csect_combined2 <- csect_combined
Idents(csect_combined2) <- "disease_group"

disease_groups_to_test <- c('HIV' = 'HIV', 'aL_HIV' = 'aL_HIV', 'VL_HIV' = 'VL_HIV')

allcells_vs_healthy_list <- lapply(disease_groups_to_test, function(x) {
  FindMarkers(csect_combined2, ident.1 = x, ident.2 = "Healthy", min.pct = 0.1, logfc.threshold = 0.25, only.pos = F) %>% mutate(gene = rownames(.), diff = pct.1 - pct.2)
})

AsympLeish_vs_SympLeish <- FindMarkers(csect_combined2, ident.1 = "aL_HIV", ident.2 = "VL_HIV", min.pct = 0.1, logfc.threshold = 0.25, only.pos = F) %>% mutate(gene = rownames(.), diff = pct.1 - pct.2)
NoLeish_vs_Leish <- FindMarkers(csect_combined2, ident.1 = c('Healthy', 'HIV'), ident.2 = c('aL_HIV', "VL_HIV"), min.pct = 0.1, logfc.threshold = 0.25, only.pos = F) %>% mutate(gene = rownames(.), diff = pct.1 - pct.2)

#################################################################
##                  Gene Set Enrichments CD4+ T                 #
#################################################################

allcd4 <- names(csect_combined@active.ident[grepl('CD4', csect_combined@active.ident)])
cd4_subset <- subset(csect_combined, cells = allcd4)

Idents(cd4_subset) <- "disease_group"

cd4_vs_healthy_list <- lapply(disease_groups_to_test, function(x) {
  FindMarkers(cd4_subset, ident.1 = x, ident.2 = "Healthy", min.pct = 0.1, logfc.threshold = 0.25, only.pos = F) %>% mutate(gene = rownames(.), diff = pct.1 - pct.2)
})

CD4_AsympLeish_vs_SympLeish <- FindMarkers(cd4_subset, ident.1 = "aL_HIV", ident.2 = "VL_HIV", min.pct = 0.1, logfc.threshold = 0.25, only.pos = F) %>% mutate(gene = rownames(.), diff = pct.1 - pct.2)
CD4_NoLeish_vs_Leish <- FindMarkers(cd4_subset, ident.1 = c('Healthy', 'HIV'), ident.2 = c('aL_HIV', "VL_HIV"), min.pct = 0.1, logfc.threshold = 0.25, only.pos = F) %>% mutate(gene = rownames(.), diff = pct.1 - pct.2)

#################################################################
##                  Gene Set Enrichments CD8+ T                 #
#################################################################

allcd8 <- names(csect_combined@active.ident[grepl('CD8', csect_combined@active.ident)])
cd8_subset <- subset(csect_combined, cells = allcd8)

Idents(cd8_subset) <- "disease_group"

cd8_vs_healthy_list <- lapply(disease_groups_to_test, function(x) {
  FindMarkers(cd8_subset, ident.1 = x, ident.2 = "Healthy", min.pct = 0.1, logfc.threshold = 0.25, only.pos = F) %>% mutate(gene = rownames(.), diff = pct.1 - pct.2)
})

CD8_AsympLeish_vs_SympLeish <- FindMarkers(cd8_subset, ident.1 = "aL_HIV", ident.2 = "VL_HIV", min.pct = 0.1, logfc.threshold = 0.25, only.pos = F) %>% mutate(gene = rownames(.), diff = pct.1 - pct.2)
CD8_NoLeish_vs_Leish <- FindMarkers(cd8_subset, ident.1 = c('Healthy', 'HIV'), ident.2 = c('aL_HIV', "VL_HIV"), min.pct = 0.1, logfc.threshold = 0.25, only.pos = F) %>% mutate(gene = rownames(.), diff = pct.1 - pct.2)

##################################################################
##                    Gene Set Enrichments APC                   #
##################################################################

allAPC <- names(csect_combined@active.ident[grepl('Monocytes', csect_combined@active.ident) | grepl('pDC', csect_combined@active.ident)])
apc_subset <- subset(csect_combined, cells = allAPC)

Idents(apc_subset) <- "disease_group"

APC_vs_healthy_list <- lapply(disease_groups_to_test, function(x) {
  FindMarkers(apc_subset, ident.1 = x, ident.2 = "Healthy", min.pct = 0.1, logfc.threshold = 0.25, only.pos = F) %>% mutate(gene = rownames(.), diff = pct.1 - pct.2)
})

APC_AsympLeish_vs_SympLeish <- FindMarkers(apc_subset, ident.1 = "aL_HIV", ident.2 = "VL_HIV", min.pct = 0.1, logfc.threshold = 0.25, only.pos = F) %>% mutate(gene = rownames(.), diff = pct.1 - pct.2)
APC_NoLeish_vs_Leish <- FindMarkers(apc_subset, ident.1 = c('Healthy', 'HIV'), ident.2 = c('aL_HIV', "VL_HIV"), min.pct = 0.1, logfc.threshold = 0.25, only.pos = F) %>% mutate(gene = rownames(.), diff = pct.1 - pct.2)

#################################################################
##                    Gene Set Enrichments NK                   #
#################################################################

allNK <- names(csect_combined@active.ident[grepl('NK', csect_combined@active.ident)])
nk_subset <- subset(csect_combined, cells = allNK)

Idents(nk_subset) <- "disease_group"

NK_vs_healthy_list <- lapply(disease_groups_to_test, function(x) {
  FindMarkers(nk_subset, ident.1 = x, ident.2 = "Healthy", min.pct = 0.1, logfc.threshold = 0.25, only.pos = F) %>% mutate(gene = rownames(.), diff = pct.1 - pct.2)
})

NK_AsympLeish_vs_SympLeish <- FindMarkers(nk_subset, ident.1 = "aL_HIV", ident.2 = "VL_HIV", min.pct = 0.1, logfc.threshold = 0.25, only.pos = F) %>% mutate(gene = rownames(.), diff = pct.1 - pct.2)
NK_NoLeish_vs_Leish <- FindMarkers(nk_subset, ident.1 = c('Healthy', 'HIV'), ident.2 = c('aL_HIV', "VL_HIV"), min.pct = 0.1, logfc.threshold = 0.25, only.pos = F) %>% mutate(gene = rownames(.), diff = pct.1 - pct.2)

##################################################################
##        Generating DE VolcanoPlots & Enrichment Networks       #
##################################################################

generate_vulcano_plots <- function(marker_list, num_genes_to_display) {
  
  vulcano_plot_list <- lapply(seq_along(marker_list), function(x) {
    tmp <- marker_list[[x]]
    
    tmp <- tmp %>%
      mutate(Trend = ifelse(p_val_adj < 0.05 & avg_log2FC > 0, "Up",
                            ifelse(p_val_adj < 0.05 & avg_log2FC < 0, "Down", "None")))
    filter <- subset(tmp, p_val_adj < 0.05)
    top_n_genes <- filter %>% top_n(n = num_genes_to_display, wt = avg_log2FC + 4*diff)
    bottom_n_genes <- filter %>% top_n(n = -num_genes_to_display, wt = avg_log2FC + 4*diff)
    
    tmp <- tmp %>%
      mutate(p_val_adj = ifelse(p_val_adj == 0, min(p_val_adj), p_val_adj)) #mutate the p-values that are lower than detection
    
    ggplot(tmp, aes(x=avg_log2FC, y=-log10(p_val_adj))) + 
      geom_point(size=0.5, color="#999999") + 
      geom_point(data=subset(tmp, Trend == "Up" | Trend == "Down"), aes(color = Trend)) + 
      ggtitle(names(marker_list)[x]) +
      labs(x = "Log2(Fold change)", y = "-Log10(P value adjusted)") +
      theme_classic(base_rect_size = 0) + 
      theme(plot.title = element_text(hjust = 0.5, family = 'Roboto Condensed', size = 20),
            text = element_text(family = 'Roboto Condensed'),
            axis.title = element_text(family = 'Roboto Condensed', size = 14),
            axis.text = element_text(family = 'Roboto Condensed'),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            panel.grid.major.x = element_line(color="#dadfe6", size=0.15),
            panel.grid.major.y = element_line(color="#dadfe6", size=0.15)) +
      scale_y_sqrt() +
      geom_vline(xintercept = 0, lty = 2) + 
      geom_hline(yintercept = 1.3, lty = 2) + 
      geom_text_repel(data=subset(tmp, gene %in% top_n_genes$gene), aes(label=gene, family = 'Roboto Condensed'), segment.size = 0.25, size=4) + 
      geom_text_repel(data=subset(tmp, gene %in% bottom_n_genes$gene & !gene %in% top_n_genes$gene), aes(label=gene, family = 'Roboto Condensed'), segment.size = 0.25, size=4) + 
      scale_color_manual(values = c("#0348A6", "#FF4B20")) + 
      guides(color = "none")
    
  })
  
  names(vulcano_plot_list) <- names(marker_list)
  return(vulcano_plot_list)
}

enrichment_plotter <- function(marker_list) {
  
  m_t2g <- msigdbr::msigdbr(species = "Homo sapiens", category = c('H')) %>% 
    dplyr::select(gs_name, entrez_gene) %>% dplyr::distinct(gs_name, entrez_gene)
  
  go_m_t2g <- msigdbr::msigdbr(species = "Homo sapiens", category = c('C5')) %>% filter(., gs_subcat != 'HPO') %>%
    dplyr::select(gs_name, entrez_gene) %>% dplyr::distinct(gs_name, entrez_gene)
  
  m_t2g <- bind_rows(m_t2g, go_m_t2g)
  
  go_jointable <- msigdbr::msigdbr(species = "Homo sapiens", category = c('C5')) %>% filter(., gs_subcat != 'HPO')
  
  enrich_result <- lapply(marker_list, FUN = function(x) {
    x <- filter(x, p_val_adj < 0.05)
    
    gene_list <- row.names(x)
    marker_expression <- x$avg_log2FC
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
      go_result_mf <- filter(go_result, ONTOLOGY == 'MF')
      go_result_cc <- filter(go_result, ONTOLOGY == 'CC')
      
      enriched@ontology <- 'BP'
      enriched@keytype <- "ENTREZID"
      enriched@organism <- "Homo sapiens"
      
      enriched@result <- go_result_bp
      enriched <- simplify(enriched)
      simplified_bp <- enriched@result
      
      enriched@result <- go_result_mf
      enriched <- simplify(enriched)
      simplified_mf <- enriched@result
      
      enriched@result <- go_result_cc
      enriched <- simplify(enriched)
      simplified_cc <- enriched@result
      
      total_simplified <- bind_rows(simplified_bp, simplified_mf, simplified_cc) %>% select(!ONTOLOGY) %>% mutate(ID = Description)
      rownames(total_simplified) <- total_simplified$ID
      total_simplified <- bind_rows(total_simplified, hallmark_result)
      
      enriched@result <- total_simplified
      
      enriched@result <- enriched@result[order(enriched@result$p.adjust),]
      
      enriched <- setReadable(enriched, org.Hs.eg.db, keyType = "ENTREZID")
      
      p <- cnetplot(enriched, showCategory = 10, foldChange=enrich_ready, circular = F, colorEdge = TRUE,
                    layout = 'kk',
                    node_label = 'gene',
                    cex_category = 1, cex_gene = 0.8, cex_label_category = 1, cex_label_gene = 0.8) + 
        guides(edge_color = 'none', size_new = 'none') +
        theme(legend.title = element_text(size=14, face='bold'),
              legend.text = element_text(size=11),
              text = element_text(family = 'Roboto Condensed')) + ggraph::set_graph_style(family = 'Roboto Condensed')
      
      p <- p + geom_node_text(aes(label = p$data[1:10,3], family = 'Roboto Condensed Bold'), data = p$data[1:10,], size = 5, repel = T)
      
      return(p)
    } else {
      return(NULL)
    }
    
  })
  
  return(enrich_result)
  
}

##################################################################
##                          TCR Analysis                         #
##################################################################

cross_sect_clonotype_axes <- DimPlot(csect_combined, group.by = 'cloneType') +
  scale_colour_manual(values = colorblind_vector(5), na.value = 'grey') + 
  theme_void() + theme(plot.title = element_blank(), legend.text = element_text(size=17), text = element_text(family = text_font))

csect_tcr_by_cluster <- expression2List(csect_combined, group = 'cluster')
csect_tcr_by_cluster <- csect_tcr_by_cluster[(grepl("T", names(csect_tcr_by_cluster)) & names(csect_tcr_by_cluster) != c('gdT'))]
csect_tcr_by_individual <- expression2List(csect_combined, group = 'orig.ident')
csect_tcr_by_disease_group <- expression2List(csect_combined, group = 'disease_group')

##----------------------------------------------------------------
##            Writing TCR data to TCRex-compliant file           -
##----------------------------------------------------------------

csect_tcrex_compliant <- data.frame(CDR3_beta = na_if(str_split(unname(csect_combined$CTaa), '_', simplify = T)[,2], ""),
                              TRBJ_gene = str_split(unname(csect_combined$CTgene), "_", simplify = T)[,2],
                              TRBV_gene = str_split(unname(csect_combined$CTgene), "_", simplify = T)[,2],
                              Barcode = names(csect_combined@active.ident),
                              Disease_group = csect_combined$disease_group,
                              Celltype = unname(csect_combined@active.ident),
                              Frequency = csect_combined$Frequency) %>%
  mutate(CDR3_beta = na_if(CDR3_beta, 'NA'),
         TRBJ_gene = na_if(TRBJ_gene, 'NA'),
         TRBV_gene = na_if(TRBV_gene, 'NA')) %>% 
  tidyr::separate_rows(., c(1:3), sep = ';') %>%
  na.omit() %>%
  mutate(TRBJ_gene = unlist(str_extract_all(TRBJ_gene, '(TRBJ\\d+-?\\d*)')),
         TRBV_gene = unlist(str_extract_all(TRBV_gene, '(TRBV\\d+-?\\d*)')))

# write.table(csect_tcrex_compliant[,1:3], file = 'cross_sect_tcrex_compliant.tsv', sep = '\t', quote = F, row.names = F)
# write.table(csect_tcrex_compliant, file = 'cross_sect_tcrex_metadata.tsv', sep = '\t', quote = F, row.names = F)

##----------------------------------------------------------------
##                      Diversity analysis                       -
##----------------------------------------------------------------

csect_diversity_by_cluster <- clonalDiversity(csect_tcr_by_cluster, cloneCall = 'gene+nt', chain = 'both', n.boots = 100, group = 'cluster')
csect_diversity_by_individual_data <- clonalDiversity(csect_tcr_by_individual, cloneCall = 'gene+nt', chain = 'both', n.boots = 100, group = 'orig.ident', exportTable = T)
csect_diversity_by_individual_data$x.axis <- NULL
csect_diversity_by_individual_data <- suppressWarnings(melt(csect_diversity_by_individual_data, id.vars = 'orig.ident'))
csect_diversity_by_individual_data <- csect_diversity_by_individual_data %>% 
  mutate(disease_group = case_when(
    orig.ident %in% levels(as.factor(csect_tcr_by_disease_group$cVL_HIV$orig.ident)) ~ 'cVL_HIV',
    orig.ident %in% levels(as.factor(csect_tcr_by_disease_group$pVL_HIV$orig.ident)) ~ 'pVL_HIV',
    orig.ident %in% levels(as.factor(csect_tcr_by_disease_group$aL_HIV$orig.ident)) ~ 'aL_HIV',
    orig.ident %in% levels(as.factor(csect_tcr_by_disease_group$HIV$orig.ident)) ~ 'HIV',
    orig.ident %in% levels(as.factor(csect_tcr_by_disease_group$Healthy$orig.ident)) ~ 'Healthy'
  ))

csect_diversity_plot <- ggplot(csect_diversity_by_individual_data,
                               aes(y = as.numeric(value), 
                                   x = factor(disease_group, levels = c('Healthy', 'HIV', 'aL_HIV', 'pVL_HIV', 'cVL_HIV')),
                                   fill = factor(disease_group, levels = c('Healthy', 'HIV', 'aL_HIV', 'pVL_HIV', 'cVL_HIV')))) +
  geom_boxplot() + facet_wrap(~variable, scales = 'free', ncol = 4) +
  geom_jitter(aes(group = disease_group)) +
  #geom_jitter(position = position_dodge(width = 0.75), aes(group = disease_group)) +
  #geom_text_repel(position = position_dodge(width = 0.75), aes(label = orig.ident), direction = 'y', size = 4.5) +
  labs(x = 'Disease groups', y = 'Diversity metric') +
  scale_fill_npg() + scale_color_npg() + theme_bw() +
  theme(axis.text.y = element_text(size = 16, family = text_font, color = 'black'),
        axis.title.y = element_text(size = 18, family = text_font),
        axis.text.x = element_text(color='black', size=16, angle = 90, vjust = 0.5, hjust = 1, family=text_font),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 18, family = text_font),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color="#dadfe6", size=0.25),
        strip.background = element_blank(),
        strip.text = element_text(size = 20, family = text_font))


# diversity_by_disease_group <- clonalDiversity(csect_tcr_by_disease_group, cloneCall = 'gene+nt', chain = 'both', n.boots = 100, group = 'disease_group')

##----------------------------------------------------------------
##                        Repertoire space                       -
##----------------------------------------------------------------

csect_combined$dg_patient <- factor(paste0(csect_combined$disease_group, '_', csect_combined$orig.ident),
                                    levels = c("Healthy_HEC1", "Healthy_HEC2", "HIV_0001M3", "HIV_0090M3", "aL_HIV_0120M3", "aL_HIV_0123M3",
                                               "VL_HIV_0114UV", "VL_HIV_0117M6", "VL_HIV_0104UV", "VL_HIV_0170M6"))

csect_clonal_occupation <- occupiedscRepertoire(csect_combined, x.axis = 'dg_patient', proportion = T) + 
  labs(y = 'Proportion') +
  scale_fill_manual(values = colorblind_vector(5)[2:5]) +
  theme(axis.text.y = element_text(size = 18, family = text_font, color = 'black'),
        axis.title.y = element_text(size = 20, family = text_font),
        axis.text.x = element_text(color='black', size=18, angle = 90, vjust = 0.5, family=text_font),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 18, family = text_font))

##----------------------------------------------------------------
##                          Gene usage                           -
##----------------------------------------------------------------

csect_TRA_V_gene_usage <- vizGenes(csect_combined, gene = 'V', chain = 'TRA', plot = 'bar', separate = 'disease_group', order = 'variance', errorbar = T)
csect_TRA_J_gene_usage <- vizGenes(csect_combined, gene = 'J', chain = 'TRA', plot = 'bar', separate = 'disease_group', order = 'variance', errorbar = T)
csect_TRB_V_gene_usage <- vizGenes(csect_combined, gene = 'V', chain = 'TRB', plot = 'bar', separate = 'disease_group', order = 'variance', errorbar = T)
csect_TRB_J_gene_usage <- vizGenes(csect_combined, gene = 'J', chain = 'TRB', plot = 'bar', separate = 'disease_group', order = 'variance', errorbar = T)

##---------------------------------------------------------------
##                    TCR Plots Publication                     -
##---------------------------------------------------------------

# layout <- "
# AAABBB
# CCCDDD
# EEEFFF
# GGGGGG
# HHHHHH
# "
#
# layout2 <- "
# AAABBB
# CCCDDD
# "
#
# csect_arranged <- csect_TRA_V_gene_usage + csect_TRB_V_gene_usage +
#   csect_TRA_J_gene_usage + csect_TRB_J_gene_usage +
#   plot_layout(design = layout2, guides = 'collect')
#
# long_arranged <- long_TRA_V_gene_usage + long_TRB_V_gene_usage +
#   long_TRA_J_gene_usage + long_TRB_J_gene_usage +
#   plot_layout(design = layout2, guides = 'collect')
#
#
# arranged_plot <- cross_sect_clonotype_axes + csect_diversity_plot +
#   longitudinal_clonotype_axes + long_diversity_plot +
#   csect_clonal_occupation + long_clonal_occupation +
#   csect_arranged + long_arranged +
#     plot_layout(design = layout) & plot_annotation(tag_levels = c('a')) & theme(plot.tag = element_text(size = 20, face = 'bold'))
#
# ggsave('CombinedTCRs_Fig7.pdf', arranged_plot, height = 30, width = 25, device = cairo_pdf)
