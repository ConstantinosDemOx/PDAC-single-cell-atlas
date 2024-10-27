# Identify recurring ecosystems of pancreatic cancer

library(tidyverse)
library(DropletUtils)
library(SpotClean)
library(rcartocolor)
library(RColorBrewer)
library(Seurat)
library(corrplot)
library(cluster)
library(ConsensusClusterPlus)


library(Seurat);library(tidyverse);library(ggpubr);library(circlize);library(ComplexHeatmap);library(factoextra);library(cluster);library(factoextra)
library(tidygraph);library(igraph);library(ggraph);library(corrr);library(viridis)

source("misc_functions.R")
source("color_functions.R")

pdac.atlas <- readRDS("pdac.atlas_final.rds")

epi.gps <- readRDS("Consensus method/Figure 2/Assets/Epithelial/epi_gps.rds")
epi.gps
epi.cells <- readRDS("Consensus method/Figure 2/Assets/Epithelial/cell.mat.split.rds")
epi.cells <- epi.cells[1:6]
names(epi.cells) <- paste0("Epi_", names(epi.cells))

pdac.atlas@meta.data <- pdac.atlas@meta.data %>% 
  mutate(GP = case_when(
    rownames(.) %in% epi.cells[[1]] ~ "Epi_KRT17",
    rownames(.) %in% epi.cells[[2]] ~ "Epi_TFF1",
    rownames(.) %in% epi.cells[[3]] ~ "Epi_HLA",
    rownames(.) %in% epi.cells[[4]] ~ "Epi_AMY2A",
    rownames(.) %in% epi.cells[[5]] ~ "Epi_FXYD2",
    rownames(.) %in% epi.cells[[6]] ~ "Epi_STMN1",
    T ~ Lvl2_Clusters))

table(pdac.atlas$Lvl1_Clusters, pdac.atlas$GP)


fib.gps <- readRDS("Consensus method/Figure 2/Assets/Fibroblast/fib_gps.rds")
fib.gps
fib.cells <- readRDS("Consensus method/Figure 2/Assets/Fibroblast/cell.mat.split.rds")

fib.cells <- fib.cells[1:8]
names(fib.cells) <- paste0("Fib_", names(fib.cells))


pdac.atlas@meta.data <- pdac.atlas@meta.data %>% 
  mutate(GP = case_when(
    rownames(.) %in% fib.cells[[1]] ~ "Fib_MMP11",
    rownames(.) %in% fib.cells[[2]] ~ "Fib_C7",
    rownames(.) %in% fib.cells[[3]] ~ "Fib_IGFBP5",
    rownames(.) %in% fib.cells[[4]] ~ "Fib_ENO1",
    rownames(.) %in% fib.cells[[5]] ~ "Fib_RGS5",
    rownames(.) %in% fib.cells[[6]] ~ "Fib_MYH11",
    rownames(.) %in% fib.cells[[7]] ~ "Fib_IGF1",
    rownames(.) %in% fib.cells[[8]] ~ "Fib_CD74",
    T ~ GP))

table(pdac.atlas$Lvl1_Clusters, pdac.atlas$GP)


ecs.gps <- readRDS("Consensus method/Figure 2/Assets/Endothelial/ecs_gps.rds")
ecs.gps
ecs.cells <- readRDS("Consensus method/Figure 2/Assets/Endothelial/cell.mat.split.rds")
ecs.cells <- ecs.cells[1:4]
names(ecs.cells) <- paste0("EC_", names(ecs.cells))


pdac.atlas@meta.data <- pdac.atlas@meta.data %>% 
  mutate(GP = case_when(
    rownames(.) %in% ecs.cells[[1]] ~ "EC_RGCC",
    rownames(.) %in% ecs.cells[[2]] ~ "EC_SEMA3G",
    rownames(.) %in% ecs.cells[[3]] ~ "EC_ACKR1",
    rownames(.) %in% ecs.cells[[4]] ~ "EC_F10",
    T ~ GP))

table(pdac.atlas$Lvl1_Clusters, pdac.atlas$GP)




cd4t.cells <- readRDS("Consensus method/Figure 2/Assets/CD4T/cell.mat.split.rds")
cd4t.cells <- cd4t.cells[1:7]
names(cd4t.cells) <- paste0("CD4T_", names(cd4t.cells))


pdac.atlas@meta.data <- pdac.atlas@meta.data %>% 
  mutate(GP = case_when(
    rownames(.) %in% cd4t.cells[[1]] ~ "CD4T_CCR7",
    rownames(.) %in% cd4t.cells[[2]] ~ "CD4T_ANXA1",
    rownames(.) %in% cd4t.cells[[3]] ~ "CD4T_FOXP3",
    rownames(.) %in% cd4t.cells[[4]] ~ "CD4T_GNLY",
    rownames(.) %in% cd4t.cells[[5]] ~ "CD4T_FOSB",
    rownames(.) %in% cd4t.cells[[6]] ~ "CD4T_STMN1",
    rownames(.) %in% cd4t.cells[[7]] ~ "CD4T_CXCL13",
    T ~ GP))

table(pdac.atlas$Lvl1_Clusters, pdac.atlas$GP)



cd8t.gps <- readRDS("Consensus method/Figure 2/Assets/CD8T/CD8T_gps.rds")
cd8t.gps
cd8t.cells <- readRDS("Consensus method/Figure 2/Assets/CD8T/cell.mat.split.rds")
cd8t.cells <- cd8t.cells[1:5]
names(cd8t.cells) <- paste0("CD8T_", names(cd8t.cells))


pdac.atlas@meta.data <- pdac.atlas@meta.data %>% 
  mutate(GP = case_when(
    rownames(.) %in% cd8t.cells[[1]] ~ "CD8T_GZMK",
    rownames(.) %in% cd8t.cells[[2]] ~"CD8T_IL7R",
    rownames(.) %in% cd8t.cells[[3]] ~ "CD8T_GZMH",
    rownames(.) %in% cd8t.cells[[4]] ~ "CD8T_KLRC1",
    rownames(.) %in% cd8t.cells[[5]] ~ "CD8T_CTLA4",
    T ~ GP))




nk.cells <- readRDS("Consensus method/Figure 2/Assets/NK/cell.mat.split.rds")
nk.cells <- nk.cells[1:3]
names(nk.cells) <- paste0("NK_", names(nk.cells))


pdac.atlas@meta.data <- pdac.atlas@meta.data %>% 
  mutate(GP = case_when(
    rownames(.) %in% nk.cells[[1]] ~ "NK_XCL1",
    rownames(.) %in% nk.cells[[2]] ~ "NK_FCGR3A",
    rownames(.) %in% nk.cells[[3]] ~ "NK_STMN1",
    T ~ GP))


b.gps <- readRDS("Consensus method/Figure 2/Assets/B/B.Pl_gps.rds")
b.gps
B.cells <- readRDS("Consensus method/Figure 2/Assets/B/cell.mat.split.rds")
B.cells <- B.cells[1:6]
names(B.cells) <- paste0("B_", names(B.cells))


pdac.atlas@meta.data <- pdac.atlas@meta.data %>% 
  mutate(GP = case_when(
    rownames(.) %in% B.cells[[1]] ~ "B_STMN1",
    rownames(.) %in% B.cells[[2]] ~ "B_CRIP1",
    rownames(.) %in% B.cells[[3]] ~ "B_TCL1A",
    rownames(.) %in% B.cells[[4]] ~ "B_CD83",
    rownames(.) %in% B.cells[[5]] ~ "B_TNF",
    rownames(.) %in% B.cells[[6]] ~ "B_ANXA1",
    T ~ GP))


MacMono_gps <- readRDS("Consensus method/Figure 2/Assets/MacMono/MacMono_gps.rds")
MacMono_gps
MacMono.cells <- readRDS("Consensus method/Figure 2/Assets/MacMono/cell.mat.split.rds")
MacMono.cells <- MacMono.cells[1:7]
names(MacMono.cells) <- paste0("MacMono_", names(MacMono.cells))


pdac.atlas@meta.data <- pdac.atlas@meta.data %>% 
  mutate(GP = case_when(
    rownames(.) %in% MacMono.cells[[1]] ~ "MacMono_SLC40A1",
    rownames(.) %in% MacMono.cells[[2]] ~ "MacMono_FCN1",
    rownames(.) %in% MacMono.cells[[3]] ~ "MacMono_APOE",
    rownames(.) %in% MacMono.cells[[4]] ~ "MacMono_MMP9",
    rownames(.) %in% MacMono.cells[[5]] ~ "MacMono_IL1A",
    rownames(.) %in% MacMono.cells[[6]] ~ "MacMono_ANXA1",
    rownames(.) %in% MacMono.cells[[7]] ~ "MacMono_STMN1",
    T ~ GP))



Neut.cells <- readRDS("Consensus method/Figure 2/Assets/Neutrophil/cell.mat.split.rds")
Neut.cells <- Neut.cells[1:3]
names(Neut.cells) <- paste0("Neut_", names(Neut.cells))


pdac.atlas@meta.data <- pdac.atlas@meta.data %>% 
  mutate(GP = case_when(
    rownames(.) %in% Neut.cells[[1]] ~ "Neut_S100A8",
    rownames(.) %in% Neut.cells[[2]] ~ "Neut_CXCR4",
    rownames(.) %in% Neut.cells[[3]] ~ "Neut_GBP1",
    T ~ GP))



DCs.cells <- readRDS("Consensus method/Figure 2/Assets/DCs/cell.mat.split.rds")
DCs.cells <- DCs.cells[1:4]
names(DCs.cells) <- paste0("DCs_", names(DCs.cells))


pdac.atlas@meta.data <- pdac.atlas@meta.data %>% 
  mutate(GP = case_when(
    rownames(.) %in% DCs.cells[[1]] ~ "DC_CD1A",
    rownames(.) %in% DCs.cells[[2]] ~ "DC_CD14",
    rownames(.) %in% DCs.cells[[3]] ~ "DC_CLEC10A",
    rownames(.) %in% DCs.cells[[4]] ~ "DC_IL3RA",
    T ~ GP))




table.df <- pdac.atlas@meta.data %>% 
  filter(str_detect(GP, "_") & Tissue.Type == "Primary Tumour") %>% 
  rownames_to_column(var = "cell.id") %>% 
  select(Lvl1_Clusters, GP) %>% distinct()


table.df <- table(pdac.atlas$Lvl1_Clusters, pdac.atlas$GP) %>% as.data.frame()

rm(epi.cells, fib.cells, ecs.cells, cd4t.cells, cd8t.cells, nk.cells, B.cells, MacMono.cells, Neut.cells, DCs.cells, Mast.cells)
rm(b.gps, cd8t.gps, ecs.gps, epi.gps, fib.gps, MacMono_gps)


ce.list <- readRDS(paste0("Consensus method/Figure 4/Assets/ce.list.17102024.rds"))
ce.list

ce.df <- stack(ce.list)

ce.df <- ce.df %>% 
  column_to_rownames(var = "values") %>% 
  dplyr::rename(CE = ind)


pdac.atlas@meta.data <- pdac.atlas@meta.data %>% 
  mutate(CE = case_when(
    GP %in% ce.list[[1]] ~ "CE1",
    GP %in% ce.list[[2]] ~ "CE2",
    GP %in% ce.list[[3]] ~ "CE3",
    GP %in% ce.list[[4]] ~ "CE4",
    T ~ "Other"
  ))

unique(pdac.atlas@meta.data$CE)

pdac.ce <- subset(pdac.atlas, subset = CE != "Other" & Tissue.Type == "Primary Tumour")

rm(pdac.atlas);gc()


Idents(pdac.ce) <- "CE"

pdac.ce


library(presto)
wilcox.res <- wilcoxauc(pdac.ce, group_by = "CE")
write.csv(wilcox.res, "Consensus method/Figure 4/Assets/ce.wilcox.res.csv")


interaction_input <- read.csv(file = 'Consensus method/Figure 4/Dataframes/interaction_input_CellChatDB.csv', row.names = 1)


all_data_ligands <- wilcox.res %>% 
  mutate(CE = group) %>% 
  filter(feature %in% unique(interaction_input$ligand)) %>%
  mutate(pct_diff = pct_in - pct_out) %>% 
  filter(pct_in > 0.2 & pct_diff > 0.1 & logFC > 0.1 & padj < 0.05) %>% 
  group_by(CE) %>% 
  arrange(desc(logFC), .by_group = T) %>% 
  mutate(ligand = feature) %>% 
  select(CE, ligand)


all_data_receptors <- wilcox.res %>% 
  mutate(CE = group) %>% 
  filter(feature %in% unique(interaction_input$receptor)) %>%
  mutate(pct_diff = pct_in - pct_out) %>% 
  filter(pct_in > 0.2 & pct_diff > 0.1 & logFC > 0.1& padj < 0.05) %>% 
  group_by(CE) %>% 
  arrange(desc(logFC), .by_group = T) %>% 
  mutate(receptor = feature) %>% 
  select(CE, receptor)


ligand_receptor_pairs <- expand.grid(all_data_ligands$ligand, all_data_receptors$receptor)
colnames(ligand_receptor_pairs) <- c("ligand", "receptor")

ligand_receptor_pairs <- ligand_receptor_pairs %>%
  left_join(all_data_ligands, by = "ligand") %>%
  left_join(all_data_receptors, by = "receptor")

# Combine ligand and receptor columns to create interaction pairs
ligand_receptor_pairs <- ligand_receptor_pairs %>%
  mutate(interaction_name = paste(ligand, receptor, sep = "_"))

unique(interaction_input$annotation)

# Merge the ligand-receptor pairs with the CellChat data based on interaction_name
merged_interactions <- interaction_input %>%
  filter(annotation %in% c("Cell-Cell Contact", "Secreted Signaling")) %>% 
  select(interaction_name) %>% 
  inner_join(ligand_receptor_pairs, by = "interaction_name") %>% 
  mutate(CE_LR = paste0(CE.x, "-", CE.y)) %>% 
  distinct()

lr_hm <- merged_interactions %>%
 dplyr::count(CE_LR) %>% 
  separate(CE_LR, sep = "-", into = c("CE.x", "CE.y")) %>% 
pivot_wider(id_cols = CE.x, names_from = CE.y, values_from = n) %>% 
  column_to_rownames(var = "CE.x") 

lr_hm[is.na(lr_hm)] <- 0
#lr_hm <- log2(lr_hm + 1)

Heatmap((lr_hm), cluster_rows =F, cluster_columns = F,
        column_order = paste0("CE", 1:4),row_title = "Ligands",
        column_title = "Receptors",border = T,
        row_order =  paste0("CE", 1:4), heatmap_legend_param = list(
          title = "L-R pairs (n)" ),
        row_names_side = "left",show_column_dend = F, show_row_dend = F,
        col  = viridis(n = 10 , option = "A"))

tiff("Consensus method/Figure 4/Assets/LR.heatmap.all.tiff", 
     width = 1200, height = 1000, res = 300)
draw(ComplexHeatmap::pheatmap(lr_hm, display_numbers = T, 
                              show_row_dend = F,
                              show_column_dend = F,color = viridis(n = 10 , option = "A"),border = T,
                              cluster_rows = F, cluster_cols = F,number_color = "black",border_color = NA,
                              column_order = paste0("CE", 1:4),
                              row_order =  paste0("CE", 1:4),
                              row_title = "Ligands",
                              column_title = "Receptors",
                              row_names_side = "left",
                              heatmap_legend_param = list(
                                title = "L-R pairs (n)" ),
                              show_rownames = T, show_colnames = T), heatmap_legend_side = "right") 
dev.off()

ce.list


# Merge the ligand-receptor pairs with the CellChat data based on interaction_name
merged_interactions <- interaction_input %>%
  #filter(annotation == "Cell-Cell Contact") %>% 
  select(interaction_name) %>% 
  inner_join(ligand_receptor_pairs, by = "interaction_name") %>% 
  mutate(CE_LR = paste0(CE.x, "-", CE.y)) %>% 
  filter(CE_LR %in% c("CE1-CE1", "CE2-CE2", "CE3-CE3", "CE4-CE4")) %>% 
  distinct()


merged_interactions_hm <- merged_interactions %>% 
 filter(str_detect(interaction_name, "HLA|IL|CXCL|CCL|LGAL|PLAU|SERP|C3|IGF|FGF|DLL|SELL|EGF|EREG|CSF|KLK|CTLA4|COL|WNT|MDK|CD22|TNF_|SIGLEC1")) %>% 
  mutate(CE_LR = factor(CE_LR, levels = c("CE1-CE1", "CE2-CE2", "CE3-CE3", "CE4-CE4"))) %>% 
  group_by(CE_LR) %>% 
  dplyr::count(interaction_name) %>% 
  slice_head(n = 50) %>% 
  pivot_wider(id_cols = CE_LR, names_from = interaction_name, values_from = n) %>% 
  column_to_rownames(var = "CE_LR")

top_ann <- interaction_input %>% 
  filter(interaction_name %in% colnames(merged_interactions_hm)) %>% 
  select(annotation) %>% t() %>% as.data.frame() %>% 
  select(colnames(merged_interactions_hm)) %>% t() %>% as.data.frame()


merged_interactions_hm <- merged_interactions_hm %>% 
  select(rownames(top_ann))

merged_interactions_hm <- ifelse(merged_interactions_hm == 1, "Yes", merged_interactions_hm)

merged_interactions_hm[is.na(merged_interactions_hm)] <- "No"

identical(colnames(merged_interactions_hm), rownames(top_ann))

rownames(merged_interactions_hm) <- paste0("CE", 1:4)


tiff("Consensus method/Figure 4/Assets/LR.heatmap.tiff", width = 4200, height = 1300, res = 300)
Heatmap(merged_interactions_hm, cluster_columns = F, col = c("grey", "#377EB8"),border = T,
        row_names_side = "left", column_title = "Enriched L-R pairs within each CE",
        top_annotation = HeatmapAnnotation(df = top_ann,
                                           simple_anno_size = unit(0.6, "cm"),
                                           show_legend = T,
                                           #  col = list(CM = CM_cols), 
                                           border = T, show_annotation_name = F), 
        heatmap_legend_param = list(title = "Enriched"),
        #row_order = c("CM1-CM1", "CM2-CM2", "CM3-CM3", "CM4-CM4", "CM5-CM5", "CM6-CM6", 
        #             "CM7-CM7", "CM8-CM8", "CM9-CM9", "CM10-CM10", "CM11-CM11", "CM12-CM12"), 
        cluster_rows = F)
dev.off()

merged_interactions <- interaction_input %>%
  #filter(annotation == "Cell-Cell Contact") %>% 
  select(interaction_name) %>% 
  inner_join(ligand_receptor_pairs, by = "interaction_name") %>% 
  mutate(CE_LR = paste0(CE.x, "-", CE.y)) %>% 
  distinct()

merged_interactions_hm <- merged_interactions %>% 
  filter(str_detect(interaction_name, "CXCL|CCL|CSF|PLAU|SERPINE1|MDK|WNT5A|C3|LGALS9|IL1|IL22")) %>% 
  group_by(CE_LR) %>% 
  dplyr::count(ligand) %>% 
  slice_head(n = 10) %>% 
  arrange(desc(n), .by_group = T)

length(unique(merged_interactions_hm$ligand))


pdac.ce <- pdac.ce %>% ScaleData(features = unique(merged_interactions_hm$ligand))

pdac.ce@meta.data <- pdac.ce@meta.data %>% 
  mutate(GP = factor(GP, levels = rownames(ce.df)))

Idents(pdac.ce) <- "GP"

ce.df <- ce.df %>% 
  mutate(Colour= case_when(
    CE == "CE1" ~ CE_cols[[1]],
    CE == "CE2" ~ CE_cols[[2]],
    CE == "CE3" ~ CE_cols[[3]],
    CE == "CE4" ~ CE_cols[[4]],
    CE == "CE5" ~ CE_cols[[5]],
    T ~ "Other"
  ))

tiff("Consensus method/Figure 4/Assets/ccl.ligands.dotplot.tiff", width = 2500, height = 2200, res = 300)
DotPlot(pdac.ce, scale = T,group.by = "GP",
        features = rev(unique(merged_interactions_hm$ligand)))+
  coord_flip()+
  xlab("")+ylab("")+
  scale_size(range = c(2,6))+
  theme_light()+
  scale_x_discrete(position = "top")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size = 13, color = ce.df$Colour),
        axis.text.y = element_text(size = 13))+
  theme(legend.position = "top")+
  scale_colour_gradient2(low = "#053061", mid = "white", high = "red")+
  ggtitle("Secreted ligands")
dev.off()

merged_interactions_hm <- merged_interactions %>% 
  filter(str_detect(interaction_name, "CXCL|CCL|CSF|PLAU|SERPINE1|MDK|WNT5A|C3|LGALS9|IL1|IL22")) %>% 
  group_by(CE_LR) %>% 
  dplyr::count(receptor) %>% 
  slice_head(n = 10) %>% 
  arrange(desc(n), .by_group = T)


tiff("Consensus method/Figure 4/Assets/ccl.receptors.dotplot.tiff", width = 2500, height = 2200, res = 300)
DotPlot(pdac.ce, scale = T,
        features = rev(unique(merged_interactions_hm$receptor)))+
  coord_flip()+
  xlab("")+ylab("")+
  scale_size(range = c(2,6))+
  theme_light()+
  #scale_x_discrete(position = "top")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size = 13, color = ce.df$Colour),
        axis.text.y = element_text(size = 13))+
  theme(legend.position = "top")+
  scale_colour_gradient2(low = "#053061", mid = "white", high = "red")
dev.off()


check.ligands <- c("PDCD1LG2", "CD274", "CD276", "PVR", "NECTIN2","LGALS9", "TNFRSF14", 
                   "FASLG", "ICOSLG", "CD80", "TNFSF18", "CD86", "TNFSF9", "CD40LG", "TNFSF4","TNFSF14", "CD70")



tiff("Consensus method/Figure 5/Assets/CC.ligands.dotplot.tiff", width = 2500, height = 1800, res = 300)
DotPlot(pdac.ce, scale = T,group.by = "GP",
        features = check.ligands)+
  coord_flip()+
  xlab("")+ylab("")+
  scale_size(range = c(2,6))+
  theme_light()+
  scale_x_discrete(position = "top")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size = 13, color = ce.df$Colour),
        axis.text.y = element_text(size = 13))+
  theme(legend.position = "top")+
  scale_colour_gradient2(low = "#053061", mid = "white", high = "red")
dev.off()


check.receptors <- c("PDCD1", "TREML2", "TIGIT", "CD226", "HAVCR2", "BTLA", "CD160","FAS", "CTLA4","LAG3", "TNFRSF18","TNFRSF9", "CD40", "TNFRSF4","TNFRSF14", "LTBR","CD27")


tiff("Consensus method/Figure 5/Assets/cc.receptors.dotplot.tiff", width = 2500, height = 1800, res = 300)
DotPlot(pdac.ce, scale = T,
        features = check.receptors)+
  coord_flip()+
  xlab("")+ylab("")+
  scale_size(range = c(2,6))+
  theme_light()+
  #scale_x_discrete(position = "top")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size = 13, color = ce.df$Colour),
        axis.text.y = element_text(size = 13))+
  theme(legend.position = "top")+
  scale_colour_gradient2(low = "#053061", mid = "white", high = "red")
dev.off()



merged_interactions_hm <- merged_interactions %>% 
  filter(ligand %in% check.ligands) %>% 
  group_by(CE_LR) %>% 
  dplyr::count(ligand) %>% 
  slice_head(n = 15)



tiff("Consensus method/Figure 5/Assets/CC.ligands.dotplot.SIGN.tiff", width = 2500, height = 1200, res = 300)
DotPlot(pdac.ce, scale = T,group.by = "GP",
        features = rev(unique(merged_interactions_hm$ligand)))+
  coord_flip()+
  xlab("")+ylab("")+
  scale_size(range = c(2,6))+
  theme_light()+
  scale_x_discrete(position = "top")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size = 13, color = ce.df$Colour),
        axis.text.y = element_text(size = 13))+
  theme(legend.position = "top")+
  scale_colour_gradient2(low = "#053061", mid = "white", high = "red")
dev.off()



merged_interactions_hm <- merged_interactions %>% 
  filter(receptor %in% check.receptors) %>% 
  group_by(CE_LR) %>% 
  dplyr::count(receptor) %>% 
  slice_head(n = 15)



tiff("Consensus method/Figure 5/Assets/cc.receptors.dotplot.SIGN.tiff", width = 2500, height = 1200, res = 300)
DotPlot(pdac.ce, scale = T,
        features = rev(unique(merged_interactions_hm$receptor)))+
  coord_flip()+
  xlab("")+ylab("")+
  scale_size(range = c(2,6))+
  theme_light()+
  #scale_x_discrete(position = "top")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size = 13, color = ce.df$Colour),
        axis.text.y = element_text(size = 13))+
  theme(legend.position = "top")+
  scale_colour_gradient2(low = "#053061", mid = "white", high = "red")
dev.off()








