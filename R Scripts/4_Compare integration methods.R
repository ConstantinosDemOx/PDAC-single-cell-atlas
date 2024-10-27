# Compare integration methods


library(Seurat)
library(tidyverse)
library(speckle)
library(limma)
library(lisi)


pdac.atlas <- readRDS("New way/Figure_1/Files/pdac.atlas_final.rds")
pdac.atlas


gc()

set.seed(1000)
sample.cells <- pdac.atlas@meta.data %>% 
  rownames_to_column(var = "Cell.ID") %>% 
  select(Cell.ID, Lvl2_Clusters) %>% 
  distinct() %>% 
  group_by(Lvl2_Clusters) %>% 
  sample_n(size = 10000, replace = T)

length(sample.cells$Cell.ID)


sample.pdac <- subset(pdac.atlas, cells = sample.cells$Cell.ID)

rm(pdac.atlas);gc()

sample.pdac <- sample.pdac %>% 
  NormalizeData() %>% FindVariableFeatures() %>% 
  ScaleData() %>% RunPCA() %>% 
  RunUMAP(reduction = "pca", dims = 1:40)


sample.pdac@reductions$umap@key <- "UMAP_"

colnames(x = sample.pdac[["umap"]]@cell.embeddings) <- paste0("UMAP_", 1:2)

sample.pdac
unique(sample.pdac$Lvl2_Clusters)

sample.pdac@meta.data <- sample.pdac@meta.data %>% 
  mutate(Lvl1_Clusters = case_when(
    str_detect(Lvl2_Clusters, "CD4T|CD8T|NK") ~ "T/NK",
    str_detect(Lvl2_Clusters, "Monocyte|Macrophage|DC|Neutrophil") ~ "Myeloid",
    str_detect(Lvl2_Clusters, "Ductal|Acinar|Hepatocyte") ~ "Epithelial\n(non-malignant)",
    str_detect(Lvl2_Clusters, "Malignant|Tuft") ~ "Epithelial\n(malignant)",
    TRUE ~ Lvl2_Clusters))

sample.pdac@meta.data <- sample.pdac@meta.data %>% 
  mutate(Lvl1_Clusters = factor(Lvl1_Clusters, levels = c(
    "Epithelial\n(malignant)", "Epithelial\n(non-malignant)", "Endocrine", "Endothelial",
    "Stellate", "Fibroblast", "Schwann", "Myeloid", "Mast", "Plasma", "B", "T/NK"
  )))

sample.pdac


library(scCustomize);library(rcartocolor)
tiff("New way/Figure_1/Integration comparisons/Unintegrated_CT.tiff",   width = 1200, height = 1200, res = 300)
DimPlot_scCustom(sample.pdac,reduction = "umap",group.by = "Lvl1_Clusters", label = F, raster = F, pt.size = 0.4, 
                 alpha= 0.1,label.size = 5,
                 aspect_ratio = 1, 
                 colors_use = carto_pal(n = 12, name = "Vivid")) +
  ggtitle("Unintegrated") + NoLegend()
dev.off()


sample.pdac@meta.data <- sample.pdac@meta.data %>% 
  mutate(Study = factor(Study , levels = c(
    "Zhou et al", "Werba et al", "Peng et al", "Zhang et al",
    "Steele et al", "Oh et al", "Lin et al", "Xue et al", "Moncada et al"
  )))

library(scCustomize);library(rcartocolor);library(RColorBrewer)
tiff("New way/Figure_1/Integration comparisons/Unintegrated_Patient.tiff",   width = 1200, height = 1200, res = 300)
DimPlot_scCustom(sample.pdac,reduction = "umap",group.by = "Study.Patient.ID", label = F, raster = F, pt.size = 0.4, 
                 alpha= 0.1,label.size = 5,
                 aspect_ratio = 1, 
                 #colors_use =brewer.pal(n = 9 , name = "Set3")
                 ) +
  ggtitle("Unintegrated") + NoLegend()
dev.off()


pdac.atlas.coord <- sample.pdac[["umap"]]@cell.embeddings


pdac.variables <- sample.pdac@meta.data %>% 
  select(Lvl1_Clusters, Study, Study.Patient.ID)

identical(rownames(pdac.atlas.coord), rownames(pdac.variables))


unintegrated_res <- compute_lisi(pdac.atlas.coord, pdac.variables, c('Lvl1_Clusters', 'Study.Patient.ID'))
mean(unintegrated_res$Lvl1_Clusters)
mean(unintegrated_res$Study.Patient.ID)

unintegrated_res <- unintegrated_res %>% 
  mutate(Method = "Unintegrated")

write.csv(unintegrated_res, "New way/Figure_1/Integration comparisons/unintegrated.lisi.csv")


# scVI integration

sample.pdac <- sample.pdac %>% 
 RunUMAP(reduction = "scVI", dims = 1:40)


sample.pdac@reductions$umap@key <- "UMAP_"

colnames(x = sample.pdac[["umap"]]@cell.embeddings) <- paste0("UMAP_", 1:2)



tiff("New way/Figure_1/Integration comparisons/scVI_CT.legend.tiff",   width = 1200, height = 1200, res = 300)
DimPlot_scCustom(sample.pdac,reduction = "umap",group.by = "Lvl1_Clusters", label = F, raster = F, pt.size = 0.4, 
                 alpha= 0.1,label.size = 5,
                 aspect_ratio = 1, 
                 colors_use = carto_pal(n = 12, name = "Vivid")) +
  ggtitle("scVI") 
dev.off()


tiff("New way/Figure_1/Integration comparisons/scVI_Study.tiff",   width = 1200, height = 1200, res = 300)
DimPlot_scCustom(sample.pdac,reduction = "umap",group.by = "Study", label = F, raster = F, pt.size = 0.4, 
                 alpha= 0.1,label.size = 5,
                 aspect_ratio = 1, 
                 colors_use =brewer.pal(n = 9 , name = "Set3")
) +
  ggtitle("scVI") + NoLegend()
dev.off()


pdac.atlas.coord <- sample.pdac[["umap"]]@cell.embeddings


pdac.variables <- sample.pdac@meta.data %>% 
  select(Lvl1_Clusters, Study.Patient.ID)

identical(rownames(pdac.atlas.coord), rownames(pdac.variables))


scvi_res <- compute_lisi(pdac.atlas.coord, pdac.variables, c('Lvl1_Clusters', 'Study.Patient.ID'))

mean(scvi_res$Lvl1_Clusters)
mean(scvi_res$Study.Patient.ID)

scvi_res <- scvi_res %>% 
  mutate(Method = "scVI")


write.csv(scvi_res, "New way/Figure_1/Integration comparisons/scvi.lisi.csv")




# scANVI integration

sample.pdac <- sample.pdac %>% 
  RunUMAP(reduction = "scANVI", dims = 1:40)



sample.pdac@reductions$umap@key <- "UMAP_"

colnames(x = sample.pdac[["umap"]]@cell.embeddings) <- paste0("UMAP_", 1:2)



tiff("New way/Figure_1/Integration comparisons/scANVI_CT.tiff",   width = 1200, height = 1200, res = 300)
DimPlot_scCustom(sample.pdac,reduction = "umap",group.by = "Lvl1_Clusters", label = F, raster = F, pt.size = 0.4, 
                 alpha= 0.1,label.size = 5,
                 aspect_ratio = 1, 
                 colors_use = carto_pal(n = 12, name = "Vivid")) +
  ggtitle("scANVI") + NoLegend()
dev.off()

tiff("New way/Figure_1/Integration comparisons/scANVI_Patient.tiff",   width = 1200, height = 1200, res = 300)
DimPlot_scCustom(sample.pdac,reduction = "umap",group.by = "Study.Patient.ID", label = F, raster = F, pt.size = 0.4, 
                 alpha= 0.1,label.size = 5,
                 aspect_ratio = 1, 
                 #colors_use =brewer.pal(n = 9 , name = "Set3")
) +
  ggtitle("scANVI") + NoLegend()
dev.off()



pdac.atlas.coord <- sample.pdac[["umap"]]@cell.embeddings


pdac.variables <- sample.pdac@meta.data %>% 
  select(Lvl1_Clusters, Study.Patient.ID)

identical(rownames(pdac.atlas.coord), rownames(pdac.variables))


scANVI_res <- compute_lisi(pdac.atlas.coord, pdac.variables, c('Lvl1_Clusters', 'Study.Patient.ID'))

mean(scANVI_res$Lvl1_Clusters)
mean(scANVI_res$Study.Patient.ID)

scANVI_res <- scANVI_res %>% 
  mutate(Method = "scANVI")


write.csv(scANVI_res, "New way/Figure_1/Integration comparisons/scANVI.lisi.csv")


# Integrate with harmony

library(harmony)

sample.pdac <- RunHarmony(sample.pdac, "Study.Patient.ID")

sample.pdac <- sample.pdac %>% 
  RunUMAP(reduction = "harmony", dims = 1:40)


sample.pdac@reductions$umap@key <- "UMAP_"

colnames(x = sample.pdac[["umap"]]@cell.embeddings) <- paste0("UMAP_", 1:2)



tiff("New way/Figure_1/Integration comparisons/Harmony_CT.tiff",   width = 1200, height = 1200, res = 300)
DimPlot_scCustom(sample.pdac,reduction = "umap",group.by = "Lvl1_Clusters", label = F, raster = F, pt.size = 0.4, 
                 alpha= 0.1,label.size = 5,
                 aspect_ratio = 1, 
                 colors_use = carto_pal(n = 12, name = "Vivid")) +
  ggtitle("Harmony") + NoLegend()
dev.off()


tiff("New way/Figure_1/Integration comparisons/Harmony_Patient.tiff",   width = 1200, height = 1200, res = 300)
DimPlot_scCustom(sample.pdac,reduction = "umap",group.by = "Study.Patient.ID", label = F, raster = F, pt.size = 0.4, 
                 alpha= 0.1,label.size = 5,
                 aspect_ratio = 1, 
                 #colors_use =brewer.pal(n = 9 , name = "Set3")
) +
  ggtitle("Harmony") + NoLegend()
dev.off()


pdac.atlas.coord <- sample.pdac[["umap"]]@cell.embeddings


pdac.variables <- sample.pdac@meta.data %>% 
  select(Lvl1_Clusters, Study.Patient.ID)

identical(rownames(pdac.atlas.coord), rownames(pdac.variables))


Harmony_res <- compute_lisi(pdac.atlas.coord, pdac.variables, c('Lvl1_Clusters', 'Study.Patient.ID'))

mean(Harmony_res$Lvl1_Clusters)
mean(Harmony_res$Study.Patient.ID)

Harmony_res <- Harmony_res %>% 
  mutate(Method = "Harmony")


write.csv(Harmony_res, "New way/Figure_1/Integration comparisons/Harmony.lisi.csv")





# Perform integration with Seurat RPCA


pdac.list <- SplitObject(sample.pdac, split.by = "Study")

# normalize and identify variable features for each dataset independently
pdac.list <- lapply(X = pdac.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})


pdac.list <- pdac.list[[-c("Steele_PT_16", "Lin_MT01", "Lin_PT02", "Lin_PT10", "Werba_M_PT_PT06", "Zhang_MET_PT4")]]
pdac.list <- pdac.list[setdiff(names(pdac.list), c("Steele_PT_P16", "Lin_MT01", "Lin_PT02", "Lin_PT10", "Werba_M_PT_PT06", "Zhang_MET_PT4"))]

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = pdac.list)
pdac.list <- lapply(X = pdac.list, FUN = function(x) {
  x <- ScaleData(x, features = features)
  x <- RunPCA(x, features = features)
})


pdac.anchors <- FindIntegrationAnchors(object.list = pdac.list, anchor.features = features, reduction = "rpca",
                                       reference = c(1, 2,3))

pdac.rpca <- IntegrateData(anchorset = pdac.anchors)


DefaultAssay(pdac.rpca) <- "integrated"

pdac.rpca <- pdac.rpca %>% 
  ScaleData() %>% RunPCA() %>% RunUMAP(dims = 1:40)


pdac.rpca@reductions$umap@key <- "UMAP_"

colnames(x = pdac.rpca[["umap"]]@cell.embeddings) <- paste0("UMAP_", 1:2)

pdac.rpca@meta.data <- pdac.rpca@meta.data %>% 
  mutate(Lvl1_Clusters = case_when(
    str_detect(Lvl2_Clusters, "CD4T|CD8T|NK") ~ "T/NK",
    str_detect(Lvl2_Clusters, "Monocyte|Macrophage|DC|Neutrophil") ~ "Myeloid",
    str_detect(Lvl2_Clusters, "Ductal|Acinar|Hepatocyte") ~ "Epithelial\n(non-malignant)",
    str_detect(Lvl2_Clusters, "Malignant|Tuft") ~ "Epithelial\n(malignant)",
    TRUE ~ Lvl2_Clusters))

pdac.rpca@meta.data <- pdac.rpca@meta.data %>% 
  mutate(Lvl1_Clusters = factor(Lvl1_Clusters, levels = c(
    "Epithelial\n(malignant)", "Epithelial\n(non-malignant)", "Endocrine", "Endothelial",
    "Stellate", "Fibroblast", "Schwann", "Myeloid", "Mast", "Plasma", "B", "T/NK"
  )))


pdac.rpca@meta.data <- pdac.rpca@meta.data %>% 
  mutate(Study = factor(Study , levels = c(
    "Zhou et al", "Werba et al", "Peng et al", "Zhang et al",
    "Steele et al", "Oh et al", "Lin et al", "Xue et al", "Moncada et al"
  )))


tiff("New way/Figure_1/Integration comparisons/Seurat_CT.tiff",   width = 1200, height = 1200, res = 300)
DimPlot_scCustom(pdac.rpca,reduction = "umap",group.by = "Lvl1_Clusters", label = F, raster = F, pt.size = 0.4, 
                 alpha= 0.1,label.size = 5,
                 aspect_ratio = 1, 
                 colors_use = carto_pal(n = 12, name = "Vivid")) +
  ggtitle("Seurat") + NoLegend()
dev.off()


tiff("New way/Figure_1/Integration comparisons/Seurat_Patient.tiff",   width = 1200, height = 1200, res = 300)
DimPlot_scCustom(pdac.rpca,reduction = "umap",group.by = "Study.Patient.ID", label = F, raster = F, pt.size = 0.4, 
                 alpha= 0.1,label.size = 5,
                 aspect_ratio = 1, 
                 #colors_use =brewer.pal(n = 9 , name = "Set3")
) +
  ggtitle("Seurat") + NoLegend()
dev.off()


pdac.atlas.coord <- pdac.rpca[["umap"]]@cell.embeddings


pdac.variables <- pdac.rpca@meta.data %>% 
  select(Lvl1_Clusters, Study.Patient.ID)

identical(rownames(pdac.atlas.coord), rownames(pdac.variables))


Seurat_res <- compute_lisi(pdac.atlas.coord, pdac.variables, c('Lvl1_Clusters', 'Study.Patient.ID'))

mean(Seurat_res$Lvl1_Clusters)
mean(Seurat_res$Study.Patient.ID)

Seurat_res <- Seurat_res %>% 
  mutate(Method = "Seurat")


write.csv(Seurat_res, "New way/Figure_1/Integration comparisons/Seurat.lisi.csv")

rm(pdac.rpca);gc()
# Perform fastMNN integration

library(batchelor)
library(SeuratWrappers)

pdac.fmnn <- RunFastMNN(object.list = pdac.list)

pdac.fmnn <- RunUMAP(pdac.fmnn, reduction = "mnn", dims = 1:40)


pdac.fmnn@reductions$umap@key <- "UMAP_"

colnames(x = pdac.fmnn[["umap"]]@cell.embeddings) <- paste0("UMAP_", 1:2)

pdac.fmnn@meta.data <- pdac.fmnn@meta.data %>% 
  mutate(Lvl1_Clusters = case_when(
    str_detect(Lvl2_Clusters, "CD4T|CD8T|NK") ~ "T/NK",
    str_detect(Lvl2_Clusters, "Monocyte|Macrophage|DC|Neutrophil") ~ "Myeloid",
    str_detect(Lvl2_Clusters, "Ductal|Acinar|Hepatocyte") ~ "Epithelial\n(non-malignant)",
    str_detect(Lvl2_Clusters, "Malignant|Tuft") ~ "Epithelial\n(malignant)",
    TRUE ~ Lvl2_Clusters))

pdac.fmnn@meta.data <- pdac.fmnn@meta.data %>% 
  mutate(Lvl1_Clusters = factor(Lvl1_Clusters, levels = c(
    "Epithelial\n(malignant)", "Epithelial\n(non-malignant)", "Endocrine", "Endothelial",
    "Stellate", "Fibroblast", "Schwann", "Myeloid", "Mast", "Plasma", "B", "T/NK"
  )))


pdac.fmnn@meta.data <- pdac.fmnn@meta.data %>% 
  mutate(Study = factor(Study , levels = c(
    "Zhou et al", "Werba et al", "Peng et al", "Zhang et al",
    "Steele et al", "Oh et al", "Lin et al", "Xue et al", "Moncada et al"
  )))


tiff("New way/Figure_1/Integration comparisons/fastMNN_CT.tiff",   width = 1200, height = 1200, res = 300)
DimPlot_scCustom(pdac.fmnn,reduction = "umap",group.by = "Lvl1_Clusters", label = F, raster = F, pt.size = 0.4, 
                 alpha= 0.1,label.size = 5,
                 aspect_ratio = 1, 
                 colors_use = carto_pal(n = 12, name = "Vivid")) +
  ggtitle("fastMNN") + NoLegend()
dev.off()


tiff("New way/Figure_1/Integration comparisons/fastMNN_Patient.tiff",   width = 1200, height = 1200, res = 300)
DimPlot_scCustom(pdac.fmnn,reduction = "umap",group.by = "Study.Patient.ID", label = F, raster = F, pt.size = 0.4, 
                 alpha= 0.1,label.size = 5,
                 aspect_ratio = 1, 
                 #colors_use =brewer.pal(n = 9 , name = "Set3")
) +
  ggtitle("FastMNN") + NoLegend()
dev.off()


pdac.atlas.coord <- pdac.fmnn[["umap"]]@cell.embeddings


pdac.variables <- pdac.fmnn@meta.data %>% 
  select(Lvl1_Clusters, Study.Patient.ID)

identical(rownames(pdac.atlas.coord), rownames(pdac.variables))


fastMNN_res <- compute_lisi(pdac.atlas.coord, pdac.variables, c('Lvl1_Clusters', 'Study.Patient.ID'))

mean(fastMNN_res$Lvl1_Clusters)
mean(fastMNN_res$Study.Patient.ID)

fastMNN_res <- fastMNN_res %>% 
  mutate(Method = "fastMNN")


write.csv(fastMNN_res, "New way/Figure_1/Integration comparisons/fastMNN.lisi.csv")


# Plot all fastMNN LISI scores



all_res <- rbind(unintegrated_res, Seurat_res, Harmony_res, fastMNN_res, scvi_res, scANVI_res)

write.csv(all_res,"New way/Figure_1/Integration comparisons/all_lisi_results_15_07_2024.csv")

all_res <- read.csv("New way/Figure_1/Integration comparisons/all_lisi_results_15_07_2024.csv")

library(ggpubr)
library(tidyverse)
ggboxplot(all_res, "Method", "Study")
ggboxplot(all_res, "Method", "Lvl1_Clusters")


all_res_df <- all_res %>% 
  group_by(Method) %>% 
  summarise(Mean_iLISI_Sc = (mean(Study) - min(Study)) / (max(Study) - min(Study)),
            Mean_cLISI_Sc = 1 - ((mean(Lvl1_Clusters) - min(Lvl1_Clusters)) / (max(Lvl1_Clusters) - min(Lvl1_Clusters))),
            Median_iLISI_Sc = (median(Study) - min(Study)) / (max(Study) - min(Study)),
            Median_cLISI_Sc = 1 - ((median(Lvl1_Clusters) - min(Lvl1_Clusters)) / (max(Lvl1_Clusters) - min(Lvl1_Clusters))))


LISI_summary_df <- list()


for (i in unique(all_res$Method)) {
  
  print(i)
  
  lisi.df <- all_res %>% 
    filter(Method == i) %>% 
    summarise(Mean_iLISI = mean(Study),
              Min_iLISI = min(Study),
              Max_iLISI = max(Study),
              Mean_cLISI = mean(Lvl1_Clusters),
              Min_cLISI = min(Lvl1_Clusters),
              Max_cLISI = max(Lvl1_Clusters))
  
  lisi.df$iLISI <- scale(lisi.df$Mean_iLISI, center = lisi.df$Min_iLISI, scale = lisi.df$Max_iLISI - lisi.df$Min_iLISI)
  lisi.df$iLISI <-(lisi.df$Mean_iLISI -lisi.df$Min_iLISI) / (lisi.df$Max_iLISI - lisi.df$Min_iLISI)
  
  
  
  
  print(min(lisi.df$iLISI))
  print(max(lisi.df$iLISI))
  
  
  lisi.df$cLISI <- 1- scale(lisi.df$Lvl1_Clusters, center = min(lisi.df$Lvl1_Clusters), scale = max(lisi.df$Lvl1_Clusters) - min(lisi.df$Lvl1_Clusters))
  
  print(min(lisi.df$cLISI))
  print(max(lisi.df$cLISI))
  
  
  lisi.df <- lisi.df %>% 
    mutate(Method = i)
  
  LISI_summary_df[[i]] <- lisi.df
  
  
  
  
  
}


LISI_summary_df <- do.call(rbind, LISI_summary_df)

colnames(LISI_summary_df)


all_res_df <- LISI_summary_df %>% 
  group_by(Method) %>% 
  summarise(Mean_iLISI = mean(iLISI),
            Mean_cLISI = mean(cLISI))

library(ggpubr)
library(ggrepel)
all_res_df %>%
  ggscatter(x = "Median_iLISI_Sc", y = "Median_cLISI_Sc",font.label = 20,
            color = "Method")+
  # geom_errorbar(aes(ymin = Mean_cLISI - SD_cLISI,ymax = Mean_cLISI + SD_cLISI, color = Integration)) + 
  #geom_errorbarh(aes(xmin = Mean_iLISI - SD_iLISI,xmax = Mean_iLISI + SD_iLISI, color = Integration)) + 
  theme_bw()+
  ylab("Bio-conservation (1- cLISI)")+
  xlab("Batch correction (iLISI)")+
 # xlim(0,1)+ ylim(0,1)+
  ggtitle("Integration comparison") +theme(legend.position = "none",
                              text = element_text(size = 10),
                              plot.title = element_text(hjust = 0.5),
                              aspect.ratio = 1)+
  #ylab("-log10(p val)") + xlab("Pearson's R")+
  # scale_color_manual(values=gp_cols)+
  geom_text_repel(aes(label = Method, color = Method))






all_res_df <- all_res %>% 
  group_by(Method) %>% 
  summarise(Mean_iLISI = mean(Study),
            Mean_cLISI =  mean(1 - Lvl1_Clusters)) %>% 
  mutate(Mean_iLISI_sc = (Mean_iLISI - min(Mean_iLISI)) / (max(Mean_iLISI) - min(Mean_iLISI))) %>% 
  mutate(Mean_cLISI_sc = (Mean_cLISI - min(Mean_cLISI)) / (max(Mean_cLISI) - min(Mean_cLISI))) %>%
  ungroup()
#  mutate(Method = factor(Method, levels = c("Unintegrated", "Seurat", "Harmony", "fastMNN", "scVI", "scANVI")))

all_res_df$Score <- rowMeans(all_res_df[,c("Mean_iLISI_sc", "Mean_cLISI_sc")])



library(ggpubr)
library(ggrepel)
library(ggsci)

int_cols <- pal_npg(palette = "nrc")(6)
names(int_cols) <- unique(all_res_df$Method)

tiff("New way/Figure_1/Integration comparisons/Integration.comparisons.tiff", width = 900, height = 900, res = 300)
all_res_df %>%
  ggscatter(x = "Mean_iLISI_sc", y = "Mean_cLISI_sc",font.label = 20,
            color = "Method", palette = int_cols)+
  # geom_errorbar(aes(ymin = Mean_cLISI - SD_cLISI,ymax = Mean_cLISI + SD_cLISI, color = Integration)) + 
  #geom_errorbarh(aes(xmin = Mean_iLISI - SD_iLISI,xmax = Mean_iLISI + SD_iLISI, color = Integration)) + 
  theme_bw()+
  ylab("Bio-conservation (1- cLISI)")+
  xlab("Batch correction (iLISI)")+
  # xlim(0,1)+ ylim(0,1)+
  ggtitle("Integration comparison") +theme(legend.position = "none",
                                           text = element_text(size = 10),
                                           plot.title = element_text(hjust = 0.5),
                                           aspect.ratio = 1)+
  #ylab("-log10(p val)") + xlab("Pearson's R")+
  # scale_color_manual(values=gp_cols)+
  geom_text_repel(aes(label = Method, color = Method))
 
dev.off()


tiff("New way/Figure_1/Integration comparisons/Integration.score.tiff", width = 1500, height = 900, res = 300)
ggbarplot(all_res_df, "Method", order = rev(c("scANVI", "scVI", "Seurat","fastMNN", "Harmony", "Unintegrated")),"Score", 
          fill = "Method", sort.val = "asc", palette = int_cols) +
  theme(legend.position = "right")+
  coord_flip()
dev.off()
























