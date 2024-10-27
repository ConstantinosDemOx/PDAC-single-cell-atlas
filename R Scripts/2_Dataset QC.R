# 1) Load libraries

library(Seurat)
library(SeuratData)
library(tidyverse)


# 1) QC Zhou et al

zhou <- readRDS("Dataframes/Raw files/Zhou/Zhou_raw_merged.rds")

zhou

meta <- zhou@meta.data

# dataset name + patient ID
zhou <- PercentageFeatureSet(zhou, pattern = "^MT-", col.name = "Percent_mito")
zhou <- PercentageFeatureSet(zhou, "^HB[^(P)]", col.name = "Percent_hb")
zhou <- PercentageFeatureSet(zhou, "^R[SP]L", col.name = "Percent_ribo")
zhou <- PercentageFeatureSet(zhou, "^HSP", col.name = "Percent_hsp")

count.df <- as.data.frame(table(zhou$orig.ident))
colnames(count.df) <- c("Sample","Frequency")

ggplot(count.df, aes(x = Sample, y = Frequency, fill = Sample)) + 
  geom_col() +
  theme_classic() +
  geom_text(aes(label = Frequency), 
            position=position_dodge(width=0.9), 
            vjust=-0.25) +
  # scale_y_continuous(breaks = seq(0,5000, by = 1000), limits = c(0,5000)) +
  ggtitle("Raw: Cells per sample") +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))

VlnPlot(zhou,
        features = c("nCount_RNA"),
        ncol = 1,
        group.by = 'orig.ident',
        pt.size = 0.1) +
  geom_hline(yintercept = 100000) + NoLegend()

VlnPlot(zhou,
        features = c("Percent_mito"),
        ncol = 1,
        group.by = 'orig.ident',
        pt.size = 0.1) +
  geom_hline(yintercept = 25) + NoLegend()


zhou_filt <- subset(zhou, subset =  Percent_mito < 25 & nCount_RNA < 100000)
zhou
zhou_filt

FeatureScatter(zhou, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(zhou_filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 

saveRDS(zhou_filt, "Dataframes/Raw files/Zhou/Zhou_raw_filt.rds")

zhou_filt <- readRDS("Dataframes/Raw files/Zhou/Zhou_raw_filt.rds")

VlnPlot(zhou_filt,
        features = c("Percent_hb"),
        ncol = 1,
        group.by = 'orig.ident',
        pt.size = 0.1) +
  geom_hline(yintercept = 100000) + NoLegend()

zhou_filt <- subset(zhou_filt, subset =  Percent_mito < 25 & nCount_RNA < 50000 & Percent_hb < 10)
zhou
zhou_filt

FeatureScatter(zhou, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(zhou_filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 

saveRDS(zhou_filt, "Dataframes/Raw files/Zhou/Zhou_raw_filt.rds")



# 1) QC Lin et al

Lin <- readRDS("Dataframes/Raw files/Lin/Lin_raw_merged.rds")

Lin

meta <- Lin@meta.data

# dataset name + patient ID
Lin <- PercentageFeatureSet(Lin, pattern = "^MT-", col.name = "Percent_mito")
Lin <- PercentageFeatureSet(Lin, "^HB[^(P)]", col.name = "Percent_hb")
Lin <- PercentageFeatureSet(Lin, "^R[SP]L", col.name = "Percent_ribo")
Lin <- PercentageFeatureSet(Lin, "^HSP", col.name = "Percent_hsp")

count.df <- as.data.frame(table(Lin$orig.ident))
colnames(count.df) <- c("Sample","Frequency")

ggplot(count.df, aes(x = Sample, y = Frequency, fill = Sample)) + 
  geom_col() +
  theme_classic() +
  geom_text(aes(label = Frequency), 
            position=position_dodge(width=0.9), 
            vjust=-0.25) +
  # scale_y_continuous(breaks = seq(0,5000, by = 1000), limits = c(0,5000)) +
  ggtitle("Raw: Cells per sample") +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))

VlnPlot(Lin,
        features = c("nCount_RNA"),
        ncol = 1,
        group.by = 'orig.ident',
        pt.size = 0.1) +
  geom_hline(yintercept = 50000) + NoLegend()

VlnPlot(Lin,
        features = c("Percent_mito"),
        ncol = 1,
        group.by = 'orig.ident',
        pt.size = 0.1) +
  geom_hline(yintercept = 25) + NoLegend()


Lin_filt <- subset(Lin, subset =  Percent_mito < 25 & nCount_RNA < 50000)
Lin
Lin_filt

FeatureScatter(Lin, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(Lin_filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 

saveRDS(Lin_filt, "Dataframes/Raw files/Lin/Lin_raw_filt.rds")

Lin_filt <- readRDS("Dataframes/Raw files/Lin/Lin_raw_filt.rds")

VlnPlot(Lin_filt,
        features = c("Percent_hb"),
        ncol = 1,
        group.by = 'orig.ident',
        pt.size = 0.1) +
  geom_hline(yintercept = 25) + NoLegend()

Lin_filt <- subset(Lin_filt, subset =  Percent_hb < 5)


# 1) QC Oh et al

Oh <- readRDS("Dataframes/Raw files/Oh/Oh_raw_merged.rds")

Oh

meta <- Oh@meta.data

# dataset name + patient ID
Oh <- PercentageFeatureSet(Oh, pattern = "^MT-", col.name = "Percent_mito")
Oh <- PercentageFeatureSet(Oh, "^HB[^(P)]", col.name = "Percent_hb")
Oh <- PercentageFeatureSet(Oh, "^R[SP]L", col.name = "Percent_ribo")
Oh <- PercentageFeatureSet(Oh, "^HSP", col.name = "Percent_hsp")

count.df <- as.data.frame(table(Oh$orig.ident))
colnames(count.df) <- c("Sample","Frequency")

ggplot(count.df, aes(x = Sample, y = Frequency, fill = Sample)) + 
  geom_col() +
  theme_classic() +
  geom_text(aes(label = Frequency), 
            position=position_dodge(width=0.9), 
            vjust=-0.25) +
  # scale_y_continuous(breaks = seq(0,5000, by = 1000), limits = c(0,5000)) +
  ggtitle("Raw: Cells per sample") +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))

VlnPlot(Oh,
        features = c("nCount_RNA"),
        ncol = 1,
        group.by = 'orig.ident',
        pt.size = 0.1) +
  geom_hline(yintercept = 50000) + NoLegend()

VlnPlot(Oh,
        features = c("Percent_mito"),
        ncol = 1,
        group.by = 'orig.ident',
        pt.size = 0.1) +
  geom_hline(yintercept = 25) + NoLegend()


Oh_filt <- subset(Oh, subset =  Percent_mito < 25 & nCount_RNA < 30000)
Oh
Oh_filt

FeatureScatter(Oh, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(Oh_filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 

saveRDS(Oh_filt, "Dataframes/Raw files/Oh/Oh_raw_filt.rds")

Oh_filt <- readRDS("Dataframes/Raw files/Oh/Oh_raw_filt.rds")
Oh_filt


VlnPlot(Oh_filt,
        features = c("Percent_hb"),
        ncol = 1,
        group.by = 'orig.ident',
        pt.size = 0.1) 

Oh_filt <- subset(Oh_filt, subset =  Percent_hb < 2.5)





# 1) QC Peng et al

Peng <- readRDS("Dataframes/Raw files/Peng/Peng_raw_merged.rds")

Peng

meta <- Peng@meta.data

# dataset name + patient ID
Peng <- PercentageFeatureSet(Peng, pattern = "^MT-", col.name = "Percent_mito")
Peng <- PercentageFeatureSet(Peng, "^HB[^(P)]", col.name = "Percent_hb")
Peng <- PercentageFeatureSet(Peng, "^R[SP]L", col.name = "Percent_ribo")
Peng <- PercentageFeatureSet(Peng, "^HSP", col.name = "Percent_hsp")

count.df <- as.data.frame(table(Peng$orig.ident))
colnames(count.df) <- c("Sample","Frequency")

ggplot(count.df, aes(x = Sample, y = Frequency, fill = Sample)) + 
  geom_col() +
  theme_classic() +
  geom_text(aes(label = Frequency), 
            position=position_dodge(width=0.9), 
            vjust=-0.25) +
  # scale_y_continuous(breaks = seq(0,5000, by = 1000), limits = c(0,5000)) +
  ggtitle("Raw: Cells per sample") +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))

VlnPlot(Peng,
        features = c("nCount_RNA"),
        ncol = 1,
        group.by = 'orig.ident',
        pt.size = 0.1) +
  geom_hline(yintercept = 50000) + NoLegend()

VlnPlot(Peng,
        features = c("Percent_mito"),
        ncol = 1,
        group.by = 'orig.ident',
        pt.size = 0.1) +
  geom_hline(yintercept = 25) + NoLegend()


Peng_filt <- subset(Peng, subset =  Percent_mito < 30 & nCount_RNA < 50000)
Peng
Peng_filt

FeatureScatter(Peng, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(Peng_filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 

saveRDS(Peng_filt, "Dataframes/Raw files/Peng/Peng_raw_filt.rds")



Peng_filt <- readRDS("Dataframes/Raw files/Peng/Peng_raw_filt.rds")

VlnPlot(Peng_filt,
        features = c("Percent_hb"),
        ncol = 1,
        group.by = 'orig.ident',
        pt.size = 0.1) 

Peng_filt <- subset(Peng_filt, subset =  Percent_hb < 5)



# 1) QC Steele et al

Steele <- readRDS("Dataframes/Raw files/Steele/Steele_raw_merged.rds")

Steele

meta <- Steele@meta.data

# dataset name + patient ID
Steele <- PercentageFeatureSet(Steele, pattern = "^MT-", col.name = "Percent_mito")
Steele <- PercentageFeatureSet(Steele, "^HB[^(P)]", col.name = "Percent_hb")
Steele <- PercentageFeatureSet(Steele, "^R[SP]L", col.name = "Percent_ribo")
Steele <- PercentageFeatureSet(Steele, "^HSP", col.name = "Percent_hsp")

count.df <- as.data.frame(table(Steele$orig.ident))
colnames(count.df) <- c("Sample","Frequency")

ggplot(count.df, aes(x = Sample, y = Frequency, fill = Sample)) + 
  geom_col() +
  theme_classic() +
  geom_text(aes(label = Frequency), 
            position=position_dodge(width=0.9), 
            vjust=-0.25) +
  # scale_y_continuous(breaks = seq(0,5000, by = 1000), limits = c(0,5000)) +
  ggtitle("Raw: Cells per sample") +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))

VlnPlot(Steele,
        features = c("nCount_RNA"),
        ncol = 1,
        group.by = 'orig.ident',
        pt.size = 0.1) +
  geom_hline(yintercept = 50000) + NoLegend()

VlnPlot(Steele,
        features = c("Percent_mito"),
        ncol = 1,
        group.by = 'orig.ident',
        pt.size = 0.1) +
  geom_hline(yintercept = 30) + NoLegend()


Steele_filt <- subset(Steele, subset =  Percent_mito < 30 & nCount_RNA < 50000)
Steele
Steele_filt

FeatureScatter(Steele, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(Steele_filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 

saveRDS(Steele_filt, "Dataframes/Raw files/Steele/Steele_raw_filt.rds")


Steele_filt <- readRDS("Dataframes/Raw files/Steele/Steele_raw_filt.rds")

VlnPlot(Steele_filt,
        features = c("Percent_hb"),
        ncol = 1,
        group.by = 'orig.ident',
        pt.size = 0.1) 

Steele_filt <- subset(Steele_filt, subset =  Percent_hb < 5)





# 1) QC Werba et al

Werba <- readRDS("Dataframes/Raw files/Werba/Werba_raw_merged.rds")

Werba

meta <- Werba@meta.data

# dataset name + patient ID
Werba <- PercentageFeatureSet(Werba, pattern = "^MT-", col.name = "Percent_mito")
Werba <- PercentageFeatureSet(Werba, "^HB[^(P)]", col.name = "Percent_hb")
Werba <- PercentageFeatureSet(Werba, "^R[SP]L", col.name = "Percent_ribo")
Werba <- PercentageFeatureSet(Werba, "^HSP", col.name = "Percent_hsp")

count.df <- as.data.frame(table(Werba$orig.ident))
colnames(count.df) <- c("Sample","Frequency")

ggplot(count.df, aes(x = Sample, y = Frequency, fill = Sample)) + 
  geom_col() +
  theme_classic() +
  geom_text(aes(label = Frequency), 
            position=position_dodge(width=0.9), 
            vjust=-0.25) +
  # scale_y_continuous(breaks = seq(0,5000, by = 1000), limits = c(0,5000)) +
  ggtitle("Raw: Cells per sample") +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))

VlnPlot(Werba,
        features = c("nCount_RNA"),
        ncol = 1,
        group.by = 'orig.ident',
        pt.size = 0.1) +
  geom_hline(yintercept = 60000) + NoLegend()

VlnPlot(Werba,
        features = c("Percent_mito"),
        ncol = 1,
        group.by = 'orig.ident',
        pt.size = 0.1) +
  geom_hline(yintercept = 30) + NoLegend()


Werba_filt <- subset(Werba, subset =  Percent_mito < 30 & nCount_RNA < 60000)
Werba
Werba_filt2

FeatureScatter(Werba, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(Werba_filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 

saveRDS(Werba_filt, "Dataframes/Raw files/Werba/Werba_raw_filt.rds")

Werba_filt <- readRDS("Dataframes/Raw files/Werba/Werba_raw_filt.rds")

VlnPlot(Werba_filt,
        features = c("Percent_hb"),
        ncol = 1,
        group.by = 'orig.ident',
        pt.size = 0.1) 

Werba_filt <- Werba_filt %>% NormalizeData()

VlnPlot(Werba_filt,
        features = c("HBB"),
        ncol = 1,
        group.by = 'orig.ident',
        pt.size = 0.1) 




Werba_filt2 <- subset(Werba_filt, subset =  HBB < 5)


saveRDS(Werba_filt2, "Dataframes/Raw files/Werba/Werba_raw_filt.rds")


# 1) QC Xue et al

Xue <- readRDS("Dataframes/Raw files/Xue/Xue_raw_merged.rds")

Xue

meta <- Xue@meta.data

# dataset name + patient ID
Xue <- PercentageFeatureSet(Xue, pattern = "^MT-", col.name = "Percent_mito")
Xue <- PercentageFeatureSet(Xue, "^HB[^(P)]", col.name = "Percent_hb")
Xue <- PercentageFeatureSet(Xue, "^R[SP]L", col.name = "Percent_ribo")
Xue <- PercentageFeatureSet(Xue, "^HSP", col.name = "Percent_hsp")

count.df <- as.data.frame(table(Xue$orig.ident))
colnames(count.df) <- c("Sample","Frequency")

ggplot(count.df, aes(x = Sample, y = Frequency, fill = Sample)) + 
  geom_col() +
  theme_classic() +
  geom_text(aes(label = Frequency), 
            position=position_dodge(width=0.9), 
            vjust=-0.25) +
  # scale_y_continuous(breaks = seq(0,5000, by = 1000), limits = c(0,5000)) +
  ggtitle("Raw: Cells per sample") +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))

VlnPlot(Xue,
        features = c("nCount_RNA"),
        ncol = 1,
        group.by = 'orig.ident',
        pt.size = 0.1) +
  geom_hline(yintercept = 50000) + NoLegend()

VlnPlot(Xue,
        features = c("Percent_mito"),
        ncol = 1,
        group.by = 'orig.ident',
        pt.size = 0.1) +
  geom_hline(yintercept = 30) + NoLegend()


Xue_filt <- subset(Xue, subset =  Percent_mito < 30 & nCount_RNA < 50000)
Xue
Xue_filt

FeatureScatter(Xue, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(Xue_filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 

saveRDS(Xue_filt, "Dataframes/Raw files/Xue/Xue_raw_filt.rds")

Xue_filt <- readRDS("Dataframes/Raw files/Xue/Xue_raw_filt.rds")

VlnPlot(Xue_filt,
        features = c("Percent_hb"),
        ncol = 1,
        group.by = 'orig.ident',
        pt.size = 0.1) 

Xue_filt <- subset(Xue_filt, subset =  Percent_hb < 5)






# 1) QC Zhang et al

Zhang <- readRDS("Dataframes/Raw files/Zhang/Zhang_raw_merged.rds")

Zhang

meta <- Zhang@meta.data

# dataset name + patient ID
Zhang <- PercentageFeatureSet(Zhang, pattern = "^MT-", col.name = "Percent_mito")
Zhang <- PercentageFeatureSet(Zhang, "^HB[^(P)]", col.name = "Percent_hb")
Zhang <- PercentageFeatureSet(Zhang, "^R[SP]L", col.name = "Percent_ribo")
Zhang <- PercentageFeatureSet(Zhang, "^HSP", col.name = "Percent_hsp")

count.df <- as.data.frame(table(Zhang$orig.ident))
colnames(count.df) <- c("Sample","Frequency")

ggplot(count.df, aes(x = Sample, y = Frequency, fill = Sample)) + 
  geom_col() +
  theme_classic() +
  geom_text(aes(label = Frequency), 
            position=position_dodge(width=0.9), 
            vjust=-0.25) +
  # scale_y_continuous(breaks = seq(0,5000, by = 1000), limits = c(0,5000)) +
  ggtitle("Raw: Cells per sample") +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))

VlnPlot(Zhang,
        features = c("nCount_RNA"),
        ncol = 1,
        group.by = 'orig.ident',
        pt.size = 0.1) +
  geom_hline(yintercept = 50000) + NoLegend()

VlnPlot(Zhang,
        features = c("Percent_mito"),
        ncol = 1,
        group.by = 'orig.ident',
        pt.size = 0.1) +
  geom_hline(yintercept = 30) + NoLegend()


Zhang_filt <- subset(Zhang, subset =  Percent_mito < 30 & nCount_RNA < 50000)
Zhang
Zhang_filt

FeatureScatter(Zhang, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(Zhang_filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 

saveRDS(Zhang_filt, "Dataframes/Raw files/Zhang/Zhang_raw_filt.rds")

Zhang_filt <- readRDS("Dataframes/Raw files/Zhang/Zhang_raw_filt.rds")

VlnPlot(Zhang_filt,
        features = c("Percent_hb"),
        ncol = 1,
        group.by = 'orig.ident',
        pt.size = 0.1) 

Zhang_filt <- subset(Zhang_filt, subset =  Percent_hb < 10)






# 1) QC Moncada et al

Moncada <- readRDS("Dataframes/Raw files/Moncada/Moncada_raw_merged.rds")

Moncada

meta <- Moncada@meta.data

# dataset name + patient ID
Moncada <- PercentageFeatureSet(Moncada, pattern = "^MT-", col.name = "Percent_mito")
Moncada <- PercentageFeatureSet(Moncada, "^HB[^(P)]", col.name = "Percent_hb")
Moncada <- PercentageFeatureSet(Moncada, "^R[SP]L", col.name = "Percent_ribo")
Moncada <- PercentageFeatureSet(Moncada, "^HSP", col.name = "Percent_hsp")

count.df <- as.data.frame(table(Moncada$orig.ident))
colnames(count.df) <- c("Sample","Frequency")

ggplot(count.df, aes(x = Sample, y = Frequency, fill = Sample)) + 
  geom_col() +
  theme_classic() +
  geom_text(aes(label = Frequency), 
            position=position_dodge(width=0.9), 
            vjust=-0.25) +
  # scale_y_continuous(breaks = seq(0,5000, by = 1000), limits = c(0,5000)) +
  ggtitle("Raw: Cells per sample") +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))

VlnPlot(Moncada,
        features = c("nCount_RNA"),
        ncol = 1,
        group.by = 'orig.ident',
        pt.size = 0.1) +
  geom_hline(yintercept = 50000) + NoLegend()

VlnPlot(Moncada,
        features = c("Percent_mito"),
        ncol = 1,
        group.by = 'orig.ident',
        pt.size = 0.1) +
  geom_hline(yintercept = 30) + NoLegend()


Moncada_filt <- subset(Moncada, subset =  Percent_mito < 30 & nCount_RNA < 10000)
Moncada
Moncada_filt

FeatureScatter(Moncada, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(Moncada_filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 

saveRDS(Moncada_filt, "Dataframes/Raw files/Moncada/Moncada_raw_filt.rds")


Moncada_filt <- readRDS("Dataframes/Raw files/Moncada/Moncada_raw_filt.rds")

VlnPlot(Moncada_filt,
        features = c("Percent_hb"),
        ncol = 1,
        group.by = 'orig.ident',
        pt.size = 0.1) 

Moncada_filt <- subset(Moncada_filt, subset =  Percent_hb < 5)























