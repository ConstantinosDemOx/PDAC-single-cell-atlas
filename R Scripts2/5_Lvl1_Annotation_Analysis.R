# Figure 1: PDAC scRNA-Seq atlas

library(Seurat);library(tidyverse);library(ggpubr);library(wesanderson);library(scCustomize);library(ggsankey)
library(nord);library(ochRe);library(rcartocolor);library(gridExtra);library(rcartocolor);library(pheatmap);library(ComplexHeatmap)
library(rstatix);library(RColorBrewer);library(scCustomize)

source("Consensus method/color_functions.R")

pdac.atlas <- readRDS("New way/Figure_1/Files/pdac.atlas_final.rds")
pdac.atlas

# Figure 1A: Tissue pie chart


clinical.df <- pdac.atlas@meta.data %>% 
  rownames_to_column(var = "Cell.ID") %>% 
  select(Study.Patient.ID, Tissue.Type) %>% distinct() %>% 
  count(Tissue.Type) 
tiff("New way/Figure_1/Assets/tissue.number.pie.tiff", width = 1200, height = 1200, res = 300   )
pie(clinical.df$n,border = "white",
    #labels = paste0(clinical.df$Tissue.Type, " (", clinical.df$n, "," ,round(100*clinical.df$n/sum(clinical.df$n)), "%)"), 
    labels = paste0(clinical.df$Tissue.Type, ",",round(100*clinical.df$n/sum(clinical.df$n)), " %"), 
    col = tissue_cols, clockwise =T,
    main = "Tissue type")+
  theme(aspect.ratio = 1)
dev.off()




clinical.df <- pdac.atlas@meta.data %>% 
  rownames_to_column(var = "Cell.ID") %>% 
  mutate(Treatment = case_when(
    Treatment == "Folfirinox" ~ "FFX",
    Treatment == "Folfirinox + Gemcitabine + nab-paclitaxel" ~ "FFX + GP",
    Treatment == "Gemcitabine + nab-paclitaxel" ~ "GP",
    Treatment == "Chemo-RT" ~ "Chemo-Radio",
    T ~ Treatment
  )) %>% 
  select(Study.Patient.ID, Treatment) %>% distinct() %>% 
  count(Treatment)


# Figure 1A: Treatment pie chart

tiff("New way/Figure_1/Assets/treatment.number.pie.tiff", width = 1500, height = 1200, res = 300   )
pie(clinical.df$n,border = "white",
    labels = paste0(clinical.df$Treatment, ",",round(100*clinical.df$n/sum(clinical.df$n)), " %"), 
    col = treatment_cols, clockwise = T,init.angle = 50,
    main = "Neoadjuvant Therapy")+
  theme(aspect.ratio = 1)
dev.off()


# Figure 1B: Clinical features

clin.df <- read.csv("Dataframes/Clinical data/Combined.clin.csv")

tnm.df <- clin.df %>% 
  filter(TNM.Stage2 != "" ) %>% 
  count(TNM.Stage2)

p1 <- ggbarplot(tnm.df, "TNM.Stage2", "n", fill = "#4393C3", width = 0.8) +
  ylab("Number of patients") +
  xlab("TNM Stage")


sex.df <- clin.df %>% 
  filter(Gender != "" ) %>% 
  count(Gender)

p2 <- ggbarplot(sex.df, "Gender", "n", fill = "#4393C3", width = 0.8) +
  ylab("Number of patients") +
  xlab("Sex")

grade.df <- clin.df %>% 
  filter(Grade != "" ) %>% 
  count(Grade)

p3 <- ggbarplot(grade.df, "Grade", "n", fill = "#4393C3", width = 0.8) +
  ylab("Number of patients") +
  xlab("Grade")

age.df <- clin.df %>% 
  filter(!is.na(Age))

p4 <-gghistogram(age.df, x = "Age",
                 add = "mean", y = "count",
                 fill = "#4393C3")+
  ylab("Number of patients") + xlab("Age (years)")


tiff("New way/Figure_1/Assets/Summary.clin.tiff", width = 2500, height = 900, res = 300)
ggarrange(p1,p2,p3,p4, nrow = 1)
dev.off()


# Figure 1C: Lvl1 Cell type UMAP

unique(pdac.atlas$Lvl1_Clusters)

pdac.atlas@meta.data <- pdac.atlas@meta.data %>% 
  mutate(Lvl1_Clusters = case_when(
    str_detect(Lvl2_Clusters, "CD4T|CD8T|NK") ~ "T/NK",
    str_detect(Lvl2_Clusters, "Monocyte|Macrophage|DC|Neutrophil") ~ "Myeloid",
    str_detect(Lvl2_Clusters, "Ductal|Acinar|Hepatocyte") ~ "Epithelial\n(non-malignant)",
    str_detect(Lvl2_Clusters, "Malignant|Tuft") ~ "Epithelial\n(malignant)",
    TRUE ~ Lvl2_Clusters))

pdac.atlas@meta.data <- pdac.atlas@meta.data %>% 
  mutate(Lvl1_Clusters = factor(Lvl1_Clusters, levels = c(
    "Epithelial\n(malignant)", "Epithelial\n(non-malignant)", "Endocrine", "Endothelial",
    "Stellate", "Fibroblast", "Schwann", "Myeloid", "Mast", "Plasma", "B", "T/NK"
  )))



tiff("New way/Figure_1/Assets/Lvl1.UMAP.tiff",   width = 2000, height = 2000, res = 300)
DimPlot_scCustom(pdac.atlas,group.by = "Lvl1_Clusters", label = T, raster = F, pt.size = 0.05, 
                 alpha= 0.1,label.size = 5,
                 aspect_ratio = 1, 
                 colors_use = carto_pal(n = 12, name = "Vivid")) +
  ggtitle("Broad cell types\n565,584 cells") + NoLegend()
dev.off()



# Figure 1D: Lvl1 Cell type Heatmap

pdac.atlas@meta.data <- pdac.atlas@meta.data %>% 
  mutate(Lvl1_Clusters = case_when(
    str_detect(Lvl2_Clusters, "CD4T|CD8T|NK") ~ "T/NK",
    str_detect(Lvl2_Clusters, "Monocyte|Macrophage|DC|Neutrophil") ~ "Myeloid",
    str_detect(Lvl2_Clusters, "Ductal|Acinar|Hepatocyte") ~ "Epithelial (non-malignant)",
    str_detect(Lvl2_Clusters, "Malignant|Tuft") ~ "Epithelial (malignant)",
    TRUE ~ Lvl2_Clusters))

pdac.atlas@meta.data <- pdac.atlas@meta.data %>% 
  mutate(Lvl1_Clusters = factor(Lvl1_Clusters, levels = c(
    "Epithelial (malignant)", "Epithelial (non-malignant)", "Endocrine", "Endothelial",
    "Stellate", "Fibroblast", "Schwann", "Myeloid", "Mast", "Plasma", "B", "T/NK"
  )))




Idents(pdac.atlas) <- "Lvl1_Clusters"

genes <- c( "KRT19", "MUC1", "TFF1", "TFF3","KRT17", "KRT18","EPCAM",
            "SLC4A4", "CFTR", "CPB1","PRSS1", "AMY2A",
            "CPA1","CHGB","CHGA","INS",
            "VWF", "PLVAP", "FLT1", "FLT4",
            "RGS5", "PDGFRB", "ACTA2", "MYH11",
            "COL1A1", "FAP", "LUM", "COL11A1", "CFD",
            "S100B", "SOX10", "CD14", "FCGR3A", "FCN1",
            "CD163", "CD68", "CSF3R", "CPA3", "KIT",
            "JCHAIN", "MZB1","IGHG1", "IGHA2", "MS4A1", "CD79A", "CD79B", "CD83",
            "IL7R", "CCL5", "CD3E", "CD3D", "GNLY", "NKG7")

average.exp <- AverageExpression(pdac.atlas,
                                 group.by = c("Lvl1_Clusters"),
                                 assays = "RNA", slot = "data", features = genes)


average.exp <- as.matrix(average.exp$RNA)

average.scaled <- average.exp %>% 
  t() %>% scale() %>% as.data.frame() 



annot.mat <- average.scaled
annot.mat <- annot.mat %>% t() %>% as.data.frame()


max_values <- apply(annot.mat, 1, max)  # Find the maximum value in each row

annot.mat$Max_Variable <- max_values

cluster_labels <- apply(annot.mat, 1, function(x) {
  cluster_id <- which.max(x)
  cluster_label <- colnames(annot.mat)[cluster_id]
  cluster_label
})

annot.mat[,"Cluster"]<- cluster_labels

annot.mat <- annot.mat %>% 
  select(Cluster)


cell_cols <- list(Cluster = lvl1_cols)


annoheatmap1 <- HeatmapAnnotation(df = annot.mat,
                                  simple_anno_size = unit(0.6, "cm"),show_legend = F,
                                  col = cell_cols, border = T, show_annotation_name = F )


library(circlize);library(viridis)

cols <- colorRamp2(colors = rev(brewer.pal(n  = 5, name = "RdBu")),
           breaks = c(-2,-1,0,1,4))

tiff("New way/Figure_1/Assets/marker.heatmap.tiff",   width = 3200, height = 1000, res = 300)
Heatmap(as.matrix(average.scaled), row_order = rownames(average.scaled), border = T,
        column_order = colnames(average.scaled), 
        col = viridis(n = 10 , option = "D"), 
        show_row_names = T , 
        top_annotation = annoheatmap1,
        heatmap_legend_param = list(
          title = "Average \nexpression"))
dev.off()



# Figure 1E: Lvl1 Cell type Heatmap


p1 <- pdac.atlas@meta.data %>% 
  select(Study,Study.Patient.ID, Tissue.Type, Lvl1_Clusters) %>% 
  dplyr::group_by(Study) %>% 
  dplyr::count(Lvl1_Clusters) %>% 
  mutate(Percentage = (n/sum(n))*100) %>% 
  ggbarplot("Study", "Percentage",
            fill = "Lvl1_Clusters", color = "black",
            palette = carto_pal(n = 12, name = "Vivid"),
            label = F, lab.col = "white", lab.pos = "in") +
  theme_pubr()+
  ylab("Cells (%)")+
  #scale_x_discrete(limits=rev)+
  xlab("")+
  #facet_grid(.~Tissue.Type, space = "free", scales = "free")+
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45,  hjust=1),
        plot.title = element_text(hjust = 0.5, size= 15))+
  labs(fill='Cell Type') +
  labs(color='Cell Type') 



p2 <- pdac.atlas@meta.data %>% 
  select(Study,Study.Patient.ID, Tissue.Type, Lvl1_Clusters) %>%
  mutate(Tissue.Type = factor(Tissue.Type, levels= c("Normal", "Primary Tumour", 
                                                     "Metastatic"))) %>% 
  dplyr::group_by(Tissue.Type) %>% 
  dplyr::count(Lvl1_Clusters) %>% 
  mutate(Percentage = (n/sum(n))*100) %>% 
  ggbarplot("Tissue.Type", "Percentage",
            fill = "Lvl1_Clusters", color = "black",
            palette = carto_pal(n = 12, name = "Vivid"),
            label = F, lab.col = "white", lab.pos = "in") +
  theme_pubr()+
  ylab("Tissue type (%)")+
  #scale_x_discrete(limits=rev)+
  xlab("")+
  #facet_grid(.~Tissue.Type, space = "free", scales = "free")+
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45,  hjust=1),
        plot.title = element_text(hjust = 0.5, size= 15))+
  labs(fill='Cell Type') +
  labs(color='Cell Type')



p3 <- pdac.atlas@meta.data %>% 
  select(Study,Study.Patient.ID, Tissue.Type, Lvl1_Clusters, Treatment) %>%
  mutate(Treatment = case_when(
    Treatment == "Folfirinox" ~ "FFX",
    Treatment == "Folfirinox + Gemcitabine + nab-paclitaxel" ~ "FFX + GP",
    Treatment == "Gemcitabine + nab-paclitaxel" ~ "GP",
    Treatment == "Chemo-RT" ~ "Chemo-Radio",
    T ~ Treatment
  )) %>%  
  mutate(Treatment = factor(Treatment, levels = c("Naive", "GP", "FFX",
                                                  "FFX + GP", "Chemo-Radio"))) %>% 
  dplyr::group_by(Treatment) %>% 
  dplyr::count(Lvl1_Clusters) %>% 
  mutate(Percentage = (n/sum(n))*100) %>% 
  ggbarplot("Treatment", "Percentage",
            fill = "Lvl1_Clusters", color = "black",
            palette = carto_pal(n = 12, name = "Vivid"),
            label = F, lab.col = "white", lab.pos = "in") +
  theme_pubr()+
  ylab("Treatment (%)")+
  #scale_x_discrete(limits=rev)+
  xlab("")+
  #facet_grid(.~Tissue.Type, space = "free", scales = "free")+
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45,  hjust=1),
        plot.title = element_text(hjust = 0.5, size= 15))+
  labs(fill='Cell Type') +
  labs(color='Cell Type')

tiff("New way/Figure_1/Assets/bar.plots.tiff", width = 2500, height = 1200, res = 300)
ggarrange(p1,p2,p3, common.legend = T, nrow = 1,widths = c(2,1.2,1.5))
dev.off()



# Figure 1F: Lvl1 Cell type vs Tissue comparisons


tissue.lvl1.props <- pdac.atlas@meta.data %>% 
  select(Study,Study.Patient.ID, Tissue.Type, Lvl1_Clusters) %>%
  mutate(Tissue.Type = factor(Tissue.Type, levels= c("Normal", "Primary Tumour", 
                                                     "Metastatic"))) %>% 
  dplyr::group_by(Study.Patient.ID,Tissue.Type) %>% 
  dplyr::count(Lvl1_Clusters)%>% 
  mutate(Percentage = (n/sum(n))*100)

tissue.wilcox.res <- list()

for (i in unique(tissue.lvl1.props$Tissue.Type)) {
 
  df <- tissue.lvl1.props
  
  broad.cell.props.res <- df %>% 
   # filter(Tissue.Type %in% c(i, "Normal")) %>% 
    mutate(Tissue.Type = ifelse(Tissue.Type == i, i, "Rest"),
           Tissue.Type = factor(Tissue.Type, levels = c(i ,"Rest"))) %>% 
    group_by(Lvl1_Clusters) %>%
    wilcox_test(Percentage ~ Tissue.Type, paired = F) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p") %>% 
    arrange(p) %>% 
    mutate(`adj_scale` = -log10(p)) 
  
  df.fc <- df %>% 
    #filter(Tissue.Type %in% c(i, "Normal")) %>% 
    mutate(Tissue.Type = ifelse(Tissue.Type == i, i, "Rest")) %>% 
    select(Lvl1_Clusters,Tissue.Type, Percentage) %>% 
    group_by(Lvl1_Clusters,Tissue.Type) %>% 
    summarise(Mean_Perc= mean(Percentage)) %>% 
    pivot_wider(id_cols = Lvl1_Clusters, names_from = Tissue.Type, values_from = Mean_Perc) %>% 
    dplyr::rename(Group1 = i) %>% 
    mutate(LogFC = log2(Group1/Rest)) 
  
  stat.fc <- merge(broad.cell.props.res, df.fc)
  
  stat.fc <- stat.fc %>% 
    mutate(Tissue.Type = i)
  
  tissue.wilcox.res[[i]] <- stat.fc
}

tissue.wilcox.res <- do.call(rbind, tissue.wilcox.res)

p1 <- tissue.wilcox.res %>% 
  mutate(Sign = case_when(
    p < 0.05 & LogFC > 0~ "Enriched (p < 0.05)",
    p < 0.05 & LogFC < 0~ "Depleted (p < 0.05)",
    T ~ "NS"
  )) %>% 
  mutate(Sign = factor(Sign , levels = c("NS", "Depleted (p < 0.05)", "Enriched (p < 0.05)"))) %>% 
  ggplot(aes(x = factor(Tissue.Type, levels = rev(c("Normal", "Primary Tumour", "Metastatic"))), 
             y = factor(Lvl1_Clusters, levels = levels(pdac.atlas)), 
             size = adj_scale, fill = Sign)) + 
  scale_size(range = c(3, 6))+
  #scale_fill_distiller(palette = "RdBu", direction = -1)+
 # scale_color_viridis_d()+
  scale_fill_viridis_d()+
  geom_point(shape = 21) +
  #ggtitle("Tissue type") +
  ylab("") + 
  xlab("") + 
  theme_bw()+
  labs(size = "-log10(p adj)")+
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.position = "right")+
  coord_flip()

p1
write.csv(tissue.wilcox.res, "New way/Figure_1/Assets/tissue.lvl1.wilcoxon.test.csv")



# Figure 1F: Lvl1 Cell type vs Treatment comparisons

treatment.lvl1.props <- pdac.atlas@meta.data %>% 
  select(Study,Study.Patient.ID, Treatment, Lvl1_Clusters) %>%
  mutate(Treatment = case_when(
    Treatment == "Folfirinox" ~ "Neoadjuvant CTx",
    Treatment == "Folfirinox + Gemcitabine + nab-paclitaxel" ~ "Neoadjuvant CTx",
    Treatment == "Gemcitabine + nab-paclitaxel" ~ "Neoadjuvant CTx",
    Treatment == "Chemo-RT" ~ "Neoadjuvant CTx",
    T ~ Treatment
  )) %>%  
  mutate(Treatment = factor(Treatment, levels = c("Naive", "Neoadjuvant CTx"))) %>% 
  dplyr::group_by(Study.Patient.ID,Treatment) %>% 
  dplyr::count(Lvl1_Clusters)%>% 
  mutate(Percentage = (n/sum(n))*100)

treatment.wilcox.res <- list()


Idents(pdac.atlas) <- "Lvl1_Clusters"


for (i in c( "Neoadjuvant CTx")) {
  
  df <- treatment.lvl1.props
  
  broad.cell.props.res <- df %>% 
    filter(Treatment %in% c(i, "Naive")) %>% 
    mutate(Treatment = ifelse(Treatment == i, i, "Rest"),
           Treatment = factor(Treatment, levels = c(i ,"Rest"))) %>% 
    group_by(Lvl1_Clusters) %>%
    wilcox_test(Percentage ~ Treatment, paired = F) %>%
    adjust_pvalue(method = "fdr") %>%
    add_significance("p") %>% 
    arrange(p) %>% 
    mutate(`adj_scale` = -log10(p)) 
  
  df.fc <- df %>% 
    filter(Treatment %in% c(i, "Naive")) %>% 
    mutate(Treatment = ifelse(Treatment == i, i, "Rest")) %>% 
    select(Lvl1_Clusters,Treatment, Percentage) %>% 
    group_by(Lvl1_Clusters,Treatment) %>% 
    summarise(Mean_Perc= mean(Percentage)) %>% 
    pivot_wider(id_cols = Lvl1_Clusters, names_from = Treatment, values_from = Mean_Perc) %>% 
    dplyr::rename(Group1 = i) %>% 
    mutate(LogFC = log2(Group1/Rest)) 
  
  stat.fc <- merge(broad.cell.props.res, df.fc)
  
  stat.fc <- stat.fc %>% 
    mutate(Treatment = i)
  
  treatment.wilcox.res[[i]] <- stat.fc
}

treatment.wilcox.res <- do.call(rbind, treatment.wilcox.res)


p2 <- treatment.wilcox.res %>% 
  mutate(Sign = case_when(
    p < 0.05 & LogFC > 0~ "Enriched (p < 0.05)",
    p < 0.05 & LogFC < 0~ "Depleted (p < 0.05)",
    T ~ "NS"
  )) %>% 
  mutate(Sign = factor(Sign , levels = c("NS", "Depleted (p < 0.05)", "Enriched (p < 0.05)"))) %>% 
  ggplot(aes(x = factor(Treatment, levels = c("Neoadjuvant CTx")), 
             y = factor(Lvl1_Clusters, levels = levels(pdac.atlas)), 
             size = adj_scale, fill = Sign)) + 
 scale_size(range = c(3, 6))+
# ggtitle("Treatment vs Naive") +
  geom_point(shape =21) +
  scale_fill_viridis_d()+
  ylab("") + 
  xlab("") + 
  theme_bw()+
  labs(size = "-log10(p adj)")+
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        #axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.position = "right") +
   coord_flip()


p2

tiff("New way/Figure_1/Assets/Treatment.sign.tiff", width = 3500, height = 800, res = 300)
ggarrange(p1,p2, common.legend = T, legend = "top", widths = c(1,1))
dev.off()



write.csv(treatment.wilcox.res, "New way/Figure_1/Assets/treatment.lvl1.wilcoxon.test.csv")
