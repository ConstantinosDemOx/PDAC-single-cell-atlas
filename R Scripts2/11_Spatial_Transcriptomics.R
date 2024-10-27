# Visium ecosystems


library(BayesPrism)
library(Seurat)
library(tidyverse)
library(InstaPrism)
library(Hmisc)
library(tidygraph)
library(ggraph)
library(corrr)
library(ggsci)
library(ggpubr)


source("Consensus method/misc_functions.R")
source("Consensus method/color_functions.R")

# read slides

ce.list <- readRDS(paste0("Consensus method/Figure 4/Assets/ce.list..rds"))
ce.list

ce.df <- stack(ce.list)

ce.df <- ce.df %>% 
  column_to_rownames(var = "values") %>% 
  dplyr::rename(CE = ind)


slide.files <- list.files(path = "Consensus method/Figure 4 Spatial/Dataframes/Visium/", pattern = "bp.rdata")

i <- slide.files[[1]]
i
cor.res <- list()
slide.bp.list <- list()


for (i in slide.files){
  
  load(paste0("Consensus method/Figure 4 Spatial/Dataframes/Visium/", i))
  
  
  slide.theta <- InstaPrism.res.initial@Post.ini.cs@theta
  dim(slide.theta)
  
  slide.bp <- slide.theta %>% t() %>% as.data.frame()
  
  slide.bp.list[[i]] <- slide.bp
  
  
  res<-rcorr(t(slide.theta), type = "pearson")
  res <- flattenCorrMatrix(res$r, res$P)
  
  res <- res %>% 
    mutate(Sample.ID = i) %>% 
    mutate(Sample.ID = gsub("prism.res.|.bp.rdata", "", Sample.ID))
  
  cor.res[[i]] <- res
  
  
    
}

cor.res.all <- do.call(rbind, cor.res)


slide.bp.all <- do.call(rbind, slide.bp.list)






cor.res.summary <- cor.res.all %>% 
  mutate(CE = case_when(
    row %in% ce.list[[1]] & column %in% ce.list[[1]] ~ "In",
    row %in% ce.list[[1]] & column %in% ce.list[[2]] ~ "Out",
    row %in% ce.list[[1]] & column %in% ce.list[[3]] ~ "Out",
    row %in% ce.list[[1]] & column %in% ce.list[[4]] ~ "Out",
  
    row %in% ce.list[[2]] & column %in% ce.list[[2]] ~ "In",
    row %in% ce.list[[2]] & column %in% ce.list[[3]] ~ "Out",
    row %in% ce.list[[2]] & column %in% ce.list[[4]] ~ "Out",
  
    row %in% ce.list[[3]] & column %in% ce.list[[3]] ~ "In",
    row %in% ce.list[[3]] & column %in% ce.list[[4]] ~ "Out",
 
    row %in% ce.list[[4]] & column %in% ce.list[[4]] ~ "In",
  
    T ~ "Other" )) %>% 
  filter(CE != "Other") %>% 
  group_by(Sample.ID, CE) %>% 
  summarise(Mean_Cor = mean(cor),
            Mean_p = mean(p)) %>% 
  mutate(Adj_scale = -log10(Mean_p + 0.0001)) %>% 
  arrange(desc(Mean_Cor))

ggpaired(cor.res.summary, x = "CE", y = "Mean_Cor",
         fill = "CE", line.color = "gray", line.size = 0.4,
         palette = "jco")+
  stat_compare_means(paired = TRUE)+
  xlab("Same CE") + ylab("Average Spearman Correlation")+
  scale_x_discrete(labels = c("Yes", "No"))


cor.res.summary <- cor.res.all %>% 
  mutate(CE = case_when(
    row %in% ce.list[[1]] & column %in% ce.list[[1]] ~ "In",
    row %in% ce.list[[1]] & column %in% ce.list[[2]] ~ "Out",
    row %in% ce.list[[1]] & column %in% ce.list[[3]] ~ "Out",
    row %in% ce.list[[1]] & column %in% ce.list[[4]] ~ "Out",
 #   row %in% ce.list[[1]] & column %in% ce.list[[5]] ~ "Out",
    T ~ "Other" )) %>% 
  filter(CE != "Other") %>% 
  group_by(Sample.ID, CE) %>% 
  summarise(Mean_Cor = mean(cor),
            Mean_p = mean(p)) %>% 
  mutate(Adj_scale = -log10(Mean_p + 0.0001)) %>% 
  arrange(desc(Mean_Cor))

p1 <- ggpaired(cor.res.summary, x = "CE", y = "Mean_Cor",
         fill = "CE", line.color = "gray", line.size = 0.4,width = 0.8,
         palette = c(CE_cols[[1]], "grey"))+
  stat_compare_means(paired = TRUE, label = "p.signif", size = 10, label.y = max(cor.res.summary$Mean_Cor), label.x = 1.5)+
  ggtitle("CE1")+
  ylim(min(cor.res.summary$Mean_Cor), max(cor.res.summary$Mean_Cor) + 0.05)+
  xlab("Same CE") + ylab("Average Pearson Correlation")+
  scale_x_discrete(labels = c("Yes", "No"))+
  theme(legend.position = "none")
p1

cor.res.summary <- cor.res.all %>% 
  mutate(CE = case_when(
    row %in% ce.list[[2]] & column %in% ce.list[[2]] ~ "In",
    row %in% ce.list[[2]] & column %in% ce.list[[1]] ~ "Out",
    row %in% ce.list[[2]] & column %in% ce.list[[3]] ~ "Out",
    row %in% ce.list[[2]] & column %in% ce.list[[4]] ~ "Out",
 #   row %in% ce.list[[2]] & column %in% ce.list[[5]] ~ "Out",
    T ~ "Other" )) %>% 
  filter(CE != "Other") %>% 
  group_by(Sample.ID, CE) %>% 
  summarise(Mean_Cor = mean(cor),
            Mean_p = mean(p)) %>% 
  mutate(Adj_scale = -log10(Mean_p + 0.0001)) %>% 
  arrange(desc(Mean_Cor))

p2 <- ggpaired(cor.res.summary, x = "CE", y = "Mean_Cor",
               fill = "CE", line.color = "gray", line.size = 0.4,width = 0.8,
               palette = c(CE_cols[[2]], "grey"))+
  stat_compare_means(paired = TRUE, label = "p.signif", size = 10, label.y = max(cor.res.summary$Mean_Cor), label.x = 1.5)+
  ggtitle("CE2")+
  xlab("Same CE") + ylab("Average Pearson Correlation")+
  ylim(min(cor.res.summary$Mean_Cor), max(cor.res.summary$Mean_Cor) + 0.05)+
  scale_x_discrete(labels = c("Yes", "No"))+
  theme(legend.position = "none")

p2



cor.res.summary <- cor.res.all %>% 
  mutate(CE = case_when(
    row %in% ce.list[[3]] & column %in% ce.list[[3]] ~ "In",
    row %in% ce.list[[3]] & column %in% ce.list[[1]] ~ "Out",
    row %in% ce.list[[3]] & column %in% ce.list[[2]] ~ "Out",
    row %in% ce.list[[3]] & column %in% ce.list[[4]] ~ "Out",
   # row %in% ce.list[[3]] & column %in% ce.list[[5]] ~ "Out",
    T ~ "Other" )) %>% 
  filter(CE != "Other") %>% 
  group_by(Sample.ID, CE) %>% 
  summarise(Mean_Cor = mean(cor),
            Mean_p = mean(p)) %>% 
  mutate(Adj_scale = -log10(Mean_p + 0.0001)) %>% 
  arrange(desc(Mean_Cor))

p3 <- ggpaired(cor.res.summary, x = "CE", y = "Mean_Cor",
               fill = "CE", line.color = "gray", line.size = 0.4,width = 0.8,
               palette = c(CE_cols[[3]], "grey"))+
  stat_compare_means(paired = TRUE, label = "p.signif", size = 10, label.y = max(cor.res.summary$Mean_Cor), label.x = 1.5)+
  ggtitle("CE3")+
  xlab("Same CE") + ylab("Average Pearson Correlation")+
  ylim(min(cor.res.summary$Mean_Cor), max(cor.res.summary$Mean_Cor) + 0.05)+
  scale_x_discrete(labels = c("Yes", "No"))+
  theme(legend.position = "none")

p3


cor.res.summary <- cor.res.all %>% 
  mutate(CE = case_when(
    row %in% ce.list[[4]] & column %in% ce.list[[4]] ~ "In",
    row %in% ce.list[[4]] & column %in% ce.list[[1]] ~ "Out",
    row %in% ce.list[[4]] & column %in% ce.list[[2]] ~ "Out",
    row %in% ce.list[[4]] & column %in% ce.list[[3]] ~ "Out",
  #  row %in% ce.list[[4]] & column %in% ce.list[[5]] ~ "Out",
    T ~ "Other" )) %>% 
  filter(CE != "Other") %>% 
  group_by(Sample.ID, CE) %>% 
  summarise(Mean_Cor = mean(cor),
            Mean_p = mean(p)) %>% 
  mutate(Adj_scale = -log10(Mean_p + 0.0001)) %>% 
  arrange(desc(Mean_Cor))

p4 <- ggpaired(cor.res.summary, x = "CE", y = "Mean_Cor",
               fill = "CE", line.color = "gray", line.size = 0.4,width = 0.8,
               palette = c(CE_cols[[4]], "grey"))+
  stat_compare_means(paired = TRUE, label = "p.signif", size = 10, label.y = max(cor.res.summary$Mean_Cor), label.x = 1.5)+
  ggtitle("CE4")+
  xlab("Same CE") + ylab("Average Pearson Correlation")+
  ylim(min(cor.res.summary$Mean_Cor), max(cor.res.summary$Mean_Cor) + 0.05)+
  scale_x_discrete(labels = c("Yes", "No"))+
  theme(legend.position = "none")

p4


p <- ggarrange(p1,p2,p3,p4, nrow = 1)


tiff("Consensus method/Figure 4/Assets/In.Out.boxplot.tiff", width = 2300, height = 1000, res = 300)
annotate_figure(p, top = text_grob("Validation of CEs in Spatial Transcriptomics", 
                                   color = "black", face = "bold", size = 14))
dev.off()




# Annotate spots



slides <- readRDS("Consensus method/Figure 4 Spatial/Dataframes/Visium/processed_seurat_files.rds")

slides_clean <- list()

for (i in names(slides)){
  slide <- slides[[i]]
  
  slide@meta.data <- slide@meta.data[,1:3]
  
  
  slides_clean[[i]] <- slide
  
  
  
}


slide.files <- list.files(path = "Consensus method/Figure 4 Spatial/Dataframes/Visium/", pattern = "bp.rdata")

slide.files <- gsub("prism.res.|.bp.rdata", "", slide.files)
slide.files


slides.anno <- list()
meta.anno <- list()


for (i in slide.files){
  
  load(paste0("Consensus method/Figure 4 Spatial/Dataframes/Visium/prism.res.", i, ".bp.rdata"))
  
  
  slide.theta <- InstaPrism.res.initial@Post.ini.cs@theta
  dim(slide.theta)
  
  slide.bp <- slide.theta %>% 
    t() %>% as.data.frame() %>% 
    rownames_to_column(var = "spot.id") %>% 
    pivot_longer(cols = -spot.id, names_to = "GP", values_to = "Percentage") %>% 
    mutate(CE = case_when(
      GP %in% ce.list[[1]]  ~ "CE1",
      GP %in% ce.list[[2]]  ~ "CE2",
      GP %in% ce.list[[3]]  ~ "CE3",
      GP %in% ce.list[[4]]  ~ "CE4",
      T ~ "Other" )) %>% 
    group_by(spot.id, CE) %>% 
   summarise(Mean_Perc = mean(Percentage)) %>% 
    pivot_wider(id_cols = spot.id, names_from = CE, values_from = Mean_Perc) %>% 
    select(-Other)
  
 
  assigned_CEs <- apply(slide.bp[,2:5], 1, function(row) {
    max_index <- which.max(row)
    
    max_index
    
    max_value <- as.vector(row[max_index])
    max_value[[1]]
    
    row[max_index] <- 0
    
    second_max_index <- which.max(row)
    second_max_value <- as.vector(row[second_max_index])
    second_max_value <- second_max_value[[1]]
    
    ratio <- second_max_value/ max_value
    ratio
    
    if (ratio < 0.85) {
      return(names(row)[max_index])
    } else {
      return("unresolved")
    }})
 
  
 slide.bp$CE <- assigned_CEs
  
  table(slide.bp$CE)
  
  
  slide.gp <- slide.theta %>% 
    t() %>% as.data.frame() %>% 
    rownames_to_column(var = "spot.id")
  
  
  SpatialDimPlot(slide, group.by = "Epi_Anno")
  
  
  tme_gene_signatures <- read.csv("Consensus method/Figure 3/Dataframes/TME gene signatures.csv", check.names = F)
  
  #tme_gene_signatures <- tme_gene_signatures %>% select(-c(EMT, `Tumour Prol`))
  tme_gene_signatures <- as.list(tme_gene_signatures)
  
  tme_gene_signatures <- lapply(tme_gene_signatures, function(x) Filter(function(y) y != "" && !is.null(y), x))
  tme_gene_signatures
  
  start <- ncol(slide@meta.data) + 1
  
  slide <- AddModuleScore(slide, features = tme_gene_signatures, name = "GS")
  
  colnames(slide@meta.data)[start:length(colnames(slide@meta.data))] <- names(tme_gene_signatures)
  
  
  
  
  slides.anno[[i]] <- slide
  meta.anno[[i]] <- slide@meta.data
  
  
}

meta.anno.sum <- do.call(rbind, meta.anno)


cor.res.top <- cor.res.all %>% 
  mutate(CE = case_when(
    row %in% ce.list[[1]] & column %in% ce.list[[1]] ~ "CE1",
    row %in% ce.list[[2]] & column %in% ce.list[[2]] ~ "CE2",
    row %in% ce.list[[3]] & column %in% ce.list[[3]] ~ "CE3",
    row %in% ce.list[[4]] & column %in% ce.list[[4]] ~ "CE4",
    #row %in% ce.list[[5]] & column %in% ce.list[[5]] ~ "CE5",
    T ~ "Other" )) %>% 
  filter(CE != "Other") %>% 
  group_by(Sample.ID, CE) %>% 
  summarise(Mean_Cor = mean(cor),
            Mean_p = mean(p)) %>% 
  mutate(Adj_scale = -log10(Mean_p + 0.0001)) %>% 
  arrange(desc(Mean_Cor), .by_group = T)



SpatialDimPlot(slides.anno[[10]], group.by = "CE", cols = CE_cols)


meta.anno.hm <- meta.anno.sum %>% 
  filter(CE != "unresolved") %>% 
  rownames_to_column(var = "spot.id") %>% 
  select(spot.id,CE, names(tme_gene_signatures)) %>% 
  pivot_longer(cols = -c(spot.id, CE), names_to = "GS", values_to = "ES") %>% 
  group_by(CE, GS) %>% 
  summarise(Mean_ES = mean(ES)) %>% 
  pivot_wider(id_cols = CE, names_from = GS, values_from = Mean_ES) %>% 
  column_to_rownames(var= "CE") %>% 
  scale() %>% as.data.frame()
colnames(meta.anno.hm)

meta.anno.hm <- meta.anno.hm %>% 
  select(Hypoxia, Matrix, Collisson_QMA,`Protumor cytokines`,
         iCAF, Endothelium, Complement, `B cells`, `Effector cells`, `Co-activation molecules`,
         `M2 signature`, `Immune Suppression by Myeloid Cells`, `Macrophage and DC traffic`)
library(ComplexHeatmap)

tiff("Consensus method/Figure 4/Assets/Spatial plots/CE.pathway.hm.tiff", width =1500, height = 1200, res = 300)
Heatmap(meta.anno.hm, col = rev(brewer.pal(n = 10 , name = "RdBu")), cluster_rows = F, cluster_columns = F, row_names_side = "left",
        heatmap_legend_param = list(title  ="Enrichment\nScore"),
        column_names_rot = 45, border = T)
dev.off()


