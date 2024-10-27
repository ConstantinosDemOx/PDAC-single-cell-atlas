# Tissue comparisons of cell states

library(Seurat);library(tidyverse);library(ggpubr);library(circlize);library(ComplexHeatmap);library(factoextra);library(cluster);library(factoextra)
library(tidygraph);library(igraph);library(ggraph);library(corrr);library(viridis)

source("Consensus method/misc_functions.R")
source("Consensus method/color_functions.R")

all_gps

pdac.atlas <- readRDS("New way/Figure_1/Files/pdac.atlas_final.rds")

epi.gps <- readRDS("Consensus method/Figure 2/Assets/Epithelial/epi_gps.rds")
epi.gps
epi.cells <- readRDS("Consensus method/Figure 2/Assets/Epithelial/cell.mat.split.rds")
epi.cells <- epi.cells[1:6]
names(epi.cells) <- paste0("Epi_", names(epi.cells))

pdac.atlas@meta.data <- pdac.atlas@meta.data %>% 
  mutate(GP = case_when(
    rownames(.) %in% epi.cells[[1]] ~ "Epi_GP1",
    rownames(.) %in% epi.cells[[2]] ~ "Epi_GP2",
    rownames(.) %in% epi.cells[[3]] ~ "Epi_GP3",
    rownames(.) %in% epi.cells[[4]] ~ "Epi_GP4",
    rownames(.) %in% epi.cells[[5]] ~ "Epi_GP5",
    rownames(.) %in% epi.cells[[6]] ~ "Epi_GP6",
    T ~ Lvl2_Clusters))

table(pdac.atlas$Lvl1_Clusters, pdac.atlas$GP)


fib.gps <- readRDS("Consensus method/Figure 2/Assets/Fibroblast/fib_gps.rds")
fib.gps
fib.cells <- readRDS("Consensus method/Figure 2/Assets/Fibroblast/cell.mat.split.rds")

fib.cells <- fib.cells[1:8]
names(fib.cells) <- paste0("Fib_", names(fib.cells))


pdac.atlas@meta.data <- pdac.atlas@meta.data %>% 
  mutate(GP = case_when(
    rownames(.) %in% fib.cells[[1]] ~ "Fib_GP1",
    rownames(.) %in% fib.cells[[2]] ~ "Fib_GP2",
    rownames(.) %in% fib.cells[[3]] ~ "Fib_GP3",
    rownames(.) %in% fib.cells[[4]] ~ "Fib_GP4",
    rownames(.) %in% fib.cells[[5]] ~ "Fib_GP5",
    rownames(.) %in% fib.cells[[6]] ~ "Fib_GP6",
    rownames(.) %in% fib.cells[[7]] ~ "Fib_GP7",
    rownames(.) %in% fib.cells[[8]] ~ "Fib_GP8",
    T ~ GP))

table(pdac.atlas$Lvl1_Clusters, pdac.atlas$GP)


ecs.gps <- readRDS("Consensus method/Figure 2/Assets/Endothelial/ecs_gps.rds")
ecs.gps
ecs.cells <- readRDS("Consensus method/Figure 2/Assets/Endothelial/cell.mat.split.rds")
ecs.cells <- ecs.cells[1:4]
names(ecs.cells) <- paste0("EC_", names(ecs.cells))


pdac.atlas@meta.data <- pdac.atlas@meta.data %>% 
  mutate(GP = case_when(
    rownames(.) %in% ecs.cells[[1]] ~ "EC_GP1",
    rownames(.) %in% ecs.cells[[2]] ~ "EC_GP2",
    rownames(.) %in% ecs.cells[[3]] ~ "EC_GP3",
    rownames(.) %in% ecs.cells[[4]] ~ "EC_GP4",
    T ~ GP))

table(pdac.atlas$Lvl1_Clusters, pdac.atlas$GP)




cd4t.cells <- readRDS("Consensus method/Figure 2/Assets/CD4T/cell.mat.split.rds")
cd4t.cells <- cd4t.cells[1:7]
names(cd4t.cells) <- paste0("CD4T_", names(cd4t.cells))


pdac.atlas@meta.data <- pdac.atlas@meta.data %>% 
  mutate(GP = case_when(
    rownames(.) %in% cd4t.cells[[1]] ~ "CD4T_GP1",
    rownames(.) %in% cd4t.cells[[2]] ~ "CD4T_GP2",
    rownames(.) %in% cd4t.cells[[3]] ~ "CD4T_GP3",
    rownames(.) %in% cd4t.cells[[4]] ~ "CD4T_GP4",
    rownames(.) %in% cd4t.cells[[5]] ~ "CD4T_GP5",
    rownames(.) %in% cd4t.cells[[6]] ~ "CD4T_GP6",
    rownames(.) %in% cd4t.cells[[7]] ~ "CD4T_GP7",
    T ~ GP))

table(pdac.atlas$Lvl1_Clusters, pdac.atlas$GP)



cd8t.gps <- readRDS("Consensus method/Figure 2/Assets/CD8T/CD8T_gps.rds")
cd8t.gps
cd8t.cells <- readRDS("Consensus method/Figure 2/Assets/CD8T/cell.mat.split.rds")
cd8t.cells <- cd8t.cells[1:5]
names(cd8t.cells) <- paste0("CD8T_", names(cd8t.cells))


pdac.atlas@meta.data <- pdac.atlas@meta.data %>% 
  mutate(GP = case_when(
    rownames(.) %in% cd8t.cells[[1]] ~ "CD8T_GP1",
    rownames(.) %in% cd8t.cells[[2]] ~"CD8T_GP2",
    rownames(.) %in% cd8t.cells[[3]] ~ "CD8T_GP3",
    rownames(.) %in% cd8t.cells[[4]] ~ "CD8T_GP4",
    rownames(.) %in% cd8t.cells[[5]] ~ "CD8T_GP5",
    T ~ GP))




nk.cells <- readRDS("Consensus method/Figure 2/Assets/NK/cell.mat.split.rds")
nk.cells <- nk.cells[1:3]
names(nk.cells) <- paste0("NK_", names(nk.cells))


pdac.atlas@meta.data <- pdac.atlas@meta.data %>% 
  mutate(GP = case_when(
    rownames(.) %in% nk.cells[[1]] ~ "NK_GP1",
    rownames(.) %in% nk.cells[[2]] ~ "NK_GP2",
    rownames(.) %in% nk.cells[[3]] ~ "NK_GP3",
    T ~ GP))


b.gps <- readRDS("Consensus method/Figure 2/Assets/B/B.Pl_gps.rds")
b.gps
B.cells <- readRDS("Consensus method/Figure 2/Assets/B/cell.mat.split.rds")
B.cells <- B.cells[1:6]
names(B.cells) <- paste0("B_", names(B.cells))


pdac.atlas@meta.data <- pdac.atlas@meta.data %>% 
  mutate(GP = case_when(
    rownames(.) %in% B.cells[[1]] ~ "B_GP1",
    rownames(.) %in% B.cells[[2]] ~ "B_GP2",
    rownames(.) %in% B.cells[[3]] ~ "B_GP3",
    rownames(.) %in% B.cells[[4]] ~ "B_GP4",
    rownames(.) %in% B.cells[[5]] ~ "B_GP5",
    rownames(.) %in% B.cells[[6]] ~ "B_GP6",
    T ~ GP))


MacMono_gps <- readRDS("Consensus method/Figure 2/Assets/MacMono/MacMono_gps.rds")
MacMono_gps
MacMono.cells <- readRDS("Consensus method/Figure 2/Assets/MacMono/cell.mat.split.rds")
MacMono.cells <- MacMono.cells[1:7]
names(MacMono.cells) <- paste0("MacMono_", names(MacMono.cells))


pdac.atlas@meta.data <- pdac.atlas@meta.data %>% 
  mutate(GP = case_when(
    rownames(.) %in% MacMono.cells[[1]] ~ "MacMono_GP1",
    rownames(.) %in% MacMono.cells[[2]] ~ "MacMono_GP2",
    rownames(.) %in% MacMono.cells[[3]] ~ "MacMono_GP3",
    rownames(.) %in% MacMono.cells[[4]] ~ "MacMono_GP4",
    rownames(.) %in% MacMono.cells[[5]] ~ "MacMono_GP5",
    rownames(.) %in% MacMono.cells[[6]] ~ "MacMono_GP6",
    rownames(.) %in% MacMono.cells[[7]] ~ "MacMono_GP7",
    T ~ GP))



Neut.cells <- readRDS("Consensus method/Figure 2/Assets/Neutrophil/cell.mat.split.rds")
Neut.cells <- Neut.cells[1:3]
names(Neut.cells) <- paste0("Neut_", names(Neut.cells))


pdac.atlas@meta.data <- pdac.atlas@meta.data %>% 
  mutate(GP = case_when(
    rownames(.) %in% Neut.cells[[1]] ~ "Neut_GP1",
    rownames(.) %in% Neut.cells[[2]] ~ "Neut_GP2",
    rownames(.) %in% Neut.cells[[3]] ~ "Neut_GP3",
    T ~ GP))



DCs.cells <- readRDS("Consensus method/Figure 2/Assets/DCs/cell.mat.split.rds")
DCs.cells <- DCs.cells[1:4]
names(DCs.cells) <- paste0("DCs_", names(DCs.cells))


pdac.atlas@meta.data <- pdac.atlas@meta.data %>% 
  mutate(GP = case_when(
    rownames(.) %in% DCs.cells[[1]] ~ "DC_GP1",
    rownames(.) %in% DCs.cells[[2]] ~ "DC_GP2",
    rownames(.) %in% DCs.cells[[3]] ~ "DC_GP3",
    rownames(.) %in% DCs.cells[[4]] ~ "DC_GP4",
    T ~ GP))


Mast.cells <- readRDS("Consensus method/Figure 2/Assets/Mast/cell.mat.split.rds")
Mast.cells <- Mast.cells[1:2]
names(Mast.cells) <- paste0("Mast_", names(Mast.cells))


pdac.atlas@meta.data <- pdac.atlas@meta.data %>% 
  mutate(GP = case_when(
    rownames(.) %in% Mast.cells[[1]] ~ "Mast_GP1",
    rownames(.) %in% Mast.cells[[2]] ~ "Mast_GP2",
    T ~ GP))


Plasma.cells <- readRDS("Consensus method/Figure 2/Assets/Plasma/cell.mat.split.rds")
Plasma.cells <- Plasma.cells[1:4]
names(Plasma.cells) <- paste0("Plasma_", names(Plasma.cells))


pdac.atlas@meta.data <- pdac.atlas@meta.data %>% 
  mutate(GP = case_when(
    rownames(.) %in% Plasma.cells[[1]] ~ "Plasma_GP1",
    rownames(.) %in% Plasma.cells[[2]] ~ "Plasma_GP2",
    rownames(.) %in% Plasma.cells[[3]] ~ "Plasma_GP3",
    rownames(.) %in% Plasma.cells[[4]] ~ "Plasma_GP4",
    T ~ GP))

dc_xcr1_cells <- WhichCells(pdac.atlas, expression = XCR1 > 0.5)


pdac.atlas@meta.data <- pdac.atlas@meta.data %>% 
  mutate(GP = case_when(
    rownames(.) %in% dc_xcr1_cells ~ "DC_GP5",
    T ~ GP))



library(rstatix)


tissue.props <- pdac.atlas@meta.data %>% 
  select(Study,Study.Patient.ID, Tissue.Type, GP) %>%
  filter(!str_detect(GP, "Plasma|Mast")) %>%
  filter(str_detect(GP, "_")) %>%
  separate(col = GP, sep = "_", into = c("CS", "GP2"), remove = F) %>% 
  mutate(Tissue.Type = factor(Tissue.Type, levels= c("Normal", "Primary Tumour", 
                                                     "Metastatic"))) %>% 
  dplyr::group_by(Study.Patient.ID,Tissue.Type,CS) %>% 
  dplyr::count(GP)%>% 
  mutate(Percentage = (n/sum(n))*100) 

write.csv(tissue.props, "Consensus method/Figure 2/Assets/tissue GP props.csv")


tissue.props.count <- table(tissue.props$Tissue.Type, tissue.props$GP) %>% as.data.frame() %>% 
  pivot_wider(id_cols = Var2, names_from = Var1, values_from = Freq) %>% 
  filter(Normal >= 3 & `Primary Tumour` >= 3 & Metastatic >= 3)


library(rstatix)
tissue.wilcox.res <- list()

for (i in unique(tissue.props$Tissue.Type)) {
  
  df <- tissue.props
  
  broad.cell.props.res <- df %>% 
    filter(GP %in% tissue.props.count$Var2) %>% 
    # filter(Tissue.Type %in% c(i, "Normal")) %>% 
    mutate(Tissue.Type = ifelse(Tissue.Type == i, i, "Rest"),
           Tissue.Type = factor(Tissue.Type, levels = c(i ,"Rest"))) %>% 
    group_by(GP) %>%
    wilcox_test(Percentage ~ Tissue.Type, paired = F) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p") %>% 
    arrange(p) %>% 
    mutate(`adj_scale` = -log10(p)) 
  
  df.fc <- df %>% 
    filter(GP %in% tissue.props.count$Var2) %>% 
    mutate(Tissue.Type = ifelse(Tissue.Type == i, i, "Rest")) %>% 
    select(GP,Tissue.Type, Percentage) %>% 
    group_by(GP,Tissue.Type) %>% 
    summarise(Mean_Perc= mean(Percentage)) %>% 
    pivot_wider(id_cols = GP, names_from = Tissue.Type, values_from = Mean_Perc) %>% 
    dplyr::rename(Group1 = i) %>% 
    mutate(LogFC = log2(Group1/Rest)) 
  
  stat.fc <- merge(broad.cell.props.res, df.fc)
  
  stat.fc <- stat.fc %>% 
    mutate(Tissue.Type = i)
  
  tissue.wilcox.res[[i]] <- stat.fc
}

tissue.wilcox.res <- do.call(rbind, tissue.wilcox.res)

write.csv(tissue.wilcox.res, "Consensus method/Figure 2/Assets/tissue wilcox GP results.csv")

met.wilcox <- tissue.wilcox.res %>% 
  filter(LogFC > 0 & p < 0.05 & Tissue.Type == "Metastatic") %>% 
  arrange(desc(LogFC)) %>% 
  separate(GP, sep = "_", into = c("CT", "GP2"), remove = F) %>% ungroup()

tiff("Consensus method/Figure 2/Assets/met.barplot.tiff", width = 1300, height = 1800, res = 300)
ggbarplot(met.wilcox, x = "GP", y = "adj_scale", sort.val = "asc", fill = "#8DA0CB")+
  geom_hline(yintercept = -log10(0.05), linetype = 2)+
  coord_flip()+
  ylab("-log10 (p-value)") +
  xlab("")+
  ggtitle("Enriched in metastatic tissue\np < 0.05")+
  theme(plot.title = element_text(hjust = 0))
dev.off()


ggbarplot(met.wilcox, x = "GP", y = "adj_scale", sort.val = "asc")+
  geom_hline(yintercept = -log10(0.05), linetype = 2)+
  # coord_flip()+
  ylab("-log10 (p-value)") +
  xlab("")+
  facet_grid(.~CT, space = "free", scales = "free")+
  ggtitle("Enriched in metastatic tissue")



normal.wilcox <- tissue.wilcox.res %>% 
  filter(LogFC > 0 & p < 0.05 & Tissue.Type == "Normal") %>% 
  arrange(desc(LogFC))

tiff("Consensus method/Figure 2/Assets/normal.barplot.tiff", width = 1300, height = 1800, res = 300)
ggbarplot(normal.wilcox, x = "GP", y = "adj_scale", sort.val = "asc", fill = "#66C2A5")+
  geom_hline(yintercept = -log10(0.05), linetype = 2)+
  coord_flip()+
  ylab("-log10 (p-value)") +
  xlab("")+
  ggtitle("Enriched in normal tissue\np < 0.05")+
  theme(plot.title = element_text(hjust = 0))
dev.off()



pt.wilcox <- tissue.wilcox.res %>% 
  filter(LogFC > 0 & p < 0.05 & Tissue.Type == "Primary Tumour") %>% 
  arrange(desc(LogFC))

tiff("Consensus method/Figure 2/Assets/normal.barplot.tiff", width = 1300, height = 1800, res = 300)
ggbarplot(normal.wilcox, x = "GP", y = "adj_scale", sort.val = "asc", fill = "#66C2A5")+
  geom_hline(yintercept = -log10(0.05), linetype = 2)+
  coord_flip()+
  ylab("-log10 (p-value)") +
  xlab("")+
  ggtitle("Enriched in normal tissue\np < 0.05")+
  theme(plot.title = element_text(hjust = 0))
dev.off()





treatment.gp.props <-pdac.atlas@meta.data %>% 
  select(Study,Study.Patient.ID, Treatment, GP) %>%
  filter(!str_detect(GP, "Plasma|Mast")) %>%
  filter(str_detect(GP, "_")) %>%
  separate(col = GP, sep = "_", into = c("CS", "GP2"), remove = F) %>% 
  mutate(Treatment = case_when(
    Treatment == "Folfirinox" ~ "Neoadjuvant CTx",
    Treatment == "Folfirinox + Gemcitabine + nab-paclitaxel" ~ "Neoadjuvant CTx",
    Treatment == "Gemcitabine + nab-paclitaxel" ~ "Neoadjuvant CTx",
    Treatment == "Chemo-RT" ~ "Neoadjuvant CTx",
    T ~ Treatment
  )) %>%  
  mutate(Treatment = factor(Treatment, levels = c("Naive", "Neoadjuvant CTx"))) %>% 
  dplyr::group_by(Study.Patient.ID,Treatment,CS) %>% 
  dplyr::count(GP)%>% 
  mutate(Percentage = (n/sum(n))*100)

treatment.wilcox.res <- list()


for (i in unique(treatment.gp.props$CS)) {
  
  df <- treatment.gp.props
  
  broad.cell.props.res <- df %>% 
    filter(CS %in% i) %>% 
    group_by(GP) %>%
    wilcox_test(Percentage ~ Treatment, paired = F) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p") %>% 
    arrange(p) %>% 
    mutate(`adj_scale` = -log10(p))  %>% 
    mutate(CS = i)
  
  df.fc <- df %>% 
    filter(CS %in% i) %>% 
    select(GP,Treatment, Percentage) %>% 
    group_by(GP,Treatment) %>% 
    summarise(Mean_Perc= mean(Percentage)) %>% 
    pivot_wider(id_cols = GP, names_from = Treatment, values_from = Mean_Perc) %>% 
    mutate(LogFC = log2(`Neoadjuvant CTx`/`Naive`)) 
  
  stat.fc <- merge(broad.cell.props.res, df.fc)
  
  treatment.wilcox.res[[i]] <- stat.fc
}

treatment.wilcox.res <- do.call(rbind, treatment.wilcox.res)



neo.wilcox <- treatment.wilcox.res %>% 
  filter(LogFC > 0 & p < 0.05 ) %>% 
  arrange(desc(LogFC)) %>% 
  separate(GP, sep = "_", into = c("CT", "GP2"), remove = F) %>% ungroup()

tiff("Consensus method/Figure 2/Assets/ctx.barplot.tiff", width = 1300, height = 1000, res = 300)
ggbarplot(neo.wilcox, x = "GP", y = "adj_scale", sort.val = "asc", fill = "#4DBBD5FF")+
  geom_hline(yintercept = -log10(0.05), linetype = 2)+
  coord_flip()+
  ylab("-log10 (p-value)") +
  xlab("")+
  ggtitle("Enriched in neo CTx samples\np < 0.05")+
  theme(plot.title = element_text(hjust = 0))
dev.off()



treatment.gp.props <-pdac.atlas@meta.data %>% 
  select(Study,Study.Patient.ID, Treatment, GP) %>%
  filter(str_detect(GP, "_")) %>%
  separate(col = GP, sep = "_", into = c("CS", "GP")) %>% 
  mutate(Treatment = case_when(
    Treatment == "Folfirinox" ~ "FFX and/or GP",
    Treatment == "Folfirinox + Gemcitabine + nab-paclitaxel" ~ "FFX and/or GP",
    Treatment == "Gemcitabine + nab-paclitaxel" ~ "FFX and/or GP",
    Treatment == "Chemo-RT" ~ "Chemo-Radio",
    T ~ Treatment
  )) %>%  
  mutate(Treatment = factor(Treatment, levels = c("Naive", "FFX and/or GP", "Chemo-Radio"))) %>% 
  dplyr::group_by(Study,Study.Patient.ID,Treatment,CS) %>% 
  dplyr::count(GP)%>% 
  mutate(Percentage = (n/sum(n))*100)


treatment.gp.samples.props <-pdac.atlas@meta.data %>% 
  select(Study,Study.Sample.ID, Treatment, GP, Study.Patient.ID) %>%
  filter(str_detect(GP, "_")) %>%
  separate(col = GP, sep = "_", into = c("CS", "GP")) %>% 
  mutate(Treatment = case_when(
    Treatment == "Folfirinox" ~ "FFX and/or GP",
    Treatment == "Folfirinox + Gemcitabine + nab-paclitaxel" ~ "FFX and/or GP",
    Treatment == "Gemcitabine + nab-paclitaxel" ~ "FFX and/or GP",
    Treatment == "Chemo-RT" ~ "Chemo-Radio",
    T ~ Treatment
  )) %>%  
  mutate(Treatment = factor(Treatment, levels = c("Naive", "FFX and/or GP", "Chemo-Radio"))) %>% 
  dplyr::group_by(Study,Study.Sample.ID,Treatment,CS) %>% 
  dplyr::count(GP)%>% 
  mutate(Percentage = (n/sum(n))*100)



p1 <- treatment.gp.props %>% 
  mutate(GP2 = paste0(CS, "_", GP)) %>% 
  #  filter(Study %in% c("Zhou et al")) %>% 
  filter(GP2 == "Fib_GP2") %>% 
  ggboxplot("Treatment", "Percentage", fill = "Treatment", palette = treatment_cols,width = 0.9) +
  ggtitle("Zhou et al and Werba et al\npatient-based") +
  stat_compare_means(comparisons = list(c("Naive", "FFX and/or GP")),
                     hide.ns = T)+
  #  facet_grid(.~Study) +
  ylab("Fib GP2 Fraction (%)")+
  xlab("") + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1),
                   plot.title = element_text(size = 10)) +
  ylim(0,100)


p2 <- treatment.gp.samples.props %>% 
  mutate(GP2 = paste0(CS, "_", GP)) %>% 
  filter(Study %in% c("Zhou et al")) %>% 
  filter(GP2 == "Fib_GP2") %>% 
  ggboxplot("Treatment", "Percentage", fill = "Treatment", palette = treatment_cols,width = 0.9) +
  ggtitle("Zhou et al\nsample-based") +
  stat_compare_means(comparisons = list(c("Naive", "FFX and/or GP")),
                     hide.ns = T)+
  #  facet_grid(.~Study) +
  ylab("Fib GP2 Fraction (%)")+
  xlab("") + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1),
                   plot.title = element_text(size = 10)) +
  ylim(0,100)





tiff("Consensus method/Figure 2/Assets/Fib.GP2.Fraction.tiff", width =1200, height = 1200, res = 300)
ggarrange(p1,p2)
dev.off()





p1 <- treatment.gp.props %>% 
  mutate(GP2 = paste0(CS, "_", GP)) %>% 
  #  filter(Study %in% c("Zhou et al")) %>% 
  filter(GP2 == "NK_GP1") %>% 
  ggboxplot("Treatment", "Percentage", fill = "Treatment", palette = treatment_cols,width = 0.8) +
  ggtitle("Zhou et al and Werba et al\npatient-based") +
  stat_compare_means(comparisons = list(c("Naive", "FFX and/or GP"), c("Naive", "Chemo-Radio")),
                     hide.ns = T)+
  #  facet_grid(.~Study) +
  ylab("NK GP1 Fraction (%)")+
  xlab("") + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1),
                   plot.title = element_text(size = 10)) +
  ylim(0,112)


p2 <- treatment.gp.samples.props %>% 
  mutate(GP2 = paste0(CS, "_", GP)) %>% 
  filter(Study %in% c("Zhou et al")) %>% 
  filter(GP2 == "NK_GP1") %>% 
  ggboxplot("Treatment", "Percentage", fill = "Treatment", palette = treatment_cols,width = 0.8) +
  ggtitle("Zhou et al\nsample-based") +
  stat_compare_means(comparisons = list(c("Naive", "FFX and/or GP"), c("Naive", "Chemo-Radio")),
                     hide.ns = T)+
  #  facet_grid(.~Study) +
  ylab("NK GP1 Fraction (%)")+
  xlab("") + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1),
                   plot.title = element_text(size = 10)) +
  ylim(0,112)



tiff("Consensus method/Figure 2/Assets/NK.GP1.Fraction.tiff", width =1200, height = 1200, res = 300)
ggarrange(p1,p2)
dev.off()



