# Discovery of ecotypes

library(tidyverse)
library(RColorBrewer)
library(ComplexHeatmap)
library(cluster)
library(factoextra)
library(viridis)
library(survminer)
library(survival)

source("Consensus method/color_functions.R")
source("Consensus method/misc_functions.R")


load("Consensus method/Figure 3/Assets/tcga.bp.rdata")

tcga.theta <- bp.res@Post.ini.cs@theta

tcga.bp <- tcga.theta %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Patient.ID") %>% 
  filter(str_detect(Patient.ID, "-01A-")) %>% 
  mutate(Patient.ID = sub("^(([^-]*-){2}[^-]*)-.*$", "\\1", Patient.ID)) %>% 
  column_to_rownames(var = "Patient.ID")



tcga.theta <- bp.res@Post.ini.ct@theta

tcga.ct <- tcga.theta %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Patient.ID") %>% 
  filter(str_detect(Patient.ID, "-01A-")) %>% 
  mutate(Patient.ID = sub("^(([^-]*-){2}[^-]*)-.*$", "\\1", Patient.ID)) %>% 
  column_to_rownames(var = "Patient.ID")


tcga.gc <- read.csv("Dataframes/Bulk Expression data/TCGA/TCGA.filt.patients.csv")

tcga.gc2 <-  tcga.gc %>% 
  mutate(Patient.ID = gsub("-01A", "", Tumor.Sample.ID))

tcga.bp.filt <- tcga.bp %>% 
  filter(rownames(.) %in% tcga.gc2$Patient.ID)

ce.list <- readRDS(paste0("Consensus method/Figure 4/Assets/ce.list.rds"))
ce.list

tcga.ce.mat <- tcga.bp.filt %>%
  as.data.frame() %>% 
  select(unlist(ce.list, use.names = F)) %>% 
  rownames_to_column(var = "Patient.ID") %>% 
  pivot_longer(cols = -Patient.ID, names_to = "GP", values_to = "Percentage") %>% 
  mutate(ECO = case_when(
    GP %in% ce.list[[1]] ~ "CE1",
    GP %in% ce.list[[2]] ~ "CE2",
    GP %in% ce.list[[3]] ~ "CE3",
    GP %in% ce.list[[4]] ~ "CE4",
     T ~ "Other")) %>% 
  filter(ECO != "Other") %>% 
  group_by(Patient.ID, ECO) %>% 
   summarise(
   IQR_Mean = mean(Percentage[Percentage >= quantile(Percentage, 0.25) & 
                                  Percentage <= quantile(Percentage, 0.75)]), 
     .groups = 'drop'
   ) %>% 
  pivot_wider(id_cols = Patient.ID, names_from = ECO, values_from = IQR_Mean) %>% 
  column_to_rownames(var = "Patient.ID") %>%  scale() %>%    as.data.frame() 



load("Consensus method/Figure 3/Assets/cptac.bp.rdata")

cptac.theta <- bp.res@Post.ini.cs@theta

cptac.bp <- cptac.theta %>% 
  as.data.frame() %>% 
  select(-matches("Normal")) %>% 
  t() %>% 
  as.data.frame()


cptac.theta <- bp.res@Post.ini.ct@theta

cptac.ct <- cptac.theta %>% 
  t() %>% 
  as.data.frame()


cptac.ce.mat <- cptac.bp %>%
  as.data.frame() %>% 
  select(unlist(ce.list, use.names = F)) %>% 
  rownames_to_column(var = "Patient.ID") %>% 
  pivot_longer(cols = -Patient.ID, names_to = "GP", values_to = "Percentage") %>% 
  mutate(ECO = case_when(
    GP %in% ce.list[[1]] ~ "CE1",
    GP %in% ce.list[[2]] ~ "CE2",
    GP %in% ce.list[[3]] ~ "CE3",
    GP %in% ce.list[[4]] ~ "CE4",
    T ~ "Other")) %>% 
  filter(ECO != "Other") %>% 
  group_by(Patient.ID, ECO) %>% 
  summarise(
    IQR_Mean = mean(Percentage[Percentage >= quantile(Percentage, 0.25) & 
                                 Percentage <= quantile(Percentage, 0.75)]), 
    .groups = 'drop'
  ) %>% 
  pivot_wider(id_cols = Patient.ID, names_from = ECO, values_from = IQR_Mean) %>% 
  column_to_rownames(var = "Patient.ID") %>%  scale() %>%    as.data.frame() 

all.ce.mat <- rbind(tcga.ce.mat, cptac.ce.mat)

max_values <- apply(all.ce.mat, 1, max)  # Find the maximum value in each row

# Define a cutoff for the fold change (or difference). Example: fold change of 1.33 (i.e., top CE must be 33% greater than second max).
cutoff <- 1

assigned_CEs <- apply(all.ce.mat, 1, function(row) {
  # Find the index and value of the maximum CE
  max_index <- which.max(row)
  max_value <- as.vector(row[max_index])
  
  # Set the maximum value to 0 temporarily to find the second maximum
  row[max_index] <- 0
  
  # Find the index and value of the second maximum CE
  second_max_index <- which.max(row)
  second_max_value <- as.vector(row[second_max_index])
  
  # Calculate the fold change between the top and second top values
  fold_change <- max_value / second_max_value
  
  # Assign the top CE only if fold change exceeds the cutoff
  if (fold_change > cutoff) {
    return(names(row)[max_index])  # Assign top CE
  } else {
    return("unresolved")  # Assign as unresolved
  }
})

# Add the assigned CEs as a new column in all.ce.mat
all.ce.mat$Ecotype <- assigned_CEs

# Print the summary of the assignments
table(all.ce.mat$Ecotype)


all.ce.mat$Ecotype <- assigned_CEs

table(all.ce.mat$Ecotype)


all.ce.mat$Max_Variable <- max_values


all.ce.mat <- all.ce.mat %>% 
  filter(!Ecotype %in% c("CE5", "unresolved"))

eco.df <- all.ce.mat %>% 
  select(Ecotype,Max_Variable) %>% 
  arrange(Ecotype, desc(Max_Variable)) %>% 
  select(Ecotype) %>% 
  mutate(Cohort = case_when(
    str_detect(rownames(.), "TCGA") ~ "TCGA-PAAD",
    T ~ "CPTAC-3"
  ))


all.hm <- all.ce.mat %>%   select(-c(Ecotype, Max_Variable)) %>% t() %>% as.data.frame() %>% 
  select(rownames(eco.df))


annoheatmap1 <- HeatmapAnnotation(df = eco.df,
                                  simple_anno_size = unit(0.6, "cm"),show_legend = T,
                                  col = list(Ecotype = CE_cols), border = T, show_annotation_name = F )

tiff("Consensus method/Figure 5/Assets/Ecotypes.tiff", width = 2000, height = 800, res = 300)
ComplexHeatmap::pheatmap(all.hm,  top_annotation = annoheatmap1, show_colnames  = F,
                         cluster_rows = F,border = T,border_color = NA,row_names_side = "left",
                         column_title = "CE-based PDAC Ecotypes",
                         heatmap_legend_param = list(title = "CE fraction\nz-score", legend_direction = "vertical"),
                         cluster_cols = F, col  = viridis(n = 10 , option = "D"),
                         breaks = c(-2,0,2))
dev.off()

# Assuming your PCA data is stored in 'all.ce.mat'
# Perform PCA
pca_result <- prcomp(t(all.hm), scale. = F)

# Extract PCA scores (principal components)
pca_scores <- as.data.frame(pca_result$x)

# Merge the PCA scores with your 'eco.df' classification
# Assuming 'eco.df' has a column "Ecotype" that you want to use for coloring
pca_scores$Ecotype <- eco.df$Ecotype
pca_scores$Ecotype <- factor(pca_scores$Ecotype)
# Create a PCA plot using ggplot2 and color by 'Ecotype'
library(ggplot2)

# Plot the first two principal components (PC1 and PC2)

tiff("Consensus method/Figure 5/Assets/PCA.tiff", width = 1100, height = 800, res = 300)
ggplot(pca_scores, aes(x = PC1, y = PC2, color = Ecotype)) +
  geom_point(size = 2) +  # Points size
  labs(title = "", x = "PC1", y = "PC2") +  # Labels
  theme_bw() +  # Minimal theme for better appearance
  scale_color_manual(values = c("#E64B35FF", "#4DBBD5FF","#00A087FF", "#3C5488FF" ))  # Use colors from ECO_cols list
dev.off()


library(corrplot)

testRes = cor.mtest(t(all.hm), conf.level = 0.95)

tiff("Consensus method/Figure 5/Assets/corr.plot.tiff", width = 1100, height = 800, res = 300)
corrplot(cor(t(all.hm)), p.mat = testRes$p, method = 'color', diag = FALSE, type = 'lower',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,tl.col = "grey20",
         col = rev(COL2("RdBu")),
         insig = 'label_sig', pch.col = 'grey20')
dev.off()




# Perform pathway analysis


tme_gene_signatures <- read.csv("Consensus method/Figure 3/Dataframes/TME gene signatures.csv", check.names = F)

#tme_gene_signatures <- tme_gene_signatures %>% select(-c(EMT, `Tumour Prol`))
tme_gene_signatures <- as.list(tme_gene_signatures)

tme_gene_signatures <- lapply(tme_gene_signatures, function(x) Filter(function(y) y != "" && !is.null(y), x))
tme_gene_signatures




tcga.exp <- read.csv("Dataframes/Bulk Expression data/TCGA/tcga.paad.pt.normalized.csv", 
                     row.names = 1, check.names = F)

set.seed(1000)
module_score(tcga.exp, y = tme_gene_signatures)

tcga.tme <- signature_scores.z %>% scale() %>% as.data.frame()



cptac.exp <- read.csv("Dataframes/Bulk Expression data/CPTAC-3/cptac.pt.normalized.csv", 
                     row.names = 1, check.names = F)

set.seed(1000)
module_score(cptac.exp, y = tme_gene_signatures)

cptac.tme <- signature_scores.z %>% scale() %>% as.data.frame()


tcga.cptac.tme <- rbind(tcga.tme, cptac.tme)


eco.tme.df <- merge.all(eco.df, tcga.cptac.tme)

eco.tme.df <- eco.tme.df %>% 
  rownames_to_column(var = "Patient.ID") %>% 
  pivot_longer(cols = names(tme_gene_signatures), values_to = "ES",names_to = "GS" ) %>% 
  group_by(Ecotype, GS) %>% 
  summarise(Mean_ES = mean(ES)) %>% 
  pivot_wider(id_cols = Ecotype, names_from = GS, values_from = Mean_ES) %>% 
  column_to_rownames(var = "Ecotype") %>% scale() %>% as.data.frame()
  
eco.tme.df <- eco.tme.df %>% 
  select(Matrix, Hypoxia, Moffitt_Basal, myCAF, Complement, iCAF, Endothelium, Stellate,
         `B cells`, `Th1 signature`, `Effector cells`, `M2 signature`, `Exhausted CD8+ T`, `Type I IFN`)

Heatmap(eco.tme.df, cluster_rows = F, row_order = c("CE1", "CE2", "CE3", "CE4"), column_names_rot = 45, cluster_columns = F, row_names_side = "left")

tiff("Consensus method/Figure 5/Assets/CE.pathway.hm.tiff", width =1500, height = 1000, res = 300)
Heatmap(eco.tme.df, col = rev(brewer.pal(n = 10 , name = "RdBu")), cluster_rows = F, cluster_columns = F, row_names_side = "left",
        heatmap_legend_param = list(title  ="Enrichment\nScore"),
        column_names_rot = 45, border = T)
dev.off()


# Clinical association of Ecotypes 

# Correlate with previous annotations


tcga.gc <- read.csv("Dataframes/Bulk Expression data/TCGA/TCGA.filt.patients.csv")


tcga.gc2 <-  tcga.gc %>% 
  mutate(Patient.ID = gsub("-01A", "", Tumor.Sample.ID)) %>% 
  select(Patient.ID,Purity.Class..high.or.low., mRNA.Bailey.Clusters..All.150.Samples..1squamous.2immunogenic.3progenitor.4ADEX, mRNA.Collisson.clusters..All.150.Samples..1classical.2exocrine.3QM,
         mRNA.Moffitt.clusters..All.150.Samples..1basal..2classical, Methylation.Clusters..All.150.Samples., miRNA.Clusters..All.150.Samples., 
         Histological.type.by.DCC, Histological.type.by.RHH, History.of.chronic.pancreatitis,
         lncRNA.Clusters..All.150.Samples., Copy.Number.Clusters..All.150.Samples., Histological.type.by.DCC, Histological.type.by.RHH, History.of.chronic.pancreatitis ) %>% 
  #na.omit() %>% 
  column_to_rownames(var = "Patient.ID")


tcga.tme.gc <- merge.all(eco.df, tcga.gc2)


table(tcga.tme.gc$Histological.type.by.DCC, tcga.tme.gc$Ecotype)

tcga.tme.gc <- tcga.tme.gc %>% 
  rownames_to_column(var = "Patient.ID") %>% 
  mutate_all(as.character)

tcga.tme.gc.anno<- tcga.tme.gc %>% 
  mutate(Bailey = case_when(
    mRNA.Bailey.Clusters..All.150.Samples..1squamous.2immunogenic.3progenitor.4ADEX == "1" ~ "Squamous",
    mRNA.Bailey.Clusters..All.150.Samples..1squamous.2immunogenic.3progenitor.4ADEX == "2" ~ "Immunogenic",
    mRNA.Bailey.Clusters..All.150.Samples..1squamous.2immunogenic.3progenitor.4ADEX == "3" ~ "Progenitor",
    T ~ "ADEX")) %>% 
  mutate(Moffitt = case_when(
    mRNA.Moffitt.clusters..All.150.Samples..1basal..2classical == "1" ~ "Basal",
    mRNA.Moffitt.clusters..All.150.Samples..1basal..2classical == "2" ~ "Classical",
    T ~ "ADEX")) %>% 
  mutate(Collisson = case_when(
    mRNA.Collisson.clusters..All.150.Samples..1classical.2exocrine.3QM == "1" ~ "Classical",
    mRNA.Collisson.clusters..All.150.Samples..1classical.2exocrine.3QM == "2" ~ "Exocrine",
    T ~ "QM"))%>% 
  mutate(Purity = case_when(
    Purity.Class..high.or.low. == "high" ~ "High",
    Purity.Class..high.or.low. == "low" ~ "Low",
    T ~ "QM")) %>% 
  mutate(miRNA = case_when(
    miRNA.Clusters..All.150.Samples. == "1" ~ "miRNA-1",
    Purity.Class..high.or.low. == "2" ~ "miRNA-2",
    T ~ "miRNA-3"))%>% 
  mutate(lncRNA = case_when(
    lncRNA.Clusters..All.150.Samples. == "1" ~ "lncRNA-1",
    lncRNA.Clusters..All.150.Samples. == "2" ~ "lncRNA-2",
    T ~ "lncRNA-3"))%>% 
  mutate(Methylation = case_when(
    Methylation.Clusters..All.150.Samples. == "1" ~ "Meth-1",
    Methylation.Clusters..All.150.Samples. == "2" ~ "Meth-2",
    T ~ "Meth-3")) %>% 
  mutate(CNA = Copy.Number.Clusters..All.150.Samples.)


chisq.test(table(tcga.tme.gc.anno$Bailey, tcga.tme.gc.anno$Ecotype))

p1 <- tcga.tme.gc.anno %>% 
  dplyr::group_by(Ecotype) %>% 
  dplyr::count(Bailey) %>% 
  mutate(Percentage = (n/sum(n))*100) %>% 
  ggbarplot("Ecotype", "Percentage",
            fill = "Bailey", color = "black",
            palette ="jama",width = 0.8,
            #label = count$Percentage, 
            lab.col = "white", lab.pos = "in") +
  theme_pubr()+
  ylab("Frequency (%)")+
  ggtitle("Bailey et al\n p-value = 1.454e-09")+
  #scale_x_discrete(limits=rev)+
  xlab("")+
  #facet_grid(.~Tissue.Type, space = "free", scales = "free")+
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45,  hjust=1),
        plot.title = element_text(hjust = 0.5, size= 15))+
  labs(fill='') +
  labs(color='')


chisq.test(table(tcga.tme.gc.anno$Moffitt, tcga.tme.gc.anno$Ecotype))

p2 <- tcga.tme.gc.anno %>% 
  dplyr::group_by(Ecotype) %>% 
  dplyr::count(Moffitt) %>% 
  mutate(Percentage = (n/sum(n))*100) %>% 
  ggbarplot("Ecotype", "Percentage",
            fill = "Moffitt", color = "black",
            palette ="d3",width = 0.8,
            #label = count$Percentage, 
            lab.col = "white", lab.pos = "in") +
  theme_pubr()+
  ylab("Frequency (%)")+
  ggtitle("Moffitt et al\n p-value = 0.0007")+
  #scale_x_discrete(limits=rev)+
  xlab("")+
  #facet_grid(.~Tissue.Type, space = "free", scales = "free")+
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45,  hjust=1),
        plot.title = element_text(hjust = 0.5, size= 15))+
  labs(fill='') +
  labs(color='')

p2
chisq.test(table(tcga.tme.gc.anno$CNA, tcga.tme.gc.anno$Ecotype))

p3 <- tcga.tme.gc.anno %>% 
  dplyr::group_by(Ecotype) %>% 
  dplyr::count(CNA) %>% 
  mutate(Percentage = (n/sum(n))*100) %>% 
  ggbarplot("Ecotype", "Percentage",
            fill = "CNA", color = "black",
            palette ="simpsons",width = 0.8,
            #label = count$Percentage, 
            lab.col = "white", lab.pos = "in") +
  theme_pubr()+
  ylab("Frequency (%)")+
  ggtitle("CNA clusters\n p-value = 0.0004")+
  #scale_x_discrete(limits=rev)+
  xlab("")+
  #facet_grid(.~Tissue.Type, space = "free", scales = "free")+
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45,  hjust=1),
        plot.title = element_text(hjust = 0.5, size= 15))+
  labs(fill='') +
  labs(color='')
p3
ggarrange(p1,p2,p3, nrow = 1)


chisq.test(table(tcga.tme.gc.anno$Collisson, tcga.tme.gc.anno$Ecotype))

p3 <- tcga.tme.gc.anno %>% 
  dplyr::group_by(Ecotype) %>% 
  dplyr::count(Collisson) %>% 
  mutate(Percentage = (n/sum(n))*100) %>% 
  ggbarplot("Ecotype", "Percentage",
            fill = "Collisson", color = "black",
            palette ="simpsons",width = 0.8,
            #label = count$Percentage, 
            lab.col = "white", lab.pos = "in") +
  theme_pubr()+
  ylab("Frequency (%)")+
  ggtitle("Collisson et al\n p-value = 0.0004")+
  #scale_x_discrete(limits=rev)+
  xlab("")+
  #facet_grid(.~Tissue.Type, space = "free", scales = "free")+
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45,  hjust=1),
        plot.title = element_text(hjust = 0.5, size= 15))+
  labs(fill='') +
  labs(color='')
p3

tiff("Consensus method/Figure 5/Assets/Previous.class.tiff" , width = 3500, height = 1000, res = 300   )
ggarrange(p1,p3,p2, nrow = 1)
dev.off()



cluster_anno.surv <- eco.df %>% 
  rownames_to_column(var = "Patient.ID")

tcga.patient <- read.delim("Dataframes/Bulk Expression data/TCGA/paad_tcga/data_clinical_patient.txt", header = F)
colnames(tcga.patient) <- tcga.patient[5,]
tcga.patient <- tcga.patient[-c(1:5),]

colnames(tcga.patient)

survival.tcga <- tcga.patient %>% 
  mutate(  DFS_Status = case_when(str_detect(DFS_STATUS, "0") ~ 0,
                                  str_detect(DFS_STATUS, "1") ~ 1,
                                  T ~ NA), DFS_Months = as.numeric(DFS_MONTHS)) %>% 
  mutate(  OS_Status = case_when(str_detect(OS_STATUS, "0") ~ 0,
                                 str_detect(OS_STATUS, "1") ~ 1,
                                 T ~ NA), OS_Months = as.numeric(OS_MONTHS)) %>% 
  mutate(  PFS_Status = case_when(str_detect(PFS_STATUS, "0") ~ 0,
                                  str_detect(PFS_STATUS, "1") ~ 1,
                                  T ~ NA), PFS_Months = as.numeric(PFS_MONTHS)) %>% 
  mutate(  DSS_Status = case_when(str_detect(DSS_STATUS, "0") ~ 0,
                                  str_detect(DSS_STATUS, "1") ~ 1,
                                  T ~ NA), DSS_Months = as.numeric(DSS_MONTHS)) %>% 
  mutate(RADIATION_THERAPY = case_when(
    RADIATION_THERAPY %in% c("No") ~ "No",
    RADIATION_THERAPY %in% c("Yes") ~ "Yes",
    T ~ NA)) %>% 
  mutate(Patient.ID = PATIENT_ID) %>% 
  filter( SUBTYPE == "PAAD" & RADIATION_THERAPY == "No") %>% 
  merge(cluster_anno.surv) %>% 
  select(Patient.ID, Ecotype, OS_Status, OS_Months)



cptac.clin <- read.csv("Dataframes/Bulk Expression data/CPTAC-3/cptac.clinical.processed.csv", row.names = 1)
colnames(cptac.clin)

table(cptac.clin$histology_diagnosis)

survival.cptac <- eco.df %>% as.data.frame() %>% rownames_to_column(var = "Patient.ID") %>% 
  merge(cptac.clin, by = "Patient.ID") 

survival.cptac <- survival.cptac %>%
  #filter(histology_diagnosis == "PDAC") %>% 
  select(Patient.ID, Ecotype, OS_Status, OS_Months)


all.survival <- rbind(survival.tcga, survival.cptac)



all.os <- ggsurvplot(survfit(Surv(OS_Months,OS_Status) ~ Ecotype, data =  survival.cptac), 
           conf.int = F,          # Add confidence interval
           legend.title = "",
           pval = T,
           risk.table = F,
           legend  = "none",
           #  legend.labs = c(paste0(gene.set, "-Low", " (", prop[[1]], ")"),
           #                paste0(gene.set, "-High", " (", prop[[2]], ")")),
           pval.coord = c(20, 0.95),
           palette =  unname(ECO_cols)[1:7],
           
           # surv.median.line = "hv",
           # title = paste0(study, ", method = ",method),
           tables.theme = clean_theme(),
           data = survival.cptac,
           xlab = "Months",   # Use xlab argument instead of x in xlab() function
           ylab = "OS (%)")
dev.off()


tiff("Consensus method/Figure 5/Assets/CPTAC.OS.tiff", width = 1000, height = 900, res = 300)
all.os$plot + ggtitle("CPTAC-3")
dev.off()


p_values <- pairwise_survdiff(Surv(OS_Months,OS_Status) ~ Ecotype, data =  survival.tcga, p.adjust.method = "none")
p_values.df <- p_values$p.value %>% t() %>% as.data.frame() %>% 
  rownames_to_column(var = "x") %>% 
  pivot_longer(cols = -x, names_to = "y", values_to = "p") %>% 
  na.omit() %>% 
#  filter(p < 0.1) %>% 
  mutate(Comparison = paste0(x, " vs ", y),
         p = round(p, digits = 3)) %>% 
  select(Comparison, p) %>% 
  dplyr::rename(`p-value` = p)

p <- ggtexttable(p_values.df, rows = NULL)

tiff("Consensus method/Figure 5/Assets/TCGA.OS.comparison.tiff", width = 1200, height = 1200, res = 300)
ggarrange(p , ncol = 2, widths = c(1.8,1))
dev.off()





cluster_anno.surv <- eco.df %>% 
  rownames_to_column(var = "Patient.ID")

tcga.patient <- read.delim("Dataframes/Bulk Expression data/TCGA/paad_tcga/data_clinical_patient.txt", header = F)
colnames(tcga.patient) <- tcga.patient[5,]
tcga.patient <- tcga.patient[-c(1:5),]

colnames(tcga.patient)

survival.tcga <- tcga.patient %>% 
  mutate(  DFS_Status = case_when(str_detect(DFS_STATUS, "0") ~ 0,
                                  str_detect(DFS_STATUS, "1") ~ 1,
                                  T ~ NA), DFS_Months = as.numeric(DFS_MONTHS)) %>% 
  mutate(  OS_Status = case_when(str_detect(OS_STATUS, "0") ~ 0,
                                 str_detect(OS_STATUS, "1") ~ 1,
                                 T ~ NA), OS_Months = as.numeric(OS_MONTHS)) %>% 
  mutate(  PFS_Status = case_when(str_detect(PFS_STATUS, "0") ~ 0,
                                  str_detect(PFS_STATUS, "1") ~ 1,
                                  T ~ NA), PFS_Months = as.numeric(PFS_MONTHS)) %>% 
  mutate(  DSS_Status = case_when(str_detect(DSS_STATUS, "0") ~ 0,
                                  str_detect(DSS_STATUS, "1") ~ 1,
                                  T ~ NA), DSS_Months = as.numeric(DSS_MONTHS)) %>% 
  mutate(RADIATION_THERAPY = case_when(
    RADIATION_THERAPY %in% c("No") ~ "No",
    RADIATION_THERAPY %in% c("Yes") ~ "Yes",
    T ~ NA)) %>% 
  mutate(Patient.ID = PATIENT_ID) %>% 
 # filter(  SUBTYPE == "PAAD" ) %>% 
  merge(cluster_anno.surv)

unique(survival.tcga$SUBTYPE)


library(survival);library(survminer)
dfs <- ggsurvplot(survfit(Surv(DFS_Months,DFS_Status) ~ Ecotype, data =  survival.tcga), 
                  conf.int = F,          # Add confidence interval
                  legend.title = "",
                  pval = T,
                  risk.table = T,
                  legend  = "none",
                  #  legend.labs = c(paste0(gene.set, "-Low", " (", prop[[1]], ")"),
                  #                paste0(gene.set, "-High", " (", prop[[2]], ")")),
                  pval.coord = c(35, 0.75),
                  palette =  unname(ECO_cols)[1:7],
                  
                  # surv.median.line = "hv",
                  # title = paste0(study, ", method = ",method),
                  tables.theme = clean_theme(),
                  data = survival.tcga,
                  xlab = "Months",   # Use xlab argument instead of x in xlab() function
                  ylab = "DFS (%)")


tiff("Consensus method/Figure 5/Assets/tcga.DFS.tiff", width = 1000, height = 900, res = 300)
dfs$plot + ggtitle("TCGA-PAAD")
dev.off()



p_values <- pairwise_survdiff(Surv(DFS_Months,DFS_Status) ~ Ecotype, data =  survival.tcga, p.adjust.method = "none")
p_values.df <- p_values$p.value %>% t() %>% as.data.frame() %>% 
  rownames_to_column(var = "x") %>% 
  pivot_longer(cols = -x, names_to = "y", values_to = "p") %>% 
  na.omit() %>% 
  filter(p < 0.1) %>% 
  mutate(Comparison = paste0(x, " vs ", y),
         p = round(p, digits = 3)) %>% 
  select(Comparison, p) %>% 
  dplyr::rename(`p-value` = p)

p <- ggtexttable(p_values.df, rows = NULL)

tiff("Consensus method/Figure 5/Assets/DFS.comparison.tiff", width = 1200, height = 1200, res = 300)
ggarrange(p , ncol = 2, widths = c(1.8,1))
dev.off()



tcga.patient <- read.delim("Dataframes/Bulk Expression data/TCGA/paad_tcga/data_clinical_patient.txt", header = F)
colnames(tcga.patient) <- tcga.patient[5,]
tcga.patient <- tcga.patient[-c(1:5),]

colnames(tcga.patient)

tcga.clin.df <- tcga.patient %>% 
  mutate( OS_Status = case_when(str_detect(OS_STATUS, "0") ~ 0,
                                str_detect(OS_STATUS, "1") ~ 1,
                                T ~ NA), 
          OS_Months = as.numeric(OS_MONTHS)) %>% 
  mutate(Stage = case_when(
    AJCC_PATHOLOGIC_TUMOR_STAGE %in% c("STAGE IA", "STAGE IB", "STAGE I") ~ "1",
    AJCC_PATHOLOGIC_TUMOR_STAGE %in% c("STAGE IIA", "STAGE IIB") ~ "2",
    AJCC_PATHOLOGIC_TUMOR_STAGE %in% c("STAGE III", "STAGE IV") ~ "3&4",
    T ~ NA)) %>% 
  mutate(RADIATION_THERAPY = case_when(
    RADIATION_THERAPY %in% c("No") ~ "No",
    RADIATION_THERAPY %in% c("Yes") ~ "Yes",
    T ~ NA)) %>% 
  mutate(Patient.ID = PATIENT_ID,
         Stage = factor(Stage, levels = c("1", "2", "3&4")),
         Age = as.numeric(AGE),
         Sex = factor(SEX),
         Therapy = RADIATION_THERAPY) %>% 
  select(Patient.ID, OS_Status, OS_Months, Stage, Sex, Age,Therapy)

cptac.clin <- read.csv("Dataframes/Bulk Expression data/CPTAC-3/cptac.clinical.processed.csv", row.names = 1)
colnames(cptac.clin)

table(cptac.clin$histology_diagnosis)


cptac.clin.df <- cptac.clin %>% 
  # filter(cause_of_death == "pancreatic carcinoma") %>% 
  mutate(Histology = case_when(histology_diagnosis == "PDAC" ~ "PDAC", T ~ "Other")) %>% 
  mutate(Stage = tumor_stage_pathological) %>% 
  mutate(Stage = case_when(
    Stage %in% c("Stage I", "Stage IA", "Stage IB") ~ "1",
    Stage %in% c("Stage II", "Stage IIA", "Stage IIB") ~ "2",
    Stage %in% c("Stage III") ~ "3&4",
    Stage %in% c("Stage IV") ~ "3&4",
    T ~ NA
  )) %>% 
  mutate(Sex = sex, Age = age) %>% 
  mutate(Stage = factor(Stage, levels = c("1", "2", "3&4")),
         Histology = factor(Histology),
         Therapy = "No",
         Age = as.numeric(Age)) %>% 
  select(Patient.ID, OS_Status, OS_Months, Stage, Sex, Age,Therapy)


all.clin.df <- rbind(tcga.clin.df, cptac.clin.df)

eco.clin.df <- all.ce.mat %>% 
  rownames_to_column(var = "Patient.ID") %>% 
  merge(all.clin.df)


library(ezcox)



cox.surv.res <- ezcox(eco.clin.df, covariates = paste0("CE", 1:4),
                       controls = c("Stage", "Sex", "Age", "Therapy"),
                       time = "OS_Months", status = "OS_Status", return_models  = T)

cox.surv.res <- cox.surv.res$res
cox.surv.res <- cox.surv.res %>% filter(contrast_level %in% paste0("CE", 1:4)) %>% mutate(Log10HR = log10(HR))

write.csv(cox.surv.res, "Consensus method/Figure 5/Assets/cox.results.csv")


cox.surv.res <- cox.surv.res %>% 
  mutate(`-log10(p val)` = -log10(p.value)
         ,p.value2 = case_when(
           p.value < 0.001 ~ "***",
           p.value < 0.01 ~ "**",
           p.value < 0.05 ~ "*",
           TRUE ~ ""
         )) %>% 
  arrange(desc(beta), .by_group = T) %>% 
  mutate(Colour = case_when(
     contrast_level == "CE1" ~ CE_cols[[1]],
     contrast_level == "CE2" ~ CE_cols[[2]],
     contrast_level == "CE3" ~ CE_cols[[3]],
     T ~ CE_cols[[4]]
   )) %>% 
  arrange(desc(beta)) %>% 
  mutate(contrast_level = factor(contrast_level, levels = contrast_level))

levels(cox.surv.res$contrast_level)
print(CE_cols)

tiff("Consensus method/Figure 5/Assets/Multivariate.COX.tiff", width = 1200, height = 1000, res = 300)

ggplot(cox.surv.res, aes(HR, Variable)) +
  geom_pointrange(aes(xmin = lower_95, xmax = upper_95, fill = contrast_level),  # Use the 'Colour' column for the fill
                  shape = 22, size = 2,
                  position = position_dodge(width = 0.5)) +
  xlab("Multivariate COX log(HR)") +
  scale_x_log10() +
  geom_vline(xintercept = 1, linetype = 2) +
  theme_pubr() +
  ylab("") +
  ggtitle("O")+
  theme(text = element_text( size = 15, color = "black"),
        legend.position = "none",
        #aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5))+
  geom_text(data=cox.surv.res, aes(x=1.6, y=contrast_level,  
                                   label=p.value2), size = 10, vjust = 0.8)+
  scale_fill_manual(values = CE_cols)  # Manually map colors to levels

dev.off()

