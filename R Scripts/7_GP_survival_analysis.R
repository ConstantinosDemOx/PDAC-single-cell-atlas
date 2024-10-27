library(tidyverse)
library(ggsci)
library(survival)
library(ezcox)
library(ggpubr)
library(patchwork)


source("Consensus method/color_functions.R")
source("Consensus method/misc_functions.R")


cox.tcga <- function(x){
  tcga.exp <- read.csv("Dataframes/Bulk Expression data/TCGA/tcga.paad.pt.normalized.csv", row.names = 1, check.names = F )
  tcga.gq <- read.csv("Dataframes/Bulk Expression data/TCGA/TCGA.good.quality.csv")
  tcga.gq <- tcga.gq %>% 
    mutate(Tumor.Sample.ID = gsub("-01A", "", Tumor.Sample.ID))
  
  tcga.exp <- tcga.exp %>% 
    select(intersect(tcga.gq$Tumor.Sample.ID, colnames(tcga.exp)))
  
  
  set.seed(1000)
  module_score(tcga.exp, x)
  
  
  tcga.patient <- read.delim("Dataframes/Bulk Expression data/TCGA/paad_tcga/data_clinical_patient.txt", header = F)
  colnames(tcga.patient) <- tcga.patient[5,]
  tcga.patient <- tcga.patient[-c(1:5),]
  
  colnames(tcga.patient)
  
  signature_scores.z <- signature_scores.z %>% as.data.frame() %>% rownames_to_column(var = "Patient.ID")
  
  multi.cox.df <- tcga.patient %>% 
    mutate( OS_Status = case_when(str_detect(OS_STATUS, "0") ~ 0,
                                  str_detect(OS_STATUS, "1") ~ 1,
                                  T ~ NA), 
            OS_Months = as.numeric(OS_MONTHS)) %>% 
    mutate(Stage = case_when(
      AJCC_PATHOLOGIC_TUMOR_STAGE %in% c("STAGE IA", "STAGE IB", "STAGE I") ~ "1",
      AJCC_PATHOLOGIC_TUMOR_STAGE %in% c("STAGE IIA", "STAGE IIB") ~ "1",
      AJCC_PATHOLOGIC_TUMOR_STAGE %in% c("STAGE III", "STAGE IV") ~ "2",
      T ~ NA)) %>% 
    mutate(RADIATION_THERAPY = case_when(
      RADIATION_THERAPY %in% c("No") ~ "No",
      RADIATION_THERAPY %in% c("Yes") ~ "Yes",
      T ~ NA)) %>% 
    filter(SUBTYPE == "PAAD") %>%
    # mutate(Age = case_when(AGE > me))
    mutate(Patient.ID = PATIENT_ID,
           Stage = factor(Stage, levels = c("1", "2")),
           Age = as.numeric(AGE),
           Sex = factor(SEX),
           Therapy = RADIATION_THERAPY) %>% 
    merge(signature_scores.z, by = "Patient.ID") %>% 
    select(Patient.ID, OS_Status, OS_Months, Stage, Sex, Age, Therapy,names(x)) %>% 
    mutate(Study = "TCGA-PAAD")#%>% 
  # filter(OS_Months > 1) 
  
  assign("TCGA.clin.df", multi.cox.df, envir=.GlobalEnv)
}
cox.cptac <- function(x){
  cptac.exp <- read.csv("Dataframes/Bulk Expression data/CPTAC-3/cptac.pt.normalized.csv", row.names = 1, check.names = F)
  
  set.seed(1000)
  module_score(cptac.exp, x)
  
  cptac.clin <- read.csv("Dataframes/Bulk Expression data/CPTAC-3/cptac.clinical.processed.csv", row.names = 1)
  colnames(cptac.clin)
  
  table(cptac.clin$histology_diagnosis)
  
  cptac.tme.os <- signature_scores.z %>% as.data.frame() %>% rownames_to_column(var = "Patient.ID") %>% 
    merge(cptac.clin, by = "Patient.ID") %>% 
    filter(histology_diagnosis == "PDAC") #%>% 
  #filter(OS_Months > 1) 
  
  unique(cptac.clin$sex)
  
  multi.cox.df <- cptac.tme.os %>% 
    # filter(cause_of_death == "pancreatic carcinoma") %>% 
    mutate(Histology = case_when(histology_diagnosis == "PDAC" ~ "PDAC", T ~ "Other")) %>% 
    mutate(Stage = tumor_stage_pathological) %>% 
    mutate(Stage = case_when(
      Stage %in% c("Stage I", "Stage IA", "Stage IB") ~ "1",
      Stage %in% c("Stage II", "Stage IIA", "Stage IIB") ~ "1",
      Stage %in% c("Stage III") ~ "2",
      Stage %in% c("Stage IV") ~ "2",
      T ~ NA
    )) %>% 
    mutate(Sex = sex, Age = age) %>% 
    mutate(Stage = factor(Stage, levels = c("1", "2")),
           Histology = factor(Histology),
           Therapy = "No",
           Age = as.numeric(Age)) %>% 
    select(Patient.ID, OS_Status, OS_Months, Stage, Sex, Age,Therapy, names(x))%>% 
    mutate(Study = "CPTAC-3")#%>% 
  # filter(OS_Months > 1)
  
  assign("CPTAC.clin.df", multi.cox.df, envir=.GlobalEnv)
  
  
}


all_gps

cox.tcga(all_gps)
cox.cptac(all_gps)

all.surv.df <- rbind(TCGA.clin.df, CPTAC.clin.df)


library(ezcox)



tcga.surv.res <- ezcox(TCGA.clin.df, covariates = names(all_gps),
                       controls = c("Stage", "Sex", "Age", "Therapy"),
                       time = "OS_Months", status = "OS_Status", return_models  = T)

tcga.surv.res <- tcga.surv.res$res
tcga.surv.hm <- tcga.surv.res %>% filter(contrast_level %in% names(all_gps)) 

write.csv(tcga.surv.hm, "Consensus method/Figure 3/Assets/tcga.cox.results.csv")


tcga.surv.hm <- tcga.surv.hm %>% 
  mutate(Adj_scale = -log10(p.value)) %>% 
  mutate(Adj_scale = case_when(beta > 0 ~ Adj_scale,
                               T ~ -Adj_scale)) %>% 
  select(Variable, Adj_scale) %>% 
  column_to_rownames(var = "Variable") %>% 
  arrange(Adj_scale) %>% 
  dplyr::rename(`Discovery\n(TCGA-PAAD)` = Adj_scale) %>% 
  filter(!str_detect(rownames(.), "STMN1|Mast|Plasma")) 



tiff("Consensus method/Figure 3/Assets/TCGA.cox.tiff", width = 1200, height = 2800, res = 300)
Heatmap(as.matrix(tcga.surv.hm),row_names_side = "left", border = T, show_row_dend = F, column_names_rot = 45,
        heatmap_legend_param = list(
          title = "Overall survival time", at = c(min(tcga.surv.hm$`Discovery\n(TCGA-PAAD)`),-1.3,0,1.3,max(tcga.surv.hm$`Discovery\n(TCGA-PAAD)`) ), 
          labels = c("0","p = 0.05", "1", "p = 0.05", "0")))
dev.off()



cptac.surv.res <- ezcox(CPTAC.clin.df, covariates = names(all_gps),
                        controls = c("Stage", "Sex", "Age"),
                        time = "OS_Months", status = "OS_Status", return_models  = T)

cptac.surv.res <- cptac.surv.res$res
cptac.surv.hm <- cptac.surv.res %>% filter(contrast_level %in%  names(all_gps))




cptac.surv.hm <- cptac.surv.hm %>% 
  mutate(Adj_scale = -log10(p.value)) %>% 
  mutate(Adj_scale = case_when(beta > 0 ~ Adj_scale,
                               T ~ -Adj_scale)) %>% 
  select(Variable, Adj_scale) %>% 
  column_to_rownames(var = "Variable") %>% 
  t() %>% as.data.frame() %>% select(rev(rownames(tcga.surv.hm))) %>% t() %>% as.data.frame() %>% 
  dplyr::rename(`CPTAC-3` = Adj_scale)%>% 
  filter(!str_detect(rownames(.), "STMN1|Mast|Plasma")) 

tiff("Consensus method/Figure 3/Assets/CPTAC.cox.tiff", width = 750, height = 2800, res = 300)
Heatmap(as.matrix(cptac.surv.hm),row_names_side = "left", border = T, show_row_dend = F, column_names_rot = 45,
        cluster_rows = F,show_row_names = F,
        heatmap_legend_param = list(
          title = "Overall survival time", at = c(min(cptac.surv.hm$`CPTAC-3`),-1.3,0,1.3,max(cptac.surv.hm$`CPTAC-3`) ), 
          labels = c("0","p = 0.05", "1", "p = 0.05", "0")))
dev.off()



concordance.df <- merge.all(tcga.surv.hm, cptac.surv.hm)

concordance.df <- concordance.df %>% 
  #  filter(!str_detect(rownames(.), "STMN1")) %>% 
  t() %>% as.data.frame() %>% select(rev(rownames(tcga.surv.hm))) %>% t() %>% as.data.frame() 


rev(brewer.pal(n = 10 , name = "RdBu"))

tiff("Consensus method/Figure 3/Assets/GP.COX.tiff", width = 1400, height = 3000, res = 300)
Heatmap(as.matrix(concordance.df),row_names_side = "left", border = T, show_row_dend = F, column_names_rot = 45,
        cluster_rows = F,show_row_names = T,cluster_columns = F, column_order = c("Discovery\n(TCGA-PAAD)", "CPTAC-3"),
        col = c("#2166AC", "white", "#B2182B"),
        heatmap_legend_param = list(
          title = "Overall survival time", at = c(min(cptac.surv.hm$`CPTAC-3`),-1.3,0,1.3,max(cptac.surv.hm$`CPTAC-3`) ), 
          labels = c("0","p = 0.05", "1", "p = 0.05", "0")))
dev.off()




concordance.summary.df <- concordance.df %>% 
  rownames_to_column(var = "GP") %>% 
  mutate(GP = factor(GP ,levels = rownames(concordance.df))) %>% 
  pivot_longer(cols = c(`Discovery\n(TCGA-PAAD)`, `CPTAC-3`), names_to = "Cohort", values_to = "Adj_Scale") %>% 
  group_by(GP) %>% 
  summarise(Mean_Adj_Scale = mean(Adj_Scale, na.rm = T)) %>% 
  mutate(Group = case_when(
    Mean_Adj_Scale > 0 ~ "Shorter",
    T ~ "Longer"
  )) %>% column_to_rownames(var = "GP") %>% 
  arrange(desc(Mean_Adj_Scale))

group.cols <- c("#2166AC", "#B2182B")
names(group.cols) <- c("Longer", "Shorter")

row_annotation <- HeatmapAnnotation(
  `-log10(p-value)` = anno_barplot(concordance.summary.df$Mean_Adj_Scale, 
                                   gp = gpar(fill = group.cols[concordance.summary.df$Group]),width =  unit(25, "mm"),
  ),which = "row",
  show_legend = T,
  
  annotation_name_rot = 0
)

concodance.hm <- concordance.df %>% t() %>% as.data.frame() %>% select(rownames(concordance.summary.df)) %>% t() %>% 
  as.data.frame()

library(circlize)

color_fun <- colorRamp2(breaks = c(min(concodance.hm), 0 , max(concodance.hm)), colors =c("#2166AC", "white", "#B2182B") )

tiff("Consensus method/Figure 3/Assets/GP.COX.tiff", width = 1600, height = 3000, res = 300)
Heatmap(as.matrix(concodance.hm),row_names_side = "left", border = T, show_row_dend = F, column_names_rot = 45,
        right_annotation = row_annotation,
        # row_order = rownames(concordance.summary.df),
        cluster_rows = F,show_row_names = T,cluster_columns = F, column_order = c("Discovery\n(TCGA-PAAD)", "CPTAC-3"),
        col = color_fun,
        column_title = "R = 0.67, p = 1.5e-07",
        column_title_gp = gpar(fontsize =20),
        heatmap_legend_param = list(
          title = "Overall\nsurvival time", at = c(min(cptac.surv.hm$`CPTAC-3`),-1.3,0,1.3,max(cptac.surv.hm$`CPTAC-3`) ), 
          labels = c("0","p = 0.05", "1", "p = 0.05", "0")))

decorate_annotation("-log10(p-value)", {
  grid.lines(x = unit(c(-1, -1), "native"), y = unit(c(0, 1), "npc"), gp = gpar(lty = 2, col = "black"))
  grid.lines(x = unit(c(1, 1), "native"), y = unit(c(0, 1), "npc"), gp = gpar(lty = 2, col = "black"))
})
dev.off()




concordance.df <- data.frame(concordance.df)
colnames(concordance.df)

concordance.df <- concordance.df %>% 
  mutate(`Cell type` = case_when(
    str_detect(rownames(concordance.df), "Epi") ~ "Epithelial",
    str_detect(rownames(concordance.df), "Fib") ~ "Fibroblast",
    str_detect(rownames(concordance.df), "MacMono") ~ "MacMono",
    str_detect(rownames(concordance.df), "Neut") ~ "Neutrophil",
    str_detect(rownames(concordance.df), "CD4T") ~ "CD4+ T",
    str_detect(rownames(concordance.df), "CD8T") ~ "CD8+ T",
    str_detect(rownames(concordance.df), "DC") ~ "DCs",
    str_detect(rownames(concordance.df), "NK_") ~ "NK",
    str_detect(rownames(concordance.df), "B_") ~ "B",
    str_detect(rownames(concordance.df), "EC_") ~ "Endothelial",
    str_detect(rownames(concordance.df), "Plasma_") ~ "Plasma",
    str_detect(rownames(concordance.df), "Mast_") ~ "Mast",
    T ~ "Other"
    
  ))



tiff("Consensus method/Figure 3/Assets/Concordance.plot.tiff", width = 1800, height = 1500, res = 300)
ggscatter(concordance.df, x = "CPTAC.3", y = "Discovery..TCGA.PAAD.",
          color = "black", shape = 21, size = 3, repel = T, fill = "Cell type",
          palette = ct_cols2,
          # label = rownames(concordance.df),
          # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson",  label.sep = "\n", size = 8)
) +theme(legend.position = "right")+
  xlab("CPTAC-3\nLonger Survival         -log10(p value)        Shorter Survival")+
  ylab("TCGA-PAAD\nLonger Survival         -log10(p value)        Shorter Survival")
dev.off()


head.gps <- head(rownames(concodance.hm), n = 10)
head.gps


tail.gps <- tail(rownames(concodance.hm), n = 10)
tail.gps


# Pathway enrichment

progeny.results <- read.csv("Consensus method/Figure 3/Progeny.results.csv", row.names = 1)


progeny.results.filt <- progeny.results %>% 
  mutate(Group = case_when(
    GP %in% c("Fibroblast_GP4","Neutrophil_GP3", "Epithelial_GP1","B_GP4",  "MacMono_GP4","CD4T_GP5", "Fibroblast_GP1",
              "MacMono_GP5", "Neutrophil_GP2", "B_GP2") ~ "Worse",
    GP %in% c( "CD8T_GP1", "CD8T_GP5","B_GP3","Epithelial_GP2", 
               "Endothelial_GP3", "CD8T_GP2", "DCs_GP1", "Endothelial_GP4", "DCs_GP4", "Fibroblast_GP2") ~ "Better",
    T ~ "Other"
  )) %>% 
  filter(Group != "Other") 

library(rstatix)

sign.pathways <- progeny.results.filt %>% 
  group_by(source) %>% 
  wilcox_test(mean ~ Group) %>% 
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p") %>% 
  filter(p < 0.1) %>% 
  arrange(p)

progeny.results.hm <- progeny.results.filt %>% 
  select(source, GP, mean) %>% 
  pivot_wider(id_cols = source, names_from = GP, values_from = mean) %>% 
  filter(source %in% sign.pathways$source) %>% 
  column_to_rownames(var = "source")

progeny.results.hm <- progeny.results.hm %>% 
  select(Fibroblast_GP4, Neutrophil_GP3, Epithelial_GP1, B_GP4, MacMono_GP4,CD4T_GP5, Fibroblast_GP1, MacMono_GP5, 
         Neutrophil_GP2,B_GP2,CD8T_GP1, CD8T_GP5,
         B_GP3,Epithelial_GP2, Endothelial_GP3, CD8T_GP2, DCs_GP1, Endothelial_GP4, Fibroblast_GP2, DCs_GP4)
head
tail.gps

colnames(progeny.results.hm) <- c(head.gps, tail.gps)

progeny.results.hm <- progeny.results.hm %>% 
  t() %>% as.data.frame() %>% select(sign.pathways$source) %>% scale()


top_anno <- concordance.summary.df %>% 
  filter(rownames(.) %in% rownames(progeny.results.hm)) %>% 
  dplyr::rename(`-log10(p-value)` = Mean_Adj_Scale) %>% 
  select(`-log10(p-value)`,Group)


group.cols <- c("#2166AC", "#B2182B")
names(group.cols) <- c("Longer", "Shorter")



anno_fun <- colorRamp2(breaks = c(min(top_anno$`-log10(p-value)`), 0 , max(top_anno$`-log10(p-value)`)), 
                       colors =c("#2166AC", "white", "#B2182B") )



column_annotation <- HeatmapAnnotation(
  df = top_anno, border = T,annotation_legend_param = 
    list(legend_direction = "horizontal"),
  col = list(Group = group.cols,
             `-log10(p-value)` = anno_fun),
  show_legend = T,
  annotation_name_rot = 0
)


color_fun <- colorRamp2(breaks = c(min(progeny.results.hm), 0 , 2), colors =viridis(n = 10 , option = "D")[c(1,5,10)] )

colnames(progeny.results.hm) <- paste0(sign.pathways$source, " p = ", round(sign.pathways$p, digits =3))


tiff("Consensus method/Figure 3/Assets/Pathway.groups.tiff", width = 2000, height = 1200, res = 300)
draw(Heatmap(t(progeny.results.hm), cluster_columns = F, top_annotation = column_annotation,
             cluster_rows = F,heatmap_legend_param = list(title = "Average\nEnrichment",legend_direction = "horizontal"),
             border = T,
             column_names_rot = 90, col = color_fun),  heatmap_legend_side = "bottom",annotation_legend_side = "bottom", merge_legend=T)
dev.off()


library(decoupleR)


net <- get_progeny(organism = 'human', top = 500)
net

tcga.exp <- read.csv("Dataframes/Bulk Expression data/TCGA/tcga.paad.pt.normalized.csv", row.names = 1, check.names = F)


sample_acts <- run_mlm(mat=tcga.exp, net=net, .source='source', .target='target',
                       .mor='weight', minsize = 5)
sample_acts

sample_acts_mat <- sample_acts %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()



progeny_km_survival_tcga <- function(gp, method, study){
  library(survminer);library(survival)
  tcga.gp <- sample_acts_mat %>% 
    scale() %>% 
    as.data.frame() %>% 
    select(gp) %>% 
    rownames_to_column(var = "Patient.ID") %>% 
    pivot_longer(cols = -Patient.ID, names_to = "gp", values_to = "Expression")
  
  
  tcga.patient <- read.delim("Dataframes/Bulk Expression data/TCGA/paad_tcga/data_clinical_patient.txt", header = F)
  colnames(tcga.patient) <- tcga.patient[5,]
  tcga.patient <- tcga.patient[-c(1:5),]
  
  colnames(tcga.patient)
  
  surv.df <- tcga.patient %>% 
    filter(SUBTYPE == "PAAD") %>% 
    dplyr::rename(Patient.ID = PATIENT_ID) %>% 
    mutate( OS_Status = case_when(str_detect(OS_STATUS, "0") ~ 0,
                                  str_detect(OS_STATUS, "1") ~ 1,
                                  T ~ NA), 
            OS_Months = as.numeric(OS_MONTHS)) %>% 
    select(Patient.ID, OS_Status, OS_Months) %>% 
    merge(tcga.gp)
  
  print(head(surv.df))
  
  
  if (method == "maxstat") {
    
    surv.df <- surv_cutpoint(
      surv.df,
      time = "OS_Months",
      event = "OS_Status",
      variables = c("Patient.ID","Expression","gp"))
    
    print(summary(surv.df))
    
    surv.df <- surv_categorize(surv.df)
    
    print(surv.df)
    
    surv.df <- data.frame(Patient.ID = as.character(surv.df$Patient.ID),
                          OS_Months = as.numeric(surv.df$OS_Months),
                          OS_Status = as.numeric(surv.df$OS_Status),
                          Expression = as.factor(surv.df$Expression),
                          gp_Name = as.character(surv.df$gp))
    
    print(surv.df)
    
    # surv.df <- surv.df %>% rename(!!gp := "setNames.as.factor.surv.df..gp.....gp.")
    
    print(surv.df)
    
    
  } else if (method == "median") {
    
    surv.df$Expression <- ifelse(surv.df$Expression >= median(surv.df$Expression), "high", "low")
  } else if (method == "quartiles") {
    
    q <- quantile(surv.df$Expression, c(0.25, 0.5, 0.75))
    
    surv.df$Expression <- ifelse(surv.df$Expression <= q[1], "low", ifelse(surv.df$Expression >= q[3], "high", "other"))
    
    surv.df <- surv.df %>% filter(surv.df$Expression %in% c("low", "high")) %>% 
      select(Patient.ID,OS_Status, OS_Months, Expression, gp)
  } else {
    stop("Invalid method selected. Please choose 'maxstat', 'median', or 'quartiles'.")
  }
  
  tiff(paste0("Consensus method/Figure 3/Assets/", gp,".",study, ".tiff"), width = 1200, height = 1200, res = 300)
  print(ggsurvplot(survfit(Surv(OS_Months, OS_Status) ~ Expression, data =  surv.df), 
                   conf.int = F,          # Add confidence interval
                   legend.title = "",
                   pval = T,
                   risk.table = F,
                   legend.labs = c(paste0(gp, "-High"),paste0(gp, "-Low") ),
                   pval.coord = c(50, 0.9),
                   palette = c("#35B779FF", "#440154FF"),
                   # surv.median.line = "hv",
                   title = paste0("Discovery (", study, ")"),
                   tables.theme = clean_theme(),
                   data = surv.df,
                   xlab = "Months",   # Use xlab argument instead of x in xlab() function
                   ylab = "Overall Survival"))
  dev.off()
  
}

progeny_km_survival_tcga(gp = "TGFb", method = "maxstat", study = "TCGA-PAAD")



cptac.exp <- read.csv("Dataframes/Bulk Expression data/CPTAC-3/cptac.pt.normalized.csv", row.names = 1, check.names = F)


sample_acts <- run_mlm(mat=cptac.exp, net=net, .source='source', .target='target',
                       .mor='weight', minsize = 5)
sample_acts

sample_acts_mat <- sample_acts %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


progeny_km_survival_cptac <- function(gp, method, study){
  
  cptac.gp <- sample_acts_mat %>% 
    scale() %>% 
    as.data.frame() %>% 
    select(gp) %>% 
    rownames_to_column(var = "Patient.ID") %>% 
    pivot_longer(cols = -Patient.ID, names_to = "gp", values_to = "Expression")
  
  
  
  cptac.patient <- read.csv("Dataframes/Bulk Expression data/CPTAC-3/cptac.clinical.processed.csv", row.names = 1)
  colnames(cptac.patient)
  
  
  surv.df <- cptac.patient %>% 
    filter(histology_diagnosis == "PDAC") %>% 
    select(Patient.ID, OS_Status, OS_Months) %>% 
    merge(cptac.gp)
  
  print(head(surv.df))
  
  
  if (method == "maxstat") {
    
    surv.df <- surv_cutpoint(
      surv.df,
      time = "OS_Months",
      event = "OS_Status",
      variables = c("Patient.ID","Expression","gp"))
    
    print(summary(surv.df))
    
    surv.df <- surv_categorize(surv.df)
    
    print(surv.df)
    
    surv.df <- data.frame(Patient.ID = as.character(surv.df$Patient.ID),
                          OS_Months = as.numeric(surv.df$OS_Months),
                          OS_Status = as.numeric(surv.df$OS_Status),
                          Expression = as.factor(surv.df$Expression),
                          gp_Name = as.character(surv.df$gp))
    
    print(surv.df)
    
    # surv.df <- surv.df %>% rename(!!gp := "setNames.as.factor.surv.df..gp.....gp.")
    
    print(surv.df)
    
    
  } else if (method == "median") {
    
    surv.df$Expression <- ifelse(surv.df$Expression >= median(surv.df$Expression), "high", "low")
  } else if (method == "quartiles") {
    
    q <- quantile(surv.df$Expression, c(0.25, 0.5, 0.75))
    
    surv.df$Expression <- ifelse(surv.df$Expression <= q[1], "low", ifelse(surv.df$Expression >= q[3], "high", "other"))
    
    surv.df <- surv.df %>% filter(surv.df$Expression %in% c("low", "high")) %>% 
      select(Patient.ID,OS_Status, OS_Months, Expression, gp)
  } else {
    stop("Invalid method selected. Please choose 'maxstat', 'median', or 'quartiles'.")
  }
  
  tiff(paste0("Consensus method/Figure 3/Assets/", gp,".",study, ".tiff"), width = 1200, height = 1200, res = 300)
  print(ggsurvplot(survfit(Surv(OS_Months, OS_Status) ~ Expression, data =  surv.df), 
                   conf.int = F,          # Add confidence interval
                   legend.title = "",
                   pval = T,
                   risk.table = F,
                   legend.labs = c(paste0(gp, "-High"),paste0(gp, "-Low") ),
                   pval.coord = c(25, 0.9),
                   palette = c("#35B779FF", "#440154FF"),
                   # surv.median.line = "hv",
                   title = paste0("Validation (", study, ")"),
                   tables.theme = clean_theme(),
                   data = surv.df,
                   xlab = "Months",   # Use xlab argument instead of x in xlab() function
                   ylab = "Overall Survival"))
  dev.off()
  
}

progeny_km_survival_cptac(gp = "Trail", method = "maxstat", study = "CPTAC-3")





tcga <- read.csv("Dataframes/Bulk Expression data/TCGA/tcga.paad.pt.normalized.csv", row.names = 1, check.names = F)

set.seed(1000)
module_score(tcga, all_gps)

tcga.gp <- signature_scores.z %>% scale() %>% as.data.frame()

GP_km_survival_tcga <- function(gp, method, study){
  library(survminer);library(survival)
  tcga.gp <- tcga.gp %>% 
    select(gp) %>% 
    rownames_to_column(var = "Patient.ID") %>% 
    pivot_longer(cols = -Patient.ID, names_to = "gp", values_to = "Expression")
  
  
  tcga.patient <- read.delim("Dataframes/Bulk Expression data/TCGA/paad_tcga/data_clinical_patient.txt", header = F)
  colnames(tcga.patient) <- tcga.patient[5,]
  tcga.patient <- tcga.patient[-c(1:5),]
  
  colnames(tcga.patient)
  
  surv.df <- tcga.patient %>% 
    filter(SUBTYPE == "PAAD") %>% 
    dplyr::rename(Patient.ID = PATIENT_ID) %>% 
    mutate( OS_Status = case_when(str_detect(OS_STATUS, "0") ~ 0,
                                  str_detect(OS_STATUS, "1") ~ 1,
                                  T ~ NA), 
            OS_Months = as.numeric(OS_MONTHS)) %>% 
    select(Patient.ID, OS_Status, OS_Months) %>% 
    merge(tcga.gp)
  
  print(head(surv.df))
  
  
  if (method == "maxstat") {
    
    surv.df <- surv_cutpoint(
      surv.df,
      time = "OS_Months",
      event = "OS_Status",
      variables = c("Patient.ID","Expression","gp"))
    
    print(summary(surv.df))
    
    surv.df <- surv_categorize(surv.df)
    
    print(surv.df)
    
    surv.df <- data.frame(Patient.ID = as.character(surv.df$Patient.ID),
                          OS_Months = as.numeric(surv.df$OS_Months),
                          OS_Status = as.numeric(surv.df$OS_Status),
                          Expression = as.factor(surv.df$Expression),
                          gp_Name = as.character(surv.df$gp))
    
    print(surv.df)
    
    # surv.df <- surv.df %>% rename(!!gp := "setNames.as.factor.surv.df..gp.....gp.")
    
    print(surv.df)
    
    
  } else if (method == "median") {
    
    surv.df$Expression <- ifelse(surv.df$Expression >= median(surv.df$Expression), "high", "low")
  } else if (method == "quartiles") {
    
    q <- quantile(surv.df$Expression, c(0.25, 0.5, 0.75))
    
    surv.df$Expression <- ifelse(surv.df$Expression <= q[1], "low", ifelse(surv.df$Expression >= q[3], "high", "other"))
    
    surv.df <- surv.df %>% filter(surv.df$Expression %in% c("low", "high")) %>% 
      select(Patient.ID,OS_Status, OS_Months, Expression, gp)
  } else {
    stop("Invalid method selected. Please choose 'maxstat', 'median', or 'quartiles'.")
  }
  
  tiff(paste0("Consensus method/Figure 3/Assets/", gp,".",study, ".tiff"), width = 1200, height = 1200, res = 300)
  print(ggsurvplot(survfit(Surv(OS_Months, OS_Status) ~ Expression, data =  surv.df), 
                   conf.int = F,          # Add confidence interval
                   legend.title = "",
                   pval = T,
                   risk.table = F,
                   legend.labs = c(paste0(gp, "-High"),paste0(gp, "-Low") ),
                   pval.coord = c(50, 0.9),
                   palette = c("#B2182B", "#2166AC"),
                   # surv.median.line = "hv",
                   title = paste0("Discovery (", study, ")"),
                   tables.theme = clean_theme(),
                   data = surv.df,
                   xlab = "Months",   # Use xlab argument instead of x in xlab() function
                   ylab = "Overall Survival"))
  dev.off()
  
}

GP_km_survival_tcga(gp = "B_CD83", method = "maxstat", study = "TCGA-PAAD")

cptac <- read.csv("Dataframes/Bulk Expression data/CPTAC-3/cptac.pt.normalized.csv", row.names = 1, check.names = F)

set.seed(1000)
module_score(cptac, all_gps)

cptac.gp <- signature_scores.z %>% scale() %>% as.data.frame()

GP_km_survival_cptac <- function(gp, method, study){
  
  cptac.gp <- cptac.gp %>% 
    select(gp) %>% 
    rownames_to_column(var = "Patient.ID") %>% 
    pivot_longer(cols = -Patient.ID, names_to = "gp", values_to = "Expression")
  
  
  
  cptac.patient <- read.csv("Dataframes/Bulk Expression data/CPTAC-3/cptac.clinical.processed.csv", row.names = 1)
  colnames(cptac.patient)
  
  
  surv.df <- cptac.patient %>% 
    filter(histology_diagnosis == "PDAC") %>% 
    select(Patient.ID, OS_Status, OS_Months) %>% 
    merge(cptac.gp)
  
  print(head(surv.df))
  
  
  if (method == "maxstat") {
    
    surv.df <- surv_cutpoint(
      surv.df,
      time = "OS_Months",
      event = "OS_Status",
      variables = c("Patient.ID","Expression","gp"))
    
    print(summary(surv.df))
    
    surv.df <- surv_categorize(surv.df)
    
    print(surv.df)
    
    surv.df <- data.frame(Patient.ID = as.character(surv.df$Patient.ID),
                          OS_Months = as.numeric(surv.df$OS_Months),
                          OS_Status = as.numeric(surv.df$OS_Status),
                          Expression = as.factor(surv.df$Expression),
                          gp_Name = as.character(surv.df$gp))
    
    print(surv.df)
    
    # surv.df <- surv.df %>% rename(!!gp := "setNames.as.factor.surv.df..gp.....gp.")
    
    print(surv.df)
    
    
  } else if (method == "median") {
    
    surv.df$Expression <- ifelse(surv.df$Expression >= median(surv.df$Expression), "high", "low")
  } else if (method == "quartiles") {
    
    q <- quantile(surv.df$Expression, c(0.25, 0.5, 0.75))
    
    surv.df$Expression <- ifelse(surv.df$Expression <= q[1], "low", ifelse(surv.df$Expression >= q[3], "high", "other"))
    
    surv.df <- surv.df %>% filter(surv.df$Expression %in% c("low", "high")) %>% 
      select(Patient.ID,OS_Status, OS_Months, Expression, gp)
  } else {
    stop("Invalid method selected. Please choose 'maxstat', 'median', or 'quartiles'.")
  }
  
  tiff(paste0("Consensus method/Figure 3/Assets/", gp,".",study, ".tiff"), width = 1200, height = 1200, res = 300)
  print(ggsurvplot(survfit(Surv(OS_Months, OS_Status) ~ Expression, data =  surv.df), 
                   conf.int = F,          # Add confidence interval
                   legend.title = "",
                   pval = T,
                   risk.table = F,
                   legend.labs = c(paste0(gp, "-High"),paste0(gp, "-Low") ),
                   pval.coord = c(25, 0.9),
                   palette = c("#B2182B", "#2166AC"),
                   # surv.median.line = "hv",
                   title = paste0("Validation (", study, ")"),
                   tables.theme = clean_theme(),
                   data = surv.df,
                   xlab = "Months",   # Use xlab argument instead of x in xlab() function
                   ylab = "Overall Survival"))
  dev.off()
  
}

GP_km_survival_cptac(gp = "B_CD83", method = "maxstat", study = "CPTAC-3")
