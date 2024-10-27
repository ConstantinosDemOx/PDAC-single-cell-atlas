# Colour functions

library(RColorBrewer);library(ggsci);library(rcartocolor);library(wesanderson)

study_cols <- brewer.pal(n = 9 , name = "Set3")
names(study_cols) <- c( "Zhou et al", "Werba et al", "Peng et al", "Zhang et al",
  "Steele et al", "Oh et al", "Lin et al", "Xue et al", "Moncada et al")


tissue_cols <- brewer.pal(n = 3, name = "Set2")
names(tissue_cols) <- c("Normal", "Primary Tumour", "Metastatic")


treatment_cols <- pal_npg(palette = "nrc")(5)
names(treatment_cols)<- c("Naive", "FFX", "GP", "FFX + GP", "Chemo-Radio")


lvl1_cols <- carto_pal(n = 12, name = "Vivid")

names(lvl1_cols) <- c("Epithelial (malignant)", "Epithelial (non-malignant)", "Endocrine", "Endothelial",
  "Stellate", "Fibroblast", "Schwann", "Myeloid", "Mast", "Plasma", "B", "T/NK"
)
cols


category_cols <- pal_jama(palette = "default")(5)
names(category_cols) <- c("Immune", "Proliferation", "Metabolism", "Signalling", "Phenotype")


gp_cols <- pal_d3(palette = "category20")(15)
names(gp_cols) <- paste0("GP", 1:15)


ct_cols <- c(pal_npg(palette = "nrc")(9), "grey")

ct_cols

names(ct_cols) <- c("B", "CD4T", "CD8T", "DC", "EC", "Epi",
                     "MacMono", "Neut", "NK", "Fib")
ct_cols



ct_cols2 <- c(pal_npg(palette = "nrc")(9), "grey", "orange", "purple")

ct_cols2

names(ct_cols2) <- c("B","Plasma" , "CD4+ T", "CD8+ T", "DCs", "Endothelial", "Epithelial",
                    "MacMono", "Neutrophil", "NK", "Fibroblast", "Mast")
ct_cols2











set.seed(100)
CE_cols <- c(wes_palette(name = "Cavalcanti1", n = 5)[c(1:2)],
             wes_palette(name = "Cavalcanti1", n = 5)[c(4:5)],
            wes_palette(name = "Darjeeling2", n = 4)[c(2:4)])

CE_cols <- pal_npg(palette = "nrc")(7)
names(CE_cols) <- paste0("CE", 1:7)


ECO_cols <- pal_npg(palette = "nrc")(7)
names(ECO_cols) <- paste0("ECO", 1:7)



set.seed(100)
CE_spatial_cols <- pal_d3("category20")(20)
names(CE_spatial_cols) <- paste0("SCM", 1:20)



