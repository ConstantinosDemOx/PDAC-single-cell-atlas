library(Seurat);library(reticulate);library(tidyverse);library(SeuratWrappers);library(sceasy)

use_condaenv("condapath/python", required = T)

scvi <- import("scvi", convert = F)
np <- import("numpy", convert = F)
os <- import("os", convert = F)
anndata <- import("anndata", convert = F)



# Filter genes: low gene counts and biotypes


pdac.merged <- readRDS("/pdac.merged.rds")

pdac.merged <- pdac.merged %>% NormalizeData()
pdac.merged

pdac.merged <- subset(pdac.merged, subset = HBB < 1.5 & Percent_mito < 15 & HBA1 < 1.5)
pdac.merged

pdac.list <- SplitObject(pdac.merged, split.by = "Study")

g2m_genes <- c("NCAPD2","ANLN","KK-TACC3","HMMR","GTSE1","NDC80","KK-AURKA",
               "TPX2","BIRC5","G2E3","CBX5","RANGAP1","CTCF","CDCA3","TTK",
               "SMC4","ECT2","CENPA","CDC20","NEK2","CENPF","TMPO","HJURP",
               "CKS2","DLGAP5","PIMREG","TOP2A","PSRC1","CDCA8","CKAP2",
               "KK-NUSAP1","KIF23","KIF11","KIF20B","CENPE","GAS2L3","KIF2C",
               "NUF2","KK-ANP32E","LBR","MKI67","CCNB2","CDC25C","HMGB2",
               "CKAP2L","BUB1","CDK1","CKS1B","UBE2C","CKAP5","AURKB","CDCA2",
               "TUBB4B","JPT1")
s_genes <- c("KK-UBR7","RFC2","RAD51","MCM2","TIPIN","MCM6","UNG","POLD3",
             "WDR76","CLSPN","CDC45","CDC6","MSH2","MCM5","POLA1","MCM4",
             "RAD51AP1","GMNN","RPA2","CASP8AP2","HELLS","E2F8","GINS2","PCNA",
             "NASP","BRIP1","DSCC1","DTL","CDCA7","CENPU","ATAD2","CHAF1B",
             "USP1","KK-SLBP","RRM1","FEN1","KK-RRM2","EXO1","CCNE2","TYMS",
             "BLM","KK-PRIM1","UHRF1")

pdac.list <- lapply(X = pdac.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- CellCycleScoring(x,g2m.features = g2m_genes,
                   s.features = s_genes)
  x <- FindVariableFeatures(x, nfeatures = 3000)
})

pdac.list


int.features <-SelectIntegrationFeatures(pdac.list, nfeatures = 3000)


mart_genes <- read.csv("/data/donc-jiang-rop/univ5173/PDAC/Raw_Seurat_files/mart_annot_genes.csv", row.names = 1)
mart_genes_pc <- mart_genes %>% filter(gene_biotype == "protein_coding")

int.features <- intersect(mart_genes_pc$external_gene_name, int.features)
int.features



# Merge normalized samples
pdac.merged <- merge(x = pdac.list[[1]],
                       y = pdac.list[2:length(pdac.list)],
                       merge.data = TRUE)
DefaultAssay(pdac.merged) <- "RNA"

# Manually set variable features of merged Seurat object
VariableFeatures(pdac.merged) <- int.features

rm(pdac.list);gc()


pdac.merged.hv <- pdac.merged[int.features]


adata <- convertFormat(pdac.merged.hv, from="seurat", 
                       to="anndata", main_layer="counts", drop_single_values=FALSE)

adata$X <- adata$X$toarray()


# run setup_anndata, use column stim for batch
scvi$model$SCVI$setup_anndata(adata, batch_key = 'Study.Sample.ID')

# create the model
model = scvi$model$SCVI(adata, n_latent=as.integer(40), gene_likelihood = "nb")

# train the model
model$train(max_epochs = as.integer(200), accelerator = "cpu", batch_size = as.integer(3000))


dir_path = "path/single_cell/scvi_output/"

model$save(dir_path, overwrite = T)

adata$write_h5ad(filename = "path/scvi_output/reference.sc.h5ad")

# get the scVI represenation
scvi_latent = model$get_latent_representation()

# put it back in our original Seurat object
scvi_latent <- as.matrix(scvi_latent)
rownames(scvi_latent) = colnames(pdac.merged)


pdac.merged[["scVI"]] <- CreateDimReducObject(embeddings = scvi_latent, key = "scVI_", assay = DefaultAssay(pdac.merged))


pdac.merged <- RunUMAP(pdac.merged, dims = 1:30, reduction = "scVI")


pdac.merged <- pdac.merged %>% NormalizeData()

VariableFeatures(pdac.merged) <- int.features

pdac.merged <- ScaleData(pdac.merged)

pdac.merged <- FindNeighbors(pdac.merged, reduction = "scVI", dims = 1:30)
pdac.merged <- FindClusters(pdac.merged, resolution = c(0.1,0.2,0.3,0.4))

saveRDS(pdac.merged, "path/scvi_output/reference_sc.rds")



# Annotate clusters using following canonical cell markers


Malignant_Markers <- c("KRT19", "KRT5", "MUC1", "KRT17", "KRT7")
Acinar_Markers <- c("CPA1","PRSS3","AMY2A", "PRSS1")
Acinar_REG_Markers <- c("REG3A", "REG3G", "REG1B")
Ductal_Markers <- c("CFTR", "ANXA4", "SLC4A4")
MUC5B_Ductal_Markers <- c("MUC5B", "CRISP3")
Tuft_Markers <- c("TRPM5", "POU2F3")
Endocrine_Markers <- c("CHGB","CHGA", "IAPP")

# Stromal 
Fibroblast_Markers <- c("CDH11","PDGFRA","PDGFRB","ACTA2")
Schwann_Markers <- c("S100B", "SOX10")
Neurons <- c("TAC1", "CHAT")
Endothelial_Markers <- c("VWF", "PECAM1","CDH5")
Stellate_Markers <- c("RGS5", "PDGFRB")
Adipocyte <- "PLIN1"
Hepatocyte <- c("APOA2", "APOC3")

# Immune: Myeloid
Macrophage_Markers <- c("CD68", "CD163", "ITGAM")
Neutrophil_Markers <- c("CSF3R")
Mast_Markers <- c("CPA3", "TPSAB1", "KIT")
Monocytes <- c("FCN1", "S100A12")
cDCs <- c("FCER1A", "CD1C")
pDCs <- c("IL3RA", "LILRA4")
RBCs <- c("HBB", "HBA2", "HBA1")


# Immune: Lymphoid
T_Cell_Markers <- c("CD3E","CD3G","CD3D")
NK_Markers <- c("NKG7", "GZMA", "KLRD1")
B_Cell_Markers <- c("CD79A","MS4A1","CD19")
Plasma_Markers <- c("CD79A","MZB1")

Dividing_Markers <- c("MKI67", "TOP2A", "STMN1")
# Perform scANVI integration in a new session


scvi <- import("scvi", convert = F)
np <- import("numpy", convert = F)
os <- import("os", convert = F)
anndata <- import("anndata", convert = F)




ref.scvi <- readRDS("/data/donc-jiang-rop/univ5173/PDAC/Integrated_output_files/reference_dataset/single_cell/reference_sc.rds")


Idents(ref.scvi) <- "RNA_snn_res.0.4"

levels(ref.scvi)

ref.scvi@meta.data <- ref.scvi@meta.data %>% 
  mutate(RNA_snn_res.0.4 = factor(RNA_snn_res.0.4, levels = 0:49))

Idents(ref.scvi) <- "RNA_snn_res.0.4"
levels(ref.scvi)


old.idents <- levels(ref.scvi)
old.idents

main_clusters <- c("T", "Fibroblasts", "Macrophages", "Malignant", "Malignant", "Malignant", "B", "Endothelial",
                   "Ductal", "Malignant", "Stellate", "NK", "T", "Monocytes", "Mast", "Neutrophils", "Acinar",
                   "Plasma", "Ductal (Atypical)", "Malignant", "Malignant", "T", "Macrophages", "Endocrine", "Acinar",
                   "Malignant", "Schwann", "Fibroblasts", "Acinar (REG+)", "Dividing", "Fibroblasts", "Remove", "Malignant",
                   "Malignant", "Malignant", "Tuft", "Hepatocytes", "Malignant (dividing)", "Malignant", "pDCs", "Malignant",
                   "Endothelial", "Macrophages", "Endothelial", "Fibroblasts", "T", "Malignant", "Acinar (REG+)", "Stellate", "Remove")

ref.scvi$Main_Clusters <- ref.scvi$RNA_snn_res.0.4

ref.scvi$Main_Clusters <- plyr::mapvalues(x = ref.scvi$Main_Clusters,
                                          from = old.idents,
                                          to = main_clusters)

Idents(ref.scvi) <- "Main_Clusters"



remove.cells <- WhichCells(ref.scvi, idents = "Macrophages", expression = KRT19 > 2)

ref.scvi.filt <- subset(ref.scvi, cells = remove.cells, invert = T)

ref.scvi.filt <-  subset(ref.scvi.filt, subset = Main_Clusters != "Remove")

ref.scvi
ref.scvi.filt

int.features <- VariableFeatures(ref.scvi)

scvi.int.hv <- ref.scvi.filt[int.features]


adata <- convertFormat(scvi.int.hv, from="seurat", 
                       to="anndata", main_layer="counts", drop_single_values=FALSE)

adata$X <- adata$X$toarray()


# run setup_anndata, use column stim for batch
scvi$model$SCANVI$setup_anndata(adata, batch_key = 'Study.Sample.ID', labels_key = "Main_Clusters",
                                unlabeled_category = "Unknown")

# create the model
model = scvi$model$SCANVI(adata, n_latent=as.integer(40), gene_likelihood = "nb")

# train the model
model$train(max_epochs = as.integer(100), accelerator = "cpu", batch_size = as.integer(2000))

dir_path = "path/scanvi_output/"

model$save(dir_path, overwrite = T)

adata$write_h5ad(filename = "path/scanvi_output/reference.sc.h5ad")

# get the scVI represenation
scanvi_latent = model$get_latent_representation()

# put it back in our original Seurat object
scanvi_latent <- as.matrix(scanvi_latent)
rownames(scanvi_latent) = colnames(ref.scvi.filt)


ref.scvi.filt[["scANVI"]] <- CreateDimReducObject(embeddings = scanvi_latent, key = "scANVI_", assay = DefaultAssay(ref.scvi.filt))


ref.scvi.filt <- RunUMAP(ref.scvi.filt, dims = 1:40, reduction = "scANVI", min.dist = 0.15)

ref.scvi.filt <- ref.scvi.filt %>% NormalizeData() %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData()

ref.scvi.filt <- FindNeighbors(ref.scvi.filt, reduction = "scANVI", dims = 1:40)
ref.scvi.filt <- FindClusters(ref.scvi.filt, resolution = c(0.3,0.4))


saveRDS(ref.scvi.filt, "path/scanvi_output/scanvi_reference_sc.rds")












