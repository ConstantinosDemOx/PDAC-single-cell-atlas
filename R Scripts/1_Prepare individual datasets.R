# 1) Load libraries

library(Seurat)
library(SeuratData)
library(tidyverse)


# 2A) Load files - Zhou et al


file.names <- list.files(path = "Dataframes/Raw 10x files/Zhou/", full.names = F, recursive = F)
file.names

Zhou <- vector("list", length = length(file.names))
names(Zhou) <- file.names
Zhou

for(i in file.names){
  Zhou[[i]] <- Read10X(data.dir = paste0("Dataframes/Raw 10x files/Zhou/",i, "/"))
}

names(Zhou)

## Create Seurat object
for(i in names(Zhou)){
  Zhou[[i]] <- CreateSeuratObject(counts = Zhou[[i]], project = names(Zhou[i]), min.cells = 10, min.features = 300)
}


Zhou_merged <- merge(x = Zhou[[1]],
                     y = Zhou[2:length(Zhou)],
                     merge.data = TRUE)


Zhou_merged


# 3) Edit meta-data

Zhou_merged$Sample.ID <- Zhou_merged$orig.ident
Zhou_merged$Study.Sample.ID <- paste0("Zhou_", Zhou_merged$Sample.ID)
Zhou_merged@meta.data <- Zhou_merged@meta.data %>% 
    mutate(Patient.ID = str_remove_all(Sample.ID, "_[1-4]"))
Zhou_merged$Study.Patient.ID <- paste0("Zhou_", Zhou_merged$Patient.ID)
Zhou_merged$Tissue.Type <- "Primary Tumour"
Zhou_merged$Study <- "Zhou et al"
Zhou_merged@meta.data <- Zhou_merged@meta.data %>% 
 mutate(Treatment = case_when(
    str_detect(Sample.ID, "FOL_") ~ "Folfirinox",
    str_detect(Sample.ID, "GP_") ~ "Gemcitabine + nab-paclitaxel",
    str_detect(Sample.ID, "M_") ~ "Folfirinox + Gemcitabine + nab-paclitaxel",
    str_detect(Sample.ID, "CRT_") ~ "Chemo-RT",
    TRUE ~ "Naive"
    ))

meta <- Zhou_merged@meta.data


saveRDS(Zhou_merged, "Dataframes/Raw 10x files/Zhou/Zhou_raw_merged.rds")


# 2A) Load files - Werba et al


file.names <- list.files(path = "Dataframes/Raw files/Werba/", full.names = F, recursive = F)
file.names
file.names <- file.names[!file.names %in% c("Meta_Data")]

Werba <- vector("list", length = length(file.names))
names(Werba) <- file.names
Werba


for(i in file.names){
  Werba[[i]] <- Read10X(data.dir = paste0("Dataframes/Raw files/Werba/",i, "/"))
}

names(Werba)

## Create Seurat object
for(i in names(Werba)){
  Werba[[i]] <- CreateSeuratObject(counts = Werba[[i]], project = names(Werba[i]), min.cells = 10, min.features = 300)
}


Werba_merged <- merge(x = Werba[[1]],
                     y = Werba[2:length(Werba)],
                     merge.data = TRUE)


Werba_merged


# 3) Edit meta-data

Werba_merged$Sample.ID <- Werba_merged$orig.ident
Werba_merged$Study.Sample.ID <- paste0("Werba_", Werba_merged$Sample.ID)
Werba_merged@meta.data <- Werba_merged@meta.data %>% 
  mutate(Patient.ID = Sample.ID)
Werba_merged$Study.Patient.ID <- paste0("Werba_", Werba_merged$Patient.ID)
Werba_merged@meta.data <- Werba_merged@meta.data %>% 
  mutate(Tissue.Type = case_when(
    str_detect(Sample.ID, "PT_") ~ "Primary Tumour",
    str_detect(Sample.ID, "MET_") ~ "Metastatic",
    TRUE ~ Sample.ID
  ))

Werba_merged$Study <- "Werba et al"
Werba_merged@meta.data <- Werba_merged@meta.data %>% 
  mutate(Treatment = case_when(
    str_detect(Sample.ID, "FOL_") ~ "Folfirinox",
    str_detect(Sample.ID, "GP_") ~ "Gemcitabine + nab-paclitaxel",
    Sample.ID == "M_PT_PT10" ~ "Chemo-RT",
    str_detect(Sample.ID, "M_") ~ "Folfirinox + Gemcitabine + nab-paclitaxel",
    Sample.ID == "MET_PT11" ~ "Folfirinox",
    TRUE ~ "Naive"
  ))

meta <- Werba_merged@meta.data

table(meta$Sample.ID, meta$Treatment)


saveRDS(Werba_merged, "Dataframes/Raw files/Werba/Werba_raw_merged.rds")



# 2C) Load files - Steele et al

file.names <- list.files(path = "Dataframes/Raw files/Steele/", full.names = F, recursive = F)
file.names
file.names <- file.names[!file.names %in% c("PT_P14")]


Steele <- vector("list", length = length(file.names))
names(Steele) <- file.names
Steele

for(i in file.names){
  Steele[[i]] <- Read10X(data.dir = paste0("Dataframes/Raw files/Steele/",i, "/"))
}

Steele$PT_P14 <- Read10X_h5(filename = "Dataframes/Raw files/Steele/PT_P14/H5/filtered_feature_bc_matrix.h5")


names(Steele)

## Create Seurat object
for(i in names(Steele)){
  Steele[[i]] <- CreateSeuratObject(counts = Steele[[i]], project = names(Steele[i]), min.cells = 3, min.features = 300)
}


Steele_merged <- merge(x = Steele[[1]],
                       y = Steele[2:length(Steele)],
                       merge.data = TRUE)


Steele_merged

# 3) Edit meta-data

Steele_merged$Sample.ID <- Steele_merged$orig.ident
Steele_merged$Study.Sample.ID <- paste0("Steele_", Steele_merged$Sample.ID)
Steele_merged@meta.data <- Steele_merged@meta.data %>% 
  mutate(Patient.ID = ifelse(Sample.ID %in% c("PT_P11A", "PT_P11B"), "PT_11", Sample.ID))

Steele_merged$Study.Patient.ID <- paste0("Steele_", Steele_merged$Patient.ID)

Steele_merged$Tissue.Type <- plyr::mapvalues(x = Steele_merged$orig.ident,
                                        from = c(file.names, "PT_P14"),
                                        to = rep(c("Normal", "Primary Tumour"), c(3, 17)))

Steele_merged$Study <- "Steele et al"
Steele_merged$Treatment <- "Naive"

meta <- Steele_merged@meta.data


saveRDS(Steele_merged, "Dataframes/Raw files/Steele/Steele_raw_merged.rds")



# 2A) Load files - Oh et al


file.names <- list.files(path = "Dataframes/Raw files/Oh/", full.names = F, recursive = F)
file.names

Oh <- vector("list", length = length(file.names))
names(Oh) <- file.names
Oh

for(i in file.names){
  Oh[[i]] <- Read10X(data.dir = paste0("Dataframes/Raw files/Oh/",i, "/"))
}

names(Oh)

## Create Seurat object
for(i in names(Oh)){
  Oh[[i]] <- CreateSeuratObject(counts = Oh[[i]], project = names(Oh[i]), min.cells = 10, min.features = 300)
}


Oh_merged <- merge(x = Oh[[1]],
                     y = Oh[2:length(Oh)],
                     merge.data = TRUE)


Oh_merged


# 3) Edit meta-data

Oh_merged$Sample.ID <- Oh_merged$orig.ident
Oh_merged$Study.Sample.ID <- paste0("Oh_", Oh_merged$Sample.ID)
Oh_merged@meta.data <- Oh_merged@meta.data %>% 
  mutate(Patient.ID = Sample.ID)
Oh_merged$Study.Patient.ID <- paste0("Oh_", Oh_merged$Patient.ID)
Oh_merged$Tissue.Type <- "Primary Tumour"
Oh_merged$Study <- "Oh et al"
Oh_merged$Treatment <- "Naive"

meta <- Oh_merged@meta.data


saveRDS(Oh_merged, "Dataframes/Raw files/Oh/Oh_raw_merged.rds")


# 2D) Load files - Xue et al


file.names <- list.files(path = "Dataframes/Raw files/Xue/", full.names = F, recursive = F)
file.names

Xue <- vector("list", length = length(file.names))
names(Xue) <- file.names



Xue

for(i in file.names){
  Xue[[i]] <- CreateSeuratObject(counts = read.table(paste0("Dataframes/Raw files/Xue/", i, "/counts.tsv"), header = TRUE, row.names = 1, sep = "\t"),
                                meta.data = read.table(paste0("Dataframes/Raw files/Xue/", i, "/cellname.txt"),header = TRUE, sep = "\t"), project = i )
}


Xue_merged <- merge(x = Xue[[1]],
                   y = Xue[2:length(Xue)],
                   merge.data = TRUE)


Xue_merged

Xue_merged$Sample.ID <- Xue_merged$orig.ident
Xue_merged$Study.Sample.ID <- paste0("Xue_", Xue_merged$Sample.ID)
Xue_merged@meta.data <- Xue_merged@meta.data %>% 
  mutate(Patient.ID = Sample.ID)
Xue_merged$Study.Patient.ID <- paste0("Xue_", Xue_merged$Patient.ID)
Xue_merged$Tissue.Type <- "Primary Tumour"
Xue_merged$Study <- "Xue et al"
Xue_merged$Treatment <- "Naive"

meta <- Xue_merged@meta.data


Xue_merged@meta.data <- Xue_merged@meta.data %>% 
  select(-c(CellName, CellIndex))

saveRDS(Xue_merged, "Dataframes/Raw files/Xue/Xue_raw_merged.rds")


# 2E) Load files - Moncada et al

## loading
Moncada_A <- read.delim("Dataframes/Raw files/Moncada/GSE111672_PDAC-A-indrop-filtered-expMat.txt", sep = "\t")
Moncada_B <- read.delim("Dataframes/Raw files/Moncada/GSE111672_PDAC-B-indrop-filtered-expMat.txt", sep = "\t")

#rownames(Moncada_A) <- make.names(Moncada_A$Genes)
Moncada_A <- distinct(Moncada_A,Genes ,.keep_all = TRUE)
rownames(Moncada_A) <- Moncada_A$Genes
Moncada_A$Genes <- NULL

#rownames(Moncada_B) <- Moncada_B$Genes
Moncada_B <- distinct(Moncada_B,Genes ,.keep_all = TRUE)
rownames(Moncada_B) <- Moncada_B$Genes
Moncada_B$Genes <- NULL

Moncada_A <- CreateSeuratObject(Moncada_A, project = "PT1", min.cells = 3, min.features = 300)
Moncada_B <- CreateSeuratObject(Moncada_B, project = "PT2", min.cells = 3, min.features = 300)
Moncada_merged <- merge(Moncada_A, Moncada_B, add.cell.ids = c("Moncada_PDAC_A","Moncada_PDAC_B"))

Moncada_merged

meta <- Moncada_merged@meta.data


Moncada_merged$Sample.ID <- Moncada_merged$orig.ident
Moncada_merged$Study.Sample.ID <- paste0("Moncada_", Moncada_merged$Sample.ID)
Moncada_merged@meta.data <- Moncada_merged@meta.data %>% 
  mutate(Patient.ID = Sample.ID)
Moncada_merged$Study.Patient.ID <- paste0("Moncada_", Moncada_merged$Patient.ID)
Moncada_merged$Tissue.Type <- "Primary Tumour"
Moncada_merged$Study <- "Moncada et al"
Moncada_merged$Treatment <- "Naive"

meta <- Moncada_merged@meta.data

saveRDS(Moncada_merged, "Dataframes/Raw files/Moncada/Moncada_raw_merged.rds")