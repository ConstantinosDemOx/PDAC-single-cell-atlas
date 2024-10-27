# Pathway enrichment analysis

library(Seurat)
library(presto)
library(tidyverse)

source("Consensus method/color_functions.R")

cell.types <- c("Epithelial", "Fibroblast", "Endothelial", "CD4T", "CD8T", "B", "NK",  "MacMono", "Neutrophil", "DCs")

seurat.files <- list.files(path = "Consensus method/Figure 2/Assets", pattern = "top2000.rds", full.names = T, recursive = T)
seurat.files

wilcox.res <- list()

for (i in unique(cell.types)){
  
  print(i)
  path <-   seurat.files[str_detect(seurat.files, paste0("Assets/", i))]
  path
  
  seurat.file <- readRDS(path)
  
  Idents(seurat.file) <- "GP"
  
  # seurat.file <- seurat.file %>% RunUMAP(reduction = "scANVI", dims = 1:20, min.dist = 0.6 )
  
  # seurat.file@reductions$umap@key <- "UMAP_"
  
  # colnames(x = seurat.file[["umap"]]@cell.embeddings) <- paste0("UMAP_", 1:2)
  
  # tiff(paste0("Consensus method/Figure 2/Assets/", i, "/top2000.UMAP.tiff"), width = 900, height = 900, res = 300)
  #  print(DimPlot(seurat.file, group.by = "GP", cols = gp_cols, label.size = 5) +
  #    ggtitle(paste0(i, " (n = ", length(unique(seurat.file$GP)), ")" )) +
  #   theme(aspect.ratio = 1) + NoLegend())
  # dev.off()
  
  
  wilcox.genes <- wilcoxauc(seurat.file, group_by = "GP")
  
  wilcox.genes <- wilcox.genes %>% 
    mutate(GP = paste0(i, "_", group))
  
  wilcox.res[[i]] <- wilcox.genes
  
  
  
  
  
  
  
  
  
  
  
  
  
}

all.wilcox.res <- do.call(rbind,wilcox.res)


# Perform progeny pathway analysis for each cell type

library(decoupleR)

net <- get_progeny(organism = 'human', top = 500)
net

progeny.res <- list()

for (i in unique(cell.types)){
  
  print(i)
  path <-   seurat.files[str_detect(seurat.files, paste0("Assets/", i))]
  path
  
  seurat.file <- readRDS(path)
  
  Idents(seurat.file) <- "GP"
  
  mat <- as.matrix(seurat.file@assays$RNA@data)
  
  # Run mlm
  acts <- run_mlm(mat=mat, net=net, .source='source', .target='target',
                  .mor='weight', minsize = 5)
  acts
  
  # Extract mlm and store it in pathwaysmlm in data
  seurat.file[['pathwaysmlm']] <- acts %>%
    pivot_wider(id_cols = 'source', names_from = 'condition',
                values_from = 'score') %>%
    column_to_rownames('source') %>%
    Seurat::CreateAssayObject(.)
  
  # Change assay
  DefaultAssay(object = seurat.file) <- "pathwaysmlm"
  
  # Scale the data
  seurat.file <- ScaleData(seurat.file)
  seurat.file@assays$pathwaysmlm@data <- seurat.file@assays$pathwaysmlm@scale.data
  
  
  df <- t(as.matrix(seurat.file@assays$pathwaysmlm@data)) %>%
    as.data.frame() %>%
    mutate(cluster = Idents(seurat.file)) %>%
    pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
    group_by(cluster, source) %>%
    summarise(mean = mean(score)) %>% 
    mutate(GP = paste0(i, "_", cluster)) %>% 
    select(-cluster)
  
  
  # Transform to wide matrix
  top_acts_mat <- df %>%
    pivot_wider(id_cols = 'GP', names_from = 'source',
                values_from = 'mean') %>%
    column_to_rownames('GP') %>%
    as.matrix()
  
  # Choose color palette
  palette_length = 100
  my_color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(palette_length)
  
  my_breaks <- c(seq(-2, 0, length.out=ceiling(palette_length/2)),
                 seq(0.05, 2, length.out=floor(palette_length/2)))
  
  
  color_fun <- colorRamp2(breaks = my_breaks, colors = my_color)
  
  
  # Plot
  print(Heatmap(top_acts_mat,  col=color_fun, heatmap_legend_param = list(title = "Pathway\nEnrichment"),
                row_names_side = "left", show_row_dend = F)) 
  
  
  progeny.res[[i]] <- df
  
  
  
  
  
  
  
  
  
  
  
  
  
}

all.progeny.res <- do.call(rbind,progeny.res)

write.csv(all.progeny.res, "Consensus method/Figure 3/Progeny.results.csv")


all.progeny.res <- read.csv("Consensus method/Figure 3/Progeny.results.csv")


for (i in unique(cell.types)){
  
  df <- all.progeny.res
  # Transform to wide matrix
  top_acts_mat <- df %>%
    filter(str_detect(GP, i)) %>% 
    pivot_wider(id_cols = 'cluster', names_from = 'source',
                values_from = 'mean') %>%
    column_to_rownames('cluster') %>%
    as.matrix() %>% 
    scale()
  n_breaks = 10
  rd_bu_colors <- rev(brewer.pal(n = n_breaks, name = "RdBu"))
  
  # Create equally spaced breaks from -max to max, with zero in the middle
  breaks <- seq(-max(top_acts_mat), max(top_acts_mat), length.out = n_breaks)
  
  # Create the color function with colorRamp2, using the breaks and reversed colors
  color_fun <- colorRamp2(breaks = breaks, colors = rd_bu_colors)
  
  # Plot
  tiff(paste0("Consensus method/Supplemental Figure 3/Assets/progeny.",i, ".tiff"), width = 1200, height = 900, res = 300)
  print(Heatmap(top_acts_mat,  col=color_fun, heatmap_legend_param = list(title = "Pathway\nEnrichment"),
                column_title = i,cluster_rows = F,border = T,cluster_columns = T,show_column_dend = F,
                row_order = paste0("GP", 1:length(rownames(top_acts_mat))),
                row_names_side = "left", show_row_dend = F)) 
  dev.off()
  
  
  
  
  
  
  
  
  
  
  
  
  
}





net <- get_collectri(organism='human', split_complexes=FALSE)
net

tf.res <- list()

for (i in unique(cell.types)){
  
  print(i)
  path <-   seurat.files[str_detect(seurat.files, paste0("Assets/", i))]
  path
  
  seurat.file <- readRDS(path)
  
  Idents(seurat.file) <- "GP"
  
  mat <- as.matrix(seurat.file@assays$RNA@data)
  
  # Run ulm
  acts <- run_ulm(mat=mat, net=net, .source='source', .target='target',
                  .mor='mor', minsize = 5)
  acts
  
  # Extract mlm and store it in pathwaysmlm in data
  seurat.file[['tfsulm']] <- acts %>%
    pivot_wider(id_cols = 'source', names_from = 'condition',
                values_from = 'score') %>%
    column_to_rownames('source') %>%
    Seurat::CreateAssayObject(.)
  
  # Change assay
  DefaultAssay(object = seurat.file) <- "tfsulm"
  
  # Scale the data
  seurat.file <- ScaleData(seurat.file)
  seurat.file@assays$tfsulm@data <- seurat.file@assays$tfsulm@scale.data
  
  
  df <- t(as.matrix(seurat.file@assays$tfsulm@data)) %>%
    as.data.frame() %>%
    mutate(cluster = Idents(seurat.file)) %>%
    pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
    group_by(cluster, source) %>%
    summarise(mean = mean(score)) %>% 
    mutate(GP = paste0(i, "_", cluster))
  
  n_tfs <- 25
  
  tfs <- df %>%
    group_by(source) %>%
    summarise(std = sd(mean)) %>%
    arrange(-abs(std)) %>%
    head(n_tfs) %>%
    pull(source)
  
  # Subset long data frame to top tfs and transform to wide matrix
  top_acts_mat <- df %>%
    filter(source %in% tfs) %>%
    pivot_wider(id_cols = 'GP', names_from = 'source',
                values_from = 'mean') %>%
    column_to_rownames('GP') %>%
    as.matrix()
  
  # Choose color palette
  palette_length = 100
  my_color = colorRampPalette(rev(brewer.pal(n = 10 , name = "RdBu")))(palette_length)
  
  
  my_breaks <- c(seq(-2, 0, length.out=ceiling(palette_length/2)),
                 seq(0.05, 2, length.out=floor(palette_length/2)))
  
  color_fun <- colorRamp2(breaks = my_breaks, colors = my_color)
  
  
  # Plot
  print(Heatmap(top_acts_mat,  col=color_fun, heatmap_legend_param = list(title = "Pathway\nEnrichment"),
                row_names_side = "left", show_row_dend = F)) 
  
  
  tf.res[[i]] <- df
  
  
  
  
  
  
  
  
  
  
  
  
  
}

all.tf.res <- do.call(rbind,tf.res)


all.tf.res2 <- all.tf.res %>% 
  mutate(Group = case_when(
    GP %in% c("Epithelial_GP1", "Fibroblast_GP4", "Fibroblast_GP1", "Neutrophil_GP3", "Neutrophil_GP2",
              "B_GP4", "CD4T_GP5", "MacMono_GP5", "B_GP2") ~ "Worse",
    GP %in% c("Epithelial_GP2", "Endothelial_GP3", "CD8T_GP2", "DCs_GP1", "Endothelial_GP4", "DCs_GP4", "Fibroblast_GP2",
              "B_GP3", "CD8T_GP5", "CD8T_GP1") ~ "Better",
    T ~ "Other"
  )) %>% 
  filter(Group != "Other") 

library(rstatix)

wilcox.tf.res <- all.tf.res2 %>% 
  group_by(source) %>% 
  wilcox_test(mean ~ Group) %>% 
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p")



all.tf.res2 %>% 
  filter(source == "GATA6") %>% 
  ggboxplot("Group", "mean", add = "jitter") +
  stat_compare_means()



write.csv(all.tf.res, "Consensus method/Figure 3/TF.results.csv")

all.tf.res <- read.csv("Consensus method/Figure 3/TF.results.csv", row.names = 1)


for (i in unique(cell.types)){
  
  df <- all.tf.res
  
  
  tfs <- df %>%
    filter(str_detect(GP, i)) %>% 
    group_by(cluster) %>%
    slice_max(n = 6, order_by = mean)
  
  
  # Subset long data frame to top tfs and transform to wide matrix
  top_acts_mat <- df %>%
    filter(str_detect(GP, i)) %>% 
    filter(source %in% tfs$source) %>%
    pivot_wider(id_cols = 'cluster', names_from = 'source',
                values_from = 'mean') %>%
    column_to_rownames('cluster') %>%
    select(tfs$source) %>% 
    as.matrix()
  # Define the number of break points (should match the number of colors in the palette)
  n_breaks <- 11  # Or any other suitable number based on the granularity you need
  
  # Generate the reversed RdBu color palette with the number of colors you specified
  rd_bu_colors <- rev(brewer.pal(n = n_breaks, name = "RdBu"))
  
  # Create equally spaced breaks from -max to max, with zero in the middle
  breaks <- seq(-max(top_acts_mat), max(top_acts_mat), length.out = n_breaks)
  
  # Create the color function with colorRamp2, using the breaks and reversed colors
  color_fun <- colorRamp2(breaks = breaks, colors = rd_bu_colors)
  
  # Plot
  
  pixel_per_column <- 65
  
  # Calculate the width based on the number of columns
  dynamic_width <- ncol(top_acts_mat) * pixel_per_column
  
  
  tiff(paste0("Consensus method/Supplemental Figure 3/Assets/TF.",i, ".tiff"), width =dynamic_width, height = 900, res = 300)
  print(Heatmap(top_acts_mat,  col=color_fun, heatmap_legend_param = list(title = "TF Activity\nScore"),
                
                column_title = i,cluster_rows = F,border = T,cluster_columns = F,show_column_dend = F,
                row_order = paste0("GP", 1:length(rownames(top_acts_mat))),
                row_names_side = "left", show_row_dend = F)) 
  dev.off()
  
  
  
  
  
  
  
  
  
  
  
  
  
}

