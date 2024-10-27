# DCs gene programme discovery

library(Seurat);library(tidyverse);library(ggpubr);library(circlize);library(ComplexHeatmap);library(factoextra);library(corrplot)


source("New way/color_functions.R")
source("New way/gene_signatures.R")



myeloid <- readRDS("New way/Figure_4/Files/myeloid.rds")

colnames(myeloid@meta.data)

DimPlot(myeloid, group.by =  "Lvl2_Clusters")

FeaturePlot(myeloid, features = "CD14")



dcs <- subset(myeloid, subset = Lvl2_Clusters %in% c("DCs"))


dcs.hq <- subset(dcs, subset = CD3E < 0.5 & CD3D < 0.5 & CD3G < 0.5 & PTPRC > 1 & CD79B < 1 & AMY2A < 1)
#dcs.hq <- subset(dcs.hq, subset = CD1A > 1 | CLEC9A > 1 | CD1C > 1 | XCR1 > 1 | LILRA4 > 1 )
dcs
dcs.hq

dcs.hq <- dcs.hq %>% NormalizeData() %>% FindVariableFeatures(nfeatures = 2000)
consensus_features <- VariableFeatures(dcs.hq)
consensus_features
consensus_features.df <- consensus_features %>% as.data.frame()



remove.genes <- read.csv("New way/Figure_2/Files/gene signatures/genes_to_remove.csv")

remove.genes <- as.list(remove.genes)

remove.genes <- lapply(remove.genes, function(x) Filter(function(y) y != "" && !is.null(y), x))

remove.genes <- unlist(remove.genes, use.names = F)

vfs_filt <- setdiff(consensus_features, remove.genes)

counts <- GetAssayData(dcs.hq, assay = "RNA")
dim(counts)

keep_rows <- which(rownames(counts) %in% vfs_filt)

counts <- counts[c(keep_rows),]
dim(counts)

dcs.hq.subset <- subset(dcs.hq, features = rownames(counts))

print(dcs.hq)
print(dcs.hq.subset)



gp_generation <- function(x) {
  
  
  cell_type_patients <- as.data.frame(table(x$Study.Patient.ID))
  cell_type_patients <- cell_type_patients %>% filter(Freq >= 100)
  
  
  x_subset <- subset(x, subset = Study.Patient.ID %in% cell_type_patients$Var1)
  
  
  
  cell_type_samples <- unique(x_subset$Study.Patient.ID)
  cell_type_samples
  
  
  resolutions <- c(0.4,0.6,0.8, 1)
  
  
  dea.res <- list()
  all.res <- list()
  
  for (sample in cell_type_samples) {
    
    print(sample)
    
    # Subset sample from mac mono data
    
    x.s <- subset(x_subset, subset = Study.Patient.ID == sample)
    
    x.s <- x.s %>% 
      NormalizeData() %>% FindVariableFeatures(nfeatures = 800) %>% ScaleData() %>% RunPCA()
    
    # Determine percent of variation associated with each PC
    pct <- x.s[["pca"]]@stdev / sum(x.s[["pca"]]@stdev) * 100
    
    # Calculate cumulative percents for each PC
    cumu <- cumsum(pct)
    
    # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
    co1 <- which(cumu > 90 & pct < 5)[1]
    
    co1
    
    # Determine the difference between variation of PC and subsequent PC
    co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
    
    # last point where change of % of variation is more than 0.1%.
    co2
    
    # Minimum of the two calculation
    pcs <- min(co1, co2)
    
    pcs
    
    
    x.s <- FindNeighbors(x.s, dims = 1:pcs)
    # Calculate clusters at different resolutions for the current sample
    
    for (res in resolutions) {
      
      x.s <- FindClusters(x.s, resolution = res)
      
      # Check and skip clusters with fewer than 20 cells
      # cluster_counts <- table(Idents(x.s))
      #small_clusters <- names(cluster_counts[cluster_counts < 20])
      #for (small_cluster in small_clusters) {
      #  x.s <- SetIdent(x.s, small_cluster, "SkipCluster")
      #}
      
      
      
    }
    
    for (resolution_id in c(paste0("RNA_snn_res.",resolutions))) {
      
      Idents(x.s) <- resolution_id
      
      cluster_counts <- table(Idents(x.s))
      cluster_counts
      keep_clusters <- names(cluster_counts[cluster_counts > 10])
      c(keep_clusters)
      
      # Check if the length of Idents() is less than 2, and skip the loop if true
      if (length(unique(keep_clusters)) == 0) {
        cat("Skipping the loop for sample:", sample, "no clusters with more than 20 cells")
        next  # Skip to the next sample
      }
      
      x.s <- subset(x.s, ident = keep_clusters)
      
      Idents(x.s) <- resolution_id
      
      # Check if the length of Idents() is less than 2, and skip the loop if true
      if (length(unique(Idents(x.s))) == 1) {
        cat("Skipping the loop for sample:", sample, "as the length of Idents() is less than 2.\n")
        next  # Skip to the next sample
      }
      
      cluster_degs <- FindAllMarkers(x.s, only.pos = T,  min.pct = 0.3, max.cells.per.ident = 2000)
      
      cluster_degs <- cluster_degs %>% 
        mutate(Group = paste0(sample, "_", resolution_id))  
      
      cluster_degs <- cluster_degs %>%
        mutate(Difference = pct.1 - pct.2) %>% 
        filter(Difference > 0) %>% 
        group_by(cluster) %>% 
        slice_max(n = 100 , order_by = avg_log2FC)
      
      dea.res[[resolution_id]] <- cluster_degs
      
      
      
      
    }
    
    gp.list <- dea.res
    
    gp.df <- do.call(rbind, gp.list)
    
    all.res[[sample]] <- gp.df
    
  }
  
  gp.df <- do.call(rbind, all.res)
  
  
  
  length(unique(gp.df$Group))
  
  gp.df2 <- gp.df %>% 
    mutate(Difference = pct.1 - pct.2) %>% 
    filter(Difference > 0) %>%  
    mutate(Study = sub("^(.*?)_.*$", "\\1", Group)) %>% 
    mutate(Patient.ID = gsub("_RNA_snn_res.(0.4|0.6|0.8|1)", "", Group)) %>% 
    mutate(Unique_GP = paste0(cluster, "_", Group)) %>% 
    group_by(Unique_GP) %>% 
    slice_max(n = 50 , order_by = avg_log2FC)
  
  
  
  gene_signatures_list <- split(x = gp.df2$gene, f = gp.df2$Unique_GP)
  
  # Filter list elements with length less than 20
  gene_signatures_list <- lapply(gene_signatures_list, function(x) {
    if (length(x) > 20) {
      return(x)
    } else {
      return(NULL)
    }
  })
  
  # Remove NULL elements from the list
  gene_signatures_list <- gene_signatures_list[sapply(gene_signatures_list, function(x) !is.null(x))]
  
  
  gp.df2 <- gp.df2 %>% 
    filter(Unique_GP %in% names(gene_signatures_list))
  
  # Keep non-redundant gene signatures (show similarity with other gene signatures (share  >75% of genes with other gene signatures))
  
  assign("gene.programmes",gp.df2,envir = .GlobalEnv)
  
  length(unique(gp.df2$Unique_GP))
 
  
  gene_signatures_list <- split(x = gp.df2$gene, f = gp.df2$Unique_GP)
  
  # Remove NULL elements from the list
  gene_signatures_list <- gene_signatures_list[sapply(gene_signatures_list, function(x) !is.null(x))]
  
  jaccard.matrix <- 1 / (outer(lengths(gene_signatures_list), lengths(gene_signatures_list), `+`) / crossprod(table(stack(gene_signatures_list))) - 1)
  jaccard.matrix[1:25,1:25]
  dim(jaccard.matrix)
  
  
  assign("jaccard.matrix",jaccard.matrix,envir = .GlobalEnv)
  
}

gp_generation(dcs.hq.subset)


dim(jaccard.matrix)
range(jaccard.matrix)

jaccard.matrix[is.infinite(jaccard.matrix)] <- 1

write.csv(gene.programmes, "Consensus method/Figure 2/Assets/DCs/gene.programmes.csv")
write.csv(jaccard.matrix, "Consensus method/Figure 2/Assets/DCs/jaccard.matrix.csv")
gene.programmes <- read.csv("Consensus method/Figure 2/Assets/DCs/gene.programmes.csv", row.names = 1)
jaccard.matrix <- read.csv("Consensus method/Figure 2/Assets/DCs/jaccard.matrix.csv", row.names = 1)



library(ConsensusClusterPlus)

rcc = ConsensusClusterPlus(as.matrix(jaccard.matrix),maxK=10,reps=500,pItem=0.8,pFeature=1,
                           seed = 1000,plot="png",
                           innerLinkage = "average", 
                           finalLinkage = "average",
                           title = "Consensus method/Figure 2/Assets/DCs/ccp_results",
                           distance="pearson",clusterAlg="pam")


saveRDS(rcc, "Consensus method/Figure 2/Assets/DCs/dcs.consensus.results.rds")
rcc <- readRDS("Consensus method/Figure 2/Assets/DCs/dcs.consensus.results.rds")


optK <-4



plot_consensus_heatmap <- function(x){
  
  
  k_clusters <- as.data.frame(rcc[[optK]][["consensusClass"]])
  k_clusters <- k_clusters %>% 
    mutate(GP = `rcc[[optK]][["consensusClass"]]`) %>% 
    select(-`rcc[[optK]][["consensusClass"]]`)  %>% 
    mutate(GP = paste0("GP", GP))
  
  
  res.mat <- rcc[[optK]][["consensusMatrix"]]
  
  colnames(res.mat) <- rownames(k_clusters)
  rownames(res.mat) <- rownames(k_clusters)
  
  cluster_anno <- k_clusters %>% 
    mutate(GP = factor(GP, levels = paste0("GP", 1:optK))) %>% 
    arrange(GP)
  
  clust <- rcc[[optK]]$consensusTree
  clust
  
  
  res.mat2 <- res.mat %>% as.data.frame() %>% select(rownames(cluster_anno)) %>% t() %>% as.data.frame() %>% select(rownames(cluster_anno))
  
  
  identical(rownames(res.mat2), rownames(cluster_anno))
  
  identical(colnames(res.mat2), rownames(cluster_anno))
  
  library(ComplexHeatmap)
  
  annoheatmap1 <- HeatmapAnnotation(df = cluster_anno,
                                    simple_anno_size = unit(0.6, "cm"),show_legend = F,
                                    col = list(GP = gp_cols), border = T, show_annotation_name = F )
  
  
  
  annoheatmap2 <- HeatmapAnnotation(df = cluster_anno,which = "row",
                                    simple_anno_size = unit(0.6, "cm"),show_legend = F,
                                    col = list(GP = gp_cols), border = T, show_annotation_name = F )
  
  
  
  dim(res.mat)
  ht <- Heatmap(res.mat2, 
                name = "my_heatmap",  # Set the heatmap name explicitly
                cluster_rows = F, 
                cluster_columns  = F, 
                top_annotation = annoheatmap1,
                
                column_title = paste0(x ," (n = ", optK, ")" ),
                heatmap_legend_param = list(title = "GP fraction\nz-score", legend_direction = "vertical"),
                # breaks = c(-0.6, 0, 0.6), 
                border = T, 
                show_column_names  = F, 
                show_row_names = F,
                col  = viridis(n = 10 , option = "D"))
  
  
  col_ord = cluster_anno$GP
  col_ord
  
  col_dup = (which(!duplicated(col_ord)) - 1)
  col_dup
  
  
  col_fract = col_dup / nrow(cluster_anno)
  col_width =  c(col_fract[-1], 1) - col_fract
  col_width
  #col_width <- col_width + 0.1
  # Rows
  
  
  row_ord = cluster_anno$GP
  row_ord
  
  row_dup = (which(!duplicated(row_ord)) - 1)
  
  
  row_dup
  
  row_fract = row_dup / nrow(cluster_anno)
  row_fract
  row_height =  c(row_fract[-1], 1) - row_fract
  row_height
  
  # Save your plot as TIFF
  tiff(paste0("Consensus method/Figure 2/Assets/", x, "/consensus.hm.tiff"), width = 800, height = 700, res = 300)
  draw(ht, heatmap_legend_side = "right")
  
  decorate_heatmap_body("my_heatmap", {
    grid.rect(x = unit(col_fract, "native"),
              y =  unit(1-row_fract, "native"), 
              width = unit(col_width, "native"),
              height = unit(row_height, "native"), 
              hjust = 0, vjust = 1,gp = gpar(col = "white", fill = NA, lwd = 3))
  })
  dev.off()
  
  
  
  
  
}


plot_consensus_heatmap("DCs")




gp.k.df <- k_clusters %>% 
  rownames_to_column(var = "Unique_GP") %>% 
  merge(gene.programmes, by = "Unique_GP")

gp.k.df <- gp.k.df %>% 
  group_by(GP) %>% 
  count(gene) %>% 
  slice_max(n = 100 , order_by = n)

gp.k.df <- gp.k.df %>% 
  pivot_wider(id_cols = GP, names_from = gene, values_from = n) %>% 
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>% 
  # mutate(across(where(is.numeric), log1p)) %>% 
  pivot_longer(cols = -GP, names_to = "gene", values_to = "n") 


library(rstatix)
fc_res <- list()
for (i in as.character(unique(gp.k.df$GP))) {
  library(rstatix)
  df <- gp.k.df
  
  df.fc <- df %>% 
    mutate(GP = ifelse(GP == i, i, "Rest")) %>% 
    select(gene,GP, n) %>% 
    group_by(gene,GP) %>% 
    summarise(Mean_Count= mean(n)) %>% 
    pivot_wider(id_cols = gene, names_from = GP, values_from = Mean_Count) %>% 
    rename(Group1 = i) %>% 
    mutate(LogFC = (Group1 + 0.001)/(Rest + 0.001))
  
  df.fc <- df.fc %>% 
    mutate(GP = i)
  
  fc_res[[i]] <- df.fc
}

# Combine marker dataframes into a single dataframe
fc_res <- do.call(rbind, fc_res)


fc_res.top <- fc_res %>% 
  group_by(GP) %>% 
  slice_max(n = 20 , order_by = LogFC) #%>% 
#ungroup() %>% 
#filter(!duplicated(gene))

DCs_gps <- split(x = fc_res.top$gene, f = fc_res.top$GP)

DCs_gps

top_20_genes <- function(genes_vector) {
  genes_vector[1:20]
}


# Apply the function to each vector in the list
DCs_gps <- lapply(DCs_gps, top_20_genes)

DCs_gps
saveRDS(DCs_gps, "Consensus method/Figure 2/Assets/DCs/DCs_gps.rds")
DCs_gps <- readRDS("Consensus method/Figure 2/Assets/DCs/DCs_gps.rds")
# Apply to DCss


colnames(dcs.hq@meta.data)
dcs.hq@meta.data <- dcs.hq@meta.data[,-c(18:length(colnames(dcs.hq@meta.data)))]

dcs.hq <- AddModuleScore(dcs.hq, features = DCs_gps, name = "GP")


dcs.gp_mat <- dcs.hq@meta.data %>% 
  select(names(DCs_gps)) %>% scale() %>% as.data.frame()


assigned_GPs <- apply(dcs.gp_mat, 1, function(row) {
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

dcs.gp_mat$GP <- assigned_GPs


table(dcs.gp_mat$GP)
colnames(dcs.gp_mat)

dcs.gp_mat <- dcs.gp_mat %>% rownames_to_column(var = "Cell.ID")
max_values <- apply(dcs.gp_mat[,2:(optK+1)], 1, max)
table(dcs.gp_mat$GP)

dcs.gp_mat$Max_Variable <- max_values

cell.mat.top2000 <- dcs.gp_mat %>%
  group_by(GP) %>% 
  slice_max(n = 2000, order_by = Max_Variable) %>% 
  filter(GP != "unresolved")

table(cell.mat.top2000$GP)



cell.mat.split <- split(x = dcs.gp_mat$Cell.ID, dcs.gp_mat$GP)


saveRDS(cell.mat.top2000, "Consensus method/Figure 2/Assets/DCs/cell.mat.top2000.csv")
saveRDS(cell.mat.split, "Consensus method/Figure 2/Assets/DCs/cell.mat.split.rds")



dcs.hq@meta.data <- dcs.hq@meta.data %>% 
  mutate(GP = case_when(
    rownames(dcs.hq@meta.data) %in% c(cell.mat.split[["GP1"]]) ~ "GP1",
    rownames(dcs.hq@meta.data) %in% c(cell.mat.split[["GP2"]]) ~ "GP2",
    rownames(dcs.hq@meta.data) %in% c(cell.mat.split[["GP3"]]) ~ "GP3",
    rownames(dcs.hq@meta.data) %in% c(cell.mat.split[["GP4"]]) ~ "GP4",
    T ~ "unresolved"))

dcs.gp_seurat <- subset(dcs.hq, cells = cell.mat.top2000$Cell.ID)

saveRDS(dcs.gp_seurat, "Consensus method/Figure 2/Assets/DCs/dcs.top2000.rds")

dcs.gp_seurat <- dcs.gp_seurat %>% RunUMAP(reduction = "scANVI", dims = 1:20, min.dist = 0.6, n.neighbors = 10 )

dcs.gp_seurat@reductions$umap@key <- "UMAP_"

colnames(x = dcs.gp_seurat[["umap"]]@cell.embeddings) <- paste0("UMAP_", 1:2)

DimPlot(dcs.gp_seurat, group.by = "GP", cols = gp_cols, label.size = 5) +
  ggtitle(paste0("DCs (n = ", optK, ")" )) +
  theme(aspect.ratio = 1) + NoLegend()

FeaturePlot(dcs.gp_seurat, features = "CLEC9A")


#DotPlot(dcs.gp_seurat, group.by = "GP", features = c( "CLEC10A", "IL3RA", "CD68", "LILRA4", "CLEC9A", "MZB1"))
DotPlot(dcs.gp_seurat, group.by = "Anno_Clusters_index", features = names(DCs_gps))


dcs.hq2 <- dcs.hq %>% subset(subset = GP != "unresolved")

dcs.hq
dcs.hq2



dcs.hq2 <- dcs.hq2 %>% RunUMAP(reduction = "scANVI", dims = 1:20, min.dist = 0.6, n.neighbors = 10 )

dcs.hq2@reductions$umap@key <- "UMAP_"

colnames(x = dcs.hq2[["umap"]]@cell.embeddings) <- paste0("UMAP_", 1:2)

tiff("Consensus method/Figure 2/Assets/DCs/dcs.UMAP.tiff", width = 900, height = 900, res = 300)
DimPlot(dcs.hq2, group.by = "GP", cols = gp_cols, label.size = 5) +
  ggtitle(paste0("DCs (n = ", optK, ")" )) +
  theme(aspect.ratio = 1) + NoLegend()
dev.off()

FeaturePlot(dcs.hq2, features = "CD14")

dcs.gp_seurat <- dcs.gp_seurat %>% RunUMAP(reduction = "scANVI", dims = 1:20, min.dist = 0.2, n.neighbors = 100 )

dcs.gp_seurat@reductions$umap@key <- "UMAP_"

colnames(x = dcs.gp_seurat[["umap"]]@cell.embeddings) <- paste0("UMAP_", 1:2)

tiff("Consensus method/Figure 2/Assets/DCs/dcs.UMAP.top2000.tiff", width = 900, height = 900, res = 300)
DimPlot(dcs.gp_seurat, group.by = "GP", cols = gp_cols, label.size = 5) +
  ggtitle(paste0("DCs (n = ", optK, ")" )) +
  theme(aspect.ratio = 1) + NoLegend()
dev.off()



plot_top3_genes <- function(x,y,z){
  library(viridis)
  
  # Function to get the first 3 genes
  top_3_genes <- function(genes_vector) {
    genes_vector[1:2]
  }
  
  
  # Apply the function to each vector in the list
  top_genes_list <- lapply(x, top_3_genes)
  
  
  # Function to get the first 3 genes
  top_20_genes <- function(genes_vector) {
    genes_vector[1:20]
  }
  
  
  # Apply the function to each vector in the list
  top_20 <- lapply(x, top_20_genes)
  
  
  
  
  matrix <- AverageExpression(y, group.by = "GP", features = unlist(top_20, use.names = F), 
                              layer = "data")
  
  
  matrix <- matrix$RNA %>% as.data.frame() %>% t() %>% scale() %>% t() %>% as.data.frame()
  
  mark_at = which(rownames(matrix) %in% unlist(top_genes_list, use.names = F))
  ha = rowAnnotation(foo = anno_mark(side = "left",at = mark_at, labels = unlist(top_genes_list, use.names = F)))
  
  top_annot_df <- data.frame(
    GP2 = paste0("GP", 1:optK)
  )
  
  top_annot_df <- top_annot_df %>% 
    mutate(GP = GP2) %>% 
    column_to_rownames(var = "GP2")
  
  annoheatmap1 <- HeatmapAnnotation(df = top_annot_df,
                                    simple_anno_size = unit(0.6, "cm"),show_legend = F,
                                    col = list(GP = gp_cols), border = T, show_annotation_name = F )
  
  row.df <- stack(DCs_gps)
  row.df <- row.df %>% 
    column_to_rownames(var = "values") %>% 
    rename(GP = ind) %>% 
    arrange(GP)
  
    ht <- Heatmap(matrix, 
                name = "my_heatmap",  # Set the heatmap name explicitly
                cluster_rows = F, 
                cluster_columns  = F, 
                top_annotation = annoheatmap1,
                width = unit(5, "cm"),
                height = unit(6, "cm"),
                left_annotation = ha,
                show_heatmap_legend = F,
                column_title = z,
                column_title_gp = gpar(fontsize = 20),
                heatmap_legend_param = list(title = "GP fraction\nz-score", legend_direction = "vertical"),
                # breaks = c(-0.6, 0, 0.6), 
                border = T, 
                show_column_names  = F, 
                show_row_names = F,
                col  = viridis(n = 10 , option = "D"))
  
  ht
  
  col_ord = top_annot_df$GP
  col_ord
  
  col_dup = (which(!duplicated(col_ord)) - 1)
  col_dup
  
  
  col_fract = col_dup / nrow(top_annot_df)
  col_width =  c(col_fract[-1], 1) - col_fract
  col_width
  #col_width <- col_width + 0.1
  # Rows
  
  
  row_ord = row.df$GP
  row_ord
  
  row_dup = (which(!duplicated(row_ord)) - 1)
  
  
  row_dup
  
  row_fract = row_dup / nrow(row.df)
  row_fract
  row_height =  c(row_fract[-1], 1) - row_fract
  row_height
  
  # Save your plot as TIFF
  tiff(paste0("Consensus method/Figure 2/Assets/", z, "/markers.hm.tiff"), width = 1000, height = 1000, res = 300)
  draw(ht, heatmap_legend_side = "right")
  
  decorate_heatmap_body("my_heatmap", {
    grid.rect(x = unit(col_fract, "native"),
              y =  unit(1-row_fract, "native"), 
              width = unit(col_width, "native"),
              height = unit(row_height, "native"), 
              hjust = 0, vjust = 1,gp = gpar(col = "white", fill = NA, lwd = 3))
  })
  dev.off()
  
  
  
  
  
  
}

plot_top3_genes(DCs_gps, dcs.gp_seurat, "DCs")

fgsea.results <- read.csv("Consensus method/Figure 2/Assets/DCs/fgsea.results.csv")
fgsea.results2 <- fgsea.results %>% 
  #filter(GP == "GP1") %>% 
  group_by(GP) %>% 
  arrange(desc(NES))

DCs_gps



DCs_markers <- list(
  `Antigen\nPresentation` = c("CD1A","HLA-DQB1", "HLA-DQB2", "HLA-DQA1"),
  #cDC2 = c( "CD1C", "FCER1A"),
  `Monocyte-like` = c("S100A9", "CD14", "VCAN"),
  cDC2 = c("CLEC10A", "CD1C"),
 pDC = c("IL3RA", "LILRA4", "CLEC4C", "JCHAIN")
  
)

left_annot_df <- data.frame(
  GP2 = paste0("GP", 1:optK)
)

left_annot_df <- left_annot_df %>% 
  mutate(GP = GP2) %>% 
  column_to_rownames(var = "GP2") %>% 
  arrange(GP)

p_y <-
  ggplot(left_annot_df) +
  geom_col(aes(y = factor(GP, levels = left_annot_df$GP), x = 1, 
               fill = unname(gp_cols[left_annot_df$GP])), width = 1) +
  scale_fill_identity() +
  coord_cartesian(expand = F) +
  theme(
    axis.text.y = element_text(colour = "black", size = 12), 
    axis.title.y = element_blank(), 
    axis.ticks.y = element_blank(), 
    panel.background = element_blank(), 
    plot.background = element_blank())+
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        panel.background = element_blank(), 
        plot.background = element_blank())+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_y_discrete(limits=rev)+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))   # Ensure no margins

Idents(dcs.gp_seurat) <- "GP"
p_dot <- DotPlot(dcs.gp_seurat, group.by = "GP",features = DCs_markers) +
  scale_y_discrete(limits=rev)+
  scale_colour_viridis_c()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_blank(), 
        axis.title.y = element_blank(),  
        plot.margin =  unit(c(0, 0, 0, 0), "cm"), 
        axis.ticks.y = element_blank()) + xlab("")+ ylab("")

p_dot
library(patchwork)
tiff(paste0("Consensus method/Figure 2/Assets/DCs/DotPlot.tiff"), width = 3800, height = 1000, res = 300)
p_y + p_dot +
  plot_layout(widths = unit(c(1, 15), "cm"))
dev.off()


# Pathway enrichment analysis

library(presto)
library(fgsea)
library(msigdbr)
library(data.table)


msigdbr_show_species()

H_df<- msigdbr(species = "Homo sapiens", category = "H")
fgsea_sets1<- H_df %>% split(x = .$gene_symbol, f = .$gs_name)

GO_df<- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
fgsea_sets2 <- GO_df %>% split(x = .$gene_symbol, f = .$gs_name)


gene_signatures <- read.csv("New way/Figure_4/Files/TME gene signatures.csv", check.names = F)
gene_signatures <- gene_signatures %>% 
  select(Hypoxia, Matrix, `Matrix remodeling`, Complement, `Antitumor cytokines`, `Protumor cytokines`, apCAF, myCAF, iCAF, Moffitt_Normal_Stroma, Moffitt_Activated_Stroma, Stellate)

gene_signatures <- as.list(gene_signatures)

gene_signatures <- lapply(gene_signatures, function(x) Filter(function(y) y != "" && !is.null(y), x))

fgsea_subtype <- c(fgsea_sets1, fgsea_sets2, gene_signatures)

Idents(dcs.hq2) <- "GP"

wilcox.genes <- wilcoxauc(dcs.hq2, group_by = "GP")
write.csv(wilcox.genes, "Consensus method/Figure 2/Assets/DCs/wilcox.results.csv")

wilcox.genes.top20 <- wilcox.genes %>% 
  group_by(group) %>% 
  slice_max(n = 20, order = auc)


gsea.gp.list <- list()

for (i in unique(wilcox.genes$group)) {
  
  cluster0.genes<- wilcox.genes %>%
    dplyr::filter(group == i) %>%
    arrange(desc(auc)) %>% 
    dplyr::select(feature, auc)
  
  
  ranks<- deframe(cluster0.genes)
  
  head(ranks)
  
  fgseaRes<- fgsea(fgsea_subtype, stats = ranks, nperm = 1000)
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES)) %>% 
    mutate(GP = i)
  
  gsea.gp.list[[i]] <- fgseaResTidy
  
  
  
  
  
  
  
}

gsea.gp.df <- do.call(rbind, gsea.gp.list)

fwrite(data.frame(gsea.gp.df), "Consensus method/Figure 2/Assets/DCs/fgsea.results.csv")

