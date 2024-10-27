
module_score <- function(x,y) {
  
  
  
  # Assign patient IDs to the row names of signature_scores
  row_names <- colnames(x)
  
  # Assign signature names to the column names of signature_scores
  col_names <- names(y)
  
  # Initialize an empty matrix to store the signature scores
  signature_scores <- matrix(0, nrow = ncol(x), ncol = length(y))
  
  dimnames(signature_scores) <- list(row_names, col_names)
  
  
  # Iterate over each cluster and compute the average expression of signature genes
  for (i in 1:length(y)) {
    signature_genes <- y[[i]]
    
    
    
    # Subset the expression matrix to include only the signature genes
    
    x1 <- data.frame(t(x))
    common_cols <- intersect(colnames(x1), signature_genes)
    
    x1 <- x1 %>% select(common_cols)
    
    x1 <- x1 %>% mutate(Row_Means = rowMeans(x1))
    
    # Calculate the average expression for each sample (row-wise mean)
    average_expression <- rowMeans(x1)
    
    # Store the average expression values as signature scores
    signature_scores[, i] <- average_expression
  }
  
  
  
  gene_average_expression <- rowMeans(x)
  
  
  # Determine the number of bins
  num_bins <- 24
  
  # Divide genes into bins based on average expression
  gene_bins <- cut(gene_average_expression, breaks = num_bins, labels = FALSE)
  
  
  # Set the number of control genes to select from each bin
  num_control_genes_per_bin <- 100
  
  # Initialize an empty vector to store the control gene IDs
  control_gene_ids <- c()
  
  # Select control genes from each bin
  for (bin in 1:num_bins) {
    bin_genes <- rownames(x)[gene_bins == bin]
    
    # Check if the number of genes in the current bin is less than the desired number of control genes
    if (length(bin_genes) < num_control_genes_per_bin) {
      # If there are fewer genes in the bin, select all available genes
      control_genes <- bin_genes
    } else {
      # Randomly select control genes from the current bin
      control_genes <- sample(bin_genes, size = num_control_genes_per_bin)
    }
    
    # Append the control gene IDs to the vector
    control_gene_ids <- c(control_gene_ids, control_genes)
    
  }
  
  
  all_sig_genes <- unlist(y)
  
  control_gene_ids <-  setdiff(control_gene_ids, all_sig_genes)
  
  
  
  # Assuming your bulk RNA-seq expression matrix is called "tcga.exp"
  # Assuming you have randomly selected control genes stored in a vector called "control_genes"
  
  # Calculate the aggregated expression of the control genes (across all samples)
  aggregated_control_expression <- colMeans(x[control_gene_ids, ])
  aggregated_control_expression
  
  # Subtract the aggregated control gene expression from the signature scores
  signature_scores.norm <- signature_scores - aggregated_control_expression
  
  
  
  # Compute the Z-Scores by scaling the signature scores within each sample
  
  # signature_scores.z <- as.data.frame(scale(signature_scores.norm))
  
  # signature_scores.z <- signature_scores.z %>% 
  # filter(TME_Subtype != "Low")
  
  
  
  assign("signature_scores.z",signature_scores.norm,envir = .GlobalEnv)
  
  
  
  
  
  
  
  
}
merge.all <- function(x, ..., by = "row.names") {
  L <- list(...)
  for (i in seq_along(L)) {
    x <- merge(x, L[[i]], by = by)
    rownames(x) <- x$Row.names
    x$Row.names <- NULL
  }
  return(x)
}
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
filter_unique_genes <- function(list_of_genes) {
  seen_genes <- c()
  unique_list <- list()
  
  for (i in names(list_of_genes)) {
    current_genes <- list_of_genes[[i]]
    # Filter out genes that have already been seen
    unique_genes <- setdiff(current_genes, seen_genes)
    # Add unique genes to the list
    unique_list[[i]] <- unique_genes
    # Update the seen genes
    seen_genes <- c(seen_genes, unique_genes)
  }
  
  return(unique_list)
}
silhouette_width <- function(k) {
  cluster_assignment <- cutree(hc, k)
  silhouette_info <- silhouette(cluster_assignment, dist = distance_matrix)
  mean(silhouette_info[, 3])
}
top_10_genes <- function(genes_vector,n) {
  genes_vector[1:10]
}


epi_gps <- readRDS("Consensus method/Figure 2/Assets/Epithelial/epi_gps.rds")

names(epi_gps) <- c("Epi_KRT17", "Epi_TFF1", "Epi_HLA", "Epi_AMY2A", "Epi_FXYD2", "Epi_STMN1")


fib_gps <- readRDS("Consensus method/Figure 2/Assets/Fibroblast/fib_gps.rds")

names(fib_gps) <- c("Fib_MMP11", "Fib_C7", "Fib_IGFBP5", "Fib_ENO1", "Fib_RGS5", "Fib_MYH11", "Fib_IGF1", "Fib_CD74")




ecs_gps <- readRDS("Consensus method/Figure 2/Assets/Endothelial/ecs_gps.rds")

names(ecs_gps) <- c("EC_RGCC", "EC_SEMA3G", "EC_ACKR1", "EC_F10")


cd4_gps <- readRDS("Consensus method/Figure 2/Assets/CD4T/cd4t_gps.rds")

names(cd4_gps) <-  c("CD4T_CCR7","CD4T_ANXA1", "CD4T_FOXP3", "CD4T_GNLY", "CD4T_FOSB", "CD4T_STMN1", "CD4T_CXCL13")



cd8_gps <- readRDS("Consensus method/Figure 2/Assets/CD8T/cd8t_gps.rds")
cd8_gps

names(cd8_gps) <-  c("CD8T_GZMK","CD8T_IL7R", "CD8T_GZMH", "CD8T_KLRC1", "CD8T_CTLA4")


nk_gps <- readRDS("Consensus method/Figure 2/Assets/NK/nk_gps.rds")

names(nk_gps) <- c("NK_XCL1","NK_FCGR3A", "NK_STMN1")



b_gps <- readRDS("Consensus method/Figure 2/Assets/B/B.Pl_gps.rds")
b_gps


names(b_gps) <- c("B_STMN1", "B_CRIP1", "B_TCL1A", "B_CD83", "B_TNF", "B_ANXA1")


macmono_gps <- readRDS("Consensus method/Figure 2/Assets/MacMono/macmono_gps.rds")
macmono_gps
names(macmono_gps) <- c("MacMono_SLC40A1", "MacMono_FCN1", "MacMono_APOE", "MacMono_MMP9",
                        "MacMono_IL1A", "MacMono_IL1RN", "MacMono_STMN1")


neut_gps <- readRDS("Consensus method/Figure 2/Assets/Neutrophil/Neutrophil_gps.rds")

names(neut_gps) <- c("Neut_S100A8", "Neut_CXCR4", "Neut_GBP1")



dcs_gps <- readRDS("Consensus method/Figure 2/Assets/DCs/DCs_gps.rds")

names(dcs_gps) <- c("DC_CD1A", "DC_CD14", "DC_CLEC10A", "DC_IL3RA")


all_gps <- c(epi_gps,fib_gps, ecs_gps, cd4_gps, cd8_gps, nk_gps, b_gps, macmono_gps, neut_gps, dcs_gps)


all_gps <- filter_unique_genes(all_gps)
