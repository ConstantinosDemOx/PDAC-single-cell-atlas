# Discovery of ecotypes

library(tidyverse)
library(RColorBrewer)
library(ComplexHeatmap)
library(cluster)
library(factoextra)
library(viridis)

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



tcga.gc <- read.csv("Dataframes/Bulk Expression data/TCGA/TCGA.good.quality.csv")

tcga.gc2 <-  tcga.gc %>% 
  mutate(Patient.ID = gsub("-01A", "", Tumor.Sample.ID))

tcga.bp.filt <- tcga.bp %>% 
  filter(rownames(.) %in% tcga.gc2$Patient.ID)



row_sd <- apply(t(tcga.bp.filt), 1, sd)
row_sd

# 2. Determine the 0.05 quantile of the standard deviations
sd_quantile <- quantile(row_sd, 0.05)
sd_quantile


# 3. Filter the matrix to retain rows with high standard deviation
filtered_mat <- tcga.bp.filt[,row_sd > sd_quantile ]

# Print the filtered matrix
filtered_mat


setdiff(colnames(tcga.bp.filt), colnames(filtered_mat))

filtered_mat_z <- filtered_mat %>% select(-c(Epi_AMY2A, Epi_FXYD2, Fib_MYH11)) %>% scale() 


library(ConsensusClusterPlus)


filtered_mat_z <- filtered_mat_z %>% as.data.frame() %>% select(unlist(ce.list, use.names = F)) %>% as.matrix()

rcc = ConsensusClusterPlus(filtered_mat_z,maxK=12,reps=1000,pItem=0.8,pFeature=1,
                           seed = 1000,plot="png",
                           innerLinkage = "ward.D2", 
                           finalLinkage = "ward.D2",
                           title = "Consensus method/Figure 4/Assets/ccp_results",
                           distance="pearson",clusterAlg="hc")


optK <-5

plot_consensus_heatmap <- function(optK, study){
  
  
  k_clusters <- as.data.frame(rcc[[optK]][["consensusClass"]])
  k_clusters <- k_clusters %>% 
    mutate(CE = `rcc[[optK]][["consensusClass"]]`) %>% 
    select(-`rcc[[optK]][["consensusClass"]]`)  %>% 
    mutate(CE = paste0("CE", CE))
  
  
  ce.list <- split(x = rownames(k_clusters), f = k_clusters$CE)
  ce.list
  
  assign(x = paste0(study, ".ce.list"), ce.list, envir = globalenv())
  
  
  res.mat <- rcc[[optK]][["consensusMatrix"]]
  
  colnames(res.mat) <- rownames(k_clusters)
  rownames(res.mat) <- rownames(k_clusters)
  
  cluster_anno <- k_clusters %>% 
    mutate(CE = factor(CE, levels = paste0("CE", 1:optK))) %>% 
    arrange(CE)
  
  clust <- rcc[[optK]]$consensusTree
  clust
  
  
  res.mat2 <- res.mat %>% as.data.frame() %>% select(rownames(cluster_anno)) %>% t() %>% as.data.frame() %>% select(rownames(cluster_anno))
  
  
  identical(rownames(res.mat2), rownames(cluster_anno))
  
  identical(colnames(res.mat2), rownames(cluster_anno))
  
  library(ComplexHeatmap)
  
  annoheatmap1 <- HeatmapAnnotation(df = cluster_anno,
                                    simple_anno_size = unit(0.6, "cm"),show_legend = F,
                                    col = list(CE = CE_cols), border = T, show_annotation_name = F )
  
  
  
  annoheatmap2 <- HeatmapAnnotation(df = cluster_anno,which = "row",
                                    simple_anno_size = unit(0.6, "cm"),show_legend = F,
                                    col = list(CE = CE_cols), border = T, show_annotation_name = F )
  
  
  
  dim(res.mat)
  ht <- Heatmap(res.mat2, 
                name = "my_heatmap",  # Set the heatmap name explicitly
                cluster_rows = F, 
                cluster_columns  = F, 
                top_annotation = annoheatmap1,
                
                column_title = paste0(study ," (n = ", optK, ")" ),
                heatmap_legend_param = list(title = "Similarity\nscore", legend_direction = "vertical"),
                # breaks = c(-0.6, 0, 0.6), 
                border = T, 
                show_column_names  = F, 
                show_row_names = F,
                col  = viridis(n = 10 , option = "D"))
  
  
  col_ord = cluster_anno$CE
  col_ord
  
  col_dup = (which(!duplicated(col_ord)) - 1)
  col_dup
  
  
  col_fract = col_dup / nrow(cluster_anno)
  col_width =  c(col_fract[-1], 1) - col_fract
  col_width
  #col_width <- col_width + 0.1
  # Rows
  
  
  row_ord = cluster_anno$CE
  row_ord
  
  row_dup = (which(!duplicated(row_ord)) - 1)
  
  
  row_dup
  
  row_fract = row_dup / nrow(cluster_anno)
  row_fract
  row_height =  c(row_fract[-1], 1) - row_fract
  row_height
  
  # Save your plot as TIFF
  tiff(paste0("Consensus method/Figure 4/Assets/", "consensus.", study, ".hm.tiff"), width = 800, height = 700, res = 300)
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


plot_consensus_heatmap(optK = 5, study = "TCGA-PAAD")




load("Consensus method/Figure 3/Assets/cptac.bp.rdata")

cptac.theta <- bp.res@Post.ini.cs@theta

cptac.bp <- cptac.theta %>% 
  t() %>% 
  as.data.frame()



row_sd <- apply(t(cptac.bp), 1, sd)
row_sd

# 2. Determine the 0.05 quantile of the standard deviations
sd_quantile <- quantile(row_sd, 0.05)
sd_quantile


# 3. Filter the matrix to retain rows with high standard deviation
filtered_mat <- cptac.bp[,row_sd > sd_quantile ]

# Print the filtered matrix
filtered_mat


setdiff(colnames(cptac.bp), colnames(filtered_mat))

filtered_mat_z <- filtered_mat %>% select(-c(Epi_AMY2A, Epi_FXYD2, Fib_MYH11)) %>% scale() 


fviz_nbclust(cor(filtered_mat_z), FUN = hcut, method = "silhouette")


library(ConsensusClusterPlus)


filtered_mat_z <- filtered_mat_z %>% as.data.frame() %>% select(unlist(ce.list, use.names = F)) %>% as.matrix()

rcc = ConsensusClusterPlus(filtered_mat_z,maxK=12,reps=1000,pItem=0.8,pFeature=1,
                           seed = 1000,plot="png",
                           innerLinkage = "ward.D2", 
                           finalLinkage = "ward.D2",
                           title = "Consensus method/Figure 4/Assets/ccp_results2",
                           distance="pearson",clusterAlg="hc")

plot_consensus_heatmap(optK = 5, study = "CPTAC-3")

names(`TCGA-PAAD.ce.list`) <- paste0("TCGA-", names(`TCGA-PAAD.ce.list`))

names(`CPTAC-3.ce.list`) <- paste0("CPTAC-", names(`CPTAC-3.ce.list`))


all.ce.list <- c(`TCGA-PAAD.ce.list`, `CPTAC-3.ce.list`)


jaccard.matrix <- 1 / (outer(lengths(all.ce.list), lengths(all.ce.list), `+`) / crossprod(table(stack(all.ce.list))) - 1)

jaccard.matrix <- jaccard.matrix %>% 
  as.data.frame() %>% select(matches("CPTAC-")) %>% 
  t() %>% as.data.frame() %>% select(!matches("CPTAC-"))


rownames(jaccard.matrix) <- gsub("CPTAC-", "", rownames(jaccard.matrix))

colnames(jaccard.matrix) <- gsub("TCGA-", "", colnames(jaccard.matrix))


tiff("Consensus method/Figure 4/Assets/TCGA.CPTAC.jaccard.tiff", 
     width = 1200, height = 1000, res = 300)
draw(ComplexHeatmap::pheatmap(jaccard.matrix, display_numbers = T, row_title = "CPTAC-3",
                              column_title = "TCGA-PAAD", show_row_dend = F,
                              show_column_dend = F,color = viridis(n = 10 , option = "D"),border = T,
                              breaks = c(0,0.6),
                              cluster_rows = F, cluster_cols = F,number_color = "black",border_color = NA,
                              heatmap_legend_param = list(title = "Similarity\nscore",
                                                          legend_direction = "vertical"),
                              show_rownames = T, show_colnames = T), heatmap_legend_side = "right")
dev.off()


ce.list.reduced <- list(
  CE1 = intersect(`TCGA-PAAD.ce.list`[[1]], `CPTAC-3.ce.list`[[1]]),
  CE2 = intersect(`TCGA-PAAD.ce.list`[[2]],  `CPTAC-3.ce.list`[[2]]),
  CE3 = intersect(`TCGA-PAAD.ce.list`[[3]],  `CPTAC-3.ce.list`[[3]]),
  CE4 = intersect(`TCGA-PAAD.ce.list`[[4]],  `CPTAC-3.ce.list`[[4]]))

ce.list.reduced


for (CE in names(ce.list.reduced)){
  library(tidygraph);library(igraph);library(ggraph);library(corrr);library(viridis)
  
  cor.graph <- all.cor %>%  
    as.data.frame() %>% 
    select(ce.list.reduced[[CE]]) %>% t() %>% as.data.frame() %>% 
    select(ce.list.reduced[[CE]]) %>% 
    rownames_to_column(var = "x") %>% 
    pivot_longer(cols = -x, names_to = "y", values_to = "r")
  
  cor.graph <- as_tbl_graph(cor.graph, directed = FALSE)
  
  # Car groups info
  group.df <- data_frame(
    name = unique(ce.list.reduced[[CE]]))
  
  group.df <- group.df %>% 
    separate(col = name , sep = "_", into = c("CT", "GP"), remove = F)
  
  
  ce.df <- stack(ce.list.reduced)
  ce.df <- ce.df %>% 
    mutate(CE = ind) %>% 
    column_to_rownames(var = "values")
  
  
  
  
  
  group.df <- ce.df %>% 
    rownames_to_column(var = "name") %>% 
    merge(group.df) %>% 
    mutate(CT = case_when(
      str_detect(CT, "Epi") ~ "Epithelial",
      str_detect(CT, "Fib") ~ "Fibroblast",
      str_detect(CT, "MacMono") ~ "MacMono",
      str_detect(CT, "Neut") ~ "Neutrophil",
      str_detect(CT, "CD4T") ~ "CD4+ T",
      str_detect(CT, "CD8T") ~ "CD8+ T",
      str_detect(CT, "DC") ~ "DCs",
      str_detect(CT, "NK") ~ "NK",
      str_detect(CT, "B") ~ "B",
      str_detect(CT, "EC") ~ "Endothelial",
      str_detect(CT, "Mast") ~ "Mast",
      T ~ "Other"
      
    ))
  
  cor.graph <- cor.graph %>%
    tidygraph::activate(nodes) %>%
    left_join(group.df, by = "name") %>%
    rename(label = GP)
  
  cor.graph <- cor.graph %>%
    activate(edges) %>%
    rename(weight = r)
  
  
  #node_positions <- layout_with_mds(cor.graph,weights = ifelse(E(cor.graph)$weight <= 0, 1, E(cor.graph)$weight))
  set.seed(1000)
  node_positions <- layout_with_fr(cor.graph,weights = ifelse(E(cor.graph)$weight <= 0, 1, E(cor.graph)$weight))
  # node_positions <- layout_with_graphopt(cor.graph)
  x_range <- range(node_positions[, 1])
  y_range <- range(node_positions[, 2])
  
  
  tiff(paste0("Consensus method/Figure 4/Assets/", CE, ".network.tiff"), width = 1000, height = 1000, res = 300)
  set.seed(1000)
  p2 <- print(cor.graph %>%
                activate(nodes) %>%
                ggraph(layout = "fr",weights = ifelse(E(cor.graph)$weight <= 0, 1, E(cor.graph)$weight)) +
                
                geom_edge_link(colour = "grey15",linemitre = 1, width = 1 ) +
                geom_node_point(aes(colour = CT),  size = 15) +
                geom_node_text(aes(label = label), repel = F,size = 4)+
                scale_color_manual(values = ct_cols2) +
                theme(legend.text = element_text(size = 20))+
                ggtitle(CE)+ 
                theme_graph()+
                ggplot2::theme(plot.title = element_text(hjust = 0.5, vjust = -5), legend.position = "none") +
                
                coord_cartesian(xlim= c(x_range[1] -0.5 , x_range[2] + 0.5 ), ylim = c(y_range[1] -0.5, y_range[2]+ 0.5 ) )+
                theme(#aspect.ratio = 1),
                  plot.margin = unit(c(0, 0, 0, 0), "cm") ))
  dev.off()
  
  
  
  
  dummy_data <- data.frame(CT = names(ct_cols2))
  dummy_data <- dummy_data %>% 
    filter(!CT %in% c("Mast", "Plasma"))
  
  
  # Create a dummy ggplot to generate the legend
  dummy_plot <- ggplot(dummy_data, aes(x = CT, y = CT, color = CT)) +
    geom_point(size = 5) +
    scale_color_manual(values = ct_cols2) +
    guides(color = guide_legend( nrow =  1)) +
    labs(color = "Cell type")+
    theme_minimal()
  
  legend <- cowplot::get_legend(dummy_plot)
  
  # Convert the legend to a ggplot object
  legend_plot <- ggplot() + 
    cowplot::draw_grob(legend) + 
    theme_void()  # Remove axes and other plot elements
  
  # Print the legend plot
  print(legend_plot)
  
  # Print the legend plot
  tiff("Consensus method/Figure 4/Assets/legend.tiff", width = 3000, height = 1600, res = 300)
  print(legend_plot)
  dev.off()
}


saveRDS(ce.list.reduced, "Consensus method/Figure 4/Assets/ce.list.rds")
