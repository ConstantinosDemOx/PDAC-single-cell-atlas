library(pROC)
library(tidyverse)
library(RColorBrewer)
library(ComplexHeatmap)
library(BayesPrism)

source("Consensus method/color_functions.R")
source("Consensus method/misc_functions.R")

load("Consensus method/Figure 6/Assets/mariathasan.bp.rdata")

imvigor.bp <- bp.res@Post.ini.cs@theta
imvigor.ct <- bp.res@Post.ini.ct@theta

imvigor.ct <- imvigor.ct %>% t() %>% as.data.frame()

ce.list <- readRDS("Consensus method/Figure 4/Assets/ce.list.17102024.rds")


mar.ce.mat <- imvigor.bp%>% t() %>% 
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


max_values <- apply(mar.ce.mat, 1, max)  # Find the maximum value in each row

# Define a cutoff for the fold change (or difference). Example: fold change of 1.33 (i.e., top CE must be 33% greater than second max).
cutoff <- 1

assigned_CEs <- apply(mar.ce.mat, 1, function(row) {
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

# Add the assigned CEs as a new column in mar.ce.mat
mar.ce.mat$Ecotype <- assigned_CEs

# Print the summary of the assignments
table(mar.ce.mat$Ecotype)


mariathan <- readRDS("Dataframes/Bulk Expression data/IMVigor210/mariathan.dataset.rds")


coldat <- as.data.frame(mariathan@colData)
exp.tpm <- as.data.frame(mariathan@assays@data$tpm)

exp.log2.tpm <- log2(exp.tpm)

exp.markers <- exp.log2.tpm %>% t() %>% as.data.frame() %>% select(CD274) 

imvigor.bp <- as.data.frame(t(imvigor.bp))

eco.df <- mar.ce.mat #%>% 
  #filter(Ecotype != "unresolved")


coldat <- coldat %>% 
  rownames_to_column(var = "remove") %>% select(-remove) %>% 
  column_to_rownames(var = "pat_id") %>% na.omit() %>% 
  merge.all(eco.df, exp.markers,imvigor.ct, imvigor.bp) %>% 
  mutate(BOR2 = case_when(BOR == "R" ~ 1,
                          T ~ 0)) %>% 
  mutate(`CE4 + TMB` = CE4 + TMB,
         `CE4 - Epi_KRT17` = CE4 - Epi_KRT17,
         `CE4 - Fib_MMP11` = CE4 - Fib_MMP11,
         `CE4 - Fib_ENO1` = CE4 - Fib_ENO1,
         `CE4 - MacMono_MMP9` = CE4 - MacMono_MMP9)
colnames(coldat)

library(boot)




# Summary of AUC results
cat("Bootstrap AUC for CD274:\n")
roc_cd274 <- roc(coldat$BOR2, glm(coldat$BOR2 ~ coldat$CD274, family = "binomial")$fitted.values, plot = T, legacy.axes = T,
                 col = "#1F78B4", percent = T, print.auc = FALSE)
roc_cd274$auc

ci_cd274 <- ci(roc_cd274) # Calculate 95% CI for AUC
ci_cd274

# Summary of AUC results
cat("Bootstrap AUC for CD8T:\n")
roc_CD8T <- roc(coldat$BOR2, glm(coldat$BOR2 ~ coldat$CD8T, family = "binomial")$fitted.values, plot = T, legacy.axes = T,
                col = "#1F78B4", percent = T, print.auc = FALSE)
roc_CD8T$auc

ci_CD8T <- ci(roc_CD8T) # Calculate 95% CI for AUC
ci_CD8T

# Summary of AUC results
cat("Bootstrap AUC for TMB:\n")
roc_TMB <- roc(coldat$BOR2, glm(coldat$BOR2 ~ coldat$TMB, family = "binomial")$fitted.values, plot = T, legacy.axes = T,
                col = "#1F78B4", percent = T, print.auc = FALSE)
roc_TMB$auc

ci_TMB <- ci(roc_TMB) # Calculate 95% CI for AUC
ci_TMB

# Summary of AUC results
cat("Bootstrap AUC for CE3:\n")
roc_CE3 <- roc(coldat$BOR2, glm(coldat$BOR2 ~ coldat$CE3, family = "binomial")$fitted.values, plot = T, legacy.axes = T,
               col = "#1F78B4", percent = T, print.auc = FALSE)
roc_CE3$auc

ci_CE3 <- ci(roc_CE3) # Calculate 95% CI for AUC
ci_CE3


# Summary of AUC results
cat("Bootstrap AUC for CE4:\n")
roc_CE4 <- roc(coldat$BOR2, glm(coldat$BOR2 ~ coldat$CE4, family = "binomial")$fitted.values, plot = T, legacy.axes = T,
               col = "#1F78B4", percent = T, print.auc = FALSE)
roc_CE4$auc

ci_CE4 <- ci(roc_CE4) # Calculate 95% CI for AUC
ci_CE4


# Summary of AUC results
cat("Bootstrap AUC for CE4:\n")
roc_CE4.KRT17 <- roc(coldat$BOR2, glm(coldat$BOR2 ~ coldat$`CE4 - Epi_KRT17`, family = "binomial")$fitted.values, plot = T, legacy.axes = T,
               col = "purple4", percent = T, print.auc = FALSE)
roc_CE4.KRT17$auc

ci_CE4.KRT17<- ci(roc_CE4.KRT17) # Calculate 95% CI for AUC
ci_CE4.KRT17


set.seed(1000)
res <- roc.test(roc_cd274, roc_CE4, alternative="less")

p.val <- round(res$p.value, digits = 2)
p.val



tiff("Consensus method/Figure 6/Assets/ROC.curve.tiff", width = 1500, height = 1500, res = 300)
plot(roc_cd274, col = "#E31A1C", percent = T, print.auc = FALSE, legacy.axes = T)
plot(roc_CD8T, col = "#33A02C", percent = T, print.auc = FALSE, legacy.axes = T, add = T)
plot(roc_CE4, col = CE_cols[[4]], percent = T, print.auc = FALSE, legacy.axes = T, add = T)
plot(roc_CE4.KRT17, col =  "black", percent = T, print.auc = FALSE, legacy.axes = T, add = T)

# Add AUC and 95% CI text for CD274

text(40, 25, "AUC (95% CI)", col = "black")
text(50, 20, "PDL1     0.73 (0.62-0.84) ", col = "#E31A1C")
text(50, 15, "CD8+ T  0.75 (0.64-0.86) ", col = "#33A02C")
#text(40, 1, "TMB     0.80 (0.68-0.89) ", col = "darkgreen")
text(50, 10, "CE4     0.81 (0.71-0.91) ", col = CE_cols[[4]])
text(50, 5, "CE4 - Epi_KRT17     0.83 (0.75-0.92) ", col = "black")

segments(x0 = 20, y0 = 20.5, x1 = 20, y1 = 10.5, lwd = 2) # Vertical line of bracket
segments(x0 = 22, y0 = 20.5, x1 = 20, y1 = 20.5, lwd = 2) # Top horizontal line
segments(x0 = 22, y0 = 10.5, x1 = 20, y1 = 10.5, lwd = 2)   # Bottom horizontal line

# Add p-value text next to the bracket
text(8, 10, paste("p =", p.val), col = "black")


dev.off()




colnames(coldat)
# List of predictor columns to iterate over (modify as needed)
predictors <- colnames(coldat)[c(2:64, 66:68)]

predictors <- setdiff(predictors, "Ecotype")

# Initialize an empty dataframe to store results
auc_results <- data.frame(Variable = character(), AUC = numeric(), stringsAsFactors = FALSE)

# Loop through each predictor and calculate AUC
for (predictor in predictors) {
  # Fit the logistic regression model
  glm_model <- glm(coldat$BOR2 ~ coldat[[predictor]], family = "binomial")
  
  # Calculate ROC
  roc_result <- roc(coldat$BOR2, glm_model$fitted.values, quiet = TRUE)
  
  # Extract AUC value
  auc_value <- auc(roc_result)
  
  # Add the result to the dataframe
  auc_results <- rbind(auc_results, data.frame(Variable = predictor, AUC = auc_value))
}

# View the resulting dataframe
print(auc_results)

library(rstatix)
wilcox_results.top10 <- coldat %>% 
  rownames_to_column(var = "Patient.ID") %>% 
  select(Patient.ID, BOR2, predictors) %>%
  pivot_longer(cols = predictors, names_to = "Variable", values_to = "Level") %>% 
  group_by(Variable) %>% 
  wilcox_test(Level~BOR2, alternative = "less") %>% 
  mutate(Adj_scale = -log10(p)) %>% 
  arrange(desc(Adj_scale)) %>% 
  mutate(Group = "response") %>% 
  slice_max(n = 10, order_by = Adj_scale)


wilcox_results.bottom10 <- coldat %>% 
  rownames_to_column(var = "Patient.ID") %>% 
  select(Patient.ID, BOR2, predictors) %>%
  pivot_longer(cols = predictors, names_to = "Variable", values_to = "Level") %>% 
  group_by(Variable) %>% 
  wilcox_test(Level~BOR2, alternative = "greater") %>% 
  mutate(Adj_scale = -log10(p)) %>% 
  arrange(desc(Adj_scale)) %>% 
  mutate(Group = "resistant") %>% 
  slice_max(n = 10, order_by = Adj_scale)


all.res <- rbind(wilcox_results.top10, wilcox_results.bottom10)
all.res <- all.res %>% 
  mutate(Adj_scale = case_when(
    Group == "resistant" ~ -Adj_scale,
    T ~ Adj_scale
  ))

tiff("Consensus method/Figure 6/Assets/Wilcox.res.tiff", width = 2500, height = 1200, res = 300)
ggbarplot(all.res, "Variable", "Adj_scale", fill = "Adj_scale", sort.val = "desc",width = 0.8) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  scale_fill_gradient2(midpoint = median(all.res$Adj_scale), 
                       low = "#B2182B", mid = "white", high = "#2166AC", 
                       space = "Lab", 
                       name = "Response", 
                       breaks = c(min(all.res$Adj_scale), median(all.res$Adj_scale), max(all.res$Adj_scale)), # Set breaks
                       labels = c("Worse", "", "Better")) +  # Custom labels for the legend
  ylab("-log10 (p-value)") +xlab("")+
  theme(legend.position = "right")
dev.off()