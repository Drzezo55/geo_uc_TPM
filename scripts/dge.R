BiocManager::install(c("GEOquery", "DESeq2", "pheatmap"), ask = FALSE)

# Load libraries
library(GEOquery)
library(DESeq2)
library(pheatmap)
library(tidyverse)
library(dplyr)
library(ggfortify)

df <- read.delim("data/GSE243625_RNAseq_colon.txt", check.names = FALSE)
df
row.names(df) <- df$`Gene symbol`
# the above will return error due to 1-mar 2-mar as they are are read by r as two similar values
rownames(df) <- make.unique(df$`Gene symbol`) # this will append .1 , .2
df$`Gene symbol` <- NULL

# create a metadata 
samples <- colnames(df)
condition <- ifelse(grepl("^N", samples), "Control",
                    ifelse(grepl("^I", samples), "Inactive_UC",
                           ifelse(grepl("^A", samples), "Active_UC", NA)))

# Filter out any columns not assigned
keep <- !is.na(condition)
df <- df[, keep]
condition <- condition[keep]

# Create metadata
coldata <- data.frame(row.names = colnames(df),
                      condition = factor(condition, levels = c("Control", "Inactive_UC", "Active_UC")))
# in the above code levels should be ordered the same as the order of ifelse loops in the condition vector
all(colnames(df) == rownames(coldata))  # should return TRUE to check if all the same 

# Transpose the data (samples as rows)
log_expr <- log2(df + 1)
pca <- prcomp(t(log_expr), scale. = FALSE)

# Plot PCA
library(ggplot2)
autoplot(pca, data = coldata, colour = "condition")

autoplot(pca, data = coldata, colour = "condition") +
  ggtitle("PCA of Gene Expression") +
  theme_minimal()

# to compare the active versus control
# Subset samples
keep <- coldata$condition %in% c("Active_UC", "Control")
df_sub <- df[, keep]
coldata_sub <- coldata[keep, , drop = FALSE]


#Run t-SNE (samples with similar gene expression locally, and outliers are outliers from one group only) 
#(distance between clusters are meaningless (local))
library(Rtsne)
set.seed(123)  # for reproducibility
tsne_out <- Rtsne(t(df), dims = 2, perplexity = 6)
tsne_df <- data.frame(tsne_out$Y, condition = coldata$condition)
colnames(tsne_df)[1:2] <- c("TSNE1", "TSNE2")
ggplot(tsne_df, aes(x = TSNE1, y = TSNE2, color = condition)) +
  geom_point(size = 3) +
  labs(title = "t-SNE on RNA-seq Samples")

#Run UMAP (global and local as the distance represnt true bioligical difference)
library(uwot)

umap_out <- umap(t(df), n_neighbors = 3)

umap_df <- data.frame(umap_out, condition = coldata$condition)
colnames(umap_df)[1:2] <- c("UMAP1", "UMAP2")

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = condition)) +
  geom_point(size = 3) +
  labs(title = "UMAP on RNA-seq Samples")
#il-6 expression by UMAP
tpm_mat <- as.matrix(df)  # if df_sub is TPM
library(uwot)
umap_out <- umap(t(tpm_mat))  # transpose: samples as rows

umap_df <- data.frame(UMAP1 = umap_out[,1],
                      UMAP2 = umap_out[,2],
                      condition = coldata$condition)

# Add IL6 expression (e.g., TPM values for IL6)
umap_df$IL6 <- tpm_mat["IL6", colnames(tpm_mat)]

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = IL6)) +
  geom_point(size = 3) +
  scale_color_viridis_c() +
  labs(title = "IL6 Expression in UMAP", color = "TPM of IL6") +
  theme_minimal()
library(pheatmap)





























# create the desing matrix
library(edgeR)
library(limma)

# Create a DGEList
dge <- DGEList(counts = df_sub)
# Use 'none' normalization because TPM/FPKM are already normalized
dge <- calcNormFactors(dge, method = "none")
# Design matrix
design <- model.matrix(~0 + coldata_sub$condition)
colnames(design) <- levels(coldata_sub$condition)
# Apply voom transformation
v <- voom(dge, design, plot = TRUE)
# Fit the linear model
fit <- lmFit(v, design)
# Define contrast: Active UC vs Control
contrast_matrix <- makeContrasts(Active_vs_Control = Active_UC - Control, levels = design)
# Fit contrast
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)
# Get results
results <- topTable(fit2, coef = "Active_vs_Control", number = Inf)
head(results)
# Filter for significant DEGs
deg_filtered <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 1, ]
library(EnhancedVolcano)
# volcano plot
EnhancedVolcano(results,
                lab = rownames(results),
                x = 'logFC',
                y = 'P.Value',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = 'Active UC vs Control')

# heatmap
top_degs <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 1, ]
top_genes <- rownames(head(top_degs[order(top_degs$adj.P.Val), ], 50))
# Extract expression of top DEGs
heatmap_data <- df_sub[top_genes, ]
# Scale (Z-score across samples)
heatmap_scaled <- t(scale(t(heatmap_data)))
# Annotation for sample groups
annotation_col <- data.frame(
  Condition = coldata_sub$condition
)
rownames(annotation_col) <- rownames(coldata_sub)
#heatmap plotting
library(pheatmap)

pheatmap(heatmap_scaled,
         annotation_col = annotation_col,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 8,
         fontsize_col = 10,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         scale = "none",  # already scaled above
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = "Top 50 DEGs: Active UC vs Control")











