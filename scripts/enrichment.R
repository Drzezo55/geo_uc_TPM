BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "DOSE"))
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)
library(tidyverse)
# Example: top DEGs from limma
deg_genes <- rownames(results[results$adj.P.Val < 0.05 & abs(results$logFC) > 1, ])
# Map gene symbols to Entrez IDs
gene_df <- bitr(deg_genes, fromType = "SYMBOL",
                toType = "ENTREZID", OrgDb = org.Hs.eg.db)
entrez_ids <- gene_df$ENTREZID
# check for the unmapped gene symbols and their percent
unmapped_genes <- setdiff(deg_genes, gene_df$SYMBOL)
# Get the top 50 gene names
top50_genes <- rownames(deg_filtered)[1:50]
# Check which unmapped genes are in the top 50
unmapped_in_top50 <- intersect(unmapped_genes, top50_genes)
unmapped_in_top50
# proceed to GO 
ego <- enrichGO(gene = entrez_ids,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "BP",       # BP: biological process
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2)

# View top terms
head(ego)
dotplot(ego, showCategory = 20) + ggtitle("GO Enrichment (BP)")
emapplot(ego)

library(enrichplot)

ego_sim <- pairwise_termsim(ego)  # compute term similarity
emapplot(ego_sim, showCategory = 30)  # plot top 30 terms

library(dplyr)
top_terms <- ego@result %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::select(Description, p.adjust) %>%
  head(20)
print(top_terms)
dotplot(ego, showCategory = 30) + ggtitle("Top GO Enrichment Terms")
barplot(ego, showCategory = 20, title = "Top GO Terms")
gene_fc_vector <- deg_filtered$logFC
names(gene_fc_vector) <- deg_filtered$ENTREZID
# Create mapping vector
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Hs.eg.db)
cnetplot(ego, categorySize = "pvalue", foldChange = gene_fc_vector)
cnetplot(ego, categorySize = "pvalue")

# proceed to KEGG
ekegg <- enrichKEGG(
  gene = entrez_ids,
  organism = 'hsa',       # 'hsa' = Homo sapiens
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)
dotplot(ekegg, showCategory = 20) + ggtitle("KEGG Pathway Enrichment")
library(enrichplot)
kegg_sim <- pairwise_termsim(ekegg)
cnetplot(ekegg, categorySize="pvalue", foldChange=gene_fc_vector)
emapplot(kegg_sim, showCategory=20)
library(org.Hs.eg.db)
library(AnnotationDbi)

gene_map <- AnnotationDbi::select(org.Hs.eg.db,
                                  keys = names(gene_fc_vector),
                                  columns = "SYMBOL",
                                  keytype = "ENTREZID")

# Replace names with symbols (optional)
names(gene_fc_vector) <- gene_map$SYMBOL[match(names(gene_fc_vector), gene_map$ENTREZID)]
library(org.Hs.eg.db)
library(AnnotationDbi)
cnetplot(ekegg, categorySize="pvalue", foldChange=gene_fc_vector)

# prepare for gsea 

# deg_filtered is your DE results data frame with 'logFC' and 'symbol'
library(clusterProfiler)
library(org.Hs.eg.db)
deg_filtered$SYMBOL <- row.names(deg_filtered)
row.names(deg_filtered) <- NULL
mapped_genes <- bitr(deg_filtered$SYMBOL, fromType="SYMBOL",
                     toType="ENTREZID", OrgDb=org.Hs.eg.db)
# Merge back to get ENTREZID and logFC
deg_filtered_mapped <- merge(deg_filtered, mapped_genes, by.x="SYMBOL", by.y="SYMBOL")
gene_list <- deg_filtered_mapped$logFC
names(gene_list) <- deg_filtered_mapped$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

library(clusterProfiler)

gsea_res <- gseKEGG(geneList = gene_list,
                    organism = "hsa",
                    nPerm = 1000,            # number of permutations
                    minGSSize = 10,
                    maxGSSize = 500,
                    pvalueCutoff = 0.05,
                    verbose = FALSE)

# Example: plot the first pathway in results
gseaplot2(gsea_res, geneSetID = gsea_res$ID[1], title = gsea_res$Description[1])
dotplot(gsea_res, showCategory=20)

emapplot(pairwise_termsim(gsea_res), showCategory=20)

ridgeplot(gsea_res, showCategory = 20)
heatplot(gsea_res, showCategory = 10, foldChange = gene_list)
gseaplot2(gsea_res, 
          geneSetID = c(gsea_res$ID[1:5]), 
          title = "Top 5 KEGG Pathways")




# specific pathway 
pathways <- as.data.frame(ekegg@result)
amoebiasis <- pathways[pathways$ID=="hsa04151",-1]
amoebiasis_ids <- amoebiasis$geneID
amoebiasis_entrez <- unlist(strsplit(amoebiasis_ids, split = "/"))
amoebiasis_genes <- bitr(amoebiasis_entrez,
                         fromType = "ENTREZID",
                         toType = "SYMBOL",
                         OrgDb = org.Hs.eg.db)

print(amoebiasis_genes)

# Optional: Save to CSV
write.csv(amoebiasis_genes, "amoebiasis_genes_hsa05146.csv", row.names = FALSE)

gene_df <- data.frame(gene = amoebiasis_genes$SYMBOL)
# If you want to do PPI only on the DEG–Amoebiasis genes
library(STRINGdb)
string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 400)
mapped <- string_db$map(gene_df, "gene", removeUnmappedRows = TRUE)
# Plot STRING network
string_db$plot_network(mapped$STRING_id)
# OR export and load into Cytoscape for more control
interactions <- string_db$get_interactions(mapped$STRING_id)
write.csv(interactions, "PPI_PI3K_Akt_pathway.csv", row.names = FALSE)



library(igraph)

# community detection
# Simplify edge list using gene names
edge_list <- merge(interactions, mapped[, c("STRING_id", "gene")], 
                   by.x = "from", by.y = "STRING_id")
edge_list <- merge(edge_list, mapped[, c("STRING_id", "gene")], 
                   by.x = "to", by.y = "STRING_id", suffixes = c("_from", "_to"))

ppi_edges <- edge_list[, c("gene_from", "gene_to", "combined_score")]
colnames(ppi_edges) <- c("from", "to", "weight")
ppi_edges$weight <- ppi_edges$weight / 1000
g <- graph_from_data_frame(ppi_edges, directed = FALSE)
# Example: Louvain community detection
communities <- cluster_louvain(g)
V(g)$community <- membership(communities)
plot(g,
     vertex.color = V(g)$community,
     vertex.label = NA,
     vertex.size = 5,
     layout = layout_with_fr(g),
     main = "STRING PPI Network with Communities")

community <- cluster_louvain(g)  # or cluster_walktrap(g)

#Walktrap (good on smaller graphs):
communities <- cluster_walktrap(g)
V(g)$community <- membership(communities)
plot(g,
     vertex.color = V(g)$community,
     vertex.label = NA,
     vertex.size = 5,
     layout = layout_with_fr(g),
     main = "STRING PPI Network with Communities")
#Infomap (information-theory based):
communities <- cluster_infomap(g)
V(g)$community <- membership(communities)
plot(g,
     vertex.color = V(g)$community,
     vertex.label = V(g)$name,
     vertex.size = 5,
     layout = layout_with_fr(g),
     main = "STRING PPI Network with Communities")


#List the Genes in Each Community
community_genes <- split(V(g)$name, V(g)$community)
print(community_genes)
# plot only one community 

# Select community number (e.g., 1)
target_comm <- 1

# Subset the graph
sub_nodes <- V(g)[V(g)$community == target_comm]
sub_g <- induced_subgraph(g, sub_nodes)

# Plot the subgraph
plot(sub_g,
     vertex.label = V(sub_g)$name,
     vertex.color = target_comm,
     vertex.size = 7,
     layout = layout_with_fr(sub_g),
     main = paste("Community", target_comm, "Subnetwork"))


# extract hub genes 
# Compute centrality measures
# Compute centrality measures
deg <- degree(g, mode = "all")
bet <- betweenness(g, directed = FALSE)

# Add to node attributes
V(g)$degree <- deg
V(g)$betweenness <- bet

#Identify Top Hub Genes per Community
# Top genes by degree
degs <- data.frame(genes = V(g)$name,
                   degrees = V(g)$degree)
degs <- degs[order(degs$degrees, decreasing = TRUE),]
top_degs <- degs[1:10,1]
# Top genes by betweenness
bets <- data.frame(genes = V(g)$name,
                   betweenness = V(g)$betweenness)
bets <- bets[order(bets$betweenness, decreasing = TRUE),]
top_bets <- bets[1:10,1]

# Genes common to both (most important!)
common_hubs <- intersect(top_degs, top_bets)
V(g)$is_top_both <- V(g)$name %in% common_hubs
plot(g,
     vertex.color = ifelse(V(g)$is_top_both, "darkred", "lightgray"),
     vertex.size = ifelse(V(g)$is_top_both, 10, 4),
     vertex.label = ifelse(V(g)$is_top_both, V(g)$name, NA),
     layout = layout_with_fr(g),
     main = "Most Influential Genes (Degree + Betweenness)")



modularity(community)
# modularity = how the network will be divided in modules (communities) it tells whether the community detection works well or not 
#usually >4 good , to find the best clustering method 


# centrality is the importane of gene 

library(igraph)

# 1. Degree Centrality
V(g)$degree <- degree(g, mode = "all")

# 2. Betweenness Centrality
V(g)$betweenness <- betweenness(g, directed = FALSE)

# 3. Closeness Centrality
V(g)$closeness <- closeness(g, normalized = TRUE)

# 4. Eigenvector Centrality
V(g)$eigenvector <- eigen_centrality(g)$vector

# 5. PageRank (optional, especially for directed graphs)
V(g)$pagerank <- page_rank(g)$vector


library(igraph)

centrality_table <- data.frame(
  Gene = V(g)$name,
  Degree = V(g)$degree,
  Betweenness = V(g)$betweenness,
  Closeness = V(g)$closeness,
  Eigenvector = V(g)$eigenvector,
  PageRank = V(g)$pagerank
)

# Save to file
write.csv(centrality_table, "Gene_Centrality_Measures.csv", row.names = FALSE)

# View top by degree
head(centrality_table[order(-centrality_table$Degree), ], 10)

















