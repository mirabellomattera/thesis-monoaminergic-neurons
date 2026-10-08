# Estrazione cluster dalla heatmap - Protein Coding
library(pheatmap)
library(RColorBrewer)

cat("Caricamento matrice...\n")
mat <- read.table(
  "/home/mirab/tesi_bioinfo/expression_matrices/protein_coding_tpm_matrix.tsv",
  header=TRUE, sep="\t", row.names=1, check.names=FALSE
)

# Stesso filtro della heatmap
gene_means <- rowMeans(mat)
mat_filtered <- mat[gene_means >= 1, ]
mat_zscore <- t(scale(t(mat_filtered)))
mat_zscore <- mat_zscore[complete.cases(mat_zscore), ]

cat("Clustering geni...\n")

# Esegui clustering gerarchico con 6 cluster
set.seed(42)
hc <- hclust(dist(mat_zscore), method="complete")
clusters <- cutree(hc, k=6)

# Crea tabella con gene e cluster di appartenenza
cluster_table <- data.frame(
  gene_name = names(clusters),
  cluster = clusters
)

# Salva tabella cluster
out_dir <- "/home/mirab/tesi_bioinfo/expression_matrices"
write.table(cluster_table,
  file.path(out_dir, "protein_coding_clusters.tsv"),
  sep="\t", row.names=FALSE, quote=FALSE
)

cat("Distribuzione geni per cluster:\n")
print(table(clusters))

# Salva un file per ogni cluster
for(cl in 1:6) {
  genes_in_cluster <- cluster_table$gene_name[cluster_table$cluster == cl]
  out_file <- file.path(out_dir, paste0("cluster_", cl, "_genes.txt"))
  writeLines(genes_in_cluster, out_file)
  cat(sprintf("Cluster %d: %d geni -> salvato in %s\n", cl, length(genes_in_cluster), out_file))
}

cat("Estrazione cluster completata!\n")
