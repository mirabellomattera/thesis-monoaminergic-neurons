# PCA Analysis - Other ncRNA
library(ggplot2)

cat("Caricamento matrice other_ncRNA...\n")
mat <- read.table(
  "/home/mirab/tesi_bioinfo/expression_matrices/other_ncRNA_tpm_matrix.tsv",
  header=TRUE, sep="\t", row.names=1, check.names=FALSE
)

# Filtra geni espressi
gene_means <- rowMeans(mat)
mat_filtered <- mat[gene_means >= 1, ]
cat("Geni dopo filtro:", nrow(mat_filtered), "\n")

# Log trasformazione
mat_log <- log2(mat_filtered + 1)

# Calcola PCA
pca <- prcomp(t(mat_log), scale=TRUE, center=TRUE)

# Varianza spiegata
var_explained <- round(summary(pca)$importance[2,] * 100, 1)
cat("Varianza spiegata PC1:", var_explained[1], "%\n")
cat("Varianza spiegata PC2:", var_explained[2], "%\n")
cat("Varianza spiegata PC3:", var_explained[3], "%\n")

# Crea dataframe per ggplot
pca_df <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  PC3 = pca$x[,3],
  Condition = c(
    "Undiff", "Undiff",
    "Undiff", "Undiff", "Undiff",
    "Undiff", "Undiff", "Undiff",
    "FGF4", "FGF4", "FGF4",
    "FGF8", "FGF8", "FGF8",
    "FGF4", "FGF4", "FGF4",
    "FGF8", "FGF8", "FGF8",
    "FGF4", "FGF4", "FGF4",
    "FGF8", "FGF8", "FGF8",
    "FGF4", "FGF4", "FGF4"
  ),
  Timepoint = c(
    "Week1", "Week1",
    "Week2", "Week2", "Week2",
    "Week2", "Week2", "Week2",
    "Week3", "Week3", "Week3",
    "Week3", "Week3", "Week3",
    "Week3", "Week3", "Week3",
    "Week4", "Week4", "Week4",
    "Week4", "Week4", "Week4",
    "Week4", "Week4", "Week4",
    "Week4", "Week4", "Week4"
  )
)

condition_colors <- c(Undiff="grey60", FGF4="steelblue", FGF8="tomato")
timepoint_shapes <- c(Week1=15, Week2=16, Week3=17, Week4=18)

# PC1 vs PC2
p1 <- ggplot(pca_df, aes(x=PC1, y=PC2, color=Condition, shape=Timepoint, label=Sample)) +
  geom_point(size=4, alpha=0.9) +
  geom_text(hjust=-0.2, vjust=0.5, size=2.5) +
  scale_color_manual(values=condition_colors) +
  scale_shape_manual(values=timepoint_shapes) +
  labs(
    title="PCA - Other ncRNA Genes",
    x=paste0("PC1 (", var_explained[1], "% variance)"),
    y=paste0("PC2 (", var_explained[2], "% variance)")
  ) +
  theme_bw() +
  theme(plot.title=element_text(size=14, face="bold"), legend.position="right")

ggsave("/home/mirab/tesi_bioinfo/figures/PCA/PCA_other_ncRNA_PC1_PC2.png",
       p1, width=10, height=8, dpi=150)

# PC1 vs PC3
p2 <- ggplot(pca_df, aes(x=PC1, y=PC3, color=Condition, shape=Timepoint, label=Sample)) +
  geom_point(size=4, alpha=0.9) +
  geom_text(hjust=-0.2, vjust=0.5, size=2.5) +
  scale_color_manual(values=condition_colors) +
  scale_shape_manual(values=timepoint_shapes) +
  labs(
    title="PCA - Other ncRNA Genes",
    x=paste0("PC1 (", var_explained[1], "% variance)"),
    y=paste0("PC3 (", var_explained[3], "% variance)")
  ) +
  theme_bw() +
  theme(plot.title=element_text(size=14, face="bold"), legend.position="right")

ggsave("/home/mirab/tesi_bioinfo/figures/PCA/PCA_other_ncRNA_PC1_PC3.png",
       p2, width=10, height=8, dpi=150)

cat("PCA other_ncRNA completata!\n")
