# ============================================================
# PCA TPM protein_coding — without outliers A1, A2 and A12
# ============================================================

library(ggplot2)

MAT_DIR <- path.expand("~/tesi_bioinfo/expression_matrices")
FIG_DIR <- path.expand("~/tesi_bioinfo/figures/PCA")
dir.create(FIG_DIR, recursive=TRUE, showWarnings=FALSE)

mat <- read.table(file.path(MAT_DIR, "protein_coding_tpm_matrix.tsv"),
                  header=TRUE, sep="\t", row.names=1, check.names=FALSE)

mat <- mat[rowMeans(mat) >= 1, ]
mat <- mat[, !colnames(mat) %in% c("FEU8234A1", "FEU8234A2", "FEU8234A12")]
cat("Samples:", ncol(mat), "\n")

mat_log <- log2(mat + 1)
pca <- prcomp(t(mat_log), scale=TRUE, center=TRUE)
var_exp <- round(summary(pca)$importance[2, 1:3] * 100, 1)

samples <- colnames(mat)

condition_map <- c(
  FEU8234A4="Undiff", FEU8234A5="Undiff", FEU8234A6="Undiff",
  FEU8234A7="Undiff", FEU8234A8="Undiff", FEU8234A9="Undiff",
  FEU8234A10="FGF4", FEU8234A11="FGF4",
  FEU8234A13="FGF8", FEU8234A14="FGF8", FEU8234A15="FGF8",
  FEU8234A16="FGF4", FEU8234A17="FGF4", FEU8234A18="FGF4",
  FEU8234A19="FGF8", FEU8234A20="FGF8", FEU8234A21="FGF8",
  FEU8234A22="FGF4", FEU8234A23="FGF4", FEU8234A24="FGF4",
  FEU8234A25="FGF8", FEU8234A26="FGF8", FEU8234A27="FGF8",
  FEU8234A28="FGF4", FEU8234A29="FGF4", FEU8234A30="FGF4"
)

timepoint_map <- c(
  FEU8234A4="Week2", FEU8234A5="Week2", FEU8234A6="Week2",
  FEU8234A7="Week2", FEU8234A8="Week2", FEU8234A9="Week2",
  FEU8234A10="Week3", FEU8234A11="Week3",
  FEU8234A13="Week3", FEU8234A14="Week3", FEU8234A15="Week3",
  FEU8234A16="Week3", FEU8234A17="Week3", FEU8234A18="Week3",
  FEU8234A19="Week4", FEU8234A20="Week4", FEU8234A21="Week4",
  FEU8234A22="Week4", FEU8234A23="Week4", FEU8234A24="Week4",
  FEU8234A25="Week4", FEU8234A26="Week4", FEU8234A27="Week4",
  FEU8234A28="Week4", FEU8234A29="Week4", FEU8234A30="Week4"
)

pca_df <- data.frame(
  Sample    = samples,
  PC1       = pca$x[, 1],
  PC2       = pca$x[, 2],
  PC3       = pca$x[, 3],
  Condition = condition_map[samples],
  Timepoint = timepoint_map[samples]
)

pca_df$Condition <- factor(pca_df$Condition, levels=c("Undiff","FGF4","FGF8"))
pca_df$Timepoint <- factor(pca_df$Timepoint, levels=c("Week2","Week3","Week4"))

condition_colors <- c(Undiff="grey60", FGF4="steelblue", FGF8="tomato")
timepoint_shapes <- c(Week2=16, Week3=17, Week4=18)

p1 <- ggplot(pca_df, aes(x=PC1, y=PC2, color=Condition,
                          shape=Timepoint, label=Sample)) +
  geom_point(size=4) +
  geom_text(size=2.5, vjust=-0.8, hjust=0.5) +
  scale_color_manual(values=condition_colors) +
  scale_shape_manual(values=timepoint_shapes) +
  labs(title="PCA TPM - Protein Coding - without outliers A1, A2 and A12",
       x=paste0("PC1 (", var_exp[1], "% variance)"),
       y=paste0("PC2 (", var_exp[2], "% variance)")) +
  theme_bw() +
  theme(legend.position="right")

p1_path <- file.path(FIG_DIR, "PCA_protein_coding_nooutliers_A1A2A12_PC1_PC2.png")
ggsave(p1_path, p1, width=10, height=7, dpi=150)
cat("Saved:", p1_path, "\n")

p2 <- ggplot(pca_df, aes(x=PC1, y=PC3, color=Condition,
                          shape=Timepoint, label=Sample)) +
  geom_point(size=4) +
  geom_text(size=2.5, vjust=-0.8, hjust=0.5) +
  scale_color_manual(values=condition_colors) +
  scale_shape_manual(values=timepoint_shapes) +
  labs(title="PCA TPM - Protein Coding - without outliers A1, A2 and A12",
       x=paste0("PC1 (", var_exp[1], "% variance)"),
       y=paste0("PC3 (", var_exp[3], "% variance)")) +
  theme_bw() +
  theme(legend.position="right")

p2_path <- file.path(FIG_DIR, "PCA_protein_coding_nooutliers_A1A2A12_PC1_PC3.png")
ggsave(p2_path, p2, width=10, height=7, dpi=150)
cat("Saved:", p2_path, "\n")

cat("=== DONE ===\n")
