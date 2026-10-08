# ============================================================
# Temporal analysis - Transcription Factors only - FGF4
# Top 50 TF by padj
# ============================================================

.libPaths(path.expand("~/R/library"))
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(tidyr)
library(dplyr)
library(gridExtra)

RES_DIR <- path.expand("~/tesi_bioinfo/results_deseq2/temporal_FGF4")
MAT_DIR <- path.expand("~/tesi_bioinfo/expression_matrices")
FIG_DIR <- path.expand("~/tesi_bioinfo/figures/DESeq2/temporal_FGF4_TF")
dir.create(FIG_DIR, recursive=TRUE, showWarnings=FALSE)

# Load TF list
tf <- read.table(path.expand("~/tesi_bioinfo/reference/Homo_sapiens_TF.txt"),
                 header=TRUE, sep="\t")
tf_clean <- tf[tf$Symbol != "", ]
cat("TF in list:", nrow(tf_clean), "\n")

# Load annotation table
annot <- read.table(path.expand("~/tesi_bioinfo/reference/annotation_table.tsv"),
                    header=TRUE, sep="\t")
annot <- unique(annot[, c("gene_id", "gene_name")])

# FGF4 samples in chronological order (without A1)
samples_FGF4 <- c("FEU8234A2",
                   "FEU8234A4","FEU8234A5","FEU8234A6",
                   "FEU8234A7","FEU8234A8","FEU8234A9",
                   "FEU8234A10","FEU8234A11","FEU8234A12",
                   "FEU8234A16","FEU8234A17","FEU8234A18",
                   "FEU8234A22","FEU8234A23","FEU8234A24",
                   "FEU8234A28","FEU8234A29","FEU8234A30")

timepoint_labels <- c("Week1",
                       "Week2","Week2","Week2",
                       "Week2","Week2","Week2",
                       "Week3","Week3","Week3",
                       "Week3","Week3","Week3",
                       "Week4","Week4","Week4",
                       "Week4","Week4","Week4")

sample_info <- data.frame(
  Timepoint = timepoint_labels,
  row.names = samples_FGF4
)

ann_colors <- list(
  Timepoint = c(Week1="#ffffe0", Week2="#ffa500", Week3="#762a83", Week4="#00441b")
)

# Load TPM matrix
tpm_path <- file.path(MAT_DIR, "protein_coding_tpm_matrix.tsv")
tpm <- read.table(tpm_path, header=TRUE, sep="\t", row.names=1, check.names=FALSE)

# Load DESeq2 temporal results
res_path <- file.path(RES_DIR, "protein_coding_temporal_FGF4.tsv")
res <- read.table(res_path, header=TRUE, sep="\t")

# Bridge gene_id -> gene_name
res <- merge(res, annot, by="gene_id", all.x=TRUE)
res$label <- ifelse(is.na(res$gene_name) | res$gene_name == "",
                    res$gene_id, res$gene_name)

# Filter significant genes
res_sig <- res[!is.na(res$padj) & res$padj < 0.05, ]

# Keep only TF
res_tf <- res_sig[res_sig$label %in% tf_clean$Symbol, ]
res_tf <- res_tf[order(res_tf$padj), ]

# Take top 50 TF
res_tf_top50 <- head(res_tf, 50)
cat("Top 50 TF selected\n")

# Get TPM for these TF
tpm_tf <- tpm[rownames(tpm) %in% res_tf_top50$label, samples_FGF4]
cat("TF found in TPM matrix:", nrow(tpm_tf), "\n")

# ------------------------------------------------------------
# 1. HEATMAP top 50 TF
# ------------------------------------------------------------
mat_z <- t(scale(t(tpm_tf)))
mat_z <- mat_z[complete.cases(mat_z), ]

heat_path <- file.path(FIG_DIR, "heatmap_temporal_FGF4_TF_top50.png")
png(heat_path, width=1600, height=1800, res=150)
pheatmap(mat_z,
         annotation_col    = sample_info,
         annotation_colors = ann_colors,
         color             = colorRampPalette(rev(brewer.pal(11,"RdBu")))(100),
         show_rownames     = TRUE,
         show_colnames     = TRUE,
         cluster_cols      = FALSE,
         cluster_rows      = TRUE,
         main              = "Top 50 Transcription Factors temporal trend - FGF4",
         fontsize_row      = 7,
         fontsize_col      = 8)
dev.off()
cat("Saved:", heat_path, "\n")

# ------------------------------------------------------------
# 2. LINE PLOT split up/down
# ------------------------------------------------------------
tpm_log <- log2(tpm_tf + 1)
tp_order <- c("Week1","Week2","Week3","Week4")

mean_by_tp <- data.frame(
  gene  = rownames(tpm_log),
  Week1 = rowMeans(tpm_log[, timepoint_labels == "Week1", drop=FALSE]),
  Week2 = rowMeans(tpm_log[, timepoint_labels == "Week2", drop=FALSE]),
  Week3 = rowMeans(tpm_log[, timepoint_labels == "Week3", drop=FALSE]),
  Week4 = rowMeans(tpm_log[, timepoint_labels == "Week4", drop=FALSE])
)

mean_by_tp$trend <- ifelse(mean_by_tp$Week4 >= mean_by_tp$Week1, "Up", "Down")
up_genes   <- mean_by_tp$gene[mean_by_tp$trend == "Up"]
down_genes <- mean_by_tp$gene[mean_by_tp$trend == "Down"]

cat("TF up in Week4:", length(up_genes), "\n")
cat("TF down in Week4:", length(down_genes), "\n")

mean_long <- pivot_longer(mean_by_tp, cols=c("Week1","Week2","Week3","Week4"),
                           names_to="Timepoint", values_to="log2TPM")
mean_long$Timepoint <- factor(mean_long$Timepoint, levels=tp_order)

p_up <- ggplot(mean_long[mean_long$gene %in% up_genes, ],
               aes(x=Timepoint, y=log2TPM, group=gene, color=gene)) +
  geom_line(linewidth=0.8) +
  geom_point(size=2) +
  labs(title="Upregulated TF - FGF4",
       x="Timepoint", y="log2(TPM + 1)", color="TF") +
  theme_bw() +
  theme(legend.position="right",
        legend.text=element_text(size=7),
        legend.key.size=unit(0.4, "cm"))

p_down <- ggplot(mean_long[mean_long$gene %in% down_genes, ],
                 aes(x=Timepoint, y=log2TPM, group=gene, color=gene)) +
  geom_line(linewidth=0.8) +
  geom_point(size=2) +
  labs(title="Downregulated TF - FGF4",
       x="Timepoint", y="log2(TPM + 1)", color="TF") +
  theme_bw() +
  theme(legend.position="right",
        legend.text=element_text(size=7),
        legend.key.size=unit(0.4, "cm"))

line_path <- file.path(FIG_DIR, "lineplot_temporal_FGF4_TF_top50_updown.png")
png(line_path, width=2400, height=800, res=120)
grid.arrange(p_up, p_down, ncol=2)
dev.off()
cat("Saved:", line_path, "\n")

cat("=== DONE ===\n")
