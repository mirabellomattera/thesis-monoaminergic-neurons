# ============================================================
# DESeq2 — Grafici: Volcano plot, MA plot, Heatmap top DEG
# FGF4 vs FGF8
# ============================================================

.libPaths(path.expand("~/R/library"))
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

CATEGORIES <- c("protein_coding", "lncRNA", "other_ncRNA")
RES_DIR  <- path.expand("~/tesi_bioinfo/results_deseq2/FGF4_vs_FGF8")
MAT_DIR  <- path.expand("~/tesi_bioinfo/expression_matrices")
FIG_DIR  <- path.expand("~/tesi_bioinfo/figures/DESeq2")
dir.create(FIG_DIR, recursive=TRUE, showWarnings=FALSE)

# Tabella annotazione per bridge gene_id -> gene_name
annot <- read.table(path.expand("~/tesi_bioinfo/reference/annotation_table.tsv"),
                    header=TRUE, sep="\t")
annot <- unique(annot[, c("gene_id", "gene_name")])

for (category in CATEGORIES) {
  cat("================================================\n")
  cat("Grafici:", category, "\n")
  cat("================================================\n")

  # Carica risultati DESeq2
  res_path <- file.path(RES_DIR, paste0(category, "_DESeq2_FGF4vsFGF8.tsv"))
  res <- read.table(res_path, header=TRUE, sep="\t")

  # Colori DEG
  res$color <- "NS"
  res$color[res$padj < 0.05 & res$log2FoldChange >  1] <- "Up in FGF8"
  res$color[res$padj < 0.05 & res$log2FoldChange < -1] <- "Up in FGF4"
  res$color[is.na(res$padj)] <- "NS"

  color_values <- c("NS"="grey70", "Up in FGF8"="tomato", "Up in FGF4"="steelblue")

  # ------------------------------------------------------------
  # 1. VOLCANO PLOT
  # ------------------------------------------------------------

  # Top 15 geni da etichettare con gene_name
  top_genes <- head(res[res$color != "NS", ], 15)
  top_genes <- merge(top_genes, annot, by="gene_id", all.x=TRUE)
  top_genes$label <- ifelse(is.na(top_genes$gene_name) | top_genes$gene_name == "",
                            top_genes$gene_id,
                            top_genes$gene_name)

  p_volcano <- ggplot(res, aes(x=log2FoldChange, y=-log10(padj), color=color)) +
    geom_point(size=1, alpha=0.6) +
    scale_color_manual(values=color_values) +
    geom_vline(xintercept=c(-1, 1), linetype="dashed", color="black", linewidth=0.4) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color="black", linewidth=0.4) +
    geom_text(data=top_genes,
              aes(label=label),
              size=2.5, hjust=-0.1, vjust=0.5, color="black") +
    labs(title=paste("Volcano plot -", category, "- FGF8 vs FGF4"),
         x="log2 Fold Change (FGF8 / FGF4)",
         y="-log10(padj)",
         color="") +
    theme_bw() +
    theme(legend.position="top")

  vol_path <- file.path(FIG_DIR, paste0("volcano_", category, ".png"))
  ggsave(vol_path, p_volcano, width=8, height=6, dpi=150)
  cat("Salvato:", vol_path, "\n")

  # ------------------------------------------------------------
  # 2. MA PLOT
  # ------------------------------------------------------------
  res$ma_color <- "NS"
  res$ma_color[res$padj < 0.05 & res$log2FoldChange >  1] <- "Up in FGF8"
  res$ma_color[res$padj < 0.05 & res$log2FoldChange < -1] <- "Up in FGF4"

  p_ma <- ggplot(res, aes(x=log10(baseMean + 1), y=log2FoldChange, color=ma_color)) +
    geom_point(size=1, alpha=0.5) +
    scale_color_manual(values=color_values) +
    geom_hline(yintercept=c(-1, 0, 1), linetype=c("dashed","solid","dashed"),
               color="black", linewidth=0.4) +
    labs(title=paste("MA plot -", category, "- FGF8 vs FGF4"),
         x="log10(Mean counts + 1)",
         y="log2 Fold Change (FGF8 / FGF4)",
         color="") +
    theme_bw() +
    theme(legend.position="top")

  ma_path <- file.path(FIG_DIR, paste0("MAplot_", category, ".png"))
  ggsave(ma_path, p_ma, width=8, height=6, dpi=150)
  cat("Salvato:", ma_path, "\n")

  # ------------------------------------------------------------
  # 3. HEATMAP TOP DEG
  # ------------------------------------------------------------
  top_deg <- res[res$color != "NS", ]
  top_deg <- top_deg[order(top_deg$padj), ]
  top_deg <- head(top_deg, 50)

  if (nrow(top_deg) == 0) {
    cat("Nessun DEG significativo per heatmap -", category, "\n\n")
    next
  }

  # Carica matrice TPM
  tpm_path <- file.path(MAT_DIR, paste0(category, "_tpm_matrix.tsv"))
  tpm <- read.table(tpm_path, header=TRUE, sep="\t", row.names=1, check.names=FALSE)

  samples_diff <- c("FEU8234A10","FEU8234A11","FEU8234A12",
                    "FEU8234A13","FEU8234A14","FEU8234A15",
                    "FEU8234A16","FEU8234A17","FEU8234A18",
                    "FEU8234A19","FEU8234A20","FEU8234A21",
                    "FEU8234A22","FEU8234A23","FEU8234A24",
                    "FEU8234A25","FEU8234A26","FEU8234A27",
                    "FEU8234A28","FEU8234A29","FEU8234A30")

  tpm_diff <- tpm[, colnames(tpm) %in% samples_diff]

  # Bridge gene_id -> gene_name
  top_deg_names <- merge(top_deg, annot, by="gene_id", all.x=TRUE)
  top_names <- top_deg_names$gene_name[!is.na(top_deg_names$gene_name)]

  tpm_top <- tpm_diff[rownames(tpm_diff) %in% top_names, ]
  cat("  Geni DEG trovati nella matrice TPM:", nrow(tpm_top), "\n")

  if (nrow(tpm_top) < 2) {
    cat("  ATTENZIONE: troppo pochi geni per heatmap - skip\n\n")
    next
  }

  # Z-score
  mat_z <- t(scale(t(tpm_top)))
  mat_z <- mat_z[complete.cases(mat_z), ]

  # Annotazione campioni
  sample_info <- data.frame(
    Condition = c("FGF4","FGF4","FGF4",
                  "FGF8","FGF8","FGF8",
                  "FGF4","FGF4","FGF4",
                  "FGF8","FGF8","FGF8",
                  "FGF4","FGF4","FGF4",
                  "FGF8","FGF8","FGF8",
                  "FGF4","FGF4","FGF4"),
    Timepoint = c("Week3","Week3","Week3",
                  "Week3","Week3","Week3",
                  "Week3","Week3","Week3",
                  "Week4","Week4","Week4",
                  "Week4","Week4","Week4",
                  "Week4","Week4","Week4",
                  "Week4","Week4","Week4"),
    row.names = samples_diff
  )

  ann_colors <- list(
    Condition = c(FGF4="steelblue", FGF8="tomato"),
    Timepoint = c(Week3="plum", Week4="darkgreen")
  )

  heat_path <- file.path(FIG_DIR, paste0("heatmap_topDEG_", category, ".png"))
  png(heat_path, width=1200, height=1400, res=150)
  pheatmap(mat_z,
           annotation_col   = sample_info,
           annotation_colors = ann_colors,
           color            = colorRampPalette(rev(brewer.pal(11,"RdBu")))(100),
           show_rownames    = TRUE,
           show_colnames    = TRUE,
           cluster_cols     = FALSE,
           cluster_rows     = TRUE,
           main             = paste("Top DEG -", category, "- FGF8 vs FGF4"),
           fontsize_row     = 7,
           fontsize_col     = 8)
  dev.off()
  cat("Salvato:", heat_path, "\n\n")
}

cat("=== DONE - tutti i grafici generati ===\n")
