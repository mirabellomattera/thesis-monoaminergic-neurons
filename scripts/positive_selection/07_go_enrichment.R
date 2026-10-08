.libPaths(c("~/R/library", .libPaths()))
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(dplyr)
library(readr)

# Carica dati v3
df <- read_tsv("~/tesi_bioinfo/positive_selection_orthodb/positive_selection_results_v3.tsv",
               show_col_types=FALSE)

# Lista di interesse: 146 geni con selezione positiva
sig_genes <- df %>% filter(n_branches_significant > 0) %>% pull(gene)

# Background: tutti i 326 geni analizzati
bg_genes <- df %>% pull(gene)

cat("Geni sotto selezione positiva:", length(sig_genes), "\n")
cat("Background (tutti analizzati):", length(bg_genes), "\n")

# Converti gene symbols in Entrez IDs
sig_entrez <- bitr(sig_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
bg_entrez  <- bitr(bg_genes,  fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

cat("Geni convertiti (interest):", nrow(sig_entrez), "\n")
cat("Geni convertiti (background):", nrow(bg_entrez), "\n")

# ---- GO enrichment ----
go_results <- enrichGO(
  gene          = sig_entrez$ENTREZID,
  universe      = bg_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",  # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

cat("\nTermini GO significativi:", nrow(go_results@result %>% filter(p.adjust < 0.05)), "\n")

# ---- KEGG enrichment ----
kegg_results <- enrichKEGG(
  gene         = sig_entrez$ENTREZID,
  universe     = bg_entrez$ENTREZID,
  organism     = "hsa",
  pAdjustMethod= "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

cat("Pathway KEGG significativi:", nrow(kegg_results@result %>% filter(p.adjust < 0.05)), "\n")

# ---- Figure GO ----
if (nrow(go_results) > 0) {
  # Dotplot GO
  p_go <- dotplot(go_results, showCategory=20, title="GO Biological Process enrichment\n(146 genes with positive selection vs 326 analyzed)") +
    theme(plot.title=element_text(face="bold", size=10))
  ggsave("~/tesi_bioinfo/positive_selection_orthodb/v3_fig_go_dotplot.png",
         p_go, width=10, height=8, dpi=300)
  cat("GO dotplot salvato\n")

  # Barplot GO
  p_go_bar <- barplot(go_results, showCategory=20,
                      title="GO Biological Process enrichment") +
    theme(plot.title=element_text(face="bold", size=10))
  ggsave("~/tesi_bioinfo/positive_selection_orthodb/v3_fig_go_barplot.png",
         p_go_bar, width=10, height=8, dpi=300)
  cat("GO barplot salvato\n")
} else {
  cat("Nessun termine GO significativo trovato\n")
}

# ---- Figure KEGG ----
if (nrow(kegg_results) > 0) {
  p_kegg <- dotplot(kegg_results, showCategory=20,
                    title="KEGG pathway enrichment\n(146 genes with positive selection vs 326 analyzed)") +
    theme(plot.title=element_text(face="bold", size=10))
  ggsave("~/tesi_bioinfo/positive_selection_orthodb/v3_fig_kegg_dotplot.png",
         p_kegg, width=10, height=8, dpi=300)
  cat("KEGG dotplot salvato\n")
} else {
  cat("Nessun pathway KEGG significativo trovato\n")
}

# Salva tabelle risultati
write_tsv(go_results@result, "~/tesi_bioinfo/positive_selection_orthodb/go_enrichment_results.tsv")
write_tsv(kegg_results@result, "~/tesi_bioinfo/positive_selection_orthodb/kegg_enrichment_results.tsv")
cat("\nTabelle salvate\n")
