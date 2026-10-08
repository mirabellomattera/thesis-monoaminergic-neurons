.libPaths(path.expand("~/R/library"))
library(jsonlite)
library(ggplot2)

meme <- fromJSON("~/tesi_bioinfo/results_deseq2/Week4_FGF4_vs_FGF8/ASPM_MEME_result.json")
content <- meme$MLE$content$`0`

df <- data.frame(
  codon = seq_len(nrow(content)),
  pvalue = content[, 7]
)
df$significant <- df$pvalue <= 0.05

p <- ggplot(df, aes(x = codon, y = pvalue, color = significant)) +
  geom_point(size = 1.2) +
  scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "red"), guide = "none") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "steelblue") +
  labs(title = "MEME site-level analysis of ASPM",
       subtitle = "4 of 3483 codon sites under episodic diversifying selection (p \u2264 0.05)",
       x = "Codon position", y = "p-value") +
  theme_minimal(base_size = 13)

ggsave("~/tesi_bioinfo/figures/meme_ASPM_full.png", p, width = 10, height = 5, dpi = 300)
cat("Fatto, immagine salvata.\n")
