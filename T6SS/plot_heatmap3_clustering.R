# Flagellar–Bacteriocin Heatmap 

library(dplyr)
library(tidyr)
library(pheatmap)


#  Read results
results <- read.delim("Flagellar_Bacteriocin_cooccurrence_results.txt",
                      sep="\t",
                      header=TRUE)

# Transform odds ratio
results <- results %>%
  mutate(log2OR = log2(odds_ratio)) %>%
  mutate(log2OR = ifelse(is.infinite(log2OR), NA, log2OR))

#  Add significance column
results <- results %>%
  mutate(significant = ifelse(FDR < 0.05, "YES", "NO"))

#  Save full dataframe
write.table(results,
            "Flagellar_Bacteriocin_heatmap_dataframe.txt",
            sep="\t",
            row.names=FALSE,
            quote=FALSE)


# Prepare heatmap matrix
heatmat_df <- results %>%
  select(flagellar_protein, bacteriocin_proteins, log2OR) %>%
  pivot_wider(names_from = bacteriocin_proteins,
              values_from = log2OR,
              values_fill = NA)

heatmat <- as.matrix(heatmat_df[,-1])
rownames(heatmat) <- heatmat_df[[1]]


# Replace NA with 0 for plotting
heatmat[is.na(heatmat)] <- 0


#  Set figure size
n_col <- ncol(heatmat)
n_row <- nrow(heatmat)

pdf("Flagellar_Bacteriocin_heatmap_Bacteriocincluster.pdf",
    width = 50*n_col/100 + 3,
    height = 30*n_row/100 + 3)


# Plot heatmap
pheatmap(
  heatmat,
  color = colorRampPalette(c("blue","white","red"))(50),
  cluster_rows = FALSE,     # keep flagellar proteins fixed
  cluster_cols = TRUE,      # cluster bacteriocin proteins
  border_color = NA,
  fontsize = 10,
  fontsize_row = 8,
  fontsize_col = 8,
  angle_col = 45,
  main = "Bacteriocin Co-occurrence Across Flagellar Proteins",
  na_col = "white",
  legend = TRUE,
  legend_title = "log2(Odds Ratio)"
)

dev.off()
