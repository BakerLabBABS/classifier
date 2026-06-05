library(dplyr)
library(tidyr)

### READ FILES
motility <- read.delim("Flagellar_motility_data.tsv", sep="\t", header=TRUE)
bacteriocin <- read.delim("bacteriocin_data.txt", sep="\t", header=TRUE)

### SPLIT PROTEIN LISTS
mot_long <- motility %>%
  separate_rows(Flagellar_motility, sep=",") %>%
  distinct(Species, Flagellar_motility) %>%   # remove duplicates
  mutate(value = 1)

bacteriocin_long <- bacteriocin %>%
  separate_rows(Bacteriocin, sep=",") %>%
  distinct(Species, Bacteriocin) %>%                 # remove duplicates
  mutate(value = 1)

### CREATE PRESENCE ABSENCE MATRICES
mot_matrix <- mot_long %>%
  pivot_wider(names_from = Flagellar_motility,
              values_from = value,
              values_fill = 0)

bacteriocin_matrix <- bacteriocin_long %>%
  pivot_wider(names_from = Bacteriocin,
              values_from = value,
              values_fill = 0)

### MERGE DATASETS
data <- full_join(mot_matrix, bacteriocin_matrix, by="Species")
data[is.na(data)] <- 0

### IDENTIFY PROTEIN COLUMNS
flagellar <- setdiff(colnames(mot_matrix), "Species")
bacteriocin_proteins <- setdiff(colnames(bacteriocin_matrix), "Species")

### CO-OCCURRENCE TEST
results <- list()
k <- 1

for (f in flagellar) {
  for (t in bacteriocin_proteins) {

    x <- data[[f]]
    y <- data[[t]]

    if(length(x) == length(y)){

      tab <- table(x, y)

      if(all(dim(tab) == c(2,2))){

        test <- fisher.test(tab)

        results[[k]] <- data.frame(
          flagellar_protein = f,
          bacteriocin_proteins = t,
          odds_ratio = test$estimate,
          p_value = test$p.value
        )

        k <- k + 1
      }
    }
  }
}

results <- bind_rows(results)

### MULTIPLE TEST CORRECTION
results$FDR <- p.adjust(results$p_value, method="fdr")

### SORT RESULTS
results <- results %>% arrange(FDR)

### SAVE RESULTS
write.table(results,
            "Flagellar_Bacteriocin_cooccurrence_results.txt",
            sep="\t",
            row.names=FALSE,
            quote=FALSE)

### VIEW RESULTS
head(results)

sig <- results %>% filter(FDR < 0.05)

write.table(sig,
            "Significant_Flagellar_Bacteriocin_associations.txt",
            sep="\t",
            row.names=FALSE,
            quote=FALSE)

