library(ggplot2)
library(reshape2)

# Read the actual data
actual <- read.table("actual_data.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(actual) <- c("Genus", "Actual_Motility")

# Read the predicted data
predicted <- read.table("predicted_data.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(predicted) <- c("Genus", "Predicted_Motility")

# Merge actual and predicted data on the 'Genus' column
merged_data <- merge(actual, predicted, by = "Genus")

# Inspect the merged data
print(merged_data)

# Create the confusion matrix
confusion_matrix <- table(Predicted = merged_data$Predicted_Motility, Actual = merged_data$Actual_Motility)
print(confusion_matrix)

# Assuming confusion matrix is already created and melted
confusion_matrix_melt <- melt(confusion_matrix)

# Plot the confusion matrix with a border
ggplot(data = confusion_matrix_melt, aes(x = Actual, y = Predicted)) +
  geom_tile(aes(fill = value), color = "black", size = 0.5) +  # Add black borders
  scale_fill_gradient(low = "white", high = "steelblue") +
  geom_text(aes(label = value), vjust = 1) +  # Add count values
  labs(title = "Confusion Matrix: Actual vs Predicted Motility", 
       x = "Actual Motility", 
       y = "Predicted Motility") +
  theme_minimal()
