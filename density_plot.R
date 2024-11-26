# Load required libraries
library(ggplot2)

# Read the TSV file into a data frame
df <- read.csv('/Users/z5439708/Downloads/histogram_data.tsv', sep='\t', header=FALSE, col.names=c('Species', 'NumGenes'))

# Set font properties 
ggplot_theme <- theme(
  text = element_text(family = "sans", size = 12),  # Changed to default font
  plot.title = element_text(size = 14, face = "bold"),
  axis.title = element_text(size = 12),
  axis.text = element_text(size = 10),
  panel.grid.major.y = element_line(color = "gray", size = 0.5)
)

# Create the density plot
plot <- ggplot(df, aes(x = NumGenes)) +
  geom_density(fill = "#2ca02c", alpha = 0.5) +  # Green fill for a nature theme
  labs(
    title = '',  # Title left blank
    x = 'Number of flagellar proteins',
    y = 'Density'
  ) +
  ggplot_theme +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  # Remove gridlines for better quality

# Save the plot as a PDF using ggsave
ggsave(filename = '/Users/z5439708/Downloads/density_plot.pdf', plot = plot, width = 10, height = 7)
