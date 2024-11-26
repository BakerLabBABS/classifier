# Load necessary libraries
library(tidyverse)
library(ape)       # To handle phylogenetic trees
library(phangorn)  # For Descendants function
library(here)

# Read the tree from the Nexus file using ape
tree <- read.nexus("~/tree.nex")

# Extract the number of tips and internal nodes
num_tips <- Ntip(tree)
num_nodes <- Nnode(tree)

# Internal nodes are those greater than the number of tips
internal_nodes <- (num_tips + 1):(num_tips + num_nodes)

# Create an empty tibble for internal nodes
d_nodes <- tibble(nodes = internal_nodes, label = NA_character_)

# Function to get tip labels for a given node using phangorn's Descendants
get_tip_labels <- function(node, tree) {
  descendants <- Descendants(tree, node, type = "tips")[[1]]  # Get descendants of the node
  tree$tip.label[descendants]  # Return the tip labels
}

# Calculate tip labels for each internal node
d_nodes <- d_nodes %>%
  rowwise() %>%
  mutate(label = paste0(get_tip_labels(nodes, tree), collapse = ' ')) %>%
  ungroup()

# Generate the commands for reconstructing ancestral states
d_nodes <- d_nodes %>%
  mutate(
    first_command = str_c('AddTag T', nodes, label, sep = ' '),
    second_command = str_c('AddNode ', nodes, ' T', nodes, sep = '')
  )

# Reshape the commands and save as a vector
to_save <- d_nodes %>%
  pivot_longer(cols = c(first_command, second_command), values_to = 'command') %>%
  pull(command)

# Save to file
write.table(to_save, here('~/tree_tag_nodes.txt'), quote = FALSE, col.names = FALSE, row.names = FALSE)

# Display the first few commands
head(to_save)
