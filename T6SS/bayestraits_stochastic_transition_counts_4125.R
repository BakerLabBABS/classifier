# BayesTraits dependent-model stochastic mapping
# ------------------------------------------------
# Uses posterior samples of the eight BayesTraits dependent-model transition
# rates to generate stochastic character maps and count realised transitions.
#
# Assumed state coding:
#   Trait 1 = flagella
#   Trait 2 = T6SS
#
# Joint states:
#   00 = flagella absent, T6SS absent
#   01 = flagella absent, T6SS present
#   10 = flagella present, T6SS absent
#   11 = flagella present, T6SS present
#
# BayesTraits rate interpretation:
#   q12: 00 -> 01  T6SS gain, flagella absent
#   q13: 00 -> 10  flagella gain, T6SS absent
#   q21: 01 -> 00  T6SS loss, flagella absent
#   q24: 01 -> 11  flagella gain, T6SS present
#   q31: 10 -> 00  flagella loss, T6SS absent
#   q34: 10 -> 11  T6SS gain, flagella present
#   q42: 11 -> 01  flagella loss, T6SS present
#   q43: 11 -> 10  T6SS loss, flagella present
#
# Tree supplied here: 4125_pruned_tree.nexus (single rooted tree; 4,125 tips)
#
# Packages:
# install.packages(c("ape", "phytools"))

library(ape)
library(phytools)

set.seed(1)

# ---------------------------
# Load data
# ---------------------------

rates_file  <- "bayestraits_rates.tsv"
tree_file   <- "4125_pruned_tree.nexus"  # supplied analysis tree
traits_file <- "4125_metadata_for_bayestraits.txt"  # #Species, FliC, T6SS

# BayesTraits trait files have headerless whitespace-delimitation:
traits_have_header <- TRUE

# Columns in traits_file. These can be numbers or column names.
species_col <- 1
flagella_col <- 2
t6ss_col <- 3

# Number of posterior rate samples to use.
# Set to Inf for the final run using all 4,700 posterior rows.
n_posterior_draws <- Inf

# Number of stochastic maps for each sampled Q.
# One map per posterior Q is usually a good way to integrate over both
# parameter uncertainty and mapping uncertainty. Increase if desired.
maps_per_draw <- 1

# Root prior used by phytools.
# "equal" gives equal prior frequencies at the root.
# "estimated" uses the stationary distribution implied by Q.
# You can also supply a named numeric vector, e.g.
# c("00"=0.25, "01"=0.25, "10"=0.25, "11"=0.25)
#
# Ideally set this to match the BayesTraits analysis.
root_prior <- "equal"

output_prefix <- "bayestraits_transition_mapping_4125"


# ---------------------------
# HELPERS
# ---------------------------

read_phylogeny <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("nex", "nexus")) {
    tr <- ape::read.nexus(path)
  } else {
    tr <- ape::read.tree(path)
  }

  # Always return a list of phylo objects.
  if (inherits(tr, "phylo")) tr <- list(tr)
  tr
}


make_Q <- function(rate_row) {
  states <- c("00", "01", "10", "11")

  Q <- matrix(
    0,
    nrow = 4,
    ncol = 4,
    dimnames = list(states, states)
  )

  Q["00", "01"] <- as.numeric(rate_row[["q12"]])
  Q["00", "10"] <- as.numeric(rate_row[["q13"]])

  Q["01", "00"] <- as.numeric(rate_row[["q21"]])
  Q["01", "11"] <- as.numeric(rate_row[["q24"]])

  Q["10", "00"] <- as.numeric(rate_row[["q31"]])
  Q["10", "11"] <- as.numeric(rate_row[["q34"]])

  Q["11", "01"] <- as.numeric(rate_row[["q42"]])
  Q["11", "10"] <- as.numeric(rate_row[["q43"]])

  diag(Q) <- -rowSums(Q)
  Q
}


make_tip_matrix <- function(species, flagella, t6ss) {
  states <- c("00", "01", "10", "11")

  flagella <- as.character(flagella)
  t6ss <- as.character(t6ss)

  ok <- flagella %in% c("0", "1") & t6ss %in% c("0", "1")
  if (!all(ok)) {
    bad <- which(!ok)
    stop(
      "Traits must be coded 0/1. Bad rows: ",
      paste(head(bad, 20), collapse = ", "),
      if (length(bad) > 20) " ..." else ""
    )
  }

  joint <- paste0(flagella, t6ss)

  X <- matrix(
    0,
    nrow = length(species),
    ncol = 4,
    dimnames = list(as.character(species), states)
  )

  X[cbind(seq_along(joint), match(joint, states))] <- 1
  X
}


summarise_columns <- function(x) {
  data.frame(
    variable = colnames(x),
    mean = vapply(x, mean, numeric(1), na.rm = TRUE),
    sd = vapply(x, sd, numeric(1), na.rm = TRUE),
    median = vapply(x, median, numeric(1), na.rm = TRUE),
    lower_95 = vapply(
      x, function(z) unname(quantile(z, 0.025, na.rm = TRUE)), numeric(1)
    ),
    upper_95 = vapply(
      x, function(z) unname(quantile(z, 0.975, na.rm = TRUE)), numeric(1)
    ),
    row.names = NULL
  )
}


# ---------------------------
# READ INPUTS
# ---------------------------

rates <- read.delim(
  rates_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_rates <- c(
  "q12", "q13", "q21", "q24",
  "q31", "q34", "q42", "q43"
)

missing_rates <- setdiff(required_rates, names(rates))
if (length(missing_rates) > 0) {
  stop("Missing rate columns: ", paste(missing_rates, collapse = ", "))
}

if (anyNA(rates[, required_rates])) {
  stop("NA values found in transition-rate columns.")
}

if (any(as.matrix(rates[, required_rates]) < 0)) {
  stop("Negative transition rates found.")
}

trees <- read_phylogeny(tree_file)

cat("Tree file:", tree_file, "\n")
cat("Trees read:", length(trees), "\n")
cat("Tips in tree 1:", length(trees[[1]]$tip.label), "\n")
cat("Total branch length, tree 1:", sum(trees[[1]]$edge.length), "\n")
cat("Posterior rows in rate file:", nrow(rates), "\n")
if ("Tree No" %in% names(rates)) {
  cat("Tree numbers used in posterior:",
      paste(sort(unique(rates[["Tree No"]])), collapse = ", "), "\n")
}

if (!file.exists(traits_file)) {
  stop(
    "Trait data file not found: '", traits_file, "'.\n",
    "The Nexus tree contains topology/taxon names but not the observed ",
    "flagella/T6SS tip states. Supply the BayesTraits data file (or an ",
    "equivalent table with species, flagella, T6SS) and set traits_file."
  )
}

traits <- read.table(
  traits_file,
  header = traits_have_header,
  stringsAsFactors = FALSE,
  comment.char = "",
  quote = ""
)

species  <- traits[[species_col]]
flagella <- traits[[flagella_col]]
t6ss     <- traits[[t6ss_col]]

if (anyDuplicated(species)) {
  stop("Duplicated species names in traits file.")
}

tip_states <- make_tip_matrix(species, flagella, t6ss)

# supplied data set should contain exactly the same taxa as the tree.
extra_trait_taxa <- setdiff(rownames(tip_states), trees[[1]]$tip.label)
missing_trait_taxa <- setdiff(trees[[1]]$tip.label, rownames(tip_states))
if (length(extra_trait_taxa) > 0 || length(missing_trait_taxa) > 0) {
  stop(
    "Tree/trait taxon mismatch: ",
    length(missing_trait_taxa), " tree tips lack trait data; ",
    length(extra_trait_taxa), " trait rows are absent from the tree."
  )
}
cat("Tree and trait taxa match exactly (", nrow(tip_states), " taxa).\n", sep = "")

cat("Observed joint tip states in trait table:\n")
print(table(paste0(as.character(flagella), as.character(t6ss))))

# Check all trees.
for (i in seq_along(trees)) {
  tr <- trees[[i]]

  if (!ape::is.rooted(tr)) {
    stop("Tree ", i, " is not rooted.")
  }

  if (is.null(tr$edge.length)) {
    stop("Tree ", i, " has no branch lengths.")
  }

  if (any(tr$edge.length < 0)) {
    stop("Tree ", i, " has negative branch lengths.")
  }

  missing_tip_data <- setdiff(tr$tip.label, rownames(tip_states))
  if (length(missing_tip_data) > 0) {
    stop(
      "Tree ", i, " contains tips missing from the traits file: ",
      paste(head(missing_tip_data, 20), collapse = ", "),
      if (length(missing_tip_data) > 20) " ..." else ""
    )
  }
}


# ---------------------------
# SELECT POSTERIOR DRAWS
# ---------------------------

if (is.infinite(n_posterior_draws)) {
  selected_rows <- seq_len(nrow(rates))
} else {
  n_use <- min(as.integer(n_posterior_draws), nrow(rates))
  selected_rows <- sample(seq_len(nrow(rates)), size = n_use, replace = FALSE)
}

# BayesTraits "Tree No" is used if present. If absent, use tree 1.
if ("Tree No" %in% names(rates)) {
  tree_numbers <- as.integer(rates[["Tree No"]])
} else {
  tree_numbers <- rep(1L, nrow(rates))
}

if (any(tree_numbers[selected_rows] < 1L) ||
    any(tree_numbers[selected_rows] > length(trees))) {
  stop(
    "A sampled 'Tree No' is outside the range of trees in tree_file. ",
    "Number of trees read = ", length(trees), "."
  )
}


# ---------------------------
# STOCHASTIC MAPPING
# ---------------------------

allowed_transitions <- data.frame(
  rate = c("q12", "q13", "q21", "q24", "q31", "q34", "q42", "q43"),
  from = c("00", "00", "01", "01", "10", "10", "11", "11"),
  to   = c("01", "10", "00", "11", "00", "11", "01", "10"),
  biological_change = c(
    "T6SS gain | flagella absent",
    "Flagella gain | T6SS absent",
    "T6SS loss | flagella absent",
    "Flagella gain | T6SS present",
    "Flagella loss | T6SS absent",
    "T6SS gain | flagella present",
    "Flagella loss | T6SS present",
    "T6SS loss | flagella present"
  ),
  stringsAsFactors = FALSE
)

n_maps_total <- length(selected_rows) * maps_per_draw

count_rows <- vector("list", n_maps_total)
time_rows  <- vector("list", n_maps_total)

out_i <- 1L

for (k in seq_along(selected_rows)) {

  rate_index <- selected_rows[k]
  rate_row <- rates[rate_index, , drop = FALSE]

  tree_no <- tree_numbers[rate_index]
  tr <- trees[[tree_no]]

  # Reorder trait matrix to this tree.
  X <- tip_states[tr$tip.label, , drop = FALSE]

  Q <- make_Q(rate_row)

  for (m in seq_len(maps_per_draw)) {

    simmap_tree <- phytools::make.simmap(
      tree = tr,
      x = X,
      nsim = 1,
      Q = Q,
      pi = root_prior,
      message = FALSE
    )

    cs <- phytools::countSimmap(
      simmap_tree,
      states = c("00", "01", "10", "11"),
      message = FALSE
    )

    Tr <- cs$Tr

    transition_counts <- setNames(
      integer(nrow(allowed_transitions)),
      allowed_transitions$rate
    )

    for (j in seq_len(nrow(allowed_transitions))) {
      transition_counts[j] <- Tr[
        allowed_transitions$from[j],
        allowed_transitions$to[j]
      ]
    }

    # Total time spent in each joint state over the whole tree.
    dwell <- colSums(simmap_tree$mapped.edge)
    dwell <- dwell[c("00", "01", "10", "11")]

    count_rows[[out_i]] <- data.frame(
      posterior_row = rate_index,
      iteration = if ("Iteration" %in% names(rates)) rate_row[["Iteration"]] else NA,
      tree_no = tree_no,
      map_within_draw = m,
      total_changes = cs$N,
      q12 = transition_counts["q12"],
      q13 = transition_counts["q13"],
      q21 = transition_counts["q21"],
      q24 = transition_counts["q24"],
      q31 = transition_counts["q31"],
      q34 = transition_counts["q34"],
      q42 = transition_counts["q42"],
      q43 = transition_counts["q43"],
      check.names = FALSE
    )

    time_rows[[out_i]] <- data.frame(
      posterior_row = rate_index,
      iteration = if ("Iteration" %in% names(rates)) rate_row[["Iteration"]] else NA,
      tree_no = tree_no,
      map_within_draw = m,
      time_00 = dwell["00"],
      time_01 = dwell["01"],
      time_10 = dwell["10"],
      time_11 = dwell["11"],
      check.names = FALSE
    )

    out_i <- out_i + 1L
  }

  if (k %% 50 == 0 || k == length(selected_rows)) {
    message(
      "Completed posterior draw ",
      k, " / ", length(selected_rows),
      " (", k * maps_per_draw, " maps)"
    )
  }
}


transition_counts <- do.call(rbind, count_rows)
dwell_times <- do.call(rbind, time_rows)

# Remove names inherited from scalar extraction.
rownames(transition_counts) <- NULL
rownames(dwell_times) <- NULL


# ---------------------------
# SUMMARIES
# ---------------------------

count_cols <- c(
  "total_changes",
  "q12", "q13", "q21", "q24",
  "q31", "q34", "q42", "q43"
)

transition_summary <- summarise_columns(
  transition_counts[, count_cols, drop = FALSE]
)

transition_summary <- merge(
  transition_summary,
  rbind(
    data.frame(
      variable = "total_changes",
      from = NA,
      to = NA,
      biological_change = "All allowed changes"
    ),
    data.frame(
      variable = allowed_transitions$rate,
      from = allowed_transitions$from,
      to = allowed_transitions$to,
      biological_change = allowed_transitions$biological_change
    )
  ),
  by = "variable",
  all.x = TRUE,
  sort = FALSE
)

# Restore intuitive ordering after merge.
transition_summary <- transition_summary[
  match(count_cols, transition_summary$variable),
]

dwell_summary <- summarise_columns(
  dwell_times[, c("time_00", "time_01", "time_10", "time_11"), drop = FALSE]
)


# Some useful combined quantities for the biological question.
derived_counts <- data.frame(
  T6SS_gains_flagella_absent  = transition_counts$q12,
  T6SS_gains_flagella_present = transition_counts$q34,
  T6SS_losses_flagella_absent  = transition_counts$q21,
  T6SS_losses_flagella_present = transition_counts$q43,
  flagella_gains_T6SS_absent   = transition_counts$q13,
  flagella_gains_T6SS_present  = transition_counts$q24,
  flagella_losses_T6SS_absent  = transition_counts$q31,
  flagella_losses_T6SS_present = transition_counts$q42
)

derived_summary <- summarise_columns(derived_counts)


# ---------------------------
# WRITE OUTPUT
# ---------------------------

write.csv(
  transition_counts,
  paste0(output_prefix, "_counts_per_map.csv"),
  row.names = FALSE
)

write.csv(
  transition_summary,
  paste0(output_prefix, "_counts_summary.csv"),
  row.names = FALSE
)

write.csv(
  dwell_times,
  paste0(output_prefix, "_dwell_times_per_map.csv"),
  row.names = FALSE
)

write.csv(
  dwell_summary,
  paste0(output_prefix, "_dwell_times_summary.csv"),
  row.names = FALSE
)

write.csv(
  derived_summary,
  paste0(output_prefix, "_biological_summary.csv"),
  row.names = FALSE
)

write.csv(
  allowed_transitions,
  paste0(output_prefix, "_transition_key.csv"),
  row.names = FALSE
)


# ---------------------------
# PRINT MAIN RESULT
# ---------------------------

cat("\nPosterior summary of realised transition counts:\n\n")
print(
  transition_summary[
    ,
    c(
      "variable", "from", "to", "biological_change",
      "mean", "median", "lower_95", "upper_95"
    )
  ],
  row.names = FALSE
)

cat("\nPosterior summary of total time spent in each joint state:\n\n")
print(dwell_summary, row.names = FALSE)

cat("\nFiles written with prefix: ", output_prefix, "\n", sep = "")

