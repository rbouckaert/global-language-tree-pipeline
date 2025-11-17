# -------------------------------------------------------------------------
# Input data:
# - input_data/global-language-trees-6636-taxa.trees [Nexus: Source trees]
#
# Output data:
# - outputs/new_trees_sliced/ [Directory: clade counts per time slice]
#   - *_clade_counts.csv [CSV: number of clades at each time slice per tree]

#   - Only uses 100 randomly selected Global trees from all_trees
#   - Cuts each tree at specified time points (250, 500, ..., 3500, 4250,
#     5000, ..., 10000 years ago)
#   - Does NOT calculate clade stats
#   - Tracks the number of clades at each time slice
# -------------------------------------------------------------------------

# Load necessary libraries
library(ape)          # for read.tree, treeSlice
library(picante)      # for node.age
library(data.table)   # optional: for data management
library(here)         # optional: for file paths
library(progress)     # optional: for progress bar
library(phytools)

# 1) Load Trees from Nexus File -------------------------------------------

# Define the path to the Nexus file
trees_file <- here("input_data", "global-language-trees-6636-taxa.trees")

# Read all trees from the Nexus file
cat("Reading trees from file:", trees_file, "\n")
all_trees <- read.nexus(trees_file)

# Number of trees
total_trees <- length(all_trees)
cat("Number of trees loaded:", total_trees, "\n")

# 2) Randomly Select 100 Trees --------------------------------------------
set.seed(123)  # For reproducibility
selected_indices <- sample(seq_along(all_trees), 100, replace = FALSE)
selected_trees <- all_trees[selected_indices]
selected_names <- names(all_trees)[selected_indices]  # Get the corresponding names
selected_names <- gsub("STATE", "tree", selected_names)  # Replace spaces with underscores

cat("Selected 100 random trees for processing.\n")

# 3) Define Time Slices ----------------------------------------------------

time_points_years <- c(seq(250, 10000, by = 250))
time_slices <- time_points_years / 1000  # convert to 'tree' units (if needed)

# 4) Initialize Progress Bar ----------------------------------------------

pb <- progress_bar$new(
  format = "[:bar] :percent ETA: :eta | Elapsed: :elapsed mins",
  total = length(selected_trees),  # Only for 100 trees
  clear = FALSE,
  width = 60
)

# 5) Process Each Selected Tree -------------------------------------------
for (i in seq_along(selected_trees)) {
  start_time <- Sys.time()  # Start timer for current tree
  
  # Get tree and its name
  tree <- selected_trees[[i]]
  tree_name <- selected_names[i] 
  
  # Define output file for this tree
  output_file <- here("outputs", "new_trees_sliced", paste0(tree_name, "_clade_counts.csv"))
  
  # Skip this tree if results already exist
  if (file.exists(output_file)) {
    cat("Skipping:", tree_name, "- already processed.\n")
    next
  }
  
  # Compute node ages
  node_ages <- picante::node.age(tree)$ages
  tot_time  <- max(node_ages)  # Root age in tree units (1000 years)
  
  # Initialize results for this tree
  tree_results <- list()
  
  # Process each time slice
  for (ts in time_slices) {
    # Calculate the cut point in tree time
    cut_point <- tot_time - ts
    
    # Slice the tree
    clades <- treeSlice(tree, cut_point, orientation = "tipwards")

    # Count number of subclades
    num_clades <- length(clades)
    
    # Measure sub-clade sizes
    sub_clade_sizes <- sapply(clades, ape::Ntip)
    var_size <- var(sub_clade_sizes, na.rm = TRUE)  # Variance of clade sizes
    median_size <- median(sub_clade_sizes, na.rm = TRUE)  # Median size
    
    # Store results
    tree_results[[length(tree_results) + 1]] <- data.frame(
      tot_time           = tot_time,
      time_slice         = ts * 1000,   # Convert back to years
      cut_point          = cut_point,
      num_clades         = num_clades,
      var_clade_size     = var_size,
      median_clade_size  = median_size
    )
  }
  
  # Combine results into a data frame
  tree_df <- data.table::rbindlist(tree_results)
  
  # Save results for this tree
  write.csv(tree_df, output_file, row.names = FALSE)
  cat("Saved results for:", tree_name, "\n")
  
  # Print elapsed time in minutes
  elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  cat(sprintf("Tree %d (%s) completed in %.2f minutes\n", i, tree_name, elapsed_time))
}

cat("\nAll selected trees processed!\n")