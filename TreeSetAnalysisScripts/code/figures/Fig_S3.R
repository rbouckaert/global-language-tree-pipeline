# ------------------------------------------------------------------
# Figure S3: Clade Counts and Variance over Time
# ------------------------------------------------------------------

library(data.table)
library(ggplot2)
library(here)
library(coda)

treeset <- "new"

# 1) Set the path and pattern based on treeset
data_path   <- here('outputs', 'new_trees_sliced')
file_pattern <- "_clade_counts\\.csv$"

# 2) List all relevant CSV files
csv_files <- list.files(
  path       = data_path,
  pattern    = file_pattern,
  full.names = TRUE
)

# 3) Create an empty list to store data from each file
all_list <- list()

# 4) Loop over files and read data
for (f in csv_files) {
  
  # Read the CSV
  tmp <- fread(f)
  
  # Safety check: skip if file is empty
  if (nrow(tmp) == 0) next
  
  # Extract a 'treeID' from the filename
  base_name <- basename(f)  
  parts <- strsplit(base_name, "_")[[1]]
  
  # Parse differently depending on old or new naming scheme
    # New: c("tree", "937", "clade", "counts.csv")
    tree_num <- parts[2]
    treeID   <- paste0("tree_", tree_num)
  
  # We keep one row per time slice
  dt_file <- data.table(
    treeID         = treeID,
    timeBP         = tmp$time_slice,      # The 'time_slice' column
    n_clades       = tmp$num_clades,      # The 'num_clades' column
    var_clade_size = tmp$var_clade_size   # The 'var_clade_size' column
  )
  
  # Append to list
  all_list[[f]] <- dt_file
}

# 5) Combine everything
df_all <- rbindlist(all_list, use.names = TRUE, fill = TRUE)

# 6) Create a summary data frame that groups by timeBP across all treeIDs
# Calculate summary statistics including 95% HPDI
df_summary <- df_all[,
                     {
                       # Convert n_clades to an MCMC object
                       clades_mcmc <- as.mcmc(n_clades)
                       
                       # Calculate 95% HPDI
                       hpd_clades <- HPDinterval(clades_mcmc, prob = 0.95)
                       
                       # Return the statistics
                       .(
                         mean_clades    = mean(n_clades, na.rm = TRUE),
                         median_clades  = median(n_clades, na.rm = TRUE),
                         sd_clades      = sd(n_clades, na.rm = TRUE),
                         hpd_lower      = hpd_clades[1],  # Lower bound of HPDI
                         hpd_upper      = hpd_clades[2],  # Upper bound of HPDI
                         mean_var       = mean(var_clade_size, na.rm = TRUE),
                         sd_var         = sd(var_clade_size, na.rm = TRUE),
                         count_trees    = .N  # Number of trees at this timeBP
                       )
                     },
                     by = timeBP]

# Sort by timeBP
setkey(df_summary, timeBP)

# ------------------------------------------------------------------
# Plotting
# ------------------------------------------------------------------

# Define a ratio for a dual y-axis
ratio <- 60000 / 1500  # example ratio

# Find min and max of timeBP for the scale
x_min <- min(df_all$timeBP, na.rm = TRUE)
x_max <- max(df_all$timeBP, na.rm = TRUE)
y_max_clades = max(df_all$n_clades, na.rm = TRUE)
  
# Build the plot
p <- ggplot() +
  # (A) Spaghetti lines for n_clades
  geom_line(
    data = df_all,
    aes(x = timeBP, y = n_clades, group = treeID),
    color = "grey40",
    alpha = 0.2
  ) +
  # (B) Spaghetti lines for var_clade_size, scaled
  geom_line(
    data = df_all,
    aes(x = timeBP, y = var_clade_size / ratio, group = treeID),
    color = "red",
    alpha = 0.2
  ) +
  # (C) Mean lines
  geom_line(
    data = df_summary,
    aes(x = timeBP, y = mean_clades),
    color = "grey20",
    size = 1
  ) +
  geom_line(
    data = df_summary,
    aes(x = timeBP, y = mean_var / ratio),
    color = "red",
    size = 1
  ) +
  # (D) Optional: vertical dashed lines at specific time points
  geom_vline(
    xintercept = c(3500, 4250, 5000),
    linetype = "dashed",
    color = "black"
  ) +
  geom_hline(
    yintercept = 0,                # Draw line at y = 0
    linetype = "solid",            # Solid line style
    color = "grey20",               # Black color
    size = 0.5,
    alpha = 0.5# Line thickness
  ) +
  # (F) Y-axes
  scale_x_continuous(
    limits = c(1000, x_max),             # Set axis limits
    breaks = c(1000,2500,3500,4250,5000,7500,10000), # Choose appropriate breaks
    name = "Years before present"          # X-axis label
  ) +
  scale_y_continuous(
    name = "Number of Clades",
    breaks = c(0, 100, 250, 500, 750, 1000, 1500),
    sec.axis = sec_axis(
      ~ . * ratio,
      name = "Variance"
    )
  ) +
  theme_classic(base_size = 14) +
  # Style secondary Y-axis title in red
  theme(
    axis.title.y.right = element_text(color = "red"), # Red label for Variance axis
    axis.text.y.right  = element_text(color = "red")   # Tick numbers in red
  ) 

print(p)

# Dynamically generate the filename based on treeset
file_prefix <- ("new_trees")

# Save the plot as a high-quality PNG
ggsave(
  
  filename = here("outputs","Fig_S3","Fig_S3.png"),
  plot = p,
  width = 15.92 / 2.54,           # Convert cm to inches
  height = 9.56 / 2.54,           # Convert cm to inches
  dpi = 600,                       # High resolution
  scale = 1.2
)

# Save the plot as a high-quality PDF
ggsave(
  filename = here("outputs","Fig_S3","Fig_S3.pdf"),
  plot = p,
  width = 15.92 / 2.54,
  height = 9.56 / 2.54,
  device = cairo_pdf,              # High-quality PDF output
  scale = 1.2
)