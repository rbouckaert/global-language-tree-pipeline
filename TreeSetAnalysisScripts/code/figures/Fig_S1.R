  ################################################################################
  # Figure S1: Language Family Diversification Rates, Imbalance, and Gamma
  ################################################################################
  # Files and folders read in this script:
  # 1. outputs (folder containing summary files)

  # Setup --------------------------------------------------------------------
  # Load required packages
  pacman::p_load(ggplot2, cowplot, ggpubr, gridExtra, RColorBrewer, ggridges, data.table, forcats, here)
  
  # Data Import --------------------------------------------------------------------
  # Read in summary file of posterior trees
  fil3 = read.csv(here("outputs", "All_Trees_Summary.csv"), 
                  stringsAsFactors = FALSE)

  # Rename and reorder real data
  fil3a <- fil3
  rm(fil3)
  fil3a <- fil3a[order(sort(fil3a$meanDR)),]
  fil3a$type = "real"

  ## Simulated data summary files
  lf1 <- list.files(here("outputs"), pattern="summary_statsSIM_2", full.names = TRUE)
  lf2 <- list.files(here("outputs"), pattern="summary_statsSIM_2", full.names = FALSE)
  lf2 <- gsub(".csv", "", lf2)
  lf2 <- gsub("SIM_", ";", lf2)
  lf3 <- read.table(text=lf2, sep=";")
  lf3$filen <- lf1
  lf3$V2 <- as.Date(lf3$V2)

  # Read in latest file (simulated data)
  fil3 <- fread(lf3[lf3$V2==(sort(lf3$V2, decreasing = TRUE)[1]), "filen"]) 
  fil3b <- fil3
  rm(fil3)
  fil3b$type = "simulated"

  # Data Processing --------------------------------------------------------------------

  # Merge real and simulated data
  fil3a$num <- NA
  common_cols <- intersect(names(fil3a), names(fil3b))
  fil3a <- fil3a[, common_cols]
  fil3b = as.data.frame(fil3b)
  fil3b <- fil3b[, common_cols]
  fil3 <- rbind(fil3a, fil3b)
  fil3 <- as.data.frame(fil3)

  # Add grouping column
  grps1 <- read.csv(here("input_data", "groups_table.csv"))
  names(grps1)[3] <- "group2"
  fil3 <- merge(fil3, grps1[, c("group", "group2")], by="group", all.x=TRUE, all.y=FALSE)

  # Check if all rows matched
  nrow(fil3[is.na(fil3$group2),])

  # Clean up language family names
  fil3$name <- gsub("_", " ", fil3$group)
  fil3$name <- gsub("Bantu (Atlantic-Congo)", "Bantu-Atlantic-Congo", fil3$name)
  fil3$name <- gsub("Central-Sudanic", "Central Sudanic", fil3$name)
  fil3$name <- gsub("North-America", "North America", fil3$name)
  fil3$name <- gsub("South-America", "South America", fil3$name)
  fil3$name <- gsub("Nuclear-Trans-New-Guinea", "Nuclear Trans New Guinea", fil3$name)
  fil3$name <- gsub("Nuclear-Torricelli", "Nuclear Torricelli", fil3$name)

# Visualizations --------------------------------------------------------------------
 # Figure S1: Imbalance and Gamma

  # Prepare data for imbalance and gamma plots
  filx <- fil3[fil3$group2 %in% c("macroarea", "all") & fil3$type=="real" ,]
  filx$gammaX <- as.numeric(filx$gammma)
  fil3b$gammaX <- as.numeric(fil3b$gammma)

  # Imbalance plot
  bal1 <- ggplot(data=filx, aes(x=(balance/richness), y=fct_reorder(name, (balance/richness)), fill=name)) +
    geom_density_ridges(bandwidth = 0.5) +
    geom_density_ridges(data=fil3[fil3$name %in% unique(filx$name) & fil3$type!="real", ], 
                        aes(x=(balance/richness), y=name, fill=name), 
                        bandwidth = 0.25, fill="light grey", alpha=0.5) +
    scale_fill_manual(values = c('#7570b350', '#e7298a50', '#00000050', '#66a61e50', '#d95f0250', '#e6ab0250')) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_cowplot(20) +
    xlab("Imbalance") +
    theme(axis.title.y = element_blank(), legend.position = "none") +
    ggtitle(label="A")

  # Gamma plot
  gamma1 <- ggplot(data=filx, aes(x=(gammaX), y=fct_reorder(name, (gammaX)), fill=name)) +
    geom_density_ridges(bandwidth = 1) +
    scale_fill_manual(values = c('#7570b350', '#e7298a50', '#00000050', '#66a61e50', '#d95f0250', '#e6ab0250')) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_cowplot(20) +
    xlab("Gamma") +
    theme(axis.title.y = element_blank(), legend.position = "none") +
    ggtitle(label="B")

  # Combine plots
  combined_plot <- grid.arrange(bal1, gamma1, nrow=1)

  # Save the combined plot
  ggsave(filename = here("outputs","Fig_S1" ,"Fig_S1.pdf"), 
        plot = combined_plot, 
        height = 5, 
        width = 12)