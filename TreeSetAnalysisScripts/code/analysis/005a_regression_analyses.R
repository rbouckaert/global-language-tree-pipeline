# ------------------------------------------------------------------------------#
#                          Data Preparation ------
# ------------------------------------------------------------------------------#
# apply INLA regression to language trees and environmental data
# D.Redding 14/07/2023

# Input data:
# - datasets_and_trees/ [Directory: Contains input CSV and tree files]
#   - *_2.csv [CSV: Language/environmental data]
#   - *_2.tre [Newick: Phylogenetic trees]

# Output data:
# - regression_results_fin_log/ [Directory: Regression results]
#   - results_*_regressions25.csv [CSV: Regression results]
#   - results_*_slopes25.csv [CSV: Slope summaries]

# ------------------------------------------------------------------------------#
#                          Data Preparation ------
# ------------------------------------------------------------------------------#

pacman::p_load("ape","nlme","ade4","dismo","picante","INLA","spdep","here",
               "progress","sp","tidyverse","data.table")

#install.packages("INLA")
e <- simpleError("test error")

# cutD function
cutD <- function(x, n) {
  cut(x,
      breaks = c(quantile(x, probs = seq(0, 1, by = 1 / n), na.rm = T)),
      include.lowest = TRUE
  )
}

# Set up folder paths
data_trees = here("datasets_and_trees_fixed_scaled_raw_scaled_logged")

# Read in CSV files ----
st0 <- list.files(data_trees, pattern = "_2.csv", full.names = TRUE)

# Tree name
st2a <- list.files(data_trees, pattern = "_2.csv", full.names = FALSE) 
st2b <- gsub(".csv", "", st2a)
st2b <- gsub("alltime", "all_time", st2b, fixed = TRUE)
st3a <- read.table(text = st2b, sep = "_", stringsAsFactors = FALSE)

# Create input file ----
files1 <- data.frame(st3a, filen = st0, uni = st2b, stringsAsFactors = FALSE)
# time column
files1[, "V7"] <- as.numeric(gsub("XXX", "", files1[, "V7"]))
# simulated trees
files1$V3[files1$V3 == "extant"] <- (-1)

# Read in trees ----
tt0 <- list.files(data_trees, pattern = "_2.tre", full.names = TRUE)
tt2a <- list.files(data_trees, pattern = "_2.tre", full.names = FALSE)
tt2a <- gsub(".tre", "", tt2a, fixed = TRUE)
tt2a <- gsub("alltime", "all_time", tt2a, fixed = TRUE)
tt3a <- read.table(text = tt2a, sep = "_", stringsAsFactors = FALSE)
tt3a$V6 <- gsub("time", "", tt3a$V6, fixed = TRUE)
tt3a$V6 <- as.numeric(tt3a$V6)
tt <- data.frame(tt3a, filen = tt0, uni = tt2a, stringsAsFactors = FALSE)
# time column
tt[, "V7"] <- as.numeric(gsub("XXX", "", tt[, "V7"]))
# simulated trees
tt$V3[tt$V3 == "extant"] <- (-1)

# Check how many trees
length(unique(tt$V3)) - 1

# Reorder sample ----
sam1 <- sample(1:nrow(files1), nrow(files1), replace = FALSE) 
files1 <- files1[sam1, ] # combine together and reorder

# ------------------------------------------------------------------------------#
#                          Regression Analysis ------
# ------------------------------------------------------------------------------#

# Initialize the progress bar
pb <- progress_bar$new(
  format = "  Processing [:bar] :percent in :elapsed, ETA: :eta",
  total = nrow(files1), clear = FALSE, width=60
)

# Initialize a vector to store iteration times
iteration_times <- numeric(nrow(files1))

# Set INLA to use 2 threads for computations
inla.setOption("num.threads", "4:1")

# Start loop ----
for (ii in 1:nrow(files1)) {
  # Progress bar and iteration time
  pb$tick()  # Update the progress bar
  start_time <- Sys.time()
  
  sam2 <- 3
  typeX <- c("nbinomial", "gaussian", "poisson")[sam2]
  
  if (paste("results_", files1[ii, "uni"], "_",
            typeX, "_regressions25.csv", sep = "") %in%
      list.files(here("regression_results_fin_log"),
                 full.names = FALSE)) {
    next
  }
  res4 <- NULL
  
  # Find file names ----
  files2 <- files1[ii, ]
  tree1 <- tt[tt$uni == files2$uni, ]
  
  # Read in files ----
  es4 <- fread(as.character(files2$filen), stringsAsFactors = FALSE)
  
  if (files1[ii, "V2"] == "simulatione") {
    es4$order <- sample(1:nrow(es4))
    es4$countup <- 1
    for (uu in 2:nrow(es4)) {
      if (es4$region[uu] == es4$region[uu - 1]) {
        es4$countup[uu] <- es4$countup[uu - 1]
      } else {
        es4$countup[uu] <- es4$countup[uu - 1] + 1
      }
    }
    es4 <- es4[order(es4$countup, es4$order), ]
  }
  
  row.names(es4) <- as.character(es4$tips)
  tr1 <- read.tree(file = as.character(tree1$filen))
  
  # Calculate node ages ----
  node_ages <- picante::node.age(tr1)$ages
  tot_time <- max(node_ages)
  
  # Regions hold out ----
  large_regions <- as.data.frame(table(es4$region))
  large_regions <- as.numeric(levels(large_regions[large_regions$Freq >= 30, "Var1"]))
  
  # Make prediction ----
  es4$pred <- NA
  
  # No hold out ----
  kf <- 1
  loos_type <- "none"
  
  # Loop for out of sample predictions ----
  for (ww in sort(as.numeric(unique(kf)))) {
    
    estr <- copy(es4)
    estr$y <- (estr$length)
    type <- "length" 
    
    # To compare against NA'd y predictions
    estr$y2 <- estr$y
    
    # Make unique names ----
    estr$names <- 1:nrow(estr) 
    
    # Make spatial ----
    esqp <- estr
    coordinates(esqp) <- ~ longitude + latitude
    rownames(esqp@data) <- esqp@data$names
    
    # K nearest neighbours ----
    nbs <- knearneigh(coordinates(esqp), k = 10, longlat = T) # k=5 nearest neighbors
    nbs <- knn2nb(nbs, sym = T) 
    mat <- nb2mat(nbs, style = "B", zero.policy = TRUE)
    rownames(mat) <- esqp$names
    colnames(mat) <- rownames(mat)
    mat <- as.matrix(mat[1:dim(mat)[1], 1:dim(mat)[1]])
    samp1 <- sample(1e10, 1)
    nb2INLA(paste0(here("cl_graph/cl_graph_"), samp1), nbs)
    H <- inla.read.graph(filename = paste0(here("cl_graph/cl_graph_"), samp1))
    
    # Make tree VCV matrix ----
    phylo_covar_mat <- ape::vcv(tr1)
    phylo_covar_mat <- phylo_covar_mat / max(phylo_covar_mat)
    phylo_covar_mat <- phylo_covar_mat / exp(determinant(phylo_covar_mat)$modulus[1] /
                                               nrow(phylo_covar_mat))
    phylo_prec_mat <- solve(phylo_covar_mat)
    
    pcprior <- list(prec = list(prior = "pc.prec", param = c(0.36, 0.1)))
    
    estr2 <- cbind(estr, data.frame(
      species = rownames(phylo_prec_mat),
      phy_id = 1:nrow(phylo_prec_mat), stringsAsFactors = FALSE
    ))
    
    b <- 15
    estr2$clade_age_g <- INLA::inla.group(estr2$clade_age, n = b)
    estr2$island2 <- estr2$island
    estr2$island2[estr2$island > 0.5] <- "islands"
    estr2$island2[estr2$island < 0.5] <- "mainland"
    estr2$subregion <- estr$names
  
    model1 <- "include_REGION"
    
    # Simplest Updated
    formx <- formula(paste0("y~  island + area_log + distancetocityyear2 + popd_log + friction_log + f(clade_age_g, model='rw2', constr = TRUE, scale.model = TRUE)", "+  f(subregion, model = 'bym2', graph = H,constr=TRUE,scale.model=TRUE) + f(phy_id, model = 'generic0', Cmatrix = phylo_prec_mat,constr = TRUE, hyper = pcprior)"))
    
    # For trycatch
    e <- simpleError("test error")
    
    # Type of model
    if (typeX == "gaussian") {
      estr2$y <- log(estr2$y)
    }
    
    # Model Fitting ----
    lm.brown <- tryCatch(
      inla(formx,
           data = estr2,
           control.fixed = list(
             mean = 0, # mean for betas
             prec = 1 / (0.3^2), # precision for intercept: sd = 0.6 --> precision =1/variance --> 1/(0.6^2) = 2.777778
             mean.intercept = 0, # mean for intercept
             prec.intercept = 1 / (0.6^2)
           ),
           family = typeX,
           control.family = list(list(link = "default")),
           control.compute = list(config = TRUE, dic = FALSE, waic = TRUE, cpo = TRUE),
           control.predictor = list(compute = TRUE, link = 1)
      ),
      error = function(e) e
    )
    
    # Check for Errors ----
    if (class(lm.brown)[1] == "simpleError") {
      next
    }
    
    # Post-Processing ----
    # Generate Predictions ----
    # Extract and store predicted values from the model
    pred1 <- lm.brown$summary.linear.predictor[1:nrow(estr), ]
    
    # predicted should go on x<- Tim Lucas
    es4$pred[kf == ww] <- pred1[kf == ww, 1]
    
    # Calculate Residuals ----
    # Compute residuals to assess model fit
    es4$resid <- estr$y2 - exp(pred1$mean)
    man1 <- data.frame(obs = NA, pvalue = NA)
    
    # Model Evaluation ----
    # Calculate metrics like pseudo R² and RMSE for model performance
    man1 <- data.frame(obs = NA, pvalue = NA)
    if (ww == max(kf)) {
      psudor2 <- cor(exp(es4$pred), estr$y, method = "pearson", use = "complete.obs")
      rmse1 <- sqrt(mean((estr$y - es4$pred)^2, na.rm = TRUE))
    } else {
      psudor2 <- NA
      rmse1 <- NA
      
      # Perform Mantel Test ----
      station.dists <- dist(cbind(es4$longitude, es4$latitude))
      ozone.dists <- dist(abs(es4$resid))
      man1 <- ade4::mantel.rtest(station.dists, ozone.dists, nrepet = 999)
    }
    
    # Summarize Regression Results ----
    # Compile fixed effects estimates and significance indicators
    res1 <- lm.brown$summary.fixed
    res1$name <- row.names(res1)
    
    res1$sig1 <- 0
    res1$sig1[(res1$`0.025quant` < 0 & res1$`0.975quant` < 0) | 
                (res1$`0.025quant` > 0 & res1$`0.975quant` > 0)] <- 1
    res1$mean_sig <- res1$mean * res1$sig1
    
    # Generate Predictions for Variables ----
    colsX <- c("island", "area_log", "distancetocityyear2", "popd_log", "friction_log")
    colsL <- c("area_log", "popd_log", "friction_log")
  
    
    for (kk in 1:length(colsX)) {
      #print(colsX[kk])
      col1 <- colsX[kk]
      if (!col1 %in% res1$name) {
        next
      }
      
      # removeL prefix from column names
      col2 = gsub("_log", "", col1, fixed = TRUE)
      
      # add term raw
      col2 = paste0(col2, "_raw")
      
      # get raw values
      range1 = estr2[,..col2]
      x_values <- seq(min(range1, na.rm = TRUE), # need to define nr of values to output for raw scale
                      max(range1, na.rm = TRUE), length.out = 1000)
      
      #Make predictions
      predscit = (res1$mean[res1$name == col1] * x_values) + res1$mean[res1$name == "(Intercept)"]
      predscit2 = (res1$mean_sig[res1$name == col1] * x_values) + res1$mean_sig[res1$name == "(Intercept)"]
      
      #Make data frame
      quads1 <- data.frame(
        col = col1,
        x = x_values,
        y = predscit,
        y2 = predscit2,
        yearbp = files1[ii, "V7"],
        tree = files1[ii, "V3"],
        kfold = ww,
        stringsAsFactors = F
      )
      
      if (kk == 1) {
        quads2 <- quads1
      } else {
        quads2 <- rbind(quads2, quads1)
      }
    }
    
    # Save Prediction Summaries ----
    write.csv(quads2, file = paste(here("regression_results_fin_log/results_"), 
                                   files1[ii, "uni"], "_slopes25.csv", sep = ""), 
              row.names = FALSE)
    
    # Finalize Regression Results ----
    res1$n <- NA # length(models[[qqq]])
    res1$model <- model1 # ref1[qqq]
    res1$model2 <- nrow(estr2) # ref2[qqq]
    res1$kfold <- ww
    res1$lastkfold <- max(kf)
    res1$AIC <- lm.brown$waic$waic
    res1$Time <- tot_time - (files1[ii, "V7"] / 1000)
    res1$Crown <- tot_time
    res1$Time2 <- files1[ii, "V7"] * -1
    res1$loop <- ii
    res1$Form <- as.character(formula(formx))[3]
    res1$type <- typeX
    res1$loos_type <- loos_type
    res1$model_name <- NA # mod_names[[ref1[qqq]]]
    res1$shapiro <- NA # round(shapiro.test(lm.brown$residuals[,1])$p.value,3)
    res1$mantcor <- man1$obs
    res1$mantel <- man1$pvalue
    res1$psudor2 <- psudor2
    res1$rmse <- rmse1
    res1$tree_type <- files1[ii, "V2"]
    res1$tree <- files1[ii, "V3"]
    res1$group <- files1[ii, "V4"]
    
    # Combine and Store Results ----
    if (is.null(res4)) {
      
      res4 <- res1
    } else {
      res4 <- rbind(res4, res1)
    }
    # Clean Up ----
    rm(res1, lm.brown)
    # End of Out-of-Sample Predictions Loop ----
  } 
  
  # Save Loop Results ----
  # After completing all folds, save the aggregated results to a CSV file
  print(ii)
  write.csv(res4, file = here("regression_results_fin_log", 
                              paste("results_", files1[ii, "uni"], "_", typeX, "_regressions25.csv", 
                                    sep = "")),
            row.names = FALSE)
  # Record the time taken
  end_time <- Sys.time()
  iteration_times[ii] <- as.numeric(difftime(end_time, start_time, units = "secs"))
  print(paste("Iteration", ii, "took", round(iteration_times[ii], 2), "seconds"))
} # End of Main Loop ----

