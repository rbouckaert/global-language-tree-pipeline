# This script runs regression models on the global trees where kfolds are performed
# based on the exclusion of specific macroregions

# Input data:
# - datasets_and_trees/ [Directory: Contains input CSV and tree files]
#   - *_2.csv [CSV: Language/environmental data]
#   - *_2.tre [Newick: Phylogenetic trees]

# Output data:
# - regression_results_regions/ [Directory: Regression results]
#   - results_*_regressions25.csv [CSV: Regression results]

library(ape) 
library(nlme) 
library(ade4) 
library(dismo) 
library(picante)
library(INLA)
library(terra)
library(INLA)
library(here)
library(fmesher)
library(spdep)
library(progress)
library(data.table)
library(sp)
library(parallel)
library(stringr)

# Check the number of cores
num_cores <- detectCores()
cat("Number of cores:", num_cores, "\n")

# Set constants
e <- simpleError("test error")
results_folder = "regression_results_regions_fixed"
#datasets_folder = "datasets_and_trees"
datasets_folder = here("datasets_and_trees_fixed_raw_scaled_logged")

# cutD function
cutD <- function(x, n) {
  cut(x,
      breaks = c(quantile(x, probs = seq(0, 1, by = 1 / n), na.rm = T)),
      include.lowest = TRUE
  )
}
# Read in CSV files ----
st0 <- list.files(datasets_folder, pattern = "_2.csv", full.names = TRUE)

# Tree name
st2a <- list.files(datasets_folder, pattern = "_2.csv", full.names = FALSE) 
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
tt0 <- list.files(datasets_folder, pattern = "_2.tre", full.names = TRUE)

# tree name
tt2a <- list.files(datasets_folder, pattern = "_2.tre", full.names = FALSE)
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

# Set INLA to use 4 threads per process
inla.setOption("num.threads", "12")
files1_copy = files1
files1=files1_copy[]

# Start loop ----
for (ii in 1:nrow(files1)) {
  
  sam2 <- 3 
  typeX <- c("nbinomial", "gaussian", "poisson")[sam2]
  
  if (paste("results_", files1[ii, "uni"], "_",
            typeX, "_regressions25.csv", sep = "") %in%
      list.files(results_folder,
                 full.names = FALSE)) {
    next
  }
  
  cat("Starting", files1[ii, "uni"])
  res4 <- NULL
  
  # Find file names ----
  files2 <- files1[ii, ]
  tree1 <- tt[tt$uni == files2$uni, ]
  
  # Read in files ----
  es4 <- fread(as.character(files2$filen), stringsAsFactors = FALSE)
  print(head(es4$y))
  
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
  
  # Macroareas ----------
  es4$macroarea = NA
  
  # Define regions and their indices
  americas_regions <- c(17, 18, 24, 3, 4, 13)
  eurasia_regions <- c(2, 6, 9, 10, 11, 12, 20, 22, 8, 26, 5, 23, 14, 19)
  africa_regions <- c(7, 15, 21, 25)
  oceania_regions <- c(16)
  
  # Update the macroarea column
  es4$macroarea[es4$SUBREGION %in% americas_regions] <- 4 #4 AMERICAS
  es4$macroarea[es4$SUBREGION %in% eurasia_regions] <- 3  #3 EURASIA
  es4$macroarea[es4$SUBREGION %in% africa_regions] <- 1   #1 AFRICA
  es4$macroarea[es4$SUBREGION %in% oceania_regions] <- 2  #2 OCEANIA
  
  kfold_macroareas = data.frame(kf = c(1,2,3,4), macroarea = c("Africa","Oceania","Eurasia","Americas"))
  
  # Update large_regions based on macroareas
  es4$region <- es4$macroarea
  large_regions <- as.data.frame(table(es4$region))
  large_regions <- as.numeric(levels(large_regions[large_regions$Freq >= 30,"Var1"]))
  
  # Make prediction ----
  es4$pred <- NA
  kf<-large_regions
  loos_type="geog"
  es4$kf = es4$region
  
  # Loop for out of sample predictions ----
  for (ww in sort(as.numeric(unique(kf)))) {
    print(ww)
    estr <- copy(es4)
    
    # unscale
    estr$y <- (estr$length)
    type <- "length" 
    
    # to compare against NA'd y predictions
    estr$y2 <- estr$y
    
    # for holdout - NA for 1 macroregion
    estr$y[estr$kf==ww]<-NA
    
    # Make unique names ----
    estr$names <- 1:nrow(estr)
    
    # Make spatial ----
    esqp <- estr
    coordinates(esqp) <- ~ longitude + latitude
    rownames(esqp@data) <- esqp@data$names
    
    # K nearest neighbours ----
    nbs <- knearneigh(coordinates(esqp), k = 10, longlat = T) # k=5 nearest neighbors
    nbs <- knn2nb(nbs, sym = T) #force symmetry!!
    mat <- nb2mat(nbs, style = "B", zero.policy = TRUE)
    rownames(mat) <- esqp$names
    colnames(mat) <- rownames(mat)
    mat <- as.matrix(mat[1:dim(mat)[1], 1:dim(mat)[1]])
    samp1 <- sample(1e10, 1)
    nb2INLA(paste0(here("cl_graph_regions/cl_graph_"), samp1), nbs)
    H <- inla.read.graph(filename = paste0(here("cl_graph_regions/cl_graph_"), samp1))
    
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
    estr2$areaL <- log(estr2$area + abs(min(estr2$area)) + 1)
    estr2$subregion <- estr$names
    
    # Cropland
    model1 <- "include_REGION"
    
    # Simplest
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
    pred1 <- lm.brown$summary.linear.predictor[1:nrow(estr), ]
    pred1$kf = es4$kf
    es4$pred[es4$kf == ww] <- pred1[pred1$kf == ww, 1]
    
    # Calculate Residuals ----
    es4$resid <- estr$y2 - exp(pred1$mean)
    
    # Model Evaluation ----
    # Calculate metrics like pseudo R² and RMSE for model performance
    man1 <- data.frame(obs = NA, pvalue = NA)
    
    if (ww == max(kf)) {
      if(!is.null(es4$y)){
        psudor2 <- cor(exp(es4$pred), estr$y, method = "pearson", use = "complete.obs")
        rmse1 <- sqrt(mean((estr$y - es4$pred)^2, na.rm = TRUE))}
      else {
        psudor2 <- cor(exp(es4$pred), es4$es_lm, method = "pearson", use = "complete.obs")
        rmse1 <- sqrt(mean((es4$es_lm - es4$pred)^2, na.rm = TRUE))
      }
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
      
      # Get x_values
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
    
    #Edit kfold to be named after macroregions
    quads2$kfold <- kfold_macroareas$macroarea[kfold_macroareas$kf==ww]
    
    # Save Prediction Summaries ----
    write.csv(quads2, file = here(results_folder, paste0("results_",files1[ii, "uni"],"_",kfold_macroareas$macroarea[kfold_macroareas$kf==ww],"_slopes25.csv")), row.names = FALSE)
    
    # Finalize Regression Results ----
    res1$n <- NA
    res1$model <- model1
    res1$model2 <- nrow(estr2)
    res1$kfold <- kfold_macroareas$macroarea[kfold_macroareas$kf==ww]
    res1$lastkfold <- max(kf)
    res1$AIC <- lm.brown$waic$waic
    res1$Time <- tot_time - (files1[ii, "V7"] / 1000)
    res1$Crown <- tot_time
    res1$Time2 <- files1[ii, "V7"] * -1
    res1$loop <- ii
    res1$Form <- as.character(formula(formx))[3]
    res1$type <- typeX
    res1$loos_type <- loos_type
    res1$model_name <- NA 
    res1$shapiro <- NA 
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
  write.csv(res4, file = here(results_folder, 
                              paste0("results_", files1[ii, "uni"], "_", typeX, "_regressions25.csv")),
            row.names = FALSE)

} # End of Main Loop ----

