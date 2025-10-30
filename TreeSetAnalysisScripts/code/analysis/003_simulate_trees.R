# ------------------------------------------------------------------------------#
#                              Tree Simulations               ----
# Input data:
# - input_data/edge6636-March-2023-no-metadata.trees [Nexus: Source trees for simulation]
# - input_data/groups_data1.csv [CSV: Group names for simulations]
# - all_trees/ [Directory: Trees for each group, used as templates]

# Output data:
# - simulationtrees/ [Directory: simulated trees and LTT outputs]
#   - simulation_trees_extant_*.tre [Tree: simulated extant trees]
#   - simulation_trees_ltt_*.tre [CSV: LTT statistics]

# Required folders:
# - input_data
# - simulationtrees
# - all_trees

# ------------------------------------------------------------------------------#
#                              Libraries                       ----
# ------------------------------------------------------------------------------#

library(moments)
library(apTreeshape)
library(caper)
library(doParallel)
library(picante)
library(data.table)
library(treespace)
library(phangorn)
library(phytools)
library(phyclust)
library(ggsci)
library(here)
library(BiocManager)
library(ggtree)
library(paleotree)
library(diversitree)
library(ape)
library(data.table)

# ------------------------------------------------------------------------------#
#                             Setup Cluster                      ----
# ------------------------------------------------------------------------------#
## Get group names
all_conts <- read.csv(file = here("input_data", "groups_data1.csv"))[, "x"]

## Make cluster
cl <- makeCluster(10)
registerDoParallel(cl)

# ------------------------------------------------------------------------------#
#                             Start Cluster                      ----
# ------------------------------------------------------------------------------#

foreach(ii = 1:length(all_conts),
                       .packages = c("base", "utils", "ape", "picante", "caper",
                                     "apTreeshape", "moments", "phytools")) %dopar% {
        kk <- ii # sample(1:length(all_conts),1)

  ## Only do if not there
  if (paste(here("simulationtrees", paste("simulation_trees_ltt_",
                                          all_conts[kk], "_1.tre", sep = "")),
                                          sep = "") %in%
      list.files(here("simulationtrees"), full.names = TRUE)) {
    return()
  }

  ## Get Trees
  tt <- list.files(here("all_trees"), pattern = all_conts[kk], full.names = TRUE)

  ### Create n random trees
  trees <- NULL
  maxt <- 50
  for (rr in 1:maxt) {
    treeN <- read.tree(tt[sample(1:length(tt), 1)])

    ntips <- length(treeN$tip.label)

    max_age <- max(node.age(treeN)$ages)

    fit1 <- fit.bd(treeN)

    ## fOr simulated trees
    ntip <- 1
    qu <- 0
    tree2 <- NULL
    while (ntip < ntips & qu < 10) {
      qu <- qu + 1
      tree2 <- pbtree(b = fit1$b, d = fit1$d, n = ntips,
                      nsim = 1, scale = max_age, extant.only = TRUE)
      if (is.null(tree2)) {
        ntip <- 1
      } else {
        ntip <- length(tree2$tip.label)
      }
    }

    if (is.null(tree2)) {
      print("fail")
      next
    }
    tr <- drop.tip(tree2, sample(1:ntip, ntip - ntips, replace = FALSE))

    ## Rename
    treeN <- ladderize(treeN)
    tr <- ladderize(tr)
    tr$tip.label <- treeN$tip.label

    rnum <- sample(1e6, 1)

    ## Write tree
    if (rr == maxt) {
      rnum <- 1
      write.tree(tr, file = here("simulationtrees", 
                                  paste("simulation_trees_extant_", 
                                  all_conts[kk], "_1.tre", sep = "")))
        } else {
      write.tree(tr, file = here("simulationtrees",
                                 paste("simulation_trees_extant_", 
                                 all_conts[kk], "_", rnum, ".tre", sep = "")))
        }

    ## Ltt
    ltt1 <- ltt(tr, plot = FALSE)

    ## Ltt times
    ltt1$times <- ltt1$times - max(ltt1$times)
    ltt1$ltt2 <- log(ltt1$ltt)

    ## Results table
    resx1 <- data.frame(tree_id_1_901 = rr, cont = all_conts[kk], ltt = ltt1$ltt,
                        logltt = ltt1$ltt2, times = ltt1$times, gamma = ltt1$gamma,
                         p = ltt1$p, max_height = max(nodeHeights(tr)), type = "simulation")

    ## Write table
    write.csv(resx1, file = here("simulationtrees",
                                paste("simulation_trees_ltt_", 
                                all_conts[kk], "_", rnum, ".tre", sep = "")))

    print(rr)
  } ## End of rr loop

} ## End ii trees

stopCluster(cl)
