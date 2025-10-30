# provide phylogenetic tree and dataframe with two columns: 
# the first comprising species names, the second comprising their associated GE2 scores (between 0 and 1)
# function returns three objects: 
# 1. dataframe with terminal branch length, GE2, ED2 and EDGE2 scores for each species
# 2. expected PD loss tree
# 3. PD and expected PD loss in MY for the clade

EDGE.2.calc <- function(tree, pext){
  require(phylobase)
  require(data.table)
  if(!class(tree) == "phylo"){
    tree <- as(tree, "phylo")
  }
  tree_dat <- data.frame(Species = as.character(unique(tree$tip.label)),
                         TBL = tree$edge.length[sapply(c(1:length(tree$tip.label)),
                                                       function(x,y) which(y==x),y=tree$edge[,2])], 
                         pext = NA, ED = NA, EDGE = NA)
  ePD.dat <- data.frame(PD = sum(tree$edge.length),ePDloss = NA)
  tree <- as(tree, "phylo4")
  names(pext) <- c("species","pext")
  for(i in 1:length(tree_dat$Species)){
    tree_dat$pext[i] <- pext$pext[pext$species == tree_dat$Species[i]]
  }
  nodes <- descendants(tree, rootNode(tree), "all")
  for(i in 1:length(nodes)){
    tips <- descendants(tree, nodes[i], "tips")
    tips <- names(tips)
    tipscores <- which(pext$species %in% tips)
    tree@edge.length[which(tree@edge[,2] == nodes[i])] <- edgeLength(tree, nodes[i])*prod(pext$pext[tipscores])
  }
  for(i in 1:length(tree_dat$Species)){
    tree_dat$EDGE[i] <- sum(tree@edge.length[which(tree@edge[,2] %in% ancestors(tree, 
                                                                                which(tipLabels(tree) == tree_dat$Species[i]), "ALL"))], na.rm=T)
    tree_dat$ED[i] <- tree_dat$EDGE[i] / tree_dat$pext[i] 
  }
  tree <- as(tree, "phylo")
  ePD.dat$ePDloss <- sum(tree$edge.length)
  edge.res <- list(tree_dat,tree,ePD.dat)
  return(edge.res)
}
