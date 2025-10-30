# Script: 013_language_analysis_bisse12b.R
# Description: Runs BiSSE (Bifurcating Speciation and Extinction) and HiSSE (Heterogeneous
# Speciation and Extinction) analyses on phylogenetic trees.

#------------------------------------------------------------------------------

# Input data:
# - code/AUX2_traitDependent_functions.R [R script: Trait-dependent functions]
# - input_data/bisse_data_in.csv [CSV: BiSSE input data]
# - all_trees/ [Directory: Phylogenetic trees]
# - input_data/tips_by_continent.csv [CSV: Tips by continent]

# Output data:
# - BISSE_runs/ [Directory: BiSSE analysis results]
#   - es4_*_BISSE_runs3.csv [CSV: BiSSE run outputs]

#------------------------------------------------------------------------------
# Load required libraries

library(caper)
library(hisse)
library(ape)
library(data.table)
library(here) # Added here package for path management

# CutD function
source(here("code", "AUX2_traitDependent_functions.R"))

# Read data file
res5 <- fread(file = here("input_data", "bisse_data_in.csv"))

# Correct ossetic
res5$Glottocode [ res5$Glottocode =="osse1243"]<-"iron1242" 

# Correct masa1300
res5$EA033_Extended[res5$Glottocode =="masa1300"]<-1

# Combine forager
res5$Forager<-res5$Forager_language_definite
res5$Forager[res5$Forager_language_possible==1]<-1

# Choose columns
cols1<-c("Marriage_Exogamous","Forager","EA033_Extended")#names(res5)[c(9,11,106,107)]

# Read in trees
tt0 <- list.files(here("all_trees"), full.names = TRUE)

# Tree name
tt2a <- list.files(here("all_trees"), full.names = FALSE)
tt2a <- gsub(".tree", "", tt2a)
tt2a <- gsub("lexicon_families", "lexiconfamilies", tt2a)
tt3a <- read.table(text = tt2a, sep = "_")


# Combine together
tt<-c(tt0)
tt2<-c(tt2a)
tt3<-rbind(tt3a)

# Only need global trees
tt<-tt[tt3$V4 %in% c("all")]
tt2<-tt2[tt3$V4 %in% c("all")]
tt3<-tt3[tt3$V4 %in% c("all"),]

# Reorder sample
sam1<-sample(1:(length(tt)),replace=FALSE)

# Combine together and reorder
tt<-tt[sam1]
tt3<-tt3[sam1,]
tt2<-tt2[sam1]

# One to drop
drop_cont<-read.csv(file=here("input_data", "tips_by_continent.csv"))

# Start loop
ii=1
for( ii in ii:length(tt)){
  res2<-NULL
  kk=1
  
  # Location
  regs1<-sample(unique(drop_cont$group),1)
  
  # Check if done
  if(paste("es4_", tt2[[ii]], "_", regs1, "_BISSE_runs3.csv", sep = "") %in% 
     list.files(here("BISSE_runs"), full.names = FALSE)) {next}
  
  # Read tree
  tree<-ape::read.tree(tt[ii])
  
  # Remove macroarea tips
  if(!regs1=="Global"){tree1<-tree; tree<-ape::drop.tip(tree1,drop_cont[drop_cont$group==regs1,"tip.labels"])} 

  # Check
  if(is.null(tree)){break}
  
  # Cut down data frame
  res6<-as.data.frame(res5[res5$Glottocode %in% tree$tip.label,])
  if(nrow(res6) != length(tree$tip.label)){print("missing data");break}
  rownames(res6)<-res6$Glottocode
  res6<-res6[order(match(res6$Glottocode,tree$tip.label)),]

  # Start loop of columns
  for (kk in kk:length(cols1)){
  
  # Cut to columns 0 is not written 1 is written
  traits<-res6[,cols1[kk]]
  states<-data.frame(name=res6$Glottocode,traits1=traits,stringsAsFactors = FALSE)
  rownames(states)<-states$name
  
  # What is missing
  states<-states[!is.na(states$traits1),]
  
  # Set missing
  fff=rep(nrow(states)/Ntip(tree),2)
  
  # Remove tips that are missing
  phy=ape::drop.tip(tree,tree$tip.label[!tree$tip.label %in% states$name])
  
  # Correct order
  states<-states[order(match(states$name,phy$tip.label)),]
  
  # Calculate DR
  eds<-caper::ed.calc(phy)$spp
  eds$DR<-1/eds$ED
  
  # Merge
  overlapping=NULL
  states2<-merge(states,eds,by.x="name",by.y="species")
  zz<-aggregate(states2$DR,by=list(states2$traits1),mean)
  names(zz)<-c("trait_value","mean_DR")
  zz2<-aggregate(states2$DR,by=list(states2$traits1),sd)
  zz3<-aggregate(states2$DR,by=list(states2$traits1),length)
  sem<-zz2$x/sqrt(zz3$x)
  zz$lower=zz$mean_DR-sem
  zz$upper=zz$mean_DR+sem
 
  # Different between 1 an 0 in mean DR
  diff1<-zz[zz$trait_value==max(zz$trait_value),"mean_DR"]-zz[zz$trait_value==min(zz$trait_value),"mean_DR"]
  
  # Test for overlap
  if(zz$lower[2]<zz$upper[1] & zz$lower[1]<zz$upper[2]){overlapping=1}else{overlapping=0}
  
  DR_test<-data.frame(tree=tt2[ii],variable=cols1[kk],name="DR_test",AIC=diff1,Parameters=overlapping,stringsAsFactors = FALSE)

  # Setting up a BiSSE model using HiSSE - (all from https://cran.r-project.org/web/packages/hisse/vignettes/hisse-new-vignette.pdf)
  
  #As with the original HiSSE implementation (see hisse.old()), the number of free parameters in the model
  #for both turnover and extinction fraction are specified as index vectors provided to the function call. Each
  #vector contains four entries that correspond to rates associated with the observed states (0 or 1) and the
  #hidden states (A or B). They are always ordered as follows for a given hidden state, i: 0i, 1i. However, in this
  #case we do not want any hidden states. But first we set up the “dull null” – i.e., turnover and extinction
  #fraction are the same for both states. Note the “f” represents the sampling fraction for each observed state
  #1 combination, which is a vector ordered in the same manner as for turnover and extinction fraction vectors:
    
  turnover <- c(1,1)
  extinction.fraction <- c(1,1)
  f <- c(1,1)
  
  # Next, we have to set up a transition matrix. There is a function provided to make this easy, and allows users
  #to customize the matrix to fit particular hypotheses. Be sure to look at the options on this function call, for
  #allowing diagonals and for customizing the matrix when you have character independent model.
  
  trans.rates.bisse <- hisse::TransMatMakerHiSSE(hidden.traits=0)
  
  # Now, we can call HiSSE and estimate the parameters under this model using the default settings:
    dull.null <- hisse::hisse(phy=phy, data=states, f=f, turnover=turnover,
                      eps=extinction.fraction, hidden.states=FALSE,
                       trans.rate=trans.rates.bisse,sann=FALSE)
  
    dn1<-data.frame(tree=tt2[ii],variable=cols1[kk],name="dull_null",dull.null$AIC,dull.null$solution,stringsAsFactors = FALSE)
    names(dn1)[4:5]<-c("AIC","Parameters")  
    
  # If you wanted to set up a true BiSSE model, where the turnover rate parameters are unlinked across the
  #observed state combinations, you would simply do the following:

  turnover <- c(1,2)
  extinction.fraction <- c(1,1)
  BiSSE <- hisse::hisse(phy=phy, data=states, f=f, turnover=turnover,
                 eps=extinction.fraction,hidden.states=FALSE,
                 trans.rate=trans.rates.bisse,sann=FALSE)
  
  
  Bis1<-data.frame(tree=tt2[ii],variable=cols1[kk],name="BiSSE",BiSSE$AIC,BiSSE$solution,stringsAsFactors = FALSE)
  names(Bis1)[4:5]<-c("AIC","Parameters")
  
  # Make overall results table
  res1<-rbind(DR_test,dn1,Bis1)
  
  # Add to last loop
  if(is.null(res2)){res2<-res1} else {res2<-rbind(res2,res1)}
  
  } # End of kk
  
  
  
  write.csv(res2, file = here("BISSE_runs", paste0("es4_", tt2[[ii]], "_", regs1, "_BISSE_runs3.csv")), row.names = FALSE)  
   
  print(ii)
} #ii
