

####Make and test EDGE scores
####threat from EGIDS and glottlog

library(ape)
library(caper)
library(rgdal)
library(sp)
library(data.table)
library(raster)
library(mice)
library(phylobase)

### create result directory
dir.create("./EDGE_scores/")


##EDGE2 code written by Ricky Gummbs
source("./EDGE2.R")

##read in all the trees
trn<-list.files("./all_trees",full.names=TRUE)


##load prepared threat data
load(file="./res5.r")

###subset needed columns
es3b<-res5[,c("glottocode","threat_glot","threat_EGIDS","L1_Users","Shap_Ar","X10")]
es3b$L1_Users<-log(es3b$L1_Users+1)
es3b$Shap_Ar<-log(es3b$Shap_Ar+1)

###impute missing values
temp1<-suppressWarnings(mice(es3b,m=5,maxit=50,method="cart",print=FALSE))
es3_1<-complete(temp1,1)
es3_2<-complete(temp1,2)
es3_3<-complete(temp1,3)
es3_4<-complete(temp1,4)
es3_5<-complete(temp1,5)

###make mean values
es3<-data.frame(glottocode=es3_1$glottocode,(es3_1[,2:6]+es3_2[,2:6]+es3_3[,2:6]+es3_4[,2:6]+es3_5[,2:6])/5)
es3$threat_EGIDS<-round(es3$threat_EGIDS,0)

##get rid of 10 - merge with 9
es3$threat_EGIDS[es3$threat_EGIDS==10]<-9

##ext_prob
names(es3)[6]<-"ext_prob"

es4<-as.data.frame(es3)#[agg2,nomatch=0]


minx<-min(es4$ext_prob[es4$ext_prob!=0])
##set any too small prob to 0
es4$ext_prob[es4$ext_prob<=0]=minx
#es4[es4$ext_prob>1,]$ext_prob=1

##rid of any extra columns
es4<-es4[,c("glottocode","ext_prob")]
names(es4) <- c("species","pext")
es4$species[es4$species=="osse1243"]<-"iron1242"

##resample order
trn<-sample(trn)

for (i in 1:(length(trn))){

filen<-gsub("all_trees/","EDGE_scores/",trn[i],fixed=TRUE)
filen<-gsub(".tree","_EDGESCORES2.csv",filen,fixed=TRUE)
filen2<-filen
filen2<-gsub(".csv","_PDLOST2.csv",filen,fixed=TRUE)

###check if done
if(filen %in% list.files("./EDGE_scores/",full.names=TRUE)){next}

##read tree
tree<-read.tree(trn[i])

##match to data
es5<-es4[match(tree$tip.label,es4$species),]

### run edge2 calculation
EDGE1<-EDGE.2.calc(tree=tree,pext=es5)

##
e1<-EDGE1[[1]]

###PD
PD0<-sum(tree$edge.length)

samp1<-na.omit(es3$glottocode[es3$threat_EGIDS>6])

##make tree with threat removed
tree2<-drop.tip(tree,tree$tip.label[tree$tip.label %in% samp1 ])

##difference in PD loosing threat
PD1<-PD0-sum(tree2$edge.length)

###drop random
tree3<-drop.tip(tree,tree$tip.label[tree$tip.label %in% sample(tree$tip.label,length(samp1))])

##difference from random
PD2<-PD0-sum(tree3$edge.length)

##make results data frame
resx1<-data.frame(tree=i,total=PD0,total_lost_threat=PD2,total_lost_random=PD2)

##write files
write.csv(resx1,file=filen2)
write.csv(EDGE1[[1]],file=filen)

print(i)
}


 