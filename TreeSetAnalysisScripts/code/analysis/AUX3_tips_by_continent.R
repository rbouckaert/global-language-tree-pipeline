
##Setting up a BiSSE model using HiSSE - (all from https://cran.r-project.org/web/packages/hisse/vignettes/hisse-new-vignette.pdf)

library(caper)
library(dismo)
library(phangorn)
library(phytools)
library(raster)
library(paleotree)
library(diversitree)
library(hisse)
library(ape)
library(data.table)
library(rgdal)
library(sp)

### setwd
setwd("C:\\Users\\xxxx\\Documents\\")

tt3a<-cbind(tt3,tt)

tt3b<-tt3a[tt3a$V2==100,]

for (k in 1:6){

tt3d<-tt3b[k,]
  
tree1<-ape::read.tree(tt3d$tt)

res1<-data.frame(tip.labels=tree1$tip.label,group=tt3d$V3)

if(k==1) {res2<-res1} else {res2<-rbind(res2,res1)}

}

table(res2$group)

res2$group[res2$group=="North-America"]<-"America"
res2$group[res2$group=="South-America"]<-"America"

write.csv(res2,file=".\\input_data\\tips_by_continent.csv")

