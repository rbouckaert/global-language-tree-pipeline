
library(moments)
library(apTreeshape)
library(ape)
library (caper)
library(RPANDA)
library(phangorn)
library(phytools)
library(doParallel)
library(picante)

##read in all the trees
####need all family trees
tt2<-list.files("./all_trees",full.names=TRUE)
tt3<-list.files("./all_trees",full.names=F)
tt3<-gsub(".tree","",tt3)

### create result directory
dir.create("./actual_trees/")

###read in family and macroareas
##read in data
res5a<-read.csv(file="./all_tips_by_year6ALL.csv")


### parallel loops
cl <- makeCluster(7)
registerDoParallel(cl)

###run loop
foreach(ii=1:length(tt2),.packages=c("base","utils","ape","picante","caper","apTreeshape","moments")) %dopar% {

  ###read in one tree per loop
  tree<-ape::read.tree(tt2[ii])
  
  ###group by
  groupbytree<-tt3[ii]
  
  #make data.frame
  res1<-data.frame(tree=ii)
  res5b<-data.frame(tree=ii,tip=tree$tip.label)
  
  ##find edge lengths
  res5b$PE<-tree$edge.length[tree$edge[,2]<= Ntip(tree)]
  
  ##mean PE
  res1$meanPE<-mean(log(res5b$PE))
  
  ##sd PE
  res1$sdPE<-sd(log(res5b$PE))
  
  ##median PE
  res1$medianPE<-median(res5b$PE)
  
  ##min max PE
  res1$minPE<-min(res5b$PE)*1000
  res1$maxPE<-max(res5b$PE)*1000
  
  ##number tips under 1000 years old
  res1$numPE1000<-length(res5b$PE[res5b$PE<1])/Ntip(tree)

    #make node ages
    nodesh<-picante::node.age(tree)$ages
    nodesh2<-max(nodesh)-nodesh
     
    ##number of splits under 1000,2000,3000,4000,5000
    res1$under500<-length(nodesh2[nodesh2<0.5])/length(nodesh2)
    res1$under1000<-length(nodesh2[nodesh2<1])/length(nodesh2)
    res1$under2000<-length(nodesh2[nodesh2<2])/length(nodesh2)
    res1$under3000<-length(nodesh2[nodesh2<3])/length(nodesh2)
    res1$under4000<-length(nodesh2[nodesh2<4])/length(nodesh2)
    res1$under5000<-length(nodesh2[nodesh2<5])/length(nodesh2)
    
    ###tree length
    res1$tree_length<-sum(tree$edge.length,na.rm=TRUE)
  
    ##crown age
    res1$crown_age<-max(nodesh)
  
    ##crown age in years
    res1$crown_age2<-res1$crown_age*-1000
  
    ##tip
    res1$richness<-length(tree$tip.label)
  
    ##gamma
    res1$gammma<-ape::gammaStat(tree)
    
    ##balance
    res1$balance<-apTreeshape::colless(as.treeshape(tree))
    
    ##tip level ED
    res5b$ED<-picante::evol.distinct(tree,type="equal.splits")$w
    
    ##tip level DR
    res5b$DR<-1/res5b$ED
    
    ##mean DR
    res1$minDR<-min(res5b$DR)
    res1$maxDR<-max(res5b$DR)
    
    ##Harmonic mean DR
    res1$HmeanDR<-(1/mean(1/res5b$DR))
    
    ##mean DR
    res1$meanDR<-mean(log(res5b$DR))
    
    ##mean DR
    res1$sdDR<-sd(log(res5b$DR))
    
    ##mean DR
    res1$medianDR<-median(res5b$DR)
    
    ##DR skew
    res1$DRskew<-moments::skewness(res5b$DR)
    
    ##DR kurtosis
    res1$DRkurtosis<-moments::kurtosis(res5b$DR)
  
    ## get details  
    temp1<-read.table(text=tt3[ii],sep="_")

    ### what group and tree
    res1$group<-temp1$V2
    res1$tree<-temp1$V1
    
    ###write results
    write.csv(res1,file=paste("./actual_trees/summary_",groupbytree,".csv",sep=""))
  
    
   
} ##end ii trees

stopCluster(cl)

###create a single summary file

fil1<-list.files("./actual_trees/",full.names=TRUE)
fil2<-list.files("./actual_trees/",full.names=F)

fil2<-gsub(".csv","",fil2)
fil2<-gsub("summary_tree_","",fil2)
fil2<-gsub("lexicon_families","lexiconfamilies",fil2)
fil2<-gsub("0_","0;",fil2)
fil2<-gsub("1_","1;",fil2)
fil2<-gsub("2_","2;",fil2)
fil2<-gsub("3_","3;",fil2)
fil2<-gsub("4_","4;",fil2)
fil2<-gsub("5_","5;",fil2)
fil2<-gsub("6_","6;",fil2)
fil2<-gsub("7_","7;",fil2)
fil2<-gsub("8_","8;",fil2)
fil2<-gsub("9_","9;",fil2)

fil2b<-read.table(text=fil2,sep=";",stringsAsFactors = F)
names(fil2b)<-c("tree","group")

#fil1<-fil1b[!fil1b %in% fil1]

for (j in 1:length(fil1)){
  
  fil2x<-read.csv(fil1[j],stringsAsFactors = FALSE)
  fil2x$tree<-fil2b[j,"tree"]
  fil2x$group<-fil2b[j,"group"]
  
    
  if(j==1){ fil3<-fil2x } else{fil3<-rbind(fil3,fil2x)}
  
  print(j)
}

###rename group
fil3$X<-fil3$group

###do HDPI
fil3<-fil3[,1:(ncol(fil3)-1)]

###summarise by group
##mean
tt1<-aggregate(fil3[,3:(ncol(fil3))],by=list(fil3$X),function (x) HPDI(x,0.89))
##67
tt2<-aggregate(fil3[,3:(ncol(fil3))],by=list(fil3$X),function (x) HPDI(x,0.67))
##97
tt3<-aggregate(fil3[,3:(ncol(fil3))],by=list(fil3$X),function (x) HPDI(x,0.97))

##median
tt4<-aggregate(fil3[,3:(ncol(fil3))],by=list(fil3$X),median,na.rm=TRUE)
##median
tt5<-aggregate(fil3[,3:(ncol(fil3))],by=list(fil3$X),min,na.rm=TRUE)
##median
tt6<-aggregate(fil3[,3:(ncol(fil3))],by=list(fil3$X),max,na.rm=TRUE)


nd2<-NULL
###reorder data frames
for (f in 2:ncol(tt4)){
  
  nd<-data.frame(type=tt1$Group.1,variable=names(tt1)[(f)],range89=tt1[(f)],range67=tt2[(f)],range97=tt3[(f)],median=tt4[,f],min=tt5[,f],max=tt6[,f])
  names(nd)[3:5]<-c("range89","range67","range97")
  
  if(is.null(nd2)){nd2<-nd} else{nd2<-rbind(nd2,nd)}
  
}

##add on group
nd3<-nd2#merge(nd2,all1,by.x="language_group",by.y="group")


#write results
write.csv(nd3,file="./summary_statsALL_global_HPDI2.csv")

