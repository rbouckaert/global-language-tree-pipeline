

library (caper)
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

e <- simpleError("test error")

### create result directory
dir.create("./BISSE_FISSE/")

###tip data
res5<-fread(file="./all_tips_by_year6ALL.csv")
res5$regions<-res5$glottolog_macroarea
res5$regions[res5$regions=="South America"]<-"Americas"
res5$regions[res5$regions=="North America"]<-"Americas"
res5$regions[res5$regions=="Papunesia"]<-"Oceania"
res5$regions[res5$regions=="Australia"]<-"Oceania"

##make spatial
res5sp<-res5
coordinates(res5sp)<-~longitude + latitude

### Read in trees
tt0<-list.files("./all_trees",full.names=TRUE)

##tree name
tt2a<-list.files("./all_trees",full.names=FALSE)
tt2a<-gsub(".tree","",tt2a)
tt3a<-read.table(text=tt2a,sep="_")

##combine together
tt<-c(tt0)
tt2<-c(tt2a)
tt3<-rbind(tt3a)

### ONLY NEED GLOBAL TREES
tt<-tt[tt3$V3=="Global"]
tt2<-tt2[tt3$V3=="Global"]
tt3<-tt3[tt3$V3=="Global",]

##reorder sample
sam1<-sample(1:(length(tt)),replace=FALSE)

##combine together and reorder
tt<-tt[sam1]
tt3<-tt3[sam1,]
tt2<-tt2[sam1]

#trees<-NULL
#res4<-NULL
ii=1

###start loop
for( ii in ii:length(tt)){
  res2<-NULL
  kk=1
  
  regs1<-"ALL"
  
  if(paste("es4_",tt2[[ii]],"_",regs1,"_BISSE_FISSE2.csv",sep="") %in% list.files("C:/Users/xxxx/Documents/BISSE_FISSE/",full.names=FALSE)) {next}
  
  
    #for simulated
  tree<-ape::read.tree(tt[ii])
  
  ##correct ossetic error
  res5$Glottocode [ res5$Glottocode =="osse1243"]<-"iron1242" 
  
  ##combine forager
  res5$Forager<-res5$`Forager language definite`
  res5$Forager[res5$`Forager language possible`==1]<-1
  
  ###Columns to run
  cols1<-c("island","written","EA033_Extended")#names(res5)[c(9,11,106,107)]

	#make new data frame
  res6<-as.data.frame(res5)
  res6<-res6[res6$regions!=regs1,]
  rownames(res6)<-res6$Glottocode
  
    ##Run loop
  for (kk in kk:length(cols1)){
  
  ##cut to columns 0 is not written 1 is written
  traits<-res6[,cols1[kk]]
  states<-data.frame(name=res6$Glottocode,traits1=traits,stringsAsFactors = FALSE)
  rownames(states)<-states$name
  
    
  ## what is missing
  states<-states[!is.na(states$traits1),]
  
  ##set missing
  fff=rep(nrow(states)/Ntip(tree),2)
  
  ##remove tips that are missing
  phy=drop.tip(tree,tree$tip.label[!tree$tip.label %in% states$name])
  
  ##correct order
  states<-states[order(match(states$name,phy$tip.label)),]
  
  ##calculate DR
  eds<-ed.calc(phy)$spp
  eds$DR<-1/eds$ED
  
  ##merge
  overlapping=NULL
  states2<-merge(states,eds,by.x="name",by.y="species")
  zz<-aggregate(states2$DR,by=list(states2$traits1),mean)
  names(zz)<-c("trait_value","mean_DR")
  zz2<-aggregate(states2$DR,by=list(states2$traits1),sd)
  zz3<-aggregate(states2$DR,by=list(states2$traits1),length)
  sem<-zz2$x/sqrt(zz3$x)
  zz$lower=zz$mean_DR-sem
  zz$upper=zz$mean_DR+sem
  
  ##different between 1 an 0 in mean DR
  diff1<-zz[zz$trait_value==max(zz$trait_value),"mean_DR"]-zz[zz$trait_value==min(zz$trait_value),"mean_DR"]
  
  #test for overlap
  if(zz$lower[2]<zz$upper[1] & zz$lower[1]<zz$upper[2]){overlapping=1}else{overlapping=0}
  
  ##run DR overlap
  DR_test<-data.frame(tree=tt2[ii],variable=cols1[kk],name="DR_test",AIC=diff1,Parameters=overlapping,stringsAsFactors = FALSE)

	###kept in for hisse
  ##set parameters
  #turnover <- c(1,1)
  #extinction.fraction <- c(1,1)
  #f <- rep(1,Ntip(phy))
  
  #Next, we have to set up a transition matrix. There is a function provided to make this easy, and allows users
  #to customize the matrix to fit particular hypotheses. Be sure to look at the options on this function call, for
  #allowing diagonals and for customizing the matrix when you have character independent model.
  
  #trans.rates.bisse <- TransMatMakerHiSSE(hidden.traits=0)
  #print(trans.rates.bisse)
  
  #Now, we can call HiSSE and estimate the parameters under this model using the default settings:
  #dull.null <- hisse(phy=phy, data=states, f=f, turnover=turnover,
   #                  eps=extinction.fraction, hidden.states=FALSE,
   #                  trans.rate=trans.rates.bisse)
  
  #dull.null <- hisse(phy=phy, data=states, hidden.states=FALSE, turnover=turnover,eps=extinction.fraction, trans.rate=trans.rates.bisse)
                    
  #If you wanted to set up a true BiSSE model, where the turnover rate parameters are unlinked across the
  
  
  #observed state combinations, you would simply do the following:
  #turnover <- c(1,2)
  #extinction.fraction <- c(1,1)
  #BiSSE <- hisse(phy=phy, data=states, turnover=turnover,
  #               eps=extinction.fraction, hidden.states=FALSE,
  #               trans.rate=trans.rates.bisse,sann=FALSE)
  
  #Bis1<-data.frame(tree=tt2[ii],variable=cols1[kk],name="BiSSE",BiSSE$AIC,BiSSE$solution,stringsAsFactors = FALSE)
  #names(Bis1)[4:5]<-c("AIC","Parameters")
  
  #Setting up a HiSSE model
  #Setting up a character-dependent HiSSE model is relatively straightforward, and relies on all the same tools
  #as above. One important thing to bear in mind, is that again, the states are ordered by state combination
  #within each hidden state. For example, if you want two hidden states, A and B, the order of the parameters
  #in the model is 0A, 1A, 0B, and 1B. So, in this case, we just need to specify the free parameters we want for
  #diversification. Here we are just going assume turnover varies across the different states:
  turnover <- c(1,2,3,4)
  extinction.fraction <- rep(1, 4)
  #f = c(1,1)
  
  #We also have to extend the transition rate matrix. This is done by specifying the number of hidden states:
  trans.rates.bisse <- TransMatMaker(hidden.states=TRUE)  #print(trans.rate.hisse)
  
  #Now, we can just plug these options into MuHiSSE:
  HiSSE <-  hisse(phy=phy, data=states, f=fff, turnover=turnover,
                  eps=extinction.fraction, hidden.states=TRUE,
                  trans.rate=trans.rates.bisse,sann=FALSE)
  
  #??"extinction.fraction"
  #extinctionfraction= extinction/speciation
    
  #turnoverrate=speciation+extinction
  
  HiSSE$solution[1,1]

  His1<-data.frame(tree=tt2[ii],variable=cols1[kk],name="HiSSE",HiSSE$AIC,HiSSE$solution,stringsAsFactors = FALSE)
  names(His1)[4:5]<-c("AIC","Parameters")
  
  #Setting up a character-independent HiSSE model (i.e., CID-2, CID-4)
  #Remember, for any character-independent model, the diversification rates must be decoupled from the
  #observed states. So, to do this, we simply set the diversification rates to be equal for all states for a given
  #hidden state. Below, I will show how to do this for a character-independent model with two rate shifts in the
  #tree, what we refer to as the CID-2 model:
  turnover <- c(1, 1, 2, 2)
  extinction.fraction <- rep(1, 4)
  f = c(1,1)
  trans.rate.cid2 <- TransMatMakerHiSSE(hidden.traits=1, make.null=TRUE)
  
  
  HiSSE.null.cid2 <- hisse(phy=phy, data=states, turnover=turnover,
                          eps=extinction.fraction, hidden.states=TRUE,
                          trans.rate=trans.rate.cid2,sann=FALSE)
  
  Hcid2<-data.frame(tree=tt2[ii],variable=cols1[kk],name="HiSSE.null.cid2",HiSSE.null.cid2$AIC,HiSSE.null.cid2$solution,stringsAsFactors = FALSE)
  names(Hcid2)[4:5]<-c("AIC","Parameters")
  
  #For the transition rate matrix, I included a setting called make.null that simply replicates the transition
  #model across the hidden states. This way the transition rates are not impacted by changes in the diversification
  #rate regime - that is, q0,A!1,A = q1,B!0,B.
  #Here is a how you would set up the diversification rates for models with three hidden state, which is equivalent
  #to CID-4 model:
  
  #turnover <- c(1, 1, 2, 2, 3, 3, 4, 4)
  #extinction.fraction <- rep(1, 8)
  #trans.rate.cid4 <- TransMatMakerHiSSE(hidden.traits=3, make.null=TRUE)
  
  
  #HiSSE.null.cid4 <- hisse(phy=phy, data=states, f=f,turnover=turnover,
   #                        eps=extinction.fraction, hidden.states=TRUE,
   #                        trans.rate=trans.rate.cid4,sann=FALSE)
  
  #Hcid4<-data.frame(tree=tt2[ii],variable=cols1[kk],name="HiSSE.null.cid4",HiSSE.null.cid4$AIC,HiSSE.null.cid4$solution,stringsAsFactors = FALSE)
  #names(Hcid4)[4:5]<-c("AIC","Parameters")
  
  ##make results
#  res1<-rbind(DR_test,Bis1,His1,Hcid2,Hcid4)#,Fisse_lambd0=resF$lambda0,Fisse_lambd1=resF$lambda1,Fisse_pval=resF$pval,Fisse_pval_2tailed=min(resF$pval, 1-resF$pval)*2)
  res1<-rbind(DR_test,His1,Hcid2)#,Fisse_lambd0=resF$lambda0,Fisse_lambd1=resF$lambda1,Fisse_pval=resF$pval,Fisse_pval_2tailed=min(resF$pval, 1-resF$pval)*2)
  
  if(is.null(res2)){res2<-res1} else {res2<-rbind(res2,res1)}
  
  }##end of kk
  
   write.csv(res2,file=paste("./BISSE_FISSE/es4_",tt2[[ii]],"_",regs1,"_BISSE_FISSE2.csv",sep=""),row.names = FALSE)  
   
  print(ii)
} ##ii
    