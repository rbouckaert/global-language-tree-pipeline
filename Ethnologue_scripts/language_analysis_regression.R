library(ape)
library(nlme)
library(ade4)
library(dismo)
library(picante)
library(INLA)
library(spdep)

##set error term
e <- simpleError("test error")

### create result directory
dir.create("./regression_results/")
dir.create("./cl_graph/")

###files to process 
st0<-list.files("./datasets_and_trees/",pattern="_4.csv",full.names=TRUE)

##trees data
st2a<-list.files("./datasets_and_trees/",pattern="_4.csv",full.names=FALSE)
st2b<-gsub(".csv","",st2a)
st3a<-read.table(text=st2b,sep="_",stringsAsFactors = FALSE)
files1<-data.frame(st3a,filen=st0,uni=st2b,stringsAsFactors = FALSE)

### Read in trees
tt0<-list.files("./datasets_and_trees/",pattern="_4.tre",full.names=TRUE)
tt0<-tt0[!tt0 %in% st0]

##tree data
tt2a<-list.files("./datasets_and_trees/",pattern="_4.tre",full.names=FALSE)
tt2a<-tt2a[!tt2a %in% st2a]
tt2a<-gsub(".tre","",tt2a,fixed=TRUE)
tt3a<-read.table(text=tt2a,sep="_",stringsAsFactors = FALSE)
tt<-data.frame(tt3a,filen=tt0,uni=tt2a,stringsAsFactors = FALSE)

##reorder sample
sam1<-sample(1:nrow(files1),nrow(files1),replace=FALSE)

##combine together and reorder
files1<-files1[sam1,]

files1[,"V5"]<-as.numeric(gsub("XXX","",files1[,"V5"]))

##start loop
for (ii in 1:nrow(files1)){
  
  ###error distribution
  typeX="poisson"
  
  ###check if done
  if(paste("results_",files1[ii,"uni"],"_",typeX,"_regressions.csv",sep="") %in% list.files("/regression_results/",full.names=FALSE)){next}  
  
  #do simulations
  if(files1[ii,"V2"]=="simulatione"){next}
  
  res4<-NULL
  ##find file names
  files2<-files1[ii,]
  tree1<-tt[tt$uni==files2$uni,]
  
  ##read in files
  es4<-read.csv(file=as.character(files2$filen),stringsAsFactors = FALSE)
  if(files1[ii,"V2"]=="simulatione"){
      es4$order=sample(1:nrow(es4))
      es4$countup=1
      for (uu in 2:nrow(es4)) {
          if(es4$region[uu]==es4$region[uu-1]){es4$countup[uu]=es4$countup[uu-1]           }else{es4$countup[uu]=es4$countup[uu-1]+1
        }
      }
      
    es4<-es4[order(es4$countup,es4$order),]
  }  
  
  row.names(es4)<-as.character(es4$tips)
  tr1<-read.tree(file=as.character(tree1$filen))
  
  ###ages
  node_ages<-picante::node.age(tr1)$ages
  tot_time<-max(node_ages)

 
  ##geographical LOOCV
  large_regions<-as.data.frame(table(es4$region))
   large_regions<-as.numeric(levels(large_regions[large_regions$Freq>=30,"Var1"]))
  
	##forage quadratic term
    es4$Forage1_2<-es4$Forage1^2
    
	###prediction vector
    es4$pred<-NA
    
    
    ##geog hold out
    kf<-large_regions;loos_type="geog"#}
    
    ##out of sample predcitions
    for (ww in sort(as.numeric(unique(kf)))){
      estr<-es4

		##make response variable
      estr$y<-(estr$length);type="length"#}
      estr$y2<-estr$y
      estr$y[kf==ww]<-NA

      ##make unique names
      estr$names<-1:nrow(estr)
	  
      ##make spatial
      esqp<-estr
      coordinates(esqp)<-~longitude+latitude
      rownames(esqp@data)<-esqp@data$names
      
      ##k nearest neightbours
      nbs<-knearneigh(coordinates(esqp), k = 10, longlat = T) #k=5 nearest neighbors
      nbs<-knn2nb(nbs,  sym = T) #force symmetry
      mat <- nb2mat(nbs, style="B",zero.policy=TRUE)
      rownames(mat)<-esqp$names
      colnames(mat) <- rownames(mat) 
      mat <- as.matrix(mat[1:dim(mat)[1], 1:dim(mat)[1]])
      samp1<-sample(1e10,1)
      nb2INLA(paste0(".//cl_graph//cl_graph_",samp1),nbs)
      H<-inla.read.graph(filename=paste0(".//cl_graph//cl_graph_",samp1))
      
      ##make tree vcv matrix
      phylo_covar_mat <- ape::vcv(tr1)
      phylo_covar_mat <- phylo_covar_mat/max(phylo_covar_mat)
      phylo_covar_mat <- phylo_covar_mat / exp(determinant(phylo_covar_mat)$modulus[1] /
                                                 nrow(phylo_covar_mat))
      phylo_prec_mat <- solve(phylo_covar_mat)
      
      pcprior <- list(prec = list(prior="pc.prec", param = c(0.36, 0.1)))
      
      estr2<-cbind(estr,data.frame(species = rownames(phylo_prec_mat),
                         phy_id = 1:nrow(phylo_prec_mat),stringsAsFactors = FALSE))
      
      
		##make subregion		
      estr2$subregion<-estr$names
     
	 ##type of model
      model1<-"include_REGION"
	  
	  ##set model formula
      formx<-formula(paste0("y~  island + I(island^2) + area + I(area^2) + distancetocityyear2 + I(distancetocityyear2^2) + popdL + I(popdL^2)  + friction + I(friction^2)  +f(clade_age, model='linear', constr = TRUE, scale.model = TRUE)","+  f(subregion, model = 'bym2', graph = H,constr=TRUE,scale.model=TRUE) + f(phy_id, model = 'generic0', Cmatrix = phylo_prec_mat,constr = TRUE, hyper = pcprior)"))
      
         
      ## Joining, by = "species"
      lm.brown <- tryCatch(inla(formx,
                    data = estr2,
                    control.fixed = list(
                    mean= 0,  #mean for betas
                    prec= 1/(0.3^2), #precision for intercept: sd = 0.6 --> precision =1/variance --> 1/(0.6^2) = 2.777778
                    mean.intercept= 0, #mean for intercept
                    prec.intercept= 1/(0.6^2)),
                    family = typeX,
                    control.family = list(list(link = 'default')),
                    control.compute=list(config=TRUE, dic=FALSE, waic= TRUE, cpo=TRUE),
                    control.predictor = list(compute = TRUE, link = 1)),
                    error=function (e) e)
  
	###if error
		if(class(lm.brown)[1]=="simpleError"){next}
   
		##get predictions
      pred1<-lm.brown$summary.linear.predictor[1:nrow(estr),]
      es4$pred[kf==ww]<-pred1[kf==ww,1]
      es4$resid<-estr$y2-exp(pred1$mean)
     
	 ##mantel test
	 man1<-data.frame(obs=NA,pvalue=NA)
      
      ##get psuedo r2 and rmse only at max kf
      if(ww==max(kf)){
        psudor2<-cor(exp(es4$pred),es4$y,method="pearson",use="complete.obs")
        rmse1<-sqrt(mean((es4$y-es4$pred)^2,na.rm=TRUE))}else{
        psudor2<-NA;rmse1<-NA
        
        actual.dists <- dist(cbind(es4$longitude, es4$latitude))
      resid.dists <- dist(abs(es4$resid))
      
      man1<-ade4::mantel.rtest(actual.dists, resid.dists, nrepet = 999)
      
        }
      
      ##create summary
      res1<-lm.brown$summary.fixed
      res1$name<-row.names(res1)
      
      res1$sig1<-0
      res1$sig1[(res1$`0.025quant`<0 & res1$`0.975quant`<0)|(res1$`0.025quant`>0 & res1$`0.975quant`>0)]<-1
      res1$mean2<-res1$mean*res1$sig1

      res1$n=NA#length(models[[qqq]])
      res1$model=model1#ref1[qqq]
      res1$model2=nrow(estr2)#ref2[qqq]
      res1$kfold=ww
      res1$lastkfold=max(kf)
      res1$AIC=lm.brown$waic$waic
      res1$Time=tot_time-(files1[ii,"V5"]/1000)
      res1$Crown=tot_time
      res1$Time2=files1[ii,"V5"]*-1
      res1$loop=ii
      res1$qqqloop=NA
      res1$Form=as.character(formula(formx))[3]
      res1$type=typeX
      res1$loos_type=loos_type
      res1$model_name=NA#mod_names[[ref1[qqq]]]
      res1$shapiro<-NA#round(shapiro.test(lm.brown$residuals[,1])$p.value,3)
      res1$mantcor<-man1$obs
      res1$mantel<-man1$pvalue
      res1$psudor2<-psudor2
      res1$rmse<-rmse1
      res1$tree_type=files1[ii,"V2"]
      res1$tree=files1[ii,"V3"]
      res1$group=files1[ii,"V4"]
      
      
      if(is.null(res4)){res4<-res1} else {res4<-rbind(res4,res1)}
      
      rm(res1,lm.brown)
      
      #
    }##end of wx forms loop
    
 
  
  print(ii)
  
  
  write.csv(res4,file=paste("./regression_results/results_",files1[ii,"uni"],"_",typeX,"_regressions.csv",sep=""),row.names = FALSE)  
  
  
}##end of ii loop

