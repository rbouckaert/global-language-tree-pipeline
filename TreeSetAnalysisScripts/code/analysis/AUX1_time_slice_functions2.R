
###functions
### get data per clade from InfTemp
get_dataenv<-function (x,env_data,age,tot_time){
  
  maxt<-(tot_time-age)*-1000#max(node.age(x)$ages)*-1000
  
  Time3<-env_data[,base::get(base::names(env_data)[base::names(env_data)=="Time3"])]
  
  nearest_time<-Time3[((Time3-maxt)^2)==min((Time3-maxt)^2)][1]
  
  glottocode<-env_data[,base::get(base::names(env_data)[base::names(env_data)=="glottocode"])]
  
  env1b<-env_data[(glottocode %in% x[["tip.label"]]) & Time3>nearest_time, ]
  
  ###
  if(nrow(env1b)==0){env3<-data.frame(env1b[1,base::names(env1b)[!base::names(env1b) %in% c("glottocode","LANG_IS")],with=FALSE]);env3[1,]<-NA;return(env3)}
  
  env2<-stats::aggregate(env1b[,base::names(env1b)[!base::names(env1b) %in% c("glottocode","LANG_IS")],with=FALSE],sum,na.rm=TRUE,by=list(env1b[,base::get(base::names(env1b)[base::names(env1b)=="glottocode"])]))
  
  return(data.frame(stats::aggregate(env2[,base::names(env2)[base::names(env2)!="Group.1"]] ,mean,na.rm=TRUE,by=list(rep(1,nrow(env2)))))[,-1])
  
  
  #subset(env2,select= -c(Group.1))
}

###regions
subs1x<-function(x,res5){
  res6<-res5[res5$glottocode %in% x$tip.label, ]
  return(stats::aggregate(res6[, c("REGION","SUBREGION")],by=list(rep(1,nrow(res6))),raster::modal,na.rm=TRUE))
}



## mean diversification rate
hm_es<-function (x) (1/mean(picante::evol.distinct(x,type="equal.splits")$w))
#es2$es_hm<-do.call(c,lapply(tr2,hm_es))
#hist(es2$es_hm)

## mean diversification rate
lm_es<-function (x) (mean(log(1/picante::evol.distinct(x,type="equal.splits")$w)))
#es2$es_lm<-do.call(c,lapply(tr2,lm_es))
#hist(es2$es_lm,main=paste("norm=",round(shapiro.test(es2$es_lm)$p.value,3)," n=",nrow(es2)," age=",round(j,3),sep=""))
#print(j)
#}
#mean vane_write? for deeper splits
#crown rate ? tree shape

#gamma 
#es2$gamma<-do.call(c,lapply(tr2,gammaStat))

##clade size
len1<-function (x) Ntip(x)
#es2$length<-do.call(c,lapply(tr2,len1))

##add length
clade_age<-function (x) max(node.age(x)$ages)*-1000
#es2$clade_age<-do.call(c,lapply(tr2, clade_age))

##skew
sq_es<-function (x) skewness(evol.distinct(x,type="equal.splits")$w)
#es2$es_sq<-do.call(c,lapply(tr2,sq_es))
#es2$es_sq[is.na(es2$es_sq)]<-0

##forager languages function
subs1<-function(x,res5){
  res6<-res5[res5$glottocode %in% x$tip.label, ]
  return(stats::aggregate(res6[, c("DPLACE subsistence","Forager_language_definite","Forager_language_possible","Shap_Ar","L1_Users","island","written","latitude","longitude" )],by=list(rep(1,nrow(res6))),mean,na.rm=TRUE))
}


##scale function
scale2<-function(x){
  
  if(!is.numeric(x)){return(x)}
  
  if(length(na.omit(x))==0){return(rep(0,length(x)))}
  
  if(var(na.omit(x))==0){return(x)}
  
  return(scale(x))
}

