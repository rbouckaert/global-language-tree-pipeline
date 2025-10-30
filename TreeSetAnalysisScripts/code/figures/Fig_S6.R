# Code to create results and plots for BISSE analyses/Extended Data Fig. 6

library(ggplot2)
library(cowplot)
library(ggpubr)
library(gridExtra)
library(RColorBrewer)
library(ggridges)
library(data.table)
library(rethinking)
library(forcats)
library(here)

# Get bisse results data
lf3 = read.csv(here("outputs", "bisse_results1.csv"))

# Load parameter names object called xxx
load(here("input_data", "bisse_param_names.r"))

# What types of test
types<-unique(lf3$name)
print(types)

# Subset and rename DR
DR1 <- lf3[lf3$name=="DR_test",]
DR1$variable <- gsub("Marriage_Endogamous", "Endogamous Marriage", DR1$variable)
DR1$variable <- gsub("Marriage_Exogamous", "Exogamous Marriage", DR1$variable)
DR1$variable <- gsub("D_Place_", "", DR1$variable)
DR1$variable <- gsub(" definite", "", DR1$variable)
DR1$variable <- gsub("island", "Island", DR1$variable)
DR1$variable <- gsub("_Pattern", "", DR1$variable)
DR1$variable <- gsub("EA033_Extended", "Political Complexity", DR1$variable)

# Plot 1 - Results of DR test --------------------------------------------------
DR1$location <- factor(DR1$location, levels = c("Global","Africa", "America","Eurasia", "Oceania"))
cols = c("black","#EE0000FF","#008B45FF","#631879FF","#008280FF")

p1 = ggplot(DR1[DR1$Parameters==0 & DR1$variable!="Endogamous Marriage",], 
            aes(x=AIC,y=location,fill=factor(location), alpha = factor(location))) +
  geom_density_ridges(bandwidth = 0.015,scale = 3.75, color = NA)+
  facet_wrap2(.~variable,scales = "free_x")+
  geom_vline(xintercept = 0,lty=2 )+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = cols) +
  scale_alpha_manual(values = c(0.9, 0.5, 0.5, 0.5, 0.5)) +
  coord_cartesian(clip = "off") +
  theme_cowplot(12)+
  labs(fill = "Macroregion excluded", alpha = "Macroregion excluded") +
  xlab("Difference in DR")+
  theme(axis.title.y = element_blank(),
        legend.position = "top",
        strip.text = element_text(face = "bold"))

print(p1)

ggsave(here("outputs","Fig_S5", "Fig_S5_ridges.png"), p1, height=3, width=8, dpi = 600)

# Figure S6 --------------------------------------------------------------------
# Updated code for faceted ridges
p5 <- ggplot(DR1,  # Use the full dataset, not just the subset
             aes(
               x = AIC, 
               y = fct_reorder(variable, AIC), 
               fill = as.factor(Parameters)  # Map fill to Parameters
             )) +
  geom_density_ridges(bandwidth = 0.015, alpha = 0.6) +  # Use alpha for transparency
  facet_wrap(. ~ variable, scales = "free") +  # Facet by variable
  geom_vline(xintercept = 0, lty = 2) +  # Vertical line at x=0
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c("#FF5733", "#3498DB"),
                    labels = c("Non-overlapping", "Overlapping"),  # Rename "0" and "1"
                    name = "95% CIs") +  # Custom fill colors
  theme_cowplot(12) +
  xlab("Difference in DR") +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    #axis.title.x = element_blank(),
    legend.position = c(0.84, 0.82)  # Add a legend for clarity
  )

print(p5)

ggsave(here("outputs", "Fig_S6.png"), p5, height=3, width=8, dpi = 600)
ggsave(here("outputs", "Fig_S6.pdf"), p5, height=3, width=8, dpi = 600)

# Continue Data Curation --------------------------------------------------
# Results table
resDR<-aggregate(DR1[,c("AIC","Parameters")],by=list(DR1$variable),mean)
names(resDR)<-c("Variable","Difference","Overlap")



# Continue Data Curation and Plot BISSE --------------------------------------------------
# Summarise and plot BISSE
bisse1<-c("dull_null","BiSSE")

# Get 0A and 1A parameters
bisp<-lf3[lf3$name %in% bisse1 ,]
bisp$param_names<-c(xxx,NA) #c("turnover0A","turnover1A",paste("scrap_",1:22))
bisp<-bisp[ bisp$param_names %in% c("turnover0A","turnover1A","eps0A","eps1A"), ]
bisp$uni<-paste(bisp$name,bisp$param_names,sep="_")
bisp$uni2<-paste(bisp$tree,bisp$variable,bisp$location,sep=";")

# Change names for clarity
bisp$param_names[bisp$param_names=="turnover1A"]<-"State 1"
bisp$param_names[bisp$param_names=="turnover0A"]<-"State 0"

# AIC results
resAIC<-aggregate(bisp$AIC,by=list(bisp$variable,bisp$name,bisp$location),mean)
names(resAIC)[ncol(resAIC)]<-"mean_AIC"

# How often parameters higher or lower
param1<-bisp[bisp$name=="BiSSE" & bisp$param_names=="State 1",]
param0<-bisp[bisp$name=="BiSSE" & bisp$param_names=="State 0",]
identical(param1$uni2,param0$uni2)

# Get mean parameter of state 1
mean1<-aggregate(param1$Parameters,by=list(param1$variable,param1$location),mean)
sd1<-aggregate(param1$Parameters,by=list(param1$variable,param1$location),sd)

# Get mean parameter of state 0
mean0<-aggregate(param0$Parameters,by=list(param0$variable,param1$location),mean)
sd0<-aggregate(param0$Parameters,by=list(param0$variable,param1$location),sd)

# Results dataframe
resParam<-mean1[,c("Group.1","Group.2")]

# Find difference
resParam$zero_over_one<-mean0$x/mean1$x
resParam$one_over_zero<-mean1$x/mean0$x

# Where state 1 is greater than state 0
param1$diff=param1$Parameters-param0$Parameters
param1$diff[param1$diff>0]<-1
param1$diff[param1$diff<0]<-0
param1$dummy<-1

# Aggregate 1 greater than 0
x<-xtabs(param1$diff~param1$variable+param1$location)
# How many runs
y<-xtabs(param1$dummy~param1$variable+param1$location)
# Standarised version of x
z<-x/y

# Diversification Conversion and Analysis --------------------------------------------------
# Wide version
bisp_wide<-dcast(bisp,uni2~uni,value.var = "Parameters")

# Calculate diversification
bisp_wide$diversification0A=(bisp_wide$BiSSE_turnover0A-(bisp_wide$BiSSE_turnover0A*bisp_wide$BiSSE_eps0A))/(1+bisp_wide$BiSSE_eps0A)
bisp_wide$diversification1A=(bisp_wide$BiSSE_turnover1A-(bisp_wide$BiSSE_turnover1A*bisp_wide$BiSSE_eps1A))/(1+bisp_wide$BiSSE_eps1A)
bisp_wide$nulldiversification0A=(bisp_wide$dull_null_turnover0A-(bisp_wide$dull_null_turnover0A*bisp_wide$dull_null_eps0A))/(1+bisp_wide$dull_null_eps0A)
bisp_wide$nulldiversification1A=(bisp_wide$dull_null_turnover1A-(bisp_wide$dull_null_turnover1A*bisp_wide$dull_null_eps1A))/(1+bisp_wide$dull_null_eps1A)

# Define what you want
bisp_wide2<-bisp_wide[,c("uni2","diversification0A","diversification1A","nulldiversification0A","nulldiversification1A")]

# Wide to long
bisp_long3<-melt(bisp_wide2,id.vars="uni2")
bisp_long4<-cbind(bisp_long3,read.table(text=bisp_long3$uni2,sep=";",stringsAsFactors = F))   
bisp_long5<-cbind(bisp_long4,read.table(text=bisp_long4$V1,sep="_",stringsAsFactors = F))  
names(bisp_long5)<-c("uni2","param_names","Parameters","uni1","variable","area",
                     "type","tree","area1","type2")
bispcols1<-read.table(text=bisp_long5$uni2,sep=";")
bisp_long5$location<-bispcols1$V3

# Create plots
# Rename variable names
bisp_long6<-bisp_long5[!is.na(bisp_long5$Parameters),]
bisp_long6$variable<-gsub("Marriage_Endogamous","Endogamous Marriage",bisp_long6$variable)
bisp_long6$variable<-gsub("Marriage_Exogamous","Exogamous Marriage",bisp_long6$variable)
bisp_long6$variable<-gsub("D_Place_","",bisp_long6$variable)
bisp_long6$variable<-gsub(" definite","",bisp_long6$variable)
bisp_long6$variable<-gsub("island","Island",bisp_long6$variable)
bisp_long6$variable<-gsub("_Pattern","",bisp_long6$variable)
bisp_long6$variable<-gsub("EA033_Extended","Political Complexity",bisp_long6$variable)

# Rename param names
bisp_long6$param_names<-gsub("nulldiversification0A","Constant Rate",as.character(bisp_long6$param_names))
bisp_long6$param_names<-gsub("nulldiversification1A","Constant Rate",as.character(bisp_long6$param_names))
bisp_long6$param_names<-gsub("diversification0A","Trait Absent",as.character(bisp_long6$param_names))
bisp_long6$param_names<-gsub("diversification1A","Trait Present",as.character(bisp_long6$param_names))
bisp_long6$Parameter<-bisp_long6$param_names

# Create unique identifier
bisp_long6$uni3<-paste0(bisp_long6$uni2,bisp_long6$Parameter)

# Remove duplicates
bisp_long6<-bisp_long6[!duplicated(bisp_long6$uni3),]

hist(bisp_long6$Parameters)

# Fig 3 Main  ------------------------------------------------------

write.csv(bisp_long6,
          file=here("outputs","BISSE_Plots3_S6_data.csv"), 
          row.names = FALSE)

# Plot results of bisse
p3 = ggplot(data=bisp_long6[bisp_long6$Parameter!="Constant Rate" & bisp_long6$location=="Global" ,],aes(x=variable,y=Parameters,col=Parameter) )+
  geom_point(fill=NA,position = position_jitterdodge(jitter.width = 0.3))+
  geom_boxplot(fill=NA,lwd=1,outlier.size = 0) +
  theme_cowplot(14)+
  #facet_wrap(.~location,scales="free")+
  #geom_hline(yintercept=0,lty=2)+
  #theme(axis.text.x=element_text(angle=45,hjust=1))+
  xlab("Trait")+
  #scale_x_discrete(guide = guide_axis(n.dodge=2))+
  ylab("Diversification Rate")+
  #ggtitle("B.")+
  scale_colour_manual(values = c("#ef8a6275","#377eb875"))#+
#theme(legend.position = c(0.05, 0.875))

ggsave(p3, filename = here("outputs", "bisse_main_figure.pdf"), width = 6, height = 4)

# Fig S6 Supplementary ---------------------------------------------
###plot results of bisse
p6 = ggplot(data=bisp_long6[bisp_long6$Parameter!="Constant Rate" & bisp_long6$location!="Global"  ,],aes(x=variable,y=Parameters,col=Parameter) )+
  geom_point(fill=NA,position = position_jitterdodge(jitter.width = 0.3))+
  geom_boxplot(fill=NA,lwd=1,outlier.size = 0) +
  theme_cowplot(14)+
  facet_wrap(.~location,ncol=2,scales="free")+
  #geom_hline(yintercept=0,lty=2)+
  #theme(axis.text.x=element_text(angle=45,hjust=1))+
  xlab("Trait")+
  #scale_x_discrete(guide = guide_axis(n.dodge=2))+
  ylab("Diversification Rate")+
  #ggtitle("B.")+
  scale_colour_manual(values = c("#ef8a6275","#377eb875"))+
  theme(legend.position = c(0.85, 0.94))

print(p6)

ggsave(p6, filename = here("outputs", "Fig_S6_square.png"), width = 10, height = 10, dpi = 600)
ggsave(p6, filename = here("outputs", "Fig_S6_narrow.pdf"), width = 5, height = 20, dpi = 600)

###comparison to NULL

##get rate variable model
lfB<-lf3[lf3$name=="BiSSE",]
lfB$param_names<-c(xxx,"no_name")
lfB$mer1<-paste0(lfB$tree,lfB$variable)

##get rate invariable modle
lfBn<-lf3[lf3$name=="dull_null",]

##check same
identical(lfB$tree,lfBn$tree)
identical(lfB$variable,lfBn$variable)

###make merge column
lfBn$mer1<-paste0(lfBn$tree,lfBn$variable)

##change names
names(lfBn)<-paste0("null_",names(lfBn))

##combine
bisse2<-merge(lfB,lfBn,by.x="mer1",by.y="null_mer1",all.x=T,all.y=F)
bisse2$delta_AIC<-bisse2$AIC-bisse2$null_AIC
bisse2$dummy=1

##create unique columns again
bisse2$uni2<-paste(bisse2$tree,bisse2$variable,bisse2$location,sep="_")
bisse2$uni3<-paste(bisse2$name,bisse2$variable,sep="_")

###make AIC
bisse2$lessAIC<-0
bisse2$lessAIC[(bisse2$delta_AIC)<=(-2)]<-1
bisse2$moreAIC<-(1-bisse2$lessAIC)


pdf(file=".\\outputs\\bisseAIC_distributions.pdf",height=8,width=12)

##make plot
ggplot(bisse2[bisse2$param_names=="turnover0A" ,], aes(x=delta_AIC,y=fct_reorder(variable,delta_AIC),fill=variable)) +
  geom_density_ridges(scale = 1)+
  facet_wrap(.~location,scales="free")+
  geom_vline(xintercept = 0,lty=2 )+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_cowplot(12)+
  xlab("Delta AIC")+
  theme(axis.title.y = element_blank(),legend.position = "none")

dev.off()

# Rename
lf3$Test<-lf3$name

# Plot differences in parameter rate estimates
pdf(file=".\\outputs\\bisse_null_comparison.pdf",height=10,width=16)

ggplot(lf3[((lf3$name=="BiSSE" |lf3$name=="dull_null") ) & lf3$variable!="Marriage_Endogamous",], aes(x=AIC,fill=Test)) +
  geom_density(adjust=2)+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_cowplot(12)+
  xlab("AIC")+
  theme(legend.position = "top")+
  facet_wrap(.~location+variable,scales="free")
dev.off()

# Main results
# Get rid of parameters for now
bisse3<-bisse2[!duplicated(bisse2$uni2) ,]

# Aggregate results
aic_diff = -2
res2<-aggregate(bisse3$dummy[bisse3$delta_AIC<=aic_diff],by=list(bisse3$variable[bisse3$delta_AIC<=aic_diff],bisse3$location[bisse3$delta_AIC<=aic_diff]),length)
res2a<-aggregate(bisse3$dummy,by=list(bisse3$variable,bisse3$location),length)
res2$prop1<-res2$x/res2a$x
names(res2)<-c("Variable","Location","Number","Proportion")
write.csv(res2,file=here("outputs", "bisse_null_AIC_proportions.csv"))
