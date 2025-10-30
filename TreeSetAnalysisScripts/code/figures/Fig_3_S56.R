# Code to create results and plots for Fig3 and FigS5/S6

# Load libraries ----------------------------------------------------------------
library(here)
library(pacman)
p_load(ggplot, cowplot, ggpubr, gridExtra, RColorBrewer, ggridges, 
       data.table, forcats, ggh4x, magrittr, here)

# Load data --------------------------------------------------------------------
# Load bisse results data
lf3 = read.csv(here("outputs", "bisse_results1.csv"))

# Load BISSE Simulation data
bisse_sim_files = list.files("outputs/BISSE_SIM", full.names = TRUE)
for (i in 1:length(bisse_sim_files)) {
  bisse_sim = read.csv(bisse_sim_files[i])
  bisse_sim$number = gsub("_Global_BISSE_runs3.csv", "", basename(bisse_sim_files[i]))
  bisse_sim$number = gsub("es4_simulationes_extant_Global_all_", "", bisse_sim$number)
  
  # Combine with existing data
  if (i == 1) {
    bisse_sim_combined = bisse_sim
  } else {
    bisse_sim_combined = rbind(bisse_sim_combined, bisse_sim)
  }
}

# Clean/relabel data --------------------------------------------------------------------
# Add metadata to BISSE Simulation data
bisse_sim_combined$location = "Global"
bisse_sim_combined$type = "simulated"
bisse_sim_combined$xxxx = "Global"
bisse_sim_combined$tree2 = "tree"

# Check for missing columns
names(lf3)[!names(lf3) %in% names(bisse_sim_combined)]

#Combine with DR restuls from non-simulated trees
DR1 = rbind(lf3, bisse_sim_combined)

# Load parameter names object called xxx
load(here("input_data", "bisse_param_names.r"))

# Specify what types of test
types<-unique(DR1$name)
print(types)

# Subset and rename DR
DR1=DR1[DR1$name=="DR_test",]
DR1$variable<-gsub("Marriage_Endogamous","Endogamous Marriage",DR1$variable)
DR1$variable<-gsub("Marriage_Exogamous","Exogamous Marriage",DR1$variable)
DR1$variable<-gsub("D_Place_","",DR1$variable)
DR1$variable<-gsub(" definite","",DR1$variable)
DR1$variable<-gsub("island","Island",DR1$variable)
DR1$variable<-gsub("_Pattern","",DR1$variable)
DR1$variable<-gsub("EA033_Extended","Political Complexity",DR1$variable)

DR1 %<>% filter(location=="Global")
DR1$type <- factor(DR1$type, levels = c("simulated", "all"))


# Figure S5 --------------------------------------------------------------------
# a) 2 categories: Updated code for faceted ridges -----
p5 <- ggplot(DR1,  # Use the full dataset, not just the subset
             aes(
               x = AIC, 
               y = fct_reorder(variable, AIC), 
               fill = as.factor(type),
               alpha = as.factor(type)# Map fill to Parameters
             )) +
  geom_density_ridges(bandwidth = 0.015, color = NA) +  # Use alpha for transparency
  facet_wrap(. ~ variable, scales = "free") +  # Facet by variable
  geom_vline(xintercept = 0, lty = 2) +  # Vertical line at x=0
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c("black","#FF5733"),
                    name = "Dataset",
                    labels = c("Simulated","Observed")) +  # Custom fill colors
  scale_alpha_manual(values = c(0.6, 0.6),
                     labels = c("Simulated","Observed")) +
  theme_cowplot(12) +
  labs(fill = "Dataset", alpha = "Dataset")+
  xlab("Difference in DR") +
  theme(legend.text = element_text(size = 9),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    legend.position = c(0.9, 0.82)  # Add a legend for clarity
  )

print(p5)

ggsave(here("outputs","Fig_S5", "Fig_S5_new.png"), p5, height=3, width=8, dpi = 600)
ggsave(here("outputs","Fig_S5", "Fig_S5_new.pdf"), p5, height=3, width=8, dpi = 600)

# Leave-one-out ridges ---------------------------------------------------------
DR1$location <- factor(DR1$location, levels = c("Global","Africa", "America","Eurasia", "Oceania"))
cols = c("black","#EE0000FF","#008B45FF","#631879FF","#008280FF")

p5r = ggplot(DR1[DR1$Parameters==0 & DR1$variable!="Endogamous Marriage",], 
            aes(x=AIC,y=factor(location),fill=factor(location), alpha = factor(location))) +
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

print(p5r)

ggsave(here("outputs","Fig_S5", "Fig_S5_ridges.png"), p5r, height=3, width=8, dpi = 600)

# Fig 3 Main --------------------------------------------------------------------
# Load BISSE plots data
bisp_long6 = read.csv(here("outputs", "BISSE_Plots3_S6_data.csv"))

# Plot results of bisse
p3 = ggplot(data=bisp_long6[bisp_long6$Parameter!="Constant Rate" & bisp_long6$location=="Global" ,],
            aes(x=variable,y=Parameters,col=Parameter) )+
  geom_point(fill=NA,position = position_jitterdodge(jitter.width = 0.3))+
  geom_boxplot(fill=NA,lwd=1,outlier.size = 0) +
  theme_cowplot(11)+
  xlab("Trait")+
  ylab("Diversification Rate")+
  scale_colour_manual(values = c("#ef8a6275","#377eb875"))+
  theme(legend.position = c(0.77, 0.875))

# Save PNG and PDF
ggsave(p3, filename = here("outputs","Fig_3", "Fig_3.pdf"), width = 5, height = 3.5, dpi = 1000)
ggsave(p3, filename = here("outputs","Fig_3", "Fig_3.png"), width = 5, height = 3.5, dpi = 1000)

# Fig S6 Supplementary --------------------------------------------------------------------
# Plot results of bisse

p6 = ggplot(data=bisp_long6[bisp_long6$Parameter!="Constant Rate" & bisp_long6$location!="Global"  ,],aes(x=variable,y=Parameters,col=Parameter) )+
  geom_point(fill=NA,position = position_jitterdodge(jitter.width = 0.3))+
  geom_boxplot(fill=NA,lwd=1,outlier.size = 0) +
  theme_cowplot(14)+
  facet_wrap(.~location,ncol=2,scales="free")+
  xlab("Trait")+
  ylab("Diversification Rate")+
  scale_colour_manual(values = c("#ef8a6275","#377eb875"))+
  theme(legend.position = c(0.85, 0.94))

# Save PNG and PDF
ggsave(p6, filename = here("outputs","Fig_S6", "Fig_S6_square.png"), width = 10, height = 10, dpi = 600)
ggsave(p6, filename = here("outputs","Fig_S6", "Fig_S6_narrow.pdf"), width = 5, height = 20, dpi = 600)

