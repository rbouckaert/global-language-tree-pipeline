# Code to generate figure 1 data plot

# Files and folders read or loaded:
# - input_data/final_env_datapoints.csv
# - input_data/taxonorder-left-to-right.txt
# - input_data/bisse_data_in.csv


library(pacman)

pacman::p_load(ggplot2, cowplot, data.table, aplot, here, patchwork, magick)


# -----------------------------------------------------------------------------|
# 1) Load data ----
# -----------------------------------------------------------------------------|

env_data<-fread(here("input_data",'fixed', "final_env_datapoints.csv"))

env_data$LANG_IS[env_data$LANG_IS=="osse1243"]<-"iron1242"

# Choose 1 date
env_data2<-env_data[env_data$Time3 ==max(env_data$Time3),]
row.names(env_data2)<-env_data2$LANG_IS

# Get regional order
ord1<-read.table(file=here("input_data","taxonorder-left-to-right.txt"),
                 stringsAsFactors = FALSE)

# Read the file, treating each line as one row
data <- fread(here("outputs", "taxonorder5.txt"), header = FALSE)

# Flatten the data, if values are space-separated in a single row
values <- unlist(strsplit(as.character(unlist(data)), "\\s+")) # Splits by any whitespace

# Create a data frame where each value is a row
result <- data.frame(V1 = values)
ord1 = result
names(ord1)<-"glottocode"
env_data6<-env_data2
setkey(env_data6,LANG_IS)

# Load up tip data
res5<-fread(file=here("input_data", "bisse_data_in.csv"))
setkey(res5,Glottocode)
res6<-env_data6[res5]
setkey(res6,LANG_IS)

# Reorder data.frame
res6<-res6[order(match(res6$LANG_IS,ord1$glottocode)),]
res6$glottocode<-factor(res6$LANG_IS,levels=res6$LANG_IS)

# -----------------------------------------------------------------------------|
# 2) Create individual plots ----
# -----------------------------------------------------------------------------|

area1<-ggplot(res6,aes(x=glottocode ,y=log(Shap_Ar+1)))+
  geom_bar(width=1,stat="identity",col='#2c7fb850')+
  theme_cowplot(14)+
  xlab("Spoken area (log km\u00B2)")+
  ylab("")+
  theme(axis.text.x =element_blank(),axis.ticks.x =element_blank(),
        legend.position = "none",axis.line = element_blank())+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))


popd1<-ggplot(res6,aes(x=glottocode ,y=log(popd+1)))+
  geom_bar(width=1,stat="identity",col='#2c7fb850')+
  theme_cowplot(14)+
  xlab("Human Pop. Density (log people per km\u00B2)")+
  ylab("")+
  theme(axis.text.x =element_blank(),axis.ticks.x =element_blank(),
        legend.position = "none",axis.line = element_blank())+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))

friction1<-ggplot(res6,aes(x=glottocode ,y=(friction)))+
  geom_bar(width=1,stat="identity",col='#2c7fb850')+
  theme_cowplot(14)+
  xlab("Landscape Friction (minutes/meter)")+
  ylab("")+
  theme(axis.text.x =element_blank(),axis.ticks.x =element_blank(),
        legend.position = "none",axis.line = element_blank())+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))+
  ylim(0,0.35)


distancetocityyear1<-ggplot(res6,aes(x=glottocode ,y=(distancetocityyear/1000000)))+
  geom_bar(width=1,stat="identity",col='#2c7fb850')+
  theme_cowplot(14)+
  xlab("Distance to City (millions of meters)")+
  ylab("")+
  theme(axis.text.x =element_blank(),axis.ticks.x =element_blank(),legend.position = "none",axis.line = element_blank())+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))

res6$dummy=1

island1 <- ggplot(res6, aes(x = as.factor(glottocode), y = dummy, 
                            fill = as.factor(ifelse(island == 1, "1", "0")))) +
  geom_bar(stat = "identity", width = 1, position=position_dodge(1)) +
  scale_fill_manual(values = c('0' = "#FFFFFA", '1' = '#2c7fb8')) + # Grey for 0, Blue for 1
  theme_cowplot(14)+
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Island Language") +
  ylab("") +
  labs(tag = "C.")+
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    axis.line = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
     plot.tag.position = c(0, 1.2),
    plot.tag = element_text(face = "bold", size = 16)
  )

# -----------------------------------------------------------------------------|
# 3) Combine subplots ----
# -----------------------------------------------------------------------------|

# Combine the plot
subplot_C <- area1 %>% 
  insert_top(distancetocityyear1) %>% insert_top(friction1) %>% 
  insert_top(popd1) %>% insert_top(island1)

# Save the plot using ggsave
# ggsave(filename = here("outputs", "Fig_1C.pdf"),  # File path
#        plot = subplot_C, height = 6, width = 14, units = "in")
# 
# ggsave(filename = here("outputs", "Fig_1C.png"),  # File path
#        plot = subplot_C, height = 6, width = 14, units = "in")


# Load the PNG image using magick
image <- image_read(here("input_data", "edge6636-densitree5cropped_cleaned.png"))

# Convert the image to a raster object
raster_image <- grid::rasterGrob(image, interpolate = TRUE)

# Create a ggplot with the image and annotations
subplot_A <- ggplot() +
  annotation_custom(raster_image, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +  # Force full-width alignment
  annotate(
    "text",
    x = 2, y = 2, label = "50",  # Custom annotation
    size = 4.5, hjust = 28, vjust=6.5, 
  ) +
  annotate(
    "text",
    x = 2, y = 2, label = "100",  # Custom annotation
    size = 4.5, hjust = 18.7, vjust=-2.9, 
  ) +
  annotate(
    "text",
    x = 2, y = 2, label = "150",  # Custom annotation
    size = 4.5, hjust = 18.7, vjust=-12.5,
  ) +
  annotate(
    "text",
    x = 2, y = 2, label = "0",  # Custom annotation
    size = 4.5, hjust = 55.5, vjust=17.4, 
  ) +
  labs(tag = "A.")+
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.tag.position = c(0.01, 1),
    plot.tag = element_text(face = "bold", size = 16, vjust = -0.5)
  
  ) 


# Use patchwork to combine the raster plot with ggplots
final_plot <- plot_spacer()+ subplot_A +
  plot_spacer()+
  island1 +
  popd1 +
  friction1 +
  distancetocityyear1 +
  area1 +
  plot_layout(
    ncol = 1,
    heights = c(0.3,3.5,0.02,0.9, 0.9, 0.9, 0.9, 0.9)  # Allocate more space to the DensiTree plot
  )

# -----------------------------------------------------------------------------|
# 4) Subplot B -------------------------------------------------------------------
# -----------------------------------------------------------------------------|
ltt1 = read.csv(here("outputs","Fig_1","LTT_Plot1_Data.csv"))

# Force 'cont' into a factor with the levels in the desired order:
ltt1$cont <- factor(ltt1$cont,
                    levels = c("Africa", "Eurasia", "Oceania",
                               "North America", "South America", "Global"))

# Plot LTT data
p4<-ggplot(ltt1[,], aes(y=logltt,x=times,col=cont,group=tree_id_1_901))+
  geom_line()+
  scale_colour_manual(values = c('#7570b3','#e7298a','#d95f02','#66a61e','#e6ab02','black'))+
  facet_wrap(.~cont,ncol=3)+
  theme_cowplot(12)+
  ggtitle(label="B.")+
  theme(legend.position = "none",
        plot.title = element_text(size=17, vjust = -3, hjust = -0.1),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(angle = 0, hjust = 0))+
  ylab("Lineages (log)")+ 
  xlab("Years Ago (KY)")
  
ggsave(here("outputs","Fig_1","LTT_Plot1_2.png"),plot=p4,height=6,width=12, dpi=600)

# -----------------------------------------------------------------------------|
# 5) Incorpotate subplot B into final_plot top right corner-----------------
# -----------------------------------------------------------------------------|

final_plot_inset <- ggdraw() +
  draw_plot(final_plot, 0, 0, 1, 1) +  # Main plot occupies full canvas
  draw_plot(p4, 0.48, 0.74, 0.5, 0.25)  # Smaller inset at top-right

# ggsave(
#   filename = here("outputs","Fig_1", paste0("figure_with_inset_cowplot_", sample(1:10000, 1), ".png")),
#   plot     = final_plot_inset,
#   height   = 13,
#   width    = 12,
#   units    = "in",
#   dpi      = 600
# )

# -----------------------------------------------------------------------------|
# Add vertical text to the final plot inset
final_plot_with_text <- final_plot_inset +
  annotate(
    "text",
    x = 0.05,               # X-coordinate 
    y = 0.8,               # Y-coordinate 
    label = "Years Ago (KY)", # The text to display
    angle = 90,            # Makes the text vertical
    size = 5,              # Adjust size as needed
    colour = "black",      # Text colour
    hjust = 0.5,           # Adjust horizontal justification
    vjust = 0.5            # Adjust vertical justification
  )

# PNG
ggsave(
  filename = here("outputs","Fig_1", paste0("Fig_1", sample(1:10000, 1), ".png")),
  plot     = final_plot_with_text,
  height   = 13,
  width    = 12,
  units    = "in",
  dpi      = 600
)

# PDF
ggsave(
  filename = here("outputs","Fig_1", paste0("Fig_1", sample(1:10000, 1), ".pdf")),
  plot     = final_plot_with_text,
  height   = 13,
  width    = 12,
  units    = "in",
  dpi      = 600
)
