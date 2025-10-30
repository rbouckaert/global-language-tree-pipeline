library(here)
library(ggplot2)
library(tidyverse)
library(viridis)
library(forcats)
library(ggridges)
library(cowplot)

# Extended Data Fig. 2 ------------------------------------------------------

# Data Import --------------------------------------------------------------------
fil3 = read.csv(here("outputs", "All_Trees_Summary.csv"), 
                stringsAsFactors = FALSE)
fil3<-as.data.frame(fil3)

# Make grouping column group2 --------------------------------------------------------------------
grps1<-read.csv(here("input_data", "groups_table.csv"), stringsAsFactors = FALSE)

names(grps1)[3]<-"group2"
fil3<-merge(fil3,grps1[,c("group","group2")],by="group",all.x=TRUE,all.y=FALSE)

# Check all matched --------------------------------------------------------------------
nrow(fil3[is.na(fil3$group2),])

# Make group more legible --------------------------------------------------------------------
fil3$name<-gsub("_"," ",fil3$group)
fil3$name<-gsub("Bantu (Atlantic-Congo)","Bantu-Atlantic-Congo",fil3$name)
fil3$name<-gsub("Central-Sudanic","Central Sudanic",fil3$name)
fil3$name<-gsub("North-America","North America",fil3$name)
fil3$name<-gsub("South-America","South America",fil3$name)
fil3$name<-gsub("Nuclear-Trans-New-Guinea","Nuclear Trans New Guinea",fil3$name)
fil3$name<-gsub("Nuclear-Torricelli","Nuclear Torricelli",fil3$name)

p = ggplot(fil3[fil3$group2=="top27families",],
           aes(x=HmeanDR,y=fct_reorder(name,HmeanDR),
               fill=name)) +
  geom_density_ridges()+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_cowplot(20)+
  xlab("Harmonic Mean\nDiversification Rate")+
  theme(axis.title.y = element_blank(),legend.position = "none")

# Save the plot as a high-quality PNG
ggsave(
  filename = here("outputs","Fig_S2","Fig_S2.png"), 
  plot = p,
  width = 12, 
  height = 8, 
  dpi = 600,  
  scale = 1
)

# Save the plot as a high-quality PDF
ggsave(
  filename = here("outputs","Fig_S2","Fig_S2.pdf" ), 
  plot = p,
  width = 12,           
  height = 8, 
  device = cairo_pdf,              
  scale = 1
)

# Alternative version - with mean-based fill ----------------------------------
# Calculate mean for each group
data_with_means <- fil3 %>%
  filter(group2 == "top27families") %>%  # Filter for the desired group
  group_by(name) %>%
  mutate(mean_value = mean(HmeanDR))     # Compute mean HmeanDR for each group

# Create the ridge plot with mean-based fill
p2 <- ggplot(data_with_means, 
            aes(x = HmeanDR, 
                y = fct_reorder(name, HmeanDR), 
                fill = mean_value)) +  # Use mean_value for fill
  geom_density_ridges() +
  scale_fill_viridis_c(option = "C", direction = -1) +  # Viridis colour scale
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_cowplot(20) +
  xlab("Harmonic Mean\nDiversification Rate") +
  theme(
    axis.title.y = element_blank(),
    legend.position = "none"  # Optional: move legend to the right
  ) +
  labs(fill = "Mean Value")  # Add legend title

# Print the plot
print(p2)


# Save the plot as a high-quality PNG
ggsave(
  filename = here("outputs","Fig_S2","Fig_S2_viridis.png"), 
  plot = p2,
  width = 12, 
  height = 8, 
  dpi = 600,  
  scale = 1
)

# Save the plot as a high-quality PNG
ggsave(
  filename = here("outputs","Fig_S2","Fig_S2_viridis.pdf"), 
  plot = p2,
  width = 12, 
  height = 8, 
  dpi = 600,  
  scale = 1
)
