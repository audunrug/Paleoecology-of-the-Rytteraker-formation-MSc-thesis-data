# loading dependencies (NB! must be installed beforehand)
library(tidyverse)
library(vegan)
library(graphics)
library(gridExtra)
library(ggimage)
library(reshape2)
library(ggtext)
library(MIAmaxent)
library(magick)
library(goeveg)
library(scales)

### loading ordination results and calculations ###
load("ordination_objects.Rdata")

### plotting orientations (as a test) ###
ggplot(all_data) +
  geom_bar(position="stack",
           aes(x=locality, fill=orientation)) +
  scale_fill_brewer(palette = "Set3")

### species proportions and abundances ###
#selecting only species abundances
only_species <- norm_species %>% select(-img_md, -img, -img_size)

prop_plots <- 0
avg_abundance <- 0
for (i in 1:length(only_species$X1)) {
  prop_plots[i] <- (sum(only_species[i,] > 0))/(length(only_species[i,]))
  avg_abundance[i] <- mean(only_species[i,][only_species[i,] > 0])
}

norm_species$prop_plots <- prop_plots
norm_species$avg_abundance <- avg_abundance
rownames(norm_species)
norm_species$prop_plots
norm_species$avg_abundance
# plotting abundance-occupancy
ggplot(data=norm_species) +
  geom_image(aes(x=prop_plots, y=avg_abundance, image=img),
             size=norm_species$img_size) +
  scale_y_continuous(trans="log10") +
  labs(x="Proportion of samples", y="Mean abundance") +
  theme_minimal()
ggsave("abundance_occupancy.png", device="png", 
       dpi=300, units = "cm", width=20, height=15)


### boxplot of abundances for each species ###
species_abundances <- melt(all_data %>% select(locality, calc_tubes:strom))

ggplot(species_abundances) +
  geom_boxplot(aes(x=variable, y=value, color=locality)) +
  scale_x_discrete(name = NULL,
                   labels = img_md) +
  theme(axis.text.x = element_markdown(color = "black", size = 11),
        legend.position = "bottom")
ggsave("boxplot.png", device="png", dpi=300, units = "cm", 
       width=20, height=14)


### violin plot of strat. height and abundances ###
species_height_T <- melt(T2_data %>% select(height_m, calc_tubes:strom), 
                         id.vars=c("height_m"))

species_height_U <- melt(U_data %>% select(height_m, calc_tubes:strom), 
                         id.vars=c("height_m"))

#calculating organism means for each strat. interval with 2+ samples
species_height_T <- species_height_T %>% 
  group_by(height_m, variable) %>%
  summarise(mean=(mean(value)))
#species_height_T <- with(species_height_T,  species_height_T[variable, ])

species_height_U <- species_height_U %>% 
  group_by(height_m, variable) %>%
  summarise(mean=(mean(value)))
#species_height_U <- with(species_height_U,  species_height_U[variable, ])


#plotting violin plot
ggplot(species_height_T) +
  geom_violin(stat="identity",
              aes(x=variable, y=height_m, violinwidth=0.03*mean), 
              fill="salmon") +
  scale_x_discrete(name = NULL,
                   labels = img_md) +
  scale_y_continuous(name="Position (m.a.b)") +
  theme(axis.text.x = element_markdown(color = "black", size = 11))
ggsave("violinTOV2.png", device="png", dpi=300, units = "cm")

ggplot(species_height_U) +
  geom_violin(stat="identity",
              aes(x=variable, y=height_m, violinwidth=0.03*mean), 
              fill="steelblue") +
  scale_x_discrete(name = NULL,
                   labels = img_md) +
  scale_y_continuous(name="Position (m.a.b)") +
  theme(axis.text.x = element_markdown(color = "black", size = 11))
#ggsave("violinULV.svg", device="svg", dpi=300, units = "cm")
ggsave("violinULV.png", device="png", dpi=300, units = "cm")


### diversity metrics ###
#making dataframe for plotting
div_data <- all_data %>% 
  select(R, E, H, locality) %>%
  melt(id.vars="locality", variable.name="metric")

# boxplot of metrics
ggplot(data=div_data) +
  geom_boxplot(aes(x=locality, y=value, color=locality)) +
  facet_wrap(~metric, scales="free_y")
ggsave("diversity_boxplot.png", device="png", dpi=300, units = "cm")

# one-way ANOVA of diversity metrics
summary(aov(R ~ locality, data = all_data))  
summary(aov(H ~ locality, data = all_data))  
summary(aov(E ~ locality, data = all_data))

all_data$img_md <- img_md

norm_species$img_md[15]
### discussion plots ###
ggplot(data=all_data[all_data$locality == "TOV2",]) +
  geom_point(aes(y=brac_thin, x=height_m), color="steelblue") +
  geom_smooth(aes(y=brac_thin, x=height_m), method="loess", color="black") +
  scale_y_continuous(
    name = norm_species$img_md[15]) +
  coord_flip() +
  theme(axis.title.x = element_markdown(color = "black", size = 11),
        legend.position = "bottom")
ggsave("thinbrach.png", device="png", dpi=300, units = "cm", 
       width=5, height=12)

ggplot(data=all_data[all_data$locality == "TOV2",]) +
  geom_point(aes(y=brac_thick, x=height_m), color="steelblue") +
  geom_smooth(aes(y=brac_thick, x=height_m), method="loess", color="black") +
  scale_y_continuous(
    name = norm_species$img_md[14]) +
  coord_flip() +
  theme(axis.title.x = element_markdown(color = "black", size = 11),
        legend.position = "bottom")
ggsave("thickbrach.png", device="png", dpi=300, units = "cm", 
       width=5, height=12)


ggplot(data=all_data[all_data$locality == "TOV2",]) +
  geom_point(aes(y=bryo_small, x=height_m), color="steelblue") +
  geom_smooth(aes(y=bryo_small, x=height_m), method="loess", color="black") +
  scale_y_continuous(
    name = norm_species$img_md[2]) +
  coord_flip() +
  theme(axis.title.x = element_markdown(color = "black", size = 11),
        legend.position = "bottom")
ggsave("smallbryo.png", device="png", dpi=300, units = "cm", 
       width=5, height=12)

ggplot(data=all_data[all_data$locality == "TOV2",]) +
  geom_point(aes(y=ostracod, x=height_m), color="steelblue") +
  geom_smooth(aes(y=ostracod, x=height_m), method="loess", color="black") +
  scale_y_continuous(
    name = norm_species$img_md[13]) +
  coord_flip() +
  theme(axis.title.x = element_markdown(color = "black", size = 11),
        legend.position = "bottom")
ggsave("ostracod.png", device="png", dpi=300, units = "cm", 
       width=5, height=12)

ggplot(data=all_data[all_data$locality == "TOV2",]) +
  geom_point(aes(y=crinoid, x=height_m), color="steelblue") +
  geom_smooth(aes(y=crinoid, x=height_m), method="loess", color="black") +
  scale_y_continuous(
    name = norm_species$img_md[17]) +
  coord_flip() +
  theme(axis.title.x = element_markdown(color = "black", size = 11),
        legend.position = "bottom")
ggsave("crinoid.png", device="png", dpi=300, units = "cm", 
       width=5, height=12)


### DISCUSSION (P4) DIAGRAMS ###
ggplot(data=all_data[all_data$locality == "TOV2",]) +
  geom_point(aes(y=brac_thick, x=height_m)) +
  geom_smooth(aes(y=brac_thick, x=height_m), method="gam")

ggplot(data=all_data[all_data$locality == "TOV2",]) +
  #geom_point(aes(y=brach_wavy, x=height_m), xlim=c(0,10)) +
  geom_point(aes(y=brach_wavy, x=height_m))

ggplot(data=all_data[all_data$locality == "TOV2",]) +
  geom_point(aes(y=bryo_small, x=height_m)) +
  geom_smooth(aes(y=bryo_small, x=height_m), method="gam")

ggplot(data=all_data[all_data$locality == "TOV2",]) +
 # geom_point(aes(y=ostracod, x=height_m)) +
  geom_smooth(aes(y=ostracod, x=height_m))

ggplot(data=all_data) +
  geom_point(aes(y=H, x=height_m, color=locality)) +
  geom_smooth(aes(y=H, x=height_m, color=locality))

ggplot(data=all_data[all_data$locality == "TOV2",]) +
  geom_point(aes(y=bryo_large, x=height_m, color=locality)) +
  geom_smooth(aes(y=bryo_large, x=height_m))

ggplot(data=all_data[all_data$locality == "TOV2",]) +
#  geom_point(aes(y=crinoid, x=height_m), xlim=c(0,10)) +
  geom_smooth(aes(y=crinoid, x=height_m))

ggplot(data=all_data[all_data$locality == "TOV2",]) +
  geom_smooth(aes(y=brach_wavy, x=height_m)) +
  geom_point(aes(y=brach_wavy, x=height_m), xlim=c(0,10))  +
  coord_flip()

ggplot(data=all_data[all_data$locality == "TOV2",]) +
  geom_smooth(aes(y=grainstone, x=height_m)) +
  geom_point(aes(y=grainstone, x=height_m), xlim=c(0,10))  +
  coord_flip()
