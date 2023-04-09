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

### loading ordination results and calculations ###
load("ordination_objects.Rdata")

### flipping axes s posteriori to fit with facet grid gradients  ###
DCA_site$DCA1[DCA_site$weighting=="Presence-absence (R = 0)"] <- 
  -DCA_site$DCA1[DCA_site$weighting=="Presence-absence (R = 0)"]

DCA_site$DCA2[DCA_site$weighting=="R = 8"|DCA_site$weighting=="R = 2"] <- 
  -DCA_site$DCA2[DCA_site$weighting=="R = 8"|DCA_site$weighting=="R = 2"]

NMDS_site$NMDS2[NMDS_site$weighting=="Presence-absence (R = 0)"] <- 
  -NMDS_site$NMDS2[NMDS_site$weighting=="Presence-absence (R = 0)"]

#species
DCA_species$DCA1[DCA_species$weighting=="Presence-absence (R = 0)"] <- 
  -DCA_species$DCA1[DCA_species$weighting=="Presence-absence (R = 0)"]

DCA_species$DCA2[DCA_species$weighting=="R = 8"|DCA_species$weighting=="R = 2"] <- 
  -DCA_species$DCA2[DCA_species$weighting=="R = 8"|DCA_species$weighting=="R = 2"]

NMDS_species$NMDS2[NMDS_species$weighting=="Presence-absence (R = 0)"] <- 
  -NMDS_species$NMDS2[NMDS_species$weighting=="Presence-absence (R = 0)"]

### facet grid of weightings and localities (axes 1 & 2) ###
# DCA
ggplot(data=DCA_site) +
  geom_point(aes(x=DCA1, y=DCA2, color=locality)) +
  facet_wrap(~weighting)
ggsave("DCA_weights.png", device="png", width=10, height=6)

# NMDS
ggplot(data=NMDS_site) +
  geom_point(aes(x=NMDS1, y=NMDS2, color=locality)) +
  facet_wrap(~weighting)
ggsave("NMDS_weights.png", device="png", width=10, height=6)

### "path" of each sample between the different weightings ###
#DCA
ggplot(data=DCA_site, aes(x=DCA1, y=DCA2, 
                          group = PMO)) +
  geom_path(color="black") +
  geom_point(aes(color=weighting, shape=locality)) +
  scale_color_brewer(palette = "PiYG")
ggsave("DCA_path.png", device="png", width=9, height=6)

#NMDS
ggplot(data=NMDS_melt, aes(x=NMDS1, y=NMDS2, 
                           group = PMO)) +
  geom_path(color="black") +
  geom_point(aes(color=weighting, shape=locality)) +
  
  scale_color_brewer(palette = "PiYG")
ggsave("NMDS_path.png", device="png", width=9, height=6)

DCA_species$taxa
### "path" of each species between the different weightings ###
#DCA
DCA_species$taxa
ggplot() +
  geom_path(data=DCA_species, 
            aes(x=DCA1, y=DCA2, group=taxa),
            color="grey") +
  geom_point(data=DCA_species, 
             aes(x=DCA1, y=DCA2, color=weighting),
             size=0.8) +
#  geom_image(data=norm_species, 
 #            aes(x=DCA1_8, y=-DCA2_8),
 #            image=img, size=img_size*0.7) +
  geom_image(data=norm_species, 
             aes(x=DCA1_r, y=DCA2_r),
             image=img, size=img_size*0.7) +
  scale_color_brewer(palette = "PiYG") +
  theme_minimal()
#  scale_color_brewer(palette = "PiYG")
ggsave("DCA_path_s.png", device="png", width=9, height=6)

#NMDS
ggplot() +
  geom_path(data=NMDS_species, 
            aes(x=NMDS1, y=NMDS2, group=taxa),
            color="grey") +
  geom_point(data=NMDS_species, 
             aes(x=NMDS1, y=NMDS2, color=weighting),
             size=0.8) +
#  geom_image(data=norm_species, 
 #            aes(x=NMDS1_pa, y=-NMDS2_pa),
  #           image=img, size=img_size) +
    geom_image(data=norm_species, 
              aes(x=NMDS1_r, y=NMDS2_r),
              image=img, size=img_size*0.7) +
  scale_color_brewer(palette = "PiYG") +
  theme_minimal()
  
  scale_color_brewer(palette = "PiYG")
ggsave("NMDS_path_s.png", device="png", width=9, height=6)


### change in stratigraphic correlation ###
ggplot(data=DCA_site) +
  geom_point(aes(x=DCA1, y=height_m, color=locality)) +
  facet_wrap(~weighting)
ggsave("DCA_weights_h.png", device="png", width=10, height=6)

ggplot(data=NMDS_melt) +
  geom_point(aes(x=NMDS1, y=height_m, color=locality)) +
  facet_wrap(~weighting)
ggsave("NMDS_weights_h.png", device="png", width=10, height=6)

### comparison of abundance_weighting ###
### procrustes comparison of ordinations ###
protest(dca.r, nmds.r)
protest(dca.16, nmds.16)
protest(dca.8, nmds.8)
protest(dca.4, nmds.4)
protest(dca.2, nmds.2)
protest(dca.pa, nmds.pa)

### kentalls t-correlation of axes###
cor.test(all_data$DCA1_r, all_data$NMDS1_r, method="k")
cor.test(all_data$DCA2_r, all_data$NMDS2_r, method="k")

cor.test(all_data$DCA1_16, all_data$NMDS1_16, method="k")
cor.test(all_data$DCA2_16, all_data$NMDS2_16, method="k")

cor.test(all_data$DCA1_8, all_data$NMDS1_8, method="k")
cor.test(all_data$DCA2_8, all_data$NMDS2_8, method="k")

cor.test(all_data$DCA1_4, all_data$NMDS1_4, method="k")
cor.test(all_data$DCA2_4, all_data$NMDS2_4, method="k")

cor.test(all_data$DCA1_2, all_data$NMDS1_2, method="k")
cor.test(all_data$DCA2_2, all_data$NMDS2_2, method="k")

cor.test(all_data$DCA1_pa, all_data$NMDS1_pa, method="k")
cor.test(all_data$DCA2_pa, all_data$NMDS2_pa, method="k")

### NMDS stress ###
nmds.r$stress
nmds.16$stress
nmds.8$stress
nmds.4$stress
nmds.2$stress
nmds.pa$stress

