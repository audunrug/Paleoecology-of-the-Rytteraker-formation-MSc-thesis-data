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
setwd("~/OneDrive - Universitetet i Oslo/master/statistics")
load("ordination_objects.Rdata")

### inspecting/plotting ordinations ###
### PCA ###
summary(pca.r) # finding axes eigenvalues
screeplot(pca.r) # displaying axes' contribution to variation

# PCA biplot
ggplot(data=all_data) +
  geom_point(aes(x=PCA1, y=-PCA2, color=locality)) +
  geom_image(data=norm_species, aes(x=PCA1, y=-PCA2),
             image=img, size=norm_species$img_size) +
  labs(x="PCA1 (16.6%)", y="PCA2 (9.0%)")
ggsave("PCA.png", device="png", dpi=300, units = "cm", 
       width=20, height=14)

# PCA1 + height
ggplot(data=all_data) +
  geom_point(aes(x=PCA1, y=height_m, color=locality)) +
  facet_wrap(~locality, scales="free_y") +
  labs(x="PCA1", y="Position (m.a.b)")
ggsave("PCA_h.png", device="png", dpi=300, units = "cm", 
       width=20, height=14)

### CA ###

summary(ca.r) # finding axes eigenvalues
screeplot(ca.r) # displaying axes' contribution to variation

# CA biplot
ggplot(data=all_data) +
  geom_point(aes(x=-CA1, y=-CA2, color=locality)) +
  geom_image(data=norm_species, aes(x=-CA1, y=-CA2),
             image=img, size=norm_species$img_size) +
  labs(x="-CA1 (18.3%)", y="-CA2 (14.8%)")
ggsave("CA.png", device="png", dpi=300, units = "cm", 
       width=20, height=14)

# CA1 + height
ggplot(data=all_data) +
  geom_point(aes(x=-CA1, y=height_m, color=locality)) +
  facet_wrap(~locality, scales="free_y") +
  labs(x="-CA1", y="Position (m.a.b)")
ggsave("CA_h.png", device="png", dpi=300, units = "cm", 
       width=20, height=14)


### DCA ###
summary(dca.r)   # NB! not statistically relevant
screeplot(dca.r) # NB! not statistically relevant

# DCA biplot
ggplot(data=all_data) +
  geom_point(aes(x=DCA1_r, y=DCA2_r, color=locality)) +
  geom_image(data=norm_species, aes(x=DCA1_r, y=DCA2_r),
             image=img, size=norm_species$img_size) +
  labs(x="DCA1", y="DCA2")
ggsave("DCA.png", device="png", dpi=300, units = "cm", 
       width=20, height=14)

# DCA1 + height
ggplot(data=all_data) +
  geom_point(aes(x=DCA1_r, y=height_m, color=locality)) +
  facet_wrap(~locality, scales="free_y") +
  labs(x="DCA1", y="Position (m.a.b)")
ggsave("DCA_h.png", device="png", dpi=300, units = "cm", 
       width=20, height=14)


### NMDS ###
nmds.r$stress #checking nmds stress
# NMDS biplot
ggplot(data=all_data) +
  geom_point(aes(x=NMDS1_r, y=NMDS2_r, color=locality)) +
  geom_image(data=norm_species, aes(x=NMDS1_r, y=NMDS2_r),
             image=img, size=img_size) +
  labs(x="NMDS1", y="-NMDS2")
ggsave("NMDS.png", device="png", dpi=300, units = "cm", 
       width=20, height=14)

# NMDS1 + height
ggplot(data=all_data) +
  geom_point(aes(x=NMDS1_r, y=height_m, color=locality)) +
  facet_wrap(~locality, scales="free_y") +
  labs(x="NMDS1", y="Position (m.a.b)")
ggsave("NMDS_h.png", device="png", dpi=300, units = "cm", 
       width=20, height=14)
all_data[all_data$locality == "TOV2"]

ggplot(data=all_data[all_data$locality == "TOV2",]) +
  geom_point(aes(y=DCA1_r, x=height_m, color=locality)) +
  geom_smooth(aes(y=DCA1_r, x=height_m)) +
  facet_wrap(~locality, scales="free") +
  labs(x="NMDS1", y="Position (m.a.b)") +
  coord_flip()


ggsave("NMDS_h.png", device="png", dpi=300, units = "cm", 
       width=20, height=14)

### comparison ###
#making facet wrap dataframe
facet_data <- all_data
facet_data$CA1 <- -facet_data$CA1

#flipping axes 2 to fit visually with the same gradients
facet_data2 <- all_data
facet_data2$CA2 <- -facet_data2$CA2
facet_data2$PCA2 <- -facet_data2$PCA2

#melting the data frames
facet_data <- melt(facet_data, 
                   measure.vars = c("PCA1", "CA1", "DCA1_r", "NMDS1_r"),
                   variable="ordination",
                   value.name="axis_1")
facet_data2 <- melt(facet_data2, 
                    measure.vars = c("PCA2", "CA2", "DCA2_r", "NMDS2_r"),
                    variable="ordination",
                    value.name="axis_2")
facet_data$axis_2 <- facet_data2$axis_2

#2D diagrams
ggplot(data=facet_data) +
  geom_point(aes(x=axis_1, y=axis_2, color=locality)) +
  scale_size_continuous(range = c(0.05, 4)) +
  facet_wrap(~ordination, scales="free",
             labeller = as_labeller(
               c(PCA1 = "PCA", CA1 = "CA",
                 DCA1_r = "DCA", NMDS1_r="NMDS"))) +
  labs(x="Axis 1", y="Axis 2")
ggsave("ord_grid.png", device="png", dpi=300, units = "cm", 
       width=20, height=16)


#mspecies score comparison
facet_data <- norm_species
facet_data$CA1 <- -norm_species$CA1

#flipping axes 2 to fit visually with the same gradients
facet_data2 <- norm_species
facet_data2$CA2 <- -norm_species$CA2
facet_data2$PCA2 <- -norm_species$PCA2

#melting the data frames
facet_data <- melt(facet_data, 
                   measure.vars = c("PCA1", "CA1", "DCA1_r", "NMDS1_r"),
                   variable="ordination",
                   value.name="axis_1")
facet_data2 <- melt(facet_data2, 
                    measure.vars = c("PCA2", "CA2", "DCA2_r", "NMDS2_r"),
                    variable="ordination",
                    value.name="axis_2")
facet_data$axis_2 <- facet_data2$axis_2

facet_data$img <- rep(img, 4)
facet_data$img_md <- rep(img_md, 4)

#2D diagrams
ggplot(data=facet_data) +
  geom_image(aes(x=axis_1, y=axis_2),
             image=img_2, size=img_size_2*2) +
  scale_size_continuous(range = c(0.05, 4)) +
  facet_wrap(~ordination, scales="free",
             labeller = as_labeller(
               c(PCA1 = "PCA", CA1 = "CA",
                 DCA1_r = "DCA", NMDS1_r="NMDS"))) +
  labs(x="Axis 1", y="Axis 2")
ggsave("ord_grid_spec.png", device="png", dpi=300, units = "cm", 
       width=20, height=16)


#height and axis 1
facet_data$locality <- facet_data$locality %>% 
  fct_collapse("ULV1+2" = c("ULV1","ULV2"))
ggplot(data=facet_data) +
  geom_point(aes(x=axis_1, y=height_m, color=locality)) +
  facet_grid(locality~ordination, scales="free",
             labeller = as_labeller(
               c(PCA1 = "PCA", CA1 = "CA",
                 DCA1_r = "DCA", NMDS1_r="NMDS",
                 TOV2="TOV2", "ULV1+2"="ULV1+2"))) +
  labs(x="Axis 1", y="Position (m.a.b)")
ggsave("ord_grid_h.png", device="png", dpi=300, units = "cm", 
       width=20, height=16)

ggplot(data=facet_data[facet_data$locality=="TOV2",]) +
  geom_point(aes(y=axis_1, x=height_m, color=locality)) +
  geom_smooth(aes(y=axis_1, x=height_m), method = "loess", color="black") +
  facet_wrap(~ordination, scales="free", nrow = 1,
             labeller = as_labeller(
               c(PCA1 = "PCA", CA1 = "CA",
                 DCA1_r = "DCA", NMDS1_r="NMDS"))) +
  coord_flip() +
  labs(y="Axis 1", x="Position (m.a.b)")
ggsave("smooth_TOV2_h.png", device="png", dpi=300, units = "cm", 
       width=20, height=10)

ggplot(data=facet_data[facet_data$locality=="ULV1+2",]) +
  geom_point(aes(y=axis_1, x=height_m, color=locality)) +
  geom_smooth(aes(y=axis_1, x=height_m), method = "loess", color="black") +
  facet_wrap(~ordination, scales="free", nrow = 1, 
             labeller = as_labeller(
              c(PCA1 = "PCA", CA1 = "CA",
               DCA1_r = "DCA", NMDS1_r="NMDS"))) +
  coord_flip() +
  labs(y="Axis 1", x="Position (m.a.b)")
ggsave("smooth_ULV_h.png", device="png", dpi=300, units = "cm", 
       width=20, height=10)


