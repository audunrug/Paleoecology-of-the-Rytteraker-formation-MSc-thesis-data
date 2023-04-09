# loading dependencies (NB: must be installed beforehand)
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

#importing data
T2_data <- read.csv("TOV2_raw_counts.csv", sep=",", dec=",", header = 1)
U1_data <- read.csv("ULV1_raw-counts.csv", sep=",", dec=",", header = 1)
U2_data <- read.csv("ULV2_raw-counts.csv", sep=",", dec=",", header = 1)
img <- c(
  "fossil_symbols/calc_tube.png", "fossil_symbols/bryo_small.png",
  "fossil_symbols/bryo_large.png", "fossil_symbols/bryo_bifoliate.png",
  "fossil_symbols/tentaculitid.png", "fossil_symbols/tabcor_favhel.png",
  "fossil_symbols/tabcor_hal.png", "fossil_symbols/tabcor_indet.png", 
  "fossil_symbols/algae.png", "fossil_symbols/rugcor.png", 
  "fossil_symbols/trilobite.png", "fossil_symbols/gastropod.png", 
  "fossil_symbols/ostracod.png", "fossil_symbols/brach_straight_thick.png", 
  "fossil_symbols/brach_straight_thin.png", "fossil_symbols/brach_wavy.png",
  "fossil_symbols/crinoid.png", "fossil_symbols/strom.png"
)

# adding image label paths
img_md <- img
a <- image_read(img_md[i])
for (i in 1:length(img_md)) {
  dims <- image_info(image_read(img_md[i]))
  if (dims$width[1] > dims$height[1]) {
    d = "width"
  } else {
    d = "height"
  }
  img_md[i] <- paste("<img src=", img_md[i], d,"='20'/>")
}
img_md

# cleaning up the minus sign formatting incompatibility ++
U2_data$height_m <- gsub("−", "-", U2_data$height_m)
U2_data$height_m <- gsub(",", ".", U2_data$height_m)
U2_data$height_m <- as.numeric(U2_data$height_m) + 70 #adding assumed missing meters to ULV2
U2_data$squares <- as.numeric(U2_data$squares)

# making a composite data frame
T2_data$locality <- "TOV2"
U1_data$locality <- "ULV1"
U2_data$locality <- "ULV2"

U_data <- full_join(U1_data, U2_data)
all_data <- full_join(T2_data, U_data)

#NB! Replacing NAs with 0
T2_data[is.na(T2_data)] <- 0
U1_data[is.na(U1_data)] <- 0
U2_data[is.na(U2_data)] <- 0
all_data[is.na(all_data)] <- 0

#remove species with fewer than 10 occurences
colSums(norm_counts!=0)
norm_countsF <- norm_counts %>% select(-calc_tubes,
                                       -tabcor_favositid,
                                       -algae, -strom)

#PCA
pca.f <- rda(norm_countsF, scale=F)
summary(pca.f)
plot(pca.f)
all_data$pcaf1 <- scores(pca.f,display="sites")[,1]

#CA
ca.f <- cca(norm_countsF, scale=T)
summary(ca.f)
plot(ca.f)
all_data$caf1 <- scores(ca.f,display="sites")[,1]

#DCA
dca.f <- decorana(norm_countsF)
summary(dca.f)
plot(dca.f)
all_data$dcaf1<-scores(dca.f,display="sites")[,1]

#NMDS
nmds.f <- metaMDS(norm_countsF, trymax=300, noshare=T, k=2,
                  autotransform=F)
nmds.f$stress
all_data$nmdsf1 <- scores(nmds.f,display="sites")[,1]

#axis 1 and height
ggplot(data=all_data) +
  geom_point(aes(x=pcaf1, y=height_m, color=locality))+
  facet_wrap(~locality, scale="free_y")

ggplot(data=all_data) +
  geom_point(aes(x=caf1, y=height_m, color=locality))+
  facet_wrap(~locality, scale="free_y")

ggplot(data=all_data) +
  geom_point(aes(x=dcaf1, y=height_m, color=locality))+
  facet_wrap(~locality, scale="free_y")

ggplot(data=all_data) +
  geom_point(aes(x=-nmdsf1, y=height_m, color=locality))+
  facet_wrap(~locality, scale="free_y")
