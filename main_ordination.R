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

### importing data from CSV ###
T2_data <- read.csv("TOV2_raw_counts.csv", sep=",", dec=",", header = 1)
U1_data <- read.csv("ULV1_raw_counts.csv", sep=",", dec=",", header = 1)
U2_data <- read.csv("ULV2_raw_counts.csv", sep=",", dec=",", header = 1)

# cleaning up the minus sign formatting incompatibility, etc.
U2_data$height_m <- gsub("−", "-", U2_data$height_m)
U2_data$height_m <- gsub(",", ".", U2_data$height_m)
U2_data$height_m <- as.numeric(U2_data$height_m) + 70 #adding assumed missing meters to ULV2
U2_data$squares <- as.numeric(U2_data$squares)

# making a composite data frame
T2_data$locality <- "TOV2"
U1_data$locality <- "ULV1"
U2_data$locality <- "ULV2"

# composite of Ulvøya localities
U_data <- full_join(U1_data, U2_data)

# composite of all localities
all_data <- full_join(T2_data, U_data)

#NB! Replacing NAs with 0
T2_data[is.na(T2_data)] <- 0
U1_data[is.na(U1_data)] <- 0
U2_data[is.na(U2_data)] <- 0
all_data[is.na(all_data)] <- 0

#plot score dataframe for ordination
norm_counts <- all_data %>% select(calc_tubes:strom)

#normalizing plot counts to 40 subplots
for (i in 1:(length(norm_counts[1,]))) {
  for (n in 1:length(norm_counts$calc_tubes)) {
    norm_counts[n,i] <- (norm_counts[n,i]/all_data$squares[n])*40
  }
}

# species dataframe for species scores
norm_species <- data.frame(t(norm_counts))
 img

### making image imformation dataframe ###
# image paths
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

# relative scaling of images
img_size <-
  c(0.05, 0.05, 0.05, 0.025, 
    0.05, 0.05, 0.05, 0.05, 0.05, 
    0.05, 0.05, 0.035, 0.05, 
    0.05, 0.05, 0.05, 0.05, 
    0.05)

# adding image Markdown label paths
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

#adding image data to species dataframe
norm_species$img <- img
norm_species$img_size <- img_size
norm_species$img_md <- img_md


### calculating diversity metrics for (raw) site data ###
# species richness (R)
all_data$R <- specnumber(norm_counts)

# Shannon's diversity index (H)
all_data$H <- diversity(norm_counts)

# species evenness (E)
all_data$E <- all_data$H/log(all_data$R)
# manually setting E in samples with only 1 species to 0
all_data$E[is.na(all_data$E)] <- 0


### abundance scale weighting of species dataframe ###
# abundance scale = 16   
norm_counts_16 <- norm_counts
for (i in 1:(length(norm_counts_16[1,]))) {
  for (n in 1:length(norm_counts_16$calc_tubes)) {
    norm_counts_16[n,i] <- (norm_counts_16[n,i])^(log(16)/log(40))
  }
}

# abundance scale = 8
norm_counts_8 <- norm_counts
for (i in 1:(length(norm_counts_8[1,]))) {
  for (n in 1:length(norm_counts_8$calc_tubes)) {
    norm_counts_8[n,i] <- (norm_counts_8[n,i])^(log(8)/log(40))
  }
}

# abundance scale = 4   
norm_counts_4 <- norm_counts
for (i in 1:(length(norm_counts_4[1,]))) {
  for (n in 1:length(norm_counts_4$calc_tubes)) {
    norm_counts_4[n,i] <- (norm_counts_4[n,i])^(log(4)/log(40))
  }
}

# abundance scale = 2
norm_counts_2 <- norm_counts
for (i in 1:(length(norm_counts_2[1,]))) {
  for (n in 1:length(norm_counts_2$calc_tubes)) {
    norm_counts_2[n,i] <- (norm_counts_2[n,i])^(log(2)/log(40))
  }
}

# abundance scale = p/a     
norm_counts_0 <- norm_counts
for (i in 1:(length(norm_counts_0[1,]))) {
  for (n in 1:length(norm_counts_0$calc_tubes)) {
    if (norm_counts_0[n,i] > 0) {
      norm_counts_0[n,i] <- 1
    }
    else {
      norm_counts_0[n,i] <- 0
    }
  }
}

### PCA ordination ###
#NB! Centering and variance-standardization set to True
pca.r <- rda(norm_counts, scale=T)

#extracting site scores
all_data$PCA1<-scores(pca.r,display="sites")[,1]
all_data$PCA2<-scores(pca.r,display="sites")[,2]

#extracting species scores
norm_species$PCA1<-scores(pca.r,display="species")[,1]
norm_species$PCA2<-scores(pca.r,display="species")[,2]


### CA ordination ###
ca.r <- cca(norm_counts)

# extracting site scores
all_data$CA1<-scores(ca.r,display="sites", hill=T)[,1]
all_data$CA2<-scores(ca.r,display="sites", hill=T)[,2]

# extracting species scores
norm_species$CA1<-scores(ca.r,display="species", hill=T)[,1]
norm_species$CA2<-scores(ca.r,display="species", hill=T)[,2]


### DCA ordinations ###
# raw data
dca.r <- decorana(norm_counts)

# site scores
all_data["DCA1_r"]<-scores(dca.r,display="sites",origin=T)[,1]
all_data["DCA2_r"]<-scores(dca.r,display="sites",origin=T)[,2]

# species scores
norm_species["DCA1_r"]<-scores(dca.r,display="species",hill=T, scaling=0)[,1]
norm_species["DCA2_r"]<-scores(dca.r,display="species",hill=T, scaling=0)[,2]

# R=16
dca.16 <- decorana(norm_counts_16)
all_data["DCA1_16"]<-scores(dca.16,display="sites",origin=T)[,1]
all_data["DCA2_16"]<-scores(dca.16,display="sites",origin=T)[,2]
norm_species["DCA1_16"]<-scores(dca.16,display="species",hill=T, scaling=0)[,1]
norm_species["DCA2_16"]<-scores(dca.16,display="species",hill=T, scaling=0)[,2]

# R=8
dca.8 <- decorana(norm_counts_8)
all_data["DCA1_8"]<-scores(dca.8,display="sites",origin=T)[,1]
all_data["DCA2_8"]<-scores(dca.8,display="sites",origin=T)[,2]
norm_species["DCA1_8"]<-scores(dca.8,display="species",hill=T, scaling=0)[,1]
norm_species["DCA2_8"]<-scores(dca.8,display="species",hill=T, scaling=0)[,2]

# R=4 ###OBS: - added in axis 2 to fit with the gradient
dca.4 <- decorana(norm_counts_4)
all_data["DCA1_4"]<-scores(dca.4,display="sites",origin=T)[,1]
all_data["DCA2_4"]<--scores(dca.4,display="sites",origin=T)[,2]
norm_species["DCA1_4"]<-scores(dca.4,display="species",hill=T, scaling=0)[,1]
norm_species["DCA2_4"]<-scores(dca.4,display="species",hill=T, scaling=0)[,2]

# R=2 ###OBS: - added to axis 2 to fit with the gradient
dca.2 <- decorana(norm_counts_2)
all_data["DCA1_2"]<-scores(dca.2,display="sites",origin=T)[,1]
all_data["DCA2_2"]<--scores(dca.2,display="sites",origin=T)[,2]
norm_species["DCA1_2"]<-scores(dca.2,display="species",hill=T, scaling=0)[,1]
norm_species["DCA2_2"]<-scores(dca.2,display="species",hill=T, scaling=0)[,2]

# R=p/a
dca.pa <- decorana(norm_counts_0)
all_data["DCA1_pa"]<-scores(dca.pa,display="sites",origin=T)[,1]
all_data["DCA2_pa"]<-scores(dca.pa,display="sites",origin=T)[,2]
norm_species["DCA1_pa"]<-scores(dca.pa,display="species",hill=T, scaling=0)[,1]
norm_species["DCA2_pa"]<-scores(dca.pa,display="species",hill=T, scaling=0)[,2]


### NMDS ordinations ###
# R=raw
nmds.r <- metaMDS(norm_counts, k=2, trymax=300,
                     noshare=T, autotransform=F)
#site scores
all_data["NMDS1_r"]<-data.frame(nmds.r$points)$MDS1
all_data["NMDS2_r"]<-data.frame(nmds.r$points)$MDS2

# species scores
norm_species["NMDS1_r"]<-data.frame(nmds.r$species)$MDS1
norm_species["NMDS2_r"]<-data.frame(nmds.r$species)$MDS2

# R=16
nmds.16 <- metaMDS(norm_counts_16, k=2, trymax=300,
                  noshare=T, autotransform=F)
nmds.16$stress
all_data$NMDS1_16 <- data.frame(nmds.16$points)$MDS1
all_data$NMDS2_16 <- data.frame(nmds.16$points)$MDS2
norm_species["NMDS1_16"]<-data.frame(nmds.16$species)$MDS1
norm_species["NMDS2_16"]<-data.frame(nmds.16$species)$MDS2

# R=8
nmds.8 <- metaMDS(norm_counts_8, k=2, trymax=300,
                   noshare=T, autotransform=F)
nmds.8$stress
all_data$NMDS1_8 <- data.frame(nmds.8$points)$MDS1
all_data$NMDS2_8 <- data.frame(nmds.8$points)$MDS2
norm_species["NMDS1_8"]<-data.frame(nmds.8$species)$MDS1
norm_species["NMDS2_8"]<-data.frame(nmds.8$species)$MDS2

# R=4
nmds.4 <- metaMDS(norm_counts_4, k=2, trymax=300,
                  noshare=T, autotransform=F)
nmds.4$stress
all_data$NMDS1_4 <- data.frame(nmds.4$points)$MDS1
all_data$NMDS2_4 <- data.frame(nmds.4$points)$MDS2
norm_species["NMDS1_4"]<-data.frame(nmds.4$species)$MDS1
norm_species["NMDS2_4"]<-data.frame(nmds.4$species)$MDS2

# R=2
nmds.2 <- metaMDS(norm_counts_2, k=2, trymax=300,
                  noshare=T, autotransform=F)
nmds.2$stress
all_data$NMDS1_2 <- data.frame(nmds.2$points)$MDS1
all_data$NMDS2_2 <- data.frame(nmds.2$points)$MDS2
norm_species["NMDS1_2"]<-data.frame(nmds.2$species)$MDS1
norm_species["NMDS2_2"]<-data.frame(nmds.2$species)$MDS2

# R=pa
nmds.pa <- metaMDS(norm_counts_0, k=2, trymax=300,
                  noshare=T, autotransform=F)
nmds.pa$stress
all_data$NMDS1_pa <- data.frame(nmds.pa$points)$MDS1
all_data$NMDS2_pa <- data.frame(nmds.pa$points)$MDS2
norm_species["NMDS1_pa"]<-data.frame(nmds.pa$species)$MDS1
norm_species["NMDS2_pa"]<-data.frame(nmds.pa$species)$MDS2


### melting the NMDS weighted data ###
NMDS_site <- melt(all_data, 
                  measure.vars = c("NMDS1_r", "NMDS1_16", "NMDS1_8", 
                                   "NMDS1_4", "NMDS1_2", "NMDS1_pa"),
                  variable.name="weighting",
                  value.name="NMDS1")

NMDS2_site <- melt(all_data, 
                   measure.vars = c("NMDS2_r", "NMDS2_16", "NMDS2_8", 
                                    "NMDS2_4", "NMDS2_2", "NMDS2_pa"),
                   value.name="NMDS2")

NMDS_site$NMDS2 <- NMDS2_site$NMDS2
levels(NMDS_site$weighting) <- c("Raw data (R = 40)", "R = 16", "R = 8",
                                 "R = 4", "R = 2", "Presence-absence (R = 0)")

### melting the DCA weighted data ###
DCA_site <- melt(all_data, 
                 measure.vars = c("DCA1_r", "DCA1_16", "DCA1_8", 
                                  "DCA1_4", "DCA1_2", "DCA1_pa"),
                 variable.name="weighting",
                 value.name="DCA1")

DCA2_site <- melt(all_data, 
                  measure.vars = c("DCA2_r", "DCA2_16", "DCA2_8", 
                                   "DCA2_4", "DCA2_2", "DCA2_pa"),
                  value.name="DCA2")

DCA_site$DCA2 <- DCA2_site$DCA2
levels(DCA_site$weighting) <- c("Raw data (R = 40)", "R = 16", "R = 8",
                                "R = 4", "R = 2", "Presence-absence (R = 0)")


### melting the NMDS species data ###
NMDS_species <- norm_species
NMDS_species$taxa <- rownames(norm_species)
NMDS_species <- melt(NMDS_species, 
                  measure.vars = c("NMDS1_r", "NMDS1_16", "NMDS1_8", 
                                   "NMDS1_4", "NMDS1_2", "NMDS1_pa"),
                  variable.name="weighting",
                  value.name="NMDS1")

NMDS2_species <- melt(norm_species, 
                   measure.vars = c("NMDS2_r", "NMDS2_16", "NMDS2_8", 
                                    "NMDS2_4", "NMDS2_2", "NMDS2_pa"),
                   value.name="NMDS2")

NMDS_species$NMDS2 <- NMDS2_species$NMDS2
levels(NMDS_species$weighting) <- c("Raw data (R = 40)", "R = 16", "R = 8",
                                 "R = 4", "R = 2", "Presence-absence (R = 0)")

### melting the DCA species data ###
DCA_species <- norm_species
DCA_species$taxa <- rownames(norm_species)
DCA_species <- melt(DCA_species, 
                 measure.vars = c("DCA1_r", "DCA1_16", "DCA1_8", 
                                  "DCA1_4", "DCA1_2", "DCA1_pa"),
                 variable.name="weighting",
                 value.name="DCA1")

DCA2_species <- melt(norm_species, 
                  measure.vars = c("DCA2_r", "DCA2_16", "DCA2_8", 
                                   "DCA2_4", "DCA2_2", "DCA2_pa"),
                  value.name="DCA2")

DCA_species$DCA2 <- DCA2_species$DCA2
levels(DCA_species$weighting) <- c("Raw data (R = 40)", "R = 16", "R = 8",
                                "R = 4", "R = 2", "Presence-absence (R = 0)")

# saving all objects as Rdata file
save(all_data, T2_data, U1_data, U2_data, img,
     img_md, img_size, norm_counts, norm_species,
     norm_counts_16, norm_counts_8, norm_counts_4,
     norm_counts_2, norm_counts_0,
     pca.r, ca.r, nmds.r, dca.r,
     nmds.16, nmds.8, nmds.4, nmds.2, nmds.pa,
     dca.16, dca.8, dca.4, dca.2, dca.pa,
     DCA_site, NMDS_site, DCA_species, NMDS_species,
     file="ordination_objects.Rdata")