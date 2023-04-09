library(tidyverse)
library(vegan)
library(graphics)
library(tidypaleo)
library(reshape2)
library(ggimage)

### importing data ###
data_raw <- read.csv("Mork_bedding-planes_data.csv", sep=",", dec=",", header = 1)

# image paths 
img <- c(
  "mork_symbols/cephalopod.png", "mork_symbols/strom.png",
  "mork_symbols/tabcor_fav.png", "mork_symbols/tabcor_hal.png", 
  "mork_symbols/tabcor_hel.png", "mork_symbols/tabcor_encrusting.png", 
  "mork_symbols/tabcor_syr.png", "mork_symbols/tabcor_indet.png", 
  "mork_symbols/rugcor.png", "mork_symbols/bryozoan.png", 
  "mork_symbols/brach_unid.png", "mork_symbols/pentamerus.png", 
  "mork_symbols/atrypid.png", "mork_symbols/gastropod.png",
  "mork_symbols/trilobite.png", "mork_symbols/crinoid.png"
)

img

#image sizes
img_size1 <-
  c(0.05, 0.05, 0.05, 0.05, 
    0.05, 0.05, 0.05, 0.05, 
    0.05, 0.05, 0.05, 0.05, 
    0.05, 0.035, 0.05, 0.05)


# Dataset where + = 1
data_1 <- data_raw
data_1[data_1 == "+"] <- 1
data_1[,3:26] <- sapply(data_1[,3:26], as.numeric)
data_1[,2] <- sapply(data_1[,2], as.factor)



# Dataset where + = 0
data_0 <- data_raw
data_0[data_0 == "+"] <- 0
data_0[,3:26] <- sapply(data_0[,3:26], as.numeric)
data_0[,2] <- sapply(data_1[,2], as.factor)


### ordination subsets ###
#morphotaxa
data_1_taxa <- data_1 %>% select(cephalopods:crinoids)
data_0_taxa <- data_0 %>% select(cephalopods:crinoids)

# species richness (R)
data_1$R <- specnumber(data_1_taxa)
data_0$R <- specnumber(data_0_taxa)

# Shannon's diversity index (H)
data_1$H <- diversity(data_1_taxa)
data_0$H <- diversity(data_0_taxa)

# species evenness (E)
data_1$E <- data_1$H/log(data_1$R)
data_0$E <- data_0$H/log(data_0$R)
# manually setting E in samples with only 1 species to 0
data_1$E[is.na(data_1$E)] <- 0
data_0$E[is.na(data_0$E)] <- 0

#species dataframe
data_1_s <- (data.frame(t(data_1_taxa)))
data_1_s$img <- img
data_1_s$img_size <- img_size
rownames(data_0_s)
data_0_s <- (data.frame(t(data_0_taxa)))
data_0_s$img <- img
data_0_s$img_size <- img_size

#lithology
data_1_lit <- data_1 %>% select(lime_finegrained:shale)
data_0_lit <- data_0 %>% select(lime_finegrained:shale)

#trace fossils
data_1_trace <- data_1 %>% select(chondrites:tracefossils_other)
data_0_trace <- data_0 %>% select(chondrites:tracefossils_other)

ggplot() +
  geom_point(aes(y=data_1$H, x=data_1$Level..7b., color="1")) +
  geom_point(aes(y=data_0$H, x=data_0$Level..7b., color="0"))

ggplot() +
  geom_point(aes(y=data_0$H, x=data_0$Level..7b.))

### DCA of +1 fossil counts ###
dca.1 <- decorana(data_1_taxa)
summary(dca.1)

data_1["DCA1"]<-scores(dca.1,display="sites",origin=T)[,1]
data_1["DCA2"]<-scores(dca.1,display="sites",origin=T)[,2]
data_1_s["DCA1"]<-scores(dca.1,display="species",origin=T)[,1]
data_1_s["DCA2"]<-scores(dca.1,display="species",origin=T)[,2]


### DCA of +0 fossil counts ###
dca.0 <- decorana(data_0_taxa)
summary(dca.0)

data_0["DCA1"]<-scores(dca.0,display="sites",origin=T)[,1]
data_0["DCA2"]<-scores(dca.0,display="sites",origin=T)[,2]
data_0_s["DCA1"]<-scores(dca.0,display="species",origin=T)[,1]
data_0_s["DCA2"]<-scores(dca.0,display="species",origin=T)[,2]


### NMDS of +1 fossil counts ###
nmds.1 <- metaMDS(data_1_taxa, noshare=T,
                  autotransform=F, trymax=300, k=2)
nmds.1$stress
data_1["NMDS1"]<-scores(nmds.1,display="sites",origin=T)[,1]
data_1["NMDS2"]<-scores(nmds.1,display="sites",origin=T)[,2]
data_1_s["NMDS1"]<-scores(nmds.1,display="species",origin=T)[,1]
data_1_s["NMDS2"]<-scores(nmds.1,display="species",origin=T)[,2]

### NMDS of +0 fossil counts ###
nmds.0 <- metaMDS(data_0_taxa, noshare=T,
                  autotransform=F, trymax=300, k=2)
nmds.0$stress
plot(nmds.0)
data_0["NMDS1"]<-scores(nmds.0,display="sites",origin=T)[,1]
data_0["NMDS2"]<-scores(nmds.0,display="sites",origin=T)[,2]
data_0_s["NMDS1"]<-scores(nmds.0,display="species",origin=T)[,1]
data_0_s["NMDS2"]<-scores(nmds.0,display="species",origin=T)[,2]

### plotting plot scores ###
ggplot(data=data_1, aes(x=DCA1, y=DCA2)) +
  geom_point(aes(color=Level..7b.)) +
  labs(col="Bedding plane") +
  theme(legend.position = "bottom")
#  scale_color_brewer(palette="Set1")
ggsave("DCA1.png", device="png", width=10, height=7.5, units="cm")

ggplot(data=data_0, aes(x=DCA1, y=DCA2)) +
  geom_point(aes(color=Level..7b.)) 
#  scale_color_brewer(palette="Set1")
ggsave("DCA0.png", device="png", width=10, height=7.5, units="cm")

ggplot(data=data_1, aes(x=-NMDS1, y=NMDS2)) +
  geom_point(aes(color=Level..7b.)) 
#  scale_color_brewer(palette="Set1")
ggsave("NMDS1.png", device="png", width=10, height=7.5, units="cm")

ggplot(data=data_0, aes(x=-NMDS1, y=NMDS2)) +
  geom_point(aes(color=Level..7b.)) 
#  scale_color_brewer(palette="Set1")
ggsave("NMDS0.png", device="png", width=10, height=7.5, units="cm")


### plotting species scores ###
ggplot(data=data_1_s, aes(x=DCA1, y=DCA2)) +
  geom_image(image=img, size=img_size*1.6)
ggsave("DCA1s.png", device="png", width=9, height=7.5, units="cm")

ggplot(data=data_0_s, aes(x=DCA1, y=DCA2)) +
  geom_image(image=img, size=img_size*1.6)
#  scale_color_brewer(palette="Set1")
ggsave("DCA0s.png", device="png", width=9, height=7.5, units="cm")

ggplot(data=data_1_s, aes(x=-NMDS1, y=NMDS2)) +
  geom_image(image=img, size=img_size*1.6)
ggsave("NMDS1s.png", device="png", width=9, height=7.5, units="cm")

ggplot(data=data_0_s, aes(x=-NMDS1, y=NMDS2)) +
  geom_image(image=img, size=img_size*1.6)
ggsave("NMDS0s.png", device="png", width=9, height=7.5, units="cm")


# correlation between axes
cor.test(data_1$DCA1,data_0$DCA1,method="k")
cor.test(data_1$DCA2,data_0$DCA2,method="k")

cor.test(data_1$NMDS1,data_0$NMDS1,method="k")
cor.test(data_1$NMDS2,data_0$NMDS2,method="k")

cor.test(data_1$DCA1,data_1$NMDS1,method="k")
cor.test(data_1$DCA2,data_1$NMDS2,method="k")

cor.test(data_0$DCA1,data_0$NMDS1,method="k")
cor.test(data_0$DCA2,data_0$NMDS2,method="k")

cor.test(data_1$DCA1,data_0$NMDS1,method="k")
cor.test(data_1$DCA2,data_0$NMDS2,method="k")

cor.test(data_1$NMDS1,data_0$DCA1,method="k")
cor.test(data_1$NMDS2,data_0$DCA2,method="k")

protest(dca.1, dca.0)
protest(nmds.0, nmds.1)
protest(dca.1, nmds.1)
protest(dca.0, nmds.0)
protest(dca.1, nmds.0)
protest(dca.0, nmds.1)

#fitting species
envfit(dca.1, data_1_taxa)
envfit(dca.0, data_0_taxa)

envfit(nmds.1, data_1_taxa)
envfit(nmds.0, data_0_taxa)

# fitting environmental variables
envfit(dca.1, data_1_lit/100)
envfit(dca.0, data_0_lit/100)

envfit(nmds.1, data_1_lit/100)
envfit(nmds.0, data_0_lit/100)

#fitting trace fossils
envfit(dca.1, data_1_trace)
envfit(dca.0, data_0_trace)
envfit(nmds.1, data_1_trace)
envfit(nmds.0, data_0_trace)


### environmental correlation ###
cor.test(data_0$DCA1, data_1_lit$lime_finegrained,method="k")


### melting data
data_0$w <- "+ = 0"
data_1$w <- "+ = 1"

data_full <- full_join(data_0, data_1)
data_full <- melt(data_full, measure.vars=c("DCA1", "NMDS1"))

ggplot(data=data_full, aes(x=DCA1, y=DCA2))


cor(data_1_trace)
symnum(cor(data_0_lit, method = "k"))


ggplot() +
  geom_point(aes(x=data_1$DCA2, y=data_0$DCA3))

cor.test()
