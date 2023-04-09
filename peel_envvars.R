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
library(moments)

### loading ordination results and calculations ###
load("ordination_objects.Rdata")

### environmental/diversity variable matrix ###
# all litological variables measured
lit_env <- all_data %>% select(bioturb_fraction:grainstone)

# all "true" continous variables
cont_env <- all_data %>% select(bioturb_fraction, R, H, E)

#all litological and faunistic variabels
all_vars <- all_data %>% select(bioturb_fraction:grainstone, 
                                R, H, E)

# kendall's between sample correlation
cor.table <- round(cor(lit_env, cont_env, method="k"), 3)

# significance tests
cor.test(cont_env$bioturb_fraction, lit_env$mudstone, method="k")
cor.test(cont_env$bioturb_fraction, lit_env$wackestone, method="k")
cor.test(cont_env$bioturb_fraction, lit_env$packstone, method="k")
cor.test(cont_env$bioturb_fraction, lit_env$grainstone, method="k")

#R
cor.test(cont_env$R, lit_env$bioturb_fraction, method="k")
cor.test(cont_env$R, lit_env$mudstone, method="k")
cor.test(cont_env$R, lit_env$wackestone, method="k")
cor.test(cont_env$R, lit_env$packstone, method="k")
cor.test(cont_env$R, lit_env$grainstone, method="k")

#H
cor.test(cont_env$H, lit_env$bioturb_fraction, method="k")
cor.test(cont_env$H, lit_env$mudstone, method="k")
cor.test(cont_env$H, lit_env$wackestone, method="k")
cor.test(cont_env$H, lit_env$packstone, method="k")
cor.test(cont_env$H, lit_env$grainstone, method="k")

#E
cor.test(cont_env$E, lit_env$bioturb_fraction, method="k")
cor.test(cont_env$E, lit_env$mudstone, method="k")
cor.test(cont_env$E, lit_env$wackestone, method="k")
cor.test(cont_env$E, lit_env$packstone, method="k")
cor.test(cont_env$E, lit_env$grainstone, method="k")



### fitting of envrionmental variables to the ordinations ###
#T-cor
PCA1_t <- round(cor(all_vars, all_data$PCA1, method="k"), 2)[,1]
CA1_t <- round(cor(all_vars, all_data$CA1, method="k"), 2)[,1]
DCA1_t <- round(cor(all_vars, all_data$DCA1_r, method="k"), 2)[,1]
NMDS1_t <- round(cor(all_vars, all_data$NMDS1_r, method="k"), 2)[,1]

axis1_varcor <- data.frame(PCA1_t, CA1_t, DCA1_t, NMDS1_t)
axis1_varcor

# p-values
cor.test(all_vars$bioturb_fraction, all_data$PCA1, method="k")
cor.test(all_vars$mudstone, all_data$PCA1, method="k")
cor.test(all_vars$wackestone, all_data$PCA1, method="k")
cor.test(all_vars$packstone, all_data$PCA1, method="k")
cor.test(all_vars$grainstone, all_data$PCA1, method="k")
cor.test(all_vars$R, all_data$PCA1, method="k")
cor.test(all_vars$H, all_data$PCA1, method="k")
cor.test(all_vars$E, all_data$PCA1, method="k")

cor.test(all_vars$bioturb_fraction, all_data$CA1, method="k")
cor.test(all_vars$mudstone, all_data$CA1, method="k")
cor.test(all_vars$wackestone, all_data$CA1, method="k")
cor.test(all_vars$packstone, all_data$CA1, method="k")
cor.test(all_vars$grainstone, all_data$CA1, method="k")
cor.test(all_vars$R, all_data$CA1, method="k")
cor.test(all_vars$H, all_data$CA1, method="k")
cor.test(all_vars$E, all_data$CA1, method="k")

cor.test(all_vars$bioturb_fraction, all_data$DCA1_r, method="k")
cor.test(all_vars$mudstone, all_data$DCA1_r, method="k")
cor.test(all_vars$wackestone, all_data$DCA1_r, method="k")
cor.test(all_vars$packstone, all_data$DCA1_r, method="k")
cor.test(all_vars$grainstone, all_data$DCA1_r, method="k")
cor.test(all_vars$R, all_data$DCA1_r, method="k")
cor.test(all_vars$H, all_data$DCA1_r, method="k")
cor.test(all_vars$E, all_data$DCA1_r, method="k")

cor.test(all_vars$bioturb_fraction, all_data$NMDS1_r, method="k")
cor.test(all_vars$mudstone, all_data$NMDS1_r, method="k")
cor.test(all_vars$wackestone, all_data$NMDS1_r, method="k")
cor.test(all_vars$packstone, all_data$NMDS1_r, method="k")
cor.test(all_vars$grainstone, all_data$NMDS1_r, method="k")
cor.test(all_vars$R, all_data$NMDS1_r, method="k")
cor.test(all_vars$H, all_data$NMDS1_r, method="k")
cor.test(all_vars$E, all_data$NMDS1_r, method="k")

#envfit
pca.fit <- envfit(pca.r, all_vars)
ca.fit <- envfit(ca.r, all_vars)
dca.fit <- envfit(dca.r, all_vars)
nmds.fit <- envfit(nmds.r, all_vars)

pca.fit
ca.fit
dca.fit
nmds.fit

png(file="nmds_fit.png")
plot(nmds.r, display=c("sites"))
plot(nmds.fit, p.max=0.001)
dev.off()

png(file="dca_fit.png")
plot(dca.r, display=c("sites"))
plot(dca.fit, p.max=0.001)
dev.off()

png(file="ca_fit.png")
plot(ca.r, display=c("sites"), xlim=c(5,-5), ylim=c(4, -3.5))
plot(ca.fit, p.max=0.001)
dev.off()

png(file="pca_fit.png")
plot(pca.r, display=c("sites"), ylim=c(2.5, -1.2), xlim=c(-2.5, 2.5))
plot(pca.fit, p.max=0.001)
dev.off()

### abundance weighted ordinations ###
all_vars <- data.frame(scale(all_vars))
dca8.fit <- envfit(dca.8, all_vars)
dca2.fit <- envfit(dca.2, all_vars)
nmds8.fit <- envfit(nmds.8, all_vars)
nmds2.fit <- envfit(nmds.2, all_vars)

DCA1_8t <- round(cor(all_vars, all_data$DCA1_8, method="k"), 2)[,1]
DCA1_2t <- round(cor(all_vars, all_data$DCA1_2, method="k"), 2)[,1]
NMDS1_8t <- round(cor(all_vars, all_data$NMDS1_8, method="k"), 2)[,1]
NMDS1_2t <- round(cor(all_vars, all_data$NMDS1_2, method="k"), 2)[,1]

axis1_varcor_w <- data.frame(DCA1_8t, DCA1_2t, NMDS1_8t, NMDS1_2t)

#p-values
cor.test(all_vars$bioturb_fraction, all_data$DCA1_8, method="k")
cor.test(all_vars$mudstone, all_data$DCA1_8, method="k")
cor.test(all_vars$wackestone, all_data$DCA1_8, method="k")
cor.test(all_vars$packstone, all_data$DCA1_8, method="k")
cor.test(all_vars$grainstone, all_data$DCA1_8, method="k")
cor.test(all_vars$R, all_data$DCA1_8, method="k")
cor.test(all_vars$H, all_data$DCA1_8, method="k")
cor.test(all_vars$E, all_data$DCA1_8, method="k")

cor.test(all_vars$bioturb_fraction, all_data$DCA1_2, method="k")
cor.test(all_vars$mudstone, all_data$DCA1_2, method="k")
cor.test(all_vars$wackestone, all_data$DCA1_2, method="k")
cor.test(all_vars$packstone, all_data$DCA1_2, method="k")
cor.test(all_vars$grainstone, all_data$DCA1_2, method="k")
cor.test(all_vars$R, all_data$DCA1_2, method="k")
cor.test(all_vars$H, all_data$DCA1_2, method="k")
cor.test(all_vars$E, all_data$DCA1_2, method="k")

cor.test(all_vars$bioturb_fraction, all_data$NMDS1_8, method="k")
cor.test(all_vars$mudstone, all_data$NMDS1_8, method="k")
cor.test(all_vars$wackestone, all_data$NMDS1_8, method="k")
cor.test(all_vars$packstone, all_data$NMDS1_8, method="k")
cor.test(all_vars$grainstone, all_data$NMDS1_8, method="k")
cor.test(all_vars$R, all_data$NMDS1_8, method="k")
cor.test(all_vars$H, all_data$NMDS1_8, method="k")
cor.test(all_vars$E, all_data$NMDS1_8, method="k")

cor.test(all_vars$bioturb_fraction, all_data$NMDS1_2, method="k")
cor.test(all_vars$mudstone, all_data$NMDS1_2, method="k")
cor.test(all_vars$wackestone, all_data$NMDS1_2, method="k")
cor.test(all_vars$packstone, all_data$NMDS1_2, method="k")
cor.test(all_vars$grainstone, all_data$NMDS1_2, method="k")
cor.test(all_vars$R, all_data$NMDS1_2, method="k")
cor.test(all_vars$H, all_data$NMDS1_2, method="k")
cor.test(all_vars$E, all_data$NMDS1_2, method="k")

# plots
png(file="nmds8_fit.png")
plot(nmds.8, display=c("sites"))
plot(nmds8.fit, p.max=0.001)
dev.off()

png(file="nmds2_fit.png")
plot(nmds.2, display=c("sites"))
plot(nmds2.fit, p.max=0.001)
dev.off()

png(file="dca8_fit.png")
plot(dca.8, display=c("sites"), ylim=c(0.5, -1))
plot(dca8.fit, p.max=0.001, ylim=c(0.5, -1.3))
dev.off()

png(file="dca2_fit.png")
plot(dca.2, display=c("sites"))
plot(dca2.fit, p.max=0.001)
dev.off()


# GAM isolines (not included in thesis)
ordisurf(nmds.r, cont_env$R)
ordisurf(pca.r, cont_env$R)
ordisurf(ca.r, cont_env$R)
ordisurf(dca.r, cont_env$R)

ordisurf(nmds.r, cont_env$H)
ordisurf(pca.r, cont_env$H)
ordisurf(ca.r, cont_env$H)
ordisurf(dca.r, cont_env$H)

ordisurf(nmds.r, cont_env$bioturb_fraction)
ordisurf(pca.r, cont_env$bioturb_fraction)
ordisurf(ca.r, cont_env$bioturb_fraction)
ordisurf(dca.r, cont_env$bioturb_fraction)