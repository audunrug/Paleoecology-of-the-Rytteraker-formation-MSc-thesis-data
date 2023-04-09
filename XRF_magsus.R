#install.packages("MIAmaxent")
library(dplyr)
library(vegan)
library(graphics)
library(ggplot2)
library(tidypaleo)
library(gridExtra)
library(moments)
library(reshape2)
library(MIAmaxent)                    
library("compare")         

### importing data ###
T2_xrf <- read.csv("toverud_xrf_summary.csv", sep=";", dec=",")
U1_magsus <- read.csv("ULV1_ms.csv", sep=";", dec=",")
U2_magsus <- read.csv("ULV2_ms.csv", sep=";", dec=",")
T2_magsus <- read.csv("TOV2_ms.csv", sep=";", dec=",")

#renaming and cleaning up
U1_magsus <- U1_magsus %>% rename("position_m" = ULV_1)
U2_magsus <- U2_magsus %>% rename("position_m" = ULV_2)
T2_magsus <- T2_magsus %>% rename("position_m" = Level..m., SI.10..6 = Mag..susc...SI.10..6.)
U1_magsus$locality <- "ULV1" #49 samples
U2_magsus$locality <- "ULV2" #61 samples
T2_magsus$locality <- "TOV2" #87 samples
U2_magsus$position_m <- U2_magsus$position_m + 70 #adding assumend missing meters to ULV2

# joining magsus data
# composite of Ulvøya localities
U_magsus <- full_join(U1_magsus, U2_magsus)

# composite of all localities
all_magsus <- full_join(T2_magsus, U_magsus)

# xrf cleanup
T2_xrf <- T2_xrf %>% select(SAMPLE, Ca, Si, Fe, Mn, Sr)
T2_xrf$"Mn/Sr" <- T2_xrf$Mn/T2_xrf$Sr
T2_xrf[,2:5] <- T2_xrf[,2:5]/1000

colnames(T2_xrf) <- c("Height_m", "Ca", "Si", "Fe","Mn", "Sr", "Mn/Sr")
T2_xrf_log <- T2_xrf %>% select(Height_m:Fe, "Mn/Sr")


### facet wrap ###
T2_xrf_melt <- melt(T2_xrf_log, id.vars="Height_m")
ggplot(T2_xrf_melt) +
  geom_point(aes(x=value, y=Height_m, color=variable)) +
  geom_lineh(aes(x=value, y=Height_m, color=variable)) +
  facet_wrap(~variable, ncol=4, scale="free_x") +
  labs(x="ppt / ratio", y="Stratigraphic position (m)",
       col="Element") +
  theme_minimal()
ggsave("XRF_l.png", device="png", dpi=300, units = "cm", 
       width=23, height=14)


#adding magnetic susceptibility readings to T2_xrf
#T2_xrf$Mag <- T2_magsus$Mag..susc...SI.10..6.[match(T2_xrf$RV, T2_magsus$Level..m.)]

T2_xrf <- data.frame(deriveVars(T2_xrf, transformtype = "L")$dvdata)

colnames(T2_xrf) <- c("Height_m", "Ca", "Si", "Fe", "Mn", "Sr", "Mn/Sr")

#"species" dataframe
T2_xrf_v <- T2_xrf %>% select(-Height_m, -"Mn/Sr") %>% t() %>% data.frame()

### environmental PCA ###
xrf_var <- T2_xrf %>% select(-Height_m, -"Mn/Sr")
pca.xrf <- rda(xrf_var, scale=T)
plot(pca.xrf)
summary(pca.xrf)
screeplot(pca.xrf)
T2_xrf$PCA1 <-scores(pca.xrf, display="sites")[,1]
T2_xrf$PCA2 <-scores(pca.xrf, display="sites")[,2]
T2_xrf_v$PCA1 <-scores(pca.xrf, display="species")[,1]
T2_xrf_v$PCA2 <-scores(pca.xrf, display="species")[,2]

ggplot(T2_xrf) +
  geom_point(aes(x=PCA1, y=Height_m), color="Coral") +
  geom_lineh(aes(x=PCA1, y=Height_m), color="Coral") +
  labs(x="PCA1 (65.1%)", y="Stratigraphic position (m)") +
  theme_minimal()
ggsave("XRF_pca1H.png", device="png", dpi=300, units = "cm", 
       width=6, height=14)

ggplot(T2_xrf) +
  geom_point(aes(x=PCA1, y=PCA2), color="Coral") +
  geom_segment(data=T2_xrf_v, 
               aes(x = 0, xend = PCA1, y = 0, yend = PCA2),
               arrow = arrow(length = unit(0.3, "cm")), 
               size=.4, colour = "black") +
  geom_text(data=T2_xrf_v, aes(x = PCA1+0.15*(PCA1/(sqrt((PCA1)^2 + (PCA2)^2))), 
                               y = PCA2+0.15*(PCA2/(sqrt((PCA1)^2 + (PCA2)^2))),
            label=rownames(T2_xrf_v)), color="black") +
  labs(x="PCA1 (65.1%)", y="PCA2 (21.6%)") +
  theme_minimal()
ggsave("XRF_pca.png", device="png", dpi=300, units = "cm", 
       width=14, height=14)

x <- na.omit(match(T2_magsus$Level..m., T2_xrf$RV))


### magsus plotting ###
all_magsus$locality2 <- all_magsus$locality %>% 
  fct_collapse("ULV1+2" = c("ULV1","ULV2"))

ggplot(all_magsus) +
  geom_point(aes(y=position_m, x=SI.10..6, color=locality2)) +
#  geom_lineh(aes(y=position_m, x=SI.10..6, color=locality2, group=locality)) +
  facet_wrap(~locality2, scale="free_y") +
  labs(x=10^-6~SI, y="Stratigraphic position (m)", col="Locality") +
  theme_minimal()
ggsave("magsus.png", device="png", dpi=300, units = "cm", 
       width=15, height=15)


### applying magsus readings to  ###
x <- na.omit(match(T2_magsus$Level..m., T2_xrf$RV))
all_data$height_m[all_data$locality=="TOV2"]
T2_magsus$position_m

#matching stratigrapgically closest values from magsus dataset to peels
ms_list <- c()
for (i in 1:length(all_data$height_m)) {
  loc <- all_data$locality[i]
  msheights <- all_magsus$position_m[all_magsus$locality==loc]
  msvals <- all_magsus$SI.10..6[all_magsus$locality==loc]

  pos <- which.min(abs(msheights - all_data$height_m[i]))
  ms_list<- append(ms_list, msvals[pos])
}

# envfit tests
envfit(nmds.r, ms_list)
plot(envfit(nmds.r, ms_list))
plot(envfit(nmds.r, all_vars))

envfit(dca.r, ms_list)
plot(dca.r)
plot(envfit(dca.r, ms_list))

envfit(ca.r, ms_list)
plot(ca.r)
plot(envfit(ca.r, ms_list))

envfit(pca.r, ms_list)
plot(pca.r)
plot(envfit(pca.r, ms_list))

#correlation tests
cor(all_vars, ms_list)
cor.test(all_data$NMDS1_r, ms_list)
cor.test(all_data$DCA1_r, ms_list)
cor.test(all_data$CA1_r, ms_list)
cor.test(all_data$PCA1, ms_list)

