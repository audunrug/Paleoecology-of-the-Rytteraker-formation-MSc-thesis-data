load("ordination_objects.Rdata")


### Relative core length of axes ###
i_max <- round(length(all_data$PMO)/10, 0)
l <- length(all_data$PMO)
all_data$pca1

outlier_table <- data.frame(Axis=c("PCA1", "CA1", "DCA1", "NMDS1"), 
                            RCL=0, 
                            Outlier_min=0, 
                            Outlier_max=0)

axis1 <- all_data %>% select(PCA1, CA1, DCA1_r, NMDS1_r)

#loop
for (n in 1:ncol(outlier_table)) {
  vect <- sort(axis1[,n], decreasing=T)
  total <- vect[1] - vect[l]
  min <- vect[1] - vect[(l-i_max)]
  for (i in 1:i_max) {
    core_l_pca <- pca_vect[i] - pca_vect[l-(i_max+i)]
    if (core_l_pca < pca_min){
      pca_min <- core_l_pca
      pca_rcl <- pca_min/total_pca
    }
    
    core_l_ca <- ca_vect[i] - ca_vect[l-(i_max+i)]
    if (core_l_ca < ca_min){
      ca_min <- core_l_ca
      ca_rcl <- ca_min/total_ca
    }
    
    core_l_dca <- dca_vect[i] - dca_vect[l-(i_max+i)]
    if (core_l_dca < dca_min){
      dca_min <- core_l_dca
      dca_rcl <- dca_min/total_dca
    }
    
    core_l_nmds <- nmds_vect[i] - nmds_vect[l-(i_max+i)]
    if (core_l_nmds < nmds_min){
      nmds_min <- core_l_nmds
      nmds_rcl <- nmds_min/total_nmds
    }
  
  }
}
#vectors
pca_vect <- sort(all_data$PCA1, decreasing=T)
pca_min <- pca_vect[1] - pca_vect[(l-i_max)]

ca_vect <- sort(all_data$CA1, decreasing=T)
ca_min <- ca_vect[1] - ca_vect[(l-i_max)]

dca_vect <- sort(all_data$DCA1_r, decreasing=T)
dca_min <- dca_vect[1] - dca_vect[(l-i_max)]

nmds_vect <- sort(all_data$NMDS1_r, decreasing=T)
nmds_min <- nmds_vect[1] - nmds_vect[(l-i_max)]

total_pca <- pca_vect[1]-pca_vect[l]
total_ca <- ca_vect[1]-ca_vect[l]
total_dca <- dca_vect[1]-dca_vect[l]
total_nmds <- nmds_vect[1]-nmds_vect[l]

#loop


#finding min and max samples
min_PCA1 <- all_data[order(-all_data$PCA1),]$PMO[1]
max_PCA1 <- all_data[order(-all_data$PCA1),]$PMO[l]

min_CA1 <- all_data[order(-all_data$CA1),]$PMO[1]
max_CA1 <- all_data[order(-all_data$CA1),]$PMO[l]

min_dCA1 <- all_data[order(-all_data$DCA1_r),]$PMO[1]
max_dCA1 <- all_data[order(-all_data$DCA1_r),]$PMO[l]

min_nmds1 <- all_data[order(-all_data$NMDS1_r),]$PMO[1]
max_nmds1 <- all_data[order(-all_data$NMDS1_r),]$PMO[l]


outlier_table$RCL <- round(c(pca_rcl, ca_rcl, 
                             dca_rcl, nmds_rcl), 2)

outlier_table$Outlier_min <- c(min_PCA1, min_CA1,
                               min_dCA1, min_nmds1)

outlier_table$Outlier_max <- c(max_PCA1, max_CA1,
                               max_dCA1, max_nmds1)

#inspecting table
outlier_table


### ab. weighted RCL ###
i_max <- round(length(all_data$PMO)/10, 0)
l <- length(all_data$PMO)

aw_outliers <- data.frame(DCA=rep(0,6),
                          NMDS=rep(0,6))

#vectors
dr_vect <- sort(all_data$PCA1, decreasing=T)
nr_vect <- sort(all_data$NMDS1_r, decreasing=T)

DCA <- all_data %>% select(DCA1_r, DCA1_16, DCA1_8,
                           DCA1_4, DCA1_2, DCA1_pa)
###DCA
for (n in 1:ncol(DCA)) {
  vect <- sort(DCA[,n], decreasing=T)
  total <- vect[1] - vect[l]
  min <- vect[1] - vect[(l-i_max)]
    for (i in 1:i_max) {
      core_l <- vect[i] - vect[l-(i_max+i)]
      if (core_l < min){
        min <- core_l
        rcl <- min/total
      }
    }
  aw_outliers$DCA[n] <- rcl
}

###NMDS
NMDS <- all_data %>% select(NMDS1_r, NMDS1_16, NMDS1_8,
                           NMDS1_4, NMDS1_2, NMDS1_pa)
#loop
for (n in 1:ncol(NMDS)) {
  vect <- sort(NMDS[,n], decreasing=T)
  total <- vect[1] - vect[l]
  min <- vect[1] - vect[(l-i_max)]
  for (i in 1:i_max) {
    core_l <- vect[i] - vect[l-(i_max+i)]
    if (core_l < min){
      min <- core_l
      rcl <- min/total
    }
  }
  aw_outliers$NMDS[n] <- rcl
}

#inspecting table
round(aw_outliers, 2)

