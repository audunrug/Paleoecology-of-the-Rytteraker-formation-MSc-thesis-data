### correlation between ordinations ###
#evenness etc.
all_data$H <- diversity(norm_counts)
all_data$R <- specnumber(norm_counts)
all_data$E <- all_data$H/log(specnumber(norm_counts))

mean(all_data$R)
sd(all_data$R)

ggplot(data=all_data)
all_data$E[is.na(all_data$E)] <- 0

cont_env <- all_data %>% select(bioturb_fraction, H, E, R)

# extracting ordination axes
axis1 <- all_data %>% select(CA1, PCA1, DCA1_r, NMDS1_r)
axis2 <- all_data %>% select(CA2, PCA2, DCA2_r, NMDS2_r)

round(cor(axis1, method="k",), 2)
round(cor(axis2, method="k"), 2)

cor.test(all_data$CA1, all_data$PCA1 , method="k")
cor.test(all_data$PCA1, all_data$DCA1_r, method="k")
cor.test(all_data$PCA1, all_data$NMDS1_r, method="k")
cor.test(all_data$CA1, all_data$DCA1_r , method="k")
cor.test(all_data$CA1, all_data$NMDS1_r, method="k")
cor.test(all_data$DCA1_r, all_data$NMDS1_r, method="k")

cor.test(all_data$CA2, all_data$PCA2 , method="k")
cor.test(all_data$PCA2, all_data$DCA2_r, method="k")
cor.test(all_data$PCA2, all_data$NMDS2_r, method="k")
cor.test(all_data$CA2, all_data$DCA2_r , method="k")
cor.test(all_data$CA2, all_data$NMDS2_r, method="k")
cor.test(all_data$DCA2_r, all_data$NMDS2_r, method="k")

protest(pca.r, ca.r)
protest(pca.r, dca.r)
protest(pca.r, nmds.r)
protest(ca.r, dca.r)
protest(ca.r, nmds.r)
protest(dca.r, nmds.r)

envfit(pca.r, cont_env)
plot(pca.r)
plot(envfit(pca.r, cont_env))
ordisurf(pca.r, cont_env$R)

envfit(nmds.full, cont_env)
plot(nmds.full)
ordisurf(nmds.full, cont_env$H)

envfit(ca.r, cont_env)
plot(ca.r)
plot(envfit(ca.r, cont_env))
ordisurf(ca.r, cont_env$H)

envfit(dca.raw, cont_env)
plot(dca.raw)
plot(envfit(dca.raw, cont_env))
ordisurf(dca.raw, cont_env$H)


# abundance weighted data
axis1_w <- all_data %>% select(DCA1_8, DCA1_2, NMDS1_8, NMDS1_2)
round(cor(axis1, method="k",), 2)



