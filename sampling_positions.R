library(arm)
library(plyr)
#library(vegan)
library(graphics)
#library(data.frame)
library(ggplot2)
library(tidypaleo)
#install.packages("svglite")
library(svglite)

all_data <- read.csv("sample_overview - PMO_peels.csv", sep=",", dec=",", header = 1)

#TOV_2
T_data <- subset(all_data[all_data$locality=="TOV_2", ])
strat_pos <- c(levels(as.factor((T_data$Height))))
strat_pos <- na.omit(
        (as.numeric(
                sub(",", ".", strat_pos), fixed=F)))
T_pos <- data.frame(1,strat_pos)

#ULV_1
U1_data <- subset(all_data[all_data$locality=="ULV_1", ])
strat_pos <- c(levels(as.factor((U1_data$Height))))
strat_pos <- na.omit(
        (as.numeric(
                sub(",", ".", strat_pos), fixed=F)))
U1_pos <- data.frame(1,strat_pos)

#ULV_2
U2_data <- subset(all_data[all_data$locality=="ULV_2", ])
strat_pos <- c(levels(as.factor((U2_data$Height))))
strat_pos <- sub("−", "-", strat_pos)
strat_pos <- na.omit(
        (as.numeric(
                sub(",", ".", strat_pos), fixed=F)))
U2_pos <- data.frame(1,strat_pos)

#U1+U2 composite
comp_pos <- append(U1_pos$strat_pos, U2_pos$strat_pos+70)
U_comp_pos <- data.frame(2, comp_pos)

#plotting the data
plot(T_pos, type = 'o', pch = '__', ylab = 'Height in the sequence (m)', xlab="", ylim=c(0,40))
plot(U1_pos, type = 'o', pch = '__', ylab = 'Height in the sequence (m)', xlab="")
plot(U2_pos, type = 'o', pch = '__', ylab = 'Height in the sequence (m)', xlab="")
plot(U_comp_pos, type = 'o', pch = '__', ylab = 'Height in the sequence (m)', xlab="")

#plotting the data in ggplot
pos <- ggplot() +
        geom_point(aes(x=X1, y=strat_pos, fill="Sample"), data=T_pos, shape=3) +
        geom_line(aes(x=X1, y=strat_pos, color="Toverud farm"), data=T_pos) +
        geom_point(aes(x=X2, y=comp_pos), data=U_comp_pos, shape=3) +
        geom_line(aes(x=X2, y=comp_pos, color="Ulvøya (composite)"), data=U_comp_pos) +
        xlim(c(0,3)) + 
        labs(x="", y="Height above formation base (m)", color="Site", fill="Legend") +
        scale_colour_brewer(palette = "Set1") +
        theme_minimal() +
        theme(axis.text.x=element_blank(),
              axis.ticks.x=element_blank())
ggsave("plot_1.png", pos, units="cm", dpi=300, height = 15, width=8)        

#flipped
pos2 <- ggplot() +
        geom_point(aes(y=-X1, x=strat_pos), data=T_pos, shape=108, size=3.5) +
        geom_line(aes(y=-X1, x=strat_pos, color="Toverud farm"), data=T_pos, size=0.5) +
        geom_point(aes(y=-X2, x=comp_pos, fill="Sample"), data=U_comp_pos, shape=108, size=3.5) +
        geom_line(aes(y=-X2, x=comp_pos, color="Ulvøya (composite)"), data=U_comp_pos) +
        ylim(c(-3,0)) + 
        labs(y="", x="Height above formation base (m)", color="Site", fill="") +
        guides(colour = guide_legend(order = 1), 
                fill = guide_legend(order = 2)) +
        scale_colour_brewer(palette = "Set1") +
        theme_minimal() +
        theme(axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank(),)
pos2
ggsave("plot_2.svg", pos2, units="cm", dpi=300, height = 5, width=15)  

#measurements (table ?)
length(T_pos$strat_pos)
length(U1_pos$strat_pos)
length(U2_pos$strat_pos)

# slabs
length(levels(as.factor(T_data$slab))) - 1 #accounting for "missing slab" level
length(levels(as.factor(U1_data$slab)))
length(levels(as.factor(U2_data$slab)))

# peels
length(levels(as.factor(T_data$PMO))) - 1 #accounting for "missing slab" level
length(levels(as.factor(U1_data$PMO)))
length(levels(as.factor(U2_data$PMO)))

knitr::kable(T_pos, format = "markdown")

