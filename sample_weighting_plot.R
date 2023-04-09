library(tidyverse)
library(vegan)
library(graphics)
library(gridExtra)
library(ggimage)
library(reshape2)

abd <- data.frame(seq(0, 40, 0.2), seq(0, 40, 0.2))
colnames(abd) <- c("x", "raw")

abd$R16 <- (abd$x^(log(16)/log(40))/16)*40
abd$R8 <- (abd$x^(log(8)/log(40))/8)*40
abd$R4 <- (abd$x^(log(4)/log(40))/4)*40
abd$R2 <- (abd$x^(log(2)/log(40))/2)*40
abd$R0 <- rep(40, 40*5+1)
abd$R0[1] <- 0
abd <- melt(abd, id.vars="x")

ggplot(abd) +
  geom_line(aes(x=x, y=value, color=variable)) +
  labs(x="Original abundance", y="Weighted abundance, standardized range",
           color="Weighting") +
  theme_minimal()
  ggsave("abundance_scores.png", device="png", dpi=300, units = "cm", 
         width=13, height=10)

