#install.packages("arm")
#library(arm)
library(tidyverse)
library(vegan)
library(graphics)
#library(data.frame)
library(ggplot2)
#library(tidypaleo)

raw_data <- read.csv("resultater - shell_measurements.csv", sep=",", header = 1, dec=",")

raw_data <- raw_data %>% select(-notes, -notes.1)
raw_data <- na.omit(raw_data)

#new data column forfraction of convex shells
raw_data$convex_frac <- raw_data$convex/raw_data$total.shells

# binomial test for the null hypothesis that the shells configurations have equal probability 
test.total <- binom.test(sum(raw_data$convex), sum(raw_data$total.shells), 0.5)
test.total

# ???
success <- 100:200
plot(success,dbinom(success,size=300,prob=.4366667),
     type='l',
     main='Binomial Distribution (n=300, p=0.436)',
     ylab='Probability',
     xlab ='# Successes',
     lwd=2)
points(success,dbinom(success,size=300,prob=.5),
     type='l',
     lwd=2)


#make df for the bar chart
fraction <- c(sum(raw_data$convex)/300, sum(raw_data$concave)/300)
orientation <- c("convex", "concave")
shell_fract <- data.frame(orientation, fraction)

# plot of the binomial test
ggplot() +
  geom_bar(data=shell_fract, aes(fill=orientation, y=fraction, x="all shells"), 
           position="stack", stat="identity") +
  theme_minimal()

#(convex, concave) ~ height
testmod <- glm (cbind(convex, concave) ~ height, family = binomial, data=raw_data)
#plot(testmod)
summary(testmod, dispersion=1)

#(convex, concave) ~ total.shells
testmod2 <- glm (cbind(convex, concave) ~ total.shells, family = binomial, data=raw_data)
summary(testmod2)
#plot(testmod2)

#(convex, concave) ~ shells.cm2
testmod3 <- glm (cbind(convex, concave) ~ shells.cm2, family = binomial, data=raw_data)
summary(testmod3)
#plot(testmod3)

testmod4 <- glm (total.shells ~ height, family = quasipoisson, data=raw_data)
summary(testmod4)

testmod5 <- lm (shells.cm2 ~ height,  data=raw_data)
summary(testmod5)

plot(testmod5)

plot(raw_data$shells.cm2, raw_data$convex_frac)
plot(raw_data$height, raw_data$convex_frac)
plot(raw_data$height, raw_data$shells.cm2)


x_height <- seq(1,37, 0.2)
y_prob <- predict(testmod, list(height = x_height),type="response", se.fit=1)

#final plot
ggplot() +
  geom_ribbon(aes(ymin = y_prob$fit-y_prob$se.fit, ymax = y_prob$fit+y_prob$se.fit, x=x_height), fill = "grey70", alpha=0.25) +
  geom_point(aes(x=height, y=convex/total.shells, size=total.shells), data=raw_data, color="steelblue", alpha=0.8) +
  scale_size(breaks = c(1, 5, 10, 20, 30)) +
  geom_line(aes(y=y_prob$fit, x=x_height), size=1, color="red",linetype="solid") +
  labs(x="Stratigraphic position (m)", y="Convex fraction", color="Legend", size="Sample size") +
  theme_minimal()

plot(convex/total.shells ~ height, data=raw_data)
lines(y_prob$fit ~ x_height)
lines((y_prob$fit+y_prob$se.fit) ~ x_height, col="red")
lines((y_prob$fit-y_prob$se.fit) ~ x_height, col="red")

ggplot(aes(succ=convex, fail=concave), data=raw_data) +
  geom_point(aes(x=height, y=convex/total.shells, size=total.shells), data=raw_data) + 
  geom_smooth(aes(x=height, y=convex/total.shells), data=raw_data, method = "glm",
              method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x) 

ggplot() +
  geom_ribbon(aes(ymin = y_prob$fit-y_prob$se.fit, ymax = y_prob$fit+y_prob$se.fit, x=x_height), fill = "grey70", alpha=0.25) +
  geom_point(aes(x=height, y=convex/total.shells, size=total.shells), data=raw_data, color="steelblue", alpha=0.8) +
  geom_line(aes(y=y_prob$fit, x=x_height), size=1, color="#F8766D",linetype="solid") +
  labs(x="Stratigraphic position (m)", y="Convex fraction", color="Legend", size="Sample size") +
  theme_minimal()

ggplot(data=raw_data) +
  geom_point(aes(y=shells.cm2, x=total.shells, colour="steelblue")) +
  theme_minimal()
  
ggplot(aes(succ=convex, fail=concave), data=raw_data) +
  geom_point(aes(x=height, y=convex/total.shells, size=total.shells), data=raw_data) + 
  geom_smooth(aes(x=height, y=convex/total.shells), data=raw_data, method = "loess") 



#(convex, concave) ~ total.shells
testmod2 <- glm (cbind(convex, concave) ~ total.shells, family = binomial, data=raw_data)
summary(testmod2)

x_height <- seq(1,34, 0.2)
y_prob <- predict(testmod2, list(total.shells = x_shells),type="response", se.fit=1)
plot(convex/total.shells ~ total.shells, data=raw_data)
lines(y_prob$fit ~ x_height)
lines((y_prob$fit+y_prob$se.fit) ~ x_height, col="red")
lines((y_prob$fit-y_prob$se.fit) ~ x_height, col="red")



x_shells <- seq(0.05,1.26, 0.01)
y_prob <- predict(testmod3, list(shells.cm2 = x_shells),type="response", se.fit=1)
plot(convex/total.shells ~ shells.cm2, data=raw_data)
lines(y_prob$fit ~ x_shells)
lines((y_prob$fit+y_prob$se.fit) ~ x_shells, col="red")
lines((y_prob$fit-y_prob$se.fit) ~ x_shells, col="red")



x_height <- seq(0,1.2, 0.01)

y_prob

raw_data$model <- predict(testmod, list(height = x_height),type="response")


plot(convex/total.shells~ height, data=raw_data)
plot(concave/total.shells ~ height, data=raw_data)

y_prob <- predict(testmod, list(height = x_height),type="response")




testmod_s <- glm (total.shells ~ height, family = poisson, data=raw_data)
plot(testmod_s)

summary(testmod_s)



exp(-0.03106)

x_height <- seq(0,40, 0.5)
y_prob <- predict(testmod, list(height = x_height),type="response")

mod <- lm(total.shells ~ height, data=raw_data)
summary(mod)

cor(raw_data$total.shells, raw_data$height)
cor(raw_data$convex/raw_data$total.shells, raw_data$height)
cor(raw_data$convex/raw_data$total.shells, raw_data$total.shells)


cor(raw_data$concave/raw_data$total.shells, raw_data$height)

linmod <- lm(convex/total.shells~ height, data=raw_data)
plot(linmod)
summary(linmod)

testmod2 <- glm (cbind(convex, concave) ~ total.shells, family = binomial, data=raw_data)
summary(testmod2)
plot(testmod2)
testmod3 <- glm (cbind(convex, concave) ~ shells.cm2, family = binomial, data=raw_data)
summary(testmod3)
plot(testmod3)

plot(x_height, y_prob)

