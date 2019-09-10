
## TO DO: metadata to add (here or separate file)
##check consistiency of data import


########################################
## Data analysis script for:
## Increased population density depresses activity but does not influence dispersal in the snail Pomatias elegans
## Maxime Dahirel, Loic Menut, Armelle Ansart
## script author: Maxime Dahirel
########################################

## Load necessary packages
library(lme4)
library(car)
library(arm)
library(merTools)
library(tidyverse)
library(tidyr)
library(cowplot)

library(matrixStats)
library(coda)
library(rstan)
library(tidybayes)
library(brms)
rstan_options(auto_write = TRUE)
options(mc.cores = 2)


## Load raw dataset
data <- read.table("./dataset_pomatias.txt",
  header = TRUE
)

morpho <- read.table("./dataset_morpho_pomatias.txt",
                   header = TRUE
)

data_discrimin <- read.table("./discrimin_pomatias.txt",
                   header = TRUE
)


########################################
## PART 0: Variable creation and experimental design checks
########################################

## creation of new variables
data$BOX <- factor(interaction(data$Date, data$Box_code))
data$Disp2 <- -0.5 + as.numeric(data$Disp == "yes") ## centered dummy variable
data$is.female <- -0.5 + as.numeric(data$Sex_true == "F") ## centered dummy variable
data$scale_Density <- scale(data$Density)

data$uniqueID <- as.numeric(str_split_fixed(data$ID,"_",3)[,2])


data <- left_join(data,morpho)

data<-data %>% gather(key="ShellMeasure",value="Shell_Mass",c("Shell_Mass1","Shell_Mass2"))


### data check: test that pretrial sex evaluation was not biased with respect to density
#### create box-level variables and create a box-level dataset
# Nmales <- by(data$Sex_true == "M", data$BOX, FUN = sum)
# Density <- by(data$Density, data$BOX, FUN = mean)
# Day <- by(data$Day, data$BOX, FUN = unique)
# boxdata <- data.frame(Density = c(Density), Nmales = c(Nmales), Day = c(Day),Week=c(Week))
# boxdata$scale_Density <- scale(boxdata$Density)
#
#### model with density as a continuous variable
# mod_dens <- glmer(cbind(Nmales,Density-Nmales) ~ scale_Density+(1|Day),
#  data = boxdata, family = binomial
# )
# Anova(mod_dens)
#
#### model with density as a factor
# mod_dens_factor <- glmer(cbind(Nmales,Density-Nmales) ~ factor(Density)+(1|Day),
#  data = boxdata, family = binomial
# )
#### note: gives a singular warning apparently because of virtually 0 among-Day variance
#### removing the random effect of Day (doing a GLM instead of a GLMM) leads to same conclusions
# Anova(mod_dens_factor)


########################################
## PART 1: Data analysis
########################################

## Dispersal model
mod_disp <- glmer(Disp ~ (scale_Density) * is.female + (1 | Date / BOX),
  data = data, family = binomial
)
Anova(mod_disp) # ; Anova(mod_disp,type="3")#virtually identical results thanks to centring/scaling

### Activity model
mod_activ <- glmer(Active ~ (scale_Density) * (Disp2 * is.female) + (1 | Date / BOX),
  data = data, family = binomial
)
Anova(mod_activ) # ; Anova(mod_activ,type="3")#virtually identical results


########################################
## PART 2: Generating figures
########################################

#### make a newdata to house predictions every 0.1 snail with everything except density averaged out
newdata <- expand.grid(
  Density = c(10:300) / 10,
  Disp2 = 0,
  is.female = 0,
  Date = data$Date[1], ### needed for predict() to work but not used (predictions on fixed effects only)
  BOX = data$BOX[1] ### same
)
newdata$scale_Density <- (newdata$Density - mean(data$Density)) / sd(data$Density) ## scale based on the *original* dataset SD


### Dispersal as a function of density
### add the predictions and their 95% confidence interval
### include.resid.var set so we have confidence intervals and not prediction intervals
PREDS_disp <- invlogit(
  predictInterval(mod_disp,
    newdata = newdata, which = "fixed",
    level = 0.95, n.sims = 10000, stat = "mean",
    type = "linear.prediction", include.resid.var = FALSE
  )
)
names(PREDS_disp) <- c("fit_disp", "ucl_disp", "lcl_disp")
newdata <- cbind(newdata, PREDS_disp)

plot_disp <- ggplot() +
  stat_summary(data = data, aes(x = Density, y = as.numeric(Disp=="yes"),group=BOX),fun.y = "mean",geom="point", col = "grey", size = 2)+
  stat_summary(data = data, aes(x = Density, y = as.numeric(Disp=="yes")),fun.y = "mean",geom="point", size = 3)+
  geom_line(data = newdata, aes(x = Density, y = fit_disp)) +
  geom_line(data = newdata, aes(x = Density, y = lcl_disp),lty=2) +
  geom_line(data = newdata, aes(x = Density, y = ucl_disp),lty=2) +
  scale_y_continuous(name = "Dispersal rate", limits = c(0, 1)) +
  scale_x_continuous(name = "", limits = c(0, 30)) + ## no axis name because present on activity subplot
  theme_classic()


### Activity as a function of density

PREDS_active <- invlogit(
  predictInterval(mod_activ,
    newdata = newdata, which = "fixed",
    level = 0.95, n.sims = 10000, stat = "mean",
    type = "linear.prediction", include.resid.var = FALSE
  )
)
names(PREDS_active) <- c("fit_activity", "ucl_activity", "lcl_activity")
newdata <- cbind(newdata, PREDS_active)

plot_activity <- ggplot() +
  stat_summary(data = data, aes(x = Density, y = Active,group=BOX),fun.y = "mean",geom="point", col = "grey", size = 2)+
  stat_summary(data = data, aes(x = Density, y = Active),fun.y = "mean",geom="point", size = 3)+
  geom_line(data = newdata, aes(x = Density, y = fit_activity)) +
  geom_line(data = newdata, aes(x = Density, y = lcl_activity),lty=2) +
  geom_line(data = newdata, aes(x = Density, y = ucl_activity),lty=2) +
  scale_y_continuous(name = "Activity probability", limits = c(0, 1)) +
  scale_x_continuous(name = "Population density (snails per box)", limits = c(0, 30)) +
  theme_classic()


### join the two subplots in a common figures using cowplot
plot_grid(plot_disp, plot_activity, labels = c("A", "B"), nrow = 2)
