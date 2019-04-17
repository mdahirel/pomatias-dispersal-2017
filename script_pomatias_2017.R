library(lme4)
library(car)
library(arm)
library(merTools)
library(tidyverse)
library(tidyr)
library(cowplot)

data <- read.table("D:/Maxime/Documents/ATER 2016-2018/RECHERCHE/projet_pomatias/dataset_pomatias.txt",
                   header=TRUE)

data$BOX <- factor(interaction(data$Date,data$Box_code))
data$Disp2 <- -0.5 + as.numeric(data$Disp == "yes")
data$is.female <- -0.5 + as.numeric(data$Sex_true == "F")
data$scale_Density <- scale(data$Density)
data$scale_Density2 <-(data$scale_Density)^2

### test that pretrial sex evaluation was not biased with respect to density
Nmales <- by(data$Sex_true == "M", data$BOX, FUN = sum)
Density <- by(data$Density, data$BOX, FUN = mean)
Day <- by(data$Day, data$BOX, FUN = unique)
Week <- by(data$Week, data$BOX, FUN = unique)
plot(Density, Nmales / Density)
plot(factor(Density), Nmales / Density)
boxdata <- data.frame(Density = c(Density), Nmales = c(Nmales), Day = c(Day),Week=c(Week))
boxdata$scale_Density <- scale(boxdata$Density)


mod1 <- glmer(cbind(Nmales,Density-Nmales) ~ scale_Density+(1|Day),
  data = boxdata, family = binomial
)

Anova(mod1)

mod1factor <- glmer(cbind(Nmales,Density-Nmales) ~ factor(Density)+(1|Day),
  data = boxdata, family = binomial
)
###gives a singular warning because no Day variance, removing the random effect of Day leads to same conclusions
Anova(mod1factor)


## Dispersal model v1
mod2 <- glmer(Disp ~ (scale_Density) * is.female + (1|Date/BOX),
  data = data, family = binomial
)

Anova(mod2)#; Anova(mod2,type="3")#virtually identical results

### activity model v1
mod3 <- glmer(Active ~ (scale_Density) * (Disp2 * is.female) + (1|Date/BOX),
  data = data, family = binomial
)

Anova(mod3)#; Anova(mod3,type="3")#virtually identical results

### do plots
### WRITE :)







###plots

###dispersal
newdata <- expand.grid(   ##for predictions
  Density = c(10:300)/10,
  Disp2 = 0,
  is.female = 0,
  Date=data$Date[1],
  BOX=data$BOX[1]
)
newdata$scale_Density=(newdata$Density-mean(data$Density))/sd(data$Density)

PREDS=invlogit(predictInterval(mod2,newdata=newdata,which="fixed",level=0.95, n.sims=10000,stat="mean",type="linear.prediction",include.resid.var = FALSE))
names(PREDS)=c("fit_disp","lcl_disp","ucl_disp")
newdata=cbind(newdata,PREDS)

test <- data[,c("Disp","Density")] %>% 
  group_by(Density) %>% 
  summarize(mean_disp=mean(Disp=="yes"))

data$meandisp=ave(data$Disp=="yes",data$BOX,FUN=mean)
plot_disp<-ggplot(data=test)+
  geom_point(data=test,aes(x = Density, y = mean_disp),size=2)+
  geom_point(data=data,aes(x = Density, y = meandisp),col="grey")+
  geom_line(data=newdata,aes(x=Density,y=fit_disp))+
  geom_line(data=newdata,aes(x=Density,y=lcl_disp))+
  geom_line(data=newdata,aes(x=Density,y=ucl_disp))+
  scale_y_continuous(name = "Dispersal rate", limits = c(0, 1)) +
  scale_x_continuous(name = "Population density (snails per box)", limits = c(0, 30)) +
  theme_classic()


###activity
newdata <- expand.grid(   ##for predictions
  Density = c(10:300)/10,
  Disp2 = 0,
  is.female = 0,
  Date=data$Date[1],
  BOX=data$BOX[1]
)
newdata$scale_Density=(newdata$Density-mean(data$Density))/sd(data$Density)

PREDS=invlogit(predictInterval(mod3,newdata=newdata,which="fixed",level=0.95, n.sims=10000,stat="mean",type="linear.prediction",include.resid.var = FALSE))
names(PREDS)=c("fit_activity","lcl_activity","ucl_activity")
newdata=cbind(newdata,PREDS)

test <- data[,c("Active","Density")] %>% 
  group_by(Density) %>% 
  summarize(mean_active=mean(Active))

data$meanactive=ave(data$Active,data$BOX,FUN=mean)

plot_activity<-ggplot(data=test)+
  geom_point(aes(x = Density, y = mean_active),size=2)+
  geom_point(data=data,aes(x = Density, y = meanactive),col="grey")+
  geom_line(data=newdata,aes(x=Density,y=fit_activity))+
  geom_line(data=newdata,aes(x=Density,y=lcl_activity))+
    geom_line(data=newdata,aes(x=Density,y=ucl_activity))+
  scale_y_continuous(name = "Activity probability", limits = c(0, 1)) +
  scale_x_continuous(name = "Population density (snails per box)", limits = c(0, 30)) +
  theme_classic()
  


plot_grid(plot_disp, plot_activity, labels = c("A", "B"), nrow = 2) ##review legend position to be sure plot 1 and 2 are aligned


