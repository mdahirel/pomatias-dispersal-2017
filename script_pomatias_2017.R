## TO DO: metadata to add (here or separate file)
##check consistiency of data import

########################################
## Data analysis script for:
## Increased population density depresses activity but does not influence dispersal in the snail Pomatias elegans
## Maxime Dahirel, Loic Menut, Armelle Ansart
## script author: Maxime Dahirel
########################################

## Load necessary packages
library(arm)
library(matrixStats)
library(tidyverse)
library(rstan)
library(tidybayes)
library(brms)
library(patchwork)
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
data<-as_tibble(data) %>% 
  mutate(BOX = factor(interaction(Date, Box_code)),
         Disp2 = -0.5 + as.numeric(Disp == "yes"), ## centered dummy variable ##for supplementary material only
         is.female = -0.5 + as.numeric(Sex_true == "F"), ## centered dummy variable
         scale_Density = scale(Density),
         uniqueID = as.numeric(str_split_fixed(ID,"_",3)[,2])
  ) %>% 
  left_join(morpho) %>% 
  mutate(behave_Disp=as.numeric(Disp=="yes"),behave_Active=Active) %>% 
  dplyr::select(Date,BOX,uniqueID,
         Density,scale_Density,
         Sex_predicted,Sex_dissection=Sex_true,is.female,
         Height,Diametre,PeristomeH,
         Shell_Mass1,Shell_Mass2,
         Disp2,behave_Active,behave_Disp,ID)

data1<-data %>% 
   dplyr::select(-c(behave_Active, behave_Disp)) %>% 
  pivot_longer(names_to="mass_measurement",
               values_to="Shell_Mass",
               cols=contains("Shell_Mass")) %>% 
  mutate(measurement_order=as.numeric(mass_measurement=="Shell_Mass2")+1)
data2<-data %>% 
   dplyr::select(-c(Shell_Mass1,Shell_Mass2)) %>% 
   pivot_longer(names_to="which_behaviour",
                values_to="behave",
                cols=contains("behave")) %>% 
  mutate(measurement_order=as.numeric(which_behaviour=="behave_Active")+1)

### narrower ranef for disp and active because binomial?
data<-left_join(data1,data2) %>% 
  mutate(scale_shell=scale(Shell_Mass),
       scale_height=scale(Height),
       scale_logshellmass=scale(log(Shell_Mass)),
       scale_logheight=scale(log(Height))) %>% 
  mutate(scale_logshellmass_NAfilled=replace_na(scale_logshellmass,0), ###needed for the subset() approach
       is.shell.valid=as.numeric(is.na(Shell_Mass)==FALSE)) ### the is.shell.valid variable ensures the filled NAs are ignored during fitting

###note, check how scaled vars are entered in tibbles

####models for supplementary

mod_S1a <- data %>% filter(mass_measurement=="Shell_Mass1") %>% 
  group_by(BOX) %>% 
  summarise(Density=mean(Density),Nfemales=sum(Sex_dissection=="F")) %>% 
  ungroup() %>% 
  brm(Nfemales|trials(Density)~1,family=binomial,
      prior=c(set_prior("normal(0,1.5)",class="Intercept")),data=., seed = 42
  )

mod_S1b <- data %>% filter(mass_measurement=="Shell_Mass1") %>% 
  group_by(BOX) %>% 
  summarise(Density=mean(Density),Nfemales=sum(Sex_dissection=="F")) %>% 
  ungroup() %>% 
  brm(Nfemales|trials(Density)~0+factor(Density),family=binomial,
      prior=c(set_prior("normal(0,1.5)",class="b")),data=., seed = 42
  )

mod_S3 = brm(behave~scale_Density+scale_logheight+is.female+Disp2+(1|BOX),family=bernoulli,
             prior=c(set_prior("normal(0,1.5)",class="Intercept"),
                     set_prior("normal(0,1)",class="b"),
                     set_prior("normal(0,1)",class="sd")),data=subset(data,which_behaviour=="behave_Active"), seed = 42
)

####main model

prior=c(
  set_prior("normal(0,1)",class="Intercept",resp=c("scalelogshellmass")),
  set_prior("normal(0,1)",class="b",resp=c("behave","scalelogshellmass")),
  set_prior("normal(0,1.5)",class="b",resp="behave",coef=c("which_behaviourbehave_Active","which_behaviourbehave_Disp")),
  set_prior("normal(0,1)",class="sd",resp=c("behave","scalelogshellmass")),
  set_prior("normal(0,1)",class="sigma",resp="scalelogshellmass"), #to relax?? hist(exp(rexp(1000,10)*0.19))
  set_prior("lkj(2)",class="cor")
)

bf_behave=bf(behave~0+which_behaviour+which_behaviour:(scale_Density+scale_logheight+is.female)+(1|p|BOX)+(1|q|ID),family=bernoulli)

bf_shell_miss=bf(scale_logshellmass|mi() ~ scale_logheight + is.female+(1|q|ID))
bf_shell_subset=bf(scale_logshellmass_NAfilled|subset(is.shell.valid) ~ scale_logheight + is.female+(1|q|ID))

mod=brm(mvbf(bf_behave+bf_shell_miss),data=data,prior=prior,iter=250, warmup=50,
        control=list(adapt_delta=0.95,max_treedepth=15),chains=4,seed=42) 

















newdata=data2 %>%
  select(behav_type) %>% 
  distinct() %>%
  expand_grid(BOX=data2$BOX[1],ID=data2$ID[1],is.female=0,scale_logheight=0,Density=c(1:30)/1,ShellMeasure="Shell_Mass1",is.shell.valid=1) %>% 
  mutate(scale_Density=(Density-attr(data$scale_Density,"scaled:center"))/attr(data$scale_Density,"scaled:scale")) %>% 
  add_fitted_draws(mod,resp="behave",re_formula = NA)

p2<-newdata %>% 
  filter(behav_type=="Active") %>% 
  ggplot()+
  stat_lineribbon(aes(x=Density,y=.value),.width=0.95,fill="lightgrey") + 
  geom_jitter(data=data %>% group_by(Density,BOX) %>% summarise(Active=mean(Active)),aes(x=Density,y=Active),col="grey")+
  geom_point(data=data %>% group_by(Density) %>% summarise(Active=mean(Active)),aes(x=Density,y=Active),cex=3)+
  scale_x_continuous("Density (snails per box)")+
  scale_y_continuous("Probability of activity")+
  theme_bw()+theme(legend.position = "none")

p1<-newdata %>% 
  filter(behav_type=="Disp") %>% 
  ggplot()+
  stat_lineribbon(aes(x=Density,y=.value),.width=0.95,fill="lightgrey") + 
  geom_jitter(data=data %>% group_by(Density,BOX) %>% summarise(Disp=mean(Disp=="yes")),aes(x=Density,y=Disp),col="grey")+
  geom_point(data=data %>% group_by(Density) %>% summarise(Disp=mean(Disp=="yes")),aes(x=Density,y=Disp),cex=3)+
  scale_x_continuous("")+
  scale_y_continuous("Dispersal rate")+
  theme_bw()+theme(legend.position = "none")

p3<-posterior_samples(mod) %>% transmute(cor=cor_ID__behave_Intercept__scalelogshell_Intercept) %>% ggplot()+
  geom_halfeyeh(aes(y=1,x=cor))+
  scale_x_continuous("individual-level correlation between behaviour probability and (relative) shell mass",limits = c(-1,1))+
  scale_y_continuous("",breaks=NULL)+
  geom_vline(xintercept=0,lty=2)+
  theme_bw()
###shell mass is a collider and cannot be included in model see DAG below
### thick is thickness (unobserved)

(p1/p2)+plot_annotation(tag_levels="A")

newdata<-data %>% select(Disp,Active,ID) %>% distinct()
newdata<-as_tibble(ranef(mod,summary = FALSE)$ID[,,"behave_Intercept"]) %>% 
  mutate(
    post_disp=select(.,all_of(subset(newdata$ID,newdata$Disp=="yes"))) %>% rowMeans(),
    post_res=select(.,all_of(subset(newdata$ID,newdata$Disp=="no"))) %>% rowMeans(),
    post_act=select(.,all_of(subset(newdata$ID,newdata$Active==1))) %>% rowMeans(),
    post_still=select(.,all_of(subset(newdata$ID,newdata$Active==0))) %>% rowMeans(),
  ) %>% 
  cbind(posterior_samples(mod) %>% 
          select(meanlogitActive=b_behave_behav_typeActive,meanlogitDisp=b_behave_behav_typeDisp)) %>% 
  mutate(post_disp=post_disp+meanlogitActive,post_res=post_res+meanlogitActive,
         post_act=post_act+meanlogitDisp,post_still=post_still+meanlogitDisp) %>% 
  select(post_disp,post_res,post_act,post_still)


p4 <- newdata %>% 
  select(post_act,post_still) %>% 
  pivot_longer(everything(),names_to="Active") %>%
  mutate(Active=fct_recode(Active,`0`="post_still",`1`="post_act")) %>% 
  ggplot()+
  geom_col(data=data %>% select(Disp,Active,ID) %>% 
             distinct() %>% 
             group_by(Active) %>% 
             summarise(Disp=mean(Disp=="yes")),aes(x=factor(Active),y=Disp),col="black",fill="white")+
  geom_eye(aes(x=Active,y=invlogit(value)),.width=c(0,0.95))+
  scale_x_discrete("Active post dispersal?",labels=c("no","yes"))+
  scale_y_continuous("Dispersal rate",limits=c(0,1))+
  theme_bw()

(newdata$post_act-newdata$post_still) %>% median_hdi()
(newdata$post_disp-newdata$post_res) %>% median_hdi()

p5 <- newdata %>% 
  select(post_disp,post_res) %>% 
  pivot_longer(everything(),names_to="Disp") %>%
  mutate(Disp=fct_recode(Disp,no="post_res",yes="post_disp")) %>% 
  ggplot()+
  geom_col(data=data %>% select(Disp,Active,ID) %>% 
             distinct() %>% 
             group_by(Disp) %>% 
             summarise(Active=mean(Active==1)),aes(x=Disp,y=Active),col="black",fill="white")+
  geom_eye(aes(x=Disp,y=invlogit(value)),.width=c(0,0.95))+
  scale_x_discrete("Dispersal status",labels=c("Resident","Disperser"))+
  scale_y_continuous("Activity probability",limits=c(0,1))+
  theme_bw()

p5+p4


newdata<-data %>% select(Disp,Active,ID) %>% distinct()
newdata<-as_tibble(ranef(mod,summary = FALSE)$ID[,,"scalelogshell_Intercept"]) %>% 
  mutate(
    post_disp=select(.,all_of(subset(newdata$ID,newdata$Disp=="yes"))) %>% rowMeans(na.rm=TRUE),
    post_res=select(.,all_of(subset(newdata$ID,newdata$Disp=="no"))) %>% rowMeans(na.rm=TRUE),
    post_act=select(.,all_of(subset(newdata$ID,newdata$Active==1))) %>% rowMeans(na.rm=TRUE),
    post_still=select(.,all_of(subset(newdata$ID,newdata$Active==0))) %>% rowMeans(na.rm=TRUE),
  ) %>% 
  cbind(posterior_samples(mod) %>% select(meanlogMass=b_scalelogshell_Intercept)) %>% 
  mutate(post_disp=post_disp+meanlogMass,post_res=post_res+meanlogMass,
         post_act=post_act+meanlogMass,post_still=post_still+meanlogMass) %>% 
  select(post_disp,post_res,post_act,post_still)


p6 <- newdata %>% 
  select(post_act,post_still) %>% 
  pivot_longer(everything(),names_to="Active") %>%
  mutate(Active=fct_recode(Active,`0`="post_still",`1`="post_act")) %>% 
  ggplot()+
  #geom_col(data=data %>% select(Disp,Active,ID) %>% 
  #           distinct() %>% 
  #           group_by(Active) %>% 
  #           summarise(Disp=mean(Disp=="yes")),aes(x=factor(Active),y=Disp),col="black",fill="white")+
  geom_eye(aes(x=fct_rev(Active),y=exp(value*attr(data$scale_logshell,"scaled:scale"))),.width=c(0,0.95))+
  scale_x_discrete("Active post dispersal?",labels=c("no","yes"))+
  scale_y_continuous("",limits=c(0.95,1.05))+
  geom_hline(yintercept=1,lty=2)+
  theme_bw()

exp((newdata$post_act-newdata$post_still)*attr(data$scale_logshell,"scaled:scale")) %>% median_hdi()
exp((newdata$post_disp-newdata$post_res)*attr(data$scale_logshell,"scaled:scale")) %>% median_hdi()

p7 <- newdata %>% 
  select(post_disp,post_res) %>% 
  pivot_longer(everything(),names_to="Disp") %>%
  mutate(Disp=fct_recode(Disp,no="post_res",yes="post_disp")) %>% 
  ggplot()+
  #geom_col(data=data %>% select(Disp,Active,ID) %>% 
  #           distinct() %>% 
  #           group_by(Disp) %>% 
  #           summarise(Active=mean(Active==1)),aes(x=Disp,y=Active),col="black",fill="white")+
  geom_eye(aes(x=fct_rev(Disp),y=exp(value*attr(data$scale_logshell,"scaled:scale"))),.width=c(0,0.95))+
  scale_x_discrete("Dispersal status",labels=c("Resident","Disperser"))+
  scale_y_continuous("Relative Shell mass",limits=c(0.95,1.05))+
  geom_hline(yintercept=1,lty=2)+
  theme_bw()

p7+p6







