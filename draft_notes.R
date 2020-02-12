### narrower ranef for disp and active because binomial?
data$scale_shell=scale(data$Shell_Mass)
data$scale_height=scale(data$Height)
data$scale_logshell=scale(log(data$Shell_Mass))
data$scale_logheight=scale(log(data$Height))
data$scale_logshell[is.na(data$scale_logshell)==TRUE]=0 ####needed for the subset approach (can't have NAs) but not included in model
data$is.shell.valid=as.numeric(is.na(data$Shell_Mass)==FALSE)


####models for supplementary

mod_S1a <- data %>% filter(ShellMeasure=="Shell_Mass1") %>% 
  group_by(BOX) %>% 
  summarise(Density=mean(Density),Nfemales=sum(Sex_true=="F")) %>% 
  ungroup() %>% 
  brm(Nfemales|trials(Density)~1,family=binomial,
             prior=c(set_prior("normal(0,1.5)",class="Intercept")),data=., seed = 42
)

mod_S1b <- data %>% filter(ShellMeasure=="Shell_Mass1") %>% 
  group_by(BOX) %>% 
  summarise(Density=mean(Density),Nfemales=sum(Sex_true=="F")) %>% 
  ungroup() %>% 
  brm(Nfemales|trials(Density)~0+factor(Density),family=binomial,
      prior=c(set_prior("normal(0,1.5)",class="b")),data=., seed = 42
  )



#mod_S1c <- data %>% filter(ShellMeasure=="Shell_Mass1") %>% 
#  brm(is.shell.valid~scale_Density+scale_logheight+is.female,family=bernoulli,
#      prior=c(set_prior("normal(0,1.5)",class="b")),data=., seed = 42
#  )







mod_S3 = brm(Active~scale_Density+scale_logheight+is.female+Disp2+(1|BOX),family=bernoulli,
             prior=c(set_prior("normal(0,1.5)",class="Intercept"),
                     set_prior("normal(0,1)",class="b"),
                     set_prior("normal(0,1)",class="sd")),data=subset(data,ShellMeasure=="Shell_Mass1"), seed = 42
)


prior=c(
  set_prior("normal(0,1)",class="Intercept",resp=c("scalelogshell")),
  set_prior("normal(0,1)",class="b",resp=c("behave","scalelogshell")),
  set_prior("normal(0,1.5)",class="b",resp="behave",coef=c("behav_typeActive","behav_typeDisp")),
  set_prior("normal(0,1)",class="sd",resp=c("behave","scalelogshell")),
  set_prior("normal(0,1)",class="sigma",resp="scalelogshell"), #to relax?? hist(exp(rexp(1000,10)*0.19))
  set_prior("lkj(2)",class="cor")
)

data2=data %>% 
  mutate(Disp=as.numeric(Disp=="yes")) %>% 
  pivot_longer(all_of(c("Disp","Active")),names_to = "behav_type",values_to="behave") %>% 
  mutate(is.activity=-0.5+1*(behav_type=="Active"))

bf_behave=bf(behave|subset(ShellMeasure=="Shell_Mass1")~0+behav_type+behav_type:(scale_Density+scale_logheight+is.female)+(1|p|BOX)+(1|q|ID),family=bernoulli)
bf_shell2=bf(scale_logshell|subset(is.shell.valid==1 & behav_type=="Active") ~ scale_logheight + is.female+(1|q|ID))



mod=brm(mvbf(bf_behave+bf_shell2),data=data2,prior=prior,iter=25000, warmup=5000,
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

library(dagitty)
g=dagitty("dag{
           sex -> size -> behavior;
           sex -> thick -> behavior;
           sex-> behavior;
           thick-> mass;
           size -> mass;
           density -> behavior
           }")



activdag=dagitty("dag{
                 Disp <- Density -> Activity
                 Disp <- Sex -> Activity
                 Disp <- size -> Activity
                 Sex -> Size
                 Disp <-> Activity
                 }")


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












