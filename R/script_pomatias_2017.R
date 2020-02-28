########################################
## Data analysis script for:
## Increased population density depresses activity but does not influence dispersal in the snail Pomatias elegans
## Maxime Dahirel, Loic Menut, Armelle Ansart
## script author: Maxime Dahirel
########################################

## Load necessary packages----
library(here)
library(arm)
library(matrixStats)
library(tidyverse)
library(rstan)
library(tidybayes)
library(bayesplot)
library(brms)
library(patchwork)
rstan_options(auto_write = TRUE)
options(mc.cores = 2)

## Load raw dataset
data0 <- read.table(here("data/dataset_main_pomatias.txt"), header = TRUE)
shell_info <- read.table(here("data/dataset_shellmass_pomatias.txt"), header = TRUE)
# data_discrimin <- read.table(here("/discrimin_pomatias.txt"), header = TRUE)
## the data_discrimin dataset contains info of the snails used for the initial LDA for sex (see manuscript), but is not used in the final analysis


######## DATA WRANGLING----

## creation of new variables
data <- as_tibble(data0) %>%
  mutate(
    BOX = factor(interaction(Date, Box_code)),
    Disp2 = -0.5 + as.numeric(Disp == "yes"), ## centered dummy variable ##for supplementary material only
    is.female = -0.5 + as.numeric(Sex_true == "F"), ## centered dummy variable
    scale_Density = scale(Density),
    uniqueID = as.numeric(str_split_fixed(ID, "_", 3)[, 2])
  ) %>%
  left_join(shell_info) %>%
  mutate(behave_Disp = as.numeric(Disp == "yes"), behave_Active = Active) %>%
  select(Date, BOX, uniqueID, ID,
         Density, scale_Density,
         Sex_predicted,
         Sex_dissection = Sex_true, is.female, ### rename Sex variable to something more "accurate"
         Shell_height = Height, Shell_width = Diametre, Aperture_height = PeristomeH, ### renamed morpho variables to something more accurate
         Shell_Mass1, Shell_Mass2,
         Disp2, behave_Active, behave_Disp, Time_sec
  )


data1 <- data %>%
  select(-c(behave_Active, behave_Disp)) %>%
  pivot_longer(
    names_prefix = "Shell_Mass",
    names_to = "measurement_order",
    values_to = "Shell_mass",
    cols = contains("Shell_Mass")
  ) %>%
  mutate(measurement_order = as.numeric(measurement_order))

data2 <- data %>%
  select(-c(Shell_Mass1, Shell_Mass2)) %>%
  pivot_longer(
    names_to = "which_behaviour",
    values_to = "behave",
    cols = contains("behave")
  ) %>%
  mutate(measurement_order = as.numeric(which_behaviour == "behave_Active") + 1)


data <- left_join(data1, data2) %>%
  mutate(which_behaviour = fct_recode(which_behaviour, Active = "behave_Active", Disp = "behave_Disp"))

data <- data %>%
  mutate(
    scale_logshellmass = scale(log(Shell_mass)),
    scale_logheight = scale(log(Shell_height))
  ) # %>%
# mutate(scale_logshellmass_NAfilled=replace_na(scale_logshellmass,0), ###needed for the subset() approach
#     is.shell.valid=as.numeric(is.na(Shell_Mass)==FALSE))  ### the is.shell.valid variable ensures the filled NAs are ignored during fitting for the subset() approach


##### MAIN MODEL----
#### first, our priors
prior <- c(
  set_prior("normal(0,1)", class = "Intercept", resp = c("scalelogshellmass")),
  set_prior("normal(0,1)", class = "b", resp = c("behave", "scalelogshellmass")),
  set_prior("normal(0,1.5)", class = "b", resp = "behave", coef = c("which_behaviourActive", "which_behaviourDisp")),
  set_prior("normal(0,1)", class = "sd", resp = c("behave", "scalelogshellmass")),
  set_prior("normal(0,1)", class = "sigma", resp = "scalelogshellmass"),
  set_prior("lkj(2)", class = "cor")
)

#### then, our submodels

bf_behave <- bf(behave ~ 0 + which_behaviour + which_behaviour:(scale_Density + scale_logheight + is.female) +
                  (0 + which_behaviour | BOX) + (1 | q | ID), family = bernoulli)

bf_shell_miss <- bf(scale_logshellmass | mi() ~ scale_logheight + is.female + (1 | q | ID))

## an alternative to missing data imputation with mi() would be to use the subset() approach, essentially ignoring the NA values, but only for the response with NA
## that is, no deletion of entire rows
### leads to sensibly the same results in our case (as expected when missing data are only in responses)
# bf_shell_subset<-bf(scale_logshellmass_NAfilled|subset(is.shell.valid) ~ scale_logheight + is.female+(1|q|ID))

mod <- brm(mvbf(bf_behave + bf_shell_miss),
           data = data, prior = prior, iter = 25000, warmup = 5000,
           control = list(adapt_delta = 0.95, max_treedepth = 15), chains = 4, seed = 42
)

#### SOME CHECKS----
summary(mod) ## everything looks good

## but summary.brmsfit gives mean and 95% quantile interval
## I would like mean and HDI instead

summary_mod <- mod %>%
  posterior_samples() %>%
  select(starts_with(c("Intercept", "b_", "sd_", "cor_", "sigma"))) %>%
  pivot_longer(everything()) %>%
  group_by(name) %>%
  mean_hdi()

summary_mod

mcmc_rank_overlay(mod, pars = summary_mod$name)
# a bit wonky at the extreme ranks for the cor_ID and sd_ID_behave, but still okay-ish hopefully
# pretty good everywhere else

pp_check(mod, resp = "behave") ## not really informative given it's a binary variable
pp_check(mod, resp = "scalelogshellmass", newdata = data %>% filter(is.na(scale_logshellmass) == FALSE))
### need to make the pp_check on a dataset without the NAs or fail
### this one is not really informative either; we have a random effect of ID so it's a given it should look good
### unless we'd spectacularly failed with the residual variation
### but still good to check
### see ?pp_check for other possibilities


###FIGURE 2----
newdata <- data %>%
  select(which_behaviour) %>%
  distinct() %>%
  expand_grid(BOX = data$BOX[1], ID = data$ID[1], is.female = 0, scale_logheight = 0, Density = c(1:30)) %>%
  mutate(scale_Density = (Density - attr(data$scale_Density, "scaled:center")) / attr(data$scale_Density, "scaled:scale")) %>%
  add_fitted_draws(mod, resp = "behave", re_formula = NA)


p1 <- newdata %>%
  filter(which_behaviour == "Disp") %>%
  ggplot() +
  stat_lineribbon(aes(x = Density, y = .value), .width = 0.95, point_interval = mean_hdi, fill = "grey90", alpha = 0.9) +
  geom_point(data = data %>% filter(which_behaviour == "Disp") %>% group_by(Density, BOX) %>% summarise(behave = mean(behave)), aes(x = Density, y = behave), col = "grey60") +
  geom_point(data = data %>% filter(which_behaviour == "Disp") %>% group_by(Density) %>% summarise(behave = mean(behave)), aes(x = Density, y = behave), cex = 3) +
  scale_x_continuous("") +
  scale_y_continuous("Dispersal rate")

p2 <- newdata %>%
  filter(which_behaviour == "Active") %>%
  ggplot() +
  stat_lineribbon(aes(x = Density, y = .value), .width = 0.95, point_interval = mean_hdi, fill = "grey90", alpha = 0.9) +
  geom_point(data = data %>% filter(which_behaviour == "Active") %>% group_by(Density, BOX) %>% summarise(behave = mean(behave)), aes(x = Density, y = behave), col = "grey60") +
  geom_point(data = data %>% filter(which_behaviour == "Active") %>% group_by(Density) %>% summarise(behave = mean(behave)), aes(x = Density, y = behave), cex = 3) +
  scale_x_continuous("Density (snails per box)") +
  scale_y_continuous("Probability of activity")

(p1 / p2) + plot_annotation(tag_levels = "A") & theme_bw() &
  theme(legend.position = "none", panel.border = element_blank(), axis.line = element_line(colour = "black"))

##### FIGURE 3----
behave_wide <- data0 %>%
  select(Disp, Active, ID) %>%
  mutate(Disp = as.numeric(Disp == "yes"))


newdata <- as_tibble(ranef(mod, summary = FALSE)$ID[, , "behave_Intercept"]) %>%
  mutate(
    post_act_ifdisp = select(., all_of(subset(behave_wide$ID, behave_wide$Disp == 1))) %>% rowMeans(),
    post_act_ifres = select(., all_of(subset(behave_wide$ID, behave_wide$Disp == 0))) %>% rowMeans(),
    post_disp_ifact = select(., all_of(subset(behave_wide$ID, behave_wide$Active == 1))) %>% rowMeans(),
    post_disp_ifstill = select(., all_of(subset(behave_wide$ID, behave_wide$Active == 0))) %>% rowMeans(),
  ) %>%
  cbind(posterior_samples(mod) %>%
          select(meanlogitActive = b_behave_which_behaviourActive, meanlogitDisp = b_behave_which_behaviourDisp)) %>%
  mutate(
    post_act_ifdisp = post_act_ifdisp + meanlogitActive, post_act_ifres = post_act_ifres + meanlogitActive,
    post_disp_ifact = post_disp_ifact + meanlogitDisp, post_disp_ifstill = post_disp_ifstill + meanlogitDisp
  ) %>%
  select(post_act_ifdisp, post_act_ifres, post_disp_ifact, post_disp_ifstill)

p3 <- newdata %>%
  select(post_act_ifdisp, post_act_ifres) %>%
  pivot_longer(everything(), names_to = "Disp") %>%
  mutate(Disp = fct_recode(Disp, `0` = "post_act_ifres", `1` = "post_act_ifdisp")) %>%
  ggplot() +
  geom_col(data = behave_wide %>%
             group_by(Disp) %>%
             summarise(Active = mean(Active == 1)), aes(x = factor(Disp), y = Active), col = "black", fill = "white") +
  stat_eye(aes(x = Disp, y = invlogit(value)),
           .width = 0.95, interval_size = 1, point_size = 3, point_interval = mean_hdi, alpha = 0.9
  ) +
  scale_x_discrete("Dispersal status", labels = c("Resident", "Disperser")) +
  scale_y_continuous("Activity probability", limits = c(0, 1))

p4 <- newdata %>%
  select(post_disp_ifact, post_disp_ifstill) %>%
  pivot_longer(everything(), names_to = "Active") %>%
  mutate(Active = fct_recode(Active, `0` = "post_disp_ifstill", `1` = "post_disp_ifact")) %>%
  ggplot() +
  geom_col(data = behave_wide %>%
             group_by(Active) %>%
             summarise(Disp = mean(Disp)), aes(x = factor(Active), y = Disp), col = "black", fill = "white") +
  stat_eye(aes(x = Active, y = invlogit(value)),
           .width = 0.95, interval_size = 1, point_size = 3, point_interval = mean_hdi, alpha = 0.9
  ) +
  scale_x_discrete("Active post dispersal?", labels = c("no", "yes")) +
  scale_y_continuous("Dispersal rate", limits = c(0, 1))

p3 + p4 &
  theme_bw() &
  theme(legend.position = "none", panel.border = element_blank(), axis.line = element_line(colour = "black"))
### it looks like there is a link

### let's confirm it!
(newdata$post_disp_ifact - newdata$post_disp_ifstill) %>% mean_hdi()
(newdata$post_act_ifdisp - newdata$post_act_ifres) %>% mean_hdi()

##### FIGURE 4----

as_tibble(ranef(mod, summary = FALSE)$ID[, , "scalelogshellmass_Intercept"]) %>%
  mutate(
    post_disp_act = select(., all_of(subset(behave_wide$ID, behave_wide$Disp == 1 & behave_wide$Active == 1))) %>% rowMeans(na.rm = TRUE),
    post_res_act = select(., all_of(subset(behave_wide$ID, behave_wide$Disp == 0 & behave_wide$Active == 1))) %>% rowMeans(na.rm = TRUE),
    post_disp_still = select(., all_of(subset(behave_wide$ID, behave_wide$Disp == 1 & behave_wide$Active == 0))) %>% rowMeans(na.rm = TRUE),
    post_res_still = select(., all_of(subset(behave_wide$ID, behave_wide$Disp == 0 & behave_wide$Active == 0))) %>% rowMeans(na.rm = TRUE),
  ) %>%
  select(post_disp_act, post_disp_still, post_res_act, post_res_still) %>%
  pivot_longer(everything()) %>%
  mutate(name = fct_recode(name,
                           `Disperser and Active` = "post_disp_act",
                           `Resident and Active` = "post_res_act",
                           `Disperser and Inactive` = "post_disp_still",
                           `Resident and Inactive` = "post_res_still"
  )) %>%
  ggplot() +
  geom_hline(yintercept = 1, lty = 2) +
  stat_eye(aes(x = name, y = exp(value * attr(data$scale_logshellmass, "scaled:scale"))),
           .width = 0.95, interval_size = 1, point_size = 3, point_interval = mean_hdi, alpha = 0.9
  ) +
  scale_x_discrete("") +
  scale_y_continuous("Average relative shell mass (posterior)", breaks = c(0.95, 1, 1.05)) +
  theme_bw() +
  theme(legend.position = "none", panel.border = element_blank(), axis.line = element_line(colour = "black"))

##### SUPPLEMENTARY MATERIALS----


mod_S1a <- data %>%
  filter(measurement_order == 1) %>% ## filter to avoid duplicates
  group_by(BOX) %>%
  summarise(Density = mean(Density), Nfemales = sum(Sex_dissection == "F")) %>%
  ungroup() %>%
  brm(Nfemales | trials(Density) ~ 1,
      family = binomial,
      prior = c(set_prior("normal(0,1.5)", class = "Intercept")), data = ., seed = 42
  )

mod_S1a %>%
  posterior_samples() %>%
  select(contains("Intercept")) %>%
  pivot_longer(everything()) %>%
  group_by(name) %>%
  mutate(value = invlogit(value)) %>%
  mean_hdi()

mod_S1b <- data %>%
  filter(measurement_order == 1) %>%
  group_by(BOX) %>%
  summarise(Density = mean(Density), Nfemales = sum(Sex_dissection == "F")) %>%
  ungroup() %>%
  brm(Nfemales | trials(Density) ~ 0 + factor(Density),
      family = binomial,
      prior = c(set_prior("normal(0,1.5)", class = "b")), data = ., seed = 42
  )

mod_S1b %>%
  posterior_samples() %>%
  select(starts_with(c("Intercept", "b_"))) %>%
  pivot_longer(everything()) %>%
  group_by(name) %>%
  mutate(value = invlogit(value)) %>%
  mean_hdi()


mod_S3 <- data %>%
  filter(which_behaviour == "Active") %>%
  brm(behave ~ scale_Density + scale_logheight + is.female + Disp2 + (1 | BOX),
      family = bernoulli,
      prior = c(
        set_prior("normal(0,1.5)", class = "Intercept"),
        set_prior("normal(0,1)", class = "b"),
        set_prior("normal(0,1)", class = "sd")
      ), data = ., seed = 42
  )

mod_S3 %>%
  posterior_samples() %>%
  select(starts_with(c("Intercept", "b_", "sd_"))) %>%
  pivot_longer(everything()) %>%
  group_by(name) %>%
  mean_hdi()
