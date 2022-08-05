## Model prevalence of infected ticks

## Doing this separately because it doesn't play nicely with the path analysis
## --Not that anything does...

library(tidyverse)
library(sjPlot)
library(patchwork)
library(MuMIn)
library(DHARMa)
library(lme4)

sem_data <- read_csv("data/path_analysis_data.csv")

## Ultimate terminal node is prevalence of infected ticks, which we believe is 
## determined primarily by tick density but may also be affected by host density.
## We have three ways to calculate prevalence of ticks infected with human pathogens:
## 1. Only human pathogens, excluding R. amblyommatis
## 2. Human pathogens, including R. amblyommatis
## 3. Human and animal pathogens, including R. amblyommatis

## We also have multiple groupings for relative host abundance
## 1. All host dung clusters
## 2. Deer and other except cottontail
## 3. Only deer

sem_data %>% 
  select(deer:total_clusters1m) %>% 
  colSums( )

sem_data <- sem_data %>%
  mutate(biomass_log = log(avg_dry_standing_gm2),
         d_since_fire_log = log(d_since_fire),
         prev_noRamb_cbnd = cbind(ticks_human_no_Ramb, tick_abundance),
         prev_human_cbnd = cbind(ticks_human, tick_abundance),
         prev_hum_ani_cbnd = cbind(ticks_human_animal, tick_abundance),
         prev_no_Ramb_beta = case_when(
           prev_human_no_Ramb == 0 ~ .0001,
           prev_human_no_Ramb == 1 ~ .9999,
           TRUE ~ prev_human_no_Ramb
         ),
         prev_noRamb_beta_log = log(prev_no_Ramb_beta)
  )

# hist(((sem_data$prev_human_no_Ramb_logit)))

# Prevalence of Infected Ticks
# Only tick abundance and host abundance were hypothesized to directly affect the prevalence of infected ticks. We therefore modeled infection prevalence outside of the path modeling framework because effects of all other factors would be indirect. We tested three models with different measures of infected tick prevalence at each plot (infected ticks / total ticks):

# 1.	Ticks infected with strictly known human pathogens
# 2.	Ticks infected with human pathogens inclusive of R. amblyomattis
# 3.	Ticks infected with human and animal pathogens, inclusive of R. amblyomattis

# Models were fit in a mixed effects framework using the glmer function from the lme4 package in R (Bates et al. 2015) with infected tick prevalence as a binomial logit-link response, tick abundance (ticks per trap) and host abundance as fixed-effects predictors, and plot ID as a random effect. 


## Prevalance of infected ticks models ####
## Strictly human pathogens, no R. amblyommatis
# all hosts
human_no_Ramb_mod_allHosts <- glmer(
  prev_noRamb_cbnd ~ ticks_per_trap + total_clusters1m + (1 | plot_id), 
  data = sem_data, 
  family = binomial(link = "logit"),
  glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=1e6)))
simulateResiduals(human_no_Ramb_mod_allHosts) %>% plot()
plot(human_no_Ramb_mod_allHosts)
plot_model(human_no_Ramb_mod_allHosts, type = "pred")
plot_model(human_no_Ramb_mod_allHosts, type = "est") +
  scale_y_continuous(limits = c(.99, 1.05))
tab_model(human_no_Ramb_mod_allHosts)
summary(human_no_Ramb_mod_allHosts)


library(glmmTMB)
m1 <- glmmTMB(prev_no_Ramb_beta ~ ticks_per_trap + total_clusters1m + (1 | plot_id), 
            data = sem_data, 
            family = beta_family(link = "logit"))
summary(m1)
plot_model(m1, type = "pred", show.data = T)
plot_model(m1, type = "est") + scale_y_continuous(limits = c(.99, 1.1))
simulateResiduals(m1) %>% plot()

## Models pretty severely violating assumptions
m2 <- lmer(prev_noRamb_beta_log ~ tpt + total_clusters1m + (1 | plot_id), 
           data = sem_data)
plot(m2)
simulateResiduals(m2) %>% plot()
summary(m2)
effectsize::standardize(m2) %>% summary()
coefs(human_no_Ramb_mod_allHosts)
coefs(m2)
## path modeling says litter cover is a missing path
m3 <- lmer(prev_noRamb_beta_log ~ tpt + total_clusters1m + logit_litter + (1 | plot_id), 
           data = sem_data)
plot(m3)
simulateResiduals(m3) %>% plot()
coefs(m3)

# hosts without rabbit
human_no_Ramb_mod_noRabbit <- glmer(
  cbind(ticks_human_no_Ramb, tick_abundance) ~ ticks_per_trap + deer_other + (1 | plot_id), 
  data = sem_data, family = binomial(link = "logit"),
  glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=1e6)))
# summary(human_no_Ramb_mod_noRabbit)

# only deer
human_no_Ramb_mod_deer <- glmer(
  cbind(ticks_human_no_Ramb, tick_abundance) ~ ticks_per_trap + deer + (1 | plot_id), 
  data = sem_data, family = binomial(link = "logit"),
  glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=1e6)))
# summary(human_no_Ramb_mod_deer)

human_no_Ramb_mod_noHosts <- glmer(
  cbind(ticks_human_no_Ramb, tick_abundance) ~ ticks_per_trap + (1 | plot_id), 
  data = sem_data, family = binomial(link = "logit"),
  glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=1e6)))
# summary(human_no_Ramb_mod_noHosts)

a <- MuMIn::AICc(human_no_Ramb_mod_allHosts, human_no_Ramb_mod_noRabbit, human_no_Ramb_mod_deer, human_no_Ramb_mod_noHosts) %>%
  arrange(AICc)
## minimal difference between model AIC
# coefTable(human_no_Ramb_mod_allHosts)
# coefTable(human_no_Ramb_mod_noRabbit)
# coefTable(human_no_Ramb_mod_deer)
# 
# plot_model(human_no_Ramb_mod_allHosts, type = "est", grid = T)/
#   plot_model(human_no_Ramb_mod_noRabbit, type = "est", grid = T)/
#   plot_model(human_no_Ramb_mod_deer, type = "est", grid = T)

p_noRamb_mods <- plot_models(human_no_Ramb_mod_noHosts, human_no_Ramb_mod_allHosts, human_no_Ramb_mod_deer, human_no_Ramb_mod_noRabbit)
# p_noRamb_mods$data
p_noRamb_mods + scale_y_continuous(limits = c(.9,1.1))

# mtab <- tab_model(human_no_Ramb_mod_allHosts, human_no_Ramb_mod_deer, human_no_Ramb_mod_noRabbit, human_no_Ramb_mod_noHosts)

## generally poor explanatory models, tiny bit of variation accounted for by 
## the plot-level random effect and tick abundance
## overall relatively few ticks infected 
sum(sem_data$ticks_human_no_Ramb)
sum(sem_data$tick_abundance)
100*(sum(sem_data$ticks_human_no_Ramb)/sum(sem_data$tick_abundance))

# simulateResiduals(human_no_Ramb_mod1) %>% plot()
# simulateResiduals(human_no_Ramb_mod2) %>% plot()
# simulateResiduals(human_no_Ramb_mod3) %>% plot()


## Human pathogens including R. amblyommatis
human_mod_allHosts <- glmer(
  cbind(ticks_human, tick_abundance) ~ ticks_per_trap + total_clusters1m + (1 | plot_id), 
  data = sem_data, family = binomial(link = "logit"),
  glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=1e6)))
# summary(human_mod_allHosts)

human_mod_noRabbit <- glmer(
  cbind(ticks_human, tick_abundance) ~ ticks_per_trap + deer_other + (1 | plot_id), 
  data = sem_data, family = binomial,
  glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=1e6)))
# summary(human_mod_noRabbit)

human_mod_deer <- glmer(
  cbind(ticks_human, tick_abundance) ~ ticks_per_trap + deer + (1 | plot_id), 
  data = sem_data, family = binomial,
  glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=1e6)))
# summary(human_mod_deer)

human_mod_noHosts <- glmer(
  cbind(ticks_human, tick_abundance) ~ ticks_per_trap + (1 | plot_id), 
  data = sem_data, family = binomial,
  glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=1e6)))
# summary(human_mod_noHosts)

b <- MuMIn::AICc(human_mod_allHosts, human_mod_noRabbit, human_mod_deer, human_mod_noHosts) %>% 
  arrange(AICc)# essentially same information
b
# simulateResiduals(human_mod1) %>% plot()
# simulateResiduals(human_mod2) %>% plot()
# simulateResiduals(human_mod3) %>% plot()
p_hum_mods <- plot_models(human_mod_allHosts, human_mod_noRabbit, human_mod_deer, human_mod_noHosts)
p_hum_mods$data <- p_hum_mods$data %>% 
  mutate(group = case_when(
    group == "cbind(ticks_human,\ntick_abundance).1" ~ "All hosts",
    group == "cbind(ticks_human,\ntick_abundance).2" ~ "No rabbit",
    group == "cbind(ticks_human,\ntick_abundance).3" ~ "Deer only",
    group == "cbind(ticks_human,\ntick_abundance).4" ~ "No hosts"
  ),
  term = case_when(
    term == "total_clusters1m" ~ "All hosts",
    term == "deer_other" ~ "No rabbit",
    term == "deer" ~ "Deer only",
    term == "ticks_per_trap" ~ "Tick abundance"
  )
  )
  
p_hum_mods + scale_y_continuous(limits = c(.95,1.05)) +
  # theme_bw() +
  theme(legend.title = element_blank())

tab_model(human_mod_allHosts, human_mod_noRabbit, human_mod_deer, human_mod_noHosts)
sum(sem_data$ticks_human)
sum(sem_data$tick_abundance)
100*(sum(sem_data$ticks_human)/sum(sem_data$tick_abundance))

## Human (incl. R. amblyommatis) and animal pathogens
human_animal_mod_allHosts <- glmer(
  cbind(ticks_human_animal, tick_abundance) ~ ticks_per_trap + total_clusters1m + (1 | inst_name/plot_id), 
  data = sem_data, family = binomial(link = "logit"),
  glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=1e6)))
# plot(simulateResiduals(human_animal_mod_allHosts))
# summary(human_animal_mod_allHosts)

human_animal_mod_noRabbit <- glmer(
  cbind(ticks_human_animal, tick_abundance) ~ ticks_per_trap + deer_other + (1 | inst_name/plot_id), 
  data = sem_data, family = binomial(link = "logit"),
  glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=1e6)))
# plot(simulateResiduals(human_animal_mod_noRabbit))
# summary(human_animal_mod_noRabbit)

human_animal_mod_deer <- glmer(
  cbind(ticks_human_animal, tick_abundance) ~ ticks_per_trap + deer + (1 | inst_name/plot_id), 
  data = sem_data, family = binomial(link = "logit"),
  glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=1e6)))
# plot(simulateResiduals(human_animal_mod_deer))
# summary(human_animal_mod_deer)

human_animal_mod_noHosts <- glmer(
  cbind(ticks_human_animal, tick_abundance) ~ ticks_per_trap + (1 | inst_name/plot_id), 
  data = sem_data, family = binomial(link = "logit"),
  glmerControl(optimizer = "Nelder_Mead",
               optCtrl=list(maxfun=1e6)))
# summary(human_animal_mod_noHosts)
# plot(simulateResiduals(human_animal_mod_noHosts))
c <- AIC(human_animal_mod_allHosts, 
         human_animal_mod_noRabbit, 
         human_animal_mod_deer, 
         human_animal_mod_noHosts)
## Very similar AIC values
# tab_model(human_animal_mod_allHosts, 
#     human_animal_mod_noRabbit, 
#     human_animal_mod_deer, 
#     human_animal_mod_noHosts)
