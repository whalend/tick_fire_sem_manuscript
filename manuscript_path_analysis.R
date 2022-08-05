## Path analysis
## standalone script

## I'm thinking of being liberal about transforming variables for the path analysis
## Then checking the direction of estimates against more generalized models.

library(tidyverse)
library(sjPlot)
library(patchwork)
library(MuMIn)
library(DHARMa)
library(lme4)
library(piecewiseSEM)

sem_data <- read_csv("data/path_analysis_data.csv")

## Ultimate terminal node is prevalence of infected ticks, which we believe is 
## determined primarily by tick density but may also be affected by host density.
## We have three ways to calculate prevalence of ticks infected with human pathogens:
## 1. Only human pathogens, excluding R. amblyommatis
## 2. Human pathogens, including R. amblyommatis
## 3. Human and animal pathogens, including R. amblyommatis
### While this would be ideal it doesn't play nicely overall with the current
### implementation of the path modeling, so it is now separated to:
### 'Rscripts/model_infected_tic_prevalence.R'

## We also have multiple groupings for relative host abundance
## 1. All host dung clusters
## 2. Deer and other except cottontail
## 3. Only deer

sem_data %>% 
  select(deer:total_clusters1m) %>% 
  colSums( ) %>% 
  knitr::kable()

sem_data <- sem_data %>%
  mutate(biomass_log = log(avg_dry_standing_gm2),
         d_since_fire_log = log(d_since_fire),
         logit_litter = car::logit(avg_pct_litter)
  )

## Make a rounded-up ticks per trap variable for Poisson and NegBin
## distribution compatibility, i.e. 0< x <1 == 1, and all other rounded to next
## larger integer
sem_data$tpt <- ceiling(sem_data$ticks_per_trap)

## Tick abundance models ####

### This "probability of ticks on traps" model doesn't make 
### much intuitive sense
# tick_mod_ah <- glmer(cbind(tick_abundance, trap_effort) ~ d_since_fire + avg_pct_litter + avg_litter_depth_all + avg_canopy_cover + avg_dry_standing_gm2 + total_clusters1m + avg_1yr_vp..Pa. + (1|inst_name/plot_id),
#   data = sem_data, family = binomial(link = "logit"),
#   na.action = "na.fail",
#   control=glmerControl(
#     optimizer = "Nelder_Mead",
#     optCtrl=list(maxfun=1e6)))
# summary(tick_mod_ah)
# simulateResiduals(tick_mod_ah) %>% plot()
# summary(effectsize::standardize(tick_mod_ah))

## This model would assume equal sampling effort - uses only 
## tick abundance, not controlling for trapping effort.
tick_abundance_nb <- glmer.nb(
  tick_abundance ~ d_since_fire_log + logit_litter + avg_litter_depth_all + avg_canopy_cover + biomass_log + total_clusters1m + avg_1yr_vp..Pa. + (1|inst_name/plot_id),
  data = sem_data, na.action = "na.fail",
  control=glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=1e6))
  )
plot(tick_abundance_nb)
simulateResiduals(tick_abundance_nb) %>% plot()
summary(tick_abundance_nb)
coefs(tick_abundance_nb)

## This model uses the rounded ticks per trap
## Outcome for relationship with time since fire is very different
tpt_nb <- glmer.nb(
  tpt ~ d_since_fire_log + logit_litter + avg_litter_depth_all + avg_canopy_cover + biomass_log + total_clusters1m + avg_1yr_vp..Pa. + (1|inst_name/plot_id),
  data = sem_data, na.action = "na.fail",
  control=glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=1e6))
  )
plot(tpt_nb)
simulateResiduals(tpt_nb) %>% plot()
effectsize::standardize(tpt_nb) %>% summary()

## unrounded ticks per trap, violates integer value assumption
## compare to rounded integer model
ticks_per_trap_nb <- glmer.nb(
  ticks_per_trap ~ d_since_fire_log + logit_litter + avg_litter_depth_all + avg_canopy_cover + biomass_log + total_clusters1m + avg_1yr_vp..Pa. + (1|inst_name/plot_id),
  data = sem_data, na.action = "na.fail",
  control=glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=1e6))
  )
plot(ticks_per_trap_nb)
simulateResiduals(ticks_per_trap_nb) %>% plot()
summary(ticks_per_trap_nb)
effectsize::standardize(ticks_per_trap_nb) %>% summary()


## Gamma ticks per trap model ####
sem_data$ticks_per_trap001 <- sem_data$ticks_per_trap + .001
ticks_per_trap_gamma <- glmer(
  ticks_per_trap001 ~ d_since_fire_log + logit_litter + avg_litter_depth_all + avg_canopy_cover + biomass_log + total_clusters1m + avg_1yr_vp..Pa. + (1|inst_name/plot_id),
  data = sem_data, na.action = "na.fail",
  family = Gamma(link = "log"),
  control=glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=1e6))
)
plot(ticks_per_trap_gamma)
simulateResiduals(ticks_per_trap_gamma) %>% plot()
summary(ticks_per_trap_gamma)
effectsize::standardize(ticks_per_trap_gamma) %>% summary()
# effectsize::standardize_parameters(ticks_per_trap_gamma)
# effectsize::standardize_posteriors(ticks_per_trap_gamma)
# coefs(ticks_per_trap_gamma)

plot_models(
  effectsize::standardize(tpt_nb), 
  effectsize::standardize(ticks_per_trap_nb), 
  effectsize::standardize(tick_abundance_nb),
  show.values = T)
tab_model(tpt_nb, ticks_per_trap_nb, tick_abundance_nb)

## Poisson ticks per trap model ####
tpt_pois <- glmer(
  tpt ~ d_since_fire_log + logit_litter + avg_litter_depth_all + avg_canopy_cover + biomass_log + total_clusters1m + avg_1yr_vp..Pa. + (1|inst_name/plot_id),
  data = sem_data, na.action = "na.fail",
  family = poisson(link = "log"),
  control=glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=1e6))
  )
plot(tpt_pois)
simulateResiduals(tpt_pois) %>% plot()
summary(tpt_pois)
effectsize::standardize(tpt_pois) %>% summary()
coefs(tpt_pois)
# plot_model(tpt_pois, show.values = T)
# plot_model(tpt_pois, type = "pred", show.data = T, grid = T)

tpt_pois2 <- glmer(
  tpt ~ d_since_fire + avg_pct_litter + avg_litter_depth_all + avg_canopy_cover + avg_dry_standing_gm2 + total_clusters1m + avg_1yr_vp..Pa. + (1|inst_name/plot_id), 
  data = sem_data, na.action = "na.fail",
  family = poisson(link = "log"),
  control=glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=1e6))
)
simulateResiduals(tpt_pois2) %>% plot()
effectsize::standardize(tpt_pois2) %>% summary()
coefs(tpt_pois2)
visreg::visreg(tpt_pois2, "d_since_fire")
plot_model(tpt_pois2, type = "pred", grid = T, auto.label = FALSE)

## Poisson tick abundance model ####
tick_abundance_pois <- glmer(
  tick_abundance ~ d_since_fire_log + logit_litter + avg_litter_depth_all + avg_canopy_cover + biomass_log + total_clusters1m + avg_1yr_vp..Pa. + (1|inst_name/plot_id),
  data = sem_data, na.action = "na.fail",
  family = poisson(link = "log"),
  control=glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=1e6))
)
simulateResiduals(tick_abundance_pois) %>% plot()
effectsize::standardize(tick_abundance_pois) %>% summary()
effectsize::standardize_parameters(tick_abundance_pois)

## Considered composite litter variable but gets really complicated because
## the litter variables are endogenous predictors (response & predictor)
# beta_cover <- coefs(tpt_pois)[2,3]
# beta_depth <- coefs(tpt_pois)[3,3]
# sem_data <- sem_data %>% 
#   mutate(litter_composite = beta_cover*logit_litter + beta_depth*avg_litter_depth_all)
# 
# tpt_pois2 <- glmer(tpt ~ d_since_fire_log + litter_composite + avg_canopy_cover + biomass_log + total_clusters1m + avg_1yr_vp..Pa. + (1|inst_name/plot_id),
#                   data = sem_data,
#                   family = poisson(link = "log"),
#                   na.action = "na.fail",
#                   control=glmerControl(
#                     optimizer = "Nelder_Mead",
#                     optCtrl=list(maxfun=1e6)))
# summary(tpt_pois2)
# effectsize::standardize(tpt_pois2) %>% summary()
# coefs(tpt_pois2)
# plot(tpt_pois2)
# simulateResiduals(tpt_pois2) %>% plot()


## Deer-only ticks model ####
# ticks_per_trap_deer_nb <- glmer.nb(
#   ticks_per_trap ~ d_since_fire_log + logit_litter + avg_litter_depth_all + avg_canopy_cover + biomass_log + deer + avg_1yr_vp..Pa. + (1|inst_name/plot_id),
#   data = sem_data, na.action = "na.fail",
#   control=glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=1e6))
#   )
# simulateResiduals(tpt_deer_nb) %>% plot()
# summary(tpt_deer_nb)
# effectsize::standardize(tpt_deer_nb) %>% summary()

tpt_deer_nb <- glmer.nb(
  tpt ~ d_since_fire_log + logit_litter + avg_litter_depth_all + avg_canopy_cover + biomass_log + deer + avg_1yr_vp..Pa. + (1|inst_name/plot_id),
  data = sem_data, na.action = "na.fail",
  control=glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=1e6))
)

tpt_deer_pois <- glmer(
  tpt ~ d_since_fire_log + logit_litter + avg_litter_depth_all + avg_canopy_cover + biomass_log + deer + avg_1yr_vp..Pa. + (1|inst_name/plot_id),
  data = sem_data, na.action = "na.fail",
  family = poisson(link = "log"),
  control=glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=1e6))
  )
# simulateResiduals(tpt_deer_pois) %>% plot()
# summary(tpt_deer_pois)
summary(tpt_deer_nb)$coefficients
summary(tpt_deer_pois)$coefficients
## Similar enough, though in the model summary the p-values are a bit smaller
## in the poisson model


## Models excluding rabbit host data ####
# ticks_per_trap_noRabb_nb <- glmer.nb(
#   ticks_per_trap ~ d_since_fire_log + logit_litter + avg_litter_depth_all + avg_canopy_cover + biomass_log + deer_other + avg_1yr_vp..Pa. + (1|inst_name/plot_id),
#   data = sem_data, na.action = "na.fail",
#   control=glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=1e6))
#   )
# simulateResiduals(ticks_per_trap_noRabb_nb) %>% plot()
# effectsize::standardize(ticks_per_trap_noRabb_nb) %>% 
  # simulateResiduals() %>% 
  # plot()
# summary(ticks_per_trap_noRabb_nb)
# effectsize::standardize(ticks_per_trap_noRabb_nb) %>% summary()

tpt_noRabb_nb <- glmer.nb(
  tpt ~ d_since_fire_log + logit_litter + avg_litter_depth_all + avg_canopy_cover + biomass_log + deer_other + avg_1yr_vp..Pa. + (1|inst_name/plot_id),
  data = sem_data, na.action = "na.fail",
  control=glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=1e6))
)

tpt_noRabb_pois <- glmer(
  tpt ~ d_since_fire_log + logit_litter + avg_litter_depth_all + avg_canopy_cover + biomass_log + deer_other + avg_1yr_vp..Pa. + (1|inst_name/plot_id),
  data = sem_data, na.action = "na.fail",
  family = poisson(link = "log"),
  control=glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=1e6))
  )

# summary(tpt_noRabb_pois)
# effectsize::standardize(tpt_noRabb_pois) %>% summary()
# effectsize::standardize(tpt_noRabb_pois) %>% 
#   simulateResiduals() %>% 
#   plot()

## Rabbit only ticks model ####
tpt_Rabb_pois <- glmer(
  tpt ~ d_since_fire_log + logit_litter + avg_litter_depth_all + avg_canopy_cover + biomass_log + cottontail + avg_1yr_vp..Pa. + (1|inst_name/plot_id),
  data = sem_data, 
  family = poisson(link = "log"),
  na.action = "na.fail",
  control=glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=1e6))
  )
# summary(tpt_Rabb_pois)


## No hosts ticks models ####
# ticks_per_trap_noHosts_nb <- glmer.nb(
#   ticks_per_trap ~ d_since_fire_log + logit_litter + avg_litter_depth_all + avg_canopy_cover + biomass_log + avg_1yr_vp..Pa. + (1|inst_name/plot_id),
#   data = sem_data, na.action = "na.fail",
#   control=glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=1e6))
#   )
# simulateResiduals(ticks_per_trap_noHosts_nb) %>% plot()
# summary(ticks_per_trap_noHosts_nb)
# effectsize::standardize(ticks_per_trap_noHosts_nb) %>% summary()


tpt_noHosts_nb <- glmer.nb(
  tpt ~ d_since_fire_log + logit_litter + avg_litter_depth_all + avg_canopy_cover + biomass_log + avg_1yr_vp..Pa. + (1|inst_name/plot_id),
  data = sem_data, 
  na.action = "na.fail",
  control=glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=1e6))
)

tpt_noHosts_pois <- glmer(
  tpt ~ d_since_fire_log + logit_litter + avg_litter_depth_all + avg_canopy_cover + biomass_log + avg_1yr_vp..Pa. + (1|inst_name/plot_id),
  data = sem_data, 
  family = poisson(link = "log"),
  na.action = "na.fail",
  control=glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=1e6))
  )
# simulateResiduals(tpt_noHosts_pois) %>% plot()
# summary(tick_mod_noHosts_gamma)
# effectsize::standardize(tpt_noHosts_pois) %>% summary()


# AICc(tpt_nb, tpt_noRabb_nb, tpt_deer_nb, tpt_noHosts_nb) %>%
  # arrange(AICc)
# AICc(tpt_pois, tpt_noRabb_pois, tpt_deer_pois, tpt_noHosts_pois) %>%
  # arrange(AICc)


## Host models: all, deer only, no rabbit ####
host_mod_all <- glmer(
  total_clusters1m ~ d_since_fire_log + avg_canopy_cover + logit_litter + avg_litter_depth_all + (1 | inst_name/plot_id), 
  data = sem_data, family = poisson(link = "log"),
  control = glmerControl(optimizer="Nelder_Mead", optCtrl=list(maxfun=1e6))
  )
coefs(host_mod_all)
simulateResiduals(host_mod_all) %>% plot()
car::vif(host_mod_all)

host_mod_all2 <- glmer(
  total_clusters1m ~ d_since_fire_log + avg_canopy_cover + logit_litter + avg_litter_depth_all + cv_30yr_fire_days + (1 | inst_name/plot_id),
  data = sem_data, family = poisson(link = "log"),
  control = glmerControl(optimizer="Nelder_Mead", optCtrl=list(maxfun=1e6))
  )
car::vif(host_mod_all2)

host_mod_all_nb <- glmer.nb(
  total_clusters1m ~ d_since_fire_log + avg_canopy_cover + biomass_log + logit_litter + avg_litter_depth_all + (1 | inst_name/plot_id),
  data = sem_data, 
  control = glmerControl(optimizer="Nelder_Mead", optCtrl=list(maxfun=1e6))
  )
simulateResiduals(host_mod_all_nb) %>% plot()


host_mod_deer <- glmer(
  deer ~ d_since_fire_log + avg_canopy_cover + logit_litter + avg_litter_depth_all + (1|inst_name/plot_id),
  data = sem_data, family = poisson(link = "log"),
  control = glmerControl(optimizer="Nelder_Mead", optCtrl=list(maxfun=1e6))
  )
coefs(host_mod_deer)
simulateResiduals(host_mod_deer) %>% plot()

host_mod_noRabbit <- glmer(
  deer_other ~ d_since_fire_log + avg_canopy_cover + logit_litter + avg_litter_depth_all + (1|inst_name/plot_id),
  data = sem_data, family = poisson(link = "log"),
  control = glmerControl(optimizer="Nelder_Mead", optCtrl=list(maxfun=1e6))
  )
coefs(host_mod_noRabbit)
simulateResiduals(host_mod_noRabbit) %>% plot()

## Litter/understory models ####
litter_cover_mod <- lmer(logit_litter ~ d_since_fire_log + avg_canopy_cover + (1|inst_name/plot_id), data = sem_data)
coefs(litter_cover_mod)

litter_depth_mod <- lmer(avg_litter_depth_all ~ d_since_fire_log + avg_canopy_cover + (1|inst_name/plot_id), data = sem_data)
coefs(litter_depth_mod)


standing_biomass_mod <- lmer(biomass_log ~ d_since_fire_log + avg_canopy_cover + (1|inst_name/plot_id), data = sem_data)
# plot(standing_biomass_mod)
simulateResiduals(standing_biomass_mod) %>% plot()
summary(standing_biomass_mod)


## Tree (overstory) canopy cover ####
# pine_mod <- lmer(logit_pct_pine ~ fri15yr_log + avg_canopy_cover +(1|inst_name/plot_id), 
#                  data = sem_data)
canopy_mod <- lmer(avg_canopy_cover ~ d_since_fire_log + fri15yr + (1|inst_name/plot_id), 
                   data = sem_data)
canopy_mod2 <- lmer(avg_canopy_cover ~ d_since_fire_log + fri15yr + avg_1yr_vp..Pa. + (1|inst_name/plot_id), 
                   data = sem_data)
# simulateResiduals(canopy_mod) %>% plot()
# summary(canopy_mod)


## Fire regime models ####
time_since_fire_mod_pred <- glmer.nb(d_since_fire ~ fri15yr + (1|inst_name/plot_id),
         data = sem_data)
simulateResiduals(time_since_fire_mod_pred) %>% plot()
plot(time_since_fire_mod_pred)
plot_model(time_since_fire_mod_pred, show.values = T)
plot_model(time_since_fire_mod_pred, type = "pred", terms = "fri15yr [all]", show.data = T)

time_since_fire_mod <- lmer(d_since_fire_log ~ fri15yr + (1|inst_name/plot_id), data = sem_data)
simulateResiduals(time_since_fire_mod) %>% plot()
plot(time_since_fire_mod)
summary(time_since_fire_mod)
plot_model(time_since_fire_mod, show.values = T)

fri_mod <- lmer(fri15yr ~ cv_30yr_fire_days + (1 | inst_name/plot_id), 
                data = sem_data)
simulateResiduals(fri_mod) %>% plot()
plot(fri_mod)

sem_data$fri15yr_int <- ceiling(sem_data$fri15yr)
fri_mod_pred <- glmer(fri15yr_int ~ cv_30yr_fire_days + (1 | inst_name/plot_id), 
                      data = sem_data,
                      family = poisson(link = "log"),
                      control = glmerControl(
                        optimizer="Nelder_Mead", optCtrl=list(maxfun=1e6))
                      )
simulateResiduals(fri_mod_pred) %>% plot()
plot(fri_mod_pred)
# plot_model(fri_mod, type = "pred", show.data = T)
# plot_model(fri_mod_pred, type = "pred", show.data = T)


## Weather ~ Climate model
avg_1yr_vp_mod <- lmer(avg_1yr_vp..Pa. ~ cv_30yr_fire_days + (1|inst_name), data = sem_data)
# simulateResiduals(avg_1yr_vp_mod) %>% plot()
# summary(avg_1yr_vp_mod)
avg_1yr_vp_mod2 <- lmer(avg_1yr_vp..Pa. ~ avg_canopy_cover + cv_30yr_fire_days + (1|inst_name), data = sem_data)
# summary(avg_1yr_vp_mod2)


## Path Models ####
# Convert to data frame class
sem_df <- as.data.frame(sem_data)

## All hosts
all_hosts_psem <- psem(
  # m2,
  tpt_pois,
  host_mod_all,
  litter_cover_mod, litter_depth_mod, canopy_mod,
  standing_biomass_mod,
  time_since_fire_mod, fri_mod,
  avg_1yr_vp_mod,
  
  data = sem_df
)
# Include human-only pathogen prevalence in the model - have to find this again
# all_hostsP_psem <- psem(
#   m2,
#   tpt_pois,
#   host_mod_all,
#   litter_cover_mod, litter_depth_mod, canopy_mod,
#   standing_biomass_mod,
#   time_since_fire_mod, fri_mod,
#   avg_1yr_vp_mod,
#   
#   data = sem_df
# )
ah_psem_summary <- summary(all_hosts_psem, .progressBar = FALSE)
ah_psem_summary
# ahP_psem_summary <- summary(all_hostsP_psem, .progressBar = FALSE)
# ahP_psem_summary
# ahP_psem_summary$AIC
ah_psem_summary$IC

tpt_mod_pred <- plot_model(tpt_pois, type = "pred", show.data = F)
tpt_mod_pred_grid <- plot_model(
  tpt_pois, type = "pred", show.data = F, grid = T
  )
tpt_mod_pred_grid + xlab(" ") + ylab("Ticks per trap")

tpt_mod_pred_grid$data <- tpt_mod_pred_grid$data %>% 
  filter(!group_col %in% c("total_clusters1m", "avg_canopy_cover","d_since_fire_log")) %>% 
  mutate(group = case_when(
    group == "logit_litter" ~ "Litter cover (logit)",
    group == "avg_litter_depth_all" ~ "Litter depth (cm)",
    group == "biomass_log" ~ "Standing biomass (log)",
    group == "avg_1yr_vp..Pa." ~ "Vapor pressure (Pa)"
    )
    )
(tpt_mod_grid_fig <- tpt_mod_pred_grid + 
  xlab(" ") + ylab("Predicted ticks per trap") + 
  theme_bw()
)
ggsave("figures/tpt_allhosts_marginal.png",tpt_mod_grid_fig, width = 5, height = 5, dpi = 600)

tpt_lcov <- tpt_mod_pred$logit_litter +
  geom_point(data = sem_df, aes(logit_litter, tpt)) +
  xlab("Litter cover (logit)") + ylab("Predicted tick count") + ggtitle(" ") +
  theme_classic()
tpt_ldep <- tpt_mod_pred$avg_litter_depth_all +
  geom_point(data = sem_df, aes(avg_litter_depth_all, tpt)) +
  xlab("Litter depth (cm)") + ylab("Predicted tick count") + ggtitle(" ") +
  theme_classic()
tpt_ccov <- tpt_mod_pred$avg_canopy_cover +
  geom_point(data = sem_df, aes(avg_canopy_cover, tpt)) +
  xlab("Canopy cover (%)") + ylab("Predicted tick count") + ggtitle(" ") +
  theme_classic()
tpt_vprs <- tpt_mod_pred$avg_1yr_vp..Pa. +
  geom_point(data = sem_df, aes(avg_1yr_vp..Pa., tpt)) +
  xlab("Vapor pressure (Pa)") + ylab("Predicted tick count") + ggtitle(" ") +
  theme_classic()
tpt_lcov + tpt_ldep + ylab(" ") + tpt_ccov + tpt_vprs + ylab(" ") &
  plot_annotation(tag_levels = 'A')
ggsave("figures/ticks_per_trap_allhosts_significant_marginal_effects.png", width = 5, height = 5, dpi = 600)


host_mod_pred <- plot_model(host_mod_all, type = "pred", show.data = F)
host_mod_pred_grid <- plot_model(
  host_mod_all, type = "pred", show.data = F, grid = T
)

host_mod_pred_grid$data <- host_mod_pred_grid$data %>% 
  mutate(group = case_when(
    group == "logit_litter" ~ "Litter cover (logit)",
    group == "avg_litter_depth_all" ~ "Litter depth (cm)",
    group == "avg_canopy_cover" ~ "Canopy cover (%)",
    group == "d_since_fire_log" ~ "Days since fire (log)"
  )
  )
(host_mod_pred_grid <- host_mod_pred_grid + 
  xlab(" ") + ylab("Predicted count of dung clusters") + 
  theme_bw()
)
ggsave("figures/hosts_allhosts_marginal.png", host_mod_pred_grid, width = 5, height = 5, dpi = 300)

host_lcov <- host_mod_pred$logit_litter +
  xlab("Litter cover (logit)") + ylab("Predicted count of dung clusters") +
  ggtitle(" ") +
  theme_classic()
host_ldepth <- host_mod_pred$avg_litter_depth_all +
  xlab("Litter depth (cm)") + ylab("Predicted count of dung clusters") + 
  ggtitle(" ") +
  theme_classic()
host_ccov <- host_mod_pred$avg_canopy_cover +
  xlab("Canopy cover (%)") + ylab("Predicted count of dung clusters") + 
  ggtitle(" ") +
  theme_classic()
host_dfire <- host_mod_pred$d_since_fire_log +
  xlab("Days since fire (log)") + ylab("Predicted count of dung clusters") + 
  ggtitle(" ") +
  theme_classic()
host_lcov + host_ldepth + ylab(" ") + host_ccov + host_dfire + ylab(" ") &
  plot_annotation(tag_levels = 'A')
ggsave("figures/dung_clusters_significant_marginal_effects.png", width = 5, height = 5, dpi = 600)


lcov_mod_pred <- plot_model(litter_cover_mod, type = "pred", show.data = F)
lcov_mod_pred_grid <- plot_model(
  litter_cover_mod, type = "pred", show.data = F, grid = T
)
lcov_mod_pred_grid$data <- lcov_mod_pred_grid$data %>% 
  mutate(group = case_when(
    group == "avg_canopy_cover" ~ "Canopy cover (%)",
    group == "d_since_fire_log" ~ "Days since fire (log)"
  )
  )
(lcov_mod_pred_grid <- lcov_mod_pred_grid + 
  xlab(" ") + ylab("Predicted litter cover (logit)") + 
  theme_bw()
)
ggsave("figures/litter_cover_allhosts_marginal.png", lcov_mod_pred_grid, width = 5, height = 5, dpi = 300)

lcov_ccov <- lcov_mod_pred$avg_canopy_cover +
  xlab("Canopy cover (%)") + ylab("Predicted litter cover (logit)") +
  ggtitle(" ") +
  theme_classic()
lcov_dfire <- lcov_mod_pred$d_since_fire_log +
  xlab("Days since fire (log)") + ylab("Predicted litter cover (logit)") + 
  ggtitle(" ") +
  theme_classic()
lcov_ccov + lcov_dfire + ylab(" ") &
  plot_annotation(tag_levels = 'A')
ggsave("figures/litter_cover_significant_marginal_effects.png", width = 5, height = 5, dpi = 600)


ldep_mod_pred <- plot_model(litter_depth_mod, type = "pred", show.data = F)
ldep_mod_pred_grid <- plot_model(
  litter_depth_mod, type = "pred", show.data = F, grid = T
)
ldep_mod_pred_grid$data <- ldep_mod_pred_grid$data %>%
  filter(group == "d_since_fire_log") %>% 
  mutate(group = "Days since fire (log)")
ldep_mod_pred_grid <- ldep_mod_pred_grid + 
  xlab(" ") + ylab("Predicted litter depth (cm)") + 
  theme_bw()

ldep_ccov <- ldep_mod_pred$avg_canopy_cover +
  xlab("Canopy cover (%)") + ylab("Predicted litter depth (cm)") +
  ggtitle(" ") +
  theme_classic()
ldep_dfire <- ldep_mod_pred$d_since_fire_log +
  xlab("Days since fire (log)") + ylab("Predicted litter depth (cm)") + 
  ggtitle(" ") +
  theme_classic()
ldep_ccov + ldep_dfire + ylab(" ") &
  plot_annotation(tag_levels = 'A')
ggsave("figures/litter_depth_marginal_effects.png", width = 5, height = 5, dpi = 600)


sbio_mod_pred <- plot_model(standing_biomass_mod, type = "pred")
sbio_mod_pred_grid <- plot_model(
  standing_biomass_mod, type = "pred", grid = T
  )
sbio_mod_pred_grid$data <- sbio_mod_pred_grid$data %>% 
  filter(group == "avg_canopy_cover") %>% 
  mutate(group = "Canopy cover (%)")
(sbio_mod_pred_grid <- sbio_mod_pred_grid +
  xlab(" ") + ylab("Predicted understory biomass (ln grams)") +
  theme_bw()
)
lcov_mod_pred_grid / (ldep_mod_pred_grid + sbio_mod_pred_grid)
ggsave("figures/veg_vars_allhosts_marginal.png", width = 5, height = 5, dpi = 300)

sbio_ccov <- sbio_mod_pred$avg_canopy_cover +
  xlab("Canopy cover (%)") + ylab("Predicted understory biomass (ln grams)") +
  ggtitle(" ") +
  theme_classic()
sbio_dfire <- sbio_mod_pred$d_since_fire_log +
  xlab("Days since fire (log)") + ylab("Predicted understory biomass (ln grams)") + 
  ggtitle(" ") +
  theme_classic()
sbio_ccov + sbio_dfire +ylab(" ") &
  plot_annotation(tag_levels = 'A')
ggsave("figures/understory_biomass_marginal_effects.png", width = 5, height = 5, dpi = 600)

canp_mod_pred <- plot_model(canopy_mod, type = "pred")
canp_mod_pred_grid <- plot_model(
  canopy_mod, type = "pred", grid = T
  )
canp_mod_pred_grid$data <- canp_mod_pred_grid$data %>% 
  filter(group == "d_since_fire_log") %>% 
  mutate(group = "Days since fire (log)")
(canp_mod_pred_grid <- canp_mod_pred_grid +
  xlab(" ") + ylab ("Predicted canopy cover (%)") +
  theme_bw()
)

canp_fri <- canp_mod_pred$fri15yr +
  xlab("15-yr fire return interval") + ylab("Predicted canopy cover (%)") +
  ggtitle(" ") +
  theme_classic()
canp_dfire <- canp_mod_pred$d_since_fire_log +
  xlab("Days since fire (log)") + ylab("Predicted canopy cover (%)") + 
  ggtitle(" ") +
  theme_classic()
canp_fri + canp_dfire + ylab(" ") &
  plot_annotation(tag_levels = 'A')
ggsave("figures/canopy_cover_marginal_effects.png", width = 5, height = 5, dpi = 600)


fire_mod_pred <- plot_model(time_since_fire_mod, type = "pred")
fire_mod_pred_grid <- plot_model(
  time_since_fire_mod, type = "pred", grid = T
)
fire_mod_pred_grid$data <- fire_mod_pred_grid$data %>% 
  mutate(group = "Fire return interval (years)")
fire_mod_pred_grid <- fire_mod_pred_grid +
  xlab(" ") + ylab("Predicted log-days since fire") +
  theme_bw()
fire_fri <- fire_mod_pred$fri15yr +
  xlab("15-yr fire return interval") + ylab("Predicted days since fire (ln)") +
  ggtitle(" ") +
  theme_classic()
fire_fri
ggsave("figures/days_since_fire_marginal_effects.png", width = 5, height = 5, dpi = 600)


fri_mod_pred <- plot_model(fri_mod, type = "pred")
fri_mod_pred_grid <- plot_model(
  fri_mod, type = "pred", grid = T
)
fri_mod_pred_grid$data <- fri_mod_pred_grid$data %>% 
  mutate(group = "30-y Average CV of Fire Days")
fri_mod_pred_grid <- fri_mod_pred_grid +
  xlab(" ") + ylab("Predicted fire return interval (years)") +
  theme_bw()
fri_cv30yr <- fri_mod_pred$cv_30yr_fire_days +
  xlab("30-y Average CV of Fire Days") + ylab("Predicted 15-y FRI (years)") +
  ggtitle(" ") +
  theme_classic()
fri_cv30yr
ggsave("figures/fire_return_interval_fire_marginal_effects.png", width = 5, height = 5, dpi = 600)

canp_mod_pred_grid + fire_mod_pred_grid + fri_mod_pred_grid
ggsave("figures/can_fire_vars_allhosts_marginal.png", width = 7, height = 3, dpi = 300)

## Treat all correlations as spurious/having a common unmeasured driver. ####
all_hosts_psem2 <- psem(
  tpt_pois,
  host_mod_all,
  litter_cover_mod, 
  litter_depth_mod,
  # litter_depth_mod2 <- lmer(avg_litter_depth_all ~ d_since_fire_log + (1 | inst_name/plot_id), data = sem_df), 
  canopy_mod,
  # canopy_mod2 <- lmer(avg_canopy_cover ~ d_since_fire_log + (1 | inst_name/plot_id), data = sem_df),
  standing_biomass_mod,
  # standing_biomass_mod2 <- lmer(biomass_log ~ avg_canopy_cover + (1 | inst_name/plot_id), data = sem_df),
  time_since_fire_mod, fri_mod,
  avg_1yr_vp_mod,
  
  avg_canopy_cover %~~% avg_1yr_vp..Pa., #should this mean something?
  total_clusters1m %~~% cv_30yr_fire_days, #strongly correlated, biologically nonsense
  total_clusters1m %~~% avg_1yr_vp..Pa., #strongly correlated, likely biologically nonsense
  ## I get this msg when I include the line above:
  # Warning message:
  # Check model convergence: log-likelihood estimates lead to negative Chi-squared!
  ## On the flip-side, the Fisher's C is acceptable. Path coefficient estimates
  ## should be the same either way. Check-in w/team on biological sense
  
  data = sem_df
)
ah_psem2_summary <- summary(all_hosts_psem2, .progressBar = FALSE)
ah_psem2_summary

# Compare models without and with correlations ####
ah_psem_summary$coefficients %>% 
  knitr::kable(digits = 3, caption = "Table xx. Coeffcient estimates all hosts path model")
ah_psem2_summary$coefficients %>% 
  knitr::kable(digits = 3, caption = "Table xx. Coeffcient estimates all hosts path model w/correlations for missing paths.")
ah_psem_summary$R2 %>% 
  knitr::kable(digits = 3, caption = "R2 all hosts path model")
ah_psem2_summary$R2 %>% 
  knitr::kable(digits = 3, caption = "R2 all hosts path model w/correlations for missing paths")
ah_psem_summary$IC
ah_psem2_summary$IC
## Same AIC
ah_psem_summary$Cstat
ah_psem2_summary$Cstat
## Change in C-statistic
# ah_psem_summary$ChiSq
# ah_psem2_summary$ChiSq
## No Chisq value because negative Chi-squared values during estimation
## stemming from the total_clusters1m %~~% avg_1yr_vp..Pa. correlation (??)


## No rabbit hosts pSEM ####
# noRabbit_psem <- psem(
#   tpt_noRabb_pois,
#   host_mod_noRabbit,
#   litter_cover_mod, litter_depth_mod, canopy_mod,
#   standing_biomass_mod,
#   time_since_fire_mod, fri_mod,
#   avg_1yr_vp_mod,
#   
#   # avg_canopy_cover %~~% avg_1yr_vp..Pa., #should this mean something?
#   # deer_other %~~% cv_30yr_fire_days, #strongly correlated, biologically nonsense
#   # deer_other %~~% avg_1yr_vp..Pa.,
#   ## Throws negative Chi-squared warning regardless
#   
#   data = sem_df
# )
# noRabbit_psem_summary <- summary(noRabbit_psem, .progressBar = FALSE)
# noRabbit_psem_summary


## Deer only pSEM ####
# deer_psem <- psem(
#   tpt_deer_pois,
#   host_mod_deer,
#   litter_cover_mod, litter_depth_mod, canopy_mod,
#   standing_biomass_mod,
#   time_since_fire_mod, fri_mod,
#   avg_1yr_vp_mod,
#   
#   # total_clusters1m %~~% cv_30yr_fire_days, #strongly correlated, biologically nonsense
#   # avg_canopy_cover %~~% avg_1yr_vp..Pa., #should this mean something?
#   
#   data = sem_df
# )
# deer_psem_summary <- summary(deer_psem, .progressBar = FALSE)
# # Warning message:
#   # Check model convergence: log-likelihood estimates lead to negative Chi-squared!
# deer_psem_summary


## No hosts pSEM ####
m <- effectsize::standardize(tpt_noHosts_pois)
noHost_psem <- psem(
  tpt_noHosts_pois,
  # m,
  litter_cover_mod, litter_depth_mod, canopy_mod,
  standing_biomass_mod,
  time_since_fire_mod, fri_mod,
  avg_1yr_vp_mod,
  
  data = sem_df
)
noHosts_psem_summary <- summary(noHost_psem, .progressBar = FALSE)
noHosts_psem_summary

rbind(ah_psem_summary$IC,
      # noRabbit_psem_summary$AIC,
      # deer_psem_summary$AIC,
      noHosts_psem_summary$IC) %>% 
  mutate(model = c("all hosts", 
                   # "no rabbit", "deer only", 
                   "no hosts")) %>% 
  knitr::kable(digits = 3)
## Model without hosts has lowest AIC

rbind(ah_psem_summary$Cstat,
      # noRabbit_psem_summary$Cstat,
      # deer_psem_summary$Cstat,
      noHosts_psem_summary$Cstat) %>% 
  mutate(model = c("all hosts", 
                   # "no rabbit", "deer only", 
                   "no hosts")) %>% 
  knitr::kable(digits = 3)

# rbind(ah_psem_summary$ChiSq,
      # noRabbit_psem_summary$ChiSq,
      # deer_psem_summary$ChiSq,
      # noHosts_psem_summary$ChiSq) %>% 
  # mutate(model = c("all hosts", 
                   # "no rabbit", "deer only", 
                   # "no hosts")) %>% 
  # knitr::kable(digits = 3)


# rbind(ah_psem_summary$coefficients[6,],
#       noRabbit_psem_summary$coefficients[6,],
#       deer_psem_summary$coefficients[6,]) %>% 
#   mutate(model = c("all hosts", 
#                    # "no rabbit", "deer only"
#                    )) %>% 
#   knitr::kable(digits = 3)


## Calculate indirect effects on tick abundance ####
## all hosts model
ah_psem_summary$coefficients
time_since_fire_i <- (.39*.52) + (.4*.22) + (.17*.4*.52) + (.17*-.37*.25)
canopy_cover_i <- (.4*.52) + (-.37*.25)
FRI_i <- (.58*.17*.4*.52) + (.58*.17*-.37*.25) + (.58*.39*.52) + (.58*.4*.22)
cv_fire_days_i <- (.1*.58*.17*.4*.52) + .1*(.58*.17*-.37*.25) + .1*(.58*.39*.52) + .1*(.58*.4*.22)

tpt_effects_allHosts <- tibble(
  variable = c("CV fire days", "FRI", "time since fire", "canopy cover", 
               "litter cover", "litter depth", "standing biomass", 
               "vapor pressure"),
  direct = c(0, 0, 0, 0, .52, .22, .25, .53),
  indirect = c(cv_fire_days_i, FRI_i, time_since_fire_i, canopy_cover_i, 0, 0, 0, 0),
  total = direct + indirect
)
tpt_effects_allHosts %>% 
  knitr::kable(digits = 3, caption = "All hosts model effects on tick abundance")

## No hosts model
noHosts_psem_summary$coefficients
time_since_fire_i <- (.39*.51) + (.4*.24) + (.17*.4*.51) + (.17*-.37*.27)
canopy_cover_i <- (.4*.51) + (-.37*.27)
FRI_i <- (.58*.17*.4*.51) + (.58*.17*-.37*.27) + (.58*.39*.51) + (.58*.4*.24)
cv_fire_days_i <- (.1*.58*.17*.4*.51) + .1*(.58*.17*-.37*.27) + .1*(.58*.39*.51) + .1*(.58*.4*.24)

tpt_effects_noHosts <- tibble(
  variable = c("CV fire days", "FRI", "time since fire", "canopy cover", 
               "litter cover", "litter depth", "standing biomass", 
               "vapor pressure"),
  direct = c(0, 0, 0, 0, .51, .24, .27, .51),
  indirect = c(cv_fire_days_i, FRI_i, time_since_fire_i, canopy_cover_i, 0, 0, 0, 0),
  total = direct + indirect
)
tpt_effects_noHosts %>% 
  knitr::kable(digits = 3, caption = "No hosts model effects on tick abundance")
