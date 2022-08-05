## Path analysis
## standalone script

## I'm thinking of being liberal about transforming variables for the path analysis
## Then checking the direction of estimates against models with fewer/no transformations

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
### While this would be ideal it doesn't play nicely overall with the 
### path modeling, so it is now separated to:
### 'Rscripts/model_infected_tick_prevalence.R'

## I also made multiple groupings for relative host abundance because rabbit
## makes up a large proportion of the total dung counts
## 1. All host dung clusters
## 2. Deer and other except cottontail
## 3. Only deer
## 4. Only rabbit

sem_data %>% 
  select(deer:total_clusters1m) %>% 
  colSums( ) %>% 
  knitr::kable(col.names = "Dung count")

sem_mod_data <- sem_data %>%
  mutate(biomass_log = log(avg_dry_standing_gm2),
         biomass_kgm2 = avg_dry_standing_gm2/1000,
         d_since_fire_log = log(d_since_fire),
         logit_litter = ifelse(avg_pct_litter==100, 99.99, avg_pct_litter),
         logit_litter = car::logit(logit_litter),
         ## Make a rounded-up ticks per trap variable for Poisson and NegBin
         ## distribution compatibility, i.e. 0< x <1 == 1, and all other rounded to next
         ## larger integer
         tpt = ceiling(ticks_per_trap),
         vp_kPa_1yr = avg_1yr_vp..Pa./1000,
         total_clusters_log = log1p(total_clusters1m),
         canopy_ratio = avg_canopy_cover/100,
         litter_ratio = avg_pct_litter/100
  ) %>% 
  select(plot_id, inst_name, visit_year, d_since_fire, fri15yr, tick_abundance, 
         trap_effort, ticks_per_trap, tpt, ticks_per_trap01, 
         pct_litter_cover = avg_pct_litter, litter_depth_cm = avg_litter_depth_all, 
         pct_canopy_cover = avg_canopy_cover, standing_biomass_gm2 = avg_dry_standing_gm2, 
         total_clusters = total_clusters1m, deer:cottontail, vp_Pa_1yr = avg_1yr_vp..Pa., 
         cv_30yr_fire_days, d_since_fire_log, fri15yr_log, biomass_log, logit_litter,
         y_since_fire, y_since_fire_log, biomass_kgm2, vp_kPa_1yr, total_clusters_log,
         canopy_ratio, litter_ratio)


## Tick abundance models ####

### This "probability of ticks on traps" model doesn't make 
### much intuitive sense
tick_mod_ah <- glmer(cbind(tick_abundance, trap_effort) ~ d_since_fire_log + pct_litter_cover + litter_depth_cm + pct_canopy_cover + biomass_log + total_clusters + vp_Pa_1yr + (1|inst_name/plot_id),
  data = sem_mod_data, family = binomial(link = "logit"),
  na.action = "na.fail",
  control=glmerControl(
    optimizer = "Nelder_Mead",
    optCtrl=list(maxfun=1e6)))
summary(tick_mod_ah)
simulateResiduals(tick_mod_ah) %>% plot()
summary(effectsize::standardize(tick_mod_ah))
plot_model(tick_mod_ah, type = "pred", grid = T)

## This model would assume equal sampling effort - uses only 
## tick abundance, not controlling for trapping effort.
tick_abundance_nb <- glmer.nb(
  tick_abundance ~ d_since_fire_log + pct_litter_cover + litter_depth_cm + pct_canopy_cover + biomass_log + total_clusters + vp_Pa_1yr + (1|inst_name/plot_id),
  data = sem_mod_data, na.action = "na.fail",
  control=glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=1e6))
  )
plot(tick_abundance_nb)
simulateResiduals(tick_abundance_nb) %>% plot()
summary(tick_abundance_nb)
coefs(tick_abundance_nb)

## This model uses the rounded ticks per trap
## Outcome for relationship with time since fire is very different
tpt_nb <- glmer.nb(
  tpt ~ d_since_fire_log + pct_litter_cover + litter_depth_cm + pct_canopy_cover + biomass_log + total_clusters + vp_Pa_1yr + (1|inst_name/plot_id),
  data = sem_mod_data, na.action = "na.fail",
  control=glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=1e6))
  )
plot(tpt_nb)
simulateResiduals(tpt_nb) %>% plot()
effectsize::standardize(tpt_nb) %>% summary()

## Unrounded ticks per trap, violates integer value assumption
## compare to rounded integer model
ticks_per_trap_nb <- glmer.nb(
  ticks_per_trap ~ d_since_fire_log + pct_litter_cover + litter_depth_cm + pct_canopy_cover + biomass_log + total_clusters + vp_Pa_1yr + (1|inst_name/plot_id),
  data = sem_mod_data, na.action = "na.fail",
  control=glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=1e6))
  )
# plot(ticks_per_trap_nb)
simulateResiduals(ticks_per_trap_nb) %>% plot()
# summary(ticks_per_trap_nb)
# effectsize::standardize(ticks_per_trap_nb) %>% summary()

tab_model(tick_mod_ah, tick_abundance_nb, tpt_nb, ticks_per_trap_nb)
tab_model(
  effectsize::standardize(tick_mod_ah), 
  effectsize::standardize(tick_abundance_nb), 
  effectsize::standardize(tpt_nb)
  )


## Gamma ticks per trap model ####
sem_mod_data$ticks_per_trap001 <- sem_mod_data$ticks_per_trap + .001
ticks_per_trap_gamma <- glmer(
  ticks_per_trap001 ~ d_since_fire_log + pct_litter_cover + litter_depth_cm + pct_canopy_cover + biomass_log + total_clusters + vp_Pa_1yr + (1|inst_name/plot_id),
  data = sem_mod_data, na.action = "na.fail",
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
  effectsize::standardize(tick_mod_ah),
  effectsize::standardize(tick_abundance_nb),
  effectsize::standardize(tpt_nb), 
  effectsize::standardize(ticks_per_trap_nb), 
  effectsize::standardize(ticks_per_trap_gamma),
  show.values = T)

tab_model(tick_mod_ah, tick_abundance_nb, tpt_nb, ticks_per_trap_nb, ticks_per_trap_gamma)


## Poisson ticks per trap model ####
tpt_pois_sem <- glmer(
  tpt ~ d_since_fire_log + logit_litter + litter_depth_cm + pct_canopy_cover + biomass_log + total_clusters + vp_Pa_1yr + (1|inst_name/plot_id),
  data = sem_mod_data, na.action = "na.fail",
  family = poisson(link = "log"),
  control=glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=1e6))
  )
plot(tpt_pois_sem)
simulateResiduals(tpt_pois_sem) %>% plot()
summary(tpt_pois_sem)
effectsize::standardize(tpt_pois_sem) %>% summary()
coefs(tpt_pois_sem)
# plot_model(tpt_pois, show.values = T)
plot_model(tpt_pois_sem, type = "pred", show.data = T, grid = T)
d <- plot_model(tpt_pois_sem, type = "re")
d[[2]]$data


p <- plot_model(tpt_pois_sem, type = "pred", pred)
p[[1]]$data

tpt_pois2 <- glmer(
  tpt ~ d_since_fire + pct_litter_cover + litter_depth_cm + pct_canopy_cover + biomass_kgm2 + total_clusters + vp_kPa_1yr + (1|inst_name/plot_id), 
  data = sem_mod_data, na.action = "na.fail",
  family = poisson(link = "log"),
  control=glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=1e6))
)
simulateResiduals(tpt_pois2) %>% plot()
effectsize::standardize(tpt_pois2) %>% summary()
coefs(tpt_pois2)
visreg(tpt_pois2, "d_since_fire")
plot_model(tpt_pois2, type = "pred", grid = T, show.data = T)
plot_model(tpt_pois2, type = "est", show.values = T) +
  scale_y_continuous(limits = c(.9, 1.8))

tab_model(tpt_pois, tpt_pois2)

plot_models(tpt_pois, tpt_pois2, show.values = T)

p <- plot_models(
  effectsize::standardize(tick_mod_ah),
  effectsize::standardize(tick_abundance_nb),
  effectsize::standardize(tpt_nb),
  effectsize::standardize(tpt_pois),
  effectsize::standardize(ticks_per_trap_nb), 
  effectsize::standardize(ticks_per_trap_gamma),
  show.values = T)
p$data <- p$data %>% 
  mutate(group = case_when(
    group=="ticks_per_trap001" ~ "gamma",
    group=="ticks_per_trap" ~ "tpt_dbl_negbin",
    group=="tpt.4" ~ "tpt_int_pois",
    group=="tpt.3" ~ "tpt_int_negbin",
    group=="tick_abundance" ~ "tick_abundance_negbin",
    TRUE ~ "tick_abundance_binom"
  ))

p1 <- plot_model(tpt_pois, type = "pred", pred.type = "re", terms = c("d_since_fire_log", "inst_name"), show.data = T)
p1$data %>% as_tibble() %>% 
  mutate(group = str_sub(group, start = 1, end =-4)) %>% 
  ggplot(., aes(x, predicted, color = group)) +
  geom_path()
p <- plot_model(tpt_pois, type = "pred", pred.type = "re", terms = "d_since_fire_log", show.data = T)
p$data %>% as_tibble() %>% 
  ggplot(., aes(x, predicted)) +
  geom_point(data = tpt_pois@frame, aes(d_since_fire_log, tpt, color = inst_name)) +
  geom_path(size = 1, linetype = "dashed") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .2) +
  geom_path(data = p1$data, aes(x, predicted, color = group))


## Poisson tick abundance model ####
tick_abundance_pois <- glmer(
  tick_abundance ~ d_since_fire_log + logit_litter + litter_depth_cm + pct_canopy_cover + biomass_log + total_clusters + vp_Pa_1yr + (1|inst_name/plot_id),
  data = sem_mod_data, na.action = "na.fail",
  family = poisson(link = "log"),
  control=glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=1e6))
)
simulateResiduals(tick_abundance_pois) %>% plot()
summary(tick_abundance_pois)
effectsize::standardize(tick_abundance_pois) %>% summary()
effectsize::standardize_parameters(tick_abundance_pois)

plot_models(tick_abundance_nb, tick_abundance_pois, ticks_per_trap_gamma)


## No hosts ticks model ####

tpt_noHosts_pois <- glmer(
  tpt ~ d_since_fire_log + logit_litter + litter_depth_cm + pct_canopy_cover + biomass_log + vp_Pa_1yr + (1|inst_name/plot_id),
  data = sem_mod_data, 
  family = poisson(link = "log"),
  na.action = "na.fail",
  control=glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=1e6))
  )
# simulateResiduals(tpt_noHosts_pois) %>% plot()
# summary(tick_mod_noHosts_gamma)
# effectsize::standardize(tpt_noHosts_pois) %>% summary()


## Host models: all, deer only, no rabbit ####
host_mod_all <- glmer(
  total_clusters ~ d_since_fire_log + logit_litter + litter_depth_cm + pct_canopy_cover  + (1 | inst_name/plot_id), 
  data = sem_mod_data, family = poisson(link = "log"),
  control = glmerControl(optimizer="Nelder_Mead", optCtrl=list(maxfun=1e6))
  )
coefs(host_mod_all)
simulateResiduals(host_mod_all, n = 1000) %>% plot()
car::vif(host_mod_all)
summary(host_mod_all)
testZeroInflation(simulateResiduals(host_mod_all, n=1000))

host_mod_all_nb <- glmer.nb(
  total_clusters ~ d_since_fire_log + logit_litter + litter_depth_cm + pct_canopy_cover  + (1 | inst_name/plot_id), 
  data = sem_mod_data,
  control = glmerControl(optimizer="Nelder_Mead", optCtrl=list(maxfun=1e6))
)
simulateResiduals(host_mod_all_nb, n = 1000) %>% plot()
summary(host_mod_all_nb)
car::vif(host_mod_all_nb)

host_mod_zip <- glmmTMB::glmmTMB(total_clusters ~ d_since_fire_log + logit_litter + litter_depth_cm + pct_canopy_cover  + (1 | inst_name/plot_id), ziformula = ~1, data = sem_mod_data, family = poisson(link = "log"))
summary(host_mod_zip)
fixef(host_mod_zip); fixef(host_mod_all)

## this model includes the fire days (climate) "missing" path
total_clusters <- glmer(
  total_clusters ~ d_since_fire_log + logit_litter + litter_depth_cm + pct_canopy_cover + cv_30yr_fire_days + (1 | inst_name/plot_id),
  data = sem_mod_data, family = poisson(link = "log"),
  control = glmerControl(optimizer="Nelder_Mead", optCtrl=list(maxfun=1e6))
  )
car::vif(host_mod_all2)


## Litter/understory models ####
litter_cover_mod <- lmer(logit_litter ~ d_since_fire_log + pct_canopy_cover + (1|inst_name/plot_id),
                         data = sem_mod_data)
# coefs(litter_cover_mod)

litter_depth_mod <- lmer(litter_depth_cm ~ d_since_fire_log + pct_canopy_cover + (1|inst_name/plot_id), 
                         data = sem_mod_data)
# coefs(litter_depth_mod)

standing_biomass_mod <- lmer(biomass_log ~ d_since_fire_log + pct_canopy_cover + (1|inst_name/plot_id), data = sem_mod_data)
# plot(standing_biomass_mod)
# simulateResiduals(standing_biomass_mod) %>% plot()
# summary(standing_biomass_mod)


## Tree (overstory) canopy cover ####
# pine_mod <- lmer(logit_pct_pine ~ fri15yr_log + avg_canopy_cover +(1|inst_name/plot_id), 
#                  data = sem_data)
canopy_mod <- lmer(pct_canopy_cover ~ d_since_fire_log + fri15yr + (1|inst_name/plot_id), 
                   data = sem_mod_data)
# simulateResiduals(canopy_mod) %>% plot()
# summary(canopy_mod)

# Model with "missing" path
canopy_mod2 <- lmer(pct_canopy_cover ~ d_since_fire_log + fri15yr + vp_Pa_1yr + (1|inst_name/plot_id), 
                   data = sem_mod_data)
coefs(canopy_mod2)


## Fire regime models ####
time_since_fire_mod_pred <- glmer.nb(d_since_fire ~ fri15yr + (1|inst_name/plot_id),
         data = sem_mod_data)
simulateResiduals(time_since_fire_mod) %>% plot()
plot(time_since_fire_mod_pred)
plot_model(time_since_fire_mod_pred, show.values = T)
plot_model(time_since_fire_mod_pred, type = "pred", terms = "fri15yr [all]", show.data = T)

time_since_fire_mod <- lmer(d_since_fire_log ~ fri15yr + (1|inst_name/plot_id), data = sem_mod_data)
# simulateResiduals(time_since_fire_mod) %>% plot()
# plot(time_since_fire_mod)
# summary(time_since_fire_mod)
# plot_model(time_since_fire_mod, show.values = T)

fri_mod <- lmer(fri15yr ~ cv_30yr_fire_days + (1 | inst_name/plot_id), 
                data = sem_mod_data)
# simulateResiduals(fri_mod) %>% plot()
# plot(fri_mod)



## Weather ~ Climate model
avg_1yr_vp_mod <- lmer(vp_Pa_1yr ~ cv_30yr_fire_days + (1|inst_name), data = sem_mod_data)
# simulateResiduals(avg_1yr_vp_mod) %>% plot()
# summary(avg_1yr_vp_mod)
avg_1yr_vp_mod2 <- lmer(vp_Pa_1yr ~ pct_canopy_cover + cv_30yr_fire_days + (1|inst_name), data = sem_mod_data)
# summary(avg_1yr_vp_mod2)


## Path Models ####
# Convert to data frame class
sem_df <- as.data.frame(sem_mod_data)

## All hosts
all_hosts_psem <- psem(
  # m2,
  tpt_pois_sem,
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
fisherC(all_hosts_psem)

tpt_mod_pred <- plot_model(tpt_pois, type = "pred", show.data = F)
tpt_mod_pred_grid <- plot_model(
  tpt_pois, type = "pred", show.data = F, grid = T
  )
tpt_mod_pred_grid + xlab(" ") + ylab("Tick count")

tpt_mod_pred_grid$data <- tpt_mod_pred_grid$data %>% 
  filter(!group_col %in% c("total_clusters1m", "avg_canopy_cover","d_since_fire_log")) %>% 
  mutate(group = case_when(
    group == "logit_litter" ~ "Litter cover (logit)",
    group == "avg_litter_depth_all" ~ "Litter depth (cm)",
    group == "biomass_log" ~ "Standing biomass (log)",
    group == "avg_1yr_vp..Pa." ~ "Vapor pressure (Pa)"
    )
    )
tpt_mod_grid_fig <- tpt_mod_pred_grid + 
  xlab(" ") + ylab("Predicted tick count") + 
  theme_bw()
ggsave("figures/tpt_allhosts_marginal.pdf",tpt_mod_grid_fig, width = 5, height = 5, dpi = 300)

tpt_lcov <- tpt_mod_pred$logit_litter +
  xlab("Litter cover (logit)") + ylab("Predicted tick count") + ggtitle(" ") +
  theme_classic()
tpt_ldep <- tpt_mod_pred$avg_litter_depth_all +
  xlab("Litter depth (cm)") + ylab("Predicted tick count") + ggtitle(" ") +
  theme_classic()
tpt_ccov <- tpt_mod_pred$avg_canopy_cover +
  xlab("Canopy cover (%)") + ylab("Predicted tick count") + ggtitle(" ") +
  theme_classic()
tpt_vprs <- tpt_mod_pred$avg_1yr_vp..Pa. +
  geom_point(data = sem_df, aes(avg_1yr_vp..Pa., tpt)) +
  xlab("Vapor pressure (Pa)") + ylab("Predicted tick count") + ggtitle(" ") +
  theme_classic()
tpt_lcov + tpt_ldep + tpt_ccov + tpt_vprs


host_mod_pred_grid <- plot_model(
  host_mod_all, type = "pred", show.data = F, grid = T
)
host_mod_pred_grid + xlab(" ") + ylab("Tick count")

host_mod_pred_grid$data <- host_mod_pred_grid$data %>% 
  mutate(group = case_when(
    group == "logit_litter" ~ "Litter cover (logit)",
    group == "avg_litter_depth_all" ~ "Litter depth (cm)",
    group == "avg_canopy_cover" ~ "Canopy cover (%)",
    group == "d_since_fire_log" ~ "Days since fire (log)"
  )
  )
host_mod_pred_grid <- host_mod_pred_grid + 
  xlab(" ") + ylab("Predicted count of dung clusters") + 
  theme_bw()
ggsave("figures/hosts_allhosts_marginal.pdf", host_mod_pred_grid, width = 5, height = 5, dpi = 300)


lcov_mod_pred_grid <- plot_model(
  litter_cover_mod, type = "pred", show.data = F, grid = T
)
lcov_mod_pred_grid + xlab(" ") + ylab("Litter cover (logit)")

lcov_mod_pred_grid$data <- lcov_mod_pred_grid$data %>% 
  mutate(group = case_when(
    group == "avg_canopy_cover" ~ "Canopy cover (%)",
    group == "d_since_fire_log" ~ "Days since fire (log)"
  )
  )
lcov_mod_pred_grid <- lcov_mod_pred_grid + 
  xlab(" ") + ylab("Predicted litter cover (logit)") + 
  theme_bw()


ldep_mod_pred_grid <- plot_model(
  litter_depth_mod, type = "pred", show.data = F, grid = T
)
ldep_mod_pred_grid$data <- ldep_mod_pred_grid$data %>%
  filter(group == "d_since_fire_log") %>% 
  mutate(group = "Days since fire (log)")
ldep_mod_pred_grid <- ldep_mod_pred_grid + 
  xlab(" ") + ylab("Predicted litter depth (cm)") + 
  theme_bw()


sbio_mod_pred_grid <- plot_model(
  standing_biomass_mod, type = "pred", grid = T
  )
sbio_mod_pred_grid$data <- sbio_mod_pred_grid$data %>% 
  filter(group == "avg_canopy_cover") %>% 
  mutate(group = "Canopy cover (%)")
sbio_mod_pred_grid <- sbio_mod_pred_grid +
  xlab(" ") + ylab("Predicted understory biomass (log)") +
  theme_bw()

lcov_mod_pred_grid / (ldep_mod_pred_grid + sbio_mod_pred_grid)
ggsave("figures/veg_vars_allhosts_marginal.pdf", width = 5, height = 5, dpi = 300)


canp_mod_pred_grid <- plot_model(
  canopy_mod, type = "pred", grid = T
  )
canp_mod_pred_grid$data <- canp_mod_pred_grid$data %>% 
  filter(group == "d_since_fire_log") %>% 
  mutate(group = "Days since fire (log)")
canp_mod_pred_grid <- canp_mod_pred_grid +
  xlab(" ") + ylab ("Predicted canopy cover (%)") +
  theme_bw()

fire_mod_pred_grid <- plot_model(
  time_since_fire_mod, type = "pred", grid = T
)
fire_mod_pred_grid$data <- fire_mod_pred_grid$data %>% 
  mutate(group = "Fire return interval (years)")
fire_mod_pred_grid <- fire_mod_pred_grid +
  xlab(" ") + ylab("Predicted log-days since fire") +
  theme_bw()

fri_mod_pred_grid <- plot_model(
  fri_mod, type = "pred", grid = T
)
fri_mod_pred_grid$data <- fri_mod_pred_grid$data %>% 
  mutate(group = "30-y Average CV of Fire Days")
fri_mod_pred_grid <- fri_mod_pred_grid +
  xlab(" ") + ylab("Predicted fire return interval (years)") +
  theme_bw()

canp_mod_pred_grid + fire_mod_pred_grid + fri_mod_pred_grid
ggsave("figures/can_fire_vars_allhosts_marginal.pdf", width = 7, height = 3, dpi = 300)

## Treat all correlations as spurious/having a common unmeasured driver.
all_hosts_psem2 <- psem(
  tpt_pois_sem,
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
  ## should be the same either way. Check-in w/team on biological sense - BFA 
  ## had no biological explanation for these relationships
  
  data = sem_df
)
ah_psem2_summary <- summary(all_hosts_psem2, .progressBar = FALSE)
ah_psem2_summary

# Compare models without and wit correlations
ah_psem_summary$coefficients %>% 
  knitr::kable(digits = 3, caption = "Table xx. Coeffcient estimates all hosts path model")
ah_psem2_summary$coefficients %>% 
  knitr::kable(digits = 3, caption = "Table xx. Coeffcient estimates all hosts path model w/correlations for missing paths.")
ah_psem_summary$R2 %>% 
  knitr::kable(digits = 3, caption = "R2 all hosts path model")
ah_psem2_summary$R2 %>% 
  knitr::kable(digits = 3, caption = "R2 all hosts path model w/correlations for missing paths")
ah_psem_summary$AIC
ah_psem2_summary$AIC
## Same AIC
ah_psem_summary$Cstat
ah_psem2_summary$Cstat
## Change in C-statistic
ah_psem_summary$ChiSq
ah_psem2_summary$ChiSq
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


rbind(ah_psem_summary$AIC,
      # noRabbit_psem_summary$AIC,
      # deer_psem_summary$AIC,
      noHosts_psem_summary$AIC) %>% 
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

rbind(ah_psem_summary$ChiSq,
      # noRabbit_psem_summary$ChiSq,
      # deer_psem_summary$ChiSq,
      noHosts_psem_summary$ChiSq) %>% 
  mutate(model = c("all hosts", 
                   # "no rabbit", "deer only", 
                   "no hosts")) %>% 
  knitr::kable(digits = 3)


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
