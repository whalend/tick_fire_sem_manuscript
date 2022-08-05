## Models to generate predicted relationships figures
## Variables are less/not transformed compared to those used in the
## path analysis where the goal was to make things play nice together

## Setup ####
library(tidyverse)
library(sjPlot)
library(patchwork)
library(MuMIn)
library(DHARMa)
library(lme4)
library(piecewiseSEM)

sem_data <- read_csv("data/path_analysis_data.csv")

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


# Ticks submodels ####
## Model used in SEM
tpt_pois_sem <- glmer(
  tpt ~ d_since_fire_log + logit_litter + litter_depth_cm + pct_canopy_cover + biomass_log + total_clusters + vp_Pa_1yr + (1|inst_name/plot_id),
  data = sem_mod_data, na.action = "na.fail",
  family = poisson(link = "log"),
  control=glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=1e6))
)

## Model without transformations
tpt_pois2 <- glmer(
  tpt ~ d_since_fire + pct_litter_cover + litter_depth_cm + pct_canopy_cover + standing_biomass_gm2 + total_clusters + vp_Pa_1yr + (1|inst_name/plot_id), 
  data = sem_mod_data, na.action = "na.fail",
  family = poisson(link = "log"),
  control=glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=1e6))
)

plot_models(tpt_pois_sem, tpt_pois2, show.values = T)
plot_model(tpt_pois_sem, type = "pred", terms = "logit_litter") |
  plot_model(tpt_pois2, type = "pred", terms = "pct_litter_cover")
plot_model(tpt_pois_sem, type = "pred", terms = "biomass_log") |
  plot_model(tpt_pois2, type = "pred", terms = "standing_biomass_gm2")
plot_model(tpt_pois_sem, type = "pred", terms = "d_since_fire_log") |
  plot_model(tpt_pois2, type = "pred", terms = "d_since_fire")

## Deer-only model #### 
tpt_deer_nb <- glmer.nb(
  tpt ~ d_since_fire_log + logit_litter + litter_depth_cm + pct_canopy_cover + biomass_log + deer + vp_Pa_1yr + (1|inst_name/plot_id),
  data = sem_mod_data, na.action = "na.fail",
  control=glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=1e6))
)
tpt_deer_pois <- glmer(
  tpt ~ d_since_fire_log + logit_litter + litter_depth_cm + pct_canopy_cover + biomass_log + deer + vp_Pa_1yr + (1|inst_name/plot_id),
  data = sem_mod_data, na.action = "na.fail",
  family = poisson(link = "log"),
  control=glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=1e6))
)
# simulateResiduals(tpt_deer_pois) %>% plot()
# summary(tpt_deer_pois)
# summary(tpt_deer_nb)$coefficients
# summary(tpt_deer_pois)$coefficients
## Similar enough, though in the model summary the p-values are a bit smaller
## in the poisson model

## Exclude rabbit data ####
tpt_noRabb_nb <- glmer.nb(
  tpt ~ d_since_fire_log + logit_litter + litter_depth_cm + pct_canopy_cover + biomass_log + deer_other + vp_Pa_1yr + (1|inst_name/plot_id),
  data = sem_mod_data, na.action = "na.fail",
  control=glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=1e6))
)

tpt_noRabb_pois <- glmer(
  tpt ~ d_since_fire_log + logit_litter + litter_depth_cm + pct_canopy_cover + biomass_log + deer_other + vp_Pa_1yr + (1|inst_name/plot_id),
  data = sem_mod_data, na.action = "na.fail",
  family = poisson(link = "log"),
  control=glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=1e6))
)
# summary(tpt_noRabb_pois)
# effectsize::standardize(tpt_noRabb_pois) %>% summary()
# effectsize::standardize(tpt_noRabb_pois) %>% 
#   simulateResiduals() %>% 
#   plot()
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

## Rabbit only ticks model ####
tpt_Rabb_pois <- glmer(
  tpt ~ d_since_fire_log + logit_litter + litter_depth_cm + pct_canopy_cover + biomass_log + cottontail + vp_Pa_1yr + (1|inst_name/plot_id),
  data = sem_mod_data, 
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
  tpt ~ d_since_fire_log + logit_litter + litter_depth_cm + pct_canopy_cover + biomass_log + vp_Pa_1yr + (1|inst_name/plot_id),
  data = sem_mod_data, 
  na.action = "na.fail",
  control=glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=1e6))
)

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

# AICc(tpt_nb, tpt_noRabb_nb, tpt_deer_nb, tpt_noHosts_nb) %>%
# arrange(AICc)
# AICc(tpt_pois, tpt_noRabb_pois, tpt_deer_pois, tpt_noHosts_pois) %>%
# arrange(AICc)
tpt_noHosts_pois2 <- glmer(
  tpt ~ d_since_fire + pct_litter_cover + litter_depth_cm + pct_canopy_cover + standing_biomass_gm2 + vp_Pa_1yr + (1|inst_name/plot_id), 
  data = sem_mod_data, na.action = "na.fail",
  family = poisson(link = "log"),
  control=glmerControl(optimizer = "Nelder_Mead", optCtrl=list(maxfun=1e6))
)

# Hosts submodels ####
## Model used in SEM
host_mod_all <- glmer(
  total_clusters ~ d_since_fire_log + logit_litter + litter_depth_cm + pct_canopy_cover + (1 | inst_name/plot_id), 
  data = sem_mod_data, family = poisson(link = "log"),
  control = glmerControl(optimizer="Nelder_Mead", optCtrl=list(maxfun=1e6))
)

host_mod_deer <- glmer(
  deer ~ d_since_fire_log + logit_litter + litter_depth_cm + pct_canopy_cover + (1|inst_name/plot_id),
  data = sem_mod_data, family = poisson(link = "log"),
  control = glmerControl(optimizer="Nelder_Mead", optCtrl=list(maxfun=1e6))
)

host_mod_rabbit <- glmer(
  cottontail ~ d_since_fire_log + logit_litter + litter_depth_cm + pct_canopy_cover + (1|inst_name/plot_id),
  data = sem_mod_data, family = poisson(link = "log"),
  control = glmerControl(optimizer="Nelder_Mead", optCtrl=list(maxfun=1e6))
)

host_mod_other <- glmer(
  deer_other ~ d_since_fire_log + logit_litter + litter_depth_cm + pct_canopy_cover + (1|inst_name/plot_id),
  data = sem_mod_data, family = poisson(link = "log"),
  control = glmerControl(optimizer="Nelder_Mead", optCtrl=list(maxfun=1e6))
)

## Model without transformation
host_mod_all2 <- glmer(
  total_clusters ~ d_since_fire + pct_litter_cover + pct_canopy_cover + litter_depth_cm + (1 | inst_name/plot_id), 
  data = sem_mod_data, family = poisson(link = "log"),
  control = glmerControl(optimizer="Nelder_Mead", optCtrl=list(maxfun=1e6))
)


## Litter/understory models ####
## Models used in SEM
litter_cover_mod <- lmer(logit_litter ~ d_since_fire_log + pct_canopy_cover + (1|inst_name/plot_id),
                         data = sem_mod_data)

litter_depth_mod <- lmer(litter_depth_cm ~ d_since_fire_log + pct_canopy_cover + (1|inst_name/plot_id), 
                         data = sem_mod_data)

standing_biomass_mod <- lmer(biomass_log ~ d_since_fire_log + pct_canopy_cover + (1|inst_name/plot_id), 
                             data = sem_mod_data)

## Models without transformation
litter_cover_mod2 <- lmer(pct_litter_cover ~ d_since_fire + pct_canopy_cover + (1|inst_name/plot_id),
                         data = sem_mod_data)

litter_depth_mod2 <- lmer(litter_depth_cm ~ d_since_fire + pct_canopy_cover + (1|inst_name/plot_id), 
                         data = sem_mod_data)

standing_biomass_mod2 <- lmer(standing_biomass_gm2 ~ d_since_fire + pct_canopy_cover + (1|inst_name/plot_id), 
                              data = sem_mod_data)


## Tree (overstory) canopy cover ####
# pine_mod <- lmer(logit_pct_pine ~ fri15yr_log + avg_canopy_cover +(1|inst_name/plot_id), 
#                  data = sem_data)
## Models used in SEM
canopy_mod <- lmer(pct_canopy_cover ~ d_since_fire_log + fri15yr + (1|inst_name/plot_id), 
                   data = sem_mod_data)

## Model without transformation
canopy_mod2 <- lmer(pct_canopy_cover ~ d_since_fire + fri15yr + (1|inst_name/plot_id), 
                   data = sem_mod_data)

# coefs(list(canopy_mod, canopy_mod2))

## Fire regime models ####
# Model used in SEM
time_since_fire_mod <- lmer(d_since_fire_log ~ fri15yr + (1|inst_name/plot_id), 
                            data = sem_mod_data)
# simulateResiduals(time_since_fire_mod) %>% plot()
# plot(time_since_fire_mod)
# summary(time_since_fire_mod)
# plot_model(time_since_fire_mod, show.values = T)
plot_model(time_since_fire_mod, type = "pred", terms = "fri15yr [all]", show.data = T)

# Model without transformation
time_since_fire_mod2 <- glmer.nb(d_since_fire ~ fri15yr + (1|inst_name/plot_id),
                                     data = sem_mod_data)
# simulateResiduals(time_since_fire_mod2) %>% plot()
# plot(time_since_fire_mod2)
# plot_model(time_since_fire_mod2, show.values = T)
plot_model(time_since_fire_mod2, type = "pred", terms = "fri15yr [n=50]", show.data = T)

# coefs(list(time_since_fire_mod, time_since_fire_mod2))
# effectsize::standardize(time_since_fire_mod2) %>% fixef()


fri_mod <- lmer(fri15yr ~ cv_30yr_fire_days + (1 | inst_name/plot_id), 
                data = sem_mod_data)
# simulateResiduals(fri_mod) %>% plot()
# plot(fri_mod)

fri_mod2 <- lmer(log(fri15yr) ~ cv_30yr_fire_days + (1 | inst_name/plot_id), 
                 data = sem_mod_data)
simulateResiduals(fri_mod2) %>% plot()
plot(fri_mod2)

sem_mod_data$fri15yr_int <- ceiling(sem_mod_data$fri15yr)
fri_mod_pred <- glmer.nb(fri15yr_int ~ cv_30yr_fire_days + (1 | inst_name/plot_id), 
                      data = sem_mod_data,
                      # family = poisson(link = "log"),
                      control = glmerControl(
                        optimizer="Nelder_Mead", optCtrl=list(maxfun=1e6))
)
# simulateResiduals(fri_mod_pred) %>% plot()
# plot(fri_mod_pred)
# plot_model(fri_mod, type = "pred", show.data = T)
# plot_model(fri_mod_pred, type = "pred", show.data = T)


## Weather ~ Climate model ####
avg_1yr_vp_mod <- lmer(vp_Pa_1yr ~ cv_30yr_fire_days + (1|inst_name), data = sem_mod_data)
# simulateResiduals(avg_1yr_vp_mod) %>% plot()
# summary(avg_1yr_vp_mod)

# Evaluate "missing" path: canopy cover
avg_1yr_vp_mod2 <- lmer(vp_Pa_1yr ~ pct_canopy_cover + cv_30yr_fire_days + (1|inst_name), data = sem_mod_data)
# summary(avg_1yr_vp_mod2)


tpt_pred <- plot_model(tpt_pois2, type = "pred", show.data = T)
# unlist(tpt_pred)

tpt_pred_data <- tibble()
n <- length(fixef(tpt_pois2))-1
sjData <- function(sjPlot_obj){
  n = length(fixef(tpt_pois2))-1
  d = sjPlot_obj
  d = tpt_pred
  for (i in 1:n) {
    j = d[[i]]$data
    tpt_pred_data = rbind(tpt_pred_data, j)
  }
  return(as_tibble(tpt_pred_data))
}
k <- sjData(tpt_pred)
unique(k$group)

tpt_est <- plot_model(tpt_pois2, show.values = T)
terms <- tpt_est$data %>% 
  filter(p.value<.05)

k <- k %>% 
  filter(group%in%terms$term)
ggplot(k, aes(x, predicted, ymin = conf.low, ymax = conf.high)) +
  geom_ribbon(alpha = .2) +
  geom_path() +
  facet_wrap(~group, scales = "free_x")

tpt_pred$pct_litter_cover
tpt_pred$litter_depth_cm
tpt_pred$standing_biomass_gm2
tpt_pred$vp_Pa_1yr

tpt_noHosts_pois2
tpt_noHosts_est <- plot_model(tpt_noHosts_pois2, show.values = T)
tpt_noHosts_pred <- plot_model(tpt_noHosts_pois2, type = "pred", show.data = T)
tpt_noHosts_pred$pct_litter_cover
tpt_noHosts_pred$litter_depth_cm
tpt_noHosts_pred$pct_canopy_cover
tpt_noHosts_pred$standing_biomass_gm2
tpt_noHosts_pred$vp_Pa_1yr

host_mod_all2
host_est <- plot_model(host_mod_all2, show.values = T)
host_pred <- plot_model(host_mod_all2, type = "pred", show.data = T)
host_pred$d_since_fire
host_pred$pct_litter_cover
host_pred$pct_canopy_cover
host_pred$litter_depth_cm


litter_cover_mod2 
litter_depth_mod2
canopy_mod2
standing_biomass_mod2
time_since_fire_mod2
fri_mod2
avg_1yr_vp_mod2
