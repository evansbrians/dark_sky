# Frequentist approach for Jared's arthropod study

# setup -------------------------------------------------------------------

library(lme4)
library(glmmTMB)
library(DHARMa)
library(AICcmodavg)
library(broom)
library(broom.mixed)
library(tidyverse)

# Read and pre-process data:

sticky <- 
  read_rds("data/trap_data.rds") %>% 
  filter(
    month(date) > 6,
    trap_type == "sticky",
  ) %>% 
  mutate(
    treatment = 
      treatment %>% 
      fct_collapse(
        dark = "control",
        other_level = "light"
      ),
    season = 
      case_when(
        month(date) == 7 ~ "early",
        .default = "late"
      ) %>% 
      factor(),
    .after = date
  ) %>% 
  filter(
    !(treatment == "light" & count == 0)
  ) %>% 
  select(
    transect_id:treatment,
    count
  )

# some summary information ------------------------------------------------

# General summary:

summary(sticky)

# How many transects sampled:

sticky %>% 
  distinct(transect_id) %>% 
  nrow()

# How many light vs. dark and early vs. late samples:

sticky %>% 
  count(treatment, season)

# Summary statistics:

sticky %>% 
  summarize(
    median = median(count),
    mean = mean(count),
    sample_size = n(),
    se = sd(count)/sqrt(sample_size),
    min = min(count),
    max = max(count),
    .by = c(treatment, season)
  ) %>% 
  arrange(treatment, season)

# How many traps are there with no arthropods?

sticky %>% 
  filter(count == 0) %>% 
  nrow()

sticky %>% 
  filter(treatment == "light") %>% 
  ggplot() +
  aes(x = count) +
  geom_histogram() +
  facet_wrap(~ season, nrow = 2)

# stats -------------------------------------------------------------------

# Make a formula list:

formulas <-
  list(
    
    # Since we are interested in a season effect, a reasonable null hypothesis
    # is that light impacts arthropod counts but season does not:
    
    light_only = count ~ treatment + (1 | transect_id),
    
    # The additive model suggests that season influence the number of arthropods
    # captured in traps, but not the influence of light vs. dark on the capture
    # rate:
    
    additive = count ~ treatment + season + (1 | transect_id),
    
    # The interaction model suggests that season influences arthropod response
    # to light.
    
    interaction = count ~ treatment * season + (1 | transect_id)
  )

## fit models: Poisson ----------------------------------------------------

# Try to fit a poisson glmer:

fits_poisson <-
  formulas %>% 
  map(
    ~ glmer(
      .x,
      data = sticky,
      family = poisson
    )
  )

# Model selection table:

aictab(fits_poisson)

# But ... check for overdispersion (c_hat > 4 represents a very poor fit):

c_hat(fits_poisson$interaction)

# Because extreme overdispersion is detected, a negative binomial is more
# appropriate!

## fit models negative binomial -------------------------------------------

fits_nb <-
  formulas %>%
  map(
    ~ glmmTMB(
      .x,
      data = sticky,
      family = nbinom2()
    )
  )

# Simulate residuals from the global model:

sim_resids <-
  fits_nb$interaction %>%
  simulateResiduals(n = 1000)

# Evaluate goodness-of-fit:

plot(sim_resids)

# If you like numbers more than visualization, this is how to interpret the Q-Q
# plot quantitatively (a p-value of > 0.05 is good -- it suggests uniformity):

testUniformity(sim_resids)

# Test residual dispersion specifically (we want a p-value of > 0.05):

testDispersion(sim_resids)

# Test for zero-inflation  (we want a p-value of > 0.05):

testZeroInflation(sim_resids)

# Plot residuals for treatment:

plotResiduals(sim_resids, form = sticky$treatment)

# Plot residuals for season:

plotResiduals(sim_resids, form = sticky$season)

# model selection and interpretation --------------------------------------

# Model selection table:

aictab(fits_nb)

# Interpretation: 
# - The data do not support our hypothesis that season influences arthropod
#   attraction to light.
# - Light treatment is significantly higher than dark (of course)
# - No evidence of a seasonal difference the number of arthropods captured
#   (additive or interaction models).
# - Although there is no formal statistical support, the data *may* indicate 
#   that with more samples the light treatment effect is greater late than 
#   early season

# Model summary

tidy(fits_nb$light_only, conf.int = TRUE) %>% 
  select(
    term,
    estimate,
    conf.low,
    conf.high
  )

# Model estimates:

new_data <-
  tibble(
    treatment = levels(sticky$treatment)
  )

# Predict on the link scale first, then back-transform:

preds <-
  fits_nb$light_only %>%
  augment(
    newdata      = new_data,
    se_fit       = TRUE,
    type.predict = "link",
    re.form      = NA
  ) %>%
  mutate(
    estimate = exp(.fitted),
    lower = 
      exp(
        .fitted - qnorm(0.975) * .se.fit
      ),
    upper = 
      exp(
        .fitted + qnorm(0.975) * .se.fit
      )
  )

# plots -------------------------------------------------------------------

# Simple boxplot:

sticky %>% 
  mutate(
    across(
      season:treatment,
      ~ str_to_title(.x)
    )
  ) %>% 
  ggplot() +
  aes(
    x = count,
    y = season,
    fill = treatment
  ) +
  geom_boxplot() +
  scale_x_continuous(
    trans  = scales::pseudo_log_trans(base = 10),
    breaks = c(0, 10, 100, 1000)
  ) +
  scale_fill_manual(
    values = 
      c(
        Dark = "#999", 
        Light = "#E8A838"
      )
  ) +
  labs(
    x = "Count",
    y = "Season",
    fill = "Treatment"
  ) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank()
  )

# Mimicking Bayesian plotting version:

sticky %>% 
  mutate(
    across(
      season:treatment,
      ~ str_to_title(.x)
    )
  ) %>% 
  ggplot() +
  aes(
    x = count,
    y = season,
    fill = treatment
  ) +
  stat_slabinterval() +
  facet_wrap(
    ~ treatment,
    ncol = 1
  ) +
  scale_x_continuous(
    trans  = scales::pseudo_log_trans(base = 10),
    breaks = c(0, 10, 100, 1000),
    labels = scales::label_comma()
  ) +
  scale_fill_manual(
    values = 
      c(
        Dark = "#999", 
        Light = "#E8A838"
      )
  ) +
  labs(
    x = "Count",
    y = "Season",
    fill = "Treatment"
  ) +
  theme_bw()

# Pattern:
# - Confidence intervals (50% and 95%) of early vs. late overlap for both dark
#   and light treatments.
# - Light treatment counts is considerable higher than dark.
# - Dark treatments seem to be fully equivalent in early and late treatments
# - Late season has a much higher (but not significant) arthropod count values
#   for the light treatment.
# - A trap in the late light treatment (but none in the early light treatment)
#   had a count of 0, which pulled the distribution down.
#
# Interpretation: Although the data do not support our hypothesis that season
# influences arthropod attraction to light, the distribution of the data *may*
# indicate that, with additional sampling, the light treatment effect is greater
# late than early season. It's too early in the sampling process to tell.

