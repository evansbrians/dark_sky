
# Bayesian workflow

# setup -------------------------------------------------------------------

library(brms)
library(ggdist)
library(tidybayes)
library(tidyverse)

# read and pre-process data -----------------------------------------------

sticky <- 
  read_rds("data/trap_data.rds") %>% 
  filter(
    month(date) > 6,
    trap_type == "sticky"
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

# stats -------------------------------------------------------------------

## make a formula list ----------------------------------------------------

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

## define priors ----------------------------------------------------------

# Priors encode our beliefs about plausible parameter values before seeing
# the data. "Weakly informative" priors constrain the parameter space enough
# to regularize estimation without strongly favoring any particular value.
#
# All fixed effects are on the log scale (the NB link function), so:
#   exp(0)  = 1.0x  (no effect)
#   exp(2)  = 7.4x  (large effect)
#   exp(-2) = 0.14x (large negative effect)
# Normal(0, 2) is wide enough to accommodate large effects while down-weighting
# implausibly extreme values.
#
# The random effect SD prior Normal(0, 1) is a half-normal in practice because
# brms constrains SD > 0. This weakly regularizes transect-level variance.
#
# The shape parameter of the negative binomial controls dispersion:
#   * Small shape = high overdispersion (variance >> mean)
#   * Large shape = low overdispersion (approaches Poisson)
# Gamma(2, 0.1) centers the prior on moderate-to-high dispersion, which is
# typical for arthropod count data.

priors <-
  c(
    prior(
      normal(0, 2),
      class = b
    ),
    prior(
      normal(0, 1), 
      class = sd
    ),
    prior(
      gamma(2, 0.1),
      class = shape
    )
  )

## fit models -------------------------------------------------------------

# Key arguments:
#   chains  = 4: run 4 independent MCMC chains (allows convergence checking)
#   iter    = 2000: total iterations per chain
#   warmup  = 1000: first 1000 iterations discarded (sampler "tuning" phase)
#   cores   = 4: run chains in parallel to reduce wall time
#   seed    = 123: ensures reproducibility of random number generation
#   save_pars = save_pars(all = TRUE): required for moment matching in LOO
#   adapt_delta = 0.99: cautious step sizes to reduce divergent transitions

fits_bayes <-
  formulas %>% 
  map(
    ~ brm(
      formula = .x,
      data = sticky,
      family = negbinomial(),
      prior = priors,
      chains = 4,
      iter = 2000,
      warmup = 1000,
      cores = 4,
      seed = 123,
      save_pars = save_pars(all = TRUE),
      control = list(adapt_delta = 0.99)
    )
  )

# mcmc diagnostics --------------------------------------------------------

# Note: mcmc stands for Markov Chain Monte Carlo

# Trace plots: chains should mix well with no trends or divergences:

plot(fits_bayes$interaction)

# Numerical summary: Rhat should be ~1.00, Bulk/Tail ESS should be > 400:

summary(fits_bayes$interaction)

## posterior predictive check (Bayesian analog of DHARMa) -----------------

# Overlays the distribution of replicated data (yrep) on observed data (y).
# A good fit shows yrep closely tracking y:

pp_check(
  fits_bayes$interaction,
  ndraws = 100
)

# Check specifically for zero-inflation:

pp_check(
  fits_bayes$interaction,
  type = "stat",
  stat = \(.x) mean(.x == 0)
)

# model selection ---------------------------------------------------------

# Model comparison via Leave-one-out (LOO) cross-validation (Bayesian analog of
# AICc).

# Add LOO criterion to each model:

fits_bayes_loo <-
  fits_bayes %>%
  map(
    ~ add_criterion(
      .x,
      criterion    = "loo",
      moment_match = TRUE
    )
  )

# Generate comparison table:

fits_bayes_loo %>%
  map(~ .x$criteria$loo) %>%
  loo_compare(x = .)

# Interpreting the output:

# * elpd stands for expected log predictive density (out-of-sample predictive
#   accuracy)
# * The top model is assigned elpd_diff = 0 with others relative to that elpd
# * se_diff is the standard error of the elpd_diff
# * If elpd_diff / se_diff is < 2 the difference is not distinguishable from 
#   noise (equivalent predictive accuracy)

# posterior predictions with credible intervals ---------------------------

# Important! I am plotting the interaction here. We wouldn't typically do so
# because this model was not supported by the data.

# Make a newdata object for predicted values:

predictions <-
  fits_bayes$interaction %>%
  epred_draws(
    newdata = 
      sticky %>% 
      distinct(treatment, season),
    re_formula = NA,
    ndraws = 4000
  ) %>%
  group_by(treatment, season) %>%
  summarize(
    estimate = median(.epred),
    lower = quantile(.epred, 0.025),
    upper = quantile(.epred, 0.975),
    .groups = "drop"
  )

## plot of estimates and credible intervals -------------------------------

predictions %>%
  ggplot() +
  aes(
    x = season,
    y = estimate,
    ymin = lower,
    ymax = upper,
    color = treatment,
    group = treatment
  ) +
  geom_pointrange(
    position = position_dodge(width = 0.3),
    size = 0.8
  ) +
  scale_y_log10(
    labels = scales::label_comma()
  ) +
  scale_color_manual(
    values = 
      c(
        dark = "#4E4E4E",
        light = "#E8A838"
      )
  ) +
  labs(
    x = "Season",
    y = "Estimated arthropod count",
    color = "Treatment"
  ) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank()
  )

## plot of estimates, credible intervals, and distributions ---------------

fits_bayes$interaction %>%
  epred_draws(
    newdata = 
      sticky %>% 
      distinct(treatment, season),
    re_formula = NA,
    ndraws = 4000
  ) %>% 
  ggplot() +
  aes(
    x = .epred,
    y = season,
    fill = treatment
  ) +
  stat_halfeye(
    alpha = 0.6,
    .width = c(0.50, 0.95),
    point_interval = median_qi,
  ) +
  scale_x_log10(
    labels = scales::label_comma()
  ) +
  scale_fill_manual(
    values =
      c(
        dark = "#4E4E4E",
        light = "#E8A838"
      )
  ) +
  facet_wrap(~ treatment, ncol = 1) +
  labs(
    x = "Estimated arthropod count (log scale)",
    y = "Season",
    fill = "Treatment"
  ) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  )

