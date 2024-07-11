# Script to simulate synthetic data for the public repository
library(tidyverse)
library(mvtnorm)
library(simsurv)
library(synthpop)
library(VAJointSurv)

# Read MLO data
# Note: file paths are omitted
data_full_mlo <- readRDS(file = "...")
data_surv_mlo <- readRDS(file = "...")

# Random seed for reproducibility
set.seed(10147298)

# Subset smaller datasets
ids <- sample(x = unique(data_surv_mlo$studieid), size = ceiling(nrow(data_surv_mlo) * 0.2))
data_full_mlo <- filter(data_full_mlo, studieid %in% ids)
data_surv_mlo <- filter(data_surv_mlo, studieid %in% ids)

# First, simulating synthetic baseline data
data_surv_mlo_sub <- data_surv_mlo |>
  select(N, t, d, fattyarea, densearea, menopause_status, age_inclusion_in_karma, bmi, hrt_status, family_bc, birth_times) |>
  mutate(d = factor(d))

# Synthesise baseline data
synth_baseline <- syn(data = data_surv_mlo_sub)
synth_surv_mlo <- synth_baseline$syn |>
  mutate(id = row_number()) |>
  select(id, age_inclusion_in_karma, bmi, hrt_status, family_bc, birth_times, N)
rm(synth_baseline)
gc()

# Now, we fit a simplified JM to derive model parameters
fa_marker <- marker_term(
  formula = sqrt(fattyarea) ~ bmi + hrt_status + family_bc,
  id = studieid,
  data = data_full_mlo,
  time_fixef = poly_term(x = age_mammogram, degree = 1, raw = TRUE), # Linear
  time_rng = poly_term(x = age_mammogram, degree = 0, intercept = TRUE, raw = TRUE) # Linear + Intercept
)
da_marker <- marker_term(
  formula = sqrt(densearea) ~ bmi + hrt_status + family_bc,
  id = studieid,
  data = data_full_mlo,
  time_fixef = poly_term(x = age_mammogram, degree = 1, raw = TRUE), # Linear
  time_rng = poly_term(x = age_mammogram, degree = 0, intercept = TRUE, raw = TRUE) # Linear + Intercept
)
ssub <- surv_term(
  Surv(single_t0, single_t, single_d) ~ bmi + hrt_status + family_bc + birth_times,
  id = studieid,
  data = data_surv_mlo,
  time_fixef = ns_term(single_t, df = 3),
  delayed = single_t0 > 0
)
comp_obj <- joint_ms_ptr(
  markers = list(fa_marker, da_marker),
  survival_terms = list(ssub),
  max_threads = parallel::detectCores()
)
start_val <- joint_ms_start_val(object = comp_obj)
opt_out <- joint_ms_opt(object = comp_obj, par = start_val)
fmt_est <- joint_ms_format(object = comp_obj, par = opt_out$par)

# Simulate random effects
# In the joint model above we assume a random intercept only:
# we now add a random slope of time, where we fix the variances to abs of the estimated
# fixed slope parameters and draw the correlations between random effects at random.
Sigma <- matrix(data = 0, ncol = 4, nrow = 4)
Sigma[1, 1] <- fmt_est$vcov$vcov_vary[1, 1]
Sigma[3, 3] <- fmt_est$vcov$vcov_vary[2, 2]
Sigma[1, 3] <- fmt_est$vcov$vcov_vary[1, 2]
Sigma[3, 1] <- fmt_est$vcov$vcov_vary[2, 1]
Sigma[2, 2] <- abs(fmt_est$markers[[1]]$fixef_vary)
Sigma[4, 4] <- abs(fmt_est$markers[[2]]$fixef_vary)
SDs <- sqrt(diag(Sigma))
Sigma <- cov2cor(Sigma)
Sigma[upper.tri(Sigma) & Sigma == 0] <- runif(n = sum(upper.tri(Sigma) & Sigma == 0), min = -1, max = 1)
Sigma[lower.tri(Sigma)] <- t(Sigma)[lower.tri(Sigma)]
Sigma <- outer(SDs, SDs) * Sigma
# Now simulate
reffs <- rmvnorm(n = length(unique(synth_surv_mlo$id)), sigma = Sigma)
reffs <- as.data.frame(reffs)
reffs$id <- synth_surv_mlo$id
reffs <- reffs |>
  rename(b0_fa = V1, b1_fa = V2, b0_da = V3, b1_da = V4)

# Simulate error terms
errs <- rmvnorm(n = length(unique(synth_surv_mlo$id)), sigma = fmt_est$vcov$vcov_marker)
errs <- as.data.frame(errs)
errs$id <- synth_surv_mlo$id
errs <- errs |>
  rename(e_fa = V1, e_da = V2)

# Create dataset with number of observations per study subject
synth_full_mlo <- map_dfr(.x = seq(nrow(synth_surv_mlo)), .f = function(i) {
  data.frame(id = synth_surv_mlo$id[i], n = seq(synth_surv_mlo$N[i]))
})
# Sample time intervals
synth_full_mlo$tdiff <- sample(x = data_full_mlo$t - data_full_mlo$t0, size = nrow(synth_full_mlo), replace = TRUE)

# Merge baseline data, random effects, error terms
synth_full_mlo <- synth_full_mlo |>
  left_join(synth_surv_mlo, by = "id") |>
  left_join(reffs, by = "id") |>
  left_join(errs, by = "id")

# Fix start-stop notation
synth_full_mlo <- synth_full_mlo |>
  mutate(tdiff = ifelse(n == 1, 0, tdiff)) |>
  mutate(tdiffcum = cumsum(tdiff), .by = "id") |>
  mutate(t0 = age_inclusion_in_karma - 40 + tdiffcum) |>
  mutate(t = lead(t0), .by = "id")

# Create fattyarea, densearea
synth_full_mlo$sqrt_fattyarea <-
  (fmt_est$markers[[1]]$fixef[1] + synth_full_mlo$b0_fa) +
  (fmt_est$markers[[1]]$fixef_vary[1] + synth_full_mlo$b1_fa) * synth_full_mlo$t0 +
  (fmt_est$markers[[1]]$fixef[2] * synth_full_mlo$bmi) +
  (fmt_est$markers[[1]]$fixef[3] * (synth_full_mlo$hrt_status == "Previous user")) +
  (fmt_est$markers[[1]]$fixef[4] * (synth_full_mlo$hrt_status == "Current user")) +
  (fmt_est$markers[[1]]$fixef[5] * (synth_full_mlo$hrt_status == "Missing")) +
  (fmt_est$markers[[1]]$fixef[6] * (synth_full_mlo$family_bc == "Yes")) +
  (fmt_est$markers[[1]]$fixef[7] * (synth_full_mlo$family_bc == "Missing")) +
  synth_full_mlo$e_fa
# Winsorise these values at the range of the observed data
synth_full_mlo$sqrt_fattyarea <- pmax(min(sqrt(data_full_mlo$fattyarea)), synth_full_mlo$sqrt_fattyarea, na.rm = TRUE)
synth_full_mlo$sqrt_fattyarea <- pmin(max(sqrt(data_full_mlo$fattyarea)), synth_full_mlo$sqrt_fattyarea, na.rm = TRUE)

#
synth_full_mlo$sqrt_densearea <-
  (fmt_est$markers[[2]]$fixef[1] + synth_full_mlo$b0_da) +
  (fmt_est$markers[[2]]$fixef_vary[1] + synth_full_mlo$b1_da) * synth_full_mlo$t0 +
  (fmt_est$markers[[2]]$fixef[2] * synth_full_mlo$bmi) +
  (fmt_est$markers[[2]]$fixef[3] * (synth_full_mlo$hrt_status == "Previous user")) +
  (fmt_est$markers[[2]]$fixef[4] * (synth_full_mlo$hrt_status == "Current user")) +
  (fmt_est$markers[[2]]$fixef[5] * (synth_full_mlo$hrt_status == "Missing")) +
  (fmt_est$markers[[2]]$fixef[6] * (synth_full_mlo$family_bc == "Yes")) +
  (fmt_est$markers[[2]]$fixef[7] * (synth_full_mlo$family_bc == "Missing")) +
  synth_full_mlo$e_da
# Winsorise these values at the range of the observed data
synth_full_mlo$sqrt_densearea <- pmax(min(sqrt(data_full_mlo$densearea)), synth_full_mlo$sqrt_densearea, na.rm = TRUE)
synth_full_mlo$sqrt_densearea <- pmin(max(sqrt(data_full_mlo$densearea)), synth_full_mlo$sqrt_densearea, na.rm = TRUE)

# Transform these values back to original scale
synth_full_mlo <- synth_full_mlo |>
  mutate(fattyarea = sqrt_fattyarea^2, densearea = sqrt_densearea^2) |>
  select(-starts_with("sqrt_"))

# Simulate time to breast cancer
# Hazard function
haz <- function(t, x, betas, ...) {
  betas[["shape"]] * (t^(betas[["shape"]] - 1)) * exp(
    betas[["betaEvent_intercept"]] +
      betas[["betaEvent_bmi"]] * x[["bmi"]] +
      betas[["betaEvent_hrt_statusPrevious"]] * (x[["hrt_status"]] == "Previous user") +
      betas[["betaEvent_hrt_statusCurrent"]] * (x[["hrt_status"]] == "Current user") +
      betas[["betaEvent_family_bcYes"]] * (x[["family_bc"]] == "Yes") +
      betas[["betaEvent_family_bcMissing"]] * (x[["family_bc"]] == "Missing") +
      betas[["betaEvent_birth_times1"]] * (x[["birth_times"]] == "1") +
      betas[["betaEvent_birth_times2"]] * (x[["birth_times"]] == "2") +
      betas[["betaEvent_birth_times3p"]] * (x[["birth_times"]] == "3+") +
      betas[["betaEvent_birth_timesMissing"]] * (x[["birth_times"]] == "Missing") +
      # Fatty
      betas[["betaEvent_assoc_fattyarea"]] * (
        betas[["betaLong_fattyarea_intercept"]] +
          betas[["betaLong_fattyarea_slope"]] * t +
          betas[["betaLong_fattyarea_bmi"]] * x[["bmi"]] +
          betas[["betaLong_fattyarea_hrt_statusPrevious"]] * (x[["hrt_status"]] == "Previous user") +
          betas[["betaLong_fattyarea_hrt_statusCurrent"]] * (x[["hrt_status"]] == "Current user") +
          betas[["betaLong_fattyarea_hrt_statusMissing"]] * (x[["hrt_status"]] == "Missing") +
          betas[["betaLong_fattyarea_family_bcYes"]] * (x[["family_bc"]] == "Yes") +
          betas[["betaLong_fattyarea_family_bcMissing"]] * (x[["family_bc"]] == "Missing")
      ) +
      # Dense
      betas[["betaEvent_assoc_densearea"]] * (
        betas[["betaLong_densearea_intercept"]] +
          betas[["betaLong_densearea_slope"]] * t +
          betas[["betaLong_densearea_bmi"]] * x[["bmi"]] +
          betas[["betaLong_densearea_hrt_statusPrevious"]] * (x[["hrt_status"]] == "Previous user") +
          betas[["betaLong_densearea_hrt_statusCurrent"]] * (x[["hrt_status"]] == "Current user") +
          betas[["betaLong_densearea_hrt_statusMissing"]] * (x[["hrt_status"]] == "Missing") +
          betas[["betaLong_densearea_family_bcYes"]] * (x[["family_bc"]] == "Yes") +
          betas[["betaLong_densearea_family_bcMissing"]] * (x[["family_bc"]] == "Missing")
      )
  )
}

# Model parameters
this_N <- nrow(synth_surv_mlo)
betas <- data.frame(
  id = synth_surv_mlo$id,
  # Weibull
  shape = rep(0.7, this_N),
  betaEvent_intercept = rep(-5.5, this_N),
  # Survival
  betaEvent_bmi = rep(fmt_est$survival[[1]]$fixef[2], this_N),
  betaEvent_hrt_statusPrevious = rep(fmt_est$survival[[1]]$fixef[3], this_N),
  betaEvent_hrt_statusCurrent = rep(fmt_est$survival[[1]]$fixef[4], this_N),
  betaEvent_hrt_statusMissing = rep(fmt_est$survival[[1]]$fixef[5], this_N),
  betaEvent_family_bcYes = rep(fmt_est$survival[[1]]$fixef[6], this_N),
  betaEvent_family_bcMissing = rep(fmt_est$survival[[1]]$fixef[7], this_N),
  betaEvent_birth_times1 = rep(fmt_est$survival[[1]]$fixef[8], this_N),
  betaEvent_birth_times2 = rep(fmt_est$survival[[1]]$fixef[9], this_N),
  betaEvent_birth_times3p = rep(fmt_est$survival[[1]]$fixef[10], this_N),
  betaEvent_birth_timesMissing = rep(fmt_est$survival[[1]]$fixef[11], this_N),
  # Association
  betaEvent_assoc_fattyarea = rep(fmt_est$survival[[1]]$associations[1], this_N),
  betaEvent_assoc_densearea = rep(fmt_est$survival[[1]]$associations[2], this_N),
  # Fatty area
  betaLong_fattyarea_intercept = rep(fmt_est$marker[[1]]$fixef[1], this_N),
  betaLong_fattyarea_slope = rep(fmt_est$markers[[1]]$fixef_vary[1], this_N),
  betaLong_fattyarea_bmi = rep(fmt_est$marker[[1]]$fixef[2], this_N),
  betaLong_fattyarea_hrt_statusPrevious = rep(fmt_est$marker[[1]]$fixef[3], this_N),
  betaLong_fattyarea_hrt_statusCurrent = rep(fmt_est$marker[[1]]$fixef[4], this_N),
  betaLong_fattyarea_hrt_statusMissing = rep(fmt_est$marker[[1]]$fixef[5], this_N),
  betaLong_fattyarea_family_bcYes = rep(fmt_est$marker[[1]]$fixef[6], this_N),
  betaLong_fattyarea_family_bcMissing = rep(fmt_est$marker[[1]]$fixef[7], this_N),
  # Dense area
  betaLong_densearea_intercept = rep(fmt_est$marker[[2]]$fixef[1], this_N),
  betaLong_densearea_slope = rep(fmt_est$markers[[2]]$fixef_vary[1], this_N),
  betaLong_densearea_bmi = rep(fmt_est$marker[[2]]$fixef[2], this_N),
  betaLong_densearea_hrt_statusPrevious = rep(fmt_est$marker[[2]]$fixef[3], this_N),
  betaLong_densearea_hrt_statusCurrent = rep(fmt_est$marker[[2]]$fixef[4], this_N),
  betaLong_densearea_hrt_statusMissing = rep(fmt_est$marker[[2]]$fixef[5], this_N),
  betaLong_densearea_family_bcYes = rep(fmt_est$marker[[2]]$fixef[6], this_N),
  betaLong_densearea_family_bcMissing = rep(fmt_est$marker[[2]]$fixef[7], this_N)
) |>
  left_join(distinct(select(synth_full_mlo, id, starts_with("b0"), starts_with("b1"))), by = "id") |>
  mutate(betaLong_fattyarea_intercept = betaLong_fattyarea_intercept + b0_fa) |>
  mutate(betaLong_densearea_intercept = betaLong_densearea_intercept + b0_da) |>
  mutate(betaLong_fattyarea_slope = betaLong_fattyarea_slope + b1_fa) |>
  mutate(betaLong_densearea_sloep = betaLong_densearea_slope + b1_da) |>
  select(-starts_with("b0"), -starts_with("b1"))
# Baseline covariates
covdat <- distinct(synth_full_mlo, id, bmi, family_bc, hrt_status, birth_times)

# Then we simulate the survival times based on the
# hazard function, covariates, and true parameter values
times <- simsurv(hazard = haz, x = covdat, betas = betas, maxt = ceiling(max(data_surv_mlo$single_t - data_surv_mlo$single_t0)))
head(times)
table(times$status)
prop.table(table(times$status))
summary(times$eventtime)
hist(times$eventtime)

# Merge this with the one-row synthetic dataset
synth_surv_mlo <- synth_surv_mlo |>
  left_join(times, by = "id") |>
  mutate(single_t0 = age_inclusion_in_karma - 40) |>
  mutate(single_t = single_t0 + eventtime) |>
  rename(single_d = status) |>
  select(id, age_inclusion_in_karma, bmi, hrt_status, family_bc, birth_times, single_t0, single_t, single_d)

# Merge this information with the long dataset
synth_full_mlo <- synth_full_mlo |>
  left_join(times, by = "id") |>
  mutate(eventtime = eventtime + age_inclusion_in_karma - 40) |>
  filter(t0 < eventtime) |>
  mutate(N = n(), .by = "id") |>
  mutate(t = ifelse(n == N, NA, t)) |>
  mutate(t = ifelse(is.na(t), eventtime, t)) |>
  mutate(d = ifelse(n == N, status, 0)) |>
  select(id, n, N, age_inclusion_in_karma, bmi, hrt_status, family_bc, birth_times, fattyarea, densearea, t0, t, d)

# Export datasets
#
saveRDS(object = synth_full_mlo, file = "data/20-synth-full-mlo.RDS")
write_csv(x = synth_full_mlo, file = "data/20-synth-full-mlo.csv")
#
saveRDS(object = synth_surv_mlo, file = "data/20-synth-surv-mlo.RDS")
write_csv(x = synth_surv_mlo, file = "data/20-synth-surv-mlo.csv")
