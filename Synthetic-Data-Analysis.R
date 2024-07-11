# Script to analyse the synthetic dataset using:
# - Univariate joint model for non-dense area, adjusting for baseline dense area;
# - Bivariate joint model for non-dense and dense area

# Required packages, they can be installed from CRAN
library(conflicted)
library(tidyverse)
library(VAJointSurv)

# Deal with functions name collision
conflicts_prefer(dplyr::filter)

# Read in synthetic data
data_full_mlo <- readRDS(file = "data/20-synth-full-mlo.RDS")
data_surv_mlo <- readRDS(file = "data/20-synth-surv-mlo.RDS")

# We add baseline dense area to the survival dataset
data_surv_mlo <- data_surv_mlo |>
  left_join(
    data_full_mlo |>
      filter(n == 1) |>
      select(id, densearea) |>
      rename(baseline_densearea = densearea),
    by = "id"
  )

# Joint model #1
# This is a univariate joint model with non-dense area as the longitudinal marker,
# and adjusting for dense area at baseline in the survival submodel.

# Marker:
m1_fa_marker <- marker_term(
  formula = sqrt(fattyarea) ~ bmi + hrt_status + family_bc,
  id = id,
  data = data_full_mlo,
  time_fixef = ns_term(x = t0, df = 3),
  time_rng = poly_term(x = t0, degree = 1, intercept = TRUE, raw = TRUE) # Random linear slope + intercept
)

# Survival:
m1_ssub <- surv_term(
  Surv(single_t0, single_t, single_d) ~ bmi + hrt_status + family_bc + birth_times + baseline_densearea,
  id = id,
  data = data_surv_mlo,
  time_fixef = ns_term(single_t, df = 3),
  delayed = single_t0 > 0
)

# Joint model:
m1_comp_obj <- joint_ms_ptr(
  markers = list(m1_fa_marker),
  survival_terms = list(m1_ssub),
  max_threads = parallel::detectCores()
)

# Initialise the joint model:
m1_start_val <- joint_ms_start_val(object = m1_comp_obj)

# Fit the joint model:
m1_opt_out <- joint_ms_opt(object = m1_comp_obj, par = m1_start_val)

# Format the estimates:
m1_fmt_est <- joint_ms_format(object = m1_comp_obj, par = m1_opt_out$par)

# Calculate the Hessian:
m1_hessian <- joint_ms_hess(object = m1_comp_obj, par = m1_opt_out$par)

# Create the variance-covariance matrix of the model parameters:
m1_vcov <- solve(m1_hessian$hessian)
m1_vcov <- as.matrix(m1_vcov)
rownames(m1_vcov) <- unlist(m1_comp_obj$param_names$param_names)
colnames(m1_vcov) <- unlist(m1_comp_obj$param_names$param_names)

# The estimated model coefficients are:
print(m1_fmt_est)

# Define a range of values to plot for:
x_range <- c(0, max(data_surv_mlo$single_t0[data_surv_mlo$single_d == 1]))

# Plot the mean marker level:
plot_marker(
  time_fixef = m1_comp_obj$markers[[1]]$time_fixef,
  time_rng = m1_comp_obj$markers[[1]]$time_rng,
  fixef_vary = m1_fmt_est$markers[[1]]$fixef_vary,
  x_range = x_range,
  vcov_vary = m1_fmt_est$vcov$vcov_vary,
  xlab = "Time",
  ylab = "Marker mean curve, sqrt(non-dense area)"
)

# Plot the baseline hazard:
plot_surv(
  time_fixef = m1_comp_obj$survival_terms[[1]]$time_fixef,
  time_rng = list(m1_comp_obj$markers[[1]]$time_rng),
  x_range = x_range,
  fixef_vary = m1_fmt_est$survival[[1]]$fixef_vary,
  vcov_vary = m1_fmt_est$vcov$vcov_vary, frailty_var = 0,
  associations = m1_fmt_est$survival[[1]]$associations,
  ders = m1_comp_obj$survival_terms[[1]]$ders,
  log_hazard_shift = m1_fmt_est$survival[[1]]$fixef[1]
)

# Joint model #2
# This is a bivariate joint model with non-dense and dense area as the longitudinal markers.

# Marker 1, non-dense area:
m2_fa_marker <- marker_term(
  formula = sqrt(fattyarea) ~ bmi + hrt_status + family_bc,
  id = id,
  data = data_full_mlo,
  time_fixef = ns_term(x = t0, df = 3),
  time_rng = poly_term(x = t0, degree = 1, intercept = TRUE, raw = TRUE) # Random linear slope + intercept
)

# Marker 2, dense area:
m2_da_marker <- marker_term(
  formula = sqrt(densearea) ~ bmi + hrt_status + family_bc,
  id = id,
  data = data_full_mlo,
  time_fixef = ns_term(x = t0, df = 3),
  time_rng = poly_term(x = t0, degree = 1, intercept = TRUE, raw = TRUE) # Random linear slope + intercept
)

# Survival:
m2_ssub <- surv_term(
  Surv(single_t0, single_t, single_d) ~ bmi + hrt_status + family_bc + birth_times,
  id = id,
  data = data_surv_mlo,
  time_fixef = ns_term(single_t, df = 3),
  delayed = single_t0 > 0
)

# Joint model:
m2_comp_obj <- joint_ms_ptr(
  markers = list(m2_fa_marker, m2_da_marker),
  survival_terms = list(m2_ssub),
  max_threads = parallel::detectCores()
)

# Initialise the joint model:
m2_start_val <- joint_ms_start_val(object = m2_comp_obj)

# Fit the joint model:
m2_opt_out <- joint_ms_opt(object = m2_comp_obj, par = m2_start_val)

# Format the estimates:
m2_fmt_est <- joint_ms_format(object = m2_comp_obj, par = m2_opt_out$par)

# Calculate the Hessian at the optimum:
m2_hessian <- joint_ms_hess(object = m2_comp_obj, par = m2_opt_out$par)

# Create the variance-covariance matrix of the model parameters:
m2_vcov <- solve(m2_hessian$hessian)
m2_vcov <- as.matrix(m2_vcov)
rownames(m2_vcov) <- unlist(m2_comp_obj$param_names$param_names)
colnames(m2_vcov) <- unlist(m2_comp_obj$param_names$param_names)

# The estimated model coefficients are:
print(m2_fmt_est)

# Plot the mean markers levels:
# Non-dense area:
plot_marker(
  time_fixef = m2_comp_obj$markers[[1]]$time_fixef,
  time_rng = m2_comp_obj$markers[[1]]$time_rng,
  fixef_vary = m2_fmt_est$markers[[1]]$fixef_vary,
  x_range = x_range,
  vcov_vary = m2_fmt_est$vcov$vcov_vary[1:2, 1:2],
  xlab = "Time",
  ylab = "Marker mean curve, sqrt(non-dense area)"
)
# Dense area:
plot_marker(
  time_fixef = m2_comp_obj$markers[[2]]$time_fixef,
  time_rng = m2_comp_obj$markers[[2]]$time_rng,
  fixef_vary = m2_fmt_est$markers[[2]]$fixef_vary,
  x_range = x_range,
  vcov_vary = m2_fmt_est$vcov$vcov_vary[3:4, 3:4],
  xlab = "Time",
  ylab = "Marker mean curve, sqrt(dense area)"
)

# Plot the baseline hazard:
plot_surv(
  time_fixef = m2_comp_obj$survival_terms[[1]]$time_fixef,
  time_rng = lapply(m2_comp_obj$markers, `[[`, "time_rng"),
  x_range = x_range,
  fixef_vary = m2_fmt_est$survival[[1]]$fixef_vary,
  vcov_vary = m2_fmt_est$vcov$vcov_vary, frailty_var = 0,
  associations = m2_fmt_est$survival[[1]]$associations,
  ders = m2_comp_obj$survival_terms[[1]]$ders,
  log_hazard_shift = m2_fmt_est$survival[[1]]$fixef[1]
)
