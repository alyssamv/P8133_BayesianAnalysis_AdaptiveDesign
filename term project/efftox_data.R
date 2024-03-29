library(trialr)

p <- efftox_solve_p(eff0 = 0.5, tox1 = 0.65, eff_star = 0.7, tox_star = 0.25)

?efftox_params
dat <- list(
  num_doses = 5,
  real_doses = c(1, 2, 4, 6.6, 10), # vector of the actual doses being administered (suppose )
  efficacy_hurdle = 0.5,
  toxicity_hurdle = 0.3,
  p_e = 0.1,
  p_t = 0.1,
  p = p,
  eff0 = 0.5,
  tox1 = 0.65,
  eff_star = 0.7,
  tox_star = 0.25,
  
  alpha_mean = -7.9593, alpha_sd = 3.5487,
  beta_mean = 1.5482, beta_sd = 3.5018,
  gamma_mean = 0.7367, gamma_sd = 2.5423,
  zeta_mean = 3.4181, zeta_sd = 2.4406,
  eta_mean = 0, eta_sd = 0.2,
  psi_mean = 0, psi_sd = 1,
  
  doses = c(),
  tox   = c(),
  eff   = c(),
  num_patients = 0
)


set.seed(1)
sims = efftox_simulate(dat, num_sims = 1, first_dose = 1, 
                       true_eff = c(0.20, 0.40, 0.60, 0.80, 0.90),
                       true_tox = c(0.05, 0.10, 0.15, 0.20, 0.40),
                       cohort_sizes = rep(3, 13))

