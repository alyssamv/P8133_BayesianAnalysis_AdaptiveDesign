#####################################################################
################# Functions for pipeline simulation #################
#####################################################################

phase2 <- function(p.true, n, p0, alpha = 0.05) {
  
  r = rbinom(n, 1, p.true)
  obs = sum(r)
  
  # 1-sample exact binomial test 
  test = binom.test(obs, n, p = p0, alternative = "greater")
  
  return(list(effect = mean(r) - p0,
              p.value = test$p.value))
  
}




pipeline <- function(p.true = NULL, # vector of DLT rates for each dose
                     n.dose = NULL, # number of test doses
                     p.response, # vector of response probabilities for each dose level
                     targetDLT = 0.1,
                     mtd.guess = 3,
                     real.doses = c(), # vector of real dose amounts (eg in mg)
                     efficacy_hurdle = 0.5, # minimum acceptable efficacy
                     toxicity_hurdle = 0.3, # maximum acceptable toxicity
                     eff0 = 0.5, # eff required when toxicity impossible
                     tox1 = 0.65, # maximum toxicity permitted when efficacy is guaranteed
                     eff_star = 0.7, # ?
                     tox_star = 0.25, # ?
                     p1.design = c("3+3", "crm", "efftox"), # design selection for phase I
                     hist, # historical response rate
                     p1.n,
                     p2.n, # phase II sample size
                     p2.alpha,
                     n.sims = 1,
                     seed = NULL
) {
  
  ### Phase I trial 
  if (p1.design %in% c("3+3", "crm")) {
    if (p1.design == "3+3") {
      phase1 = UBCRM::sim3p3(truerate = p.true,
                             seed = seed)
      n1 = sum(phase1$data$npt)
      rp2d = phase1$mtd
    } else if (p1.design == "crm") {
      prior = dfcrm::getprior(halfwidth = 0.05, 
                              target = targetDLT, 
                              nu = mtd.guess, 
                              nlevel = n.dose)
      firstdose = 3
      phase1 = dfcrm::crmsim(PI = p.true, 
                             prior = prior, 
                             target = targetDLT, 
                             n = p1.n, 
                             x0 = firstdose, 
                             nsim = 1,
                             seed = NULL)
      n1 = p1.n
      rp2d = phase1$MTD
    }
    
    if (!(rp2d %in% 0:n.dose)) {
      return(list("p1.n" = n1, 
                  "rp2d" = "Dose not selected", 
                  "rp2d.eff" = NA, 
                  "phase2" = NA))
    } else {
      # Efficacy for RP2D
      p = p.response[[rp2d]]
      
      ### Phase II SAT
      phase2 = phase2(p.true = p, 
                      n = p2.n, 
                      p0 = hist, 
                      alpha = p2.alpha)
      
      
      return(list("p1.n" = n1, 
                  "rp2d" = rp2d, 
                  "rp2d.eff" = p, 
                  "phase2" = phase2))
    }
    
    
    
  } else if (p1.design == "efftox") {
    dat <- trialr::efftox_parameters_demo()
    #dat$real_doses = real.doses
    #dat$num_doses = n.dose
    
    p1 = trialr::efftox_simulate(dat, num_sims = n.sims, first_dose = 1, 
                                 true_eff = p.response,
                                 true_tox = p.true,
                                 cohort_sizes = rep(3, p1.n/3))
    
    
    results = vector(mode = "list", length = n.sims)
    
    for (kk in 1:n.sims) {
      rp2d = p1$recommended_dose[[kk]]
      n1 = length(p1$efficacies[[kk]])
      
      if (is.na(rp2d)) {
        results[[kk]] = list("p1.n" = n1, 
                             "rp2d" = "Dose not selected", 
                             "rp2d.eff" = NA, 
                             "phase2" = NA)
      } else {
        
        p = p.response[[rp2d]]
        
        phase2 = phase2(p.true = p, 
                        n = p2.n, 
                        p0 = hist, 
                        alpha = p2.alpha)
        
        results[[kk]] = list("p1.n" = n1, 
                             "rp2d" = rp2d, 
                             "rp2d.eff" = p, 
                             "phase2" = phase2)
      }
    }
    return(results)
  }
  

}
