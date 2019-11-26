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
                     p1.design = c("3+3", "crm", "efftox"), # design selection for phase I
                     hist, # historical response rate
                     p1.n,
                     p2.n, # phase II sample size
                     p2.alpha,
                     seed = NULL
                     ) {
  
  ### Phase I trial 
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
                           nsim = 1)
    n1 = p1.n
    rp2d = phase1$MTD
  }
  
  if (!(rp2d %in% 1:n.dose)) {
    return(c("p1.n" = NA, 
             "rp2d" = rp2d, 
             "rp2d.eff" = NA, 
             NA))
  } else {
    # Efficacy for RP2D
    p = p.response[[rp2d]]
    
    ### Phase II SAT
    phase2 = phase2(p.true = p, 
                    n = p2.n, 
                    p0 = hist, 
                    alpha = p2.alpha)
    
    
    return(c("p1.n" = n1, 
             "rp2d" = rp2d, 
             "rp2d.eff" = p, 
             phase2))
  }
  
}
