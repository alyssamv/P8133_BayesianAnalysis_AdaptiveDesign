#####################################################################
################# Functions for pipeline simulation #################
#####################################################################


# p.true = vector of true DLT probabilities for each test dose
# n.dose = number (integer) of test doses
# N = number of simulations to run
sim_3plus3 <- function(p.true, n.dose){
  results = vector(mode = "list", length = N)
  
  pt = 1 # patient on dose k
  n = 1 # patient on trial
  dose = 1 # dose k
  esc = 1 # escalation yes/no
  mtd = 0 # MTD
  
  dlts = vector(mode = 'list', length = n.dose)
  
  while (esc == 1 && dose <= n.dose && pt <= 6) {
    dlts[[dose]][pt] = rbinom(1, 1, p.true[dose])
    
    # If <3 pts enrolled on current dose and <2 DLTs, enroll next pt on same dose
    if (sum(dlts[[dose]]) <= 1 && length(dlts[[dose]]) < 3) {
      pt = pt + 1
      n = n + 1
      dose = dose
      # If 1/3 pts on current dose observe DLT, enroll next pt on same dose
    } else if (sum(dlts[[dose]]) == 1 && length(dlts[[dose]]) == 3) {
      pt = pt + 1
      n = n + 1
      dose = dose
      # If 0/3 pts on current dose observe DLT, escalate dose
    } else if (sum(dlts[[dose]]) == 0 && length(dlts[[dose]]) == 3) {
      pt = 1
      dose = dose + 1
      n = n + 1
      # If 2 patients out of 2-3 observe DLT, stop trial and declare the previous dose MTD
    } else if (sum(dlts[[dose]]) == 2 && length(dlts[[dose]]) <= 3) {
      esc = 0
      mtd = dose - 1
      # If <=1 patient out of 4-6 observe DLT, stay at same dose
    } else if (sum(dlts[[dose]]) == 1 && length(dlts[[dose]]) > 3 && length(dlts[[dose]]) < 6) {
      pt = pt + 1
      n = n + 1
      dose = dose
      # If we observe 1/6 DLTs, stop escalation and declare this dose level the MTD
    } else if (sum(dlts[[dose]]) == 1 && length(dlts[[dose]]) == 6) {
      esc = 0 # MTD
      mtd = dose
      # If we observe 2 DLTs for 4-6 pts, stop the trial and declare the previous dose the MTD
    } else if (sum(dlts[[dose]]) == 2 && length(dlts[[dose]]) > 3 && length(dlts[[dose]]) <= 6) {
      esc = 0
      mtd = dose - 1
      # In case of error
    } else {
      dlts[[dose]][pt] = NA # reassign entry as NA in case of error
    }
  }
  
  results = list(dlts = unlist(lapply(dlts, sum)),
                 nperdose = unlist(lapply(dlts, length)),
                 n = n,
                 mtd = mtd)
  return(results)
}




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
                     p2.alpha
                     ) {
  
  if (design == "3+3") {
    phase1 = sim_3plus3(p.true, n.dose)
    rp2d = phase1$mtd
  } else if (design == "crm") {
    prior = dfcrm::getprior(halfwidth = 0.05, 
                            target = targetDLT, 
                            nu = mtd.guess, 
                            nlevel = n.dose)
    firstdose = 3
    phase1 = crmsim(PI = p.true, 
                    prior = prior, 
                    target = targetDLT, 
                    n = p1.n, 
                    x0 = firstdose, 
                    nsim = 1)
    rp2d = phase1$MTD
  }
  
  
  p = p.response[[rp2d]]
  
  phase2 = phase2(p.true = p, 
                  n = p2.n, 
                  p0 = hist, 
                  alpha = p2.alpha)
  
  return(phase2)
  
}