source(file.path(getwd(), "term project/fns.R"))

target = 0.25

n.p1 = 33
n.p2 = 48

dose.amount = c(1, 2, 4, 6.6, 10)
dlt.true = c(0.011, 0.08, 0.15, 0.25, 0.4)  # c(0.011, 0.064, 0.195, 0.330, 0.446) # From Sumithra
n.dose = length(dlt.true)
mtd.true = tail(which(dlt.true <= target), n = 1)

# Response curves when drug is NOT effective
response.null = list(
  c(0.1, 0.2, 0.25, 0.35, 0.4), # correct dose = 4
  rep(0.45, n.dose), # CD = 1
  c(0.1, 0.2, 0.35, 0.35, 0.35), # CD = 4
  c(0.45, 0.35, 0.3, 0.25, 0.1), # CD = 1
  c(0.3, 0.45, 0.35, 0.25, 0.1) # CD = 3
)

# Response curves when drug is effective
response.alt = list(
  c(0.1, 0.2, 0.45, 0.6, 0.7), # correct dose = 4
  rep(0.5, n.dose), # CD = 1
  c(0.2, 0.5, 0.8, 0.8, 0.8), # CD = 4
  c(0.8, 0.6, 0.4, 0.3, 0.2), # CD = 1
  c(0.3, 0.8, 0.65, 0.45, 0.3) # CD = 3
)

N = 1e3 # number of simulations per efficacy scenario
p.hist = 0.4 # historical response rate 

correct.dose = c(4, 1, 3, 1, 2)


############################################ Simulations ############################################

############ 3 + 3 #############

three.sims.null = lapply(response.null, function(i){
  lapply(1:N, function(j) {
    pipeline(p.true = dlt.true,
             n.dose = n.dose,
             p.response = i,
             p1.design = "3+3",
             hist = p.hist,
             p1.n = n.p1,
             p2.n = n.p2,
             p2.alpha = 0.05)
  })
})


############# CRM ##############

crm.sims.null = lapply(response.null, function(i){
  lapply(1:N, function(j) {
    pipeline(p.true = dlt.true,
             n.dose = n.dose,
             p.response = i,
             p1.design = "crm",
             hist = p.hist,
             targetDLT = target,
             p1.n = n.p1,
             p2.n = n.p2,
             p2.alpha = 0.05)
  })
})


############ EffTox #############
library(foreach)
library(doParallel)

numCores = detectCores()
registerDoParallel(numCores)

start = Sys.time()
efftox.sims.null = foreach(yyy = 1:length(response.null)) %dopar% {
  pipeline(p.true = dlt.true,
           n.dose = n.dose,
           p.response = response.null[[yyy]],
           p1.design = "efftox",
           real.doses = dose.amount, ## need to get a reasonable vector from paper
           hist = p.hist,
           targetDLT = target,
           p1.n = n.p1,
           p2.n = n.p2,
           p2.alpha = 0.05,
           n.sims = N)
}
end = Sys.time()
end - start
  
save.image(file = file.path(getwd(), "term project/sims_null.RData"))

