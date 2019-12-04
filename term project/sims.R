source(file.path(getwd(), "term project/fns.R"))

target.dlt = 0.25

dlt.true = c(0.011, 0.10, 0.15, 0.25, 0.4)  # c(0.011, 0.064, 0.195, 0.330, 0.446) # From Sumithra
n.dose = length(dlt.true)

# from Sumithra
response = list(
  c(0.1, 0.2, 0.3, 0.4, 0.7), # correct dose = 4
  rep(0.5, n.dose), # CD = 1
  c(0.2, 0.5, 0.8, 0.8, 0.8), # CD = 4
  c(0.8, 0.6, 0.4, 0.3, 0.2), # CD = 1
  c(0.3, 0.8, 0.65, 0.45, 0.3) # CD = 3
)

N = 1e3 # number of simulations per efficacy scenario
p.hist = 0.4 # historical response rate 

mtd.true = tail(which(dlt.true <= target.dlt), n = 1)
correct.dose = c(4, 1, 3, 1, 2)

# Figure 1
par(mfrow = c(3, 2))
for (i in 1:length(response)) {
  plot(x = 1:n.dose, y = dlt.true, type = 'b',
       ylim = c(0, 1),
       col = "red",
       xlab = "",
       ylab = "")
  lines(x = 1:n.dose, y = response[[i]], type = 'b',
        ylim = c(0, 1),
        col = "green")
  abline(h = target.dlt, lty = 2)
  title(main = paste("Scenario", i, sep = " "),
        xlab = "Dose",
        ylab = "Probability")
  points(x = mtd.true, y = dlt.true[[mtd.true]], pch = 17, cex = 1.2) # MTD
  points(x = correct.dose[[i]], response[[i]][correct.dose[[i]]], pch = 8, cex = 1.2) # correct dose
}
plot(1, type = "n", xlab = "", ylab = "", xlim = c(0, 10), ylim = c(0, 10), axes = FALSE,
     xaxt = "n", yaxt = "n")
legend(x = 0, y = 10, legend = c("Dose-toxicity", "Dose-efficacy", "Target toxicity", "MTD", "Correct dose"),
       col = c("red", "green", "black", "black", "black"),
       lty = c(1, 1, 2, NA, NA), pch = c(1, 1, NA, 17, 8),
       cex = 1.2, bty = "n")


############################################ Simulations ############################################

############ 3 + 3 #############

three.sims = lapply(response, function(i){
  lapply(1:N, function(j) {
    pipeline(p.true = dlt.true,
             n.dose = n.dose,
             p.response = i,
             p1.design = "3+3",
             hist = p.hist,
             p1.n = 33,
             p2.n = 40,
             p2.alpha = 0.05)
  })
})


############# CRM ##############

crm.sims = lapply(response, function(i){
  lapply(1:N, function(j) {
    pipeline(p.true = dlt.true,
             n.dose = n.dose,
             p.response = i,
             p1.design = "crm",
             hist = p.hist,
             targetDLT = 0.2,
             p1.n = 33,
             p2.n = 40,
             p2.alpha = 0.05)
  })
})


############ EffTox #############
library(foreach)
library(doParallel)

numCores = detectCores()
registerDoParallel(numCores)

start = Sys.time()
#vector(mode = "list", length = length(response))
efftox.sims = foreach (yyy = 1:length(response)) %dopar% {
  pipeline(p.true = dlt.true,
           n.dose = n.dose,
           p.response = response[[yyy]],
           p1.design = "efftox",
           real.doses = c(1, 2, 4, 6.6, 10), ## need to get a reasonable vector from paper
           hist = p.hist,
           targetDLT = 0.2,
           p1.n = 33,
           p2.n = 40,
           p2.alpha = 0.05,
           n.sims = N)
}
end = Sys.time()
end - start
  
save.image(file = file.path(getwd(), "term project/sims.RData"))

