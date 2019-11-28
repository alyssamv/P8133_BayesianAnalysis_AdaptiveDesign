source(file.path(getwd(), "term project/fns.R"))

dlt.true = c(0.008, 0.011, 0.064, 0.195, 0.330, 0.446)
n.dose = 6

response = list(
  c(0.1, 0.2, 0.3, 0.4, 0.7, 0.9), # correct dose = 4
  rep(0.5, n.dose), # CD = 1
  c(0.2, 0.5, 0.8, 0.8, 0.8, 0.8), # CD = 4
  c(0.8, 0.6, 0.4, 0.3, 0.2, 0.1), # CD = 1
  c(0.3, 0.6, 0.8, 0.6, 0.3, 0.1) # CD = 3
)

N = 1e3 # number of simulations per efficacy scenario
p.hist = 0.4 # historical response rate 

target.dlt = 0.2
mtd.true = 4
correct.dose = c(4, 1, 4, 1, 3)

# plot dose-tox-eff
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


### Simulations 
three = lapply(1:N, function(i) { 
  pipeline(p.true = dlt.true, 
           n.dose = n.dose, 
           p.response = response[[1]], 
           p1.design = "3+3",
           hist = p.hist, 
           p1.n = 33, 
           p2.n = 40, 
           p2.alpha = 0.05) 
}
)

crm = lapply(1:N, function(i) { 
  pipeline(p.true = dlt.true, 
           n.dose = n.dose, 
           p.response = response[[1]], 
           p1.design = "crm",
           hist = p.hist, 
           targetDLT = 0.2,
           p1.n = 33, 
           p2.n = 40, 
           p2.alpha = 0.05) 
}
)

efftox = lapply(1:N, function(i) { 
  pipeline(p.true = dlt.true, 
           n.dose = n.dose, 
           p.response = response[[1]], 
           p1.design = "efftox",
           hist = p.hist, 
           targetDLT = 0.2,
           p1.n = 33, 
           p2.n = 40, 
           p2.alpha = 0.05) 
}
)
