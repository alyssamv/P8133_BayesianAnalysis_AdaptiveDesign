source(file.path(getwd(), "term project/fns.R"))

dlt.true = c(0.017, 0.043, 0.10, 0.22, 0.41)
response = c(0, 0.1, 0.3, 0.35, 0.35)
n.dose = 5
N = 1e3
p = 0.4

plot(x = 1:n.dose, y = dlt.true, type = 'b',
     ylim = c(0, 0.6),
     col = "red")
lines(x = 1:n.dose, y = response, type = 'b',
     ylim = c(0, 0.6),
     col = "green")

sims = sapply(1:N, function(i) {pipeline(p.true = dlt.true, 
         n.dose = n.dose, 
         p.response = response, 
         p1.design = "3+3",
         hist = 0.2, 
         p1.n = 33, 
         p2.n = 40, 
         p2.alpha = 0.05) })


