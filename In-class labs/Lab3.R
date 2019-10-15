library(dfcrm)
p = 0.25
p0 = c(0.08, 0.14, 0.25, 0.37, 0.52) # skeleton
x = c(3, 4, 2, 2, 3, 2)
y = c(0, 1, 0, 0, 1, 0)
crmfit = crm(prior = p0, target = p, tox = y, level = x)
crmfit$mtd
crmfit$ptox # posterior dose-tox curve
plot(crmfit, ask = TRUE)


### 2a
PI = c(0.002, 0.01, 0.051, 0.23, 0.62)
simfit = crmsim(PI, prior = p0, target = 0.25, n = 24, x0 = 3, nsim = 50)
plot(simfit)
simfit$MTD


### 2b
p0 = c(0.05, 0.15, 0.25, 0.50, 0.70)
simfit2 = crmsim(PI, prior = p0, target = 0.25, n = 24, x0 = 3, nsim = 50)
plot(simfit2)
simfit2$MTD


### 2c
p0 = c(0.05, 0.25, 0.35, 0.50, 0.70)
simfit3 = crmsim(PI, prior = p0, target = 0.25, n = 24, x0 = which(p0 == p), nsim = 50)
plot(simfit3)
simfit3$MTD

p0 = c(0.05, 0.10, 0.15, 0.25, 0.70)
simfit4 = crmsim(PI, prior = p0, target = 0.25, n = 24, x0 = which(p0 == p), nsim = 50)
plot(simfit4)
simfit4$MTD

# when the skeleton assigns the MTD to a smaller dose, it seems to perform better

### 2d
p0 = c(0.08, 0.14, 0.25, 0.37, 0.52) # c(0.05, 0.25, 0.35, 0.50, 0.70)
PI = c(0.002, 0.010, 0.05, 0.12, 0.25)
simfit = crmsim(PI, prior = p0, target = 0.25, n = 24, x0 = which(p0 == p), nsim = 50)
plot(simfit)
simfit$MTD

# changing the truth affects whether the specified skeleton gives good performance


#### 3
p0 = getprior(halfwidth = 0.06, target = 0.25, nu = 3, nlevel = 5)

PI = c(0.002, 0.01, 0.051, 0.23, 0.62)
simfit = crmsim(PI, prior = p0, target = 0.25, n = 24, x0 = which(p0 == p), nsim = 50)
plot(simfit)
simfit$MTD

PI = c(0.002, 0.010, 0.05, 0.12, 0.25)
simfit = crmsim(PI, prior = p0, target = 0.25, n = 24, x0 = which(p0 == p), nsim = 50)
plot(simfit)
simfit$MTD


### 3b

p0 = getprior(halfwidth = 0.02, target = 0.25, nu = 3, nlevel = 5)
plot(1:5, p0)

PI = c(0.002, 0.01, 0.051, 0.23, 0.62)
simfit = crmsim(PI, prior = p0, target = 0.25, n = 24, x0 = which(p0 == p), nsim = 50)
plot(simfit)
simfit$MTD

PI = c(0.002, 0.010, 0.05, 0.12, 0.25)
simfit = crmsim(PI, prior = p0, target = 0.25, n = 24, x0 = which(p0 == p), nsim = 50)
plot(simfit)
simfit$MTD


p0 = getprior(halfwidth = 0.1, target = 0.25, nu = 3, nlevel = 5)
plot(1:5, p0)

PI = c(0.002, 0.01, 0.051, 0.23, 0.62)
simfit = crmsim(PI, prior = p0, target = 0.25, n = 24, x0 = which(p0 == p), nsim = 50)
plot(simfit)
simfit$MTD

PI = c(0.002, 0.010, 0.05, 0.12, 0.25)
simfit = crmsim(PI, prior = p0, target = 0.25, n = 24, x0 = which(p0 == p), nsim = 50)
plot(simfit)
simfit$MTD

