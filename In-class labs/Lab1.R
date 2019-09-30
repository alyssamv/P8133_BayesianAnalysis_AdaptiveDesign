### 1. Learn beta

N = 10000
a = 0.5
b = 0.5 # as 
p = rbeta(N, a, b)
mean(p); var(p)
hist(p, prob = T, main = "pdf of beta(a, b)")

# The prior probability Pr(p > 0.25):
a = 0.5
b = 1.5
1 - pbeta(0.25, a, b) # 0.391


### 2. Posterior computation
a = 0.5
b = 1.5
n = 62
s = 22

# Posterior probability Pr(p > 0.25 | data)
a_cond = a + s
b_cond = b + n - s
1 - pbeta(0.25, a_cond, b_cond) # 0.961


### 3. Uncertainty
a.s = 25
b.s = 75
d = 0.15
alpha = 0.5

a_cond = a + s
b_cond = b + n - s

ps = rbeta(N, a.s, b.s)
pe = rbeta(N, a_cond, b_cond)

par(mfrow = c(2, 1))
hist(ps, xlim = c(0, 1))
hist(pe, xlim = c(0, 1))

t = sum(pe > (ps + d)) / N; t # ~0.25 < alpha so conclude no-go

# Find d and alpha s.t. 22/62 is the go/no-go boundary
s = 22
a_cond = a + s
b_cond = b + n - s
pe = rbeta(N, a_cond, b_cond)

alpha = seq(0.01, 0.8, 0.01)

t = sum(pe > (ps + d)) / N; t > alpha # 0.25 -> no-go
# take only those values of alpha for which s=22 lends a "go" decision

alpha = alpha[(which((t > alpha) == TRUE))]; alpha 
t1 = sum(pe > (ps + d)) / N; t1 > alpha # 0.25 -> no-go

s = 21
a_cond = a + s
b_cond = b + n - s

pe = rbeta(N, a_cond, b_cond)
t2 = sum(pe > (ps + d)) / N; t2 > alpha # 0.32


## Verify that the selected alpha makes the decision boundary s = 22
index = intersect(which((t1 > alpha) == TRUE), which((t2 > alpha) == FALSE))
alpha.test = alpha[index] # new alphas

# This should be a "go" decision
s = 22
a_cond = a + s
b_cond = b + n - s
pe = rbeta(N, a_cond, b_cond)
t1 = sum(pe > (ps + d)) / N; t1 > alpha.test # TRUE -> go

# This should be a "no-go" decision
s = 21
a_cond = a + s
b_cond = b + n - s
pe = rbeta(N, a_cond, b_cond)
t2 = sum(pe > (ps + d)) / N; t2 > alpha.test # FALSE -> no-go

