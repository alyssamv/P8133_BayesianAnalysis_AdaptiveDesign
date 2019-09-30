

N = 10000
a = 0.5
b = 1.5 # as 
p = rbeta(N, a, b)
mean(p); var(p)
hist(p, prob = T, main = "pdf of beta(a, b)")


post = function(x2, x1 = 1, a = 0.5, b = 1.5) {
  gamma(a + b + x1)/(gamma(a + x1)*gamma(1 - x1 + b)) * (gamma(x2 + x1 + a)*gamma(b - x1 - x2 + 2))/gamma(a + b + 2)
}

post(x2 = 1)
post(x2 = 0)

###### 2

x1 = 1
a = 0.5 + x1
b = 1.5
pE = rbeta(N, a, b)
x2 = rbinom(N, 1, pE)
table(x2)


###### 3
n = 20
x10 = 10
a = 0.5 + x10
b = 1.5 + n - x10

pE = rbeta(N, a, b)

mean(rbinom(N, 62 - n, pE) >= (22 - x10))

