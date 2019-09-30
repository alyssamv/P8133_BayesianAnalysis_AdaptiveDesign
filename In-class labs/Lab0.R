### In-class example

p0 = 0.25
p1 = 0.4
alpha = 0.05
power = 0.8

z1 = 6
n1 = 20
z2 = 24
ntotal = 71

## Type I error rate?
(1 - pbinom(z1 - 1, n1, p0)) * (1 - pbinom(z2 - 1, ntotal, p0)) # 0.023

## Power?
(1 - pbinom(z1 - 1, n1, p1)) * (1 - pbinom(z2 - 1, ntotal, p1)) # 0.772

## Expected sample size under the null
(pbinom(z1 - 2, n1, p0) * n1) + ((1 - pbinom(z1 - 2, n1, p0)) * ntotal) # 40




# ### Interim
# 
# n = 20
# 
# ## Type I error
# i = 1
# t1 = 1
# while (t1 > alpha) {
#   z = i
#   t1 = 1 - pbinom(z - 1, n, p0)
#   i = i + 1
# }
# i - 1
# 
# ## Power
# j = 1
# pow = 0
# while (pow < power) {
#   z = j
#   pow = pbinom(z - 1, n, p1)
#   j = j + 1
# }
# j - 1
# 
# ### Total
# n = 51
# 
# i = 1
# t1 = 1
# while (t1 > alpha) {
#   z = i
#   t1 = 1 - pbinom(z - 1, n, p0)
#   i = i + 1
# }
# i - 1