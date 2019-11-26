#### Functions for Bayesian analysis and adaptive designs

################################################################
# go.2stg = function to calculate probability of a "go" decision in a two-stage adaptive design
# n1 = stage 1 sample size
# s1 = minimum number of patients required to continue to stage 2
# n2 = stage 2 sample size
# s2 = minimum number of patients (both stages) to conclude efficacy
# p = "true" response probability
################################################################

go.2stg <- function(n1, s1, n2, s2, p) {
  f = c()
  for (i in s1:n1) {
    f[i] = dbinom(i, n1, p) * (1 - pbinom(s2 - i - 1, n2, p))
  }
  pgo = sum(f, na.rm = T)
  return(pgo)
}


################################################################
# nogo.Bayes = function to calculate Bayesian formulation of a "no-go" decision
# pe = simulated exp response probabilities
# ps = simulated control response probabilities
# delta = effect margin
# N = number of simulated probabilities 
################################################################

nogo.Bayes <- function(pe, ps, delta, N) {
  sum(pe > (ps + d)) / N
}
