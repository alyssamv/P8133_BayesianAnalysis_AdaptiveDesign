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
  t = c()
  for (i in s1:n1) {
    e1 = dbinom(i, n1, p) # Pr(s1 = i | p = p0)
    e2 = 0
    if (s2 - i > 0) {
      for (j in (s2 - i):(n2)) {
        e2 = e2 + dbinom(j, n2, p) # sum {Pr(s1 >= s2-i | p = p0)}
      }
    } else {
      for (j in 0:(ntotal - n1)) {
        e2 = e2 + dbinom(j, n2, p)
      }
    }
    t[i] = e1*e2
  }
  pgo.2stg = sum(t, na.rm = TRUE)
  return(pgo.2stg)
}


################################################################
# go.2stg = function to calculate probability of a "go" decision in a two-stage adaptive design
# n1 = stage 1 sample size
# s1 = minimum number of patients required to continue to stage 2
# n2 = stage 2 sample size
# s2 = minimum number of patients (both stages) to conclude efficacy
# p = "true" response probability
################################################################


