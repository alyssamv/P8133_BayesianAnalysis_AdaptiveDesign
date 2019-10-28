data {
  int<lower=0> N; // Total number of observations for all drugs
  int y[N]; // Observed outcome of response
  int<lower=1> D; // Number of drugs
  int<lower=1> M; // Number of biomarkers
  int<lower=1, upper=M> m[N]; // Index for biomarkers 
  row_vector[D] X[N]; // Design matrix
  int S[N]; // Number of trials
}
parameters {
  real mu[D];
  real<lower=0> sigma[D];
  vector[D] beta[M];
}
model {
  for (d in 1:D) {
    mu[d] ~ normal(0, 100);
    for (i in 1:M)
      beta[i,d] ~ normal(mu[d], sigma[d]);
  }
  for (n in 1:N)
    y[n] ~ binomial(S[n], inv_logit(X[n] * beta[m[n]]));
}
