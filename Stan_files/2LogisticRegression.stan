data {
  int<lower=0> N; // Total number of observations for all drugs
  int y[N]; // Observed outcome of response
  int D; // Dimension of design matrix (Number of drugs)
  matrix [N, D] X; // Design matrix X
  int S[N]; // Number of trials
}
parameters {
  vector[D] beta;
}
transformed parameters {
  vector[D] eta;
  vector[D] theta;
  eta = X * beta;
  for(i in 1:D)
  theta[i] = exp(eta[i])/(1+exp(eta[i]));
}
model {
  //Prior
    beta[2] ~ normal(1, 1); // non-informative prior
  //Likelihood
  y ~ binomial_logit(S, eta);
}
