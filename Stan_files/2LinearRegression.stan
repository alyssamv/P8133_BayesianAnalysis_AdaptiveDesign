data {
  int<lower=0> N;
  int<lower=1> K; // Dimension of design matrix
  matrix[N, K] x; // Design matrix
  vector[N] y;
}
parameters {
  vector[K] beta; 
  real<lower=0> sigma;
}
model {
  // Prior
  beta[1] ~ normal(10, 0.1);
  beta[2] ~ normal(6, 0.1);
  
  //Likelihood
  y ~ normal(x * beta, sigma);
}
