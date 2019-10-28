data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}
parameters {
  real beta0;
  real beta1;
  real<lower=0> sigma;
}
model {
  y ~ normal(beta0 + beta1 * x, sigma);
}
