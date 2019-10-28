data {
  int<lower=0> N; // Total number of observations for all drugs
  int y[N]; // Observed outcome of response
  int S[N]; // Number of trials
}
parameters {
  vector<lower=0, upper=1>[N] theta;
}
model {
  y ~ binomial(S, theta);
}
