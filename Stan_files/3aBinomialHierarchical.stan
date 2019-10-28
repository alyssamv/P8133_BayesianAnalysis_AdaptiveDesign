data {
  int<lower=0> N; // Total number of observations for all drugs
  int y[N]; // Observed outcome of response
  int S[N]; // Number of trials
}
parameters {
  vector<lower=0, upper=1>[N] theta;
  real<lower=1> a;
  real<lower=1> b;
}
model {
  //Prior distirbution
  theta ~ beta(a, b);
  
  //Hyperprior distribution
  a ~ exponential(1);
  b ~ exponential(1);
  
  //Likehood function
  for (n in 1:N){
    y[n] ~ binomial(S[n], theta[n]);
  }
}
