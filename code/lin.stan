data {
  int<lower=0> N;
  vector[N] att;
  vector[N] afd;
}
parameters {
  vector[2] beta;
  real<lower=0> sigma;
}
model {
  att ~ normal(beta[1] + beta[2] * afd, sigma);
}

