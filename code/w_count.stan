data {
  int<lower=0> N;
  int<lower=0> N_edges;
  int<lower=1, upper=N> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=N> node2[N_edges];  // and node1[i] < node2[i]
  int att[N];
  vector[N] afd;
}
parameters {
  vector[N] phi;
  real<lower=0,upper=1> lambda;  // could be used to restrict lambda
  // real<lower=-1,upper=1> lambda;
  vector[1] beta; // vector[2] beta for 2 betas... 
}
model {
target += -0.5 * dot_self(phi[node1] - phi[node2]); // soft sum-to-zero constraint on phi  
sum(phi) ~ normal(0, 0.01 * N);  // equivalent to mean(phi) ~ normal(0,0.01);
att ~ poisson_log(beta[1] * afd + lambda * phi[N]);  // lambda estimates spatial influence 
// lambda ~ ... // choose priors or diffuse default priors are used 
}
