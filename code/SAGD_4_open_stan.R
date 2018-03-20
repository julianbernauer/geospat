###############################################################
#CDSS Course "Spatial Analysis of Geographic Data", Spring 2018
#Julian Bernauer, Data and Methods Unit, MZES
#Session 4, March 20, 2018 
#Open-ended spatial analysis in R 

###
# Going Bayesian - Working up to a spatial multilevel count model in Stan 

library(maptools) 
library(spdep)
library(sp) 
library(rgdal)
library(RColorBrewer) 
library(classInt)
library(rjags)
library(R2WinBUGS)
library(rstan)
library(raster)


#setwd("...")
#setwd("U:/Service MZES/Geodata/workshop cdss/code")
setwd("C:/Users/Dr. J/Desktop/code")

# data with variable for arson / assaults already included 
load("distr_btw17_arsas.Rdata")


###
# a few variables 

# arson and assaults 
att <- distr_btw17_arsas$attacks_arsas 

#afd
afd <- (distr_btw17_arsas$afd_zweit17/distr_btw17_arsas$gültige_zweit17)*10

# postcom 
postcom <- distr_btw17_arsas$postcom

# N's and id for states 
idl <- distr_btw17_arsas$land_id
N <- length(afd)
NC <- max(idl)


### 
# neighbourhood structure for electoral districts to be fed into the model 

# list of physical neighbours - queens logic 
nb <- poly2nb(distr_btw17_arsas, row.names=distr_btw17_arsas$WKR_NR)

#A list of vectors into a single vector 
nb2 <- unlist(nb)

#Count of the number of neighbors per district 
num <- card(nb)

#Adding weights - 1, for example 
weight <- rep(1, times=length(nb2))

# style="B" to avoid row-standardization 
nbw_raw <- nb2listw(nb, zero.policy=TRUE,style="B")


####
# intro to Stan 
# http://mc-stan.org/
# Stan: named after Stanislaw Ulam, pioneer of Monte Carlo methods 
# or "Sampling Through Adaptive Neighbourhoods" 
# -> but still Stan, not STAN, just like Stata  

# quick start, examples: https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started

# Models from Gelman's book on Bayesian Analysis translated to Stan
# https://github.com/stan-dev/example-models/wiki/ARM-Models-Sorted-by-Type

# a linear model - usually saved in text file with extension .stan 
# which is recognized by RStudio 

# contents of linear model file lin.stan 

lin.stan <- "
data {
  int<lower=0> N;
  vector[N] att;
  vector[N] afd;
}
parameters{
  vector[2] beta;
  real<lower=0> sigma;
}
model {
  att ~ normal(beta[1] + beta[2] * afd, sigma);
}
"

# always pattern of data / parameters / model 
# data: declare all variables, make restrictions 
# parameters: unobserved parameters to be estimated, 
# often featuring transformations 
# model: the usual, priors are uniform if not declared otherwise 

# "stan_model" compiles code to very general programming language C++ 
#-> tons of flexibility 
# if it takes time, silence please: 
# Stan is "figuring out the gradient functions for the 
# Hamiltonian dynamics"
# (Kruschke 2015: 408)
# good for testing the model 

# a model needs to be compiled only once in a session 
# and will be reused if unchanged 
linstan <- stan_model('lin.stan')

# fitting a Stan model using lin.stan 
# (compiles as well if not already done so)
linfit <- stan(file = 'lin.stan', data = list(N,att,afd), 
               iter = 1000, warmup = 500, chains = 4)

# does 1000 iterations for each chain, discarding 500 as "burn-in"
# 4 chains because it's possible - usually, three is fine 

# print results 
linfit

# plot results 
plot(linfit)

#traceplot 
traceplot(linfit, pars=c("beta[2]"))
# mixing nicely, no autocorrelation 

# posterior density of parameter pairs  
pairs(linfit)


###
# brms package 
# translates models form the "arm" package to Stan - thanks! 
library(brms)
standata <- data.frame(att,afd,idl)

# example: zero-inflated poisson multilevel model 

# generating Stan code with "make_stancode" 
infl_count.stan <- make_stancode(att ~ afd + (1|idl),
                                 data = standata, 
                                 family = "zero_inflated_poisson")

# code is often a bit fancy, well-commented -> nice learning feature 
infl_count.stan

# direct estimation via brms 
infl_count.mod <- brm(att ~ afd + (1|idl),
                      data =  standata, 
                      family = "zero_inflated_poisson")

# seems to use 1000 iterations for burn-in/inference as default
# and 4 chains - why 3 when you can have 4, right? 

# nice amount of explanation
# zi parameter shows that zero inflation seems appropriate 
print(infl_count.mod)
# rather fancy default graphs - parameters as distributions! 
plot(infl_count.mod)

# web interface for further diagnostics 
# launch_shinystan(count.mod)


############
# exercise #
# estimate the model "lm (att ~ afd)"
# compare the results to the "lin.stan" model 
# translate "lm (att ~ afd)" using the brms package  
# compare the Stan code to "lin.stan"
# add the explanatory variable "postcom" 
# to "lin.stan" and its brms version, estimate the models,
# and plot the results 
###



#### 
# a spatial Bayesian model 
# implementation of ICAR in Stan 
# following the case study provided by Mitzi Morris under
# http://mc-stan.org/users/documentation/case-studies/icar_stan.html
# -> visit site 

# ICAR: spatial random effects  

# starting with an "empty" linear spatial model 
# using German electoral district neighbourhood structures 
# only priors for spatial effects, no data! 

# model: saved with .stan extension (simple_iar.stan )
# contents: 

simple_iar.stan <- "
data {
  int<lower=0> N;
  int<lower=0> N_edges;
  int<lower=1, upper=N> node1[N_edges];  
  // node1[i] adjacent to node2[i]
  int<lower=1, upper=N> node2[N_edges];  
  // and node1[i] < node2[i]
}
parameters {
  vector[N] phi; // spatial random effects 
}
model {
  target += -0.5 * dot_self(phi[node1] - phi[node2]); 
  // target += from C++, adds term to the log posterior, 
  // same as assigning distribution   
  // dot_self computes the dot product of a vector with itself
  // soft sum-to-zero constraint on phi
  sum(phi) ~ normal(0, 0.01 * N);  
  // equivalent to mean(phi) ~ normal(0,0.01) 
  // -> needed for identification
}
"

# Uses and transforms WinBUGS-style neighbourhood information: 
# "We have written a helper function mungeCARdata4stan 
# which can transform the fields data$adj and data$num into a list 
# structure with fields N, N_edges, node1, and node2 
# which correspond to the inputs required by the Stan model"

mungeCARdata4stan = function(adjBUGS,numBUGS) {
  N = length(numBUGS);
  nn = numBUGS;
  N_edges = length(adjBUGS) / 2;
  node1 = vector(mode="numeric", length=N_edges);
  node2 = vector(mode="numeric", length=N_edges);
  iAdj = 0;
  iEdge = 0;
  for (i in 1:N) {
    for (j in 1:nn[i]) {
      iAdj = iAdj + 1;
      if (i < adjBUGS[iAdj]) {
        iEdge = iEdge + 1;
        node1[iEdge] = i;
        node2[iEdge] = adjBUGS[iAdj];
      }
    }
  }
  return (list("N"=N,"N_edges"=N_edges,"node1"=node1,"node2"=node2));
}

# from neighbours list to Stan specification, or:
# "To efficiently store the neighborhood structure,
# we encode the spatial adjacency matrix as an array of edges 
# of an undirected graph instead of using a large square matrix."

nbs = mungeCARdata4stan(nb2, num)
node1 = nbs$node1
node2 = nbs$node2
N_edges = nbs$N_edges

# compile code to C++
iar_stan = stan_model("simple_iar.stan")

# compiling and estimating in one step 
fit_stan <- stan(file = 'simple_iar.stan', 
                 data = list(N,N_edges,node1,node2), 
                 iter = 1000, warmup = 500, chains = 3)

# you get 299 estimated parameters for the spatial priors per district
print(fit_stan)

stan_plot(fit_stan)

# diagnostics 
stan_diag(fit_stan)

# traceplot 
traceplot(fit_stan)



###
# adding count model of attacks 
# and afd as an explanatory variable 
# (so far only prior information and neighbourhood structure!) 

# file: w_count.stan

w_count_stan = stan_model('w_count.stan')

fit_stan <- stan(file = 'w_count.stan', 
                 data = list(N,N_edges,node1,node2,att,afd),
                 iter = 1000, warmup = 500, chains = 3,
                 control = list(adapt_delta=.9))

print(fit_stan)

# lambda: spatial dependency parameter 
# phi[158]: spatial random effect for "Sächsische Schweiz - Osterzgebirge"
stan_plot(fit_stan, pars=c("beta","lambda","phi[158]"), show_density=T, fill_color="pink")

stan_trace(fit_stan, pars=c("beta","lambda","phi[158]"))

stan_dens(fit_stan, pars=c("beta","lambda","phi[158]"))

# some things saved in Stan object 
fit_stan@model_pars

dim(fit_stan)

dimnames(fit_stan)

# get simulations, summarize spatial random effects 
postm <- get_posterior_mean(fit_stan, pars=c("phi"))

summary(postm)

# fourth column averages the three chains 
summary(postm[,4])

# plot of spatial random effects 
hist(postm[,4])


###
# next steps: integrate multilevel, zero-inflated poisson... 
###
