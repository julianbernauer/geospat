###############################################################
#CDSS Course "Spatial Analysis of Geographic Data", Spring 2018
#Julian Bernauer, Data and Methods Unit, MZES
#Session 4, March 20, 2018 
#Open-ended spatial analysis in R 

###
# Going Bayesian - time and space 

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

###
# tried to use the daily attacks data, but too rare for modelling 
# daily time series 
# after a look at the daily data, 
# we'll walk through the model assuming 
# that we have time series data 
# with less 0s / NAs 

# data requirements: 
# neighbourhood information on districts as before 
# but dependent variable = districts-day (e.g.) dyads! 
# new ids for district (WKR_NR) and time (date) needed 
# as well as N of districts -> additional context level 


###
# a look at the daily data 
attacks <- readOGR(dsn = "vorfaelle.geojson")
# using all attacks to maximize data input 
date <- attacks$date
date 
# could be worked on with lubridate package 

# Connecting districts and attacks via coordinates 
# - with "over", but not count per district as before 
load("distrbtw17.Rdata")
atts <- over(attacks,distr_btw17[,"WKR_NR"])  

# returns WKR_NR per attack 
atts_wkrnr <- as.numeric(unlist(atts))

# combine with date data, generate data frame
daydata <- data.frame(atts_wkrnr,date)

# leaves us with relatively little daily data over a long period 
length(daydata$date)
# about 3300 of about 210'000 district-day dyads, or 1.6 % 

# we could now add zeros for all other district-days, 
# generate district and time ids, 
# and try to model the AR(1) process,
# but this is not promising... 


###
# options
# 1) abandon time series approach 
# 2) use all attacks, at the price of heterogeneity 
# 3) switch to monthly analysis 
# -> probably still sparse data, but might be worth a try  


# Examples: Sächsische Schweiz - Osterzgebirge 
daydata[daydata$atts_wkrnr==158,]
# Mannheim 
daydata[daydata$atts_wkrnr==275,]
# some district 
daydata[daydata$atts_wkrnr==99,]


### discuss 
# maybe go for monthly data? 
# how useful is daily data? 


###
# for now: walking through the model parts 
###


###
# Parts for the goal: 
# Spatial logistic model with AR(1) component within districts 
# (starting from Lunn et al. 2013: 257ff)
# Multilevel components for time within in districts in addition 
# to districts within states 
# using idd: district id; ND: number of districts 

##
# needs combination of this: AR(1) time series for a single unit 
# NT: points in time observed 

model.winbugs.logit.AR = "model{

# outcome - bernoulli distribution with logistic link function 
for(t in 1:NT){
att_daily[t]  ~ dbern(p.att[t])
p.att[t] <- 1/(1+exp(-z[t]))
}

# note the index starting from 2!
for(t in 2:NT){
z[t]     <- b.0 + theta*z[t-1]
eps[t]   <- attd[t] - z[t]
}

# for t=1! 
z[1]     <- attd[1] - eps[1] 

eps[1] ~ dnorm(0, .001)
b.0 ~ dnorm(0, .001)
theta ~ dnorm(0, .001) 

}"


##
# and this: logistic spatial multilevel model across districts... 
# NC: state id
# N: districts 

model.winbugs.logit = "model{

for(i in 1:N){
att[i]  ~ dbern(p.att[i])
p.att[i] <- 1/(1+exp(-z[i]))
z[i] <- b.0[[idl[i]]] + b.afd*afd[i] + v.att[i] 
}

b.afd ~ dnorm(0, .001) 

v.att[1:N] ~ car.normal(nb2[], weight[], num[], vtau.att)
vtau.att <- pow(vsigma.att, -2)
vsigma.att ~ dunif(0, 100)

for(j in 1:NC){
b.0[j] ~ dnorm(mu.0,tau.0)
}

mu.0 ~ dunif(0,100)
tau.0 <- pow(sigma.0, -2)
sigma.0 ~ dunif(0, 100)

}"


###
# will be shown in the next session... 
###
