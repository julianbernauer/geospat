###############################################################
#CDSS Course "Spatial Analysis of Geographic Data", Spring 2018
#Julian Bernauer, Data and Methods Unit, MZES
#Session 4, March 20, 2018 
#Open-ended spatial analysis in R 

###
# Going Bayesian - Spatial Interpolation 

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
# sunshine hours - data 
load("btw_sunstat.Rdata")

#Sun hours in electoral districts data 
colors <- brewer.pal(9, "YlOrRd") 
brks<-classIntervals(btw_sunstat$Jahr, n=9, style="quantile")
brks<- brks$brks 
plot(btw_sunstat, col=colors[findInterval(btw_sunstat$Jahr, brks,all.inside=TRUE)])
# already warns about missing data 

# Using WinBUGS - missing data is more difficult to deal with in Stan 
# But it's possible too! 

# linear spatial multilevel model for exposition - working up 


###
# neighbourhood information 

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


###
# sunshine hours 
sun <- btw_sunstat$Jahr
summary(sun)
# 119 missing values 

# districts and N 
idl <- distr_btw17_arsas$land_id
N <- length(afd)
NC <- max(idl)



###
# Prediciton from prior in linear model 

model.winbugs = "model{

# outcome 
for(k in 1:N){
sun[k] ~ dnorm(mu[k],tau)
mu[k] <- b.0   
}

# diffuse priors 
b.0 ~ dnorm(0, .001)

tau <- pow(sigma, -2)
sigma ~ dunif(0, 100)

}"


# send to Winbugs model file 
writeLines(model.winbugs,con="model.winbugs.txt")

# the data as a list 
winbugs.data <- list("N", "sun")

# initial values - not the priors 
winbugs.inits <- function() {list(b.0 = 1000)} 

# parameters to be kept 
# for prediction, just monitor mu!
winbugs.parameters <- c("b.0","sigma","mu")

# actual estimation, calls external Windows application which tends to become unresponsive if you touch anything - just go and get a coffee...  
# always use debug=TRUE for even more fun 
winbugs.linear <- bugs(data=c(winbugs.data), winbugs.inits, winbugs.parameters, "model.winbugs.txt", 
                        n.chains=3, bugs.directory = "c:/program files/winbugs14/", 
                        n.iter=10000, n.burnin=5000, n.thin=50, debug=TRUE)

# you can inspect the chains for convergence, normality, non-autocorrelation 

# get the results into R 
attach.bugs(winbugs.linear)

# a quick look at the results 
winbugs.linear 
plot(winbugs.linear)

ls(winbugs.linear)
winbugs.linear$mean$mu

#get predicted sunshine hours 
pred_linear <- winbugs.linear$mean$mu

linpred_data <- data.frame(sun,pred_linear)
linpred_data$linpred <- linpred_data$sun
linpred_data$linpred[is.na(linpred_data$linpred)] <- linpred_data$pred_linear[is.na(linpred_data$linpred)]
lpred <- linpred_data$linpred

colors <- brewer.pal(9, "YlOrRd") 
brks<-classIntervals(lpred, n=9, style="quantile")
brks<- brks$brks 
plot(btw_sunstat, col=colors[findInterval(lpred, brks,all.inside=TRUE)])


### 
# Borrowing strength in a multilevel model 

model.winbugs = "model{

# outcome 
for(k in 1:N){
sun[k] ~ dnorm(mu[k],tau)
mu[k] <- b.0[idl[k]] 
}

# no prior needed here -> multilevel, b.0 varies across countries 

tau <- pow(sigma, -2)
sigma ~ dunif(0, 100)

# random intercepts for states 
for(j in 1:NC){
b.0[j] ~ dnorm(mu.0,tau.0)
}

mu.0 ~ dunif(0,3000)
tau.0 <- pow(sigma.0, -2)
sigma.0 ~ dunif(0, 100)

}"


writeLines(model.winbugs,con="model.winbugs.txt")

winbugs.data <- list("N", "NC", "idl", "sun")

winbugs.inits <- function() {list(mu.0 = 1500)} 

winbugs.parameters <- c("mu.0","b.0","sigma","sigma.0","mu")

winbugs.multilevel <- bugs(data=c(winbugs.data), winbugs.inits, winbugs.parameters, "model.winbugs.txt", 
                       n.chains=3, bugs.directory = "c:/program files/winbugs14/", 
                       n.iter=10000, n.burnin=5000, n.thin=50, debug=TRUE)

attach.bugs(winbugs.multilevel)
winbugs.multilevel 
plot(winbugs.multilevel)

#get predicted sunshine hours 
pred_multi <- winbugs.multilevel$mean$mu

multipred_data <- data.frame(sun,pred_multi)
multipred_data$multipred <- multipred_data$sun
multipred_data$multipred[is.na(multipred_data$multipred)] <- multipred_data$pred_multi[is.na(multipred_data$multipred)]
mpred <- multipred_data$multipred

colors <- brewer.pal(9, "YlOrRd") 
brks<-classIntervals(mpred, n=9, style="quantile")
brks<- brks$brks 
plot(btw_sunstat, col=colors[findInterval(mpred, brks,all.inside=TRUE)])


###
# Adding a spatial component 

model.winbugs = "model{

# outcome 
for(k in 1:N){
sun[k] ~ dnorm(mu[k],tau)
mu[k] <- b.0[idl[k]] + v[k]  
}

tau <- pow(sigma, -2)
sigma ~ dunif(0, 100)

# spatially autocorrelated errors -> car.normal implements conditional autoregressive model (CAR) 
v[1:N] ~ car.normal(nb2[], weight[], num[], vtau)
vtau <- pow(vsigma, -2)
vsigma ~ dunif(0, 100)

# random intercepts for states 
for(j in 1:NC){
b.0[j] ~ dnorm(mu.0,tau.0)
}

mu.0 ~ dunif(0,3000)
tau.0 <- pow(sigma.0, -2)
sigma.0 ~ dunif(0, 100)

}"


writeLines(model.winbugs,con="model.winbugs.txt")

winbugs.data <- list("N", "NC", "idl", "nb2", "num", "weight", "sun")

winbugs.inits <- function() {list(mu.0 = 1500)} 

winbugs.parameters <- c("b.0","sigma","vsigma","sigma.0","v","mu.0","mu")

winbugs.spatial <- bugs(data=c(winbugs.data), winbugs.inits, winbugs.parameters, "model.winbugs.txt", 
                       n.chains=3, bugs.directory = "c:/program files/winbugs14/", 
                       n.iter=10000, n.burnin=5000, n.thin=50, debug=TRUE)


attach.bugs(winbugs.spatial)
winbugs.spatial 
plot(winbugs.spatial)

#get predicted sunshine hours 
pred_spat <- winbugs.spatial$mean$mu

spatpred_data <- data.frame(sun,pred_spat)
spatpred_data$spatpred <- spatpred_data$sun
spatpred_data$spatpred[is.na(spatpred_data$spatpred)] <- spatpred_data$pred_spat[is.na(spatpred_data$spatpred)]
spred <- spatpred_data$spatpred

colors <- brewer.pal(9, "YlOrRd") 
brks<-classIntervals(spred, n=9, style="quantile")
brks<- brks$brks 
plot(btw_sunstat, col=colors[findInterval(spred, brks,all.inside=TRUE)])


### you could also add covariates... see Selb and Munzert (2011)


### that's all folks! 
