###############################################################
#CDSS Course "Spatial Analysis of Geographic Data", Spring 2018
#Julian Bernauer, Data and Methods Unit, MZES
#Session 3, March 6, 2018 
#Implementing spatial analysis in R 

###
# Going Bayesian 
# today: WinBUGS to estimate a count model 
library(rjags)
library(R2WinBUGS)

#setwd("...")
setwd("U:/Service MZES/Geodata/workshop cdss/code")

load("distrbtw17.Rdata")

### 
# neighbourhood structure for electoral districts to be fed into the model 

# list of physical neighbours - queens logic 
nb <- poly2nb(distr_btw17, row.names=distr_btw17$WKR_NR)

#A list of vectors into a single vector 
nb2 <- unlist(nb)

#Count of the number of neighbors per district 
num <- card(nb)

#Adding weights - 1, for example -> Here's the reason why W is not a weighting matrix 
weight <- rep(1, times=length(nb2))

# style="B" to avoid row-standardization 
nbw_raw <- nb2listw(nb, zero.policy=TRUE,style="B")


###
# a few variables 
# attacks (arson/assaults) without log -> we'll go for a count model! 
# attacks <- readOGR(dsn = "vorfaelle.geojson")
attacks_demo <- attacks[attacks$title=="Kundgebung/Demo",]
attacks_other <- attacks[attacks$title=="Sonstige Angriffe auf Unterkünfte",]
# arson or assault - 221 + 505 = 726 cases 
attacks_help <- attacks[attacks$title!="Sonstige Angriffe auf Unterkünfte",]
attacks_arsas <- attacks_help[attacks_help$title!="Kundgebung/Demo",]
att_ars <- over(attacks_arsas,distr_btw17[,"WKR_NR"])  
att_ars[1:10,]
attack_count <- as.numeric(table(att_ars))
attack_count
length(attack_count)
attack_wkr <- as.numeric(names(table(att_ars)))
att_data <- cbind(attack_wkr,attack_count)
add_data <- read.csv("adddata_arsas.csv",header=T,sep=";")
att_data <- rbind(att_data,add_data)
length(att_data$attack_wkr)
att_data <- att_data[order(att_data$attack_wkr),]
distr_btw17$attacks_arsas <- att_data$attack_count
att <- distr_btw17$attacks_arsas 

# dummy attacks or not for logistic model 
attd <- att
attd[att_d>0] <- 1
table(attd)

#afd
afd <- (distr_btw17$afd_zweit17/distr_btw17$g�ltige_zweit17)*100

# postcom 
postcom <- distr_btw17$postcom

# N's and id for states 
idl <- distr_btw17$land_id
N <- length(afd)
NC <- max(idl)



###
# comparison SEM spdep 
afdatt_spater2 <- errorsarlm(att ~ afd + postcom, ,nbw_raw)
summary(afdatt_spater2)


###
# code for the Bayesian model 
# linear variant - error model 
# WinBUGS / GeoBUGS -> requires external installation 

model.winbugs = "model{

# outcome 
for(k in 1:N){
att[k] ~ dnorm(mu.att[k],tau.att)
mu.att[k] <- b.0[idl[k]] + b.afd*afd[k] + b.postcom*postcom[k] + v.att[k]  
}

# diffuse priors 
b.afd ~ dnorm(0, .001) 
b.postcom ~ dnorm(0, .001)

tau.att <- pow(sigma.att, -2)
sigma.att ~ dunif(0, 100)

# spatially autocorrelated errors -> car.normal implements conditional autoregressive model (CAR) 
v.att[1:N] ~ car.normal(nb2[], weight[], num[], vtau.att)
vtau.att <- pow(vsigma.att, -2)
vsigma.att ~ dunif(0, 100)

# random intercepts for states 
for(j in 1:NC){
b.0[j] ~ dnorm(mu.0,tau.0)
}

mu.0 ~ dunif(0,100)
tau.0 <- pow(sigma.0, -2)
sigma.0 ~ dunif(0, 100)

}"


# send to winbugs model file 
writeLines(model.winbugs,con="model.winbugs.txt")

# the data as a list 
winbugs.data <- list("N", "NC", "idl", "nb2", "num", "weight", "att", "afd", "postcom")

# initial values - not the priors 
winbugs.inits <- function() {list(mu.0 = .15, b.postcom=0, b.afd=0)} 

# parameters to be kept 
winbugs.parameters <- c("b.afd","b.postcom","tau.att", "sigma.att","vtau.att", "vsigma.att","mu.0","sigma.0","b.0","v.att")

# actual estimation, calls external Windows application which tends to become unresponsive if you touch anything - just go and get a coffee...  
# always use debug=TRUE or even more fun 
winbugs.attacks <- bugs(data=c(winbugs.data), winbugs.inits, winbugs.parameters, "model.winbugs.txt", n.chains=3, bugs.directory = "c:/program files/winbugs14/", n.iter=10000, n.burnin=5000, n.thin=50, debug=TRUE)

# you can inspect the chains for convergence, normality, non-autocorrelation 

# get the results into R 
attach.bugs(winbugs.attacks)

# a quick look at the results 
winbugs.attacks 



# logistic regression - using att_d

model.winbugs.logit = "model{

# outcome - bernoulli distribution with logistic link function 
for(k in 1:N){
attd[k] ~ dbern(p.att[k])
p.att[k] <- 1/(1+exp(-z[k]))
z[k] <- b.0[idl[k]] + b.afd*afd[k] + b.postcom*postcom[k] + v.att[k]  
}

b.afd ~ dnorm(0, .001) 
b.postcom ~ dnorm(0, .001)

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


writeLines(model.winbugs.logit,con="model.winbugs.txt")
winbugs.data <- list("N", "NC", "idl", "nb2", "num", "weight", "attd", "afd", "postcom")
winbugs.inits <- function() {list(mu.0 = .15, b.postcom=0, b.afd=0)} 
winbugs.parameters <- c("b.afd","b.postcom","vtau.att", "vsigma.att","mu.0","sigma.0","b.0")
winbugs.attacks <- bugs(data=c(winbugs.data), winbugs.inits, winbugs.parameters, "model.winbugs.txt", n.chains=3, bugs.directory = "c:/program files/winbugs14/", n.iter=10000, n.burnin=5000, n.thin=50, debug=TRUE)

# does not converge 

attach.bugs(winbugs.attacks)
winbugs.attacks 


###
# poisson 

# without spatial structure 

model.winbugs.poisson = "model{

# outcome - poisson distribution with exponential link function 
for(k in 1:N){
att[k] ~ dpois(mu.att[k])
log(mu.att[k]) <- b.0[idl[k]] + b.afd*afd[k] + b.postcom*postcom[k]   
}

b.afd ~ dnorm(0, .001) 
b.postcom ~ dnorm(0, .001)

for(j in 1:NC){
b.0[j] ~ dnorm(mu.0,tau.0)
}

mu.0 ~ dunif(0,100)
tau.0 <- pow(sigma.0, -2)
sigma.0 ~ dunif(0, 100)

}"


writeLines(model.winbugs.poisson,con="model.winbugs.txt")
winbugs.data <- list("N", "NC", "idl", "att", "afd", "postcom")
winbugs.inits <- function() {list(mu.0 = .15, b.postcom=0, b.afd=0)} 
winbugs.parameters <- c("b.afd","b.postcom","mu.0","sigma.0","b.0")
winbugs.attacks <- bugs(data=c(winbugs.data), winbugs.inits, winbugs.parameters, "model.winbugs.txt", n.chains=3, bugs.directory = "c:/program files/winbugs14/", n.iter=10000, n.burnin=5000, n.thin=50, debug=TRUE)
attach.bugs(winbugs.attacks)
winbugs.attacks 


###
# poisson with spatial structure - simplified (no state random effects, no postcom control)    

model.winbugs.poisson = "model{
    
    # outcome - poisson distribution with exponential link function 
    for(k in 1:N){
    att[k] ~ dpois(mu.att[k])
    log(mu.att[k]) <- b.0 + b.afd*afd[k] + v.att[k]  
    }
    
    # diffuse priors 
    b.0  ~ dnorm(0, .001) 
    b.afd ~ dnorm(0, .001) 

    # spatially autocorrelated errors -> car.normal implements conditional autoregressive model (CAR) 
    v.att[1:N] ~ car.normal(nb2[], weight[], num[], vtau.att)
    vtau.att <- pow(vsigma.att, -2)
    vsigma.att ~ dunif(0, 100)
    
  }"


writeLines(model.winbugs.poisson,con="model.winbugs.txt")
winbugs.data <- list("N", "nb2", "num", "weight", "att", "afd")
winbugs.inits <- function() {list(b.0 = .15, b.afd=0)} 
winbugs.parameters <- c("b.afd","vtau.att", "vsigma.att","b.0")
winbugs.attacks <- bugs(data=c(winbugs.data), winbugs.inits, winbugs.parameters, "model.winbugs.txt", n.chains=3, bugs.directory = "c:/program files/winbugs14/", n.iter=10000, n.burnin=5000, n.thin=50, debug=TRUE)
attach.bugs(winbugs.attacks)
winbugs.attacks 



# more on WinBUGS / GeoBUGS    
# https://www.mrc-bsu.cam.ac.uk/wp-content/uploads/geobugs12manual.pdf     


###
# Outlook 
# I will look into options for spatial lags and zero-inflation  
# and explore JAGS/Stan spatial modelling 
# Interpolation (e.g. WinBUGS), time and space 
# Putting it all together for the running example: 
# zero-inflated spatio-temporal multilevel model 
