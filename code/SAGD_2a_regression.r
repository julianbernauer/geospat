###############################################################
#CDSS Course "Spatial Analysis of Geographic Data", Spring 2018
#Julian Bernauer, Data and Methods Unit, MZES
#Session 2, February 27, 2018 

#################
#Part II: MODELS#
#################

# CONTENTS:
# Spatial regression
# Spatial correlation  
# Illustrations using the running example: AfD vote and attacks on refugees 
# Starting your own project 

# set your working directory 
setwd("...")
#setwd("U:/Service MZES/Geodata/workshop cdss/code")


###
# Packages 

# maptools to read in the polygon data 
library(maptools) 

# spdep and sp for spatial analysis 
library(spdep)
library(sp) 

# rgdal helps dealing with more recent formats such as geoJSON 
library(rgdal)

# to make nice choropleth maps 
library(RColorBrewer) 
library(classInt)


###
# load object with 2017 BTW results and data on attacks on 
# refugees (see SAGD_1a_geodata.r for details)
load("distrbtw17.Rdata")
ls()


###
# repetition: generating W, the connectivity matrix 

# list of physical neighbours - queens logic 
nb <- poly2nb(distr_btw17, row.names=distr_btw17$WKR_NR)
nb

# generating W from a list of neighbours - two variants for comparison 
# default: row-standardized, usually wrong according to Neumayer and Plümper 2016
nbw <- nb2listw(nb, zero.policy=TRUE, style="W")
nbw 

# style="B" to avoid row-standardization 
nbw_raw <- nb2listw(nb, zero.policy=TRUE,style="B")
nbw_raw
# we'll work with this variant 


###
# spatial regression 

###
#a few variable for further analysis  
# afd in per cent 
afd <- (distr_btw17$afd_zweit17/distr_btw17$gültige_zweit17)*100
# attacks using the log - involves assumption! 
plot(density(distr_btw17$attacks))
att <- log(distr_btw17$attacks+1)
plot(density(att))
# slightly bimodal 
# count data, theoretical difference between 0 and 1 or more, ponder other models... 
postcom <- distr_btw17$postcom
npd <- distr_btw17$npd_zweit17/distr_btw17$gültige_zweit17
npd[is.na(npd)] <- 0


###
# OLS as baseline 
afdatt_ols1 <- lm(att ~ afd) 
summary(afdatt_ols1)


###
# testing SEM vs. SLM using a Lagrange Multiplier Test 
# starting from OLS model and checking residuals for spatial patterns  
res <- lm.LMtests(afdatt_ols1, nbw_raw, test=c("LMerr", "LMlag"))
summary(res)
# both are at the margins of conventional levels of statistical significance 
# tendency of favouring SEM 
# HARD TO TELL  


###
# actual spatial regression model 
# from the spdep package: SLM - spatial lag model
# maximum likelihood - should (!) consider endogeneity 
# tests contagion: attacks from the surrounding areas affect attack levels 
afdatt_spatlag1 <- lagsarlm(att ~ afd, ,nbw_raw ,type="lag")
summary(afdatt_spatlag1)
# afd vote is positively associated (causality unclear, time gap...) with attacks 
# rho: attacks tend to be slightly more likely in areas where attacks occur 


###
# from the spdep package: SEM - spatial error model 
# tests unobserved spatial confounders 
afdatt_spater1 <- errorsarlm(att ~ afd, ,nbw_raw)
summary(afdatt_spater1)
#lambda: unobserved spatial factors slightly affect the attack levels  

# -> models reflect the Lagrange Multiplier Tests 
# -> which way would you go? 


###
# adding controls - usually from the start

###
# postcom in OLS
afdatt_ols2 <- lm(att ~ afd + postcom) 
summary(afdatt_ols2)

# plot of residuals actually asks us to specify W! 
plot(afdatt_ols2)
# not really helpful to learn about spatial structure -> see Lagrange Multiplier Test above


###
# quick exercise: test the interaction between afd and postcom 
# what idea are we testing here? 
###


###
# SLM 
# controlling for postcom 
afdatt_spatlag2 <- lagsarlm(att ~ afd + postcom, ,nbw_raw ,type="lag")
summary(afdatt_spatlag2)
# what happens here? 


###
# SEM 
afdatt_spater2 <- errorsarlm(att ~ afd + postcom, ,nbw_raw)
summary(afdatt_spater2)


###
# effects of raw vs. row-standardized W 
# example of SEM 

# with raw W
summary(afdatt_spater1)
# with row-standardized W 
afdatt_spater3 <- errorsarlm(att ~ afd, ,nbw)
summary(afdatt_spater3)
# affects lambda - props to Neumaye and Plümper 2016! 


### 
# SAR - SEM and SLM in one 
# in spdep: "SAC"
afdatt_sar <- sacsarlm(att ~ afd, ,nbw_raw)
summary(afdatt_sar)
# a bit too much? 


###
# Durbin model: adding spatially lagged independent variables to the lag model 
afdatt_durb <- lagsarlm(att ~ afd , ,nbw ,type="mixed")
summary(afdatt_durb)
# spatial lag of AfD vote has no systematic positive coefficent 
# but affects rho, which now looks bigger! 
# afd vote itself is more strongly clustered geographically than attacks (see below)



#####################

###
# spatial correlation 


###
# Global correlation: Moran's I
# attacks and also afd vote 
moran.test(att,nbw_raw)
moran.test(afd,nbw_raw)


###
# afd correlogram - correlation for different spatial lags 
sp.cor <- sp.correlogram(nb, afd, order=10, method="I", randomisation=FALSE)
sp.cor
plot(sp.cor,main="Spatial structure AfD share")


# attacks 
sp.cor2 <- sp.correlogram(nb, att, order=10, method="I", randomisation=FALSE)
plot(sp.cor2,main="Spatial structure cumulated attacks")
# a bit different - gap after two lags 


###
# alternative measure: C 
sp.cor2 <- sp.correlogram(nb, att, order=10, method="C", randomisation=FALSE)
plot(sp.cor2,main="Spatial structure cumulated attacks")


###
# local spatial statistic: variants of Moran's I and Geary's C 
# here: local Moran's I 
res_locI <- localmoran(att, nbw_raw)
# gives us matrix with statistic on local dependency and a few auxilliary values  
colnames(res_locI)
# local coefficients 
res_locI[,1]
hist(res_locI[,1])
# p-values 
hist(res_locI[,5])
#local Moran's I for the electoral districts 
printCoefmat(data.frame(res_locI[distr_btw17$WKR_NR,1], row.names=distr_btw17$WKR_NAME[distr_btw17$WKR_NR]), check.names=FALSE)


###
# quick exercise 
# put these local coefficients into a choropleth map 
###


###
# Moran plot - similarity of an observation with its surroundings 
moran.plot(att, nbw_raw, labels=distr_btw17$WKR_NAME, xlim=c(-2,5))
# lower right-hand side: more attacks than expected in surrounding area 
moran.plot(afd, nbw_raw, labels=distr_btw17$WKR_NAME)



###################################################

###
###exercise - running example attack data hackathon 

# first, we split up the dependent variable to just use incidences of arson and assault 

attacks <- readOGR(dsn = "vorfaelle.geojson")
plot(attacks)
summary(attacks$title)

# demos only - 345 cases 
attacks_demo <- attacks[attacks$title=="Kundgebung/Demo",]
summary(attacks_demo$title)
plot(attacks_demo)

#other only - 2240 cases!
attacks_other <- attacks[attacks$title=="Sonstige Angriffe auf UnterkÃ¼nfte",]
summary(attacks_other$title)
plot(attacks_other)

# arson or assault - 221 + 505 = 726 cases 
attacks_help <- attacks[attacks$title!="Sonstige Angriffe auf UnterkÃ¼nfte",]
attacks_arsas <- attacks_help[attacks_help$title!="Kundgebung/Demo",]
summary(attacks_arsas$title)
plot(attacks_arsas)

# count of arson/assaults added to map data 
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

# distribution 
plot(density(distr_btw17$attacks_arsas))

# choropleth  
colors <- brewer.pal(9, "YlOrRd") 
brks<-classIntervals(distr_btw17$attacks_arsas, n=9, style="quantile")
brks<- brks$brks 
plot(distr_btw17, col=colors[findInterval(distr_btw17$attacks_arsas, brks,all.inside=TRUE)], axes=F)

# in map of afd vote 
colors <- brewer.pal(9, "Blues") 
brks<-classIntervals(distr_btw17$afd_zweit17/distr_btw17$gültige_zweit17, n=9, style="quantile")
brks<- brks$brks 
plot(distr_btw17, col=colors[findInterval(distr_btw17$afd_zweit17/distr_btw17$gültige_zweit17, brks,all.inside=TRUE)], axes=F)
plot(attacks_arsas, add=TRUE)


### 
# now, get some structural data for the 2017 BTW election - like unemployment or population density... 
# option 1: from the Bundeswahlleiter with distric ids for merging, https://bundeswahlleiter.de/bundestagswahlen/2017/strukturdaten.html
# option 2: from a open geodata repository with geographic information, https://hub.arcgis.com/datasets/esri-de::wahlkreise-2017-mit-strukturdaten
# ...and play around with a few model specifications such as: 
# lags and/or errors, control for other parties (npd...), include structural variables (unemployment...), interactions... 
###



##################################

### 
# last part: starting your project 
# first idea for RQ 
# look for data, for instance starting here: http://gisgeography.com/best-free-gis-data-sources-raster-vector/
# think about spatial mechanisms -> W, errors, lags... 
# maybe do some first exploration 
# round of reporting and discussing 
###
