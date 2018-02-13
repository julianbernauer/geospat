#CDSS Course "Spatial Analysis of Geographic Data", Spring 2018
#Julian Bernauer, Data and Methods Unit, MZES
#Session 1, February 20, 2018 

##################
#Part I: CONCEPTS#
##################

#CONTENTS:
#Running example: The AfD vote in the 2017 Bundestag election 
#Reading geodata into R: polygons and points 
#From points to area data 
#A few maps 
#Generating W
#Quick spatial correlation and regression
#Exercise: solutions will be posted as file "SAGB_1b"

setwd("U:/Forschung MZES/Geography of Conflict/btw17_geometrie_wahlkreise_vg250_geo_shp")

#A few packages which are useful for the analysis of geodata 
#Install them first if not done yet, then load 
library(spdep) 
library(sp) 
library(rgdal)
library(maptools) 
library(RColorBrewer) 
library(classInt)
#Details of the packages will be discussed 


################
#Introducing the running example: The AfD vote in the 2017 Bundestag election 
################
#Reproducing the map on the course plan / in the slides 


###
#Geography 
#Source for electoral districts: https://www.bundeswahlleiter.de/bundestagswahlen/2017/wahlkreiseinteilung/downloads.html
# -> Of main interest is the .shp "shapefile", contains the information about the polygons
# -> Comes in a folder with further files such as .dbf, which contains for instance area attribute information such as identification (here: "WKR_NR")
# -> Does not yet include further data such as election results 

#German electoral districts -> polygon data imported into an object using the readShapePoly command from the maptools package 
districts <- readShapePoly("Geometrie_Wahlkreise_19DBT_VG250_geo.shp") 

#atypical object: different levels of variables - click on view 

summary(districts)

#A first map with the simple plot command - automatically produces map 
plot(districts)


###
#AfD results of 2017 BTW: https://www.bundeswahlleiter.de/bundestagswahlen/2017/ergebnisse.html 
#Election results have been edited manually, smaller parties excluded, identification variable named "WKR_NR" -> file "btw17_kerg_edited.csv"
#Just attribute data for the electoral districts with area identification -> can be merged with the geographic information retrieved above 

#Read data into R 
btw17 <- read.csv("btw17_kerg_edited.csv",header=TRUE,sep=";",na.strings=c("\t"))

#removing some aggregate summary data 
btw17 <- btw17[btw17$land_id!=0,]
btw17 <- btw17[btw17$land_id!=99,]

#merge election results with geodata based on the common identification variable "WKR_NR"
distr_btw17 <- merge(districts,btw17,by="WKR_NR")
#a quick check -> plenty of election results now available, including the Garten and the Hip Hop parties! 
summary(distr_btw17)


#Congratulations! That's been the hardest part... 


###
#Chloropleth map: color shades for intervals 
 
#defining colors using the RColorBrewer package 
colors <- brewer.pal(9, "Blues") 

#creating 9 shades using the classInt package 
brks<-classIntervals(distr_btw17$afd_zweit17/distr_btw17$gültige_zweit17, n=9, style="quantile")
brks<- brks$brks 

#AfD results plot with shades of blue
plot(distr_btw17, col=colors[findInterval(distr_btw17$afd_zweit17/distr_btw17$gültige_zweit17, brks,all.inside=TRUE)], axes=F)


###
#attacks on refugees from https://www.mut-gegen-rechte-gewalt.de/chronik-karte - download in February 2017 
#attacks between 2015 and early 2017 recored (more recent data available, but contains only a few more cases) 
#geoJSON point data 
#geoJSON is a relatively new format, see https://macwright.org/2015/03/23/geojson-second-bite  
#as shown, geoJSON files transform into an interactive map if loaded on GitHub 
#read using the readOGR command from the rgdal package 
attacks <- readOGR(dsn = "vorfaelle.geojson")
summary(attacks)
plot(attacks) 


###
#reproducing the map on the course program 
#combine plots to show coincidence of attacks and AfD vote share 
#beware of different points in time/periods covered: attacks before vote! 
#replot AfD vote
plot(distr_btw17, col=colors[findInterval(distr_btw17$afd_zweit17/distr_btw17$gültige_zweit17, brks,all.inside=TRUE)], axes=F)
plot(attacks, add=TRUE)


###
#integrating data polygon and point data: election results and attacks  
#SpatialPolygonsDataFrame and SpatialPointsDataFrame
#coordinates of attacks to count of attacks per district 

#adapting CRS - Coordinate Reference System - don't ask 
attacks@proj4string
distr_btw17@proj4string
CRS.new <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
proj4string(distr_btw17) <- CRS.new 
distr_btw17@proj4string

#Creating Variable assigning attacks to districts with "over" command 
afd_att <- over(attacks,distr_btw17[,"WKR_NR"])  
afd_att[1:10,]

#Attacks per district 
attack_count <- as.numeric(table(afd_att))
attack_count
#158 Sächsische Schweiz: 112 attacks

#Not all districts covered, 17 without attacks 
length(attack_count)
#get district names from variable assigning attacks to districts 
attack_wkr <- as.numeric(names(table(afd_att)))
#combine with the attack count variable 
att_data <- cbind(attack_wkr,attack_count)
#add districts without attacks -> manual data 
add_data <- read.csv("adddata.csv",header=T,sep=";")
att_data <- rbind(att_data,add_data)
#now 299 districts 
length(att_data$attack_wkr)
#sort by district number 
att_data <- att_data[order(att_data$attack_wkr),]

#finally, add count of attacks to the polygon data with afd shares 
distr_btw17$attacks <- att_data$attack_count

#plot attack levels as chloropleth map, this time in red 
colors <- brewer.pal(9, "YlOrRd") 
brks<-classIntervals(distr_btw17$attacks, n=9, style="quantile")
brks<- brks$brks 
plot(distr_btw17, col=colors[findInterval(distr_btw17$attacks, brks,all.inside=TRUE)], axes=F)

#because it's possible: add attack locations
plot(attacks, add=TRUE)


#a few variable for further analysis  
afd <- distr_btw17$afd_zweit17/distr_btw17$gültige_zweit17
att <- distr_btw17$attacks
wkr <- att_data$attack_wkr
idl <- btw17$land_id
postcom <- btw17$postcom
npd <- distr_btw17$npd_zweit17/distr_btw17$gültige_zweit17
npd[is.na(npd)] <- 0


#Have a break?! 


#############
#Generating W 
#############

#W - the connectivity matrix 

#Generating W for the German districts - step-by-step 

#List of neighbors - queens logic (alternative: distance bands)
nb <- poly2nb(distr_btw17, row.names=distr_btw17$WKR_NR)
nb

### *** W! *** ###
#Generating w from a list of neighbors
nbw <- nb2listw(nb, zero.policy=TRUE)
nbw

#Count of the number of neighbors per district 
num <- card(nb)
summary(num)
num

#A list of vectors into a single vector 
nb2 <- unlist(nb)
#1638 pairs of districts with common borders 

#Adding weights - 1, for example -> Here's the reason why W is not a weighting matrix 
weight <- rep(1, times=length(nb2))
#Here, a matrix with different weights for pairs of districts could be used... 


###
#A quick spatial correlation 

#Moran's I
moran.test(afd,nbw)

#Correlogram - correlation for different spatial lags 
sp.cor <- sp.correlogram(nb, afd, order=10, method="I", randomisation=FALSE)
sp.cor

print(sp.cor)
plot(sp.cor,main="Spatial structure AfD share")


###
#A few quick spatial regressions -  mind the time gap! rather complete causal nonsense at the moment! 

#Autoregressive error model 
spaterr_spdep <- errorsarlm(att~afd,,nbw)
summary(spaterr_spdep)
#lambda: unobserved spatial factors may affect the attack levels 
#afd vote is positively associated (causality unclear, time gap...) with attacks 


#From the spdep package: autoregressive lag model, attacks from the surrounding areas affect attack levels 
spatlag1_spdep <- lagsarlm(att~afd,,nbw,type="lag")
summary(spatlag1_spdep)
#rho: attacks are more likely in areas where attacks occur 


#Durbin model: adding spatially lagged independent variables 
spatlag2_spdep <- lagsarlm(att~afd,,nbw,type="mixed")
summary(spatlag2_spdep)
#spatial lag of AfD vote has no systematic positive coefficent 


###
#Exercise working with the Bundestag election data: 
###
#0) Starting from an empty workspace, load all packages and the data 
#1) Choose a party other than the AfD and make a chloropleth map of its vote share with a suiting color 
#   Hints: 
#   - make life easier by defining a variable containing the party's votre share in the districts 
#   - to choose a color palette: display.brewer.all()
#2) Add the data on attacks on refugee shelters to the map 
#3) Copy job: integrate attack and voting data and create W
#4) Plot the spatial correlation for the party chosen
#5) If you got this far, play around with the spatial regression models: 
#   add your party's vote share as a control variable in spatial lag, error or Durbin models...  

##############the end#########################################
  
  

