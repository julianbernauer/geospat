#Guess what you see 

library(rgdal)

####
#Example with rgdal package: Location of "???" in Mannheim 
#Files retrieved 6 Oktober 2017 from https://mannheim.opendatasoft.com/explore/dataset/standorte_hts_gesamt/export/

#load the shapefile from the working directory into a R object with the rgdal package (alternative: maptools)
x1 <- readOGR(dsn = ".")

#contents of shapefile 
summary(x1)

#plot shapefile
plot(x1) 
