###############################################################
#CDSS Course "Spatial Analysis of Geographic Data", Spring 2018
#Julian Bernauer, Data and Methods Unit, MZES
#Session 3, March 6, 2018 
#Implementing spatial analysis in R 

###
# contents
# getting fancy with geodata in R: more packages and maps 
# plotting in R - plotly, ggplot2 
# exercise: replicating Economist plot 
# varied outcomes: categorical and count data -> Bayesian framework, second file 

####
library(maptools) 
library(spdep)
library(sp) 
library(rgdal)
library(RColorBrewer) 
library(classInt)

#setwd("...")
setwd("U:/Service MZES/Geodata/workshop cdss/code")

###
# a different example data set starting from coordinates (e.g. given by GPS recording): 
# climate data Germany from Deutscher Wetter Dienst 
# yearly/monthly sunshine hours between 1981 and 2010 

# ftp://ftp-cdc.dwd.de/pub/CDC/observations_germany/climate/multi_annual/mean_81-10/
# stations with coordinates   
stations_deu <- read.csv("Sonnenscheindauer_1981-2010_Stationsliste_festerStandort.txt",header=TRUE,sep=";",na.strings=c("\t"))
View(stations_deu)

# merge sunshine hours 
sun_deu <- read.csv("Sonnenscheindauer_1981-2010_festerStandort.txt",header=TRUE,sep=";",na.strings=c("\t"))
sunstat <- merge(stations_deu,sun_deu,by="Stations_id")
sunstat$X.x <- NULL

# use coordinates to create SpatialPointsDataFrame
# "structure" command to create data frame and rename variables  
coord <- structure(list(longitude = sunstat$geogr..Laenge, latitude = sunstat$geogr..Breite), .Names = c("longitude","latitude"), class = "data.frame", row.names = c(NA, -260L))
# or quick data frame  
#coord2 <- data.frame(sunstat$geogr..Laenge, sunstat$geogr..Breite)
statdeu_spdf <- SpatialPointsDataFrame(coords = coord, data = sunstat, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
plot(statdeu_spdf)

# from points to polygons: connecting sunshine hours and electoral districts 
# or: not everything has to make sense, right? 
load("distrbtw17.Rdata")

# creating variable assigning stations to districts with "over" command 
stations_electoral <- over(statdeu_spdf,distr_btw17[,"WKR_NR"])  
stations_electoral[1:10,]
sunstat$WKR_NR <- as.numeric(unlist(stations_electoral))

# stations per district 
stat_count <- as.numeric(table(stations_electoral))
stat_count

# collapse by district 
sunstat_sub <- subset(sunstat, select=c(Jahr,Feb.,WKR_NR))
sunstat_d <- aggregate(sunstat_sub, by=list(sunstat_sub$WKR_NR), FUN=mean, na.rm=TRUE)

# merge with district map data 
btw_sunstat <- merge(distr_btw17,sunstat_d,by="WKR_NR")


#Sonnenstunden in den Wahlkreisen 
colors <- brewer.pal(9, "YlOrRd") 
brks<-classIntervals(btw_sunstat$Jahr, n=9, style="quantile")
brks<- brks$brks 
plot(btw_sunstat, col=colors[findInterval(btw_sunstat$Jahr, brks,all.inside=TRUE)])
#Wetterstationen dazu 
plot(statdeu_spdf, add=T)


# by the way: labels:
saar <- distr_btw17[distr_btw17$LAND_NAME=="Saarland",]
plot(saar)
# finds the middle of the polygon 
invisible(text(getSpPPolygonsLabptSlots(saar), labels=as.character(saar$WKR_NR), cex=1.2))


###
#Dunkeldeutschland? 
attacks <- readOGR(dsn = "vorfaelle.geojson")
plot(btw_sunstat, col=colors[findInterval(btw_sunstat$Jahr, brks,all.inside=TRUE)])
plot(attacks, add=T)


###
# -> such data can also be used for interpolating missing values 
# -> idea of small area estimation... 
### 


###
# alternative: grid data germany with sunshine 
# ftp://ftp-cdc.dwd.de/pub/CDC/grids_germany/multi_annual/

# usually, packages should be loaded at the start, but to indicate they are needed here... 
library(SDMTools)
deu_sun_map <- read.asc.gz("grids_germany_multi_annual_sunshine_duration_1981-2010_17.asc.gz")
deu_sun_map

library(raster)
deu_sun_map_rast <- raster(deu_sun_map)
plot(deu_sun_map_rast, col=heat.colors(100))
plot(deu_sun_map_rast, col=terrain.colors(10))

 
  
###
# visualizing geodata: from random mapping in R to a more systematic approach


### 
# maps package has some, well, maps 
# a world map 
library(maps)
data(worldMapEnv)
# map command 
map('world', fill = TRUE, col = 1:10, wrap=c(-180,180))

# US states 
data(stateMapEnv)
map('state', fill = TRUE,col = palette())


###
# dismo package 
# gmaps - pictures of maps, not actual map data  
# also has an API with access to more data 
# -> we'll have a look at open street map (osm) fpr such functions 
library(dismo)
germ <- gmap("Deutschland")
plot(germ)

# you can type a search term and it looks for a fitting map! 
# sometimes it finds a map it didn't find before 
ice <- gmap("Iceland", type="satellite")
plot(ice)

ma <- gmap("Mannheim A5,6")
plot (ma)

zug <- gmap("Zugspitze", type="satellite")
plot(zug)

# find atlantis 
atl <- gmap("Atlantis")
plot(atl)

# zoom out by setting exp to higher values 
# (below 1 for zooming in)
mar <- gmap("Marianegraben", exp=100)
plot(mar)

mar <- gmap("Marianegraben", exp=1000)
plot(mar)


###
# quick exercise 
# get a satellite picture of Berlin and zoom in to find the Berghain 


###
# package ggmap 
# map own data on gmaps material 
library(ggmap)

# find a location via single coordinate point and zoom 
germd <- get_googlemap(c(lon=10.4475, lat=51.163), zoom=6)
ggmap(germd)

# add for instance sun hours data via tis coordinates 
ggmap(germd) + geom_point(aes(x=longitude, y=latitude), data=coord, colour="red", size=2)


###
# open street map - osm 
# has a lot of crowdsourced data 
# application programming interface - API 
# package "osmar"
library(osmar)

# define a box: left, bottom, right, top   
# lat and lon from https://www.openstreetmap.org/export#map=16/49.4871/8.4693
box <- corner_bbox(8.450, 49.485, 8.465, 49.490)
box 

# alternative: center point, heigth and width in meters 
cbox <- center_bbox(8.458772, 49.487384, 100, 100)
cbox

# extract map data from osm for box using the api 
api <- osmsource_api(url = "http://api.openstreetmap.org/api/0.6/")
MA <- get_osm(cbox, source = api)
MA

# gives us nodes, ways, and relations 

# all
plot(MA)

# nodes 
plot_nodes(MA)

# ways 
plot_ways(MA)

# contents - plenty of information 
summary(MA$nodes)
summary(MA$ways)
summary(MA$relations)

# looking for a tree? 
trees <- find(MA, node(tags(v == "tree")))
# five should be around the corner 
trees
treedata <- subset(MA, node_ids = trees)

# plot the trees 
plot(MA)
plot_nodes(treedata, add=TRUE, col="red")

#convert to sp polygon data for further analysis - makes sense especially for buildings... 
ma_poly <- as_sp(MA, "polygons")
plot(ma_poly)

# see for more, e.g. using osm to program a navigator:
# https://journal.r-project.org/archive/2013-1/eugster-schlesinger.pdf


###
# quick exercise 
# get a piece of osm map data: Bundeskanzerlamt 
# subset it looking for trees close to it
# highlight the trees on the map and find the closest tree 



###
# a step back - more systematically 
# plotting in R 
# book by Wickham and Grolemund (2017): R for Data Science 
# -> sort of pimp my R, 2010s instead of 1990s 
# -> integrates well with RStudio (Wickham works there) 
# R in terms of tidyverse, ggplot2 and plotly...


###
# tidyverse: sort of a "tidy" environment within R  
# loads ggplot2, tibble, dplyr...
library(tidyverse)

# data as tibble: neat data frame which gets rid of common problems such as factors... 
wkr_nr <- distr_btw17$WKR_NR
wkr_name <- distr_btw17$WKR_NAME
land_id <- distr_btw17$land_id
land_name <- distr_btw17$LAND_NAME
postcom <- distr_btw17$postcom
afd <- distr_btw17$afd_zweit17/distr_btw17$gültige_zweit17
att <- distr_btw17$attacks

ataf_help  <- data.frame(wkr_nr,wkr_name,land_id,land_name,postcom,afd,att)
ataf_help
ataf <- as.tibble(ataf_help)
ataf
# immediately nicer presentation: cuts off earlier, useful comments 
ataf_help$att
ataf$att


###
# quickly a few other nice to haves with dplyr - transformation 
filter(ataf, land_name=="Sachsen")

arrange(ataf, -afd)

select(ataf, att,afd)

# tidyverse also has features for strings, time/dates, pipe operators, spreading/gathering, merging... see book 


###
# ggplot2: gives you more control over plots
# might involve more code for simpler and less code for complex plots 
ggplot(data = ataf) +
  geom_point(mapping = aes(x = afd, y = att))

ggplot(data = ataf) +
  geom_point(mapping = aes(x = afd, y = att, color = postcom))

# or alphs, or shape... 

# only something going on in sachsen-anhalt, brandenburg and mecklenburg-vorpommern? 
ggplot(data = ataf) +
  geom_point(mapping = aes(x = afd, y = att)) +
  facet_wrap(~ land_name, nrow = 6)

# combine as you wish... 
afd_ggplot <- ggplot(data = ataf, mapping = aes(x = afd, y = att)) +
  geom_point(mapping = aes(color = postcom)) + 
  geom_smooth()
afd_ggplot


###
# ggplot2 works with maps, too (p. 32)!
# e.g. on map from the maps package (map_data... but seems to fail), or other polygon data 

#polygons with readOGR works 
setwd("U:/Forschung MZES/Geography of Conflict/btw17_geometrie_wahlkreise_vg250_geo_shp")
districts <- readOGR("Geometrie_Wahlkreise_19DBT_VG250_geo.shp", layer="Geometrie_Wahlkreise_19DBT_VG250_geo") 
btw17 <- read.csv("btw17_kerg_edited.csv",header=TRUE,sep=";",na.strings=c("\t"))
btw17 <- btw17[btw17$land_id!=0,]
btw17 <- btw17[btw17$land_id!=99,]
distrbtw <- merge(districts,btw17,by="WKR_NR")
distrbtw$afd <- (distrbtw$afd_zweit17/distrbtw$gültige_zweit17)*100
setwd("U:/Service MZES/Geodata/workshop cdss/code")

# "fortify" - a lot of data 
library(plyr)
library(dplyr)
distrbtw@data$id = rownames(distrbtw@data)
distrbtw.points = fortify(distrbtw, region="id")
distrbtw.df = join(distrbtw.points, distrbtw@data, by="id")

# raw 
ggplot() +  geom_polygon(data=distrbtw.df, aes(x=long, y=lat, group=group))

# coord_quickmap to get projection right 
ggplot(distrbtw.df) + 
  aes(long,lat,group=group,fill=afd) + 
  geom_polygon() +
  geom_path(color="white") +
  coord_quickmap() 


# another example with choropleth map: http://ggplot2.tidyverse.org/reference/map_data.html


###
# plotly: intractive maps for the web 
library(plotly)

# example from  https://plot.ly/r/#maps -> https://plot.ly/r/choropleth-maps/

df <- read.csv('https://raw.githubusercontent.com/plotly/datasets/master/2014_world_gdp_with_codes.csv')

# light grey boundaries
l <- list(color = toRGB("grey"), width = 0.5)

# specify map projection/options
g <- list(
  showframe = FALSE,
  showcoastlines = FALSE,
  projection = list(type = 'Mercator')
)

p <- plot_geo(df) %>%
  add_trace(
    z = ~GDP..BILLIONS., color = ~GDP..BILLIONS., colors = 'Blues',
    text = ~COUNTRY, locations = ~CODE, marker = list(line = l)
  ) %>%
  colorbar(title = 'GDP Billions US$', tickprefix = '$') %>%
  layout(
    title = '2014 Global GDP<br>Source:<a href="https://www.cia.gov/library/publications/the-world-factbook/fields/2195.html">CIA World Factbook</a>',
    geo = g
  )

p

# Create a shareable link to your chart
# Set up API credentials: https://plot.ly/r/getting-started
# chart_link = plotly_POST(p, filename="choropleth/world")
# chart_link

# see https://cran.r-project.org/web/packages/plotly/plotly.pdf 


###
# ggplotly: from ggplot2 to plotly
ggplotly(afd_ggplot)

ggplotly(
  ggplot(distrbtw.df) + 
  aes(long,lat,group=group,fill=afd) + 
  geom_polygon() +
  geom_path(color="white") +
  coord_equal() 
)


###
# exercise Economist plot replication 
# teams of two 
# -> look in The Economists for maps,  
# their sources and try to replicate as close as possible 
# example of such an exercise using ggplot2: http://tutorials.iq.harvard.edu/R/Rgraphics/Rgraphics.html
# also see for R Economist theme: http://latticeextra.r-forge.r-project.org/man/theEconomist.theme.html

