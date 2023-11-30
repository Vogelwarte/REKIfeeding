###########################################################################################################################
###### PROJECTING ANALYSIS TO IDENTIFY ANTHROPOGENIC FEEDING SITES OF RED KITES ACROSS SWITZERLAND ################
###########################################################################################################################
# original idea by Nathalie Heiniger (MSc thesis 2020)
# uses output provided by "Feeding_analysis.r"
# created by steffen.oppel@vogelwarte.ch in July 2023
# building and forest layers provided by Jerome Guelat

library(tidyverse)
library(dplyr, warn.conflicts = FALSE)
options(dplyr.summarise.inform = FALSE)
library(sf)
library(lubridate)
library(leaflet)
library(adehabitatHR)
library(stars)
library(rworldmap)
library(trip)
library(sqldf)
data(countriesLow)
filter<-dplyr::filter
sf_use_s2(FALSE)


## set root folder for project
try(setwd("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/REKIfeeding"),silent=T)
try(setwd("C:/STEFFEN/OneDrive - Vogelwarte/REKI/Analysis/REKIfeeding"),silent=T)



###### LOADING DATA THAT PREDICTED FEEDING PROBABILITY ------------------------------------------------------------

CHgrid<-readRDS("output/REKI_feeding_grid.rds")
range(CHgrid$FEEDER_predicted)





########## READ IN RAW TRACKING DATA ########################
trackingdata<-readRDS(file = "data/REKI_trackingdata_raw.rds")



# DATA PREPARATION -------------------------------------------------------------

if("long_wgs" %in% names(trackingdata)){
  trackingdata<-trackingdata %>% rename(long=long_wgs, lat=lat_wgs)
}
trackingdata <- trackingdata %>%
  # filter(long>5.9) %>%
  # filter(lat>45.8) %>%
  # filter(long<10.6) %>%
  # filter(lat<48) %>%
  filter(!is.na(timestamp)) %>%
  filter(!is.na(long)) %>%
  filter(!is.na(year_id))

dim(trackingdata)


# converting to metric CRS prior to estimating distances
track_sf <- trackingdata %>% 
  st_as_sf(coords = c("long", "lat"))
st_crs(track_sf) <- 4326

track_sf <- track_sf %>%
  st_transform(crs = 3035) %>%
  dplyr::mutate(long_eea = sf::st_coordinates(.)[,1],
                lat_eea = sf::st_coordinates(.)[,2])

head(track_sf)


### REMOVE DUPLICATE TIME STAMPS AND THOSE WITH TOO FEW LOCATIONS ####
toofew<-track_sf %>% mutate(count=1) %>% group_by(year_id) %>%
  summarise(N=sum(count)) %>%
  filter(N<10)
track_sf<-track_sf %>% filter(!(year_id %in% toofew$year_id))
duplicates<-track_sf %>% mutate(count=1) %>% group_by(year_id, t_) %>%
  summarise(N=sum(count)) %>%
  filter(N>1)

track_sf[c(1526,1529),]
track_sf %>% filter(year_id=="2016_10") %>% filter(t_==ymd_hms("2016-06-20 21:00:00"))




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# QUANTIFY TIME IN AREA OF EACH SWISS GRID CELL
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


REKI_trips<-trip(track_sf, TORnames=c("t_","year_id"), correct_all=TRUE)			### switch to "DateTime" when using the raw locations
trg <- tripGrid(EV_all_trips, grid=grd,method="pixellate")						### this will provide the number of bird seconds spent in each grid cell
spplot(trg)			## plots the trips with a legend
proj4string(turbSPDF)<-proj4string(EVSP_all)


### CONVERT SPATIAL GRID TO SOMETHING WE CAN PLOT
spdf <- SpatialPixelsDataFrame(points=trg, data=trg@data)
HOTSPOTS<-data.frame(lat=spdf@coords[,2],long=spdf@coords[,1],time=spdf@data$z)
HOTSPOTS$time<-HOTSPOTS$time/(3600*24)								### this converts the seconds into bird days
summary(HOTSPOTS)

HOTSPOTS<-HOTSPOTS[HOTSPOTS$time>10,]

##### CONVERT TO SPATIAL POLYGONS FOR OVERLAY ###
ras<- raster(spdf)		# converts the SpatialPixelDataFrame into a raster
spoldf <- rasterToPolygons(ras, n=4) # converts the raster into quadratic polygons
proj4string(spoldf)<-proj4string(EVSP_all)


### PRODUCE NICE AND SHINY MAP WITH THE MOST IMPORTANT TEMPORARY CONGREGATION SITES ###

#MAP <- get_map(EGVUbox, source="google", zoom=4, color = "bw")		### retrieves a map from Google (requires live internet connection)
MAP <- get_map(location = c(lon = mean(long), lat = mean(lat)), source="google", zoom=4, color = "bw")		### retrieves a map from Google (requires live internet connection)

pdf("EGVU_MIGRATION_HOTSPOTS.pdf", width=12, height=11)
ggmap(MAP)+geom_tile(data=HOTSPOTS, aes(x=long,y=lat, fill = time)) +
  scale_fill_gradient(name = 'N bird days', low="white", high="red", na.value = 'transparent', guide = "colourbar", limits=c(10, 30))+
  theme(axis.ticks = element_blank(),axis.text = element_blank(),axis.title = element_blank())+
  theme(strip.text.y = element_text(size = 20, colour = "black"), strip.text.x = element_text(size = 15, colour = "black"))+
  geom_point(data=turbs, aes(x=Longitude, y=Latitude), pch=16, col='darkolivegreen', size=0.5)
dev.off()

