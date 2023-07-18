###########################################################################################################################
###### USE TRACK2KBA TO IDENTIFY ANTHROPOGENIC FEEDING SITES OF RED KITES BASED ON GPS TRACKING DATA ################
###########################################################################################################################
# uses csv table instead of Movebank download to save time
# created by steffen.oppel@vogelwarte.ch in June 2023
# inspired by thesis of Nathalie Heiniger

## RUNNING THE ENTIRE SCRIPT TAKES >12 HRS

library(tidyverse)
library(track2KBA)
library(sf)
library(dplyr, warn.conflicts = FALSE)
options(dplyr.summarise.inform = FALSE)
library(recurse)
library(lubridate)
library(leaflet)
sf_use_s2(FALSE) # deactivating spherical geometry s2


## set root folder for project
setwd("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/Feeding")



# LOADING DATA -----------------------------------------------------------------
### LOAD THE ANNOTATED TRACKING DATA
trackingdata<-fread("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/REKI_annotated_tracking_data.csv") %>%
  mutate(Date=as.Date(timestamp)) %>%
  mutate(time=format(timestamp, format= "%H:%M:%S")) %>%
  filter(long_wgs>5.9) %>%
  filter(lat_wgs>45.8) %>%
  filter(long_wgs<10.6) %>%
  filter(lat_wgs<48)
head(trackingdata)

# nestdata need to read in coordinates of nests

colony <- trackingdata %>% group_by(year_id) %>%
  summarise(Latitude=first(nest_lat),Longitude=first(nest_long)) %>%
  filter(!is.na(Latitude)) %>%
  rename(ID=year_id)

dataGroup <- formatFields(
  dataGroup = trackingdata[trackingdata$year_id %in% colony$ID,], 
  fieldID   = "year_id", 
  fieldDateTime = "timestamp", 
  fieldLon  = "long_wgs", 
  fieldLat  = "lat_wgs"
)
head(dataGroup)

# ### FIX STUPID MIDNIGHT TIME ERROR
# dataGroup %>% filter(is.na(DateTime))
# 
# test<-dataGroup %>% filter(event_id %in% c(4831857242,4831857284))
# test<-dataGroup %>% filter(ID=="2016_1")

dataGroup <- dataGroup %>%
  mutate(DateTime=if_else(hour(DateTime) == minute(DateTime) &
                          minute(DateTime) == second(DateTime) &
                          second(DateTime) ==0,DateTime+second(1),DateTime))

### SPLIT INTO TRIPS BY EXCLUDING NEST LOCATIONS
trips <- tripSplit(
  dataGroup  = dataGroup,
  colony     = colony,
  innerBuff  = 0.4,      # 0.4 kilometers
  returnBuff = 1,
  duration   = 2.5,      # hours set to 2.5 because otherwise a single location is a trip
  rmNonTrip  = TRUE,
  nests=T
)
trips <- subset(trips, trips$Returns == "Yes" )


### tripSummary function takes a very long time to run (thousands of trips)
sumTrips <- tripSummary(trips = trips, colony = colony, nests=T)
#sumTrips
fwrite(sumTrips,"output/REKI_trip_summaries.csv")


### REPROJECT DATA 
tracks <- projectTracks(dataGroup = trips, projType = 'azim', custom=TRUE )
class(tracks)

#### NEED TO RUN SCALE FUNCTION TO INFORM KDE
findScale(
  tracks,
  scaleARS = TRUE,
  res = NULL,
  sumTrips = sumTrips,
  scalesFPT = NULL,
  peakWidth = 1,
  peakMethod = "first"
)




# ESTIMATE UDs
tracks <- tracks[tracks$ColDist > 0.3, ] # remove trip start and end points near nest
head(tracks)


### MEMORY LIMITATIONS - NEED TO SPLIT BY YEAR
outlist<-list()

for (y in unique(tracks$year)){
  outlist[[y-2015]]<- estSpaceUse(
    tracks = tracks[tracks$year==y,], 
    scale = 0.4, 
    levelUD = 50,
    res=0.5,
    polyOut = TRUE)
}


# ASSESS REPRESENTATIVITY
## this takes a long time and is not very informative since we deal with territorial birds, hence overlap is expected to be low!
outrep<-data.frame(year=unique(tracks$year),rep=0)

for (y in unique(tracks$year)){
  outrep[outrep$year==y,2] <- repAssess(
  tracks    = tracks[tracks$year==y,], 
  KDE       = outlist[[y-2015]]$KDE.Surface, 
  levelUD   = 50,
  iteration = 50, 
  bootTable = FALSE)
}
fwrite(outrep,"output/REKI_rep_assessment.csv")


# FIND SITES 
outsite<-list()

for (y in unique(tracks$year)){
outsite[[y-2015]] <- findSite(
  KDE = outlist[[y-2015]]$KDE.Surface,
  represent = 100,
  levelUD = 50,
  popSize = 500,     # 500 individual kites
  polyOut = TRUE
)
}

class(Site)




##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
########## COMBINE THE DATA FOR DIFFERENT YEARS  #############
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################

# READ IN FEEDING LOCATIONS
plot_feeders<-fread("data/Private_Feeders/private_feeders_upd2022.csv") %>% 
  filter(!is.na(coordX)) %>%
  st_as_sf(coords = c("coordX", "coordY"), crs=21781) %>%
  st_transform(crs = 4326) 


## split the multipolygon into separate polygons
cong_site_polygons <- st_cast(outsite[[1]], "POLYGON") %>%
  mutate(Year=2016)

for (y in unique(tracks$year)[-1]){
  cong_site_polygons <- st_cast(outsite[[y-2015]], "POLYGON") %>%
    mutate(Year=y) %>%
    bind_rows(cong_site_polygons)

}

cong_site_polygons


# CALCULATE DISTANCE TO FEEDING SITES
feed_site_distances<-st_distance(cong_site_polygons,plot_feeders, by_element=F,tolerance=5) 
cong_site_polygons <- cong_site_polygons %>%
  mutate(dist_feeder=apply(feed_site_distances,1,min)/1000, ### distance in km
         area=as.numeric(st_area(cong_site_polygons))/10000) %>%  ### area in ha
  arrange(desc(dist_feeder)) %>%
  filter(N_IND>0) %>%  ### exclude 0 ind polygons
  filter(area<100)   ### exclude massive box polygons of >100 ha


# SELECT SITES OF HIGH CONCENTRATION BUT NO NEARBY FEEDER
focal_sites <- cong_site_polygons %>%
  filter(N_IND>2) %>%
  filter(dist_feeder>1)





##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
########## PLOT THE MAP FOR POTENTIAL CONGREGATION SITES   #############
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################

########## CREATE A LEAFLET MAP OF PREDICTED FEEDING LOCATIONS ########################


## need to specify color palette 
# If you want to set your own colors manually:
pal <- colorNumeric(c("cornflowerblue","orange","firebrick"), c(4,5,6))
year.pal <- colorFactor(topo.colors(7), cong_site_polygons$Year)


m <- leaflet(options = leafletOptions(zoomControl = F)) %>% #changes position of zoom symbol
  htmlwidgets::onRender("function(el, x) {L.control.zoom({ 
                           position: 'bottomright' }).addTo(this)}"
  ) %>% #Esri.WorldTopoMap #Stamen.Terrain #OpenTopoMap #Esri.WorldImagery
  addProviderTiles("Esri.WorldImagery", group = "Satellite",
                   options = providerTileOptions(opacity = 0.6, attribution = F)) %>%
  addProviderTiles("OpenTopoMap", group = "Roadmap", options = providerTileOptions(attribution = F)) %>%  
  addLayersControl(baseGroups = c("Satellite", "Roadmap")) %>%  
  
  addCircleMarkers(
    data=plot_feeders,
    radius = 10,
    stroke = TRUE, color = "blue", weight = 1,
    fillColor = "green", fillOpacity = 0.8
  ) %>%
  
  addPolygons(
    data=cong_site_polygons,
    stroke = TRUE, color = ~year.pal(Year), weight = 1,
    fillColor = ~year.pal(Year), fillOpacity = 0.3
  ) %>%
  
  addPolygons(
    data=focal_sites,
    stroke = TRUE, color = "red", weight = 1,
    fillColor = "red", fillOpacity = 0.5
  ) %>%
  
  addLegend(     # legend for date (viridis scale)
    position = "topleft",
    pal = year.pal,
    values = cong_site_polygons$Year,
    opacity = 1,
    title = "Year"
  ) %>%
  
addScaleBar(position = "bottomright", options = scaleBarOptions(imperial = F))

m


htmltools::save_html(html = m, file = "C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/Feeding/output/track2kba_congregation_sites.html")

st_write(cong_site_polygons,"output/REKI_track2kba_congregation_areas.kml")



save.image("output/REKI_track2KBA_output.RData")




