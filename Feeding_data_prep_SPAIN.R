###########################################################################################################################
###### DATA PREPARATION TO ANTHROPOGENIC FEEDING MODELS OF RED KITES TO SPAIN ################
###########################################################################################################################
# modified by steffen.oppel@vogelwarte.ch in October 2023
# relies on pre-prepared data stes based on the following scripts:
# REKI_Eurokite_to_hourly.r
# Feeding_data_prep.R

library(tidyverse)
library(rnaturalearth)
library(sf)
library(amt)
library(dplyr, warn.conflicts = FALSE)
options(dplyr.summarise.inform = FALSE)
library(recurse)
library(lubridate)
library(leaflet)
library(leafgl)  ## needed to plot 1million points
#library(move)
library(data.table); setDTthreads(percent = 65)
sf_use_s2(FALSE) # deactivating spherical geometry s2


#define different coordinate systems
CH_LV95_coords = "+init=epsg:2056"
EU_coords = "+init=epsg:5635"
WGS84_coords ="+init=epsg:4326"
CH_LV03_coords ="+init=epsg:21781"


## set root folder for project
setwd("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/REKIFeeding")


### LOAD THE TRACKING DATA AND INDIVIDUAL SEASON SUMMARIES
trackingdata<-readRDS(file = "data/REKI_trackingdata_raw.rds")
#trackingdata<-fread("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/NestTool/REKI/output/02_preprocessing/03_milvus_combined.csv")
indseasondata<-fread("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/NestTool/REKI/output/01_validation/03_validation_combined.csv")
trackingdata<-trackingdata %>% left_join(indseasondata, by="year_id") %>%
  select(year_id,timestamp,long_wgs,lat_wgs,sex,age_cy) %>%
  rename(long=long_wgs, lat=lat_wgs)
trackingdata<-readRDS(file = "C:/Users/sop/OneDrive - Vogelwarte/gps_data/EUROKITE_all_tracks.rds")  %>%
  rename(sex=sex_of_bird) %>%
  mutate(age_cy=ifelse(is.na(as.numeric(age_of_bird_code)),ifelse(age_of_bird_code %in% c("J","F"),1,3),
                       as.numeric(age_of_bird_code))) %>%
  select(year_id,timestamp,long,lat,sex,age_cy) %>%
  bind_rows(trackingdata)

unique(trackingdata$sex)
head(trackingdata)

# DATA PREPARATION -------------------------------------------------------------
# keeping only the information of relevant locations in Spain

trackingdata <- trackingdata %>%
  filter(long<3.5) %>%
  filter(lat<43.5) %>%
  filter(!is.na(timestamp)) %>%
  filter(!is.na(long)) %>%
  filter(!is.na(year_id)) %>%
  filter(sex %in% c("m","f","male","female", "probably male","probably female")) %>%
  mutate(sex = ifelse(sex %in% c("m","male","probably male"),"m","f"))

dim(trackingdata)


# converting to metric CRS prior to estimating distances
track_sf <- trackingdata %>% 
  st_as_sf(coords = c("long", "lat"))
st_crs(track_sf) <- 4326

track_sf <- track_sf %>%
  st_transform(crs = 5635) %>%   ### LAEA EUROPE for Spain
  dplyr::mutate(long_eea = sf::st_coordinates(.)[,1],
              lat_eea = sf::st_coordinates(.)[,2])

head(track_sf)


# Creating a track to exclude nocturnal locations (later on)
# 3 mins
track_amt <- track_sf %>%
  mk_track(
    .x = long_eea,
    .y = lat_eea,
    .t = timestamp,
    sex= sex,
    age_cy=age_cy,
    id = year_id,
    crs = 5635
  ) %>%
  time_of_day(include.crepuscule = T) %>% # if F, crepuscule is considered as night
  arrange(id, t_)


### CALCULATE OTHER METRICS
track_amt$step_length<-amt::step_lengths(track_amt)       # include step lengths
track_amt$turning_angle<-amt::direction_abs(track_amt,append_last=T)      # include turning angles
track_amt$speed<-amt::speed(track_amt)      # include speed

head(track_amt)


# RECURSIONS FOR EACH LOCATION 50 m BUFFER--------------------------------------------------
# splitting track into a list with each single id grouped to an element
track_amt <- as.data.frame(track_amt)
track_amt_list <- split(track_amt, track_amt$id)

# calculating recursions
rm(trackingdata,track_sf)  ### clean up workspace and memory
gc()

track_amt_recurse <- lapply(track_amt_list, function(x)
  getRecursions(x = x[1:4], radius = 50, timeunits = "hours"))

# allocating recurse information to track data frame (15 mins)
track_amt$revisits <- NA
track_amt$residence_time <- NA
for (i in 1:length(track_amt_recurse)) {
  track_amt$revisits[track_amt$id == unique(track_amt$id)[i]] <-
    track_amt_recurse[[i]]$revisits
  track_amt$residence_time[track_amt$id == unique(track_amt$id)[i]] <-
    track_amt_recurse[[i]]$residenceTime
  
# CALCULATING FIRST AND LAST REVISIT AND DURATION AND TEMPORAL PERSISTENCE OF REVISITS -----------------------------------------------------------------
  tempout<-
    track_amt_recurse[[i]]$revisitStats %>%
    mutate(jday.ent=yday(entranceTime),jday.ex=yday(exitTime)) %>%
    group_by(coordIdx) %>%
    summarise(first=min(entranceTime, na.rm=T),
              last=max(exitTime, na.rm=T),
              meanFreqVisit=mean(timeSinceLastVisit, na.rm=T),
              n_days=length(unique(c(jday.ent,jday.ex)))) %>%
    mutate(TimeSpan=as.numeric(difftime(last,first,unit="days"))) %>%
    mutate(TempEven=n_days/TimeSpan) %>%
    mutate(meanFreqVisit=ifelse(is.na(meanFreqVisit),0,meanFreqVisit)) %>%   ## set the frequency of visit to 0 for locations never revisited
    mutate(TempEven=ifelse(n_days==1,1,TempEven)) %>%   ## set the evenness to 1 for locations never revisited on more than a single day
    select(meanFreqVisit,n_days,TimeSpan,TempEven)
  track_amt$meanFreqVisit[track_amt$id == unique(track_amt$id)[i]] <-tempout$meanFreqVisit
  track_amt$n_days[track_amt$id == unique(track_amt$id)[i]] <-tempout$n_days
  track_amt$TimeSpan[track_amt$id == unique(track_amt$id)[i]] <-tempout$TimeSpan
  track_amt$TempEven[track_amt$id == unique(track_amt$id)[i]] <-tempout$TempEven
}


# CALCULATING MOVING AVERAGE FOR SPEED AND ANGLE -----------------------------------------------------------------
rm(track_amt_recurse,track_amt_list,trackingdata)
track_amt <- track_amt %>% 
  arrange(id, t_) %>%
  group_by(id) %>%
  mutate(mean_speed=frollmean(speed,n=3,na.rm=T, align="center")) %>%
  mutate(mean_angle=frollmean(abs(turning_angle),n=5,na.rm=T, align="center")) %>%
  mutate(mean_speed=ifelse(is.na(mean_speed),speed,mean_speed)) %>%
  mutate(mean_angle=ifelse(is.na(mean_angle),turning_angle,mean_angle)) 
head(track_amt)



### READ IN SHAPEFILES OF RUBBISH DUMPS
## add layers that Jaume sent and from 

dumps <- st_read("data/Spain/Vertederos.shp") %>% 
  st_drop_geometry() %>%   ## for some reason the CRS is weird and wrong
  st_as_sf(coords = c("LONGITUD", "LATITUD"))
st_crs(dumps) <- 4326





### remove locations outside of study area
track_sf <- track_amt %>% 
  st_as_sf(coords = c("x_", "y_"))
st_crs(track_sf) <- 5635
head(track_sf)



### reformat data
track_sf <- track_sf %>% 
  ungroup() %>%
  st_transform(crs = 4326) %>%
  dplyr::mutate(long = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2]) %>%
  # mutate(long = unlist(map(track_sf$geometry,1)),
  #        lat = unlist(map(track_sf$geometry,2)))
  #st_drop_geometry() %>%
  rename(year_id=id)


# fwrite(as.data.frame(track_sf),"data/REKI_annotated_feeding.csv")
saveRDS(track_sf, file = "data/REKI_trackingdata_Spain.rds")
# head(track_sf)
# dim(track_sf)



#### look at a simple map ###

mF <- leaflet(options = leafletOptions(zoomControl = F)) %>% #changes position of zoom symbol
  setView(lng = mean(track_sf$long), lat = mean(track_sf$lat), zoom = 5) %>%
  htmlwidgets::onRender("function(el, x) {L.control.zoom({ 
                           position: 'bottomright' }).addTo(this)}"
  ) %>% #Esri.WorldTopoMap #Stamen.Terrain #OpenTopoMap #Esri.WorldImagery
  addProviderTiles("Esri.WorldImagery", group = "Satellite",
                   options = providerTileOptions(opacity = 0.4, attribution = F,minZoom = 5, maxZoom = 20)) %>%
  # addCircleMarkers(
  #   data=track_sf,
  #   radius = 1,
  #   stroke = TRUE, color = "red", weight = 0.8,
  #   fillColor = "red", fillOpacity = 0.8
  # ) %>%
  
  addGlPoints(data = track_sf, fillColor = "red", fillOpacity = 0.8) %>%
  
  addCircleMarkers(
    data=dumps,
    radius = 10,
    stroke = TRUE, color = "green", weight = 0.8,
    fillColor = "green", fillOpacity = 0.8
  ) %>%
  
  addScaleBar(position = "bottomright", options = scaleBarOptions(imperial = F))

mF




