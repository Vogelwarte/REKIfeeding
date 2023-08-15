###########################################################################################################################
###### DATA PREPARATION TO IDENTIFY ANTHROPOGENIC FEEDING SITES OF RED KITES BASED ON GPS TRACKING DATA ################
###########################################################################################################################
# original script by Nathalie Heiniger (MSc thesis 2020)
# uses csv table instead of Movebank download to save time
# modified by steffen.oppel@vogelwarte.ch in June 2023

### DATA INPUT IS A PERSISTENT HEADACHE
### originally tried Movebank download, but failed, then thinned the data but recursions fail with std::bad_alloc (probably memory overload)
### reverted to using hourly data that were pre-prepared by Ursin Beeli, as no way forward using other data


### FINALISED on 30 JUNE 2023 with hourly data


library(tidyverse)
library(rnaturalearth)
library(sf)
library(amt)
library(dplyr, warn.conflicts = FALSE)
options(dplyr.summarise.inform = FALSE)
library(recurse)
library(lubridate)
#library(move)
library(data.table); setDTthreads(percent = 65)
sf_use_s2(FALSE) # deactivating spherical geometry s2


## set root folder for project
setwd("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/REKIFeeding")




###  LOADING DATA -----------------------------------------------------------------
# ATTEMPT TO LOAD DATA FROM MOVEBANK ABANDONED ON 20 June 2023 because of memory limits
# saved Movebank data to csv file and read in csv file on 28 June 2023
# instead using data given to Ursin Beeli for NestTool analysis
movemil_GSM <- read.csv("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/NestTool/REKI/output/02_preprocessing/01_milvus_gsm.csv")
movemil_Milsar <- read.csv("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/NestTool/REKI/output/02_preprocessing/02_milvus_milsar.csv")
dat_GSM <- as.data.frame(movemil_GSM) %>%
  dplyr::select(individual.local.identifier,timestamp,year_id,year,month,day,long_wgs,lat_wgs,long_eea,lat_eea) %>%
  mutate(timestamp=as.POSIXct(timestamp,format ="%Y-%m-%d %H:%M:%S", tz = "UTC")) %>%
  dplyr::filter(!is.na(timestamp))
myDF <- as.data.frame(movemil_Milsar) %>%
  dplyr::select(individual.local.identifier,timestamp,year_id,year,month,day,long_wgs,lat_wgs,long_eea,lat_eea) %>%
  mutate(timestamp=as.POSIXct(timestamp,format ="%Y-%m-%d %H:%M:%S", tz = "UTC")) %>%
  dplyr::filter(!is.na(timestamp)) %>%
  bind_rows(dat_GSM) %>%
  arrange(year_id,timestamp)

############# THIS CHUNK LAST RUN ON 28 JUNE 2023 ##############
# movemil_Milsar <- fread("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/Feeding/data/Milvusmilvus_Milsar_SOI_final.csv")    ### insufficient memory to read in those data
# dim(movemil_Milsar)
# 
# 
# #read in your data directly from movebank
# username <- "Steffen"
# password <- "xxxxx"
# curl <- movebankLogin(username=username,  password=password)
# #
# # #read in your data directly from movebank
# # getMovebankID(study="Milvusmilvus_Milsar_SOI_final", login=curl)
# # unique(getMovebankSensors(study=1356790386, login=curl)$sensor_type_id)
# # getMovebankSensors(login=curl)
# movemil_GSM <- getMovebankLocationData(study="Milvusmilvus_GSM_SOI", sensorID=653, login=curl)
# # movemil_Milsar <- getMovebankLocationData(study="Milvusmilvus_Milsar_SOI_final", sensorID=653, login=curl)  ## does not work due to Individual(s): '423' do(es) not have data for one or more of the selected sensor(s).
# 
# 
# dat_GSM <- as.data.frame(movemil_GSM) %>%
#   dplyr::select(individual.local.identifier,timestamp,location.long,location.lat) %>%
#   rename(id=individual.local.identifier, lat=location.lat,long=location.long) %>%
#   mutate(year=year(timestamp),
#           month=month(timestamp),
#           day=day(timestamp),
#           yday=yday(timestamp),
#          year_id=paste(year,id,sep="_")) %>%
#   dplyr::filter(!is.na(timestamp))
# 
# 
# myDF <- as.data.frame(movemil_Milsar) %>%
#   dplyr::select(`individual-local-identifier`,timestamp,`location-long`,`location-lat`) %>%
#   rename(id=`individual-local-identifier`, lat=`location-lat`,long=`location-long`) %>%
#   mutate(year=year(timestamp),
#          month=month(timestamp),
#          day=day(timestamp),
#          yday=yday(timestamp),
#          year_id=paste(year,id,sep="_")) %>%
#   dplyr::filter(!is.na(timestamp)) %>%
#   bind_rows(dat_GSM) %>%
#   arrange(year_id,timestamp)
# dim(myDF)
trackingdata <- myDF[!duplicated(paste0(myDF$timestamp,myDF$id)),] ## this is to exclude duplicated timestamps (if present)
saveRDS(trackingdata, file = "data/REKI_trackingdata_raw.rds", version=3)
dim(trackingdata)
rm(myDF,dat_GSM,movemil_GSM,movemil_Milsar)

### LOAD THE TRACKING DATA AND INDIVIDUAL SEASON SUMMARIES
trackingdata<-readRDS(file = "data/REKI_trackingdata_raw.rds")



# DATA PREPARATION -------------------------------------------------------------
# keeping only the information of relevant locations in Switzerland

if("long_wgs" %in% names(trackingdata)){
  trackingdata<-trackingdata %>% rename(long=long_wgs, lat=lat_wgs)
}
trackingdata <- trackingdata %>%
  filter(long>5.9) %>%
  filter(lat>45.8) %>%
  filter(long<10.6) %>%
  filter(lat<48) %>%
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


# Creating a track to exclude nocturnal locations (later on)
# 3 mins
track_amt <- track_sf %>%
  mk_track(
    .x = long_eea,
    .y = lat_eea,
    .t = timestamp,
    id = year_id,
    crs = 3035
  ) %>%
  time_of_day(include.crepuscule = T) %>% # if F, crepuscule is considered as night
  arrange(id, t_)



#### SWITCHED OFF ON 30 JUNE AFTER REVERTING TO HOURLY DATASET
# turn into a 30 min resolution track - this doesn't work with multiple animals as 'negative time differences occur'
# amt_resampled<-tibble()
# for(i in unique(track_amt$id)){
#   singletrack<-track_amt %>% filter(id==i)
#   single_amt <- track_resample(singletrack, rate = minutes(60), tolerance = minutes(10), start = 1)
#   amt_resampled<-amt_resampled %>% bind_rows(single_amt)
#   ### troubleshooting ###
#   # t_ <- as.numeric(singletrack$t_)
#   # summary(diff(t_))
#   # neg.time<-which(diff(t_) < 0)
#   # track_amt[c(neg.time[1:3]),] %>% arrange(id, t_)
#   # test<-singletrack[2500:2510,]
# }
# 
# if(dim(track_amt)[1]-dim(amt_resampled)[1]>0){
#   track_amt <- amt_resampled %>%
#     mk_track(
#       .x = x_,
#       .y = y_,
#       .t = t_,
#       id = id,
#       crs = 3035
#     ) %>%
#     time_of_day(include.crepuscule = T) %>% # if F, crepuscule is considered as night
#     arrange(id, t_)
#   rm(amt_resampled,singletrack,single_amt)
# }

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
# rm(trackingdata,track_sf,buildings,forest)  ### clean up workspace and memory
# gc()

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


# COMBINING DATA WITH FORESTS AND BUILDINGS -----------------------------------------------------------------



#trackingdata<-fread("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/NestTool/REKI/output/02_preprocessing/03_milvus_combined.csv")
indseasondata<-fread("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/NestTool/REKI/output/01_validation/03_validation_combined.csv")
nestdata<-fread("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/NestTool/REKI/data/Basic_nest_list_2015_2022.csv")

#define different coordinate systems
CH_LV95_coords = "+init=epsg:2056"
EU_coords = "+init=epsg:3035"
WGS84_coords ="+init=epsg:4326"
CH_LV03_coords ="+init=epsg:21781"



### READ IN SHAPEFILES OF FOREST AND BUILDINGS
buildings <- st_read("data/Buildings/tlm_buildings_studyareaExtra_size65_buff_50m_sf_singlepoly.shp", stringsAsFactors=FALSE) %>%
  st_transform(crs = 3035) 
forest <- st_read("data/Forest/vec25_forest_buff_20m_studyarea.shp", stringsAsFactors=FALSE) %>%
  st_transform(crs = 3035) %>%
  select(AREA)


### READ IN FEEDERS
FEEDERS<- fread("data/Private_Feeders/private_feeders_upd2022.csv") %>% 
  filter(!is.na(coordX)) %>%
  st_as_sf(coords = c("coordX", "coordY"))
st_crs(FEEDERS) <- 21781
FEEDER_buff<- FEEDERS %>% st_transform(crs = 3035) %>%
  st_buffer(dist=50) %>%
  select(ID,type_of_food,frequency) %>%
  rename(feeder_id=ID)


### EXPERIMENTAL FEEDING STATIONS
EXPFEEDERS<- fread("data/experimental_feeding.csv") %>%
  mutate(long=as.numeric(lon)) %>%
  # group_by(nest_name) %>%
  # summarise(lat=mean(lat,na.rm=T),long=mean(long, na.rm=T)) %>%
  filter(!is.na(long)) %>%
  filter(!is.na(lat)) %>%
  st_as_sf(coords = c("long", "lat"))
st_crs(EXPFEEDERS) <- 4326
EXPFEEDERS_buff<- EXPFEEDERS %>% st_transform(crs = 3035) %>%
  st_buffer(dist=50) %>%
  select(nest_name,Timepoint,food_placed) %>%
  rename(feeder_id=nest_name)


### FEEDING PLATFORMS
PLATFORMS<- fread("data/feeding_platforms_15_16.csv") %>%
  filter(!is.na(x)) %>%
  mutate(start_date=dmy(start_date), end_date=dmy(end_date)) %>%
  st_as_sf(coords = c("x", "y"))
st_crs(PLATFORMS) <- 21781
PLATFORMS_buff<- PLATFORMS %>% st_transform(crs = 3035) %>%
  st_buffer(dist=50) %>%
  select(Name,year,n_event,start_date,end_date) %>%
  rename(feeder_id=Name)
PLATFORMS_buff


### SUMMARY OF NUMBER OF FEEDERS
## this makes only sense if the experimental feeders have been grouped
dim(FEEDERS)[1] +
dim(EXPFEEDERS)[1] +
dim(PLATFORMS)[1]


### NESTS
nests<- nestdata %>% dplyr::select(nest_name,tree_spec,latitude,longitude) %>%
  st_as_sf(coords = c("longitude", "latitude"))
st_crs(nests) <- 4326
nests<-nests %>% st_transform(crs = 3035)
nests_buff<- nests %>%
  st_buffer(dist=50)
nests_buff



### remove locations outside of study area
track_sf <- track_amt %>% 
  st_as_sf(coords = c("x_", "y_"))
st_crs(track_sf) <- 3035
head(track_sf)
track_sf <- track_sf %>% st_crop(x=track_sf, y=st_bbox(forest))



### spatial joins with forests, buildings and feeders
head(track_sf)
track_sf <- track_sf %>%
  st_join(forest,
          join = st_intersects, 
          left = TRUE) %>%
  st_join(buildings,
          join = st_intersects, 
          left = TRUE) %>%
  st_join(FEEDER_buff,
          join = st_intersects, 
          left = TRUE) %>%
  st_join(nests_buff,
          join = st_intersects, 
          left = TRUE) %>%
  mutate(BUILD=ifelse(is.na(build_id),0,1),
         NEST=ifelse(is.na(nest_name),0,1),
         FOREST=ifelse(is.na(AREA),0,1),
         FEEDER=ifelse(is.na(feeder_id),"NO","YES")) %>%   ### FOR PUBLIC FEEDERS SPATIAL OVERLAP IS YES OR NO
  mutate(FEED_ID=feeder_id) %>%
  select(-feeder_id) %>%
  st_join(EXPFEEDERS_buff,
          join = st_intersects, 
          left = TRUE) %>%
  mutate(FEEDER=ifelse(is.na(feeder_id),FEEDER,
                       ifelse(t_ %within% interval(Timepoint,Timepoint+days(30)),"YES",FEEDER))) %>%   ### FOR EXPERIMENTAL FEEDERS NEED TEMPORAL OVERLAP TO WITHIN A MONTH OF Timepoint
  mutate(FEED_ID=ifelse(is.na(feeder_id),FEED_ID,feeder_id)) %>%
  select(-feeder_id) %>%
  
  st_join(PLATFORMS_buff,
          join = st_intersects, 
          left = TRUE) %>%
  mutate(FEEDER=ifelse(is.na(feeder_id),FEEDER,
                       ifelse(t_ %within% interval(start_date,end_date+days(30)),"YES",FEEDER))) %>%   ### FOR FEEDING PLATFORMS NEED TEMPORAL OVERLAP TO WITHIN A MONTH OF end time
  mutate(FEED_ID=ifelse(is.na(feeder_id),FEED_ID,feeder_id)) %>%
  select(-year,-n_event,-start_date,-end_date,-Timepoint,-food_placed,-nest_name,-tree_spec,-feeder_id) %>%
  rename(forest_size=AREA)

head(track_sf)

st_bbox(FEEDER_buff)
st_bbox(forest)
st_bbox(buildings)

table(track_sf$FEEDER)
table(track_sf$FOREST)
table(track_sf$BUILD)


### calculate distance to nearest nest
nest_site_distances<-st_distance(track_sf,nests) 
track_sf <- track_sf %>%
  mutate(dist_nest=apply(nest_site_distances,1,min)/1000)  ### distance in km




### join data with individual info
track_sf <- track_sf %>% 
  st_transform(crs = 4326) %>%
  dplyr::mutate(long = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2]) %>%
  # mutate(long = unlist(map(track_sf$geometry,1)),
  #        lat = unlist(map(track_sf$geometry,2)))
  st_drop_geometry() %>%
  rename(year_id=id) %>%
  left_join(indseasondata, by="year_id")


fwrite(as.data.frame(track_sf),"data/REKI_annotated_feeding.csv")
saveRDS(track_sf, file = "data/REKI_trackingdata_annotated.rds")
head(track_sf)
dim(track_sf)
