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

### REVISED on 28 May 2024 with newly downloaded tracking data
### revision in different branch to work with HIGH RES data - not at hourly resolution

## ADDED individual life history data from all birds on 3 June 2024

### REVISED ON 16 Jan 2025 to include UHF DATA 
## saved as new file to include survey feeders as well

## CHANGED ON 20 JANUARY 2025 to read in building and forest layer across ALL of Switzerland (used to be just study area)
## spatial overlay failed after 48 hrs - script never completed running
## tried new strategy on 22 Jan 2025 to reduce tracking data

rm(list=ls())
library(tidyverse)
library(rnaturalearth)
library(sf)
library(amt)
library(suncalc)
library(dplyr, warn.conflicts = FALSE)
options(dplyr.summarise.inform = FALSE)
library(recurse)
library(readxl)
library(lubridate)
#library(move)
library(data.table); setDTthreads(percent = 65)
sf_use_s2(FALSE) # deactivating spherical geometry s2
library(tictoc)

## set root folder for project
#setwd("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/REKIfeeding")
#setwd("C:/STEFFEN/OneDrive - Vogelwarte/REKI/Analysis/REKIfeeding")




# ###  LOADING DATA -----------------------------------------------------------------
## outsourced to script REKI_data_download_high_res.r
## updated on 15 Jan 2025 with script ANALYSIS/DataPrep/REKI_move2bank_download_filter.r

### LOAD THE TRACKING DATA AND INDIVIDUAL SEASON SUMMARIES
#trackingdata<-readRDS(file = "data/REKI_trackingdata_raw.rds")
#trackingdata<-readRDS(file = "data/REKI_trackingdata_raw2024.rds") %>%
trackingdata<-readRDS(file = "C:/Users/sop/OneDrive - Vogelwarte/General/DATA/REKI_filtered_15min_ALLData_15Jan2025.rds")

### LOAD INDIVIDUAL LIFE HISTORIES
#trackingdata<-fread("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/NestTool/REKI/output/02_preprocessing/03_milvus_combined.csv")
#indseasondata<-fread("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/NestTool/REKI/output/01_validation/03_validation_combined.csv")
indseasondata<-read_excel("C:/Users/sop/OneDrive - Vogelwarte/General/DATA/Individual_life_history_2015-2023.xlsx", sheet="Individual_life_history_2015-20") %>% # updated on 3 June 2024 to include birds from 2022
  dplyr::select(bird_id,ring_number,tag_year,sex_compiled, age, hatch_year) %>%
  rename(ring_id=ring_number) %>%
  mutate(hatch_year=if_else(is.na(as.numeric(hatch_year)),tag_year-3,as.numeric(hatch_year)))
nestdata<-fread("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/NestTool/REKI/data/Basic_nest_list_2015_2022.csv")



# DATA PREPARATION -------------------------------------------------------------
# converting to metric CRS prior to estimating distances
track_sf <- trackingdata %>%
  select(-deployment_id) %>%
  filter(long>5.9) %>%
  filter(lat>45.8) %>%
  filter(long<10.6) %>%
  filter(lat<48) %>%
  filter(!is.na(timestamp)) %>%
  filter(!is.na(long)) %>%
  mutate(year_id=paste(year(timestamp),bird_id, sep="_")) %>%
  filter(!is.na(year_id)) %>%
  # left_join(indseasondata, by=c("bird_id")) %>%
  # mutate(age_cy=year(timestamp)-hatch_year) %>%
  # rename(sex=sex_compiled) %>%
  dplyr::select(year_id,sex,age_cy,timestamp,geometry) %>%
  st_transform(crs = 3035) %>%
  ungroup() %>%
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
    sex=sex,
    age_cy=age_cy,
    crs = 3035,
    order_by_ts=TRUE
  ) %>%
  time_of_day(include.crepuscule = T) %>% # if F, crepuscule is considered as night
  arrange(id, t_)

## clean up workspace
rm(trackingdata,track_sf)
gc()


### FILTER OUT DUPLICATE LOCATIONS WITHIN 60 seconds
## this step did not complete in 24 hours on 8 Jan 2025 - commented out because the input data are already filter to 15 min
# track_amt<-track_amt %>%
# 	mutate(DOP=1) %>%
# 	flag_duplicates(.,gamma=seconds(60),time_unit="sec",DOP="DOP")


### CALCULATE OTHER METRICS
track_amt$step_length<-amt::step_lengths(track_amt)       # include step lengths
track_amt$turning_angle<-abs(amt::direction_rel(track_amt,append_last=T, full_circle = FALSE,lonlat = FALSE))      # include RELATIVE turning angles - changed from absolute
track_amt$turning_angle<-ifelse(is.na(track_amt$turning_angle),0,track_amt$turning_angle)  ## replace NA turning angles (caused by 0 distance) with 0
track_amt$speed<-amt::speed(track_amt)      # include speed
track_amt$locid<-seq_along(track_amt$t_)

## check how speed is calculated by hand
#track_amt %>% mutate(dt=as.numeric(difftime(dplyr::lead(t_),t_, units="sec"))) %>%
#	mutate(speed2=step_length/dt)

### ABOVE METRICS NEED TO BE SET TO NA FOR THE LAST POSITION OF EACH INDIVIDUAL - THEY ARE NOT CALCULATED FOR YEAR_ID and therefore calculate distances between different animals
lastlocs<- track_amt %>% st_drop_geometry() %>%
	arrange(id,t_) %>%
	group_by(id) %>%
	summarise(last=max(t_), lastloc=max(locid))
track_amt$speed[track_amt$locid %in% lastlocs$lastloc]<-NA
track_amt$step_length[track_amt$locid %in% lastlocs$lastloc]<-NA
track_amt$turning_angle[track_amt$locid %in% lastlocs$lastloc]<-NA

head(track_amt)
hist(track_amt$turning_angle*(180/pi))
range(track_amt$speed, na.rm=T)


# RECURSIONS FOR EACH LOCATION 50 m BUFFER--------------------------------------------------
# splitting track into a list with each single id grouped to an element------
# this causes a memory limit error, so try and do it in a loop

track_amt <- as.data.frame(track_amt)
#track_amt_list <- split(track_amt, track_amt$id)



### RECURSIONS IN A LOOP - takes 2 hrs --------------------------------
# calculating recursions
rm(trackingdata,track_sf,buildings,forest)  ### clean up workspace and memory
gc()

# track_amt_recurse <- lapply(track_amt_list, function(x)
#   getRecursions(x = x[1:4], radius = 50, timeunits = "hours"))
tic()
track_amt$revisits <- NA
track_amt$residence_time <- NA
for (i in unique(track_amt$id)){
  x<-track_amt %>% filter(id==i)
  xr<-getRecursions(x = x[1:4], radius = 50, timeunits = "hours")
  
  track_amt$revisits[track_amt$id == i] <- xr$revisits
  track_amt$residence_time[track_amt$id == i] <- xr$residenceTime
  
  # CALCULATING FIRST AND LAST REVISIT AND DURATION AND TEMPORAL PERSISTENCE OF REVISITS -----------------------------------------------------------------
  tempout<-
    xr$revisitStats %>%
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
  track_amt$meanFreqVisit[track_amt$id == i] <-tempout$meanFreqVisit
  track_amt$n_days[track_amt$id == i] <-tempout$n_days
  track_amt$TimeSpan[track_amt$id == i] <-tempout$TimeSpan
  track_amt$TempEven[track_amt$id == i] <-tempout$TempEven
  rm(tempout,x,xr)
}
toc()



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
dim(track_amt)

# INTERMEDIATE SAVE TO NOT RERUN RECURSIONS UPON RESTART -----------------------------------------------------------------
saveRDS(track_amt,"data/track_amt2025.rds")
#track_amt<-readRDS("data/track_amt2025.rds")

# COMBINING DATA WITH FORESTS AND BUILDINGS -----------------------------------------------------------------

#define different coordinate systems
CH_LV95_coords = "+init=epsg:2056"
EU_coords = "+init=epsg:3035"
WGS84_coords ="+init=epsg:4326"
CH_LV03_coords ="+init=epsg:21781"



### READ IN SHAPEFILES OF FOREST AND BUILDINGS - INCLUDE ALL OF SWITZERLAND OUTSIDE OF STUDY AREA
buildings <- st_read("data/Buildings/tlm_buildings_size65_buff_50m_dissolved.shp", stringsAsFactors=FALSE) %>%
  st_transform(crs = 3035) 
forest <- st_read("data/Forest/vec25_forest_buff_20m_CH.shp", stringsAsFactors=FALSE) %>%
  st_transform(crs = 3035) %>%
  select(AREA)


### READ IN FEEDERS
FEEDERS<- fread("data/Private_Feeders/private_feeders_upd2022.csv") %>% 
  filter(!is.na(coordX)) %>%
  st_as_sf(coords = c("coordX", "coordY"))

### ADD SURVEY DATA FROM EVA CEREGHETTI AND FIONA PELLET  #############
SURVEYS<- fread("C:/Users/sop/OneDrive - Vogelwarte/General/MANUSCRIPTS/AnthropFeeding/DataArchive/REKI_validation_feeders.csv") %>%
  select(-geometry) %>%
  st_as_sf(coords = c("long", "lat"), crs=4326) %>%
  filter(FEEDER_surveyed==1) %>%
  mutate(Type="surveys") %>%
  select(ID,Type,geometry)
SURVEYS_buff<- SURVEYS %>% st_transform(crs = 3035) %>%
  st_buffer(dist=50) %>%
  select(ID,Type,geometry) %>%
  mutate(frequency="unknown") %>%
  rename(feeder_id=ID,type_of_food=Type)

st_crs(FEEDERS) <- 21781
FEEDER_buff<- FEEDERS %>% st_transform(crs = 3035) %>%
  st_buffer(dist=50) %>%
  select(ID,type_of_food,frequency) %>%
  rename(feeder_id=ID) %>%
  mutate(feeder_id=as.character(feeder_id)) %>%
  bind_rows(SURVEYS_buff)


### EXPERIMENTAL FEEDING STATIONS
### THE PROBLEM IS THE SPATIAL DUPLICATION AT MULTIPLE TIME POINTS; SO WE TAKE THE MEDIAN FOR EACH NEST NAME
EXPFEEDERS_buff<- fread("data/experimental_feeding.csv") %>%
  dplyr::filter(!event_id %in% c(5610,4413)) %>%  ## remove 2 offending stations with duplicate coordinates that cause problems when overlaying with locations
  mutate(long=as.numeric(lon)) %>%
  filter(!is.na(long)) %>%
  filter(!is.na(lat)) %>%
  group_by(nest_name) %>%
  summarise(start=min(Timepoint), end=max(Timepoint), lat=mean(lat,na.rm=T),long=mean(long, na.rm=T), food_placed=median(food_placed, na.rm=T)) %>%
  st_as_sf(coords = c("long", "lat"), crs=4326) %>%
  ungroup() %>%
  st_transform(crs = 3035) %>%
  st_buffer(dist=50) %>%
  select(nest_name,start, end,food_placed) %>%
  rename(feeder_id=nest_name)
st_difference(EXPFEEDERS_buff)


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
dim(EXPFEEDERS_buff)[1] +
dim(PLATFORMS)[1] +
dim(SURVEYS)[1]


### NESTS
nests<- nestdata %>% dplyr::select(nest_name,tree_spec,latitude,longitude) %>%
  st_as_sf(coords = c("longitude", "latitude"))
st_crs(nests) <- 4326
nests<-nests %>% st_transform(crs = 3035)
nests_buff<- nests %>%
  st_buffer(dist=50)
nests_buff



### recreate a spatial feature
track_sf <- track_amt %>% 
  st_as_sf(coords = c("x_", "y_"))  #
st_crs(track_sf) <- 3035
head(track_sf)
#track_sf <- track_sf %>% st_crop(x=track_sf, y=st_bbox(forest)) ## this reduces n locs from 10 mio to 8 mio
dim(track_sf)
st_difference(FEEDER_buff)

### spatial joins with forests, buildings and feeders
## this operation bizarrely ADDs duplicate rows where polygons overlap
## need to include st_difference() for all layers to prevent this: https://gis.stackexchange.com/questions/351429/sf-st-intersection-returning-duplicate-features
head(track_sf)
track_sf <- track_sf %>%
  st_join(st_difference(forest),
          join = st_intersects,
          left = TRUE) %>%
  st_join(st_difference(buildings),
          join = st_intersects,
          left = TRUE) %>%
  st_join(st_difference(FEEDER_buff),
          join = st_intersects,
          left = TRUE) %>%
  st_join(st_difference(nests_buff),
          join = st_intersects,
          left = TRUE) %>%
  mutate(BUILD=ifelse(is.na(build_id),0,1),
         NEST=ifelse(is.na(nest_name),0,1),
         FOREST=ifelse(is.na(AREA),0,1),
         FEEDER=ifelse(is.na(feeder_id),"NO","YES")) %>%   ### FOR PUBLIC FEEDERS SPATIAL OVERLAP IS YES OR NO
  mutate(FEED_ID=feeder_id) %>%
  select(-feeder_id) %>%
  # st_join(st_difference(SURVEYS_buff),
  #         join = st_intersects,
  #         left = TRUE) %>%
  # mutate(FEEDER=ifelse(is.na(feeder_id),FEEDER,"YES")) %>%   
  # select(-feeder_id) %>%
  st_join(st_difference(EXPFEEDERS_buff),
          join = st_intersects, 
          left = TRUE) %>%
  mutate(FEEDER=ifelse(is.na(feeder_id),FEEDER,
                       ifelse(t_ %within% interval(start,end+days(30)),"YES",FEEDER))) %>%   ### FOR EXPERIMENTAL FEEDERS NEED TEMPORAL OVERLAP TO WITHIN A MONTH OF Timepoint - changed to start and end of interval
  mutate(FEED_ID=ifelse(is.na(feeder_id),FEED_ID,feeder_id)) %>%
  select(-feeder_id) %>%
  st_join(st_difference(PLATFORMS_buff),
          join = st_intersects, 
          left = TRUE) %>%
  mutate(FEEDER=ifelse(is.na(feeder_id),FEEDER,
                       ifelse(t_ %within% interval(start_date,end_date+days(30)),"YES",FEEDER))) %>%   ### FOR FEEDING PLATFORMS NEED TEMPORAL OVERLAP TO WITHIN A MONTH OF end time
  mutate(FEED_ID=ifelse(is.na(feeder_id),FEED_ID,feeder_id)) %>%
  select(-year,-n_event,-start_date,-end_date,-start,-end,-food_placed,-nest_name,-tree_spec,-feeder_id) %>%
  rename(forest_size=AREA)

head(track_sf)
dim(track_sf)
st_bbox(FEEDER_buff)
st_bbox(forest)
st_bbox(buildings)

table(track_sf$FEEDER)
table(track_sf$FOREST)
table(track_sf$BUILD)


### calculate distance to nearest nest
## this causes memory allocation error, so need to do it in a loop, which takes 3278.25 seconds
## NEED TO FACTOR IN YEAR AS WELL ??
rm(track_amt,EXPFEEDERS)
gc()

# nest_site_distances<-st_distance(track_sf,nests) 
# track_sf <- track_sf %>%
#   mutate(dist_nest=apply(nest_site_distances,1,min)/1000)  ### distance in km

tic()
track_sf$dist_nest <- NA

for (i in unique(track_sf$id)){
  x<-track_sf %>% filter(id==i)
  x_distances<-st_distance(x,nests)
  track_sf$dist_nest[track_sf$id==i] <- apply(x_distances,1,min)/1000  ### distance in km
}
toc()


# EXPORT THE ANNOTATED DATA FOR ALL OF SUI -----------------------------------------------------------------



### export data in two different formats
track_out <- track_sf %>% 
  st_transform(crs = 4326) %>%
  dplyr::mutate(long = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2]) %>%
  # mutate(long = unlist(map(track_sf$geometry,1)),
  #        lat = unlist(map(track_sf$geometry,2)))
  st_drop_geometry() %>%
  rename(year_id=id) #%>%
  #left_join(indseasondata, by="year_id")

dim(track_sf)
fwrite(as.data.frame(track_out),"data/REKI_annotated_feeding2025_CH.csv")
saveRDS(track_sf, file = "data/REKI_trackingdata_annotated2025_CH.rds")


