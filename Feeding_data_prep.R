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
#setwd("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/REKIFeeding")




# ###  LOADING DATA -----------------------------------------------------------------
## outsourced to script REKI_data_download_high_res.r
## updated on 15 Jan 2025 with script ANALYSIS/DataPrep/REKI_move2bank_download_filter.r

### LOAD THE TRACKING DATA AND INDIVIDUAL SEASON SUMMARIES
#trackingdata<-readRDS(file = "data/REKI_trackingdata_raw.rds")
#trackingdata<-readRDS(file = "data/REKI_trackingdata_raw2024.rds") %>%
trackingdata<-readRDS(file = "C:/Users/sop/OneDrive - Vogelwarte/General/DATA/REKI_filtered_15min_ALLData_15Jan2025.rds")
# uhfdata<-fread("C:/Users/sop/OneDrive - Vogelwarte/General/DATA/UHF_data/REKI_UHF_data_locations.csv") %>%
#   filter(!is.na(Latitude)) %>%
#   select(bird_id,timestamp,Latitude, Longitude) %>%
#   rename(long=Longitude, lat=Latitude)
# ## check whether uhfdata are in trackingdata - NO THEY ARE NOT
# # trackingdata %>% filter(ymo_id=="101_2016_06") %>% filter(day(timestamp)==22) %>% filter(hour(timestamp)==7)

# trackingdata<-trackingdata %>%
#   ungroup() %>%
#   select(-deployment_id) %>%
#   # st_transform(4326) %>%
#   # dplyr::mutate(long = sf::st_coordinates(.)[,1],
#   #               lat = sf::st_coordinates(.)[,2]) %>%
#   st_drop_geometry() %>%
#   #rename(bird_id=id) %>%
#   select(bird_id,timestamp,lat,long) %>%
#   #bind_rows(uhfdata) %>%
#   arrange(bird_id, timestamp)
# # rm(uhfdata)
# # gc()
### LOAD INDIVIDUAL LIFE HISTORIES
#trackingdata<-fread("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/NestTool/REKI/output/02_preprocessing/03_milvus_combined.csv")
#indseasondata<-fread("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/NestTool/REKI/output/01_validation/03_validation_combined.csv")
indseasondata<-read_excel("C:/Users/sop/OneDrive - Vogelwarte/General/DATA/Individual_life_history_2015-2023.xlsx", sheet="Individual_life_history_2015-20") %>% # updated on 3 June 2024 to include birds from 2022
  dplyr::select(bird_id,ring_number,tag_year,sex_compiled, age, hatch_year) %>%
  rename(ring_id=ring_number) %>%
  mutate(hatch_year=if_else(is.na(as.numeric(hatch_year)),tag_year-3,as.numeric(hatch_year)))
nestdata<-fread("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/NestTool/REKI/data/Basic_nest_list_2015_2022.csv")



# DATA PREPARATION -------------------------------------------------------------
# # keeping only the information of relevant locations in Switzerland
# 
# # if("long_wgs" %in% names(trackingdata)){
# #   trackingdata<-trackingdata %>% rename(long=long_wgs, lat=lat_wgs)
# # }
# trackingdata <- trackingdata %>%
#   filter(long>5.9) %>%
#   filter(lat>45.8) %>%
#   filter(long<10.6) %>%
#   filter(lat<48) %>%
#   filter(!is.na(timestamp)) %>%
#   filter(!is.na(long)) %>%
#   mutate(year_id=paste(year(timestamp),bird_id, sep="_")) %>%
#   filter(!is.na(year_id))
# 
# dim(trackingdata)
# 
# 
# # filling in gaps in age and sex and making age uniform (in years)
# 
# trackingdata <- trackingdata %>%
#   left_join(indseasondata, by=c("bird_id")) %>%
#   mutate(age_cy=year(timestamp)-hatch_year) %>%
#   rename(sex=sex_compiled)
# 
# dim(trackingdata)
# summary(trackingdata$age_cy)
# head(trackingdata)

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


########### THIS SECTION COMMENTED OUT ON 16 JAN 2025 because data were already speed filtered before ############################

# track_amt %>% filter(speed<0)
# track_amt %>% filter(id=="2015_139") %>% filter(yday(t_) > yday(ymd("2015-10-27"))) %>% st_drop_geometry() %>% print(n=30)
# track_amt %>% filter(speed>50)
# track_amt %>% filter(locid %in% seq(21245,21258,1))
# 
# 
# ### FILTER OUT CRAZY SPEED LOCATIONS
# dim(track_amt) 
# crazyspeedlocs<- track_amt %>% filter(speed>50) %>% st_drop_geometry()
# 
# #l=21248
# #l=1067771
# #l=226069
# #l=6473339
# for (l in crazyspeedlocs$locid) {
# 	chunk<-track_amt %>% filter(locid %in% seq(l-5,l+5,1)) %>% st_drop_geometry() %>% select(-age_cy,-sex,-tod_) %>%
# 	mutate(prev_t=dplyr::lag(t_), post_t=dplyr::lead(t_), prev_dist=dplyr::lag(step_length),post_dist=dplyr::lead(step_length)) %>%
# 	mutate(pre_dt=as.numeric(difftime(t_,prev_t, units="sec")),post_dt=as.numeric(difftime(post_t,t_, units="sec"))) %>%
# 	rowwise() %>%
# 	mutate(dt=min(pre_dt,post_dt, na.rm=T)) %>%
# 	ungroup()
# 
# 	## filter the location that is (1) among the top 2 step_lengths, (2) among the top 2 turning angles, (3) min of (pre_dt and post_dt) among the bottom 2 dts
# 	sel1a<-slice_max(chunk,step_length,n=2)
# 	sel1b<-slice_max(chunk,speed,n=2)
# 	sel2<-slice_max(chunk,turning_angle,n=2)
# 	sel3<-slice_min(chunk,dt,n=2)
# 	sel3b<-slice_min(chunk,post_dt,n=1)
# 	badid<-Reduce(intersect, list(unique(sel1a$locid,sel1b$locid),sel2$locid,sel3$locid))
# 	if(length(badid)<1){badid<-Reduce(intersect, list(sel1b$locid,sel3b$locid))}
# 	if(length(badid)<1){badid<-Reduce(intersect, list(unique(sel1a$locid,sel1b$locid),sel2$locid))}
# 	if(length(badid)==1){track_amt<-track_amt %>% filter(locid!=badid)
# 	rm(sel1,sel2,sel3,chunk,badid)}
# }
# dim(track_amt)
# 
# ### RECALCULATE METRICS AFTER HAVING ELIMINATED CRAZY SPEED LOCS
# track_amt$locid<-seq_along(track_amt$t_)
# track_amt$step_length<-amt::step_lengths(track_amt)       # include step lengths
# track_amt$turning_angle<-abs(amt::direction_rel(track_amt,append_last=T, full_circle = FALSE,lonlat = FALSE))      # include RELATIVE turning angles - changed from absolute
# track_amt$turning_angle<-ifelse(is.na(track_amt$turning_angle),0,track_amt$turning_angle)  ## replace NA turning angles (caused by 0 distance) with 0
# track_amt$speed<-amt::speed(track_amt)      # include speed
# lastlocs<- track_amt %>% st_drop_geometry() %>%
# 	arrange(id,t_) %>%
# 	group_by(id) %>%
# 	summarise(last=max(t_), lastloc=max(locid))
# track_amt$speed[track_amt$locid %in% lastlocs$lastloc]<-NA
# track_amt$step_length[track_amt$locid %in% lastlocs$lastloc]<-NA
# track_amt$turning_angle[track_amt$locid %in% lastlocs$lastloc]<-NA
# dim(track_amt)
# 
# ## check whether the speeds are now better
# crazyspeedlocs<- track_amt %>% filter(speed>50) %>% st_drop_geometry()
# crazyspeedlocs




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



# COMBINING DATA WITH FORESTS AND BUILDINGS -----------------------------------------------------------------

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
## there are two duplicates that cause errors when using st_difference because they lead to identical buffers
# EXPFEEDERS<- fread("data/experimental_feeding.csv") %>%
#   mutate(long=as.numeric(lon)) %>%
#   filter(!is.na(long)) %>%
#   filter(!is.na(lat)) %>%
#   st_as_sf(coords = c("long", "lat"))
# ### troubleshoot the duplicates
# # EXPFEEDERS %>% st_transform(crs = 3035) %>% group_by(geometry) %>% summarise(N=length(nest_name),event=min(event_id),event2=max(event_id)) %>% filter(N>1)
# # EXPFEEDERS %>% filter(event_id %in% c(5610,5611,4413,4414))
# # EXPFEEDERS %>% filter(nest_name %in% c("Alterswil_Grabach"))
# st_crs(EXPFEEDERS) <- 4326
# 
# EXPFEEDERS_buff<- EXPFEEDERS %>%
#   #dplyr::filter(!event_id %in% c(5610,4413)) %>%  ## remove 2 offending stations with duplicate coordinates that cause problems when overlaying with locations
#   #dplyr::filter(!event_id %in% c(206,218)) %>%  ## remove 2 offending stations with duplicate coordinates that cause problems when overlaying with locations
#   st_transform(crs = 3035) %>%
#   st_buffer(dist=50) %>%
#   select(nest_name,Timepoint,food_placed) %>%
#   rename(feeder_id=nest_name)

# plot(EXPFEEDERS_buff[which(st_coordinates(EXPFEEDERS_buff)[,2]==2638298.634522798),])
# st_make_valid(EXPFEEDERS_buff[which(st_coordinates(EXPFEEDERS_buff)[,2]==2638298.634522798),])
# st_difference(st_make_valid(EXPFEEDERS_buff))

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
## this causes memory allocation error, so need to do it in a loop, which takes 45 min
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


# EXPORT THE ANNOTATED DATA FOR ALL OF SUI AND STUDY AREA -----------------------------------------------------------------



### join data with individual info
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

#### SIMPLE SUMMARY FOR MANUSCRIPT
track_out %>%  mutate(extra=year_id) %>%
  separate_wider_delim(extra, delim="_", names=c("year","bird_id")) %>%
  group_by(bird_id) %>%
  summarise(N= length(unique(year_id)))


##### CROP TO STUDY AREA ONLY ###########

track_sf <- track_sf %>% st_crop(x=track_sf, y=st_bbox(forest)) ## this reduces n locs from 10 mio to 8 mio

track_out <- track_sf %>% 
  st_transform(crs = 4326) %>%
  dplyr::mutate(long = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2]) %>%
  st_drop_geometry() %>%
  rename(year_id=id)

fwrite(as.data.frame(track_out),"data/REKI_annotated_feeding2025_study_area.csv")
saveRDS(track_sf, file = "data/REKI_trackingdata_annotated2025_study_area.rds")

head(track_sf)
dim(track_sf)


