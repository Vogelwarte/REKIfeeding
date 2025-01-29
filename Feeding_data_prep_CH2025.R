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
## reverted to model building in study area - did not work across all of Switzerland

## ADDED NEW FEEDERS REPORTED AT MAT on 25/26 Jan 2025


### KILLED PROCESS AFTER 24 HRS on 29 JAN 2025 - cannot figure out why spatial overlay takes so long....

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
track_amt<-readRDS("data/track_amt2025.rds")

### LOAD INDIVIDUAL LIFE HISTORIES
#trackingdata<-fread("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/NestTool/REKI/output/02_preprocessing/03_milvus_combined.csv")
#indseasondata<-fread("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/NestTool/REKI/output/01_validation/03_validation_combined.csv")
indseasondata<-read_excel("C:/Users/sop/OneDrive - Vogelwarte/General/DATA/Individual_life_history_2015-2023.xlsx", sheet="Individual_life_history_2015-20") %>% # updated on 3 June 2024 to include birds from 2022
  dplyr::select(bird_id,ring_number,tag_year,sex_compiled, age, hatch_year) %>%
  rename(ring_id=ring_number) %>%
  mutate(hatch_year=if_else(is.na(as.numeric(hatch_year)),tag_year-3,as.numeric(hatch_year)))
nestdata<-fread("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/NestTool/REKI/data/Basic_nest_list_2015_2022.csv")


# COMBINING DATA WITH FORESTS AND BUILDINGS -----------------------------------------------------------------

#define different coordinate systems
CH_LV95_coords = "+init=epsg:2056"
EU_coords = "+init=epsg:3035"
WGS84_coords ="+init=epsg:4326"
CH_LV03_coords ="+init=epsg:21781"



### READ IN SHAPEFILES OF FOREST AND BUILDINGS
st_layers("data/buildings.gpkg")
buildings<- st_read("data/buildings.gpkg", "buildings") %>%
  st_transform(crs = 3035) %>%
  rename(build_id=id)

st_layers("data/wald.gpkg")
forest<- st_read("data/wald.gpkg", "wald") %>%
  st_transform(crs = 3035)




### READ IN FEEDERS
FEEDERS<- fread("data/Private_Feeders/private_feeders_upd2022.csv") %>% 
  filter(!is.na(coordX)) %>%
  st_as_sf(coords = c("coordX", "coordY"))

FEEDERS<-fread("data/Private_Feeders/reki_feeding_stations_MATadditions2025.csv") %>% 
  st_as_sf(coords = c("long", "lat")) %>%
  bind_rows(FEEDERS)
  


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
# track_sf <- track_sf %>% st_crop(x=track_sf, y=st_bbox(forest)) ## this reduces n locs from 8 mio to 5.5 mio
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


