###########################################################################################################################
###### ANALYSIS TO QUANTIFY USAGE OF ANTHROPOGENIC FEEDING SITES OF RED KITES ACROSS SWITZERLAND ################
###########################################################################################################################
# original idea by Nathalie Heiniger (MSc thesis 2020)
# uses output provided by "Feeding_projection.r"
# created by steffen.oppel@vogelwarte.ch in November 2023
# inspired by REKI Team meeting on 14 Nov 2023 - include metrics of usage in manuscript

## updated on 2 Oct 2024 to align year from 1 Sept to 31 Aug

## completely revised in January 2025 to use more recent tracking data
rm(list=ls())
library(tidyverse)
library(dplyr, warn.conflicts = FALSE)
options(dplyr.summarise.inform = FALSE)
library(dtplyr)
library(sf)
library(lubridate)
library(stars)
library(sqldf)
library(readxl)
filter<-dplyr::filter
select<-dplyr::select
sf_use_s2(FALSE)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# READ IN REGULARISED TRACKING DATA 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### script outsourced to "REKI_track_interpolation.r"
setwd("C:/Users/sop/OneDrive - Vogelwarte/General/DATA")
track_sf<-readRDS("REKI_regular_15min_tracking_data_January2025_projected.rds")

## set root folder for project
setwd("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/REKIfeeding")
#setwd("C:/STEFFEN/OneDrive - Vogelwarte/REKI/Analysis/REKIfeeding")

### read in Switzerland map
# SUI<-st_read("S:/rasters/outline_maps/swiss_map_overview/layers.gpkg") %>% filter(country=="Switzerland")
# saveRDS(SUI,"data/Swiss_border.rds")
SUI<-readRDS("data/Swiss_border.rds") %>%
  st_transform(3035)
plot(SUI)


###### LOADING DATA THAT PREDICTED FEEDING PROBABILITY ------------------------------------------------------------
CHgrid<-readRDS("output/REKI_feeding_grid2024.rds")
range(CHgrid$FEEDER_predicted)
dim(track_sf)
plot(CHgrid)

### READ IN FEEDERS
feeders<-fread("data/Private_Feeders/private_feeders_upd2022.csv") %>% 
  filter(!is.na(coordX)) %>%
  st_as_sf(coords = c("coordX", "coordY"), crs=21781) %>%
  st_transform(crs = 4326) %>%
  mutate(Type="Private") %>%
  select(Type)

feeders<-fread("data/Private_Feeders/reki_feeding_stations_MATadditions2025.csv") %>% 
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
  mutate(Type="Private") %>%
  select(Type) %>%
  bind_rows(feeders)

### EXPERIMENTAL FEEDING STATIONS
feeders2<- fread("data/experimental_feeding.csv") %>%
  filter(!is.na(lon)) %>%
  filter(!is.na(lat)) %>%
  st_as_sf(coords = c("lon", "lat"), crs=4326)%>%
  mutate(Type="Experimental") %>%
  select(Type)

### FEEDING PLATFORMS - are also experimental
feeders3<- fread("data/feeding_platforms_15_16.csv") %>%
  filter(!is.na(x)) %>%
  mutate(start_date=dmy(start_date), end_date=dmy(end_date)) %>%
  st_as_sf(coords = c("x", "y"), crs=21781) %>%
  st_transform(crs = 4326) %>%
  mutate(Type="Platforms") %>%
  select(Type)

### ADD SURVEY DATA FROM EVA CEREGHETTI AND FIONA PELLET  #############
## read in survey data from Eva Cereghetti
feeders4<- fread("C:/Users/sop/OneDrive - Vogelwarte/General/MANUSCRIPTS/AnthropFeeding/DataArchive/REKI_validation_feeders.csv") %>%
  select(-geometry) %>%
  filter(FEEDER_surveyed==1) %>%
  st_as_sf(coords = c("long", "lat"), crs=4326) %>%
  mutate(Type="Private") %>%
  select(Type)

feeders<-rbind(feeders,feeders2,feeders3, feeders4) %>%
  st_transform(3035)


###### LOADING INDIVIDUAL DATA THAT PREDICTED FEEDING PROBABILITY ------------------------------------------------------------

inddat<-read_excel("C:/Users/sop/OneDrive - Vogelwarte/General/DATA/Individual_life_history_2015-2023.xlsx", sheet="Individual_life_history_2015-20") %>% # updated on 3 June 2024 to include birds from 2022
  dplyr::select(bird_id,ring_number,tag_year,sex_compiled, age, hatch_year) %>%
  rename(ring_id=ring_number,sex=sex_compiled) %>%
  mutate(bird_id=as.character(bird_id)) %>%
  mutate(hatch_year=if_else(is.na(as.numeric(hatch_year)),tag_year-3,as.numeric(hatch_year)))

rm(feeders2,feeders3,feeders4,SUI)
gc()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# JOIN FEEDERS WITH GRID TO UPDATE PRED_FEED TO 1 where we know that there are feeders
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

CHgrid<-CHgrid %>%
  sf::st_join(.,feeders) %>%
  mutate(Type=ifelse(is.na(Type),"None",Type)) %>%
  mutate(FEEDER_predicted=ifelse(Type=="Private",0.8,FEEDER_predicted))

unique(CHgrid$Type)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# JOIN TRACKING DATA WITH PREDICTION GRID AND SUMMARISE HOURS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
allids<-unique(track_sf$season_id)
length(allids)

## check what has been done already
done<-list.files("./output/ind")
doneids<-as.character()
for(f in 1:length(done)){
	doneids[f]<-substr(done[f],17,nchar(done[f])-4)
}
length(doneids)

## complete the remaining ids
allids<-allids[!(allids %in% doneids)]
length(allids)

for(i in allids) {
FEED_TRACK<-track_sf %>%
  filter(season_id==i) %>%
  #st_crop(SUI) %>%  ## remove all locations outside - should not be needed because raw data are already cropped to Switzerland
  sf::st_join(.,CHgrid) %>%
  mutate(FEEDER_predicted=ifelse(is.na(FEEDER_predicted),0,FEEDER_predicted)) %>%
  mutate(FEEDER_predicted=ifelse(Type=="Experimental" & startsWith(season_id,"B_2020"),0.5,FEEDER_predicted)) %>%
  mutate(FEEDER_predicted=ifelse(Type=="Platforms" & startsWith(season_id,"B_2016"),0.5,FEEDER_predicted)) %>%
  separate(season_id, into=c("season","year","bird_id"), sep="_") %>%
  #mutate(season = ifelse(month(date) %in% c(3,4,5,6,7,8) ,"B","NB")) %>%   ## 1 SEPT AS CUT OFF
  st_transform(4326) %>%
  dplyr::mutate(long = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2]) %>%
  st_drop_geometry() %>%
  left_join(inddat, by='bird_id') %>%
  select(-n, -N_ind, -N_feed_points, -N_feed_ind, -prop_feed, -prop_pts, -gridid) %>%
  mutate(season_id = paste(season,year,bird_id,sep="_")) %>%
  mutate(age_cy=(year(date)-tag_year)+1)
saveRDS(FEED_TRACK,sprintf("output/ind/REKI_feed_index_%s.rds",i))
}




### FIX WRONG YEAR IN NB SEASON FILES
# year is accidentally year+1 and not year-1 for all NB season files
# read in - relabel, and re-save
# completed once and then moved files into 'old' folder
# 
# done<-list.files("./output/ind")
# wrong<-done[which(startsWith(done,"REKI_feed_index_NB"))]
# for(f in 1:length(wrong)){
#   X <- readRDS(paste0("./output/ind/",wrong[f])) %>%
#     separate(season_id, into=c("season","year","bird_id"), sep="_") %>%
#     mutate(year=as.numeric(year), bird_id=as.numeric(bird_id)) %>%
#     mutate(year=year-1) %>%
#     mutate(season_id = paste(season,year,bird_id,sep="_"))
#   new_seas<-unique(X$season_id)
#   saveRDS(X,sprintf("output/ind/new/REKI_feed_index_%s.rds",new_seas))
# }
# 




### MERGE ALL DATA INTO SINGLE FILE OF SEASONAL PROPORTIONAL TIME SPENT FEEDING

FEED_DATA <- data.frame()
done<-list.files("./output/ind")
done<-done[-1]  ## remove the folder for old wrong files
for(f in 2544:length(done)){
 FEED_DATA <- readRDS(paste0("./output/ind/",done[f])) %>%
  mutate(hrs=0.25) %>%
  mutate(FEED_hrs=FEEDER_predicted*hrs) %>%
  group_by(season_id) %>%
  summarise(FEEDsum=sum(FEED_hrs, na.rm=T),trackhrs=sum(hrs)) %>%
  mutate(FEED=FEEDsum/trackhrs) %>%
  bind_rows(FEED_DATA)
saveRDS(FEED_DATA,"output/REKI_seasonal_feed_index2025.rds")
}
  

