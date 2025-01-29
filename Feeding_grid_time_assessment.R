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
library(adehabitatHR)
library(stars)
library(trip)
library(sqldf)
library(amt)
library(readxl)
library(adehabitatLT)
filter<-dplyr::filter
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






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# JOIN TRACKING DATA WITH PREDICTION GRID AND SUMMARISE HOURS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm('all_dat15min')
gc()

FEED_TRACK<-track_sf %>%
  #st_crop(SUI) %>%  ## remove all locations outside - should not be needed because raw data are already 
  st_join(CHgrid) %>%
  mutate(FEEDER_predicted=ifelse(is.na(FEEDER_predicted),0,FEEDER_predicted)) %>%
  separate(season_id, into=c("year","bird_id"), sep="_")

# rm('track_sf')
# gc()
# 
# SUMMARY<-FEED_TRACK  %>%
#   st_drop_geometry() %>%
#   #mutate(season=ifelse(substr(season_id,1,1)=="B","breeding","non-breeding")) %>%
#   mutate(season_id = ifelse(yday(date)<70 ,
#                             paste(season,year(date)-1,bird_id,sep="_"),
#                             paste(season,year(date),bird_id,sep="_"))) %>%
#   group_by(season_id, season) %>%
#   summarise(FEED=mean(FEEDER_predicted))
# 
# 
# ggplot(SUMMARY, aes(x=season,y=FEED)) + geom_boxplot()





###### LOADING INDIVIDUAL DATA THAT PREDICTED FEEDING PROBABILITY ------------------------------------------------------------



inddat<-read_excel("C:/Users/sop/OneDrive - Vogelwarte/General/DATA/Individual_life_history_2015-2023.xlsx", sheet="Individual_life_history_2015-20") %>% # updated on 3 June 2024 to include birds from 2022
  dplyr::select(bird_id,ring_number,tag_year,sex_compiled, age, hatch_year) %>%
  rename(ring_id=ring_number) %>%
  rename(sex=sex_compiled) %>%
  mutate(hatch_year=if_else(is.na(as.numeric(hatch_year)),tag_year-3,as.numeric(hatch_year)))


### MERGE tracks with INDIVIDUAL DATA

FEED_DATA <- FEED_TRACK %>%
  select(-bird_id) %>%
  rename(bird_id=year, year=season) %>%
  mutate(bird_id = as.integer(bird_id), year=as.integer(year)) %>%
  mutate(season = ifelse(month(date) %in% c(3,4,5,6,7,8) ,"B","NB")) %>%   ## 1 SEPT AS CUT OFF
  st_transform(4326) %>%
  dplyr::mutate(long = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2]) %>%
  st_drop_geometry() %>%
  left_join(inddat, by='bird_id') %>%
  select(-n, -N_ind, -N_feed_points, -N_feed_ind, -prop_feed, -prop_pts, -gridid) %>%
  mutate(season_id = paste(season,year,bird_id,sep="_")) %>%
  #mutate(year=year(date)) %>%
  mutate(age_cy=(year-tag_year)+1)
  # mutate(HR=ifelse(home_range_id>0,"settled","not settled")) %>%
  # mutate(BR=ifelse(nest_id>0,"breeding","not breeding"))
  
saveRDS(FEED_DATA,"output/REKI_food_supplementation_index2025.rds")
