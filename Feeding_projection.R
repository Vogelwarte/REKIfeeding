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
library(ranger)
library(caret)
library(randomForest)
library(sf)
library(lubridate)
library(leaflet)
library(adehabitatHR)
library(stars)
filter<-dplyr::filter
sf_use_s2(FALSE)
library(pROC)
library(amt)
library(recurse)

## set root folder for project
setwd("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/Feeding")



###### LOADING DATA -----------------------------------------------------------------

load("Feeding_analysis.RData")

#define different coordinate systems
CH_LV95_coords = "+init=epsg:2056"
EU_coords = "+init=epsg:3035"
WGS84_coords ="+init=epsg:4326"
CH_LV03_coords ="+init=epsg:21781"

### READ IN SHAPEFILES OF FOREST AND BUILDINGS
st_layers("data/buildings.gpkg")
buildings_CH <- st_read("data/buildings.gpkg", "buildings") %>%
  st_transform(crs = 3035) %>%
  rename(build_id=id)

st_layers("data/wald.gpkg")
forest_CH <- st_read("data/wald.gpkg", "wald") %>%
  st_transform(crs = 3035)



###### REGENERATE TRACKING DATA FOR ALL OF SWITZERLAND -----------------------------------------------------------------

### LOAD THE TRACKING DATA AND INDIVIDUAL SEASON SUMMARIES
trackingdata<-readRDS(file = "data/REKI_trackingdata_raw.rds")

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

### CALCULATE OTHER METRICS
track_amt$step_length<-amt::step_lengths(track_amt)       # include step lengths
track_amt$turning_angle<-amt::direction_abs(track_amt,append_last=T)      # include turning angles
track_amt$speed<-amt::speed(track_amt)      # include speed

head(track_amt)


# RECURSIONS FOR EACH LOCATION 50 m BUFFER--------------------------------------------------
# splitting track into a list with each single id grouped to an element
track_amt <- as.data.frame(track_amt)
track_amt_list <- split(track_amt, track_amt$id)


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

#trackingdata<-fread("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/NestTool/REKI/output/02_preprocessing/03_milvus_combined.csv")
indseasondata<-fread("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/NestTool/REKI/output/01_validation/03_validation_combined.csv")
nestdata<-fread("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/NestTool/REKI/data/Basic_nest_list_2015_2022.csv")

### NESTS
nests<- nestdata %>% dplyr::select(nest_name,tree_spec,latitude,longitude) %>%
  st_as_sf(coords = c("longitude", "latitude"))
st_crs(nests) <- 4326
nests<-nests %>% st_transform(crs = 3035)
nests_buff<- nests %>%
  st_buffer(dist=50)
nests_buff

### remove locations outside of Switzerland
dim(track_amt)
track_sf <- track_amt %>% 
  st_as_sf(coords = c("x_", "y_"))
st_crs(track_sf) <- 3035
head(track_sf)
track_sf <- track_sf %>% st_crop(x=track_sf, y=st_bbox(forest_CH))
dim(track_sf)


### spatial joins with forests, buildings and feeders
head(track_sf)
track_sf <- track_sf %>%
  st_join(forest_CH,
          join = st_intersects, 
          left = TRUE) %>%
  st_join(buildings_CH,
          join = st_intersects, 
          left = TRUE) %>%
  st_join(nests_buff,
          join = st_intersects, 
          left = TRUE) %>%
  mutate(BUILD=ifelse(is.na(build_id),0,1),
         NEST=ifelse(is.na(nest_name),0,1),
         FOREST=ifelse(is.na(objektart),0,1)) %>%   
  select(-objektart,-nest_name,-tree_spec,-build_id) %>%
  rename(building_size=area)

head(track_sf)

st_bbox(forest)
st_bbox(buildings)

### calculate distance to nearest KNOWN nest
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

saveRDS(track_sf, file = "data/REKI_trackingdata_annotated_CH.rds")
head(track_sf)


### LOAD THE TRACKING DATA AND INDIVIDUAL SEASON SUMMARIES
track_sf<-readRDS(file = "data/REKI_trackingdata_annotated_CH.rds")




###### APPLY MODELS ACROSS EXTENT OF TRACKING DATA  -----------------------------------------------------------------


################ CREATE VARIOUS SUBSETS TO IMPROVE RATIO OF CLASS MEMBERSHIP
track_day<-track_sf %>% dplyr::filter(tod_=="day")
track_nofor<-track_sf %>% filter(FOREST==0)
track_nofor_day<-track_sf %>% filter(FOREST==0) %>% filter(tod_=="day")

track_nofor_day_build<-track_nofor_day %>%
  filter(!(BUILD==0)) 


################ PREPARE DATA FOR PREDICTION

DATA <- track_nofor_day_build %>%
  mutate(YDAY=yday(t_), hour=hour(t_), month=month(t_)) %>%
  filter(!is.na(step_length)) %>%
  filter(!is.na(turning_angle)) %>%
  filter(!is.na(speed)) %>%
  filter(!is.na(mean_speed)) %>%
  filter(!is.na(mean_angle)) %>%
  select(-tod_,-FOREST) %>%
  mutate(point_id=seq_along(t_))
head(DATA)



# CREATING BASELINE MAP OF SAMPLING INTENSITY -----------------------------------------------------------------

grid_CH <- forest_CH %>% 
  st_make_grid(cellsize = 500, what = "polygons",
               square = FALSE) # This statements leads to hexagons
track_sf <-   st_as_sf(track_sf, coords = c("long", "lat"), crs = 4326) %>%
  st_transform(3035)
tab <- st_intersects(grid_CH, track_sf)
lengths(tab)
countgrid <- st_sf(n = lengths(tab), geometry = st_cast(grid_CH, "MULTIPOLYGON")) %>%
  st_transform(3035)
summary(log(countgrid$n+1))
pal <- colorNumeric(c("cornflowerblue","firebrick"), seq(0,10))
leaflet(options = leafletOptions(zoomControl = F)) %>% #changes position of zoom symbol
  htmlwidgets::onRender("function(el, x) {L.control.zoom({ 
                           position: 'bottomright' }).addTo(this)}"
  ) %>% #Esri.WorldTopoMap #Stamen.Terrain #OpenTopoMap #Esri.WorldImagery
  addProviderTiles("Esri.WorldImagery", group = "Satellite",
                   options = providerTileOptions(opacity = 0.6, attribution = F)) %>%
  addProviderTiles("OpenTopoMap", group = "Roadmap", options = providerTileOptions(attribution = F)) %>%  
  addLayersControl(baseGroups = c("Satellite", "Roadmap")) %>%  
  
  addPolygons(
    data=countgrid %>%
      st_transform(4326),
    stroke = TRUE, color = ~pal(log(n+1)), weight = 1,
    fillColor = ~pal(log(n+1)), fillOpacity = 0.5
  )



# ##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
# ########## PREDICTING INDIVIDUAL LOCATIONS THAT MATCH FEEDING LOCATION PATTERN   #############
# ##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################


#### apply RF2 model across all data

PRED<-stats::predict(RF2,data=DATA, type = "response")

DATA <- DATA %>%
  dplyr::bind_cols(PRED$predictions) %>%
  dplyr::rename(no_feed_prob = NO, feed_prob = YES) %>%
  dplyr::mutate(FEEDER_predicted = as.factor(dplyr::case_when(feed_prob > prevalence ~ "YES",
                                                              feed_prob < prevalence ~ "NO")))



##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
########## COUNT POINTS AND INDIVIDUALS IN GRID AND PREDICT FEEDING GRIDS   #############
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################

#### FIRST COUNT N INDIVIDUALS PER GRID CELL

for(c in 1:length(tab)){
  countgrid$N_ind[c]<-length(unique(track_sf$year_id[tab[c][[1]]]))
}

#### SECOND COUNT N INDIVIDUALS AND PREDICTED FEEDING LOCS PER GRID CELL
OUT_sf<-DATA %>%
  filter(FEEDER_predicted=="YES") %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
  st_transform(3035)
tab2 <- st_intersects(grid_CH, OUT_sf)
countgrid$N_feed_points<-lengths(tab2)
for(c in 1:length(tab2)){
  countgrid$N_feed_ind[c]<-length(unique(OUT_sf$year_id[tab2[c][[1]]]))
}



### MANIPULATE COUNTED ANIMALS INTO PROPORTIONS

countgrid<-countgrid %>%
  mutate(prop_feed=N_feed_ind/N_ind) %>%
  mutate(prop_pts=N_feed_points/n)

hist(countgrid$prop_feed)
hist(countgrid$prop_pts)



##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
########## PREDICT SECOND RANDOM FOREST MODEL TO PREDICT FEEDING SITE AT GRID CELL LEVEL   #############
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################


PRED_GRID<-countgrid %>% 
  mutate(gridid=seq_along(n)) %>%
  filter(n>10) %>%
  st_drop_geometry()

#### classification success of training data

PRED<-stats::predict(RF3,data=PRED_GRID, type = "response")
PRED_GRID <- PRED_GRID %>%
  dplyr::mutate(FEEDER_predicted=PRED$predictions[,2])
dim(PRED$predictions)
dim(PRED_GRID)



#### BASIC STATISTICS OF PREDICTIONS
dim(PRED_GRID %>% filter(FEEDER_observed==1  & n>10))
dim(PRED_GRID %>% filter(FEEDER_observed==1 & FEEDER_predicted>0.5))
dim(PRED_GRID %>% filter(FEEDER_observed==1 & n>10 & FEEDER_predicted>mean(PRED_GRID$FEEDER_observed)))/dim(PRED_GRID %>% filter(FEEDER_observed==1  & n>10))
hist(PRED_GRID$FEEDER_predicted)
mean(PRED_GRID$FEEDER_observed)


########## CREATE OUTPUT GRID WITH PREDICTED FEEDING LOCATIONS ########################

OUTgrid<-countgrid %>% select(-FEEDER) %>%
  mutate(gridid=seq_along(n)) %>%
  left_join(PRED_GRID, by=c("gridid","n","N_ind","N_feed_points","N_feed_ind","prop_feed","prop_pts")) %>%
  filter(!is.na(FEEDER_predicted))



##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
########## VALIDATE PREDICTIONS WITH SURVEY DATA FROM EVA CEREGHETTI  #############
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
## read in survey data from Eva Cereghetti
## Q1 is the question whether they feed or not
setwd("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/Feeding")
feed_surveys<-fread("data/survey.final.csv") %>% #filter(Q1=="Ja") %>%
  mutate(FEEDER_surveyed=ifelse(Q1=="Ja",1,0)) %>%
  select(nr,coord_x,coord_y,square,random,area,building.type,FEEDER_surveyed) %>%
  st_as_sf(coords = c("coord_x", "coord_y"), crs=21781) %>%
  st_transform(crs = 3035) 

validat <- st_intersection(OUTgrid, feed_surveys)
AUC_TEST<-auc(roc(data=validat,response=FEEDER_surveyed,predictor=FEEDER_predicted))

### summarise the predicted sites
validat<-validat %>% filter(FEEDER_surveyed==1)
hist(validat$FEEDER_predicted)
mean(validat$FEEDER_predicted)
length(validat[validat$FEEDER_predicted>0.5,])/dim(validat)[1]






##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
########## PLOT THE MAP FOR KNOWN AND OBSERVED FEEDING STATIONS   #############
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################

########## CREATE A LEAFLET MAP OF PREDICTED FEEDING LOCATIONS ########################

## need to specify color palette 
# If you want to set your own colors manually:
pred.pal <- colorNumeric(c("cornflowerblue","firebrick"), seq(0,1))
pal <- colorNumeric(c("cornflowerblue","firebrick"), seq(0,15))
#year.pal <- colorFactor(topo.colors(7), KDE_sf$id)
feed.pal <- colorFactor(c("darkgreen","lightgreen"), unique(plot_feeders$Type))
m2 <- leaflet(options = leafletOptions(zoomControl = F)) %>% #changes position of zoom symbol
  setView(lng = mean(st_coordinates(plot_feeders)[,1]), lat = mean(st_coordinates(plot_feeders)[,2]), zoom = 11) %>%
  htmlwidgets::onRender("function(el, x) {L.control.zoom({ 
                           position: 'bottomright' }).addTo(this)}"
  ) %>% #Esri.WorldTopoMap #Stamen.Terrain #OpenTopoMap #Esri.WorldImagery
  addProviderTiles("Esri.WorldImagery", group = "Satellite",
                   options = providerTileOptions(opacity = 0.6, attribution = F,minZoom = 5, maxZoom = 20)) %>%
  addProviderTiles("OpenTopoMap", group = "Roadmap", options = providerTileOptions(attribution = F,minZoom = 5, maxZoom = 15)) %>%  
  addLayersControl(baseGroups = c("Satellite", "Roadmap")) %>%  
  
  addCircleMarkers(
    data=plot_OUT,
    radius = 2,
    stroke = TRUE, color = "black", weight = 0.5,
    fillColor = "grey75", fillOpacity = 0.5,
    popup = ~ paste0("year ID: ", plot_OUT$year_id, "<br>", plot_OUT$timestamp)
  ) %>%
  addPolygons(
    data=OUTgrid %>%
      st_transform(4326),
    stroke = TRUE, color = ~pred.pal(FEEDER_predicted), weight = 1,
    fillColor = ~pred.pal(FEEDER_predicted), fillOpacity = 0.5
  ) %>%
  
  addCircleMarkers(
    data=plot_feeders,
    radius = 3,
    stroke = TRUE, color = ~feed.pal(Type), weight = 1,
    fillColor = ~feed.pal(Type), fillOpacity = 0.2
  ) %>%
  
  addCircleMarkers(
    data=validat %>%
      st_transform(4326),
    radius = 5,
    stroke = TRUE, color = "red", weight = 1,
    fillColor = "red", fillOpacity = 0.5
  ) %>%
  

  addLegend(     # legend for predicted prob of feeding
    position = "topleft",
    pal = pred.pal,
    values = OUTgrid$FEEDER_predicted,
    opacity = 1,
    title = "Predicted probability of </br>anthropogenic feeding"
  ) %>%
  addLegend(     # legend for known feeding sites
    position = "topleft",
    pal = feed.pal,
    values = plot_feeders$Type,
    opacity = 1,
    title = "Type of feeding station"
  ) %>%
  
  
  
  addScaleBar(position = "bottomright", options = scaleBarOptions(imperial = F))

m2


htmltools::save_html(html = m2, file = "C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/Feeding/output/Potential_feeding_grids.html")
mapview::mapshot(m2, url = "C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/Feeding/output/Potential_feeding_grids.html")
st_write(OUTgrid,"output/REKI_predicted_anthropogenic_feeding_areas.kml",append=FALSE)




##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
########## COMPILE THE OUTPUT REPORT   #############
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
library(knitr)
library(markdown)
library(rmarkdown)

### create HTML report for overall summary report
# Sys.setenv(RSTUDIO_PANDOC="C:/Users/Inge Oppel/AppData/Local/Pandoc")
# Sys.setenv(RSTUDIO_PANDOC="C:/Program Files/RStudio/bin/pandoc")
#
rmarkdown::render('C:\\Users\\sop\\OneDrive - Vogelwarte\\REKI\\Analysis\\Feeding\\Feeding_Analysis_report2023_v2.Rmd',
                  output_file = "Feeding_Analysis_Progress2023.html",
                  output_dir = 'C:\\Users\\sop\\OneDrive - Vogelwarte\\REKI\\Analysis\\Feeding\\output')

