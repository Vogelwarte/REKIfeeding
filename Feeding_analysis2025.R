###########################################################################################################################
###### DATA ANALYSIS TO IDENTIFY ANTHROPOGENIC FEEDING SITES OF RED KITES BASED ON GPS TRACKING DATA ################
###########################################################################################################################
# original idea by Nathalie Heiniger (MSc thesis 2020)
# uses csv table prepared by "Feeding_data_prep.r"
# created by steffen.oppel@vogelwarte.ch in June 2023
# 15 Aug: updated graphs to make nice ppt presentation

## 11 December 2023: included new evaluation method - Boyce index https://www.r-bloggers.com/2022/05/model-evaluation-with-presence-points-and-raster-predictions/
## 22 Jan 2024: finalised revisions of comments by Martin Grüebler

## RE-RUN ANALYSIS WITH 2024 data at very high resolution

## 11 June 2024: changed the count of ID (from year_id to bird_id)

## 30 June 2024: curtailed whole analysis to study area

## 16 JANUARY 2025 - re-ran analysis with all data up to 2025 (including UHF data)

#library(stringi)
rm(list=ls())
library(tidyverse)
library(dplyr, warn.conflicts = FALSE)
options(dplyr.summarise.inform = FALSE)
library(ranger)
library(caret)
library(readxl)
library(randomForest)
library(sf)
library(lubridate)
library(leaflet)
library(adehabitatHR)
library(stars)
library(tmap)
filter<-dplyr::filter
select<-dplyr::select
sf_use_s2(FALSE)
library(pROC)
library(modEvA)  ### for calculating the Boyce index
library(pdp) ### for creating partial dependence plots

## set root folder for project
setwd("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/REKIfeeding")
#setwd("C:/STEFFEN/OneDrive - Vogelwarte/REKI/Analysis/REKIfeeding")

# LOADING DATA -----------------------------------------------------------------
### LOAD THE TRACKING DATA AND INDIVIDUAL SEASON SUMMARIES
#track_sf<-readRDS(file = "data/REKI_trackingdata_annotated2025_CH.rds")
track_sf<-fread("data/REKI_annotated_feeding2025_CH.csv") %>%
  mutate(extra=year_id) %>%
  separate_wider_delim(extra, delim="_", names=c("year","bird_id"))

#track_sf<-load(file = "data/REKI_trackingdata_annotated2024_sf.rds")

#define different coordinate systems
CH_LV95_coords = "+init=epsg:2056"
EU_coords = "+init=epsg:3035"
WGS84_coords ="+init=epsg:4326"
CH_LV03_coords ="+init=epsg:21781"


### read in Switzerland map
# SUI<-st_read("S:/rasters/outline_maps/swiss_map_overview/layers.gpkg") %>% filter(country=="Switzerland")
# saveRDS(SUI,"data/Swiss_border.rds")
SUI<-readRDS("data/Swiss_border.rds")
plot(SUI)

# STUDY_AREA<-st_read("C:/STEFFEN/OneDrive - Vogelwarte/General/DATA/REKI_study_area_sm.kml")
# saveRDS(STUDY_AREA,"data/REKI_study_area.rds")
STUDY_AREA<-readRDS("data/REKI_study_area.rds")
plot(STUDY_AREA)
st_area(STUDY_AREA)/1000000

 ##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
 ########## CALCULATE DISTANCE FROM PREDICTED FEEDING SITE TO NEAREST KNOWN SITE   #############
 ##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
 
 # READ IN FEEDING LOCATIONS
 plot_feeders<-fread("data/Private_Feeders/private_feeders_upd2022.csv") %>% 
   filter(!is.na(coordX)) %>%
   st_as_sf(coords = c("coordX", "coordY"), crs=21781) %>%
   st_transform(crs = 4326) %>%
   mutate(Type="Private") %>%
   select(Type)

plot_feeders<-fread("data/Private_Feeders/reki_feeding_stations_MATadditions2025.csv") %>% 
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
  mutate(Type="Private") %>%
  select(Type) %>%
  bind_rows(plot_feeders)
 
 ### EXPERIMENTAL FEEDING STATIONS
 plot_feeders2<- fread("data/experimental_feeding.csv") %>%
   filter(!is.na(lon)) %>%
   filter(!is.na(lat)) %>%
   st_as_sf(coords = c("lon", "lat"), crs=4326)%>%
   mutate(Type="Experimental") %>%
   select(Type)
 
 ### FEEDING PLATFORMS - are also experimental
 plot_feeders3<- fread("data/feeding_platforms_15_16.csv") %>%
   filter(!is.na(x)) %>%
   mutate(start_date=dmy(start_date), end_date=dmy(end_date)) %>%
   st_as_sf(coords = c("x", "y"), crs=21781) %>%
   st_transform(crs = 4326) %>%
   mutate(Type="Experimental") %>%
   select(Type)
 
 ### ADD SURVEY DATA FROM EVA CEREGHETTI AND FIONA PELLET  #############
 ## read in survey data from Eva Cereghetti
 plot_feeders4<- fread("C:/Users/sop/OneDrive - Vogelwarte/General/MANUSCRIPTS/AnthropFeeding/DataArchive/REKI_validation_feeders.csv") %>%
   select(-geometry) %>%
   filter(FEEDER_surveyed==1) %>%
   st_as_sf(coords = c("long", "lat"), crs=4326) %>%
   mutate(Type="Private") %>%
   select(Type)
 
 plot_feeders<-rbind(plot_feeders,plot_feeders2,plot_feeders3, plot_feeders4) %>% dplyr::mutate(long = sf::st_coordinates(.)[,1],
                                                                                                        lat = sf::st_coordinates(.)[,2])
 
# #fwrite(plot_feeders,"C:/Users/sop/MAT/REKI_feeding_stations.csv")
# 
# 
# ########## PLOT FEEDERS AND CHECK EXTENT OF STUDY AREA   #############
# 
# tmap_mode("view")
# tm_basemap(server="OpenStreetMap") +
#   tm_shape(SUI)  +
#   tm_polygons(col = "firebrick", fill="grey95", alpha=0.1) +
#   tm_shape(STUDY_AREA)  +
#   tm_polygons(col = "firebrick", fill="firebrick", alpha=0.2) +
#   tm_shape(plot_feeders)+
#   tm_symbols(col = 'green', size = 0.0025)
# 
# 
# 
# 
# ### READ IN SHAPEFILES OF FOREST AND BUILDINGS
# buildings <- st_read("data/Buildings/tlm_buildings_studyareaExtra_size65_buff_50m_sf_singlepoly.shp", stringsAsFactors=FALSE) %>%
#   st_transform(crs = 3035) 
# forest <- st_read("data/Forest/vec25_forest_buff_20m_studyarea.shp", stringsAsFactors=FALSE) %>%
#   st_transform(crs = 3035) %>%
#   select(AREA)
# 


################ CREATE VARIOUS SUBSETS TO IMPROVE RATIO OF CLASS MEMBERSHIP
# track_day<-track_sf %>% dplyr::filter(tod_=="day") %>%
#   st_intersection(,SUI) %>% filter(!is.na(admin))
# track_nofor<-track_sf %>% filter(FOREST==0)
track_nofor_day_build<-track_sf %>%
  filter(FOREST==0) %>%
  filter(tod_=="day") %>%
  filter(!(BUILD==0 & FEEDER =="NO")) %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
  st_transform(2056) %>%
  st_intersection(.,SUI) %>% filter(!is.na(country)) %>% ### remove all data outside of Switzerland
  st_transform(4326) %>%
  #st_intersection(.,STUDY_AREA) %>% filter(!is.na(Name)) %>% ### remove all data outside of study area   --------- BIG DECISION WHETHER TO BUILD MODEL FOR STUDY AREA ONLY OR ALL OF SUI
  dplyr::mutate(long = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2]) %>%
  st_drop_geometry()
 

table(track_nofor_day_build$FEEDER)


################ PREPARE DATA FOR PREDICTION

DATA <- track_nofor_day_build %>%
  mutate(YDAY=yday(t_), hour=hour(t_), month=month(t_)) %>%
  # mutate(extra=year_id) %>%
  # separate_wider_delim(extra, delim="_", names=c("year","bird_id")) %>%
  filter(!is.na(step_length)) %>%
  filter(!is.na(turning_angle)) %>%
  filter(!is.na(speed)) %>%
  filter(!is.na(mean_speed)) %>%
  filter(!is.na(mean_angle)) %>%
  select(-tod_,-FOREST,-forest_size,-build_id,-type_of_food,-frequency) %>%
  mutate(point_id=seq_along(t_)) %>%
  mutate(year=as.numeric(year), bird_id=as.factor(bird_id)) %>%
  filter(!is.na(age_cy)) ##### TEMPORARY APPROACH UNTIL I CAN GET DATA FROM 2022 onwards
head(DATA)



# CREATING BASELINE MAP OF SAMPLING INTENSITY -----------------------------------------------------------------

grid <- SUI %>% st_transform(3035) %>%
  st_make_grid(cellsize = 500, what = "polygons",
               square = FALSE) %>% # This statements leads to hexagons
  st_transform(2056) %>%
  st_intersection(.,SUI) %>% ### remove all data outside of Switzerland
  st_transform(3035)

sum(st_area(grid))/1000000  ## size of study area in sq km
track_sf <-   st_as_sf(track_sf, coords = c("long", "lat"), crs = 4326) %>%
  st_transform(3035)
tab <- st_intersects(grid, track_sf)
lengths(tab)
countgrid <- st_sf(n = lengths(tab), geometry = st_cast(grid, "MULTIPOLYGON")) %>%
  st_transform(3035)
summary(log(countgrid$n+1))

#### COUNT N INDIVIDUALS PER GRID CELL
## this is tab and not tab 2 because it is independent of behaviour
for(c in 1:length(tab)){
  #countgrid$N_ind[c]<-length(unique(track_sf$year_id[tab[c][[1]]]))
  countgrid$N_ind[c]<-length(unique(track_sf$bird_id[tab[c][[1]]]))
}




# view the map
tmap_mode("view")
tm_basemap(server="OpenStreetMap") +
  #tm_shape(countgrid)  +
  #tm_polygons(col = 'N_ind', fill='N_ind', alpha=0.4) +
  tm_shape(track_nofor_day_build)  +
  tm_dots(col = 'FEEDER')




##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
########## REDUCE WORKSPACE TO MAKE SAVING OUTPUT EASIER  #############
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
ls()
rm(forest,buildings,track_sf,track_nofor_day_build,plot_feeders2,plot_feeders3,plot_feeders4,tab,grid)
gc()



##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
########## ANALYSING DATA IN RANDOM FOREST  #############
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
### completely random subsetting for model evaluation - removed to use maximum of information for country-wide projection

##### RUN MODEL ##########
## takes about 45 minutes with 2024 dataset
#DATA<-DATA %>% 
  #dplyr::mutate(long = sf::st_coordinates(.)[,1],
  #              lat = sf::st_coordinates(.)[,2]) %>%
  #st_drop_geometry() %>%
#  dplyr::select(-Description)
RF2 <- ranger::ranger(as.factor(FEEDER) ~ sex + age_cy + YDAY + hour + month +                            ## basic variables such as age, sex, and time
                        revisits + residence_time + meanFreqVisit + n_days + TimeSpan + TempEven +        ## several variables dealing with the temporal revisitation pattern
                        step_length + turning_angle + speed +                                             ## several variables dealing with the current movement characteristics
                        mean_speed + mean_angle +                                                         ## variables dealing with the average movement characteristics
                        BUILD + NEST + dist_nest,     #                                                      ## variables dealing with distance to structures and nests
                      data=DATA, mtry=2, num.trees=2500, replace=F, importance="none", oob.error=T, write.forest=T, probability=T)
saveRDS(RF2, "output/feed_site_RF_model2025_CH.rds")

#### classification success of training data

PRED<-stats::predict(RF2,data=DATA, type = "response")
prevalence<-table(DATA$FEEDER)[2]/sum(table(DATA$FEEDER))
hist(PRED$predictions)
OUT <- DATA %>%
  dplyr::mutate(FEEDER_observed = factor(FEEDER,levels=c("NO","YES"))) %>%
  dplyr::bind_cols(PRED$predictions) %>%
  dplyr::rename(no_feed_prob = NO, feed_prob = YES) %>%
  dplyr::mutate(feed_obs_num = ifelse(FEEDER=="NO",0,1))
  

suppressWarnings({testmat<-caret::confusionMatrix(data = OUT$FEEDER_predicted, reference = OUT$FEEDER_observed, positive="YES")})


### DEFINE PREDICTION THRESHOLD FOR FEEDING LOCATIONS AND CONVERT TO BINARY OUTCOME###
## we cannot predict correct ABSENCE of feeding locations - even if one household does not feed their neighbours may and the prediction is therefore useless (and falsifying accuracy)
## but we use ROC curve to identify threshold
ROC_val<-pROC::roc(data=OUT,
                   response=feed_obs_num,
                   predictor=feed_prob)
AUC_VAL<-pROC::auc(ROC_val)
THRESH_pts<-pROC::coords(ROC_val, "best", "threshold")$threshold
OUT <- OUT %>%
  dplyr::mutate(FEEDER_predicted = ifelse(feed_prob < THRESH_pts, "NO","YES")) 


### CALCULATE BOYCE INDEX FOR VALIDATION DATA ###
BI<-Boyce(obs = OUT$feed_obs_num, pred = OUT$feed_prob)$Boyce
BI






##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
########## SECOND LEVEL PREDICTION: COUNT POINTS AND INDIVIDUALS IN GRID   #############
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################

#### ASSESS NUMBER OF AVAILABLE LOCATIONS
table(OUT$FEEDER_predicted)
OUT %>% dplyr::filter(FEEDER_predicted=="YES") %>% summarise(n=length(unique(year_id)))

#### CREATE A COUNT GRID SIMPLE FEATURE
OUT_sf<-OUT %>%
  filter(FEEDER_predicted=="YES") %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
  st_transform(3035)
tab2 <- st_intersects(countgrid, OUT_sf)
countgrid$N_feed_points<-lengths(tab2)


#### SECOND COUNT N INDIVIDUALS AND PREDICTED FEEDING LOCS PER GRID CELL

for(c in 1:length(tab2)){
  #countgrid$N_feed_ind[c]<-length(unique(OUT_sf$year_id[tab2[c][[1]]]))
  countgrid$N_feed_ind[c]<-length(unique(OUT_sf$bird_id[tab2[c][[1]]]))
}



### MANIPULATE COUNTED ANIMALS INTO PROPORTIONS

countgrid<-countgrid %>%
  mutate(prop_feed=N_feed_ind/N_ind) %>%
  mutate(prop_pts=N_feed_points/n)

hist(countgrid$prop_feed)
hist(countgrid$prop_pts)
summary(countgrid)

countgrid %>% filter(is.na(prop_feed))




##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
########## SUMMARISE AND PLOT GRID CELLS WITH FEEDERS   #############
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################

### overlay feeder data with countgrid
feed_grd<-st_intersects(countgrid,(plot_feeders %>% st_transform(3035)))
countgrid$FEEDER<-lengths(feed_grd) 


PRED_GRID<-countgrid %>% 
  mutate(gridid=seq_along(n)) %>%
  filter(n>10) %>%
  st_transform(4326) %>%
  #st_intersection(.,STUDY_AREA) %>% filter(!is.na(Name)) %>% ### remove all data outside of study area   --------- BIG DECISION WHETHER TO BUILD MODEL FOR STUDY AREA ONLY OR ALL OF SUI
  st_drop_geometry() %>%
  mutate(FEEDER=ifelse(FEEDER==0,0,1))
dim(PRED_GRID)

##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
########## FIT SECOND RANDOM FOREST MODEL TO PREDICT FEEDING SITE AT GRID CELL LEVEL   #############
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################

RF3 <- ranger::ranger(FEEDER ~ n+N_ind+N_feed_points+N_feed_ind+prop_feed+prop_pts,     
                      data=PRED_GRID, mtry=2, num.trees=2500, replace=F, importance="none", oob.error=T, write.forest=T, probability=T)
saveRDS(RF3, "output/feed_grid_RF_model2025_CH.rds")


#### classification success of training data

PRED<-stats::predict(RF3,data=PRED_GRID, type = "response")
PRED_GRID <- PRED_GRID %>%
  dplyr::mutate(FEEDER_observed = FEEDER) %>%
  dplyr::mutate(FEEDER_predicted=PRED$predictions[,2]) %>%
  dplyr::mutate(feed_obs_num = ifelse(FEEDER=="NO",0,1))
hist(PRED$predictions[,2])
dim(PRED_GRID)

ROC_grid<-pROC::roc(data=PRED_GRID,response=FEEDER_observed,predictor=FEEDER_predicted)
AUC<-pROC::auc(ROC_grid)
THRESH_grd<-pROC::coords(ROC_grid, "best", "threshold")$threshold
THRESH<-table(PRED_GRID$FEEDER)[2]/dim(PRED_GRID)[1]


#### BASIC STATISTICS OF PREDICTIONS
hist(PRED_GRID$FEEDER_predicted)
range(PRED_GRID$FEEDER_predicted)
range(PRED_GRID$prop_feed)
range(PRED_GRID$prop_pts)
names(PRED_GRID)



########## CREATE OUTPUT GRID WITH PREDICTED FEEDING LOCATIONS ########################

OUTgrid<-countgrid %>% select(-FEEDER) %>%
  mutate(gridid=seq_along(n)) %>%
  left_join(PRED_GRID, by=c("gridid","n","N_ind","N_feed_points","N_feed_ind","prop_feed","prop_pts")) %>%
  filter(!is.na(FEEDER_predicted))
  


BI2<-Boyce(obs = PRED_GRID$FEEDER_observed , pred = PRED_GRID$FEEDER_predicted, n.bins=10)$Boyce
BI2



#### SHOW FINAL OUTPUT OF PREDICTIONS

tmap_mode("view")
tm_basemap(server="OpenStreetMap") +
  tm_shape(OUTgrid)  +
  tm_polygons(col = 'FEEDER_predicted', fill='FEEDER_predicted', alpha=0.4)


save.image("output/Feeding_analysis2025.RData")  



##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
########## PLOT THE MAP FOR KNOWN AND OBSERVED FEEDING STATIONS   #############
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################

########## CREATE A LEAFLET MAP OF PREDICTED FEEDING LOCATIONS #####################

icons <- awesomeIcons(
  icon = "fast-food-outline",
  iconColor = 'black',
  library = 'ion',
  markerColor = 'red'
)
 
# If you want to set your own colors manually:
pred.pal <- colorNumeric(c("cornflowerblue","firebrick"), seq(0,1))

m4 <- leaflet(options = leafletOptions(zoomControl = F)) %>% #changes position of zoom symbol
  setView(lng = mean(st_coordinates(OUT_sf %>%  st_transform(4326))[,1]),
          lat = mean(st_coordinates(OUT_sf %>%  st_transform(4326))[,2]),
          zoom = 8) %>%
  htmlwidgets::onRender("function(el, x) {L.control.zoom({ 
                           position: 'bottomright' }).addTo(this)}"
  ) %>% #Esri.WorldTopoMap #Stamen.Terrain #OpenTopoMap #Esri.WorldImagery
  addProviderTiles("Esri.WorldImagery", group = "Satellite",
                   options = providerTileOptions(opacity = 0.3, attribution = F,minZoom = 5, maxZoom = 14)) %>%
  addProviderTiles("OpenTopoMap", group = "Roadmap", options = providerTileOptions(attribution = F,minZoom = 5, maxZoom = 14)) %>%
  addLayersControl(baseGroups = c("Satellite", "Roadmap")) %>%

  addPolygons(
    data=SUI %>%
      st_transform(4326),
    stroke = TRUE, color = "black", weight = 3,
    fillColor = NULL, fillOpacity = 0
  ) %>%

  addAwesomeMarkers(data=plot_feeders %>%
      st_transform(4326),icon=icons, label=~as.character(Type), size=0.5, opacity = 0.5) %>%

  addPolygons(
    data=OUTgrid %>%
      st_transform(4326),
    stroke = TRUE,
	color = ~pred.pal(FEEDER_predicted), weight = 1.5,
    fillColor = ~pred.pal(FEEDER_predicted), fillOpacity = 0.5,
    popup = ~as.character(paste(round(FEEDER_predicted,3),"/ N_ind=",N_ind,"/ Prop feed pts=",round(prop_feed,3), sep=" ")),
    label = ~as.character(round(FEEDER_predicted,3))
  ) %>%
  
  addLegend(     # legend for predicted prob of feeding
    position = "topleft",
    pal = pred.pal,
    values = OUTgrid$FEEDER_predicted,
    opacity = 1,
    title = "Predicted probability of </br>anthropogenic feeding"
  ) %>%
  
  addScaleBar(position = "bottomright", options = scaleBarOptions(imperial = F))

m4

htmltools::save_html(html = m4, file = "C:/Users/sop/MAT/SUI_REKI_feed_projection2025.html")
mapview::mapshot(m4, url = "C:/Users/sop/MAT/SUI_REKI_feed_projection.html")
htmlwidgets::saveWidget(m4,"C:/Users/sop/MAT/SUI_REKI_feed_projection.html", selfcontained = T)


save.image("output/Feeding_grid_CH_2025.RData")

