###########################################################################################################################
###### PROJECTING ANALYSIS TO IDENTIFY ANTHROPOGENIC FEEDING SITES OF RED KITES ACROSS SWITZERLAND ################
###########################################################################################################################
# original idea by Nathalie Heiniger (MSc thesis 2020)
# uses output provided by "Feeding_analysis.r"
# created by steffen.oppel@vogelwarte.ch in July 2023
# building and forest layers provided by Jerome Guelat

# adapted to 2025 data including UHF data to make projections for all animals

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
library("rnaturalearth")
library("rnaturalearthdata")


### read in Switzerland map
# SUI<-st_read("S:/rasters/outline_maps/swiss_map_overview/layers.gpkg") %>% filter(country=="Switzerland")
# saveRDS(SUI,"data/Swiss_border.rds")
SUI<-readRDS("data/Swiss_border.rds")
plot(SUI)

# STUDY_AREA<-st_read("C:/STEFFEN/OneDrive - Vogelwarte/General/DATA/REKI_study_area_sm.kml")
# saveRDS(STUDY_AREA,"data/REKI_study_area.rds")
STUDY_AREA<-readRDS("data/REKI_study_area.rds")
plot(STUDY_AREA)

## set root folder for project
try(setwd("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/REKIfeeding"),silent=T)
try(setwd("C:/STEFFEN/OneDrive - Vogelwarte/REKI/Analysis/REKIfeeding"),silent=T)


###### LOADING DATA -----------------------------------------------------------------

load("output/Feeding_analysis2025.RData")

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



### LOAD THE TRACKING DATA AND INDIVIDUAL SEASON SUMMARIES 
track_sf<-fread("data/REKI_annotated_feeding2025_CH.csv") %>%
  mutate(extra=year_id) %>%
  separate_wider_delim(extra, delim="_", names=c("year","bird_id"))
  # st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
  # st_transform(3035)
dim(track_sf)



###### APPLY MODELS ACROSS FULL EXTENT OF TRACKING DATA  -----------------------------------------------------------------


################ CREATE SUBSETS TO APPLY MODEL TO DAYTIME LOCATIONS OUTSIDE OF FOREST
track_nofor_day_build<-track_sf %>%
  filter(FOREST==0) %>%
  filter(tod_=="day") %>%
  filter(!(BUILD==0)) 


################ PREPARE DATA FOR PREDICTION

DATA_CH <- track_nofor_day_build %>%
  mutate(YDAY=yday(t_), hour=hour(t_), month=month(t_)) %>%
  filter(!is.na(step_length)) %>%
  filter(!is.na(turning_angle)) %>%
  filter(!is.na(speed)) %>%
  filter(!is.na(mean_speed)) %>%
  filter(!is.na(mean_angle)) %>%
  select(-tod_,-FOREST) %>%
  mutate(point_id=seq_along(t_))
head(DATA_CH)
dim(DATA_CH)


# CREATING BASELINE MAP OF SAMPLING INTENSITY -----------------------------------------------------------------

track_sf<-track_sf %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
  st_transform(3035)

grid_CH <- forest_CH %>% 
  st_make_grid(cellsize = 500, what = "polygons",
               square = FALSE) # This statements leads to hexagons

tab <- st_intersects(grid_CH, track_sf)
lengths(tab)
countgrid <- st_sf(n = lengths(tab), geometry = st_cast(grid_CH, "MULTIPOLYGON")) %>%
  st_transform(3035) %>%
  filter(n>20)

dim(grid_CH)
dim(countgrid)


# ##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
# ########## PREDICTING INDIVIDUAL LOCATIONS THAT MATCH FEEDING LOCATION PATTERN   #############
# ##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################


#### apply RF2 model across all data
RF2<-readRDS("output/feed_site_RF_model2025.rds")
PRED<-stats::predict(RF2,data=DATA_CH, type = "response")
hist(PRED$predictions)
DATA_CH <- DATA_CH %>%
  dplyr::bind_cols(PRED$predictions) %>%
  dplyr::rename(no_feed_prob = NO, feed_prob = YES) %>%
  dplyr::mutate(FEEDER_predicted = as.factor(dplyr::case_when(feed_prob > 0.5 ~ "YES",   ### replaced THRESH_pts with 0.5 because histogram is pretty bimodal
                                                              feed_prob < 0.5 ~ "NO")))



##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
########## COUNT POINTS AND INDIVIDUALS IN GRID AND PREDICT FEEDING GRIDS   #############
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################

#### COUNT N INDIVIDUALS PER GRID CELL
tabind <- st_intersects(countgrid, track_sf)
lengths(tabind)
for(c in 1:length(tabind)){
  countgrid$N_ind[c]<-length(unique(track_sf$bird_id[tabind[c][[1]]]))
  #countgrid$N_ind[c]<-length(unique(track_sf$year_id[tabind[c][[1]]]))
}

#### COUNT N INDIVIDUALS AND PREDICTED FEEDING LOCS PER GRID CELL
OUT_sf<-DATA_CH %>%
  filter(FEEDER_predicted=="YES") %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
  st_transform(3035)
tabfeed <- st_intersects(countgrid, OUT_sf)
countgrid$N_feed_points<-lengths(tabfeed)
for(c in 1:length(tabfeed)){
  #countgrid$N_feed_ind[c]<-length(unique(OUT_sf$year_id[tabfeed[c][[1]]]))
  countgrid$N_feed_ind[c]<-length(unique(OUT_sf$bird_id[tabfeed[c][[1]]]))
}


### MANIPULATE COUNTED ANIMALS INTO PROPORTIONS

countgrid<-countgrid %>%
  #filter(n>10) %>%
  mutate(prop_feed=N_feed_ind/N_ind) %>%
  mutate(prop_pts=N_feed_points/n)

hist(countgrid$prop_feed)
hist(countgrid$prop_pts)



##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
########## PREDICT SECOND RANDOM FOREST MODEL TO PREDICT FEEDING SITE AT GRID CELL LEVEL   #############
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################

PROJ_GRID<-countgrid %>% 
  mutate(gridid=seq_along(n)) %>%
  filter(n>10) %>%
  filter(prop_feed>-1) %>% ## to exclude NaN that occur when no points occur
  st_drop_geometry()
str(PROJ_GRID)
PROJ_GRID %>% filter(is.na(prop_feed))

#### classification success of training data
RF3<-readRDS("output/feed_grid_RF_model2025.rds")
PRED<-stats::predict(RF3,data=PROJ_GRID, type = "response")
PROJ_GRID <- PROJ_GRID %>%
  dplyr::mutate(FEEDER_predicted=PRED$predictions[,2])
dim(PRED$predictions)
dim(PROJ_GRID)


hist(PROJ_GRID$FEEDER_predicted)


########## CREATE OUTPUT GRID WITH PREDICTED FEEDING LOCATIONS ########################
THRESH<-table(PRED_GRID$FEEDER)[2]/dim(PRED_GRID)[1]
CHgrid<-countgrid %>%
  mutate(gridid=seq_along(n)) %>%
  left_join(PROJ_GRID, by=c("gridid","n","N_ind","N_feed_points","N_feed_ind","prop_feed","prop_pts")) %>%
  filter(!is.na(FEEDER_predicted))

range(CHgrid$FEEDER_predicted)
hist(CHgrid$FEEDER_predicted)

FEEDgrid<-CHgrid %>%
  filter(FEEDER_predicted>THRESH)

dim(FEEDgrid)[1]/dim(CHgrid)[1]


saveRDS(CHgrid,"output/REKI_feeding_grid2025.rds")
CHgrid<-readRDS("output/REKI_feeding_grid2025.rds")


export<-CHgrid %>% select(N_ind,N_feed_points,N_feed_ind,prop_feed,prop_pts,n,FEEDER_predicted, gridid, geometry)
st_write(export, "C:/Users/sop/MAT/REKI_feeding_probability_2025.gpkg")



##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
########## PLOT THE MAP FOR KNOWN AND OBSERVED FEEDING STATIONS   #############
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################

########## CREATE A LEAFLET MAP OF PREDICTED FEEDING LOCATIONS ########################
SUIpol <- ne_countries(scale = "large", returnclass = "sf") %>%
  filter(admin=="Switzerland") %>%
  st_transform(3035)
CHgrid<-CHgrid %>% st_intersection(.,SUIpol) %>% filter(!is.na(admin))
FEEDgrid<-FEEDgrid %>% st_intersection(.,SUIpol) %>% filter(!is.na(admin))
## need to specify color palette 
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
    data=CHgrid %>%
      st_transform(4326),
    stroke = TRUE, color = ~pred.pal(FEEDER_predicted), weight = 1,
    fillColor = ~pred.pal(FEEDER_predicted), fillOpacity = 0.5
  ) %>%
  addPolygons(
    data=FEEDgrid %>%
      st_transform(4326),
    stroke = TRUE, color = "red", weight = 1.5,
    fillColor = ~pred.pal(FEEDER_predicted), fillOpacity = 0.5,
    popup = ~as.character(paste(round(FEEDER_predicted,3),"/ N_ind=",N_ind,"/ Prop feed pts=",round(prop_feed,3), sep=" ")),
    label = ~as.character(round(FEEDER_predicted,3))
  ) %>%
  addPolygons(
    data=SUI %>%
      st_transform(4326),
    stroke = TRUE, color = "black", weight = 3,
    fillColor = NULL, fillOpacity = 0
  ) %>%
  
  addLegend(     # legend for predicted prob of feeding
    position = "topleft",
    pal = pred.pal,
    values = CHgrid$FEEDER_predicted,
    opacity = 1,
    title = "Predicted probability of </br>anthropogenic feeding"
  ) %>%
  
  addScaleBar(position = "bottomright", options = scaleBarOptions(imperial = F))

m4

htmltools::save_html(html = m4, file = "C:/Users/sop/MAT/SUI_REKI_feed_projection.html")
mapview::mapshot(m4, url = "C:/Users/sop/MAT/SUI_REKI_feed_projection.html")
htmlwidgets::saveWidget(m4,"C:/Users/sop/MAT/SUI_REKI_feed_projection.html", selfcontained = T)


save.image("output/Feeding_grid_CH_2025.RData")
