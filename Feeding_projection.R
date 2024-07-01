###########################################################################################################################
###### PROJECTING ANALYSIS TO IDENTIFY ANTHROPOGENIC FEEDING SITES OF RED KITES ACROSS SWITZERLAND ################
###########################################################################################################################
# original idea by Nathalie Heiniger (MSc thesis 2020)
# uses output provided by "Feeding_analysis.r"
# created by steffen.oppel@vogelwarte.ch in July 2023
# building and forest layers provided by Jerome Guelat

# fully revised in June 2024 to include high resolution data

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

load("output/Feeding_analysis2024.RData")

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
track_sf<-fread("data/REKI_annotated_feeding2024_CH.csv") %>%
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
dim(countgrid)/359869
1-(dim(countgrid)[1]/359869) ## reported in results - fraction eliminated


# summary(log(countgrid$n+1))
# pal <- colorNumeric(c("cornflowerblue","firebrick"), seq(0,12))
# leaflet(options = leafletOptions(zoomControl = F)) %>% #changes position of zoom symbol
#   htmlwidgets::onRender("function(el, x) {L.control.zoom({ 
#                            position: 'bottomright' }).addTo(this)}"
#   ) %>% #Esri.WorldTopoMap #Stamen.Terrain #OpenTopoMap #Esri.WorldImagery
#   addProviderTiles("Esri.WorldImagery", group = "Satellite",
#                    options = providerTileOptions(opacity = 0.6, attribution = F)) %>%
#   addProviderTiles("OpenTopoMap", group = "Roadmap", options = providerTileOptions(attribution = F)) %>%  
#   addLayersControl(baseGroups = c("Satellite", "Roadmap")) %>%  
#   
#   addCircleMarkers(
#     data=track_sf %>% st_transform(4326) %>% slice_sample(n=1000),
#     radius = 1,
#     stroke = TRUE, color = "goldenrod", weight = 0.8,
#     fillColor = "goldenrod", fillOpacity = 0.8
#   ) %>%
#   
#   addPolygons(
#     data=st_as_sfc(st_bbox(track_sf %>% st_transform(4326))),
#     stroke = TRUE, color = "red", weight = 1,
#     fillOpacity = 0.5
#   )
#   # 
#   # addPolygons(
#   #   data=grid_CH %>%
#   #     st_transform(4326),
#   #   stroke = TRUE, color = "green", weight = 1,
#   #   fillColor = "grey17", fillOpacity = 0.5
#   # )
# 



# ##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
# ########## PREDICTING INDIVIDUAL LOCATIONS THAT MATCH FEEDING LOCATION PATTERN   #############
# ##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################


#### apply RF2 model across all data

PRED<-stats::predict(RF2,data=DATA_CH, type = "response")

DATA_CH <- DATA_CH %>%
  dplyr::bind_cols(PRED$predictions) %>%
  dplyr::rename(no_feed_prob = NO, feed_prob = YES) %>%
  dplyr::mutate(FEEDER_predicted = as.factor(dplyr::case_when(feed_prob > THRESH_pts ~ "YES",
                                                              feed_prob < THRESH_pts ~ "NO")))



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

PRED<-stats::predict(RF3,data=PROJ_GRID, type = "response")
PROJ_GRID <- PROJ_GRID %>%
  dplyr::mutate(FEEDER_predicted=PRED$predictions[,2])
dim(PRED$predictions)
dim(PROJ_GRID)


hist(PROJ_GRID$FEEDER_predicted)


########## CREATE OUTPUT GRID WITH PREDICTED FEEDING LOCATIONS ########################

CHgrid<-countgrid %>%
  mutate(gridid=seq_along(n)) %>%
  left_join(PROJ_GRID, by=c("gridid","n","N_ind","N_feed_points","N_feed_ind","prop_feed","prop_pts")) %>%
  filter(!is.na(FEEDER_predicted))

range(CHgrid$FEEDER_predicted)

FEEDgrid<-CHgrid %>%
  filter(FEEDER_predicted>THRESH)

dim(FEEDgrid)[1]/dim(CHgrid)[1]


saveRDS(CHgrid,"output/REKI_feeding_grid2024.rds")



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


# htmltools::save_html(html = m4, file = "C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/REKIfeeding/output/Potential_feeding_grids_CH.html")
# mapview::mapshot(m4, url = "C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/REKIfeeding/output/Potential_feeding_grids_CH.html")
# st_write(CHgrid,"output/REKI_predicted_anthropogenic_feeding_areas_CH.kml",append=FALSE)
# saveRDS(CHgrid, "output/pred_anthro_feed_grid.rds")
htmltools::save_html(html = m4, file = "C:/Users/sop/OneDrive - Vogelwarte/General/MANUSCRIPTS/AnthropFeeding/Figure_2.html")
mapview::mapshot(m4, url = "C:/Users/sop/OneDrive - Vogelwarte/General/MANUSCRIPTS/AnthropFeeding/Figure_2.html")
htmlwidgets::saveWidget(m4,"C:/Users/sop/OneDrive - Vogelwarte/General/MANUSCRIPTS/AnthropFeeding/Figure_2.html", selfcontained = T)


save.image("output/Feeding_grid_CH.RData")  
load("output/Feeding_grid_CH.RData")




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
rmarkdown::render('C:\\Users\\sop\\OneDrive - Vogelwarte\\REKI\\Analysis\\REKIfeeding\\Feeding_Analysis_report2023_v2.Rmd',
                  output_file = "Feeding_Analysis_report2023_v2.html",
                  output_dir = 'C:\\Users\\sop\\OneDrive - Vogelwarte\\REKI\\Analysis\\REKIfeeding\\output')

