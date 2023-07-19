###########################################################################################################################
###### DATA ANALYSIS TO IDENTIFY ANTHROPOGENIC FEEDING SITES OF RED KITES BASED ON GPS TRACKING DATA ################
###########################################################################################################################
# original idea by Nathalie Heiniger (MSc thesis 2020)
# uses csv table prepared by "Feeding_data_prep.r"
# created by steffen.oppel@vogelwarte.ch in June 2023

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

## set root folder for project
setwd("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/REKIfeeding")

### TO DO :
## create subset of known feeding locations for training the model and retain some for testing (from similar point intensity locations)


# LOADING DATA -----------------------------------------------------------------
### LOAD THE TRACKING DATA AND INDIVIDUAL SEASON SUMMARIES 
track_sf<-fread("data/REKI_annotated_feeding.csv")
#track_sf<-load(file = "data/REKI_trackingdata_annotated.rds")

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



################ CREATE VARIOUS SUBSETS TO IMPROVE RATIO OF CLASS MEMBERSHIP
track_day<-track_sf %>% dplyr::filter(tod_=="day")
track_nofor<-track_sf %>% filter(FOREST==0)
track_nofor_day<-track_sf %>% filter(FOREST==0) %>% filter(tod_=="day")
table(track_nofor_day$FEEDER)

track_nofor_day_build<-track_nofor_day %>%
  filter(!(BUILD==0 & FEEDER =="NO")) 

table(track_nofor_day_build$FEEDER)


################ PREPARE DATA FOR PREDICTION

DATA <- track_nofor_day_build %>%
  mutate(YDAY=yday(t_), hour=hour(t_), month=month(t_)) %>%
  filter(!is.na(step_length)) %>%
  filter(!is.na(turning_angle)) %>%
  filter(!is.na(speed)) %>%
  filter(!is.na(mean_speed)) %>%
  filter(!is.na(mean_angle)) %>%
  select(-tod_,-FOREST,-forest_size,-build_id,-type_of_food,-frequency) %>%
  mutate(point_id=seq_along(t_))
head(DATA)


### basic stats
length(unique(DATA$year_id))
length(unique(DATA$bird_id))
length(unique(DATA$FEED_ID))


### random feeder subsetting
## results in poor transferability
# nfeed<-length(unique(DATA$FEED_ID,na.rm=T))
# selectids<-sample(unique(DATA$FEED_ID),nfeed*0.67)
# 
# DATA_TRAIN<- DATA %>% filter(FEED_ID %in% selectids)
# DATA_TEST<- DATA %>% filter(!(FEED_ID %in% selectids))


# nind<-length(unique(DATA$bird_id))
# selectids<-sample(unique(DATA$bird_id),nind*0.67)
# 
# DATA_TRAIN<- DATA %>% filter(bird_id %in% selectids)
# DATA_TEST<- DATA %>% filter(!(bird_id %in% selectids))

### completely random subsetting
DATA_TRAIN<- DATA %>% slice_sample(prop=0.67, by = FEEDER, replace = FALSE)
DATA_TEST<- DATA %>% filter(!(point_id %in% DATA_TRAIN$point_id))

table(DATA_TRAIN$FEEDER)
table(DATA_TEST$FEEDER)



# DATA %>% filter(FEED_ID=="") %>% filter(FEEDER=="YES")
# DATA %>% filter(!(FEED_ID=="")) %>% filter(FEEDER=="NO")   ### these are feeders that were only temporarily active, so feeder=NO, despite location in spatial proximity

# CREATING BASELINE MAP OF SAMPLING INTENSITY -----------------------------------------------------------------

grid <- forest %>% 
  st_make_grid(cellsize = 500, what = "polygons",
               square = FALSE) # This statements leads to hexagons
track_sf <-   st_as_sf(track_sf, coords = c("long", "lat"), crs = 4326) %>%
  st_transform(3035)
tab <- st_intersects(grid, track_sf)
lengths(tab)
countgrid <- st_sf(n = lengths(tab), geometry = st_cast(grid, "MULTIPOLYGON")) %>%
  st_transform(3035)
summary(log(countgrid$n+1))
pal <- colorNumeric(c("cornflowerblue","firebrick"), seq(0,15))
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

# ANALYSING DATA -----------------------------------------------------------------



############### LOOP OVER TUNING SETTINGS TO IMPROVE PERFORMANCE ##################
### this takes too long >3 hrs
# tuning.out<-expand.grid(m=c(4,8,14),t=c(500,750,1000,2000)) %>%
#   dplyr::mutate(oob.error=0)
# for (m in c(4,8,14)) {
#   for(t in c(500,750,1000,2000)){
#     RFtest<-ranger::ranger(as.factor(FEEDER) ~ sex + revisits + residence_time + age_cy +
#                              step_length + turning_angle + speed + 
#                            mean_speed + mean_angle + YDAY + hour + month + BUILD,
#                            data = DATA, mtry = m, num.trees = t, replace = T, importance = "permutation", oob.error=T, write.forest=F)
#     tuning.out[tuning.out$t==t & tuning.out$m==m,3] <-RFtest$prediction.error
#   }
# }
# 
# tuning.out<-tuning.out %>% dplyr::arrange(oob.error)


##### RUN MODEL ##########


RF2 <- ranger::ranger(as.factor(FEEDER) ~ sex + age_cy + YDAY + hour + month +                            ## basic variables such as age, sex, and time
                        revisits + residence_time + meanFreqVisit + n_days + TimeSpan + TempEven +        ## several variables dealing with the temporal revisitation pattern
                        step_length + turning_angle + speed +                                             ## several variables dealing with the current movement characteristics
                        mean_speed + mean_angle +                                                         ## variables dealing with the average movement characteristics
                        BUILD + NEST + dist_nest,     #                                                      ## variables dealing with distance to structures and nests
                      data=DATA_TRAIN, mtry=2, num.trees=2500, replace=F, importance="permutation", oob.error=T, write.forest=T, probability=T)
#saveRDS(RF2, "output/feed_site_RF_model.rds")
IMP<-as.data.frame(RF2$variable.importance) %>%
  dplyr::mutate(variable=names(RF2$variable.importance)) %>%
  dplyr::rename(red.accuracy=`RF2$variable.importance`) %>%
  dplyr::arrange(dplyr::desc(red.accuracy)) %>%
  dplyr::mutate(rel.imp=(red.accuracy/max(red.accuracy))*100) %>%
  dplyr::select(variable,red.accuracy,rel.imp)


#### classification success of training data

PRED<-stats::predict(RF2,data=DATA_TRAIN, type = "response")
prevalence<-table(DATA_TRAIN$FEEDER)[2]/sum(table(DATA_TRAIN$FEEDER))
DATA_TRAIN <- DATA_TRAIN %>%
  dplyr::mutate(FEEDER_observed = factor(FEEDER,levels=c("NO","YES"))) %>%
  dplyr::bind_cols(PRED$predictions) %>%
  dplyr::rename(no_feed_prob = NO, feed_prob = YES) %>%
  dplyr::mutate(FEEDER_predicted = as.factor(dplyr::case_when(feed_prob > prevalence ~ "YES",
                                                            feed_prob < prevalence ~ "NO")))

suppressWarnings({trainmat<-caret::confusionMatrix(data = DATA_TRAIN$FEEDER_observed, reference = DATA_TRAIN$FEEDER_predicted, positive="YES")})

#### classification success of test data

PRED<-stats::predict(RF2,data=DATA_TEST, type = "response")
prevalence<-table(DATA_TEST$FEEDER)[2]/sum(table(DATA_TEST$FEEDER))
DATA_TEST <- DATA_TEST %>%
  dplyr::mutate(FEEDER_observed = factor(FEEDER,levels=c("NO","YES"))) %>%
  dplyr::bind_cols(PRED$predictions) %>%
  dplyr::rename(no_feed_prob = NO, feed_prob = YES) %>%
  dplyr::mutate(FEEDER_predicted = as.factor(dplyr::case_when(feed_prob > prevalence ~ "YES",
                                                              feed_prob < prevalence ~ "NO")))

suppressWarnings({testmat<-caret::confusionMatrix(data = DATA_TEST$FEEDER_observed, reference = DATA_TEST$FEEDER_predicted, positive="YES")})

## export data for further use in outcome prediction
OUT<-dplyr::bind_rows(DATA_TEST, DATA_TRAIN)


#### CREATE PLOT FOR VARIABLE IMPORTANCE
  mylevels<-IMP$variable[10:1]
  impplot<-IMP[10:1,] %>%
    dplyr::mutate(variable=forcats::fct_relevel(variable,mylevels)) %>%
    ggplot2::ggplot(ggplot2::aes(x=variable, y=rel.imp)) +
    ggplot2::geom_bar(stat='identity', fill='lightblue') +
    ggplot2::coord_flip()+
    ggplot2::ylab("Variable importance (%)") +
    ggplot2::xlab("Explanatory variable") +
    ggplot2::scale_y_continuous(limits=c(-5,105), breaks=seq(0,100,20), labels=seq(0,100,20))+
    ggplot2::annotate("text",x=2,y=80,label=paste("Accuracy = ",round(testmat$byClass[3],3)),size=8) +
    ggplot2::theme(panel.background=ggplot2::element_rect(fill="white", colour="black"), 
                   axis.text.x=ggplot2::element_text(size=18, color="black"),
                   axis.text.y=ggplot2::element_text(size=16, color="black"), 
                   axis.title=ggplot2::element_text(size=20), 
                   panel.grid.major = ggplot2::element_blank(), 
                   panel.grid.minor = ggplot2::element_blank(), 
                   panel.border = ggplot2::element_blank())
  print(impplot)




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
  st_transform(crs = 4326)%>%
  mutate(Type="Experimental") %>%
  select(Type)

plot_feeders<-rbind(plot_feeders,plot_feeders2,plot_feeders3)

## FIRST, split the multipolygon into separate polygons:
# feed_site_polygons <- st_cast(feed_site_sf, "POLYGON")
# #nearneighb<-st_nearest_feature(feed_site_polygons,plot_feeders)
# feed_site_distances<-st_distance(feed_site_polygons,plot_feeders) 
# feed_site_polygons <- feed_site_polygons %>%
#   mutate(dist_feeder=apply(feed_site_distances,1,min)/1000) %>%  ### distance in km
#   arrange(desc(dist_feeder))




# 
# ##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
# ########## PLOT THE MAP FOR KNOWN AND OBSERVED FEEDING STATIONS   #############
# ##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
#   
#   ########## CREATE A LEAFLET MAP OF PREDICTED FEEDING LOCATIONS ########################
#  
# head(OUT)    
#   
# # create plotting frame
#   
plot_OUT<-OUT %>%
  filter(FEEDER_predicted=="YES") %>%
    st_as_sf(coords = c("long", "lat"), crs = 4326)



##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
########## ALTERNATIVE APPROACH TO SIMPLY COUNT POINTS AND INDIVIDUALS IN GRID   #############
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################

#### FIRST COUNT N INDIVIDUALS PER GRID CELL

for(c in 1:length(tab)){
  countgrid$N_ind[c]<-length(unique(track_sf$year_id[tab[c][[1]]]))
}

#### SECOND COUNT N INDIVIDUALS AND PREDICTED FEEDING LOCS PER GRID CELL
OUT_sf<-OUT %>%
  filter(FEEDER_predicted=="YES") %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
  st_transform(3035)
tab2 <- st_intersects(grid, OUT_sf)
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
########## FIT SECOND RANDOM FOREST MODEL TO PREDICT FEEDING SITE AT GRID CELL LEVEL   #############
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################

### overlay feeder data with countgrid
feed_grd<-st_intersects(countgrid,(plot_feeders %>% st_transform(3035)))
countgrid$FEEDER<-lengths(feed_grd) 


PRED_GRID<-countgrid %>% 
  mutate(gridid=seq_along(n)) %>%
  filter(n>10) %>%
  st_drop_geometry() %>%
  mutate(FEEDER=ifelse(FEEDER==0,0,1))

RF3 <- ranger::ranger(FEEDER ~ n+N_ind+N_feed_points+N_feed_ind+prop_feed+prop_pts,     
                      data=PRED_GRID, mtry=2, num.trees=2500, replace=F, importance="permutation", oob.error=T, write.forest=T, probability=T)
saveRDS(RF3, "output/feed_grid_RF_model.rds")
IMP2<-as.data.frame(RF3$variable.importance) %>%
  dplyr::mutate(variable=names(RF3$variable.importance)) %>%
  dplyr::rename(red.accuracy=`RF3$variable.importance`) %>%
  dplyr::arrange(dplyr::desc(red.accuracy)) %>%
  dplyr::mutate(rel.imp=(red.accuracy/max(red.accuracy))*100) %>%
  dplyr::select(variable,red.accuracy,rel.imp)


#### classification success of training data

PRED<-stats::predict(RF3,data=PRED_GRID, type = "response")
prevalence<-table(DATA_TRAIN$FEEDER)[2]/sum(table(DATA_TRAIN$FEEDER))
PRED_GRID <- PRED_GRID %>%
  dplyr::mutate(FEEDER_observed = FEEDER) %>%
  dplyr::mutate(FEEDER_predicted=PRED$predictions[,2])
dim(PRED$predictions)
dim(PRED_GRID)

AUC<-auc(roc(data=PRED_GRID,response=FEEDER_observed,predictor=FEEDER_predicted))
#### CREATE PLOT FOR VARIABLE IMPORTANCE
mylevels2<-IMP2$variable[6:1]
impplot2<-IMP2[6:1,] %>%
  dplyr::mutate(variable=forcats::fct_relevel(variable,mylevels2)) %>%
  ggplot2::ggplot(ggplot2::aes(x=variable, y=rel.imp)) +
  ggplot2::geom_bar(stat='identity', fill='lightblue') +
  ggplot2::coord_flip()+
  ggplot2::ylab("Variable importance (%)") +
  ggplot2::xlab("Explanatory variable") +
  ggplot2::scale_y_continuous(limits=c(-5,105), breaks=seq(0,100,20), labels=seq(0,100,20))+
  ggplot2::annotate("text",x=2,y=80,label=paste("AUC = ",round(AUC,3)),size=8) +
  ggplot2::theme(panel.background=ggplot2::element_rect(fill="white", colour="black"), 
                 axis.text.x=ggplot2::element_text(size=18, color="black"),
                 axis.text.y=ggplot2::element_text(size=16, color="black"), 
                 axis.title=ggplot2::element_text(size=20), 
                 panel.grid.major = ggplot2::element_blank(), 
                 panel.grid.minor = ggplot2::element_blank(), 
                 panel.border = ggplot2::element_blank())
print(impplot2)




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
setwd("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/REKIfeeding")
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
########## VALIDATE PREDICTIONS WITH INTERVIEW DATA FROM FIONA PELLE  #############
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
## survey provided with addresses only
## Jerome Guelat provided R script to convert addresses to coordinates
source("C:/Users/sop/OneDrive - Vogelwarte/General/ANALYSES/DataPrep/swisstopo_address_lookup.r")
library(stringi)

## Feeding is the question whether they feed or not
setwd("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/REKIfeeding")
feed_surveys2<-read_csv("data/FeedersFionaPelle.csv", locale = locale(encoding = "UTF-8")) %>%
#<-fread("data/FeedersFionaPelle.csv", encoding = 'UTF-8') %>% 
  mutate(FEEDER_surveyed=ifelse(Feeding=="Yes",1,0))

## generate coordinates from addresses
feed_surveys2_locs<-swissgeocode(address=as.character(feed_surveys2$Address), nresults=3)

feed_surv2_sf<-feed_surveys2_locs %>% rename(Address=address_origin) %>%
  left_join(feed_surveys2, by="Address",relationship ="many-to-many") %>%
  filter(!is.na(x)) %>%
  filter(!is.na(FEEDER_surveyed)) %>%
  st_as_sf(coords = c("lon", "lat"), crs=4326) %>%
  st_transform(crs = 3035) 
feed_surv2_sf$FEEDER_surveyed

validat2 <- st_intersection(OUTgrid, feed_surv2_sf)
validat2$FEEDER_predicted
AUC_TEST2<-auc(roc(data=validat2,response=FEEDER_surveyed,predictor=FEEDER_predicted))

### summarise the predicted sites
validat2<-validat2 %>% filter(FEEDER_surveyed==1)
hist(validat2$FEEDER_predicted)
mean(validat2$FEEDER_predicted)
length(validat2[validat2$FEEDER_predicted>0.5,])/dim(validat2)[1]




##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
########## COMBINE VALIDATION DATA  #############
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
VAL_DAT<-bind_rows(validat %>% select(n,N_ind,N_feed_points,N_feed_ind,prop_feed,prop_pts,FEEDER_observed,FEEDER_predicted),
validat2 %>% select(n,N_ind,N_feed_points,N_feed_ind,prop_feed,prop_pts,FEEDER_observed,FEEDER_predicted)) %>%
  mutate(Classification=ifelse(FEEDER_predicted>0.5,"correct","missed"))


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
val.pal <- colorFactor(c("green","red"), unique(VAL_DAT$Classification))
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
    fillColor = ~pred.pal(FEEDER_predicted), fillOpacity = 0.5,
    popup = ~as.character(paste("N_pts=",n,"N_ind=",N_ind,"Prop feed pts=",round(prop_feed,3), sep=" / ")),
    label = ~as.character(round(FEEDER_predicted,3))
  ) %>%
  
  addCircleMarkers(
    data=plot_feeders,
    radius = 3,
    stroke = TRUE, color = ~feed.pal(Type), weight = 1,
    fillColor = ~feed.pal(Type), fillOpacity = 0.2
  ) %>%
  
  addCircleMarkers(
    data=VAL_DAT %>%
      st_transform(4326),
    radius = 5,
    stroke = TRUE, color = ~val.pal(Classification), weight = 1,
    fillColor = ~val.pal(Classification), fillOpacity = 1,
    popup = ~as.character(paste(round(FEEDER_predicted,3),"N_ind=",N_ind,"Prop feed pts=",round(prop_feed,3), sep=" / ")),
    label = ~as.character(round(FEEDER_predicted,3))
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
  addLegend(     # legend for known feeding sites
    position = "topleft",
    pal = val.pal,
    values = VAL_DAT$Classification,
    opacity = 1,
    title = "Validation (interviews)"
  ) %>%
  
  
  
  addScaleBar(position = "bottomright", options = scaleBarOptions(imperial = F))

m2


htmltools::save_html(html = m2, file = "C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/REKIfeeding/output/Potential_feeding_grids.html")
mapview::mapshot(m2, url = "C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/REKIfeeding/output/Potential_feeding_grids.html")
st_write(OUTgrid,"output/REKI_predicted_anthropogenic_feeding_areas.kml",append=FALSE)



##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
########## IDENTIFY AREAS WHERE PREDICTION DOES NOT WORK   #############
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################


### MISSED PREDICTION - FIGURE OUT WHY NOT PREDICTED AS FEEDER

### LOOK AT THE GRID CELLS THAT ARE NOT USED FOR RF DUE TO FEW LOCATIONS
## does not contain the missed predictions
# MISS_GRID<-countgrid %>% 
#   mutate(gridid=seq_along(n)) %>%
#   filter(n<11) 
# m2 %>% addPolygons(
#   data=MISS_GRID %>%
#     st_transform(4326),
#   stroke = TRUE, color = "red", weight = 3,
#   fillColor = "red", fillOpacity = 0.2)

### LOOK AT THE GRID CELLS THAT ARE NOT INCLUDED IN PLOT
# does not contain the missed predictions
MISS_GRID<-countgrid %>%
  mutate(gridid=seq_along(n)) %>%
  filter(!(gridid %in% PRED_GRID$gridid))
m2 %>% addPolygons(
  data=MISS_GRID %>%
    st_transform(4326),
  stroke = TRUE, color = "red", weight = 3,
  fillColor = "red", fillOpacity = 0.2)


ERROR_GRID<-OUTgrid %>% filter(FEEDER_observed==1) %>%
  filter(FEEDER_predicted<0.1)

m2 %>% addPolygons(
  data=ERROR_GRID %>%
    st_transform(4326),
  stroke = TRUE, color = "red", weight = 3,
  fillColor = "red", fillOpacity = 0.2)



### PREDICTION OF FEEDERS WHERE NONE EXISTS
ERROR_GRID<-OUTgrid %>% filter(FEEDER_observed==0) %>%
  filter(FEEDER_predicted>0.5)

m2 %>% addPolygons(
  data=ERROR_GRID %>%
    st_transform(4326),
  stroke = TRUE, color = "red", weight = 3,
  fillColor = "red", fillOpacity = 0.2)



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


save.image("output/Feeding_analysis.RData")  
load("output/Feeding_analysis.RData")














###########~~~~~~~~~~~~~     ABANDONED CODE       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
## abandoned kernelUD in early July 2023 because prediction with RF in grid cells was better than HR estimates and arbitrary polygons

# ##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
# ########## CALCULATE 50% kernelUD for positive predicted locations   #############
# ##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
### REMOVED on 3 July 2023 as simpler grid approach is faster and more efficient
# OUT_sp<-OUT %>%
#     filter(FEEDER_predicted=="YES") %>%
#     st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
#     st_transform(3035) %>%
#     select(year) %>%   ### everything below fails for year_id
#     as("Spatial")
#   
# plot(OUT_sp)  
#   ### FILTER OUT year_ID with too few locations
#   # exclude<-as_tibble(OUT_sp) %>% group_by(year_id) %>% summarise(n=length(coords.x1)) %>% dplyr::filter(n<10)
#   # OUT_sp<- OUT_sp[!(OUT_sp$year_id %in% exclude$year_id),] 
#   
#   
#   ### CREATE CUSTOM GRID for kernelUD (instead of same4all=TRUE) --------------
#   extendX <- max(diff(range(coordinates(OUT_sp)[, 1])) * 0.05, 50 * 100)
#   extendY <- max(diff(range(coordinates(OUT_sp)[, 2])) * 0.05, 50 * 100)
#   
#   minX <- min(coordinates(OUT_sp)[, 1]) - extendX
#   maxX <- max(coordinates(OUT_sp)[, 1]) + extendX
#   minY <- min(coordinates(OUT_sp)[, 2]) - extendY
#   maxY <- max(coordinates(OUT_sp)[, 2]) + extendY
#   res <- (max(abs(minX - maxX) / 1000, abs(minY - maxY) / 1000)) / 1000
#   xrange <- seq(minX, maxX, by = res * 1000)
#   yrange <- seq(minY, maxY, by = res * 1000)
#   grid.locs <- expand.grid(x = xrange, y = yrange)
#   INPUTgrid <- SpatialPixels(
#     SpatialPoints(grid.locs), proj4string = proj4string(OUT_sp)
#   )
#   
# #### convert grid to INPUT grid
# ## is too small for kernels  
# # INPUTgrid <- countgrid  %>%
# #   st_transform(3035) %>% st_rasterize() %>%
# #   as("Spatial")
# 
#   
# #### THIS FAILED FOR THE year_id level - can only be done at year level
#   plot(INPUTgrid,  col='grey87')
#   plot(OUT_sp, add=T, col='red')
# 
# ### THIS RUNS FOR >30 MINUTES
# KDE.Surface <-   adehabitatHR::kernelUD(OUT_sp,h = 100, grid = INPUTgrid, same4all = FALSE)
# #image(KDE.Surface)
# #KDE.Surface <-   adehabitatHR::kernelUD(OUT_sp,h = "href", grid=10000, same4all = TRUE)
# KDE_sp <- adehabitatHR::getverticeshr(KDE.Surface, percent=75, unin = "m", unout = "ha")
# 
# KDE_sf <- st_as_sf(KDE_sp) %>%
#   st_transform(4326)  
# 
# 
# 
# 
#     
# ### COUNT OVERLAPPING CORE AREAS
# KDE <- adehabitatHR::estUDm2spixdf(KDE.Surface)
# pixArea <- KDE@grid@cellsize[[1]]
# 
# ## output reported by kernelUD is intensity / m2. This intensity is multiplied
# # by pixel area and sums to 1 for each individual (exceptions near borders)
# ## we sort this calc cumulative sum --> this is effectively the "% UD"
# KDE@data <- KDE@data %>%
#   mutate(rowname = seq_len(nrow(KDE@data))) %>%
#   tidyr::pivot_longer(!rowname, names_to = "ID", values_to = "UD") %>%
#   mutate(usage = UD * (pixArea^2)) %>%
#   arrange(ID, desc(usage)) %>%
#   group_by(ID) %>%
#   mutate(cumulUD = cumsum(usage)) %>%
#   dplyr::select(rowname, ID, cumulUD) %>%
#   arrange(rowname) %>%
#   tidyr::pivot_wider(names_from = ID, values_from = cumulUD) %>%
#   dplyr::select(-rowname)
# 
# Noverlaps <- KDE
# 
# Noverlaps@data <- as.data.frame(
#   ifelse(Noverlaps@data < 0.5, 1, 0)
# ) %>%
#   mutate(N_YEAR = rowSums(.))
# 
# cols <- colnames(Noverlaps@data)[-which(colnames(Noverlaps@data) == "N_YEAR")]
# 
# Noverlaps@data <- Noverlaps@data %>% dplyr::select(N_YEAR)
# 
# ### THIS OPERATION IS VERY SLOW
# # OUTMAP <- raster::aggregate(
# #   as(Noverlaps, "SpatialPolygonsDataFrame"),
# #   c("N_YEAR")
# # )
# 
# ### CONVERT INTO SIMPLE FEATURE AS OUTPUT AND FOR PLOTTING
# feed_site_sf <- sf::st_as_sf(as(Noverlaps, "SpatialPolygonsDataFrame")) %>%
#   dplyr::filter(N_YEAR>1) %>%
#   sf::st_union(by_feature = TRUE) %>%
#   sf::st_transform(4326) %>%
#   arrange(N_YEAR)
# 
# length(feed_site_sf$geometry[[1]])
# 
# 
# 
  