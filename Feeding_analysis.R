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
sf_use_s2(FALSE)
library(pROC)
library(modEvA)  ### for calculating the Boyce index
library(pdp) ### for creating partial dependence plots

## set root folder for project
#setwd("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/REKIfeeding")

# LOADING DATA -----------------------------------------------------------------
### LOAD THE TRACKING DATA AND INDIVIDUAL SEASON SUMMARIES 
track_sf<-fread("data/REKI_annotated_feeding2024_CH.csv") %>%
  mutate(extra=year_id) %>%
  separate_wider_delim(extra, delim="_", names=c("year","bird_id"))

#track_sf<-load(file = "data/REKI_trackingdata_annotated2024_sf.rds")

#define different coordinate systems
CH_LV95_coords = "+init=epsg:2056"
EU_coords = "+init=epsg:3035"
WGS84_coords ="+init=epsg:4326"
CH_LV03_coords ="+init=epsg:21781"



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


##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
########## PROVIDE SUMMARY INFO FOR MANUSCRIPT   #############
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################

dim(track_sf)
hist(track_sf$turning_angle*(180/pi))

## calculate and inspect the average time between subsequent locations
# abandoned here and moved to 'REKI_data_download_high_res.r'
# pointless to do this after filtering because long gaps when birds left SUI!!!


dim(DATA)
length(unique(DATA$year_id))
range(DATA$age_cy)
range(DATA$year)
length(unique(DATA$bird_id))
length(unique(DATA$FEED_ID))
DATA %>% group_by(bird_id) %>%
  summarise(N= length(unique(year_id))) %>%
  ungroup() %>%
  summarise(avg=mean(N), max=max(N))

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
dim(DATA_TRAIN)
dim(DATA_TEST)

table(DATA_TRAIN$FEEDER)
table(DATA_TEST$FEEDER)

table(DATA_TRAIN$NEST)
table(DATA_TEST$NEST)

# DATA %>% filter(FEED_ID=="") %>% filter(FEEDER=="YES")
# DATA %>% filter(!(FEED_ID=="")) %>% filter(FEEDER=="NO")   ### these are feeders that were only temporarily active, so feeder=NO, despite location in spatial proximity

# CREATING BASELINE MAP OF SAMPLING INTENSITY -----------------------------------------------------------------

grid <- forest %>% 
  st_make_grid(cellsize = 500, what = "polygons",
               square = FALSE) # This statements leads to hexagons
sum(st_area(grid))/1000000  ## size of study area in sq km
track_sf <-   st_as_sf(track_sf, coords = c("long", "lat"), crs = 4326) %>%
  st_transform(3035)
tab <- st_intersects(grid, track_sf)
lengths(tab)
countgrid <- st_sf(n = lengths(tab), geometry = st_cast(grid, "MULTIPOLYGON")) %>%
  st_transform(3035)
summary(log(countgrid$n+1))
pal <- colorNumeric(c("cornflowerblue","firebrick"), seq(0,12))
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



########## CREATE A LEAFLET MAP TO SHOW APPROACH OF OVERLAYING FEEDING AND NON-FEEDING LOCATIONS ########################

### NEED TO FIND GOOD EXAMPLE BIRD FOR DEMO ###
sort(unique(DATA$FEED_ID))
DATA %>% filter (FEEDER=="YES") %>% filter(revisits>2000) %>% arrange(revisits)

## need to specify color palette 
# If you want to set your own colors manually:
res.pal <- colorNumeric(c("grey15","red"), seq(0,3300))
mF <- leaflet(options = leafletOptions(zoomControl = F)) %>% #changes position of zoom symbol
  setView(lng = mean(st_coordinates(plot_feeders)[,1]), lat = mean(st_coordinates(plot_feeders)[,2]), zoom = 11) %>%
  htmlwidgets::onRender("function(el, x) {L.control.zoom({ 
                           position: 'bottomright' }).addTo(this)}"
  ) %>% #Esri.WorldTopoMap #Stamen.Terrain #OpenTopoMap #Esri.WorldImagery
  addProviderTiles("Esri.WorldImagery", group = "Satellite",
                   options = providerTileOptions(opacity = 0.4, attribution = F,minZoom = 5, maxZoom = 20)) %>%
  # addPolylines(
  #   data=DATA %>% filter(year_id=="2018_459"),
  #   lng = ~long, lat = ~lat, group = ~year_id,
  #   color = "red", weight = 1
  #  ) %>%
  addCircleMarkers(
    data=track_sf %>% filter(year_id=="2019_481") %>% st_transform(4326),
    radius = 5,
    stroke = TRUE, color = ~res.pal(revisits), weight = 0.8,
    fillColor = ~res.pal(revisits), fillOpacity = 0.8
  ) %>%
  addCircleMarkers(
    data=plot_feeders,
    radius = 10,
    stroke = TRUE, color = "cornflowerblue", weight = 1,   ###~feed.pal(Type)
    fillColor = "cornflowerblue", fillOpacity = 0.8   ## ~feed.pal(Type)
  ) %>%

  
  addLegend(     # legend for predicted prob of feeding
    position = "topleft",
    pal = res.pal,
    values = DATA$revisits,
    opacity = 1,
    title = "Number of </br>repeat visits"
  ) %>%
  
  addScaleBar(position = "bottomright", options = scaleBarOptions(imperial = F))

mF




##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
########## ANALYSING DATA IN RANDOM FOREST  #############
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################



##### RUN MODEL ##########
## takes about 45 minutes with 2024 dataset

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
write.table(IMP,"clipboard", sep="\t")

#### classification success of training data

PRED<-stats::predict(RF2,data=DATA_TRAIN, type = "response")
prevalence<-table(DATA_TRAIN$FEEDER)[2]/sum(table(DATA_TRAIN$FEEDER))
DATA_TRAIN <- DATA_TRAIN %>%
  dplyr::mutate(FEEDER_observed = factor(FEEDER,levels=c("NO","YES"))) %>%
  dplyr::bind_cols(PRED$predictions) %>%
  dplyr::rename(no_feed_prob = NO, feed_prob = YES) %>%
  dplyr::mutate(FEEDER_predicted = as.factor(dplyr::case_when(feed_prob > prevalence ~ "YES",
                                                            feed_prob < prevalence ~ "NO")))

suppressWarnings({trainmat<-caret::confusionMatrix(data = DATA_TRAIN$FEEDER_predicted, reference = DATA_TRAIN$FEEDER_observed, positive="YES")})

#### classification success of test data

PRED<-stats::predict(RF2,data=DATA_TEST, type = "response")
prevalence<-table(DATA_TEST$FEEDER)[2]/sum(table(DATA_TEST$FEEDER))
DATA_TEST <- DATA_TEST %>%
  dplyr::mutate(FEEDER_observed = factor(FEEDER,levels=c("NO","YES"))) %>%
  dplyr::bind_cols(PRED$predictions) %>%
  dplyr::rename(no_feed_prob = NO, feed_prob = YES) %>%
  dplyr::mutate(FEEDER_predicted = as.factor(dplyr::case_when(feed_prob > prevalence ~ "YES",
                                                              feed_prob < prevalence ~ "NO")))

suppressWarnings({testmat<-caret::confusionMatrix(data = DATA_TEST$FEEDER_predicted, reference = DATA_TEST$FEEDER_observed, positive="YES")})

# hist(DATA_TEST$feed_prob, breaks=c(0,0.005,0.01,0.015,0.02,0.025,0.03,0.35,0.1,0.5,1))
# max(DATA_TEST$feed_prob)


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
    ggplot2::scale_x_discrete(name="",limit = mylevels,
                     labels = c("step length [m]",
                                "building present [y/n]",
                                "age [years]",                              
                                "distance to nest [m]",
                                "attendance pattern [visits/day]",
                                "time span of visitation [days]",
                                "frequency of visits",
                                "N of days",
                                "attendance duration [hrs]",
                                "N of repeat visits"))  +
    ggplot2::annotate("text",x=2,y=80,label=paste("Accuracy = ",round(testmat$byClass[3],3)),size=8) +
    ggplot2::theme(panel.background=ggplot2::element_rect(fill="white", colour="black"), 
                   axis.text.x=ggplot2::element_text(size=18, color="black"),
                   axis.text.y=ggplot2::element_text(size=16, color="black"), 
                   axis.title=ggplot2::element_text(size=20), 
                   panel.grid.major = ggplot2::element_blank(), 
                   panel.grid.minor = ggplot2::element_blank(), 
                   panel.border = ggplot2::element_blank())
  print(impplot)

#  ggsave("output/REKI_feed_ind_loc_variable_importance.jpg", height=7, width=11)
  
trainmat
testmat

# save.image("output/Feeding_analysis2024.RData")  
# load("output/Feeding_analysis2024.RData")
  
  
    
########### PARTIAL DEPENDENCE PLOTS ###############
#tic()
## this seems to take forever if parallel=T, otherwise 15 min
  
# ndpart<-partial(RF2, pred.var="dist_nest", trim.outliers=F, progress=T, prob=T)
# autoplot(ndpart)
# range(DATA_TRAIN$dist_nest)
# range(DATA_TEST$dist_nest)
# hist(OUT$dist_nest)
#   
# agepart<-partial(RF2, pred.var="age_cy", progress=T, trim.outliers=T, prob=T)  
# autoplot(agepart)    
   
  
  #  ggsave("output/REKI_feed_ind_loc_variable_importance.jpg", height=7, width=11)
 #toc() 
  
  ##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
  ########## VALIDATE PREDICTIONS WITH SURVEY DATA FROM EVA CEREGHETTI AND FIONA PELLET  #############
  ##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
  ## read in survey data from Eva Cereghetti
  ## Q1 is the question whether they feed or not
  setwd("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/REKIfeeding")
  feed_surveys<-fread("data/survey.final.csv") %>% #filter(Q1=="Ja") %>%
    mutate(FEEDER_surveyed=ifelse(Q1=="Ja",1,0)) %>%
    mutate(ID=paste0("Eva",nr)) %>%
    select(ID,coord_x,coord_y,square,random,area,building.type,FEEDER_surveyed) %>%
    st_as_sf(coords = c("coord_x", "coord_y"), crs=21781) %>%
    st_transform(crs = 3035) 
  
  ## read in survey from Fiona Pellet provided with addresses only
  ## Jerome Guelat provided R script to convert addresses to coordinates
  source("C:/Users/sop/OneDrive - Vogelwarte/General/ANALYSES/DataPrep/swisstopo_address_lookup.r")
  library(stringi)
  
  ## Feeding is the question whether they feed or not
  feed_surveys2<-read_csv("data/FeedersFionaPelle.csv", locale = locale(encoding = "UTF-8")) %>%
    #<-fread("data/FeedersFionaPelle.csv", encoding = 'UTF-8') %>% 
    mutate(FEEDER_surveyed=ifelse(Feeding=="Yes",1,0))
  
  ## generate coordinates from addresses
  feed_surveys2_locs<-swissgeocode(address=as.character(feed_surveys2$Address), nresults=3)
  
  feed_surv2_sf<-feed_surveys2_locs %>% rename(Address=address_origin) %>%
    left_join(feed_surveys2, by="Address",relationship ="many-to-many") %>%
    filter(!is.na(x)) %>%
    filter(!is.na(FEEDER_surveyed)) %>%
    mutate(ID=paste0("Fiona",ID)) %>%
    group_by(ID,lon,lat) %>%
    summarise(FEEDER_surveyed=max(FEEDER_surveyed)) %>%
    st_as_sf(coords = c("lon", "lat"), crs=4326) %>%
    st_transform(crs = 3035) %>%
    bind_rows(feed_surveys) %>% distinct()
  feed_surv2_sf$FEEDER_surveyed

## REMOVE DUPLICATE LOCATIONS
mtx_distance <- st_distance(feed_surv2_sf[feed_surv2_sf$FEEDER_surveyed==1,], feed_surv2_sf[feed_surv2_sf$FEEDER_surveyed==1,])
mtx_distance<-as_tibble(mtx_distance)
names(mtx_distance)<-feed_surv2_sf$ID[feed_surv2_sf$FEEDER_surveyed==1]
eliminate<-mtx_distance %>%
	mutate(ID=feed_surv2_sf$ID[feed_surv2_sf$FEEDER_surveyed==1]) %>%
	gather(key="ID2", value="dist",-ID) %>%
	mutate(dist=as.numeric(dist)) %>%
	filter(!(ID==ID2)) %>%
	filter(dist<40) %>%
	arrange(desc(dist))
dim(feed_surv2_sf)
for (l in 1:dim(eliminate)[1]) {
	#gridids<-length(unique(feed_surv2_sf$gridid[feed_surv2_sf$ID %in% c(eliminate$ID[l],eliminate$ID2[l])]))  ## check whether the two points are in the same grid cell
	pair1<-eliminate$ID2[l]
	if((eliminate$ID[l] %in% feed_surv2_sf$ID)){
		feed_surv2_sf<-feed_surv2_sf %>% filter(!(ID==pair1))
		}
}
dim(feed_surv2_sf)

## create 50 m buffer around feeders  
VAL_FEED_BUFF <-feed_surv2_sf %>%
  st_buffer(dist=50) %>%
  select(FEEDER_surveyed)

VAL_DAT<-DATA_TEST  %>%
  st_as_sf(coords = c("long", "lat"), crs=4326) %>%
  st_transform(crs = 3035) %>%
  st_join(VAL_FEED_BUFF,
          join = st_intersects, 
          left = TRUE) %>%
  mutate(FEEDER_surveyed=ifelse(is.na(FEEDER_surveyed),0,FEEDER_surveyed)) %>%
  dplyr::mutate(FEEDER_surveyed = as.factor(dplyr::case_when(FEEDER_surveyed==1 ~ "YES",
                                                             FEEDER_surveyed==0 ~ "NO"))) %>%
  mutate(FEEDER_observed=as.factor(dplyr::if_else(as.character(FEEDER_surveyed)=="YES",FEEDER_surveyed,FEEDER_observed))) %>%
  mutate(feed_obs_num=as.numeric(FEEDER_observed)-1)
str(VAL_DAT)  
suppressWarnings({valmat<-caret::confusionMatrix(data = VAL_DAT$FEEDER_predicted, reference = VAL_DAT$FEEDER_observed, positive="YES")})
valmat
5829/(5829+1025)
30973/(30973+ 247031)

  
## we cannot predict correct ABSENCE of feeding locations - even if one household does not feed their neighbours may and the prediction is therefore useless (and falsifying accuracy)
## but we use ROC curve to identify threshold
ROC_val<-pROC::roc(data=VAL_DAT,
             response=feed_obs_num,
             predictor=feed_prob)
AUC_VAL<-pROC::auc(ROC_val)

  

### DEFINE PREDICTION THRESHOLD FOR FEEDING LOCATIONS ###
THRESH_pts<-pROC::coords(ROC_val, "best", "threshold")$threshold


### CALCULATE BOYCE INDEX FOR VALIDATION DATA ###
BI<-Boyce(obs = VAL_DAT$feed_obs_num, pred = VAL_DAT$feed_prob)$Boyce
BI




# ##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
# ########## PLOT THE MAP FOR KNOWN AND OBSERVED FEEDING STATIONS   #############
# ##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
#   
#  
# head(OUT)    
#   
# # create plotting frame
#   
plot_OUT<-OUT %>%
  filter(feed_prob>THRESH_pts) %>%
  #filter(FEEDER_predicted=="YES") %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326)



########## CREATE A LEAFLET MAP TO SHOW APPROACH OF OVERLAYING FEEDING AND NON-FEEDING LOCATIONS ########################

## need to specify color palette 
# If you want to set your own colors manually:
pred.pal <- colorNumeric(c("cornflowerblue","firebrick"), seq(0,1))
m1 <- leaflet(options = leafletOptions(zoomControl = F)) %>% #changes position of zoom symbol
  setView(lng = mean(st_coordinates(plot_feeders)[,1]), lat = mean(st_coordinates(plot_feeders)[,2]), zoom = 11) %>%
  htmlwidgets::onRender("function(el, x) {L.control.zoom({ 
                           position: 'bottomright' }).addTo(this)}"
  ) %>% #Esri.WorldTopoMap #Stamen.Terrain #OpenTopoMap #Esri.WorldImagery
  addProviderTiles("Esri.WorldImagery", group = "Satellite",
                   options = providerTileOptions(opacity = 0.4, attribution = F,minZoom = 5, maxZoom = 20)) %>%
  
  addCircleMarkers(
    data=track_nofor_day_build,
    radius = 2,
    stroke = TRUE, color = "grey15", weight = 0.6,
    fillColor = "grey15", fillOpacity = 0.6
  ) %>%
  addCircleMarkers(
    data=plot_OUT,
    radius = 5,
    stroke = TRUE, color = ~pred.pal(feed_prob), weight = 0.8,
    fillColor = ~pred.pal(feed_prob), fillOpacity = 0.8,
    popup = ~ paste0("year ID: ", plot_OUT$year_id, "<br>", plot_OUT$timestamp)
  ) %>%
  addPolygons(
    data=grid %>%
      st_transform(4326),
    stroke = TRUE, color = "black", weight = 1,
    fillColor = NA, fillOpacity = 0
  ) %>%

  addLegend(     # legend for predicted prob of feeding
    position = "topleft",
    pal = pred.pal,
    values = plot_OUT$feed_prob,
    opacity = 1,
    title = "Probability of </br>anthropogenic feeding"
  ) %>%

  addScaleBar(position = "bottomright", options = scaleBarOptions(imperial = F))

m1



##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
########## ASSESS REASONS FOR FAILED PREDICTIONS (FIG S1) #############
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
head(OUT)
### remove outliers to make plot more informative
## https://stackoverflow.com/questions/59140960/remove-outliers-and-reduce-ylim-appropriately-for-each-facet-in-ggplot2

# filtering function - turns outliers into NAs to be removed
filter_lims <- function(x){
  l <- boxplot.stats(x)$stats[1]
  u <- boxplot.stats(x)$stats[5]
  
  for (i in 1:length(x)){
    x[i] <- ifelse(x[i]>l & x[i]<u, x[i], NA)
  }
  return(x)
}


# New facet label names for predictor variable
var.labs <- c("Step length", "Distance to nest", "Mean visit frequency","N days", "Time at location", "N revisits", "Evenness of attendance", "Time span of attendance")


DATA_TEST %>%
  mutate(Classification=ifelse(FEEDER_predicted==FEEDER_observed,"correct",
                               ifelse(FEEDER_predicted=="YES","false pred","missed pred"))) %>%
  # mutate(Classification=ifelse(Classification!="correct",Classification,
  #                              ifelse(FEEDER_predicted=="YES","correct pres","correct abs"))) %>%
  select(Classification,step_length,speed,revisits,residence_time,meanFreqVisit,n_days,TimeSpan,TempEven,dist_nest) %>%
  gather(key=variable, value=value,-Classification) %>%
  dplyr::mutate(variable=forcats::fct_relevel(variable,mylevels[c(1,4,7,8,9,10,6,5)])) %>%
  #group_by(variable) %>%
  mutate(value2 = filter_lims(value)) %>%  # new variable (value2) so as not to displace first one)
  mutate(value2 = filter_lims(value2)) %>%  # new variable (value2) so as not to displace first one)
  filter(!is.na(value2)) %>%
  
  ggplot(aes(x=Classification, y=value2)) +
  geom_boxplot(na.rm = TRUE) +
  facet_wrap(~variable, ncol=3, scales="free_y") +
  #facet_wrap(~variable, ncol=3, scales="free_y",labeller = labeller(variable = var.labs)) +
  labs(x="Prediction of locations near known anthropogenic feeding",
       y="Value of respective variable") +
  ggplot2::theme(panel.background=ggplot2::element_rect(fill="white", colour="black"), 
                 axis.text=ggplot2::element_text(size=16, color="black"),
                 strip.background=ggplot2::element_rect(fill="white", colour="black"), 
                 strip.text=ggplot2::element_text(size=16, color="black"),
                 axis.title=ggplot2::element_text(size=20), 
                 panel.grid.major = ggplot2::element_blank(), 
                 panel.grid.minor = ggplot2::element_blank(), 
                 panel.border = ggplot2::element_blank())

ggsave("C:/Users/sop/OneDrive - Vogelwarte/General/MANUSCRIPTS/AnthropFeeding/Figure_S1_2024.jpg", width=12, height=11, dpi=300)



pred_succ_var_summary<-OUT %>%
  mutate(Classification=ifelse(FEEDER_predicted==FEEDER_observed,"correct",ifelse(FEEDER_predicted=="YES","false pred","missed pred"))) %>%
  select(Classification,step_length,speed,revisits,residence_time,meanFreqVisit,n_days,TimeSpan,TempEven,dist_nest) %>%
  gather(key=variable, value=value,-Classification) %>%
  group_by(variable, Classification) %>%
  summarise(mean=median(value), min=min(value), max=max(value))

fwrite(pred_succ_var_summary,"output/RF1_pred_succ_var_summary_2024.csv")



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
tab2 <- st_intersects(grid, OUT_sf)
countgrid$N_feed_points<-lengths(tab2)


#### FIRST COUNT N INDIVIDUALS PER GRID CELL

for(c in 1:length(tab)){
  #countgrid$N_ind[c]<-length(unique(track_sf$year_id[tab[c][[1]]]))
  countgrid$N_ind[c]<-length(unique(track_sf$bird_id[tab[c][[1]]]))
}

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
  st_drop_geometry() %>%
  mutate(FEEDER=ifelse(FEEDER==0,0,1))

# reporting for manuscript
table(PRED_GRID$N_feed_points)[1]
dim(PRED_GRID)
summary(PRED_GRID)

# tmap_mode("view")
# tm_basemap(server="OpenStreetMap") +
#   tm_shape(OUT_sf)+
#   tm_symbols(col = 'feed_prob', size = 0.01) +
#   tm_shape(plot_feeders)+
#   tm_symbols(col = 'green', size = 0.05) +
#   tm_shape(countgrid)  +
#   tm_polygons(col = "firebrick", fill="N_feed_points", alpha=0.2)



##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
########## FIT SECOND RANDOM FOREST MODEL TO PREDICT FEEDING SITE AT GRID CELL LEVEL   #############
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################

RF3 <- ranger::ranger(FEEDER ~ n+N_ind+N_feed_points+N_feed_ind+prop_feed+prop_pts,     
                      data=PRED_GRID, mtry=2, num.trees=2500, replace=F, importance="permutation", oob.error=T, write.forest=T, probability=T)
#saveRDS(RF3, "output/feed_grid_RF_model.rds")
IMP2<-as.data.frame(RF3$variable.importance) %>%
  dplyr::mutate(variable=names(RF3$variable.importance)) %>%
  dplyr::rename(red.accuracy=`RF3$variable.importance`) %>%
  dplyr::arrange(dplyr::desc(red.accuracy)) %>%
  dplyr::mutate(rel.imp=(red.accuracy/max(red.accuracy))*100) %>%
  dplyr::select(variable,red.accuracy,rel.imp)
#write.table(IMP2,"clipboard", sep="\t")

#### classification success of training data

PRED<-stats::predict(RF3,data=PRED_GRID, type = "response")
PRED_GRID <- PRED_GRID %>%
  dplyr::mutate(FEEDER_observed = FEEDER) %>%
  dplyr::mutate(FEEDER_predicted=PRED$predictions[,2])
dim(PRED$predictions)
dim(PRED_GRID)

# ROC_train<-pROC::roc(data=PRED_GRID,response=FEEDER_observed,predictor=FEEDER_predicted)
# AUC<-pROC::auc(ROC_train)
# AUC
# THRESH<-pROC::coords(ROC_train, "best", "threshold")$threshold
THRESH<-table(PRED_GRID$FEEDER)[2]/dim(PRED_GRID)[1]


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
  ggplot2::scale_x_discrete(name="",limit = mylevels2,
                            labels = c("N individuals",
                                       "N feeding individuals",
                                       "N feeding locations",
                                       "Total N of GPS locations",
                                       "Proportion of feeding locations",
                                       "Proportion of feeding individuals"))  +
  #ggplot2::annotate("text",x=2,y=80,label=paste("AUC = ",round(AUC,3)),size=8) +
  ggplot2::theme(panel.background=ggplot2::element_rect(fill="white", colour="black"), 
                 axis.text.x=ggplot2::element_text(size=18, color="black"),
                 axis.text.y=ggplot2::element_text(size=16, color="black"), 
                 axis.title=ggplot2::element_text(size=20), 
                 panel.grid.major = ggplot2::element_blank(), 
                 panel.grid.minor = ggplot2::element_blank(), 
                 panel.border = ggplot2::element_blank())
print(impplot2)

#ggsave("output/REKI_feed_congregation_variable_importance2024.jpg", height=7, width=11)



#### BASIC STATISTICS OF PREDICTIONS
dim(PRED_GRID %>% filter(FEEDER_observed==1  & n>10))
dim(PRED_GRID %>% filter(FEEDER_observed==1 & FEEDER_predicted>0.5))
dim(PRED_GRID %>% filter(FEEDER_observed==1 & n>10 & FEEDER_predicted>mean(PRED_GRID$FEEDER_observed)))/dim(PRED_GRID %>% filter(FEEDER_observed==1  & n>10))
hist(PRED_GRID$FEEDER_predicted)
mean(PRED_GRID$FEEDER_observed)
range(PRED_GRID$prop_feed)
range(PRED_GRID$prop_pts)
names(PRED_GRID)



########## CREATE OUTPUT GRID WITH PREDICTED FEEDING LOCATIONS ########################

OUTgrid<-countgrid %>% select(-FEEDER) %>%
  mutate(gridid=seq_along(n)) %>%
  left_join(PRED_GRID, by=c("gridid","n","N_ind","N_feed_points","N_feed_ind","prop_feed","prop_pts")) %>%
  filter(!is.na(FEEDER_predicted))



########### PARTIAL DEPENDENCE PLOTS ###############
tic()
## this seems to take forever if parallel=T, otherwise 15 min

ndpart<-partial(RF3, pred.var="N_feed_points", trim.outliers=F, progress=T, prob=T)
autoplot(ndpart)
range(PRED_GRID$N_feed_points)

agepart<-partial(RF3, pred.var="prop_feed", progress=T, trim.outliers=T, prob=T)
autoplot(agepart)

ptspart<-partial(RF3, pred.var="prop_pts", progress=T, trim.outliers=T, prob=T)
autoplot(ptspart)

toc()


##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
########## VALIDATE PREDICTIONS WITH SURVEY DATA FROM EVA CEREGHETTI AND FIONA PELLET  #############
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
### validation data already read in above
feed_surv2_sf$FEEDER_surveyed

#validat <- st_intersection(OUTgrid, feed_surv2_sf) %>%
validat <- st_intersection(feed_surv2_sf,OUTgrid) %>%
  filter(!is.na(n))  ## excludes bullshit addresses outside of study area
validat$FEEDER_surveyed
summary(validat$FEEDER_predicted)

## we cannot predict correct ABSENCE of feeding locations - even if one household does not feed their neighbours may and the prediction is therefore useless (and falsifying accuracy)
## but you cannot fit an ROC curve if the response is only 1 and not 0
# validat<-validat %>% filter(FEEDER_surveyed==1)
# ROC_val<-pROC::roc(data=validat,response=FEEDER_surveyed,predictor=FEEDER_predicted)
# AUC_TEST<-pROC::auc(ROC_val)
# pROC::coords(ROC_val, "best", "threshold")



##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
########## SUMMARISE VALIDATION DATA  #############
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
validat

### find and remove duplicates
# duplfeed<-validat %>% 
#   st_drop_geometry() %>%
#   distinct() %>%
#   group_by(n,N_ind,N_feed_points,N_feed_ind,prop_feed,prop_pts,gridid) %>%
#   summarise(nv=length(FEEDER_surveyed)) %>% filter(nv>1)


VAL_DAT2<-validat %>% #st_drop_geometry() %>%
  select(-square,-random,-area,-building.type) %>% distinct() %>% #filter(gridid==890) 
  #filter(!(is.na(nr) & (gridid %in% duplfeed$gridid))) %>%
  filter(FEEDER_surveyed==1) %>%
  select(ID, n,N_ind,N_feed_points,N_feed_ind,prop_feed,prop_pts,FEEDER_predicted) %>%
  mutate(Classification=ifelse(FEEDER_predicted>THRESH,"correct","missed"))

### summarise the predicted sites
table(VAL_DAT2$Classification)
summary(VAL_DAT2$FEEDER_predicted)
mean(VAL_DAT2$FEEDER_predicted)
table(VAL_DAT2$Classification)[1]/dim(VAL_DAT2)[1]
min(VAL_DAT2$n[VAL_DAT2$Classification=="correct"])
min(VAL_DAT2$N_ind[VAL_DAT2$Classification=="correct"])
min(VAL_DAT2$N_feed_points[VAL_DAT2$Classification=="correct"])
min(VAL_DAT2$N_feed_ind[VAL_DAT2$Classification=="correct"])
VAL_DAT2 %>% filter(N_feed_points==0)

### CALCULATE BOYCE INDEX FOR VALIDATION DATA ###
BI2<-Boyce(obs = validat$FEEDER_surveyed , pred = validat$FEEDER_predicted, n.bins=10)$Boyce
BI2



##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
########## PLOT THE MAP FOR KNOWN AND OBSERVED FEEDING STATIONS   #############
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################

########## CREATE A LEAFLET MAP OF PREDICTED FEEDING LOCATIONS ########################

## need to specify color palette 
# If you want to set your own colors manually:
pred.pal <- colorNumeric(c("cornflowerblue","firebrick"), seq(0,1))
feed.pal <- colorFactor(c("darkgreen","lightgreen"), unique(plot_feeders$Type))
val.pal <- colorFactor(c("green","red"), unique(VAL_DAT$Classification))
m2 <- leaflet(options = leafletOptions(zoomControl = F)) %>% #changes position of zoom symbol
  setView(lng = mean(st_coordinates(plot_feeders)[,1]), lat = mean(st_coordinates(plot_feeders)[,2]), zoom = 11) %>%
  htmlwidgets::onRender("function(el, x) {L.control.zoom({ 
                           position: 'bottomright' }).addTo(this)}"
  ) %>% #Esri.WorldTopoMap #Stamen.Terrain #OpenTopoMap #Esri.WorldImagery
  addProviderTiles("Esri.WorldImagery", group = "Satellite",
                   options = providerTileOptions(opacity = 0.6, attribution = F,minZoom = 5, maxZoom = 22)) %>%
  addProviderTiles("OpenTopoMap", group = "Roadmap", options = providerTileOptions(attribution = F,minZoom = 5, maxZoom = 22)) %>%  
  addLayersControl(baseGroups = c("Satellite", "Roadmap")) %>%  
  
  # addCircleMarkers(
  #   data=plot_OUT,
  #   radius = 2,
  #   stroke = TRUE, color = "black", weight = 0.5,
  #   fillColor = "grey75", fillOpacity = 0.5,
  #   popup = ~ paste0("year ID: ", plot_OUT$year_id, "<br>", plot_OUT$timestamp)
  # ) %>%
  addPolygons(
    data=OUTgrid %>%
      st_transform(4326),
    stroke = TRUE, color = ~pred.pal(FEEDER_predicted), weight = 1,
    fillColor = ~pred.pal(FEEDER_predicted), fillOpacity = 0.5,
    popup = ~as.character(paste("N_pts=",n,"/ N_ind=",N_ind,"/ Prop feed pts=",round(prop_feed,3), sep=" ")),
    label = ~as.character(round(FEEDER_predicted,3))
  ) %>%
  
  addCircleMarkers(
    data=plot_feeders,
    radius = 2,
    stroke = TRUE, color = "black", weight = 1,   ###~feed.pal(Type)
    fillColor = "grey25", fillOpacity = 0.5   ## ~feed.pal(Type)
  ) %>%
  
  # addCircleMarkers(
  #   data=VAL_DAT2 %>%
  #     st_transform(4326) %>% filter(N_feed_points==0),
  #   radius = 8,
  #   stroke = TRUE, color = "goldenrod", weight = 1,
  #   fillColor = "goldenrod"
  # ) %>%
  

  addCircleMarkers(
    data=VAL_DAT2 %>%
      st_transform(4326),
    radius = 5,
    stroke = TRUE, color = ~val.pal(Classification), weight = 1,
    fillColor = ~val.pal(Classification), fillOpacity = 1,
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
  # addLegend(     # legend for known feeding sites
  #   position = "topleft",
  #   pal = feed.pal,
  #   values = plot_feeders$Type,
  #   opacity = 1,
  #   title = "Type of feeding station"
  # ) %>%
  addLegend(     # legend for known feeding sites
    position = "topleft",
    pal = val.pal,
    values = VAL_DAT$Classification,
    opacity = 1,
    title = "Validation (interviews)"
  ) %>%
  
  
  
  addScaleBar(position = "bottomright", options = scaleBarOptions(imperial = F))

m2


htmltools::save_html(html = m2, file = "C:/Users/sop/OneDrive - Vogelwarte/General/MANUSCRIPTS/AnthropFeeding/Figure_1_2024.html")
mapview::mapshot(m2, url = "C:/Users/sop/OneDrive - Vogelwarte/General/MANUSCRIPTS/AnthropFeeding/Figure_1_2024.html")
st_write(OUTgrid,"output/REKI_predicted_anthropogenic_feeding_areas_2024.kml",append=FALSE)



##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
########## COMPARE DATA BETWEEN SUCCESSFUL AND FAILED PREDICTIONS   #############
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
### remove outliers to make plot more informative
## https://stackoverflow.com/questions/59140960/remove-outliers-and-reduce-ylim-appropriately-for-each-facet-in-ggplot2

# filtering function - turns outliers into NAs to be removed
filter_lims <- function(x){
  l <- boxplot.stats(x)$stats[1]
  u <- boxplot.stats(x)$stats[5]
  
  for (i in 1:length(x)){
    x[i] <- ifelse(x[i]>l & x[i]<u, x[i], NA)
  }
  return(x)
}


# New facet label names for predictor variable
var.labs <- c("Total N of GPS locations","N individuals","N feeding locations","N feeding individuals","Proportion of feeding individuals",
              "Proportion of feeding locations")
names(var.labs) <- names(VAL_DAT2)[2:7]


VAL_DAT2 %>%
  st_drop_geometry() %>%
  select(-FEEDER_predicted,-ID) %>%
  gather(key=variable, value=value,-Classification) %>% 
  filter(!(variable=="n" & value>45000)) %>%  ### remove a single outlier value
  filter(!(variable=="N_feed_ind" & value>30)) %>%  ### remove a single outlier value
  filter(!(variable=="N_feed_points" & value>350)) %>%  ### remove a single outlier value
  filter(!(variable=="prop_pts" & value>0.3)) %>%  #filter(variable=='n') %>% arrange(desc(value))### remove a single outlier value
  dplyr::mutate(variable=forcats::fct_relevel(variable,mylevels2)) %>%
  #mutate(value2 = filter_lims(value)) %>%  # new variable (value2) so as not to displace first one)
  
  ggplot(aes(x=Classification, y=value)) +
  geom_boxplot(na.rm = TRUE, coef = 5) +
  facet_wrap(~variable, ncol=3, scales="free_y",labeller = labeller(variable = var.labs)) +
  labs(x="Prediction of grid cells with known anthropogenic feeding",
       y="Value of respective variable") +
  ggplot2::theme(panel.background=ggplot2::element_rect(fill="white", colour="black"), 
                 axis.text=ggplot2::element_text(size=16, color="black"),
                 strip.background=ggplot2::element_rect(fill="white", colour="black"), 
                 strip.text=ggplot2::element_text(size=16, color="black"),
                 axis.title=ggplot2::element_text(size=20), 
                 panel.grid.major = ggplot2::element_blank(), 
                 panel.grid.minor = ggplot2::element_blank(), 
                 panel.border = ggplot2::element_blank())
  

ggsave("C:/Users/sop/OneDrive - Vogelwarte/General/MANUSCRIPTS/AnthropFeeding/Figure_S2_2024.jpg", width=10, height=8, dpi=300)






pred2_succ_var_summary<-VAL_DAT2 %>%
  st_drop_geometry() %>%
  select(-FEEDER_predicted,-ID) %>%
  gather(key=variable, value=value,-Classification) %>%
  filter(!(variable=="n" & value>45000)) %>%  ### remove a single outlier value
  # filter(!(variable=="N_feed_ind" & value>30)) %>%  ### remove a single outlier value
  # filter(!(variable=="N_feed_points" & value>500)) %>%  ### remove a single outlier value
  group_by(variable, Classification) %>%
  summarise(mean=median(value), min=min(value), max=max(value))

fwrite(pred2_succ_var_summary,"output/RF2_pred_succ_var_summary2024.csv")



##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
########## COMPILE THE OUTPUT REPORT   #############
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
library(knitr)
library(markdown)
library(rmarkdown)

### create HTML report for overall summary report
# Sys.setenv(RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools")
#
# rmarkdown::render('C:\\Users\\sop\\OneDrive - Vogelwarte\\REKI\\Analysis\\Feeding\\Feeding_Analysis_report2023_v2.Rmd',
#                   output_file = "Feeding_Analysis_Progress2023.html",
#                   output_dir = 'C:\\Users\\sop\\OneDrive - Vogelwarte\\REKI\\Analysis\\Feeding\\output')


save.image("output/Feeding_analysis2024.RData")  
load("output/Feeding_analysis2024.RData")

