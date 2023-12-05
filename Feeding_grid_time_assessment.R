###########################################################################################################################
###### ANALYSIS TO QUANTIFY USAGE OF ANTHROPOGENIC FEEDING SITES OF RED KITES ACROSS SWITZERLAND ################
###########################################################################################################################
# original idea by Nathalie Heiniger (MSc thesis 2020)
# uses output provided by "Feeding_projection.r"
# created by steffen.oppel@vogelwarte.ch in November 2023
# inspired by REKI Team meeting on 14 Nov 2023 - include metrics of usage in manuscript

library(tidyverse)
library(dplyr, warn.conflicts = FALSE)
options(dplyr.summarise.inform = FALSE)
library(sf)
library(lubridate)
library(adehabitatHR)
library(stars)
library(trip)
library(sqldf)
library(amt)
library(adehabitatLT)
filter<-dplyr::filter
sf_use_s2(FALSE)



########## READ IN RAW TRACKING DATA ########################
tracks<-readRDS(file = "C:/Users/sop/OneDrive - Vogelwarte/General/DATA/REKI_raw_tracking_data_30Nov2023.rds")



# DATA PREPARATION -------------------------------------------------------------


tracks <- tracks %>%
  filter(!is.na(long)) %>%
  filter(!is.na(lat)) %>%
  filter(long<11) %>%
  filter(long>-10) %>%
  filter(lat<49) %>%
  filter(lat>35) %>%
  filter(!is.na(timestamp)) %>%
  mutate(year=year(timestamp), 
         year_id=paste(year,id,sep="_")) %>%
  rename(bird_id = id) %>%
  filter(!is.na(year_id)) %>%
  mutate(season1 = ifelse(yday(timestamp)>70 ,"B","NB")) %>%
  mutate(season2 = ifelse(yday(timestamp)>175 ,"NB","B")) %>%
  mutate(season_id = ifelse(yday(timestamp)<70 ,
                              paste(season1,year-1,bird_id,sep="_"),
                              paste(season2,year,bird_id,sep="_"))
    ) %>%
  select(-season1,-season2) %>%
  arrange(season_id,timestamp)
dim(tracks)
head(tracks)


##### eliminate bullshit birds with <100 locations
exclude<-tracks %>% group_by(year_id) %>%
  summarise(n = length(unique(yday(timestamp)))) %>%
  filter(n<10)
tracks<-tracks %>% filter(!(year_id %in% exclude$year_id))
dim(tracks)



  
### REMOVE DUPLICATE TIME STAMPS AND THOSE WITH TOO FEW LOCATIONS ####
dupes<-duplicated(tracks,by=c("timestamp","year_id","year","season_id","bird_id"))
which(dupes==TRUE)
tracks<-tracks[!dupes==T,]
dim(tracks)





########### clean up workspace #######################
ls()
rm('data','exclude','dupes','trackingdata','CHgrid')
gc()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# INTERPOLATE ALL TRACKING DATA TO 15 MIN INTERVALS TO ENSURE THAT NO MIGRATION IS MISSED THROUGH A CELL
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Convert to LTRAJ TO INTERPOLATE DATA
traj<-as.ltraj(xy=tracks[,3:4],date=tracks[,2],id=tracks[,7],infolocs=tracks[,c(1:7)],typeII=TRUE)

## Rediscretization every 900 seconds
tr <- redisltraj(traj, 900, type="time")

## Convert output into a data frame
all_dat15min<-data.frame()
for (l in 1:length(unique(tracks$season_id))){
  out<-tr[[l]]
  out$season_id<-as.character(attributes(tr[[l]])[4])				#### extracts the MigID from the attribute 'id'
  all_dat15min<-rbind(all_dat15min,out)				#### combines all data
}

setwd("C:/Users/sop/OneDrive - Vogelwarte/General/DATA")
saveRDS(all_dat15min,"REKI_regular_15min_tracking_data_30Nov2023.rds")

### remove unnecessary data
rm('tracks','tr','traj','out')
gc()



# converting to metric CRS prior to curtailing data to extent of CHgrid
track_sf <- all_dat15min %>% 
  select(x,y,date,season_id) %>%
  st_as_sf(coords = c("x", "y"), crs = 4326) %>%
  st_transform(crs = 3035)
head(track_sf)
dim(track_sf)

setwd("C:/Users/sop/OneDrive - Vogelwarte/General/DATA")
saveRDS(track_sf,"REKI_regular_15min_tracking_data_30Nov2023_projected.rds")





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# JOIN TRACKING DATA WITH PREDICTION GRID AND SUMMARISE HOURS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm('all_dat15min')
gc()

###### LOADING DATA THAT PREDICTED FEEDING PROBABILITY ------------------------------------------------------------
## set root folder for project
try(setwd("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/REKIfeeding"),silent=T)
# try(setwd("C:/STEFFEN/OneDrive - Vogelwarte/REKI/Analysis/REKIfeeding"),silent=T)
track_sf<-readRDS("C:/Users/sop/OneDrive - Vogelwarte/General/DATA/REKI_regular_15min_tracking_data_30Nov2023_projected.rds")
# track_sf<-readRDS("C:/STEFFEN/OneDrive - Vogelwarte/General/DATA/REKI_regular_15min_tracking_data_30Nov2023_projected.rds")
CHgrid<-readRDS("output/REKI_feeding_grid.rds")
range(CHgrid$FEEDER_predicted)


FEED_TRACK<-track_sf %>% st_join(CHgrid) %>%
  mutate(FEEDER_predicted=ifelse(is.na(FEEDER_predicted),0,FEEDER_predicted)) %>%
  separate(season_id, into=c("season","year","bird_id"), sep="_")

rm('track_sf')
gc()

SUMMARY<-FEED_TRACK  %>%
  mutate(season=ifelse(substr(season_id,1,1)=="B","breeding","non-breeding")) %>%
  group_by(season_id, season) %>%
  summarise(FEED=mean(FEEDER_predicted))


ggplot(SUMMARY, aes(x=season,y=FEED)) + geom_boxplot()





###### LOADING INDIVIDUAL DATA THAT PREDICTED FEEDING PROBABILITY ------------------------------------------------------------


### TO DO #############
## link foraging area use to individuals
## create 4-panel plot with 
## x-axis: age
## colour= sex
## panels: breeding vs non-breeding season
## ask Steph to exclude migrants OR use simple latitudinal cutoff

inddat<-fread("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/NestTool/REKI/data/Individual_life_history_2015-2022.csv") %>%
  select(bird_id,sex_compiled,hatch_year, tag_year) %>%
  rename(sex=sex_compiled) %>%
  mutate(hatch_year= as.numeric(ifelse(hatch_year=="unknown",(tag_year-4),hatch_year)))

# inddat<-fread("C:/STEFFEN/OneDrive - Vogelwarte/REKI/Analysis/NestTool/REKI/data/Individual_life_history_2015-2022.csv") %>%
#   select(bird_id,sex_compiled,hatch_year, tag_year) %>%
#   rename(sex=sex_compiled) %>%
#   mutate(hatch_year= as.numeric(ifelse(hatch_year=="unknown",(tag_year-4),hatch_year)))

### MERGE tracks with INDIVIDUAL DATA

FEED_DATA <- as_tibble(FEED_TRACK) %>% #mutate(year_id = substr(season_id,3,nchar(season_id))) %>%
  mutate(bird_id = as.integer(bird_id), year=as.integer(year)) %>%
  left_join(inddat, by='bird_id') %>%
  select(-geometry,-n, -N_ind, -N_feed_points, -N_feed_ind, -prop_feed, -prop_pts, -gridid) %>%
  mutate(season_id = paste(season,year,bird_id,sep="_")) %>%
  #mutate(year=year(date)) %>%
  mutate(age_cy=(year-tag_year)+1)
  # mutate(HR=ifelse(home_range_id>0,"settled","not settled")) %>%
  # mutate(BR=ifelse(nest_id>0,"breeding","not breeding"))
  
saveRDS(FEED_DATA,"output/REKI_feed_data.rds")


########### clean up workspace #######################
rm('FEED_TRACK')
gc()


### SUMMARISE BY SEX AGE AND BREEDING STATUS

INDSUMMARY<-FEED_DATA %>%
  #select(-geometry,-n, -N_ind, -N_feed_points, -N_feed_ind, -prop_feed, -prop_pts, -gridid) %>%
  mutate(FEED_hrs=FEEDER_predicted*0.25) %>%
  mutate(DAY=yday(date)) %>%
  group_by(season_id, season, sex, age_cy,DAY) %>%
  summarise(FEED=sum(FEED_hrs)) %>%
  ungroup() %>%
  group_by(season_id, season, sex, age_cy) %>%
  summarise(FEED=mean(FEED, na.rm=T))
INDSUMMARY %>% filter(is.na(sex))

PLOTDAT<- INDSUMMARY %>%
  filter(sex %in% c("m","f")) %>% ### remove the unassigned and unknown sexes (3 birds)
  ungroup() %>%
  group_by(season, sex, age_cy) %>%
  summarise(n=length(unique(season_id)), median=quantile(FEED,0.5),lcl=quantile(FEED,0.15),ucl=quantile(FEED,0.85))


### CREATE A GRAPH OF USAGE BY AGE, SEX, and BREEDING SEASON
range(PLOTDAT$ucl)
PLOTDAT %>%  
  mutate(age_cy=ifelse(sex=="m", age_cy-0.2, age_cy+0.2)) %>%
  mutate(season=ifelse(season=="B", "breeding season (Mar - Jul)", "non-breeding season (Aug - Feb)")) %>%
  
  ggplot() +
  geom_point(aes(y=median, x=age_cy, colour=sex),size=1.5)+
  geom_errorbar(aes(x=age_cy, ymin=lcl, ymax=ucl, colour=sex), width=0.1)+
  scale_y_continuous(name="Mean daily usage (hrs) of anthropogenic feeding", limits=c(0,24), breaks=seq(0,24,3)) +
  scale_x_continuous(name="Age (in years)", limits=c(0.5,8.5), breaks=seq(1,8,1)) +
  facet_wrap(~season, ncol=1, scales="free_y", dir="v") +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text.y=element_text(size=16, color="black"), 
        axis.text.x=element_text(size=13, color="black"), 
        axis.title=element_text(size=18), 
        strip.text=element_text(size=18, color="black"),
        legend.text=element_text(size=14, color="black"),
        legend.title=element_text(size=18, color="black"),
        legend.key=element_blank(),
        legend.position=c(0.08,0.92),
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())

ggsave("output/REKI_feeding_site_usage.jpg", height=9, width=11)






# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # ABANDONED BECAUSE TIME IN GRID WAS VERY HARD
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# # Points in SWISS GRID
# track_sf_CH <- st_filter(track_sf,CHgrid)
# duplicated(track_sf_CH)
# dim(track_sf_CH)
# 
# 
# length(unique(track_sf_CH$season_id))
# toofew<-track_sf_CH %>% mutate(count=1) %>% group_by(season_id) %>%
#   summarise(N=sum(count)) %>%
#   filter(N<10)
# track_sf_CH<-track_sf_CH %>% filter(!(season_id %in% toofew$season_id))
# dim(track_sf_CH)
# 
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # QUANTIFY TIME IN AREA OF EACH SWISS GRID CELL - THIS CALCULATES THE OVERALL POPULATION LEVEL
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CHspdf<-as(CHgrid, "Spatial")
# grdSUI<-makeGridTopology(obj=CHspdf)
# REKI_trips<-trip(track_sf_CH, TORnames=c("timestamp","season_id"), correct_all=TRUE)		### switch to "DateTime" when using the raw locations
# TIME_GRID <- tripGrid(REKI_trips, grid=grdSUI,method="pixellate")						### this will provide the number of bird seconds spent in each grid cell
# spplot(TIME_GRID)			## plots the trips with a legend
# proj4string(turbSPDF)<-proj4string(EVSP_all)
# 
# 
# ### CONVERT SPATIAL GRID TO SOMETHING WE CAN PLOT
# spdf <- SpatialPixelsDataFrame(points=trg, data=trg@data)
# HOTSPOTS<-data.frame(lat=spdf@coords[,2],long=spdf@coords[,1],time=spdf@data$z)
# HOTSPOTS$time<-HOTSPOTS$time/(3600*24)								### this converts the seconds into bird days
# summary(HOTSPOTS)
# 
# HOTSPOTS<-HOTSPOTS[HOTSPOTS$time>10,]
# 
# ##### CONVERT TO SPATIAL POLYGONS FOR OVERLAY ###
# ras<- raster(spdf)		# converts the SpatialPixelDataFrame into a raster
# spoldf <- rasterToPolygons(ras, n=4) # converts the raster into quadratic polygons
# proj4string(spoldf)<-proj4string(EVSP_all)
# 
# 
# ### PRODUCE NICE AND SHINY MAP WITH THE MOST IMPORTANT TEMPORARY CONGREGATION SITES ###
# 
# #MAP <- get_map(EGVUbox, source="google", zoom=4, color = "bw")		### retrieves a map from Google (requires live internet connection)
# MAP <- get_map(location = c(lon = mean(long), lat = mean(lat)), source="google", zoom=4, color = "bw")		### retrieves a map from Google (requires live internet connection)
# 
# pdf("EGVU_MIGRATION_HOTSPOTS.pdf", width=12, height=11)
# ggmap(MAP)+geom_tile(data=HOTSPOTS, aes(x=long,y=lat, fill = time)) +
#   scale_fill_gradient(name = 'N bird days', low="white", high="red", na.value = 'transparent', guide = "colourbar", limits=c(10, 30))+
#   theme(axis.ticks = element_blank(),axis.text = element_blank(),axis.title = element_blank())+
#   theme(strip.text.y = element_text(size = 20, colour = "black"), strip.text.x = element_text(size = 15, colour = "black"))+
#   geom_point(data=turbs, aes(x=Longitude, y=Latitude), pch=16, col='darkolivegreen', size=0.5)
# dev.off()
# 
# 
# 
# 
# 
# 
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # TRYING TO USE RECURSE WITH AN OLD INSTALLATION OF RECURSE AND RGEOS
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# require(devtools)
# install_version("recurse", version = "0.x.x", repos = "http://cran.us.r-project.org")
# 
# 
# 
