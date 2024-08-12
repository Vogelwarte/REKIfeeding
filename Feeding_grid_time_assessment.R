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


### read in Switzerland map
# SUI<-st_read("S:/rasters/outline_maps/swiss_map_overview/layers.gpkg") %>% filter(country=="Switzerland")
# saveRDS(SUI,"data/Swiss_border.rds")
SUI<-readRDS("data/Swiss_border.rds")
plot(SUI)


########## READ IN RAW TRACKING DATA ########################
#tracks<-readRDS(file = "C:/Users/sop/OneDrive - Vogelwarte/General/DATA/REKI_raw_tracking_data_30Nov2023.rds")
tracks<-readRDS(file = "data/REKI_trackingdata_raw2024.rds") %>%
  st_transform(2056) %>%
  st_intersection(.,SUI) %>% filter(!is.na(country)) %>% ### remove all data outside of Switzerland
  st_transform(4326) %>%
  dplyr::mutate(long = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2]) %>% st_drop_geometry()


# DATA PREPARATION -------------------------------------------------------------

tracks <- tracks %>%
  filter(!is.na(long)) %>%
  filter(!is.na(lat)) %>%
  filter(!is.na(timestamp)) %>%
  mutate(year=year(timestamp), 
         year_id=paste(year,bird_id,sep="_")) %>%
  filter(!is.na(year_id)) %>%
  mutate(season1 = ifelse(yday(timestamp)>70 ,"B","NB")) %>%   ## 10 MARCH
  mutate(season2 = ifelse(yday(timestamp)>175 ,"NB","B")) %>%  ## 20 June
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
rm('data','exclude','dupes','trackingdata','CHgrid','SUI')
gc()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# INTERPOLATE ALL TRACKING DATA TO 15 MIN INTERVALS TO ENSURE THAT NO MIGRATION IS MISSED THROUGH A CELL
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Convert to LTRAJ TO INTERPOLATE DATA
trajtracks<-as.data.frame(tracks)  ### LTRAJ only works on data.frames but not on tibbles!
traj<-as.ltraj(xy=trajtracks[,7:8],date=trajtracks[,5],id=trajtracks[,11],infolocs=trajtracks[,c(1:11)],typeII=TRUE)

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
saveRDS(all_dat15min,"REKI_regular_15min_tracking_data_June2024.rds")

### remove unnecessary data
rm('tracks','tr','traj','out','trajtracks')
gc()



# converting to metric CRS prior to curtailing data to extent of CHgrid
track_sf <- all_dat15min %>% 
  select(x,y,date,season_id) %>%
  st_as_sf(coords = c("x", "y"), crs = 4326) %>%
  st_transform(crs = 3035)
head(track_sf)
dim(track_sf)

setwd("C:/Users/sop/OneDrive - Vogelwarte/General/DATA")
saveRDS(track_sf,"REKI_regular_15min_tracking_data_June2024_projected.rds")





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# JOIN TRACKING DATA WITH PREDICTION GRID AND SUMMARISE HOURS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm('all_dat15min')
gc()

###### LOADING DATA THAT PREDICTED FEEDING PROBABILITY ------------------------------------------------------------
## set root folder for project
try(setwd("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/REKIfeeding"),silent=T)
# try(setwd("C:/STEFFEN/OneDrive - Vogelwarte/REKI/Analysis/REKIfeeding"),silent=T)
track_sf<-readRDS("C:/Users/sop/OneDrive - Vogelwarte/General/DATA/REKI_regular_15min_tracking_data_13Jun2024_projected.rds")
# track_sf<-readRDS("C:/STEFFEN/OneDrive - Vogelwarte/General/DATA/REKI_regular_15min_tracking_data_13Jun2024_projected.rds")
CHgrid<-readRDS("output/REKI_feeding_grid.rds")
range(CHgrid$FEEDER_predicted)


FEED_TRACK<-track_sf %>% st_join(CHgrid) %>%
  mutate(FEEDER_predicted=ifelse(is.na(FEEDER_predicted),0,FEEDER_predicted)) %>%
  separate(season_id, into=c("season","year","bird_id"), sep="_")

rm('track_sf')
gc()

SUMMARY<-FEED_TRACK  %>%
  st_drop_geometry() %>%
  #mutate(season=ifelse(substr(season_id,1,1)=="B","breeding","non-breeding")) %>%
  mutate(season_id = ifelse(yday(date)<70 ,
                            paste(season,year(date)-1,bird_id,sep="_"),
                            paste(season,year(date),bird_id,sep="_"))) %>%
  group_by(season_id, season) %>%
  summarise(FEED=mean(FEEDER_predicted))


ggplot(SUMMARY, aes(x=season,y=FEED)) + geom_boxplot()





###### LOADING INDIVIDUAL DATA THAT PREDICTED FEEDING PROBABILITY ------------------------------------------------------------


library(readxl)
inddat<-read_excel("C:/Users/sop/OneDrive - Vogelwarte/General/DATA/Individual_life_history_2015-2023.xlsx", sheet="Individual_life_history_2015-20") %>% # updated on 3 June 2024 to include birds from 2022
  dplyr::select(bird_id,ring_number,tag_year,sex_compiled, age, hatch_year) %>%
  rename(ring_id=ring_number) %>%
  rename(sex=sex_compiled) %>%
  mutate(hatch_year=if_else(is.na(as.numeric(hatch_year)),tag_year-3,as.numeric(hatch_year)))

inddat<-read_excel("C:/STEFFEN/OneDrive - Vogelwarte/General/DATA/Individual_life_history_2015-2023.xlsx", sheet="Individual_life_history_2015-20") %>% # updated on 3 June 2024 to include birds from 2022
  dplyr::select(bird_id,ring_number,tag_year,sex_compiled, age, hatch_year) %>%
  rename(ring_id=ring_number) %>%
  rename(sex=sex_compiled) %>%
  mutate(hatch_year=if_else(is.na(as.numeric(hatch_year)),tag_year-3,as.numeric(hatch_year)))

### MERGE tracks with INDIVIDUAL DATA

FEED_DATA <- FEED_TRACK %>% #mutate(year_id = substr(season_id,3,nchar(season_id))) %>%
  st_drop_geometry() %>%
  mutate(bird_id = as.integer(bird_id), year=as.integer(year)) %>%
  left_join(inddat, by='bird_id') %>%
  select(-n, -N_ind, -N_feed_points, -N_feed_ind, -prop_feed, -prop_pts, -gridid) %>%
  mutate(season_id = paste(season,year,bird_id,sep="_")) %>%
  #mutate(year=year(date)) %>%
  mutate(age_cy=(year-tag_year)+1)
  # mutate(HR=ifelse(home_range_id>0,"settled","not settled")) %>%
  # mutate(BR=ifelse(nest_id>0,"breeding","not breeding"))
  
# saveRDS(FEED_DATA,"output/REKI_food_supplementation_index.rds")
# 
# 
# ########### clean up workspace #######################
# FEED_DATA <- readRDS("output/REKI_feed_data2024.rds")
rm('FEED_TRACK')
gc()


SEASONSUMMARY<-FEED_DATA %>%
  mutate(FEED_hrs=FEEDER_predicted*0.25) %>%
  mutate(DAY=yday(date)) %>%
  group_by(season_id, year, bird_id, season, sex, age_cy) %>%
  summarise(FEED=sum(FEED_hrs)) %>%
  ungroup()
saveRDS(SEASONSUMMARY,"C:\\Users\\sop\\OneDrive - Vogelwarte\\General\\ANALYSES\\REKIpopmod/data/CH_popmodel/REKI_food_supplementation_index.rds")


### SUMMARISE BY SEX AGE AND BREEDING STATUS

INDSUMMARY<-FEED_DATA %>%
  #select(-geometry,-n, -N_ind, -N_feed_points, -N_feed_ind, -prop_feed, -prop_pts, -gridid) %>%
  mutate(FEED_hrs=FEEDER_predicted*0.25) %>%
  mutate(DAY=yday(date)) %>%
  group_by(season_id, year, bird_id, season, sex, age_cy,DAY) %>%
  summarise(FEED=sum(FEED_hrs)) %>%
  ungroup() %>%
  group_by(season_id, year, bird_id, season, sex, age_cy) %>%
  summarise(FEED=mean(FEED, na.rm=T))
INDSUMMARY %>% filter(is.na(sex))

PLOTDAT<- INDSUMMARY %>%
  filter(sex %in% c("m","f")) %>% ### remove the unassigned and unknown sexes (3 birds)
  ungroup() %>%
  group_by(season, sex, age_cy) %>%
  summarise(n=length(unique(season_id)), median=quantile(FEED,0.5),lcl=quantile(FEED,0.025),ucl=quantile(FEED,0.975))


### CREATE A GRAPH OF USAGE BY AGE, SEX, and BREEDING SEASON
range(PLOTDAT$ucl)
PLOTDAT %>%  
  mutate(age_cy=ifelse(sex=="m", age_cy-0.2, age_cy+0.2)) %>%
  mutate(season=ifelse(season=="B", "breeding season", "non-breeding season")) %>%
  
  ggplot() +
  geom_point(aes(y=median, x=age_cy, colour=sex),size=2)+
  geom_errorbar(aes(x=age_cy, ymin=lcl, ymax=ucl, colour=sex), width=0.2)+
  scale_y_continuous(name="Mean daily usage (hrs) of anthropogenic feeding", limits=c(0,20), breaks=seq(0,20,4)) +
  scale_x_continuous(name="Age (in years)", limits=c(0.5,8.5), breaks=seq(1,8,1)) +
  facet_wrap(~season, ncol=1, scales="free_y", dir="v") +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=16, color="black"), 
        axis.title=element_text(size=18), 
        strip.text=element_text(size=18, color="black"),
        legend.text=element_text(size=16, color="black"),
        legend.title=element_text(size=18, color="black"),
        legend.key=element_blank(),
        legend.position=c(0.08,0.39),
        legend.background=element_blank(), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())

ggsave("C:/Users/sop/OneDrive - Vogelwarte/General/MANUSCRIPTS/AnthropFeeding/Figure_3.jpg", width=11, height=9, dpi=300)



### REPORT OUTPUT NUMBERS ####

INDSUMMARY %>%
  filter(sex %in% c("m","f")) %>% arrange(desc(FEED))






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
