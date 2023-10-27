###########################################################################################################################
###### DATA PREPARATION TO IDENTIFY ANTHROPOGENIC FEEDING SITES OF RED KITES BASED ON GPS TRACKING DATA ################
###########################################################################################################################
# original script by Nathalie Heiniger (MSc thesis 2020)
# uses csv table instead of Movebank download to save time
# modified by steffen.oppel@vogelwarte.ch in June 2023
# ALTERNATIVE APPROACH ATTEMPT to use aniMotum for gamma estimates - NEED TO WEED OUT DATA WITH FEW LOCATIONS FIRST

### DATA INPUT IS A PERSISTENT HEADACHE
### originally tried Movebank download, but failed, then thinned the data but recursions fail with std::bad_alloc (probably memory overload)
### reverted to using hourly data that were pre-prepared by Ursin Beeli, as no way forward using other data


### FINALISED on 30 JUNE 2023 with hourly data

library(aniMotum)
library(tidyverse)
library(rnaturalearth)
library(sf)
library(amt)
library(dplyr, warn.conflicts = FALSE)
options(dplyr.summarise.inform = FALSE)
library(recurse)
library(lubridate)
#library(move)
library(data.table); setDTthreads(percent = 65)
sf_use_s2(FALSE) # deactivating spherical geometry s2


## set root folder for project
setwd("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/REKIFeeding")




###  LOADING DATA -----------------------------------------------------------------
# ATTEMPT TO LOAD DATA FROM MOVEBANK ABANDONED ON 20 June 2023 because of memory limits
# saved Movebank data to csv file and read in csv file on 28 June 2023
# instead using data given to Ursin Beeli for NestTool analysis
movemil_GSM <- read.csv("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/NestTool/REKI/output/02_preprocessing/01_milvus_gsm.csv")
movemil_Milsar <- read.csv("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/NestTool/REKI/output/02_preprocessing/02_milvus_milsar.csv")
dat_GSM <- as.data.frame(movemil_GSM) %>%
  dplyr::select(individual.local.identifier,timestamp,year_id,year,month,day,long_wgs,lat_wgs,long_eea,lat_eea) %>%
  mutate(timestamp=as.POSIXct(timestamp,format ="%Y-%m-%d %H:%M:%S", tz = "UTC")) %>%
  dplyr::filter(!is.na(timestamp))
myDF <- as.data.frame(movemil_Milsar) %>%
  dplyr::select(individual.local.identifier,timestamp,year_id,year,month,day,long_wgs,lat_wgs,long_eea,lat_eea) %>%
  mutate(timestamp=as.POSIXct(timestamp,format ="%Y-%m-%d %H:%M:%S", tz = "UTC")) %>%
  dplyr::filter(!is.na(timestamp)) %>%
  bind_rows(dat_GSM) %>%
  arrange(year_id,timestamp)

############# THIS CHUNK LAST RUN ON 28 JUNE 2023 ##############
# movemil_Milsar <- fread("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/Feeding/data/Milvusmilvus_Milsar_SOI_final.csv")    ### insufficient memory to read in those data
# dim(movemil_Milsar)
# 
# 
# #read in your data directly from movebank
# username <- "Steffen"
# password <- "xxxxx"
# curl <- movebankLogin(username=username,  password=password)
# #
# # #read in your data directly from movebank
# # getMovebankID(study="Milvusmilvus_Milsar_SOI_final", login=curl)
# # unique(getMovebankSensors(study=1356790386, login=curl)$sensor_type_id)
# # getMovebankSensors(login=curl)
# movemil_GSM <- getMovebankLocationData(study="Milvusmilvus_GSM_SOI", sensorID=653, login=curl)
# # movemil_Milsar <- getMovebankLocationData(study="Milvusmilvus_Milsar_SOI_final", sensorID=653, login=curl)  ## does not work due to Individual(s): '423' do(es) not have data for one or more of the selected sensor(s).
# 
# 
# dat_GSM <- as.data.frame(movemil_GSM) %>%
#   dplyr::select(individual.local.identifier,timestamp,location.long,location.lat) %>%
#   rename(id=individual.local.identifier, lat=location.lat,long=location.long) %>%
#   mutate(year=year(timestamp),
#           month=month(timestamp),
#           day=day(timestamp),
#           yday=yday(timestamp),
#          year_id=paste(year,id,sep="_")) %>%
#   dplyr::filter(!is.na(timestamp))
# 
# 
# myDF <- as.data.frame(movemil_Milsar) %>%
#   dplyr::select(`individual-local-identifier`,timestamp,`location-long`,`location-lat`) %>%
#   rename(id=`individual-local-identifier`, lat=`location-lat`,long=`location-long`) %>%
#   mutate(year=year(timestamp),
#          month=month(timestamp),
#          day=day(timestamp),
#          yday=yday(timestamp),
#          year_id=paste(year,id,sep="_")) %>%
#   dplyr::filter(!is.na(timestamp)) %>%
#   bind_rows(dat_GSM) %>%
#   arrange(year_id,timestamp)
# dim(myDF)
trackingdata <- myDF[!duplicated(paste0(myDF$timestamp,myDF$id)),] ## this is to exclude duplicated timestamps (if present)
saveRDS(trackingdata, file = "data/REKI_trackingdata_raw.rds", version=3)
dim(trackingdata)
rm(myDF,dat_GSM,movemil_GSM,movemil_Milsar)

### LOAD THE TRACKING DATA AND INDIVIDUAL SEASON SUMMARIES
trackingdata<-readRDS(file = "data/REKI_trackingdata_raw.rds")



# DATA PREPARATION -------------------------------------------------------------
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
  filter(!is.na(year_id)) #%>%
  #filter(!(year_id %in% c("2017_164","2017_228")))  # shit data only 10 locations between 14 Feb and 10 Aug



#### FILTER OUT SHIT SEASONS
eliminate<-trackingdata %>% group_by(year_id) %>%
	summarise(N=length(unique(timestamp))) %>%
	arrange(N) %>%
	filter(N<100)

dim(trackingdata)

trackingdata<-trackingdata %>% filter(!(year_id %in% eliminate$year_id))
dim(trackingdata)


# ATTEMPT TO USE SSM TO GET GAMMA ESTIMATES--------------------------------------------------
head(sese)
# GPS data should have 5 columns in the following order: id, date, lc, lon, lat. 
SSMdat<-trackingdata %>% rename(date=timestamp, lon=long,id=year_id) %>%
	mutate(lc="G") %>%
	select(id, date, lc, lon, lat)
dim(SSMdat)
fit <- fit_ssm(SSMdat, 
               vmax= 50, 
               model = "mp", 
               time.step = 1,
               control = ssm_control(verbose = 0))

#plot(fit, type = 3, pages = 1, ncol = 2)




##### EXTRACT THE DATA FROM THE LIST OF LISTS ####

ssm_out<-fit$ssm[[1]]$predicted
for (l in 2:length(unique(SSMdat$id))) {
 ssm_out<-bind_rows(ssm_out,fit$ssm[[l]]$predicted)
}
dim(ssm_out)
saveRDS(ssm_out, file = "data/REKI_SSM_AniMotum_out.rds")


ssm_outfit<-fit$ssm[[1]]$fitted
for (l in 2:length(unique(SSMdat$id))) {
 ssm_outfit<-bind_rows(ssm_out,fit$ssm[[l]]$fitted)
}
dim(ssm_outfit)
saveRDS(ssm_outfit, file = "data/REKI_SSM_AniMotum_outfitted.rds")





##### CHECK FOR SOME IND WHY THEY HAVE MORE LOCS AFTER SSM THAN BEFORE ####
SSMdat %>% filter(id=="2018_271")
ssm_out %>% filter(id=="2018_271")