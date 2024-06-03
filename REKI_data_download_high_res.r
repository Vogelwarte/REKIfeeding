##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
########## DOWNLOAD TRACKING DATA OF RED KITES FROM MOVEBANK  #############
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######################################
### streamlined script using move2 package
### download data for feeding analysis from Movebank

####### LIBRARIES REQUIRED------------------
library(tidyverse)
library(sf)
library(move2)
library(dplyr, warn.conflicts = FALSE)
library(lubridate)
library(data.table); setDTthreads(percent = 65)
library(leaflet)
library(units)
library(foreach)
library(geosphere)
library(keyring)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SET UP DOWNLOAD OF TRACKING DATA FROM MOVEBANK --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

####### SPECIFY THE MOVEBANK ID OF THE STUDY FOR WHICH DATA ARE SHARED
MYSTUDY<-c(230545451,1356790386)
MYUSERNAME<-"Steffen"
movebank_store_credentials(username=MYUSERNAME, key_name = getOption("move2_movebank_key_name"), force = TRUE)
movebank_download_study_info(study_id=MYSTUDY[2])$sensor_type_ids

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DOWNLOAD MOVEBANK DATA AND ANIMAL INFO ----------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

birds1<-movebank_retrieve(study_id=MYSTUDY[1], entity_type="individual") %>%
  dplyr::rename(individual_id=id,bird_id=local_identifier) %>%
  dplyr::select(individual_id,comments, bird_id,ring_id,sex,latest_date_born) 
birds2<-movebank_retrieve(study_id=MYSTUDY[2], entity_type="individual") %>%
  dplyr::rename(individual_id=id,bird_id=local_identifier) %>%
  dplyr::select(individual_id,comments, bird_id,ring_id,sex,latest_date_born)

birds<-bind_rows(birds1,birds2) %>%
  mutate(bird_id=as.numeric(as.character(bird_id)))

### this returns the move2 format which is critical for the filter functions - but they don't work
# locs1<-movebank_download_study(study_id=MYSTUDY[1],
#                                attributes=NULL,
#                          sensor_type_id="gps",
#                          progress=T)
# locs2<-movebank_download_study(study_id=MYSTUDY[2],
#                                attributes=NULL,
#                                sensor_type_id="gps",
#                                progress=T)

locs1<-movebank_retrieve(study_id=MYSTUDY[1],
                         entity_type="event",
                         sensor_type_id="gps",
                         timestamp_start=ymd_hms(paste(year(Sys.time())-10,"-02-15 12:00:00",sep="")),
                         timestamp_end=Sys.time(),
                         progress=T)
locs2<-movebank_retrieve(study_id=MYSTUDY[2],
                         entity_type="event",
                         sensor_type_id="gps",
                         timestamp_start=ymd_hms(paste(year(Sys.time())-10,"-02-15 12:00:00",sep="")),
                         timestamp_end=Sys.time(),
                         progress=T)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FILTER AND COMBINE DATA ----------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### these functions failed on 28 May 2024 with Error: Not all locations are non empty points. Called from: assert_that(mt_is_time_ordered_non_empty_points(x))
## CANNOT BE REINSTATED because these functions require sf objects (download changed to tibble only for speed)
# locs1<-rmOutlier(locs1,maxspeed=35)
# locs1<-speedfilter(locs1,thrspeed=35)
# locs2<-rmOutlier(locs2,maxspeed=35)
# locs2<-speedfilter(locs2,thrspeed=35)


### MANUAL filter function to remove all locations that are within a certain time window
maxtimelag<-60  ## maximum duration of 60 seconds to consider subsequent locations separate
filterlocs<-bind_rows(locs1,locs2) %>%
  group_by(individual_id) %>%
  mutate(prev_t=dplyr::lag(timestamp), prev_id=dplyr::lag(individual_id), prev_lat=dplyr::lag(location_lat),prev_long=dplyr::lag(location_long)) %>%
  mutate(dt=as.numeric(difftime(timestamp,prev_t, units="sec"))) %>%
  mutate(dt=dplyr::if_else(prev_id==individual_id & round(prev_lat,3)==round(location_lat,3) & round(prev_long,3)==round(location_long,3),dt,maxtimelag*2)) %>%
  mutate(dt=dplyr::if_else(is.na(dt),maxtimelag*2,dt)) %>%
  filter(dt>maxtimelag) %>%
  ungroup() %>%
  dplyr::select(timestamp,location_lat,location_long,individual_id)
dim(filterlocs)


#### turned into function to contribute to move2: https://gitlab.com/bartk/move2/-/issues/60


#' @param locs tibble. Tracking data downloaded from Movebank using \code{move2::movebank_retrieve}
#' @param mintimelag numeric. Time lag (in seconds) that must have elapsed between two subsequent locations of the same individual for both locations to be retained.
#' @param loc_res integer. Precision of global coordinates (in EPSG:4326) to which subsequent latitudes and longitudes are rounded when assessing whether identical locations should be filtered.
# mt_filter_lag<-function(locs,maxtimelag=60,loc_res=3){
#   out<-locs %>%
#     group_by(individual_id) %>%
#     mutate(prev_t=dplyr::lag(timestamp), prev_id=dplyr::lag(individual_id), prev_lat=dplyr::lag(location_lat),prev_long=dplyr::lag(location_long)) %>%
#     mutate(dt=as.numeric(difftime(timestamp,prev_t, units="sec"))) %>%
#     mutate(dt=dplyr::if_else(prev_id==individual_id & round(prev_lat,loc_res)==round(location_lat,loc_res) & round(prev_long,loc_res)==round(location_long,loc_res),dt,maxtimelag*2)) %>%
#     mutate(dt=dplyr::if_else(is.na(dt),maxtimelag*2,dt)) %>%
#     filter(dt>maxtimelag) %>%
#     ungroup() %>%
#     dplyr::select(timestamp,location_lat,location_long,individual_id)
#   return(out)
# }
# 
# mt_filter_lag(locs1,maxtimelag=120,loc_res=3)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FILTER AND COMBINE DATA ----------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


locs<-filterlocs %>%
  left_join(birds, by="individual_id") %>%
  mutate(tag_year=as.numeric(comments)) %>%
  mutate(age_cy=as.integer((timestamp-latest_date_born)/365)) %>%
  dplyr::select(bird_id,ring_id,sex,age_cy,timestamp,location_lat,location_long) %>%
  rename(long_wgs=location_long,lat_wgs=location_lat) %>%
  dplyr::filter(!is.na(timestamp)) %>%
  dplyr::filter(!is.na(lat_wgs)) %>%
  dplyr::filter(!is.na(bird_id)) %>%
  dplyr::filter(!is.na(long_wgs)) %>%
  filter(long_wgs<15) %>%
  filter(long_wgs>-10) %>%
  filter(lat_wgs<54) %>%
  filter(lat_wgs>35) %>%
  st_as_sf(coords = c("long_wgs", "lat_wgs"), crs=4326)
head(locs)
dim(locs)


### SAVE THE TRACKING DATA
saveRDS(locs, file = "data/REKI_trackingdata_raw2024.rds")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FILTER DATA BY OUTLIERS AND SPEED ----------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## outlier removal copied from MoveApps: https://github.com/movestore/RemoveOutliers/blob/master/RFunction.R
## speedfilter: https://github.com/movestore/SegmentData-bySpeed/blob/master/RFunction.R

rmOutlier <- function(data, maxspeed=NULL, MBremove=TRUE, FUTUREremove=TRUE, accuracy_var=NULL, minaccuracy=NULL)
{
  Sys.setenv(tz="UTC") 
  
  if (!is.null(accuracy_var))
  {
    if (accuracy_var %in% names(data)==FALSE)
    {
      print("Your defined accuracy variable name does not exist in the data. Please double check. Here it is set to NULL, leading to no removal of high error locations.")
      accuracy_var <- NULL
    }
  }
  
  if (is.null(maxspeed) & MBremove==FALSE & FUTUREremove==FALSE & (is.null(accuracy_var) | is.null(minaccuracy))) print("No maximum speed provided, no accuracy variable/minimum accuracy defined and required to leave Movebank marked Outliers and future timestamp locations in. Return input data set.")
  
  if (!is.null(maxspeed)) print(paste0("Remove positions with maximum speed > ", maxspeed,"m/s")) else print("No maximum speed provided, so no filtering by it.")
  if (MBremove==TRUE) print("In Movebank marked outliers will be removed.") else print("In Movebank marked outliers will be retained.")
  if (FUTUREremove==TRUE) print("Locations with future timestamps will be removed.") else print("Locations with future timestamps will be retained.")
  if (!is.null(accuracy_var) & !is.null(minaccuracy)) print(paste("Remove positions with high location error:",accuracy_var,">",minaccuracy)) else print("Data will not be filtered for location error.")
  
  #take out unrealistic coordinates
  ixNN <- which (st_coordinates(data)[,1]<(-180) | st_coordinates(data)[,1]>180 | st_coordinates(data)[,2]<(-90) | st_coordinates(data)[,2]>90)
  if (length(ixNN)>0)
  {
    print(paste(length(ixNN),"locations have longitude/latitude outside of the usual ranges [-180,180],[-90,90]. Those locations are removed from the data set"))
    data <- data[-ixNN] #if one complete animal is taken out, no problem with moveStack structure :)
  }
  
  data.split <- split(data,mt_track_id(data))
  clean <- lapply(data.split, function(datai) {
    print(unique(mt_track_id(datai)))
    if (MBremove==TRUE) 
    {
      ix <- which(as.logical(datai$visible)==FALSE) #all marked outliers go together in "visible"
      print(paste("For this animal",length(ix),"in Movebank marked outliers were removed."))
      if (length(ix)>0) datai <- datai[-ix,]
    }
    if (!is.null(accuracy_var) & !is.null(minaccuracy)) 
    {
      ixA <- which(as.numeric(as.data.frame(datai)[,accuracy_var])>=as.numeric(minaccuracy))
      print(paste("For this animal",length(ixA),"positions with high errors are removed:",accuracy_var,">",minaccuracy))
      if (length(ixA)>0) datai <- datai[-ixA,] 
    }
    if (!is.null(maxspeed) & nrow(datai)>0) #here changed to while loop (Jan2024)
    {
      len0 <- nrow(datai)
      while (any(units::set_units(mt_speed(datai),m/s)[-nrow(datai)]>units::set_units(maxspeed,m/s)))
      {
        if (length(datai)>1) ixS <- which(units::set_units(mt_speed(datai),m/s)>units::set_units(maxspeed,m/s)) else ixS <- numeric()  #fix for tracks with 1 locations
        if (length(ixS)>0) datai <- datai[-ixS,]
      }
      print(paste("For this animal",len0-nrow(datai),"positions are removed due to between location speeds >",maxspeed,"m/s"))
    }
    datai
  })
  names(clean) <- names(data.split) #clean is still list of move objects
  if (length(clean)>1) clean_move2 <- mt_stack(clean,.track_combine = "rename") else clean_move2 <- clean[[1]]
  
  if (FUTUREremove==TRUE & nrow(clean_move2)>0)
  {
    time_now <- Sys.time()
    clean_nofuture <- lapply(clean, function(cleani) {
      print(unique(mt_track_id(cleani)))
      if (any(mt_time(cleani)>time_now)) 
      {
        ix_future <- which(mt_time(cleani)>time_now)
        print(paste("Warning! Data of the animal",unique(mt_track_id(cleani)),"contain",length(ix_future),"timestamps in the future. They are removed here."))
        cleani[-ix_future] 
      } else 
      {
        print("There are no locations with timestamps in the future.")
        cleani
      }
    })
  } else clean_nofuture <- clean
  
  if (length(clean_nofuture)>1) result <- mt_stack(clean_nofuture,.track_combine="rename") else result <- clean_nofuture[[1]]
  
  if (nrow(result)==0)
  {
    print("Your output file contains no positions. Return NULL.")
    result <- NULL
  }
  
  print(nrow(result))
  
  result
}




########## function to annotate and remove speeds exceeding a threshold

speedfilter <- function(data, speedoption="step", thrspeed=NULL, direc="above")
{
  
  #speedx <- function(x) #input move object
  #{
  #  N <- length(x)
  #  distVincentyEllipsoid(coordinates(x))/as.numeric(difftime(timestamps(x)[-1],timestamps(x)[-N],units="secs"))
  #}
  
  if (speedoption=="ground") 
  {
    if (any(names(data) == "ground.speed" | any(names(data)=="ground_speed"))) print("You have selected to use ground speed for you data selection. This variable existis in your data. However, in case there are locations where ground.speed is NA, distance based speed is estimated (averaged speed from previous location and speed to next location) for the respective steps.") else 
    {
      print("You have selected to use ground speed for your data selection/annotation. However, this variable does not existis in your input data set. Therefore, the calculations are performed using distance based speed estimates (averaged speed from previous location and speed to next location).")
      speedoption <- "step" #
      print(speedoption)
    }
    
  } else print("You have selected to use distance based speed (distance to previous or next location/duration from previous or next location) for your data selection. Note that ground speed at the locations can differ, especially if data resolution is low.")
  
  if (is.null(thrspeed)) 
  {
    print("You have not selected a threshold speed. Please change. Here returning full data set.")
    result <- data
  } else
  {
    if(direc != "annotate") print(paste("You have selected to filter for positions with speed", direc, thrspeed,"m/s")) else print("You have selected to annotate your data with the attribute `speed_class`, indicating if the location is `high` (speed above threshold) or `low` (speed below/equal to threshold).")
    
    thrspeed <- units::set_units(thrspeed,m/s)
    
    data.split <- split(data,mt_track_id(data))
    
    if (direc=="above")
    {
      segm <- foreach(datai = data.split) %do% {
        print(unique(mt_track_id(datai)))
        if (speedoption=="step")
        {
          ix <- which(units::set_units(mt_speed(datai),m/s)>thrspeed)
          dataix <- datai[sort(unique(c(ix,ix+1))),]
        } else
        {
          if (any(names(data) == "ground.speed")) gsi <- units::set_units(datai$ground.speed,m/s) else gsi <- units::set_units(datai$ground_speed,m/s) #if none of those names exist speedoption has been set to "step" above
          
          if (any(is.na(gsi)))
          {
            ixna <- which(is.na(gsi))
            if (1 %in% ixna) 
            {
              gsi[1] <- units::set_units(mt_speed(datai),m/s)[1]
              ixna <- ixna[-1]
            }
            leni <- nrow(datai)
            if (leni %in% ixna)
            {
              gsi[leni] <- units::set_units(mt_speed(datai),m/s)[leni-1]
              ixna <- ixna[-nrow(ixna)]
            }
            if (nrow(ixna)>0)
            {
              gsi[ixna] <- (units::set_units(mt_speed(datai),m/s)[ixna-1]+units::set_units(mt_speed(datai),m/s)[ixna])/2 #average speed of before and after movement
            }
          }
          dataix <- datai[which(gsi>thrspeed),]
        }
        return(dataix)
      }
      names (segm) <- names(data.split)
      
      result <- mt_stack(segm,.track_combine="rename")
      if (dim(result)[1]== 0)
      {
        print("Your output file contains no positions. Return NULL.")
        result <- NULL
      }
      
    } else if (direc=="below")
    {
      segm <- foreach(datai = data.split) %do% {
        print(unique(mt_track_id(datai)))
        if (speedoption=="step")
        {
          ix <- which(units::set_units(mt_speed(datai),m/s)<=thrspeed)
          dataix <- datai[sort(unique(c(ix,ix+1))),]
        } else
        {
          if (any(names(data) == "ground.speed")) gsi <- units::set_units(datai$ground.speed,m/s) else gsi <- units::set_units(datai$ground_speed,m/s)
          if (any(is.na(gsi)))
          {
            ixna <- which(is.na(gsi))
            if (1 %in% ixna) 
            {
              gsi[1] <- units::set_units(mt_speed(datai),m/s)[1]
              ixna <- ixna[-1]
            }
            leni <- nrow(datai)
            if (leni %in% ixna)
            {
              gsi[leni] <- units::set_units(mt_speed(datai),m/s)[leni-1]
              ixna <- ixna[-length(ixna)]
            }
            if (length(ixna)>0)
            {
              gsi[ixna] <- (units::set_units(mt_speed(datai),m/s)[ixna-1]+units::set_units(mt_speed(datai),m/s)[ixna])/2 #average speed of before and after movement
            }
          }
          dataix <- datai[which(gsi<=thrspeed),]
        }
        return(dataix)
      }
      names (segm) <- names(data.split)
      
      result <- mt_stack(segm,.track_combine="rename")
      if (dim(result)[1]== 0)
      {
        print("Your output file contains no positions. Return NULL.")
        result <- NULL
      }
      
    } else if (direc=="annotate")
    {
      if (speedoption=="step")
      {
        ix <- which(units::set_units(mt_speed(data),m/s)<=thrspeed)
        data$speed_class <- NA
        data$speed_class[ix] <- "low"
        data$speed_class[-ix] <- "high"
      } else
      {
        if (any(names(data) == "ground.speed")) gsi <- units::set_units(data$ground.speed,m/s) else gsi <- units::set_units(data$ground_speed,m/s)
        ix <- which (gsi<=thrspeed)
        data$speed_class <- NA
        data$speed_class[ix] <- "low"
        data$speed_class[-ix] <- "high"
      }
      result <- data
      print("Your full data set has been annotated with `speed_class`.")
      
    } else
    {
      print("Your indication of direction was not correct. Returning full data set.")
      result <- data
    }
  }
  
  #Artefakt, plot speed histogram with cut-off
  data.split.nn <- data.split[unlist(lapply(data.split,nrow)>1)] # take out individuals with only one position, else speed error
  
  if (speedoption=="step") 
  {
    hist.tab <- foreach(datai = data.split.nn, .combine=rbind) %do% {
      data.frame("speed"=units::set_units(mt_speed(datai),m/s),"id"=unique(mt_track_id(datai)))
    }
  } else
  {
    hist.tab <- foreach(datai = data.split.nn, .combine=rbind) %do% {
      if (any(names(datai) == "ground.speed")) data.frame("speed"=datai$ground.speed,"id"=unique(mt_track_id(datai))) else data.frame("speed"=datai$ground_speed,"id"=unique(mt_track_id(datai)))
    }
  }
  
  if (!is.null(hist.tab))
  {
    speed.plot <- ggplot(hist.tab, aes(x = speed, fill = id)) +
      geom_histogram(position = "identity", alpha = 0.2, bins = 100) +
      geom_vline(xintercept = thrspeed,lty=2) +
      ggtitle("Histogram of the speeds with selected threshold (unit m/s)")
    
    pdf(appArtifactPath("speed_artefakt.pdf"))
    print(speed.plot)
    dev.off() 
  }
  
  return(result)
}






############ CREATE LEAFLET MAP------------------------


pred.pal <- colorNumeric(c("cornflowerblue","firebrick"), seq(0,1))
leaflet(options = leafletOptions(zoomControl = F)) %>% #changes position of zoom symbol
  setView(lng = median(st_coordinates(locs)[,1]), lat = median(st_coordinates(locs)[,2]), zoom = 10) %>%
  htmlwidgets::onRender("function(el, x) {L.control.zoom({ 
                           position: 'bottomright' }).addTo(this)}"
  ) %>% #Esri.WorldTopoMap #Stamen.Terrain #OpenTopoMap #Esri.WorldImagery
  addProviderTiles("Esri.WorldImagery", group = "Satellite",
                   options = providerTileOptions(opacity = 0.6, attribution = F,minZoom = 5, maxZoom = 20)) %>%
  addProviderTiles("OpenTopoMap", group = "Roadmap", options = providerTileOptions(attribution = F,minZoom = 5, maxZoom = 15)) %>%  
  addLayersControl(baseGroups = c("Satellite", "Roadmap")) %>%  
  
  addCircleMarkers(
    data=locs,
    radius = 3,
    stroke = TRUE, color = "white", weight = 0.5,
    fillColor = "firebrick", fillOpacity = 0.8,
    popup = ~ paste0("bird ID: ", locs$bird_id, "<br>", locs$timestamp)
  ) %>%
  
  addScaleBar(position = "bottomright", options = scaleBarOptions(imperial = F))



