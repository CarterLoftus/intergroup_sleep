

library( ctmm )

library(amt) 

library(ezknitr)
library(ggmap)
library(jpeg)
library(knitr)
library(leaflet)
library(maptools)
library(move)
#library(plotKML)
library(plotly)
library(plotrix)

library(spacetime)
library(stats)
library(tibble)



library(data.table)
library(stringr)
library(lubridate)
library(plyr)
library(sp)
library(adehabitatHR)
library(rgdal)
library(raster)

which_spec <- "baboon"


############## Separate and clean the data ##############

# # enter Movebank login information
# login <- movebankLogin()
# 
# # load the data from Movebank
# testdata <- getMovebankData( study = 'Leopards, vervets, and baboons in Laikipia, Kenya', login = login, removeDuplicatedTimestamps = T )
# 
# # extract the names of the individuals in the data
# IDs <- testdata@idData$local_identifier
# 
# # remove the data. We are going to load the data back in one individual at a time so that we can add their name and their species to their dataframe
# rm( testdata )
# 
# # create an empty vector that we will turn into a list of dataframes with each element corresponding to the data of one individual
# df_list <- c()
# 
# # add each individual's data to the list
# for(i in 1:length( IDs ) ){ # for each individaul...
# 
#   # load their data from Movebank
#   temp <- getMovebankData( study = 'Leopards, vervets, and baboons in Laikipia, Kenya', animalName = IDs[ i ], login = login, removeDuplicatedTimestamps = T )
# 
#   # add the individual's species to its data
#   temp@data$individual_taxon_canonical_name <- temp@idData$taxon_canonical_name
# 
#   # add the individual's name to its data
#   temp@data$individual_local_identifier <- temp@idData$local_identifier
# 
#   # turn its data into a dataframe
#   temp_df <- data.frame( temp@data )
# 
#   # add its dataframe to the ongoing list of individual dataframes
#   df_list[ i ] <- list( temp_df )
# }
# 
# # turn the list of dataframes into one large data frame
# df <- do.call( rbind, df_list )
# 
# # df <- fread( 'Leopards, vervets, and baboons in Laikipia, Kenya.csv' ) ## reads in the data downloaded from Movebank
# 
# # df <- as.data.frame( df ) ## turns the data.table into a data.frame
# 
# df$timestamp <- as.POSIXct(x= df$timestamp, format=c("%Y-%m-%d %H:%M:%S"), tz='UTC') ## turns the timestamp into a POSIX element
# 
# df$local_timestamp <- df$timestamp + 3*60*60  ## this makes the timestamp into local time by adding three hours. But don't be confused by the fact that it is still labeled as UTC. The timestamp is now local Kenyan time. I prefer to keep all timestamps in UTC regardless of their actual time zone
# 
# df <- df[ , c( 'individual_local_identifier', 'local_timestamp' , 'location_long' , 'location_lat' , 'ground_speed' , 'heading' , 'individual_taxon_canonical_name' ) ] ## keep only the necessary columns
# 
# dir.create( paste0( getwd(), '/DATA' ) ) # create the DATA directory
# 
# write.csv(df, "DATA/Leopards, vervets, and baboons in Laikipia, Kenya_trim.csv", row.names = F) ## save this trimmed data.frame for future use
# 
# 
# ## in case the above steps have already been done, one can start here:
# df <- fread( "DATA/Leopards, vervets, and baboons in Laikipia, Kenya_trim.csv" )
# 
# ## turn df from a data.table into a data.frame
# df <- as.data.frame( df )
# 
# ## shorten the column names for convenience
# names( df ) <- c( 'id', 'local_timestamp', 'lon', 'lat', 'speed', 'heading', 'species' )
# 
# df$id <- gsub( ' ', '_', df$id )
# 
# ## remove rows with unsuccessful GPS fixes
# df <- df[ !is.na( df$lat ), ]
# 
# ## Make the local_timestamp column into a POSIX element
# df$local_timestamp <- as.POSIXct( df$local_timestamp, tz = 'UTC')
# 
# # add the UTM coordinates
# crs_longlat <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") # define the lat-lon coordinate reference system
# 
# crs_utm <- CRS("+proj=utm +zone=37 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") # define the UTM coordinate reference system for the study area
# 
# df_sp <- SpatialPointsDataFrame(coords = df[, c('lon','lat')], proj4string = crs_longlat, data = df) # turn the dataframe into a spatial points dataframe
# 
# df_sp <- spTransform(df_sp, crs_utm) # transform the coordinate reference system from lat-lon to UTM
# 
# df <- as.data.frame(df_sp) # turn the spatial points dataframe back into a normal dataframe
# 
# # rename the new columns appropriately
# names(df)[ names(df) == c('lon.1') ] <- 'x'
# names(df)[ names(df) == c('lat.1') ] <- 'y'
# 
# rm( df_sp ) # remove the spatial points dataframe
# 
# ## split the dataframe into a list of dataframes. Each element of the list is all of the data for one species
# split_df <- split( df, df$species )
# 
# ## add a day and time column for each dataframe. The day column is the day of the study period. Day 1 is the first day that any individual of the given species began collecting successful fixes
# 
# for(i in 1:length( split_df ) ){
#   ## add the day column
#   split_df[[ i ]]$day <- as.numeric( as.Date( split_df[[ i ]]$local_timestamp) - min( as.Date( split_df[[ i ]]$local_timestamp ) ) + 1)
# 
#   ## add the time column. This is the same as the local_timestamp, but with the date removed
#   split_df[[ i ]]$time <- str_split_fixed( split_df[[ i ]]$local_timestamp, " ", 2 )[ ,2 ]
# }
# 
# ## separate out the data.frames of the different species
# ## save the data of each species in the alphabetical order of their individuals' names. This is not essential, but is nice for organization.
# 
# verv_df <- split_df$Chlorocebus
# verv_df <- verv_df[order(as.numeric(as.factor(verv_df$id))),]
# 
# leo_df <- split_df$`Panthera pardus`
# leo_df <- leo_df[order(as.numeric(as.factor(leo_df$id))),]
# 
# bab_df <- split_df$Papio
# bab_df <- bab_df[order(as.numeric(as.factor(bab_df$id))),]
# 
# ## write the csv for each species separately
# write.csv(verv_df,"DATA/verv_gps.csv", row.names = F)
# write.csv(leo_df,"DATA/leo_gps.csv", row.names = F)
# write.csv(bab_df,"DATA/bab_gps.csv", row.names = F)

############## Standardizing local_timestamps to determine synchronous locations ############## 


### The analysis proceeds from here using only one species. For this particular analysis, we focus on the vervets

## First, specify which species we are working with

which_spec <- "baboon"  ## options are "baboon", "vervet", or "leopard"

## read in the dataframe

if( which_spec == "baboon" ){
  spec_df <- read.csv( "DATA/bab_gps.csv" )
}else{
  if( which_spec == "vervet" ){
    spec_df <- read.csv( "DATA/verv_gps.csv" )
  }else{
    if( which_spec == "leopard" ){
      spec_df <- read.csv( "DATA/leo_gps.csv" )
    }
  }
}

## number of baboon days of data: 1644
nrow( unique( spec_df[ , c( 'id', 'day' )] ) )

## Make the local_timestamp into a POSIX element. Remember, even though it is labeled UTC, it is actually EAT (East Africa Time)
spec_df$local_timestamp <- as.POSIXct( x= spec_df$local_timestamp, tz = "UTC" )

## Make the time into a character, not a factor
spec_df$time <- as.character( spec_df$time )

## Make a column to label the group that the individual belongs to. Note that this is only relevant for the primates

if( which_spec != "leopard" ){
  spec_df$group <- str_split_fixed( spec_df$id, "_", 3)[ ,2 ]
}

## Visualize the data collection period for each animal
plot( as.numeric( as.factor( spec_df$id ) ) ~ spec_df$local_timestamp, cex=0.3, pch=16, main = "When does each collar collect data?", xlab="", xaxt='n', yaxt='n', ylab="ID") ## Make the plot. Breaks in the solid dots represent breaks in successful data collection for a given individual. 

axis(2,at=1:length(unique(as.factor(spec_df$id))),labels=sort(unique(as.factor(spec_df$id))),las=1,cex=0.3) ## Label the y-axis

axis.POSIXct(1,at=seq(min(spec_df$local_timestamp),max(spec_df$local_timestamp),by="10 day"), labels = as.Date(seq(min(spec_df$local_timestamp),max(spec_df$local_timestamp),by="10 day")), las=2,cex.axis = 0.5) ## Label the x-axis

## make table for supplemental materials

supp_mat_mins <- aggregate( spec_df$local_timestamp, by = list( spec_df$id, spec_df$group ), FUN = min )

names( supp_mat_mins ) <- c( 'id', 'group', 'GPS_start' )

supp_mat_mins$GPS_start <- as.Date( supp_mat_mins$GPS_start, tz = 'UTC' )

supp_mat_mins

supp_mat_maxs <- aggregate( spec_df$local_timestamp, by = list( spec_df$id, spec_df$group ), FUN = max )

names( supp_mat_maxs ) <- c( 'id', 'group', 'GPS_end' )

supp_mat_maxs$GPS_end <- as.Date( supp_mat_maxs$GPS_end, tz = 'UTC' )

supp_mat_maxs


supp_fig_1_temp <- merge( x = supp_mat_mins, y = supp_mat_maxs, by = c( 'id', 'group' ), all = T, sort = F )

dir.create( paste0( getwd(), '/RESULTS' ) )

write.csv( supp_fig_1_temp, 'RESULTS/supp_fig_1_temp.csv', row.names = F )



# first interpolate locations at the quarter hour marks when there is a fix shortly before or shortly after the quarter hour (within 7.5 minutes -- these would have been rounded to the nearest quarter hour anyway using the old method above)

# the loop below is sensitive to the order of the data frame, so make sure the rows are ordered by id and timestamp
spec_df <- spec_df[ order( spec_df$local_timestamp ), ]
spec_df <- spec_df[ order( spec_df$id ), ]

new_spec <- spec_df[ 0 , ] # save an empty dataframe with the same column names as spec_df. We will add each individual's data to this as we loop through each individual

tag_names <- as.character( unique( spec_df$id ) ) # save the name of each individual

# recreate each individual's data with interpolated points at the desired sampling times. Add these recreated dataframes to the new_spec dataframe

for( tag in tag_names ){ # for each individual...
  
  id_dat <- spec_df[ spec_df$id == tag, ] # subset the data to this individual's data
  
  start_time <- floor_date( min( id_dat$local_timestamp ), '15 mins') # find the earliest timestamp of this individual's data, and round down to the nearest 15 minutes
  
  end_time <- ceiling_date( max( id_dat$local_timestamp ), '15 mins') # find the latest timestamp of this individual's data, and round up to the nearest 15 mins
  
  new_id_dat <- data.frame( id = tag, local_timestamp = seq( start_time, end_time, '15 mins'), x = NA, y = NA ) # create a new dataframe for this individual with every 15 minutes accounted for from the start to the end of the period for which this individual has data
  
  counter <- 0 # set the counter to 0
  
  for( i in 1:nrow( new_id_dat ) ){ # for each timestamp of the ideal sampling schedule (each row of the new data frame)...
    
    # save the current timestamp
    curr_timestamp <- new_id_dat$local_timestamp[ i ]
    
    # subset the data to the data occurring prior to this timestamp
    before <- id_dat[ id_dat$local_timestamp <= curr_timestamp, ]
    
    if( nrow( before ) != 0 ){ # if there are data with timestamps occurring prior to this timestamp
      
      before_time <- before$local_timestamp[ nrow( before ) ] # find the maximum time of data occurring prior to this timestamp (i.e. the closest time to this current timestamp in the before data) (because the full data is ordered by timestamp, this subsetting will find us the maximum time)
      
      before_diff <- abs( as.numeric( curr_timestamp - before_time, unit = 'mins' ) ) # find the time difference between the current timestamp of the ideal sampling schedule and the most recent GPS fix
      
      after <- id_dat[ id_dat$local_timestamp >= curr_timestamp, ] # subset the data to the data occurring after this timestamp
      
      if( nrow( after ) != 0 ){ # if there are data with timestamps occurring after this timestamp
        
        after_time <- after$local_timestamp[ 1 ] # find the minimum time of data occurring after this timestamp (i.e. the closest time to this current timestamp in the after data) (because the full data is ordered by timestamp, this subsetting will find us the minimum time)
        
        after_diff <- abs( as.numeric( curr_timestamp - after_time, unit = 'mins' ) ) # find the time difference between the current timestamp of the ideal sampling schedule and the GPS fix following closest in time
        
        total_diff <- abs( as.numeric( before_time - after_time, unit = 'mins' ) ) # find the total time difference between the previous GPS fix and the next GPS fix
        
        if( before_diff < 7.5 | after_diff < 7.5 ){ # if either the previous GPS fix or the following GPS fix occurred close to the current timestamp (within the range that would have rounded to this timestamp if we had just done rounding to the nearest quarter hour)
          
          counter <- counter + 1 # advance the counter. This counter will tell us the proportion of the idealized sampling schedule that we could interpolate a fix for using only very nearby fixes
          
          # interpolate the y UTM and x UTM, by calculating a weighted average of the previous y UTM and x UTM and next y UTM and x UTM, respectively, and weighting by the time difference to the next lon/lat and the previous lon/lat respectively 
          new_id_dat$x[ i ] <- stats::weighted.mean( x = c( before$x[ nrow( before ) ], after$x[ 1 ] ), w = ( c( after_diff, before_diff ) ) )
          
          new_id_dat$y[ i ] <- stats::weighted.mean( x = c( before$y[ nrow( before ) ], after$y[ 1 ] ), w = ( c( after_diff, before_diff ) ) )
          
        }
      }
    }
  }
  
  new_spec <- rbind( new_spec, new_id_dat ) # add this individual's interpolated data to the ongoing dataframe
  
  print( paste0( 'success rate = ', counter/nrow( new_id_dat ) ) ) # print the success rate of this individual's collar
}


## now we will interpolate missing data (i.e. data for timestamps that would not have had any data if we simply rounded to the nearest 15 minutes). We will interpolate all data as long as there GPS fix within an hour before and after the given timestamp, with no more than a total of 75 minutes being between the GPS fix before the timestamp and the GPS fix after the timestamp

new_spec$interp_x <- NA
new_spec$interp_y <- NA

for( tag in tag_names ){ # for each individual...
  
  id_dat <- new_spec[ new_spec$id == tag, ] # subset the data to this individual's data
  
  counter <- 0 # set the counter to 0
  
  for( i in 1:nrow( id_dat ) ){ # for each timestamp...
    
    # save the current timestamp
    curr_timestamp <- id_dat$local_timestamp[ i ]
    
    # subset the data to the data occurring prior to this timestamp (when there is actually data, i.e. id_dat$x and id_dat$y are not simply NA)
    before <- id_dat[ id_dat$local_timestamp < curr_timestamp & !is.na( id_dat$x ), ]
    
    if( nrow( before ) != 0 ){ # if there are data with timestamps occurring prior to this timestamp
      
      before_time <- before$local_timestamp[ nrow( before ) ] # find the maximum time of data occurring prior to this timestamp (because the full data is ordered by timestamp, this subsetting will find us the maximum time)
      
      before_diff <- abs( as.numeric( curr_timestamp - before_time, unit = 'mins' ) ) # find the time difference between the current timestamp and the previous timestamp with data
      
      after <- id_dat[ id_dat$local_timestamp > curr_timestamp & !is.na( id_dat$x ), ] # subset the data to the data occurring after this timestamp
      
      if( nrow( after ) != 0 ){ # if there are data with timestamps occurring after this timestamp
        
        after_time <- after$local_timestamp[ 1 ] # find the minimum time of data occurring after this timestamp (because the full data is ordered by timestamp, this subsetting will find us the minimum time)
        
        after_diff <- abs( as.numeric( curr_timestamp - after_time, unit = 'mins' ) ) # find the time difference between the current timestamp of the ideal sampling schedule and the GPS fix following closest in time
        
        total_diff <- abs( as.numeric( before_time - after_time, unit = 'mins' ) ) # find the total time difference between the previous GPS fix and the next GPS fix
        
        if( before_diff <= 60 & after_diff <= 60 & total_diff <= 75 ){ # if the conditions are met to allow for interpolation of a data point
          
          counter <- counter + 1 # advance the counter
          
          #print( stats::weighted.mean( x = c( before_time, after_time ), w = ( c( after_diff, before_diff ) / sum( after_diff, before_diff ) ) ) )
          
          id_dat$interp_x[ i ] <- stats::weighted.mean( x = c( before$x[ nrow( before ) ], after$x[ 1 ] ), w = ( c( after_diff, before_diff ) ) )
          
          id_dat$interp_y[ i ] <- stats::weighted.mean( x = c( before$y[ nrow( before ) ], after$y[ 1 ] ), w = ( c( after_diff, before_diff ) ) )
          
        }
      }
    }
  }
  
  new_spec[ new_spec$id == tag, ] <- id_dat # insert the individual's updated data back into the full data frame
  
  print( paste0( 'success rate = ', counter/nrow( id_dat ) ) ) # print the success rate of interpolation
}


new_spec$interp_dist <- sqrt( ( new_spec$x - new_spec$interp_x )**2  + ( new_spec$y - new_spec$interp_y )**2 ) # find the distance between the observed GPS fixes and the interpolated GPS fixes

mean( new_spec$interp_dist, na.rm = T) # mean distance between a location and its interpolated counterpart is 38.1 m

median( new_spec$interp_dist, na.rm = T) # median distance between a location and its interpolated counterpart is 15.6 m

sd( new_spec$interp_dist, na.rm = T ) / sum( !is.na( new_spec$interp_dist ) ) # standard error of the distance between a location and its interpolated counterpart is 0.00036 m

interp_df <- new_spec # copy the datframe to a new dataframe where we will edit x and y coordinates if there is missing data or large error

interp_df$interp_missing <- 0 # make an empty column. This column will contain binary values indicating whether the x UTM and y UTM used are interpolated due to missing observed location at that timestamp

interp_df$interp_error <- 0 # make an empty column. This column will contain binary values indicating whether the x UTM and y UTM used are interpolated due to large discrepency betweeen the observed location and the interpolated location

interp_df[ ( is.na( interp_df$x ) & !is.na( interp_df$interp_x ) ), 'interp_missing' ] <- 1 # mark the interpolated observations in the interp column

interp_df[ ( is.na( interp_df$x ) & !is.na( interp_df$interp_x ) ), c( 'x', 'y' ) ] <- interp_df[ ( is.na( interp_df$x ) & !is.na( interp_df$interp_x ) ), c( 'interp_x', 'interp_y' ) ] # insert the interpolated GPS fixes into the GPS fix columns (x and y UTM locations)

interp_df <- interp_df[ !is.na( interp_df$x ), ] # remove rows where we don't have a fix and weren't able to interpolate one

# use the interpolated location for GPS fixes that are greater than the 99th percentile away from their interpolated counterpart

interp_df[ which( interp_df$interp_dist > quantile( interp_df$interp_dist, 0.99, na.rm = T ) ), 'interp_error' ] <- 1 # mark the locations as interpolated

interp_df[ which( interp_df$interp_dist > quantile( interp_df$interp_dist, 0.99, na.rm = T ) ), c( 'x', 'y' ) ] <-  interp_df[ which( interp_df$interp_dist > quantile( interp_df$interp_dist, 0.99, na.rm = T ) ), c( 'interp_x', 'interp_y' ) ] # insert the interpolated GPS fixes in the GPS columns (x and y UTM locations, overwriting the information that was previously here)

interp_df$interp <- interp_df$interp_missing + interp_df$interp_error # make a new column that says generally whether the location was interpolated

sum( interp_df$interp_missing, na.rm = T ) # number of missing fixes that were replaced with an interpolation
sum( interp_df$interp_error, na.rm = T ) # number of anomalous fixes that were replace with an interpolation

interp_df <- interp_df[order(interp_df$local_timestamp),]
interp_df <- interp_df[order(as.numeric(as.factor(interp_df$id))),]

## add a day and time column for each dataframe. The day column is the day of the study period. Day 1 is the first day that any individual of the given species began collecting successful fixes

## add the day column
interp_df$day <- as.numeric(as.Date(interp_df$local_timestamp) - min(as.Date(interp_df$local_timestamp)) + 1)

## add the time column. This is the same as the local_timestamp, but with the date removed
interp_df$time <- str_split_fixed(interp_df$local_timestamp," ",2)[,2]

if( which_spec != 'leopard' ){
  
  ## Make a column to label the group that the individual belongs to
  interp_df$group <- str_split_fixed(interp_df$id,"_",3)[,2]
  
}

# add the lon lat data back to the dataframe

crs_longlat <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") # define the lat-lon coordinate reference system

crs_utm <- CRS("+proj=utm +zone=37 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") # define the UTM coordinate reference system for the study area

interp_df_sp <- SpatialPointsDataFrame( coords = interp_df[, c('x','y')], proj4string = crs_utm, data = interp_df ) # turn the dataframe into a spatial points dataframe

interp_df_sp <- spTransform( interp_df_sp, crs_longlat ) # transform the coordinate reference system from lat-lon to UTM

interp_df_lonlat <- as.data.frame( interp_df_sp ) # turn the spatial points dataframe back into a normal dataframe

# rename the new columns appropriately
names(interp_df_lonlat)[ names(interp_df_lonlat) == c('x.1') ] <- 'lon'
names(interp_df_lonlat)[ names(interp_df_lonlat) == c('y.1') ] <- 'lat'

rm( interp_df_sp ) # remove the spatial points dataframe

## write the new dataframe into a csv

if(which_spec == "baboon"){
  write.csv(interp_df_lonlat, "DATA/bab_df.csv", row.names = F)
}else{
  if(which_spec == "vervet"){
    write.csv(interp_df_lonlat, "DATA/verv_df.csv", row.names = F)
  }else{
    if(which_spec == "leopard"){
      write.csv(interp_df_lonlat, "DATA/leo_df.csv", row.names = F)
    }
  }
}





