

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


which_spec <- 'baboon'

if(which_spec == "baboon"){
  spec_df <- read.csv( "DATA/bab_dyad.csv" )
}else{
  if(which_spec == "vervet"){
    spec_df <- read.csv( "DATA/verv_dyad.csv" )
  }else{
    if(which_spec == "leopard"){
      spec_df <- read.csv( "DATA/leo_dyad.csv" )
    }
  }
}



## find the distance between group-mates for the Methods section+

groupmate_df <- spec_df[ ( spec_df$group1 == spec_df$group2 ) & ( spec_df$id1 != spec_df$id2 ), ]

groupmate_df <- groupmate_df[ !duplicated( groupmate_df[ c( 'group1', 'group2', 'local_timestamp' ) ] ), ]

mean( groupmate_df$dyadDist, na.rm = T )

sd( groupmate_df$dyadDist, na.rm = T ) / sqrt( sum( !is.na( groupmate_df$dyadDist ) ) )

# remove the rows that give the dyadic distances between groupmates (or between an individual and itself, which are currently entered as NAs). We don't need these rows for identifying intergroup interactions
spec_df <- spec_df[ spec_df$group1 != spec_df$group2, ]

# remove the rows for which we don't have dyadic distance data
spec_df_nonas <- spec_df[ !is.na( spec_df$dyadDist ), ]

# only keep the necessary columns
spec_df_trim_nonas <- spec_df_nonas[ , c( 'group1', 'group2', 'local_timestamp', 'dyadDist' ) ]

head( spec_df_nonas )

# make the timestamp a POSIX character
spec_df_trim_nonas$local_timestamp <- as.POSIXct( spec_df_trim_nonas$local_timestamp, tz = 'UTC' )

## the following function will provide a unique name to each unique group dyad
dy_name <- function(vec) {
  
  temp <- sort( c( vec[ 'group1' ], vec[ 'group2' ] ) )
  
  return( paste( temp[1], temp[2], sep = '-') )
}

if(which_spec != "leopard"){
  
  spec_df_trim_nonas$group1 <- as.character( spec_df_trim_nonas$group1 )
  spec_df_trim_nonas$group2 <- as.character( spec_df_trim_nonas$group2 )
  
  ## use function from above to assign dyad ID name
  spec_df_trim_nonas$dyadID <- apply( spec_df_trim_nonas, 1, FUN = dy_name)
  
}


## for the primates, we want to know encounters at the group level, not the individual level. So we will remove rows from the dyadic distance dataframe that are redundant. (i.e. if group 1 has individuals A and B and group 2 has individual C, we don't need to know about an encounter between A and C and an encounter between B and C if they are simultaneous). We will keep the rows that minimize the dyadic distance so we don't miss any encounters that may have potentially happened.

agg_dat <- aggregate( spec_df_trim_nonas$dyadDist, by = list( spec_df_trim_nonas$local_timestamp, spec_df_trim_nonas$dyadID ), FUN = min )

names( agg_dat ) <- c( 'local_timestamp', 'dyadID', 'dyadDist')


# # for a particular dyad, we want group1 to always be group1 and group2 to always be group2, for use in the dyadDurat code below. So let's use the sorted dyad names to make this happen
agg_dat$group1 <- str_split_fixed( agg_dat$dyadID, '-', 2 )[ , 1 ]

agg_dat$group2 <- str_split_fixed( agg_dat$dyadID, '-', 2 )[ , 2 ]

# reorder the dataframe by id1, id2, and timestamp
agg_dat <- agg_dat[ order( agg_dat$group2 ), ]
agg_dat <- agg_dat[ order( agg_dat$local_timestamp ), ]
agg_dat <- agg_dat[ order( agg_dat$group1 ), ]



############## Use the dyadic distances to pull out encounters ##############


## set the distance threshold for what we consider a potential encounter
dist_thres <- 600 # change this if runnning the leopard interactions

## set a time threshold. Encounters that happen within this number of minutes will be merged with each other rather than considered as distinct encounters
time_thres <- 75 ## interactions that occur within this number of minutes of each other will be merged

## set the sampling interval of the study
samp_int <- 15

sum(duplicated(agg_dat[,c('local_timestamp','dyadID')])) ## this should be 0 to confirm that we have removed redundant information


## create a dataframe with 0 rows that we will add rows to
dyadDurat<-data.frame(id1=character(), id2=character(), start_local_timestamp=character(), end_local_timestamp=character(), duration=character())

## save the names of the individuals in the study
group_names <- unique( c( agg_dat$group1, agg_dat$group2 ) )

for(a in 1:(length(group_names)-1)){
  
  for(b in (a+1):(length(group_names))){
    
    # for each dyad, subset the dyadic distance data to only the data for the dyad members, and to observations in which they are closer to each other than the distance threshold hyperparameter
    tempDF <- agg_dat[(agg_dat$group1 == group_names[ a ] & agg_dat$group2 == group_names[ b ] & agg_dat$dyadDist < dist_thres) | (agg_dat$group1 == group_names[ b ] & agg_dat$group2==group_names[ a ] & agg_dat$dyadDist < dist_thres) , ]
    
    if(nrow(tempDF) != 0){ # if the dyad comes within this distance of each other...
      # find the time differences between observations of potential encounters
      diff_min <- diff(tempDF$local_timestamp)
      
      # find out how many of these are consecutive GPS fixes
      tsig <- c(F,(abs(diff_min) == samp_int))
      
      # find the indices that represent the start of a potential encounter (i.e. more than the sampling interval after the previous time the dyad was in close proximity)
      startIndex <- which(tsig == F)
      
      # find the indices associated with the ends of potential encounters (when the next time of close proximity is more than the sampling interval after an observation of proximity)
      endIndex <- c((startIndex[-1] - 1), nrow(tempDF))
      
      # find the duration of each potential encounter
      dur_interact <- as.numeric(tempDF$local_timestamp[endIndex] - tempDF$local_timestamp[startIndex], units = 'mins')
      
      if(length(startIndex) >= 1){ # if there is at least one potential encounter (which there should be because we only went into this if statement for dyads that do come within the threshold for proximity)
        
        # make a temporary data frame documenting this potential encounter
        tempDF2 <- data.frame( group1 = rep(group_names[a], each = length( startIndex )) , group2 = rep(group_names[b], each = length( startIndex )), start_local_timestamp = tempDF$local_timestamp[startIndex], end_local_timestamp = tempDF$local_timestamp[endIndex], duration= dur_interact )
        
        # add this temporary dataframe to the running dataframe of all potential encounter events in the data
        dyadDurat<-rbind(dyadDurat,tempDF2)
        
      }
    }
  }
}


## Add the group identity and the dyad name to the encounter data
if(which_spec != 'leopard'){
  
  # add dyadID
  dyadDurat$dyadID <- apply(dyadDurat, 1, FUN = dy_name)
}else{
  
  #for leopards, the dyad ID will just be the individual-level dyad, rather than the group-level dyad used for the monkeys
  dyadDurat$dyadID <- paste(dyadDurat$group1, dyadDurat$group2, sep = '_')
}

par( bg = 'black' )

# Plot a histogram of encounter durations of all encounters
hist( dyadDurat$duration, xlab='Duration (mins)', main = '', col = 'white', col.axis = 'white', col.lab = 'white' )
axis( 1, labels = F, col = 'white' )
axis( 2, labels = F, col = 'white' )

# Plot a historgram of the dyadic distances of all dyads (except intra-group dyads)
hist(agg_dat$dyadDist,xlab='Distance (m)',main='Histogram of dyadic distances over the study duration (excluding groupmates)')

head(dyadDurat)

## merge interactions that occur within the time threshold of each other. For the monkeys, this merge will occur across members of the same group. So it will merge all potential encounters between two groups within time threshold, not just two individuals.

## create a dataframe with 0 rows. After we merged encounter data to it, this will be the final data of potential encounters
finalDurat <- data.frame(id1=character(), id2=character(), start_local_timestamp=character(), end_local_timestamp=character(), duration=character())

## save the unique group dyads
dyad_names <- unique(dyadDurat$dyadID)

for(i in 1:length(dyad_names)){ # for each unique dyad...
  
  # subset the encounter data to just this dyad
  dyadDF <- dyadDurat[dyadDurat$dyadID == dyad_names[i],]
  
  if(nrow(dyadDF) != 0){ # if this dyad does have encounters...
    
    # order the encounters in consecutive order by the time the encounter started
    dyadDF <- dyadDF[order(dyadDF$start_local_timestamp),]
    
    # find the time between the end of one interaction and the start of the next
    time_between <- as.numeric(dyadDF$start_local_timestamp[-1] - dyadDF$end_local_timestamp[- nrow(dyadDF)], units = 'mins')
    
    # determine which of these are under the time threshold required for them to be considered two distinct encounters, and render a vector of booleans telling whether they are NOT distinct (i.e. True means they should be merged)
    bool <- c(F, time_between < time_thres)
    
    # find when the distinct encounters (of those that require merging) end (indicated by a -1)
    inds <- diff( c( bool, F ) )
    
    # in the row where an encounter starts, replace the original end time with the end time that should result from merging encounters within the threshold time
    dyadDF$end_local_timestamp[inds == 1] <- dyadDF$end_local_timestamp[inds == -1]
    
    # only keep the rows that have the start of a distinct encounter
    dyadDF <- dyadDF[ bool == F, ]
    
    # add this dyad's ID to the total running dataframe containing all of the merged encounters
    finalDurat <- rbind(finalDurat, dyadDF)
  }
  
}


# recalcute the durations of the encounters
finalDurat$duration <- as.numeric(finalDurat$end_local_timestamp - finalDurat$start_local_timestamp, units = 'mins')

finalDurat <- finalDurat[ order( finalDurat$group2 ), ]
finalDurat <- finalDurat[ order( finalDurat$group1 ), ]

# not sure if we actually need this
spec_df$local_timestamp <- as.POSIXct( spec_df$local_timestamp, tz = 'UTC' )


# write the csv of the data of the potential encounters
dir.create( paste0( getwd(), "/DATA/DyadDurats/" ) )

rownames( finalDurat ) <- NULL

write.csv( finalDurat, paste("DATA/DyadDurats/", which_spec, "_dyadDurat.csv", sep = "" ), row.names = F)


### creating KMLs of encounters #####
finalDurat <- read.csv( paste( "DATA/DyadDurats/", which_spec, "_dyadDurat.csv", sep = "" ) )

finalDurat$start_local_timestamp <- as.POSIXct( finalDurat$start_local_timestamp, tz = 'UTC' )
finalDurat$end_local_timestamp <- as.POSIXct( finalDurat$end_local_timestamp, tz = 'UTC' )

# read in the GPS data that we will use for the KML
if(which_spec == "baboon"){
  gps_df <- read.csv( "DATA/bab_complete.csv" )
}else{
  if(which_spec == "vervet"){
    gps_df <- read.csv( "DATA/verv_complete.csv" )
  }else{
    if(which_spec == "leopard"){
      gps_df <- read.csv( "DATA/leo_complete.csv" )
    }
  }
}


head( gps_df )

gps_df <- gps_df[ !is.na( gps_df$x ), ]

gps_df$local_timestamp <- as.POSIXct( gps_df$local_timestamp, tz = 'UTC' )

#### Split id1 and id2 so they can have different colours

### Make sure to set correct directory, all kml-files will be stored here!

dir.create( paste0( getwd( ), '/DATA/interaction_KMLs' ) )

setwd( paste0( getwd( ), '/DATA/interaction_KMLs' ) )

idData <- data.frame( matrix( nrow = nrow( finalDurat ), ncol = 4 ) )

names(idData)<-c( 'interaction_id', 'original_row', 'group1', 'group2' )

library(plotKML)
library(sp)

crs_utm <- CRS("+proj=utm +zone=37 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") # define the UTM coordinate reference system for the study area

# Creating kml's

for(i in 1:nrow( finalDurat ) ){ ### for each interaction
  
  start.time <- finalDurat$start_local_timestamp[i] - 60*60## save the time an hour before the interaction begins
  
  end.time <- finalDurat$end_local_timestamp[i] + 60*60 ## save the time an hour after the interaction ends
  
  # susbset the data to just the time of the interaction
  gps_sub <- gps_df[ gps_df$local_timestamp >= start.time & gps_df$local_timestamp <= end.time, ]
  
  # save the names the individuals who have data at this time
  present_tags <- as.character( unique( gps_sub$id ) )
  
  # plot the KML for each individual
  for( t in 1:length( present_tags ) ){ #for each individual...
    
    id_gps_df <- gps_sub[ gps_sub$id == present_tags[ t ], ] # subset to just this individual's data
    
    id_gps_df <- id_gps_df[ order( id_gps_df$local_timestamp ), ] # order the data by timestamp
    
    group <- as.character( unique( id_gps_df$group ) ) # save the group identity that this individual is a part of
    
    times <- id_gps_df$local_timestamp # save the timestamps in their own vector
    
    coordinates( id_gps_df )<-c("x","y") # make x and y the coordinates
    
    proj4string( id_gps_df )<- crs_utm # set the projection
    
    id_gps_df_sp <-as( id_gps_df ,"SpatialPoints") # turn the individual's data into a spatial points dataframe
    
    id_ST<- STIDF( id_gps_df_sp, times, data=data.frame( time = id_gps_df$local_timestamp ) ) # turn the data into a spacetime object
    
    # plot the KML, with the color depending on the identity of the group - group1 will be colored yellow to magenta, group2 will be colored 
    if( group == finalDurat$group1[ i ] ){
      
      plotKML( id_ST, points_names = paste0( 'group1_', t ), colour_scale = c( "lightyellow", "darkgreen" ),
               file.name = paste( i, '_', gsub( ':', '', gsub( '-', '', gsub( ' ', '_', finalDurat$start_local_timestamp[i] ) ) ), '_', t, ".kml",sep = "" ), fixed = TRUE)
      
    }else{
      
      if( group == finalDurat$group2[ i ] ){
        
        plotKML( id_ST, points_names = paste0( 'group2_', t ), colour_scale = c( "lightblue", "darkblue" ),
                 file.name = paste( i, '_', gsub( ':', '', gsub( '-', '', gsub( ' ', '_', finalDurat$start_local_timestamp[i] ) ) ), '_', t, ".kml",sep = "" ), fixed = TRUE)
        
      }else{
        
        plotKML( id_ST, points_names = t, colour_scale = c( "lightgrey", "darkgrey" ),
                 file.name = paste( i, '_', gsub( ':', '', gsub( '-', '', gsub( ' ', '_', finalDurat$start_local_timestamp[i] ) ) ), '_', t, ".kml",sep = "" ), fixed = TRUE)
        
        
      }
    }
  }
}

# number in front of the interaction KML files corresponds to rows of finalDurat

