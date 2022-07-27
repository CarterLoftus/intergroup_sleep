

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

library( hms )

library(data.table)
library(stringr)
library(lubridate)
library(plyr)
library(sp)
library(adehabitatHR)
library(rgdal)
library(raster)
library(rgeos)


which_spec <- 'baboon'

if(which_spec == "baboon"){
  spec_df <- read.csv("DATA/bab_df.csv")
}else{
  if(which_spec == "vervet"){
    spec_df <- read.csv("DATA/verv_df.csv")
  }else{
    if(which_spec == "leopard"){
      spec_df <- read.csv("DATA/leo_df.csv")
    }
  }
}

# save the coordinate reference systems
crs_longlat <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") # define the lat-lon coordinate reference system

crs_utm <- CRS("+proj=utm +zone=37 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") # define the UTM coordinate reference system for the study area


# Make local_timestamp POSIX class 
spec_df$local_timestamp<-as.POSIXct( spec_df$local_timestamp, format="%Y-%m-%d %H:%M:%OS", tz="UTC" )

# make a column for HH:MM:SS time
spec_df$time <- as.character( spec_df$time )

# remove the rows that don't have GPS fixes
spec_df <- spec_df[ !is.na( spec_df$x ), ]

# save the names of the individuals
tag_names <- as.character( unique( spec_df$id ) )


# how many points are interpolated

sum( spec_df$interp_missing ) 

sum( spec_df$interp_missing ) / nrow( spec_df )

sum( spec_df$interp_error ) 

sum( spec_df$interp_error ) / nrow( spec_df )




### Spatially discretize to account for GPS error ###

# make a column for the step length from the previous timestamp to the current timestamp
spec_df$spat.disc.dist <- NA

# make a column for the step length from the current timestamp to the next timestamp (this is a bit redundant with that from above, but I just want to make sure this variable is lined up correctly because I will use it in statistical analysis)
spec_df$spat.disc.dist_next <- NA

# make a column for the spatially discretized x coordinate
spec_df$x_disc <- spec_df$x

# make a column for the spatially discretized y coordinate
spec_df$y_disc <- spec_df$y

n <- 30 ### distance threshold to spatially discretize at

# discretize each individual's data and then add it back to the full dataframe

for(j in unique(tag_names)){ # for each individual...
  
  idDat <- spec_df[spec_df$id == j,] ## subset the data for the given individual
  
  temp.y<-idDat$y[1] ## set their first location to the first temp location
  temp.x<-idDat$x[1]
  
  temp.count <- 1
  
  dist <- 0
  
  i <- 0 ## set counter to 0
  
  totalDist <- 0 ## Set step length distance to 0
  
  while(i <= nrow(idDat)){
    
    while(dist < n){
      
      if(i == nrow(idDat)){
        break
        
      }
      
      if(i == 0 || i == temp.count){
        i <- i+1
        
        dist <- sqrt((idDat$y[i]-temp.y)**2 + (idDat$x[i]-temp.x)**2)
        
      }else{
        
        idDat$spat.disc.dist[i]<-0
        idDat$y_disc[i] <- idDat$y[temp.count]
        idDat$x_disc[i] <- idDat$x[temp.count]
        i<-i+1
        dist<-sqrt((idDat$y[i]-temp.y)**2 + (idDat$x[i]-temp.x)**2)
      }
    }
    if(dist < n){
      idDat$spat.disc.dist[i] <- 0
      break
      
    }else{
      
      idDat$spat.disc.dist[i] <- dist
      temp.y<-idDat$y[i]
      temp.x<-idDat$x[i]
      temp.count <- i
      dist <- 0
      
    }
  }
  
  idDat$spat.disc.dist_next <- c( idDat$spat.disc.dist[ - 1 ], NA )
  
  spec_df[spec_df$id == j,] <- idDat
  
}


sum( spec_df$spat.disc.dist_next !=  c( spec_df$spat.disc.dist[ - 1 ], NA ), na.rm = T )

## finding heading
## current heading is calculated as the heading that the baboon took to get from the previous location to the current location. Be careful not to confuse this when you get to rel_heading_dyad later, which is the heading that the baboon takes in the step following that fix timestamp

spec_df$heading <- NA

for(g in tag_names){
  idDat <- spec_df[spec_df$id == g, ]
  
  vec <- c( NA )
  
  for(i in 2:(nrow(idDat))){
    
    y_disc <- idDat$y_disc[i]
    x_disc <- idDat$x_disc[i]
    prev.y_disc <- idDat$y_disc[i-1]
    prev.x_disc <- idDat$x_disc[i-1]
    north.y_disc <- prev.y_disc + 100
    north.x_disc <- prev.x_disc
    if(x_disc == prev.x_disc){
      if(y_disc == prev.y_disc){
        angle <- NA
      }else{
        if(y_disc > prev.y_disc){
          angle <- 0
        }else{
          angle <- pi
        }
      }
    }else{
      a <- sqrt((y_disc-prev.y_disc)**2+(x_disc-prev.x_disc)**2)
      b <- 100
      c <- sqrt((y_disc-north.y_disc)**2+(x_disc-north.x_disc)**2)
      angle <- acos((a**2 + b**2 - c**2) / (2 * a * b))
    }
    if(x_disc < prev.x_disc){
      angle <- 2*pi - angle
    }
    vec <- c(vec, angle)
  }
  spec_df[spec_df$id == g, 'heading'] <- vec
}


## finding turning angle

## the turning angle at a given location is the turn that is about to be made (or perhaps more accurately, the turn that is in progress). The current location is the vertex of the two vectors whose angle with each other represent the turning angle
# subtract new heading from previous heading. If the difference is more than 180, then subtract 360. If the difference is less than -180, add 360.

spec_df$turn_angle <- NA

s_days <- as.character(unique(spec_df$day))

for( id in tag_names ){
  
  id_df <-  spec_df[spec_df$id == id, ]
  
  for( day in s_days ){
    
    day_df <-  id_df[ id_df$day == day, ]
    
    if(nrow(day_df) != 0){
      
      for(i in 1:( length(day_df$heading)-1) ){
        
        if( !is.na(day_df$heading[ i + 1 ]) ){
          
          till_now <- day_df$heading[ 1 : i ]
          
          inds <- which( !is.na( till_now ) )
          
          last_ind <- inds[ length(inds) ]
          
          if( length(last_ind) != 0){
            
            day_df$turn_angle[i] <- day_df$heading[ i + 1 ] - day_df$heading[ last_ind ]
            
          }
        }
      }
    }
    id_df[ id_df$day == day, ] <- day_df
  }
  spec_df[spec_df$id == id, ] <- id_df
}

spec_df$turn_angle[ spec_df$turn_angle > pi & !is.na(spec_df$turn_angle)] <- spec_df$turn_angle[ spec_df$turn_angle > pi & !is.na(spec_df$turn_angle)] - 2*pi

spec_df$turn_angle[ spec_df$turn_angle < - pi & !is.na(spec_df$turn_angle)] <- spec_df$turn_angle[ spec_df$turn_angle < - pi & !is.na(spec_df$turn_angle)] + 2*pi


############### Identify cosleeping ############### 

## save the x and y coordinates of the location where each baboon slept each night. I will calculate this location as the mean x and y coordinates of their location from 20:00 to 04:00

## set the minimum amount of fixes that need to be taken over night to for use to calculate the sleeping location
min_sleep_fixes <- 3

spec_df$sleep_x <- NA
spec_df$sleep_y <- NA

for(i in 1:length(tag_names)){ ## for each individual...
  
  # save this individual's data
  id_df <- spec_df[spec_df$id == tag_names[i],]
  
  # find the start and end day of the individual's data
  start <- min(id_df$day)
  end <- max(id_df$day)
  
  for(d in start:end){ # for each day on which the individual had data (each day between start and end of their data) 
    
    if( nrow( id_df[ id_df$day == d, ]) != 0 ){ ## if there is data for the day...
      
      # subset the data to the data between to 20:00 and 04:00
      sub <- id_df[ (id_df$day==d & id_df$time >= '20:00:00') | (id_df$day == (d+1) & id_df$time <= '04:00:00') ,]
      
      # if there are enough fixes during this night subset to robustly determine the sleep location
      if( nrow( sub ) >= min_sleep_fixes ){
        
        # take the mean of the x coordinate of the locations within this subset and save it as the x coordinate of the sleep site 
        id_df[ id_df$day == d, 'sleep_x' ] <- mean( sub$x, na.rm=T )
        
        # take the mean of the y coordinate of the locations within this subset and save it as the y coordinate of the sleep site
        id_df[ id_df$day == d, 'sleep_y'] <- mean( sub$y, na.rm=T )
        
      }
    }
  }
  
  # fill the total data frame with this individual's complete data
  spec_df[spec_df$id == tag_names[i], ] <- id_df
}

nas <- spec_df[ is.na( spec_df$sleep_x ), ]

unique( nas[, c( 'id', 'day' ) ] ) ### prints out the days we don't have sleep site info (but do have other info) for each baboon


## Use cluster analysis to find DISTINCT sleep sites ##
par(mfrow = c(1,1))

sleep_sites <- unique( spec_df[,c('sleep_x','sleep_y','id')] )

sleep_sites <- sleep_sites[!is.na(sleep_sites$sleep_x), ]

sleep_sites <- sleep_sites[order(sleep_sites$sleep_y), ] ### Orders the sleep site by latitudes so cluster number roughly refers to latitude later

plot(sleep_sites[, 1:2], col = as.factor(sleep_sites$id), asp=1, xlab = '',ylab = '',main='Sleep sites plotted by id') ## the points here are colored by id

legend('bottomright',legend = unique(sleep_sites$id), col = as.factor(unique(sleep_sites$id)),pch=1)

dist_matrix <- dist(sleep_sites[, 1 : 2])

clus_tree <- hclust(dist_matrix, method = 'complete')

cluster_membership <- cutree( tree = clus_tree, h = 300 )

sleep_sites$cluster <- cluster_membership

plot(sleep_sites$sleep_x,sleep_sites$sleep_y, col=(sleep_sites$cluster), asp=1, main="Sleep sites plotted by cluster", xlab='', ylab='')

#text( sleep_sites$sleep_x,sleep_sites$sleep_y, labels = sleep_sites$cluster)

### Determining the locations of the cluster. Using the mean y and x of all members of a cluster
cluster_site_y <- tapply(sleep_sites$sleep_y,INDEX = sleep_sites$cluster, FUN=mean)
cluster_site_x <- tapply(sleep_sites$sleep_x,INDEX = sleep_sites$cluster, FUN=mean)
uniq_sites <- data.frame(cluster = 1:length(cluster_site_y),y = cluster_site_y, x = cluster_site_x) ### This is a dataframe of each unique sleep site and its location

#### Use look-up tables to get this information back to spec_df ####
## first put the info into the sleep_sites dataframe, and then use this as the ultimate look-up table for spec_df
sleep_sites$clus.y <- sapply(sleep_sites$cluster, function(x) uniq_sites$y[match(x, uniq_sites$cluster)])
sleep_sites$clus.x <- sapply(sleep_sites$cluster, function(x) uniq_sites$x[match(x, uniq_sites$cluster)])

spec_df$clus_sleep_y <- sapply(spec_df$sleep_y,  function(x) sleep_sites$clus.y[match(x, sleep_sites$sleep_y)])
spec_df$clus_sleep_x <- sapply(spec_df$sleep_x,  function(x) sleep_sites$clus.x[match(x, sleep_sites$sleep_x)])
spec_df$sleep_clus <- sapply(spec_df$sleep_x,  function(x) sleep_sites$cluster[match(x, sleep_sites$sleep_x)])

# check to make sure that both individuals with the same group slept in the same sleep site every night
same_sleep_check <- aggregate( spec_df$sleep_clus, by = list( spec_df$group, spec_df$day ), FUN = function( x ) sum( !is.na( unique( x ) ) ) )

names( same_sleep_check ) <- c( 'group', 'day', 'num_sleep_sites')

same_sleep_check[ same_sleep_check$num_sleep_sites != 1, ]


####### time of arrival and departure time from sleep site ########

## find the first time (in the afternoon) that baboons entered within 200 m of the sleep location where the are going to sleep in -- this is the arrival time. Find the first time in the day that they exit a 200 m radius of their sleep location -- this is the departure time
spec_df$arrive_sleep_site <- NA
spec_df$leave_sleep_site <- NA

# set the distance 
dist_thresh <- 200

sum( spec_df$sleep_x < 100000, na.rm = T )

num_assigned_noon <- 0 

for( tag in tag_names ){
  
  idDat <- spec_df[spec_df$id == tag, ]
  
  if( sum( idDat$sleep_x < 100000, na.rm = T ) != 0 ) stop( 'arg' )
  
  for( day in unique( idDat$day ) ){
    
    dayDat <- idDat[ idDat$day == day, ]
    
    if( sum( dayDat$sleep_x < 100000, na.rm = T ) != 0 ) stop( 'arggg' )
    
    # I can't figure out why, but the y locations are becoming characters upon subsetting to dayDat but the x locations aren't (and the y locations aren't character in idDat or spec_df...??)
    dayDat$sleep_x <- as.numeric( dayDat$sleep_x )
    dayDat$sleep_y <- as.numeric( dayDat$sleep_y )
    
    idDat$sleep_x <- as.numeric( idDat$sleep_x )
    idDat$sleep_y <- as.numeric( idDat$sleep_y )
    
    if( sum( idDat$sleep_x < 100000, na.rm = T ) != 0 | sum( dayDat$sleep_x < 100000, na.rm = T ) != 0 ) stop( 'argggg' )
    
    if( sum( idDat[idDat$day == (day - 1), 'sleep_x'] < 100000, na.rm = T ) != 0 ) stop( 'argggggggg' )
    
    prev.sleep_x <- idDat[idDat$day == (day - 1), 'sleep_x'][1]
    prev.sleep_y <- idDat[idDat$day == (day - 1), 'sleep_y'][1]
    
    
    # if we know where they slept last night, determine when they leave the sleep site in the morning
    if( length( prev.sleep_x ) != 0 & !is.na( prev.sleep_x ) ){
      
      # find the distance of each location from the previous sleep location
      dist_from_prev <- sqrt( ( prev.sleep_y - dayDat$y_disc )**2 + ( prev.sleep_x - dayDat$x_disc )**2 )
      
      # if they never leave the distance threshold radius from their sleep location and that is not just because their data cuts out, declare their leave time 12:00:00
      if( sum( dist_from_prev > dist_thresh ) == 0 ){ # if they never leave the distance threshold
        
        # declare that they never left
        print( 'never left sleep site' )
        
        # if they actually had data from the middle of the day, when the baboons should be out foraging and stuff, then give them a departure time of 12:00:00, as the latest possible departure time. If they don't actually have data from the middle of the day, we will leave the departure time as NA
        if( nrow( dayDat[ dayDat$time > '11:00:00' & dayDat$time < '15:00:00', ] ) != 0 ){
          
          print( 'really never left sleep site' )
          
          # set the departure time as noon
          dayDat$leave_sleep_site <-  '12:00:00'
          
          num_assigned_noon <- num_assigned_noon + 1
          
        }
        
        # if the do leave the distance threshold radius of their sleep location as expected, save the first time in which they do this as the departure time  
      }else{
        
        # find the first instance in which they are greater than the distance threshold away from their previous sleep location
        ind_of_leave <- min( which( dist_from_prev > dist_thresh ) )
        
        # fill in the time of departure column with the time at which they left the previous sleep location
        dayDat$leave_sleep_site <- dayDat$time[ ind_of_leave ]
        
      }
      
    }
    
    # if we know where they are sleeping the coming night, determine the time of arrival at that sleep site
    if( sum( !is.na( dayDat$sleep_x ) ) != 0 ){
      
      # find the distance from each location to the sleep location for the coming night
      dist_from_next <- sqrt( ( dayDat$sleep_y - dayDat$y_disc )**2 + ( dayDat$sleep_x - dayDat$x_disc )**2 )
      
      # find the first instance in which they are closer than the distance threshold away from their next sleep location, conditioned on it being afternoon (because they often wake up in the morning where they will sleep the next night, and if we don't condition on it being in the afternoon, it will look like they arrived at the sleep site when they were actually waking up and leaving it in the morning). If they never leave the sleep site (or they arrive at it before noon), the time of arrival will be marked as noon
      ind_of_arrive <- min( which( dist_from_next < dist_thresh & dayDat$time > "12:00:00" ) )
      
      # fill in the time of arrival column with the time at which they arrived
      dayDat$arrive_sleep_site <- dayDat$time[ ind_of_arrive ]
      
      
    }
    
    
    # put this day of data back into the individual's data
    idDat[idDat$day == day, ] <- dayDat
  }
  
  # put this individual's data back into the full data
  spec_df[spec_df$id == tag, ] <- idDat
  
}

dev.off()

# make a dataframe with one row per day, so we have one observation per individual per day of when they left the previous sleep site and arrived at the next one
daily_dat <- unique( spec_df[ , c( 'id', 'day', 'arrive_sleep_site', 'leave_sleep_site' ) ] )

num_assigned_noon

sum( !is.na( daily_dat$leave_sleep_site ) )

as_hms( median( as_hms(  daily_dat$arrive_sleep_site ), na.rm = T ) )

as_hms( median( as_hms(  daily_dat$leave_sleep_site ), na.rm = T ) )

# make a histogram of the time of arrivals at the sleep site
x_min <- as.numeric( as_hms( '12:00:00' ) ) 
x_max <- as.numeric( as_hms( '24:00:00' ) ) 

axis_ticks <- seq( x_min, x_max, by = 60*60 )

hist( as.numeric( as_hms( daily_dat$arrive_sleep_site ) ), xlim = c( x_min, x_max ), xaxt = 'n', main = 'Time of arrival at sleep site (by definition cannot be earlier than noon)', xlab = 'Time of arrival' )

axis( 1, at = axis_ticks, labels = as_hms( axis_ticks ) )

# make a histogram of the time of departures from the sleep site
x_min <- as.numeric( as_hms( '00:00:00' ) ) 
x_max <- as.numeric( as_hms( '18:00:00' ) ) 

axis_ticks <- seq( x_min, x_max, by = 60*60 )

hist( as.numeric( as_hms( daily_dat$leave_sleep_site ) ), xlim = c( x_min, x_max ), xaxt = 'n', main = 'Time of departure from sleep site (if never departed, set to noon)', xlab = 'Time of departure' )

axis( 1, at = axis_ticks, labels = as_hms( axis_ticks ) )


##### more on cosleeping ####

# remove LI TH's data when she is dying. The day chosen here is the day on which she starts sleeping separate from the rest of her group
spec_df <- spec_df[ !( spec_df$id == 'Pa_LI_TH' & spec_df$local_timestamp >= as.POSIXct( '2014-05-28 00:00:00', tz = 'UTC') ), ]

## how many sleep sites are there? How many are shared?

uniq_sleeps <- unique( spec_df[ c( 'group', 'sleep_clus' ) ] )

uniq_sleeps <- uniq_sleeps[ !is.na( uniq_sleeps$sleep_clus ), ]

length( unique( uniq_sleeps$sleep_clus ) )

uniq_sleeps$sleep_clus[ duplicated( uniq_sleeps$sleep_clus ) ]

length( uniq_sleeps$sleep_clus[ duplicated( uniq_sleeps$sleep_clus ) ] )

## on how many occasions do groups co-sleep?

uniq_day_sleeps <- unique( spec_df[ c( 'group', 'day', 'sleep_clus' ) ] )

uniq_day_sleeps <- uniq_day_sleeps[ !is.na( uniq_day_sleeps$sleep_clus ), ]

cosleep_events <- uniq_day_sleeps[ duplicated( uniq_day_sleeps[ , c( 'day', 'sleep_clus' ) ] ) | duplicated( uniq_day_sleeps[ , c( 'day', 'sleep_clus' ) ], fromLast = T ) , ]

cosleep_events <- cosleep_events[ order( cosleep_events$day ), ]

nrow( unique( cosleep_events[ , c( 'day', 'sleep_clus' ) ] ) )

# make a matrix of how many times each group cosleeps with each other

cosleep_arr <- array( NA, dim = c( length( unique( spec_df$group ) ), length( unique( spec_df$group ) ) ), dimnames = list( sort( unique( spec_df$group ) ), sort( unique( spec_df$group ) ) ) )

group_names <- as.character( sort( unique( spec_df$group ) ) ) 

for( a in 1:( length( group_names ) - 1 ) ){
  
  for( b in ( a + 1 ):length( group_names ) ){
    
    pair_uniq_day_sleeps <- uniq_day_sleeps[ uniq_day_sleeps$group %in% group_names[ c( a, b ) ], ]
    
    cosleep_arr[ group_names[ a ], group_names[ b ] ] <- sum( duplicated( pair_uniq_day_sleeps[ , c( 'day', 'sleep_clus' ) ] ) )
  }
}

cosleep_arr

sum( cosleep_arr, na.rm = T )



## write the complete data into a dataframe
if(which_spec == "baboon"){
  write.csv(spec_df, "DATA/bab_complete.csv", row.names = F)
}else{
  if(which_spec == "vervet"){
    write.csv(spec_df, "DATA/verv_complete.csv", row.names = F)
  }else{
    if(which_spec == "leopard"){
      write.csv(spec_df, "DATA/leo_complete.csv", row.names = F)
    }
  }
}


