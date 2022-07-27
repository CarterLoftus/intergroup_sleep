


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
library( dtw )


########### Relative angles and distances to other groups ############

which_spec <- 'baboon'

if(which_spec == "baboon"){
  spec_df <- read.csv( "DATA/bab_complete.csv" )
}else{
  if(which_spec == "vervet"){
    spec_df <- read.csv( "DATA/verv_complete.csv" )
  }else{
    if(which_spec == "leopard"){
      spec_df <- read.csv( "DATA/leo_complete.csv" )
    }
  }
}


## make the timestamp into a POSIXct element
spec_df$local_timestamp <- as.POSIXct( spec_df$local_timestamp, tz = 'UTC' )

## save the names of the individuals
tag_names <- as.character( unique( spec_df$id ) ) 

# make a wide dataframe with a column for each individual's x location and y location
coors <- as.data.frame( ( data.table::dcast( as.data.table( spec_df ), local_timestamp ~ id, value.var = c( 'x_disc', 'y_disc'), crop = F ) ) )

head( coors )

## for the analysis below, it is important that all times are consecutive. So I am going to add timestamps for which there are no fixes, but that are during the study period, back into the dataframe
missing_times <- seq( min( spec_df$local_timestamp ), max( spec_df$local_timestamp ) , by = '15 mins' )[ ( ! seq( min( spec_df$local_timestamp ), max( spec_df$local_timestamp ) , by = '15 mins' ) %in% spec_df$local_timestamp ) ]

# make an empty dataframe that we will fill with the timestamps and empty entries where the coordinates would be
times_to_add <- data.frame( matrix( nrow = length( missing_times ), ncol = ncol( coors ) ) )

# give the dataframe the same column names as the wide format coordinates dataframe
names( times_to_add ) <- names( coors )

# add the missing times to the timestamp column of the dataframe 
times_to_add$local_timestamp <- missing_times

# add the dataframe with the missing times to the full dataframe
coors_full <- rbind( coors, times_to_add )

# reorder the dataframe by timestamp
coors_full <- coors_full[ order( coors_full$local_timestamp ), ]

# make sure that all the timestamps are in the dataframe
#print( sum( diff( coors_full$local_timestamp ) != 15 ) )

## make sure we have all the timestamps accounted for (again)
seq( min( spec_df$local_timestamp ), max( spec_df$local_timestamp ) , by = '15 mins' )[ ( ! seq( min( spec_df$local_timestamp ), max( spec_df$local_timestamp ) , by = '15 mins' ) %in% coors_full$local_timestamp ) ]

## create a dataframe that is going to store our dyadic data. We are going to repeat each individual's data 6 times, so that at each timestamp for which it has data it has one row per other individual in the data (actually including itself) in which it can store dyadic metrics at this timestamp
spec_df_full <- spec_df[ rep( seq_len( nrow( spec_df ) ), each = length( tag_names ) ), ]

## order the dataframe by the timestamp
spec_df_full <- spec_df_full[ order( spec_df_full$local_timestamp) , ]

## change the column names of id and group to be id1 and group1. This will be our focal individual
names( spec_df_full )[ names( spec_df_full ) == 'id' ] <- 'id1'
names( spec_df_full )[ names( spec_df_full ) == 'group' ] <- 'group1'

## create a column for the names of the non-focal dyad member
spec_df_full$id2 <- tag_names

## the distance changes saved here represent how the distances change in the next step. In other words, if two groups are close now, how does the focal group respond (i.e. what does it do next)? I figured that was more interesting than looking at how their distances changed in the previous step -- that would be more along the lines of, what did the groups do to come into contact, rather than what did the groups do in response to coming into contact
info_to_extract <- c( 'dyadDist', 'traj_dist', 'dtw_dist' )

# create an empty dataframe to add to the full data
cols_to_add <- as.data.frame( matrix( NA, nrow = nrow( spec_df_full ), ncol = length( info_to_extract ) ) )

# add the names of the columns we will add to the full data
names( cols_to_add ) <- info_to_extract  

# add these empty columns to the full data
spec_df_fuller <- cbind( spec_df_full, cols_to_add )

dtw_window <- 5

# calculate the information above (see info_to_extract) for each focal animal - non-focal animal dyad pair for the entire dataset
for(a in 1:length(tag_names)){ # for each pair of individuals
  
  for(b in 1:length(tag_names)){
    
    if( a != b ) { # if the dyad members are not the same individual...
      
      ## pull out the x and y location columns associated with each dyad member
      id1_x_col <- which( grepl( tag_names[ a ], names( coors_full ) ) &  grepl( 'x', names( coors_full ) ) )
      
      id1_y_col <- which( grepl( tag_names[ a ], names( coors_full ) ) &  grepl( 'y', names( coors_full ) ) )
      
      id2_x_col <- which( grepl( tag_names[ b ], names( coors_full ) ) &  grepl( 'x', names( coors_full ) ) )
      
      id2_y_col <- which( grepl( tag_names[ b ], names( coors_full ) ) &  grepl( 'y', names( coors_full ) ) )
      
      ## save the time-unshifted coordinates of each individual
      x_1_1 <- coors_full[ , id1_x_col ]
      
      x_2_1 <- coors_full[ , id2_x_col ]
      
      y_1_1 <- coors_full[ , id1_y_col ]
      
      y_2_1 <- coors_full[ , id2_y_col ]
      
      
      ## save the trajectories of each individual (i.e. the change in x and the change in y coordinates) (we will use these to see when cohesive movement starts and ends so we want the timestamp to match up with the step that is to come -- this will help us identify when the first similar step between two individuals is taken)
      x_change_1 <- c( diff( x_1_1 ), NA )
      x_change_2 <- c( diff( x_2_1 ), NA )
      
      y_change_1 <- c( diff( y_1_1 ), NA )
      y_change_2 <- c( diff( y_2_1 ), NA )
      
      ## send up a flag in case there are two individuals that have the same x or y coordinates. It just seems suspicious if they do
      if( sum( x_1_1 == x_2_1, na.rm = T ) != 0 | sum( y_1_1 == y_2_1, na.rm = T ) != 0 ) print( paste( tag_names[ a ], tag_names[ b ], coors_full$local_timestamp[ c( which(  y_1_1 == y_2_1 ), which(  y_1_1 == y_2_1 ) ) ] ) )
      
      ## calculate the dyadic distance at this point in time
      dyad_dist <- sqrt( ( x_1_1 - x_2_1 )**2 + ( y_1_1 - y_2_1 )**2 )
      
      ## find the timestamps associated with the dyadic distances that aren't NAs
      times_to_match <- coors_full$local_timestamp[ !is.na( dyad_dist ) ]
      
      ## insert the dyadic distances that aren't NAs into the full dataframe by lining up the timestamps
      spec_df_fuller[ spec_df_fuller$id1 == tag_names[ a ] & spec_df_fuller$id2 == tag_names[ b ], 'dyadDist' ][ match( times_to_match, spec_df_fuller[ spec_df_fuller$id1 == tag_names[ a ] & spec_df_fuller$id2 == tag_names[ b ], 'local_timestamp' ] ) ] <- dyad_dist[ !is.na( dyad_dist ) ]
      
      ## just a double check to insure that the subsetting above is lining up the timestamps correctly between the data we want to add to the dataframe and the dataframe itself
      sum( spec_df_fuller[ spec_df_fuller$id1 == tag_names[ a ] & spec_df_fuller$id2 == tag_names[ b ], 'local_timestamp'][ match( times_to_match, spec_df_fuller[ spec_df_fuller$id1 == tag_names[ a ] & spec_df_fuller$id2 == tag_names[ b ], 'local_timestamp' ] ) ] != times_to_match )
      
      
      
      # calculate the trajectory similarity (well technically dissimilarity, because more similar trajectories are indicated by lower numbers). This is the Euclidean distance between the dyad member's change in x coordinates and change in y coordinates. This will be a symmetrical measure for the dyad
      traj_dissim <- sqrt( ( x_change_1 - x_change_2 )**2 + ( y_change_1 - y_change_2 )**2 )
      
      # match the timestamps to insert the data into the dataframe
      times_to_match_5 <- coors_full$local_timestamp[ !is.na( traj_dissim ) ]
      
      # insert the data into the dataframe
      spec_df_fuller[ spec_df_fuller$id1 == tag_names[ a ] & spec_df_fuller$id2 == tag_names[ b ], 'traj_dist' ][ match( times_to_match_5, spec_df_fuller[ spec_df_fuller$id1 == tag_names[ a ] & spec_df_fuller$id2 == tag_names[ b ], 'local_timestamp' ] ) ] <- traj_dissim[ !is.na( traj_dissim ) ]
      
      
      ### finding the dynamic time warping distance of the movement trajectories with a sliding time window
      
      dtw_distance <- c()
      
      traj_1 <- matrix( c( x_change_1, y_change_1 ), ncol = 2, byrow = F )
      traj_2 <- matrix( c( x_change_2, y_change_2 ), ncol = 2, byrow = F )
      
      for( inde in 1:( nrow( coors_full ) - dtw_window ) ){
        
        ## if both individuals have data to calculate the dtw distance between their trajectories over the dtw sliding window
        if( sum( is.na( x_change_1[ inde : ( inde + dtw_window ) ] ) ) == 0 & sum( is.na( x_change_2[ inde : ( inde + dtw_window ) ] ) ) == 0 ){
          
          dtw_obj <- dtw( traj_1[ inde : ( inde + dtw_window ), ], traj_2[ inde : ( inde + dtw_window ), ] )
          
          dtw_distance <- c( dtw_distance, dtw_obj$distance )
          
        }else{
          
          dtw_distance <- c( dtw_distance, NA )
          
        }
        
      } 
      
      dtw_distance <- c( dtw_distance, rep( NA, length = dtw_window ) )
      
      # match the timestamps to insert the data into the dataframe
      times_to_match_6 <- coors_full$local_timestamp[ !is.na( dtw_distance ) ]
      
      # insert the data into the dataframe
      spec_df_fuller[ spec_df_fuller$id1 == tag_names[ a ] & spec_df_fuller$id2 == tag_names[ b ], 'dtw_dist' ][ match( times_to_match_6, spec_df_fuller[ spec_df_fuller$id1 == tag_names[ a ] & spec_df_fuller$id2 == tag_names[ b ], 'local_timestamp' ] ) ] <- dtw_distance[ !is.na( dtw_distance ) ]
      
    }
  }
}

### positive values for the changes and move columns imply that an individual is moving away from the other individual, and negative values imply that the focal individual is moving towards the other individual (id2)


# save the name of the group of the other dyad member
spec_df_fuller$group2 <- str_split_fixed( spec_df_fuller$id2, '_', 3 )[ , 2 ]


## write the new dataframe into a csv

if(which_spec == "baboon"){
  write.csv(spec_df_fuller, "DATA/bab_dyad.csv", row.names = F)
}else{
  if(which_spec == "vervet"){
    write.csv(spec_df_fuller, "DATA/verv_dyad.csv", row.names = F)
  }else{
    if(which_spec == "leopard"){
      write.csv(spec_df_fuller, "DATA/leo_dyad.csv", row.names = F)
    }
  }
}
