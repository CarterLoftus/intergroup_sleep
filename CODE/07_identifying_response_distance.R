
library( stringr )
library( data.table )
library( tidyr )


transp <- function(col, alpha=.5){
  res <- apply(col2rgb(col),2, function(c) rgb(c[1]/255, c[2]/255, c[3]/255, alpha))
  return(res)
}

##############  using permutations to find detection distance ##################

which_spec <- 'baboon'

# read in the data
if(which_spec == "baboon"){
  spec_df <- read.csv("DATA/bab_complete.csv")
}else{
  if(which_spec == "vervet"){
    spec_df <- read.csv("DATA/verv_complete.csv")
  }else{
    if(which_spec == "leopard"){
      spec_df <- read.csv("DATA/leo_complete.csv")
    }
  }
}

na_sum <- function( vec ){
  
  if( sum( !is.na( vec ) ) == 0){
    
    return( NA )
    
  }else{
    
    return( sum( vec, na.rm = T) )
  }
}


# make the timestamp into a POSIX element
spec_df$local_timestamp <- as.POSIXct( spec_df$local_timestamp, tz = 'UTC' )

# only keep the relevent columns
spec_df <- spec_df[ , ! (names( spec_df ) %in% c( 'lon', 'lat', 'group', 'spat.disc.dist', 'heading', 'turn_angle', 'UD', 'sleep_x', 'sleep_y', 'clus_sleep_y', 'clus_sleep_x', 'sleep_clus', 'rel.sleep.heading') ) ]

head( spec_df )

# set the number of permutations
n <- 1000

# set the maximum number of days that two paired days of data can be apart after permuting (chunk size)
perm_day_thresh <- 30

# set the proportion of days within a chunk of data that need to be there to be allow permutation of that data. We will be very conservative here and say there needs to be data every day 
min_to_perm <- 1

# set the size of the distance bins that we will test to see in which bin a group can first detect another group
bin_size <- 100

# set the maximum bin, beyond which we don't need to test for the possibility of detection
max_bin <- 7000

# save the time as its own column. This will be useful when creating fake timestamps during the permutations
spec_df$time <- str_split_fixed( spec_df$local_timestamp, ' ', 2)[ , 2 ]

# subset the data to just the data during the day
spec_df <- spec_df[ spec_df$time > '09:00:00' & spec_df$time < '17:00:00', ]

# remove LI TH and AI WG because they are already represented by their group-mates, whose data last longer. For this analysis, it is important that it remains only one individual's data and not a mixture of the two, because we are going to be looking at the probability of moving towards, away, and doing nothing, and so we need one individual's time series or it may look like a group is moving towards/away from another group only because it switched to tracking the other individual in the group
spec_df <- spec_df[ !spec_df$id %in% c( 'Pa_LI_TH', 'Pa_AI_WG' ), ]

tag_names <- as.character( unique( spec_df$id ) ) # save the tag names

# create empty vectors. These vectors will be filled with entries declaring which individuals need to be compared in the permutations
vec_a <- c()
vec_b <- c()

for( a in 1: length( tag_names ) ){
  
  for( b in 1: length( tag_names ) ){
    
    if( a != b ){
      
      # create vectors that represent the unique combinations that can be made of individuals/groups
      vec_a <- c( vec_a, tag_names[ a ] )
      
      vec_b <- c( vec_b, tag_names[ b ] )
      
    }
    
  }
}

# make a vector of the upper bound of each bin. The upper bound is exclusive, while the lower bound of each bin (which is just the upper bound - bin_size) is inclusive
dists <- seq( bin_size, max_bin, by = bin_size )

# set up the final dataframe that will be filled out to give the empirical probabilities of different actions at each distance bin between each dyad
emp_detect_dist <- data.frame( id_a = rep( vec_a, each = length( dists ) ), id_b = rep( vec_b, each = length( dists ) ), dist_bin = rep( dists, times = length( vec_a ) ), num_move_away = NA, num_move_towards = NA, num_do_nothing = NA, num_samples = NA, stringsAsFactors = F)

# set up the final dataframe that will be filled out to give the probabilities of different actions at each distance bin between each dyad produced by the permutations
perm_detect_dist <- emp_detect_dist[ rep( 1:nrow( emp_detect_dist ), times = n ), ]

# add a column to record the number of the permutation
perm_detect_dist$perm_n <- rep( 1:n, each = nrow( emp_detect_dist ) )

# create a dataframe that will contain the metadata for the permutations. I don't know if this is necessary, but I'll keep it
meta_spec_df <- data.frame( id_a = vec_a, id_b = vec_b, start_perm_at = rep( NA, length( vec_a ) ), end_perm_at =  rep( NA, length( vec_a ) ), stringsAsFactors = F )

# saves the first day of the study
start_date <- as.Date( min( spec_df$local_timestamp ) )

# fill in both the empirical and permutation data frames, one dyad at a time

set.seed( 500 )

for( a in 1:( length( tag_names ) - 1 ) ){
  
  for( b in ( a + 1 ): length( tag_names ) ){
    
    # subset the full gps dataframe to just including the dyad's data
    pair_spec_df <- spec_df[ spec_df$id %in% c( tag_names[ a ], tag_names[ b ] ), ]
    
    # just making sure 'id' is a character and not a factor. It messes things up if it is a factor
    pair_spec_df$id <- as.character( pair_spec_df$id ) 
    
    # make a column for the day since the start of the study
    pair_spec_df$day <- as.numeric( as.Date( pair_spec_df$local_timestamp )  -  start_date + 1 , units = 'days' )
    
    # trims the dataframe so it starts on the first day that both members of the dyad have data. We don't want to analyze anything before this
    pair_spec_df <- pair_spec_df[ pair_spec_df$day >=  max( aggregate( pair_spec_df$day , by=list( pair_spec_df$id ), FUN = min )[ , 2 ] ) , ] 
    
    # makes a column with both of their day columns starting at 1 on the first day when they both have the data. Important for chunking up the data in the next line
    pair_spec_df$temp_day <- as.numeric( pair_spec_df$day - min( pair_spec_df$day ) + 1, units = 'days' )
    
    # assign each row of data to a subset so that days of the data are only permuted within an range determined by perm_day_thresh
    pair_spec_df$chunk_num <- ceiling( pair_spec_df$temp_day / perm_day_thresh ) 
    
    # split up the data into the subsets created in the line above. Now we have a list of dataframes, which each dataframe corresponding to one time chunk within which permutation is allowable
    chunked_spec_df <- split( pair_spec_df, f = pair_spec_df$chunk_num )
    
    # create an empty dataframe that will represent a dyadic distance at every simultaneous fix during the current study period
    spec_df_real <- pair_spec_df[ 0, ]
    
    # perform n permutations of the dyad's dataset, building a permuted set of data each time and measuring the number of times an individual/group moved towards, away from, or did nothing with respect to the other dyad member
    for( i in 1:n ){
      
      if( i %% 100 == 0 ) print( i ) # progress tracker
      
      # create an empty dataframe that will represent a dyadic distance at every derived simultaneous fix of the study period that results after the permutation
      spec_df_perm <- data.frame( local_timestamp = character(), dist = numeric() )
      
      # loop through each chunk and permute the data within the chunks
      for( w in 1:length( chunked_spec_df ) ){ # For each chunk of data
        
        # save the chunk of data to the dataframe "chunk"
        chunk <- chunked_spec_df[[w]]
        
        # save the number of days of data that individuals a and b have in this chunk
        num_unique_days_a <- length( unique( chunk[ chunk$id == tag_names[ a ], 'day' ] ) ) 
        num_unique_days_b <- length( unique( chunk[ chunk$id == tag_names[ b ], 'day' ] ) ) 
        
        # if one of the individuals has no data for this chunk, skip the rest of the body of the loop and move to the next chunk
        if( num_unique_days_a == 0 || num_unique_days_b == 0 ){
          next
        }
        
        # if either individual a or b don't have the necessary number of days worth of data as determined by the min_to_perm parameter, skip the rest of the body of the loop
        if( num_unique_days_a < (min_to_perm * perm_day_thresh)  || num_unique_days_b < (min_to_perm * perm_day_thresh) ){ 
          next
        }
        
        # The first time through, we will make the empirical dataframe of dyadic distances and use that to determine when an individual moved towards/away from/did nothing with respect to the other dyad member. We only need to do this once
        if( i == 1 ){
          
          # add this to the running dataframe of dyadic distances over the whole study for this dyad (we are recreating the original dataframe chunk by chunk so that if chunks get left out of the permutation, they also aren't included in the empirical calculations)
          spec_df_real <- rbind( spec_df_real, chunk )
          
          # save the latest time that will be successfully permuted. For the first run through for each dyad, this will get updated every time the conditions above are surpassed such that a chunk of data is successfully permuted. Eventually this will serve to mark the end of successful permutations
          end_of_perm <- max(chunked_spec_df[[w]]$local_timestamp)
        }
        
        ### this is where you can decide whether to shift by days or shift by multiples of 15 mins
        
        # determine how many days to shift ID b's data by for the permutation
        ID_b_shifts <- sample( 1:( perm_day_thresh - 1 ), 1, replace = TRUE )
        
        # shift ID b's data by the amount determined above
        new_days_b <- chunk[ chunk$id == tag_names[ b ], 'day'] + ID_b_shifts
        
        # complicated line of code, but all it does is wrap the end of b's data back around to match the beginning of a's data. So if we shifted b's data by 3, a's data will still be 1, 2, 3, 4, 5, 6, 7; and b's data will be 5, 6, 7, 1, 2, 3, 4.
        new_days_b[ new_days_b > max( chunk[ chunk$id == tag_names[ b ], 'day' ] ) ] <- new_days_b[ new_days_b > max( chunk [ chunk$id == tag_names[ b ], 'day' ] ) ] - max( chunk[ chunk$id == tag_names[ b ], 'day' ] ) + min( chunk[ chunk$id == tag_names[ b ], 'day'] ) - 1
        
        # replace b's day data with these 'fake' shifted days
        chunk[ chunk$id == tag_names[ b ], 'day'] <- new_days_b
        
        # creates a fake date based on the fake day and the known start date of the study
        chunk$fake_date <- start_date + chunk$day - 1
        
        # creates a fake timestamp by appending the real time onto the fake date
        chunk$local_timestamp <- as.POSIXct(paste(chunk$fake_date, chunk$time, sep = ' '), tz = 'UTC')
        
        # add this to the running dataframe of permuted data over the whole study
        spec_df_perm <- rbind( spec_df_perm, chunk )
        
      }
      
      # if we are on the first run through the loop, save the empirical CA of the dyad
      if( i == 1 ){
        
        # make a wide format dataframe. Each row will be a timestamp and the columns will contain the x and y coordinates for both members of the dyad
        coords <- as.data.frame( ( data.table::dcast( as.data.table( spec_df_real ), local_timestamp ~ id, value.var = c( 'x_disc', 'y_disc'), crop = F ) ) )
        
        # save the index of the column with the x coordinates for id1
        id1_x_col <- which( grepl( tag_names[ a ], names( coords ) ) &  grepl( 'x', names( coords ) ) )
        
        # save the index of the column with the y coordinates for id1
        id1_y_col <- which( grepl( tag_names[ a ], names( coords ) ) &  grepl( 'y', names( coords ) ) )
        
        # save the index of the column with the x coordinates for id2
        id2_x_col <- which( grepl( tag_names[ b ], names( coords ) ) &  grepl( 'x', names( coords ) ) )
        
        # save the index of the column with the y coordinates for id2
        id2_y_col <- which( grepl( tag_names[ b ], names( coords ) ) &  grepl( 'y', names( coords ) ) )
        
        # save the x coordinates of id1 at time point 1 (unshifted time)
        x_1_1 <- coords[ , id1_x_col ]
        
        # save the x coordinates of id2 at time point 1 (unshifted time)
        x_2_1 <- coords[ , id2_x_col ]
        
        # save the y coordinates of id1 at time point 1 (unshifted time)
        y_1_1 <- coords[ , id1_y_col ]
        
        # save the y coordinates of id2 at time point 1 (unshifted time)
        y_2_1 <- coords[ , id2_y_col ]
        
        # calculate the dyadic distances at each timestamp
        coords$dy_dist <- sqrt( ( x_1_1 - x_2_1 )**2 + ( y_1_1 - y_2_1 )**2 )
        
        # go distance bin by distance bin to determine the number of each action taken at each distance bin
        for( dist in dists ){ # for each distance bin...
          
          coords_sub <- coords[ !is.na( coords$dy_dist ) & coords$dy_dist >= ( dist - bin_size ) & coords$dy_dist < dist , ] # subset to the just the data that occurs within this distance bin
          
          if( nrow( coords_sub ) != 0 ){ # if there is any data in this distance bin...
            
            # record the number of times that id1 was in this distance bin
            emp_detect_dist[ emp_detect_dist$id_a == tag_names[ a ] & emp_detect_dist$id_b == tag_names[ b ] & emp_detect_dist$dist_bin == dist, 'num_samples' ] <- nrow( coords_sub )
            
            # record the number of times that id2 was in this distance bin (should match id1's entry for this column above)
            emp_detect_dist[ emp_detect_dist$id_a == tag_names[ b ] & emp_detect_dist$id_b == tag_names[ a ] & emp_detect_dist$dist_bin == dist, 'num_samples' ] <- nrow( coords_sub )
          }
          
        }
        
        # add the metadata for this dyad to the metadata dataframe
        meta_spec_df[ meta_spec_df$id_a == tag_names[ a ] & meta_spec_df$id_b == tag_names[ b ], c( 'start_perm_at', 'end_perm_at' ) ] <- c( min( chunked_spec_df [[ 1 ]]$local_timestamp ), end_of_perm )
        
      }
      
      
      ## repeat the same as above but for the randomization
      
      # make a wide format dataframe. Each row will be a timestamp and the columns will contain the x and y coordinates for both members of the dyad
      coords <- as.data.frame( ( data.table::dcast( as.data.table( spec_df_perm ), local_timestamp ~ id, value.var = c( 'x_disc', 'y_disc'), crop = F ) ) )
      
      # save the index of the column with the x coordinates for id1
      id1_x_col <- which( grepl( tag_names[ a ], names( coords ) ) &  grepl( 'x', names( coords ) ) )
      
      # save the index of the column with the y coordinates for id1
      id1_y_col <- which( grepl( tag_names[ a ], names( coords ) ) &  grepl( 'y', names( coords ) ) )
      
      # save the index of the column with the x coordinates for id2
      id2_x_col <- which( grepl( tag_names[ b ], names( coords ) ) &  grepl( 'x', names( coords ) ) )
      
      # save the index of the column with the y coordinates for id2
      id2_y_col <- which( grepl( tag_names[ b ], names( coords ) ) &  grepl( 'y', names( coords ) ) )
      
      # save the x coordinates of id1 at time point 1 (unshifted time)
      x_1_1 <- coords[ , id1_x_col ]
      
      # save the x coordinates of id2 at time point 1 (unshifted time)
      x_2_1 <- coords[ , id2_x_col ]
      
      # save the y coordinates of id1 at time point 1 (unshifted time)
      y_1_1 <- coords[ , id1_y_col ]
      
      # save the y coordinates of id2 at time point 1 (unshifted time)
      y_2_1 <- coords[ , id2_y_col ]
      
      # calculate the dyadic distances at each timestamp
      coords$dy_dist <- sqrt( ( x_1_1 - x_2_1 )**2 + ( y_1_1 - y_2_1 )**2 )
      
      # go distance bin by distance bin to determine the number of each action taken at each distance bin
      for( dist in dists ){ # for each distance bin...
        
        coords_sub <- coords[ !is.na( coords$dy_dist ) & coords$dy_dist >= ( dist - bin_size ) & coords$dy_dist < dist , ] # subset to the just the data that occurs within this distance bin
        
        if( nrow( coords_sub ) != 0 ){ # if there is any data in this distance bin...
          
          # record the number of times that id1 was in this distance bin
          perm_detect_dist[ perm_detect_dist$id_a == tag_names[ a ] & perm_detect_dist$id_b == tag_names[ b ] & perm_detect_dist$dist_bin == dist & perm_detect_dist$perm_n == i, 'num_samples' ] <- nrow( coords_sub )
          
          # record the number of times that id2 was in this distance bin (should match id1's entry for this column above)
          perm_detect_dist[ perm_detect_dist$id_a == tag_names[ a ] & perm_detect_dist$id_b == tag_names[ b ] & perm_detect_dist$dist_bin == dist & perm_detect_dist$perm_n == i, 'num_samples' ] <- nrow( coords_sub )
          
        }
        
      }
      
    }
    
    print( paste( tag_names[ a ], tag_names[ b ] ) ) # progress tracker
    
  }
}


write.csv( emp_detect_dist, 'DATA/emp_detect_dist_100m_7000_seed_500.csv', row.names = F )

write.csv( perm_detect_dist, 'DATA/perm_detect_dist_100m_7000_seed_500.csv', row.names = F )


emp_detect_dist <- read.csv( 'DATA/emp_detect_dist_100m_7000_seed_500.csv' )
perm_detect_dist <- read.csv( 'DATA/perm_detect_dist_100m_7000_seed_500.csv' )


dev.off() 

emp_detect_dist$dyad_name <- apply( emp_detect_dist[, c( 'id_a', 'id_b' ) ], 1, function( x ) paste( sort( x ), collapse = '-' ) )
emp_detect_dist <- emp_detect_dist[ !duplicated( emp_detect_dist[ c( 'dist_bin', 'dyad_name') ] ), ]

perm_detect_dist$dyad_name <- paste( perm_detect_dist$id_a, perm_detect_dist$id_b, sep = '-' )

perm_detect_dist <- perm_detect_dist[ !is.na( perm_detect_dist$num_samples ), ]


emp_agg <- aggregate( emp_detect_dist$num_samples, by = list( emp_detect_dist$dist_bin ), FUN = sum, na.rm = T )
names( emp_agg ) <- c( 'dist_bin', 'num_samples' )

perm_agg <- aggregate( perm_detect_dist$num_samples, by = list( perm_detect_dist$dist_bin, perm_detect_dist$perm_n ), FUN = sum, na.rm = T )
names( perm_agg ) <- c( 'dist_bin', 'perm_n', 'num_samples' )



plot( perm_agg$dist_bin, perm_agg$num_samples, col = transp( 'black' , 0.01 ), pch = 16, cex = 1, main = '', bty = 'l', xlab = 'Distance bin', ylab = 'Number of samples within distance bin' ) 

points( emp_agg$dist_bin, emp_agg$num_samples, col = 'red', pch = 16, cex = 1.2 )



perm_agger_high <- aggregate( perm_agg$num_samples, by = list( perm_agg$dist_bin ), function( x ) quantile( x, 0.95 ) )
names( perm_agger_high ) <- c( 'dist_bin', 'upper_quant' )

perm_agger_low <- aggregate( perm_agg$num_samples, by = list( perm_agg$dist_bin ), function( x ) quantile( x, 0.05 ) )
names( perm_agger_low ) <- c( 'dist_bin', 'lower_quant' )


summary_stats_temp <- merge( x = emp_agg, y = perm_agger_high, by = 'dist_bin', all.x = T, all.y = T, sort = F )
summary_stats <- merge( x = summary_stats_temp, y = perm_agger_low, by = 'dist_bin', all.x = T, all.y = T, sort = F )

segments( summary_stats$dist_bin - 50, summary_stats$lower_quant, summary_stats$dist_bin + 50, summary_stats$lower_quant, col = '#4D4DFF', lwd = 3 )
segments( summary_stats$dist_bin - 50, summary_stats$upper_quant, summary_stats$dist_bin + 50, summary_stats$upper_quant, col = '#4D4DFF', lwd = 3 )

axis( 1, at = 600 )

abline( v = 650, lty = 2 )



summary_stats$sig <- as.numeric( summary_stats$upper_quant < summary_stats$num_samples)









