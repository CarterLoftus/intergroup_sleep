


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


#### na_sum is an alternative to the sum function. It is a function that returns the sum of a vector, but it will return NA (not 0) if all of the entries in the vector are NA
# vec is a vector of numbers

na_sum <- function( vec ){
  
  if( sum( !is.na( vec ) ) == 0){
    
    return( NA )
    
  }else{
    
    return( sum( vec, na.rm = T) )
  }
}



which_spec <- "baboon"

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

# make the timestamp into a POSIX element
spec_df$local_timestamp <- as.POSIXct( spec_df$local_timestamp, tz = 'UTC' )


# set the data to use for the permutations. This can either be 'full', 'day', or 'night'
data_to_use <- 'day'

# determine the number of permutations to do
n <- 1000

# set the distance threshold for what is considered an encounter
dist_thresh <- 600

# set the threshold of the number of days beyond which we won't swap data (due to potential effects of seasonality)
perm_day_thresh <- 30

# set the proportion of days that need to be represented in a chunk of data to allow its permutation. We will be strict here and start with 1
min_to_perm <- 1

# declare whether the permutation shift should occur by units of a day or simply by shifting the rows a permom number of times (so here that would be units of 15 mins)
shift_by_days <- FALSE

# set a spatial threshold for use when evaluating when groups moved towards or away from other group (or changed the overall dyadic distance). Any change below this distance will not be considered a real change
spat_thresh <- 30

# set the sampling interval
samp_int <- 15*60


## trim the data to the day or night, according to the data_to_use input above
if( data_to_use == 'day' ){
  
  spec_df <- spec_df[ spec_df$time > "09:00:00" & spec_df$time < "17:00:00", ]
  
}

if( data_to_use == 'night' ){
  
  spec_df <- spec_df[ spec_df$time < "06:00:00" & spec_df$time < "18:45:00", ]
  
  spec_df$local_timestamp <- spec_df$local_timestamp - 12*60*60 # subtract 12 hours so that during the a night is not split in the middle when it is wrapped back to the front of the data after having a shift by a permom number of days applied
}

# remove the data for LI TH and AI WG, because we will be calculating metrics such as headings of movement towards other dyad member, probability of moving closer or moving farther away, and for these it is important that are measuring only from one individual
spec_df <- spec_df[ !spec_df$id %in% c( 'Pa_LI_TH', 'Pa_AI_WG' ), ]

# only keep the relevent columns
spec_df <- spec_df[ , ! (names( spec_df ) %in% c( 'lon', 'lat', 'UD', 'interp_dist', 'sleep_x', 'sleep_y', 'clus_sleep_y', 'clus_sleep_x', 'sleep_clus', 'rel.sleep.heading', 'interp_x', 'interp_y', 'interp_missing', 'interp_error', 'interp', 'arrive_sleep_site', 'leave_sleep_site', 'dist_from_HR_edge' ) ) ]

# save the tag names
tag_names <- as.character( unique( spec_df$id ) )

# create empty vectors. These vectors will be filled with entries that eventually be put into the final dataframe that declares which individuals need to be compared
vec_a <- c()
vec_b <- c()


for( a in 1:( length( tag_names ) - 1 ) ){
  
  for( b in ( a + 1 ): length( tag_names ) ){
    
    # create vectors that represent the unique combinations that can be made of individuals in this category
    vec_a <- c( vec_a, tag_names[ a ] )
    
    vec_b <- c( vec_b, tag_names[ b ] )
    
  }
}


# set up the final dataframe that will be filled out to give the empirical metrics of interaction between the dyads
real_final <- data.frame( id_a = vec_a, id_b = vec_b, sum_dyad_dist = NA, num_within_dist_thresh = NA, total_dyad_dist_samples = NA, sum_dyad_dist_within_dist_thresh = NA, sum_heading_diff_in_dist_thresh = NA, length_heading_diff_in_dist_thresh = NA, sum_step_length_diff_in_dist_thresh = NA, length_step_length_diff_in_dist_thresh = NA,   sum_dtw_dist_in_dist_thresh = NA, length_dtw_dist_in_dist_thresh = NA, stringsAsFactors = F )

# set up the final dataframe that will be filled out to give metrics of interaction between the dyads that are produced by the permutations
perm_final <- real_final[ rep( 1:nrow( real_final ), times = n ), ]

perm_final$perm_n <- rep( 1:n, each = nrow( real_final ) )

# create a dataframe that will contain the metadata for the permutations
meta_spec_df <- data.frame( id_a = vec_a, id_b = vec_b, start_perm_at = rep( NA, length( vec_a ) ), end_perm_at =  rep( NA, length( vec_a ) ), stringsAsFactors = F )

# saves the first day of the study
start_date <- as.Date( min( spec_df$local_timestamp ) )

# save the time as its own column. This will be useful when creating fake timestamps during the permutations
spec_df$time <- str_split_fixed( spec_df$local_timestamp, ' ', 2)[ , 2 ]

for( row in 1:nrow( real_final) ){

  print( row / nrow(real_final) )

  # subset the full gps dataframe to just including the individual dyad's data during the correct period. This appropriate combination is determined by the set up of the dataframes above
  pair_spec_df <- spec_df[ spec_df$id %in% c( real_final[ row, c('id_a', 'id_b') ] ), ]

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

  # split up the data into the subsets created in the line above. Now we have a list of dataframes, which each dataframe corresponding to one time chunk within permutation is allowable
  chunked_spec_df <- split( pair_spec_df, f = pair_spec_df$chunk_num )

  # create an empty dataframe that will represent the real locations and data for the individuals during the current study period (and will match the permutations in which time chunks are included/excluded)
  spec_df_real <- pair_spec_df[ 0, ]

  # perform n permutations of the dyad's dataset
  for( i in 1:n ){

    if( i %% 50 == 0 ){ print( i ) } # progress tracker

    # create an empty dataframe that will represent a dyadic distance at every derived simultaneous fix of the study period that results after the permutation
    spec_df_perm <- pair_spec_df[ 0, ]

    # loop through each chunk and permute the data within the chunks
    for( w in 1:length( chunked_spec_df ) ){ # For each chunk of data

      # save the chunk of data to the dataframe "chunk"
      chunk <- chunked_spec_df[[w]]

      # save the number of days of data that individuals a and b have in this chunk
      num_unique_days_a <- length( unique( chunk[ chunk$id == real_final$id_a[ row ], 'day' ] ) )
      num_unique_days_b <- length( unique( chunk[ chunk$id == real_final$id_b[ row ], 'day' ] ) )

      # if one of the individuals has no data for this chunk, skip the rest of the body of the loop and move to the next chunk
      if( num_unique_days_a == 0 || num_unique_days_b == 0 ){
        next
      }

      # if either individual a or b don't have the necessary number of days worth of data as determined by the min_to_perm parameter, skip the rest of the body of the loop
      if( num_unique_days_a < (min_to_perm * perm_day_thresh)  || num_unique_days_b < (min_to_perm * perm_day_thresh) ){
        next
      }

      # The first time through, we will make the empirical dataframe of dyadic distances. We only need to do this once
      if( i == 1 ){

        # add this to the running dataframe of dyadic distances over the whole study for this dyad (we are recreating the original dataframe chunk by chunk so that if chunks get left out of the permutation, they also aren't included in the empirical calculations)
        spec_df_real <- rbind( spec_df_real, chunk )

        # save the latest time that will be successfully permuted. For the first run through for each individual dyad, this will get updated every time the conditions above are surpassed such that a chunk of data is successfully permuted. Eventually this will serve to mark the end of successful permutations
        end_of_perm <- max(chunked_spec_df[[w]]$local_timestamp)

      }

      if( shift_by_days == T ){

        # determine how many days to shift ID b's data by for the permutation
        ID_b_shifts <- sample( 1:( perm_day_thresh - 1 ), 1, replace = TRUE )

        # shift ID b's data by the amount determined above
        new_days_b <- chunk[ chunk$id == real_final$id_b[ row ], 'day'] + ID_b_shifts

        # complicated line of code, but all it does is wrap the end of b's data back around to match the beginning of a's data. So if we shifted b's data by 3, a's data will still be 1, 2, 3, 4, 5, 6, 7; and b's data will be 5, 6, 7, 1, 2, 3, 4.
        new_days_b[ new_days_b > max( chunk[ chunk$id == real_final$id_b[ row ], 'day' ] ) ] <- new_days_b[ new_days_b > max( chunk [ chunk$id == real_final$id_b[ row ], 'day' ] ) ] - max( chunk[ chunk$id == real_final$id_b[ row ], 'day' ] ) + min( chunk[ chunk$id == real_final$id_b[ row ], 'day'] ) - 1

        # replace b's day data with these 'fake' shifted days
        chunk[ chunk$id == real_final$id_b[ row ], 'day'] <- new_days_b

        # creates a fake date based on the fake day and the known start date of the study
        chunk$fake_date <- start_date + chunk$day - 1

        # creates a fake timestamp by appending the real time onto the fake date
        chunk$local_timestamp <- as.POSIXct(paste( chunk$fake_date, chunk$time, sep = ' '), tz = 'UTC')


      }else{

        chunk_b <- chunk[ chunk$id == real_final$id_b[ row ], ]

        ID_b_shifts <- sample( 1: ( nrow( chunk_b ) - 1), 1, replace = TRUE )

        new_row_inds <- seq_len( nrow( chunk_b ) ) + ID_b_shifts

        new_row_inds [ new_row_inds > nrow( chunk_b ) ] <- new_row_inds [ new_row_inds > nrow( chunk_b ) ] - nrow( chunk_b )

        chunk_b$local_timestamp <- chunk_b$local_timestamp[ new_row_inds ]

        chunk[ chunk$id == real_final$id_b[ row ], ] <- chunk_b

      }

      # add this to the running dataframe of permuted data over the whole study
      spec_df_perm <- rbind( spec_df_perm, chunk )

    }

    # if we are on the first run through the loop, save the empirical CA of the dyad
    if( i == 1 ){

           # use the dyad_dist function to calculate the dyadic distance for this pair of individuals
      wide_dat <- as.data.frame( dcast( as.data.table( spec_df_real ), local_timestamp ~ id, value.var = c( 'x_disc', 'y_disc', 'spat.disc.dist_next', 'heading' ), drop = F ) )

      # save the index of the column with the x coordinates for id1
      id1_x_col <- which( grepl( real_final$id_a[ row ], names( wide_dat ) ) &  grepl( 'x_disc', names( wide_dat ) ) )

      # save the index of the column with the y coordinates for id1
      id1_y_col <- which( grepl( real_final$id_a[ row ], names( wide_dat ) ) &  grepl( 'y_disc', names( wide_dat ) ) )

      # save the index of the column with the x coordinates for id2
      id2_x_col <- which( grepl( real_final$id_b[ row ], names( wide_dat ) ) &  grepl( 'x_disc', names( wide_dat ) ) )

      # save the index of the column with the y coordinates for id2
      id2_y_col <- which( grepl( real_final$id_b[ row ], names( wide_dat ) ) &  grepl( 'y_disc', names( wide_dat ) ) )

      # save the x coordinates of id1 at time point 1 (unshifted time)
      x_1_1 <- wide_dat[ , id1_x_col ]

      # save the x coordinates of id2 at time point 1 (unshifted time)
      x_2_1 <- wide_dat[ , id2_x_col ]

      # save the y coordinates of id1 at time point 1 (unshifted time)
      y_1_1 <- wide_dat[ , id1_y_col ]

      # save the y coordinates of id2 at time point 1 (unshifted time)
      y_2_1 <- wide_dat[ , id2_y_col ]

      # calculate the dyadic distances at each timestamp
      wide_dat$dy_dist <- sqrt( ( x_1_1 - x_2_1 )**2 + ( y_1_1 - y_2_1 )**2 )

      # add the sum total of the dyadic distances to the final dataframe
      real_final[ row, 'sum_dyad_dist' ] <- sum( wide_dat$dy_dist, na.rm = T )

      # add the number of instances within the distance threshold to the dataframe
      real_final[ row, 'num_within_dist_thresh' ] <- sum( wide_dat$dy_dist <= dist_thresh, na.rm = T )

      # add the total number of dyadic distance samples present for this dyad to the dataframe
      real_final[ row, 'total_dyad_dist_samples' ] <- sum( !is.na( wide_dat$dy_dist ) )

      # add the sum of the dyadic distances within the distance threshold to the final dataframe
      real_final[ row, 'sum_dyad_dist_within_dist_thresh' ] <- na_sum( wide_dat$dy_dist[ which( wide_dat$dy_dist <= dist_thresh ) ] )

      # save the index of the column indicating id_a's heading
      id_a_heading_col <- which( grepl( real_final$id_a[ row ], names( wide_dat ) ) &  grepl( 'heading', names( wide_dat ) ) )

      # save the index of the column indicating id_b's heading
      id_b_heading_col <- which( grepl( real_final$id_b[ row ], names( wide_dat ) ) &  grepl( 'heading', names( wide_dat ) ) )

      # calculate the difference in headings between id_a and id_b while in the distance threshold
      angle_diff <- wide_dat[ , id_a_heading_col ][ which( wide_dat$dy_dist <= dist_thresh ) ] - wide_dat[ , id_b_heading_col ][ which( wide_dat$dy_dist <= dist_thresh ) ]

      # correct for the fact that 0 and 2 pi are the same
      angle_diff[ angle_diff > pi & !is.na(angle_diff) ] <- angle_diff[ angle_diff > pi & !is.na(angle_diff)] - 2*pi

      angle_diff[ angle_diff < - pi & !is.na(angle_diff) ] <- angle_diff[ angle_diff < - pi & !is.na(angle_diff)] + 2*pi

      # save the sum of the difference in headings within the distance threshold in the final dataframe
      real_final[ row, 'sum_heading_diff_in_dist_thresh' ] <- na_sum( abs( angle_diff ) )

      # save the number of samples used to calculate the sum of the difference in headings within the distance threshold in the final dataframe
      real_final[ row, 'length_heading_diff_in_dist_thresh' ] <- sum( !is.na( abs( angle_diff ) ) )

      # save the index of the column indicating id_a's heading
      id_a_step_col <- which( grepl( real_final$id_a[ row ], names( wide_dat ) ) &  grepl( 'spat.disc.dist', names( wide_dat ) ) )

      # save the index of the column indicating id_b's heading
      id_b_step_col <- which( grepl( real_final$id_b[ row ], names( wide_dat ) ) &  grepl( 'spat.disc.dist', names( wide_dat ) ) )

      # extract id_a's step lengths within the distance threshold
      id_a_step_lengths <- wide_dat[ , id_a_step_col ][ which( wide_dat$dy_dist <= dist_thresh ) ]

      # extract id_b's step lengths within the distance threshold
      id_b_step_lengths <- wide_dat[ , id_b_step_col ][ which( wide_dat$dy_dist <= dist_thresh ) ]

      # save the sum of the difference in step lengths between id_a and id_b while they are in the distance threshold to the final dataframe
      real_final[ row, 'sum_step_length_diff_in_dist_thresh' ] <- na_sum( abs( id_a_step_lengths - id_b_step_lengths ) )

      # save the number of sampls used to calculate the sum of the difference in step lengths between id_a and id_b while they are in the distance threshold to the final dataframe
      real_final[ row, 'length_step_length_diff_in_dist_thresh' ] <- sum( !is.na( abs( id_a_step_lengths - id_b_step_lengths ) ) )

      # calculate the change in the x coordinate at each timestamp for id_a (on a row it will represent the change in the x coordinate to get to the next row)
      x_change_A <- c( diff( wide_dat[ , id1_x_col ] ), NA )

      # calculate the change in the y coordinate at each timestamp for id_a (on a row it will represent the change in the y coordinate to get to the next row)
      y_change_A <- c( diff( wide_dat[ , id1_y_col ] ), NA )

      # calculate the change in the x coordinate at each timestamp for id_b (on a row it will represent the change in the x coordinate to get to the next row)
      x_change_B <- c( diff( wide_dat[ , id2_x_col ] ), NA )

      # calculate the change in the y coordinate at each timestamp for id_b (on a row it will represent the change in the y coordinate to get to the next row)
      y_change_B <- c( diff( wide_dat[ , id2_y_col ] ), NA )

      # if there's only one instance within the distance threshold, the below code will make a vector instead of a matrix of coordinates unless drop = F in the subsetting

      # make a matrix representing A's trajectory while within the distance threshold
      A_traj <- matrix( c( x_change_A, y_change_A ), ncol = 2, byrow = F )[ which( wide_dat$dy_dist <= dist_thresh ), , drop = F ]

      # make a matrix representing B's trajectory while within the distance threshold
      B_traj <- matrix( c( x_change_B, y_change_B ), ncol = 2, byrow = F )[ which( wide_dat$dy_dist <= dist_thresh ), , drop = F ]

      # remove NA's from A's trajectory
      A_traj_trim <- A_traj[ complete.cases( A_traj ) & complete.cases( B_traj ), , drop = F ]

      # remove NA's from A's trajectory
      B_traj_trim <- B_traj[ complete.cases( A_traj ) & complete.cases( B_traj ), , drop = F ]

      if( nrow( A_traj_trim  ) != 0 ){

        # calculate the dtw distance between the two trajectories
        dtw_traj <- dtw( x = A_traj_trim, y = B_traj_trim )

        # save the dtw distance between the trajectories while within the distance threshold to the final dataframe
        real_final[ row, 'sum_dtw_dist_in_dist_thresh' ] <- dtw_traj$distance

      }

      # save the length of the trajectory for which we calculated the dtw distance to the final dataframe
      real_final[ row, 'length_dtw_dist_in_dist_thresh' ] <- nrow( A_traj_trim )

      meta_spec_df[ row, c( 'start_perm_at', 'end_perm_at' ) ] <- as.character( c( min( chunked_spec_df[[ 1 ]]$local_timestamp ), end_of_perm ) )

    }
    
    
    # use the dyad_dist function to calculate the dyadic distance for this pair of individuals
    wide_dat <- as.data.frame( dcast( as.data.table( spec_df_perm ), local_timestamp ~ id, value.var = c( 'x_disc', 'y_disc', 'spat.disc.dist_next', 'heading' ), drop = F ) )

    # save the index of the column with the x coordinates for id1
    id1_x_col <- which( grepl( perm_final$id_a[ perm_final$id_a == real_final$id_a[ row ] & perm_final$id_b == real_final$id_b[ row ] & perm_final$perm_n == i ], names( wide_dat ) ) &  grepl( 'x_disc', names( wide_dat ) ) )

    # save the index of the column with the y coordinates for id1
    id1_y_col <- which( grepl( perm_final$id_a[ perm_final$id_a == real_final$id_a[ row ] & perm_final$id_b == real_final$id_b[ row ] & perm_final$perm_n == i ], names( wide_dat ) ) &  grepl( 'y_disc', names( wide_dat ) ) )

    # save the index of the column with the x coordinates for id2
    id2_x_col <- which( grepl( perm_final$id_b[ perm_final$id_a == real_final$id_a[ row ] & perm_final$id_b == real_final$id_b[ row ] & perm_final$perm_n == i ], names( wide_dat ) ) &  grepl( 'x_disc', names( wide_dat ) ) )

    # save the index of the column with the y coordinates for id2
    id2_y_col <- which( grepl( perm_final$id_b[ perm_final$id_a == real_final$id_a[ row ] & perm_final$id_b == real_final$id_b[ row ] & perm_final$perm_n == i ], names( wide_dat ) ) &  grepl( 'y_disc', names( wide_dat ) ) )

    # save the x coordinates of id1 at time point 1 (unshifted time)
    x_1_1 <- wide_dat[ , id1_x_col ]

    # save the x coordinates of id2 at time point 1 (unshifted time)
    x_2_1 <- wide_dat[ , id2_x_col ]

    # save the y coordinates of id1 at time point 1 (unshifted time)
    y_1_1 <- wide_dat[ , id1_y_col ]

    # save the y coordinates of id2 at time point 1 (unshifted time)
    y_2_1 <- wide_dat[ , id2_y_col ]

    # calculate the dyadic distances at each timestamp
    wide_dat$dy_dist <- sqrt( ( x_1_1 - x_2_1 )**2 + ( y_1_1 - y_2_1 )**2 )

    # add the sum total of the dyadic distances to the final dataframe
    perm_final[ perm_final$id_a == real_final$id_a[ row ] & perm_final$id_b == real_final$id_b[ row ] & perm_final$perm_n == i, 'sum_dyad_dist' ] <- sum( wide_dat$dy_dist, na.rm = T )

    # add the number of instances within the distance threshold to the dataframe
    perm_final[ perm_final$id_a == real_final$id_a[ row ] & perm_final$id_b == real_final$id_b[ row ] & perm_final$perm_n == i, 'num_within_dist_thresh' ] <- sum( wide_dat$dy_dist <= dist_thresh, na.rm = T )

    # add the total number of dyadic distance samples present for this dyad to the dataframe
    perm_final[ perm_final$id_a == real_final$id_a[ row ] & perm_final$id_b == real_final$id_b[ row ] & perm_final$perm_n == i, 'total_dyad_dist_samples' ] <- sum( !is.na( wide_dat$dy_dist ) )

    # add the sum of the dyadic distances within the distance threshold to the final dataframe
    perm_final[ perm_final$id_a == real_final$id_a[ row ] & perm_final$id_b == real_final$id_b[ row ] & perm_final$perm_n == i, 'sum_dyad_dist_within_dist_thresh' ] <- na_sum( wide_dat$dy_dist[ which( wide_dat$dy_dist <= dist_thresh ) ] )

    # save the index of the column indicating id_a's heading
    id_a_heading_col <- which( grepl( perm_final$id_a[ perm_final$id_a == real_final$id_a[ row ] & perm_final$id_b == real_final$id_b[ row ] & perm_final$perm_n == i ], names( wide_dat ) ) &  grepl( 'heading', names( wide_dat ) ) )

    # save the index of the column indicating id_b's heading
    id_b_heading_col <- which( grepl( perm_final$id_b[ perm_final$id_a == real_final$id_a[ row ] & perm_final$id_b == real_final$id_b[ row ] & perm_final$perm_n == i ], names( wide_dat ) ) &  grepl( 'heading', names( wide_dat ) ) )

    # calculate the difference in headings between id_a and id_b while in the distance threshold
    angle_diff <- wide_dat[ , id_a_heading_col ][ which( wide_dat$dy_dist <= dist_thresh ) ] - wide_dat[ , id_b_heading_col ][ which( wide_dat$dy_dist <= dist_thresh ) ]

    # correct for the fact that 0 and 2 pi are the same
    angle_diff[ angle_diff > pi & !is.na(angle_diff) ] <- angle_diff[ angle_diff > pi & !is.na(angle_diff)] - 2*pi

    angle_diff[ angle_diff < - pi & !is.na(angle_diff) ] <- angle_diff[ angle_diff < - pi & !is.na(angle_diff)] + 2*pi

    # save the sum of the difference in headings within the distance threshold in the final dataframe
    perm_final[ perm_final$id_a == real_final$id_a[ row ] & perm_final$id_b == real_final$id_b[ row ] & perm_final$perm_n == i, 'sum_heading_diff_in_dist_thresh' ] <- na_sum( abs( angle_diff ) )

    # save the number of samples used to calculate the sum of the difference in headings within the distance threshold in the final dataframe
    perm_final[ perm_final$id_a == real_final$id_a[ row ] & perm_final$id_b == real_final$id_b[ row ] & perm_final$perm_n == i, 'length_heading_diff_in_dist_thresh' ] <- sum( !is.na( abs( angle_diff ) ) )

    # save the index of the column indicating id_a's heading
    id_a_step_col <- which( grepl( perm_final$id_a[ perm_final$id_a == real_final$id_a[ row ] & perm_final$id_b == real_final$id_b[ row ] & perm_final$perm_n == i ], names( wide_dat ) ) &  grepl( 'spat.disc.dist', names( wide_dat ) ) )

    # save the index of the column indicating id_b's heading
    id_b_step_col <- which( grepl( perm_final$id_b[ perm_final$id_a == real_final$id_a[ row ] & perm_final$id_b == real_final$id_b[ row ] & perm_final$perm_n == i ], names( wide_dat ) ) &  grepl( 'spat.disc.dist', names( wide_dat ) ) )

    # extract id_a's step lengths within the distance threshold
    id_a_step_lengths <- wide_dat[ , id_a_step_col ][ which( wide_dat$dy_dist <= dist_thresh ) ]

    # extract id_b's step lengths within the distance threshold
    id_b_step_lengths <- wide_dat[ , id_b_step_col ][ which( wide_dat$dy_dist <= dist_thresh ) ]

    # save the sum of the difference in step lengths between id_a and id_b while they are in the distance threshold to the final dataframe
    perm_final[ perm_final$id_a == real_final$id_a[ row ] & perm_final$id_b == real_final$id_b[ row ] & perm_final$perm_n == i, 'sum_step_length_diff_in_dist_thresh' ] <- na_sum( abs( id_a_step_lengths - id_b_step_lengths ) )

    # save the number of sampls used to calculate the sum of the difference in step lengths between id_a and id_b while they are in the distance threshold to the final dataframe
    perm_final[ perm_final$id_a == real_final$id_a[ row ] & perm_final$id_b == real_final$id_b[ row ] & perm_final$perm_n == i, 'length_step_length_diff_in_dist_thresh' ] <- sum( !is.na( abs( id_a_step_lengths - id_b_step_lengths ) ) )

    # calculate the change in the x coordinate at each timestamp for id_a (on a perm_final$id_a == real_final$id_a[ row ] & perm_final$id_b == real_final$id_b[ row ] & perm_final$perm_n == i it will represent the change in the x coordinate to get to the next perm_final$id_a == real_final$id_a[ row ] & perm_final$id_b == real_final$id_b[ row ] & perm_final$perm_n == i)
    x_change_A <- c( diff( wide_dat[ , id1_x_col ] ), NA )

    # calculate the change in the y coordinate at each timestamp for id_a (on a perm_final$id_a == real_final$id_a[ row ] & perm_final$id_b == real_final$id_b[ row ] & perm_final$perm_n == i it will represent the change in the x coordinate to get to the next perm_final$id_a == real_final$id_a[ row ] & perm_final$id_b == real_final$id_b[ row ] & perm_final$perm_n == i)
    y_change_A <- c( diff( wide_dat[ , id1_y_col ] ), NA )

    # calculate the change in the x coordinate at each timestamp for id_b (on a perm_final$id_a == real_final$id_a[ row ] & perm_final$id_b == real_final$id_b[ row ] & perm_final$perm_n == i it will represent the change in the x coordinate to get to the next perm_final$id_a == real_final$id_a[ row ] & perm_final$id_b == real_final$id_b[ row ] & perm_final$perm_n == i)
    x_change_B <- c( diff( wide_dat[ , id2_x_col ] ), NA )

    # calculate the change in the y coordinate at each timestamp for id_b (on a perm_final$id_a == real_final$id_a[ row ] & perm_final$id_b == real_final$id_b[ row ] & perm_final$perm_n == i it will represent the change in the x coordinate to get to the next perm_final$id_a == real_final$id_a[ row ] & perm_final$id_b == real_final$id_b[ row ] & perm_final$perm_n == i)
    y_change_B <- c( diff( wide_dat[ , id2_y_col ] ), NA )

    # make a matrix representing A's trajectory while within the distance threshold
    A_traj <- as.data.frame( matrix( c( x_change_A, y_change_A ), ncol = 2, byrow = F )[ which( wide_dat$dy_dist <= dist_thresh ), , drop = F ] )

    # make a matrix representing B's trajectory while within the distance threshold
    B_traj <- as.data.frame( matrix( c( x_change_B, y_change_B ), ncol = 2, byrow = F )[ which( wide_dat$dy_dist <= dist_thresh ), , drop = F ] )

    # remove NA's from A's trajectory
    A_traj_trim <- A_traj[ complete.cases( A_traj ) & complete.cases( B_traj ), , drop = F ]

    # remove NA's from A's trajectory
    B_traj_trim <- B_traj[ complete.cases( A_traj ) & complete.cases( B_traj ), , drop = F ]

    if( nrow( A_traj_trim ) != 0 ){

      # calculate the dtw distance between the two trajectories
      dtw_traj <- dtw( x = A_traj_trim, y = B_traj_trim )

      # save the dtw distance between the trajectories while within the distance threshold to the final dataframe
      perm_final[ perm_final$id_a == real_final$id_a[ row ] & perm_final$id_b == real_final$id_b[ row ] & perm_final$perm_n == i, 'sum_dtw_dist_in_dist_thresh' ] <- dtw_traj$distance

    }

    # save the length of the trajectory for which we calculated the dtw distance to the final dataframe
    perm_final[ perm_final$id_a == real_final$id_a[ row ] & perm_final$id_b == real_final$id_b[ row ] & perm_final$perm_n == i, 'length_dtw_dist_in_dist_thresh' ] <- nrow( A_traj_trim )

    
  }
  
}

meta_spec_df$start_perm_at <- as.POSIXct( meta_spec_df$start_perm_at, tz = 'UTC', origin ='1970-01-01 00:00:00' )

meta_spec_df$end_perm_at <- as.POSIXct( meta_spec_df$end_perm_at, tz = 'UTC', origin ='1970-01-01 00:00:00' )

dir.create( paste0( getwd(), '/RESULTS/movement_path_permutation_results' ) )

write.csv( real_final, paste0('RESULTS/movement_path_permutation_results/real_final_', data_to_use, '_by_day_', shift_by_days, '.csv' ), row.names = F )

write.csv( perm_final, paste0('RESULTS/movement_path_permutation_results/perm_final_', data_to_use, '_by_day_', shift_by_days, '.csv' ), row.names = F )

real_final <- read.csv( paste0('RESULTS/movement_path_permutation_results/real_final_', data_to_use, '_by_day_', shift_by_days, '.csv' ) )

perm_final <- read.csv( paste0('RESULTS/movement_path_permutation_results/perm_final_', data_to_use, '_by_day_', shift_by_days, '.csv' ) )


### mean dyadic distance 

perm_agg <- aggregate( perm_final[ , - c( 1, 2 ) ], by = list( perm_final$perm_n ), FUN = sum, na.rm = T )

names( perm_agg )[ names ( perm_agg ) == "Group.1" ] <- 'perm_n'

real_final$for_agg <- 1

real_agg <- aggregate( real_final[ , - c(1, 2) ], by = list( real_final$for_agg ) ,FUN = sum, na.rm = T )

dev.off()

real_mean_dy_dist <- real_agg$sum_dyad_dist / real_agg$total_dyad_dist_samples

perm_mean_dy_dist <- perm_agg$sum_dyad_dist / perm_agg$total_dyad_dist_samples

hist( perm_mean_dy_dist, xlab = 'Mean dyadic distance (m)', xlim = range( c( real_mean_dy_dist, perm_mean_dy_dist ) ), main = "All Dyads" )

abline( v = real_mean_dy_dist, col = 'red' )



dens <- density( perm_mean_dy_dist )

dens_x <- dens$x

dens_y <- dens$y

plot( dens_x, dens_y, type = 'l', bty = 'l', xlab = paste( 'Mean dyadic distance (m)' ), xlim = range( c( real_mean_dy_dist, perm_mean_dy_dist ) ), ylab = "Probability density" )

abline( v = real_mean_dy_dist, col = 'red' )

abline( v = quantile( perm_mean_dy_dist, c( 0.025, 0.975 ) ), col = 'black', lty = 3)


par( bg = 'black' )

dens <- density( perm_mean_dy_dist )

dens_x <- dens$x

dens_y <- dens$y

plot( dens_x, dens_y, type = 'l', bty = 'l', xlab = paste( 'Mean dyadic distance (m)' ), xlim = range( c( real_mean_dy_dist, perm_mean_dy_dist ) ), ylab = "Probability density", col = 'white', col.axis = 'white', col.lab = 'white', col.main = 'white' )

abline( v = real_mean_dy_dist, col = 'red' )

abline( v = quantile( perm_mean_dy_dist, c( 0.025, 0.975 ) ), col = 'white', lty = 3)

axis( 1, col = 'white', at = seq( 5790, 5850, by = 10 ), labels = rep( '', length( seq( 5790, 5850, by = 10 ) ) ) )
axis( 2, col = 'white', at = seq( 0, 0.04, by = 0.01 ), labels = rep( '', length( seq( 0, 0.04, by = 0.01 ) ) ) )



p <- sum( real_mean_dy_dist <= perm_mean_dy_dist ) / max( perm_agg$perm_n )

final_p <- ifelse( p > 0.5, sum( real_mean_dy_dist >= perm_mean_dy_dist ) / max( perm_agg$perm_n ), p )

#mtext( paste( 'p =', final_p ), side = 4 )


### proportion of time within the distance threshold 

dev.off()

real_for_comparison <- sum( real_final$num_within_dist_thresh ) / sum( real_final$total_dyad_dist_samples )

perm_for_comparison <- perm_agg$num_within_dist_thresh / perm_agg$total_dyad_dist_samples

hist( perm_for_comparison, xlab = paste( 'Proportion of time within', dist_thresh, '(m)' ), xlim = range( c( real_for_comparison, perm_for_comparison ) ), main = "" )


dens <- density( perm_for_comparison )

dens_x <- dens$x

dens_y <- dens$y

plot( dens_x, dens_y, type = 'l', bty = 'l', xlab = paste( 'Proportion of time engaged in an encounter' ), xlim = range( c( real_for_comparison, perm_for_comparison ) ), ylab = "Probability density" )

abline( v = real_for_comparison, col = 'red' )

abline( v = quantile( perm_for_comparison, c( 0.025, 0.975 ) ), col = 'black', lty = 3)


par( bg = 'black' )

dens <- density( perm_for_comparison )

dens_x <- dens$x

dens_y <- dens$y

plot( dens_x, dens_y, type = 'l', bty = 'l', xlab = paste( 'Proportion of time within the detection distance' ), xlim = range( c( real_for_comparison, perm_for_comparison ) ), ylab = "Probability density", col.axis = 'white', col.lab = 'white', col.main = 'white', col = 'white' )

abline( v = real_for_comparison, col = 'red' )

abline( v = quantile( perm_for_comparison, c( 0.025, 0.975 ) ), col = 'white', lty = 3)


axis( 1, col = 'white', at = seq( 0.005, 0.02, by = 0.005 ), labels = rep( '', length( seq( 0.005, 0.02, by = 0.005 ) ) ) )
axis( 2, col = 'white', at = seq( 0, 300, by = 50 ), labels = rep( '', length( seq( 0, 300, by = 50 ) ) ) )


p <- sum( real_for_comparison <= perm_for_comparison ) / max( perm_sub$perm_n )

final_p <- ifelse( p > 0.5, sum( real_for_comparison >= perm_for_comparison ) / max( perm_sub$perm_n ), p )

#mtext( paste( 'p =', final_p ), side = 4 )


### Mean dyadic distance while within the distance threshold

dev.off()

real_for_comparison <- real_agg$sum_dyad_dist_within_dist_thresh / real_agg$num_within_dist_thresh

if( is.na( real_for_comparison ) ){
  
  plot( 1, type = 'n', bty = 'n', main = paste( tag_a, tag_b ) )
}else{
  
  perm_for_comparison <- perm_agg$sum_dyad_dist_within_dist_thresh / perm_agg$num_within_dist_thresh
  
  perm_for_comparison <- perm_for_comparison[ !is.na( perm_for_comparison ) ]
  
  hist( perm_for_comparison, xlab = 'Mean dyadic distance within the distance threshold', xlim = range( c( real_for_comparison, perm_for_comparison ) ), main = "All dyads" )
  
  abline( v = real_for_comparison, col = 'red' )
  
  p <- sum( real_for_comparison <= perm_for_comparison ) / length( perm_for_comparison )
  
  final_p <- ifelse( p > 0.5, sum( real_for_comparison >= perm_for_comparison ) / length( perm_for_comparison ), p )
  
  mtext( paste( 'p =', final_p ), side = 4 )
  
}



dens <- density( perm_for_comparison )

dens_x <- dens$x

dens_y <- dens$y

plot( dens_x, dens_y, type = 'l', bty = 'l', xlab = "Mean dyadic distance during encounters (m)", xlim = range( c( real_for_comparison, perm_for_comparison ) ), ylab = "Probability density" )

abline( v = real_for_comparison, col = 'red' )

abline( v = quantile( perm_for_comparison, c( 0.025, 0.975 ) ), col = 'black', lty = 3)




par( bg = 'black' )

dens <- density( perm_for_comparison )

dens_x <- dens$x

dens_y <- dens$y

plot( dens_x, dens_y, type = 'l', bty = 'l', xlab = paste( 'Mean dyadic distance when within the detection distance (m)' ), xlim = range( c( real_for_comparison, perm_for_comparison ) ), ylab = "Probability density", col.axis = 'white', col.lab = 'white', col.main = 'white', col = 'white' )

abline( v = real_for_comparison, col = 'red' )

abline( v = quantile( perm_for_comparison, c( 0.025, 0.975 ) ), col = 'white', lty = 3)


axis( 1, col = 'white', at = seq( 275, 425, by = 25 ), labels = rep( '', length( seq( 275, 425, by = 25 ) ) ) )
axis( 2, col = 'white', at = seq( 0, 0.025, by = 0.005 ), labels = rep( '', length( seq( 0, 0.025, by = 0.005 ) ) ) )



### average difference in headings of the dyad members within the distance threshold

dev.off()

real_for_comparison <- real_agg$sum_heading_diff_in_dist_thresh / real_agg$length_heading_diff_in_dist_thresh

if( is.na( real_for_comparison ) ){
  
  plot( 1, type = 'n', bty = 'n', main = paste( tag_a, tag_b ) )
}else{
  
  perm_for_comparison <- perm_agg$sum_heading_diff_in_dist_thresh / perm_agg$length_heading_diff_in_dist_thresh
  
  perm_for_comparison <- perm_for_comparison[ !is.na( perm_for_comparison ) ]
  
  hist( perm_for_comparison, xlab = 'Mean difference in headings within the distance threshold', xlim = range( c( real_for_comparison, perm_for_comparison ) ), main = "All dyads" )
  
  abline( v = real_for_comparison, col = 'red' )
  
  p <- sum( real_for_comparison <= perm_for_comparison ) / length( perm_for_comparison )
  
  final_p <- ifelse( p > 0.5, sum( real_for_comparison >= perm_for_comparison ) / length( perm_for_comparison ), p )
  
  mtext( paste( 'p =', final_p ), side = 4 )
}


### average difference in step lengths of the dyad members within the distance threshold

dev.off()

real_for_comparison <- real_agg$sum_step_length_diff_in_dist_thresh / real_agg$length_step_length_diff_in_dist_thresh

if( is.na( real_for_comparison ) ){
  
  plot( 1, type = 'n', bty = 'n', main = paste( tag_a, tag_b ) )
}else{
  
  perm_for_comparison <- perm_agg$sum_step_length_diff_in_dist_thresh / perm_agg$length_step_length_diff_in_dist_thresh
  
  perm_for_comparison <- perm_for_comparison[ !is.na( perm_for_comparison ) ]
  
  hist( perm_for_comparison, xlab = 'Mean difference in step lengths within the distance threshold', xlim = range( c( real_for_comparison, perm_for_comparison ) ), main = "All dyads" )
  
  abline( v = real_for_comparison, col = 'red' )
  
  p <- sum( real_for_comparison <= perm_for_comparison ) / length( perm_for_comparison )
  
  final_p <- ifelse( p > 0.5, sum( real_for_comparison >= perm_for_comparison ) / length( perm_for_comparison ), p )
  
  mtext( paste( 'p =', final_p ), side = 4 )
}


### mean of dynamic time warping distance between the trajectories within the distance threshold 

dev.off()


real_for_comparison <- real_agg$sum_dtw_dist_in_dist_thresh / real_agg$length_dtw_dist_in_dist_thresh

if( is.na( real_for_comparison ) ){
  
  plot( 1, type = 'n', bty = 'n', main = paste( tag_a, tag_b ) )
}else{
  
  perm_for_comparison <- perm_agg$sum_dtw_dist_in_dist_thresh / perm_agg$length_dtw_dist_in_dist_thresh
  
  perm_for_comparison <- perm_for_comparison[ !is.na( perm_for_comparison ) ]
  
  hist( perm_for_comparison, xlab = "Mean dtw distance between the trajectories in the distance threshold", xlim = range( c( real_for_comparison, perm_for_comparison ) ), main = 'All dyads' )
  
  abline( v = real_for_comparison, col = 'red' )
  
  p <- sum( real_for_comparison <= perm_for_comparison ) / length( perm_for_comparison )
  
  final_p <- ifelse( p > 0.5, sum( real_for_comparison >= perm_for_comparison ) / length( perm_for_comparison ), p )
  
  mtext( paste( 'p =', final_p ), side = 4 )
}



dens <- density( perm_for_comparison )

dens_x <- dens$x

dens_y <- dens$y

plot( dens_x, dens_y, type = 'l', bty = 'l', xlab = "Mean DTW distance between trajectories during encounters", xlim = range( c( real_for_comparison, perm_for_comparison ) ), ylab = "Probability density" )

abline( v = real_for_comparison, col = 'red' )

abline( v = quantile( perm_for_comparison, c( 0.025, 0.975 ) ), col = 'black', lty = 3)




par( bg = 'black' )

dens <- density( perm_for_comparison )

dens_x <- dens$x

dens_y <- dens$y

plot( dens_x, dens_y, type = 'l', bty = 'l', xlab = paste( "Mean DTW distance between trajectories during encounters" ), xlim = range( c( real_for_comparison, perm_for_comparison ) ), ylab = "Probability density", col.axis = 'white', col.lab = 'white', col.main = 'white', col = 'white' )

abline( v = real_for_comparison, col = 'red' )

abline( v = quantile( perm_for_comparison, c( 0.025, 0.975 ) ), col = 'white', lty = 3)


axis( 1, col = 'white', at = seq( 160, 280, by = 20 ), labels = rep( '', length( seq( 160, 280, by = 20 ) ) ) )
axis( 2, col = 'white', at = seq( 0, 0.02, by = 0.005 ), labels = rep( '', length( seq( 0, 0.02, by = 0.005 ) ) ) )




