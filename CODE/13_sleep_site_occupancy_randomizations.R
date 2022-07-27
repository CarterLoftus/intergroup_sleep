


#### sleep site randomizations ######

library( infotheo )
library( hms )
library( brms )
library( boot )
library( data.table )

## function for normalizing a vector
normalize_func <- function( x ) return( (x - mean( x, na.rm = T ) )/ sd( x, na.rm = T ) )


which_spec <- 'baboon'

if( which_spec == "baboon" ){
  spec_df <- read.csv( "DATA/bab_complete.csv" )
}else{
  if( which_spec == "vervet" ){
    spec_df <- read.csv( "DATA/verv_complete.csv" )
  }else{
    if( which_spec == "leopard" ){
      spec_df <- read.csv( "DATA/leo_complete.csv" )
    }
  }
}


## just to confirm that each group has a maximum of one sleep site assigned each night 
same_sleep_check <- aggregate( spec_df$sleep_clus, by = list( spec_df$group, spec_df$day ), FUN = function( x ) sum( !is.na( unique( x ) ) ) )

names( same_sleep_check ) <- c( 'group', 'day', 'num_sleep_sites')

same_sleep_check[ same_sleep_check$num_sleep_sites != 1, ]


# make a dataframe with one row per group per night, stating when they left the sleep site that morning and then they arrived at their sleep site in the evening (averaged across the individuals in the group when there is more than one)
cosleep_dat <- aggregate( spec_df[ , c( 'arrive_sleep_site', 'leave_sleep_site' ) ], by = list( spec_df$group, spec_df$day, spec_df$sleep_clus ), FUN = function( x ) as_hms( mean( as.numeric( as_hms( as.character( x ) ) ), na.rm = T ) ) )

# rename the columns of the dataframe
names( cosleep_dat )[ 1:3 ] <- c( 'group', 'day', 'sleep_clus' )

# reorder the dataframe
cosleep_dat <- cosleep_dat[ order( cosleep_dat$day ), ]
cosleep_dat <- cosleep_dat[ order( cosleep_dat$group ), ]

# make a column that will declare whether the group coslept with another group that night. Instantiate the column with all 0s
cosleep_dat$cosleep <- 0

# fill in the cosleep columns with 1's on nights when the group slept at the same site as another group
cosleep_dat[ duplicated( cosleep_dat[ , c( 'day', 'sleep_clus' ) ] ) | duplicated( cosleep_dat[ , c( 'day', 'sleep_clus' ) ], fromLast = T ), 'cosleep' ] <- 1


#### how many nights do they cosleep on?
sum( duplicated( cosleep_dat[ , c( 'day', 'sleep_clus' ) ] ) )



#### rand parameters ####

# make a dataframe with one row per group per night, stating when they left the sleep site that morning and then they arrived at their sleep site in the evening (averaged across the individuals in the group when there is more than one)
cosleep_dat <- aggregate( spec_df[ , c( 'arrive_sleep_site', 'leave_sleep_site' ) ], by = list( spec_df$group, spec_df$day, spec_df$sleep_clus ), FUN = function( x ) as_hms( mean( as.numeric( as_hms( as.character( x ) ) ), na.rm = T ) ) )

# rename the columns of the dataframe
names( cosleep_dat )[ 1:3 ] <- c( 'group', 'day', 'sleep_clus' )

cosleep_dat$id <- cosleep_dat$group
cosleep_dat$func_ind <- cosleep_dat$group


tag_names <- as.character( unique( cosleep_dat$id ) )

n <- 1000

rand_day_thresh <- 30

min_to_rand <- 1

# create empty vectors. These vectors will be filled with entries that eventually be put into the final dataframe that declares which functuals need to be compared
vec_a <- c()
vec_b <- c()


for( a in 1:( length( tag_names ) - 1 ) ){
  
  for( b in ( a + 1 ): length( tag_names ) ){
    
    # create vectors that represent the unique combinations that can be made of functuals in this category
    vec_a <- c( vec_a, tag_names[ a ] )
    
    vec_b <- c( vec_b, tag_names[ b ] )
    
  }
}

vec_cat_a <- rep( 'baboon', length( vec_a ) )

vec_cat_b <- vec_cat_a

# set up the final dataframe that will be filled out to give the empirical coefficient of associations between each dyad of functuals
real_final_cosleep <- data.frame( func_a = vec_a, func_b = vec_b, num = rep( NA, times = length( vec_a ) ), denom = rep( NA, times = length( vec_a ) ), MI = NA, cat_a = vec_cat_a, cat_b = vec_cat_b , stringsAsFactors = F)

# set up the final dataframe that will be filled out to give the coefficient of associations between each dyad of functuals produced by the randomizations
rand_final_cosleep <- data.frame( func_a = rep( vec_a, times = n ), func_b = rep( vec_b, times = n ), rand_n = rep( 1:n , each = length( vec_a ) ), num = rep( NA, times = length( rep( vec_a, times = n ) ) ), denom = rep( NA, times = length( rep( vec_a, times = n ) ) ), MI = NA, cat_a = rep( vec_cat_a, times = n ), cat_b = rep( vec_cat_b, times = n ), stringsAsFactors = F)

# create a dataframe that will contain the metadata for the randomizations
meta_cosleep_dat <- data.frame( func_a = vec_a, func_b = vec_b, start_rand_at = rep( NA, length( vec_a ) ), end_rand_at =  rep( NA, length( vec_a ) ), cat_a = vec_cat_a, cat_b = vec_cat_b, stringsAsFactors = F )

cosleep_dat$day <- as.numeric( cosleep_dat$day )

# saves the first day of the study
start_date <- min( cosleep_dat$day )

cosleep_sites_real <- data.frame( func_a = character(), func_b = character(), shared_site = integer() )
cosleep_sites_rand <- data.frame( func_a = character(), func_b = character(), rand_num = integer(), shared_site = integer() )

set.seed( 111 )

for( row in 1:nrow( real_final_cosleep) ){
  
  print( row / nrow(real_final_cosleep) )
  
  # subset the full gps dataframe to just including the functual dyad's data during the correct period. This appropriate combination is determined by the set up of the dataframes above
  pair_cosleep_dat <- cosleep_dat[ cosleep_dat$id %in% c( real_final_cosleep[ row, c('func_a', 'func_b') ]), ]
  
  # just making sure 'func_ind' is a character and not a factor. It messes things up if it is a factor
  pair_cosleep_dat$func_ind <- as.character( pair_cosleep_dat$func_ind ) 
  
  # trims the dataframe so it starts on the first day that both members of the dyad have data. We don't want to analyze anything before this
  pair_cosleep_dat <- pair_cosleep_dat[ pair_cosleep_dat$day >=  max( aggregate( pair_cosleep_dat$day , by=list( pair_cosleep_dat$func_ind ), FUN = min )[ , 2 ] ) , ] 
  
  # makes a column with both of their day columns starting at 1 on the first day when they both have the data. Important for chunking up the data in the next line
  pair_cosleep_dat$temp_day <- as.numeric( pair_cosleep_dat$day - min( pair_cosleep_dat$day ) + 1, units = 'days' )
  
  # assign each row of data to a subset so that days of the data are only randomized within an range determined by rand_day_thresh
  pair_cosleep_dat$chunk_num <- ceiling( pair_cosleep_dat$temp_day / rand_day_thresh ) 
  
  # split up the data into the subsets created in the line above. Now we have a list of dataframes, which each dataframe corresponding to one time chunk within randomizaation is allowable
  chunked_cosleep_dat <- split( pair_cosleep_dat, f = pair_cosleep_dat$chunk_num )
  
  # create an empty dataframe that will represent a dyadic distance at every simultaneous fix during the current study period
  total_real <- data.frame( day = integer(), sleep_site_a = integer(), sleep_site_b = integer() )
  
  # perform n permutations of the dyad's dataset
  for( i in 1:n ){
    
    # create an empty dataframe that will represent a dyadic distance at every derived simultaneous fix of the study period that results after the randomization
    total_rand <- data.frame( day = integer(), sleep_site_a = integer(), sleep_site_b = integer() )
    
    # loop through each chunk and permute the data within the chunks
    for( w in 1:length( chunked_cosleep_dat ) ){ # For each chunk of data
      
      # save the chunk of data to the dataframe "chunk"
      chunk <- chunked_cosleep_dat[[w]]
      
      # save the number of days of data that individuals a and b have in this chunk
      num_unique_days_a <- length( unique( chunk[ chunk$func_ind == real_final_cosleep$func_a[ row ], 'day' ] ) ) 
      num_unique_days_b <- length( unique( chunk[ chunk$func_ind == real_final_cosleep$func_b[ row ], 'day' ] ) ) 
      
      # if one of the functuals has no data for this chunk, skip the rest of the body of the loop and move to the next chunk
      if( num_unique_days_a == 0 || num_unique_days_b == 0 ){
        next
      }
      
      # if either individual a or b don't have the necessary number of days worth of data as determined by the min_to_rand parameter, skip the rest of the body of the loop
      if( num_unique_days_a < (min_to_rand * rand_day_thresh)  || num_unique_days_b < (min_to_rand * rand_day_thresh) ){ 
        next
      }
      
      # The first time through, we will make the empirical dataframe of dyadic distances. We only need to do this once
      if( i == 1 ){
        
        # use the dyad_dist function to calculate the dyadic distance for this pair of functuals
        real_sub <- as.data.frame( dcast( as.data.table( chunk ), day ~ id, value.var = 'sleep_clus', drop = F ) ) 
        
        names( real_sub ) <- c( 'day', 'sleep_site_a', 'sleep_site_b' )
        
        # add this to the running dataframe of dyadic distances over the whole study for this dyad
        total_real <- rbind( total_real, real_sub )
        
        # save the latest time that will be successfully randomized. For the first run through for each functual dyad, this will get updated every time the conditions above are surpassed such that a chunk of data is successfully randomized. Eventually this will serve to mark the end of successful randomizations
        end_of_rand <- max(chunked_cosleep_dat[[w]]$day)
        
      }
      
      # determine how many days to shift ID b's data by for the randomization
      ID_b_shifts <- sample( 1:( rand_day_thresh - 1 ), 1, replace = TRUE )
      
      # shift ID b's data by the amount determined above
      new_days_b <- chunk[ chunk$func_ind == real_final_cosleep$func_b[ row ], 'day'] + ID_b_shifts
      
      # complicated line of code, but all it does is wrap the end of b's data back around to match the beginning of a's data. So if we shifted b's data by 3, a's data will still be 1, 2, 3, 4, 5, 6, 7; and b's data will be 5, 6, 7, 1, 2, 3, 4.
      new_days_b[ new_days_b > max( chunk[ chunk$func_ind == real_final_cosleep$func_b[ row ], 'day' ] ) ] <- new_days_b[ new_days_b > max( chunk [ chunk$func_ind == real_final_cosleep$func_b[ row ], 'day' ] ) ] - max( chunk[ chunk$func_ind == real_final_cosleep$func_b[ row ], 'day' ] ) + min( chunk[ chunk$func_ind == real_final_cosleep$func_b[ row ], 'day'] ) - 1
      
      # replace b's day data with these 'fake' shifted days
      chunk[ chunk$func_ind == real_final_cosleep$func_b[ row ], 'day'] <- new_days_b
      
      rand_sub <- as.data.frame( dcast( as.data.table( chunk ), day ~ id, value.var = 'sleep_clus', drop = F ) ) 
      
      names( rand_sub ) <- c( 'day', 'sleep_site_a', 'sleep_site_b' )
      
      # add this to the running dataframe of derived dyadic distances for the randomized data over the whole study
      total_rand <- rbind( total_rand, rand_sub )
      
    }
    
    # if we are on the first run through the loop, save the empirical CA of the dyad
    if( i == 1 ){
      
      # adds the empirical coefficient of association for this dyad
      real_final_cosleep[ row , 'num' ] <- sum( total_real$sleep_site_a == total_real$sleep_site_b, na.rm = T ) 
      
      real_final_cosleep[ row , 'denom' ] <- sum( !is.na( total_real$sleep_site_a == total_real$sleep_site_b ) )
      
      real_final_cosleep[ row , 'MI' ] <- mutinformation( as.character( total_real$sleep_site_a ) , as.character( total_real$sleep_site_b ) )
      
      meta_cosleep_dat[ row, c( 'start_rand_at', 'end_rand_at' ) ] <- c( min( chunked_cosleep_dat [[ 1 ]]$day ), end_of_rand )
      
      sites_real <- total_real$sleep_site_a[ total_real$sleep_site_a == total_real$sleep_site_b & !is.na( total_real$sleep_site_a == total_real$sleep_site_b ) ] 
      
      if( length( sites_real ) != 0 ){
        
        sites_real_sub <- data.frame( func_a = real_final_cosleep$func_a[ row ], func_b = real_final_cosleep$func_b[ row ], shared_site = sites_real )
        
        cosleep_sites_real <- rbind( cosleep_sites_real, sites_real_sub )
        
      }
      
    }
    
    # save the randomized CA of the functual dyad for this randomization
    rand_final_cosleep[ ( row + nrow( real_final_cosleep ) * ( i - 1 ) ), 'num' ] <- sum( total_rand$sleep_site_a == total_rand$sleep_site_b, na.rm = T ) 
    
    rand_final_cosleep[ ( row + nrow( real_final_cosleep ) * ( i - 1 ) ), 'denom' ] <- sum( !is.na( total_rand$sleep_site_a == total_rand$sleep_site_b ) )
    
    rand_final_cosleep[ ( row + nrow( real_final_cosleep ) * ( i - 1 ) ) , 'MI' ] <- mutinformation( as.character( total_rand$sleep_site_a ), as.character( total_rand$sleep_site_b ) )
    
    sites_rand <- total_rand$sleep_site_a[ total_rand$sleep_site_a == total_rand$sleep_site_b & !is.na( total_rand$sleep_site_a == total_rand$sleep_site_b ) ]
    
    if( length( sites_rand ) != 0 ){
      
      sites_rand_sub <- data.frame( func_a = real_final_cosleep$func_a[ row ], func_b = real_final_cosleep$func_b[ row ], rand_num = i, shared_site = sites_rand )
      
      cosleep_sites_rand <- rbind( cosleep_sites_rand, sites_rand_sub )
      
    }
  }
}




real_agg <- aggregate( real_final_cosleep[ , c('num', 'denom') ], by = list( real_final_cosleep$cat_a, real_final_cosleep$cat_b ), FUN = sum, na.rm = T)

names( real_agg ) <- c( 'cat_a', 'cat_b', 'num', 'denom' )

real_agg$prop <- real_agg$num / real_agg$denom

rand_agg <- aggregate( rand_final_cosleep[ , c('num', 'denom') ], by = list( rand_final_cosleep$cat_a, rand_final_cosleep$cat_b, rand_final_cosleep$rand_n ), FUN = sum, na.rm = T)

names( rand_agg ) <- c( 'cat_a', 'cat_b', 'rand_n', 'num', 'denom' )

rand_agg$prop <- rand_agg$num / rand_agg$denom

p <- sum( real_agg$prop <= rand_agg$prop ) / max( rand_agg$rand_n )

final_p <- ifelse( p > 0.5, sum( real_agg$prop >= rand_agg$prop ) / max( rand_agg$rand_n ), p )



par( bg = 'white' )

dens_rand <- density( rand_agg$prop )

plot( dens_rand$x, dens_rand$y, col = 'black', col.axis = 'black', col.lab = 'black', col.main = 'black', xlab = 'Probability of two groups sleeping at same site', type = 'l', main = '', ylab = 'Probability density', bty = 'l' )

axis(1, col = 'black', tick = T, labels = F  )

axis(2, col = 'black', tick = T, labels = F  )

abline( v = quantile( rand_agg$prop, c( 0.025, 0.975 ) ), col = 'black', lty = 3)

abline( v = real_agg$prop, col = 'red', lty = 1 )

print( paste( 'p-value = ', round( final_p, 4) ), side = 4 )




par( bg = 'black' )

dens_rand <- density( rand_agg$prop )

plot( dens_rand$x, dens_rand$y, col = 'white', col.axis = 'white', col.lab = 'white', col.main = 'white', xlab = 'Probability of two groups sleeping at same site', type = 'l', main = '', ylab = 'Probability density', bty = 'l' )

axis(1, col = 'white', tick = T, labels = F  )

axis(2, col = 'white', tick = T, labels = F  )

abline( v = quantile( rand_agg$prop, c( 0.025, 0.975 ) ), col = 'white', lty = 3)

abline( v = real_agg$prop, col = 'red', lty = 1 )

print( paste( 'p-value = ', round( final_p, 4) ), side = 4 )



