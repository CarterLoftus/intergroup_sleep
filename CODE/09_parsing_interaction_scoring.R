

library( stringr )
library( hms )

#### parsing interactions


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



# remove the rows of spec_df that do not have location data
spec_df <- spec_df[ !is.na( spec_df$x ), ]

# make the timestamp into a POSIX element
spec_df$local_timestamp <- as.POSIXct( spec_df$local_timestamp, tz = 'UTC' )

# turn the day column into a date rather than day from the start of the study period (makes it easier to merge information across dataframes and between day and night dataframes)
spec_df$day <- as.Date( spec_df$local_timestamp, tz = 'UTC' )

# remove LI TH's data when she is dying. The day chosen here is the day on which she starts sleeping separate from the rest of her group
spec_df <- spec_df[ !( spec_df$id == 'Pa_LI_TH' & spec_df$local_timestamp >= as.POSIXct( '2014-05-28 00:00:00') ), ]

# make a dataframe containing a row for each group on each day of data it has
temp_dyad_interact_df <- unique( spec_df[ , c( 'group', 'day' ) ] ) 

# save the names of the groups
group_names <- as.character( unique( temp_dyad_interact_df$group ) )

# repeat this dataframe by the number of groups that there are, because this is going to be a dyadic dataframe (one row for each dyad on each day on which the focal group has data)
dyad_interact_df <- temp_dyad_interact_df[ rep( seq_len( nrow( temp_dyad_interact_df ) ), each = length( group_names ) ), ]

## order the dataframe by the timestamp
dyad_interact_df <- dyad_interact_df[ order( dyad_interact_df$group) , ]
dyad_interact_df <- dyad_interact_df[ order( dyad_interact_df$day) , ]

dyad_interact_df$day <- as.Date( dyad_interact_df$day, tz = 'UTC' )

## change the column names of id and group to be id1 and group1. This will be our focal individual
names( dyad_interact_df )[ names( dyad_interact_df ) == 'group' ] <- 'group1'

## create a column for the names of the non-focal dyad member
dyad_interact_df$group2 <- group_names

# remove the rows where both dyad members are the same group
dyad_interact_df <- dyad_interact_df[ dyad_interact_df$group1 != dyad_interact_df$group2, ]

# make a column for the interaction that the dyad had on that day. We will instantiate this with "None" for no interaction. Then we will go through and add the interactions that did occur. After that, we will go through the days on which there wasn't an interaction between the dyad and insert NAs in columns where one of the dyad members didn't have enough data on that day to actually say that there was no interaction
dyad_interact_df$interaction <- 'None'


# now read in the interaction dataframes so that we can fill in the dyad_interact_df$interaction column with the interactions that did occur
full_interact <- read.csv( 'RESULTS/dyadDurat_interaction_scoring.csv' ) # this is a dataframe that is produced by taking "baboon_dyadDurat.csv" produced in 07_extracting_encounters.R and adding columns manually based on the manual scoring of the interactions (from the KMLs, also produced in 07_extracting_encounters.R). It is then saved manually in excel and used as an input here

# extract the final outcome of the interaction
full_interact$final_outcome <- ifelse( grepl( 'cosleep', full_interact$outcome ) & grepl( 'comove', full_interact$outcome ), "co-sleep with co-move", ifelse( grepl( 'cosleep', full_interact$outcome ), "co-sleep",  ifelse( grepl( 'comove', full_interact$outcome ) | grepl( 'milling', full_interact$context ), "co-move", ifelse( ( grepl( 'avoid', full_interact$context ) &  grepl( 'sleep_site', full_interact$context ) ) |  ( grepl( 'displace', full_interact$context ) &  grepl( 'sleep_site', full_interact$context ) ), "avoid or displace at sleep site", ifelse( grepl( "avoid", full_interact$context ) | grepl( 'displace', full_interact$context ) | grepl( 'follow', full_interact$context ) , 'avoid, displace, and pursue', ifelse( grepl( "crossing_path", full_interact$context ), 'crossing paths', NA ) ) ) ) ) )


full_interact[ is.na( full_interact$final_outcome ), ] # there should be no rows with NA for final_outcome


barplot( table( full_interact$final_outcome ), las = 1, ylab = 'frequency', xlab = 'outcome', main = 'Outcome of encounters' )

uniq_dyads <- unique( full_interact[ , c( 'group1', 'group2' ) ] )

for( i in 1:nrow( uniq_dyads ) ){
  
  dy_dat <- full_interact[ full_interact$group1 == uniq_dyads$group1[ i ] & full_interact$group2 == uniq_dyads$group2[ i ], ]
  
  barplot( table( dy_dat$final_outcome ), las = 1, ylab = 'frequency', main = paste(uniq_dyads$group1[ i ], uniq_dyads$group2[ i ] ) )
  
}

write.csv( full_interact, 'RESULTS/dyadDurat_interaction_scoring_with_final_outcome.csv', row.names = F )


full_interact$start_local_timestamp <- as.POSIXct( full_interact$start_local_timestamp, tz = 'UTC' )
full_interact$end_local_timestamp <- as.POSIXct( full_interact$end_local_timestamp, tz = 'UTC' )

# reorder the full_interact dataframe by timestamp of the start of the interaction
full_interact <- full_interact[ order( full_interact$start_local_timestamp ), ]

### percent of comove ###

sum( grepl( 'comove', full_interact$outcome ) )/ nrow( full_interact )

sum( full_interact$sleep_site_related )
sum( grepl( 'cosleep', full_interact$outcome ) )/ nrow( full_interact )
sum( grepl( 'cosleep', full_interact$outcome ) )
sum( grepl( 'cosleep', full_interact$context ) )/ nrow( full_interact )

full_interact[ full_interact$sleep_site_related == 0 &  grepl( 'sleep', full_interact$outcome ) , 'sleep_site_related' ] <- 1
sum( full_interact$sleep_site_related )

full_interact[ full_interact$sleep_site_related == 1 &  !grepl( 'cosleep', full_interact$outcome ) , ]


# add the interactions into the dyad_interact_df dataframe
for( i in 1:nrow( full_interact ) ){
  
  start_date <- as.Date( full_interact$start_local_timestamp[ i ], tz = 'UTC' )
  end_date <- as.Date( full_interact$end_local_timestamp[ i ], tz = 'UTC' )
  
  if( start_date == end_date ){
    
    if( full_interact$final_outcome[ i ] == 'crossing paths' ){
      
      dyad_interact_df[ ( ( dyad_interact_df$group1 == full_interact$group1[ i ] & dyad_interact_df$group2 == full_interact$group2[ i ] ) | ( dyad_interact_df$group1 == full_interact$group2[ i ] & dyad_interact_df$group2 == full_interact$group1[ i ] ) ) & dyad_interact_df$day >= start_date & dyad_interact_df$day <= end_date, 'interaction' ] <- 'crossing_paths'
      
    }
    
    if( grepl( 'avoid', full_interact$context[ i ] ) ){
      
      if( grepl( 'mutual', full_interact$context[ i ] ) ){
        
        dyad_interact_df[ ( ( dyad_interact_df$group1 == full_interact$group1[ i ] & dyad_interact_df$group2 == full_interact$group2[ i ] ) | ( dyad_interact_df$group1 == full_interact$group2[ i ] & dyad_interact_df$group2 == full_interact$group1[ i ] ) ) & dyad_interact_df$day == start_date, 'interaction' ] <- 'mutual_avoid'
        
      }else{
        
        if( grepl( 'sleep', full_interact$context[ i ] ) ){
          
          dyad_interact_df[ dyad_interact_df$group1 == full_interact[ i , str_split_fixed( full_interact$context[ i ], "_", 5 )[ , 1 ] ] & dyad_interact_df$group2 == full_interact[ i , str_split_fixed( full_interact$context[ i ], "_", 5 )[ , 3 ] ] & dyad_interact_df$day == start_date, 'interaction' ] <- 'sleep_site_avoider'
          
          dyad_interact_df[ dyad_interact_df$group1 == full_interact[ i , str_split_fixed( full_interact$context[ i ], "_", 5 )[ , 3 ] ] & dyad_interact_df$group2 == full_interact[ i , str_split_fixed( full_interact$context[ i ], "_", 5 )[ , 1 ] ] & dyad_interact_df$day == start_date, 'interaction' ] <- 'sleep_site_avoided'
          
        }else{
          
          dyad_interact_df[ dyad_interact_df$group1 == full_interact[ i , str_split_fixed( full_interact$context[ i ], "_avoid_", 2 )[ , 1 ] ] & dyad_interact_df$group2 == full_interact[ i , str_split_fixed( full_interact$context[ i ], "_avoid_", 2 )[ , 2 ] ] & dyad_interact_df$day == start_date, 'interaction' ] <- 'avoider'
          
          dyad_interact_df[ dyad_interact_df$group1 == full_interact[ i , str_split_fixed( full_interact$context[ i ], "_avoid_", 2 )[ , 2 ] ] & dyad_interact_df$group2 == full_interact[ i , str_split_fixed( full_interact$context[ i ], "_avoid_", 2 )[ , 1 ] ] & dyad_interact_df$day == start_date, 'interaction' ] <- 'avoided'
          
        }
        
      }
      
    }
    
    if( grepl( 'follow', full_interact$context[ i ] ) ){
      
      # if(  dyad_interact_df[ dyad_interact_df$group == full_interact[ i , str_split_fixed( full_interact$context[ i ], "_follow_", 2 )[ , 1 ] ] & dyad_interact_df$day == start_date, 'interaction' ] != 'None' ) stopppppppp
      
      dyad_interact_df[ dyad_interact_df$group1 == full_interact[ i , str_split_fixed( full_interact$context[ i ], "_follow_", 2 )[ , 1 ] ] & dyad_interact_df$group2 == full_interact[ i , str_split_fixed( full_interact$context[ i ], "_follow_", 2 )[ , 2 ] ] & dyad_interact_df$day == start_date, 'interaction' ] <- 'follower'
      
      dyad_interact_df[ dyad_interact_df$group1 == full_interact[ i , str_split_fixed( full_interact$context[ i ], "_follow_", 2 )[ , 2 ] ] & dyad_interact_df$group2 == full_interact[ i , str_split_fixed( full_interact$context[ i ], "_follow_", 2 )[ , 1 ] ] & dyad_interact_df$day == start_date, 'interaction' ] <- 'followed'
      
    }
    
    if( grepl( 'displace', full_interact$context[ i ] ) ){
      
      if( grepl( 'sleep', full_interact$context[ i ] ) ){
        
        dyad_interact_df[ dyad_interact_df$group1 == full_interact[ i , str_split_fixed( full_interact$context[ i ], "_", 5 )[ , 1 ] ] & dyad_interact_df$group2 == full_interact[ i , str_split_fixed( full_interact$context[ i ], "_", 5 )[ , 3 ] ] & dyad_interact_df$day == start_date, 'interaction' ] <- 'sleep_site_displacer'
        
        dyad_interact_df[ dyad_interact_df$group1 == full_interact[ i , str_split_fixed( full_interact$context[ i ], "_", 5 )[ , 3 ] ] & dyad_interact_df$group2 == full_interact[ i , str_split_fixed( full_interact$context[ i ], "_", 5 )[ , 1 ] ] & dyad_interact_df$day == start_date, 'interaction' ] <- 'sleep_site_displaced'
        
      }else{
        
        dyad_interact_df[ dyad_interact_df$group1 == full_interact[ i , str_split_fixed( full_interact$context[ i ], "_displace_", 2 )[ , 1 ] ] & dyad_interact_df$group2 == full_interact[ i , str_split_fixed( full_interact$context[ i ], "_displace_", 2 )[ , 2 ] ] & dyad_interact_df$day == start_date, 'interaction' ] <- 'displacer'
        
        dyad_interact_df[ dyad_interact_df$group1 == full_interact[ i , str_split_fixed( full_interact$context[ i ], "_displace_", 2 )[ , 2 ] ] & dyad_interact_df$group2 == full_interact[ i , str_split_fixed( full_interact$context[ i ], "_displace_", 2 )[ , 1 ] ] & dyad_interact_df$day == start_date, 'interaction' ] <- 'displaced'
        
      }
    }
    
    if( full_interact$final_outcome[ i ] == 'co-move' ){
      
      dyad_interact_df[ ( ( dyad_interact_df$group1 == full_interact$group1[ i ] & dyad_interact_df$group2 == full_interact$group2[ i ] ) | ( dyad_interact_df$group1 == full_interact$group2[ i ] & dyad_interact_df$group2 == full_interact$group1[ i ] ) ) & dyad_interact_df$day == start_date, 'interaction' ]  <- 'comove'
      
    }
    
    if( full_interact$final_outcome[ i ] == 'co-sleep' ){
      
      #if there is also comove in there, then this is wrong. Let's just check and make sure that there isn't
      if( grepl( 'comove', full_interact$outcome[ i ] ) ){
        
        stop( 'co-sleep with co-move falsely labeled as co-sleep' )
        
      }
      
      ## if the cosleep starts in the morning that it is actually a cosleep for the previous night and not the current night
      start_time <- str_split_fixed( full_interact$start_local_timestamp[ i ], " ", 2 )[ , 2 ]
      
      if( start_time < '12:00:00' ){ #if the coleep interaction starts in the morning, then the cosleep occurred on the previous night
        
        dyad_interact_df[ ( ( dyad_interact_df$group1 == full_interact$group1[ i ] & dyad_interact_df$group2 == full_interact$group2[ i ] ) | ( dyad_interact_df$group1 == full_interact$group2[ i ] & dyad_interact_df$group2 == full_interact$group1[ i ] ) ) & dyad_interact_df$day == ( start_date - 1 ), 'interaction' ] <- 'cosleep'
        
      }else{ # if the cosleep interaction happens after noon, however, then the cosleep occurs that night
        
        dyad_interact_df[ ( ( dyad_interact_df$group1 == full_interact$group1[ i ] & dyad_interact_df$group2 == full_interact$group2[ i ] ) | ( dyad_interact_df$group1 == full_interact$group2[ i ] & dyad_interact_df$group2 == full_interact$group1[ i ] ) ) & dyad_interact_df$day == ( start_date ), 'interaction' ] <- 'cosleep'
        
      }
      
    }
    
    if( full_interact$final_outcome[ i ] == 'co-sleep with co-move' ){
      
      order_of_events <- trimws( str_split( full_interact$outcome[ i ], ',', simplify = T ) )
      
      if( order_of_events[ 1 ] == 'cosleep' ){ ## I spot checked this so that I know it will work
        
        # then on the night before, mark them down for a cosleep followed by a comove
        dyad_interact_df[ ( ( dyad_interact_df$group1 == full_interact$group1[ i ] & dyad_interact_df$group2 == full_interact$group2[ i ] ) | ( dyad_interact_df$group1 == full_interact$group2[ i ] & dyad_interact_df$group2 == full_interact$group1[ i ] ) ) & dyad_interact_df$day == ( start_date - 1 ), 'interaction' ] <- 'cosleep_then_comove'
        
      }else{
        
        dyad_interact_df[ ( ( dyad_interact_df$group1 == full_interact$group1[ i ] & dyad_interact_df$group2 == full_interact$group2[ i ] ) | ( dyad_interact_df$group1 == full_interact$group2[ i ] & dyad_interact_df$group2 == full_interact$group1[ i ] ) ) & dyad_interact_df$day == start_date, 'interaction' ] <- 'comove_then_cosleep'
        
      }
      
    }
    
    
  }else{ ## if the start_date does not equal the end date (aka it is a multiday interaction)
    
    if( full_interact$final_outcome[ i ] == 'co-sleep' ){
      
      #if there is also comove in there, then this is wrong. Let's just check and make sure that there isn't
      if( grepl( 'comove', full_interact$outcome[ i ] ) ){
        
        stop( 'co-sleep with co-move falsely labeled as co-sleep' )
        
      }
      
      dyad_interact_df[ ( ( dyad_interact_df$group1 == full_interact$group1[ i ] & dyad_interact_df$group2 == full_interact$group2[ i ] ) | ( dyad_interact_df$group1 == full_interact$group2[ i ] & dyad_interact_df$group2 == full_interact$group1[ i ] ) ) & dyad_interact_df$day == ( start_date ), 'interaction' ] <- 'cosleep'
      
    }
    
    
    if( full_interact$final_outcome[ i ] == 'co-move' ){ ## this one is just hard-coded
      
      dyad_interact_df[ dyad_interact_df$group1 == full_interact[ i , 'group1' ] & dyad_interact_df$group2 == full_interact[ i , 'group2' ] & dyad_interact_df$day == ( start_date ), 'interaction' ] <- 'sleep_site_avoider'
      
      dyad_interact_df[ dyad_interact_df$group1 == full_interact[ i , 'group2' ] & dyad_interact_df$group2 == full_interact[ i , 'group1' ] & dyad_interact_df$day == ( start_date ), 'interaction' ] <- 'sleep_site_avoided'
      
      
      dyad_interact_df[ ( ( dyad_interact_df$group1 == full_interact$group1[ i ] & dyad_interact_df$group2 == full_interact$group2[ i ] ) | ( dyad_interact_df$group1 == full_interact$group2[ i ] & dyad_interact_df$group2 == full_interact$group1[ i ] ) ) & dyad_interact_df$day == ( start_date + 1 ), 'interaction' ] <- 'comove'
      
    }
    
    if( full_interact$final_outcome[ i ] == 'co-sleep with co-move' ){
      
      order_of_events <- trimws( str_split( full_interact$outcome[ i ], ',', simplify = T ) )
      
      sleeps <- which( grepl( 'sleep', order_of_events ) )
      
      if( length( sleeps ) == 1 ){
        
        comoves <- which( grepl( 'comove', order_of_events ) )
        
        comove_before_sleep <- ifelse( sum( comoves < sleeps ) > 0, 1, 0 )
        
        comove_after_sleep <- ifelse( sum( comoves > sleeps ) > 0, 1, 0 )
        
        if( comove_before_sleep == 1 & comove_after_sleep == 0 ){
          
          dyad_interact_df[ ( ( dyad_interact_df$group1 == full_interact$group1[ i ] & dyad_interact_df$group2 == full_interact$group2[ i ] ) | ( dyad_interact_df$group1 == full_interact$group2[ i ] & dyad_interact_df$group2 == full_interact$group1[ i ] ) ) & dyad_interact_df$day == ( start_date ), 'interaction' ] <- 'comove_then_cosleep'
          
        }
        
        if( comove_before_sleep == 0 & comove_after_sleep == 1 ){
          
          dyad_interact_df[ ( ( dyad_interact_df$group1 == full_interact$group1[ i ] & dyad_interact_df$group2 == full_interact$group2[ i ] ) | ( dyad_interact_df$group1 == full_interact$group2[ i ] & dyad_interact_df$group2 == full_interact$group1[ i ] ) ) & dyad_interact_df$day == ( start_date ), 'interaction' ] <- 'cosleep_then_comove'
          
          dyad_interact_df[ ( ( dyad_interact_df$group1 == full_interact$group1[ i ] & dyad_interact_df$group2 == full_interact$group2[ i ] ) | ( dyad_interact_df$group1 == full_interact$group2[ i ] & dyad_interact_df$group2 == full_interact$group1[ i ] ) ) & dyad_interact_df$day == ( start_date + 1 ), 'interaction' ] <- 'comove'
        }
        
        if( comove_before_sleep == 1 & comove_after_sleep == 1 ){
          
          dyad_interact_df[ ( ( dyad_interact_df$group1 == full_interact$group1[ i ] & dyad_interact_df$group2 == full_interact$group2[ i ] ) | ( dyad_interact_df$group1 == full_interact$group2[ i ] & dyad_interact_df$group2 == full_interact$group1[ i ] ) ) & dyad_interact_df$day == ( start_date ), 'interaction' ] <- 'comove_then_cosleep_then_comove'
          
          dyad_interact_df[ ( ( dyad_interact_df$group1 == full_interact$group1[ i ] & dyad_interact_df$group2 == full_interact$group2[ i ] ) | ( dyad_interact_df$group1 == full_interact$group2[ i ] & dyad_interact_df$group2 == full_interact$group1[ i ] ) ) & dyad_interact_df$day == ( start_date + 1 ), 'interaction' ] <- 'comove'
          
        }
        
      }else{ ### this is just hard-coded because there are only 2 cases of multi-day interactions and they are the same order of interactions
        
        dyad_interact_df[ ( ( dyad_interact_df$group1 == full_interact$group1[ i ] & dyad_interact_df$group2 == full_interact$group2[ i ] ) | ( dyad_interact_df$group1 == full_interact$group2[ i ] & dyad_interact_df$group2 == full_interact$group1[ i ] ) ) & dyad_interact_df$day == ( start_date ), 'interaction' ] <- 'cosleep_then_comove'
        
        dyad_interact_df[ ( ( dyad_interact_df$group1 == full_interact$group1[ i ] & dyad_interact_df$group2 == full_interact$group2[ i ] ) | ( dyad_interact_df$group1 == full_interact$group2[ i ] & dyad_interact_df$group2 == full_interact$group1[ i ] ) ) & dyad_interact_df$day == ( start_date + 1 ), 'interaction' ] <- 'comove_then_cosleep_then_comove'
        
        dyad_interact_df[ ( ( dyad_interact_df$group1 == full_interact$group1[ i ] & dyad_interact_df$group2 == full_interact$group2[ i ] ) | ( dyad_interact_df$group1 == full_interact$group2[ i ] & dyad_interact_df$group2 == full_interact$group1[ i ] ) ) & dyad_interact_df$day == ( start_date + 2 ), 'interaction' ] <- 'comove'
        
      }
    }
  }
}


# View( dyad_interact_df[ dyad_interact_df$interaction != 'None', ] )
# some spot checking suggests that this dataframe has been created correctly

# Now go through the days on which a dyad doesn't have an interaction, and determine whether we actually have enough data from both individuals on that day to say whether they interacted or not

## This threshold is the proportion of successful fixes that must have occurred from 06:00 AM to 06:00 AM the next day to say whether we know whether the group did not have an interaction on days that there are no interaction recorded (interactions that happen before 06:00 would only be cosleeps that would count towards the previous day)
data_sufficiency_threshold <- 0.7


fixes_in_12_hours <- length( seq( as.numeric( as_hms( "00:00:00" ) ), as.numeric( as_hms( "12:00:00" ) ), by =  15*60 ) )

fixes_in_24_hours <- ( fixes_in_12_hours - 1 ) * 2

no_interact_inds <- which( dyad_interact_df$interaction == 'None' )

for( i in no_interact_inds ){
  
  
  ## can also do the unique by group and local timestamp so that if the two group members together have enough coverage, then it would be sufficient
  id1_data <- spec_df[ spec_df$group == dyad_interact_df$group1[ i ] & spec_df$local_timestamp >= as.POSIXct( paste( as.Date( dyad_interact_df$day[ i ], origin = '1970-01-01' ), '06:00:00' ), tz = 'UTC' ) & spec_df$local_timestamp <= as.POSIXct( paste( as.Date( dyad_interact_df$day[ i ] + 1 , origin = '1970-01-01' ), '06:00:00' ), tz = 'UTC' ), ]
  
  id1_data <- id1_data[ !duplicated( id1_data$local_timestamp ), ]
  
  id2_data <- spec_df[ spec_df$group == dyad_interact_df$group2[ i ] & spec_df$local_timestamp >= as.POSIXct( paste( as.Date( dyad_interact_df$day[ i ], origin = '1970-01-01' ), '06:00:00' ), tz = 'UTC' ) & spec_df$local_timestamp <= as.POSIXct( paste( as.Date( dyad_interact_df$day[ i ] + 1 , origin = '1970-01-01' ), '06:00:00' ), tz = 'UTC' ), ]
  
  id2_data <- id2_data[ !duplicated( id2_data$local_timestamp ), ]
  
  
  if( ( nrow( id1_data ) < data_sufficiency_threshold*fixes_in_24_hours ) | ( nrow( id2_data ) < data_sufficiency_threshold*fixes_in_24_hours  ) ){
    
    dyad_interact_df$interaction [ i ] <- NA
    
  }
  
}


write.csv( dyad_interact_df, "DATA/dyad_interact_df.csv", row.names =  F )
