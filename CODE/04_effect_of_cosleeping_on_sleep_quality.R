




library(sjPlot)
library(insight)
library(httr)
library(brms)
library( zoo )
library( hms )
library( data.table )
library( stringr )
library( lubridate )
library( lmerTest )
library( plotrix )
library( suncalc )
library( LaplacesDemon )
library( dplyr )
library( purrr )
library( HDInterval )
library(multcomp)
library( nlme )
library(tidyr) 
library(lmerTest)
library( sp )
library( rgdal )
library( stats )
library(rgeos)
library( entropy )
library( reshape2 )
library( plyr )
library(rstan)
library( brms )
library(fitdistrplus)
library( gamm4 )
library(glmmTMB)
library( mgcv )
library( rstudioapi )

################# Functions #########################

## function for normalizing a vector
normalize_func <- function( x ) return( (x - mean( x, na.rm = T ) )/ sd( x, na.rm = T ) )


## function for setting transparency of a color while plotting
transp <- function(col, alpha=.5){
  res <- apply(col2rgb(col),2, function(c) rgb(c[1]/255, c[2]/255, c[3]/255, alpha))
  return(res)
}

## function for plotting times from noon to noon. It will make correct 12:00 - 24:00 to be 00:00 - 12:00 and 00:00 - 12:00 to be 12:00 to 24:00

ts_func <- function( time_vec ){
  
  num_time <- as.numeric( as_hms( time_vec ) )
  
  corr_time <- ifelse( num_time < 12*60*60, num_time + 12*60*60, num_time - 12*60*60 )
  
  return( corr_time )
  
}

#### sleep period function

## missing_mins: this is the maximum total number of minutes of data that can be missing from a day and still have that day included in the analysis

## time_gap: this is the maximum allowable time gap between two accelerometer bursts (in seconds) that can exist in a day without removing this day from the data

## move_window: this is the size of the moving window (in minutes) used in calculating the rolling median 

## percentile: this is the percentile threshold of the variable (within the noon-to-noon time period) used to classify activity vs. inactivity 

## multiplier: this is the multiplier of the threshold value determined by the percentile above. Values below the threshold value multiplied by the multiplier are considered inactive and those equal to or above are considered active

## block_size: duration in minutes of the blocks of continuous inactivity that will be considered sleep

## gap_size: maximum duration between sleep blocks that will be merged

## title: if title == T, the plot will include the tag number and the night number at the top of the plot

## x_axis: if x_axis == F, the plot will be plotted without an x axis
## waso_percentile: this is the percentile threshold of the variable (within the noon-to-noon time period) used to classify activity vs. inactivity within the sleep period

## waso_multiplier: this is the multiplier of the threshold value determined by the percentile above. Values WITHIN the sleep period below the threshold value multiplied by the multiplier are considered inactive and those equal to or above are considered active. This threshold value may be different than the one above even when waso_percentile = percentile and waso_multiplier = multiplier, because the value of the given percentile of this variable depends on the smoothing window (see waso_window below; aka waso window might not be equal to mov_window)

## waso_window: this is the size of the moving window (in minutes) used in calculating the rolling median that will be used to find periods of wake after sleeping. A waso_window of 1 is the same as using the raw data without a rolling median

## waso_block: this is the number of consecutive minutes of inactivity needed to classify a period as sleep within the sleep period. A waso_block of 1 means that anytime the value is below the threshold, the baboon in considered sleeping and anytime the value is above the threshold the baboon is considered awake

## las: las sets las in the plotting window


sleep_per_func <- function( tag, night_num, missing_mins = 45, time_gap = 20*60, mov_window = 9, percentile = 0.10, multiplier = 1.125, block_size = 30, gap_size = 45, title = F, x_axis = T, plot_waso = F, waso_window = 1, waso_block = 1, waso_percentile = 0.10, waso_multiplier = 1.125, las = 1, ... ){
  
  ## save a variable denoting the total number of minutes in the day
  mins_in_day <- 60*24
  
  ## subset the data to the given tag on the given night
  night_dat <- d1[ d1$tag == tag & d1$night_num == night_num, ]
  
  ## sort the timestamps (they are probably already sorted)
  sorted_times <- sort( night_dat$local_time )
  
  ## find the time difference in seconds between each burst
  time_diffs <- as.numeric( diff( as_hms( sorted_times ) ), units = 'secs' )
  
  ## if there is more than a single burst...
  if(length(time_diffs) != 0){
    
    ## if the number of bursts exceed the minimum required number of bursts in a night (determined by missing mins) and if the gaps in the data are within the allowable gap size (determined by time_gap)...
    if( nrow( night_dat ) > ( mins_in_day - missing_mins ) & max( time_diffs ) < time_gap ){
      
      ## take the rolling median of the log VeDBA and save it as a column
      night_dat$roll_log_vedba <- rollmedian( night_dat$log_vedba, mov_window, fill = NA, align = 'center' )
      
      ## determine the threshold activity vs. inactivity threshold based on the percentile, multiplier, and the rolling median just produced
      thresh <- quantile( night_dat$roll_log_vedba, percentile, na.rm = T) * multiplier
      
      ## put the rows of the dataframe in order from noon to noon (they should already be in this order, so this should be redundant)
      night_dat <- night_dat[ order( ts_func( night_dat$local_time ) ), ]
      
      ## turn the times into numerical elements for plotting
      ts_time <- ts_func( night_dat$local_time )
      
      if( title == F ){
        ## plot the log VeDBA
        #plot( ts_time, night_dat$log_vedba, type = 'l', xlab = 'Time', ylab = '', xaxt = 'n', las = las )
        
        plot( ts_time, night_dat$log_vedba, type = 'l', xlab = 'Time', ylab = '', xaxt = 'n', las = las, ylim = c( 1.9, 7 ), ... )
        
      }else{
        ## plot the log VeDBA
        plot( ts_time, night_dat$log_vedba, type = 'l', xlab = 'Time', ylab = '', main = paste( tag, night_num ), xaxt = 'n', las = las, ... )
        
        
      }
      
      if( x_axis == T ){
        
        axis( 1, at = seq( 0, 60*24*60, 60*60), labels = c( as_hms( seq( 12*60*60, 60*23*60, 60*60) ), as_hms( seq( 0, 60*12*60, 60*60) ) ) ) 
        
      }
      
      title( ylab = 'log VeDBA', line = 3.9 )
      ## plot the rolling median of the log VeDBA
      lines( ts_time, rollmedian( night_dat$log_vedba, mov_window, fill = NA, align = 'center' ), col = 'red')
      
      ## plot the threshold of the log VeDBA
      abline( h = thresh, col = 'blue' )
      
      ### find blocks of continuous inactivity
      
      ## find the run length encoding of periods above and below the threshold
      temp <- rle(as.numeric( night_dat$roll_log_vedba < thresh ) ) 
      
      ## mark the rows that are part of runs (i.e. part of chunks that are greater than the block_size of either continuous activity or continuous inactivity )
      night_dat$runs <- as.numeric( rep( temp$lengths > block_size, times = temp$lengths ) )
      
      ## mark the rows corresponding to sleep bouts. These sleep bouts are runs of inactivity
      night_dat$sleep_bouts <- as.numeric( night_dat$roll_log_vedba < thresh & night_dat$runs == 1 )
      
      ## find when sleep bouts start and end
      diffs <- diff( c(0, night_dat$sleep_bouts ) )
      starts <- which( diffs == 1 ) [ -1 ]
      ends <- which( diffs == -1 )
      
      ## if there are any sleep bouts...
      if( length( which( diffs == 1 ) ) != 0){
        
        ## find the duration of the gaps between each sleep bout (the end of one sleep bout and the start of the next)
        gaps <- as.numeric( night_dat$local_timestamp [ starts ] - night_dat$local_timestamp [ ends[ 1: length( starts ) ] ], units = 'mins' )
        
        ## sleep bouts separated by gaps that are shorter than that specified by gap_size will be merged. Note which of these gaps are shorter than the gap_size
        inds_to_remove <- which( gaps < gap_size ) 
        
        ## if there are NO gaps between sleep bouts that are to be removed...
        if( length( inds_to_remove ) == 0 ){
          
          ## set sleep onset index to be the start of sleep bouts
          onset <- which( diffs == 1 ) 
          
          ## set waking index to be the end of sleep bouts
          wake <- ends
          
        }else{ ## if there ARE gaps between sleep bouts that are to be removed...
          
          ## set sleep onset index to be the start of sleep bouts that do not correspond to the gaps to be removed (because these will be within sleep periods, not a start of a new bout)
          onset <- which( diffs == 1 ) [ - (inds_to_remove + 1) ]
          
          ## set waking index to be the end of sleep bouts that do not correspond to the gaps to be removed
          wake <- ends [ - inds_to_remove ]
          
        }
        
        ## determine which sleep period is the longest
        per_ind <- which.max( as.numeric( night_dat$local_timestamp[ wake ] - night_dat$local_timestamp[ onset ], units = 'secs' ) )
        
        ## plot the sleep onset time and waking time on the log VeDBA plot
        abline( v = c( ts_time[ onset[ per_ind ] ], ts_time[ wake[ per_ind ] ] ), col = 'orange', lty = 3, lwd = 4 )
        
        ## if you also want to plot WASO
        if( plot_waso == T ){
          
          ## calculate the rolling median of the log VeDBA using the waso_window
          night_dat$SPT_roll_log_vedba <- rollmedian( night_dat$log_vedba, waso_window, fill = NA, align = 'center' )
          
          ## plot this rolling median
          #lines( ts_time, night_dat$SPT_roll_log_vedba, col = 'red', lty = 1 )
          
          ## calculate the threshold for sleeping and waking within the sleep period
          SPT_thresh <- quantile( night_dat$SPT_roll_log_vedba, waso_percentile, na.rm = T) * waso_multiplier
          
          ## plot the threshold
          abline( h = SPT_thresh, col = 'blue', lty = 1, lwd = 2 )
          
          ## subset the night's data to only the sleep period
          trim_night <- night_dat[ night_dat$local_timestamp >= night_dat$local_timestamp[ onset[ per_ind ] ] & night_dat$local_timestamp <= night_dat$local_timestamp[ wake[ per_ind ] ] , ]
          
          ### find blocks of continuous inactivity
          
          ## calcuate the run length encoding
          temp <- rle(as.numeric( trim_night$SPT_roll_log_vedba < SPT_thresh ) ) 
          
          ## mark the runs of activity or inactivity
          trim_night$night_runs <- as.numeric( rep( temp$lengths >= waso_block, times = temp$lengths ) )
          
          ## mark the runs of inactivity as sleep bouts
          trim_night$night_sleep_bouts <- as.numeric( trim_night$SPT_roll_log_vedba < SPT_thresh & trim_night$night_runs == 1 )
          
          ## find the starts and ends of waking bouts
          diffs <- diff( c(1, trim_night$night_sleep_bouts ) )
          
          starts <- which( diffs == -1 )
          
          ## add back in a "- 1" at the end of this line if you wish for the start and end times to be accurate. Now I just want to make it so the polygons show up even without a border
          ends <- which( diffs == 1 )
          
          ## if there are waking bouts
          if( length( which( diffs == -1 ) ) != 0){
            
            ## if the last waking bout never ends...
            if( length( starts ) != length( ends ) ){
              
              ## make it end at the end of the sleep period
              ends <- c(ends, length( diffs) )
              
            }
            
            ## save the start times and end times of waking bouts
            starts <- ts_func( trim_night$local_time[ starts ] )
            ends <- ts_func( trim_night$local_time[ ends ] )
            
            ## plot a polygon for each distinct waking bout
            for( n in 1:length( starts ) ){
              
              polygon( x = c(starts[ n ], ends[ n ], ends[ n ], starts[ n ], starts[ n ] ), y = c( 0, 0, 10, 10, 0), col = transp('blue', .25), border = NA )
              
            }
          }
        }
        ## fill in the sleep period data frame with the sleep onset and waking time associated with the longest sleep period in the day (noon to noon)
        return( c( night_dat$local_timestamp[ onset[ per_ind ] ], night_dat$local_timestamp[ wake[ per_ind ] ] ) )
      } 
    }
  } 
}




################## Read in the d1 (accelerometer burst) data ###################

## d1 is a dataframe with a row for each minute for each baboon. Each row contains the raw (or interpolated) GPS acc burst, and several different measures calculated from bursts (like VeDBA)
d1 <- fread( "DATA/sleep_analysis/full_night_and_day_data.csv" )

## turn the data table into a dataframe
d1 <- as.data.frame( d1 )

## turn timestamp into POSIX element and time into character
d1$local_timestamp <- as.POSIXct( d1$local_timestamp, tz = 'UTC' )

# round local_timestamps to the nearest minute
d1$local_timestamp <- floor_date( d1$local_timestamp, unit = 'minute' )

# recreate the time column with the rounded times
d1$time <- str_split_fixed( d1$local_timestamp, ' ', 2 )[ , 2 ]

## specify that the time is local_time
names( d1 )[ names( d1 ) == 'time' ] <- 'local_time'

## assign each minute of data to a given night. A night lasts from noon to noon. First, apply a time shift so that each night is a unit, and not each day
time_shift <- d1$local_timestamp - 12*60*60

## save the date of the first night of the study (the date of the night is always the date of the evening at the beginning of that night; so the first night of the study is 2012-07-31, although the data starts on 2012-08-01, because the data on that first morning is still technically part of the data for the previous night, as a night is noon to noon)
start_date <- as.Date( min( d1$local_timestamp )- 12*60*60 )

## assign night as number of nights from the start of the study, with all data before the first noon representing night 1
d1$night_num <- as.numeric( as.Date( time_shift ) - start_date + 1 )

## show how many baboon-nights there are
nrow( unique( d1[ , c( 'tag', 'night' ) ] ) )

## check where the night changes from one night to the next to see if it is at noon
d1[ ( diff( c( d1$night_num ) ) == 1),]
d1[ ( diff( c( as.Date( d1$night ) ) ) == 1),]
identical( d1[ ( diff( c( d1$night_num ) ) == 1),], d1[ ( diff( c( as.Date( d1$night ) ) ) == 1),] )

## this next line has already been done in the previous script, so it should just be redundant here
# remove LI TH's data when she is dying
d1 <- d1[ !( d1$tag == 'Pa LI TH' & d1$local_timestamp >= as.POSIXct( '2014-05-28 00:00:00', tz = 'UTC' ) ), ]

nrow( unique( d1[ , c( 'tag', 'night' ) ] ) )


## view the dataframe and it's summary
head(d1)

summary(d1)

sum( duplicated( d1[ , c( 'tag', 'local_timestamp' ) ] ) ) # just confirms that there are no duplicate timestamps even after the floor_date rounding above

# #################### Read in the GPS data #######################

bdat <- fread( "DATA/bab_complete.csv", fill = T )

## turning bdat into a dataframe from a data table
bdat <- as.data.frame( bdat )

## make the local_timestamp into a POSIXct element
bdat$local_timestamp <- as.POSIXct(bdat$local_timestamp, tz = 'UTC' )

## this next line is taken care of by a previous script so it is redundant
# remove LI TH's data when she is dying
bdat <- bdat[ ! ( bdat$id == 'Pa_LI_TH' & bdat$local_timestamp >= as.POSIXct( '2014-05-28 00:00:00', tz = 'UTC') ), ]

## make a column named day that represents the day from the start of the study period, with the first day being day 1
bdat$day <- as.Date( bdat$local_timestamp, tz = 'UTC' )

# check where the day changes from one day to the next
bdat[ ( diff( c( as.Date( bdat$day ) ) ) == 1),]

## make a column named time that will have just the time component, and not the date, of the local_timestamp
bdat$local_time <- str_split_fixed( bdat$local_timestamp, " ",2)[,2]

## remove rows of the dataframe that don't have a fix
bdat <- bdat[!is.na( bdat$lat ),]

## open the first rows of the dataframe
head( bdat )

## find the average longitude and latitude of the study GPS points. These will be used to find the times of sunset and sunrise
ave_lon <- mean( bdat$lon )
ave_lat <- mean( bdat$lat )

## look at how many 'baboon days' there are in the data
nrow( unique( bdat[ , c( 'id', 'day' ) ] ) )


############## Make a dataframe of sunrise and sunset times ###############

## make an empty dataframe with each date of the study. For each date, we will fill in when sunset occurred, when the dark period began (the end of evening astronomical twilight), when the dark period ended the following morning (the beginning of morning astronomical twilight), and when sunrise the following morning occured
sun_dat <- data.frame( date = c( ( min( as.Date( d1$local_timestamp ) ) - 2:1 ), unique( as.Date( d1$local_timestamp ) ) ), sunset = NA, night_start = NA, night_end = NA, sunrise = NA )

## fill in sunset time and night start time (dark period start time) on each date with the getsunlighttimes function
sun_dat[, c( 'sunset', 'night_start' ) ] <- getSunlightTimes( date = sun_dat$date, lon = ave_lon, lat = ave_lat )[, c( 'sunset', 'night' ) ]

## fill in rise time and night end time (dark period end time) on each date with those from the following date with the getsunlighttimes function. The reason we are using the following date is because we want to know when a night ended, which happens on the date following the start of that night
sun_dat[, c( 'sunrise', 'night_end' ) ] <- getSunlightTimes( date = ( sun_dat$date + 1 ), lon = ave_lon, lat = ave_lat )[, c( 'sunrise', 'nightEnd' ) ]

## put sun data in local time
sun_dat[ , 2:5 ] <- sun_dat[ , 2:5 ] + 3*60*60

############## Make a dataframe of moon phases ###############

## save the unique dates associated with the local timestamps in d1
dates <- c( ( min( as.Date( d1$local_timestamp ) ) - 2:1 ), unique( as.Date( d1$local_timestamp ) ) )

## using the getmoonillumination function, get information about the moon on the study dates
moon_dat <- getMoonIllumination( date = c( min( dates ) - 1, dates, max( dates ) + 1 ) )

####################### Determining sleep periods with modification of Van Hees et al. 2018 method ####################################


## save a variable denoting the total number of minutes in the day
mins_in_day <- 60*24

missing_mins <- 120 ## this is the maximum total number of minutes of data that can be missing from a day and still have that day included in the analysis (for sleep period time and sleep based analyses; i.e. not ave_vedba)

night_missing_mins <- 45 ## this is the maximum total number of minutes of data that can be missing from a dark period and still have that day included in the analysis

time_gap <- 20*60 ## this is the maximum allowable time gap between two accelerometer bursts (in seconds) that can exist in a noon-to-noon period without removing this noon-to-noon period from the data

mov_window <- 9 ## this is the size of the moving window (in minutes) used in calculating the rolling median of the average VeDBA ## I had this at 17 due to visual inspection but if we want it to match the validated algorithm, we need to keep it at 9

percentile <- 0.10 ## this is the percentile threshold of the variable (within the noon-to-noon time period) used to classify activity vs. inactivity (VeDBA here)

multiplier <- 1.125 ## this is the multiplier of the threshold value determined by the percentile above. Values below the threshold value multiplied by the multiplier are considered inactive and those equal to or above are considered active

block_size <- 30 ## duration in minutes of the blocks of continuous inactivity that will be considered sleep

gap_size <- 45 ## maximum duration between sleep blocks that will be merged ## I had this at 60 due to visual inspection but if we want it to match the validated algorithm, we need to keep it at 45

waso_percentile <- 0.10 ## this is the percentile threshold of the variable (within the noon-to-noon time period) used to classify activity vs. inactivity within the sleep period

waso_multiplier <- 1.125 ## this is the multiplier of the threshold value determined by the percentile above. Values WITHIN the sleep period below the threshold value multiplied by the multiplier are considered inactive and those equal to or above are considered active. This threshold value may be different than the one above even when waso_percentile = percentile and waso_multiplier = multiplier, because the value of the given percentile of this variable depends on the smoothing window (see waso_window below; aka waso window might not be equal to mov_window)

waso_window <- 1 ## this is the size of the moving window (in minutes) used in calculating the rolling median that will be used to find periods of wake after sleeping. A waso_window of 1 is the same as using the raw data without a rolling median

waso_block <- 3 ## this is the number of consecutive minutes of inactivity needed to classify a period as sleep within the sleep period. A waso_block of 1 means that anytime the value is below the threshold, the baboon in considered sleeping and anytime the value is above the threshold the baboon is considered awake

nap_block <- waso_block ## this is the number of consecutive minutes of inactivity needed to classify a period as sleep during the waking period. A nap_block of 1 means that anytime the value is below the threshold, the baboon in considered sleeping and anytime the value is above the threshold the baboon is considered awake

frag_block <- 2 ## this is the number of minutes of waking that need to be consecutive to be considered a wake bout during the night (other epochs of wake that do not meet this criterion will still be considered wake for WASO and wake_bouts, but not frag_wake_bouts)

hist( as.numeric( as_hms( bdat$arrive_sleep_site ) ), xaxt = 'n' )
axis( 1, at = seq( 0, 60*24*60, 60*60), labels = as_hms( seq( 0, 60*24*60, 60*60) ) )
abline( v = median( as.numeric( as_hms( bdat$arrive_sleep_site ) ), na.rm = T ) )
as_hms( median( as.numeric( as_hms( bdat$arrive_sleep_site ) ), na.rm = T ) )

hist( as.numeric( as_hms( bdat$leave_sleep_site ) ), xaxt = 'n' )
axis( 1, at = seq( 0, 60*24*60, 60*60), labels = as_hms( seq( 0, 60*24*60, 60*60) ) )
abline( v = median( as.numeric( as_hms( bdat$leave_sleep_site ) ), na.rm = T ) )
as_hms( median( as.numeric( as_hms( bdat$leave_sleep_site ) ), na.rm = T ) )


median_sleep_site_arrive_time <- as.character( as_hms( median( as.numeric( as_hms( bdat$arrive_sleep_site ) ), na.rm = T ) ) )

median_sleep_site_leave_time <- as.character( as_hms( median( as.numeric( as_hms( bdat$leave_sleep_site ) ), na.rm = T ) ) )


sep_day_night <- F ## this determines whether sleep periods and non-sleep periods are separated before finding runs of inactivity to consider as sleep

## shows the time (as well as one previous time and one later time) where a minute is skipped
sort( unique( d1$local_time ) ) [ which( diff( as_hms( sort( unique( d1$local_time ) ) ) ) != as_hms( '00:01:00' ) ) + -1:1 ]

## again confirms that every minute is represented in the data except for one (can tell this by comparing this number to the minutes_in_day variable above)
length( unique(d1$local_time) )

## create a vector containing the names of each baboon
tag_names <- unique( d1$tag )

## make a copy of d1. We will fill in this new dataframe with information about if the baboon was asleep in each epoch
full_dat <- d1

full_dat$sleep_per <- NA ## binary indicating whether a row belongs to the sleep period window
full_dat$sleep_bouts <- NA ## binary indicating whether the row is considered sleep, based on the waso or nap requirements
full_dat$n_bursts <- NA ## the number of bursts collected in a given noon-to-noon period (to be compared to the total number of minutes in the day). This column will indicate whether the data for a given night is insufficient to calculate the sleep period (and thus: onset, waking, SPT, sleep_eff, TST, sleep_bouts -- because this depends on a percentile of bursts' log vedba, WASO, wake_bouts, summed_VeDBA, night_VeDBA_corr, dark_TST, prev_naps, prev_day_sleep)
full_dat$max_time_diff <- NA ## the maximum difference between consecutive fixes in a given noon-to-noon period. With the previous column, this column will indicate whether the data is insufficient to calculate the sleep period (and thus: onset, waking, SPT, sleep_eff, TST, WASO, wake_bouts, summed_VeDBA, night_VeDBA_corr, prev_naps )
full_dat$dark <- NA ## a binary indicating whether a row belongs to the period of darkness (between the end of astrological evening twilight and the beginning of astrological morning twilight)
full_dat$poss_dark_bursts <- NA ## the number of potential bursts in the dark period on this night, i.e. the number of minutes between the start and end of the night
full_dat$n_dark_bursts <- NA ## the total number of bursts actually taken during the dark period on this night. This will be compared to the previous column to determine whether the data during the dark period is sufficient to calculate the dark_TST and ave_vedba
full_dat$max_dark_time_diff <- NA ## the maximum difference between consecutive fixes in a given dark period. With the previous column, this column will indicate whether the data is insufficient to calculate the dark_TST and ave_vedba
full_dat$poss_day_bursts <- NA ## the number of potential bursts in preceding day (light) period, i.e. the number of minutes between the end of morning astrological twilight (maybe this should be changed to sunrise?) and the start of evening astrological twilight (maybe this should be changed to sunset?)
full_dat$n_day_bursts <- NA ## the total number of bursts actually taken during the preceding day (light) period. This will be compared to the previous column to determine whether the data during the day period is sufficient to calculate the prev_day_sleep and prev_day_ave_vedba.
full_dat$max_day_time_diff <- NA ## the maximum difference between consecutive fixes in a given day (light) period. With the previous column, this column will indicate whether the data is insufficient to calculate the prev_day_sleep and prev_day_ave_vedba
## prev_naps depends on having sufficient data to calculate both the current SPT as well as the previous day's SPT


## create a vector containing the names of each baboon
tag_names <- unique( d1$tag )

## for each individual...
for( tag in tag_names ){
  
  ## subset the data to just this individual's data
  id_dat <- d1[ d1$tag == tag, ]
  
  ## create a vector the nights for which this individual has data
  nights <- unique( id_dat$night )
  
  ## for each night on which this individual has data
  for( night in nights ){
    
    ## subset this individual's data to just that night
    night_dat <- id_dat[ id_dat$night == night, ]
    
    ## create empty columns for the sleep period and sleep bout binary variables
    night_dat$sleep_per <- NA
    night_dat$sleep_bouts <- NA
    
    
    ## save a column of the total number of bursts for that day. This will also make it easier to remove these days from the dataframe later
    night_dat$n_bursts <- nrow( night_dat )
    
    ## sort the timestamps (they are probably already sorted)
    sorted_times <- c( '00:00:00', sort( night_dat$local_time ), '23:59:00' )
    
    ## find the time difference in seconds between each burst
    time_diffs <- as.numeric( diff( as_hms( sorted_times ) ), units = 'secs' )
    
    if( length( time_diffs ) > 0 ){ ### There is one night for one baboon with only one single burst, which is why this if statement is needed
      
      ## save a column of the maximum time difference between burst for that day (this will make it easier to pull out days with insufficient data later)
      night_dat$max_time_diff <- max( time_diffs )
      
    }else{
      
      night_dat$max_time_diff <- NA
      
    }
    
    ## take the rolling median of the log VeDBA and save it as a column
    roll_log_vedba <- rollmedian( night_dat$log_vedba, mov_window, fill = NA, align = 'center' )
    
    ## determine the threshold activity vs. inactivity threshold based on the percentile, multiplier, and the rolling median just produced
    thresh <- quantile( roll_log_vedba, percentile, na.rm = T) * multiplier
    
    ### find blocks of continuous inactivity
    
    ## find the run length encoding of periods above and below the threshold
    temp <- rle(as.numeric( roll_log_vedba < thresh ) )
    
    ## mark the rows that are part of runs (i.e. part of chunks that are greater than the block_size of either continuous activity or continuous inactivity )
    sleep_per_runs <- as.numeric( rep( temp$lengths > block_size, times = temp$lengths ) )
    
    ## mark the rows corresponding to sleep bouts. These sleep bouts are runs of inactivity
    sleep_per_sleep_bouts <- as.numeric( roll_log_vedba < thresh & sleep_per_runs == 1 )
    
    ## find when sleep bouts start and end
    diffs <- diff( c(0, sleep_per_sleep_bouts ) )
    starts <- which( diffs == 1 ) [ -1 ]
    ends <- which( diffs == -1 )
    
    ## if there are any sleep bouts...
    if( length( which( diffs == 1 ) ) != 0 ){
      
      ## find the duration of the gaps between each sleep bout (the end of one sleep bout and the start of the next)
      gaps <- as.numeric( night_dat$local_timestamp [ starts ] - night_dat$local_timestamp [ ends[ 1: length( starts ) ] ], units = 'mins' )
      
      ## sleep bouts separated by gaps that are shorter than that specified by gap_size will be merged. Note which of these gaps are shorter than the gap_size
      inds_to_remove <- which( gaps < gap_size )
      
      ## if there are NO gaps between sleep bouts that are to be removed...
      if( length( inds_to_remove ) == 0 ){
        
        ## set sleep onset index to be the start of sleep bouts
        onset <- which( diffs == 1 )
        
        ## set waking index to be the end of sleep bouts
        wake <- ends
        
      }else{ ## if there ARE gaps between sleep bouts that are to be removed...
        
        ## set sleep onset index to be the start of sleep bouts that do not correspond to the gaps to be removed (because these will be within sleep periods, not a start of a new bout)
        onset <- which( diffs == 1 ) [ - (inds_to_remove + 1) ]
        
        ## set waking index to be the end of sleep bouts that do not correspond to the gaps to be removed
        wake <- ends [ - inds_to_remove ]
        
      }
      
      ## determine which sleep period is the longest
      per_ind <- which.max( as.numeric( night_dat$local_timestamp[ wake ] - night_dat$local_timestamp[ onset ], units = 'secs' ) )
      
      ## fill in the sleep period data frame with the sleep onset and waking time associated with the longest sleep period in the day (noon to noon)
      
      night_dat$sleep_per <- as.numeric( night_dat$local_timestamp >= night_dat$local_timestamp[ onset[ per_ind ] ] & night_dat$local_timestamp <= night_dat$local_timestamp[ wake[ per_ind ] ] )
      
    }else{ ## if there aren't any sleep bouts, record all rows as a 0 in the sleep_per column
      
      night_dat$sleep_per <- 0
      
    }
    
    ## take the rolling median of the log VeDBA
    night_dat$roll_log_vedba <- rollmedian( night_dat$log_vedba, waso_window, fill = NA, align = 'center' )
    
    ### find blocks of continuous inactivity
    
    if( sep_day_night == T ){
      
      # first during the night. Splitting this up by day and night will allow us to use different threshold dor the duration of inactivity required to be considered sleep between day and night
      
      ## subset the night data to only the sleep period window
      trim_night <- night_dat[ night_dat$sleep_per == 1, ]
      
      ## find the run length encoding of periods above and below the threshold
      temp <- rle(as.numeric( trim_night$roll_log_vedba < thresh ) )
      
      ## mark which rows are part of runs of activity vs. inactivity, as determined by waso_block
      runs <- as.numeric( rep( temp$lengths >= waso_block, times = temp$lengths ) )
      
      sleep_bouts_night <- as.numeric( trim_night$roll_log_vedba < thresh & runs == 1 )
      
      ## mark which rows are part of runs of inactivity. These are the periods of sleep within the sleep period
      trim_night$sleep_bouts <- sleep_bouts_night
      
      ## return this sleep period data back into night_dat
      night_dat[ night_dat$sleep_per == 1, ] <- trim_night
      
      if( sum( night_dat$sleep_per ) > 0 ){ ## if there are sleep periods
        
        ## subset the night data to only the times outside the sleep period window
        trim_day_before <- night_dat[ night_dat$sleep_per != 1 & night_dat$local_timestamp < min( night_dat$local_timestamp[ night_dat$sleep_per == 1 ] ), ]
        
        ## find the run length encoding of periods above and below the threshold
        temp <- rle(as.numeric( trim_day_before$roll_log_vedba < thresh ) )
        
        ## mark which rows are part of runs of activity vs. inactivity, as determined by nap_block
        runs <- as.numeric( rep( temp$lengths >= nap_block, times = temp$lengths ) )
        
        
        sleep_bouts_day_before <- as.numeric( trim_day_before$roll_log_vedba < thresh & runs == 1 )
        
        ## mark which rows are part of runs of inactivity. These are the periods of sleep within the sleep period
        trim_day_before$sleep_bouts <- sleep_bouts_day_before
        
        ## return this data from outside the sleep period window back to night_data
        night_dat[ night_dat$sleep_per != 1 & night_dat$local_timestamp < min( night_dat$local_timestamp[ night_dat$sleep_per == 1 ] ), ] <- trim_day_before
        
        
        ## subset the night data to only the times outside the sleep period window
        trim_day_after <- night_dat[ night_dat$sleep_per != 1 & night_dat$local_timestamp > max( night_dat$local_timestamp[ night_dat$sleep_per == 1 ] ), ]
        
        ## find the run length encoding of periods above and below the threshold
        temp <- rle(as.numeric( trim_day_after$roll_log_vedba < thresh ) )
        
        ## mark which rows are part of runs of activity vs. inactivity, as determined by nap_block
        runs <- as.numeric( rep( temp$lengths >= nap_block, times = temp$lengths ) )
        
        sleep_bouts_day_after <- as.numeric( trim_day_after$roll_log_vedba < thresh & runs == 1 )
        
        ## mark which rows are part of runs of inactivity. These are the periods of sleep within the sleep period
        trim_day_after$sleep_bouts <- sleep_bouts_day_after
        
        ## return this data from outside the sleep period window back to night_data
        night_dat[ night_dat$sleep_per != 1 & night_dat$local_timestamp > max( night_dat$local_timestamp[ night_dat$sleep_per == 1 ] ), ] <- trim_day_after
        
        
      }else{
        
        ## subset the night data to only the times outside the sleep period window
        trim_day <- night_dat[ night_dat$sleep_per != 1, ]
        
        ## find the run length encoding of periods above and below the threshold
        temp <- rle(as.numeric( trim_day$roll_log_vedba < thresh ) )
        
        ## mark which rows are part of runs of activity vs. inactivity, as determined by nap_block
        runs <- as.numeric( rep( temp$lengths >= nap_block, times = temp$lengths ) )
        
        sleep_bouts_day <- as.numeric( trim_day$roll_log_vedba < thresh & runs == 1 )
        
        ## mark which rows are part of runs of inactivity. These are the periods of sleep within the sleep period
        trim_day$sleep_bouts <- sleep_bouts_day
        
        
        ## return this data from outside the sleep period window back to night_data
        night_dat[ night_dat$sleep_per != 1, ] <- trim_day
        
      }
      
    }else{
      
      ## find the run length encoding of periods above and below the threshold
      temp <- rle(as.numeric( night_dat$roll_log_vedba < thresh ) )
      
      ## mark which rows are part of runs of activity vs. inactivity, as determined by waso_block
      runs <- as.numeric( rep( temp$lengths >= waso_block, times = temp$lengths ) )
      
      sleep_bouts <- as.numeric( night_dat$roll_log_vedba < thresh & runs == 1 )
      
      ## make which rows are part of runs of inactivity. These are the periods of sleep within the sleep period
      night_dat$sleep_bouts <- sleep_bouts
      
    }
    
    
    
    dark_start <- sun_dat$night_start[ sun_dat$date == night ]
    
    dark_end <- sun_dat$night_end[ sun_dat$date == night ]
    
    night_dat$dark <- as.numeric( night_dat$local_timestamp >= dark_start & night_dat$local_timestamp <= dark_end )
    
    night_dat$poss_dark_bursts <- floor( as.numeric( dark_end - dark_start, units = 'mins' ) )
    
    night_dat$n_dark_bursts <- sum( night_dat$local_timestamp >= dark_start & night_dat$local_timestamp <= dark_end )
    
    ## sort the timestamps (they are probably already sorted)
    sorted_times <- c( dark_start, sort( night_dat[ night_dat$dark == 1, 'local_timestamp'] ), dark_end )
    
    ## find the time difference in seconds between each burst
    time_diffs <- as.numeric( diff( sorted_times ), units = 'secs' )
    
    if( length( time_diffs ) != 0 ){
      ## save the maximum time difference during the night
      night_dat$max_dark_time_diff <- max( time_diffs )  ## so if this variable is NA, that means that there are no bursts taken from the nighttime period on this day for this individual
    }else{
      
      night_dat$max_dark_time_diff <- NA
    }
    
    
    
    first_light <- sun_dat$night_end[ sun_dat$date == ( night - 1 ) ]
    
    night_dat$poss_day_bursts <- floor( as.numeric( dark_start - first_light, units = 'mins' ) )
    
    day_dat <- id_dat[ id_dat$local_timestamp > first_light & id_dat$local_timestamp < dark_start, ]
    
    night_dat$n_day_bursts <- nrow( day_dat )
    
    
    ## sort the timestamps (they are probably already sorted)
    sorted_times <- c( first_light, sort( day_dat$local_timestamp ), dark_start )
    
    ## find the time difference in seconds between each burst
    time_diffs <- as.numeric( diff( sorted_times ), units = 'secs' )
    
    if( length( time_diffs ) != 0 ){
      ## save the maximum time difference during the night
      night_dat$max_day_time_diff <- max( time_diffs )  ## so if this variable is NA, that means that there are no bursts taken from the nighttime period on this day for this individual
    }else{
      
      night_dat$max_day_time_diff <- NA
    }
    
    
    
    
    
    night_dat <- night_dat[ ,  names( night_dat ) != 'roll_log_vedba' ]
    
    
    ### put the night data back into full_dat
    full_dat[ full_dat$tag == tag & full_dat$night == night, ] <- night_dat
  }
  
}

pre_clean_full <- full_dat

study_nights <- as.Date( min( d1$night ):max( d1$night ) )

day_lim_start_time <- "07:30:00"

day_lim_end_time <- "17:30:00"

sleep_per <- data.frame( tag = rep( unique( d1$tag ), each = length( study_nights ) ), night = rep( study_nights, times = length( tag ) ), onset = NA, waking = NA, SPT = NA, WASO = NA, TST = NA, TST_at_sleep_site = NA, sleep_eff = NA, wake_bouts = NA, frag_wake_bouts = NA, frag_wake_bouts_in_median_SPT = NA, frag_wake_bouts_at_sleep_site = NA, summed_VeDBA = NA, night_VeDBA_corr = NA, ave_vedba = NA, dark_TST = NA, light_TST = NA, dark_sleep_eff = NA, light_sleep_eff = NA, prev_naps = NA, prev_day_sleep = NA, prev_day_sleep_lim = NA, prev_day_ave_vedba = NA, poss_dark_bursts = NA, n_dark_bursts = NA, poss_day_bursts = NA, n_day_bursts = NA, max_time_diff = NA, n_bursts= NA, max_dark_time_diff = NA, max_day_time_diff = NA )


## create empty vectors for the durations of sleep and wake bouts. We will fill these in to see if the distributions of the durations of these bouts later
sleep_durs <- c()
wake_durs <- c()


## for each individual...
for( tag in tag_names ){
  
  ## subset the data to just this individual's data
  id_dat <- full_dat[ full_dat$tag == tag, ]
  
  ## create a vector the nights for which this individual has data
  nights <- unique( id_dat$night )
  
  ## for each night on which this individual has data
  for( night in nights ){
    
    ## subset this individual's data to just that night
    night_dat <- id_dat[ id_dat$night == night, ]
    
    ## should already be in order, but just in case
    night_dat <- night_dat[ order( night_dat$local_timestamp ), ]
    
    sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$n_bursts <- unique( night_dat$n_bursts )
    
    sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$max_time_diff <- unique( night_dat$max_time_diff )
    
    sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$poss_dark_bursts <- unique( night_dat$poss_dark_bursts )
    
    sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$n_dark_bursts <- unique( night_dat$n_dark_bursts )
    
    sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$max_dark_time_diff <- unique( night_dat$max_dark_time_diff )
    
    sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$poss_day_bursts <- unique( night_dat$poss_day_bursts )
    
    sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$n_day_bursts <- unique( night_dat$n_day_bursts )
    
    sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$max_day_time_diff <- unique( night_dat$max_day_time_diff )
    
    
    SPT_dat <- night_dat[ night_dat$sleep_per == 1, ]
    
    if( nrow( SPT_dat ) > 0 ){
      
      onset <- min( SPT_dat$local_timestamp )
      
      sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$onset <- onset
      
      waking <- max( SPT_dat$local_timestamp )
      
      sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$waking <- waking
      
      SPT <- as.numeric( waking - onset, units = 'mins' ) + 1
      
      sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$SPT <- SPT
      
      WASO <- sum( SPT_dat$sleep_bouts == 0 )
      
      sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$WASO <- WASO
      
      TST <- sum( SPT_dat$sleep_bouts == 1 )
      
      sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$TST <- TST
      
      sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$sleep_eff <- TST/ nrow( SPT_dat )
      
      sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$summed_VeDBA <- sum( SPT_dat$log_vedba )
      
      sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$night_VeDBA_corr <- sum( SPT_dat$log_vedba ) / SPT
      
      temp <- rle( SPT_dat$sleep_bouts )
      
      runs <- as.numeric( rep( temp$lengths >= frag_block, times = temp$lengths ) )
      
      frag_wake_bouts <- as.numeric( SPT_dat$sleep_bouts == 0 & runs == 1 )
      
      diffs <- diff( c( 1, frag_wake_bouts ) )
      
      sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$frag_wake_bouts <- sum( diffs == 1 )
      
      ## find the distinct sleep bouts (i.e. epochs of sleep separated by waking)
      diffs <- diff( c( 0, SPT_dat$sleep_bouts ) )
      
      ## save the number of distinct wake bouts
      sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$wake_bouts <- sum( diffs == -1 )
      
      ## find durations of sleep and wake bouts
      temp <- rle( SPT_dat$sleep_bouts )
      
      ## add the duration of sleep bouts to the sleep bout duration vector
      sleep_durs <- c( sleep_durs, temp$lengths[ temp$values == 1 ] )
      ## add the duration of wake bouts to the wake bout duration vector
      wake_durs <- c( wake_durs, temp$lengths[ temp$values == 0 ] )
      
      
      waking_dat <- id_dat[ id_dat$local_timestamp > sleep_per[ sleep_per$tag == tag & sleep_per$night == ( night - 1 ), ]$waking & id_dat$local_timestamp <  onset   , ]
      
      sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$prev_naps <- sum( waking_dat$sleep_bouts )
      
      
    }
    
    
    if( !is.null( median_sleep_site_arrive_time ) ){
      
      median_sleep_site_dat <- night_dat[ night_dat$local_time >= median_sleep_site_arrive_time | night_dat$local_time <= median_sleep_site_leave_time, ]
      
      TST_at_sleep_site <- sum( median_sleep_site_dat$sleep_bouts == 1 )
      
      sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$TST_at_sleep_site <- TST_at_sleep_site
      
      # finding fragmentation within the median sleep period
      temp <- rle( median_sleep_site_dat$sleep_bouts )
      
      runs <- as.numeric( rep( temp$lengths >= frag_block, times = temp$lengths ) )
      
      frag_wake_bouts_in_med_sleep_site <- as.numeric( median_sleep_site_dat$sleep_bouts == 0 & runs == 1 )
      
      diffs <- diff( c( 1, frag_wake_bouts_in_med_sleep_site ) )
      
      sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$frag_wake_bouts_at_sleep_site <- sum( diffs == 1 )
      
    }
    
    dark_dat <- night_dat[ night_dat$dark == 1, ]
    
    sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$ave_vedba <- mean( dark_dat$log_vedba )
    
    sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$dark_TST <- sum( dark_dat$sleep_bouts == 1 )
    
    sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$dark_sleep_eff <- sum( dark_dat$sleep_bouts == 1 ) / nrow( dark_dat )
    
    light_dat <- night_dat[ night_dat$dark == 0 & night_dat$sleep_per == 1, ]
    
    sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$light_TST <- sum( light_dat$sleep_bouts == 1 )
    
    sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$light_sleep_eff <- sum( light_dat$sleep_bouts == 1 ) / nrow( dark_dat )
    
    
    first_light <- sun_dat$night_end[ sun_dat$date == ( night - 1 ) ]
    
    last_light <- sun_dat$night_start[ sun_dat$date == night ]
    
    day_dat <- id_dat[ id_dat$local_timestamp > first_light & id_dat$local_timestamp < last_light, ]
    
    sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$prev_day_sleep <- sum( day_dat$sleep_bouts )
    
    sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$prev_day_ave_vedba <- mean( day_dat$log_vedba )
    
    day_lim_start <- as.POSIXct( paste( str_split_fixed( first_light, " ", 2 )[ , 1 ], day_lim_start_time ), tz = 'UTC' )
    
    day_lim_end <- as.POSIXct( paste( str_split_fixed( last_light, " ", 2 )[ , 1 ], day_lim_end_time ), tz = 'UTC' )
    
    day_lim <- id_dat[ id_dat$local_timestamp > day_lim_start & id_dat$local_timestamp < day_lim_end, ]
    
    sleep_per[ sleep_per$tag == tag & sleep_per$night == night, ]$prev_day_sleep_lim <- sum( day_lim$sleep_bouts )
    
  }
}


sleep_per

# make a column for the number of nights since the beginning of the study, just for easy notation and subsetting later
sleep_per$night_num <- as.numeric( sleep_per$night - min( sleep_per$night ) + 1 )

## reformat sleep timestamp
sleep_per$onset <- as.POSIXct( sleep_per$onset, origin = "1970-01-01 00:00:00", tz = "UTC" )

## reformat waking timestamp
sleep_per$waking <- as.POSIXct( sleep_per$waking, origin = "1970-01-01 00:00:00", tz = "UTC" )


########### Add in moon phase data ###########

sleep_per <- merge( x = sleep_per, y = moon_dat, by.x = 'night', by.y = 'date', all.x = T, all.y = F, sort = F )

sum( !is.na( sleep_per$SPT ) )

### check number of nights for which sleep period was calculated and inspect those for which no sleep period was calculated ###

sleep_per_nona <- sleep_per[ !is.na( sleep_per$SPT ), ]

# nrow( sleep_per_nona )
#
# left_out <- unique( d1[ , c( 'tag', 'night' ) ] )[ !paste( unique( d1[ , c( 'tag', 'night' ) ] )$tag, unique( d1[ , c( 'tag', 'night' ) ] )$night ) %in% paste( sleep_per_nona$tag, sleep_per_nona$night ), ]
#
#
# for( i in 1:nrow( left_out ) ){
#
#   tag_night_dat <- d1[ d1$tag == left_out$tag[ i ] & d1$night == left_out$night[ i ], ]
#
#   plot( tag_night_dat$local_timestamp, tag_night_dat$log_vedba )
#
# }
#
# nrow( left_out )
#
# tag_night_dat <- d1[ d1$tag == tag_names[ 1 ] & d1$night_num == 5, ]
#
# plot( tag_night_dat$local_timestamp, tag_night_dat$log_vedba )
#


############# Cleaning the dataframes of data on nights with insufficient data ################

## for days on the later side of a noon-to-noon period with a lot of missing data, we can't have a reliable threshold for what was considered sleep in the morning, and we may be missing a lot of epochs. Therefore, remove: sleep time during the day (day defined by SPT (wake time to onset time)), sleep time during the day (day defined as prev dark period end to dark period start), sleep time during the day (defined as 0730 to 1730)

sleep_per[ paste( sleep_per$tag, sleep_per$night, sep = '_' ) %in% paste( sleep_per[ sleep_per$n_bursts < ( mins_in_day - missing_mins ) & !is.na( sleep_per$n_bursts ), 'tag' ], ( sleep_per[ sleep_per$n_bursts < ( mins_in_day - missing_mins ) & !is.na( sleep_per$n_bursts ) , 'night' ] + 1 ), sep = '_' ), c( 'prev_naps', 'prev_day_sleep', 'prev_day_sleep_lim' ) ] <- NA

## for days on the later side of a noon-to-noon period with large gaps of missing data, we can't have a reliable waking time. Therefore, remove: sleep time during the day (day defined by SPT (wake time to onset time))
sleep_per[ paste( sleep_per$tag, sleep_per$night, sep = '_' ) %in% paste( sleep_per[ sleep_per$max_time_diff > time_gap & !is.na( sleep_per$max_time_diff ), 'tag' ], ( sleep_per[ sleep_per$max_time_diff > time_gap & !is.na( sleep_per$max_time_diff ) , 'night' ] + 1 ), sep = '_' ), c( 'prev_naps' ) ] <- NA

## remove all these variable from the night, and from the days on the early side of the noon-to-noon period if the noon-to-noon period is missing a lot of data (because then we might not be able to reliably calculate the sleep VeDBA threshold, and a lot of epochs might be missing, which would skew TST and such)
sleep_per[ sleep_per$n_bursts < ( mins_in_day - missing_mins ) & !is.na( sleep_per$n_bursts ), c( 'onset', 'waking', 'SPT', 'sleep_eff', 'TST', 'WASO', 'wake_bouts', 'summed_VeDBA', 'night_VeDBA_corr', 'dark_TST', 'dark_sleep_eff', 'light_TST', 'light_sleep_eff', 'prev_naps', 'prev_day_sleep', 'prev_day_sleep_lim', 'TST_at_sleep_site', 'frag_wake_bouts', 'frag_wake_bouts_in_median_SPT', 'frag_wake_bouts_at_sleep_site' ) ] <- NA

sleep_per_nona <- sleep_per[ !is.na( sleep_per$SPT ), ]

nrow( sleep_per_nona )


## remove all these variable from the night, and from the days on the early side of the noon-to-noon period (only for those depending on SPT) if the noon-to-noon period has large gaps of missing data (because then we can't reliably calculate the SPT)
sleep_per[ sleep_per$max_time_diff > time_gap & !is.na( sleep_per$max_time_diff ), c( 'onset', 'waking', 'SPT', 'sleep_eff', 'TST', 'WASO', 'wake_bouts', 'summed_VeDBA', 'light_TST', 'light_sleep_eff', 'night_VeDBA_corr', 'prev_naps', 'TST_at_sleep_site', 'frag_wake_bouts', 'frag_wake_bouts_in_median_SPT', 'frag_wake_bouts_at_sleep_site' ) ] <- NA

sleep_per_nona <- sleep_per[ !is.na( sleep_per$SPT ), ]

nrow( sleep_per_nona )

## remove data for sleep period and sleep bouts on days when there is a lot of missing data, because we cannot reliably calculate the sleep VeDBA threshold and there may be a lot of missing epochs
full_dat[ full_dat$n_bursts < ( mins_in_day - missing_mins ), c( 'sleep_per', 'sleep_bouts' ) ] <- NA

## remove data for sleep period on days when there are large gaps of missing data, becasue we can't reliably calculate the SPT with gaps in the data
full_dat[ full_dat$max_time_diff > time_gap, 'sleep_per'  ] <- NA


## make columns for just the time part of the sleep onset and waking timestamps
sleep_per$onset_time <- as_hms( sleep_per$onset )
sleep_per$waking_time <- as_hms( sleep_per$waking )


write.csv( sleep_per, 'DATA/sleep_analysis/sleep_per.csv', row.names = F )
write.csv( full_dat, 'DATA/sleep_analysis/full_dat.csv', row.names = F )

hist( sleep_per$TST, breaks = 100, xlim = c( 0, 900 ) )

for( tag in unique( sleep_per$tag ) ) {
  
  tag_dat <- sleep_per[ sleep_per$tag == tag, ]
  
  hist( tag_dat$TST, breaks = 100, xlim = c( 0, 900 ), main = tag )
}


double_check <- which( sleep_per$TST < 350 )

# for( ni in double_check ){
#
#     sleep_per_func( tag = sleep_per$tag[ ni ], night_num = sleep_per$night_num[ ni ], title = T, plot_waso = T, x_axis = T, las = 1, missing_mins = missing_mins, time_gap = time_gap, mov_window = mov_window, percentile = percentile, multiplier = multiplier, block_size = block_size, gap_size = gap_size, waso_window = waso_window, waso_block = waso_block, waso_percentile = waso_percentile, waso_multiplier = waso_multiplier )
#
# }




################ Statistical analysis ##############

# sleep_per <- read.csv( 'DATA/sleep_analysis/sleep_per.csv' )
#
# full_dat <- read.csv( 'DATA/sleep_analysis/full_dat.csv' )

sleep_per$night <- as.Date( sleep_per$night, tz = 'UTC' )
sleep_per$group <- str_split_fixed( sleep_per$tag, ' ', 3 )[ , 2 ]
sleep_per$onset_time <- as_hms( sleep_per$onset_time )
sleep_per$waking_time <- as_hms( sleep_per$waking_time )

as_hms( as.numeric( mean( sleep_per$onset_time, na.rm = T ) ) )
as.numeric( std.error( sleep_per$onset_time, na.rm = T ) ) / 60

as_hms( as.numeric( mean( sleep_per$waking_time, na.rm = T ) ) )
as.numeric( std.error( sleep_per$waking_time, na.rm = T ) ) / 60

mean( sleep_per$SPT / 60, na.rm = T )
std.error( sleep_per$SPT / 60, na.rm = T)

mean( sleep_per$TST / 60, na.rm = T )
std.error( sleep_per$TST / 60 , na.rm = T)

mean( sleep_per$sleep_eff, na.rm = T )
std.error( sleep_per$sleep_eff, na.rm = T)

#### Modeling the effects of cosleeping ####

head( bdat )

# produced a dataframe with one row for each day for each group which says which sleep site the group slept in on that night
ss_dat <- unique( bdat[ , c( 'group', 'day', 'sleep_clus' ) ] )

# remove the rows for which we don't have data on where the group slept
ss_dat <- ss_dat[ !is.na( ss_dat$sleep_clus ), ]

ss_dat$day <- as.Date( ss_dat$day, tz = 'UTC' )

ss_dat[ duplicated( ss_dat[ , c( 'group', 'day' ) ] ), ] # each group has a single sleep site that we can assign it to for each night

# make a column that will declare whether the group coslept with another group that night. Instantiate the column with all 0s
ss_dat$cosleep <- 0

# fill in the cosleep columns with 1's on nights when the group slept at the same site as another group
ss_dat[ duplicated( ss_dat[ , c( 'day', 'sleep_clus' ) ] ) | duplicated( ss_dat[ , c( 'day', 'sleep_clus' ) ], fromLast = T ), 'cosleep' ] <- 1

# make a column in the sleep_per dataframe for group (this is a duplicate of above, but shouldn't cause any problems)
sleep_per$group <- str_split_fixed( sleep_per$tag, " ", 3 )[ , 2 ]

# merge the dataframe
sleep_per_cosleep <- merge( x = sleep_per, y = ss_dat, by.x = c( 'group', 'night' ), by.y = c( 'group', 'day' ), all.x = T, all.y = T, sort = F )

sleep_per_cosleep$cosleep <- as.factor( sleep_per_cosleep$cosleep )

sleep_per_cosleep$sleep_clus <- as.factor( sleep_per_cosleep$sleep_clus )

# save the last date for which the first baboon group runs out of data run. We will create some models limited to nights before this group's last night of data, because we don't know for sure whether groups coslept after that or not (they may have coslept with the group who no longer has tracking data)
max_night <- min( aggregate( bdat$day, by = list( bdat$group ), FUN = max )[ , 2 ] )

library( boot )


#### prep the data for modeling ######

sleep_per_cosleep$onset_time_num <- as.numeric( as_hms( as.character( sleep_per_cosleep$onset_time ) ) )

hist( sleep_per_cosleep$onset_time_num, xaxt = 'n', breaks = 40 )
axis( 1, at = seq( 0, 60*24*60, 60*60), labels = as_hms( seq( 0, 60*24*60, 60*60) ) )

min( sleep_per_cosleep$onset_time_num, na.rm = T )

## we need to make the morning onset times later than the evening onset times, so we will add the numeric equivalent to 24 hours to them

noon_num <- as.numeric( as_hms( as.character( '12:00:00' ) ) )

sleep_per_cosleep$onset_time_num[ which( sleep_per_cosleep$onset_time_num < noon_num ) ] <- sleep_per_cosleep$onset_time_num[ which( sleep_per_cosleep$onset_time_num < noon_num ) ] + 24*60*60

## checking that it is all better now
hist( sleep_per_cosleep$onset_time_num, xaxt = 'n', breaks = 40 )
axis( 1, at = seq( 0, 60*24*60, 60*60), labels = as_hms( seq( 0, 60*24*60, 60*60) ) )

abline( v = median( sleep_per_cosleep$onset_time_num, na.rm = T ) )

# find median sleep onset time
as_hms( median( sleep_per_cosleep$onset_time_num, na.rm = T ) )


sleep_per_cosleep$waking_time_num <- as.numeric( as_hms( as.character( sleep_per_cosleep$waking_time ) ) )

hist( sleep_per_cosleep$waking_time_num, xaxt = 'n', breaks = 40 )
axis( 1, at = seq( 0, 60*24*60, 60*60), labels = as_hms( seq( 0, 60*24*60, 60*60) ) )

## we need to make the night-time waking times earlier than the morning waking times, so we will subtract the numeric equivalent of 24 hours from them

sleep_per_cosleep$waking_time_num[ which( sleep_per_cosleep$waking_time_num > noon_num ) ] <- sleep_per_cosleep$waking_time_num[ which( sleep_per_cosleep$waking_time_num > noon_num ) ] - 24*60*60

## checking that it is all better now
hist( sleep_per_cosleep$waking_time_num, xaxt = 'n', breaks = 40 )
axis( 1, at = seq( 0, 60*24*60, 60*60), labels = as_hms( seq( 0, 60*24*60, 60*60) ) )

abline( v = median( sleep_per_cosleep$waking_time_num, na.rm = T ) )

as_hms( median( sleep_per_cosleep$waking_time_num, na.rm = T ) )

sleep_per_cosleep$onset_time_num_std <- normalize_func( sleep_per_cosleep$onset_time_num )
sleep_per_cosleep$waking_time_num_std <- normalize_func( sleep_per_cosleep$waking_time_num )


sleep_per_cosleep$TST_std <- normalize_func( sleep_per_cosleep$TST )

sleep_per_cosleep$frag_wake_bouts[ is.na( sleep_per_cosleep$TST ) ] <- NA

sleep_per_cosleep$fragG2 <- sleep_per_cosleep$frag_wake_bouts / ( sleep_per_cosleep$TST / 60 )

sleep_per_cosleep$fragG2_std <- normalize_func( sleep_per_cosleep$fragG2 )

sleep_per_cosleep$SPT_std <- normalize_func( sleep_per_cosleep$SPT )

sleep_per_cosleep$ave_vedba_std <- normalize_func( sleep_per_cosleep$ave_vedba )

sleep_per_cosleep$sleep_eff_std <- normalize_func( sleep_per_cosleep$sleep_eff )

sleep_per_cosleep$logit_sleep_eff <- logit( sleep_per_cosleep$sleep_eff )

sleep_per_cosleep$logit_sleep_eff_std <- normalize_func( sleep_per_cosleep$logit_sleep_eff )

sleep_per_cosleep$TST_at_sleep_site_std <- normalize_func( sleep_per_cosleep$TST_at_sleep_site )

sleep_per_cosleep$tag <- as.factor( sleep_per_cosleep$tag )
sleep_per_cosleep$group <- as.factor( sleep_per_cosleep$group )
sleep_per_cosleep$night <- as.factor( sleep_per_cosleep$night )
sleep_per_cosleep$cosleep <- as.factor( sleep_per_cosleep$cosleep )
sleep_per_cosleep$sleep_clus <- as.factor( sleep_per_cosleep$sleep_clus )


# need to remove rows of incomplete data for the brm function below
sleep_sub <- sleep_per_cosleep[ !is.na( sleep_per_cosleep$TST ) & !is.na( sleep_per_cosleep$cosleep ) & !is.na( sleep_per_cosleep$sleep_clus ), ]


write.csv( sleep_per_cosleep, "DATA/sleep_per_cosleep.csv", row.names = F )

sleep_per_cosleep <- read.csv( "DATA/sleep_per_cosleep.csv" )



## how long to do baboons sleep at the sleep site?
median( sleep_per_cosleep$TST_at_sleep_site / 60 , na.rm = T )

sd( sleep_per_cosleep$TST_at_sleep_site / 60 , na.rm = T ) / sqrt( sum( !is.na( sleep_per_cosleep$TST_at_sleep_site / 60 ) ) )


## when does sleep onset occur?
hist( sleep_per_cosleep$onset_time_num, xaxt = 'n' )
axis( 1, at = seq( 0, 60*24*60, 60*60), labels = as_hms( seq( 0, 60*24*60, 60*60) ) )
as_hms( median( sleep_per_cosleep$onset_time_num, na.rm = T ) )
sd( sleep_per_cosleep$onset_time_num, na.rm = T )/sqrt( sum( !is.na( sleep_per_cosleep$onset_time_num ) ) ) / 60

## when does waking occur
hist( sleep_per_cosleep$waking_time_num, xaxt = 'n' )
axis( 1, at = seq( 0, 60*24*60, 60*60), labels = as_hms( seq( 0, 60*24*60, 60*60) ) )
as_hms( median( sleep_per_cosleep$waking_time_num, na.rm = T ) )
sd( sleep_per_cosleep$waking_time_num, na.rm = T )/sqrt( sum( !is.na( sleep_per_cosleep$waking_time_num ) ) ) / 60

## what was the efficiency within the sleep period?
median( sleep_per_cosleep$sleep_eff, na.rm = T )
sd( sleep_per_cosleep$sleep_eff, na.rm = T ) / sqrt( sum( !is.na( sleep_per_cosleep$sleep_eff ) ) )

median( sleep_per_cosleep$fragG2, na.rm = T )
sd( sleep_per_cosleep$fragG2, na.rm = T ) / sqrt( sum( !is.na( sleep_per_cosleep$fragG2 ) ) )



##### model the effect of cosleeping on total sleep time 

# plot the histogram of total sleep time
hist( sleep_per_cosleep$TST_at_sleep_site )

# need to remove rows of incomplete data 
sleep_sub <- sleep_per_cosleep[ !is.na( sleep_per_cosleep$TST ) & !is.na( sleep_per_cosleep$cosleep ), ]

class( sleep_sub$night )

### Bayesian model ###



library( brms )

options( mc.cores = parallel::detectCores() )

priors <- c( prior(normal(0,2), class = "b", coef = cosleep1 ),
             prior(student_t(3, 0, 2.5), class = "Intercept" ) )


options( mc.cores = parallel::detectCores() )

# with group, tag, and night as random effects
#TST_at_sleep_site_mod <- brm( TST_at_sleep_site_std ~ cosleep + ( 1 | sleep_clus ) + ( 1 | group ) + ( 1 | group:tag ) + ( 1 | night ), data = sleep_sub, family= "student", iter = 5000, prior = priors, control = list(max_treedepth = 15, adapt_delta = .999999999999 ), chains = 4, save_pars = save_pars( all = T ) )


#saveRDS( object = TST_at_sleep_site_mod, file = "RESULTS/models/TST_at_sleep_site_student_mod_std.rds" )

TST_at_sleep_site_student_mod <- readRDS( file = "RESULTS/models/TST_at_sleep_site_student_mod_std.rds" )

summary( TST_at_sleep_site_student_mod, prob = 0.9 )

#tab_model( TST_at_sleep_site_student_mod, show.ci = 0.9 )

pp_check( TST_at_sleep_site_student_mod )

library( ggplot2 )

plot( conditional_effects( TST_at_sleep_site_student_mod, effects = c( 'cosleep' ), spaghetti = F, re_formula = NULL ),  line_args = list( colour = "white", size = 2 ), errorbar_args = list( colour = "white" ), cat_args = list( colour = "white" ), theme = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_rect(fill = "black", colour = "black"), panel.background = element_rect(fill = "black", colour = "black"), axis.line = element_line(colour = "white"), axis.text = element_text(size = rel(0.8), colour = "white"), axis.ticks = element_line(colour = "white")  ) )

sleep_sub$cosleep_jitter <- jitter( as.numeric( sleep_sub$cosleep ) )

plot( conditional_effects( TST_at_sleep_site_student_mod, effects = c( 'cosleep' ), spaghetti = F ), line_args = list( colour = "blue", size = 2 ), errorbar_args = list( colour = "blue", size = 2 ), cat_args = list( colour = "black" ), theme = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) )[[ 1 ]] + geom_point( data = sleep_sub, aes( x = cosleep_jitter, y = TST ), alpha = 0.1, inherit.aes = F )



posterior <- posterior_samples( TST_at_sleep_site_student_mod )

mean( -posterior$b_cosleep1*sd( sleep_sub$TST_at_sleep_site ) )

-0.52*sd( sleep_sub$TST_at_sleep_site )

coef_dens <- density( -posterior$b_cosleep1*sd( sleep_sub$TST_at_sleep_site ) )

transp <- function(col, alpha=.5){
  res <- apply(col2rgb(col),2, function(c) rgb(c[1]/255, c[2]/255, c[3]/255, alpha))
  return(res)
}


plot( coef_dens$x, coef_dens$y, type = 'l', xlab = 'Minutes of sleep lost when sharing the sleep site with another group', bty = 'l', las = 1, ylab = 'Posterior probability density', col = 'white', col.axis = 'white', col.lab = 'white', col.main = 'white' )
axis( 1, labels = F, col = 'white' )
axis( 2, labels = F, col = 'white' )


rel_quantiles <- quantile( -posterior$b_cosleep1*sd( sleep_sub$TST_at_sleep_site ), c( 0.05, 0.95 ) )

abline( v = rel_quantiles, lty = 2, col = 'white' )

poly_x <- c( rel_quantiles[ 1 ], coef_dens$x[ coef_dens$x > rel_quantiles[ 1 ] & coef_dens$x < rel_quantiles[ 2 ] ], rel_quantiles[ 2 ] )
poly_y <- c( -0.01, coef_dens$y[ coef_dens$x > rel_quantiles[ 1 ] & coef_dens$x < rel_quantiles[ 2 ] ], -0.01 )

polygon( x = poly_x, y = poly_y, col = transp( 'white', 0.3), border = F )



## white for non-presentation

posterior <- posterior_samples( TST_at_sleep_site_student_mod )

mean( -posterior$b_cosleep1*sd( sleep_sub$TST_at_sleep_site ) )

-0.52*sd( sleep_sub$TST_at_sleep_site )

coef_dens <- density( -posterior$b_cosleep1*sd( sleep_sub$TST_at_sleep_site ) )

transp <- function(col, alpha=.5){
  res <- apply(col2rgb(col),2, function(c) rgb(c[1]/255, c[2]/255, c[3]/255, alpha))
  return(res)
}


plot( coef_dens$x, coef_dens$y, type = 'l', xlab = 'Sleep deficit when sharing the sleep site (min)', bty = 'l', las = 1, ylab = 'Posterior probability density' )

rel_quantiles <- quantile( -posterior$b_cosleep1*sd( sleep_sub$TST_at_sleep_site ), c( 0.05, 0.95 ) )

abline( v = rel_quantiles, lty = 2 )

poly_x <- c( rel_quantiles[ 1 ], coef_dens$x[ coef_dens$x > rel_quantiles[ 1 ] & coef_dens$x < rel_quantiles[ 2 ] ], rel_quantiles[ 2 ] )
poly_y <- c( -0.01, coef_dens$y[ coef_dens$x > rel_quantiles[ 1 ] & coef_dens$x < rel_quantiles[ 2 ] ], -0.01 )

polygon( x = poly_x, y = poly_y, col = transp( 'grey', 0.3), border = F )







## sleep efficiency model

options( mc.cores = parallel::detectCores() )

priors <- c( prior(normal(0,2), class = "b", coef = cosleep1 ),
             prior(student_t(3, 0, 2.5), class = "Intercept" ) )

sleep_sub$logit_sleep_eff <- log( sleep_sub$sleep_eff / ( 1 - sleep_sub$sleep_eff ) )

sleep_sub$logit_sleep_eff_std <- normalize_func( sleep_sub$logit_sleep_eff )

class( sleep_sub$sleep_clus )
class( sleep_sub$cosleep )

# with group and tag as nested random effects
#sleep_eff_mod_std <- brm( logit_sleep_eff_std ~ cosleep + ( 1 | sleep_clus ) + ( 1 | group ) + ( 1 | group:tag )  + ( 1 | night ), data = sleep_sub, family= "gaussian", iter = 5000, prior = priors, control = list(max_treedepth = 15, adapt_delta = .999999999999 ), chains = 4  )


# saveRDS( object = sleep_eff_mod_std, file = "RESULTS/models/sleep_eff_mod_std_sleep_clus_re.rds" )

sleep_eff_mod_std <- readRDS( file = "RESULTS/models/sleep_eff_mod_std_sleep_clus_re.rds" )

summary( sleep_eff_mod_std )

# tab_model( sleep_eff_mod_std, show.ci = 0.90 )

pp_check( sleep_eff_mod_std )

effects <- fixef( sleep_eff_mod_std )[ -1, ]

par( bg = 'black' )

plot( 0, type = 'n', xlim = c( -3, max( effects )), ylim = c( 0, nrow( effects  ) ), bty = 'n', xlab = '', ylab = '', yaxt = 'n', xaxt = 'n', col.axis = 'white', col.lab = 'white', col.main = 'white', col = 'white'  )

axis( 1, seq( -1.5, 1, by = .5 ), col = 'white', col.ticks = 'white', col.axis = 'white')
abline( v = 0, lty = 2, col = 'white' )

for( i in 1:nrow( effects ) ){
  
  points( x = effects[ i, 1 ], y = i, pch = 16, col = 'white' )
  
  segments( x0 = effects[ i, 3 ], x1 = effects[ i, 4 ], y0 = i, y1 = i, col = 'white'  )
  
  text( x = -1.8, y = i, labels = rownames( effects )[ i ], col = 'white'  )
  
}



plot( conditional_effects( sleep_eff_mod_std ) )


posterior <- posterior_samples( sleep_eff_mod_std )


# intercept is 
invlogit( mean( sleep_sub$logit_sleep_eff ) + 0.03*sd( sleep_sub$logit_sleep_eff ) )

invlogit( ( mean( sleep_sub$logit_sleep_eff ) + 0.03*sd( sleep_sub$logit_sleep_eff ) ) - 0.20*sd( sleep_sub$logit_sleep_eff ) )

invlogit( mean( sleep_sub$logit_sleep_eff ) + 0.03*sd( sleep_sub$logit_sleep_eff ) ) - invlogit( ( mean( sleep_sub$logit_sleep_eff ) + 0.03*sd( sleep_sub$logit_sleep_eff ) ) - 0.20*sd( sleep_sub$logit_sleep_eff ) )


changes <- - apply( posterior[ , c( 'b_Intercept', 'b_cosleep1' ) ], 1, function( x )  ( ( invlogit( mean( sleep_sub$logit_sleep_eff ) + x[ 1 ]*sd( sleep_sub$logit_sleep_eff ) ) - ( invlogit( ( mean( sleep_sub$logit_sleep_eff ) + x[ 1 ]*sd( sleep_sub$logit_sleep_eff ) ) + x[ 2 ]*sd( sleep_sub$logit_sleep_eff ) ) ) ) ) ) *100



#changes <- apply( posterior[ , c( 'b_Intercept', 'b_cosleep1' ) ], 1, function( x )  ( ( invlogit( mean( sleep_sub$logit_sleep_eff ) ) - ( invlogit( ( mean( sleep_sub$logit_sleep_eff ) + x[ 2 ]*sd( sleep_sub$logit_sleep_eff ) ) ) ) ) ) )

dev.off()

coef_dens <- density( changes )

transp <- function(col, alpha=.5){
  res <- apply(col2rgb(col),2, function(c) rgb(c[1]/255, c[2]/255, c[3]/255, alpha))
  return(res)
}


plot( coef_dens$x, coef_dens$y, type = 'l', xlab = 'Change in sleep efficiency when sharing the sleep site (%)', bty = 'l', las = 1, ylab = 'Posterior probability density' )

rel_quantiles <- quantile( changes, c( 0.05, 0.95 ) )

abline( v = rel_quantiles, lty = 2 )

poly_x <- c( rel_quantiles[ 1 ], coef_dens$x[ coef_dens$x > rel_quantiles[ 1 ] & coef_dens$x < rel_quantiles[ 2 ] ], rel_quantiles[ 2 ] )
poly_y <- c( -.03, coef_dens$y[ coef_dens$x > rel_quantiles[ 1 ] & coef_dens$x < rel_quantiles[ 2 ] ], -.03 )

polygon( x = poly_x, y = poly_y, col = transp( 'grey', 0.3), border = F )


par( bg = 'black' )

plot( coef_dens$x, coef_dens$y, type = 'l', xlab = 'Sleep efficieny lost when sharing the sleep site with another group (%)', bty = 'l', las = 1, ylab = 'Posterior probability density', col = 'white', col.axis = 'white', col.lab = 'white', col.main = 'white' )
axis( 1, labels = F, col = 'white' )
axis( 2, labels = F, col = 'white' )


rel_quantiles <- quantile( changes, c( 0.05, 0.95 ) )

abline( v = rel_quantiles, lty = 2 , col = 'white' )

poly_x <- c( rel_quantiles[ 1 ], coef_dens$x[ coef_dens$x > rel_quantiles[ 1 ] & coef_dens$x < rel_quantiles[ 2 ] ], rel_quantiles[ 2 ] )
poly_y <- c( -.03, coef_dens$y[ coef_dens$x > rel_quantiles[ 1 ] & coef_dens$x < rel_quantiles[ 2 ] ], -.03 )

polygon( x = poly_x, y = poly_y, col = transp( 'white', 0.3), border = F )



## sleep fragmentation model

priors <- c( prior(normal(0,2), class = "b", coef = cosleep1 ),
             prior(student_t(3, 0, 2.5), class = "Intercept" ) )


hist( sleep_per_cosleep$fragG2 )

#fragG2_mod_std <- brm( fragG2_std ~ cosleep + ( 1 | sleep_clus ) + ( 1 | group ) + ( 1 | group:tag ) + ( 1 | night ), data = sleep_sub, family= "gaussian", iter = 5000, prior = priors, control = list(max_treedepth = 15, adapt_delta = .999999999999 ), chains = 4  )
# 
#saveRDS( object = fragG2_mod_std, file = "RESULTS/models/fragG2_mod_std_sleep_clus_re.rds" )


fragG2_mod_std <- readRDS( file = "RESULTS/models/fragG2_mod_std_sleep_clus_re.rds" )

summary( fragG2_mod_std )

# tab_model( fragG2_mod_std, show.ci = 0.90 )

pp_check( fragG2_mod_std )



posterior <- posterior_samples( fragG2_mod_std )

mean( posterior$b_cosleep1*sd( sleep_sub$fragG2 ) )

#0.2*sd( sleep_sub$fragG2 )

coef_dens <- density( posterior$b_cosleep1*sd( sleep_sub$fragG2 ) )

transp <- function(col, alpha=.5){
  res <- apply(col2rgb(col),2, function(c) rgb(c[1]/255, c[2]/255, c[3]/255, alpha))
  return(res)
}


plot( coef_dens$x, coef_dens$y, type = 'l', xlab = 'Change in sleep fragmentation when sharing the sleep site\n(awakenings/hour)', bty = 'l', las = 1, ylab = 'Posterior probability density', col = 'white', col.axis = 'white', col.lab = 'white', col.main = 'white' )
axis( 1, labels = F, col = 'white' )
axis( 2, labels = F, col = 'white' )

rel_quantiles <- quantile( posterior$b_cosleep1*sd( sleep_sub$fragG2 ), c( 0.05, 0.95 ) )

abline( v = rel_quantiles, lty = 2, col = 'white' )

poly_x <- c( rel_quantiles[ 1 ], coef_dens$x[ coef_dens$x > rel_quantiles[ 1 ] & coef_dens$x < rel_quantiles[ 2 ] ], rel_quantiles[ 2 ] )
poly_y <- c( -0.4, coef_dens$y[ coef_dens$x > rel_quantiles[ 1 ] & coef_dens$x < rel_quantiles[ 2 ] ], -0.4 )

polygon( x = poly_x, y = poly_y, col = transp( 'white', 0.3), border = F )






plot( coef_dens$x, coef_dens$y, type = 'l', xlab = 'Change in sleep fragmentation when sharing the sleep site\n(awakenings/hour)', bty = 'l', las = 1, ylab = 'Posterior probability density', col = 'black', col.axis = 'black', col.lab = 'black', col.main = 'black' )

rel_quantiles <- quantile( posterior$b_cosleep1*sd( sleep_sub$fragG2 ), c( 0.05, 0.95 ) )

abline( v = rel_quantiles, lty = 2, col = 'black' )

poly_x <- c( rel_quantiles[ 1 ], coef_dens$x[ coef_dens$x > rel_quantiles[ 1 ] & coef_dens$x < rel_quantiles[ 2 ] ], rel_quantiles[ 2 ] )
poly_y <- c( -0.4, coef_dens$y[ coef_dens$x > rel_quantiles[ 1 ] & coef_dens$x < rel_quantiles[ 2 ] ], -0.4 )

polygon( x = poly_x, y = poly_y, col = transp( 'grey', 0.3), border = F )






effects <- fixef( fragG2_mod_std )[ -1,  ]

par( bg = 'black' )

plot( 0, type = 'n', xlim = c( -3, max( effects )), ylim = c( 0, nrow( effects  ) ), bty = 'n', xlab = '', ylab = '', yaxt = 'n', xaxt = 'n', col.axis = 'white', col.lab = 'white', col.main = 'white', col = 'white'  )

axis( 1, seq( -1.5, 1, by = .5 ), col = 'white', col.ticks = 'white', col.axis = 'white')
abline( v = 0, lty = 2, col = 'white' )

for( i in 1:nrow( effects ) ){
  
  points( x = effects[ i, 1 ], y = i, pch = 16, col = 'white' )
  
  segments( x0 = effects[ i, 3 ], x1 = effects[ i, 4 ], y0 = i, y1 = i, col = 'white'  )
  
  text( x = -1.8, y = i, labels = rownames( effects )[ i ], col = 'white'  )
  
}


plot( conditional_effects( fragG2_mod_std ) )








