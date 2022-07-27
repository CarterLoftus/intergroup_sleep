



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
library( boot )

##### The effect of sleeping together on the probability of moving together the next day ######

## function for normalizing a vector
normalize_func <- function( x ) return( (x - mean( x, na.rm = T ) )/ sd( x, na.rm = T ) )

na_min <- function( vec ){
  
  if( sum( !is.na( vec ) ) == 0 ){
    
    return( NA )
    
  }else{
    
    return( min( vec, na.rm = T ) )
  }
}

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


## fix threshold is the proportion of successful fixes that must have occurred from 6 PM to 6 AM and from 6 AM to 6 PM to say that we know whether the group coslept the night before or comoved with another group during the following day, respectively
cosleep_fix_threshold <- 0.5
comove_fix_threshold <- 0.85


# remove the rows of spec_df that do not have location data
spec_df <- spec_df[ !is.na( spec_df$x ), ]

# make the timestamp into a POSIX element
spec_df$local_timestamp <- as.POSIXct( spec_df$local_timestamp, tz = 'UTC' )

# turn the day column into a date rather than day from the start of the study period (makes it easier to merge information across dataframes and between day and night dataframes)
spec_df$day <- as.Date( spec_df$local_timestamp, tz = 'UTC' )


### read in the interaction dataframe
dyad_interact_df <- read.csv( "DATA/dyad_interact_df.csv" )

# change the date to a real date element
dyad_interact_df$day <- as.Date( dyad_interact_df$day, tz = 'UTC' )

dyad_interact_df$interaction[ dyad_interact_df$interaction == 'cosleep_then_comove' ] <- 'cosleep'


fixes_in_12_hours <- length( seq( as.numeric( as_hms( "00:00:00" ) ), as.numeric( as_hms( "12:00:00" ) ), by =  15*60 ) )


# make a dataframe containing a row for each group on each day of data it has
temp_info_df <- unique( spec_df[ , c( 'group', 'day' ) ] ) 

# save the names of the groups
group_names <- as.character( unique( temp_info_df$group ) )

# repeat this dataframe by the number of groups that there are, because this is going to be a dyadic dataframe (one row for each dyad on each day on which the focal group has data)
soc_info_df <- temp_info_df[ rep( seq_len( nrow( temp_info_df ) ), each = length( group_names ) ), ]

## order the dataframe by the timestamp
soc_info_df <- soc_info_df[ order( soc_info_df$group) , ]
soc_info_df <- soc_info_df[ order( soc_info_df$day) , ]

## change the column names of id and group to be id1 and group1. This will be our focal individual
names( soc_info_df )[ names( soc_info_df ) == 'group' ] <- 'group1'

## create a column for the names of the non-focal dyad member
soc_info_df$group2 <- group_names

# remove the rows where both dyad members are the same group
soc_info_df <- soc_info_df[ soc_info_df$group1 != soc_info_df$group2, ]

soc_info_df$comove_today <- NA

## fill in the columns indicating whether the groups are comoving on the current day. Don't fill out this information for days on which at least one of the dyad members does not have sufficient data -- instead, put an NA in the column on those days
for( i in 1:nrow( soc_info_df ) ){
  
  day_after_group1 <- spec_df[ spec_df$group == soc_info_df$group1[ i ] & spec_df$local_timestamp >= as.POSIXct( paste( as.Date( soc_info_df$day[ i ], origin = '1970-01-01' ), '06:00:00' ), tz = 'UTC' ) & spec_df$local_timestamp <= as.POSIXct( paste( as.Date( soc_info_df$day[ i ], origin = '1970-01-01' ), '18:00:00' ), tz = 'UTC' ), ]
  
  day_after_group1 <- day_after_group1[ !duplicated( day_after_group1$local_timestamp ), ]
  
  day_after_group2 <- spec_df[ spec_df$group == soc_info_df$group2[ i ] & spec_df$local_timestamp >= as.POSIXct( paste( as.Date( soc_info_df$day[ i ], origin = '1970-01-01' ), '06:00:00' ), tz = 'UTC' ) & spec_df$local_timestamp <= as.POSIXct( paste( as.Date( soc_info_df$day[ i ], origin = '1970-01-01' ), '18:00:00' ), tz = 'UTC' ), ]
  
  day_after_group2 <- day_after_group2[ !duplicated( day_after_group2$local_timestamp ), ]
  
  if( ( nrow( day_after_group1 ) > comove_fix_threshold*fixes_in_12_hours ) & ( nrow( day_after_group2 ) > comove_fix_threshold*fixes_in_12_hours ) ){
    
    today_interaction <- dyad_interact_df[ dyad_interact_df$group1 == soc_info_df$group1[ i ] & dyad_interact_df$group2 == soc_info_df$group2[ i ] & dyad_interact_df$day == as.Date( soc_info_df$day[ i ], origin = '1970-01-01' ), 'interaction' ]
    
    soc_info_df[ i, 'comove_today'] <- ifelse( grepl( 'comove', today_interaction ), 1, 0 )
    
    soc_info_df[ i , 'full_interaction_today' ] <- today_interaction
    
  }
}

soc_info_df$comove_today <- as.factor( soc_info_df$comove_today )

### check coslept last night ###

soc_info_df$auto_cosleep <- NA

for( i in 1:nrow( soc_info_df ) ){
  
  sleep_site_1 <- unique( spec_df[ spec_df$group == soc_info_df$group1[ i ] & spec_df$day == as.Date( soc_info_df$day[ i ] - 1, origin = '1970-01-01' ), 'sleep_clus' ] )
  
  sleep_site_1 <- sleep_site_1[ !is.na( sleep_site_1 ) ]
  
  sleep_site_2 <- unique( spec_df[ spec_df$group == soc_info_df$group2[ i ] & spec_df$day == as.Date( soc_info_df$day[ i ] - 1, origin = '1970-01-01' ), 'sleep_clus' ] )
  
  sleep_site_2 <- sleep_site_2[ !is.na( sleep_site_2 ) ]
  
  if( length( sleep_site_1 ) > 1 | length( sleep_site_2 ) > 1 ) stop( 'more than one sleep site for this day' )
  
  ## if we know where both groups slept on the previous night
  if( length( sleep_site_1 ) == 1 & length( sleep_site_2 ) == 1  ){
    
    # add whether they slept at the same site or not as a binary variable to the dataframe
    soc_info_df$auto_cosleep[ i ] <- as.numeric( sleep_site_1 == sleep_site_2 )
    
  }
}



########## models ###########

soc_info_df$any_interaction <- ifelse( soc_info_df$full_interaction_today == 'None', 0, ifelse( !is.na( soc_info_df$full_interaction_today ), 1, NA ) )

soc_info_df$dy_name <- apply( soc_info_df[ , c( 'group1', 'group2' ) ], 1, function( x ) paste( sort( x ), collapse = '_' ) )

soc_info_df$auto_cosleep <- as.factor( as.character( soc_info_df$auto_cosleep ) )


soc_sub <- soc_info_df[ !duplicated( soc_info_df[ , c( 'dy_name', 'day' ) ] ), ]


## test the effect of cosleeping on having an interaction the following day

soc_subber <- soc_sub[ !is.na( soc_sub$auto_cosleep ) & !is.na( soc_sub$any_interaction ), ]

options( mc.cores = parallel::detectCores() )

priors <- c( prior(normal(0,2), class = "b", coef = auto_cosleep1 ),
             prior(student_t(3, 0, 2.5), class = "Intercept" ) )


interaction_vs_cosleep_mod <- brm( any_interaction ~ auto_cosleep + ( 1| group1 ) + ( 1| group2 ) + ( 1 | dy_name ), data = soc_subber, family= "bernoulli", iter = 5000, prior = priors, control = list(max_treedepth = 15, adapt_delta = .999999999999 ) )

dir.create( paste0( getwd( ), '/RESULTS/models' ), recursive = T )

saveRDS( interaction_vs_cosleep_mod, 'RESULTS/models/interaction_vs_cosleep_mod.rds' )

interaction_vs_cosleep_mod <- readRDS( 'RESULTS/models/interaction_vs_cosleep_mod.rds' )

summary( interaction_vs_cosleep_mod, prob = 0.9 )

# tab_model( interaction_vs_cosleep_mod, show.ci = 0.90, transform = NULL )

pp_check( interaction_vs_cosleep_mod )

## plotting model predictions
cosleep_inds <- which( soc_subber$auto_cosleep == 1 )
non_cosleep_inds <- which( soc_subber$auto_cosleep == 0 )


post_pred <- posterior_linpred( interaction_vs_cosleep_mod, transform = T, newdata = soc_subber, re_formula = NULL )

cosleep_vec <- as.numeric( post_pred[  , cosleep_inds ] )

non_cosleep_vec <- as.numeric( post_pred[  , non_cosleep_inds ] )


plot( x = c( 0, 1 ), y = c( mean( non_cosleep_vec ), mean( cosleep_vec ) ), xlim = c( -0.5, 1.5 ), ylim = c( 0, 1 ), pch = 16, cex = 2, las = 1, bty = 'l', xaxt = 'n', xlab = '', ylab = 'Probability of engaging in an interaction with another group', cex.lab = 1 )

axis(1, at = c( 0, 1 ), labels = c( '', '' ) )
axis(1, at = c( 0, 1 ), labels = c('After sleeping as\na lone group' , 'After sharing sleep site\nwith the other group'), line = 1, tick = F )


CIs <- quantile( non_cosleep_vec, c( 0.05, 0.95 ) )
segments( 0, CIs[ 1 ], 0, CIs[ 2 ] )

segments( -0.2, CIs[ 1 ], 0.2, CIs[ 1 ] )
segments( -0.2, CIs[ 2 ], 0.2, CIs[ 2 ] )


CIs <- quantile( cosleep_vec, c( 0.05, 0.95 ) )
segments( 1, CIs[ 1 ], 1, CIs[ 2 ] )

segments( 0.8, CIs[ 1 ], 1.2, CIs[ 1 ] )
segments(  0.8, CIs[ 2 ], 1.2, CIs[ 2 ] )



### black for presentation

cosleep_inds <- which( soc_subber$auto_cosleep == 1 )
non_cosleep_inds <- which( soc_subber$auto_cosleep == 0 )


post_pred <- posterior_linpred( interaction_vs_cosleep_mod, transform = T, newdata = soc_subber, re_formula = NULL )

cosleep_vec <- as.numeric( post_pred[  , cosleep_inds ] )

non_cosleep_vec <- as.numeric( post_pred[  , non_cosleep_inds ] )


plot( x = c( 0, 1 ), y = c( mean( non_cosleep_vec ), mean( cosleep_vec ) ), xlim = c( -0.5, 1.5 ), ylim = c( 0, 1 ), pch = 16, cex = 2, las = 1, bty = 'l', xaxt = 'n', xlab = '', ylab = 'Probability of engaging in an interaction with another group', cex.lab = 1, col = 'white', col.axis = 'white', col.lab = 'white', col.main = 'white', )

axis(2, tick = T, labels = F,  col = 'white' )

axis(1, at = c( 0, 1 ), labels = c( '', '' ), col = 'white' )
axis(1, at = c( 0, 1 ), labels = c('After sleeping as\na lone group' , 'After sharing sleep site\nwith the other group'), line = 1, tick = F, col = 'white' )


CIs <- quantile( non_cosleep_vec, c( 0.05, 0.95 ) )
segments( 0, CIs[ 1 ], 0, CIs[ 2 ], col = 'white' )

segments( -0.2, CIs[ 1 ], 0.2, CIs[ 1 ], col = 'white' )
segments( -0.2, CIs[ 2 ], 0.2, CIs[ 2 ], col = 'white' )


CIs <- quantile( cosleep_vec, c( 0.05, 0.95 ) )
segments( 1, CIs[ 1 ], 1, CIs[ 2 ], col = 'white' )

segments( 0.8, CIs[ 1 ], 1.2, CIs[ 1 ], col = 'white' )
segments(  0.8, CIs[ 2 ], 1.2, CIs[ 2 ], col = 'white' )



## test the effect of cosleeping on an interaction involving comoving the following day
soc_subber <- soc_sub[ !is.na( soc_sub$auto_cosleep ) & !is.na( soc_sub$comove_today ) & soc_sub$full_interaction_today != 'None' & !is.na( soc_sub$full_interaction_today ), ]


priors <- c( prior(normal(0,2), class = "b", coef = auto_cosleep1 ),
             prior(student_t(3, 0, 2.5), class = "Intercept" ) )

comove_vs_cosleep_mod_condition_on_interaction <- brm( comove_today ~ auto_cosleep + ( 1| group1 ) + ( 1| group2 ) + ( 1 | dy_name ), data = soc_subber, family= "bernoulli", iter = 5000, prior = priors, control = list(max_treedepth = 15, adapt_delta = .999999999999 ) )

saveRDS( comove_vs_cosleep_mod_condition_on_interaction, 'RESULTS/models/comove_vs_cosleep_mod_condition_on_interaction.rds' )

comove_vs_cosleep_mod_condition_on_interaction <- readRDS( 'RESULTS/models/comove_vs_cosleep_mod_condition_on_interaction.rds' )

summary( comove_vs_cosleep_mod_condition_on_interaction, prob = 0.9 )

# tab_model( comove_vs_cosleep_mod_condition_on_interaction, show.ci = 0.90, transform = NULL )

pp_check( comove_vs_cosleep_mod_condition_on_interaction )

## plotting model predictions
cosleep_inds <- which( soc_subber$auto_cosleep == 1 )
non_cosleep_inds <- which( soc_subber$auto_cosleep == 0 )


post_pred <- posterior_linpred( comove_vs_cosleep_mod_condition_on_interaction, transform = T, newdata = soc_subber, re_formula = NULL )

cosleep_vec <- as.numeric( post_pred[  , cosleep_inds ] )

non_cosleep_vec <- as.numeric( post_pred[  , non_cosleep_inds ] )


plot( x = c( 0, 1 ), y = c( mean( non_cosleep_vec ), mean( cosleep_vec ) ), xlim = c( -0.5, 1.5 ), ylim = c( 0, 1 ), pch = 16, cex = 2, las = 1, bty = 'l', xaxt = 'n', xlab = '', ylab = 'Probability of an interaction involving cohesive movement', cex.lab = 1 )

axis(1, at = c( 0, 1 ), labels = c( '', '' ) )
axis(1, at = c( 0, 1 ), labels = c('After sleeping as\na lone group' , 'After sharing sleep site\nwith the interacting group'), line = 1, tick = F )


CIs <- quantile( non_cosleep_vec, c( 0.05, 0.95 ) )
segments( 0, CIs[ 1 ], 0, CIs[ 2 ] )

segments( -0.2, CIs[ 1 ], 0.2, CIs[ 1 ] )
segments( -0.2, CIs[ 2 ], 0.2, CIs[ 2 ] )


CIs <- quantile( cosleep_vec, c( 0.05, 0.95 ) )
segments( 1, CIs[ 1 ], 1, CIs[ 2 ] )

segments( 0.8, CIs[ 1 ], 1.2, CIs[ 1 ] )
segments(  0.8, CIs[ 2 ], 1.2, CIs[ 2 ] )




### black for presentation

cosleep_inds <- which( soc_subber$auto_cosleep == 1 )
non_cosleep_inds <- which( soc_subber$auto_cosleep == 0 )


post_pred <- posterior_linpred( comove_vs_cosleep_mod_condition_on_interaction, transform = T, newdata = soc_subber, re_formula = NULL )

cosleep_vec <- as.numeric( post_pred[  , cosleep_inds ] )

non_cosleep_vec <- as.numeric( post_pred[  , non_cosleep_inds ] )


plot( x = c( 0, 1 ), y = c( mean( non_cosleep_vec ), mean( cosleep_vec ) ), xlim = c( -0.5, 1.5 ), ylim = c( 0, 1 ), pch = 16, cex = 2, las = 1, bty = 'l', xaxt = 'n', xlab = '', ylab = 'Probability of an interaction involving cohesive movement', cex.lab = 1, col = 'white', col.axis = 'white', col.lab = 'white', col.main = 'white' )

axis(2, labels = F, col = 'white' )
axis(1, at = c( 0, 1 ), labels = c( '', '' ), col = 'white' )
axis(1, at = c( 0, 1 ), labels = c('After sleeping as\na lone group' , 'After sharing sleep site\nwith the interacting group'), line = 1, tick = F, col = 'white' )


CIs <- quantile( non_cosleep_vec, c( 0.05, 0.95 ) )
segments( 0, CIs[ 1 ], 0, CIs[ 2 ], col = 'white' )

segments( -0.2, CIs[ 1 ], 0.2, CIs[ 1 ], col = 'white' )
segments( -0.2, CIs[ 2 ], 0.2, CIs[ 2 ], col = 'white' )


CIs <- quantile( cosleep_vec, c( 0.05, 0.95 ) )
segments( 1, CIs[ 1 ], 1, CIs[ 2 ], col = 'white' )

segments( 0.8, CIs[ 1 ], 1.2, CIs[ 1 ], col = 'white' )
segments(  0.8, CIs[ 2 ], 1.2, CIs[ 2 ], col = 'white' )

