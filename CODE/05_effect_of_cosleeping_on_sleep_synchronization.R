

library( data.table )
library( stringr )
library( hms )
library( brms ) 
library( sjPlot )

## function for normalizing a vector
normalize_func <- function( x ) return( (x - mean( x, na.rm = T ) )/ sd( x, na.rm = T ) )


transp <- function(col, alpha=.5){
  res <- apply(col2rgb(col),2, function(c) rgb(c[1]/255, c[2]/255, c[3]/255, alpha))
  return(res)
}

na_mean <- function( vec_1, vec_2 ){
  
  if( sum( !is.na( vec_1 ) ) == 0 | sum( !is.na( vec_2 ) ) == 0 ){
    
    return( NA )
    
  }else{
    
    return( mean( vec_1 == vec_2, na.rm = T ) )
    
  }
}


bdat <- fread( "DATA/bab_complete.csv", fill = T )

## turning bdat into a dataframe from a data table
bdat <- as.data.frame( bdat )

## make the local_timestamp into a POSIXct element
bdat$local_timestamp <- as.POSIXct(bdat$local_timestamp, tz = 'UTC' )

## make a column named day that reprsents the day from the start of the study period, with the first day being day 1
bdat$day <- as.Date( bdat$local_timestamp, tz = 'UTC' )

# check where the day changes from one day to the next
bdat[ ( diff( c( as.Date( bdat$day ) ) ) == 1),]

## make a column named time that will have just the time component, and not the date, of the local_timestamp
bdat$local_time <- str_split_fixed( bdat$local_timestamp, " ",2)[,2]

## remove rows of the dataframe that don't have a fix
bdat <- bdat[!is.na( bdat$lat ),]

## open the first rows of the dataframe
head( bdat )

## limit the data to the time between the median sleep site arrival time and median sleep site departure time
median_sleep_site_arrive_time <- as.character( as_hms( median( as.numeric( as_hms( bdat$arrive_sleep_site ) ), na.rm = T ) ) )

median_sleep_site_leave_time <- as.character( as_hms( median( as.numeric( as_hms( bdat$leave_sleep_site ) ), na.rm = T ) ) ) 


full_dat <- as.data.frame( fread( 'DATA/sleep_analysis/full_dat.csv' ) )

head( full_dat )


night_dat <- full_dat[ full_dat$local_time > median_sleep_site_arrive_time | full_dat$local_time < median_sleep_site_leave_time , ]

tag_names <- as.character( unique( night_dat$tag ) )

vec_a <- c()
vec_b <- c()

for( a in 1:( length( tag_names ) - 1 ) ){
  
  for( b in ( a + 1 ): length( tag_names ) ){
    
    vec_a <- c( vec_a, tag_names[ a ] )
    
    vec_b <- c( vec_b, tag_names[ b ] )
  } 
}

nights <- as.character( unique( night_dat$night ) )

sleep_sync_dat <- data.frame( tag_a = rep( vec_a, times = length( nights ) ), tag_b = rep( vec_b, times = length( nights ) ), night = rep( nights, each = length( vec_a ) ), sync = NA )


night_wide_dat <- as.data.frame( ( data.table::dcast( as.data.table( night_dat ), local_timestamp + night ~ tag, value.var = c( 'sleep_bouts'), crop = F ) ) )


#### find the synchronization score for each group dyad on each night ####
for( i in 1:length( vec_a ) ){
  
  dy_dat <- night_wide_dat[ , c( 'night', vec_a[ i ], vec_b[ i ] ) ]
  
  for( night in nights ){
    
    dy_night_dat <- dy_dat[ dy_dat$night == night, ]
    
    sleep_sync_dat[ sleep_sync_dat$tag_a == vec_a[ i ] & sleep_sync_dat$tag_b == vec_b[ i ] & sleep_sync_dat$night == night, 'sync' ] <- na_mean( dy_night_dat[ , 2 ], dy_night_dat[ , 3 ] )
    
  }
  
}

#### determine on which nights co-sleeping occurred between groups ####
sleep_site_dat <- unique( bdat[ , c( 'id', 'day', 'sleep_clus' ) ] )

sleep_site_dat$day <- as.character( sleep_site_dat$day )

## need to add underscores to the id's in sleep_sync_dat
sleep_sync_dat$tag_a <- gsub( ' ', '_', sleep_sync_dat$tag_a )
sleep_sync_dat$tag_b <- gsub( ' ', '_', sleep_sync_dat$tag_b )

sync_site_temp <- merge( x = sleep_sync_dat, y = sleep_site_dat, by.x = c( 'tag_a', 'night' ), by.y = c( 'id', 'day' ), all.x = T, all.y = F, sort = F )

names( sync_site_temp )[ names( sync_site_temp ) == 'sleep_clus' ] <- 'sleep_clus_a'

sync_site_dat <- merge( x = sync_site_temp, y = sleep_site_dat, by.x = c( 'tag_b', 'night' ), by.y = c( 'id', 'day' ), all.x = T, all.y = F, sort = F )

names( sync_site_dat )[ names( sync_site_dat ) == 'sleep_clus' ] <- 'sleep_clus_b'


sync_site_dat$group_a <- str_split_fixed( sync_site_dat$tag_a, '_', 3 )[ , 2 ]
sync_site_dat$group_b <- str_split_fixed( sync_site_dat$tag_b, '_', 3 )[ , 2 ]

sync_site_dat$same_site <- as.factor( as.numeric( sync_site_dat$sleep_clus_a == sync_site_dat$sleep_clus_b ) )
sync_site_dat$same_group <- as.factor( as.numeric( sync_site_dat$group_a == sync_site_dat$group_b ) )

sync_site_dat$dyad_name <- paste( sync_site_dat$tag_a, sync_site_dat$tag_b, sep = '__' )

sync_site_dat$cat <- ifelse( sync_site_dat$same_group == 0 & sync_site_dat$same_site == 0, 'diff_group_diff_site', ifelse( sync_site_dat$same_group == 0 & sync_site_dat$same_site == 1, 'diff_group_same_site', 'same_group_same_site' ) )

sync_site_dat$cat <- factor( sync_site_dat$cat, levels = c( 'diff_group_diff_site', 'diff_group_same_site', 'same_group_same_site' ) )

hist( sync_site_dat$sync, breaks = 100 )

sync_sub <- sync_site_dat[ !is.na( sync_site_dat$sync ) & !is.na( sync_site_dat$same_site ) & sync_site_dat$same_group != 1, ]


options( mc.cores = parallel::detectCores() )

sync_sub$sync_std <- normalize_func( sync_sub$sync )

priors <- c( prior(normal(0,2), class = "b", coef = same_site1 ) )

sleep_sync_mod <- brm( sync_std ~ same_site + ( 1| tag_a ) + ( 1| tag_b ) + ( 1 | dyad_name ), data = sync_sub, family= "gaussian", iter = 5000, prior = priors, control = list(max_treedepth = 15, adapt_delta = .999999999999 ) )

saveRDS( sleep_sync_mod, 'RESULTS/models/sleep_sync_mod.rds' )

sleep_sync_mod <- readRDS( 'RESULTS/models/sleep_sync_mod.rds' )

tab_model( sleep_sync_mod, show.ci = .90 )

pp_check( sleep_sync_mod )




posterior <- posterior_samples( sleep_sync_mod )

mean( posterior$b_same_site1*sd( sync_sub$sync ) )

coef_dens <- density( posterior$b_same_site1*sd( sync_sub$sync ) )

par( bg = 'black' )

plot( coef_dens$x, coef_dens$y, type = 'l', xlab = 'Increase in synchronization when sleeping together', bty = 'l', las = 1, ylab = 'Posterior probability density', col = 'white', col.axis = 'white', col.lab = 'white', col.main = 'white' )
axis( 1, labels = F, col = 'white' )
axis( 2, labels = F, col = 'white' )

rel_quantiles <- quantile( posterior$b_same_site1*sd( sync_sub$sync ), c( 0.05, 0.95 ) )

abline( v = rel_quantiles, lty = 2, col = 'white' )

poly_x <- c( rel_quantiles[ 1 ], coef_dens$x[ coef_dens$x > rel_quantiles[ 1 ] & coef_dens$x < rel_quantiles[ 2 ] ], rel_quantiles[ 2 ] )
poly_y <- c( -5, coef_dens$y[ coef_dens$x > rel_quantiles[ 1 ] & coef_dens$x < rel_quantiles[ 2 ] ], -5 )

polygon( x = poly_x, y = poly_y, col = transp( 'white', 0.3), border = F )


dev.off()

h1 = hypothesis(sleep_sync_mod, c( "same_site1 > 0", "same_site1 < 0" ) )

h1                         







