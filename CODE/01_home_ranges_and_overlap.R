

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
library(rgeos)


############## Determine home ranges and home range overlap zones ##############

which_spec <- 'baboon'

# read in the data
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


# Make local_timestamp POSIX class 
spec_df$local_timestamp<-as.POSIXct(spec_df$local_timestamp, format="%Y-%m-%d %H:%M:%OS", tz="UTC")

### We want a track and home range at the group level, not the individual level. So remove any duplicate local_timestamps at the group level. We are going to treat the monkey groups as individuals here, so temporarily overwrite the id column with the information in the group column, and then remove duplicates. This is okay to have one group represented by different monkeys in this code, because we don't need to maintain the movement characteristics exactly to calculate home ranges

if(which_spec != "leopard"){
  spec_df$group <- str_split_fixed(spec_df$id, "_", 3)[,2]
  spec_df$id <- spec_df$group
}

spec_df <- spec_df[ !duplicated( spec_df[, c( 'id', 'local_timestamp' ) ] ), ]

# now that we have removed some of locations that were duplicates when analyzing at the group level, we need to reorder the dataframe by timestamp and id
spec_df <- spec_df[ order( spec_df$local_timestamp ), ]
spec_df <- spec_df[ order( spec_df$id ), ]



### finding home ranges using the adehabitatHR package

crs_longlat <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

crs_utm <- CRS("+proj=utm +zone=37 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

## create a SpatialPointsDataFrame to use for finding the home ranges with adeHabitat
new_sp <- SpatialPointsDataFrame(coords = spec_df[, c('x','y')], proj4string = crs_utm, data = spec_df) 

## First, estimate the utilization density
kud <- kernelUD(new_sp[, 'id'],
                h = "href", 
                same4all = F, 
                grid = 100,
                hlim = c(0.1, 2.0),
                kern = c("bivnorm"),
                boundary = NULL,
                extent = 0.2)

## Calculating home ranges from kernel utilization densities

## specify that we want the 95% KDE as our home range
level <- 95

homerange <- getverticeshr(kud,
                           percent = level,
                           unin = 'm',
                           unout = 'm2') # get the home ranges for each group

dir.create( paste0( getwd(), "/RESULTS/" ) )

dir.create( paste0( getwd(), "/RESULTS/home_ranges/" ) )

dir.create( paste0( getwd(), "/RESULTS/home_ranges/adeHabitat" ) )

dir.create( paste0( getwd(), "/RESULTS/home_ranges/adeHabitat/", which_spec ) )

dir.create( paste0( getwd(), "/RESULTS/home_ranges/adeHabitat/", which_spec, '/home_ranges/' ) )

### write the SpatialPolygonsDataFrame as a shapefile
writeOGR( homerange, paste( "RESULTS/home_ranges/adeHabitat/", which_spec, '/home_ranges/full_data', sep = ""), layer = paste("home_ranges_", level, sep = ''), driver="ESRI Shapefile")

## now just using the 99% KDE home ranges to find the home range overlap zones

#for each dyad, find the polygon of intersection between the dyad's home range polygons
for(i in 1 : (nrow(homerange)-1) ){
  
  for(j in (i+1) : nrow(homerange) ){
    
    #for each dyad, find the polygon of intersection between the dyad's home range polygons
    intersection <- raster::intersect(homerange[i,]
                                      , homerange[j,])
    
    if(!is.null(intersection)){
      
      ## if there is an intersection between the polygons, write the home range overlap zone polygon as a shapefile
      writeOGR(intersection, paste0( getwd(), "/RESULTS/home_ranges/adeHabitat/", which_spec, "/HR_overlap", sep = ""), layer = paste("HR_overlap_", homerange$id[i], '_', homerange$id[j], sep = ''), driver="ESRI Shapefile")
    }
  }
}


### calculating home range overlap

which_spec <- 'baboon'

dat <- 'full'

level <- 95

crs_longlat <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

crs_utm <- CRS("+proj=utm +zone=37 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

# read in the homerange data

homerange <- readOGR( paste0( "RESULTS/home_ranges/adeHabitat/", which_spec, '/home_ranges/', dat, '_data' ), layer = paste("home_ranges_", level, sep = '') )


# read in the GPS data
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


# Make local_timestamp POSIX class 
spec_df$local_timestamp <- as.POSIXct( spec_df$local_timestamp, format="%Y-%m-%d %H:%M:%OS", tz="UTC" )

### We want a track and home range at the group level, not the individual level. So remove any duplicate local_timestamps at the group level. We are going to treat the monkey groups as individuals here, so temporarily overwrite the id column with the information in the group column, and then remove duplicates. This is okay to have one group represented by different monkeys in this code, because we don't need to maintain the movement characteristics exactly to calculate home ranges

if(which_spec != "leopard"){
  spec_df$group <- str_split_fixed(spec_df$id, "_", 3)[,2]
  spec_df$id <- spec_df$group
}

spec_df <- spec_df[ !duplicated( spec_df[, c( 'id', 'local_timestamp' ) ] ), ]

# now that we have removed some of locations that were duplicates when analyzing at the group level, we need to reorder the dataframe by timestamp and id
spec_df <- spec_df[ order( spec_df$local_timestamp ), ]
spec_df <- spec_df[ order( spec_df$id ), ]


#### diving into home range overlap ######

## create a SpatialPointsDataFrame to use for finding the home ranges with adeHabitat
new_sp <- SpatialPointsDataFrame(coords = spec_df[, c('x','y')], proj4string = crs_utm, data = spec_df) 

## transparency function

transp <- function( col, alpha = 0.5 ){
  
  res <- apply( col2rgb( col ), 2, function( c ) rgb( c[ 1 ]/255, c[ 2 ]/255, c[ 3 ]/255, alpha ) )
  
  return( res )
}

plot( homerange, col = transp( as.numeric( as.factor( homerange$id ) ) ) )

HRO_mat <- kerneloverlap( new_sp[ , 'id' ], meth="HR", percent = level, conditional = FALSE, grid = 400 ) # calculate the homerange overlap. This is the proportion of one focal group's (row IDs) homerange that overlaps with another group's (column IDs) homerange

BA_mat <- kerneloverlap( new_sp[ , 'id' ], meth="BA", percent = level, conditional = FALSE, grid = 400 ) # similar to the volume of intersection but takes both UDs into account at all locations (with a joint probability) and not just the minimum UD at any given point. The best metric of the similarity of two UDs




