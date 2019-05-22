## ******************************************************************** ##
## fral_link_patch2pop.R
##
## Author: Matthew Aiello-Lammens
## Date Created: 2013-12-09
##
## Purpose:
## Link the patch structure to the population size. Patch structure is
## determined by applying RAMAS Spatial/Patch to the SDM projections and
## the Population Size is determined by the results of a RAMAS Metapop
## simulation run.
## ******************************************************************** ##

#mp_file <- 'ramas/fral_denseff_add_popd_LDD499.mp' # Low LOSS example
mp_file <- 'ramas/fral_denseff_add_popd_LDD265.mp' # High Sensitivity example
patch_raster <- 'ramas/fral_patch.ASC'

## ******************************************************************** ##
## Source requisite packages and files
## ******************************************************************** ##
## Source the fral_demog_setup.R script to setup required functions
source('R/fral_demog_setup.R')

# R Packages
require(raster)


## ******************************************************************** ##
## Read in files
## ******************************************************************** ##

# Read mp files
mp <- mp.read.results(mpFile=mp_file)
mp <- mp$mp.file
mp_res <- mp$results
  
# Read patch raster
patches <- raster(patch_raster)

## ******************************************************************** ##
## Reclassify raster values:
## Replace patch number values in raster with final population size 
## values
## ******************************************************************** ##

pop_size_increment <- c( 1, 20, 40, 60, 80, 100 )

for ( t in pop_size_increment ) {
  # Get population size values
  final_popVals <- mp_res$PopInd[ t ,1,]
  # Take the log of the population size
  final_popVals <- log(final_popVals)
  # Set any population with a 0 (-Inf) value to -9999
  final_popVals[ final_popVals==-Inf ] <- -9999
  
  # Get the patch ids from the raster
  patch_ids <- as.data.frame(patches)
  length(which(patch_ids>0))
  # ### TEMP - remove some of the patches because of the difference between
  # ### TEMP - the most recent patch structure and that used in the mp run
  # final_popVals <- c(final_popVals,rep(0,15))
  
  patch_finalPop <- subs(patches,data.frame(id=1:3423,v=final_popVals) ,
  #filename=paste( "data_gis/fral_", t, "_popSize_lowLoss.asc", sep=""),
  filename=paste( "data_gis/fral_", t, "_popSize_highRes.asc", sep=""),
  overwrite=TRUE)
  
  #plot( patch_finalPop )
}


## -------------------------------------------------------------------- ##
## Make a carrying capacity patch example
# Get population size values
patch_k <- mp$PopData_df$K

# Get the patch ids from the raster
patch_ids <- as.data.frame(patches)
length(which(patch_ids>0))

patch_k_map <- subs(patches,data.frame(id=1:3423,v=patch_k) ,
                    #filename=paste( "data_gis/fral_", t, "_popSize_lowLoss.asc", sep=""),
                    filename=paste( "data_gis/fral_k.asc", sep=""),
                    overwrite=TRUE)
