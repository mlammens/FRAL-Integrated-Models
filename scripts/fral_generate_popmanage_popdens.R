## ******************************************************************** ##
## fral_generate_popmange_popdens.R
##
## Author: Matthew Aiello-Lammens
## Date Created: 2013-12-17
##
## Purpose:
## Generate random population managment actions to apply to a set of 
## RAMAS Metapop models. In this script, movement of species is weighted
## by the population density of grid cells (Populations)
##
## ******************************************************************** ##

## Goal
## ****
## I want to generate a set of popution management actions that can be
## used to parameterize a set of *mp files. These actions need to be 
## generated such that the are completely replicable, so they can be 
## applied to pair simulation models.
##
## In this script, as opposed to `fral_generate_popmange.R`, random 
## transplants are weighted by human population density in the grid cells
## (patches) the species could potentially occupy.
##
## -------------------------------------------------------------------- ##
## Setup required packages and scripts
source('R/fral_demog_setup.R')
source('~/Dropbox/F-Alnus/Chapter-4/R/spatial_data_setup.R')

# Read in patch raster file
patch_raster <- 'ramas/fral_patch.ASC'
patches <- raster( patch_raster )

# Get a list of all of the human population density HYDE layers
hyde.popd.layers <- Sys.glob(file.path(FRAL_HYDE_LAYERS,"popd*laea.asc"))
# Remove the first five, which are for 1850 to 1890
hyde.popd.layers <- hyde.popd.layers[ -(1:5) ]
# Remove year 2005, use only 2000 for last ten years
hyde.popd.layers <- hyde.popd.layers[ -12 ]

## ******************************************************************** ##
## Read in hyde layers as raster stack. Extract values for only
## FRAL patches
## ******************************************************************** ##
pop_dens <- stack(hyde.popd.layers)

pop_dens_temp <- as.data.frame(pop_dens[[1]])
plot(pop_dens)

# Make a two dimensional matrix with population densities for each of the
# raster cells as rows and each popd raster as a column
pop_dens_mat <- c()
for ( pop_yr in 1:nlayers(pop_dens) ){
  pop_dens_temp <- as.data.frame(pop_dens[[pop_yr]])
  pop_dens_temp <- pop_dens_temp[[1]]
  # Add to matrix
  pop_dens_mat <- cbind(pop_dens_mat, pop_dens_temp)
}

# Get the row IDs for the FRAL patches
fral_patches <- as.data.frame(patches)
fral_patches <- which( fral_patches$fral_patch > 0 )

# Remove rows in pop_dens_mat not asssociated with a FRAL patch
pop_dens_mat <- pop_dens_mat[ fral_patches, ]

# Because of the conversion from WGS 1984 to LAEA, 
# some pop density values are less than 0. Because I am only
# interested in the relative density of patchs, I'm going to
# add the absolute value of the minimum + 0.1 to the whole
# matrix. The 0.1 value is to make it such that it is not
# impossible to select any patches
pop_dens_mat <-
  pop_dens_mat + ( abs( min(pop_dens_mat,na.rm=TRUE) ) + 0.1 )

# NB: There are 19 patches that have no population density
# information. I am going to set the pop density values for
# these populations to be the median values for that year
pop_dens_mat[ is.na(pop_dens_mat[,1]), ] <-
  apply(pop_dens_mat,MARGIN=2,FUN=median,na.rm=TRUE)

## ******************************************************************** ##
## Read in existing PopManage files - these were created by the 
## fral_generate_popmanage_popdens.R script.
## ******************************************************************** ##

# PopManage <- read.csv('ramas/PopManage.csv')
# Model_PopManage_Num <- read.csv('ramas/Model_PopManage_Num.csv')
PopManage <- read.csv('ramas/PopManage_2.csv')
Model_PopManage_Num <- read.csv('ramas/Model_PopManage_Num_2.csv')
Model_PopManage_Num <- Model_PopManage_Num$x

## ******************************************************************** ##
## Choose new target populations based on population density
## ******************************************************************** ##

# Set PopManage_Saved to TRUE if you have already generated the new 
# target populations.
PopManage_Saved <- TRUE

if( !PopManage_Saved ){
  # PopManage <- read.csv('ramas/PopManage.csv')
  # Model_PopManage_Num <- read.csv('ramas/Model_PopManage_Num.csv')
  PopManage <- read.csv('ramas/PopManage_2.csv')
  Model_PopManage_Num <- read.csv('ramas/Model_PopManage_Num_2.csv')
  Model_PopManage_Num <- Model_PopManage_Num$x
  
  ## Set preliminary values
  ## ----------------------
  # Initial pops that can be "pulled" from: 
  # These populations are always started with a non-zero initial abundance
  initial_pops <- c(2559,2820,1990,2145,1626,404,1308,1991,1827,1910,1992)
  
  # Total number of populations
  total_pops <- 3423
  
  # Get possible target populations
  target_pops <- 1:total_pops
  target_pops <- target_pops[-initial_pops]
  
  # Remove the source populations from `pop_dens_mat` as well
  pop_dens_mat <- pop_dens_mat[ -initial_pops, ]
  
  ## -------------------------------------------------------------------- ##
  # Setup the number of time steps. Use the modulus operator to seperate
  # out time steps into 10 year increments that match the hyde layers
  
  # Number of time steps
  time_steps <- 100
  decade_steps <- unique( ( 1:time_steps %/% 10 ) + 1 )
  
  # Make a decades steps column for the PopManage data.frame
  PopManage$decade_step <- ( PopManage$begin.time %/% 10 ) + 1
  
  # Loop through the decade steps and find the tranlocations for
  # each decade
  for ( dec in 1:length(decade_steps) ){
    # Find rows of PopManage data.frame associated with the 
    # current decade step
    PopManage_dec_rows <- which(PopManage$decade_step==dec)
    # Generate a new set of target poplation ideas, using the
    # population density to weight selection
    new_target_pops <- sample( target_pops,
                               size=length(PopManage_dec_rows),
                               replace=TRUE,
                               prob=pop_dens_mat[ ,dec] )
    # Assign these new values to the PopManage data.frame
    PopManage$to.pop[ PopManage_dec_rows ] <- new_target_pops    
  }

  ## Save the PopMange data.frame
  #write.csv(PopManage,file='ramas/PopManage.csv',row.names=FALSE,quote=FALSE)
  # UPDATE: Make a new PopManage set for next 500 models
  write.csv(PopManage,file='ramas/PopManage_2_popdens.csv',row.names=FALSE,quote=FALSE)
} else {
#   PopManage <- read.csv('ramas/PopManage.csv')
#   Model_PopManage_Num <- read.csv('ramas/Model_PopManage_Num.csv')
  PopManage <- read.csv('ramas/PopManage_2_popdens.csv')
  Model_PopManage_Num <- read.csv('ramas/Model_PopManage_Num_2.csv')
  Model_PopManage_Num <- Model_PopManage_Num$x
}

## ******************************************************************** ##
## Assign managment scenarios to population models
## ******************************************************************** ##

# Set indexes for which rows of PopManage should be associated
# with each of 500 simulations
Model_PopManage_Index_End <- cumsum(Model_PopManage_Num)
Model_PopManage_Index_Start <- Model_PopManage_Index_End+1
# Remove the last entry
Model_PopManage_Index_Start <- 
  Model_PopManage_Index_Start[-length(Model_PopManage_Index_Start)] 
# Begin at row 1
Model_PopManage_Index_Start <-
  c(1,Model_PopManage_Index_Start)


## -------------------------------------------------------------------- ##
# 
## On 2014-02-01 I created 500 models using the fral_denseff_add set as
## base models, and used the new target population values based on 
## human population density.

setwd('ramas/')

mp_files <- list.files(pattern='fral_denseff_add[0-9].*mp$')
for( ind in 1:length(mp_files)){
  # Read in temporary mp file
  mp <- mp_files[ind]
  mp_temp <- mp.read(mpFile=mp)
  mp_temp <- mp_temp$mp.file
  
  # Change the number of population management actions
  mp_temp$NPopManage <- Model_PopManage_Num[ind]
  mp_temp$PopManageProp <- 
    PopManage[Model_PopManage_Index_Start[ind]:Model_PopManage_Index_End[ind],]
  
  # Write new file
  mp_new <- sub(pattern='add',replacement='add_popd_LDD',mp)
  mp.write(mp.new=mp_temp,version="51",mp.new.file=mp_new)
  system( paste(sens.base.dir,'nix2win.sh ',mp_new, sep="") )
}
