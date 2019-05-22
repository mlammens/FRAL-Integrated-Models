## ******************************************************************** ##
## fral_generate_popmange.R
##
## Author: Matthew Aiello-Lammens
## Date Created: 2013-12-17
##
## Purpose:
## Generate random population managment actions to apply to a set of 
## RAMAS Metapop models
##
## ******************************************************************** ##

## Goal
## ****
## I want to generate a set of popution management actions that can be
## used to parameterize a set of *mp files. These actions need to be 
## generated such that the are completely replicable, so they can be 
## applied to pair simulation models.
##
## I will do this by generating a very large set of population managment
## actions, saving this matrix, then loading it and partitioning it out
## to mp files in the same order for paired simulations.

## -------------------------------------------------------------------- ##
## Setup required packages and scripts
source('R/fral_demog_setup.R')

## -------------------------------------------------------------------- ##

PopManage_Saved <- TRUE

if( !PopManage_Saved ){
  ## -------------------------------------------------------------------- ##
  ## Read in a set of population management actions from one of the MP files
  ## used the GSA
  mp_popman <- mp.read('ramas/frallo.MP')
  mp_popman_ver <- mp_popman$version
  mp_popman <- mp_popman$mp.file
  
  # Get Pop Managment data.frame
  PopManage_base <- mp_popman$PopManageProp
  
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
  
  # Number of time steps
  time_steps <- 100
  
  ## Make a new data.frame of 250,000 PopMangeProp lines
  ## This number allows for a maximum of 500 translocations per simulation.
  ## 500 translocations over 100 years would equate to approximately 5 per year,
  ## provided the years are randomized
  popmanage_num <- 250000
  PopManage <- PopManage_base[rep(1,popmanage_num),]
  
  # Randomize source populations
  PopManage$from.pop <- sample(initial_pops,size=popmanage_num,replace=TRUE)
  
  # Randomize target population
  PopManage$to.pop <- sample(target_pops,size=popmanage_num,replace=TRUE)
  
  # Randomize **when** action happens
  popmanage_when <- sample(1:time_steps,size=popmanage_num,replace=TRUE)
  PopManage$begin.time <- popmanage_when
  PopManage$end.time <- popmanage_when
  
  ## Generate the number of tranlocations to occur in each simulation
  ## This should be a number between 0 and 500.
  Model_PopManage_Num <- sample(1:500,size=500,replace=TRUE)
  
  ## Save the PopMange data.frame
  #write.csv(PopManage,file='ramas/PopManage.csv',row.names=FALSE,quote=FALSE)
  # UPDATE: Make a new PopManage set for next 500 models
  write.csv(PopManage,file='ramas/PopManage_2.csv',row.names=FALSE,quote=FALSE)
  ## Save the Model_PopMange_Num 
  #write.csv(Model_PopManage_Num,file='ramas/Model_PopManage_Num.csv',row.names=FALSE,quote=FALSE)
  write.csv(Model_PopManage_Num,file='ramas/Model_PopManage_Num_2.csv',row.names=FALSE,quote=FALSE)
} else {
#   PopManage <- read.csv('ramas/PopManage.csv')
#   Model_PopManage_Num <- read.csv('ramas/Model_PopManage_Num.csv')
  PopManage <- read.csv('ramas/PopManage_2.csv')
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

# ## -------------------------------------------------------------------- ##
# 
# ## On 2013-12-24 I ran 500 models (fral_denseff1.mp to fral_denseff500.mp) using
# ## a new PopManage set
# 
# setwd('ramas/')
# 
# mp_files <- list.files(pattern='fral_denseff[0-9].*mp$')
# for( ind in 1:length(mp_files)){
#   # Read in temporary mp file
#   mp <- mp_files[ind]
#   mp_temp <- mp.read(mpFile=mp)
#   mp_temp <- mp_temp$mp.file
#   
#   # Change the number of population management actions
#   mp_temp$NPopManage <- Model_PopManage_Num[ind]
#   mp_temp$PopManageProp <- 
#     PopManage[Model_PopManage_Index_Start[ind]:Model_PopManage_Index_End[ind],]
#   
#   # Write new file
#   mp.write(mp.new=mp_temp,version="51",mp.new.file=mp)
#   system( paste(sens.base.dir,'nix2win.sh ',mp, sep="") )
# }
# 
# ## -------------------------------------------------------------------- ##
# 
# ## On 2013-12-24 I ran 500 models (fral_1.mp to fral_denseff500.mp) using
# ## a new PopManage set
# 
# #setwd('ramas/')
# 
# mp_files <- list.files(pattern='fral_ceil[0-9].*mp$')
# for( ind in 1:length(mp_files)){
#   # Read in temporary mp file
#   mp <- mp_files[ind]
#   mp_temp <- mp.read(mpFile=mp)
#   mp_temp <- mp_temp$mp.file
#   
#   # Change the number of population management actions
#   mp_temp$NPopManage <- Model_PopManage_Num[ind]
#   mp_temp$PopManageProp <- 
#     PopManage[Model_PopManage_Index_Start[ind]:Model_PopManage_Index_End[ind],]
#   
#   # Write new file
#   mp.write(mp.new=mp_temp,version="51",mp.new.file=mp)
#   system( paste(sens.base.dir,'nix2win.sh ',mp, sep="") )
# }

## -------------------------------------------------------------------- ##
# 
# ## On 2014-01-18 I ran 500 models (fral_denseff_add1.mp to fral_denseff_add500.mp) using
# ## a new PopManage set
# 
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
  mp.write(mp.new=mp_temp,version="51",mp.new.file=mp)
  system( paste(sens.base.dir,'nix2win.sh ',mp, sep="") )
}
