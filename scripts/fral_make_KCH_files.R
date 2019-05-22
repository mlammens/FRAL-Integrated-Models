## ******************************************************************** ##
## fral_make_KCH_files.R
##
## Author: Matthew Aiello-Lammens
## Date Created: 2013-12-23
##
## Purpose:
## Make KCH files for MED and LOW scenarios for user defined density
## dependence (UDD) models, and all three HIGH, MED, and LOW scenarios for
## Ceiling Carrying Capacity models.
## The original models constructed were for the HIGH UDD scenario
##
## Also make KCH files for No Land Change scenario. In this case, I am
## using the previously created KCH files, but each KCH file is a repeat
## of the first value
##
## ******************************************************************** ##

# source 'fral_deomg_setup.R' script
source('R/fral_demog_setup.R')

## ******************************************************************** ##
## Get list of KCH files
## ******************************************************************** ##

kch_files <- list.files(path="ramas",pattern="frall.*KCH",full.names=TRUE)

# Check that there are an equal number of KCH files as there are populaitons
if( length(kch_files) != 3423 ){
  stop("ERROR: Incorrect number of KCH files found by 'list.files' function.")
}

# Read all kch files
kch <- mclapply(X=kch_files,FUN=read.table,header=FALSE)

for (k in 1:length(kch_files)){
  # Rename KCH high
  kch_udd_high <- floor(kch[[k]])
  # Write new file
  kch_temp <- sub(pattern='frall',replacement='fral_udd_high_',kch_files[k])
  write.table(kch_udd_high,file=kch_temp,quote=FALSE,row.names=FALSE,col.names=FALSE)
  system( paste(sens.base.dir,'nix2win.sh ',kch_temp, sep="") )
  
  # Calculate kch values for MED and LOW UDD scenarios
  kch_udd_med <- floor(kch[[k]]/2)
  # Write new file
  kch_temp <- sub(pattern='frall',replacement='fral_udd_med_',kch_files[k])
  write.table(kch_udd_med,file=kch_temp,quote=FALSE,row.names=FALSE,col.names=FALSE)
  system( paste(sens.base.dir,'nix2win.sh ',kch_temp, sep="") )
  
  kch_udd_low <- floor(kch[[k]]/4)
  # Write new file
  kch_temp <- sub(pattern='frall',replacement='fral_udd_low_',kch_files[k])
  write.table(kch_udd_low,file=kch_temp,quote=FALSE,row.names=FALSE,col.names=FALSE)
  system( paste(sens.base.dir,'nix2win.sh ',kch_temp, sep="") )
  
  # Calculate kch values for High, Med, Low Ceiling scenarios
  kch_ceil_high <- floor(kch[[k]])*40
  # Limit maximum carrying capacity to be below max usable by RAMAS
  kch_ceil_high$V1 <- ifelse(kch_ceil_high$V1>2.1e9,yes=2100000000,no=kch_ceil_high$V1)
  # Write new file
  kch_temp <- sub(pattern='frall',replacement='fral_ceil_high_',kch_files[k])
  write.table(kch_ceil_high,file=kch_temp,quote=FALSE,row.names=FALSE,col.names=FALSE)
  system( paste(sens.base.dir,'nix2win.sh ',kch_temp, sep="") )
 
  kch_ceil_med <- floor(kch[[k]]/2)*40
  # Write new file
  kch_temp <- sub(pattern='frall',replacement='fral_ceil_med_',kch_files[k])
  write.table(kch_ceil_med,file=kch_temp,quote=FALSE,row.names=FALSE,col.names=FALSE)
  system( paste(sens.base.dir,'nix2win.sh ',kch_temp, sep="") )

  kch_ceil_low <- floor(kch[[k]]/4)*40
  # Write new file
  kch_temp <- sub(pattern='frall',replacement='fral_ceil_low_',kch_files[k])
  write.table(kch_ceil_low,file=kch_temp,quote=FALSE,row.names=FALSE,col.names=FALSE)
  system( paste(sens.base.dir,'nix2win.sh ',kch_temp, sep="") )
  
}

## ******************************************************************** ##
## Make No Change KCH files
## ------------------------
## Here I am taking advatange of the fact that RAMAS will repeat the K
## value for the last timestep in the KCH file for the duration of the
## simulation. Therefore, each file only needs to contain a single
## numeric value (K in year 1910) 
## ******************************************************************** ##

## Get all kch files
kch_files_all <- list.files(path='ramas',pattern='fral_[ceil|udd].*KCH',full.names=TRUE)

## Read in all KCH files
kch_all <- mclapply(X=kch_files_all,FUN=read.table,header=FALSE)

for (k in 1:length(kch_all)){
  # Get K for timestep 1
  k_first <- kch_all[[k]]$V1[1]
  # Write a new kch file
  kch_temp <- sub(pattern='fral',replacement='fral_nochng',kch_files_all[k])
  write.table(k_first,file=kch_temp,quote=FALSE,row.names=FALSE,col.names=FALSE)
  system( paste(sens.base.dir,'nix2win.sh ',kch_temp, sep="") )
}
