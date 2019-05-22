## ******************************************************************** ##
## fral_edit_MP_files.R
##
## Author: Matthew Aiello-Lammens
## Date Created: 2013-12-15
##
## Purpose:
## This script adapted from fral_misc_demog.R. It is written to create 
## the base MP files which are used in the sensitivity analysis runs.
## ******************************************************************** ##

##**** WARNING ****#
## Re-running any of the code in this file may result in overwriting
## files. The SA scripts are written to *NOT* overwrite files, but 
## other parts of the code below do not have the same FAIL SAFES
##**** WARNING ****#

## ******************************************************************** ##
## Source requisite packages and files
## ******************************************************************** ##

source('R/fral_demog_setup.R')

## ******************************************************************** ##
## Setup iniital population sizes
## -------------------------------------------------------------------- ##
## Need to determine which populations are occupied by 1910 then set
## initial populaitons sizes
## ******************************************************************** ##

set_initial_abundance <- FALSE

if( set_initial_abundance) {
  ## Read in patch id maps
  patch_raster <- '/Users/mlammens/Dropbox/F-alnus/Chapter-5-FRAL-Demography/ramas/fral_patch.ASC'
  patch_ids <- raster(patch_raster)
  
  ## Get csv of occurence records through time
  occurences <- read.csv('ramas/FRAL_Occurence_LAEA.csv')
  occurences_pre1910 <- occurences[ !occurences$YEAR>1910, ]
  
  ## Get populations that were occupied pre-1910
  pops_pre1910 <- extract(patch_ids,y=cbind(occurences_pre1910$LON_LAEA,occurences_pre1910$LAT_LAEA))
  pops_pre1910
  
  ## Assign the 0 patch value to 2559
  pops_pre1910[1] <- 2559
  # Only want unique pops
  pops_pre1910 <- unique(pops_pre1910)
  
  ## Read in the mp file
  mp_file <- 'ramas/frallo.MP'
  mp <- mp.read(mpFile=mp_file)
  mp_ver <- mp$version
  mp <- mp$mp.file
  
  ## Set initial abundance for these populations to 10,000
  mp$PopData_df[pops_pre1910,4] <- 10000
  
  ## Write the new mp file
  mp.write(mp.new=mp,version=mp_ver,mp.new.file='ramas/frallo.mp')
  print( 'Converting *.mp file to Windows format' )
  system( paste(sens.base.dir,'nix2win.sh ','ramas/frallo.mp', sep="") )
  
}

## ******************************************************************** ##
## Read mp file with user defined DLL
## ******************************************************************** ##

# ## Read in the mp file
# mp_file <- 'ramas/fral_sp_DensEff.mp'
# mp <- mp.read(mpFile=mp_file)
# mp_ver <- mp$version
# mp <- mp$mp.file
# 
# View(mp$PopData_df)

## ******************************************************************** ##
## Sensisivity analysis: Setting upper and lower bounds
## ******************************************************************** ##

calc_param_limits <- FALSE

if (calc_param_limits){
  ## WARNING:
  ## ********
  ## Before running the code below, run all of the code chunks in the 
  ## Fral_Demog_Methods.Rmd script.
  
  ## -------------------------------------------------------------------- ##
  ## FECUNDITY
  ## -------------------------------------------------------------------- ##
  
  # Set establishment probability - see manuscript for explinations
  # Estab_prob = Seeds per fruit * germination rate * within season survival
  establishment_prob_lo <- 2*.05*0.0177
  establishment_prob_hi <- 3*0.90*0.9311
  
  ## Function copied from Fral_Demog_Methods.Rmd, and modified 
  ## Fecundity 
  ## -- This function assumes that I'm using the ANCOVA
  ## -- results using DAH super class as a categorical predictor
  
  f.yx.mod <- function(x,establishment_prob,ancova.fit,upper=TRUE) {
    # Get standard deviations of the coefficients
    coef_sd <- unname(sqrt(diag(vcov(ancova.fit))))
    if (upper){
      fruit.by.class <- coefficients(ancova.fit)[1] + coefficients(ancova.fit)[2:3]
      fruit.by.class <- fruit.by.class + coef_sd[1]
      fruit.by.class <- fruit.by.class + coef_sd[2:3]
      fruit.by.class <- c((coefficients(ancova.fit)[1]+coef_sd[1]),fruit.by.class)
    } else {
      fruit.by.class <- coefficients(ancova.fit)[1] + coefficients(ancova.fit)[2:3]
      #fruit.by.class <- fruit.by.class - coef_sd[1]
      #fruit.by.class <- fruit.by.class - coef_sd[2:3]
      #fruit.by.class <- c((coefficients(ancova.fit)[1]-coef_sd[1]),fruit.by.class)
      fruit.by.class <- c((coefficients(ancova.fit)[1]),fruit.by.class)
      
    }
    fruit.by.class <- exp(fruit.by.class)
    fruit.by.class <- c(0,fruit.by.class)
    fruit.by.class <- unname(fruit.by.class)
    
    # Assign density independent fruit values to each size class
    frt <- x*0
    frt[x<0.5] <- fruit.by.class[1]
    frt[(x>=0.5 & x<2)] <- fruit.by.class[2]
    frt[(x>=2 & x<3.5)] <- fruit.by.class[3]
    frt[(x>=3.5)] <- fruit.by.class[4]
    
    # Apply establishment probability
    fec <- frt * establishment_prob
    
    # Return fecundity vector
    return(fec)
  }
  
  ## Build high matrix
  ## -------------------------------------------------------------------- ##
  
  ## Calculate stage specific survivals
  # Assume low survival - i.e., all missing are mortalities
  params$surv.int=coefficients(surv.reg)[1]
  params$surv.slope=coefficients(surv.reg)[2]
  # Calcualte stage specific survivals
  S <- s.x(y,params=params)
  # Set S[1] = same survival as stage 1 -> stage 2. 
  # This will only be used for the first couple of row elements,
  # while the rest of the row repsents
  # establishment rate, which is accounted for in the Fecundity
  # function
  S[1] <- seedling.betweenYr.surv
  # Combine surivial and growth
  P <- G
  for( i in 1:n) { P[,i]<-G[,i]*S[i] }
  # Set some of the matrix values to 0, as they're 
  # not biologically reasonable.
  P[which(P<1e-5)] <- 0
  # Also set elements of the first Row associated with feconudity 
  # to 0
  P[1,which(y>0.5)] <- 0
  
  ## Calculate fecundity
  Fec.row1.hi <- 
    f.yx.mod(y,establishment_prob_hi,fruit.by.dens.ancova,upper=TRUE)
  # Make a fecundity matrix to add to the Growth-Survival matrix
  Fec <- P*0
  Fec[1,] <- Fec.row1.hi
  # Combine Surivial-Growth and Fecundity kernels
  K.hi=P+Fec
  eigen(K.hi)$values[1]
  
  
  ## Build low matrix
  ## -------------------------------------------------------------------- ##
  
  ## Calculate stage specific survivals
  # Assume low survival - i.e., all missing are mortalities
  params$surv.int=coefficients(survLow.reg)[1]
  params$surv.slope=coefficients(survLow.reg)[2]
  # Calcualte stage specific survivals
  S <- s.x(y,params=params)
  # Set S[1] = same survival as stage 1 -> stage 2. 
  # This will only be used for the first couple of row elements,
  # while the rest of the row repsents
  # establishment rate, which is accounted for in the Fecundity
  # function
  S[1] <- seedling.betweenYr.surv
  # Combine surivial and growth
  P <- G
  for( i in 1:n) { P[,i]<-G[,i]*S[i] }
  # Set some of the matrix values to 0, as they're 
  # not biologically reasonable.
  P[which(P<1e-5)] <- 0
  # Also set elements of the first Row associated with feconudity 
  # to 0
  P[1,which(y>0.5)] <- 0
  
  ## Calculate fecundity
  Fec.row1.lo <-
    f.yx.mod(y,establishment_prob_lo,fruit.by.dens.ancova,upper=FALSE)
  # Make a fecundity matrix to add to the Growth-Survival matrix
  Fec <- P*0
  Fec[1,] <- Fec.row1.lo
  # Combine Surivial-Growth and Fecundity kernels
  K.lo=P+Fec
  eigen(K.lo)$values[1]
  
  ## -------------------------------------------------------------------- ##
  ## Estimate low and high slopes for DLL
  dll_slopes <- unname(coefficients(fruit.by.dens.ancova)[4:6])
  dll_slopes_sd <- unname(sqrt(diag(vcov(fruit.by.dens.ancova))))[4:6]
  dll_slopes_lo <- dll_slopes - dll_slopes_sd
  dll_slopes_hi <- dll_slopes + dll_slopes_sd
  # Threshold the hi values at 0 (can't have positive density dependence)
  dll_slopes_hi <- ifelse(dll_slopes_hi<0,yes=dll_slopes_hi,no=0)
  
}

## ******************************************************************** ##
## Set parameters for low and high mp files

# Set version for all ("51")
mp_ver <- "51"

## Low
## ***
# Ceiling DD
mp_lo <- mp.read('ramas/frallo.MP')
mp_lo <- mp_lo$mp.file
# DLL DD
mp_lo_denseff <- mp.read('ramas/frallo_DensEff.MP')
mp_lo_denseff <- mp_lo_denseff$mp.file

# Set new transition matrix
mp_lo$StMatr[[1]]$Matr <- K.lo
mp_lo_denseff$StMatr[[1]]$Matr <- K.lo

# Set initial abundance values
mp_lo$PopData_df[pops_pre1910,4] <- 1000
mp_lo_denseff$PopData_df[pops_pre1910,4] <- 1000

# Set DLL pars - ONLY DO THIS FOR DENSEFF MODELS
mp_lo_denseff$PopData_df[,c('udd_2','udd_3','udd_4')] <- 
  matrix(rep(dll_slopes_lo,nrow(mp_lo_denseff$PopData_df)),nrow=nrow(mp_lo_denseff$PopData_df),byrow=TRUE)

# Set KCH Scenarios
mp_lo$PopData_df$KchangeSt <- 
  sub(pattern='FRALL',replacement='fral_ceil_low_',mp_lo$PopData_df$KchangeSt)
mp_lo_denseff$PopData_df$KchangeSt <- 
  sub(pattern='FRALL',replacement='fral_udd_low_',mp_lo_denseff$PopData_df$KchangeSt)

## -------------------------------------------------------------------- ##
## Write the new mp files to disk
mp.write(mp.new=mp_lo,version=mp_ver,mp.new.file='ramas/frallo.MP')
mp.write(mp.new=mp_lo_denseff,version=mp_ver,mp.new.file='ramas/frallo_DensEff.MP')

# Convert to windows format
system( paste(sens.base.dir,'nix2win.sh ','ramas/frallo.MP', sep="") )
system( paste(sens.base.dir,'nix2win.sh ','ramas/frallo_DensEff.MP', sep="") )

## -------------------------------------------------------------------- ##
## High
## ****
# Ceiling DD
mp_hi <- mp.read('ramas/fralhi.MP')
mp_hi <- mp_hi$mp.file
# DLL DD
mp_hi_denseff <- mp.read('ramas/fralhi_DensEff.MP')
mp_hi_denseff <- mp_hi_denseff$mp.file

# Set new transition matrix
mp_hi$StMatr[[1]]$Matr <- K.hi
mp_hi_denseff$StMatr[[1]]$Matr <- K.hi

# Set initial abundance values
mp_hi$PopData_df[pops_pre1910,4] <- 10000
mp_hi_denseff$PopData_df[pops_pre1910,4] <- 10000

# Set DLL pars - ONLY DO THIS FOR DENSEFF MODELS
mp_hi_denseff$PopData_df[,c('udd_2','udd_3','udd_4')] <- 
  matrix(rep(dll_slopes_hi,nrow(mp_hi_denseff$PopData_df)),nrow=nrow(mp_hi_denseff$PopData_df),byrow=TRUE)

# Set KCH Scenarios
mp_hi$PopData_df$KchangeSt <- 
  sub(pattern='FRALL',replacement='fral_ceil_high_',mp_hi$PopData_df$KchangeSt)
mp_hi_denseff$PopData_df$KchangeSt <- 
  sub(pattern='FRALL',replacement='fral_udd_high_',mp_hi_denseff$PopData_df$KchangeSt)

## -------------------------------------------------------------------- ##
## Write the new mp files to disk
mp.write(mp.new=mp_hi,version=mp_ver,mp.new.file='ramas/fralhi.MP')
mp.write(mp.new=mp_hi_denseff,version=mp_ver,mp.new.file='ramas/fralhi_DensEff.MP')

# Convert to windows format
system( paste(sens.base.dir,'nix2win.sh ','ramas/fralhi.MP', sep="") )
system( paste(sens.base.dir,'nix2win.sh ','ramas/fralhi_DensEff.MP', sep="") )

## -------------------------------------------------------------------- ##
## Make medium scenarios using the low estimates. The only MED values
## are the KCH values
mp_med <- mp_lo
mp_med_denseff <- mp_lo_denseff

# Set KCH Scenarios
mp_med$PopData_df$KchangeSt <- 
  sub(pattern='FRALL',replacement='fral_ceil_med_',mp_med$PopData_df$KchangeSt)
mp_med_denseff$PopData_df$KchangeSt <- 
  sub(pattern='fral_udd_low_',replacement='fral_udd_med_',mp_med_denseff$PopData_df$KchangeSt)

## -------------------------------------------------------------------- ##
## Write the new mp files to disk
mp.write(mp.new=mp_med,version=mp_ver,mp.new.file='ramas/fralmed.MP')
mp.write(mp.new=mp_med_denseff,version=mp_ver,mp.new.file='ramas/fralmed_DensEff.MP')

# Convert to windows format
system( paste(sens.base.dir,'nix2win.sh ','ramas/fralmed.MP', sep="") )
system( paste(sens.base.dir,'nix2win.sh ','ramas/fralmed_DensEff.MP', sep="") )

## ******************************************************************** ##
## ******************************************************************** ##
## Run sensitivty analysis
## ******************************************************************** ##
## ******************************************************************** ##

run_sa <- FALSE

if (run_sa){
  setwd('ramas')
  
  # ## Start sensitivity analysis for fralA set
  # sensitivity('sens.config_fral.txt')
  # 
  # ## Start sensitivity analsysi for fralB set
  # sensitivity('sens.config_fralB.txt')
  
  ## Start sensitivity analysis for fral_denseff set
  sensitivity('sens.config_fral_denseff.txt')
  
  ## Start sensitivity analysis for fral_denseff set
  sensitivity('sens.config_fral_ceiling.txt')
  
  ## Start sensitivity analysis for fral_denseff_2 set
  ## --This set was run after the others, but if everything
  ## --is to be re-run, it could go here. Alternatively, the other
  ## --sensitivity analyses could be re-run with 1000 replicates
  ## Start sensitivity analysis for fral_denseff set
  sensitivity('sens.config_fral_denseff_2.txt')
  
}

## ******************************************************************** ##
## Make fral_nochng MP files
## ******************************************************************** ##

##**** WARNING ****#
## This script uses regular expressions to make substitutions in the KCH
## file names. Make sure you don't overwrite any MP files with non-existant 
## KCH files
##**** WARNING ****#

make_fral_nochng <- FALSE

if (make_fral_nochng){
  # Get list of models
  fral_mpFiles <- list.files(path='ramas/',pattern='fral_[ceil|dens].*[0-9].*mp',full.names=TRUE)
  
  ## Setup a cluster to change the mp files
  require(foreach)
  require(doSNOW)
  cl <- makeCluster(3,"SOCK")
  registerDoSNOW(cl)
  
  foreach (f_ind = 1:length(fral_mpFiles),.packages=c('lhs','fields')) %dopar% { #length(fral_mpFiles)
    mp <- mp.read(mpFile=fral_mpFiles[f_ind])
    mp <- mp$mp.file
    
    ## Change KCH file names
    mp$PopData_df$KchangeSt <-
      sub(pattern='FRAL',replacement='fral_nochng',mp$PopData_df$KchangeSt)
    
    ## Write "no change" model
    ## -------------------------------------------------------------------- ##
    # Make a new model name
    fral_mp_temp <- sub(pattern='fral',replacement='fral_nochng',fral_mpFiles[f_ind])
    # Write the model
    mp.write(mp.new=mp,version='51',mp.new.file=fral_mp_temp)
    # If on unix system, convert to Windows (comment out otherwise)
    system( paste(sens.base.dir,'nix2win.sh ',mp.new.file=fral_mp_temp, sep="") )
  }
  
  # Stop cluster
  stopCluster(cl)
}

## -------------------------------------------------------------------- ##
## Edit batch files

# Get batch files
batch_files <- list.files(path='ramas/',pattern='fral_batch_.*[ceil|denseff][0-9].*bat',full.names=TRUE)
# Remove Rerun files
batch_files <- batch_files[!grepl(pattern='Rerun',batch_files)]
# Read in all batch files
batch_all <- mclapply(X=batch_files,FUN=readLines)
# Change 'fral' to 'fral_nochng'
batch_all_new <- lapply(batch_all,function(x){sub(pattern='fral',replacement='fral_nochng',x)})
# Write these files
for ( x in 1:length(batch_files) ){
  # Make a new batch file name
  batch_new <- sub(pattern='fral',replacement='fral_nochng',batch_files[x])
  # Write the new file
  write.table(batch_all_new[[x]],file=batch_new,
              quote=FALSE,row.names=FALSE,col.names=FALSE)
}