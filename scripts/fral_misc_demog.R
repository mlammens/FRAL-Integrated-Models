## ******************************************************************** ##
## fral_misc_demog.R
##
## Author: Matthew Aiello-Lammens
## Date Created: 2013-12-15
##
## Purpose:
## Carryout miscellaneous tasks for demographic modeling of F. alnus
## ******************************************************************** ##

## ******************************************************************** ##
## Source requisite packages and files
## ******************************************************************** ##

source('R/fral_demog_setup.R')

## ******************************************************************** ##
## Determining dispersal rates
## ******************************************************************** ##

# Make dispersal-distance function
disp <- function(D,a,b,c){
  M <- a*exp(-(D^c)/b)
  return(M)
}

D <- c(0,20,28,40,44.7,56.5,60,63.2)#,72.1,84.8)
cnt <- c(1,4,4,4,8,4,4,8)#,8,4)
M_lo <- disp(D,a=0.04,b=5,c=0.75) #Low
M_hi <- disp(D,a=0.07,b=5,c=0.6) #High
#M <- disp(D,a=0.06,b=5,c=0.4) #High

# ## Testing uncertainty
# dpars_l <- c(.04,5,.75)
# dpars_h <- c(.07,5,.6)
# ((dpars_h-dpars_l)*runif(1))+dpars_l
# M <- disp(D,a=0.07,b=5,c=0.6)
# M <- disp(D,a=0.068,b=5,c=0.61)
# sum(M*cnt)

View(data.frame(dist=D,disp=M_lo))
View(data.frame(dist=D,disp=M_hi))
sum(M_lo*cnt)
sum(M_hi*cnt)
plot(D,M_hi,ylim=c(0,0.3))
points(D,M_lo,pch=19)

# Total dispersal of low *outside* of plot
sum(M_lo*cnt) - 0.04
# Total dispersal of high *outside* of plot
sum(M_hi*cnt) - 0.07

## ******************************************************************** ##
## Modify a template mp file used in the spatial analysis of FRAL
## ******************************************************************** ##

mp_file <- 'ramas/fral_ceiling.mp'
mp <- mp.read(mpFile=mp_file)
mp_ver <- mp$version
mp <- mp$mp.file

## Change StMig so that only first stage migrates / disperses
mp$StMig <- c(1,rep(0,49))

## Set replicates to 10
mp$MaxRep <- 10

## Set Dispersal depends on targe K
mp$DispersalDependsOnTargetPopK <- "Yes"

## -------------------------------------------------------------------- ##
## Write the new mp file
mp.write(mp.new=mp,version=mp_ver,mp.new.file='ramas/fral_ceiling.mp')

print( 'Converting *.mp file to Windows format' )
system( paste(sens.base.dir,'nix2win.sh ','ramas/fral_ceiling.mp', sep="") )

## ******************************************************************** ##
## Setup iniital population sizes
## -------------------------------------------------------------------- ##
## Need to determine which populations are occupied by 1910 then set
## initial populaitons sizes
## ******************************************************************** ##

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
mp_file <- 'ramas/fral.MP'
mp <- mp.read(mpFile=mp_file)
mp_ver <- mp$version
mp <- mp$mp.file

## Set initial abundance for these populations to 10,000
mp$PopData_df[pops_pre1910,4] <- 10000

## Write the new mp file
mp.write(mp.new=mp,version=mp_ver,mp.new.file='ramas/fral_InitAb.mp')
print( 'Converting *.mp file to Windows format' )
system( paste(sens.base.dir,'nix2win.sh ','ramas/fral_InitAb.mp', sep="") )

## ******************************************************************** ##
## Read mp file with user defined DLL
## ******************************************************************** ##

## Read in the mp file
mp_file <- 'ramas/fral_sp_DensEff.mp'
mp <- mp.read(mpFile=mp_file)
mp_ver <- mp$version
mp <- mp$mp.file

View(mp$PopData_df)

## Write the new mp file
mp.write(mp.new=mp,version=mp_ver,mp.new.file='ramas/TEMP.mp')
print( 'Converting *.mp file to Windows format' )
system( paste(sens.base.dir,'nix2win.sh ','ramas/TEMP.mp', sep="") )

## ******************************************************************** ##
## Sensisivity analysis: Setting upper and lower bounds
## ******************************************************************** ##

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

## -------------------------------------------------------------------- ##
## Write the new mp files to disk
mp.write(mp.new=mp_hi,version=mp_ver,mp.new.file='ramas/fralhi.MP')
mp.write(mp.new=mp_hi_denseff,version=mp_ver,mp.new.file='ramas/fralhi_DensEff.MP')

# Convert to windows format
system( paste(sens.base.dir,'nix2win.sh ','ramas/fralhi.MP', sep="") )
system( paste(sens.base.dir,'nix2win.sh ','ramas/fralhi_DensEff.MP', sep="") )

## ******************************************************************** ##
## -------------------------------------------------------------------- ##
## Run sensitivty analysis

# setwd('ramas')
# 
# ## Start sensitivity analysis for fralA set
# sensitivity('sens.config_fral.txt')
# 
# ## Start sensitivity analsysi for fralB set
# sensitivity('sens.config_fralB.txt')

## ******************************************************************** ##
## ******************************************************************** ##
## Fix fral_ceil mps using HIGH K values
## -------------------------------------
## In these models, some populations exceed RAMAS maximum carrying
## capacity allowed (slightly greater than 2.1e9)
## ******************************************************************** ##
## ******************************************************************** ##

## Determine which models did not get run
fral_ceil_mpFiles <- list.files(path='ramas/',pattern='fral_ceil[0-9].*mp',full.names=TRUE)
fral_ceil_mpFiles_run <- list.files(path='ramas/',pattern='fral_ceil[0-9].*SCL',full.names=TRUE)
fral_ceil_mpFiles_run <- sub(pattern='SCL',replacement='mp',x=fral_ceil_mpFiles_run)
fral_ceil_mpFiles_ToRun <- setdiff(fral_ceil_mpFiles,fral_ceil_mpFiles_run)

## Setup a cluster to change the mp files
require(foreach)
require(doSNOW)
cl <- makeCluster(3,"SOCK")
registerDoSNOW(cl)

foreach (fc_ind = 1:length(fral_ceil_mpFiles_ToRun),.packages=c('lhs','fields')) %dopar% { #
  mp <- mp.read(mpFile=fral_ceil_mpFiles_ToRun[fc_ind])
  mp <- mp$mp.file
  
  # Fix K in PopData_df
  mp$PopData_df$K <- ifelse(mp$PopData_df$K>2100000000,yes=2100000000,no=mp$PopData_df$K)
  
  # Write fixed model
  mp.write(mp.new=mp,version='51',mp.new.file=fral_ceil_mpFiles_ToRun[fc_ind])
  system( paste(sens.base.dir,'nix2win.sh ',mp.new.file=fral_ceil_mpFiles_ToRun[fc_ind], sep="") )
}

# Stop cluster
stopCluster(cl)

# Write new *.bat files to run these models
## I did this by hand after exporting the list of fral_ceil_mpFiles_ToRun

## ******************************************************************** ##
## Fix fral_ceil mps using MED K values
## -------------------------------------
## In these models, I failed to properly set the fralmed.MP file to call
## the Medium KCH files. Instead they call the low kch files.
## ******************************************************************** ##

# First fix the template
mp <- mp.read(mpFile='ramas/fralmed.MP')
mp <- mp$mp.file

## Change KCH file names
mp$PopData_df$KchangeSt <-
  sub(pattern='low',replacement='med',mp$PopData_df$KchangeSt)

# Write fixed model
mp.write(mp.new=mp,version='51',mp.new.file='ramas/fralmed.MP')
system( paste(sens.base.dir,'nix2win.sh ',mp.new.file='ramas/fralmed.MP', sep="") )


## -------------------------------------------------------------------- ##
## Fix all other files

mp_files_toReRun <- Sim_Results$mp.file[grep('MED',Sim_Results$kch.file[501:1000])]

setwd('ramas')

## Setup a cluster to change the mp files
require(foreach)
require(doSNOW)
cl <- makeCluster(3,"SOCK")
registerDoSNOW(cl)

foreach (fc_ind = 1:length(mp_files_toReRun),.packages=c('lhs','fields','popbio')) %dopar% { #
  mp <- mp.read(mpFile=mp_files_toReRun[fc_ind])
  mp <- mp$mp.file
  
  ## Change KCH file names
  mp$PopData_df$KchangeSt <-
    sub(pattern='LOW',replacement='med',mp$PopData_df$KchangeSt)
  
  # Write fixed model
  mp.write(mp.new=mp,version='51',mp.new.file=mp_files_toReRun[fc_ind])
  system( paste(sens.base.dir,'nix2win.sh ',mp.new.file=mp_files_toReRun[fc_ind], sep="") )
}

# Stop cluster
stopCluster(cl)


## -------------------------------------------------------------------- ##
## Need to do the same for the no_change scenarios as well :(
mp_files_nochng_toReRun <- sub(pattern='fral',replacement='fral_nochng',x=mp_files_toReRun)
mp_files_nochng_toReRun

## Setup a cluster to change the mp files
require(foreach)
require(doSNOW)
cl <- makeCluster(3,"SOCK")
registerDoSNOW(cl)

foreach (fc_ind = 1:length(mp_files_nochng_toReRun),.packages=c('lhs','fields','popbio')) %dopar% { #
  mp <- mp.read(mpFile=mp_files_nochng_toReRun[fc_ind])
  mp <- mp$mp.file
  
  ## Change KCH file names
  mp$PopData_df$KchangeSt <-
    sub(pattern='LOW',replacement='med',mp$PopData_df$KchangeSt)
  
  # Write fixed model
  mp.write(mp.new=mp,version='51',mp.new.file=mp_files_nochng_toReRun[fc_ind])
  system( paste(sens.base.dir,'nix2win.sh ',mp.new.file=mp_files_nochng_toReRun[fc_ind], sep="") )
}

# Stop cluster
stopCluster(cl)

## Make new batch files
## Determine which models did not get run
fral_nochng_ceil_mpFiles <- list.files(pattern='fral_nochng_ceil[0-9].*mp',full.names=TRUE)
fral_nochng_ceil_mpFiles_run <- list.files(pattern='fral_nochng_ceil[0-9].*SCL',full.names=TRUE)
fral_nochng_ceil_mpFiles_run <- sub(pattern='SCL',replacement='mp',x=fral_nochng_ceil_mpFiles_run)
fral_nochng_ceil_mpFiles_ToRun <- setdiff(fral_nochng_ceil_mpFiles,fral_nochng_ceil_mpFiles_run)

# Write these to a batch file
write.table(fral_nochng_ceil_mpFiles_ToRun,file='fral_nochng_batch_ceil_Rerun0.bat',
            quote=FALSE,row.names=FALSE,col.names=FALSE)
## Mannually changed after this writing in TextWrangler

## ******************************************************************** ##
## 2014-02-09 Combine the different cum_occupancy_mpMult extractions
## ******************************************************************** ##

## Read in saved cum_occupancy_mpMult files

## Load saved 'cum_occupancy_mpMult' as needed
load('ramas/cum_occupancy_mpMult_1.RData')
#load('ramas/cum_occupancy_mpMult_1000.RData')
#load('ramas/cum_occupancy_mpMult_2000.RData')

## Save cum_occupancy_mpMult as another different df temporalily
cum_occ_data <- cum_occupancy_mpMult

## Load saved additional data extractions
load('ramas/cum_occupancy_mpMult_add_1.RData')
#load('ramas/cum_occupancy_mpMult_add_1000.RData')
#load('ramas/cum_occupancy_mpMult_add_2000.RData')

## Save cum_occupancy_mpMult as a seperate df
cum_occ_data_add <- cum_occupancy_mpMult

## Load saved popd_LDD data extractions
load('ramas/cum_occupancy_mpMult_add_popd_1.RData')
#load('ramas/cum_occupancy_mpMult_add_popd_1000.RData')
#load('ramas/cum_occupancy_mpMult_add_popd_2000.RData')

## Use combined cumulative occupancy mpMult data sets
cum_occupancy_mpMult <- c(cum_occ_data,cum_occ_data_add,cum_occupancy_mpMult)

rm(cum_occ_data, cum_occ_data_add)

## Save the new combined data sets
save(cum_occupancy_mpMult,file='ramas/cum_occupancy_mpMult_all_1.RData',compress=TRUE)
#save(cum_occupancy_mpMult,file='ramas/cum_occupancy_mpMult_all_1000.RData',compress=TRUE)
#save(cum_occupancy_mpMult,file='ramas/cum_occupancy_mpMult_all_2000.RData',compress=TRUE)


