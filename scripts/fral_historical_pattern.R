## ******************************************************************** ##
## fral_historical_pattern.R
##
## Author: Matthew Aiello-Lammens
## Date Created: 2013-12-22
##
## Purpose:
## Create data sets that represent the historical pattern of spread of
## F. alnus (sensu Chapter 3).
## ******************************************************************** ##

# Source fral_demog_setup.R file
source('R/fral_demog_setup.R')

## ******************************************************************** ##
## Read in patch raster and historical occurence records
## -----------------------------------------------------
## This map provides a unique number to each patch (population) in the 
## metapopulation. 
## The historical occurence records are from collected herbarium records
## (Chapter 3). The occurence records file read here also includes all
## field observations as well.
## ******************************************************************** ##

## Read in patch id maps
# patch_raster <- 'ramas/fral_patch.ASC'
patch_raster <- "/Volumes/Garage/Projects/F-alnus/Chapter-5_ramas_folder/ramas/fral_patch.ASC"
patch_ids <- raster(patch_raster)
patch_total <- max(as.matrix(patch_ids))

## Get csv of occurence records through time. NB: using records in
## Lambert Equal Area projection.
occurences <- read.csv('/Volumes/Garage/Projects/F-alnus/Chapter-5_ramas_folder/ramas/FRAL_Occurence_LAEA.csv')
# Remove any occurences with no 'YEAR' information
occurences <- occurences[ !is.na(occurences$YEAR), ]
# Extract patch values for all occurence records
occurences$PatchID <- extract(patch_ids,cbind(occurences$LON_LAEA,occurences$LAT_LAEA))
# Remove patchs IDed as NA
occur_na <- occurences[ is.na(occurences$PatchID), ]
occurences <- occurences[ !is.na(occurences$PatchID), ]
# Seperate points that fall outside of populations
occur_outside <- occurences[ occurences$PatchID==0, ]
occurences <- occurences[ occurences$PatchID!=0, ]

## What percentage of points are being ommitted
nrow(occur_outside)/(nrow(occur_outside)+nrow(occurences)) # Omitting > 10%

# How many unique patches are occupied during all time?
length(unique(occurences$PatchID)) # 490 - This includes 2011 and 2012!
# what about only prior to 2011
length(unique(occurences[occurences$YEAR<=2010,]$PatchID)) # 459 
## Two years incraese this by 31 grid cells.

# What years are missing
#View(occur_outside[order(occur_outside$YEAR),])
## I am eliminating some early 20th century records in removing these 
## points, but an alternative approach would require IDing each of
## these by eye. A number of them are points georeferenced to the county
## level.

## ******************************************************************** ##
## Create a data.frame of temporal trends of occupancy for all patches
## ******************************************************************** ##

# Initiate cumulative occupancy data.frame
cum_occupancy_hist <- data.frame(Pop=1:patch_total)

for ( yr in 1910:2010 ){
  # Identify populations that are occupied by year 'yr'
  pop_occ <- 
    unique(occurences$PatchID[ occurences$YEAR <= yr ])
  # Make a vector of zeroes of length 'patch_total'
  pop_temp <- rep(0,patch_total)
  # Set patches found in 'pop_occ' equal to 1
  pop_temp[pop_occ] <- 1
  # Add this vector to the data.frame
  cum_occupancy_hist <- 
    cbind(cum_occupancy_hist,pop_temp)
}

## Set population 2559 to being occupied for all years
## This was one of the original populations in the mp models
cum_occupancy_hist[2559,2:102] <- 1

# Rename column names
yr_names <- paste('YR.',1910:2010,sep='')
names(cum_occupancy_hist) <- c('Pop',yr_names)


## ******************************************************************** ##
## Look at cummulative occupancy pattern across **all** patches
## ******************************************************************** ##
cum_occ_hist_allPops <- colSums(cum_occupancy_hist)
cum_occ_hist_allPops <- cum_occ_hist_allPops[-1]
cum_occ_hist_allPops_df <- data.frame(Year=1910:2010,CumOcc=cum_occ_hist_allPops)

cum_occ_lm2 <- lm(cum_occ_hist_allPops_df$CumOcc~cum_occ_hist_allPops_df$Year+
                    I(cum_occ_hist_allPops_df$Year^2))

ggplot(cum_occ_hist_allPops_df,aes(x=Year,y=sqrt(CumOcc))) +
  geom_point() +
  #stat_smooth(method='lm') +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), colour = "red") +
  theme_bw()

## This pattern appears to qualitavitly match what I found in Chapter 3. In 
## particular, the best fit linear model is a quadratic, concave up curve,
## indicating that the rate at which cells become occupied accelerates.

## ******************************************************************** ##
## Read in historical **ASSOCIATED SPECIES** records
## ******************************************************************** ##
## NB: The associated data used in chapter 3 were converted from WGS84 to
## LAEA in QGIS. Here I am reading in the LAEA data.
assoc_spec <- read.csv('/Volumes/Garage/Projects/F-alnus/Chapter-5_ramas_folder/ramas/Assoc_Occurence_LAEA.csv',as.is=TRUE)
# Make Collection year a numeric value - this removes some vales
assoc_spec$CllctnY <- as.numeric(assoc_spec$CllctnY)
# Remove records with no Observation Data
assoc_spec <- assoc_spec[!is.na(assoc_spec$CllctnY),]

# Get 20x20 km patches that associated data are located in
assoc_spec$PatchID <- extract(patch_ids,cbind(assoc_spec$Lon_LAEA,assoc_spec$Lat_LAEA))
# Remove records that did not occur in one of the patches
assoc_spec <- assoc_spec[!is.na(assoc_spec$PatchID),]
# Remove patches labeled as zero
assoc_spec <- assoc_spec[assoc_spec$PatchID!=0,]

# How many unique patches are filled?
length(unique(assoc_spec$PatchID))
## Surprisingly, only ~700 patches are occupied by the associated species 

# # Remove records with uncertainty greater than 20 km
# assoc_spec <- assoc_spec[!assoc_spec$CrdUncr>20000,]

## ******************************************************************** ##
## Create a data.frame of temporal trends of occupancy for that ASSOCIATED
## speices for all patches
## ******************************************************************** ##

# Initiate cumulative occupancy data.frame
cum_occ_assoc_hist <- data.frame(Pop=1:patch_total)

for ( yr in 1910:2010 ){
  # Identify populations that are occupied by year 'yr'
  pop_occ <- 
    unique(assoc_spec$PatchID[ assoc_spec$CllctnY <= yr ])
  # Make a vector of zeroes of length 'patch_total'
  pop_temp <- rep(0,patch_total)
  # Set patches found in 'pop_occ' equal to 1
  pop_temp[pop_occ] <- 1
  # Add this vector to the data.frame
  cum_occ_assoc_hist <- 
    cbind(cum_occ_assoc_hist,pop_temp)
}

# Rename column names
yr_names <- paste('YR.',1910:2010,sep='')
names(cum_occ_assoc_hist) <- c('Pop',yr_names)

## ******************************************************************** ##
## Look at cummulative occupancy pattern across **all** patches
## ******************************************************************** ##
cum_occ_assoc_hist_allPops <- colSums(cum_occ_assoc_hist)
cum_occ_assoc_hist_allPops <- cum_occ_assoc_hist_allPops[-1]
cum_occ_assoc_hist_allPops_df <- data.frame(Year=1910:2010,CumOcc=cum_occ_assoc_hist_allPops)

ggplot(cum_occ_assoc_hist_allPops_df,aes(x=Year,y=sqrt(CumOcc))) +
  geom_point() +
  geom_smooth(method='lm') +
  #geom_smooth(method = "lm", formula = y ~ poly(x, 2), colour = "red") +
  theme_bw()

## ******************************************************************** ##
## Look at cummulative occupancy patterns for both FRAL and associated spec
## ******************************************************************** ##
cum_occ_all <- rbind(cum_occ_hist_allPops_df,cum_occ_assoc_hist_allPops_df)
cum_occ_all$Spec <- c(rep('fral',nrow(cum_occ_hist_allPops_df)),
                      rep('assoc',nrow(cum_occ_assoc_hist_allPops_df)))

ggplot(cum_occ_all,aes(x=Year,y=sqrt(CumOcc),colour=Spec)) +
  geom_point() +
  theme_bw()

## ******************************************************************** ##
## Make a dataset that represents all patches that were sampled as of
## a given year. That is, patches were either historical FRAL or 
## associated species were observed.
## ******************************************************************** ##

# Start with past FRAL occurences as template
hist_pops_sampled <- cum_occupancy_hist
# Add past associated species occurences
hist_pops_sampled[2:102] <- 
  hist_pops_sampled[2:102] + cum_occ_assoc_hist[2:102]
# Threshold these values to be 1 if >0
hist_pops_sampled[2:102] <-
  ifelse(hist_pops_sampled[2:102]>0,yes=1,no=0)

# As a check, calculate the cumulative populations sampled
cum_hist_pops_sampled <- colSums(hist_pops_sampled)
cum_hist_pops_sampled <- cum_hist_pops_sampled[-1]
cum_hist_pops_sampled <- data.frame(Year=1910:2010,CumOcc=cum_hist_pops_sampled)

ggplot(cum_hist_pops_sampled,aes(x=Year,y=sqrt(CumOcc))) +
  geom_point()

## ******************************************************************** ##
## Look at cummulative occupancy patterns for FRAL, associated spec,
## and the combined occurences
## ******************************************************************** ##
cum_occ_all <- rbind(cum_occ_hist_allPops_df,
                     cum_occ_assoc_hist_allPops_df,
                     cum_hist_pops_sampled)
cum_occ_all$Spec <- c(rep('fral',nrow(cum_occ_hist_allPops_df)),
                      rep('assoc',nrow(cum_occ_assoc_hist_allPops_df)),
                      rep('comb',nrow(cum_hist_pops_sampled)))

ggplot(cum_occ_all,aes(x=Year,y=sqrt(CumOcc),colour=Spec)) +
  geom_point() +
  theme_bw()
