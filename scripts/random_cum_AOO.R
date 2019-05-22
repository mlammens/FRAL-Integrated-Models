## ******************************************************************** ##
## random_cum_AOO.R
##
## Author: Matthew Aiello-Lammens
## Date Created: 2014-02-08
##
## Purpose:
## Generate random cumulative AOO curves and calculate the loss function
## between these curves and the historical FRAL cum AOO
##
## ******************************************************************** ##

## ******************************************************************** ##
## Call necessary functions and packages
## ******************************************************************** ##

## Source the fral_demog_setup.R script to setup required functions
source('R/fral_demog_setup.R')

## Source the fral_historical_pattern.R script to compare historical
## occurences to simulated occurences
source('R/fral_historical_pattern.R')

## -------------------------------------------------------------------- ##
## FUNCTION: make_cum_aoo
## ***
## Generate a random cumulative AOO curve
## -------------------------------------------------------------------- ##
make_cum_aoo <- function() {
  # Set up time step one
  #cum_aoo_temp <- sample( 0:cum_hist_pops_sampled$CumOcc[1], size=1 )
  cum_aoo_temp <- 11 # In the simulation I make, all start at 11
  for ( x in 2:100){
    #cum_ts <- sample( cum_aoo_temp[x-1]:cum_hist_pops_sampled$CumOcc[x], size=1 )
    cum_ts <- sample( 0:cum_hist_pops_sampled$CumOcc[x], size=1 )
    cum_ts <- max( cum_ts, cum_aoo_temp[x-1] )
    cum_aoo_temp <- c(cum_aoo_temp, cum_ts)
  }
  return( cum_aoo_temp )
}
## -------------------------------------------------------------------- ##


## ******************************************************************** ##
## Make a data set of simulations
## ******************************************************************** ##

# Monitor time to complete loop
proc_time_start <- proc.time()

## Set the number of simulations
sims_num <- 100000
## Make simulation names
sim_names <- paste( 'sim_', 1:sims_num, sep='' )
## Use the make_cum_aoo function
sim_aoo <- c()
for ( x in 1:sims_num ){
  sim_aoo <- cbind( sim_aoo, make_cum_aoo() )
}
## Assign names to data.frame
sim_aoo <- as.data.frame(sim_aoo)
names(sim_aoo) <- sim_names

# Determine the amount of time taken to run loop
proc_time_finish <- proc.time()
print(proc_time_finish - proc_time_start)

## Make a histogram of the final cumulative AOO
hist(as.numeric(sim_aoo[100,]))
