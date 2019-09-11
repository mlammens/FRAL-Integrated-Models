## ******************************************************************** ##
## fral_compile_results.R
##
## Author: Matthew Aiello-Lammens
## Date Created: 2013-12-21
##
## Purpose:
## Compile results from running F. alnus simulation models
## ******************************************************************** ##

## ******************************************************************** ##
## Call necessary functions and packages
## ******************************************************************** ##

## Source the fral_demog_setup.R script to setup required functions
source('scripts/fral_demog_setup.R')

## Source the fral_historical_pattern.R script to compare historical
## occurences to simulated occurences
source('scripts/fral_historical_pattern.R')

## -------------------------------------------------------------------- ##
## FUNCTION: Calculate simulation sensitivity
sens_sim <- function(x){
  length(which(x==4)) / (length(which(x==4))+length(which(x==3)))
}
## -------------------------------------------------------------------- ##
## FUNCTION: Calculate historical data sensitivity
sens_hist <- function(x){
  length(which(x==4)) / (length(which(x==4))+length(which(x==1)))
}
## -------------------------------------------------------------------- ##
## FUNCTION: Calculate OVERLAP between historical and simulation occurence
sim_hist_overlap <- function(x){
  length(which(x==4))
}
## -------------------------------------------------------------------- ##
## FUNCTION: Return historical ONLY occurences
hist_occur <- function(x){
  length(which(x==1))
}
## -------------------------------------------------------------------- ##
## FUNCTION: Return simulation ONLY occurences
sim_occur <- function(x){
  length(which(x==3))
}
## -------------------------------------------------------------------- ##
## FUNCTION: Compare simulation cumulative presence through time with 
## historical presence through time
comp_sim_hist <- function(sim_pres,hist_pres,hist_occ_only=FALSE){
  # Multiply sim_pres by 3, to give it a different index numbber than
  # the hisotical occurence data
  sim_pres[,2:101] <- sim_pres[,2:101]*3
  
  # Add sim_pres and hist_pres data_frames together.
  # In the new data.frame:
  # 1 = hist_pres only
  # 3 = sim_pres only
  # 4 = pres in both sim and hist (i.e., area of overlap)
  hist_sim_occur <- sim_pres[,2:101]+hist_pres[,2:101]
  
  ## If using only grid cells that were sampled in the past (either for 
  ## FRAL or assocaited speices), then multiply the hist_sim_occur by
  ## the hist_pops_sampled matrix
  if (hist_occ_only){
    hist_sim_occur <- hist_sim_occur * hist_pops_sampled[3:102]
  }
  
  ## Area of overlap / Simulation presences - sensitivity-simulation
  sens_sim_all <- 
    length(which(hist_sim_occur==4)) / (length(which(hist_sim_occur==4))+length(which(hist_sim_occur==3)))
  
  ## Area of overlap / Historical presences - sensitivity-historical
  sens_hist_all <-
    length(which(hist_sim_occur==4)) / (length(which(hist_sim_occur==4))+length(which(hist_sim_occur==1)))
  
  sens_list <-
    list(sens_hist_ann=apply(hist_sim_occur,MARGIN=2,FUN=sens_hist),
         sens_sim_ann=apply(hist_sim_occur,MARGIN=2,FUN=sens_sim),
         grid_overlap=apply(hist_sim_occur,MARGIN=2,FUN=sim_hist_overlap),
         hist_occur_only=apply(hist_sim_occur,MARGIN=2,FUN=hist_occur),
         sim_occur_only=apply(hist_sim_occur,MARGIN=2,FUN=sim_occur),
         sens_hist_all=sens_hist_all,
         sens_sim_all=sens_sim_all)
  
  return(sens_list)  
}

## ******************************************************************** ##
## Read MP files:
## --------------
## Search for files ending in *.SCL, as these indicate MP simulations
## that have been completed.
## ******************************************************************** ##

# Get list of completed fral simulations 
# fral_files <- list.files(path='ramas/',pattern='fral_.*SCL',full.names=TRUE)
fral_files <- list.files(path='~/Dropbox/F-alnus/RAMAS_Models_Storage/',
                         pattern='fral_.*SCL',full.names=TRUE)
fral_files <- sub(pattern='SCL',replacement='mp',x=fral_files)

## Seperate additional (i.e., 'add') files out to maintain order
## ***
# Get the additional files
fral_files_add <- fral_files[grepl(pattern='add',fral_files)]
# Seperate out the popd ones
fral_files_add_popd_LDD <- fral_files_add[grepl(pattern='popd',fral_files_add)]
fral_files_add <- fral_files_add[!grepl(pattern='popd',fral_files_add)]
# Remove the 'add' sims from the fral_files list
fral_files <- fral_files[!grepl(pattern='add',fral_files)]
# Re-combine file names
fral_files <- c(fral_files,fral_files_add,fral_files_add_popd_LDD)


run_extraction <- FALSE

## Set the population size threshold. This number is the populations
## size for a patch to be consdidered occupied
## Thresholds used: 0 , 1000, 2000
pop_size_thresh <- 1000

if (run_extraction) {
  ## Setup a cluster to run the results extraction and compilation
  require(foreach)
  require(doSNOW)
  cl <- makeCluster(7,"SOCK") ### WARNING: Change number of cores
  registerDoSNOW(cl)
  
  # Monitor time to complete loop
  proc_time_start <- proc.time()
  
  # Reset 'list' object to store extracted results
  cum_occupancy_mpMult <- c()
  
  cum_occupancy_mpMult <- foreach  (f_ind = 1:length(fral_files),.packages=c('lhs','fields','popbio')) %dopar% { # length(fral_files)
    mp <- mp.read.results(mpFile=fral_files[f_ind])
    mp <- mp$mp.file
    mp_res <- mp$results
    
    ## Extract simulation results metrics
    fral_results <- mp.results(mpFile=fral_files[f_ind],mpList=mp,spatial=TRUE,mac=TRUE)
    # Add the mp.file name
    fral_results$mp.file <- basename(row.names(fral_results))
    # Add the first kch file name
    fral_results$kch.file <- mp$PopData_df$KchangeSt[1]
    # Add the number of translocations to the fral_results
    fral_results$NumTranslocs <- nrow(mp$PopManageProp)
    
    ## Extract the stage matrix
    fral_StMatr <- mp$StMatr[[1]]$Matr
    
    ## -------------------------------------------------------------------- ##
    ## Make a cumulative occupancy data.frame for the mp results
    
    # Get population size through time
    pop_size_time <- mp_res$PopInd[,1,]
        
    # Get the cumulative maximum population size through time
    pop_size_time_cumMax <- t(apply(pop_size_time,MARGIN=2,FUN=cummax))
    
    # Threshold population size to make binary matrix
    pop_pres_time <- ifelse(pop_size_time_cumMax>pop_size_thresh,yes=1,no=0) ## Thresholds used: 0 , 1000, 2000
    pop_pres_time <- as.data.frame(pop_pres_time)
    
    pop_pres_time <- cbind(1:3423,pop_pres_time)
    names(pop_pres_time) <- c('Pop',paste('YR.',1911:2010,sep=''))
    
    ## -------------------------------------------------------------------- ##
    ## Look at cummulative occupancy pattern across **all** patches
    
    cum_occupancy_allPops <- colSums(pop_pres_time)
    cum_occupancy_allPops <- cum_occupancy_allPops[-1]
    
    ## -------------------------------------------------------------------- ##
    ## Look at cummulative occupancy pattern across patches **sampled** in 
    ## a give year
        
    pop_pres_time_histSamp <- pop_pres_time[2:101] * hist_pops_sampled[3:102]
    cum_occupancy_allPops_histSamp <- colSums(pop_pres_time_histSamp)

    ## -------------------------------------------------------------------- ##
    ## Return a list object with extracted results
    list(MP_file=fral_files[f_ind],
         Patch_Pres=pop_pres_time,
         Cum_Occ=cum_occupancy_allPops,
         Cum_Occ_histSamp=cum_occupancy_allPops_histSamp,
         MP_res=fral_results,
         StMatr=fral_StMatr)
  }
  
  # Stop cluster
  stopCluster(cl)
  
  # Determine the amount of time taken to run loop
  proc_time_finish <- proc.time()
  print(proc_time_finish - proc_time_start)
  
  ## Save the results read in from the mp files incase of RStudio crash
  #save(cum_occupancy_mpMult,file='ramas/cum_occupancy_mpMult.RData')
  #save(cum_occupancy_mpMult,file='ramas/cum_occupancy_mpMult_1000.RData')
  #save(cum_occupancy_mpMult,file='ramas/cum_occupancy_mpMult_2000.RData')
  ## Save the 500 additional simulation results
  #save(cum_occupancy_mpMult,file='ramas/cum_occupancy_mpMult_add.RData')
  #save(cum_occupancy_mpMult,file='ramas/cum_occupancy_mpMult_add_1000.RData')
  #save(cum_occupancy_mpMult,file='ramas/cum_occupancy_mpMult_add_2000.RData')
  ## Save the 500 additional simulation results that used popd migration
  #save(cum_occupancy_mpMult,file='ramas/cum_occupancy_mpMult_add_popd_1.RData',compress=TRUE)
  #save(cum_occupancy_mpMult,file='ramas/cum_occupancy_mpMult_add_popd_1000.RData',compress=TRUE)
  #save(cum_occupancy_mpMult,file='ramas/cum_occupancy_mpMult_add_popd_2000.RData',compress=TRUE)
  
}

## ******************************************************************** ##
## Extract final population size informaiton from the simulations
## ******************************************************************** ##

extract_pop_sizes <- FALSE

if (extract_pop_sizes){
  ## Setup a cluster to run the results extraction and compilation
  require(foreach)
  require(doSNOW)
  cl <- makeCluster(3,"SOCK") ### WARNING: Change number of cores
  registerDoSNOW(cl)
  
  # Monitor time to complete loop
  proc_time_start <- proc.time()
  
  mp_pop_size <- c()
  
  pop_sizes_mpMult <- foreach  (f_ind = 1:length(fral_files),
                                .packages=c('lhs','fields','popbio'),
                                .combine=rbind) %dopar% { # length(fral_files)
    mp <- mp.read.results(mpFile=fral_files[f_ind])
    mp <- mp$mp.file
    mp_res <- mp$results

    # Get population size through time
    pop_size_time <- mp_res$PopInd[,1,]
    
    # Get the mp.file name
    mp_pop_size$mp.file <- basename(fral_files[f_ind])
    
    # Get mean and SD population size for the last time step
    # -- add these to the fral_results data.frame
    mp_pop_size$pop_size_t100_mean <- mean( pop_size_time[100,] )
    mp_pop_size$pop_size_t100_sd <- sd( pop_size_time[100,] )
    
    mp_pop_size
    
                                }
    
  # Stop cluster
  stopCluster(cl)
  
  # Determine the amount of time taken to run loop
  proc_time_finish <- proc.time()
  print(proc_time_finish - proc_time_start)
  
  pop_sizes_mpMult_df <- data.frame( mp.file=as.character(pop_sizes_mpMult[,1]),
                                     pop_size_t100_mean=as.numeric(pop_sizes_mpMult[,2]),
                                     pop_size_t100_sd=as.numeric(pop_sizes_mpMult[,3])
  )
  
  # Write this information to file
  write.csv(pop_sizes_mpMult_df,'outputs/pop_sizes_mpMult.csv',row.names=FALSE,quote=FALSE)
  
}

## ******************************************************************** ##
## Read in saved cum_occupmancy_mpMult files
## ******************************************************************** ##

## Load saved 'cum_occupancy_mpMult' as needed
# load('~/Dropbox/F-Alnus/Chapter-5_ramas_folder/ramas/cum_occupancy_mpMult_all_1.RData')
load('~/Dropbox/F-Alnus/Chapter-5_ramas_folder/ramas/cum_occupancy_mpMult_all_1000.RData')
# load('~/Dropbox/F-Alnus/Chapter-5_ramas_folder/ramas/cum_occupancy_mpMult_all_2000.RData')

## ******************************************************************** ##
## Caclulate population occurence overlap
## ******************************************************************** ##
## Compare occurence patterns of simulation with historical occurences
## In order to make this comparison, first drop occurences from 1910, since
## this comparison is not meaningful (the simulation occurence is set to be
## equal to the observed occurence record for 1910).
pop_hist_time <- cum_occupancy_hist
pop_hist_time$YR.1910 <- NULL

## ******************************************************************** ##
## Compile simulation results
## ******************************************************************** ##

## -------------------------------------------------------------------- ##
## Make a new data.frame of only the 'Patch_Pres' results
Patch_Pres_Sim <- 
  lapply(cum_occupancy_mpMult,function(x){x$Patch_Pres}) 

## Get senstivity values for historical and simulation occurences
Hist_Sim_Sens <- 
  #lapply(Patch_Pres_Sim,function(x){comp_sim_hist(sim_pres=x,hist_pres=pop_hist_time)})
  mclapply(Patch_Pres_Sim,function(x){comp_sim_hist(sim_pres=x,hist_pres=pop_hist_time)})

## Get senstivity values for historical and simulation occurences, for
## grid cells that were observed in the past
Hist_Sim_Sens_ObsGrid <- 
  #lapply(Patch_Pres_Sim,function(x){comp_sim_hist(sim_pres=x,hist_pres=pop_hist_time,hist_occ_only=TRUE)})
  mclapply(Patch_Pres_Sim,function(x){comp_sim_hist(sim_pres=x,hist_pres=pop_hist_time,hist_occ_only=TRUE)})

## -------------------------------------------------------------------- ##
## For each simulation, calculate the mean Historical and 
## Simulation Sensitivity
Hist_Sim_Sens_All_Mean <-
  ldply(lapply(Hist_Sim_Sens,function(x){data.frame(hist_sens=mean(x$sens_hist_ann),
                                                    sim_sens=mean(x$sens_sim_ann) )}))

Hist_Sim_Sens_All_Mean_ObsGrid <-
  ldply(lapply(Hist_Sim_Sens_ObsGrid,function(x){data.frame(hist_sens_obsgrd=mean(x$sens_hist_ann),
                                                            sim_sens_obsgrd=mean(x$sens_sim_ann) )}))

## -------------------------------------------------------------------- ##
## Get the simulation results for all simulations as a data.frame

Sim_Results <- 
  ldply(lapply(cum_occupancy_mpMult,function(x){x$MP_res}))

## -------------------------------------------------------------------- ##
## Calculate mean fecundity and mean survival values
Sim_Results$Fec.Mean <- 
  ldply(lapply(cum_occupancy_mpMult,function(x){mean(x$StMatr[1,5:50])}))$V1

Sim_Results$Surv.Mean <- 
  #ldply(lapply(cum_occupancy_mpMult,function(x){mean(x$StMatr[2:50,])}))$V1
  ldply(lapply(cum_occupancy_mpMult,
               function(x){ mean( colSums( x$StMatr[2:50,] ) ) }
               ) )$V1

## -------------------------------------------------------------------- ##
## Add the Hist_Sim_Sens_All_Mean values to this data.frame
Sim_Results <- cbind(Sim_Results,
                     Hist_Sim_Sens_All_Mean, 
                     Hist_Sim_Sens_All_Mean_ObsGrid)

## ******************************************************************** ##
## Calculate sensitivity and extent values for historical occurence
## and simulations For all years.
## ******************************************************************** ##

# Calculate year 1 values
Hist_Sim_Sens_All <-
  ldply(lapply(Hist_Sim_Sens,function(x){data.frame(Sim_Sens=x$sens_sim_ann[1],
                                                    Hist_Sens=x$sens_hist_ann[1],
                                                    Grid_Overlap=x$grid_overlap[1],
                                                    Hist_Occur=x$hist_occur_only[1],
                                                    Sim_Occur=x$sim_occur_only[1])}))

# Make vector of years - based on timepoints in simulation, NOT actual years
years <- seq(from=2,to=100)

for ( y_ind in 1:length(years) ){
  Hist_Sim_Sens_TEMP <-
    ldply(lapply(Hist_Sim_Sens,function(x){data.frame(Sim_Sens=x$sens_sim_ann[years[y_ind]],
                                                      Hist_Sens=x$sens_hist_ann[years[y_ind]],
                                                      Grid_Overlap=x$grid_overlap[years[y_ind]],
                                                      Hist_Occur=x$hist_occur_only[years[y_ind]],
                                                      Sim_Occur=x$sim_occur_only[years[y_ind]])}))
  Hist_Sim_Sens_All <- data.frame(Hist_Sim_Sens_All,Hist_Sim_Sens_TEMP)
}

# Add the names for these new variables
hist_sim_sens_Names <- 
  paste( rep(c('sens_sim','sens_hist','grid_overlap','hist_only','sim_only'),length(years)+1),
         rep(seq(1911,2010),each=5), 
         sep='.' )

names(Hist_Sim_Sens_All) <- hist_sim_sens_Names

# Add the simulation name
Hist_Sim_Sens_All <-
  cbind(Sim_Results$mp.file,Hist_Sim_Sens_All)
names(Hist_Sim_Sens_All)[1] <- 'mp.file'

## DEPRECIATED SECTION - THESE PLOTS NOW CONSTRUCTED IN fral_make_cumulative_AOO_figures.R
# # Plot trends through time
# Hist_Sim_Sens_All_m <-
#   melt(Hist_Sim_Sens_All)
# # Add years column
# Hist_Sim_Sens_All_m$Year <- as.numeric(sub(pattern='.*\\.',replacement='',Hist_Sim_Sens_All_m$variable))
# 
# ggplot(Hist_Sim_Sens_All_m[grepl(pattern='sens_hist',Hist_Sim_Sens_All_m$variable),],
#        aes(x=Year,y=value,colour=mp.file)) +
#   geom_line() +
#   ylab('Sensitivity') +
#   guides(colour=FALSE) +
#   theme_bw()

# ## Plot only popd dense effect simulations
# Hist_Sim_Sens_Popd_Denseff <- Hist_Sim_Sens_All[ 2501:3000, ]
# 
# # Plot trends through time
# Hist_Sim_Sens_Popd_Denseff_m <-
#   melt(Hist_Sim_Sens_Popd_Denseff)
# # Add years column
# Hist_Sim_Sens_Popd_Denseff_m$Year <- 
#   as.numeric(sub(pattern='.*\\.',replacement='',Hist_Sim_Sens_Popd_Denseff_m$variable))
# 
# ggplot(Hist_Sim_Sens_Popd_Denseff_m[grepl(pattern='sens_hist',Hist_Sim_Sens_Popd_Denseff_m$variable),],
#        aes(x=Year,y=value,colour=mp.file)) +
#   geom_line() +
#   ylab('Sensitivity') +
#   guides(colour=FALSE) +
#   geom_hline(yintercept=0.50) +
#   theme_bw()
# ggsave(filename='figures/Sensitivity_PopD_DensEff_Models.pdf',width=7.5,height=7.5,units='in')

## ******************************************************************** ##
## Calculate sensitivity and extent values for historical occurence
## and simulations For all years, using observations in historical 
## occerence grids only.
## ******************************************************************** ##

# Calculate year 1 values
Hist_Sim_Sens_ObsGrid_All <-
  ldply(lapply(Hist_Sim_Sens_ObsGrid,function(x){data.frame(Sim_Sens=x$sens_sim_ann[1],
                                                            Hist_Sens=x$sens_hist_ann[1],
                                                            Grid_Overlap=x$grid_overlap[1],
                                                            Hist_Occur=x$hist_occur_only[1],
                                                            Sim_Occur=x$sim_occur_only[1])}))

# Make vector of years - based on timepoints in simulation, NOT actual years
years <- seq(from=2,to=100)

for ( y_ind in 1:length(years) ){
  Hist_Sim_Sens_TEMP <-
    ldply(lapply(Hist_Sim_Sens_ObsGrid,function(x){data.frame(Sim_Sens=x$sens_sim_ann[years[y_ind]],
                                                      Hist_Sens=x$sens_hist_ann[years[y_ind]],
                                                      Grid_Overlap=x$grid_overlap[years[y_ind]],
                                                      Hist_Occur=x$hist_occur_only[years[y_ind]],
                                                      Sim_Occur=x$sim_occur_only[years[y_ind]])}))
  Hist_Sim_Sens_ObsGrid_All <- data.frame(Hist_Sim_Sens_ObsGrid_All,Hist_Sim_Sens_TEMP)
}

# Add the names for these new variables
hist_sim_sens_Names <- 
  paste( rep(c('sens_sim','sens_hist','grid_overlap','hist_only','sim_only'),length(years)+1),
         rep(seq(1911,2010),each=5), 
         sep='.' )

names(Hist_Sim_Sens_ObsGrid_All) <- hist_sim_sens_Names

# Add the simulation name
Hist_Sim_Sens_ObsGrid_All <-
  cbind(Sim_Results$mp.file,Hist_Sim_Sens_ObsGrid_All)
names(Hist_Sim_Sens_ObsGrid_All)[1] <- 'mp.file'


## ******************************************************************** ##
## ******************************************************************** ##
## Look at cumulative area of occupancy curves
## ******************************************************************** ##
## ******************************************************************** ##

cum_occ_mpMult_df <- make_cum_occ_df( cum_occ_list=cum_occupancy_mpMult, histSamp=TRUE )
cum_occ_mpMult_df_ALL <- make_cum_occ_df( cum_occ_list=cum_occupancy_mpMult, histSamp=FALSE )

# Grab the historical pattern and remove the first line
hist_temp <- cum_occ_hist_allPops_df[-1,]

## -------------------------------------------------------------------- ##
## Calculate a loss function comparing the two cumulative distribution
## curves

## Melt the cum_occ_mpMult_df data.frame
cum_occ_mpMult_df_m <- melt(cum_occ_mpMult_df,id.vars='Year')

## Calculate the Loss for each model
cum_AOO_comp_Loss <- apply(cum_occ_mpMult_df[2:ncol(cum_occ_mpMult_df)],
                           MARGIN=2,
                           FUN=function(col){sum(abs(hist_temp$CumOcc-col))})
hist(cum_AOO_comp_Loss,breaks=40)

## Make a data.frame of loss functoin values and model names
cum_AOO_comp_Loss_df <- data.frame(mp_mod=colnames(cum_occ_mpMult_df[-1]),
                                   loss=cum_AOO_comp_Loss)

## Add the value of the loss function as an indicator to the melted data.frame
cum_occ_mpMult_df_m <- merge(cum_occ_mpMult_df_m,cum_AOO_comp_Loss_df,
                             by.x='variable',by.y='mp_mod')

## -------------------------------------------------------------------- ##
## Plot the cumulative area of occupancy curves
ggplot() + 
  geom_line(data=cum_occ_mpMult_df_m,
            aes(x=Year,y=sqrt(value), group=variable,colour=loss)) +
  #scale_colour_gradient(limits=c(2500,23000),low="blue",high="red") +
  scale_colour_gradient(limits=c(0,3000),low="blue",high="red") +
  ylab('Sqrt(Cumulative 20x20 km Grid Cells)') +
  geom_line(data=hist_temp,aes(x=Year,y=sqrt(CumOcc)),size=2) +
  theme_bw()

ggplot() + 
  geom_line(data=subset(cum_occ_mpMult_df_m,subset=loss<=5000),
            aes(x=Year,y=sqrt(value), group=variable,colour=loss)) +
  #scale_colour_gradient(limits=c(2500,23000),low="blue",high="red") +
  scale_colour_gradient(limits=c(0,3000),low="blue",high="red") +
  ylab('Sqrt(Cumulative 20x20 km Grid Cells)') +
  geom_line(data=hist_temp,aes(x=Year,y=sqrt(CumOcc)),size=2) +
  theme_bw()

#ggsave(filename='figures/Cumulative_AOO_Curves_LOSS.pdf',width=14,height=14,units='in')

## ******************************************************************** ##
## Add values to Sim_Results
## ******************************************************************** ##

## -------------------------------------------------------------------- ##
## Make a KCH type variable
Sim_Results$kch.type <- NA
# Low values
Sim_Results$kch.type[ grepl(pattern='LOW',Sim_Results$kch.file) ] <- 'LOW'
Sim_Results$kch.type[ grepl(pattern='MED',Sim_Results$kch.file) ] <- 'MED'
Sim_Results$kch.type[ grepl(pattern='HIGH',Sim_Results$kch.file) ] <- 'HIGH'
Sim_Results$kch.type <- as.factor(Sim_Results$kch.type)

## -------------------------------------------------------------------- ##
## Make a variable indicating whether there was habitat suitability 
## change through time or not
Sim_Results$hs.chng <- TRUE
Sim_Results$hs.chng[ grepl(pattern='nochng',Sim_Results$mp.file) ] <- FALSE

## -------------------------------------------------------------------- ##
## Change the variable indicating whether effective density or ceiling
## type density dependence are used - just to make the factors more
## understandable

## First check the current levels and make sure they are correct w/r
## to model runs
if (!all( levels(Sim_Results$dd.type) == c("CE", "UD Frangula\\FrangulaEDss.dll") )){
  print('!!!WARNING!!! Levels for density dependence type appear incorrect')
}

## Next rename the levels
Sim_Results$dd.type <-
  revalue( Sim_Results$dd.type, c("CE"="ceil","UD Frangula\\FrangulaEDss.dll"="eff_dens") )

## -------------------------------------------------------------------- ##
## Add sensitivity for key years
Sim_Results$sens_hist.1950 <- Hist_Sim_Sens_All$sens_hist.1950
Sim_Results$sens_hist.2010 <- Hist_Sim_Sens_All$sens_hist.2010

## Add simulation sensitivity (positive predictive power)
Sim_Results$sens_sim_obsgrd.1950 <- Hist_Sim_Sens_ObsGrid_All$sens_sim.1950
Sim_Results$sens_sim_obsgrd.2010 <- Hist_Sim_Sens_ObsGrid_All$sens_sim.2010

## -------------------------------------------------------------------- ##
## Add loss function values
Sim_Results$cum_AOO_loss <- cum_AOO_comp_Loss

## -------------------------------------------------------------------- ##
## Add the number of populations that are occupied by the end of the 
## simulation
Sim_Results$occupied_pops_all <-
  as.numeric( cum_occ_mpMult_df_ALL[100, 2:ncol(cum_occ_mpMult_df_ALL)] ) 
Sim_Results$occupied_pops_obsgrd <-
  as.numeric( cum_occ_mpMult_df[100, 2:ncol(cum_occ_mpMult_df)] ) 

## -------------------------------------------------------------------- ##
## Save Sim_Results data.frame
# write.csv(Sim_Results, 'outputs/sim_results_1.csv',row.names=FALSE,quote=FALSE)
# write.csv(Sim_Results, 'outputs/sim_results_1000.csv',row.names=FALSE,quote=FALSE)
# write.csv(Sim_Results, 'outputs/sim_results_2000.csv',row.names=FALSE,quote=FALSE)

## -------------------------------------------------------------------- ##
## Save Hist_Sim_Sens* grids
# write.csv(Hist_Sim_Sens_All,'outputs/hist_sim_sens_all_1.csv',row.names=FALSE,quote=FALSE)
# write.csv(Hist_Sim_Sens_ObsGrid_All,'outputs/hist_sim_sens_obsgrd_1.csv',row.names=FALSE,quote=FALSE)
# write.csv(Hist_Sim_Sens_All,'outputs/hist_sim_sens_all_1000.csv',row.names=FALSE,quote=FALSE)
# write.csv(Hist_Sim_Sens_ObsGrid_All,'outputs/hist_sim_sens_obsgrd_1000.csv',row.names=FALSE,quote=FALSE)
# write.csv(Hist_Sim_Sens_All,'outputs/hist_sim_sens_all_2000.csv',row.names=FALSE,quote=FALSE)
# write.csv(Hist_Sim_Sens_ObsGrid_All,'outputs/hist_sim_sens_obsgrd_2000.csv',row.names=FALSE,quote=FALSE)

## -------------------------------------------------------------------- ##
## Save the cumulative occupancy through time data.frames
# write.csv(cum_occ_mpMult_df_ALL,'outputs/cumulative_occupancy_all_1.csv',row.names=FALSE,quote=FALSE)
# write.csv(cum_occ_mpMult_df,'outputs/cumulative_occupancy_obsgrid_1.csv',row.names=FALSE,quote=FALSE)
# write.csv(cum_occ_mpMult_df_ALL,'outputs/cumulative_occupancy_all_1000.csv',row.names=FALSE,quote=FALSE)
# write.csv(cum_occ_mpMult_df,'outputs/cumulative_occupancy_obsgrid_1000.csv',row.names=FALSE,quote=FALSE)
# write.csv(cum_occ_mpMult_df_ALL,'outputs/cumulative_occupancy_all_2000.csv',row.names=FALSE,quote=FALSE)
# write.csv(cum_occ_mpMult_df,'outputs/cumulative_occupancy_obsgrid_2000.csv',row.names=FALSE,quote=FALSE)

## -------------------------------------------------------------------- ##
## Save workspace image
#save.image('R_Session_BKUP_1000.RData',compress=TRUE)
#save.image('R_Session_BKUP_1.RData',compress=TRUE)
#save.image('R_Session_BKUP_2000.RData',compress=TRUE)