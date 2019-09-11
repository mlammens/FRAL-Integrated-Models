## ******************************************************************** ##
## fral_demog_setup.R
##
## Author: Matthew Aiello-Lammens
## Date Created: 2013-12-22
##
## Purpose:
## Load required packages and source user defined functions for F. alnus
## demographic analysis (Chapter 4 of Aiello-Lammens PhD Dissertation)
## ******************************************************************** ##

## Required packages and user defined functions
require(raster)
require(ggplot2)
require(reshape2)
require(plyr)
require(lhs)
require(fields)
require(GGally)
require(gridExtra)

# if(tolower(Sys.info()['sysname'])!='windows'){
#   library("multicore")
# }

## Call packages for BRT analysis
require(gbm)
require(dismo)

## External functions (i.e., Sensitivity Analysis functions)

# # Get node name of computer
# system_nodename <- tolower(Sys.info()['nodename'])
# 
# # Change directories based on system node name
# if (system_nodename=='chloe-pc'){
#   source('/Users/plover/Dropbox/RA-Sensitivity/SACode_Sandbox/sensitivity.setup.r')
#   sensitivity.setup('/Users/plover/Dropbox/RA-Sensitivity/SACode_Sandbox/')
# } else if (system_nodename=='jessie-labpc') {
#   ## Directories to use on Jessies-Lab computer
#   source('/Users/jessie/Dropbox/RA-Sensitivity/SACode_Sandbox/sensitivity.setup.r')
#   sensitivity.setup('/Users/jessie/Dropbox/RA-Sensitivity/SACode_Sandbox/')
# } else if (system_nodename=='akcakaya-lab') {
#   source('/Users/Matthew/My Documents/My Dropbox/RA-Sensitivity/SACode_Sandbox/sensitivity.setup.r')
#   sensitivity.setup('/Users/Matthew/My Documents/My Dropbox/RA-Sensitivity/SACode_Sandbox/')
# } else if (system_nodename=='protea') {
#   source('/home/punctata/Dropbox/RA-Sensitivity/SACode_Sandbox/sensitivity.setup.r')
#   sensitivity.setup('/home/punctata/Dropbox/RA-Sensitivity/SACode_Sandbox/')
# } else { # Assume on MA-L MacBook Pro
#   source('~/Dropbox/RA-Sensitivity/SACode_Sandbox/sensitivity.setup.r')
#   sensitivity.setup('~/Dropbox/RA-Sensitivity/SACode_Sandbox/')
# }

source('~/Dropbox/RA-Sensitivity/SACode_Sandbox/sensitivity.setup.r')
sensitivity.setup('~/Dropbox/RA-Sensitivity/SACode_Sandbox/')


## ******************************************************************** ##
## Custom functions
## ******************************************************************** ##

## -------------------------------------------------------------------- ##
## FUNCTION: make_cum_occ_df
## Short Discription: Calculate cumulative occupancy data.frames
## 
## Args:
##   cum_occ_list: List object with cumulative occupancy information for
##     each simulation
##
##   histSamp: (TRUE/FALSE) Use only patches that were sampled historically
##
## Returns:
##   cum_occ_df: Cumulative occupancy data.frame
##
## -------------------------------------------
## Make a data.frame of the cumulative area of occupancy curves for each 
## simulation. These curves can then be plots along with the pattern from the 
## historic data (sensu Chapter 3)
## -------------------------------------------------------------------- ##

make_cum_occ_df <- function( cum_occ_list, histSamp=TRUE ){
  # Get data from list object
  if ( histSamp ){
    cum_occ_df <- ldply(lapply(cum_occ_list,function(x){x$Cum_Occ_histSamp}))
  } else if ( !histSamp ){
    cum_occ_df <- ldply(lapply(cum_occ_list,function(x){x$Cum_Occ}))
  } else {
    stop("ERROR: Incorrect setting for argument 'histSamp'")
  }
  # Reformat the data.frame and add years
  cum_occ_df <- as.data.frame(t(cum_occ_df))
  cum_occ_df <- cbind(1911:2010,cum_occ_df)
  # Get the names of the simulations
  cum_occ_sim_names <- ldply( lapply( cum_occ_list, function(x){ basename(x$MP_file) } ) )
  cum_occ_sim_names <- cum_occ_sim_names$V1
  # Add names to make column names for final data.frame
  names(cum_occ_df) <- c('Year',cum_occ_sim_names)
  # Return cum_occ_df
  return( cum_occ_df )
}

## -------------------------------------------------------------------- ##
## FUNCTION: calc_loss
## Short Discription: Calculate the LOSS values between historical FRAL
##   sampling and simulated spread
## 
## Args:
##   cum_occ_df: data.frame object with cumulative occupancy information for
##     each simulation. Each column represents a unique simulation.
##
##   hist_temp: data.frame object containing historical cumulative AOO 
##     pattern
##
## Returns:
##   cum_occ_df_m: Melted cumulative AOO data.frame, including loss values
##
## -------------------------------------------
## Melt the cum_occ_df and compare with the historical cumulativa area of
## occurence curve to calculate loss
## -------------------------------------------------------------------- ##

calc_loss <- function( cum_occ_df, hist_temp ){
  ## Melt the cum_occ_mpMult_df data.frame
  cum_occ_df_m <- melt(cum_occ_df,id.vars='Year')
  
  ## Calculate the Loss for each model
  cum_AOO_comp_Loss <- apply(cum_occ_df[2:ncol(cum_occ_df)],
                             MARGIN=2,
                             FUN=function(col){sum(abs(hist_temp$CumOcc-col))})
  #hist(cum_AOO_comp_Loss,breaks=40)
  
  ## Make a data.frame of loss functoin values and model names
  cum_AOO_comp_Loss_df <- data.frame(mp_mod=colnames(cum_occ_df[-1]),
                                     loss=cum_AOO_comp_Loss)
  
  ## Add the value of the loss function as an indicator to the melted data.frame
  cum_occ_df_m <- merge(cum_occ_df_m,cum_AOO_comp_Loss_df,
                               by.x='variable',by.y='mp_mod')
  return( cum_occ_df_m )

}

## -------------------------------------------------------------------- ##
## FUNCTION: make_gbm_interaction_plots
## ***
## Make four interaction plots with labels
make_gbm_interaction_plots <- function( gbm_model, zrange, fig_file){
  # Plot the main interaction between Fecundity and Number of Translocations
  graphics.off()
  pdf(file=fig_file,width=7.5,height=7.5)
  par(mfrow=c(2,2),mar=c(1,1,1,1))
  gbm.perspec(gbm_model,x=7,y=6,z.range=zrange,
              x.label='Mean fecundity',
              y.label='Numer of LDD events',
              z.label='Mean annual sensitivity')
  # Look at other interactions
  gbm.perspec(gbm_model,x=7,y=4,z.range=zrange,
              x.label='Mean fecundity',
              y.label='Mean dispersal',
              z.label='')
  gbm.perspec(gbm_model,x=8,y=6,z.range=zrange,
              x.label='Mean survival-growth',
              y.label='Numer of LDD events',
              z.label='Mean annual sensitivity')
  gbm.perspec(gbm_model,x=7,y=8,z.range=zrange,
              x.label='Mean fecundity',
              y.label='Mean survival-growth',
              z.label='')
  dev.off()
}
## -------------------------------------------------------------------- ##


## -------------------------------------------------------------------- ##
## FUNCTION: gbm.plot.cust
##
## Modified gbm.plot function 
gbm.plot.cust <- function (gbm.object, variable.no = 0, smooth = FALSE, rug = TRUE, 
                           n.plots = length(pred.names), common.scale = TRUE, write.title = TRUE, 
                           y.label = "fitted function", x.label = NULL, show.contrib = TRUE, 
                           plot.layout = c(3, 4), ...) 
{
  if (!require(gbm)) {
    stop("you need to install the gbm package to run this function")
  }
  if (!require(splines)) {
    stop("you need to install the splines package to run this function")
  }
  gbm.call <- gbm.object$gbm.call
  gbm.x <- gbm.call$gbm.x
  pred.names <- gbm.call$predictor.names
  response.name <- gbm.call$response.name
  dataframe.name <- gbm.call$dataframe
  data <- eval(parse(text = dataframe.name))
  max.plots <- plot.layout[1] * plot.layout[2]
  plot.count <- 0
  n.pages <- 1
  if (length(variable.no) > 1) {
    stop("only one response variable can be plotted at a time")
  }
  if (variable.no > 0) {
    n.plots <- 1
  }
  max.vars <- length(gbm.object$contributions$var)
  if (n.plots > max.vars) {
    n.plots <- max.vars
    warning("reducing no of plotted predictors to maximum available (", 
            max.vars, ")")
  }
  predictors <- list(rep(NA, n.plots))
  responses <- list(rep(NA, n.plots))
  for (j in c(1:n.plots)) {
    if (n.plots == 1) {
      k <- variable.no
    }
    else k <- match(gbm.object$contributions$var[j], pred.names)
    if (is.null(x.label)) 
      var.name <- gbm.call$predictor.names[k]
    else var.name <- x.label
    pred.data <- data[, gbm.call$gbm.x[k]]
    response.matrix <- plot.gbm(gbm.object, k, return.grid = TRUE)
    predictors[[j]] <- response.matrix[, 1]
    if (is.factor(data[, gbm.call$gbm.x[k]])) {
      predictors[[j]] <- factor(predictors[[j]], levels = levels(data[, 
                                                                      gbm.call$gbm.x[k]]))
    }
    responses[[j]] <- response.matrix[, 2] #- mean(response.matrix[, 
                                            #                      2])
    if (j == 1) {
      ymin = min(responses[[j]])
      ymax = max(responses[[j]])
    }
    else {
      ymin = min(ymin, min(responses[[j]]))
      ymax = max(ymax, max(responses[[j]]))
    }
  }
  op <- par(no.readonly = TRUE)
  par(mfrow = plot.layout)
  for (j in c(1:n.plots)) {
    if (plot.count == max.plots) {
      plot.count = 0
      n.pages <- n.pages + 1
    }
    plot.count <- plot.count + 1
    if (n.plots == 1) {
      k <- match(pred.names[variable.no], gbm.object$contributions$var)
      if (show.contrib) {
        x.label <- paste(var.name, "  (", round(gbm.object$contributions[k, 
                                                                         2], 1), "%)", sep = "")
      }
    }
    else {
      k <- match(gbm.object$contributions$var[j], pred.names)
      var.name <- gbm.call$predictor.names[k]
      if (show.contrib) {
        x.label <- paste(var.name, "  (", round(gbm.object$contributions[j, 
                                                                         2], 1), "%)", sep = "")
      }
      else x.label <- var.name
    }
    if (common.scale) {
      plot(predictors[[j]], responses[[j]], ylim = c(ymin, 
                                                     ymax), type = "l", xlab = x.label, ylab = y.label, 
           ...)
    }
    else {
      plot(predictors[[j]], responses[[j]], type = "l", 
           xlab = x.label, ylab = y.label, ...)
    }
    if (smooth & is.vector(predictors[[j]])) {
      temp.lo <- loess(responses[[j]] ~ predictors[[j]], 
                       span = 0.3)
      lines(predictors[[j]], fitted(temp.lo), lty = 2, 
            col = 2)
    }
    if (plot.count == 1) {
      if (write.title) {
        title(paste(response.name, " - page ", n.pages, 
                    sep = ""))
      }
      if (rug & is.vector(data[, gbm.call$gbm.x[variable.no]])) {
        rug(quantile(data[, gbm.call$gbm.x[variable.no]], 
                     probs = seq(0, 1, 0.1), na.rm = TRUE))
      }
    }
    else {
      if (write.title & j == 1) {
        title(response.name)
      }
      if (rug & is.vector(data[, gbm.call$gbm.x[k]])) {
        rug(quantile(data[, gbm.call$gbm.x[k]], probs = seq(0, 
                                                            1, 0.1), na.rm = TRUE))
      }
    }
  }
  par(op)
}
