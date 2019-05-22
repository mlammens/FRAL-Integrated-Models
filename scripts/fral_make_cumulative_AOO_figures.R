## ******************************************************************** ##
## fral_make_cumulative_AOO_curves.R
##
## Author: Matthew Aiello-Lammens
## Date Created: 2014-02-16
##
## Purpose:
## Make cumulative AOO curve figures for manuscript based on compiled
## data from `fral_compile_results.R`
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

## ******************************************************************** ##
## ******************************************************************** ##
## Cumulative AOO through time
## ******************************************************************** ##
## ******************************************************************** ##

## ******************************************************************** ##
## Read in data
## ******************************************************************** ##

## -------------------------------------------------------------------- ##
## Look at cumulative occupancy for **all** patches

# Read in three data.frames of cumulative AOO values
cum_occupancy_all <- lapply( X=list.files(path='outputs/', pattern='cumulative_occupancy_all*', full.names=TRUE),
                             FUN=read.csv )
# Combine into a single data.frame
cum_occupancy_all_df <- ldply( cum_occupancy_all )
# Add a factor indicating which occupancy threshold was used
cum_occupancy_all_df$occ_threshold <- NA
cum_occupancy_all_df$occ_threshold[ 1:100 ] <- 1
cum_occupancy_all_df$occ_threshold[ 101:200 ] <- 1000
cum_occupancy_all_df$occ_threshold[ 201:300 ] <- 2000

# Melt the data.frame
cum_occupancy_all_df_m <-
  melt( cum_occupancy_all_df, id.vars=c('Year','occ_threshold') )

# Make a grouping variable based on mp_name and occ_threshold
cum_occupancy_all_df_m$grouping <-
  paste( cum_occupancy_all_df_m$variable, cum_occupancy_all_df_m$occ_threshold, sep="_" )

# Make a model_type variable
cum_occupancy_all_df_m$model_type <-
  sub(pattern='[0-9].*',replacement='',as.character(cum_occupancy_all_df_m$variable))

ggplot( cum_occupancy_all_df_m,
        aes(x=Year,y=value) ) +
  geom_line( aes( group=grouping, colour=factor(occ_threshold) ) )


## -------------------------------------------------------------------- ##
## Look at cumultive occupancy for obsgrd only

# Read in three data.frames of cumulative AOO values
cum_occupancy_obsgrd <- lapply( X=list.files(path='outputs/', pattern='cumulative_occupancy_obsgrd*', full.names=TRUE),
                             FUN=read.csv )
# Combine into a single data.frame
cum_occupancy_obsgrd_df <- ldply( cum_occupancy_obsgrd )
# Add a factor indicating which occupancy threshold was used
cum_occupancy_obsgrd_df$occ_threshold <- NA
cum_occupancy_obsgrd_df$occ_threshold[ 1:100 ] <- 1
cum_occupancy_obsgrd_df$occ_threshold[ 101:200 ] <- 1000
cum_occupancy_obsgrd_df$occ_threshold[ 201:300 ] <- 2000

# Melt the data.frame
cum_occupancy_obsgrd_df_m <-
  melt( cum_occupancy_obsgrd_df, id.vars=c('Year','occ_threshold') )

# Make a grouping variable based on mp_name and occ_threshold
cum_occupancy_obsgrd_df_m$grouping <-
  paste( cum_occupancy_obsgrd_df_m$variable, cum_occupancy_obsgrd_df_m$occ_threshold, sep="_" )

# Make a model_type variable
cum_occupancy_obsgrd_df_m$model_type <-
  sub(pattern='[0-9].*',replacement='',as.character(cum_occupancy_obsgrd_df_m$variable))

## Make a plot with all lines
ggplot() +
  geom_line( data=subset( cum_occupancy_obsgrd_df_m, subset=model_type=="fral_denseff_add_popd_LDD" ),
             aes(x=Year, y=sqrt(value), group=grouping, colour=factor(occ_threshold) ),
             alpha=0.4 ) +
  ylab('Sqrt(Cumulative 20x20 km Grid Cells)') +
  geom_line(data=cum_occ_hist_allPops_df[-1,],
            aes(x=Year,y=sqrt(CumOcc)),size=2) +
  scale_colour_discrete(name="Occupancy\nThreshold") +
  theme_bw() +
  theme(text=element_text(size=12,family="Times",face="bold")) 
#ggsave(filename='figures/Cumulative_AOO_DensEff_PopD_allsims.pdf',width=7.5,height=7.5,units='in')
ggsave(filename='figures/Diss_Fig_4_7.jpg',width=6.5,height=6.5,units='in')

## -------------------------------------------------------------------- ##
## Make plot for the defense talk
ggplot() +
  geom_line( data=subset( cum_occupancy_obsgrd_df_m, 
                          subset= ( model_type=="fral_denseff_add_popd_LDD" & occ_threshold==1000 ) ),
             aes(x=Year, y=sqrt(value), group=grouping, colour=factor(occ_threshold) ),
             alpha=0.4 ) +
  ylab('Sqrt(Cumulative 20x20 km grid cells)') +
  ylim( c(-1,32) ) +
  geom_line(data=cum_occ_hist_allPops_df[-1,],
            aes(x=Year,y=sqrt(CumOcc)),size=2) +
  scale_colour_manual(name="Occupancy\nThreshold", values=c("dodgerblue1") ) +
  guides( colour = FALSE ) +
  theme_bw() +
  theme(text=element_text( size=20, face="bold"),
        axis.title.y=element_text(vjust=0.3) ) 
# ggsave(filename='figures/Pres_Fig_Cumulative_AOO_DensEff_PopD_OccThresh_1.pdf',
#        width=9.5, height=6, units='in')
# ggsave(filename='figures/Pres_Fig_Cumulative_AOO_DensEff_PopD_OccThresh_1.jpg',
#        width=9.5, height=6, units='in')
ggsave(filename='~/Desktop/temp_cumoo.jpg',
       width=9.5, height=6, units='in')


ggplot() +
  geom_line( data=subset( cum_occupancy_obsgrd_df_m, 
                          subset= ( model_type=="fral_denseff_add_popd_LDD" & 
                                      ( occ_threshold==1 | occ_threshold==1000) ) ),
             aes(x=Year, y=sqrt(value), group=grouping, colour=factor(occ_threshold) ),
             alpha=0.4 ) +
  ylab('Sqrt(Cumulative 20x20 km grid cells)') +
  ylim( c(-1,32) ) +
  geom_line(data=cum_occ_hist_allPops_df[-1,],
            aes(x=Year,y=sqrt(CumOcc)),size=2) +
#  scale_colour_discrete(name="Occupancy\nThreshold") +
  scale_colour_manual(name="Occupancy\nThreshold", values=c("dodgerblue1", "blue") ) +
  theme_bw() +
  theme(text=element_text( size=20, face="bold"),
        axis.title.y=element_text(vjust=0.3) ) 
# ggsave(filename='figures/Pres_Fig_Cumulative_AOO_DensEff_PopD_OccThresh_1_1000.pdf',
#        width=9.5, height=6, units='in')
ggsave(filename='figures/Pres_Fig_Cumulative_AOO_DensEff_PopD_OccThresh_1_1000.jpg',
       width=9.5, height=6, units='in')

ggplot() +
  geom_line( data=subset( cum_occupancy_obsgrd_df_m, subset=model_type=="fral_denseff_add_popd_LDD" ),
             aes(x=Year, y=sqrt(value), group=grouping, colour=factor(occ_threshold) ),
             alpha=0.4 ) +
  ylab('Sqrt(Cumulative 20x20 km grid cells)') +
  ylim( c(-1,32) ) +
  geom_line(data=cum_occ_hist_allPops_df[-1,],
            aes(x=Year,y=sqrt(CumOcc)),size=2) +
#  scale_colour_discrete(name="Occupancy\nThreshold") +
  scale_colour_manual(name="Occupancy\nThreshold", values=c("dodgerblue1", "blue", "navyblue") ) +
  theme_bw() +
  theme(text=element_text( size=20, face="bold"),
        axis.title.y=element_text(vjust=0.3) ) 
# ggsave(filename='figures/Pres_Fig_Cumulative_AOO_DensEff_PopD_allsims.pdf',
#        width=9.5, height=6, units='in')
ggsave(filename='figures/Pres_Fig_Cumulative_AOO_DensEff_PopD_allsims.jpg',
       width=9.5, height=6, units='in')


## Mkae a plot comparing land use change and no land use change
ggplot() +
#   geom_line( data=subset( cum_occupancy_obsgrd_df_m, 
#                           subset= ( (model_type=="fral_denseff_add_popd_LDD" | model_type=="fral_denseff_add")
#                                     & occ_threshold==1 ) ),
#              aes(x=Year, y=sqrt(value), group=grouping, colour=factor(model_type) ),
#              alpha=0.4 ) +
  ylab('Sqrt(Cumulative 20x20 km grid cells)') +
  ylim( c(-1,32) ) +
  geom_line(data=cum_occ_hist_allPops_df[-1,],
            aes(x=Year,y=sqrt(CumOcc)),size=2) +
  #scale_colour_manual(name="Occupancy\nThreshold", values=c("dodgerblue1") ) +
  theme_bw() +
  theme(text=element_text( size=20, face="bold"),
        axis.title.y=element_text(vjust=0.3) ) 
ggsave(filename='figures/Pres_Fig_Cumulative_AOO_Hist.jpg',
       width=9.5, height=6, units='in')



## -------------------------------------------------------------------- ##

## Make a plot that summarizes over simulations
ggplot() +
  geom_smooth( data=subset( cum_occupancy_obsgrd_df_m, subset=model_type=="fral_denseff_add_popd_LDD" ),
               aes(x=Year, y=sqrt(value), colour=factor(occ_threshold) ), 
               level=(1-1e-15) ) +
  ylab('Sqrt(Cumulative AOO)') +
  geom_line(data=cum_occ_hist_allPops_df[-1,],
            aes(x=Year,y=sqrt(CumOcc)),size=2) +
  scale_colour_discrete(name="Occupancy\nThreshold") +
  theme_bw() +
  theme(text=element_text(size=12,family="Times",face="bold")) 

## Make a summary data.frame
cum_occupancy_obsgrd_summary <-
  ddply( cum_occupancy_obsgrd_df_m, 
         .variables=c('Year','occ_threshold','model_type'),
         mean_occupancy=mean(value),
         sd_occupancy=sd(value),
         max_occupancy=max(value),
         min_occupancy=min(value),
         summarize )


## Make a grouping varialbe
cum_occupancy_obsgrd_summary$grouping <-
  paste(cum_occupancy_obsgrd_summary$model_type, cum_occupancy_obsgrd_summary$occ_threshold, sep="_")

## Plot summary lines
temp_data <- subset( cum_occupancy_obsgrd_summary, subset=model_type=="fral_denseff_add_popd_LDD" )

# temp_data <- subset( cum_occupancy_obsgrd_summary, 
#                      subset = (model_type=="fral_denseff_add_popd_LDD" |
#                                  model_type=="fral_denseff_add") )

ggplot() +
  geom_line(data=temp_data,
            aes(x=Year,y=sqrt(mean_occupancy),group=grouping,colour=factor(occ_threshold) ),
            size=1.5 ) +
  geom_line(data=temp_data,
            aes(x=Year,y=(sqrt(mean_occupancy + sd_occupancy)),
                group=grouping,colour=factor(occ_threshold) ),
            linetype='dashed') +
  geom_line(data=temp_data,
            aes(x=Year,y=(sqrt(mean_occupancy - sd_occupancy)),
                group=grouping,colour=factor(occ_threshold) ),
            linetype='dashed') +
#   geom_line(data=temp_data,
#             aes(x=Year,y=sqrt(max_occupancy),group=grouping,colour=factor(occ_threshold) ),
#             linetype='dotted') +
#   geom_line(data=temp_data,
#             aes(x=Year,y=sqrt(min_occupancy),group=grouping,colour=factor(occ_threshold) ),
#             linetype='dotted') +
  ylab('Sqrt(Cumulative AOO)') +
  geom_line(data=cum_occ_hist_allPops_df[-1,],
            aes(x=Year,y=sqrt(CumOcc)),size=2) +
  scale_colour_discrete(name="Occupancy\nThreshold") +
  theme_bw() +
  theme(text=element_text(size=12,family="Times",face="bold")) 
ggsave(filename='figures/Cumulative_AOO_DensEff_PopD_Summary.pdf',width=7.5,height=7.5,units='in')


## -------------------------------------------------------------------- ##
## Make figure for presentation
ggplot() +
  geom_line(data=temp_data,
            aes(x=Year,y=sqrt(mean_occupancy),group=grouping,colour=factor(occ_threshold) ),
            size=1.5 ) +
  geom_line(data=temp_data,
            aes(x=Year,y=(sqrt(mean_occupancy + sd_occupancy)),
                group=grouping,colour=factor(occ_threshold) ),
            linetype='dashed') +
  geom_line(data=temp_data,
            aes(x=Year,y=(sqrt(mean_occupancy - sd_occupancy)),
                group=grouping,colour=factor(occ_threshold) ),
            linetype='dashed') +
  #   geom_line(data=temp_data,
  #             aes(x=Year,y=sqrt(max_occupancy),group=grouping,colour=factor(occ_threshold) ),
  #             linetype='dotted') +
  #   geom_line(data=temp_data,
  #             aes(x=Year,y=sqrt(min_occupancy),group=grouping,colour=factor(occ_threshold) ),
  #             linetype='dotted') +
  ylab('Sqrt(Cumulative 20x20 km grid cells)') +
  geom_line(data=cum_occ_hist_allPops_df[-1,],
            aes(x=Year,y=sqrt(CumOcc)),size=2) +
  scale_colour_manual(name="Occupancy\nThreshold", values=c("dodgerblue1", "blue", "navyblue") ) +
  theme_bw() +
  theme(text=element_text(size=20, face="bold"),
        axis.title.y=element_text(vjust=0.3) ) 
ggsave(filename='figures/Pres_Fig_Cumulative_AOO_DensEff_PopD_summarysims.pdf',
       width=9.5, height=6, units='in')





## ******************************************************************** ##
## Calculate cumulative AOO Loss functions
## -------------------------------------------------------------------- ##
## Most of these calculations are duplicated from fral_compile_results.R
## code.
## ******************************************************************** ##

#
cum_occupancy_obsgrd_wide <- 
  reshape(data=cum_occupancy_obsgrd_df, direction='wide',timevar='occ_threshold', idvar='Year' )

# Grab the historical pattern and remove the first line
hist_temp <- cum_occ_hist_allPops_df[-1,]
# Simplify the data to get only the vector of cumulative AOO values
hist_temp <- hist_temp$CumOcc

## Calculate the Loss for each model
cum_AOO_comp_Loss <- apply(cum_occupancy_obsgrd_wide[ 2:ncol(cum_occupancy_obsgrd_wide) ],
                           MARGIN=2,
                           FUN=function(col){sum(abs(hist_temp-col))})
hist(cum_AOO_comp_Loss,breaks=40)

## Make a data.frame of loss functoin values and model names

# First get the model names
mod_names <- colnames(cum_occupancy_obsgrd_wide)[ 2:ncol(cum_occupancy_obsgrd_wide) ]
# Replace "mp.[occ_thresh]" with "mp_[occ_thresh]"
mod_names <- sub( pattern="mp\\.", replacement="mp_", mod_names )
# Make this into a data.frame
cum_AOO_comp_Loss_df <- data.frame(mp_mod=mod_names,
                                   loss=cum_AOO_comp_Loss)

## Add the value of the loss function as an indicator to the melted data.frame
cum_occupancy_obsgrd_df_m <- merge(cum_occupancy_obsgrd_df_m,cum_AOO_comp_Loss_df,
                             by.x='grouping',by.y='mp_mod')


## Code to plot a figure of cumulative AOO with lines colour coded by 
## loss values
# # Take a subset of the data - model_type = "fral_denseff_add_popd_LDD"
# cum_occupancy_temp <- subset(cum_occupancy_obsgrd_df_m, 
#                              subset=model_type=="fral_denseff_add_popd_LDD" )
# 
# ggplot( cum_occupancy_temp ) +
#   geom_line( aes(x=Year,y=value,group=grouping,colour=loss ) ) +
#   #scale_fill_discrete(name="Occupancy\nThreshold") +
#   xlab('Year') +
#   #ylab('Cumulative AOO Loss') +
#   theme(text=element_text(size=12,family="Times",face="bold")) +
#   theme_bw()

## -------------------------------------------------------------------- ##
## Make a box plot for cumulative AOO loss

# Make an indicator variable for model_type
cum_AOO_comp_Loss_df$model_type <-
  sub(pattern='[0-9].*',replacement='',as.character(cum_AOO_comp_Loss_df$mp_mod))
# Make an indicator variable for occupancy threshold
cum_AOO_comp_Loss_df$occ_threshold <-
  sub(pattern='.*_',replacement='',as.character(cum_AOO_comp_Loss_df$mp_mod))

## Box plots with all model types
ggplot( cum_AOO_comp_Loss_df,
        aes(x=occ_threshold,y=loss,fill=model_type) ) +
  geom_boxplot()

## Box plots with only denseff_popd_LDD
ggplot( subset( cum_AOO_comp_Loss_df, subset=model_type=="fral_denseff_add_popd_LDD" ),
        aes(x=occ_threshold,y=loss) ) +
  geom_boxplot() +
  xlab('Occupancy Threshold') +
  ylab('Cumulative AOO Loss') +
  theme_bw() +
  theme(text=element_text(size=12,family="Times",face="bold")) +
  annotate("text", x = 1, y = 34000, label="1" ) +
  annotate("text", x = 2, y = 21000, label="2" ) +
  annotate("text", x = 3, y = 21000, label="2" ) 
#ggsave(filename='figures/Cumulative_AOO_Loss_DensEff_PopD_BoxPlot.pdf',width=7.5,height=7.5,units='in')
ggsave(filename='figures/Diss_Fig_4_8.pdf',width=6.5,height=6.5,units='in')


## -------------------------------------------------------------------- ##
## Make this figure for presentation
ggplot( subset( cum_AOO_comp_Loss_df, subset=model_type=="fral_denseff_add_popd_LDD" ),
        aes(x=occ_threshold,y=loss) ) +
  geom_boxplot() +
  xlab('Occupancy Threshold') +
  ylab('Loss function for cumulative occupancy') +
  theme_bw() +
  theme(text=element_text(size=20,face="bold")) +
  annotate("text", x = 1, y = 34000, label="1" ) +
  annotate("text", x = 2, y = 21000, label="2" ) +
  annotate("text", x = 3, y = 21000, label="2" ) 
ggsave(filename='figures/Pres_Fig_Cumulative_AOO_Loss_DensEff_PopD_BoxPlot.pdf',
       width=7.5,height=6.5,units='in')


## Simple stats test to say whether there are significant differences 
## between Cumulative AOO loss values
cum_aoo_temp <- subset( cum_AOO_comp_Loss_df, subset=model_type=="fral_denseff_add_popd_LDD" )

head(cum_aoo_temp)
anova( lm( data=cum_aoo_temp, formula=loss~as.factor(occ_threshold) ) )
TukeyHSD( aov( data=cum_aoo_temp, formula=loss~as.factor(occ_threshold) ) )

## ******************************************************************** ##
## ******************************************************************** ##
## Sensitivity through time
## ******************************************************************** ##
## ******************************************************************** ##

## ******************************************************************** ##
## Read in data
## ******************************************************************** ##

## -------------------------------------------------------------------- ##
## Sensitivity and positve predictive power considering all patches

# Read in the three data.frames of sensitivity valeus
hist_sim_all <- lapply( X=list.files(path='outputs/', pattern='hist_sim_sens_all*',full.names=TRUE),
                        FUN=read.csv )
# Combine into a single data.frame
hist_sim_all_df <- ldply( hist_sim_all )
# Add a factor indicating which occupancy threshold was used
hist_sim_all_df$occ_threshold <- NA
hist_sim_all_df$occ_threshold[ 1:3000 ] <- 1
hist_sim_all_df$occ_threshold[ 3001:6000 ] <- 1000
hist_sim_all_df$occ_threshold[ 6001:9000 ] <- 2000

# Make a model type vector
hist_sim_all_df$model_type <- sub(pattern='[0-9].*',replacement='',as.character(hist_sim_all_df$mp.file))

## -------------------------------------------------------------------- ##
## Sensitivity and positve predictive power considering patches in 
## group of observed historical patches (i.e., Fral or associated species)

# Read in the three data.frames of sensitivity valeus
hist_sim_obsgrd <- lapply( X=list.files(path='outputs/', pattern='hist_sim_sens_obsgrd*',full.names=TRUE),
                        FUN=read.csv )
# Combine into a single data.frame
hist_sim_obsgrd_df <- ldply( hist_sim_obsgrd )
# Add a factor indicating which occupancy threshold was used
hist_sim_obsgrd_df$occ_threshold <- NA
hist_sim_obsgrd_df$occ_threshold[ 1:3000 ] <- 1
hist_sim_obsgrd_df$occ_threshold[ 3001:6000 ] <- 1000
hist_sim_obsgrd_df$occ_threshold[ 6001:9000 ] <- 2000

# Make a model type vector
hist_sim_obsgrd_df$model_type <- sub(pattern='[0-9].*',replacement='',as.character(hist_sim_obsgrd_df$mp.file))

## ******************************************************************** ##
## Make sensitivity through time figure
## ******************************************************************** ##

## Select a model type to look at
unique(hist_sim_all_df$model_type)
hist_sim_all_df_short <- hist_sim_all_df[ hist_sim_all_df$model_type=='fral_denseff_add_popd_LDD', ]
## Make a grouping variable to uniquely ID mp.file and occupancy threshold
hist_sim_all_df_short$grouping <- 
  paste( hist_sim_all_df_short$mp.file, hist_sim_all_df_short$occ_threshold, sep='')

## Remove model type column
hist_sim_all_df_short$model_type <- NULL
## Remove mp.file
hist_sim_all_df_short$mp.file <- NULL

# Plot trends through time
hist_sim_all_df_m <-
  melt(hist_sim_all_df_short,id.vars=c('grouping','occ_threshold'))
# Add years column
hist_sim_all_df_m$Year <- as.numeric(sub(pattern='.*\\.',replacement='',hist_sim_all_df_m$variable))

## Plots the GAM fitted curves, which show overall model fits
ggplot(hist_sim_all_df_m[grepl(pattern='sens_hist',hist_sim_all_df_m$variable),],
       aes(x=Year,y=value)) + #,group=grouping,colour=factor(occ_threshold) )) +
  #geom_line(aes(group=grouping,colour=factor(occ_threshold))) +
  #geom_point(aes(group=grouping,colour=factor(occ_threshold))) +
  geom_smooth(aes(colour=factor(occ_threshold)),level=(1-1e-15) ) +
  ylab('Sensitivity') +
  #guides(colour=FALSE) +
  #facet_grid(.~occ_threshold) +
  theme_bw()

## ******************************************************************** ##
## Make box plot figure for a few key years
## ******************************************************************** ##

hist_sim_all_df_box <- hist_sim_all_df[ c('mp.file','occ_threshold','model_type',
                                          'sens_hist.1925','sens_hist.1950',
                                          'sens_hist.1975','sens_hist.2000') ]
names(hist_sim_all_df_box) <- c('mp.file','occ_threshold','model_type',
                                '1925','1950','1975','2000')
hist_sim_all_df_box_m <- 
  melt(hist_sim_all_df_box, id.vars=c('mp.file','occ_threshold','model_type') )

ggplot( subset(hist_sim_all_df_box_m, subset=model_type=='fral_denseff_add_popd_LDD') ) +
  geom_boxplot( aes(x=variable,y=value, fill=factor(occ_threshold) ) ) +
  xlab('Year') +
  ylab('Sensitivity') +
  scale_fill_discrete(name="Occupancy\nThreshold") +
  theme_bw() +
  theme(text=element_text(size=12,family="Times",face="bold")) 
        #axis.title.y=element_text(size=12,family="Times",face="bold")) +
#ggsave(filename='figures/Sensitivity_DensEff_PopD_BoxPlot.pdf',width=7.5,height=7.5,units='in')
ggsave(filename='figures/Diss_Fig_4_4.pdf',width=6.5,height=6.5,units='in')

## -------------------------------------------------------------------- ##
## Make figure for presentation

ggplot( subset(hist_sim_all_df_box_m, subset=model_type=='fral_denseff_add_popd_LDD') ) +
  geom_boxplot( aes(x=variable,y=value, fill=factor(occ_threshold) ) ) +
  xlab('Year') +
  ylab('Sensitivity') +
  scale_fill_discrete(name="Occupancy\nThreshold") +
  theme_bw() +
  theme( text=element_text( size=20, face="bold") ) 
#axis.title.y=element_text(size=12,family="Times",face="bold")) +
ggsave( filename='figures/Pres_Fig_Sensitivity_DensEff_PopD_BoxPlot.pdf',
        width=7.5, height=6.5, units='in')


## ******************************************************************** ##
## Make a box plot figure for the number of grid cells occupied
## -------------------------------------------------------------------- ##
## Note that here I am using only the gridcells observed as either 
## having FRAL or an associated species: 974 patches
## ******************************************************************** ##

# First make the dataset
occupied_df <- hist_sim_obsgrd_df[ c('mp.file','occ_threshold','model_type',
                                  'grid_overlap.1925','sim_only.1925',
                                  'grid_overlap.1950','sim_only.1950',
                                  'grid_overlap.1975','sim_only.1975',
                                  'grid_overlap.2000','sim_only.2000' ) ]
occupied_df$'1925' <- occupied_df$grid_overlap.1925 + occupied_df$sim_only.1925
occupied_df$'1950' <- occupied_df$grid_overlap.1950 + occupied_df$sim_only.1950
occupied_df$'1975' <- occupied_df$grid_overlap.1975 + occupied_df$sim_only.1975
occupied_df$'2000' <- occupied_df$grid_overlap.2000 + occupied_df$sim_only.2000

occupied_df <- occupied_df[ c('mp.file','occ_threshold','model_type',
                              '1925','1950','1975','2000') ]

occupied_df_m <- 
  melt(occupied_df, id.vars=c('mp.file','occ_threshold','model_type') )

ggplot( subset(occupied_df_m, subset=model_type=='fral_denseff_add_popd_LDD') ) +
  geom_boxplot( aes(x=variable,y=value, fill=factor(occ_threshold) ) ) +
  xlab('Year') +
  ylab('Cumulative occupied patchs') +
  scale_fill_discrete(name="Occupancy\nThreshold") +
  geom_hline(yintercept=974,colour="red",size=1.5) +
  annotate("text", x=1, y=934,label="Max patches") +
#   geom_hline(yintercept=459,colour="black",size=1.5) +
#   annotate("text", x=1, y=419,label="Total F. alnus patches") +
  theme_bw() +
  theme(text=element_text(size=12,family="Times",face="bold")) 

#ggsave(filename='figures/Occupancy_DensEff_PopD_BoxPlot.pdf',width=7.5,height=7.5,units='in')
ggsave(filename='figures/Diss_Fig_4_5.pdf',width=6.5,height=6.5,units='in')

## -------------------------------------------------------------------- ##
## Make figure for presentaiton

ggplot( subset(occupied_df_m, subset=model_type=='fral_denseff_add_popd_LDD') ) +
  geom_boxplot( aes(x=variable,y=value, fill=factor(occ_threshold) ) ) +
  xlab('Year') +
  ylab('Cumulative occupied patchs') +
  scale_fill_discrete(name="Occupancy\nThreshold") +
  geom_hline(yintercept=974,colour="red",size=1.5) +
  annotate("text", x=1, y=934,label="Max patches") +
  #   geom_hline(yintercept=459,colour="black",size=1.5) +
  #   annotate("text", x=1, y=419,label="Total F. alnus patches") +
  theme_bw() +
  theme(text=element_text(size=20, face="bold")) 
ggsave(filename='figures/Pres_Fig_Occupancy_DensEff_PopD_BoxPlot.pdf',
       width=7.5, height=6.5,units='in')


## Get number of patches occupied in time
cum_occ_hist_allPops_df[ which( cum_occ_assoc_hist_allPops_df$Year %in% c(1925,1950,1975,2000,2010) ), ]