## ******************************************************************** ##
## fral_lagphase_compare.R
##
## Author: Matthew Aiello-Lammens
## Date Created: 2014-02-01
##
## Purpose:
## Compare lag phase calculations base on observations made by herbarium
## versus all patches
##
## ******************************************************************** ##


## ******************************************************************** ##
## Look at cumulative area of occupancy curves

cum_occ_mpMult_histSamp_df <- make_cum_occ_df( cum_occ_list=cum_occupancy_mpMult, histSamp=TRUE )
cum_occ_mpMult_df <- make_cum_occ_df( cum_occ_list=cum_occupancy_mpMult, histSamp=FALSE )

## Only use *denseff* simulations
sims2use <- which( grepl( pattern='denseff', names(cum_occ_mpMult_df) ) )
## Add first column back in (Year)
sims2use <- c(1,sims2use)
cum_occ_mpMult_df <- cum_occ_mpMult_df[ sims2use ]
cum_occ_mpMult_histSamp_df <- cum_occ_mpMult_histSamp_df[ sims2use ]

## Remove *nochng* simulations
sims2use <- which( !grepl( pattern='nochng', names(cum_occ_mpMult_df) ) )
cum_occ_mpMult_df <- cum_occ_mpMult_df[ sims2use ]
cum_occ_mpMult_histSamp_df <- cum_occ_mpMult_histSamp_df[ sims2use ]


# Grab the historical pattern and remove the first line
hist_temp <- cum_occ_hist_allPops_df[-1,]

## -------------------------------------------------------------------- ##
## Calculate LOSS between the simulation and historic cumulative AOO
## curves
cum_occ_mpMult_df_m <- 
  calc_loss( cum_occ_df=cum_occ_mpMult_df, hist_temp=hist_temp )
cum_occ_mpMult_histSamp_df_m <- 
  calc_loss( cum_occ_df=cum_occ_mpMult_histSamp_df, hist_temp=hist_temp )

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
  geom_line(data=cum_occ_mpMult_histSamp_df_m,
            aes(x=Year,y=sqrt(value), group=variable,colour=loss)) +
  #scale_colour_gradient(limits=c(2500,23000),low="blue",high="red") +
  scale_colour_gradient(limits=c(0,3000),low="blue",high="red") +
  ylab('Sqrt(Cumulative 20x20 km Grid Cells)') +
  geom_line(data=hist_temp,aes(x=Year,y=sqrt(CumOcc)),size=2) +
  theme_bw()


## -------------------------------------------------------------------- ##
## Plot both sets of data on the same plot
cum_occ_mpMult_df_m$SampleGrids <- 'all'
cum_occ_mpMult_histSamp_df_m$SampleGrids <- 'hist'
# Combine these
cum_occ_mpMult_df_m_COMBINED <-
  rbind( cum_occ_mpMult_df_m, cum_occ_mpMult_histSamp_df_m )

# Make a new 'variable'
cum_occ_mpMult_df_m_COMBINED$mp_samp <- 
  paste( cum_occ_mpMult_df_m_COMBINED$variable, cum_occ_mpMult_df_m_COMBINED$SampleGrids, sep='_' )

ggplot() + 
  geom_line(data=cum_occ_mpMult_df_m_COMBINED,
            aes(x=Year,y=sqrt(value), group=mp_samp,colour=SampleGrids), alpha=0.4) +
  ylab('Sqrt(Cumulative 20x20 km Grid Cells)') +
  geom_line(data=hist_temp,aes(x=Year,y=sqrt(CumOcc)),size=2) +
  #geom_hline(yintercept=sqrt(3423)) + # Max line
  theme_bw()

## ******************************************************************** ##
## Calculate the slopes of the cumulative AOO curves for two time 
## periods. 
## -------------------------------------------------------------------- ##
## I'm going to use a few time periods. When Resit and I first 
## discussed this, we thought of using 1910 to 1925 versus 1925 to 1940.
## Based on the results of Larkin 2011 I would also like to use up to 1940
## versus 1940 to 1970.
## ******************************************************************** ##

## In order to do this analysis, I have to subset my cumulative occupancy
## data.frame
cum_occ_1911_1925 <- cum_occ_mpMult_df_m_COMBINED[ (cum_occ_mpMult_df_m_COMBINED$Year %in% 1911:1925), ]
cum_occ_1926_1940 <- cum_occ_mpMult_df_m_COMBINED[ (cum_occ_mpMult_df_m_COMBINED$Year %in% 1926:1940), ]

## Calculate slope and R^2 values
lm_fit_results_1911_1925 <- 
  ddply( .data=cum_occ_1911_1925, .variables=c('variable','SampleGrids'),
         lm_slope=summary( lm( value ~ Year ) )$coefficients[2,1],
         lm_R2=summary( lm( value ~ Year ) )$r.squared,
         summarize )
lm_fit_results_1926_1940 <- 
  ddply( .data=cum_occ_1926_1940, .variables=c('variable','SampleGrids'),
         lm_slope=summary( lm( value ~ Year ) )$coefficients[2,1],
         lm_R2=summary( lm( value ~ Year ) )$r.squared,
         summarize )
## Merge the two data.frames
lm_fit_results_1911_1940_comb <-
  cbind( lm_fit_results_1911_1925, lm_fit_results_1926_1940[ c('lm_slope','lm_R2') ] )
names(lm_fit_results_1911_1940_comb) <- 
  c('mp_sim','sample_grids','slope_1911_1925','r2_1911_1925','slope_1926_1940','r2_1926_1940')
## Calculate change in slopes
lm_fit_results_1911_1940_comb$slope_delta <- 
  lm_fit_results_1911_1940_comb$slope_1926_1940 - lm_fit_results_1911_1940_comb$slope_1911_1925

## Plot histogram of slopes
ggplot( lm_fit_results_1911_1940_comb, aes(x=slope_delta,fill=sample_grids) ) +
  #geom_density()
  geom_histogram(position="dodge")

## Make a data.frame to calculate the paired differences between restricting ourselves
## to the historic sample grid versus not
slope_diffs_1911_1940 <-
  dcast(data=lm_fit_results_1911_1940_comb,formula=mp_sim~sample_grids,value.var='slope_delta')
slope_diffs_1911_1940$paired_diff <-
  slope_diffs_1911_1940$all - slope_diffs_1911_1940$hist

## Melt this data.frame to plot it
slope_diffs_1911_1940_m <-
  melt(slope_diffs_1911_1940)

ggplot( slope_diffs_1911_1940_m, aes(x=value,colour=variable) ) +
  geom_density()

## -------------------------------------------------------------------- ##
## Look at the second set of years 1911-1940 and 1941-1970
## In order to do this analysis, I have to subset my cumulative occupancy
## data.frame
cum_occ_1911_1940 <- cum_occ_mpMult_df_m_COMBINED[ (cum_occ_mpMult_df_m_COMBINED$Year %in% 1911:1940), ]
cum_occ_1941_1970 <- cum_occ_mpMult_df_m_COMBINED[ (cum_occ_mpMult_df_m_COMBINED$Year %in% 1941:1970), ]

## Calculate slope and R^2 values
lm_fit_results_1911_1940 <- 
  ddply( .data=cum_occ_1911_1940, .variables=c('variable','SampleGrids'),
         lm_slope=summary( lm( value ~ Year ) )$coefficients[2,1],
         lm_R2=summary( lm( value ~ Year ) )$r.squared,
         summarize )
lm_fit_results_1941_1970 <- 
  ddply( .data=cum_occ_1941_1970, .variables=c('variable','SampleGrids'),
         lm_slope=summary( lm( value ~ Year ) )$coefficients[2,1],
         lm_R2=summary( lm( value ~ Year ) )$r.squared,
         summarize )
## Merge the two data.frames
lm_fit_results_1911_1970 <-
  cbind( lm_fit_results_1911_1940, lm_fit_results_1941_1970[ c('lm_slope','lm_R2') ] )
names(lm_fit_results_1911_1970) <- 
  c('mp_sim','sample_grids','slope_1911_1940','r2_1911_1940','slope_1941_1970','r2_1941_1970')
## Calculate change in slopes
lm_fit_results_1911_1970$slope_delta <- 
  lm_fit_results_1911_1970$slope_1941_1970 - lm_fit_results_1911_1970$slope_1911_1940

## Plot histogram of slopes
ggplot( lm_fit_results_1911_1970, aes(x=slope_delta,colour=sample_grids) ) +
  geom_density()

## Make a data.frame to calculate the paired differences between restricting ourselves
## to the historic sample grid versus not
slope_diffs_1911_1970 <-
  dcast(data=lm_fit_results_1911_1970,formula=mp_sim~sample_grids,value.var='slope_delta')
slope_diffs_1911_1970$paired_diff <-
  slope_diffs_1911_1970$all - slope_diffs_1911_1970$hist

## Melt this data.frame to plot it
slope_diffs_1911_1970_m <-
  melt(slope_diffs_1911_1970)

ggplot( slope_diffs_1911_1970_m, aes(x=value,colour=variable) ) +
  geom_density() +
  xlab('Change in slope from 1911-1940 to 1940-1970')
