## ******************************************************************** ##
## fral_analyze_occthresh_2000_results.R
##
## Author: Matthew Aiello-Lammens
## Date Created: 2014-01-01
##
## Purpose:
## Carryout a Boosted Regression Tree analysis on compiled results to 
## determine the relationship between input value and model outups AND
## to examine parameter importance on model outcomes.
##
## USE DATA SET ASSUMING OCCURENCE THRESHOLD = 1
##
## In running this script I am assuming that 'fral_compile_results.R' 
## has alredy been run
## ******************************************************************** ##

## ******************************************************************** ##
## Call necessary functions and packages
## ******************************************************************** ##

## Source the fral_demog_setup.R script to setup required functions
source('R/fral_demog_setup.R')

## Read in simulation results
Sim_Results_2000 <- read.csv('outputs/sim_results_2000.csv')

## -------------------------------------------------------------------- ##
## Read in mean population size data.frame
Sim_PopSizes <- read.csv('outputs/pop_sizes_mpMult.csv')
# Make a model type vector
Sim_PopSizes$model_type <- sub(pattern='[0-9].*',replacement='',as.character(Sim_PopSizes$mp.file))

ggplot(Sim_PopSizes, aes(x=model_type, y=pop_size_t100_mean)) +
  geom_boxplot() +
  coord_cartesian(ylim=c(0,1e7))

ggplot(Sim_PopSizes, aes(x=pop_size_t100_mean)) +
  geom_density(aes(colour=model_type))

## Recall that there are 1e8 2x2 m grid cells per patch.

## ******************************************************************** ##
## Set predictor and response variables
## ******************************************************************** ##

## Set predictors
predictors <- c(1,7,8,47,60,63,64,65)
pred_cat <- c(70)
names(Sim_Results_2000)[c(predictors,pred_cat)]

## Define response variables
resp_hist_sens_obsgrd <- grep( "hist_sens_obsgrd", names(Sim_Results_2000) ) # Mean sensitivity obsgrd
resp_cum_AOO_loss <- grep( "cum_AOO_loss", names(Sim_Results_2000) ) # cumulative AOO loss function values
resp_sens_1950 <- grep( "sens_hist.1950", names(Sim_Results_2000) ) # 1950 sensitivity
resp_sens_2010 <- grep( "sens_hist.2010", names(Sim_Results_2000) ) # 2010 sensitivity

# Choose reponse to use in further analyses
response <- resp_cum_AOO_loss
## -------------------------------------------------------------------- ##

Sim_Results_2000_Small <- Sim_Results_2000[ c(predictors,pred_cat,response) ]
Sim_Results_2000_Small_m <- melt(Sim_Results_2000_Small,measure.vars=names(Sim_Results_2000)[predictors])

ggplot(Sim_Results_2000_Small_m,aes(x=value,y=cum_AOO_loss,colour=kch.type)) + ### #y=hist_sens
  geom_point(alpha=0.5) +
  facet_wrap(~variable,scales='free')

## This is a very messy plot, but helps to demonstrate that there are no
## obvious differences between the two dd.type scenarios. There are also
## no obvious differences between the kch.types. 

## ggsave(filename='figures/Scatter_Pred_vs_MeanAnnSens.pdf',width=7.5,height=7.5,units='in')

## ******************************************************************** ##
## Compare different model structures
## ******************************************************************** ##
## Guide to Different Scenarios/ Structures
## 
## Models     | Dens Type | Hab Type  | Paired
## -----------|-----------|-----------|-------
## 1 - 500    | Ceiling   | Dynamic   | TRUE
## 501 - 1000 | Eff Dens  | Dynamic   | TRUE
## 1001-1500  | Ceiling   | Static    | TRUE
## 1501-2000  | Eff Dens  | Static    | TRUE
## 2001-2500  | Eff Dens  | Dynamic   | w/ 2501-3000 
## 2501-3000  | Eff Dens  | Dynamic   | w/ 2001-2500
##
## ** LDD random for 1 - 2500
## ** LDD based on population density for 2501-3000
## ******************************************************************** ##

## Check on the order of models by looking at the mp.file name for the 
## first model in each section of 500
Sim_Results_2000$mp.file[ c(1,501,1001,1501,2001,2501) ]
## -------------------------------------------------------------------- ##

## ******************************************************************** ##
## Subset density dependence types
## ******************************************************************** ##

# Sensitivity
Sim_Results_2000_Diff_HistSens <- Sim_Results_2000$hist_sens[1:500] - Sim_Results_2000$hist_sens[501:1000]
mean(Sim_Results_2000_Diff_HistSens)
qplot(x=Sim_Results_2000_Diff_HistSens,geom="histogram",binwidth=0.001)
t.test(x=Sim_Results_2000$hist_sens[1:500],
       y=Sim_Results_2000$hist_sens[501:1000],
       paired=TRUE)

# Cumulative AOO Loss
Sim_Results_2000_Diff_loss <- Sim_Results_2000$cum_AOO_loss[1:500] - Sim_Results_2000$cum_AOO_loss[501:1000]
mean(Sim_Results_2000_Diff_loss)
qplot(x=Sim_Results_2000_Diff_loss, geom="histogram")
t.test(x=Sim_Results_2000$cum_AOO_loss[1:500],
       y=Sim_Results_2000$cum_AOO_loss[501:1000],
       paired=TRUE)

# Positive predictive power
Sim_Results_2000_Diff_ppp <- Sim_Results_2000$sim_sens[1:500] - Sim_Results_2000$sim_sens[501:1000]
mean(Sim_Results_2000_Diff_ppp)
qplot(x=Sim_Results_2000_Diff_ppp, geom="histogram")
t.test(x=Sim_Results_2000$sim_sens[1:500],
       y=Sim_Results_2000$sim_sens[501:1000],
       paired=TRUE)

# Delta EMA
Sim_Results_2000_Diff_EMA <- Sim_Results_2000$exp.min.n[1:500] - Sim_Results_2000$exp.min.n[501:1000]
mean(Sim_Results_2000_Diff_EMA)
qplot(x=Sim_Results_2000_Diff_EMA,geom="histogram")
t.test(x=Sim_Results_2000$exp.min.n[1:500],
       y=Sim_Results_2000$exp.min.n[501:1000],
       paired=TRUE)

# Delta N_Mean
Sim_Results_2000_Diff_N <- Sim_Results_2000$n.mean[1:500] - Sim_Results_2000$n.mean[501:1000]
mean(Sim_Results_2000_Diff_N)
median(Sim_Results_2000_Diff_N)
qplot(x=Sim_Results_2000_Diff_N,geom="histogram")
t.test(x=Sim_Results_2000$n.mean[1:500],
       y=Sim_Results_2000$n.mean[501:1000],
       paired=TRUE)

Sim_PopSizes_Diff_N <- Sim_PopSizes$pop_size_t100_mean[1:500] - Sim_PopSizes$pop_size_t100_mean[501:1000]
mean(Sim_PopSizes_Diff_N)
qplot(Sim_PopSizes_Diff_N)
t.test(x=Sim_PopSizes$pop_size_t100_mean[1:500],
       y=Sim_PopSizes$pop_size_t100_mean[501:1000],
       paired=TRUE)

## There is some evidence that effective density results in higher 
## total meta-population sizes, but this distribution still straddles
## zero. The largest differences have very high fecundity, which makes
## sense because if you can survive past seedling stage, your survival
## is relatively high in this model.

## ******************************************************************** ##
## Subset HS change types - EFFECTIVE POPULATION DENSITY DEPENDENCE
## ******************************************************************** ##

# Difference in cumulative AOO loss
Sim_Results_2000_hsDiff_EffDens <-
  Sim_Results_2000$cum_AOO_loss[501:1000] - Sim_Results_2000$cum_AOO_loss[1501:2000]
mean(Sim_Results_2000_hsDiff_EffDens)
qplot(x=Sim_Results_2000_hsDiff_EffDens,geom="histogram")
t.test(x=Sim_Results_2000$cum_AOO_loss[501:1000],
       y=Sim_Results_2000$cum_AOO_loss[1501:2000],
       paired=TRUE)

# Difference in sensitivity
Sim_Results_2000_hsDiff_EffDens_Sens <-
  Sim_Results_2000$hist_sens[501:1000] - Sim_Results_2000$hist_sens[1501:2000]
mean(Sim_Results_2000_hsDiff_EffDens_Sens)
qplot(x=Sim_Results_2000_hsDiff_EffDens_Sens,geom="histogram")
t.test(x=Sim_Results_2000$hist_sens[501:1000],
       y=Sim_Results_2000$hist_sens[1501:2000],
       paired=TRUE)

# PPP
Sim_Results_2000_hsDiff_EffDens_PPP <-
  Sim_Results_2000$sim_sens[501:1000] - Sim_Results_2000$sim_sens[1501:2000]
mean(Sim_Results_2000_hsDiff_EffDens_PPP,na.rm=TRUE)
qplot(x=Sim_Results_2000_hsDiff_EffDens_PPP,geom="histogram")
t.test(x=Sim_Results_2000$sim_sens[501:1000],
       y=Sim_Results_2000$sim_sens[1501:2000],
       paired=TRUE)

# Difference in EMA
Sim_Results_2000_hsDiff_EffDens_EMA <-
  Sim_Results_2000$exp.min.n[501:1000] - Sim_Results_2000$exp.min.n[1501:2000]
mean(Sim_Results_2000_hsDiff_EffDens_EMA)
ggplot(data=NULL,aes(x=Sim_Results_2000_hsDiff_EffDens_EMA)) +
  geom_histogram() +
  geom_vline(x=0,col='red') +
  geom_vline(x=mean(Sim_Results_2000_hsDiff_EffDens_EMA),col='blue')
t.test(x=Sim_Results_2000$exp.min.n[501:1000],
       y=Sim_Results_2000$exp.min.n[1501:2000],
       paired=TRUE)

# Difference in final N
# Delta N_Mean
Sim_Results_2000_hsDiff_EffDens_N <- 
  Sim_Results_2000$n.mean[501:1000] - Sim_Results_2000$n.mean[1501:2000]
mean(Sim_Results_2000_hsDiff_EffDens_N)
median(Sim_Results_2000_hsDiff_EffDens_N)
qplot(x=Sim_Results_2000_hsDiff_EffDens_N,geom="histogram")
t.test(x=Sim_Results_2000$n.mean[501:1000],
       y=Sim_Results_2000$n.mean[1501:2000],
       paired=TRUE)

## Here there is some evidence that population size decreases due to 
## habitat change, which makes sense because there is an overall decrease
## in carrying capacity across the landscape through time when the 
## dynamic landscape is considered. However, for the most part,
## there are no compelling differences between these two scenarios.

## ******************************************************************** ##
## Subset HS change types - CEILING TYPE DENSITY DEPENDENCE
## ******************************************************************** ##

# Difference in cumulative AOO loss
Sim_Results_2000_hsDiff_Ceil <-
  Sim_Results_2000$cum_AOO_loss[1:500] - Sim_Results_2000$cum_AOO_loss[1001:1500]
mean(Sim_Results_2000_hsDiff_Ceil)
qplot(x=Sim_Results_2000_hsDiff_Ceil,geom="histogram")
t.test(x=Sim_Results_2000$cum_AOO_loss[1:500],
       y=Sim_Results_2000$cum_AOO_loss[1001:1500],
       paired=TRUE)

# Difference in sensitivity
Sim_Results_2000_hsDiff_Ceil_Sens <-
  Sim_Results_2000$hist_sens[1:500] - Sim_Results_2000$hist_sens[1001:1500]
mean(Sim_Results_2000_hsDiff_Ceil_Sens)
qplot(x=Sim_Results_2000_hsDiff_Ceil_Sens,geom="histogram")
t.test(x=Sim_Results_2000$hist_sens[1:500],
       y=Sim_Results_2000$hist_sens[1001:1500],
       paired=TRUE)

# Difference in EMA
Sim_Results_2000_hsDiff_Ceil_EMA <-
  Sim_Results_2000$exp.min.n[1:500] - Sim_Results_2000$exp.min.n[1001:1500]
mean(Sim_Results_2000_hsDiff_Ceil_EMA)
ggplot(data=NULL,aes(x=Sim_Results_2000_hsDiff_Ceil_EMA)) +
  geom_histogram() +
  geom_vline(x=0,col='red') +
  geom_vline(x=mean(Sim_Results_2000_hsDiff_Ceil_EMA),col='blue')
t.test(x=Sim_Results_2000$exp.min.n[1:500],
       y=Sim_Results_2000$exp.min.n[1001:1500],
       paired=TRUE)

# Difference in final N
# Delta N_Mean
Sim_Results_2000_hsDiff_Ceil_N <- 
  Sim_Results_2000$n.mean[1:500] - Sim_Results_2000$n.mean[1001:1500]
mean(Sim_Results_2000_hsDiff_Ceil_N)
median(Sim_Results_2000_hsDiff_Ceil_N)
qplot(x=Sim_Results_2000_hsDiff_Ceil_N,geom="histogram")
t.test(x=Sim_Results_2000$n.mean[1:500],
       y=Sim_Results_2000$n.mean[1001:1500],
       paired=TRUE)


## ******************************************************************** ##
## ******************************************************************** ##
## Look at only denseff models WITH HS change
## -------------------------------------------------------------------- ##
## By further examining these sets of models, we are looking at only
## 1500 simulations. 1000 simulations that use random long-distance
## dispersal (LDD) and 500 simulations that use LDD based on human
## population density through time.
## ******************************************************************** ##
## ******************************************************************** ##

Sim_Results_2000_DensEff <- Sim_Results_2000[Sim_Results_2000$dd.type=='eff_dens',]
Sim_Results_2000_DensEff <- Sim_Results_2000_DensEff[!grepl(pattern='nochng',Sim_Results_2000_DensEff$mp.file),]

# Make an indicator for the simulations that used a LDD model based on population density
Sim_Results_2000_DensEff$popd <- 'rand_LDD'
Sim_Results_2000_DensEff$popd[ grepl( 'popd', Sim_Results_2000_DensEff$mp.file ) ] <- 'popd_LDD'

## Add an indicator of cum_AOO_loss threshold value
Sim_Results_2000_DensEff$loss_cat <- NA
Sim_Results_2000_DensEff$loss_cat[Sim_Results_2000_DensEff$cum_AOO_loss<=2500] <- 'A'
Sim_Results_2000_DensEff$loss_cat[Sim_Results_2000_DensEff$cum_AOO_loss>2500 & 
                               Sim_Results_2000_DensEff$cum_AOO_loss<=5000] <- 'B'
Sim_Results_2000_DensEff$loss_cat[Sim_Results_2000_DensEff$cum_AOO_loss>5000 & 
                               Sim_Results_2000_DensEff$cum_AOO_loss<=10000] <- 'C'
Sim_Results_2000_DensEff$loss_cat[Sim_Results_2000_DensEff$cum_AOO_loss>10000]  <- 'D'

## Add a binary indicator of cum_AOO_loss threshold value
Sim_Results_2000_DensEff$loss_bin <- NA
Sim_Results_2000_DensEff$loss_bin[Sim_Results_2000_DensEff$cum_AOO_loss<=3000] <- TRUE #'lt3000'
Sim_Results_2000_DensEff$loss_bin[Sim_Results_2000_DensEff$cum_AOO_loss>3000 ] <- FALSE #'gt3000'

## Add a binary indicator for sensitivity values
Sim_Results_2000_DensEff$sens_bin <- NA
Sim_Results_2000_DensEff$sens_bin[Sim_Results_2000_DensEff$hist_sens<0.5] <- 'lt_0.5'
Sim_Results_2000_DensEff$sens_bin[Sim_Results_2000_DensEff$hist_sens>=0.5 ] <- 'ge_0.5'

## Calcualte the difference between sensitivity and PPP
Sim_Results_2000_DensEff$sens_diff <-
  Sim_Results_2000_DensEff$hist_sens_obsgrd - Sim_Results_2000_DensEff$sim_sens_obsgrd

## Add a binary indicator for optimum sensitivity vs PPP values
Sim_Results_2000_DensEff$sens_ppp_opt <- 0
Sim_Results_2000_DensEff$sens_ppp_opt[(Sim_Results_2000_DensEff$sens_diff<=0.1 & 
                                         Sim_Results_2000_DensEff$hist_sens>=0.5)] <- 1

# Seperate out these two sets into seperate data.frames
Sim_Results_2000_DensEff_rand_LDD <- Sim_Results_2000_DensEff[ Sim_Results_2000_DensEff$popd=='rand_LDD', ]
Sim_Results_2000_DensEff_popd_LDD <- Sim_Results_2000_DensEff[ Sim_Results_2000_DensEff$popd=='popd_LDD', ]

## Examine distribution of loss values for these models
mean(Sim_Results_2000_DensEff_rand_LDD$cum_AOO_loss)
median(Sim_Results_2000_DensEff_rand_LDD$cum_AOO_loss)
quantile(Sim_Results_2000_DensEff_rand_LDD$cum_AOO_loss,probs=c(0.05,0.10,0.15,0.20))
quantile(Sim_Results_2000_DensEff_popd_LDD$cum_AOO_loss,probs=c(0.05,0.10,0.15,0.20))

## -------------------------------------------------------------------- ##
## Plot histogram of Cumulative AOO Loss values for the two different
## LDD scenarios
ggplot(Sim_Results_2000_DensEff,aes(x=cum_AOO_loss)) +
  geom_histogram() +
  #geom_vline(xintercept=mean(Sim_Results_2000_DensEff$cum_AOO_loss),linetype='longdash',colour='grey') +
  #geom_vline(xintercept=median(Sim_Results_2000_DensEff$cum_AOO_loss),colour='red') +
  facet_grid(popd~.,scales="free") +
  xlab("Cumulative AOO Loss") +
  theme_bw()

## -------------------------------------------------------------------- ##
## Plot histograms of predictor variables with NO classifications
ggplot(melt(Sim_Results_2000_DensEff_popd_LDD,measure.vars=names(Sim_Results_2000)[ predictors ]),
       aes(x=value)) + #,fill=popd)) +
  geom_histogram(position="dodge") + #,aes(y=..density..)) +
  #geom_density(size=1) +
  facet_wrap(~variable,scales='free_x') +
  theme_bw()
# ggsave(filename='figures/Histogram_Predictor_Variables_NoClass.pdf',width=7.5,height=7.5)

## Make a pairs plot for these data
ggpairs( Sim_Results_2000_DensEff_popd_LDD[ c(predictors,pred_cat) ] )

## -------------------------------------------------------------------- ##
## Plot histograms of predictor variables WITH classifications

## Four classifications for cumulative AOO loss
ggplot(melt(Sim_Results_2000_DensEff_rand_LDD,measure.vars=names(Sim_Results_2000)[predictors]),
       aes(x=value,colour=loss_cat)) +
  geom_histogram(position="dodge",alpha=0.4,aes(y=..density..,fill=loss_cat)) +
  geom_density(size=1) +
  facet_wrap(~variable,scales='free') +
  theme_bw()

## Binary classification for cumulative AOO loss
ggplot(melt(Sim_Results_2000_DensEff_rand_LDD,measure.vars=names(Sim_Results_2000)[predictors]),
       aes(x=value,colour=loss_bin)) +
  geom_histogram(position="dodge",alpha=0.4,aes(y=..density..,fill=loss_bin)) +
  geom_density(size=1) +
  facet_wrap(~variable,scales='free') +
  scale_colour_discrete(guide=FALSE) +
  scale_fill_discrete(name="Cumulative\nAOO Loss",
                      breaks=c(FALSE,TRUE),
                      labels=c("gt3000","lt300") ) +
  theme( legend.position=c(1,0), legend.justification=c(1,0) ) +
  theme_bw()
# ggsave(filename='figures/Histogram_Predictor_Variables_CumAOOLoss_BinClass.pdf',width=7.5,height=7.5)

## ******************************************************************** ##
## Make density plots of all predictor values and only those that matched
## well the simulation fit metrics
## -------------------------------------------------------------------- ##
## These plots are similar to prior - posterior plots
## ******************************************************************** ##

## sensitivity - positive predictive power optimum
ggplot() +
  geom_density( data=melt(Sim_Results_2000_DensEff_popd_LDD,measure.vars=names(Sim_Results_2000)[ predictors ]),
                aes(x=value), fill='grey40', alpha=0.4) +
  geom_density( data=melt(Sim_Results_2000_DensEff_popd_LDD[ Sim_Results_2000_DensEff_popd_LDD$sens_ppp_opt==1, ]
                          ,measure.vars=names(Sim_Results_2000)[ predictors ]),
                aes(x=value), fill='grey80', alpha=0.4) +
  facet_wrap(~variable,scales='free') +
  theme_bw()
# ggsave(filename='figures/Prior_posterior_density_2000_sens_ppp_opt.pdf',width=7.5,height=7.5)


## Cumulative AOO loss less than or equal to 5000 - as is, this is 189 of 500 simulations
ggplot() +
  geom_density( data=melt(Sim_Results_2000_DensEff_popd_LDD,measure.vars=names(Sim_Results_2000)[ predictors ]),
                aes(x=value), fill='grey40', alpha=0.4) +
  geom_density( data=melt(Sim_Results_2000_DensEff_popd_LDD[ Sim_Results_2000_DensEff_popd_LDD$cum_AOO_loss<=5000, ],
                          measure.vars=names(Sim_Results_2000)[ predictors ]),
                aes(x=value), fill='grey80', alpha=0.4) +
  facet_wrap(~variable,scales='free') +
  theme_bw()
# ggsave(filename='figures/Prior_posterior_density_2000_aoo_loss.pdf',width=7.5,height=7.5)

## sensitivity values greater than 0.5
ggplot() +
  geom_density( data=melt(Sim_Results_2000_DensEff_popd_LDD,measure.vars=names(Sim_Results_2000)[ predictors ]),
                aes(x=value), fill='grey40', alpha=0.4) +
  geom_density( data=melt(Sim_Results_2000_DensEff_popd_LDD[ Sim_Results_2000_DensEff_popd_LDD$hist_sens>=0.5, ]
                          ,measure.vars=names(Sim_Results_2000)[ predictors ]),
                aes(x=value), fill='grey80', alpha=0.4) +
  facet_wrap(~variable,scales='free') +
  theme_bw()
# ggsave(filename='figures/Prior_posterior_density_2000_hist_sens.pdf',width=7.5,height=7.5)


## ******************************************************************** ##
## Plot xy-plots
## ******************************************************************** ##

## Predictors versus cumulative AOO Loss response
ggplot(melt(Sim_Results_2000_DensEff,measure.vars=names(Sim_Results_2000)[predictors]),
       aes(x=value,y=cum_AOO_loss,colour=popd)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~variable,scales='free')

## Predictors versus number of occupied patches within the set of historically
## observed patches
ggplot(melt(Sim_Results_2000_DensEff,measure.vars=names(Sim_Results_2000)[predictors]),
       aes(x=value,y=occupied_pops_obsgrd,colour=popd)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~variable,scales='free')

## Predictors versus number sensitivity
ggplot(melt(Sim_Results_2000_DensEff,measure.vars=names(Sim_Results_2000)[predictors]),
       aes(x=value,y=hist_sens,colour=popd)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~variable,scales='free')

## -------------------------------------------------------------------- ##
## Plot various endpoint measures against each other
ggpairs( Sim_Results_2000_DensEff[ c(66,69,76:78,79) ], colour='popd' )

ggplot(Sim_Results_2000_DensEff,
       aes(x=occupied_pops_obsgrd,y=cum_AOO_loss,colour=popd)) +
  geom_point()

ggplot(Sim_Results_2000_DensEff,
       aes(x=hist_sens,y=sens_hist.2010,colour=popd)) +
  geom_point()

ggplot(Sim_Results_2000_DensEff,aes(x=sens_hist.2010,fill=popd)) +
  geom_histogram(aes(y=..density..),position="dodge")
  
## -------------------------------------------------------------------- ##
## Look at the relationship with hist_sens and the predictors, given 
## categorization based on loss
ggplot(melt(Sim_Results_2000_DensEff,measure.vars=names(Sim_Results_2000)[predictors]),
       aes(x=value,y=hist_sens_obsgrd,colour=loss_cat)) +#,colour=dd.type)) +
  geom_point() +
  facet_wrap(~variable,scales='free')

## -------------------------------------------------------------------- ##
## Make pair plots
## ***
## This plot takes a very long time to construct, so only do this if needed
Sim_Results_2000_DensEff_short <- Sim_Results_2000_DensEff_popd_LDD[c(predictors,pred_cat,80)]
### ggpairs(Sim_Results_2000_DensEff_short,colour='loss_cat')

# ## Plot only loss less than 3000
Sim_Results_2000_DensEff_LossLT3000 <- Sim_Results_2000_DensEff_rand_LDD[ Sim_Results_2000_DensEff_rand_LDD$cum_AOO_loss <=3000, ]
### ggpairs(Sim_Results_2000_DensEff_LossLT3000[c(predictors,pred_cat)])


## ******************************************************************** ##
## Look at cum aoo with other predictors
## ******************************************************************** ##

## Merge cum_occ_mpMult_df_m with Sim_Results data.frame 
# cum_occ_mpMult_df_m_FULL <- merge(cum_occ_mpMult_df_m,Sim_Results_2000_DensEff,
#                                   by.x='variable',by.y='mp.file')
# 
# cum_occ_temp <- cum_occ_mpMult_df_m_FULL
# # cum_occ_temp <- cum_occ_mpMult_df_m_FULL[ (cum_occ_mpMult_df_m_FULL$dd.type=='eff_dens' &
# #                                              cum_occ_mpMult_df_m_FULL$hs.chng), ]
# 
# ggplot() + 
#   geom_line(data=cum_occ_temp,
#             aes(x=Year,y=sqrt(value), group=variable,
#                 colour=cum_AOO_loss) # colour=hist_sens_obsgrd)
#             ,alpha=0.5) + 
#   scale_colour_gradient(low="blue",high="red") +
#   #scale_colour_gradient(limits=c(0,3000),low="blue",high="red") +
#   ylab('Sqrt(Cumulative 20x20 km Grid Cells)') +
#   geom_line(data=hist_temp,aes(x=Year,y=sqrt(CumOcc)),size=2) +
#   facet_grid(popd~.) +
#   theme_bw()
# # ggsave(filename='figures/Cumulative_AOO_Loss_curves.pdf',width=10.5,height=7.5)

## ******************************************************************** ##
## Compare population density based LDD with random LDD
## -------------------------------------------------------------------- ##
## These are paired simulations for 500 of the random LDD simulations
## of the total of 1000 random LDD simulations.
## ******************************************************************** ##

## Check mp.file names
Sim_Results_2000_DensEff$mp.file[ c(1, 501, 1001) ]

## This confirms that 501-1000 is paired with 1001-1500

## Cumulative AOO Loss
LDD_Scenario_Diff_2000_cum_aoo_loss <- 
  Sim_Results_2000_DensEff$cum_AOO_loss[ 1001:1500 ] - Sim_Results_2000_DensEff$cum_AOO_loss[ 501:1000 ]
mean( LDD_Scenario_Diff_2000_cum_aoo_loss )
qplot( LDD_Scenario_Diff_2000_cum_aoo_loss, geom="histogram" )
t.test( x=Sim_Results_2000_DensEff$cum_AOO_loss[ 1001:1500 ],
        y=Sim_Results_2000_DensEff$cum_AOO_loss[ 501:1000 ],
        paired=TRUE )

## While the distibution overlaps 0, it is heavily left tailed, indicating
## that cumulative AOO loss is *lower* using the popd_LDD

## Sensitivity
LDD_Scenario_Diff_2000_sens <- 
  Sim_Results_2000_DensEff$hist_sens[ 1001:1500 ] - Sim_Results_2000_DensEff$hist_sens[ 501:1000 ]
mean( LDD_Scenario_Diff_2000_sens )
qplot( LDD_Scenario_Diff_2000_sens, geom="histogram" )
t.test( x=Sim_Results_2000_DensEff$hist_sens[ 1001:1500 ],
        y=Sim_Results_2000_DensEff$hist_sens[ 501:1000 ],
        paired=TRUE )

## Historic sensitivity shows an increase in the popd_LDD scenario

## PPP
LDD_Scenario_Diff_2000_ppp <- 
  Sim_Results_2000_DensEff$sim_sens_obsgrd[ 1001:1500 ] - Sim_Results_2000_DensEff$sim_sens_obsgrd[ 501:1000 ]
mean( LDD_Scenario_Diff_2000_ppp )
qplot( LDD_Scenario_Diff_2000_ppp, geom="histogram" )
t.test( x=Sim_Results_2000_DensEff$sim_sens_obsgrd[ 1001:1500 ],
        y=Sim_Results_2000_DensEff$sim_sens_obsgrd[ 501:1000 ],
        paired=TRUE )

## The effect on positive predictive rate is less compelling

## ******************************************************************** ##
## ******************************************************************** ##
## Run BRT analysis
## ******************************************************************** ##
## ******************************************************************** ##

## -------------------------------------------------------------------- ##
## Cumulative AOO (continuous)
fral_denseff_2000_BRT_cum_aoo_loss <- 
  gbm.step(data=Sim_Results_2000_DensEff_rand_LDD,gbm.x=c(predictors,pred_cat),gbm.y=response, #Sim_Results_2000_DensEff_popd_LDD
           family="gaussian",tree.complexity=4,learning.rate=0.01,
           tolerance.method="fixed",tolerance=.01,
           bag.fraction=0.5)

fral_denseff_2000_BRT_cum_aoo_loss_summary <- summary(fral_denseff_2000_BRT_cum_aoo_loss)
fral_denseff_2000_BRT_cum_aoo_loss_inter <- gbm.interactions(gbm.object=fral_denseff_2000_BRT_cum_aoo_loss)

pdf(file='figures/BRT_Response_Curves_2000_cum_AOO_loss.pdf',width=7.5,height=7.5)
gbm.plot.cust(fral_denseff_2000_BRT_cum_aoo_loss,n.plots=9,write.title=TRUE,plot.layout=c(3,3))#,common.scale=FALSE)
dev.off()

## Make interaction plots
make_gbm_interaction_plots(gbm_model=fral_denseff_2000_BRT_cum_aoo_loss,zrange=c(0,20000),
                           fig_file='figures/BRT_Interactions_2000_cum_AOO_loss.pdf')

## Fecundity by metapop.initab plotted seperately
gbm.perspec( fral_denseff_2000_BRT_cum_aoo_loss, 
             x=7, y=1, z.range=c(0,20000) )

## -------------------------------------------------------------------- ##
## Thresholded cumulative AOO loss (categorical/ binary)
fral_denseff_2000_BRT_loss_bin <- 
  gbm.step(data=Sim_Results_2000_DensEff_rand_LDD,gbm.x=c(predictors,pred_cat),gbm.y=81,
           family="bernoulli",tree.complexity=4,learning.rate=0.001,
           #tolerance.method="fixed",tolerance=.01,
           bag.fraction=0.5)

fral_denseff_2000_BRT_loss_bin_summary <- summary(fral_denseff_2000_BRT_loss_bin)
fral_denseff_2000_BRT_loss_bin_inter <- gbm.interactions(gbm.object=fral_denseff_2000_BRT_loss_bin)

pdf(file='figures/BRT_Response_Curves_2000_loss_bin.pdf',width=7.5,height=7.5)
gbm.plot.cust(fral_denseff_2000_BRT_loss_bin,n.plots=9,write.title=TRUE,plot.layout=c(3,3))#,common.scale=FALSE)
dev.off()

## -------------------------------------------------------------------- ##
## Historical sensitivity (continuous)
fral_denseff_2000_BRT_histSens <- 
  gbm.step(data=Sim_Results_2000_DensEff_rand_LDD,gbm.x=c(predictors,pred_cat),gbm.y=66,
           family="gaussian",tree.complexity=4,learning.rate=0.001,
           #tolerance.method="fixed",tolerance=.01,
           bag.fraction=0.5)

fral_denseff_2000_BRT_histSens_summary <- summary(fral_denseff_2000_BRT_histSens)
fral_denseff_2000_BRT_histSens_inter <- gbm.interactions(gbm.object=fral_denseff_2000_BRT_histSens)

pdf(file='figures/BRT_Response_Curves_2000_hist_sens.pdf',width=7.5,height=7.5)
gbm.plot.cust(fral_denseff_2000_BRT_histSens,n.plots=9,write.title=TRUE,plot.layout=c(3,3))#,common.scale=FALSE)
dev.off()

## Make interaction plots
make_gbm_interaction_plots(gbm_model=fral_denseff_2000_BRT_histSens,zrange=c(0,1),
                           fig_file='figures/BRT_Interactions_2000_histSens.pdf')


## -------------------------------------------------------------------- ##
## Trade-off between sensitivity and PPP (categorical/ binary)

# First look at how many simulations meet this threshold
sum(Sim_Results_2000_DensEff_rand_LDD$sens_ppp_opt) # 94 / 1000

fral_denseff_2000_BRT_sens_ppp_opt <- 
  gbm.step(data=Sim_Results_2000_DensEff_rand_LDD,gbm.x=c(predictors,pred_cat),gbm.y=84,
           family="bernoulli",tree.complexity=4,learning.rate=0.001,
           #tolerance.method="fixed",tolerance=.01,
           bag.fraction=0.5)

fral_denseff_2000_BRT_sens_ppp_opt_summary <- summary(fral_denseff_2000_BRT_sens_ppp_opt)
fral_denseff_2000_BRT_sens_ppp_opt_inter <- gbm.interactions(gbm.object=fral_denseff_2000_BRT_sens_ppp_opt)

pdf(file='figures/BRT_Response_Curves_2000_sens_ppp_opt.pdf',width=7.5,height=7.5)
gbm.plot.cust(fral_denseff_2000_BRT_sens_ppp_opt,n.plots=9,write.title=TRUE,plot.layout=c(3,3))#,common.scale=FALSE)
dev.off()

## Make interaction plots
make_gbm_interaction_plots(gbm_model=fral_denseff_2000_BRT_sens_ppp_opt,zrange=c(0,1),
                           fig_file='figures/BRT_Interactions_2000_sens_ppp_opt.pdf')

gbm.perspec( fral_denseff_2000_BRT_sens_ppp_opt, 
             x=7, y=1 )

## ******************************************************************** ##
## Look at the popd_LDD responses
## ******************************************************************** ##

## -------------------------------------------------------------------- ##
## Cumulative AOO (continuous)
fral_denseff_popd_2000_BRT_cum_aoo_loss <- 
  gbm.step(data=Sim_Results_2000_DensEff_rand_LDD,gbm.x=c(predictors,pred_cat),gbm.y=response, #Sim_Results_2000_DensEff_popd_LDD
           family="gaussian",tree.complexity=4,learning.rate=0.01,
           tolerance.method="fixed",tolerance=.01,
           bag.fraction=0.5)

fral_denseff_popd_2000_BRT_cum_aoo_loss_summary <- summary(fral_denseff_popd_2000_BRT_cum_aoo_loss)
fral_denseff_popd_2000_BRT_cum_aoo_loss_inter <- gbm.interactions(gbm.object=fral_denseff_popd_2000_BRT_cum_aoo_loss)

pdf(file='figures/BRT_Response_Curves_popd_2000_cum_AOO_loss.pdf',width=7.5,height=7.5)
gbm.plot.cust(fral_denseff_popd_2000_BRT_cum_aoo_loss,n.plots=9,write.title=TRUE,plot.layout=c(3,3))#,common.scale=FALSE)
dev.off()

## Make interaction plots
make_gbm_interaction_plots(gbm_model=fral_denseff_popd_2000_BRT_cum_aoo_loss,zrange=c(0,20000),
                           fig_file='figures/BRT_Interactions_popd_2000_cum_AOO_loss.pdf')

## Fecundity by metapop.initab plotted seperately
gbm.perspec( fral_denseff_popd_2000_BRT_cum_aoo_loss, 
             x=7, y=1, z.range=c(0,20000) )

## -------------------------------------------------------------------- ##
## Historical sensitivity (continuous)
fral_denseff_popd_2000_BRT_histSens <- 
  gbm.step(data=Sim_Results_2000_DensEff_popd_LDD,gbm.x=c(predictors,pred_cat),gbm.y=66,
           family="gaussian",tree.complexity=4,learning.rate=0.001,
           #tolerance.method="fixed",tolerance=.01,
           bag.fraction=0.5)

fral_denseff_popd_2000_BRT_histSens_summary <- summary(fral_denseff_popd_2000_BRT_histSens)
fral_denseff_popd_2000_BRT_histSens_inter <- gbm.interactions(gbm.object=fral_denseff_popd_2000_BRT_histSens)

pdf(file='figures/BRT_Response_Curves_popd_2000_hist_sens.pdf',width=7.5,height=7.5)
gbm.plot.cust(fral_denseff_popd_2000_BRT_histSens,n.plots=9,write.title=TRUE,plot.layout=c(3,3))#,common.scale=FALSE)
dev.off()

## Make interaction plots
make_gbm_interaction_plots(gbm_model=fral_denseff_popd_2000_BRT_histSens,zrange=c(0,1),
                           fig_file='figures/BRT_Interactions_popd_2000_histSens.pdf')


## -------------------------------------------------------------------- ##
## Trade-off between sensitivity and PPP (categorical/ binary)

# First look at how many simulations meet this threshold
sum(Sim_Results_2000_DensEff_popd_LDD$sens_ppp_opt) # 97 / 500

fral_denseff_popd_2000_BRT_sens_ppp_opt <- 
  gbm.step(data=Sim_Results_2000_DensEff_popd_LDD,gbm.x=c(predictors,pred_cat),gbm.y=84,
           family="bernoulli",tree.complexity=4,learning.rate=0.001,
           #tolerance.method="fixed",tolerance=.01,
           bag.fraction=0.5)

fral_denseff_popd_2000_BRT_sens_ppp_opt_summary <- summary(fral_denseff_popd_2000_BRT_sens_ppp_opt)
fral_denseff_popd_2000_BRT_sens_ppp_opt_inter <- gbm.interactions(gbm.object=fral_denseff_popd_2000_BRT_sens_ppp_opt)

pdf(file='figures/BRT_Response_Curves_popd_2000_sens_ppp_opt.pdf',width=7.5,height=7.5)
gbm.plot.cust(fral_denseff_popd_2000_BRT_sens_ppp_opt,n.plots=9,write.title=TRUE,plot.layout=c(3,3))#,common.scale=FALSE)
dev.off()

## Make interaction plots
make_gbm_interaction_plots(gbm_model=fral_denseff_popd_2000_BRT_sens_ppp_opt,zrange=c(0,1),
                           fig_file='figures/BRT_Interactions_popd_2000_sens_ppp_opt.pdf')

gbm.perspec( fral_denseff_popd_2000_BRT_sens_ppp_opt, 
             x=7, y=1 )


## ******************************************************************** ##
## ******************************************************************** ##
## Calculate some population size values
## ******************************************************************** ##
## ******************************************************************** ##

## Maximum Fral density based on maximum-mean patch population size
max(Sim_PopSizes$pop_size_t100_mean)/(20000^2) # 0.09 plants per m^2

## Mean Fral densities for each simulation type
Mean_Pop_Sizes <- 
  ddply(.data=Sim_PopSizes,.variables='model_type',
        Mean_Pop_Size=mean(pop_size_t100_mean),
        Max_Pop_Size=max(pop_size_t100_mean),
        summarize)
## Calculate mean densities
Mean_Pop_Sizes$Mean_Density <-
  Mean_Pop_Sizes$Mean_Pop_Size/(20000^2)

Mean_Pop_Sizes$Max_Density <-
  Mean_Pop_Sizes$Max_Pop_Size/(20000^2)

print( Mean_Pop_Sizes )

## ******************************************************************** ##
## ******************************************************************** ##
### Development Work = Principal Component Analysis ###
## ******************************************************************** ##
## ******************************************************************** ##

## Multivariate Analyses
sim_res_denseff_short <- Sim_Results_DensEff_rand_LDD[predictors]
sim_res_denseff_short <- scale(sim_res_denseff_short)

require(vegan)
fral_pca <- rda(sim_res_denseff_short)
fral_pca$CA$v
biplot(fral_pca,type=c('text','points'))
biplot(fral_pca,type=c('text','points'),choices=c(1,3))
biplot(fral_pca,type=c('text','points'),choices=c(2,3))


sim_res_denseff_pca <- cbind(Sim_Results_DensEff_rand_LDD,fral_pca$CA$u)

ggplot(sim_res_denseff_pca,aes(x=PC1,y=PC2,colour=factor(sens_ppp_opt))) +
  geom_point()

## There is no evidence that a PCA reduces the axes of variation in 
## this data.set. Nor do the PC axes help seperate out measures of
## simulation fit to the historic pattern.
