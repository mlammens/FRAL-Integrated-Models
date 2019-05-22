## dev_work.R

sim_res <- rbind(Sim_Results_DensEff_All, Sim_Results_DensEff)
dim(sim_res)

# Make and indicator if popd or not
sim_res$popd <- FALSE
sim_res$popd[ grepl( 'popd', sim_res$mp.file ) ] <- TRUE

## Mean historical sensitivity versus mean simulation sensitivity
ggplot(sim_res,aes(x=hist_sens_obsgrd,y=sim_sens_obsgrd,colour=loss_cat_char,shape=popd)) + 
  geom_point(size=3) +
  theme_bw()

# Check for match between popd and no-popd scenarios
which( sim_res$NumTranslocs[501:1000] != sim_res$NumTranslocs[1001:1500] )

cum_aoo_diff <- sim_res$cum_AOO_loss[1001:1500] - sim_res$cum_AOO_loss[501:1000]
qplot(cum_aoo_diff, geom="histogram")

hist_sens_diff <- sim_res$hist_sens_obsgrd[1001:1500] - sim_res$hist_sens_obsgrd[501:1000]
qplot(hist_sens_diff, geom="histogram")

ggplot(data=NULL,aes(x=cum_aoo_diff,y=hist_sens_diff)) +
  geom_point() +
  geom_smooth(method='lm')
summary(lm(hist_sens_diff~cum_aoo_diff)) # No predictive power

ggplot(data=NULL, aes(x=sim_res$hist_sens_obsgrd[1001:1500], y=sim_res$hist_sens_obsgrd[501:1000]) ) +
  geom_point() +
  geom_abline(slope=1,intercept=0,colour='red')

ggplot(data=NULL, aes(x=sim_res$cum_AOO_loss[1001:1500], y=sim_res$cum_AOO_loss[501:1000]) ) +
  geom_point() +
  geom_abline(slope=1,intercept=0,colour='red') +
  geom_vline(xintercept=3000,colour='blue')


ggplot(melt(sim_res[1001:1500,],measure.vars=names(Sim_Results)[predictors]),
       aes(x=value,colour=loss_cat)) +
  geom_histogram(position="identity",alpha=0.4,aes(y=..density..)) +
  geom_density(size=1) +
  facet_wrap(~variable,scales='free') +
  theme_bw()



## Mean historical sensitivity versus mean simulation sensitivity
ggplot(sim_res,aes(x=hist_sens_obsgrd,y=sim_sens_obsgrd,colour=loss_cat_char)) + 
  geom_point(size=3) +
  theme_bw()

ggplot(sim_res,aes(x=hist_sens_obsgrd,y=sens_hist.2010,colour=loss_cat_char,shape=popd)) + 
  geom_point(size=3) +
  theme_bw()
