## ******************************************************************** ##
## fral_make_figures.R
##
## Author: Matthew Aiello-Lammens
## Date Created: 2014-01-26
##
## Purpose:
## This is a script to keep most of the figure making code in one 
## place.
##
## First exectue fral_compile_results.R and fral_analyze_results.R
##
## ******************************************************************** ##

## Define a function for extracting a ggplot figure legend
## This code is lifted from: http://stackoverflow.com/questions/12539348/ggplot-separate-legend-and-plot
## and: http://stackoverflow.com/questions/12041042/how-to-plot-just-the-legends-in-ggplot2
#Extract Legend 
g_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 

## ******************************************************************** ##
## Make scatter plots with marginal density plots
## ******************************************************************** ##

## Calculate a convex hull for the points in the lt3000 group
#sim_res_temp <- Sim_Results_DensEff[ Sim_Results_DensEff$loss_cat=='lt3000', ]
#hull <- sim_res_temp[ chull( sim_res_temp$Fec.Mean, sim_res_temp$NumTranslocs ), ]

## Plot Fecundity versus LDD (NumTranslocs)
scatter <- ggplot() +
  geom_point(data=Sim_Results_DensEff_popd_LDD,
             aes(x=Fecundity,y=LDD,shape=factor(sens_ppp_opt),colour=factor(sens_ppp_opt)), #colour=loss_bin),
             size=4) +
  #scale_shape_manual(values=c(20,17)) +
  #geom_polygon( data=hull, aes(x=Fec.Mean,y=NumTranslocs), alpha=1, colour='black', fill=NA ) +
  theme_bw() +
  theme(legend.position='none') 
fec_dens <- ggplot() +
  geom_density(data=Sim_Results_DensEff_popd_LDD,aes(x=Fecundity,colour=factor(sens_ppp_opt))) + #colour=loss_bin)) +
  theme_bw() +
  theme(legend.position='none',
        axis.text.x=element_blank(), axis.title.x=element_blank()) 
ldd_dens <- ggplot() +
  geom_density(data=Sim_Results_DensEff_popd_LDD,aes(x=LDD,colour=factor(sens_ppp_opt))) + #colour=loss_bin)) +
  coord_flip() +
  theme_bw() +
  theme(legend.position='none',
        axis.text.y=element_blank(), axis.title.y=element_blank()) 
# empty <- ggplot()+geom_point(aes(1,1), colour="white")+
#   theme(axis.ticks=element_blank(), 
#        panel.background=element_blank(), 
#        axis.text.x=element_blank(), axis.text.y=element_blank(),           
#        axis.title.x=element_blank(), axis.title.y=element_blank())

# Get legend
fec_dens_wLegend <- ggplot() +
  geom_density(data=Sim_Results_DensEff_popd_LDD,aes(x=Fecundity,colour=factor(sens_ppp_opt))) + #colour=loss_bin)) +
  theme_bw()
legend <- g_legend(fec_dens_wLegend)

grid.arrange(fec_dens,legend,scatter,ldd_dens,ncol=2,nrow=2,widths=c(4,1),heights=c(1,4))

# ## Save plot
# pdf('figures/Fec_by_LDD.pdf',width=7.5,height=7.5)
# grid.arrange(fec_dens,legend,scatter,ldd_dens,ncol=2,nrow=2,widths=c(4,1),heights=c(1,4))
# dev.off()

### Temp stuff
print(scatter)
print(fec_dens_wLegend)
print(ldd_dens)

## ******************************************************************** ##
## Make scatter plots with marginal density plots - Fec x mean.t0.disp.rate
## ******************************************************************** ##

## Plot Fecundity versus LDD (mean.t0.disp.rate)
scatter <- ggplot() +
  geom_point(data=Sim_Results_DensEff_popd_LDD,
             aes(x=Fec.Mean,y=mean.t0.disp.rate,shape=loss_bin,colour=loss_bin),
             size=4) +
  #scale_shape_manual(values=c(20,17)) +
  #geom_polygon( data=hull, aes(x=Fec.Mean,y=mean.t0.disp.rate), alpha=1, colour='black', fill=NA ) +
  theme_bw() +
  theme(legend.position='none') 
fec_dens <- ggplot() +
  geom_density(data=Sim_Results_DensEff_popd_LDD,aes(x=Fec.Mean,colour=loss_bin)) +
  theme_bw() +
  theme(legend.position='none',
        axis.text.x=element_blank(), axis.title.x=element_blank()) 
ldd_dens <- ggplot() +
  geom_density(data=Sim_Results_DensEff_popd_LDD,aes(x=mean.t0.disp.rate,colour=loss_bin)) +
  coord_flip() +
  theme_bw() +
  theme(legend.position='none',
        axis.text.y=element_blank(), axis.title.y=element_blank()) 

# Get legend
fec_dens_wLegend <- ggplot() +
  geom_density(data=Sim_Results_DensEff_popd_LDD,aes(x=Fec.Mean,colour=loss_bin)) +
  theme_bw()
legend <- g_legend(fec_dens_wLegend)

grid.arrange(fec_dens,legend,scatter,ldd_dens,ncol=2,nrow=2,widths=c(4,1),heights=c(1,4))

## ******************************************************************** ##
## Repeat above plot, but use sensisitivity instead
## ******************************************************************** ##
## Plot Fecundity versus LDD (NumTranslocs)
scatter <- ggplot() +
  geom_point(data=Sim_Results_DensEff_popd_LDD,
             aes(x=Fec.Mean,y=NumTranslocs,shape=factor(sens_ppp_opt),colour=factor(sens_ppp_opt)),
             size=4) +
  #scale_shape_manual(values=c(20,17)) +
  #geom_polygon( data=hull, aes(x=Fec.Mean,y=NumTranslocs), alpha=1, colour='black', fill=NA ) +
  theme_bw() +
  theme(legend.position='none') 
fec_dens <- ggplot() +
  geom_density(data=Sim_Results_DensEff_popd_LDD,aes(x=Fec.Mean,colour=factor(sens_ppp_opt))) +
  theme_bw() +
  theme(legend.position='none',
        axis.text.x=element_blank(), axis.title.x=element_blank()) 
ldd_dens <- ggplot() +
  geom_density(data=Sim_Results_DensEff_popd_LDD,aes(x=NumTranslocs,colour=factor(sens_ppp_opt))) +
  coord_flip() +
  theme_bw() +
  theme(legend.position='none',
        axis.text.y=element_blank(), axis.title.y=element_blank()) 
# empty <- ggplot()+geom_point(aes(1,1), colour="white")+
#   theme(axis.ticks=element_blank(), 
#        panel.background=element_blank(), 
#        axis.text.x=element_blank(), axis.text.y=element_blank(),           
#        axis.title.x=element_blank(), axis.title.y=element_blank())

# Get legend
fec_dens_wLegend <- ggplot() +
  geom_density(data=Sim_Results_DensEff_popd_LDD,aes(x=Fec.Mean,colour=factor(sens_ppp_opt))) +
  theme_bw()
legend <- g_legend(fec_dens_wLegend)

grid.arrange(fec_dens,legend,scatter,ldd_dens,ncol=2,nrow=2,widths=c(4,1),heights=c(1,4))


## LDD (NumTranslocs) versus survival-growth
scatter <- ggplot() +
  geom_point(data=Sim_Results_DensEff_popd_LDD,
             aes(x=Surv.Mean,y=NumTranslocs,shape=factor(sens_ppp_opt),colour=factor(sens_ppp_opt)),
             size=4) +
  #scale_shape_manual(values=c(20,17)) +
  #geom_polygon( data=hull, aes(x=Surv.Mean,y=NumTranslocs), alpha=1, colour='black', fill=NA ) +
  theme_bw() +
  theme(legend.position='none') 
Surv_dens <- ggplot() +
  geom_density(data=Sim_Results_DensEff_popd_LDD,aes(x=Surv.Mean,colour=factor(sens_ppp_opt))) +
  theme_bw() +
  theme(legend.position='none',
        axis.text.x=element_blank(), axis.title.x=element_blank()) 
ldd_dens <- ggplot() +
  geom_density(data=Sim_Results_DensEff_popd_LDD,aes(x=NumTranslocs,colour=factor(sens_ppp_opt))) +
  coord_flip() +
  theme_bw() +
  theme(legend.position='none',
        axis.text.y=element_blank(), axis.title.y=element_blank()) 
# empty <- ggplot()+geom_point(aes(1,1), colour="white")+
#   theme(axis.ticks=element_blank(), 
#        panel.background=element_blank(), 
#        axis.text.x=element_blank(), axis.text.y=element_blank(),           
#        axis.title.x=element_blank(), axis.title.y=element_blank())

# Get legend
Surv_dens_wLegend <- ggplot() +
  geom_density(data=Sim_Results_DensEff_popd_LDD,aes(x=Surv.Mean,colour=factor(sens_ppp_opt))) +
  theme_bw()
legend <- g_legend(Surv_dens_wLegend)

grid.arrange(Surv_dens,legend,scatter,ldd_dens,ncol=2,nrow=2,widths=c(4,1),heights=c(1,4))


## ******************************************************************** ##
## Make density plots of all predictor values and only those that matched
## well the simulation fit metrics
## -------------------------------------------------------------------- ##
## These plots are similar to prior - posterior plots
## ******************************************************************** ##

## sensitivity - positive predictive power optimum
temp <- melt(Sim_Results_DensEff_popd_LDD,measure.vars=names(Sim_Results)[ predictors ])
temp$variable <- factor( temp$variable, fral_denseff_popd_1k_BRT_sens_ppp_opt_summary$var )

ggplot() +
  geom_density( data=temp,
                aes(x=value), fill='grey40', alpha=0.4) +
  geom_density( data=melt(Sim_Results_DensEff_popd_LDD[ Sim_Results_DensEff_popd_LDD$sens_ppp_opt==1, ]
                          ,measure.vars=names(Sim_Results)[ predictors ]),
                aes(x=value), fill='grey80', alpha=0.4) +
  facet_wrap(~variable,scales='free') +
  theme_bw() +
  theme( text=element_text( size=12, family="Times"),
         axis.ticks.y=element_blank(),
         axis.text.y=element_blank(),
         axis.title.y=element_blank() )
#ggsave(filename='figures/Prior_posterior_density_1k_sens_ppp_opt.pdf',width=7.5,height=7.5)
# ggsave(filename='figures/Diss_Fig_4_13.pdf',width=6.5,height=6.5)

## -------------------------------------------------------------------- ##
## Make figure for presentation

## Only Fecundity for starters
ggplot() +
  geom_density( data = Sim_Results_DensEff_popd_LDD, 
                aes( x=Fecundity ), fill="grey40", alpha=0.4 ) +
  geom_density( data = Sim_Results_DensEff_popd_LDD[ Sim_Results_DensEff_popd_LDD$sens_ppp_opt==1, ], 
                aes( x=Fecundity ), fill="grey80", alpha=0.4 ) +
  theme_bw() +
  theme( text=element_text( size=20, face="bold" ),
         axis.ticks.y=element_blank(),
         axis.text.y=element_blank(),
         axis.title.y=element_blank() )
# ggsave(filename='figures/Pres_Fig_Prior_posterior_density_1k_sens_ppp_opt_Fecundity.pdf',width=5.5,height=5.5)




ggplot() +
  geom_density( data=temp,
                aes(x=value), fill='grey40', alpha=0.4) +
  geom_density( data=melt(Sim_Results_DensEff_popd_LDD[ Sim_Results_DensEff_popd_LDD$sens_ppp_opt==1, ]
                          ,measure.vars=names(Sim_Results)[ predictors ]),
                aes(x=value), fill='grey80', alpha=0.4) +
  facet_wrap(~variable,scales='free') +
  theme_bw() +
  theme( text=element_text( size=18 ),
         axis.ticks.y=element_blank(),
         axis.text.y=element_blank(),
         axis.title.y=element_blank() )
# ggsave(filename='figures/Pres_Fig_Prior_posterior_density_1k_sens_ppp_opt.pdf',width=7.5,height=7.5)


## Also make a plot of the differences between kch scenarios
Sim_Results_DensEff_popd_LDD$Carrying_Capacity <- 
  factor( Sim_Results_DensEff_popd_LDD$Carrying_Capacity, c( "LOW", "MED", "HIGH" ) )

ggplot() +
  geom_histogram( data=Sim_Results_DensEff_popd_LDD, 
                  aes( x=Carrying_Capacity, fill=factor( sens_ppp_opt) ),
                  position="dodge") +
  scale_fill_manual( values = c( "grey40", "grey80" ) ) +
  theme_bw() +
  theme( text=element_text( size=18 ),
         axis.title.x=element_blank(),
         legend.position = "none" )
# ggsave(filename='figures/Pres_Fig_Prior_posterior_density_1k_sens_ppp_opt_K.pdf',width=2.5,height=2.5)

  

rm(temp)

## Cumulative AOO loss less than or equal to 5000 - as is, this is 189 of 500 simulations
temp <- melt(Sim_Results_DensEff_popd_LDD,measure.vars=names(Sim_Results)[ predictors ])
temp$variable <- factor( temp$variable, fral_denseff_popd_1k_BRT_cum_aoo_loss_summary$var )

ggplot() +
  geom_density( data=temp,
                aes(x=value), fill='grey40', alpha=0.4) +
  geom_density( data=melt(Sim_Results_DensEff_popd_LDD[ Sim_Results_DensEff_popd_LDD$cum_AOO_loss<=5000, ],
                          measure.vars=names(Sim_Results)[ predictors ]),
                aes(x=value), fill='grey80', alpha=0.4) +
  facet_wrap(~variable,scales='free') +
  theme_bw() +
  theme( text=element_text( size=12, family="Times"),
         axis.ticks.y=element_blank(),
         axis.text.y=element_blank(),
         axis.title.y=element_blank() )

#ggsave(filename='figures/Prior_posterior_density_1k_aoo_loss.pdf',width=7.5,height=7.5)
ggsave(filename='figures/Diss_Fig_4_10.pdf',width=6.5,height=6.5)

rm(temp)

## sensitivity values greater than 0.5
ggplot() +
  geom_density( data=melt(Sim_Results_DensEff_popd_LDD,measure.vars=names(Sim_Results)[ predictors ]),
                aes(x=value), fill='grey40', alpha=0.4) +
  geom_density( data=melt(Sim_Results_DensEff_popd_LDD[ Sim_Results_DensEff_popd_LDD$hist_sens>=0.5, ]
                          ,measure.vars=names(Sim_Results)[ predictors ]),
                aes(x=value), fill='grey80', alpha=0.4) +
  facet_wrap(~variable,scales='free') +
  theme_bw()
ggsave(filename='figures/Prior_posterior_density_1k_hist_sens.pdf',width=7.5,height=7.5)

                      
## ******************************************************************** ##
## Make plots of predictors versus sensitivity measures
## ******************************************************************** ##

## Plot historgrams of two sensitivity measures
ggplot(melt(Hist_Sim_Sens_All_Mean),aes(x=value)) +
  geom_histogram() +
  facet_grid(variable~.)

## Plot xy scatter plot of two sensitivity measures
ggplot(Hist_Sim_Sens_All_Mean,aes(x=hist_sens,y=sim_sens)) +
  geom_point() +
  ggtitle('Mean annual sensitivity for\n historical and simulation data')

ggplot(Sim_Results_DensEff,aes(x=hist_sens,y=sim_sens_obsgrd)) +
  geom_point() +
  ggtitle('Mean annual sensitivity for\n historical and simulation data')

## Look at Population growth rate (eigen value) versus sensitivity values
ggplot(Sim_Results_DensEff,aes(x=EigenVal,y=hist_sens)) +
  geom_point()

ggplot(Sim_Results_DensEff,aes(x=EigenVal,y=sim_sens_obsgrd)) +
  geom_point()

## Look at fec.cv.avg - limit the x-range
ggplot(Sim_Results,aes(x=fec.cv.avg,y=hist_sens)) + 
  geom_point() +
  coord_cartesian(xlim=c(0,125))

## Look at histograms of final pop sizes
ggplot(Sim_Results,aes(x=exp.min.n)) +
  geom_histogram()

ggplot(Sim_Results,aes(x=n.mean)) +
  geom_histogram()

## Set predictors
predictors <- c(1,7,8,47,60,63,64,65)
names(Sim_Results_DensEff)[predictors]
## Set response
response <- 68

Sim_Results_Small <- Sim_Results[ c(predictors,response) ]
Sim_Results_Small_m <- melt(Sim_Results_Small,measure.vars=names(Sim_Results)[predictors])

ggplot(Sim_Results_Small_m,aes(x=value,y=hist_sens_obsgrd)) +
  geom_point() +
  facet_wrap(~variable,scales='free')



## ******************************************************************** ##
## Misc plots
## ******************************************************************** ##

## -------------------------------------------------------------------- ##
## create empty cumulative AOO loss category vector
loss_cat_char <- c()
## write function to create loss_cat_char
loss_cat_labeller <- function( loss_cat, loss_cat_char ){
  loss_cat_char[ loss_cat=="A" ] <- "< 2500"
  loss_cat_char[ loss_cat=="B" ] <- "2500-5000"
  loss_cat_char[ loss_cat=="C" ] <- "5000-10000"
  loss_cat_char[ loss_cat=="D" ] <- "> 10000"
  return(loss_cat_char)
}

Sim_Results_DensEff$loss_cat_char <- loss_cat_labeller( Sim_Results_DensEff$loss_cat, loss_cat_char )
Sim_Results_DensEff$loss_cat_char <- factor( Sim_Results_DensEff$loss_cat_char,
                                             levels=c("< 2500","2500-5000","5000-10000","> 10000") )

## -------------------------------------------------------------------- ##

## Mean historical sensitivity versus mean simulation sensitivity
ggplot(subset( Sim_Results_DensEff, subset=popd=='popd_LDD' ),
       aes(x=hist_sens_obsgrd,y=sim_sens_obsgrd,color=as.factor(sens_ppp_opt))) + 
  geom_point(size=2) +
  #geom_hline(yintercept=0.5) +
  #geom_vline(xintercept=0.5) +
  #geom_abline(intercept=(-0.1),slope=1) +
  xlab("Sensitivity") +
  ylab("Positive predictive power") +
  scale_color_manual(name="Combined metric\nvalue", values = c("grey", "black")) +
  theme_bw() +
  #annotate("text", x=0.65, y=0.7, label="Combined metric = 1", family="Times") +
  coord_fixed(ratio = 1) +
  theme(text=element_text(size=12,family="Times",face="bold"))   
#ggsave(filename='figures/Sensitivity_vs_PPP_popd_1k.pdf',width=7.5,height=7.5,units='in')
#ggsave(filename='figures/Diss_Fig_4_6.pdf',width=6.5,height=6.5,units='in')
ggsave(filename='figures/Figure_PPP_Sens.pdf',width=6.5,height=6.5,units='in')

## -------------------------------------------------------------------- ##
## Make figure for presentation
## ***
## Mean historical sensitivity versus mean simulation sensitivity
ggplot(subset( Sim_Results_DensEff, subset=popd=='popd_LDD' ),
       aes(x=hist_sens_obsgrd,y=sim_sens_obsgrd,colour=loss_cat_char)) + 
  geom_point(size=3) +
  #geom_hline(yintercept=0.5) +
  geom_vline(xintercept=0.5) +
  geom_abline(intercept=(-0.1),slope=1) +
  xlab("Sensitivity") +
  ylab("Positive predictive power") +
  scale_colour_discrete(name="Loss function\nvalue") +
  theme_bw() +
  #annotate("text", x=0.65, y=0.7, label="Combined metric = 1") +
  theme(text=element_text(size=29,face="bold"))   


## 2010 sensitivity versus mean historical sensitivity
ggplot(Sim_Results_DensEff,aes(x=hist_sens_obsgrd,y=sens_hist.2010,colour=loss_cat_char)) + 
  geom_point(size=3) +
  theme_bw()

## 2010 sensitivity versus positive predictive power
ggplot(Sim_Results_DensEff,aes(x=sens_hist.2010,y=sens_sim_obsgrd.2010,colour=loss_cat_char)) + 
  geom_point(size=3) +
  geom_hline(yintercept=0.47) +
  theme_bw()


ggplot(Sim_Results_DensEff, aes(x=sim_sens,y=hist_sens)) + 
  geom_point()

ggplot(Sim_Results_DensEff, aes(x=cum_AOO_loss,y=log(n.mean),colour=kch.type)) + 
  geom_point()

ggplot(Sim_Results_DensEff, aes(kch.type,n.mean)) + 
  geom_boxplot()

ggplot(Sim_Results_DensEff, aes(x=sim_sens_obsgrd,y=cum_AOO_loss)) + 
  geom_point()

ggplot(Sim_Results_DensEff, aes(x=hist_sens_obsgrd,y=cum_AOO_loss)) + 
  geom_point()



ggplot(Sim_Results_DensEff, aes(x=sens_hist.2010,y=cum_AOO_loss)) + 
  geom_point()

ggplot(Sim_Results_DensEff, aes(x=sim_sens,y=sim_sens_obsgrd)) + 
  geom_point()

ggplot(Sim_Results_DensEff, aes(x=NumTranslocs,y=cum_AOO_loss)) + 
  geom_point()

ggplot(Sim_Results,aes(x=sim_sens,y=sim_sens_obsgrd)) +
  geom_point() +
  geom_abline(intercept=0,slope=1,col='red')

ggplot(Sim_Results,aes(x=hist_sens,y=hist_sens_obsgrd)) +
  geom_point() +
  geom_abline(intercept=0,slope=1,col='red')
