###################### Table Of Contents #######################
# Search for figs.XX to jump to the corresponding plot section #
# In each Chapter there are sometimes several plots for the    #
# different variables etc. Some of the Parts are quite long    #
#                                                              #
# Figs.01: Count Series and Station Locations                  #
# Figs.02: plot sample time series for lon/lat                 #
# Figs.03: Time Series for indices (full period)               #
# Figs.04: Sample Year Plots                                   #
# Figs.05: Spread-Error Ratio Analysis                         #
# Figs.06: Talagrant Diagram                                   #
# Figs.07: Maps of Variable-Spread                             #
# Figs.08: Average Update Plots                                #
# Figs.09: Validation Correlation Maps                         #
# Figs.10: Correlation Difference Maps                         #
# Figs.11: Validation Bias Maps                                #
# Figs.12: Reduction of Error Maps                             #
# Figs.13: RootMeanSquareError Maps                            #
# Figs.14: Validation Indices TS                               #
# Figs.15: Validation Indices Correlations                     #
### Figs.16: Plot Distance Weight                                #
# Appendix: Commented Plot parts moved to here                 #
################################################################
rm(list=ls())

# syrtot and eyrtot are only used for the total indices time series (e.g. ENH.temp2, Global temp2 etc)
#syrtot=1903 #set to the same syr and eyr of the prepplots script (default 1602)
#eyrtot=1904 #(default 2000) 

syr=1902 #validation period: syr>=1902, eyr<2004. Syr should be the later of the two 1902 and syr in prepplots
eyr=2000

user <- system("echo $USER",intern=T)
print(paste('User:',user))
if (user=="veronika") {
  # workdir('/scratch/veronika/rerun/r_code')
  workdir ='/scratch3/veronika/reuse/reuse_git/' # where are the scripts from github
} else if (user=="joerg") {
  workdir='/scratch3/joerg/projects/reuse/git/'
} else {
  stop("Unknown user!")
  
}
dataextdir='/mnt/climstor/giub/EKF400/'
dataintdir=paste0(workdir,'../data/')
setwd(workdir)


source('EnSRF_switches.R')
source('EnSRF_functions.R')



if (monthly_out) {
  prepplotdir=paste0("../data/prepplot/EKF400_",version,"_",expname,"/prepplot_monthly/")
  nseas=12
  s.plot=12
} else if (pseudo_prox) {
  prepplotdir=paste0("../data/prepplot/EKF400_",version,"_",expname,"/prepplot_annual/") 
  nseas=1
  s.plot=1
} else {
  prepplotdir=paste0("../data/prepplot/EKF400_",version,"_",expname,"/prepplot_seasonal/") 
  nseas=2
  s.plot=2
}
figpath=paste0('../figures/EKF400_',version,'_',expname,'_',syr,'-',eyr) #format(Sys.time(), "%Y%m%d_%H-%M_")
dir.create(figpath)

pnames <- paste(c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k'), ')', sep='')


################################################################################
################## Figs.01: Count Series and Station Locations #################
################################################################################



################ Count Series Plot ####################

if (countseries & !pseudo_prox) {
  # count number of assimilation data of each type and add to plot
  nmxd <- ntrw <- nprox <- ndoc <-ninstslp <-ninsttemp <- ninst <- rep(NA,length(syr:eyr))
  i <- 1
  for (cyr in syr:eyr) {
    #print(cyr)
    if (every2grid) {
      load(file=paste0(prepplotdir,'analysis_',cyr,'_2ndgrid.Rdata')) 
    } else {
      load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata')) 
    }
    ninst[i] <- length(which(calibrate$sour=="inst"))
    ninsttemp[i] <- length(which(calibrate$sour=="inst"&calibrate$names=="temp2"))
    ninstslp[i] <- length(which(calibrate$sour=="inst"&calibrate$names=="slp"))
    ndoc[i] <- length(which(calibrate$sour=="doc"))
    nprox[i] <- length(which(calibrate$sour=="prox"))-length(which(apply(is.na(calibrate$mr[which(calibrate$sour=="prox"),]),1,all)))-length(which(apply(is.na(calibrate$data[which(calibrate$sour=="prox"),]),1,all)&apply(!is.na(calibrate$mr[which(calibrate$sour=="prox"),]),1,all)))
    i <- i+1
  }
  nrecords <- cbind((syr:eyr),ninst,ninsttemp,ninstslp,ndoc,nprox,rep(0,length(syr:eyr)))
  #nrecords[41,] <- nrecords[40,] <- nrecords[39,] <- nrecords[38,] # check input data 1639-41!
  pdf(paste(figpath,'record_no.pdf',sep='/'), width=9, height=3.5, paper='special')
  par(oma=c(0,0,0,0),mar=c(4,4,0.4,0.4))
  
  plot(nrecords[,1],nrecords[,7],ty="l",col='white',ylim=c(0,5000),xlab="year",ylab="No of records")
  polygon(c(nrecords[,1],rev(nrecords[,1])),c(nrecords[,6],rev(nrecords[,7])),col="seagreen3")
  polygon(c(nrecords[,1],rev(nrecords[,1])),c((nrecords[,5]+nrecords[,6]),rev(nrecords[,6])),col="plum2")
  polygon(c(nrecords[,1],rev(nrecords[,1])), c((nrecords[,3]+nrecords[,5]+nrecords[,6]),rev((nrecords[,5]+nrecords[,6]))),col="firebrick1")
  polygon(c(nrecords[,1],rev(nrecords[,1])), c((nrecords[,4]+nrecords[,3]+nrecords[,5]+nrecords[,6]),rev((nrecords[,3]+nrecords[,5]+nrecords[,6]))),col="goldenrod1")
  legend("topleft", c("Instr. SLP","Instr. temp.", "Docum. temp.",'Proxies'), 
         pch=rep(15,6), col=c("goldenrod1","firebrick1", "plum2", "seagreen3"), pt.cex=1, pt.lwd=1, 
         inset=0.005, bg='transparent', box.col='transparent', cex=1)
  dev.off()
  
  ################ Station Locations Plot ########################
  
  if (every2grid) {
    load(file=paste0(prepplotdir,'analysis_',station_yr,'_2ndgrid.Rdata')) 
  } else {
    load(file=paste0(prepplotdir,'analysis_',station_yr,'.Rdata')) 
  }
  dat <- echam
  dat$data <- array(echam$ensmean[1:4608,]*0, c(4608, 1, 2))
  dat$lat <- dat$lat[1:dim(dat$data)[1]]
  dat$lon <- dat$lon[1:dim(dat$data)[1]]
  dat$names <- dat$names[1:dim(dat$data)[1]]
  pdf(paste0(figpath,'/stat_locations_',station_yr,'.pdf'), width=9, height=3.5, paper='special')
  layout(matrix(c(1,2), 1, 2, byrow = TRUE))
  par(oma=c(0,0,0,0))
  levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0.05,0.2,0.4,0.6,0.8,1)
  plot_echam(dat,type='data',cex.pt=1.5,names=pnames[1:dim(dat$data)[3]],lev=levs,colorbar=F,
             stations=calibrate,st.col=NULL,statpanel=c(1,2),add=T,wcol='darkgrey')
  legend("bottomleft",c("kept","excluded"),pch=c(3,3),pt.lwd=1.5,col=c(rgb(0,10,0,8,maxColorValue=10),"firebrick4"))
  dev.off()
}




################################################################################
############################## Validation Plots ################################
################################################################################

for (v in validation_set){
  # loop over all validation datasets in validation_set, specified in the SWITCHES script.
  # 1. reads validation statistics file with everything related to the 20th century vali period
  # 2. reads the indices_tot file for plotting of time series over the entire e.g. 400-yr period
  if (every2grid) {
    if (monthly_out) {
      load(file=paste0("../data/image/EKF400_",version,"_",expname,"/prepplot_",v,"_calc_vali_stat_image_",
                       syr,"-",eyr,"_monthly_2ndgrid.Rdata"))
#      load(file=paste0("../data/image/EKF400_",version,"_",expname,"/indices_tot_",syrtot,"-",eyrtot,
#                       "_monthly_2ndgrid.Rdata"))
    } else if (pseudo_prox) {
      load(file=paste0("../data/image/EKF400_",version,"_",expname,"/prepplot_",v,"_calc_vali_stat_image_",
                       syr,"-",eyr,"_annual_2ndgrid.Rdata"))
#      load(file=paste0("../data/image/EKF400_",version,"_",expname,"/indices_tot_",syrtot,"-",eyrtot,
#                       "_annual_2ndgrid.Rdata"))
    } else {
      load(file=paste0("../data/image/EKF400_",version,"_",expname,"/prepplot_",v,"_calc_vali_stat_image_",
                       syr,"-",eyr,"_seasonal_2ndgrid.Rdata"))
#      load(file=paste0("../data/image/EKF400_",version,"_",expname,"/indices_tot_",syrtot,"-",eyrtot,
#                       "_seasonal_2ndgrid.Rdata"))
    }  
  } else {
    if (monthly_out) {
      load(file=paste("../data/image/EKF400_",version,"_",expname,"/prepplot_",v,"_calc_vali_stat_image_",
                      syr,"-",eyr,"_monthly.Rdata",sep=""))
#      load(file=paste0("../data/image/EKF400_",version,"_",expname,"/indices_tot_",syrtot,"-",eyrtot,
#                       "_monthly.Rdata"))
    } else if (pseudo_prox) {
      load(file=paste("../data/image/EKF400_",version,"_",expname,"/prepplot_",v,"_calc_vali_stat_image_",
                      syr,"-",eyr,"_annual.Rdata",sep=""))
#      load(file=paste0("../data/image/EKF400_",version,"_",expname,"/indices_tot_",syrtot,"-",eyrtot,
#                       "_annual.Rdata"))
    } else {
      load(file=paste("../data/image/EKF400_",version,"_",expname,"/prepplot_",v,"_calc_vali_stat_image_",
                      syr,"-",eyr,"_seasonal.Rdata",sep=""))
#      load(file=paste0("../data/image/EKF400_",version,"_",expname,"/indices_tot_",syrtot,"-",eyrtot,
#                       "_seasonal.Rdata"))
    }
  }
  #if (monthly_out) {
  #  stop("insert spos and epos for monthly plots! Not coded yet")
  #} else {
  #  spos <- which(vind[[1]]$time == syr)
  #  epos <- which(vind[[1]]$time == paste0(eyr,'.5'))
  #}
  validate.init <- validate
  validate.anom.init <- validate.anom
  vind.init <- vind
  vind.anom.init <- vind.anom
#  vind.tot.init <- vind.tot
#  vind.anom.tot.init <- vind.anom.tot
  
  validate <- validate.init[[v]]
  validate.anom<-validate.anom.init[[v]]
  vind<-vind.init[[v]]
  vind.anom<-vind.anom.init[[v]]
#  vind.tot<-vind.tot.init[[v]]
#  vind.anom.tot<-vind.anom.tot.init[[v]]
  
  data.dim <- dim(echam$data)[c(1,2,2,3)]
  data.dim[2:3] <- c(s.plot,data.dim[3]/s.plot)
  ens.dim <- c(nrow(echam$ensmean), s.plot, ncol(echam$ensmean)/s.plot)
  print(paste("monthly_out:",monthly_out))
  
  
  
  ###############################################################
  ####### figs.02: plot sample time series for lon/lat: #########
  ###############################################################
  
  if (v=="cru_vali"){
    valilegend<-"CRU"
  }else if(v=="twentycr_vali"){
    valilegend<-"20CR"
  }
  
  lonlatnr<-c(592,676,278,632,1460) #,2734)
  #look at plot_example_ts in function script.
  #codeline beneath finds gridbox with lowest correlation (was used because of strange bad results)
  #which(corr$ensmean[,2]==min(corr$ensmean[which(corr$lat<(-13)&corr$lat>(-18)&corr$lon<(170)&corr$lon>(165)&corr$names=="temp2"),2],na.rm=T))
  # vanatu -15.4492793 167.595411
  for (nr in lonlatnr){
    plot_example_ts(validate,analysis,echam,nr)
    plot_example_ts(validate.anom,analysis.anom,echam.anom,nr,type="anomaly")
  }
  
  ###############################################################
  
  
  # 
  # ###############################################################
  # ####### figs.03: Time Series for indices (full period) ########
  # ###############################################################
  # if (!pseudo_prox) {
  # inds <- c('ENH.temp2', 'EU.temp2', 'NEU.temp2', 'GLO.temp2')
  # indices <- c('ENH.temp2','NAM.temp2','SAM.temp2','AFR.temp2',
  #              'ASI.temp2','AUS.temp2','ARC.temp2','ANT.temp2',
  #              'NEU.temp2','MED.temp2',
  #              'GLO.temp2','NAM.precip','SAM.precip','AFR.precip',
  #              'ASI.precip','AUS.precip','ARC.precip','ANT.precip',
  #              'NEU.precip','MED.precip',
  #              'SH.temp2','NAM.slp','SAM.slp','AFR.slp',
  #              'ASI.slp','AUS.slp','ARC.slp','ANT.slp',
  #              'NEU.slp','MED.slp',
  #              'NH.temp2', 'EU.temp2', 'EU.precip', 'EU.slp')
  # vind.tot$names<-indices
  # for(ind in inds){
  #   pdf(paste(figpath,'/timeseries_',v,'_',ind,'.pdf',sep=''), width=15, height=4.5, paper='special')
  #     winter<-seq(1,length(eind.tot$time),by=2)
  #     summer<-seq(2,length(eind.tot$time), by=2)
  #     valiwinter<-seq(1,length(vind.tot$time),by=2)
  #     valisummer<-seq(2,length(vind.tot$time),by=2)
  #     # winter
  #     ylimmin<-min(cbind(c(vind.tot$ensmean[which(vind.tot$names==ind),valiwinter],rep(NA,syr-syrtot)), 
  #                        eind.tot$ensmean[which(eind.tot$names==ind),winter],
  #                        aind.tot$ensmean[which(aind.tot$names==ind),winter]),na.rm=T)-1
  #     ylimmax<-max(cbind(c(vind.tot$ensmean[which(vind.tot$names==ind),valiwinter],rep(NA,syr-syrtot)), 
  #                        eind.tot$ensmean[which(eind.tot$names==ind),winter],
  #                        aind.tot$ensmean[which(aind.tot$names==ind),winter]),na.rm=T)+1
  #     par(oma=c(0,0,0,0),mar=c(4,4,2,0.5),mfrow=c(1,2))
  #     plot(eind.tot$time[winter],eind.tot$ensmean[which(eind.tot$names==ind),winter], ylim=c(ylimmin,ylimmax),ty='l',col="black",main="Oct.-Mar.",
  #          xlab='year',ylab=paste(ind,'[°C]'))
  #     polygon(c(eind.tot$time[winter],rev(eind.tot$time[winter])),c(apply(eind.tot$data[which(eind.tot$names==ind),winter,],1,max),
  #             rev(apply(eind.tot$data[which(eind.tot$names==ind),winter,],1,min))),density=NA, col=rgb(1,1,1,3,maxColorValue=10))
  #     lines(aind.tot$time[winter],aind.tot$ensmean[which(aind.tot$names==ind),winter],col="red")
  #     polygon(c(aind.tot$time[winter],rev(aind.tot$time[winter])),c(apply(aind.tot$data[which(aind.tot$names==ind),winter,],1,max),
  #             rev(apply(aind.tot$data[which(aind.tot$names==ind),winter,],1,min))),density=NA, col=rgb(10,0,0,3,maxColorValue=10))
  #     lines(vind.tot$time[valiwinter], vind.tot$ensmean[which(vind.tot$names==ind),valiwinter],col="blue")
  #     legend("bottomleft", c(valilegend, 'CCC400', "EFK400"),col=c("blue", "black", "red"),
  #            lty=c(1,1,1),lwd=c(1,1,1),pt.cex=1, pt.lwd=1,inset=0.005, bg='transparent',
  #            box.col='transparent', cex=1)
  #     # summer
  #     ylimmin<-min(cbind(c(vind.tot$ensmean[which(vind.tot$names==ind),valisummer],rep(NA,syr-syrtot)), 
  #                        eind.tot$ensmean[which(eind.tot$names==ind),summer],
  #                        aind.tot$ensmean[which(aind.tot$names==ind),summer]),na.rm=T)-1
  #     ylimmax<-max(cbind(c(vind.tot$ensmean[which(vind.tot$names==ind),valisummer],rep(NA,syr-syrtot)), 
  #                        eind.tot$ensmean[which(eind.tot$names==ind),summer],
  #                        aind.tot$ensmean[which(aind.tot$names==ind),summer]),na.rm=T)+1
  #     plot(eind.tot$time[summer],eind.tot$ensmean[which(eind.tot$names==ind),summer], ylim=c(ylimmin,ylimmax),
  #          ty='l',col="black",main="Apr.-Sept.",xlab='year',ylab='')
  #     polygon(c(eind.tot$time[summer],rev(eind.tot$time[summer])),c(apply(eind.tot$data[which(eind.tot$names==ind),summer,],1,max),
  #             rev(apply(eind.tot$data[which(eind.tot$names==ind),summer,],1,min))),density=NA, col=rgb(1,1,1,3,maxColorValue=10))
  #     lines(aind.tot$time[summer],aind.tot$ensmean[which(aind.tot$names==ind),summer],col="red")
  #     polygon(c(aind.tot$time[summer],rev(aind.tot$time[summer])),c(apply(aind.tot$data[which(aind.tot$names==ind),summer,],1,max),
  #             rev(apply(aind.tot$data[which(aind.tot$names==ind),summer,],1,min))),density=NA, col=rgb(10,0,0,3,maxColorValue=10))
  #     lines(vind.tot$time[valisummer], vind.tot$ensmean[which(vind.tot$names==ind),valisummer],col="blue")
  #   dev.off()
  #   
  #   # indices timeserie whole periode smoothed by runningmean wdow = window of running mean
  #   wdow<-11
  #   vind.run<-runmean(vind.tot$ensmean[which(vind.tot$names==ind),], 2*wdow,endrule="NA")
  #   aind.run<-runmean(aind.tot$ensmean[which(aind.tot$names==ind),], 2*wdow,endrule="NA")
  #   eind.run<-runmean(eind.tot$ensmean[which(eind.tot$names==ind),], 2*wdow,endrule="NA")
  #   ylimmin<-min(rbind(vind.run,eind.run,aind.run),na.rm=T)-1
  #   ylimmax<-max(rbind(vind.run,eind.run,aind.run),na.rm=T)+1
  #   pdf(paste(figpath,'/timeseries_',v,'_',ind,'_smoothed.pdf',sep=''), width=15, height=4.5, paper='special')
  #     plot(eind.tot$time,eind.run, ylim=c(ylimmin,ylimmax),ty='l',col="black",
  #          xlab='year',ylab=paste(ind,'[°C]'), main=paste(wdow,'year running mean'))
  #     lines(vind.tot$time, vind.run, col="blue")
  #     lines(aind.tot$time, aind.run, col="red")
  #     legend("bottomleft", c(valilegend, 'CCC400', "EFK400"),col=c("blue", "black", "red"),
  #            lty=c(1,1,1),lwd=c(1,1,1),pt.cex=1, pt.lwd=1,inset=0.005, bg='transparent',
  #            box.col='transparent', cex=1)
  #   dev.off()
  # }
  # 
  # # same but with anomalies
  # vind.anom.tot$names<-indices
  # for(ind in inds){
  #   pdf(paste(figpath,'/timeseries_',v,'_',ind,'_anom.pdf',sep=''), width=15, height=4.5, paper='special')
  #     winter<-seq(1,length(eind.anom.tot$time),by=2)
  #     summer<-seq(2,length(eind.anom.tot$time), by=2)
  #     valiwinter<-seq(1,length(vind.anom.tot$time),by=2)
  #     valisummer<-seq(2,length(vind.anom.tot$time),by=2)
  #     ylimmin<-min(cbind(c(vind.anom.tot$ensmean[which(vind.anom.tot$names==ind),valiwinter],rep(NA,syr-syrtot)), 
  #                        eind.anom.tot$ensmean[which(eind.anom.tot$names==ind),winter],
  #                        aind.anom.tot$ensmean[which(aind.anom.tot$names==ind),winter]),na.rm=T)-1
  #     ylimmax<-max(cbind(c(vind.anom.tot$ensmean[which(vind.anom.tot$names==ind),valiwinter],rep(NA,syr-syrtot)), 
  #                        eind.anom.tot$ensmean[which(eind.anom.tot$names==ind),winter],
  #                        aind.anom.tot$ensmean[which(aind.anom.tot$names==ind),winter]),na.rm=T)+1
  #     par(oma=c(0,0,0,0),mar=c(4,4,2,0.5),mfrow=c(1,2))
  #     # winter
  #     plot(eind.anom.tot$time[winter],eind.anom.tot$ensmean[which(eind.anom.tot$names==ind),winter], 
  #          ylim=c(ylimmin,ylimmax),ty='l',col="black",main="Oct.-Mar.",xlab='year',ylab=paste(ind,'[°C]'))
  #     polygon(c(eind.anom.tot$time[winter],rev(eind.anom.tot$time[winter])),
  #             c(apply(eind.anom.tot$data[which(eind.anom.tot$names==ind),winter,],1,max),
  #             rev(apply(eind.anom.tot$data[which(eind.anom.tot$names==ind),winter,],1,min))),
  #             density=NA, col=rgb(1,1,1,3,maxColorValue=10))
  #     lines(aind.anom.tot$time[winter],aind.anom.tot$ensmean[which(aind.anom.tot$names==ind),winter],col="red")
  #     polygon(c(aind.anom.tot$time[winter],rev(aind.anom.tot$time[winter])),
  #             c(apply(aind.anom.tot$data[which(aind.anom.tot$names==ind),winter,],1,max),
  #             rev(apply(aind.anom.tot$data[which(aind.anom.tot$names==ind),winter,],1,min))),
  #             density=NA, col=rgb(10,0,0,3,maxColorValue=10))
  #     lines(vind.anom.tot$time[valiwinter], vind.anom.tot$ensmean[which(vind.anom.tot$names==ind),valiwinter],col="blue")
  #     legend("bottomleft", c(valilegend, 'CCC400', "EFK400"),col=c("blue", "black", "red"),
  #            lty=c(1,1,1),lwd=c(1,1,1),pt.cex=1, pt.lwd=1,inset=0.005, bg='transparent',
  #            box.col='transparent', cex=1)
  #     # summer
  #     ylimmin<-min(cbind(c(vind.anom.tot$ensmean[which(vind.anom.tot$names==ind),valisummer],
  #                  rep(NA,syr-syrtot)), eind.anom.tot$ensmean[which(eind.anom.tot$names==ind),summer],
  #                  aind.anom.tot$ensmean[which(aind.anom.tot$names==ind),summer]),na.rm=T)-1
  #     ylimmax<-max(cbind(c(vind.anom.tot$ensmean[which(vind.anom.tot$names==ind),valisummer],
  #                  rep(NA,syr-syrtot)), eind.anom.tot$ensmean[which(eind.anom.tot$names==ind),summer],
  #                  aind.anom.tot$ensmean[which(aind.anom.tot$names==ind),summer]),na.rm=T)+1
  #     plot(eind.anom.tot$time[summer],eind.anom.tot$ensmean[which(eind.anom.tot$names==ind),summer], ylim=c(ylimmin,ylimmax),ty='l',col="black",main="Apr.-Sept.",
  #          xlab='year',ylab='')
  #     polygon(c(eind.anom.tot$time[summer],rev(eind.anom.tot$time[summer])),
  #             c(apply(eind.anom.tot$data[which(eind.anom.tot$names==ind),summer,],1,max),
  #             rev(apply(eind.anom.tot$data[which(eind.anom.tot$names==ind),summer,],1,min))),
  #             density=NA,col=rgb(1,1,1,3,maxColorValue=10))
  #     lines(aind.anom.tot$time[summer],aind.anom.tot$ensmean[which(aind.anom.tot$names==ind),summer],col="red")
  #     polygon(c(aind.anom.tot$time[summer],rev(aind.anom.tot$time[summer])),
  #             c(apply(aind.anom.tot$data[which(aind.anom.tot$names==ind),summer,],1,max),
  #             rev(apply(aind.anom.tot$data[which(aind.anom.tot$names==ind),summer,],1,min))),
  #             density=NA,col=rgb(10,0,0,3,maxColorValue=10))
  #     lines(vind.anom.tot$time[valisummer],vind.anom.tot$ensmean[which(vind.anom.tot$names==ind),valisummer],
  #           col="blue")
  #   dev.off()
  #   
  #   
  #   # indices timeserie whole periode smoothed by runningmean wdow = window of running mean
  #   wdow<-11
  #   vind.run<-runmean(vind.anom.tot$ensmean[which(vind.anom.tot$names==ind),], 2*wdow,endrule="NA")
  #   aind.run<-runmean(aind.anom.tot$ensmean[which(aind.anom.tot$names==ind),], 2*wdow,endrule="NA")
  #   eind.run<-runmean(eind.anom.tot$ensmean[which(eind.anom.tot$names==ind),], 2*wdow,endrule="NA")
  #   ylimmin<-min(rbind(vind.run, eind.run,aind.run),na.rm=T)-1
  #   ylimmax<-max(rbind(vind.run, eind.run,aind.run),na.rm=T)+1
  #   pdf(paste(figpath,'/timeseries_',v,'_',ind,'_anom_smoothed.pdf',sep=''), width=15, height=4.5, paper='special')
  #     plot(eind.anom.tot$time,eind.run, ylim=c(ylimmin,ylimmax),ty='l',col="black",
  #          xlab='year',ylab=paste(ind,'[°C]'), main=paste(wdow,'year running mean'))
  #     lines(vind.anom.tot$time, vind.run, col="blue")
  #     lines(aind.anom.tot$time, aind.run, col="red")
  #     legend("bottomleft", c(valilegend, 'CCC400', "EFK400"),col=c("blue", "black", "red"),
  #            lty=c(1,1,1),lwd=c(1,1,1),pt.cex=1, pt.lwd=1,inset=0.005, bg='transparent',
  #            box.col='transparent', cex=1)
  #   dev.off()
  # }
  # ###############################################################
  # } # end if not pseudo_prox
  # 
  
  ###############################################################
  ################# figs.04: Sample Year Plots ##################
  ###############################################################
  
  #chose timestep for sample year plots
  t=2
  if (!monthly_out&!yearly_out) {
    # Nevin: in the Appendix there is a code part to make a for loop over several yrs which includes all plots in one pdf
    plotdata=echam
    plotdata$data <- array(c(echam.anom$ensmean[,t],analysis.anom$ensmean[,t]), c(nrow(echam.anom$ensmean),1,2))
    #plotdata$names <- c(echam.anom$names,analysis.anom$names)
    #plotdata$data <- array(c(echam$ensmean[,t],analysis$ensmean[,t]), c(nrow(echam$ensmean),1,2))
    
    pdf(paste(figpath,'/example_year_',plotdata$time[t],'_temp_anom.pdf',sep=''), width=9, height=4.5, paper='special')
    layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
    par(oma=c(0,0,0,0))
    levs <- c(-Inf, seq(-1.2,1.2,0.2), Inf)
    plot_echam(plotdata, varname='temp2', type='data', cex.pt=1.5, names=pnames[1:dim(plotdata$data)[3]], lev=levs, st.col=NULL, stations=NULL, add=T)
    dev.off()
    
    pdf(paste(figpath,'/example_year_',plotdata$time[t],'_precip_anom.pdf',sep=''), width=9, height=4.5, paper='special')
    layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
    par(oma=c(0,0,0,0))
    levs <- c(-Inf, seq(-16,16,2), Inf)
    plot_echam(plotdata, varname='precip', type='data', cex.pt=1.5, names=pnames[1:dim(plotdata$data)[3]], lev=levs, st.col=NULL, stations=NULL, add=T)
    dev.off()
    
    pdf(paste(figpath,'/example_year_',plotdata$time[t],'_slp_anom.pdf',sep=''), width=9, height=4.5, paper='special')
    layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
    par(oma=c(0,0,0,0))
    levs <- c(-Inf, seq(-1.2,1.2,0.2), Inf)
    plot_echam(plotdata, varname='slp', type='data', cex.pt=1.5, names=pnames[1:dim(plotdata$data)[3]], lev=levs, st.col=NULL, stations=NULL, add=T)
    dev.off()
  }
  
  ###############################################################
  
  
  
  ###############################################################
  ########### figs.05: Spread-Error Ratio #######################
  ###############################################################
  
  
  #################### SER Analysis ##############################
  
  if (!monthly_out){
    plotdata=echam
    plotdata$data <- array(cbind(sprerr.win,sprerr.sum,sprerr.c.win,sprerr.c.sum), 
                           c(length(sprerr.c.win), 1, 4))
    plotdata$ensmean <- array(cbind(sprerr.win,sprerr.sum,sprerr.c.win,sprerr.c.sum), 
                              c(length(sprerr.c.win), 4))
    plotdata$names <- plotdata$names[echam$names=="temp2"]
    plotdata$lon <- plotdata$lon[echam$names=="temp2"]
    plotdata$lat <- plotdata$lat[echam$names=="temp2"]
    
    levs <- c(0,0.33,0.5,0.67,0.83,1.2,1.5,2,3,Inf)
    
    plot_echam4(plotdata, varname='temp2', cex.pt=1.5, names=pnames[1:dim(plotdata$data)[3]], lev=levs , 
                type='ensmean', st.col=NULL, stations=plstat,NHseason=NULL,
                plotname=paste0('spread_error_ratio_anom_tmean_analysis.pdf'), paper='special')
    
    
    ##################### SER Echam ###########################
    plotdata=echam
    plotdata$data <- array(cbind(ech.sprerr.win,ech.sprerr.sum,ech.sprerr.c.win,ech.sprerr.c.sum), 
                           c(length(ech.sprerr.c.win), 1, 4))
    plotdata$ensmean <- array(cbind(ech.sprerr.win,ech.sprerr.sum,ech.sprerr.c.win,ech.sprerr.c.sum), 
                              c(length(ech.sprerr.c.win), 4))
    plotdata$names <- plotdata$names[echam$names=="temp2"]
    plotdata$lon <- plotdata$lon[echam$names=="temp2"]
    plotdata$lat <- plotdata$lat[echam$names=="temp2"]
    levs <- c(0,0.33,0.5,0.67,0.83,1.2,1.5,2,3,Inf)
    plot_echam4(plotdata, varname='temp2', cex.pt=1.5, names=pnames[1:dim(plotdata$data)[3]], lev=levs , 
                type='ensmean', st.col=NULL, stations=plstat,NHseason=NULL,
                plotname=paste0('spread_error_ratio_anom_tmean_echam.pdf'), paper='special')
  }
  ###############################################################
  
  
  
  ###############################################################
  ################# figs.06: Talagrant Diagram ##################
  ###############################################################
  if (vali) {
    if (!recon_vali) {
      if (anomaly_assim) {
        ereliable.summer <- ereliable.anom.summer
        areliable.summer <- areliable.anom.summer
        # erel_obserr <- erel_obserr.anom
        # arel_obserr <- arel_obserr.anom
      }
      if (length(names(ereliable$temp2)) == 31) {
        pdf(paste(figpath,'talagrant_temp_summer.pdf',sep='/'), width=6, height=6, paper='special') 
        par(mfrow=c(2,1), mar=c(3,3,1,1), oma=c(0,0,0,0))
        names(ereliable.summer$temp2) <- seq(0,30)
        barplot(ereliable.summer$temp2)
        names(areliable.summer$temp2) <- seq(0,30)
        barplot(areliable.summer$temp2)
        dev.off()
      }
      if (length(names(ereliable$precip)) == 31) {
        pdf(paste(figpath,'talagrant_precip_summer.pdf',sep='/'), width=6, height=6, paper='special') 
        par(mfrow=c(2,1), mar=c(3,3,1,1), oma=c(0,0,0,0))
        names(ereliable.summer$precip) <- seq(0,30)
        barplot(ereliable.summer$precip)
        names(areliable.summer$precip) <- seq(0,30)
        barplot(areliable.summer$precip)
        dev.off()
      }
      if (length(names(ereliable$slp)) == 31) {
        pdf(paste(figpath,'talagrant_slp_summer.pdf',sep='/'), width=6, height=6, paper='special') 
        par(mfrow=c(2,1), mar=c(3,3,1,1), oma=c(0,0,0,0))
        names(ereliable.summer$slp) <- seq(0,30)
        barplot(ereliable.summer$slp)
        names(areliable.summer$slp) <- seq(0,30)
        barplot(areliable.summer$slp)
        dev.off()
      }
      # if (length(names(erel_obserr$temp2)) == 31) {
      #   pdf(paste(figpath,'talagrant_obserr_temp.pdf',sep='/'), width=6, height=6, paper='special') 
      #   par(mfrow=c(2,1), mar=c(3,3,1,1), oma=c(0,0,0,0))
      #   names(erel_obserr$temp2) <- seq(0,30)
      #   barplot(erel_obserr$temp2)
      #   names(arel_obserr$temp2) <- seq(0,30)
      #   barplot(arel_obserr$temp2)
      #   dev.off()
      # }
    }
  }
  
  if (vali & !monthly_out) {
    if (!recon_vali) {
      if (anomaly_assim) {
        ereliable <- ereliable.anom
        areliable <- areliable.anom
        erel_obserr <- erel_obserr.anom
        arel_obserr <- arel_obserr.anom
      }
      
      
      # ereliable <- tapply(apply(echam$data[1:(dim(validate$data)[1]),,] > 
      #                             as.vector(validate$data),1:2,mean), rep(echam$names,
      #                             length=length(validate$data)), table)
      # areliable <- tapply(apply(analysis$data[1:(dim(validate$data)[1]),,] > 
      #                             as.vector(validate$data),1:2,mean), rep(analysis$names,
      #                             length=length(validate$data)), table)
      if (length(names(ereliable$temp2)) == 31) {
        pdf(paste(figpath,'talagrant_temp.pdf',sep='/'), width=6, height=6, paper='special') 
        par(mfrow=c(2,1), mar=c(3,3,1,1), oma=c(0,0,0,0))
        names(ereliable$temp2) <- seq(0,30)
        barplot(ereliable$temp2)
        names(areliable$temp2) <- seq(0,30)
        barplot(areliable$temp2)
        dev.off()
      }
      if (length(names(ereliable$precip)) == 31) {
        pdf(paste(figpath,'talagrant_precip.pdf',sep='/'), width=6, height=6, paper='special') 
        par(mfrow=c(2,1), mar=c(3,3,1,1), oma=c(0,0,0,0))
        names(ereliable$precip) <- seq(0,30)
        barplot(ereliable$precip)
        names(areliable$precip) <- seq(0,30)
        barplot(areliable$precip)
        dev.off()
      }
      if (length(names(ereliable$slp)) == 31) {
        pdf(paste(figpath,'talagrant_slp.pdf',sep='/'), width=6, height=6, paper='special') 
        par(mfrow=c(2,1), mar=c(3,3,1,1), oma=c(0,0,0,0))
        names(ereliable$slp) <- seq(0,30)
        barplot(ereliable$slp)
        names(areliable$slp) <- seq(0,30)
        barplot(areliable$slp)
        dev.off()
      }
      if (length(names(erel_obserr$temp2)) == 31) {
        pdf(paste(figpath,'talagrant_obserr_temp.pdf',sep='/'), width=6, height=6, paper='special') 
        par(mfrow=c(2,1), mar=c(3,3,1,1), oma=c(0,0,0,0))
        names(erel_obserr$temp2) <- seq(0,30)
        barplot(erel_obserr$temp2)
        names(arel_obserr$temp2) <- seq(0,30)
        barplot(arel_obserr$temp2)
        dev.off()
      }
      # next lines are commented because we don't have uncertainty estimates for instr. prec. and slp
      # if (length(names(erel_obserr$precip)) == 31) {
      #   pdf(paste(figpath,'talagrant_obserr_precip.pdf',sep='/'), width=6, height=6, paper='special') 
      #   par(mfrow=c(2,1), mar=c(3,3,1,1), oma=c(0,0,0,0))
      #   names(erel_obserr$precip) <- seq(0,30)
      #   barplot(erel_obserr$precip)
      #   names(arel_obserr$precip) <- seq(0,30)
      #   barplot(arel_obserr$precip)
      #   dev.off()
      # }
      # if (length(names(erel_obserr$slp)) == 31) {
      #   pdf(paste(figpath,'talagrant_obserr_slp.pdf',sep='/'), width=6, height=6, paper='special') 
      #   par(mfrow=c(2,1), mar=c(3,3,1,1), oma=c(0,0,0,0))
      #   names(erel_obserr$slp) <- seq(0,30)
      #   barplot(erel_obserr$slp)
      #   names(arel_obserr$slp) <- seq(0,30)
      #   barplot(arel_obserr$slp)
      #   dev.off()
      # }
    }
  }
  ###############################################################
  
  
  
  ###############################################################
  ############### figs.07: Maps of Variable-Spread ##############
  ###############################################################
  
  
  ############## Temperature Spread ##############
  
  if(!monthly_out){
    espread <- echam
    espread$data <- array(ech.spread, c(nrow(ech.spread), 1, ncol(ech.spread)))
    espread$ensmean <- array(ech.spread, c(nrow(ech.spread),ncol(ech.spread)))
    
    aspread <- echam
#    if (pseudoproxy) { 
#      ana.spread.bkp <- ana.spread
#      ana.spread <- ana.spread.bkp[4] 
#    }
    #if ((instrumental) || (inst_at_proxy) || (real_proxies)) {ana.spread <- ana.spread.bkp[2] }
    #aspread$data <- array(unlist(ana.spread)/as.vector(ech.spread)*100, c(nrow(ech.spread), 1, ncol(ech.spread)*(length(ana.spread))))
    aspread$data <- array(unlist(ana.spread)/as.vector(ech.spread)*100, c(nrow(ech.spread), 1, ncol(ech.spread)))
    aspread$ensmean <- array(unlist(ana.spread)/as.vector(ech.spread)*100, c(nrow(ech.spread), ncol(ech.spread)))
    
    ## Luca: couldn't convert this part into plot_echam4 yet
    
    pdf(paste(figpath,'spread_temp.pdf',sep='/'), width=9, height=6, paper='special')
    oldpar <- par(no.readonly=TRUE)
    layout(matrix(c(1,2,3,3,4,5,6,6), 4, 2, byrow = TRUE), height=c(3,1,3,1))
    par(oma=c(0,0,0,0))
    plot_echam(espread, varname='temp2', names=pnames[1:dim(espread$data)[3]], symmetric=FALSE, 
               cex.pt=1.3, add=TRUE)
    plot_echam(aspread, varname='temp2', names=pnames[dim(espread$data)[3] + 1:dim(aspread$data)[3]], 
               lev=c(seq(0,90,10),101), cex.pt=1.3, add=TRUE, st.col=NULL, stations=plstat) #calibrate)
    
    dev.off()
  }
  
  
  
  ############## Precipitation Spread ##############
  
  ## Luca: couldn't convert this to plot_echam4 yet
  if (!monthly_out){
    
    pdf(paste(figpath,'spread_precip.pdf',sep='/'), width=9, height=6, paper='special')
    oldpar <- par(no.readonly=TRUE)
    layout(matrix(c(1,2,3,3,4,5,6,6), 4, 2, byrow = TRUE), height=c(3,1,3,1))
    par(oma=c(0,0,0,0))
    plot_echam(espread, varname='precip', names=pnames[1:dim(espread$data)[3]], symmetric=FALSE, 
               cex.pt=1.3, add=TRUE, levs=c(0,1,2,3,5,8,13,22,36,60,100))
    plot_echam(aspread, varname='precip', names=pnames[dim(espread$data)[3] + 1:dim(aspread$data)[3]], 
               lev=seq(0,100,10), cex.pt=1.3, add=TRUE, st.col=NULL, stations=plstat)
    dev.off()
  }
  
  ############## Sea Level Pressure ##############
  
  ## Luca: couldn't convert this to plot_echam4 yet
  if (!monthly_out){
    pdf(paste(figpath,'spread_slp.pdf',sep='/'), width=9, height=6, paper='special')
    oldpar <- par(no.readonly=TRUE)
    layout(matrix(c(1,2,3,3,4,5,6,6), 4, 2, byrow = TRUE), height=c(3,1,3,1))
    par(oma=c(0,0,0,0))
    plot_echam(espread, varname='slp', names=pnames[1:dim(espread$data)[3]], symmetric=FALSE, cex.pt=1.3, add=TRUE)
    plot_echam(aspread, varname='slp', names=pnames[dim(espread$data)[3] + 1:dim(aspread$data)[3]], 
               lev=c(seq(0,90,10),101), cex.pt=1.3, add=TRUE, st.col=NULL, stations=plstat)
    dev.off()
  }
  ###############################################################
  
  
  
  ###############################################################
  ############### figs.08: Average Update Plots #################
  ###############################################################
  
  
  ############## Temperature Update ##############
  
  upd <- apply(array(analysis$data - echam$data, data.dim), 1:2, mean, na.rm=T)
  update <- echam
  update$data <- array(upd,c(dim(upd)[1],1,dim(upd)[2]))
  update$ensmean <- array(update$data,dim=c(dim(update$data)[1],dim(update$data)[3]))
  tmax=round(max(abs(update$ensmean[which(update$names=="temp2"),])),digits=2)
  slpmax=round(max(abs(update$ensmean[which(update$names=="slp"),])),digits=2)
  prmax=round(max(abs(update$ensmean[which(update$names=="precip"),])),digits=1)
  if (monthly_out) {
    plot_echam4(update, varname='temp2', cex.pt=1.3, names=pnames[1:dim(update$data)[3]], type='ensmean', 
                st.col=NULL, stations=plstat,NHseason="summer",
                plotname=paste0('avg_upd_temp_summer.pdf'), paper='special')
    plot_echam4(update, varname='temp2', cex.pt=1.3, names=pnames[1:dim(update$data)[3]], type='ensmean', 
                st.col=NULL, stations=plstat,NHseason="winter",
                plotname=paste0('avg_upd_temp_winter.pdf'), paper='special')
  } else {
    plot_echam4(update, varname='temp2', cex.pt=1.3, lev=round(tmax/9*c(-Inf,-7.5:-0.5,0.5:7.5,Inf),digits=3), 
                names=pnames[1:dim(update$data)[3]], type='ensmean', st.col=NULL, stations=plstat,
                NHseason=NULL, plotname=paste0('avg_upd_temp.pdf'), paper='special')
  }
  
  ############## Precipitation Update ##############
  
  if (monthly_out) {
    plot_echam4(update, varname='precip', cex.pt=1.3, names=pnames[1:dim(update$data)[3]], type='ensmean', 
                st.col=NULL, stations=plstat, NHseason="summer",
                plotname=paste0('avg_upd_precip_summer.pdf'), paper='special')
    plot_echam4(update, varname='precip', cex.pt=1.3, names=pnames[1:dim(update$data)[3]], type='ensmean', 
                st.col=NULL, stations=plstat, NHseason="winter",
                plotname=paste0('avg_upd_precip_winter.pdf'), paper='special')
  } else {
    plot_echam4(update, varname='precip', cex.pt=1.3,lev=round(prmax/9*c(-Inf,-7.5:-0.5,0.5:7.5,Inf),digits=2), 
                names=pnames[1:dim(update$data)[3]], type='ensmean', st.col=NULL, stations=plstat,NHseason=NULL,
                plotname=paste0('avg_upd_precip.pdf'), paper='special')
  }
  
  ############## Sea Level Pressure ##############
  
  if (monthly_out) {
    plot_echam4(update, varname='slp', cex.pt=1.3, names=pnames[1:dim(update$data)[3]] , type='ensmean', 
                st.col=NULL, stations=plstat, NHseason="summer",
                plotname=paste0('avg_upd_slp_summer.pdf'), paper='special')
    plot_echam4(update, varname='slp', cex.pt=1.3, names=pnames[1:dim(update$data)[3]] , type='ensmean', 
                st.col=NULL, stations=plstat, NHseason="winter",
                plotname=paste0('avg_upd_slp_winter.pdf'), paper='special')
  } else {
    plot_echam4(update, varname='slp', cex.pt=1.3,lev=round(slpmax/9*c(-Inf,-7.5:-0.5,0.5:7.5,Inf),digits=3), 
                names=pnames[1:dim(update$data)[3]], type='ensmean', st.col=NULL, stations=plstat,
                NHseason=NULL,plotname=paste0('avg_upd_slp.pdf'), paper='special')
  }
  rm(upd,update)
  ###############################################################
  
  
  
  ###############################################################
  ########### figs.09: Validation Correlation Maps  #############
  ###############################################################
  
  
  
  ###################### Temperature Correlation ################
  
  corr.tot <- echam
  #corr.tot$data <- array(cbind(corr.ech$ensmean,corr$ensmean),c(dim(corr$ensmean)[1],1,dim(corr$ensmean)[2]*2)) #[,,1:2,drop=F]
  ## it doesn't make sense to put $ensmean into $data
  ## and the choosing type='data'
  ## different approach: put $ensmean into $ensmean and 
  ## choose type='ensmean'
  #corr.tot$ensmean <- array(cbind(corr.ech$ensmean,corr$ensmean),c(dim(corr$ensmean)[1],dim(corr$ensmean)[2]*2)) #[,,1:2,drop=F]
  corr.tot$ensmean <- array(cbind(corr.ech$ensmean,corr$ensmean),c(dim(corr$ensmean)[1],s.plot*2)) #[,,1:2,drop=F]
  
  levs <- c(-1,seq(-0.9,0.9,0.2),1)
  if (monthly_out) {
    plot_echam4(corr.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(corr.tot$data)[3]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat, NHseason="summer",
                plotname=paste0('corr_echam_anal-',v,'_temp_summer.pdf'), paper='special')
    plot_echam4(corr.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(corr.tot$data)[3]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat, NHseason="winter",
                plotname=paste0('corr_echam_anal-',v,'_temp_winter.pdf'), paper='special')
  } else {
    plot_echam4(corr.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(corr.tot$data)[3]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat, NHseason=NULL,
                plotname=paste0('corr_echam_anal-',v,'_temp.pdf'), paper='special')
  }
  
  #################### Precipitation Correlation ##################
  
  levs <- c(-1,seq(-0.9,0.9,0.2),1)
  if (monthly_out){
    plot_echam4(corr.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(corr.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason="summer",plotname=paste0('corr_echam_anal-',v,'_precip_summer.pdf'), paper='special')
    plot_echam4(corr.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(corr.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason="winter",plotname=paste0('corr_echam_anal-',v,'_precip_winter.pdf'), paper='special')
  }else{
    plot_echam4(corr.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(corr.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason=NULL,plotname=paste0('corr_echam_anal-',v,'_precip.pdf'), paper='special')
  }
  
  
  ################# Sea Level Pressure Correlation #################
  
  levs <- c(-1,seq(-0.9,0.9,0.2),1)
  
  if (monthly_out){
    plot_echam4(corr.tot, varname='slp', cex.pt=1.5, names=pnames[1:dim(corr.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason="summer",plotname=paste0('corr_echam_anal-',v,'_slp_summer.pdf'), paper='special')
    plot_echam4(corr.tot, varname='slp', cex.pt=1.5, names=pnames[1:dim(corr.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason="winter",plotname=paste0('corr_echam_anal-',v,'_slp_winter.pdf'), paper='special')
  }else{
    plot_echam4(corr.tot, varname='slp', cex.pt=1.5, names=pnames[1:dim(corr.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason=NULL,plotname=paste0('corr_echam_anal-',v,'_slp.pdf'), paper='special')
  }
  
  
  ###############################################################
  
  
  
  ###############################################################
  ########### figs.10: Correlation Difference Maps  #############
  ###############################################################
  
  corr.diff <- corr.tot
  if (pseudo_prox) {
    corr.diff$ensmean <- array(corr.diff$ensmean[,2]-corr.diff$ensmean[,1],dim=c(dim(corr.diff$ensmean)[1],1))
  } else {
    corr.diff$data <- corr.diff$data[,,(ncol(corr.tot$ensmean)/2+1):ncol(corr.tot$ensmean)]-corr.diff$data[,,1:(ncol(corr.tot$ensmean)/2)]
    corr.diff$ensmean <- corr.diff$ensmean[,(ncol(corr.tot$ensmean)/2+1):ncol(corr.tot$ensmean)]-corr.diff$ensmean[,1:(ncol(corr.tot$ensmean)/2)]
  }
  
  ###################### Temperature CorrDiff ###################
  
  if (monthly_out) {
    plot_echam4(corr.diff, varname='temp2', cex.pt=1.5, names=pnames[1:dim(corr.diff$data)[2]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat,NHseason="summer",
                plotname=paste0('corr.diff_echam_anal-',v,'_temp_summer.pdf'), paper='special')
    plot_echam4(corr.diff, varname='temp2', cex.pt=1.5, names=pnames[1:dim(corr.diff$data)[2]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat,NHseason="winter",
                plotname=paste0('corr.diff_echam_anal-',v,'_temp_winter.pdf'), paper='special')
  } else if (pseudo_prox) {
    plot_echam4(corr.diff, varname='temp2', cex.pt=1.5, names=pnames[1], lev=levs, type='ensmean', 
                st.col=NULL, stations=plstat,NHseason=NULL,
                plotname=paste0('corr.diff_echam_anal-',v,'_temp.pdf'), paper='special')
  } else {
    plot_echam4(corr.diff, varname='temp2', cex.pt=1.5, names=pnames[1:dim(corr.diff$data)[2]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat, NHseason=NULL,
                plotname=paste0('corr.diff_echam_anal-',v,'_temp.pdf'), paper='special')
  }
  
  
  ###################### Precipitation CorrDiff ################
  
  if (monthly_out) {
    plot_echam4(corr.diff, varname='precip', cex.pt=1.5, names=pnames[1:dim(corr.diff$data)[2]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat,NHseason="summer",
                plotname=paste0('corr.diff_echam_anal-',v,'_precip_summer.pdf'), paper='special')
    plot_echam4(corr.diff, varname='precip', cex.pt=1.5, names=pnames[1:dim(corr.diff$data)[2]], 
                lev=levs, type='ensmean', st.col=NULL, stations=plstat,NHseason="winter",
                plotname=paste0('corr.diff_echam_anal-',v,'_precip_winter.pdf'), paper='special')
  } else if (pseudo_prox) {
    plot_echam4(corr.diff, varname='precip', cex.pt=1.5, names=pnames[1], lev=levs, type='ensmean', 
                st.col=NULL, stations=plstat,NHseason=NULL,
                plotname=paste0('corr.diff_echam_anal-',v,'_precip.pdf'), paper='special')
  } else {
    plot_echam4(corr.diff, varname='precip', cex.pt=1.5, names=pnames[1:dim(corr.diff$data)[2]], 
                lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason=NULL,
                plotname=paste0('corr.diff_echam_anal-',v,'_precip.pdf'), paper='special')
  }
  
  ################# Sea Level Pressure CorrDiff #################
  
  if (monthly_out) {
    plot_echam4(corr.diff, varname='slp', cex.pt=1.5, names=pnames[1:dim(corr.diff$data)[2]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat,NHseason="summer",
                plotname=paste0('corr.diff_echam_anal-',v,'_slp_summer.pdf'), paper='special')
    plot_echam4(corr.diff, varname='slp', cex.pt=1.5, names=pnames[1:dim(corr.diff$data)[2]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat,NHseason="winter",
                plotname=paste0('corr.diff_echam_anal-',v,'_slp_winter.pdf'), paper='special')
  } else if (pseudo_prox) {
    plot_echam4(corr.diff, varname='slp', cex.pt=1.5, names=pnames[1], lev=levs, type='ensmean', 
                st.col=NULL, stations=plstat,NHseason=NULL,
                plotname=paste0('corr.diff_echam_anal-',v,'_slp.pdf'), paper='special')
  } else {
    plot_echam4(corr.diff, varname='slp', cex.pt=1.5, names=pnames[1:dim(corr.diff$data)[2]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat,NHseason=NULL,
                plotname=paste0('corr.diff_echam_anal-',v,'_slp.pdf'), paper='special')
  }
  
  ###############################################################
  
  
  
  ###############################################################
  ############# figs.11: Validation Bias Maps  ##################
  ###############################################################
  
  
  ##################### Temperature Bias ########################
  
  if (!monthly_out){
    bias.tot <- echam
    bias.tot$data <- array(cbind(bias.ech$ensmean,bias$ensmean),c(dim(bias$ensmean)[1],1,dim(bias$ensmean)[2]*2))
    bias.tot$ensmean <- array(cbind(bias.ech$ensmean,bias$ensmean),c(dim(bias$ensmean)[1],dim(bias$ensmean)[2]*2))
    
    levs <- c(-Inf,-10,-5,-2,-1,0,1,2,5,10,Inf)
    
    plot_echam4(bias.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(bias.tot$data)[3]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat, NHseason=NULL,
                plotname=paste0('bias_echam_anal-',v,'_temp.pdf'), paper='special')
    
    
    ##################### Precipitation Bias ########################
    
    
    levs <- c(-Inf,-40,-20,-10,-5,0,5,10,20,40,Inf)
    
    plot_echam4(bias.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(bias.tot$data)[3]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat, NHseason=NULL,
                plotname=paste0('bias_echam_anal-',v,'_precip.pdf'), paper='special')
    
    
    ##################### Sea Level Pressure Bias ########################
    
    
    levs <- c(-Inf,-10,-5,-2,-1,0,1,2,5,10,Inf)
    
    plot_echam4(bias.tot, varname='slp', cex.pt=1.5, names=pnames[1:dim(bias.tot$data)[3]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat, NHseason=NULL,
                plotname=paste0('bias_echam_anal-',v,'_slp.pdf'), paper='special')
    
  }
  
  ###############################################################
  
  
  
  ###############################################################
  ############ figs.12: Reduction of Error Maps  ################
  ###############################################################
  
  
  ################### Temperature RE Anomaly ###################
  
  RE.tot <- echam
  RE.tot$data <- array(RE.anom$ensmean,c(dim(RE$ensmean)[1], 1, dim(RE$ensmean)[2]))
  RE.tot$ensmean <- array(RE.anom$ensmean,c(dim(RE$ensmean)[1], dim(RE$ensmean)[2]))
  
  levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0.05,0.2,0.4,0.6,0.8,1)
  
  if (monthly_out){
    plot_echam4(RE.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat, NHseason="summer",
                plotname=paste0('re_anom_echam_anal-',v,'_temp_summer.pdf'), paper='special')
    plot_echam4(RE.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat, NHseason="winter",
                plotname=paste0('re_anom_echam_anal-',v,'_temp_winter.pdf'), paper='special')
  }else{
    plot_echam4(RE.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat, NHseason=NULL,
                plotname=paste0('re_anom_echam_anal-',v,'_temp.pdf'), paper='special')
  }
  
  #################### Precipitation RE Anomaly ######################
  
  if (monthly_out){
    plot_echam4(RE.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat, NHseason="summer",
                plotname=paste0('re_anom_echam_anal-',v,'_precip_summer.pdf'), paper='special')
    plot_echam4(RE.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat, NHseason="winter",
                plotname=paste0('re_anom_echam_anal-',v,'_precip_winter.pdf'), paper='special')
  }else{
    plot_echam4(RE.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat, NHseason=NULL,
                plotname=paste0('re_anom_echam_anal-',v,'_precip.pdf'), paper='special')
  }
  
  ##################### Sea Level Pressure RE Anomaly ########################
  
  levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0,0.2,0.4,0.6,0.8,1)
  
  if (monthly_out){
    plot_echam4(RE.tot, varname='slp', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat, NHseason="summer",
                plotname=paste0('re_anom_echam_anal-',v,'_slp_summer.pdf'), paper='special')
    plot_echam4(RE.tot, varname='slp', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat, NHseason="winter",
                plotname='re_anom_echam_anal-',v,'_slp_winter.pdf', paper='special')
  }else{
    plot_echam4(RE.tot, varname='slp', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat, NHseason=NULL,
                plotname=paste0('re_anom_echam_anal-',v,'_slp.pdf'), paper='special')
  }
  
  
  ##################### Temperature RE Absolute ########################
  
  RE.tot <- echam
  RE.tot$data <- array(RE$ensmean,c(dim(RE$ensmean)[1], 1, dim(RE$ensmean)[2]))
  RE.tot$ensmean <- array(RE$ensmean,c(dim(RE$ensmean)[1], dim(RE$ensmean)[2]))
  
  
  levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0.05,0.2,0.4,0.6,0.8,1)
  if (monthly_out){
    plot_echam4(RE.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat, NHseason="summer",
                plotname=paste0('re_echam_anal-',v,'_temp_summer.pdf'), paper='special')
    plot_echam4(RE.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat, NHseason="winter",
                plotname=paste0('re_echam_anal-',v,'_temp_winter.pdf'), paper='special')
  }else{
    plot_echam4(RE.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat, NHseason=NULL,
                plotname=paste0('re_echam_anal-',v,'_temp.pdf'), paper='special')
  }
  
  ##################### Precipitation RE Absolute ########################
  
  if (monthly_out){
    plot_echam4(RE.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat, NHseason="summer",
                plotname=paste0('re_echam_anal-',v,'_precip_summer.pdf'), paper='special')
    plot_echam4(RE.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat, NHseason="winter",
                plotname=paste0('re_echam_anal-',v,'_precip_winter.pdf'), paper='special')
  }else{
    plot_echam4(RE.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat, NHseason=NULL,
                plotname=paste0('re_echam_anal-',v,'_precip.pdf'), paper='special')
  }
  
  ##################### Sea Level Pressure RE Absolute ########################
  
  levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0,0.2,0.4,0.6,0.8,1)
  
  if (monthly_out){
    plot_echam4(RE.tot, varname='slp', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat, NHseason="summer",
                plotname=paste0('re_echam_anal-',v,'_slp_summer.pdf'), paper='special')
    plot_echam4(RE.tot, varname='slp', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs,
                type='ensmean', st.col=NULL, stations=plstat, NHseason="winter",
                plotname=paste0('re_echam_anal-',v,'_slp_winter.pdf'), paper='special')
  }else{
    plot_echam4(RE.tot, varname='slp', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat, NHseason=NULL,
                plotname=paste0('re_echam_anal-',v,'_slp.pdf'), paper='special')
  }
  
  
  ##################### Temperature RE Climatology Absolute ########################
  
  RE.tot <- echam
  RE.tot$data <- array(RE.echam.clim$ensmean,c(dim(RE$ensmean)[1], 1, dim(RE$ensmean)[2]))
  RE.tot$ensmean <- array(RE.echam.clim$ensmean,c(dim(RE$ensmean)[1],dim(RE$ensmean)[2]))
  
  levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0.05,0.2,0.4,0.6,0.8,1)
  
  if (monthly_out){
    plot_echam4(RE.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat, NHseason="summer",
                plotname=paste0('re_clim_echam_anal-',v,'_temp_summer.pdf'), paper='special')
    plot_echam4(RE.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat, NHseason="winter",
                plotname=paste0('re_clim_echam_anal-',v,'_temp_winter.pdf'), paper='special')
  }else{
    plot_echam4(RE.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat, NHseason=NULL,
                plotname=paste0('re_clim_echam_anal-',v,'_temp.pdf'), paper='special')
  }
  
  ##################### Precipitation RE Climatology Absolute ########################
  
  if (monthly_out){
    plot_echam4(RE.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat, NHseason="summer",
                plotname=paste0('re_clim_echam_anal-',v,'_precip_summer.pdf'), paper='special')
    plot_echam4(RE.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat, NHseason="winter",
                plotname=paste0('re_clim_echam_anal-',v,'_precip_winter.pdf'), paper='special')
  }else{
    plot_echam4(RE.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat, NHseason=NULL,
                plotname=paste0('re_clim_echam_anal-',v,'_precip.pdf'), paper='special')
  }
  
  ##################### Sea Level Pressure RE Climatology Absolute ########################
  
  levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0,0.2,0.4,0.6,0.8,1)
  
  if (monthly_out){
    plot_echam4(RE.tot, varname='slp', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat, NHseason="summer",
                plotname=paste0('re_clim_echam_anal-',v,'_slp_summer.pdf'), paper='special')
    plot_echam4(RE.tot, varname='slp', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat, NHseason="winter",
                plotname=paste0('re_clim_echam_anal-',v,'_slp_winter.pdf'), paper='special')
  }else{
    plot_echam4(RE.tot, varname='slp', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat, NHseason=NULL,
                plotname=paste0('re_clim_echam_anal-',v,'_slp.pdf'), paper='special')
  }
  
  ##################### Temperature RE Climatology Anomaly ########################
  
  RE.tot <- echam
  RE.tot$data <- array(RE.echam.clim.anom$ensmean,c(dim(RE$ensmean)[1], 1, dim(RE$ensmean)[2]))
  RE.tot$ensmean <- array(RE.echam.clim.anom$ensmean,c(dim(RE$ensmean)[1], dim(RE$ensmean)[2]))
  
  levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0.05,0.2,0.4,0.6,0.8,1)
  
  if (monthly_out){
    plot_echam4(RE.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat, NHseason="summer",
                plotname=paste0('re_clim_anom_echam_anal-',v,'_temp_summer.pdf'), paper='special')
    plot_echam4(RE.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat, NHseason="winter",
                plotname=paste0('re_clim_anom_echam_anal-',v,'_temp_winter.pdf'), paper='special')
  }else{
    plot_echam4(RE.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat, NHseason=NULL,
                plotname=paste0('re_clim_anom_echam_anal-',v,'_temp.pdf'), paper='special')
  }
  
  ##################### Precipitation RE Climatology Anomaly ########################
  
  if (monthly_out){
    plot_echam4(RE.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat, NHseason="summer",
                plotname=paste0('re_clim_anom_echam_anal-',v,'_precip_summer.pdf'), paper='special')
    plot_echam4(RE.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat, NHseason="winter",
                plotname=paste0('re_clim_anom_echam_anal-',v,'_precip_winter.pdf'), paper='special')
  }else{
    plot_echam4(RE.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat, NHseason=NULL,
                plotname=paste0('re_clim_anom_echam_anal-',v,'_precip.pdf'), paper='special')
  }
  
  ##################### Sea Level Pressure RE Climatology Anomaly ########################
  
  levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0,0.2,0.4,0.6,0.8,1)
  
  if (monthly_out){
    plot_echam4(RE.tot, varname='slp', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat, NHseason="summer",
                plotname=paste0('re_clim_anom_echam_anal-',v,'_slp_summer.pdf'), paper='special')
    plot_echam4(RE.tot, varname='slp', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat, NHseason="winter",
                plotname=paste0('re_clim_anom_echam_anal-',v,'_slp_winter.pdf'), paper='special')
  }else{
    plot_echam4(RE.tot, varname='slp', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, 
                type='ensmean', st.col=NULL, stations=plstat, NHseason=NULL,
                plotname=paste0('re_clim_anom_echam_anal-',v,'_slp.pdf'), paper='special')
  }
  
  
  ###############################################################
  
  
  
  ###############################################################
  ############ figs.13: RootMeanSquareError Maps  ###############
  ###############################################################
  
  # Diff: Analysis-Echam-> the smaller the value the better
  # here there is a for loop for all variables 
  rm(levs)
  rmse.diff<-rmse
  rmse.diff$data<-rmse$data-rmse.ech$data
  rmse.diff$ensmean<-rmse$ensmean-rmse.ech$ensmean
  rmse.anom.diff<-rmse.anom
  rmse.anom.diff$data<-rmse.anom$data-rmse.ech.anom$data
  rmse.anom.diff$ensmean<-rmse.anom$ensmean-rmse.ech.anom$ensmean
  for (vari in unique(rmse$names)){
    
    ###################################################################################
    # Fig.: rmse analysis
    levs<-c(pretty(rmse$ensmean[which(rmse$names==vari)],n=11),Inf)
    if (monthly_out){
      plot_echam4(rmse, varname=vari, cex.pt=1.5, names=pnames[1:dim(rmse$data)[2]],  type='ensmean', 
                  st.col=NULL, stations=plstat, NHseason="summer",
                  plotname=paste0('rmse_anal-',v,'_',vari,'_summer.pdf'), paper='special')
      plot_echam4(rmse, varname=vari, cex.pt=1.5, names=pnames[1:dim(rmse$data)[2]],  type='ensmean', 
                  st.col=NULL, stations=plstat, NHseason="winter",
                  plotname=paste0('rmse_anal-',v,'_',vari,'_winter.pdf'), paper='special')
    }else{
      plot_echam4(rmse, varname=vari, cex.pt=1.5, names=pnames[1:dim(rmse$data)[2]],lev=levs,
                  cols=c('lightblue','royalblue','navy'), type='ensmean', st.col=NULL, 
                  stations=plstat,NHseason=NULL,plotname=paste0('rmse_anal-',v,'_',vari,'.pdf'), paper='special')
      
    }
    
    ###################################################################################
    # Fig.: rmse analysis anom
    levs<-c(pretty(rmse.anom$ensmean[which(rmse.anom$names==vari)],n=11),Inf)
    if (monthly_out){
      plot_echam4(rmse.anom, varname=vari, cex.pt=1.5, names=pnames[1:dim(rmse.anom$data)[2]], 
                  symmetric=F, type='ensmean', st.col=NULL, stations=plstat, NHseason="summer",
                  plotname=paste0('rmse_anom_anal-',v,'_',vari,'_summer.pdf'), paper='special')
      plot_echam4(rmse.anom, varname=vari, cex.pt=1.5, names=pnames[1:dim(rmse.anom$data)[2]], 
                  symmetric=F, type='ensmean', st.col=NULL, stations=plstat,NHseason="winter",
                  plotname=paste0('rmse_anom_anal-',v,'_',vari,'_winter.pdf'), paper='special')
    }else{
      plot_echam4(rmse.anom, varname=vari, cex.pt=1.5, names=pnames[1:dim(rmse.anom$data)[2]],
                  lev=levs,cols=c('lightblue','royalblue','navy'), type='ensmean', st.col=NULL, 
                  stations=plstat,NHseason=NULL,plotname=paste0('rmse_anom_anal-',v,'_',vari,'.pdf'), paper='special')
    }
    
    ###################################################################################
    # Fig.: rmse echam
    levs<-pretty(rmse.ech$ensmean[which(rmse.ech$names==vari)],n=11)
    if (monthly_out){
      plot_echam4(rmse.ech, varname=vari, cex.pt=1.5, names=pnames[1:dim(rmse.ech$data)[2]], 
                  symmetric=F, type='ensmean', st.col=NULL, stations=plstat,NHseason="summer",
                  plotname=paste0('rmse_echam-',v,'_',vari,'_summer.pdf'), paper='special')
      plot_echam4(rmse.ech, varname=vari, cex.pt=1.5, names=pnames[1:dim(rmse.ech$data)[2]], 
                  symmetric=F, type='ensmean', st.col=NULL, stations=plstat,NHseason="winter",
                  plotname=paste0('rmse_echam-',v,'_',vari,'_winter.pdf'), paper='special')
    }else{
      plot_echam4(rmse.ech, varname=vari, cex.pt=1.5, names=pnames[1:dim(rmse.ech$data)[2]],
                  lev=levs,cols=c('lightblue','royalblue','navy'), type='ensmean', st.col=NULL, 
                  stations=plstat,NHseason=NULL,plotname=paste0('rmse_echam-',v,'_',vari,'.pdf'), paper='special')
    }
    
    ###################################################################################
    # Fig.: rmse echam anom
    levs<-c(pretty(rmse.ech.anom$ensmean[which(rmse.ech.anom$names==vari)],n=11),Inf)
    if (monthly_out){
      plot_echam4(rmse.ech.anom, varname=vari, cex.pt=1.5, names=pnames[1:dim(rmse.ech.anom$data)[2]], 
                  symmetric=F, type='ensmean', st.col=NULL, stations=plstat, NHseason="summer",
                  plotname=paste0('rmse_anom_echam-',v,'_',vari,'_summer.pdf'), paper='special')
      plot_echam4(rmse.ech.anom, varname=vari, cex.pt=1.5, names=pnames[1:dim(rmse.ech.anom$data)[2]], 
                  symmetric=F, type='ensmean', st.col=NULL, stations=plstat, NHseason="winter",
                  plotname=paste0('rmse_anom_echam-',v,'_',vari,'_winter.pdf'), paper='special')
    }else{
      plot_echam4(rmse.ech.anom, varname=vari, cex.pt=1.5, names=pnames[1:dim(rmse.ech.anom$data)[2]],
                  lev=levs,cols=c('lightblue','royalblue','navy'), type='ensmean', st.col=NULL, 
                  stations=plstat,NHseason=NULL,plotname=paste0('rmse_anom_echam-',v,'_',vari,'.pdf'), paper='special')
    }
    
    ###################################################################################
    #Fig.: rmse diff
    vmax<-max(rmse.diff$ensmean[which(rmse.diff$names==vari),],na.rm=T)
    vmin<-min(rmse.diff$ensmean[which(rmse.diff$names==vari),],na.rm=T)
    vmaxmax<-max(vmax,-1*vmin)
    levs<-c(-Inf,signif(seq(-vmaxmax,vmaxmax,by=2*vmaxmax/9),2),Inf)
    if (monthly_out){
      plot_echam4(rmse.diff, varname=vari, cex.pt=1.5, names=pnames[1:dim(rmse.diff$data)[2]], 
                  symmetric=F, type='ensmean', st.col=NULL, stations=plstat, NHseason="summer",
                  plotname=paste0('rmse_diff-',v,'_',vari,'_summer.pdf'), paper='special')
      plot_echam4(rmse.diff, varname=vari, cex.pt=1.5, names=pnames[1:dim(rmse.diff$data)[2]], 
                  symmetric=F, type='ensmean', st.col=NULL, stations=plstat, NHseason="winter",
                  plotname=paste0('rmse_diff-',v,'_',vari,'_winter.pdf'), paper='special')
    }else{
      plot_echam4(rmse.diff, varname=vari, cex.pt=1.5, names=pnames[1:dim(rmse.diff$data)[2]],
                  lev=levs,cols=c('yellow','white','navy'), type='ensmean', st.col=NULL, 
                  stations=plstat,NHseason=NULL,plotname=paste0('rmse_diff-',v,'_',vari,'.pdf'), paper='special')
    }
    
    ###################################################################################
    #Fig.: rmse.anom diff
    vmax<-max(rmse.anom.diff$ensmean[which(rmse.anom.diff$names==vari),],na.rm=T)
    vmin<-min(rmse.anom.diff$ensmean[which(rmse.anom.diff$names==vari),],na.rm=T)
    vmaxmax<-max(vmax,-1*vmin)
    levs<-c(-Inf,signif(seq(-vmaxmax,vmaxmax,by=2*vmaxmax/9),2),Inf)
    if (monthly_out){
      plot_echam4(rmse.anom.diff, varname=vari, cex.pt=1.5, names=pnames[1:dim(rmse.anom.diff$data)[2]], 
                  symmetric=F, type='ensmean', st.col=NULL, stations=plstat, NHseason="summer",
                  plotname=paste0('rmse_diff_anom-',v,'_',vari,'_summer.pdf'), paper='special')
      plot_echam4(rmse.anom.diff, varname=vari, cex.pt=1.5, names=pnames[1:dim(rmse.anom.diff$data)[2]],
                  symmetric=F,  type='ensmean', st.col=NULL, stations=plstat, NHseason="winter",
                  plotname=paste0('rmse_diff_anom-',v,'_',vari,'_winter.pdf'), paper='special')
    }else{
      plot_echam4(rmse.anom.diff, varname=vari, cex.pt=1.5, names=pnames[1:dim(rmse.anom.diff$data)[2]],
                  lev=levs,cols=c('yellow','white','navy'), type='ensmean', st.col=NULL, stations=plstat,
                  NHseason=NULL,plotname=paste0('rmse_diff_anom-',v,'_',vari,'.pdf'), paper='special')
    }
  } #end rmse plots for loop
  
  
  ###############################################################
  
  
  
  ###############################################################
  ############### figs.14: Validation Indices TS  ###############
  ###############################################################
  
  #if (ind_ECHAM&!monthly_out) {
  if (!monthly_out) {
    ################# Absolute values ##########################
    pdf(paste(figpath,'/vali_ind.pdf',sep='/'), width=9, height=6, paper='special')
      par(mfrow=c(4,2), cex.axis=1.4, cex.lab=1.4, mar=c(3,5,1,1), oma=c(0,0,0,0))
      vind2=vind
      vind2$data=array(vind$data,c(dim(vind$data)[1],2,dim(vind$data)[2]/2))
      eind2=eind
      eind2$ensmean=array(eind$ensmean,c(dim(eind$ensmean)[1],2,dim(eind$ensmean)[2]/2))
      eind2$min=array(apply(eind$data,1:2,min),c(dim(eind$ensmean)[1],2,dim(eind$ensmean)[2]/2))
      eind2$max=array(apply(eind$data,1:2,max),c(dim(eind$ensmean)[1],2,dim(eind$ensmean)[2]/2))
      aind2=aind
      aind2$ensmean=array(aind$ensmean,c(dim(aind$ensmean)[1],2,dim(aind$ensmean)[2]/2))
      aind2$min=array(apply(aind$data,1:2,min),c(dim(aind$ensmean)[1],2,dim(aind$ensmean)[2]/2))
      aind2$max=array(apply(aind$data,1:2,max),c(dim(aind$ensmean)[1],2,dim(aind$ensmean)[2]/2))
      if(!tps_only){
        indices <- c('ENH.temp2','NAM.temp2','SAM.temp2','AFR.temp2',
                     'ASI.temp2','AUS.temp2','ARC.temp2','ANT.temp2',
                     'NEU.temp2','MED.temp2',
                     'GLO.temp2','NAM.precip','SAM.precip','AFR.precip',
                     'ASI.precip','AUS.precip','ARC.precip','ANT.precip',
                     'NEU.precip','MED.precip',
                     'SH.temp2','NAM.slp','SAM.slp','AFR.slp',
                     'ASI.slp','AUS.slp','ARC.slp','ANT.slp',
                     'NEU.slp','MED.slp',
                     'NH.temp2', 'EU.temp2', 'EU.precip', 'EU.slp',
                     'HC.calc', 'SJ_u200.calc', 'SJ_slp.calc', 'PV.calc',
                     'PWC.calc', 'DIMI.calc') #,'ITCZ.calc', 'NAO.calc', 'PNA.calc')
        vind2$names[which(vind2$names=="ind_recon_dimi")]<-"DIMI.calc"
        vind2$names[which(vind2$names=="ind_recon_z100")]<-"PV.calc"
        #vind2$names[which(vind2$names=="ind_recon_z300")]<-midlattitude circulation
        vind2$names[which(vind2$names=="ind_recon_pwc")]<-"PWC.calc"
        vind2$names[which(vind2$names=="ind_recon_hc")]<-"HC.calc"
        vind2$names[which(vind2$names=="ind_recon_sj")]<-"SJ_slp.calc"
        units<-c(rep('[ºC]',11),rep('[mm]',9),'[ºC]',rep('[hPa]',9),rep('[ºC]',2),'[mm]','[hPa]',rep('',6))
      } else {
        indices <- c('ENH.temp2','NAM.temp2','SAM.temp2','AFR.temp2',
                     'ASI.temp2','AUS.temp2','ARC.temp2','ANT.temp2',
                     'NEU.temp2','MED.temp2',
                     'GLO.temp2','NAM.precip','SAM.precip','AFR.precip',
                     'ASI.precip','AUS.precip','ARC.precip','ANT.precip',
                     'NEU.precip','MED.precip',
                     'SH.temp2','NAM.slp','SAM.slp','AFR.slp',
                     'ASI.slp','AUS.slp','ARC.slp','ANT.slp',
                     'NEU.slp','MED.slp',
                     'NH.temp2', 'EU.temp2', 'EU.precip', 'EU.slp')
        units<-c(rep('[ºC]',11),rep('[mm]',9),'[ºC]',rep('[hPa]',9),rep('[ºC]',2),'[mm]','[hPa]')
      }
      ti=seq(1,length(vind$time),by=2)
      i=1
      for (ind in indices) {
        mainname=ind
        for (seas in c(1,2)) {
          if(all(is.nan(eind2$ensmean[which(eind2$names==ind),,]))){
            next()
          }else{
            if (ind=="SJ_u200.calc"){ # only one SJ in vind but 2 different calculations ind aind and eind
              if (seas == 1) {color='darkblue'; color3='blue'; color2='cyan'}
              if (seas == 2) {color='darkred'; color3='red'; color2='orange'}
              if (units[i]=='[ºC]'){
                ymin=min(eind2$min[which(eind2$names==ind),seas,])-5
                ymax=max(eind2$max[which(eind2$names==ind),seas,])+5
              } else if (units[i]=='[mm]') {
                ymin=min(eind2$min[which(eind2$names==ind),seas,])-20
                ymax=max(eind2$max[which(eind2$names==ind),seas,])+20
              } else if (units[i]=='[hPa]') {
                ymin=min(eind2$min[which(eind2$names==ind),seas,])-2
                ymax=max(eind2$max[which(eind2$names==ind),seas,])+2
              } else {
                ymin=min(eind2$min[which(eind2$names==ind),seas,])
                ymax=max(eind2$max[which(eind2$names==ind),seas,])
              }
              plot(vind2$time[ti],vind2$data[which(vind2$names=='SJ_slp.calc'),seas,],ty='l',
                   col=color,lty=1,ylim=c(ymin,ymax),main=mainname,ylab=units[i],xlab='')
              lines(eind2$time[ti],eind2$ensmean[which(eind2$names==ind),seas,],ty='l',col=color2,lwd=2,lty=2,main='')
              lines(eind2$time[ti],eind2$min[which(eind2$names==ind),seas,],ty='l',col=color2,lty=2,lwd=1,main='')
              lines(eind2$time[ti],eind2$max[which(eind2$names==ind),seas,],ty='l',col=color2,lty=2,lwd=1,main='')
              lines(aind2$time[ti],aind2$ensmean[which(aind2$names==ind),seas,],ty='l',col=color3,lwd=2,lty=3,main='')
              lines(aind2$time[ti],aind2$min[which(aind2$names==ind),seas,],ty='l',col=color3,lty=3,lwd=1,main='')
              lines(aind2$time[ti],aind2$max[which(aind2$names==ind),seas,],ty='l',col=color3,lty=3,lwd=1,main='')
            } else {
              if (seas == 1) {color='darkblue'; color3='blue'; color2='cyan'}
              if (seas == 2) {color='darkred'; color3='red'; color2='orange'}
              if (units[i]=='[ºC]'){
                ymin=min(eind2$min[which(eind2$names==ind),seas,])-5
                ymax=max(eind2$max[which(eind2$names==ind),seas,])+5
              } else if (units[i]=='[mm]') {
                ymin=min(eind2$min[which(eind2$names==ind),seas,])-20
                ymax=max(eind2$max[which(eind2$names==ind),seas,])+20
              } else if (units[i]=='[hPa]') {
                ymin=min(eind2$min[which(eind2$names==ind),seas,])-2
                ymax=max(eind2$max[which(eind2$names==ind),seas,])+2
              } else {
                ymin=min(eind2$min[which(eind2$names==ind),seas,],na.rm=T)-2
                ymax=max(eind2$max[which(eind2$names==ind),seas,],na.rm=T)+2
              }
              plot(vind2$time[ti],vind2$data[which(vind2$names==ind),seas,],ty='l',col=color,
                   lty=1,ylim=c(ymin,ymax),main=mainname,ylab=units[i],xlab='')
              lines(eind2$time[ti],eind2$ensmean[which(eind2$names==ind),seas,],ty='l',col=color2,lwd=2,lty=2,main='')
              lines(eind2$time[ti],eind2$min[which(eind2$names==ind),seas,],ty='l',col=color2,lty=2,lwd=1,main='')
              lines(eind2$time[ti],eind2$max[which(eind2$names==ind),seas,],ty='l',col=color2,lty=2,lwd=1,main='')
              lines(aind2$time[ti],aind2$ensmean[which(aind2$names==ind),seas,],ty='l',col=color3,lwd=2,lty=3,main='')
              lines(aind2$time[ti],aind2$min[which(aind2$names==ind),seas,],ty='l',col=color3,lty=3,lwd=1,main='')
              lines(aind2$time[ti],aind2$max[which(aind2$names==ind),seas,],ty='l',col=color3,lty=3,lwd=1,main='')
            }
          }
        }
        i=i+1
      }
    dev.off()
    
    
    #################### Anomaly Indices ##########################
    pdf(paste(figpath,'/vali_ind_anom.pdf',sep='/'), width=9, height=6, paper='special')
      par(mfrow=c(4,2), cex.axis=1.4, cex.lab=1.4, mar=c(3,5,1,1), oma=c(0,0,0,0))
      vind2=vind.anom
      vind2$data=array(vind.anom$data,c(dim(vind.anom$data)[1],2,dim(vind.anom$data)[2]/2))
      eind2=eind.anom
      eind2$ensmean=array(eind.anom$ensmean,c(dim(eind.anom$ensmean)[1],2,dim(eind.anom$ensmean)[2]/2))
      eind2$min=array(apply(eind.anom$data,1:2,min),c(dim(eind.anom$ensmean)[1],2,dim(eind.anom$ensmean)[2]/2))
      eind2$max=array(apply(eind.anom$data,1:2,max),c(dim(eind.anom$ensmean)[1],2,dim(eind.anom$ensmean)[2]/2))
      aind2=aind.anom
      aind2$ensmean=array(aind.anom$ensmean,c(dim(aind.anom$ensmean)[1],2,dim(aind.anom$ensmean)[2]/2))
      aind2$min=array(apply(aind.anom$data,1:2,min),c(dim(aind.anom$ensmean)[1],2,dim(aind.anom$ensmean)[2]/2))
      aind2$max=array(apply(aind.anom$data,1:2,max),c(dim(aind.anom$ensmean)[1],2,dim(aind.anom$ensmean)[2]/2))
      if(!tps_only){
        indices <- c('ENH.temp2','NAM.temp2','SAM.temp2','AFR.temp2',
                     'ASI.temp2','AUS.temp2','ARC.temp2','ANT.temp2',
                     'NEU.temp2','MED.temp2',
                     'GLO.temp2','NAM.precip','SAM.precip','AFR.precip',
                     'ASI.precip','AUS.precip','ARC.precip','ANT.precip',
                     'NEU.precip','MED.precip',
                     'SH.temp2','NAM.slp','SAM.slp','AFR.slp',
                     'ASI.slp','AUS.slp','ARC.slp','ANT.slp',
                     'NEU.slp','MED.slp',
                     'NH.temp2', 'EU.temp2', 'EU.precip', 'EU.slp',
                     'HC.calc', 'SJ_u200.calc', 'SJ_slp.calc', 'PV.calc',
                     'PWC.calc', 'DIMI.calc') #,'ITCZ.calc', 'NAO.calc', 'PNA.calc')
        vind2$names[which(vind2$names=="ind_recon_dimi")]<-"DIMI.calc"
        vind2$names[which(vind2$names=="ind_recon_z100")]<-"PV.calc"
        #vind2$names[which(vind2$names=="ind_recon_z300")]<-midlattitude circulation
        vind2$names[which(vind2$names=="ind_recon_pwc")]<-"PWC.calc"
        vind2$names[which(vind2$names=="ind_recon_hc")]<-"HC.calc"
        vind2$names[which(vind2$names=="ind_recon_sj")]<-"SJ_slp.calc"
        units<-c(rep('[ºC]',11),rep('[mm]',9),'[ºC]',rep('[hPa]',9),rep('[ºC]',2),'[mm]','[hPa]',rep('',6))
      } else {
        indices <- c('ENH.temp2','NAM.temp2','SAM.temp2','AFR.temp2',
                     'ASI.temp2','AUS.temp2','ARC.temp2','ANT.temp2',
                     'NEU.temp2','MED.temp2',
                     'GLO.temp2','NAM.precip','SAM.precip','AFR.precip',
                     'ASI.precip','AUS.precip','ARC.precip','ANT.precip',
                     'NEU.precip','MED.precip',
                     'SH.temp2','NAM.slp','SAM.slp','AFR.slp',
                     'ASI.slp','AUS.slp','ARC.slp','ANT.slp',
                     'NEU.slp','MED.slp',
                     'NH.temp2', 'EU.temp2', 'EU.precip', 'EU.slp')
        units<-c(rep('[ºC]',11),rep('[mm]',9),'[ºC]',rep('[hPa]',9),rep('[ºC]',2),'[mm]','[hPa]')
      }
      ti=seq(1,length(vind$time),by=2)
      i=1
      for (ind in indices) {
        mainname=ind
        for (seas in c(1,2)) {
          if (all(is.nan(eind2$ensmean[which(eind2$names==ind),,]))) {
            next()
          } else {
            if (ind=="SJ_u200.calc"){ # only one SJ in vind but 2 different calculations ind aind and eind
              if (seas == 1) {color='darkblue'; color3='blue'; color2='cyan'}
              if (seas == 2) {color='darkred'; color3='red'; color2='orange'}
              if (units[i]=='[ºC]') {
                ymin=min(eind2$min[which(eind2$names==ind),seas,])
                ymax=max(eind2$max[which(eind2$names==ind),seas,])
              } else if (units[i]=='[mm]') {
                ymin=min(eind2$min[which(eind2$names==ind),seas,])
                ymax=max(eind2$max[which(eind2$names==ind),seas,])
              } else if (units[i]=='[hPa]') {
                ymin=min(eind2$min[which(eind2$names==ind),seas,])
                ymax=max(eind2$max[which(eind2$names==ind),seas,])
              } else {
                ymin=min(eind2$min[which(eind2$names==ind),seas,])
                ymax=max(eind2$max[which(eind2$names==ind),seas,])
              }
              plot(vind2$time[ti],vind2$data[which(vind2$names=='SJ_slp.calc'),seas,],ty='l',
                   col=color,lty=1,ylim=c(ymin,ymax),main=mainname,ylab=units[i],xlab='')
              lines(eind2$time[ti],eind2$ensmean[which(eind2$names==ind),seas,],ty='l',col=color2,lwd=2,lty=2,main='')
              lines(eind2$time[ti],eind2$min[which(eind2$names==ind),seas,],ty='l',col=color2,lty=2,lwd=1,main='')
              lines(eind2$time[ti],eind2$max[which(eind2$names==ind),seas,],ty='l',col=color2,lty=2,lwd=1,main='')
              lines(aind2$time[ti],aind2$ensmean[which(aind2$names==ind),seas,],ty='l',col=color3,lwd=2,lty=3,main='')
              lines(aind2$time[ti],aind2$min[which(aind2$names==ind),seas,],ty='l',col=color3,lty=3,lwd=1,main='')
              lines(aind2$time[ti],aind2$max[which(aind2$names==ind),seas,],ty='l',col=color3,lty=3,lwd=1,main='')
            } else {
              if (seas == 1) {color='darkblue'; color3='blue'; color2='cyan'}
              if (seas == 2) {color='darkred'; color3='red'; color2='orange'}
              if (units[i]=='[ºC]') {
                ymin=min(eind2$min[which(eind2$names==ind),seas,])-5
                ymax=max(eind2$max[which(eind2$names==ind),seas,])+5
              } else if (units[i]=='[mm]') {
                ymin=min(eind2$min[which(eind2$names==ind),seas,])-20
                ymax=max(eind2$max[which(eind2$names==ind),seas,])+20
              } else if (units[i]=='[hPa]') {
                ymin=min(eind2$min[which(eind2$names==ind),seas,])-2
                ymax=max(eind2$max[which(eind2$names==ind),seas,])+2
              } else {
                ymin=min(eind2$min[which(eind2$names==ind),seas,],na.rm=T)-2
                ymax=max(eind2$max[which(eind2$names==ind),seas,],na.rm=T)+2
              }
              plot(vind2$time[ti],vind2$data[which(vind2$names==ind),seas,],ty='l',col=color,lty=1,
                   ylim=c(ymin,ymax),main=mainname,ylab=units[i],xlab='')
              lines(eind2$time[ti],eind2$ensmean[which(eind2$names==ind),seas,],ty='l',col=color2,lwd=2,lty=2,main='')
              lines(eind2$time[ti],eind2$min[which(eind2$names==ind),seas,],ty='l',col=color2,lty=2,lwd=1,main='')
              lines(eind2$time[ti],eind2$max[which(eind2$names==ind),seas,],ty='l',col=color2,lty=2,lwd=1,main='')
              lines(aind2$time[ti],aind2$ensmean[which(aind2$names==ind),seas,],ty='l',col=color3,lwd=2,lty=3,main='')
              lines(aind2$time[ti],aind2$min[which(aind2$names==ind),seas,],ty='l',col=color3,lty=3,lwd=1,main='')
              lines(aind2$time[ti],aind2$max[which(aind2$names==ind),seas,],ty='l',col=color3,lty=3,lwd=1,main='')
            }
          }
        }
        i=i+1
      }
    dev.off()
    ###############################################################
    
    
    
    
    
    ###############################################################
    ######### figs.15: Validation Indices Correlations  ###########
    ###############################################################
    
    if (!pseudo_prox) {
    inds <- c('ENH.temp2', 'EU.temp2', 'EU.precip', 'EU.slp')
    indices <- c('ENH.temp2','NAM.temp2','SAM.temp2','AFR.temp2',
                 'ASI.temp2','AUS.temp2','ARC.temp2','ANT.temp2',
                 'NEU.temp2','MED.temp2',
                 'GLO.temp2','NAM.precip','SAM.precip','AFR.precip',
                 'ASI.precip','AUS.precip','ARC.precip','ANT.precip',
                 'NEU.precip','MED.precip',
                 'SH.temp2','NAM.slp','SAM.slp','AFR.slp',
                 'ASI.slp','AUS.slp','ARC.slp','ANT.slp',
                 'NEU.slp','MED.slp',
                 'NH.temp2', 'EU.temp2', 'EU.precip', 'EU.slp')
    RE.mat <- array(NA, c(length(inds), 2, 31))
    corr.mat <- array(NA, c(length(inds)*2, 2, 31))
    ind.names <- list()
    for (se in 1:2){
      ind.names[[se]] <- inds
      for (k in seq(inds)){
        RE.mat[k,se,1:30] <- RE.ind$data[RE.ind$names == inds[k],se,]
        RE.mat[k,se,31] <- RE.ind$ensmean[RE.ind$names == inds[k],se]
        corr.mat[k*2,se,1:30] <- acorr.ind$data[acorr.ind$names == inds[k],se,]
        corr.mat[k*2,se,31] <- acorr.ind$ensmean[acorr.ind$names == inds[k],se]
        corr.mat[k*2-1,se,1:30] <- ecorr.ind$data[ecorr.ind$names == inds[k],se,]
        corr.mat[k*2-1,se,31] <- ecorr.ind$ensmean[ecorr.ind$names == inds[k],se]
      }
    }
    
    pdf(paste(figpath,'/index_correlation.pdf',sep='/'), width=9, height=6, paper='special')
      par(mfrow=c(2,2), mar=c(1,1,0,0), oma=c(0,2,1,0), cex.axis=1.4, cex.lab=1.4)
      ind.col <- list()
      ind.col2 <- list()
      ind.col[[1]] <- hcl(c(0,0,120,240,270,300), c=40, l=50)
      ind.col[[2]] <- hcl(c(0,0,120,240,270,200), c=40, l=50)
      ind.col2[[1]] <- hcl(rep(c(0,0,120,240,270,300), each=2), c=40, l=c(90,50))
      ind.col2[[2]] <- hcl(rep(c(0,0,120,240,270,200), each=2), c=40, l=c(90,50))
      for (se in 1:2){
        plot(0, type='n', xaxt='n', ylim=c(-1.1,1.1), xlim=c(0.5,length(inds) + .5), xlab='', ylab='', yaxt=if (se == 2) 'n' else 's')
        polygon(c(-1,-1,9,9), c(0,-2,-2,0), border=NA, col=grey(0.9))
#        if (pseudoproxy){
#          boxplot(t(RE.mat[,se,1:29]), at=seq(inds), col=hcl(rep(seq(indices)/length(indices)*360, each=1), c=45, l=75),add=T, xaxt='n', yaxt='n', xlab='', ylab='', lwd=1, boxwex=0.5, range=0, lty=1)
#          points(seq(inds) - 0.35, RE.mat[,se,30], pch='>', lwd=3, cex=1.4)
#        } else {
        boxplot(t(RE.mat[,se,1:30]), at=seq(inds), col=hcl(rep(seq(indices)/length(indices)*360, each=1), c=45, l=75),add=T, xaxt='n', yaxt='n', xlab='', ylab='', lwd=1, boxwex=0.5, range=0, lty=1)
        points(seq(inds) - 0.35, RE.mat[,se,31], pch='>', lwd=3, cex=1.4)
 #       }
        text(seq(inds), rep(-0.9,3), ind.names[[se]], cex=1.2, adj=c(0.5,1))
        text(0.5,1.1,paste(pnames[se], ' Skill (RE) in ', c('October to March', 'April to September')[se], sep='') , cex=1.4, adj=c(0,1))
      }
      indind <- c(1,2,4,5,7,8,10,11)
      for (se in 1:2){
        plot(0, type='n', xaxt='n', ylim=c(-1.1,1.1), xlim=c(0,max(indind)+1), xlab='', ylab='', yaxt=if (se == 2) 'n' else 's')
        polygon(c(-1,-1,max(indind)+c(2,2)), c(0,-2,-2,0), border=NA, col=grey(0.9))
        boxplot(t(corr.mat[,se,1:30]), at=indind, col=hcl(rep(seq(indices)/length(indices)*360, each=2), c=c(30,60), l=c(95,50)),
                add=T, xaxt='n', yaxt='n', xlab='', ylab='', lwd=1, range=0, lty=1)
        points(seq(1.5,max(indind),3) - 1.2, corr.mat[seq(1,nrow(corr.mat),2),se,31], pch='>', lwd=3, cex=1.4)
        points(seq(1.5,max(indind),3) + 1.2, corr.mat[seq(2,nrow(corr.mat),2),se,31], pch='<', lwd=3, cex=1.4)
        text(seq(1.5,max(indind),3), rep(-0.9,3), ind.names[[se]], cex=1.2, adj=c(0.5,1))
        text(0,1.1,paste(pnames[se + 2], ' Correlation in ', c('October to March', 'April to September')[se], sep=''), cex=1.4, adj=c(0,1))
      }
    dev.off()
    } # end not pseudo_prox
  } # end if(echam_ind)
  ###############################################################
  
  
  
  
   
  # ###############################################################
  # ############ figs.16: Plot Distance Weight  ###################
  # ###############################################################
  # 
  # # plot distance weights
  # if (plot_dweights) {
  #   library(raster)
  #   pdf(paste(figpath,'distance_weights.pdf',sep='/'), width=4.5, height=4.5, 
  #       paper='special')
  #   image(d.weights_all)
  #   dev.off()
  # }
  
  
  ###############################################################
  
  
  
  ###############################################################
  ##### figs.15: Continuous Ranked Probability Score Maps  ######
  ###############################################################
  
    ### plot CRPS score difference ech-ana (positive values are good) with plot_echam: Temp
  if(CRPS){
    pdf(paste(figpath,'CRPS_new.pdf',sep='/'), width=9, height=6, paper='special')
    crps <- crps.ech-crps.ana
    crps.tot <- echam
    crps.tot$ensmean <- crps
    
    levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0.01,0.02,0.04,0.08,0.1,Inf)
    
    if (monthly_out){
      plot_echam4(crps.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(crps.tot$ensmean)[2]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason="summer",plotname=paste0('crps_anom_echam_anal-',v,'_temp_summer.pdf'), paper='special')
      plot_echam4(crps.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(crps.tot$ensmean)[2]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason="winter",plotname=paste0('crps_anom_echam_anal-',v,'_temp_winter.pdf'), paper='special')
    }else{
      plot_echam4(crps.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(crps.tot$ensmean)[2]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason=NULL,plotname=paste0('crps_anom_echam_anal-',v,'_temp.pdf'), paper='special')
    }
    
    
    
    
    ### plot CRPS score difference ech-ana (positive values are good) with plot_echam: Precip
    
    levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0.05,0.2,0.4,0.6,0.8,Inf)
    if (monthly_out) {
      plot_echam4(crps.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(crps.tot$ensmean)[2]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason="summer",plotname=paste0('crps_anom_echam_anal-',v,'_precip_summer.pdf'), paper='special')
      plot_echam4(crps.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(crps.tot$ensmean)[2]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason="winter",plotname=paste0('crps_anom_echam_anal-',v,'_precip_winter.pdf'), paper='special')
    }else{
      plot_echam4(crps.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(crps.tot$ensmean)[2]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason=NULL,plotname=paste0('crps_anom_echam_anal-',v,'_precip.pdf'), paper='special')
    }
    dev.off()
  }

# plot(density(RE.anom$ensmean[which(RE.anom$names  == varname),i],na.rm=T,weights=area_wgt[!is.na(RE.anom$ensmean[which(RE.anom$names  == varname),i])]))


  validate <- validate.init
} #end validation_set for loop
  
  #############################################################################
  
  ################################################################################
  ################## Appendix: Commented Plot parts moved to here ################
  ################################################################################
  
  
  if (loo) {  # quick plots to check loo
    #Nein, damit meine ich dass die Standardabweichung der Analyse an der ausgelassenen 
    #Gitterzelle dem mittleren quadrierten Fehler (korrigiert um den Beobachtungsfehler) 
    #entsprechen soll. Dies beinhaltet also eine Mittelung entweder über die Zeit für 
    #Aussagen pro Gitterbox, oder über den Raum (pro Zeitschritt). Ich denke man braucht
    #beide Varianten. Erstere zum Abschätzen ob die Analyse überall ungefähr gleich 
    #verlässlich ist, und Letztere zum Illustrieren des Effekts der sich ändernden 
    #Beobachtungsdichte.
    
    #Ja, oder andersrum: 
    #rmse ~= spread ; sqrt((y_loo – mean(x_loo,i))^2 – sigma_obs^2) = sqrt(var(x_loo,i)), 
    # wobei x_loo,i das ite ensemble member der LOO Analyse und  y_loo die ausgelassene 
    # Beobachtung ist.
    for (cyr in syr2:eyr){
      #cyr=1811
      load(paste0(dataintdir,'loo/loo_results_',cyr,'.Rdata'))
      print(dim(looensmeanres))
      if (cyr==syr2) {
        loodata <- cbind(paste(loodata[,1],loodata[,2],sep="-"),loodatares)
        looensmean <- cbind(paste(looensmean[,1],looensmean[,2],sep="-"),looensmeanres)
      } else {
        tmp <- merge(looensmean,looensmeanres,by=intersect(looensmean[,1],looensmeanres[,1]))
      }
      #lev <- pretty(looensmeanres[,3:14],11)
      lev <- c(-Inf,-2,-1,-0.5,0.5,1,2,Inf) #pretty(looensmeanres[,3:14],11)
      br <- length(lev)
      colpal <- two.colors(n=br,start="darkblue", end="darkred", alpha=0.5)
      colpal[4:5] <- "#FFFFFFFF"
      pdf(file=paste0(figpath,'/loo_',cyr,'.pdf'),width=5,height=36)
      set.panel(12,1)
      par(oma=c(.1,.1,.1,.1))
      for (i in 3:14){
        datcol <- colpal[as.numeric(cut(looensmeanres[,i],breaks=br))]
        #plot(NULL,xlim=c(-170,40),ylim=c(20,80),xlab="lon",ylab="lon")
        plot(NULL,xlim=c(-180,180),ylim=c(-60,80),xlab="lon",ylab="lon")
        points(looensmeanres[,1], looensmeanres[,2],cex=0.5,pch=15,col=datcol,
               xlab="Longitude",ylab="Latitude")                # point fill color
        map("world",interior=F,add=T,ty='l',col='black',
            xlim=c(-180,180),ylim=c(-60,80))
        #xlim=c(-170,40),ylim=c(20,80))
        legend("bottomleft", inset=0.01, as.character(lev),fill=colpal,bty="o",
               bg="white",box.col="white",box.lwd=0,cex=0.7)
      }
      dev.off()
    }
  }
  
  
  # simple usage of jonas plotting functions:
  # plot_echam(bias.ech, varname='temp2', type='ensmean', cex.pt=1.5, ti=1)
  # plot_echam(bias$'Localized analysis', varname='temp2', type='ensmean', cex.pt=2, ti=1)
  # plot_echam(echam.anom, varname='temp2', type='ensmean', cex.pt=2, ti=1)
  # par(mar=c(0,0,0,0))
  # map("world",interior=F,col='darkgrey') 
  # map.axes()
  # title("")
  # #map.scale(-180, -50, relwidth = 0.15, metric=T, ratio=F,cex=0.5) # ACHTUNG, wenn Projektion nicht längentreu
  # points(calibrate$lon, calibrate$lat, main='1943 station locations',pch='.',col='red')
  # points(calibrate$lon[which(calibrate$sour=="prox")], 
  #       calibrate$lat[which(calibrate$sour=="prox")],pch=1,col='red',cex=0.5)
  # points(analysis$lon[H.i[,1]],analysis$lat[H.i[,1]],col='blue',cex=0.5,pch=20) # nur die proxies sind in H.i
  
  #######################SampleYrPlotsCombined#######################
  # for (t in 1:60) {
  # plotdata=echam
  # plotdata$data <- array(c(echam.anom$ensmean[,t],analysis.anom$ensmean[,t]), c(nrow(echam.anom$ensmean),1,2))
  # plotdata$names <- c(echam.anom$names,analysis.anom$names)
  # plotdata$data <- array(c(echam$ensmean[,t],analysis$ensmean[,t]), c(nrow(echam$ensmean),1,2))
  # 
  # pdf(paste(figpath,'/example_year_',plotdata$time[t],'_temp_anom.pdf',sep=''), width=9, height=4.5, paper='special')
  # layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
  # par(oma=c(0,0,0,0))
  # levs <- c(-Inf, seq(-1.2,1.2,0.2), Inf)
  # plot_echam4(plotdata, varname='temp2', type='data', cex.pt=1.5, names=pnames[1:dim(plotdata$data)[3]], lev=levs, st.col=NULL, stations=NULL, add=T,plotname = paste0('example_year_',plotdata$time[t],'_temp_anom.pdf'),paper = 'special')
  # dev.off()
  # }
  


