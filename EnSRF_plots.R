# use EnSRF_prepplot_v???.R to load the data

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


rm(list=ls())


syrtot=1901
eyrtot=2000

syr=1902
eyr=2000

user <- system("echo $USER",intern=T)
print(paste('User:',user))
if (user=="veronika") {
  # workdir('/scratch/veronika/rerun/r_code')
  workdir ='/scratch3/veronika/reuse/reuse_git/' # where are the scripts from github
} else if (user=="joerg") {
  workdir='/scratch3/joerg/projects/reuse/reuse_git/'
} else if (user=="lucaf") {
  workdir='/scratch3/lucaf/reuse/reuse_git/'
} else if (user == "nevin"){
  workdir = '/scratch3/nevin/reuse/reuse_git/'
} else {
  stop("Unknown user!")

}
dataextdir='/mnt/climstor/giub/EKF400/'
dataintdir=paste0(workdir,'../data/')
setwd(workdir)


source('EnSRF_switches.R')
source('EnSRF_functions.R')

statyr=1905

if (monthly_out) {
  if(calc_prepplot){
    dir.create(paste0("../data/prepplot/",expname))
    dir.create(paste0("../data/prepplot/",expname,'/prepplot_monthly/'))
  }
  prepplotdir=paste0("../data/prepplot/",expname,'/prepplot_monthly/')
  nseas=12
} else {
  if(calc_prepplot){
    dir.create(paste0("../data/prepplot/",expname))
    dir.create(paste0("../data/prepplot/",expname,'/prepplot_seasonal/'))
  }
  prepplotdir=paste0("../data/prepplot/",expname,'/prepplot_seasonal/') 
  nseas=2
}
figpath=paste0('../figures/',expname,'_',syr,'-',eyr) #format(Sys.time(), "%Y%m%d_%H-%M_")
dir.create(figpath)


pnames <- paste(c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k'), ')', sep='')
# the plots will be reordered in future
################################################################################
# Fig. 0: station locations 
#
if (countseries) {
  # count number of assimilation data of each type and add to plot
  nmxd <- ntrw <- nprox <- ndoc <-ninstslp <-ninsttemp <- ninst <- rep(NA,length(syr:eyr))
  i <- 1
  for (cyr in syr:eyr) {
    print(cyr)
    if (every2grid) {
      load(file=paste0(prepplotdir,'analysis_',cyr,'_2ndgrid.Rdata')) 
    } else {
      load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata')) 
    }
    ninst[i] <- length(which(calibrate$sour=="inst"))
    ninsttemp[i] <- length(which(calibrate$sour=="inst"&calibrate$names=="temp2"))
    ninstslp[i] <- length(which(calibrate$sour=="inst"&calibrate$names=="slp"))
    ndoc[i] <- length(which(calibrate$sour=="doc"))
    nprox[i] <- length(which(calibrate$sour=="prox"))-length(which(apply(is.na(calibrate$data[which(calibrate$sour=="prox"),]),1,all)))
    ctrw <- 0
    cmxd <- 0
    if (!is.null(nrow(calibrate$mr))) {
      for (j in 1:nrow(calibrate$mr)) {
        #print(c(j,length(which(!is.na(calibrate$mr[j,])))))
        if (length(which(!is.na(calibrate$mr[j,])))==5) cmxd <- cmxd+1
        if (length(which(!is.na(calibrate$mr[j,])))==8) ctrw <- ctrw+1
      }
      ntrw[i] <- ctrw
      nmxd[i] <- cmxd
    } else {
      ntrw[i] <- 0
      nmxd[i] <- 0
    }
    i <- i+1
  }
  nrecords <- cbind((syr:eyr),ninst,ninsttemp,ninstslp,ndoc,nprox,ntrw,nmxd,rep(0,length(syr:eyr)))
  #nrecords[41,] <- nrecords[40,] <- nrecords[39,] <- nrecords[38,] # check input data 1639-41!
  pdf(paste(figpath,'record_no.pdf',sep='/'), width=9, height=3.5, paper='special')
  par(oma=c(0,0,0,0),mar=c(4,4,0.4,0.4))
  
  plot(nrecords[,1],nrecords[,9],ty="l",col='white',ylim=c(0,5000),xlab="year",ylab="No of records")
  polygon(c(nrecords[,1],rev(nrecords[,1])),c(nrecords[,6],rev(nrecords[,9])),col="seagreen3")
  polygon(c(nrecords[,1],rev(nrecords[,1])),c((nrecords[,7]+nrecords[,6]),rev((nrecords[,6]))),col="salmon4")
  polygon(c(nrecords[,1],rev(nrecords[,1])),c((nrecords[,8]+nrecords[,7]+nrecords[,6]),rev((nrecords[,7]+nrecords[,6]))),col="royalblue")
  
  polygon(c(nrecords[,1],rev(nrecords[,1])),c((nrecords[,5]+nrecords[,8]+nrecords[,7]+nrecords[,6]),rev((nrecords[,8]+nrecords[,7]+nrecords[,6]))),col="plum2")
  polygon(c(nrecords[,1],rev(nrecords[,1])), c((nrecords[,3]+nrecords[,5]+nrecords[,8]+nrecords[,7]+nrecords[,6]),rev((nrecords[,5]+nrecords[,8]+nrecords[,7]+nrecords[,6]))),col="firebrick1")
  polygon(c(nrecords[,1],rev(nrecords[,1])), c((nrecords[,4]+nrecords[,3]+nrecords[,5]+nrecords[,8]+nrecords[,7]+nrecords[,6]),rev((nrecords[,3]+nrecords[,5]+nrecords[,8]+nrecords[,7]+nrecords[,6]))),col="goldenrod1")
  legend("topleft", c("Instr. SLP","Instr. temp.", "Docum. temp.", 'MXD','TRW','Proxies'), 
         pch=rep(15,6), col=c("goldenrod1","firebrick1", "plum2", "royalblue", "salmon4", "seagreen3"), pt.cex=1, pt.lwd=1, 
         inset=0.005, bg='transparent', box.col='transparent', cex=1)
  dev.off()
  
  
  #statsyr=1790
  #stateyr=1804
  if (every2grid) {
    load(file=paste0(prepplotdir,'analysis_',statyr,'_2ndgrid.Rdata')) 
  } else {
    load(file=paste0(prepplotdir,'analysis_',statyr,'.Rdata')) 
  }
  #statll <- unique(cbind(calibrate$lon,calibrate$lat))
  #for (statyr in seq(from=(statsyr+1),to=stateyr)) {
  #  load(paste0(datadir,"EnSRF_analysis/prepplot_v3_seasonal/analysis_",statyr,".Rdata"))
  #  statll <- unique(rbind(statll,unique(cbind(calibrate$lon,calibrate$lat))))
  #}
  dat <- echam
  #dat$data <- array(corr$ensmean*0, c(dim(corr$ensmean)[1], 1, dim(corr$ensmean)[2]))
  dat$data <- array(echam$ensmean[1:4608,]*0, c(4608, 1, 2))
  dat$lat <- dat$lat[1:dim(dat$data)[1]]
  dat$lon <- dat$lon[1:dim(dat$data)[1]]
  dat$names <- dat$names[1:dim(dat$data)[1]]
  pdf(paste0(figpath,'/stat_locations_',statyr,'.pdf'), width=9, height=3.5, paper='special')
  #layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
  layout(matrix(c(1,2), 1, 2, byrow = TRUE))
  par(oma=c(0,0,0,0))
  levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0.05,0.2,0.4,0.6,0.8,1)
  plot_echam(dat,type='data',cex.pt=1.5,names=pnames[1:dim(dat$data)[3]],lev=levs,colorbar=F,
             stations=calibrate,st.col=NULL,statpanel=c(1,2),add=T,wcol='darkgrey')
  dev.off()
}


################VALIDATION PLOTS################


for (v in validation_set){
  
  ## CHOOSE THE VALIDATION DATA SET IN SWITCHES (validation_set), PLOTS SCRIPT HAS TO BE RUN FOR VALIDATION SETS SEPARATELY
  if (every2grid) {
    if (monthly_out) {
      load(file=paste0("../data/image/",expname,"/prepplot_",v,"_calc_vali_stat_image_",
                       syr,"-",eyr,"_monthly_2ndgrid.Rdata"))
      load(file=paste0("../data/image/",expname,"/indices_tot_",syrtot,"-",eyrtot,"_monthly_2ndgrid.Rdata"))
    } else {
      load(file=paste0("../data/image/",expname,"/prepplot_",v,"_calc_vali_stat_image_",
                       syr,"-",eyr,"_seasonal_2ndgrid.Rdata"))
      load(file=paste0("../data/image/",expname,"/indices_tot_",syrtot,"-",eyrtot,"_seasonal_2ndgrid.Rdata"))
    }  
  } else {
    if (monthly_out) {
      load(file=paste("../data/image/",expname,"/prepplot_",v,"_calc_vali_stat_image_",
                      syr,"-",eyr,"_monthly.Rdata",sep=""))
      load(file=paste0("../data/image/",expname,"/indices_tot_",syrtot,"-",eyrtot,"_monthly.Rdata"))
    } else {
      load(file=paste("../data/image/",expname,"/prepplot_",v,"_calc_vali_stat_image_",
                      syr,"-",eyr,"_seasonal.Rdata",sep=""))
      load(file=paste0("../data/image/",expname,"/indices_tot_",syrtot,"-",eyrtot,"_seasonal.Rdata"))
    }
  }
  validate.init <- validate
  vind.init<-vind
  vind.tot.init<-vind.tot
  validate.anom.init<-validate.anom
  
  validate <- validate.init[[v]]
  vind<-vind.init[[v]]
  vind.tot<-vind.tot.init[[v]]
  validate.anom<-validate.anom.init[[v]]  
  
  data.dim <- dim(echam$data)[c(1,2,2,3)]
  data.dim[2:3] <- c(s.plot,data.dim[3]/s.plot)
  ens.dim <- c(nrow(echam$ensmean), s.plot, ncol(echam$ensmean)/s.plot)
  print(paste0("monthly_out: ",monthly_out))
  #####################################################################################
  #if (recalc){
  #  source('EnSRF_data.R')
  #  source('EnSRF_prepplots.R')
  #} 
  #if (reload){
  #  load(paste("../data/EnSRF_",syr,"-",eyr,".Rdata",sep=""))
  #}
  
  if ((monthly)&(!sixmonstatevector)) { 
    # just plot jan und jul
    ana.spread <- ana.spread[,c(1,7)]
    ech.spread <- ech.spread[,c(1,7)]
    corr.ech$ensmean <- corr.ech$ensmean[,c(1,7)]
    corr$ensmean <- corr$ensmean[,c(1,7)]
    bias.ech$ensmean <- bias.ech$ensmean[,c(1,7)]
    bias$ensmean <- bias$ensmean[,c(1,7)]
    RE.anom$ensmean <- RE.anom$ensmean[,c(1,7)]
    RE$ensmean <- RE$ensmean[,c(1,7)]
  }  
  
  if ((monthly)&(sixmonstatevector)&(monthly_out)) { 
    # just plot jan und jul
    if (!landcorrected){
      ana.spread <- ana.spread[,c(4,10)]
      ech.spread <- ech.spread[,c(4,10)]
    }
    corr.ech$ensmean <- corr.ech$ensmean[,c(4,10)]
    corr$ensmean <- corr$ensmean[,c(4,10)]
    bias.ech$ensmean <- bias.ech$ensmean[,c(4,10)]
    bias$ensmean <- bias$ensmean[,c(4,10)]
    RE.anom$ensmean <- RE.anom$ensmean[,c(4,10)]
    RE$ensmean <- RE$ensmean[,c(4,10)]
  }  
  
  
  
  ################################################################################
  # plot example year anomalies
  #plot_echam(validate, varname='temp2', type='ensmean', ti=33, lty=3, add=F)
  #plot_echam(analysis[[2]], varname='temp2', type='ensmean', ti=33, lty=3, add=F)
  #plot_echam(analysis.anom, varname='temp2', type='ensmean', ti=32, lty=3, add=F)
  #plot_echam(echam, varname='temp2', type='ensmean', ti=32, lty=3, add=F)
  #plot_echam(echam.anom, varname='temp2', type='ensmean', ti=33, lty=3, add=F)
  #plotdata <- echam
  
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
  
  
  
  
  
  
  
  # plot sample time series for Norway lon/lat:
  # plot sample time series for Siberia lon/lat:
  
  # 
  # stat.pos1 <- which(calibrate$lon>50&calibrate$lon<70&calibrate$lat>60&calibrate$lat<70) #sibiria
  # 
  # for (i in 1:length(stat.pos1)){
  #   if (monthly_out) {
  #     pdf(paste(figpath,'/example_timeseries_sibiria_',v,'_mon',i,'.pdf',sep=''), 
  #         width=4.5, height=6, paper='special')
  #   }else if(yearly_out){
  #     pdf(paste(figpath,'/example_timeseries_sibiria_',v,'_yrly',i,'.pdf',sep=''), 
  #         width=4.5, height=6, paper='special')
  #   } else {
  #     pdf(paste(figpath,'/example_timeseries_sibiria_',v,'_',i,'.pdf',sep=''), 
  #         width=7, height=6, paper='special')
  #   }
  #   par(oma=c(0,0,0,0),mar=c(2,4,2,0.5),mfrow=c(3,1))
  #   stat.pos<-stat.pos1[i]  
  #   pos <- getgridboxnum(calibrate,echam) [stat.pos]
  #   
  #   
  #   par(oma=c(0,0,0,0),mar=c(2,4,2,0.5),mfrow=c(3,1))
  #   
  #   
  #   # stat.pos <- which(calibrate$lon>5&calibrate$lon<30&calibrate$lat>56&calibrate$lat<70) #scandinavia
  #   if (monthly_out) {
  #     #3yr monthly temp
  #     pos2 <- seq(1,36,1)
  #     period <- as.Date(paste(rep(seq(from=1950,to=1952),each=12),rep(seq(1,12),3),rep(15,36),sep='-'))
  #   }else if (yearly_out){
  #     pos2<-seq(1,ncol(validate$data),1)
  #     period<-validate$time
  #   } else {
  #     #30yr summer temp
  #     pos2 <- seq(2,ncol(validate$data),2)
  #     period <- validate$time[pos2]
  #   }
  #   plot(period,validate$data[pos,pos2],ylim=c(min(cbind(echam$data[pos,pos2,],validate$data[pos,pos2]),na.rm=T),
  #                                              max(cbind(echam$data[pos,pos2,],validate$data[pos,pos2]),na.rm=T)),ty='l',col=rgb(0,0,10,10,maxColorValue=10),
  #        main="Temperature",xlab='',ylab='ºC',xaxt='n')
  #   polygon(c(period,rev(period)),c(apply(echam$data[pos,pos2,],1,max),
  #                                   rev(apply(echam$data[pos,pos2,],1,min))),density=NA,
  #           col=rgb(1,1,1,3,maxColorValue=10))
  #   lines(period,echam$ensmean[pos,pos2],col=rgb(0,0,0,10,maxColorValue=10))
  #   polygon(c(period,rev(period)),c(apply(analysis$data[pos,pos2,],1,max),
  #                                   rev(apply(analysis$data[pos,pos2,],1,min))),density=NA,
  #           col=rgb(10,0,0,3,maxColorValue=10))
  #   lines(period,analysis$ensmean[pos,pos2],col=rgb(10,0,0,10,maxColorValue=10))
  #   lines(period,validate$data[pos,pos2],ylim=c(min(echam$data[pos,pos2,],na.rm=T),
  #                                               max(echam$data[pos,pos2,],na.rm=T)),col=rgb(0,0,10,10,maxColorValue=10))
  #   
  #   # plot(stat$time[seq(2,48,2)],stat$data[pos,][201,seq(2,48,2)],ylim=c(min(stat$data[201,pos2]),
  #   #  max(stat$data[201,pos2])),ty='l',col=rgb(0,0,10,10,maxColorValue=10),
  #   #      main="Temperature",xlab='',ylab='ºC',xaxt='n')
  #   # lines(stat$time[seq(2,48,2)],stat$data[pos,][216,seq(2,48,2)],col=rgb(10,0,0,10,maxColorValue=10))
  #   # lines(stat$time[seq(2,48,2)],stat$data[pos,][219,seq(2,48,2)],col=rgb(0,10,0,10,maxColorValue=10))
  #   # lines(analysis.anom$time[seq(2,48,2)],analysis.anom$ensmean[676,seq(2,48,2)],
  #   #       col=rgb(10,0,0,10,maxColorValue=10),ty='l')
  #   #  legend("topleft", c('Instrumental CRU TS3', 'CCC400', "EFK400"),col=c("blue", "black", "red"), 
  #   #         lty=c(1,1,1),lwd=c(1,1,1),pt.cex=1, pt.lwd=1,inset=0.005, bg='white', 
  #   #         box.col='white', cex=1)
  #   #precip
  #   plot(period,validate$data[4608+pos,pos2],ylim=c(min(echam$data[4608+pos,pos2,],na.rm=T),
  #                                                   max(echam$data[4608+pos,pos2,],na.rm=T)),ty='l',col=rgb(0,0,10,10,maxColorValue=10),
  #        main="Precipitation",xlab='',ylab='mm',xaxt='n')
  #   polygon(c(period,rev(period)),c(apply(echam$data[4608+pos,pos2,],1,max),
  #                                   rev(apply(echam$data[4608+pos,pos2,],1,min))),density=NA,
  #           col=rgb(1,1,1,3,maxColorValue=10))
  #   lines(period,echam$ensmean[4608+pos,pos2],col=rgb(0,0,0,10,maxColorValue=10))
  #   polygon(c(period,rev(period)),c(apply(analysis$data[4608+pos,pos2,],1,max),
  #                                   rev(apply(analysis$data[4608+pos,pos2,],1,min))),density=NA,
  #           col=rgb(10,0,0,3,maxColorValue=10))
  #   lines(period,analysis$ensmean[4608+pos,pos2],col=rgb(10,0,0,10,maxColorValue=10))
  #   lines(period,validate$data[4608+pos,pos2],ylim=c(min(echam$data[pos,pos2,],na.rm=T),
  #                                                    max(echam$data[4608+pos,pos2,],na.rm=T)),col=rgb(0,0,10,10,maxColorValue=10))
  #   lines(period,calibrate$data[stat.pos,pos2]+74,col="yellow")
  #   # lines(period,apply(calibrate$data[stat.pos,pos2],2,mean)+74,col="yellow")
  #   legend("topleft", c('CRU TS3', 'CCC400', "EFK400"),col=c("blue", "black", "red"), 
  #          lty=c(1,1,1),lwd=c(1,1,1),pt.cex=1, pt.lwd=1,inset=0.005, bg='white', 
  #          box.col='white', cex=1)
  #   #slp
  #   plot(period,validate$data[2*4608+pos,pos2],ylim=c(min(echam$data[2*4608+pos,pos2,],na.rm=T),
  #                                                     max(echam$data[2*4608+pos,pos2,],na.rm=T)),ty='l',col=rgb(0,0,10,10,maxColorValue=10),
  #        main="SLP",xlab='',ylab='hPa')
  #   polygon(c(period,rev(period)),c(apply(echam$data[2*4608+pos,pos2,],1,max),
  #                                   rev(apply(echam$data[2*4608+pos,pos2,],1,min))),density=NA,
  #           col=rgb(1,1,1,3,maxColorValue=10))
  #   lines(period,echam$ensmean[2*4608+pos,pos2],col=rgb(0,0,0,10,maxColorValue=10))
  #   polygon(c(period,rev(period)),c(apply(analysis$data[2*4608+pos,pos2,],1,max),
  #                                   rev(apply(analysis$data[2*4608+pos,pos2,],1,min))),density=NA,
  #           col=rgb(10,0,0,3,maxColorValue=10))
  #   lines(period,analysis$ensmean[2*4608+pos,pos2],col=rgb(10,0,0,10,maxColorValue=10))
  #   lines(period,validate$data[2*4608+pos,pos2],ylim=c(min(echam$data[pos,pos2,],na.rm=T),
  #                                                      max(echam$data[2*4608+pos,pos2,],na.rm=T)),col=rgb(0,0,10,10,maxColorValue=10))
  #   # legend("bottomleft", c('Instrumental CRU TS3', 'CCC400', "EFK400"),col=c("blue", "black", "red"), 
  #   #   lty=c(1,1,1),lwd=c(1,1,1),pt.cex=1, pt.lwd=1,inset=0.005, bg='white', 
  #   #   box.col='white', cex=1)
  #   dev.off()
  # }
  
  
  ###############################################################
  ############ plot sample time series for lon/lat: #############
  ###############################################################
  
  #codeline beneath finds gridbox with lowest correlation (was used because of strange bad results)
  #which(corr$ensmean[,2]==min(corr$ensmean[which(corr$lat<34&corr$lat>24&corr$lon<(90)&corr$lon>(68)&corr$names=="temp2"),2],na.rm=T))
  if (v=="cru_vali"){
    valilegend<-"CRU"
  }else if(v=="twentycr_vali"){
    valilegend<-"20CR"
  }
  
  if (monthly_out) {
    #3yr monthly temp
    pos2 <- seq(1,36,1)
    period <- as.Date(paste(rep(seq(from=1950,to=1952),each=12),rep(seq(1,12),3),rep(15,36),sep='-'))
  }else if (yearly_out){
    pos2<-seq(1,ncol(validate$data),1)
    period<-validate$time
  } else {
    #30yr summer temp
    pos2 <- seq(2,ncol(validate$data),2)
    period <- validate$time[pos2]
  }
  
  
  ################### absolute value s#####################
  
  # Siberia
  paste(validate$lon[592],validate$lat[592])
  if (monthly_out) {
    pdf(paste(figpath,'/example_timeseries_sibiria_',v,'_mon.pdf',sep=''), 
        width=4.5, height=6, paper='special')
  }else if(yearly_out){
    pdf(paste(figpath,'/example_timeseries_sibiria_',v,'_yrly.pdf',sep=''), 
        width=7, height=6, paper='special')
  } else {
    pdf(paste(figpath,'/example_timeseries_sibiria_',v,'.pdf',sep=''), 
        width=7, height=6, paper='special')
  }
  par(oma=c(0,0,0,0),mar=c(2,4,2,2),mfrow=c(3,1))
  #50yr summer temp
  plot(echam$time[pos2],echam$ensmean[592,pos2],ylim=c(0,9),ty='l',col="black",main="Temperature",
       xlab='',ylab='[ºC]',xaxt='n',axes=F)
  
  polygon(c(echam$time[pos2],rev(echam$time[pos2])),c(apply(echam$data[592,pos2,],1,max),
                                                      rev(apply(echam$data[592,pos2,],1,min))),density=NA, col=rgb(1,1,1,3,maxColorValue=10))
  lines(analysis$time[pos2],analysis$ensmean[592,pos2],col="red")
  polygon(c(analysis$time[pos2],rev(analysis$time[pos2])),c(apply(analysis$data[592,pos2,],1,max),
                                                            rev(apply(analysis$data[592,pos2,],1,min))),density=NA, col=rgb(10,0,0,3,maxColorValue=10))
  
  axis(2,col="red")
  par(new=T)
  plot(validate$data[592,pos2],ylim=c(1.5,10.5),ty='l',col="blue",main="",
       xlab='',ylab='[ºC]',xaxt='n',axes=F)
  axis(4, col="blue")
  box()
  #precip
  
  plot(echam$time[pos2],echam$ensmean[4608+592,pos2],ylim=c(0,100),ty='l',col="black",main="Precipitation",
       xlab='',ylab='[mm]',xaxt='n',axes=F)
  
  polygon(c(echam$time[pos2],rev(echam$time[pos2])),c(apply(echam$data[4608+592,pos2,],1,max),
                                                      rev(apply(echam$data[4608+592,pos2,],1,min))),density=NA, col=rgb(1,1,1,3,maxColorValue=10))
  lines(analysis$time[pos2],analysis$ensmean[4608+592,pos2],col="red")
  polygon(c(analysis$time[pos2],rev(analysis$time[pos2])),c(apply(analysis$data[4608+592,pos2,],1,max),
                                                            rev(apply(analysis$data[4608+592,pos2,],1,min))),density=NA, col=rgb(10,0,0,3,maxColorValue=10))
  
  axis(2,col="red")
  par(new=T)
  plot(validate$data[4608+592,pos2],ylim=c(-5,95),ty='l',col="blue",main="",
       xlab='',ylab='[mm]',xaxt='n',axes=F)
  axis(4, col="blue")
  box()
  #slp
  
  plot(echam$time[pos2],echam$ensmean[2*4608+592,pos2],ylim=c(1006,1020),ty='l',col="black",main="Sea level pressure",
       xlab='',ylab='[hPa]',xaxt='n',axes=F)
  
  polygon(c(echam$time[pos2],rev(echam$time[pos2])),c(apply(echam$data[2*4608+592,pos2,],1,max),
                                                      rev(apply(echam$data[2*4608+592,pos2,],1,min))),density=NA, col=rgb(1,1,1,3,maxColorValue=10))
  lines(analysis$time[pos2],analysis$ensmean[2*4608+592,pos2],col="red")
  polygon(c(analysis$time[pos2],rev(analysis$time[pos2])),c(apply(analysis$data[2*4608+592,pos2,],1,max),
                                                            rev(apply(analysis$data[2*4608+592,pos2,],1,min))),density=NA, col=rgb(10,0,0,3,maxColorValue=10))
  
  axis(2,col="red")
  par(new=T)
  plot(validate$time[pos2], validate$data[2*4608+592,pos2],ylim=c(1005,1019),ty='l',col="blue",main="",
       xlab='',ylab='[hPa]',xaxt='n',axes=F)
  axis(4, col="blue")
  axis(1)
  box()
  
  legend("bottomleft", c(valilegend, 'CCC400', "EFK400"),col=c("blue", "black", "red"),
         lty=c(1,1,1),lwd=c(1,1,1),pt.cex=1, pt.lwd=1,inset=0.005, bg='transparent',
         box.col='transparent', cex=1)
  dev.off()
  
  # Norway
  paste(validate$lon[676],validate$lat[676])
  
  if (monthly_out) {
    pdf(paste(figpath,'/example_timeseries_norway_',v,'_mon.pdf',sep=''), 
        width=4.5, height=6, paper='special')
  }else if(yearly_out){
    pdf(paste(figpath,'/example_timeseries_norway_',v,'_yrly.pdf',sep=''), 
        width=7, height=6, paper='special')
  } else {
    pdf(paste(figpath,'/example_timeseries_norway_',v,'.pdf',sep=''), 
        width=7, height=6, paper='special')
  }
  par(oma=c(0,0,0,0),mar=c(2,4,2,2),mfrow=c(3,1))
  #50yr summer temp
  plot(echam$time[pos2],echam$ensmean[676,pos2],ylim=c(1,10),ty='l',col="black",main="Temperature",
       xlab='',ylab='[ºC]',xaxt='n',axes=F)
  
  polygon(c(echam$time[pos2],rev(echam$time[pos2])),c(apply(echam$data[676,pos2,],1,max),
                                                      rev(apply(echam$data[676,pos2,],1,min))),density=NA, col=rgb(1,1,1,3,maxColorValue=10))
  lines(analysis$time[pos2],analysis$ensmean[676,pos2],col="red")
  polygon(c(analysis$time[pos2],rev(analysis$time[pos2])),c(apply(analysis$data[676,pos2,],1,max),
                                                            rev(apply(analysis$data[676,pos2,],1,min))),density=NA, col=rgb(10,0,0,3,maxColorValue=10))
  
  axis(2,col="red")
  par(new=T)
  plot(validate$data[676,pos2],ylim=c(3,12),ty='l',col="blue",main="",
       xlab='',ylab='[ºC]',xaxt='n',axes=F)
  axis(4, col="blue")
  box()
  #precip
  
  plot(echam$time[pos2],echam$ensmean[4608+676,pos2],ylim=c(30,150),ty='l',col="black",main="Precipitation",
       xlab='',ylab='[mm]',xaxt='n',axes=F)
  
  polygon(c(echam$time[pos2],rev(echam$time[pos2])),c(apply(echam$data[4608+676,pos2,],1,max),
                                                      rev(apply(echam$data[4608+676,pos2,],1,min))),density=NA, col=rgb(1,1,1,3,maxColorValue=10))
  lines(analysis$time[pos2],analysis$ensmean[4608+676,pos2],col="red")
  polygon(c(analysis$time[pos2],rev(analysis$time[pos2])),c(apply(analysis$data[4608+676,pos2,],1,max),
                                                            rev(apply(analysis$data[4608+676,pos2,],1,min))),density=NA, col=rgb(10,0,0,3,maxColorValue=10))
  
  axis(2,col="red")
  par(new=T)
  plot(validate$data[4608+676,pos2],ylim=c(21.5,141.5),ty='l',col="blue",main="",
       xlab='',ylab='[mm]',xaxt='n',axes=F)
  axis(4, col="blue")
  box()
  #slp
  
  plot(echam$time[pos2],echam$ensmean[2*4608+676,pos2],ylim=c(1005,1018),ty='l',col="black",main="Sea level pressure",
       xlab='',ylab='[hPa]',xaxt='n',axes=F)
  
  polygon(c(echam$time[pos2],rev(echam$time[pos2])),c(apply(echam$data[2*4608+676,pos2,],1,max),
                                                      rev(apply(echam$data[2*4608+676,pos2,],1,min))),density=NA, col=rgb(1,1,1,3,maxColorValue=10))
  lines(analysis$time[pos2],analysis$ensmean[2*4608+676,pos2],col="red")
  polygon(c(analysis$time[pos2],rev(analysis$time[pos2])),c(apply(analysis$data[2*4608+676,pos2,],1,max),
                                                            rev(apply(analysis$data[2*4608+676,pos2,],1,min))),density=NA, col=rgb(10,0,0,3,maxColorValue=10))
  
  axis(2,col="red")
  par(new=T)
  plot(validate$time[pos2], validate$data[2*4608+676,pos2],ylim=c(1005,1018),ty='l',col="blue",main="",
       xlab='',ylab='[hPa]',xaxt='n',axes=F)
  axis(4, col="blue")
  axis(1)
  box()
  
  legend("bottomleft", c(valilegend, 'CCC400', "EFK400"),col=c("blue", "black", "red"),
         lty=c(1,1,1),lwd=c(1,1,1),pt.cex=1, pt.lwd=1,inset=0.005, bg='transparent',
         box.col='transparent', cex=1)
  dev.off()
  
  
  # Greenland
  paste(validate$lon[278],validate$lat[278])
  if (monthly_out) {
    pdf(paste(figpath,'/example_timeseries_greenland_',v,'_mon.pdf',sep=''), 
        width=4.5, height=6, paper='special')
  }else if(yearly_out){
    pdf(paste(figpath,'/example_timeseries_greenland_',v,'_yrly.pdf',sep=''), 
        width=7, height=6, paper='special')
  } else {
    pdf(paste(figpath,'/example_timeseries_greenland_',v,'.pdf',sep=''), 
        width=7, height=6, paper='special')
  }
  par(oma=c(0,0,0,0),mar=c(2,4,2,2),mfrow=c(3,1))
  #50yr summer temp
  plot(echam$time[pos2],echam$ensmean[278,pos2],ylim=c(-19,-8),ty='l',col="black",main="Temperature",
       xlab='',ylab='[ºC]',xaxt='n',axes=F)
  
  polygon(c(echam$time[pos2],rev(echam$time[pos2])),c(apply(echam$data[278,pos2,],1,max),
                                                      rev(apply(echam$data[278,pos2,],1,min))),density=NA, col=rgb(1,1,1,3,maxColorValue=10))
  lines(analysis$time[pos2],analysis$ensmean[278,pos2],col="red")
  polygon(c(analysis$time[pos2],rev(analysis$time[pos2])),c(apply(analysis$data[278,pos2,],1,max),
                                                            rev(apply(analysis$data[278,pos2,],1,min))),density=NA, col=rgb(10,0,0,3,maxColorValue=10))
  
  axis(2,col="red")
  par(new=T)
  plot(validate$data[278,pos2],ylim=c(-21.5,-10.5),ty='l',col="blue",main="",
       xlab='',ylab='[ºC]',xaxt='n',axes=F)
  axis(4, col="blue")
  box()
  #precip
  
  plot(echam$time[pos2],echam$ensmean[4608+278,pos2],ylim=c(0,30),ty='l',col="black",main="Precipitation",
       xlab='',ylab='[mm]',xaxt='n',axes=F)
  
  polygon(c(echam$time[pos2],rev(echam$time[pos2])),c(apply(echam$data[4608+278,pos2,],1,max),
                                                      rev(apply(echam$data[4608+278,pos2,],1,min))),density=NA, col=rgb(1,1,1,3,maxColorValue=10))
  lines(analysis$time[pos2],analysis$ensmean[4608+278,pos2],col="red")
  polygon(c(analysis$time[pos2],rev(analysis$time[pos2])),c(apply(analysis$data[4608+278,pos2,],1,max),
                                                            rev(apply(analysis$data[4608+278,pos2,],1,min))),density=NA, col=rgb(10,0,0,3,maxColorValue=10))
  
  axis(2,col="red")
  par(new=T)
  plot(validate$data[4608+278,pos2],ylim=c(6,36),ty='l',col="blue",main="",
       xlab='',ylab='[mm]',xaxt='n',axes=F)
  axis(4, col="blue")
  box()
  #slp
  
  plot(echam$time[pos2],echam$ensmean[2*4608+278,pos2],ylim=c(1010,1023),ty='l',col="black",main="Sea level pressure",
       xlab='',ylab='[hPa]',xaxt='n',axes=F)
  
  polygon(c(echam$time[pos2],rev(echam$time[pos2])),c(apply(echam$data[2*4608+278,pos2,],1,max),
                                                      rev(apply(echam$data[2*4608+278,pos2,],1,min))),density=NA, col=rgb(1,1,1,3,maxColorValue=10))
  lines(analysis$time[pos2],analysis$ensmean[2*4608+278,pos2],col="red")
  polygon(c(analysis$time[pos2],rev(analysis$time[pos2])),c(apply(analysis$data[2*4608+278,pos2,],1,max),
                                                            rev(apply(analysis$data[2*4608+278,pos2,],1,min))),density=NA, col=rgb(10,0,0,3,maxColorValue=10))
  
  axis(2,col="red")
  par(new=T)
  plot(validate$time[pos2], validate$data[2*4608+278,pos2],ylim=c(1008,1021),ty='l',col="blue",main="",
       xlab='',ylab='[hPa]',xaxt='n',axes=F)
  axis(4, col="blue")
  axis(1)
  box()
  
  legend("bottomleft", c(valilegend, 'CCC400', "EFK400"),col=c("blue", "black", "red"),
         lty=c(1,1,1),lwd=c(1,1,1),pt.cex=1, pt.lwd=1,inset=0.005, bg='transparent',
         box.col='transparent', cex=1)
  dev.off()
  
  
  # Alaska
  # line beneath was used to determine closest lon/lat to validate lat/lon 632
  #calipos<-which(abs(calibrate$lon-validate$lon[632])+abs(calibrate$lat-validate$lat[632])==min(abs(calibrate$lon-validate$lon[632])+abs(calibrate$lat-validate$lat[632])))
  
  paste(validate$lon[632],validate$lat[632])
  if (monthly_out) {
    pdf(paste(figpath,'/example_timeseries_alaska_',v,'_mon.pdf',sep=''), 
        width=4.5, height=6, paper='special')
  }else if(yearly_out){
    pdf(paste(figpath,'/example_timeseries_alaska_',v,'_yrly.pdf',sep=''), 
        width=7, height=6, paper='special')
  } else {
    pdf(paste(figpath,'/example_timeseries_alaska_',v,'.pdf',sep=''), 
        width=7, height=6, paper='special')
  }
  
  par(oma=c(0,0,0,0),mar=c(2,4,2,2),mfrow=c(3,1))
  #50yr summer temp
  plot(echam$time[pos2],echam$ensmean[632,pos2],ylim=c(-3,9),ty='l',col="black",main="Temperature",
       xlab='',ylab='[ºC]',xaxt='n',axes=F)
  
  polygon(c(echam$time[pos2],rev(echam$time[pos2])),c(apply(echam$data[632,pos2,],1,max),
                                                      rev(apply(echam$data[632,pos2,],1,min))),density=NA, col=rgb(1,1,1,3,maxColorValue=10))
  lines(analysis$time[pos2],analysis$ensmean[632,pos2],col="red")
  polygon(c(analysis$time[pos2],rev(analysis$time[pos2])),c(apply(analysis$data[632,pos2,],1,max),
                                                            rev(apply(analysis$data[632,pos2,],1,min))),density=NA, col=rgb(10,0,0,3,maxColorValue=10))
  
  axis(2,col="red")
  par(new=T)
  plot(validate$data[632,pos2],ylim=c(-1.5,10.5),ty='l',col="blue",main="",
       xlab='',ylab='[ºC]',xaxt='n',axes=F)
  axis(4, col="blue")
  box()
  #precip
  
  plot(echam$time[pos2],echam$ensmean[4608+632,pos2],ylim=c(10,130),ty='l',col="black",main="Precipitation",
       xlab='',ylab='[mm]',xaxt='n',axes=F)
  
  polygon(c(echam$time[pos2],rev(echam$time[pos2])),c(apply(echam$data[4608+632,pos2,],1,max),
                                                      rev(apply(echam$data[4608+632,pos2,],1,min))),density=NA, col=rgb(1,1,1,3,maxColorValue=10))
  lines(analysis$time[pos2],analysis$ensmean[4608+632,pos2],col="red")
  polygon(c(analysis$time[pos2],rev(analysis$time[pos2])),c(apply(analysis$data[4608+632,pos2,],1,max),
                                                            rev(apply(analysis$data[4608+632,pos2,],1,min))),density=NA, col=rgb(10,0,0,3,maxColorValue=10))
  
  axis(2,col="red")
  par(new=T)
  plot(validate$data[4608+632,pos2],ylim=c(-28.5,91.5),ty='l',col="blue",main="",
       xlab='',ylab='[mm]',xaxt='n',axes=F)
  axis(4, col="blue")
  box()
  
  #slp
  
  plot(echam$time[pos2],echam$ensmean[2*4608+632,pos2],ylim=c(1005,1017),ty='l',col="black",main="Sea level pressure",
       xlab='',ylab='[hPa]',xaxt='n',axes=F)
  
  polygon(c(echam$time[pos2],rev(echam$time[pos2])),c(apply(echam$data[2*4608+632,pos2,],1,max),
                                                      rev(apply(echam$data[2*4608+632,pos2,],1,min))),density=NA, col=rgb(1,1,1,3,maxColorValue=10))
  lines(analysis$time[pos2],analysis$ensmean[2*4608+632,pos2],col="red")
  polygon(c(analysis$time[pos2],rev(analysis$time[pos2])),c(apply(analysis$data[2*4608+632,pos2,],1,max),
                                                            rev(apply(analysis$data[2*4608+632,pos2,],1,min))),density=NA, col=rgb(10,0,0,3,maxColorValue=10))
  
  axis(2,col="red")
  par(new=T)
  plot(validate$time[pos2], validate$data[2*4608+632,pos2],ylim=c(1006,1018),ty='l',col="blue",main="",
       xlab='',ylab='[hPa]',xaxt='n',axes=F)
  axis(4, col="blue")
  axis(1)
  box()
  
  legend("bottomleft", c(valilegend, 'CCC400', "EFK400"),col=c("blue", "black", "red"),
         lty=c(1,1,1),lwd=c(1,1,1),pt.cex=1, pt.lwd=1,inset=0.005, bg='transparent',
         box.col='transparent', cex=1)
  dev.off()
  
  # Somewhere near Himalaya (Pakistan)
  paste(validate$lon[1460],validate$lat[1460])
  if (monthly_out) {
    pdf(paste(figpath,'/example_timeseries_pakistan_',v,'_mon.pdf',sep=''), 
        width=4.5, height=6, paper='special')
  }else if(yearly_out){
    pdf(paste(figpath,'/example_timeseries_pakistan_',v,'_yrly.pdf',sep=''), 
        width=7, height=6, paper='special')
  } else {
    pdf(paste(figpath,'/example_timeseries_pakistan_',v,'.pdf',sep=''), 
        width=7, height=6, paper='special')
  }
  par(oma=c(0,0,0,0),mar=c(2,4,2,2),mfrow=c(3,1))
  #50yr summer temp
  plot(echam$time[pos2],echam$ensmean[1460,pos2],ylim=c(29,36),ty='l',col="black",main="Temperature",
       xlab='',ylab='ºC',xaxt='n',axes=F)
  
  polygon(c(echam$time[pos2],rev(echam$time[pos2])),c(apply(echam$data[1460,pos2,],1,max),
                                                      rev(apply(echam$data[1460,pos2,],1,min))),density=NA, col=rgb(1,1,1,3,maxColorValue=10))
  lines(analysis$time[pos2],analysis$ensmean[1460,pos2],col="red")
  polygon(c(analysis$time[pos2],rev(analysis$time[pos2])),c(apply(analysis$data[1460,pos2,],1,max),
                                                            rev(apply(analysis$data[1460,pos2,],1,min))),density=NA, col=rgb(10,0,0,3,maxColorValue=10))
  
  axis(2,col="red")
  par(new=T)
  plot(validate$data[1460,pos2],ylim=c(25,32),ty='l',col="blue",main="",
       xlab='',ylab='[ºC]',xaxt='n',axes=F)
  axis(4, col="blue")
  box()
  #precip
  
  plot(echam$time[pos2],echam$ensmean[4608+1460,pos2],ylim=c(-20,80),ty='l',col="black",main="Precipitation",
       xlab='',ylab='[mm]',xaxt='n',axes=F)
  
  polygon(c(echam$time[pos2],rev(echam$time[pos2])),c(apply(echam$data[4608+1460,pos2,],1,max),
                                                      rev(apply(echam$data[4608+1460,pos2,],1,min))),density=NA, col=rgb(1,1,1,3,maxColorValue=10))
  lines(analysis$time[pos2],analysis$ensmean[4608+1460,pos2],col="red")
  polygon(c(analysis$time[pos2],rev(analysis$time[pos2])),c(apply(analysis$data[4608+1460,pos2,],1,max),
                                                            rev(apply(analysis$data[4608+1460,pos2,],1,min))),density=NA, col=rgb(10,0,0,3,maxColorValue=10))
  
  axis(2,col="red")
  par(new=T)
  plot(validate$data[4608+1460,pos2],ylim=c(15,115),ty='l',col="blue",main="",
       xlab='',ylab='[mm]',xaxt='n',axes=F)
  axis(4, col="blue")
  box()
  
  #slp
  
  plot(echam$time[pos2],echam$ensmean[2*4608+1460,pos2],ylim=c(997,1004),ty='l',col="black",main="Sea level pressure",
       xlab='',ylab='[hPa]',xaxt='n',axes=F)
  
  polygon(c(echam$time[pos2],rev(echam$time[pos2])),c(apply(echam$data[2*4608+1460,pos2,],1,max),
                                                      rev(apply(echam$data[2*4608+1460,pos2,],1,min))),density=NA, col=rgb(1,1,1,3,maxColorValue=10))
  lines(analysis$time[pos2],analysis$ensmean[2*4608+1460,pos2],col="red")
  polygon(c(analysis$time[pos2],rev(analysis$time[pos2])),c(apply(analysis$data[2*4608+1460,pos2,],1,max),
                                                            rev(apply(analysis$data[2*4608+1460,pos2,],1,min))),density=NA, col=rgb(10,0,0,3,maxColorValue=10))
  
  axis(2,col="red")
  par(new=T)
  plot(validate$time[pos2], validate$data[2*4608+1460,pos2],ylim=c(998.5,1005.5),ty='l',col="blue",main="",
       xlab='',ylab='[hPa]',xaxt='n',axes=F)
  axis(4, col="blue")
  axis(1)
  box()
  
  legend("bottomleft", c(valilegend, 'CCC400', "EFK400"),col=c("blue", "black", "red"),
         lty=c(1,1,1),lwd=c(1,1,1),pt.cex=1, pt.lwd=1,inset=0.005, bg='transparent',
         box.col='transparent', cex=1)
  dev.off()
  
  
  ########################### anomaly values##############################
  
  # Siberia
  paste(validate.anom$lon[592],validate.anom$lat[592])
  if (monthly_out) {
    pdf(paste(figpath,'/example_timeseries_sibiria_anom_',v,'_mon.pdf',sep=''), 
        width=4.5, height=6, paper='special')
  }else if(yearly_out){
    pdf(paste(figpath,'/example_timeseries_sibiria_anom_',v,'_yrly.pdf',sep=''), 
        width=7, height=6, paper='special')
  } else {
    pdf(paste(figpath,'/example_timeseries_sibiria_anom_',v,'.pdf',sep=''), 
        width=7, height=6, paper='special')
  }
  par(oma=c(0,0,0,0),mar=c(2,4,2,2),mfrow=c(3,1))
  #50yr summer temp
  plot(validate.anom$time[pos2],validate.anom$ensmean[592,pos2],ty='l',col="blue",main="Temperature",
       xlab='',ylab='[ºC]',xaxt='n',ylim=c(-5,5))
  lines(echam.anom$time[pos2],echam.anom$ensmean[592,pos2],col="black")
  polygon(c(echam.anom$time[pos2],rev(echam.anom$time[pos2])),c(apply(echam.anom$data[592,pos2,],1,max),
                                                                rev(apply(echam.anom$data[592,pos2,],1,min))),density=NA, col=rgb(1,1,1,3,maxColorValue=10))
  lines(analysis.anom$time[pos2],analysis.anom$ensmean[592,pos2],col="red")
  polygon(c(analysis.anom$time[pos2],rev(analysis.anom$time[pos2])),c(apply(analysis.anom$data[592,pos2,],1,max),
                                                                      rev(apply(analysis.anom$data[592,pos2,],1,min))),density=NA, col=rgb(10,0,0,3,maxColorValue=10))
  
  #precip
  
  plot(validate.anom$time[pos2],validate.anom$ensmean[4608+592,pos2],ty='l',col="blue",main="Precipitation",
       xlab='',ylab='[mm]',xaxt='n',ylim=c(-50,50))
  lines(echam.anom$time[pos2],echam.anom$ensmean[4608+592,pos2],col="black")
  polygon(c(echam.anom$time[pos2],rev(echam.anom$time[pos2])),c(apply(echam.anom$data[4608+592,pos2,],1,max),
                                                                rev(apply(echam.anom$data[4608+592,pos2,],1,min))),density=NA, col=rgb(1,1,1,3,maxColorValue=10))
  lines(analysis.anom$time[pos2],analysis.anom$ensmean[4608+592,pos2],col="red")
  polygon(c(analysis.anom$time[pos2],rev(analysis.anom$time[pos2])),c(apply(analysis.anom$data[4608+592,pos2,],1,max),
                                                                      rev(apply(analysis.anom$data[4608+592,pos2,],1,min))),density=NA, col=rgb(10,0,0,3,maxColorValue=10))
  
  #slp
  
  plot(validate.anom$time[pos2],validate.anom$ensmean[2*4608+592,pos2],ty='l',col="blue",main="Sea level pressure",
       xlab='',ylab='[hPa]',ylim=c(-6,6))
  lines(echam.anom$time[pos2],echam.anom$ensmean[2*4608+592,pos2],col="black")
  polygon(c(echam.anom$time[pos2],rev(echam.anom$time[pos2])),c(apply(echam.anom$data[2*4608+592,pos2,],1,max),
                                                                rev(apply(echam.anom$data[2*4608+592,pos2,],1,min))),density=NA, col=rgb(1,1,1,3,maxColorValue=10))
  lines(analysis.anom$time[pos2],analysis.anom$ensmean[2*4608+592,pos2],col="red")
  polygon(c(analysis.anom$time[pos2],rev(analysis.anom$time[pos2])),c(apply(analysis.anom$data[2*4608+592,pos2,],1,max),
                                                                      rev(apply(analysis.anom$data[2*4608+592,pos2,],1,min))),density=NA, col=rgb(10,0,0,3,maxColorValue=10))
  
  legend("bottomleft", c(valilegend, 'CCC400', "EFK400"),col=c("blue", "black", "red"),
         lty=c(1,1,1),lwd=c(1,1,1),pt.cex=1, pt.lwd=1,inset=0.005, bg='transparent',
         box.col='transparent', cex=1)
  dev.off()
  
  # Norway
  paste(validate.anom$lon[676],validate.anom$lat[676])
  if (monthly_out) {
    pdf(paste(figpath,'/example_timeseries_norway_anom_',v,'_mon.pdf',sep=''), 
        width=4.5, height=6, paper='special')
  }else if(yearly_out){
    pdf(paste(figpath,'/example_timeseries_norway_anom_',v,'_yrly.pdf',sep=''), 
        width=7, height=6, paper='special')
  } else {
    pdf(paste(figpath,'/example_timeseries_norway_anom_',v,'.pdf',sep=''), 
        width=7, height=6, paper='special')
  }
  par(oma=c(0,0,0,0),mar=c(2,4,2,2),mfrow=c(3,1))
  #50yr summer temp
  plot(validate.anom$time[pos2],validate.anom$ensmean[676,pos2],ty='l',col="blue",main="Temperature",
       xlab='',ylab='[ºC]',xaxt='n',ylim=c(-4,4))
  lines(echam.anom$time[pos2],echam.anom$ensmean[676,pos2],col="black")
  polygon(c(echam.anom$time[pos2],rev(echam.anom$time[pos2])),c(apply(echam.anom$data[676,pos2,],1,max),
                                                                rev(apply(echam.anom$data[676,pos2,],1,min))),density=NA, col=rgb(1,1,1,3,maxColorValue=10))
  lines(analysis.anom$time[pos2],analysis.anom$ensmean[676,pos2],col="red")
  polygon(c(analysis.anom$time[pos2],rev(analysis.anom$time[pos2])),c(apply(analysis.anom$data[676,pos2,],1,max),
                                                                      rev(apply(analysis.anom$data[676,pos2,],1,min))),density=NA, col=rgb(10,0,0,3,maxColorValue=10))
  
  #precip
  
  plot(validate.anom$time[pos2],validate.anom$ensmean[4608+676,pos2],ty='l',col="blue",main="Precipitation",
       xlab='',ylab='[mm]',xaxt='n',ylim=c(-50,50))
  lines(echam.anom$time[pos2],echam.anom$ensmean[4608+676,pos2],col="black")
  polygon(c(echam.anom$time[pos2],rev(echam.anom$time[pos2])),c(apply(echam.anom$data[4608+676,pos2,],1,max),
                                                                rev(apply(echam.anom$data[4608+676,pos2,],1,min))),density=NA, col=rgb(1,1,1,3,maxColorValue=10))
  lines(analysis.anom$time[pos2],analysis.anom$ensmean[4608+676,pos2],col="red")
  polygon(c(analysis.anom$time[pos2],rev(analysis.anom$time[pos2])),c(apply(analysis.anom$data[4608+676,pos2,],1,max),
                                                                      rev(apply(analysis.anom$data[4608+676,pos2,],1,min))),density=NA, col=rgb(10,0,0,3,maxColorValue=10))
  
  #slp
  
  plot(validate.anom$time[pos2],validate.anom$ensmean[2*4608+676,pos2],ty='l',col="blue",main="Sea level pressure",
       xlab='',ylab='[hPa]',ylim=c(-7,7))
  lines(echam.anom$time[pos2],echam.anom$ensmean[2*4608+676,pos2],col="black")
  polygon(c(echam.anom$time[pos2],rev(echam.anom$time[pos2])),c(apply(echam.anom$data[2*4608+676,pos2,],1,max),
                                                                rev(apply(echam.anom$data[2*4608+676,pos2,],1,min))),density=NA, col=rgb(1,1,1,3,maxColorValue=10))
  lines(analysis.anom$time[pos2],analysis.anom$ensmean[2*4608+676,pos2],col="red")
  polygon(c(analysis.anom$time[pos2],rev(analysis.anom$time[pos2])),c(apply(analysis.anom$data[2*4608+676,pos2,],1,max),
                                                                      rev(apply(analysis.anom$data[2*4608+676,pos2,],1,min))),density=NA, col=rgb(10,0,0,3,maxColorValue=10))
  
  legend("bottomleft", c(valilegend, 'CCC400', "EFK400"),col=c("blue", "black", "red"),
         lty=c(1,1,1),lwd=c(1,1,1),pt.cex=1, pt.lwd=1,inset=0.005, bg='transparent',
         box.col='transparent', cex=1)
  dev.off()
  
  
  # Greenland
  paste(validate.anom$lon[278],validate.anom$lat[278])
  if (monthly_out) {
    pdf(paste(figpath,'/example_timeseries_greenland_anom_',v,'_mon.pdf',sep=''), 
        width=4.5, height=6, paper='special')
  }else if(yearly_out){
    pdf(paste(figpath,'/example_timeseries_greenland_anom_',v,'_yrly.pdf',sep=''), 
        width=7, height=6, paper='special')
  } else {
    pdf(paste(figpath,'/example_timeseries_greenland_anom_',v,'.pdf',sep=''), 
        width=7, height=6, paper='special')
  }
  par(oma=c(0,0,0,0),mar=c(2,4,2,2),mfrow=c(3,1))
  #50yr summer temp
  plot(validate.anom$time[pos2],validate.anom$ensmean[278,pos2],ty='l',col="blue",main="Temperature",
       xlab='',ylab='[ºC]',xaxt='n',ylim=c(-4,4))
  lines(echam.anom$time[pos2],echam.anom$ensmean[278,pos2],col="black")
  polygon(c(echam.anom$time[pos2],rev(echam.anom$time[pos2])),c(apply(echam.anom$data[278,pos2,],1,max),
                                                                rev(apply(echam.anom$data[278,pos2,],1,min))),density=NA, col=rgb(1,1,1,3,maxColorValue=10))
  lines(analysis.anom$time[pos2],analysis.anom$ensmean[278,pos2],col="red")
  polygon(c(analysis.anom$time[pos2],rev(analysis.anom$time[pos2])),c(apply(analysis.anom$data[278,pos2,],1,max),
                                                                      rev(apply(analysis.anom$data[278,pos2,],1,min))),density=NA, col=rgb(10,0,0,3,maxColorValue=10))
  
  #precip
  
  plot(validate.anom$time[pos2],validate.anom$ensmean[4608+278,pos2],ty='l',col="blue",main="Precipitation",
       xlab='',ylab='[mm]',xaxt='n',ylim=c(-20,20))
  lines(echam.anom$time[pos2],echam.anom$ensmean[4608+278,pos2],col="black")
  polygon(c(echam.anom$time[pos2],rev(echam.anom$time[pos2])),c(apply(echam.anom$data[4608+278,pos2,],1,max),
                                                                rev(apply(echam.anom$data[4608+278,pos2,],1,min))),density=NA, col=rgb(1,1,1,3,maxColorValue=10))
  lines(analysis.anom$time[pos2],analysis.anom$ensmean[4608+278,pos2],col="red")
  polygon(c(analysis.anom$time[pos2],rev(analysis.anom$time[pos2])),c(apply(analysis.anom$data[4608+278,pos2,],1,max),
                                                                      rev(apply(analysis.anom$data[4608+278,pos2,],1,min))),density=NA, col=rgb(10,0,0,3,maxColorValue=10))
  
  #slp
  
  plot(validate.anom$time[pos2],validate.anom$ensmean[2*4608+278,pos2],ty='l',col="blue",main="Sea level pressure",
       xlab='',ylab='[hPa]',ylim=c(-6,6))
  lines(echam.anom$time[pos2],echam.anom$ensmean[2*4608+278,pos2],col="black")
  polygon(c(echam.anom$time[pos2],rev(echam.anom$time[pos2])),c(apply(echam.anom$data[2*4608+278,pos2,],1,max),
                                                                rev(apply(echam.anom$data[2*4608+278,pos2,],1,min))),density=NA, col=rgb(1,1,1,3,maxColorValue=10))
  lines(analysis.anom$time[pos2],analysis.anom$ensmean[2*4608+278,pos2],col="red")
  polygon(c(analysis.anom$time[pos2],rev(analysis.anom$time[pos2])),c(apply(analysis.anom$data[2*4608+278,pos2,],1,max),
                                                                      rev(apply(analysis.anom$data[2*4608+278,pos2,],1,min))),density=NA, col=rgb(10,0,0,3,maxColorValue=10))
  
  legend("bottomleft", c(valilegend, 'CCC400', "EFK400"),col=c("blue", "black", "red"),
         lty=c(1,1,1),lwd=c(1,1,1),pt.cex=1, pt.lwd=1,inset=0.005, bg='transparent',
         box.col='transparent', cex=1)
  dev.off()
  
  
  # Alaska
  # line beneath was used to determine closest lon/lat to validate.anom lat/lon 632
  #calipos<-which(abs(calibrate$lon-validate.anom$lon[632])+abs(calibrate$lat-validate.anom$lat[632])==min(abs(calibrate$lon-validate.anom$lon[632])+abs(calibrate$lat-validate.anom$lat[632])))
  
  paste(validate.anom$lon[632],validate.anom$lat[632])
  if (monthly_out) {
    pdf(paste(figpath,'/example_timeseries_alaska_anom_',v,'_mon.pdf',sep=''), 
        width=4.5, height=6, paper='special')
  }else if(yearly_out){
    pdf(paste(figpath,'/example_timeseries_alaska_anom_',v,'_yrly.pdf',sep=''), 
        width=7, height=6, paper='special')
  } else {
    pdf(paste(figpath,'/example_timeseries_alaska_anom_',v,'.pdf',sep=''), 
        width=7, height=6, paper='special')
  }
  
  par(oma=c(0,0,0,0),mar=c(2,4,2,2),mfrow=c(3,1))
  #50yr summer temp
  plot(validate.anom$time[pos2],validate.anom$ensmean[632,pos2],ty='l',col="blue",main="Temperature",
       xlab='',ylab='[ºC]',xaxt='n',ylim=c(-5,5))
  lines(echam.anom$time[pos2],echam.anom$ensmean[632,pos2],col="black")
  polygon(c(echam.anom$time[pos2],rev(echam.anom$time[pos2])),c(apply(echam.anom$data[632,pos2,],1,max),
                                                                rev(apply(echam.anom$data[632,pos2,],1,min))),density=NA, col=rgb(1,1,1,3,maxColorValue=10))
  lines(analysis.anom$time[pos2],analysis.anom$ensmean[632,pos2],col="red")
  polygon(c(analysis.anom$time[pos2],rev(analysis.anom$time[pos2])),c(apply(analysis.anom$data[632,pos2,],1,max),
                                                                      rev(apply(analysis.anom$data[632,pos2,],1,min))),density=NA, col=rgb(10,0,0,3,maxColorValue=10))
  
  #precip
  
  plot(validate.anom$time[pos2],validate.anom$ensmean[4608+632,pos2],ty='l',col="blue",main="Precipitation",
       xlab='',ylab='[mm]',xaxt='n',ylim=c(-50,50))
  lines(echam.anom$time[pos2],echam.anom$ensmean[4608+632,pos2],col="black")
  polygon(c(echam.anom$time[pos2],rev(echam.anom$time[pos2])),c(apply(echam.anom$data[4608+632,pos2,],1,max),
                                                                rev(apply(echam.anom$data[4608+632,pos2,],1,min))),density=NA, col=rgb(1,1,1,3,maxColorValue=10))
  lines(analysis.anom$time[pos2],analysis.anom$ensmean[4608+632,pos2],col="red")
  polygon(c(analysis.anom$time[pos2],rev(analysis.anom$time[pos2])),c(apply(analysis.anom$data[4608+632,pos2,],1,max),
                                                                      rev(apply(analysis.anom$data[4608+632,pos2,],1,min))),density=NA, col=rgb(10,0,0,3,maxColorValue=10))
  
  #slp
  
  plot(validate.anom$time[pos2],validate.anom$ensmean[2*4608+632,pos2],ty='l',col="blue",main="Sea level pressure",
       xlab='',ylab='[hPa]',ylim=c(-6,6))
  lines(echam.anom$time[pos2],echam.anom$ensmean[2*4608+632,pos2],col="black")
  polygon(c(echam.anom$time[pos2],rev(echam.anom$time[pos2])),c(apply(echam.anom$data[2*4608+632,pos2,],1,max),
                                                                rev(apply(echam.anom$data[2*4608+632,pos2,],1,min))),density=NA, col=rgb(1,1,1,3,maxColorValue=10))
  lines(analysis.anom$time[pos2],analysis.anom$ensmean[2*4608+632,pos2],col="red")
  polygon(c(analysis.anom$time[pos2],rev(analysis.anom$time[pos2])),c(apply(analysis.anom$data[2*4608+632,pos2,],1,max),
                                                                      rev(apply(analysis.anom$data[2*4608+632,pos2,],1,min))),density=NA, col=rgb(10,0,0,3,maxColorValue=10))
  
  legend("bottomleft", c(valilegend, 'CCC400', "EFK400"),col=c("blue", "black", "red"),
         lty=c(1,1,1),lwd=c(1,1,1),pt.cex=1, pt.lwd=1,inset=0.005, bg='transparent',
         box.col='transparent', cex=1)
  dev.off()
  
  # Somewhere near Himalaya (Pakistan)
  paste(validate.anom$lon[1460],validate.anom$lat[1460])
  if (monthly_out) {
    pdf(paste(figpath,'/example_timeseries_pakistan_anom_',v,'_mon.pdf',sep=''), 
        width=4.5, height=6, paper='special')
  }else if(yearly_out){
    pdf(paste(figpath,'/example_timeseries_pakistan_anom_',v,'_yrly.pdf',sep=''), 
        width=7, height=6, paper='special')
  } else {
    pdf(paste(figpath,'/example_timeseries_pakistan_anom_',v,'.pdf',sep=''), 
        width=7, height=6, paper='special')
  }
  par(oma=c(0,0,0,0),mar=c(2,4,2,2),mfrow=c(3,1))
  #50yr summer temp
  plot(validate.anom$time[pos2],validate.anom$ensmean[1460,pos2],ty='l',col="blue",main="Temperature",
       xlab='',ylab='[ºC]',xaxt='n',ylim=c(-4,4))
  lines(echam.anom$time[pos2],echam.anom$ensmean[1460,pos2],col="black")
  polygon(c(echam.anom$time[pos2],rev(echam.anom$time[pos2])),c(apply(echam.anom$data[1460,pos2,],1,max),
                                                                rev(apply(echam.anom$data[1460,pos2,],1,min))),density=NA, col=rgb(1,1,1,3,maxColorValue=10))
  lines(analysis.anom$time[pos2],analysis.anom$ensmean[1460,pos2],col="red")
  polygon(c(analysis.anom$time[pos2],rev(analysis.anom$time[pos2])),c(apply(analysis.anom$data[1460,pos2,],1,max),
                                                                      rev(apply(analysis.anom$data[1460,pos2,],1,min))),density=NA, col=rgb(10,0,0,3,maxColorValue=10))
  
  #precip
  
  plot(validate.anom$time[pos2],validate.anom$ensmean[4608+1460,pos2],ty='l',col="blue",main="Precipitation",
       xlab='',ylab='[mm]',xaxt='n',ylim=c(-40,60))
  lines(echam.anom$time[pos2],echam.anom$ensmean[4608+1460,pos2],col="black")
  polygon(c(echam.anom$time[pos2],rev(echam.anom$time[pos2])),c(apply(echam.anom$data[4608+1460,pos2,],1,max),
                                                                rev(apply(echam.anom$data[4608+1460,pos2,],1,min))),density=NA, col=rgb(1,1,1,3,maxColorValue=10))
  lines(analysis.anom$time[pos2],analysis.anom$ensmean[4608+1460,pos2],col="red")
  polygon(c(analysis.anom$time[pos2],rev(analysis.anom$time[pos2])),c(apply(analysis.anom$data[4608+1460,pos2,],1,max),
                                                                      rev(apply(analysis.anom$data[4608+1460,pos2,],1,min))),density=NA, col=rgb(10,0,0,3,maxColorValue=10))
  
  #slp
  
  plot(validate.anom$time[pos2],validate.anom$ensmean[2*4608+1460,pos2],ty='l',col="blue",main="Sea level pressure",
       xlab='',ylab='[hPa]',ylim=c(-4,4))
  lines(echam.anom$time[pos2],echam.anom$ensmean[2*4608+1460,pos2],col="black")
  polygon(c(echam.anom$time[pos2],rev(echam.anom$time[pos2])),c(apply(echam.anom$data[2*4608+1460,pos2,],1,max),
                                                                rev(apply(echam.anom$data[2*4608+1460,pos2,],1,min))),density=NA, col=rgb(1,1,1,3,maxColorValue=10))
  lines(analysis.anom$time[pos2],analysis.anom$ensmean[2*4608+1460,pos2],col="red")
  polygon(c(analysis.anom$time[pos2],rev(analysis.anom$time[pos2])),c(apply(analysis.anom$data[2*4608+1460,pos2,],1,max),
                                                                      rev(apply(analysis.anom$data[2*4608+1460,pos2,],1,min))),density=NA, col=rgb(10,0,0,3,maxColorValue=10))
  
  legend("bottomleft", c(valilegend, 'CCC400', "EFK400"),col=c("blue", "black", "red"),
         lty=c(1,1,1),lwd=c(1,1,1),pt.cex=1, pt.lwd=1,inset=0.005, bg='transparent',
         box.col='transparent', cex=1)
  dev.off()
  
  
  
  
  ########################################################################
  ############timeseries plots for indices of the whole period############
  ########################################################################
  inds <- c('ENH.temp2', 'EU.temp2', 'NEU.temp2', 'GLO.temp2')
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
  vind.tot$names<-indices
  for(ind in inds){
    pdf(paste(figpath,'/timeseries_',v,'_',ind,'.pdf',sep=''), width=15, height=4.5, paper='special')
    
    winter<-seq(1,length(eind.tot$time),by=2)
    summer<-seq(2,length(eind.tot$time), by=2)
    
    valiwinter<-seq(1,length(vind.tot$time),by=2)
    valisummer<-seq(2,length(vind.tot$time),by=2)
    
    ylimmin<-min(cbind(c(vind.tot$ensmean[which(vind.tot$names==ind),valiwinter],rep(NA,syr-syrtot)), eind.tot$ensmean[which(eind.tot$names==ind),winter],
                       aind.tot$ensmean[which(aind.tot$names==ind),winter]),na.rm=T)-1
    
    ylimmax<-max(cbind(c(vind.tot$ensmean[which(vind.tot$names==ind),valiwinter],rep(NA,syr-syrtot)), eind.tot$ensmean[which(eind.tot$names==ind),winter],
                       aind.tot$ensmean[which(aind.tot$names==ind),winter]),na.rm=T)+1
    
    par(oma=c(0,0,0,0),mar=c(4,4,2,0.5),mfrow=c(1,2))
    plot(eind.tot$time[winter],eind.tot$ensmean[which(eind.tot$names==ind),winter], ylim=c(ylimmin,ylimmax),ty='l',col="black",main="Oct.-Mar.",
         xlab='year',ylab=paste(ind,'[°C]'))
    polygon(c(eind.tot$time[winter],rev(eind.tot$time[winter])),c(apply(eind.tot$data[which(eind.tot$names==ind),winter,],1,max),
                                                                  rev(apply(eind.tot$data[which(eind.tot$names==ind),winter,],1,min))),density=NA, col=rgb(1,1,1,3,maxColorValue=10))
    
    
    lines(aind.tot$time[winter],aind.tot$ensmean[which(aind.tot$names==ind),winter],col="red")
    polygon(c(aind.tot$time[winter],rev(aind.tot$time[winter])),c(apply(aind.tot$data[which(aind.tot$names==ind),winter,],1,max),
                                                                  rev(apply(aind.tot$data[which(aind.tot$names==ind),winter,],1,min))),density=NA, col=rgb(10,0,0,3,maxColorValue=10))
    lines(vind.tot$time[valiwinter], vind.tot$ensmean[which(vind.tot$names==ind),valiwinter],col="blue")
    legend("bottomleft", c(valilegend, 'CCC400', "EFK400"),col=c("blue", "black", "red"),
           lty=c(1,1,1),lwd=c(1,1,1),pt.cex=1, pt.lwd=1,inset=0.005, bg='transparent',
           box.col='transparent', cex=1)
    
    ylimmin<-min(cbind(c(vind.tot$ensmean[which(vind.tot$names==ind),valisummer],rep(NA,syr-syrtot)), eind.tot$ensmean[which(eind.tot$names==ind),summer],
                       aind.tot$ensmean[which(aind.tot$names==ind),summer]),na.rm=T)-1
    
    ylimmax<-max(cbind(c(vind.tot$ensmean[which(vind.tot$names==ind),valisummer],rep(NA,syr-syrtot)), eind.tot$ensmean[which(eind.tot$names==ind),summer],
                       aind.tot$ensmean[which(aind.tot$names==ind),summer]),na.rm=T)+1
    
    plot(eind.tot$time[winter],eind.tot$ensmean[which(eind.tot$names==ind),summer], ylim=c(ylimmin,ylimmax),ty='l',col="black",main="Apr.-Sept.",
         xlab='year',ylab='')
    polygon(c(eind.tot$time[summer],rev(eind.tot$time[winter])),c(apply(eind.tot$data[which(eind.tot$names==ind),summer,],1,max),
                                                                  rev(apply(eind.tot$data[which(eind.tot$names==ind),summer,],1,min))),density=NA, col=rgb(1,1,1,3,maxColorValue=10))
    
    
    lines(aind.tot$time[winter],aind.tot$ensmean[which(aind.tot$names==ind),summer],col="red")
    polygon(c(aind.tot$time[winter],rev(aind.tot$time[winter])),c(apply(aind.tot$data[which(aind.tot$names==ind),summer,],1,max),
                                                                  rev(apply(aind.tot$data[which(aind.tot$names==ind),summer,],1,min))),density=NA, col=rgb(10,0,0,3,maxColorValue=10))
    lines(vind.tot$time[valiwinter], vind.tot$ensmean[which(vind.tot$names==ind),valisummer],col="blue")
    dev.off()
    
    
    # indices timeserie whole periode smoothed by runningmean wdow = window of running mean
    wdow<-11
    vind.run<-runmean(vind.tot$ensmean[which(vind.tot$names==ind),], 2*wdow,endrule="NA")
    aind.run<-runmean(aind.tot$ensmean[which(aind.tot$names==ind),], 2*wdow,endrule="NA")
    eind.run<-runmean(eind.tot$ensmean[which(eind.tot$names==ind),], 2*wdow,endrule="NA")
    
    ylimmin<-min(rbind(vind.run, eind.run,
                       aind.run),na.rm=T)-1
    
    ylimmax<-max(rbind(vind.run, eind.run,
                       aind.run),na.rm=T)+1
    
    pdf(paste(figpath,'/timeseries_',v,'_',ind,'_smoothed.pdf',sep=''), width=15, height=4.5, paper='special')
    plot(eind.tot$time,eind.run, ylim=c(ylimmin,ylimmax),ty='l',col="black",
         xlab='year',ylab=paste(ind,'[°C]'), main=paste(wdow,'year running mean'))
    lines(vind.tot$time, vind.run, col="blue")
    lines(aind.tot$time, aind.run, col="red")
    legend("bottomleft", c(valilegend, 'CCC400', "EFK400"),col=c("blue", "black", "red"),
           lty=c(1,1,1),lwd=c(1,1,1),pt.cex=1, pt.lwd=1,inset=0.005, bg='transparent',
           box.col='transparent', cex=1)
    dev.off()
  }
  
  #chose timestep for sample year plots
  t=2
  if (!monthly_out&!yearly_out) {
    #for (t in 1:60) {
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
    #}
    
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
  
  ###############################################################################
  # Fig. xx: spread-error ratio analysis
  if (!monthly_out){
    plotdata=echam
    #plotdata$ensmean <- cbind(ech_spr_err_ratio,ana_spr_err_ratio)
    #plotdata$ensmean <- cbind(sprerr.win,sprerr.sum)
    plotdata$data <- array(cbind(sprerr.win,sprerr.sum,sprerr.c.win,sprerr.c.sum), 
                           c(length(sprerr.c.win), 1, 4))
    plotdata$ensmean <- array(cbind(sprerr.win,sprerr.sum,sprerr.c.win,sprerr.c.sum), 
                              c(length(sprerr.c.win), 4))
    plotdata$names <- plotdata$names[echam$names=="temp2"]
    plotdata$lon <- plotdata$lon[echam$names=="temp2"]
    plotdata$lat <- plotdata$lat[echam$names=="temp2"]
    # v3: corrected for uncertainties in instrumental data
    # v4: only 1901-1980 because of validation data error after 1980
    # v5: obs error not taken into account (top), taken in account (bottom)
    
    # pdf(paste0(figpath,'/spread_error_ratio_anom_tmean.pdf'), width=9, height=7, paper='special')
    # layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE), height=c(3,3,1))
    # #layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
    # par(oma=c(0,2,3,0))
    # #levs <- c(0,0.33,0.5,0.57,0.67,0.77,0.83,0.91,1.1,1.2,1.3,1.5,1.75,2,3,Inf)
    levs <- c(0,0.33,0.5,0.67,0.83,1.2,1.5,2,3,Inf)
    # plot_echam(plotdata, varname='temp2', type='data', cex.pt=1.5, names=pnames[1:dim(plotdata$data)[3]], 
    #            lev=levs, st.col=NULL, stations=NULL, add=T, #ti=(1:2),
    #            colnames=c("Oct.-Apr.","May-Sep."),rownames=c("w/o obs. error","w/ obs. error"))
    # dev.off()
    plot_echam4(plotdata, varname='temp2', cex.pt=1.5, names=pnames[1:dim(plotdata$data)[3]], lev=levs , 
                type='ensmean', st.col=NULL, stations=plstat,NHseason=NULL,
                plotname=paste0('spread_error_ratio_anom_tmean_analysis.pdf'), paper='special')
    
    
    # Fig. xx: spread-error ratio echam
    plotdata=echam
    plotdata$data <- array(cbind(ech.sprerr.win,ech.sprerr.sum,ech.sprerr.c.win,ech.sprerr.c.sum), 
                           c(length(ech.sprerr.c.win), 1, 4))
    plotdata$ensmean <- array(cbind(ech.sprerr.win,ech.sprerr.sum,ech.sprerr.c.win,ech.sprerr.c.sum), 
                              c(length(ech.sprerr.c.win), 4))
    plotdata$names <- plotdata$names[echam$names=="temp2"]
    plotdata$lon <- plotdata$lon[echam$names=="temp2"]
    plotdata$lat <- plotdata$lat[echam$names=="temp2"]
    # v3: corrected for uncertainties in instrumental data
    # v4: only 1901-1980 because of validation data error after 1980
    # v5: obs error not taken into account (top), taken in account (bottom)
    
    #layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
    
    #levs <- c(0,0.33,0.5,0.57,0.67,0.77,0.83,0.91,1.1,1.2,1.3,1.5,1.75,2,3,Inf)
    levs <- c(0,0.33,0.5,0.67,0.83,1.2,1.5,2,3,Inf)
    plot_echam4(plotdata, varname='temp2', cex.pt=1.5, names=pnames[1:dim(plotdata$data)[3]], lev=levs , 
                type='ensmean', st.col=NULL, stations=plstat,NHseason=NULL,
                plotname=paste0('spread_error_ratio_anom_tmean_echam.pdf'), paper='special')
    
    
    
    # plotdata$data <- array(mes_obserr,c(nrow(mes_obserr),1,2))
    # pdf(paste(figpath,'/spread_error_ratio_obserr.pdf',sep=''), width=9, height=4.5, paper='special')
    # layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
    # par(oma=c(0,0,0,0))
    # levs <- c(-Inf, seq(0,2,0.25), Inf)
    # plot_echam(plotdata, varname='temp2', type='data', cex.pt=1.5, names=pnames[1:dim(plotdata$data)[3]], lev=levs, st.col=NULL, stations=NULL, add=T)
    # dev.off()
    
  }
  ################################################################################
  # Fig. xx: Talagrant diagram
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
  
  # if (!real_proxies){
  # ################################################################################
  # # Fig. 1a: maps of proxy locations and correlation with cru
  # pdf(paste(figpath,'corr_prox-cru_temp.pdf',sep='/'), width=9, height=4.5, paper='special')
  # levs <- c(-1,seq(-0.8,0.8,0.2),1)
  # #levs <- c(-1,seq(0.3,0.7,0.1),1)
  # #if (real_proxies) {levs <- c(-1,seq(-0.4,0.4,0.2),1)}
  # cor.col <- array(as.numeric(cut(prox.cru.corr.temp, levs)), dim(prox.cru.corr.temp))
  # cols <- rbfun(length(levs) - 1)
  # #cornames <- c('< 0.3', paste(seq(0.3,0.6,0.1), seq(0.4,0.7,0.1), sep='-'), '> 0.7')
  # cornames <- c('< -0.8', paste(c(-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6), c(-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8), sep='-'), '> 0.8')
  # #if (real_proxies) {cornames <- c('< -0.4', paste(seq(-0.4,0.2,0.2), seq(-0.2,0.4,0.2), sep='-'), '> 0.7')}
  #   
  # par(mar=rep(0.5,4))
  # map(interior=F)
  # map(add=T, region=c('Caspian', 'Great Lakes'))
  # points(calibrate$lon[calibrate$name=='temp2'], calibrate$lat[calibrate$name=='temp2'], col=cols[cor.col[,1]], bg=cols[cor.col[,2]], pch=21, lwd=2, cex=2)
  # legend('bottomleft', c(cornames, 'winter', 'summer'), pch=c(rep(19, length(cols)), 1, 19), col=c(cols, cols[1], cols[1]), pt.cex=1, pt.lwd=1, inset=0.005, bg='white', title='Correlation', box.col='white', cex=1.2)
  # box()
  # dev.off()
  # #\caption{Correlation of pseudoproxies with reference time series in winter (October to March, open circles) and summer (April to September, filled dots)}
  # 
  # 
  # ################################################################################
  # # Fig. 1aa: maps of proxy locations and correlation with cru
  # pdf(paste(figpath,'corr_prox-cru_precip.pdf',sep='/'), width=9, height=4.5, paper='special')
  # levs <- c(-1,seq(-0.8,0.8,0.2),1)
  # #levs <- c(-1,seq(0.3,0.7,0.1),1)
  # #if (real_proxies) {levs <- c(-1,seq(-0.4,0.4,0.2),1)}
  # cor.col <- array(as.numeric(cut(prox.cru.corr.precip, levs)), dim(prox.cru.corr.precip))
  # cols <- rbfun(length(levs) - 1)
  # #cornames <- c('< 0.3', paste(seq(0.3,0.6,0.1), seq(0.4,0.7,0.1), sep='-'), '> 0.7')
  # cornames <- c('< -0.8', paste(c(-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6), c(-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8), sep='-'), '> 0.8')
  # #if (real_proxies) {cornames <- c('< -0.4', paste(seq(-0.4,0.2,0.2), seq(-0.2,0.4,0.2), sep='-'), '> 0.7')}
  # 
  # par(mar=rep(0.5,4))
  # map(interior=F)
  # map(add=T, region=c('Caspian', 'Great Lakes'))
  # points(calibrate$lon[calibrate$name=='precip'], calibrate$lat[calibrate$name=='precip'], col=cols[cor.col[,1]], bg=cols[cor.col[,2]], pch=21, lwd=2, cex=2)
  # legend('bottomleft', c(cornames, 'winter', 'summer'), pch=c(rep(19, length(cols)), 1, 19), col=c(cols, cols[1], cols[1]), pt.cex=1, pt.lwd=1, inset=0.005, bg='white', title='Correlation', box.col='white', cex=1.2)
  # box()
  # dev.off()
  # #\caption{Correlation of pseudoproxies with reference time series in winter (October to March, open circles) and summer (April to September, filled dots)}
  # 
  # ################################################################################
  # # Fig. 1aaa: maps of proxy locations and correlation with cru
  # pdf(paste(figpath,'corr_prox-cru_slp.pdf',sep='/'), width=9, height=4.5, paper='special')
  # levs <- c(-1,seq(-0.8,0.8,0.2),1)
  # #levs <- c(-1,seq(0.3,0.7,0.1),1)
  # #if (real_proxies) {levs <- c(-1,seq(-0.4,0.4,0.2),1)}
  # cor.col <- array(as.numeric(cut(prox.cru.corr.slp, levs)), dim(prox.cru.corr.slp))
  # cols <- rbfun(length(levs) - 1)
  # #cornames <- c('< 0.3', paste(seq(0.3,0.6,0.1), seq(0.4,0.7,0.1), sep='-'), '> 0.7')
  # cornames <- c('< -0.8', paste(c(-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6), c(-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8), sep='-'), '> 0.8')
  # #if (real_proxies) {cornames <- c('< -0.4', paste(seq(-0.4,0.2,0.2), seq(-0.2,0.4,0.2), sep='-'), '> 0.7')}
  # 
  # par(mar=rep(0.5,4))
  # map(interior=F)
  # map(add=T, region=c('Caspian', 'Great Lakes'))
  # points(calibrate$lon[calibrate$name=='slp'], calibrate$lat[calibrate$name=='slp'], col=cols[cor.col[,1]], bg=cols[cor.col[,2]], pch=21, lwd=2, cex=2)
  # legend('bottomleft', c(cornames, 'winter', 'summer'), pch=c(rep(19, length(cols)), 1, 19), col=c(cols, cols[1], cols[1]), pt.cex=1, pt.lwd=1, inset=0.005, bg='white', title='Correlation', box.col='white', cex=1.2)
  # box()
  # dev.off()
  # #\caption{Correlation of pseudoproxies with reference time series in winter (October to March, open circles) and summer (April to September, filled dots)}
  # 
  # 
  # 
  # 
  # 
  # 
  # ################################################################################
  # # Fig. 1b: maps of proxy locations and bias wrt cru
  # pdf(paste(figpath,'bias_prox-cru.pdf',sep='/'), width=9, height=4.5, paper='special')
  # if (inst_pseudoproxy || pseudoproxy) {
  #   blevs <- c(-10, -3:3,10)
  #   bias.col <- array(as.numeric(cut(apply(array(dist.arr, c(length(prox.cru.corr), length(dist.arr)/length(prox.cru.corr))), 1, mean), blevs)), dim(prox.cru.corr))
  # } else {
  #   blevs <- c(-Inf, -3:3,Inf)
  #   bias.col <- array(as.numeric(cut(prox.cru.bias, blevs)), dim(prox.cru.bias))
  # }
  #   bcols <- rbfun(length(blevs) - 1)
  #   bcornames <- c('< -3', paste(seq(-3,2), seq(-2, 3), sep=' - '), '> 3')
  #   map(interior=F)
  #   map(add=T, region=c('Caspian', 'Great Lakes'))
  #   points(calibrate$lon, calibrate$lat, col=bcols[bias.col[,1]], bg=bcols[bias.col[,2]], pch=21, lwd=2, cex=1.4)
  #   legend('bottomleft', c(bcornames, 'winter', 'summer'), pch=c(rep(19, length(bcols)), 1, 19), col=c(bcols, bcols[1], bcols[1]), pt.cex=1, pt.lwd=1, inset=0.005, bg='white', title='Bias', box.col='white')
  #   box()
  # # }
  # dev.off()
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # ################################################################################
  # # Fig. 1c: maps of proxy locations and correlation with echam ensmean
  # pdf(paste(figpath,'corr_prox-echam.pdf',sep='/'), width=9, height=4.5, paper='special')
  # levs <- c(-1,seq(-0.6,0.6,0.2),1)
  # #levs <- c(-1,seq(0.3,0.7,0.1),1)
  # cor.col <- array(as.numeric(cut(prox.echam.corr, levs)), dim(prox.echam.corr))
  # cols <- rbfun(length(levs) - 1)
  # #cornames <- c('< 0.3', paste(seq(0.3,0.6,0.1), seq(0.4,0.7,0.1), sep='-'), '> 0.7')
  # cornames <- c('< -0.6', paste(c(-0.6,-0.4,-0.2,0.0,0.2,0.4), c(-0.4,-0.2,0.0,0.2,0.4), sep='-'), '> 0.6')
  # 
  # par(mar=rep(0.5,4))
  # map(interior=F)
  # map(add=T, region=c('Caspian', 'Great Lakes'))
  # points(calibrate$lon, calibrate$lat, col=cols[cor.col[,1]], bg=cols[cor.col[,2]], pch=21, lwd=2, cex=2)
  # legend('bottomleft', c(cornames, 'winter', 'summer'), pch=c(rep(19, length(cols)), 1, 19), col=c(cols, cols[1], cols[1]), pt.cex=1, pt.lwd=1, inset=0.005, bg='white', title='Correlation', box.col='white', cex=1.2)
  # box()
  # dev.off()
  # #\caption{Correlation of pseudoproxies with reference time series in winter (October to March, open circles) and summer (April to September, filled dots)}
  # 
  # 
  # 
  # ################################################################################
  # # Fig. 1d: maps of proxy locations and bias wrt echam ensmean
  # pdf(paste(figpath,'bias_prox-echam.pdf',sep='/'), width=9, height=4.5, paper='special')
  # if (inst_pseudoproxy || pseudoproxy) {
  #   blevs <- c(-10, -3:3,10)
  #   bias.col <- array(as.numeric(cut(apply(array(dist.arr, c(length(prox.echam.corr), length(dist.arr)/length(prox.echam.corr))), 1, mean), blevs)), dim(prox.echam.corr))
  # } else {
  #   blevs <- c(-Inf, -3:3,Inf)
  #   bias.col <- array(as.numeric(cut(prox.echam.bias, blevs)), dim(prox.echam.bias))
  # }
  #   bcols <- rbfun(length(blevs) - 1)
  #   bcornames <- c('< -3', paste(seq(-3,2), seq(-2, 3), sep=' - '), '> 3')
  #   map(interior=F)
  #   map(add=T, region=c('Caspian', 'Great Lakes'))
  #   points(calibrate$lon, calibrate$lat, col=bcols[bias.col[,1]], bg=bcols[bias.col[,2]], pch=21, lwd=2, cex=1.4)
  #   legend('bottomleft', c(bcornames, 'winter', 'summer'), pch=c(rep(19, length(bcols)), 1, 19), col=c(bcols, bcols[1], bcols[1]), pt.cex=1, pt.lwd=1, inset=0.005, bg='white', title='Bias', box.col='white')
  #   box()
  # # }
  # dev.off()
  # 
  # 
  # } # end if real_proxies
  
  
  
  ################################################################################
  # Fig. 2a: maps of temperature spread
  # spread calc now done in EnSRF_data.R
  #data.dim=c((694*3),2,29,30)
  #ech.spread <- apply(sqrt(apply(array(echam$data[1:(694*3)] - as.vector(echam$ensmean[1:(694*3)]), data.dim)**2, 1:3, mean,na.rm=T)), 1:2, mean, na.rm=T)
  #ana.spread <- lapply(analysis, function(x) apply(sqrt(apply(array(x$data - as.vector(x$ensmean), data.dim)**2, 1:3, mean, na.rm=T)), 1:2, mean, na.rm=T))
  if(!monthly_out){
    espread <- echam
    espread$data <- array(ech.spread, c(nrow(ech.spread), 1, ncol(ech.spread)))
    espread$ensmean <- array(ech.spread, c(nrow(ech.spread),ncol(ech.spread)))
    
    aspread <- echam
    if (pseudoproxy) { 
      ana.spread.bkp <- ana.spread
      ana.spread <- ana.spread.bkp[4] 
    }
    #if ((instrumental) || (inst_at_proxy) || (real_proxies)) {ana.spread <- ana.spread.bkp[2] }
    #aspread$data <- array(unlist(ana.spread)/as.vector(ech.spread)*100, c(nrow(ech.spread), 1, ncol(ech.spread)*(length(ana.spread))))
    aspread$data <- array(unlist(ana.spread)/as.vector(ech.spread)*100, c(nrow(ech.spread), 1, ncol(ech.spread)))
    aspread$ensmean <- array(unlist(ana.spread)/as.vector(ech.spread)*100, c(nrow(ech.spread), ncol(ech.spread)))
    
    ## Luca: couldn't convert this part into plot_echam4 yet
    
    
    #pdf('figures/inst/spread.pdf', width=9, height=6, paper='special')
    pdf(paste(figpath,'spread_temp.pdf',sep='/'), width=9, height=6, paper='special')
    oldpar <- par(no.readonly=TRUE)
    #layout(matrix(c(1:3,3,3+ seq(1,length(ana.spread)*2), rep(4,2) + length(ana.spread)*2),length(ana.spread)+3, 2, byrow=T), height=c(5,lcm(2), rep(5, length(ana.spread)),lcm(2)))
    #layout(matrix(c(1:3,3,3+ seq(1,2), rep(4,2) + 2),3, 2, byrow=T), height=c(5,lcm(2), rep(5, 1),lcm(2)))
    layout(matrix(c(1,2,3,3,4,5,6,6), 4, 2, byrow = TRUE), height=c(3,1,3,1))
    par(oma=c(0,0,0,0))
    plot_echam(espread, varname='temp2', names=pnames[1:dim(espread$data)[3]], symmetric=FALSE, cex.pt=1.3, add=TRUE)
    plot_echam(aspread, varname='temp2', names=pnames[dim(espread$data)[3] + 1:dim(aspread$data)[3]], lev=c(seq(0,90,10),101), cex.pt=1.3, add=TRUE, st.col=NULL, stations=plstat) #calibrate)
    #par(oldpar)
    #plot_echam(espread, varname='temp2', names=pnames[1:dim(espread$data)[3]], symmetric=FALSE, cex.pt=1.3, add=TRUE)
    #plot_echam(aspread, varname='temp2', names=pnames[dim(espread$data)[3] + 1:dim(aspread$data)[3]], lev=seq(0,100,10), cex.pt=1.3, add=TRUE, st.col=1, stations=calibrate) #calibrate)
    #par(oldpar)
    # \caption{Average temperature spread of the ECHAM ensemble in a) winter (October to March) and b) summer (April to September). Percentage of spread in the analysis ensembles with respect to the ECHAM ensemble for the EnSRF analysis with perfect proxies in c) and d), with perfect proxies and localization in e) and f), with pseudoproxies in g) and h), and with pseudoproxies and localization in i) and j). }
    dev.off()
  }
  ################################################################################
  # Fig. 2b: maps of precipitation spread
  
  ## Luca: couldn't convert this to plot_echam4 yet
  if (!monthly_out){
    
    pdf(paste(figpath,'spread_precip.pdf',sep='/'), width=9, height=6, paper='special')
    #<<label=spread_precip, echo=FALSE, fig=TRUE, width=8, height=11, results=hide, eps=FALSE>>=
    oldpar <- par(no.readonly=TRUE)
    layout(matrix(c(1,2,3,3,4,5,6,6), 4, 2, byrow = TRUE), height=c(3,1,3,1))
    par(oma=c(0,0,0,0))
    plot_echam(espread, varname='precip', names=pnames[1:dim(espread$data)[3]], symmetric=FALSE, cex.pt=1.3, add=TRUE, levs=c(0,1,2,3,5,8,13,22,36,60,100))
    plot_echam(aspread, varname='precip', names=pnames[dim(espread$data)[3] + 1:dim(aspread$data)[3]], lev=seq(0,100,10), cex.pt=1.3, add=TRUE, st.col=NULL, stations=plstat)
    #layout(matrix(c(1:3,3,3+ seq(1,length(ana.spread)*2), rep(4,2) + length(ana.spread)*2), length(ana.spread)+3, 2, byrow=T), height=c(5,lcm(2), rep(5, length(ana.spread)),lcm(2)))
    #par(oma=c(0,0,0,0))
    #plot_echam(espread, varname='precip', names=pnames[1:dim(espread$data)[3]], symmetric=FALSE, cex.pt=1.3, add=TRUE, levs=c(0,1,2,3,5,8,13,22,36,60,100))
    #plot_echam(aspread, varname='precip', cols=rbfun(10), names=pnames[dim(espread$data)[3] + 1:dim(aspread$data)[3]], levs=seq(0,100,10), cex.pt=1.3, add=TRUE, st.col=1, stations=plstat) #calibrate)
    #par(oldpar)
    dev.off()
    # \caption{According to figure \ref{fig:spread} but for precipitation in winter (left) and summer (right). Please note the quasi-logarithmic shading in panels a) and b).}
    
  }
  
  ################################################################################
  # Fig. 2c: maps of SLP spread
  
  ## Luca: couldn't convert this to plot_echam4 yet
  if (!monthly_out){
    pdf(paste(figpath,'spread_slp.pdf',sep='/'), width=9, height=6, paper='special')
    #<<label=spread_slp, echo=FALSE, fig=TRUE, width=8, height=11, results=hide, eps=FALSE>>=
    oldpar <- par(no.readonly=TRUE)
    layout(matrix(c(1,2,3,3,4,5,6,6), 4, 2, byrow = TRUE), height=c(3,1,3,1))
    par(oma=c(0,0,0,0))
    plot_echam(espread, varname='slp', names=pnames[1:dim(espread$data)[3]], symmetric=FALSE, cex.pt=1.3, add=TRUE)
    plot_echam(aspread, varname='slp', names=pnames[dim(espread$data)[3] + 1:dim(aspread$data)[3]], lev=c(seq(0,90,10),101), cex.pt=1.3, add=TRUE, st.col=NULL, stations=plstat)
    #oldpar <- par(no.readonly=TRUE)
    #layout(matrix(c(1:3,3,3+ seq(1,length(ana.spread)*2), rep(4,2) + length(ana.spread)*2),
    #              length(ana.spread)+3, 2, byrow=T), height=c(5,lcm(2), rep(5, length(ana.spread)),lcm(2)))
    #par(oma=c(0,0,0,0))
    #plot_echam(espread, varname='slp', names=pnames[1:dim(espread$data)[3]], symmetric=FALSE, cex.pt=1.3, add=TRUE)
    #plot_echam(aspread, varname='slp', cols=rbfun(10), names=pnames[dim(espread$data)[3] + 1:dim(aspread$data)[3]], levs=seq(0,100,10), cex.pt=1.3, add=TRUE, st.col=1, stations=plstat) #calibrate)
    #par(oldpar)
    dev.off()
    #    \caption{According to figure \ref{fig:spread} but for mean sea level pressure in winter (left) and summer (right). }
  }
  
  
  ################################################################################
  # Fig. 3a: average temperature update
  
  #<<label=update, echo=FALSE, fig=TRUE, width=8, height=9, results=hide, eps=FALSE>>=
  #pdf(paste(figpath,'avg_upd_temp.pdf',sep='/'), width=4, height=9, paper='special')
  #upd <- lapply(analysis, function(x) apply(array(x$data - echam$data, data.dim), 1:2, mean, na.rm=T))
  upd <- apply(array(analysis$data - echam$data, data.dim), 1:2, mean, na.rm=T)
  
  update <- echam
  #update$data <- array(unlist(upd), c(nrow(echam$data), 1, ncol(upd[[1]])*length(upd)))
  update$data <- array(upd,c(dim(upd)[1],1,dim(upd)[2]))
  
  update$ensmean <- array(update$data,dim=c(dim(update$data)[1],dim(update$data)[3]))
  
  tmax=round(max(abs(update$ensmean[which(update$names=="temp2")])),digits=2)
  
  slpmax=round(max(abs(update$ensmean[which(update$names=="slp")])),digits=1)
  
  prmax=round(max(abs(update$ensmean[which(update$names=="precip")])),digits=1)
  
  
  if (monthly_out) {
    
    plot_echam4(update, varname='temp2', cex.pt=1.3, names=pnames[1:dim(update$data)[3]], type='ensmean', st.col=NULL, stations=plstat,NHseason="summer",plotname=paste0('avg_upd_temp_summer.pdf'), paper='special')
    
    plot_echam4(update, varname='temp2', cex.pt=1.3, names=pnames[1:dim(update$data)[3]], type='ensmean', st.col=NULL, stations=plstat,NHseason="winter",plotname=paste0('avg_upd_temp_winter.pdf'), paper='special')
  } else {
    plot_echam4(update, varname='temp2', cex.pt=1.3, lev=round(tmax/9*c(-Inf,-7.5:-0.5,0.5:7.5,Inf),digits=3), names=pnames[1:dim(update$data)[3]], type='ensmean', st.col=NULL, stations=plstat,NHseason=NULL,plotname=paste0('avg_upd_temp.pdf'), paper='special')
  }
  
  
  #    \caption{Average update of the echam ensemble members according to EnSRF with perfect proxies in a) and b), localized perfect proxies in c) and d), pseudoproxies in e) and f), and localized pseudoproxies in g) and h). Results for winter (October to March) are shown in the left column, results for summer (April to October) in the right column. In case of instrumental data there are no pseudoproxies: a) NO localization winter, b) NO localization summer, c) WITH localization winter, d WITH localization summer}
  
  
  
  ################################################################################
  # Fig. 3b: average precipitation update
  
  
  if (monthly_out) {
    plot_echam4(update, varname='precip', cex.pt=1.3, names=pnames[1:dim(update$data)[3]], type='ensmean', st.col=NULL, stations=plstat,NHseason="summer",plotname=paste0('avg_upd_precip_summer.pdf'), paper='special')
    plot_echam4(update, varname='precip', cex.pt=1.3, names=pnames[1:dim(update$data)[3]], type='ensmean', st.col=NULL, stations=plstat,NHseason="winter",plotname=paste0('avg_upd_precip_winter.pdf'), paper='special')
  } else {
    plot_echam4(update, varname='precip', cex.pt=1.3,lev=round(prmax/9*c(-Inf,-7.5:-0.5,0.5:7.5,Inf),digits=2), names=pnames[1:dim(update$data)[3]], type='ensmean', st.col=NULL, stations=plstat,NHseason=NULL,plotname=paste0('avg_upd_precip.pdf'), paper='special')
  }
  #    \caption{According to figure \ref{fig:update} but for precipitation in winter (left) and summer (right)}
  
  
  
  ################################################################################
  # Fig. 3c: average SLP update
  
  #<<update_slp, echo=FALSE, fig=TRUE, width=8, height=9, results=hide, eps=FALSE>>=
  
  #    \caption{According to figure \ref{fig:update} but for mean sea level pressure in winter (left) and summer (right)}
  
  if (monthly_out) {
    plot_echam4(update, varname='slp', cex.pt=1.3, names=pnames[1:dim(update$data)[3]] , type='ensmean', st.col=NULL, stations=plstat,NHseason="summer",plotname=paste0('avg_upd_slp_summer.pdf'), paper='special')
    plot_echam4(update, varname='slp', cex.pt=1.3, names=pnames[1:dim(update$data)[3]] , type='ensmean', st.col=NULL, stations=plstat,NHseason="winter",plotname=paste0('avg_upd_slp_winter.pdf'), paper='special')
  } else {
    plot_echam4(update, varname='slp', cex.pt=1.3,lev=round(slpmax/9*c(-Inf,-7.5:-0.5,0.5:7.5,Inf),digits=2), names=pnames[1:dim(update$data)[3]], type='ensmean', st.col=NULL, stations=plstat,NHseason=NULL,plotname=paste0('avg_upd_slp.pdf'), paper='special')
  }
  rm(upd,update)
  
  ################################################################################
  # Fig. 4a: average temperature corr echam and analysis vs validation wrt ens. mean 
  #<<label=rmse, echo=FALSE, fig=TRUE, width=8, height=9, results=hide, eps=FALSE>>=
  
  corr.tot <- echam
  corr.tot$data <- array(cbind(corr.ech$ensmean,corr$ensmean),c(dim(corr$ensmean)[1],1,dim(corr$ensmean)[2]*2)) #[,,1:2,drop=F]
  ## it doesn't make sense to put $ensmean into $data
  ## and the choosing type='data'
  ## different approach: put $ensmean into $ensmean and 
  ## choose type='ensmean'
  corr.tot$ensmean <- array(cbind(corr.ech$ensmean,corr$ensmean),c(dim(corr$ensmean)[1],dim(corr$ensmean)[2]*2)) #[,,1:2,drop=F]
  
  corr.diff <- corr.tot
  corr.diff$data <- corr.diff$data[,,(ncol(corr.tot$ensmean)/2+1):ncol(corr.tot$ensmean)]-corr.diff$data[,,1:(ncol(corr.tot$ensmean)/2)]
  corr.diff$ensmean <- corr.diff$ensmean[,(ncol(corr.tot$ensmean)/2+1):ncol(corr.tot$ensmean)]-corr.diff$ensmean[,1:(ncol(corr.tot$ensmean)/2)]
  
  #corr.ech <- corr               # analysis vs. validate
  #corr.ech$Analysis <- corr.ech  # echam vs. validate
  #names(corr.ech) <- c('echam', 'analysis_localized')
  #corr.tot$data <- array(sapply(corr.ech, function(x) x$ensmean), c(nrow(corr.ech[[1]]$ensmean), 1, ncol(corr.ech[[1]]$ensmean)*length(corr.ech)))
  
  #corr.tot$data <- array(corr.tot$data[,,5:8], c(nrow(corr[[1]]$ensmean), 1, 4))
  #corr.tot$data <- array(corr.tot$data[,,], c(nrow(corr[[1]]$ensmean), 1, 4))
  
  # pdf(paste(figpath,'corr_echam_anal-cru_temp.pdf',sep='/'), width=9, height=5, paper='special')
  #pdf(paste(figpath,'corr_echam_anal-cru_temp.pdf',sep='/'), width=9, height=12, paper='special')
  # layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE), height=c(3,3,1))
  
  # layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,13), 7, 2, byrow = TRUE), height=c(3,3,3,3,3,3,1))
  # par(oma=c(0,0,0,0))
  levs <- c(-1,seq(-0.9,0.9,0.2),1)
  #levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0,0.1,0.2,0.3, 0.5, 0.7, 1)
  
  if (monthly_out) {
    plot_echam4(corr.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(corr.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason="summer",plotname=paste0('corr_echam_anal-',v,'_temp_summer.pdf'), paper='special')
    plot_echam4(corr.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(corr.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason="winter",plotname=paste0('corr_echam_anal-',v,'_temp_winter.pdf'), paper='special')
  } else {
    plot_echam4(corr.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(corr.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason=NULL,plotname=paste0('corr_echam_anal-',v,'_temp.pdf'), paper='special')
  }
  
  #plot_echam(corr.tot, cex.pt=1.5, names=pnames[1:dim(corr.tot$data)[3]], lev=levs, st.col=1, stations=plstat,add=TRUE)
  # dev.off()
  #    \caption{temp corr echam ens mean with cru validation (top) and analysis ens mean with cru validation (bottom). Winter left, summer right.}
  
  ################################################################################
  # Fig.4.1a: correlation difference map (analysis-echam) TEMP
  
  
  if (monthly_out) {
    plot_echam4(corr.diff, varname='temp2', cex.pt=1.5, names=pnames[1:dim(corr.diff$data)[2]], lev=levs, type='ensmean', st.col=NULL, stations=plstat,NHseason="summer",plotname=paste0('corr.diff_echam_anal-',v,'_temp_summer.pdf'), paper='special')
    plot_echam4(corr.diff, varname='temp2', cex.pt=1.5, names=pnames[1:dim(corr.diff$data)[2]], lev=levs, type='ensmean', st.col=NULL, stations=plstat,NHseason="winter",plotname=paste0('corr.diff_echam_anal-',v,'_temp_winter.pdf'), paper='special')
  } else {
    plot_echam4(corr.diff, varname='temp2', cex.pt=1.5, names=pnames[1:dim(corr.diff$data)[2]], lev=levs, type='ensmean', st.col=NULL, stations=plstat,NHseason=NULL,plotname=paste0('corr.diff_echam_anal-',v,'_temp.pdf'), paper='special')
  }
  
  ################################################################################
  # Fig. 4b: average precip corr echam and analysis vs validation wrt ens. mean 
  
  # pdf(paste(figpath,'corr_echam_anal-cru_precip.pdf',sep='/'), width=9, height=6, paper='special')
  # layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE), height=c(3,3,1))
  # par(oma=c(0,0,0,0))
  levs <- c(-1,seq(-0.9,0.9,0.2),1)
  
  if (monthly_out){
    plot_echam4(corr.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(corr.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason="summer",plotname=paste0('corr_echam_anal-',v,'_precip_summer.pdf'), paper='special')
    plot_echam4(corr.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(corr.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason="winter",plotname=paste0('corr_echam_anal-',v,'_precip_winter.pdf'), paper='special')
  }else{
    plot_echam4(corr.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(corr.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason=NULL,plotname=paste0('corr_echam_anal-',v,'_precip.pdf'), paper='special')
  }
  # dev.off()
  #    \caption{precip corr echam ens mean with cru validation (top) and analysis ens mean with cru validation (bottom). Winter left, summer right.}
  
  ################################################################################
  # Fig.4.1b: correlation difference map (analysis-echam) PRECIP
  
  if (monthly_out) {
    plot_echam4(corr.diff, varname='precip', cex.pt=1.5, names=pnames[1:dim(corr.diff$data)[2]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason="summer",plotname=paste0('corr.diff_echam_anal-',v,'_precip_summer.pdf'), paper='special')
    plot_echam4(corr.diff, varname='precip', cex.pt=1.5, names=pnames[1:dim(corr.diff$data)[2]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason="winter",plotname=paste0('corr.diff_echam_anal-',v,'_precip_winter.pdf'), paper='special')
  } else {
    plot_echam4(corr.diff, varname='precip', cex.pt=1.5, names=pnames[1:dim(corr.diff$data)[2]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason=NULL,plotname=paste0('corr.diff_echam_anal-',v,'_precip.pdf'), paper='special')
  }
  
  ################################################################################
  # Fig. 4c: average slp corr echam and analysis vs validation wrt ens. mean 
  
  levs <- c(-1,seq(-0.9,0.9,0.2),1)
  
  if (monthly_out){
    plot_echam4(corr.tot, varname='slp', cex.pt=1.5, names=pnames[1:dim(corr.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason="summer",plotname=paste0('corr_echam_anal-',v,'_slp_summer.pdf'), paper='special')
    plot_echam4(corr.tot, varname='slp', cex.pt=1.5, names=pnames[1:dim(corr.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason="winter",plotname=paste0('corr_echam_anal-',v,'_slp_winter.pdf'), paper='special')
  }else{
    plot_echam4(corr.tot, varname='slp', cex.pt=1.5, names=pnames[1:dim(corr.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason=NULL,plotname=paste0('corr_echam_anal-',v,'_slp.pdf'), paper='special')
  }
  
  
  #    \caption{slp corr echam ens mean with cru validation (top) and analysis ens mean with cru validation (bottom). Winter left, summer right.}
  
  ################################################################################
  # Fig.4.1c: correlation difference map (analysis-echam) SLP
  
  if (monthly_out) {
    plot_echam4(corr.diff, varname='slp', cex.pt=1.5, names=pnames[1:dim(corr.diff$data)[2]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason="summer",plotname=paste0('corr.diff_echam_anal-',v,'_slp_summer.pdf'), paper='special')
    plot_echam4(corr.diff, varname='slp', cex.pt=1.5, names=pnames[1:dim(corr.diff$data)[2]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason="winter",plotname=paste0('corr.diff_echam_anal-',v,'_slp_winter.pdf'), paper='special')
  } else {
    plot_echam4(corr.diff, varname='slp', cex.pt=1.5, names=pnames[1:dim(corr.diff$data)[2]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason=NULL,plotname=paste0('corr.diff_echam_anal-',v,'_slp.pdf'), paper='special')
  }
  
  ################################################################################
  # Fig. 5a: average temperature bias echam and analysis vs. validate wrt ens. mean 
  #<<label=rmse, echo=FALSE, fig=TRUE, width=8, height=9, results=hide, eps=FALSE>>=
  if (!monthly_out){
    bias.tot <- echam
    bias.tot$data <- array(cbind(bias.ech$ensmean,bias$ensmean),c(dim(bias$ensmean)[1],1,dim(bias$ensmean)[2]*2))
    bias.tot$ensmean <- array(cbind(bias.ech$ensmean,bias$ensmean),c(dim(bias$ensmean)[1],dim(bias$ensmean)[2]*2))
    
    #bias.ech <- bias
    #bias.ech$Analysis <- bias.ech
    #names(bias.ech) <- c('echam', 'analysis_localized')
    #bias.tot$data <- array(sapply(bias.ech, function(x) x$ensmean), c(nrow(bias.ech[[1]]$ensmean), 1, ncol(bias.ech[[1]]$ensmean)*length(bias.ech)))
    
    #bias.tot$data <- array(bias.tot$data[,,5:8], c(nrow(bias[[1]]$ensmean), 1, 4))
    #bias.tot$data <- array(bias.tot$data[,,], c(nrow(bias[[1]]$ensmean), 1, 4))
    
    
    #    \caption{}
    
    levs <- c(-Inf,-10,-5,-2,-1,0,1,2,5,10,Inf)
    
    plot_echam4(bias.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(bias.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason=NULL,plotname=paste0('bias_echam_anal-',v,'_temp.pdf'), paper='special')
    
    
    ################################################################################
    # Fig. 5b: average precip bias echam and analysis vs. validate wrt ens. mean 
    
    
    
    levs <- c(-Inf,-40,-20,-10,-5,0,5,10,20,40,Inf)
    
    plot_echam4(bias.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(bias.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason=NULL,plotname=paste0('bias_echam_anal-',v,'_precip.pdf'), paper='special')
    
    ################################################################################
    # Fig. 5c: average precip bias echam and analysis vs. validate wrt ens. mean 
    
    
    levs <- c(-Inf,-10,-5,-2,-1,0,1,2,5,10,Inf)
    
    plot_echam4(bias.tot, varname='slp', cex.pt=1.5, names=pnames[1:dim(bias.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason=NULL,plotname=paste0('bias_echam_anal-',v,'_slp.pdf'), paper='special')
    
  }
  ################################################################################
  # Fig. 6a: average temperature RE wrt ens. mean 
  # NEW RE BASED ON ANOMALIES BECAUSE WITH BIASES BETWEEN MODEL AND VALIDATION TARGET WE COULD NOT REACH 1
  #<<label=rmse, echo=FALSE, fig=TRUE, width=8, height=9, results=hide, eps=FALSE>>=
  #rmse <- echam
  #rmse$data <- array(unlist(RE), c(nrow(rmse.ech), 1, ncol(rmse.ech)*length(RE)))
  RE.tot <- echam
  #RE.bkp <- RE
  #RE <- RE.bkp
  #anomalies and only 1 analysis with distance weighting
  RE.tot$data <- array(RE.anom$ensmean,c(dim(RE$ensmean)[1], 1, dim(RE$ensmean)[2]))
  RE.tot$ensmean <- array(RE.anom$ensmean,c(dim(RE$ensmean)[1], dim(RE$ensmean)[2]))
  
  levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0.05,0.2,0.4,0.6,0.8,1)
  
  if (monthly_out){
    plot_echam4(RE.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason="summer",plotname=paste0('re_anom_echam_anal-',v,'_temp_summer.pdf'), paper='special')
    plot_echam4(RE.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason="winter",plotname=paste0('re_anom_echam_anal-',v,'_temp_winter.pdf'), paper='special')
  }else{
    plot_echam4(RE.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason=NULL,plotname=paste0('re_anom_echam_anal-',v,'_temp.pdf'), paper='special')
  }
  
  #RE.tot$data <- array(cbind(RE.anom$ensmean,RE.ech.clim.anom$ensmean), 
  #                     c(dim(RE$ensmean)[1], 1, (dim(RE$ensmean)[2]*2)))
  #RE.tot$data <- array((rmse.anom$ensmean-rmse.clim.anom$ensmean), 
  #                     c(dim(RE$ensmean)[1], 1, (dim(RE$ensmean)[2])))
  
  #if (real_proxies) {
  #RE.tot$data <- array(RE$ensmean, c(dim(RE$ensmean)[1], 1, dim(RE$ensmean)[2]))
  #}
  #RE.tot$data <- array(sapply(RE, function(x) x$ensmean), c(nrow(RE[[1]]$ensmean), 1, ncol(RE[[1]]$ensmean)*length(RE)))
  
  # # only pseudoproxies including noise
  # #RE.tot$data <- array(RE.tot$data[,,5:8], c(nrow(RE[[1]]$ensmean), 1, 4))
  # # instrumental data (only 4 fields (winter,summer,w/o and w localisation) because no noise added)
  # #RE.tot$data <- array(RE.tot$data[,,1:4], c(nrow(RE[[1]]$ensmean), 1, 4))
  # pdf(paste(figpath,'/re_anom_echam_anal-cru_temp.pdf',sep='/'), width=9, height=3.5, paper='special')
  # #pdf(paste(figpath,'/re_clim_echam-cru_temp.pdf',sep='/'), width=9, height=6, paper='special')
  # #pdf(paste(figpath,'/re_echam_anal-cru_temp_12mon.pdf',sep='/'), width=9, height=12, paper='special')
  # layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
  # #layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE), height=c(3,3,1))
  # #layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,13), 7, 2, byrow = TRUE), height=c(3,3,3,3,3,3,1))
  # par(oma=c(0,0,0,0))
  # levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0.05,0.2,0.4,0.6,0.8,1)
  # #levs <- c(-Inf,seq(-1,1,0.2),Inf)
  # plot_echam(RE.tot, cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, st.col=NULL,
  #            stations=calibrate,add=TRUE)
  # #plot_echam(RE.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, st.col=1, stations=plstat,add=TRUE)
  # dev.off()
  # #    \caption{According to figure \ref{fig:update} but for average reduction of error of the analysis ensemble mean with respect to the ECHAM5 ensemble mean. Positive values indicate that the analysis is closer to the validation measurement than the unconstrained ensemble mean in the respective season.}
  
  
  
  ################################################################################
  # Fig. 6b: average precip RE wrt ens. mean 
  #<<label=rmse, echo=FALSE, fig=TRUE, width=8, height=9, results=hide, eps=FALSE>>=
  if (monthly_out){
    plot_echam4(RE.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason="summer",plotname=paste0('re_anom_echam_anal-',v,'_precip_summer.pdf'), paper='special')
    plot_echam4(RE.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason="winter",plotname=paste0('re_anom_echam_anal-',v,'_precip_winter.pdf'), paper='special')
  }else{
    plot_echam4(RE.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason=NULL,plotname=paste0('re_anom_echam_anal-',v,'_precip.pdf'), paper='special')
  }
  
  # 
  # pdf(paste(figpath,'re_anom_echam_anal-cru_precip.pdf',sep='/'), width=9, height=3.5, paper='special')
  # layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
  # #layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE), height=c(3,3,1))
  # par(oma=c(0,0,0,0))
  # levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0,0.2,0.4,0.6,0.8,1)
  # plot_echam(RE.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, st.col=NULL, stations=plstat,add=TRUE)
  # dev.off()
  # #    \caption{According to figure \ref{fig:update} but for average reduction of error of the analysis ensemble mean with respect to the ECHAM5 ensemble mean. Positive values indicate that the analysis is closer to the validation measurement than the unconstrained ensemble mean in the respective season.}
  
  
  ################################################################################
  # Fig. 6c: average slp RE wrt ens. mean 
  
  
  #    \caption{According to figure \ref{fig:update} but for average reduction of error of the analysis ensemble mean with respect to the ECHAM5 ensemble mean. Positive values indicate that the analysis is closer to the validation measurement than the unconstrained ensemble mean in the respective season.}
  
  levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0,0.2,0.4,0.6,0.8,1)
  
  
  if (monthly_out){
    plot_echam4(RE.tot, varname='slp', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason="summer",plotname=paste0('re_anom_echam_anal-',v,'_slp_summer.pdf'), paper='special')
    plot_echam4(RE.tot, varname='slp', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason="winter",plotname='re_anom_echam_anal-',v,'_slp_winter.pdf', paper='special')
  }else{
    plot_echam4(RE.tot, varname='slp', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason=NULL,plotname=paste0('re_anom_echam_anal-',v,'_slp.pdf'), paper='special')
  }
  
  ################################################################################
  # Fig. 6a: average temperature RE wrt ens. mean 
  # OLD version without anomalies
  # RE BASED ON ANOMALIES BECAUSE WITH BIASES BETWEEN MODEL AND VALIDATION TARGET WE COULD NOT REACH 1
  #<<label=rmse, echo=FALSE, fig=TRUE, width=8, height=9, results=hide, eps=FALSE>>=
  #rmse <- echam
  #rmse$data <- array(unlist(RE), c(nrow(rmse.ech), 1, ncol(rmse.ech)*length(RE)))
  RE.tot <- echam
  #RE.bkp <- RE
  #RE <- RE.bkp
  #anomalies and only 1 analysis with distance weighting
  RE.tot$data <- array(RE$ensmean,c(dim(RE$ensmean)[1], 1, dim(RE$ensmean)[2]))
  RE.tot$ensmean <- array(RE$ensmean,c(dim(RE$ensmean)[1], dim(RE$ensmean)[2]))
  
  
  levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0.05,0.2,0.4,0.6,0.8,1)
  if (monthly_out){
    plot_echam4(RE.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason="summer",plotname=paste0('re_echam_anal-',v,'_temp_summer.pdf'), paper='special')
    plot_echam4(RE.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason="winter",plotname=paste0('re_echam_anal-',v,'_temp_winter.pdf'), paper='special')
  }else{
    plot_echam4(RE.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason=NULL,plotname=paste0('re_echam_anal-',v,'_temp.pdf'), paper='special')
  }
  #RE.tot$data <- array(cbind(RE.anom$ensmean,RE.ech.clim.anom$ensmean), 
  #                     c(dim(RE$ensmean)[1], 1, (dim(RE$ensmean)[2]*2)))
  #RE.tot$data <- array((rmse.anom$ensmean-rmse.clim.anom$ensmean), 
  #                     c(dim(RE$ensmean)[1], 1, (dim(RE$ensmean)[2])))
  
  #if (real_proxies) {
  #RE.tot$data <- array(RE$ensmean, c(dim(RE$ensmean)[1], 1, dim(RE$ensmean)[2]))
  #}
  #RE.tot$data <- array(sapply(RE, function(x) x$ensmean), c(nrow(RE[[1]]$ensmean), 1, ncol(RE[[1]]$ensmean)*length(RE)))
  
  # # only pseudoproxies including noise 
  # #RE.tot$data <- array(RE.tot$data[,,5:8], c(nrow(RE[[1]]$ensmean), 1, 4)) 
  # # instrumental data (only 4 fields (winter,summer,w/o and w localisation) because no noise added)
  # #RE.tot$data <- array(RE.tot$data[,,1:4], c(nrow(RE[[1]]$ensmean), 1, 4))
  # pdf(paste(figpath,'re_echam_anal-cru_temp.pdf',sep='/'), width=9, height=3.5, paper='special')
  # #pdf(paste(figpath,'re_clim_echam-cru_temp.pdf',sep='/'), width=9, height=6, paper='special')
  # #pdf(paste(figpath,'re_echam_anal-cru_temp_12mon.pdf',sep='/'), width=9, height=12, paper='special')
  # layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
  # #layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE), height=c(3,3,1))
  # #layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,13), 7, 2, byrow = TRUE), height=c(3,3,3,3,3,3,1))
  # par(oma=c(0,0,0,0))
  # levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0.05,0.2,0.4,0.6,0.8,1)
  # #levs <- c(-Inf,seq(-1,1,0.2),Inf)
  # plot_echam(RE.tot, cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, st.col=NULL, 
  #            stations=calibrate,add=TRUE)
  # #plot_echam(RE.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, st.col=1, stations=plstat,add=TRUE)
  # dev.off()
  # #    \caption{According to figure \ref{fig:update} but for average reduction of error of the analysis ensemble mean with respect to the ECHAM5 ensemble mean. Positive values indicate that the analysis is closer to the validation measurement than the unconstrained ensemble mean in the respective season.}
  
  
  
  ################################################################################
  # Fig. 6b: average precip RE wrt ens. mean 
  #<<label=rmse, echo=FALSE, fig=TRUE, width=8, height=9, results=hide, eps=FALSE>>=
  
  if (monthly_out){
    plot_echam4(RE.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason="summer",plotname=paste0('re_echam_anal-',v,'_precip_summer.pdf'), paper='special')
    plot_echam4(RE.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason="winter",plotname=paste0('re_echam_anal-',v,'_precip_winter.pdf'), paper='special')
  }else{
    plot_echam4(RE.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason=NULL,plotname=paste0('re_echam_anal-',v,'_precip.pdf'), paper='special')
  }
  # 
  # 
  # pdf(paste(figpath,'re_echam_anal-cru_precip.pdf',sep='/'), width=9, height=3.5, paper='special')
  # layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
  # #layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE), height=c(3,3,1))
  # par(oma=c(0,0,0,0))
  # levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0,0.2,0.4,0.6,0.8,1)
  # plot_echam(RE.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, st.col=NULL, stations=plstat,add=TRUE)
  # dev.off()
  # #    \caption{According to figure \ref{fig:update} but for average reduction of error of the analysis ensemble mean with respect to the ECHAM5 ensemble mean. Positive values indicate that the analysis is closer to the validation measurement than the unconstrained ensemble mean in the respective season.}
  
  
  ################################################################################
  # Fig. 6c: average slp RE wrt ens. mean 
  
  
  #    \caption{According to figure \ref{fig:update} but for average reduction of error of the analysis ensemble mean with respect to the ECHAM5 ensemble mean. Positive values indicate that the analysis is closer to the validation measurement than the unconstrained ensemble mean in the respective season.}
  
  levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0,0.2,0.4,0.6,0.8,1)
  
  if (monthly_out){
    plot_echam4(RE.tot, varname='slp', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason="summer",plotname=paste0('re_echam_anal-',v,'_slp_summer.pdf'), paper='special')
    plot_echam4(RE.tot, varname='slp', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason="winter",plotname=paste0('re_echam_anal-',v,'_slp_winter.pdf'), paper='special')
  }else{
    plot_echam4(RE.tot, varname='slp', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason=NULL,plotname=paste0('re_echam_anal-',v,'_slp.pdf'), paper='special')
  }
  ################################################################################
  # Fig. 6a: average temperature RE wrt ens. mean 
  # NEW RE BASED ON ANOMALIES BECAUSE WITH BIASES BETWEEN MODEL AND VALIDATION TARGET WE COULD NOT REACH 1
  #<<label=rmse, echo=FALSE, fig=TRUE, width=8, height=9, results=hide, eps=FALSE>>=
  #rmse <- echam
  #rmse$data <- array(unlist(RE), c(nrow(rmse.ech), 1, ncol(rmse.ech)*length(RE)))
  RE.tot <- echam
  #RE.bkp <- RE
  #RE <- RE.bkp
  #anomalies and only 1 analysis with distance weighting
  RE.tot$data <- array(RE.echam.clim$ensmean,c(dim(RE$ensmean)[1], 1, dim(RE$ensmean)[2]))
  RE.tot$ensmean <- array(RE.echam.clim$ensmean,c(dim(RE$ensmean)[1],dim(RE$ensmean)[2]))
  
  
  levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0.05,0.2,0.4,0.6,0.8,1)
  
  if (monthly_out){
    plot_echam4(RE.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason="summer",plotname=paste0('re_clim_echam_anal-',v,'_temp_summer.pdf'), paper='special')
    plot_echam4(RE.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason="winter",plotname=paste0('re_clim_echam_anal-',v,'_temp_winter.pdf'), paper='special')
  }else{
    plot_echam4(RE.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason=NULL,plotname=paste0('re_clim_echam_anal-',v,'_temp.pdf'), paper='special')
  }
  #RE.tot$data <- array(cbind(RE.anom$ensmean,RE.ech.clim.anom$ensmean), 
  #                     c(dim(RE$ensmean)[1], 1, (dim(RE$ensmean)[2]*2)))
  #RE.tot$data <- array((rmse.anom$ensmean-rmse.clim.anom$ensmean), 
  #                     c(dim(RE$ensmean)[1], 1, (dim(RE$ensmean)[2])))
  
  #if (real_proxies) {
  #RE.tot$data <- array(RE$ensmean, c(dim(RE$ensmean)[1], 1, dim(RE$ensmean)[2]))
  #}
  #RE.tot$data <- array(sapply(RE, function(x) x$ensmean), c(nrow(RE[[1]]$ensmean), 1, ncol(RE[[1]]$ensmean)*length(RE)))
  
  # only pseudoproxies including noise 
  #RE.tot$data <- array(RE.tot$data[,,5:8], c(nrow(RE[[1]]$ensmean), 1, 4)) 
  # instrumental data (only 4 fields (winter,summer,w/o and w localisation) because no noise added)
  #RE.tot$data <- array(RE.tot$data[,,1:4], c(nrow(RE[[1]]$ensmean), 1, 4))
  # pdf(paste(figpath,'re_clim_echam_anal-cru_temp.pdf',sep='/'), width=9, height=3.5, paper='special')
  # #pdf(paste(figpath,'re_clim_echam-cru_temp.pdf',sep='/'), width=9, height=6, paper='special')
  # #pdf(paste(figpath,'re_echam_anal-cru_temp_12mon.pdf',sep='/'), width=9, height=12, paper='special')
  # layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
  # #layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE), height=c(3,3,1))
  # #layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,13), 7, 2, byrow = TRUE), height=c(3,3,3,3,3,3,1))
  # par(oma=c(0,0,0,0))
  # levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0.05,0.2,0.4,0.6,0.8,1)
  # #levs <- c(-Inf,seq(-1,1,0.2),Inf)
  # plot_echam(RE.tot, cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, st.col=NULL, 
  #            stations=calibrate,add=TRUE)
  # #plot_echam(RE.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, st.col=1, stations=plstat,add=TRUE)
  # dev.off()
  # #    \caption{According to figure \ref{fig:update} but for average reduction of error of the analysis ensemble mean with respect to the ECHAM5 ensemble mean. Positive values indicate that the analysis is closer to the validation measurement than the unconstrained ensemble mean in the respective season.}
  
  
  
  ################################################################################
  # Fig. 6b: average precip RE wrt ens. mean 
  #<<label=rmse, echo=FALSE, fig=TRUE, width=8, height=9, results=hide, eps=FALSE>>=
  
  if (monthly_out){
    plot_echam4(RE.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason="summer",plotname=paste0('re_clim_echam_anal-',v,'_precip_summer.pdf'), paper='special')
    plot_echam4(RE.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason="winter",plotname=paste0('re_clim_echam_anal-',v,'_precip_winter.pdf'), paper='special')
  }else{
    plot_echam4(RE.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason=NULL,plotname=paste0('re_clim_echam_anal-',v,'_precip.pdf'), paper='special')
  }
  
  # 
  # pdf(paste(figpath,'re_clim_echam_anal-cru_precip.pdf',sep='/'), width=9, height=3.5, paper='special')
  # layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
  # #layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE), height=c(3,3,1))
  # par(oma=c(0,0,0,0))
  # levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0,0.2,0.4,0.6,0.8,1)
  # plot_echam(RE.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, st.col=NULL, stations=plstat,add=TRUE)
  # dev.off()
  #    \caption{According to figure \ref{fig:update} but for average reduction of error of the analysis ensemble mean with respect to the ECHAM5 ensemble mean. Positive values indicate that the analysis is closer to the validation measurement than the unconstrained ensemble mean in the respective season.}
  
  
  ################################################################################
  # Fig. 6c: average slp RE wrt ens. mean 
  
  #    \caption{According to figure \ref{fig:update} but for average reduction of error of the analysis ensemble mean with respect to the ECHAM5 ensemble mean. Positive values indicate that the analysis is closer to the validation measurement than the unconstrained ensemble mean in the respective season.}
  
  levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0,0.2,0.4,0.6,0.8,1)
  
  if (monthly_out){
    plot_echam4(RE.tot, varname='slp', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason="summer",plotname=paste0('re_clim_echam_anal-',v,'_slp_summer.pdf'), paper='special')
    plot_echam4(RE.tot, varname='slp', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason="winter",plotname=paste0('re_clim_echam_anal-',v,'_slp_winter.pdf'), paper='special')
  }else{
    plot_echam4(RE.tot, varname='slp', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason=NULL,plotname=paste0('re_clim_echam_anal-',v,'_slp.pdf'), paper='special')
  }
  #################################################################################
  
  RE.tot <- echam
  #RE.bkp <- RE
  #RE <- RE.bkp
  #anomalies and only 1 analysis with distance weighting
  RE.tot$data <- array(RE.echam.clim.anom$ensmean,c(dim(RE$ensmean)[1], 1, dim(RE$ensmean)[2]))
  RE.tot$ensmean <- array(RE.echam.clim.anom$ensmean,c(dim(RE$ensmean)[1], dim(RE$ensmean)[2]))
  
  levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0.05,0.2,0.4,0.6,0.8,1)
  
  if (monthly_out){
    plot_echam4(RE.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason="summer",plotname=paste0('re_clim_anom_echam_anal-',v,'_temp_summer.pdf'), paper='special')
    plot_echam4(RE.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason="winter",plotname=paste0('re_clim_anom_echam_anal-',v,'_temp_winter.pdf'), paper='special')
  }else{
    plot_echam4(RE.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason=NULL,plotname=paste0('re_clim_anom_echam_anal-',v,'_temp.pdf'), paper='special')
  }
  #RE.tot$data <- array(cbind(RE.anom$ensmean,RE.ech.clim.anom$ensmean), 
  #                     c(dim(RE$ensmean)[1], 1, (dim(RE$ensmean)[2]*2)))
  #RE.tot$data <- array((rmse.anom$ensmean-rmse.clim.anom$ensmean), 
  #                     c(dim(RE$ensmean)[1], 1, (dim(RE$ensmean)[2])))
  
  #if (real_proxies) {
  #RE.tot$data <- array(RE$ensmean, c(dim(RE$ensmean)[1], 1, dim(RE$ensmean)[2]))
  #}
  #RE.tot$data <- array(sapply(RE, function(x) x$ensmean), c(nrow(RE[[1]]$ensmean), 1, ncol(RE[[1]]$ensmean)*length(RE)))
  
  # # only pseudoproxies including noise 
  # #RE.tot$data <- array(RE.tot$data[,,5:8], c(nrow(RE[[1]]$ensmean), 1, 4)) 
  # # instrumental data (only 4 fields (winter,summer,w/o and w localisation) because no noise added)
  # #RE.tot$data <- array(RE.tot$data[,,1:4], c(nrow(RE[[1]]$ensmean), 1, 4))
  # pdf(paste(figpath,'/re_clim_anom_echam_anal-cru_temp.pdf',sep='/'), width=9, height=3.5, paper='special')
  # #pdf(paste(figpath,'/re_clim_echam-cru_temp.pdf',sep='/'), width=9, height=6, paper='special')
  # #pdf(paste(figpath,'/re_echam_anal-cru_temp_12mon.pdf',sep='/'), width=9, height=12, paper='special')
  # layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
  # #layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE), height=c(3,3,1))
  # #layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,13), 7, 2, byrow = TRUE), height=c(3,3,3,3,3,3,1))
  # par(oma=c(0,0,0,0))
  # levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0.05,0.2,0.4,0.6,0.8,1)
  # #levs <- c(-Inf,seq(-1,1,0.2),Inf)
  # plot_echam(RE.tot, cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, st.col=NULL, 
  #            stations=calibrate,add=TRUE)
  # #plot_echam(RE.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, st.col=1, stations=plstat,add=TRUE)
  # dev.off()
  # #    \caption{According to figure \ref{fig:update} but for average reduction of error of the analysis ensemble mean with respect to the ECHAM5 ensemble mean. Positive values indicate that the analysis is closer to the validation measurement than the unconstrained ensemble mean in the respective season.}
  
  
  
  ################################################################################
  # Fig. 6b: average precip RE wrt ens. mean 
  #<<label=rmse, echo=FALSE, fig=TRUE, width=8, height=9, results=hide, eps=FALSE>>=
  
  if (monthly_out){
    plot_echam4(RE.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason="summer",plotname=paste0('re_clim_anom_echam_anal-',v,'_precip_summer.pdf'), paper='special')
    plot_echam4(RE.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason="winter",plotname=paste0('re_clim_anom_echam_anal-',v,'_precip_winter.pdf'), paper='special')
  }else{
    plot_echam4(RE.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason=NULL,plotname=paste0('re_clim_anom_echam_anal-',v,'_precip.pdf'), paper='special')
  }
  
  # 
  # pdf(paste(figpath,'re_clim_anom_echam_anal-cru_precip.pdf',sep='/'), width=9, height=3.5, paper='special')
  # layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
  # #layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE), height=c(3,3,1))
  # par(oma=c(0,0,0,0))
  # levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0,0.2,0.4,0.6,0.8,1)
  # plot_echam(RE.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, st.col=NULL, stations=plstat,add=TRUE)
  # dev.off()
  # #    \caption{According to figure \ref{fig:update} but for average reduction of error of the analysis ensemble mean with respect to the ECHAM5 ensemble mean. Positive values indicate that the analysis is closer to the validation measurement than the unconstrained ensemble mean in the respective season.}
  
  
  ################################################################################
  # Fig. 6c: average slp RE wrt ens. mean 
  
  
  #    \caption{According to figure \ref{fig:update} but for average reduction of error of the analysis ensemble mean with respect to the ECHAM5 ensemble mean. Positive values indicate that the analysis is closer to the validation measurement than the unconstrained ensemble mean in the respective season.}
  
  levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0,0.2,0.4,0.6,0.8,1)
  
  if (monthly_out){
    plot_echam4(RE.tot, varname='slp', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason="summer",plotname=paste0('re_clim_anom_echam_anal-',v,'_slp_summer.pdf'), paper='special')
    plot_echam4(RE.tot, varname='slp', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason="winter",plotname=paste0('re_clim_anom_echam_anal-',v,'_slp_winter.pdf'), paper='special')
  }else{
    plot_echam4(RE.tot, varname='slp', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs , type='ensmean', st.col=NULL, stations=plstat,NHseason=NULL,plotname=paste0('re_clim_anom_echam_anal-',v,'_slp.pdf'), paper='special')
  }
  
  
  
  
  
  #####################################################################################
  ################################# rmse PLOTS ########################################
  #####################################################################################
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
      plot_echam4(rmse, varname=vari, cex.pt=1.5, names=pnames[1:dim(rmse$data)[2]],  type='ensmean', st.col=NULL, stations=plstat,NHseason="summer",plotname=paste0('rmse_anal-',v,'_',vari,'_summer.pdf'), paper='special')
      plot_echam4(rmse, varname=vari, cex.pt=1.5, names=pnames[1:dim(rmse$data)[2]],  type='ensmean', st.col=NULL, stations=plstat,NHseason="winter",plotname=paste0('rmse_anal-',v,'_',vari,'_winter.pdf'), paper='special')
    }else{
      plot_echam4(rmse, varname=vari, cex.pt=1.5, names=pnames[1:dim(rmse$data)[2]],lev=levs,cols=c('lightblue','royalblue','navy'), type='ensmean', st.col=NULL, stations=plstat,NHseason=NULL,plotname=paste0('rmse_anal-',v,'_',vari,'.pdf'), paper='special')
    }
    
    ###################################################################################
    # Fig.: rmse analysis anom
    levs<-c(pretty(rmse.anom$ensmean[which(rmse.anom$names==vari)],n=11),Inf)
    if (monthly_out){
      plot_echam4(rmse.anom, varname=vari, cex.pt=1.5, names=pnames[1:dim(rmse.anom$data)[2]], symmetric=F, type='ensmean', st.col=NULL, stations=plstat,NHseason="summer",plotname=paste0('rmse_anom_anal-',v,'_',vari,'_summer.pdf'), paper='special')
      plot_echam4(rmse.anom, varname=vari, cex.pt=1.5, names=pnames[1:dim(rmse.anom$data)[2]], symmetric=F, type='ensmean', st.col=NULL, stations=plstat,NHseason="winter",plotname=paste0('rmse_anom_anal-',v,'_',vari,'_winter.pdf'), paper='special')
    }else{
      plot_echam4(rmse.anom, varname=vari, cex.pt=1.5, names=pnames[1:dim(rmse.anom$data)[2]],lev=levs,cols=c('lightblue','royalblue','navy'), type='ensmean', st.col=NULL, stations=plstat,NHseason=NULL,plotname=paste0('rmse_anom_anal-',v,'_',vari,'.pdf'), paper='special')
    }
    
    ###################################################################################
    # Fig.: rmse echam
    levs<-pretty(rmse.ech$ensmean[which(rmse.ech$names==vari)],n=11)
    if (monthly_out){
      plot_echam4(rmse.ech, varname=vari, cex.pt=1.5, names=pnames[1:dim(rmse.ech$data)[2]], symmetric=F, type='ensmean', st.col=NULL, stations=plstat,NHseason="summer",plotname=paste0('rmse_echam-',v,'_',vari,'_summer.pdf'), paper='special')
      plot_echam4(rmse.ech, varname=vari, cex.pt=1.5, names=pnames[1:dim(rmse.ech$data)[2]], symmetric=F, type='ensmean', st.col=NULL, stations=plstat,NHseason="winter",plotname=paste0('rmse_echam-',v,'_',vari,'_winter.pdf'), paper='special')
    }else{
      plot_echam4(rmse.ech, varname=vari, cex.pt=1.5, names=pnames[1:dim(rmse.ech$data)[2]],lev=levs,cols=c('lightblue','royalblue','navy'), type='ensmean', st.col=NULL, stations=plstat,NHseason=NULL,plotname=paste0('rmse_echam-',v,'_',vari,'.pdf'), paper='special')
    }
    
    ###################################################################################
    # Fig.: rmse echam anom
    levs<-c(pretty(rmse.ech.anom$ensmean[which(rmse.ech.anom$names==vari)],n=11),Inf)
    if (monthly_out){
      plot_echam4(rmse.ech.anom, varname=vari, cex.pt=1.5, names=pnames[1:dim(rmse.ech.anom$data)[2]], symmetric=F, type='ensmean', st.col=NULL, stations=plstat,NHseason="summer",plotname=paste0('rmse_anom_echam-',v,'_',vari,'_summer.pdf'), paper='special')
      plot_echam4(rmse.ech.anom, varname=vari, cex.pt=1.5, names=pnames[1:dim(rmse.ech.anom$data)[2]], symmetric=F, type='ensmean', st.col=NULL, stations=plstat,NHseason="winter",plotname=paste0('rmse_anom_echam-',v,'_',vari,'_winter.pdf'), paper='special')
    }else{
      plot_echam4(rmse.ech.anom, varname=vari, cex.pt=1.5, names=pnames[1:dim(rmse.ech.anom$data)[2]],lev=levs,cols=c('lightblue','royalblue','navy'), type='ensmean', st.col=NULL, stations=plstat,NHseason=NULL,plotname=paste0('rmse_anom_echam-',v,'_',vari,'.pdf'), paper='special')
    }
    
    ###################################################################################
    #Fig.: rmse diff
    levs<-c(-Inf,pretty(rmse.diff$ensmean[which(rmse.diff$names==vari)],n=10),Inf)
    if (monthly_out){
      plot_echam4(rmse.diff, varname=vari, cex.pt=1.5, names=pnames[1:dim(rmse.diff$data)[2]], symmetric=F, type='ensmean', st.col=NULL, stations=plstat,NHseason="summer",plotname=paste0('rmse_diff-',v,'_',vari,'_summer.pdf'), paper='special')
      plot_echam4(rmse.diff, varname=vari, cex.pt=1.5, names=pnames[1:dim(rmse.diff$data)[2]], symmetric=F, type='ensmean', st.col=NULL, stations=plstat,NHseason="winter",plotname=paste0('rmse_diff-',v,'_',vari,'_winter.pdf'), paper='special')
    }else{
      plot_echam4(rmse.diff, varname=vari, cex.pt=1.5, names=pnames[1:dim(rmse.diff$data)[2]],lev=levs,cols=c('lightblue','royalblue','navy'), type='ensmean', st.col=NULL, stations=plstat,NHseason=NULL,plotname=paste0('rmse_diff-',v,'_',vari,'.pdf'), paper='special')
    }
    
    ###################################################################################
    #Fig.: rmse.anom diff
    levs<-c(-Inf,pretty(rmse.anom.diff$ensmean[which(rmse.anom.diff$names==vari)],n=10),Inf)
    if (monthly_out){
      plot_echam4(rmse.anom.diff, varname=vari, cex.pt=1.5, names=pnames[1:dim(rmse.anom.diff$data)[2]], symmetric=F, type='ensmean', st.col=NULL, stations=plstat,NHseason="summer",plotname=paste0('rmse_diff_anom-',v,'_',vari,'_summer.pdf'), paper='special')
      plot_echam4(rmse.anom.diff, varname=vari, cex.pt=1.5, names=pnames[1:dim(rmse.anom.diff$data)[2]],symmetric=F,  type='ensmean', st.col=NULL, stations=plstat,NHseason="winter",plotname=paste0('rmse_diff_anom-',v,'_',vari,'_winter.pdf'), paper='special')
    }else{
      plot_echam4(rmse.anom.diff, varname=vari, cex.pt=1.5, names=pnames[1:dim(rmse.anom.diff$data)[2]],lev=levs,cols=c('lightblue','royalblue','navy'), type='ensmean', st.col=NULL, stations=plstat,NHseason=NULL,plotname=paste0('rmse_diff_anom-',v,'_',vari,'.pdf'), paper='special')
    }
  } #end rmse plots for loop
  
  
  
  
  
  
  
  
  
  
  
  
  
  # 
  # 
  # ################################################################################
  # # Fig. 7a: validation time series 
  # 
  # #<<label=indices, echo=FALSE, fig=TRUE, width=8, height=6, results=hide, eps=FALSE>>=
  # #
  # pdf(paste(figpath,'/vali_ind_1.pdf',sep='/'), width=9, height=6, paper='special')
  # indices <- setdiff(unique(echam[['names']]), c('temp2', 'precip', 'slp','bias'))
  # par(mfrow=c(ceiling(length(indices)/4), 2), cex.axis=1.4, cex.lab=1.4, mar=c(3,5,1,1), oma=c(0,0,0,0))
  # # create normalized echam and analysis because cru is normalized (only mean = 0)
  # if ((instrumental) || (inst_at_proxy) || (real_proxy)) {
  #   echam.ind=echam
  #   l1=array(NA,c(length(indices),dim(echam$ensmean)[2]))
  #   l2=array(NA,c(length(indices),dim(echam$data)[2:3]))
  #   echam.ind$ensmean=l1
  #   echam.ind$data=l2
  #   echam.ind$lon=0
  #   echam.ind$lat=0
  #   echam.ind$height=0
  #   echam.ind$names=indices[1:length(indices)]
  #   #echam.ind=list(l1,l2)
  #   #names(echam.ind)=c('ensmean','data')
  #   rownames(echam.ind$ensmean)=indices[1:length(indices)]
  #   rownames(echam.ind$data)=indices[1:length(indices)]
  #   for (i in indices[1:length(indices)]){
  #     # split into summer and winter for separate means
  #     tmp <- array(echam$ensmean[echam$names == i,], c(2,dim(echam$ensmean)[2]/2))
  #     tmpmean=rowMeans(tmp)
  #     tmpano=tmp #-tmpmean
  #     tmp2 <- array(echam$data[echam$names == i,,], c(2,dim(echam$data)[2]/2, dim(echam$data)[3]))
  #     tmp2mean=(apply(tmp2,c(1,3),mean))
  #     tmp2ano=array(NA,dim(echam$data)[2:3])
  #     for (j in 1:dim(echam$data)[3]){
  #       tmp2ano[,j]=tmp2[,,j] #-tmp2mean[,j]
  #     }
  #     echam.ind$ensmean[i,]=tmpano
  #     echam.ind$data[i,,]=tmp2ano
  #   }
  # 
  #   analysis.ind=analysis$'Localized analysis'
  #   l1=array(NA,c(length(indices),dim(analysis$'Localized analysis'$ensmean)[2]))
  #   l2=array(NA,c(length(indices),dim(analysis$'Localized analysis'$data)[2:3]))
  #   analysis.ind$ensmean=l1
  #   analysis.ind$data=l2
  #   analysis.ind$lon=0
  #   analysis.ind$lat=0
  #   analysis.ind$height=0
  #   analysis.ind$names=indices[1:length(indices)]
  #   #analysis.ind=list(l1,l2)
  #   #names(analysis.ind)=c('ensmean','data')
  #   rownames(analysis.ind$ensmean)=indices[1:length(indices)]
  #   rownames(analysis.ind$data)=indices[1:length(indices)]
  #   for (i in indices[1:length(indices)]){
  #     # split into summer and winter for separate means
  #     tmp <- array(analysis$'Localized analysis'$ensmean[analysis$'Localized analysis'$names == i,], c(2,dim(analysis$'Localized analysis'$ensmean)[2]/2))
  #     tmpmean=rowMeans(tmp)
  #     tmpano=tmp #-tmpmean
  #     tmp2 <- array(analysis$'Localized analysis'$data[analysis$'Localized analysis'$names == i,,], c(2,dim(analysis$'Localized analysis'$data)[2]/2, dim(analysis$'Localized analysis'$data)[3]))
  #     tmp2mean=(apply(tmp2,c(1,3),mean))
  #     tmp2ano=array(NA,dim(analysis$'Localized analysis'$data)[2:3])
  #     for (j in 1:dim(analysis$'Localized analysis'$data)[3]){
  #       tmp2ano[,j]=tmp2[,,j] #-tmp2mean[,j]
  #     }
  #     analysis.ind$ensmean[i,]=tmpano
  #     analysis.ind$data[i,,]=tmp2ano
  #   }
  # 
  #   for (i in indices[1:floor(length(indices)/2)]){
  #     plot_echam(analysis.ind, varname=i, type='data', lty=3, add=F)
  #     plot_echam(analysis.ind, varname=i, type='ensmean', lty=3, add=T)
  #     plot_echam(echam.ind, varname=i, type='data', lty=2, add=T)
  #     plot_echam(echam.ind, varname=i, type='ensmean', lty=2, add=T)
  #     plot_echam(validate, varname=i, type='data', lty=1, add=T)
  #   }
  # } else {
  #   for (i in indices[1:(length(indices)/2)]){
  #     plot_echam(echam, varname=i, type='data', lty=2)
  #     plot_echam(analysis.ind[[2]], varname=i, type='data', lty=3, add=T)
  #     plot_echam(echam, varname=i, type='ensmean', lty=2, add=T)
  #     plot_echam(analysis[[2]], varname=i, type='ensmean', lty=3, add=T)
  #     plot_echam(validate, varname=i, type='data', lty=1, add=T)
  #   }
  # }
  # dev.off()
  # #    \caption{Validation time series (solid), ECHAM5 ensemble (thin dashed lines) and ensemble mean (thick dashed lines), and analysis with perfect proxies (thin dotted lines) and analysis mean (thick dotted lines) for different indices (\Sexpr{paste(pnames[1:length(indices)], indices, sep=' ', collapse=', ')}). Blue denotes winter (October to March) and summer (April to September). }
  
  
  # 
  # 
  # ################################################################################
  # # Fig. 7b: validation time series, aggregated
  #
  # #<<label=indices2, echo=FALSE, fig=TRUE, width=8, height=6, results=hide, eps=FALSE>>=
  # pdf(paste(figpath,'/vali_ind_2.pdf',sep='/'), width=9, height=6, paper='special')
  # par(mfrow=c(floor(length(indices)/4), 2), cex.axis=1.4, cex.lab=1.4, mar=c(3,5,1,1), oma=c(0,0,0,0))
  #
  # #for (i in indices[(length(indices):(length(indices)/2)]){
  # for (i in indices[(length(indices)-1):(floor(length(indices)/2))]){
  #   plot_echam(analysis.ind, varname=i, type='data', lty=3, add=F)
  #   plot_echam(echam.ind, varname=i, type='data', lty=2, add=T)
  #   plot_echam(echam.ind, varname=i, type='ensmean', lty=2, add=T)
  #   plot_echam(analysis.ind, varname=i, type='ensmean', lty=3, add=T)
  #   plot_echam(validate, varname=i, type='data', lty=1, add=T)
  # }
  # dev.off()
  # #    \caption{As in figure \ref{fig:indices} but for aggregated temperature, precipitation, mean sea level pressure }
  if (ind_ECHAM) {
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
    
    for (ind in c(1,2,3,4,6)) {
      if (ind == 1) {mainname='NHt2m'}
      if (ind == 2) {mainname='NAM.temp2'}
      if (ind == 3) {mainname='SAM.temp2'}
      if (ind == 4) {mainname='AFR.temp2'}
      if (ind == 6) {mainname='AUS.temp2'}
      for (seas in c(1,2)) {
        if (seas == 1) {color='darkblue'; color3='blue'; color2='cyan'}
        if (seas == 2) {color='darkred'; color3='red'; color2='orange'}
        ymin=min(eind2$min[ind,seas,])
        ymax=max(eind2$max[ind,seas,])
        plot(vind2$data[ind,seas,],ty='l',col=color,lty=1,ylim=c(ymin,ymax),main=mainname)
        lines(eind2$ensmean[ind,seas,],ty='l',col=color2,lwd=2,lty=2,main='')
        lines(eind2$min[ind,seas,],ty='l',col=color2,lty=2,lwd=1,main='')
        lines(eind2$max[ind,seas,],ty='l',col=color2,lty=2,lwd=1,main='')
        lines(aind2$ensmean[ind,seas,],ty='l',col=color3,lwd=2,lty=3,main='')
        lines(aind2$min[ind,seas,],ty='l',col=color3,lty=3,lwd=1,main='')
        lines(aind2$max[ind,seas,],ty='l',col=color3,lty=3,lwd=1,main='')
      }
    }
    dev.off()
    
    ################################################################################
    # Fig. 8:  index correlation
    # new like in the paper
    #<<label=indices_corr, echo=FALSE, fig=TRUE, width=12, height=7, results=hide, eps=FALSE>>=
    
    inds <- c('ENH.temp2', 'NAM.temp2', 'SAM.temp2', 'AFR.temp2')
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
    
    if (pseudoproxy) {
      
      RE.mat <- array(NA, c(length(inds), 2, 30))
      corr.mat <- array(NA, c(length(inds)*2, 2, 30))
      ind.names <- list()
      for (se in 1:2){
        #        if (se == 2){
        #          inds <- c('NH.temp2', 'NEU.temp2', 'NEU.precip', 'HC', 'SJ', 'Z100', 'Z300', 'PWC')
        #        }
        ind.names[[se]] <- inds
        for (k in seq(inds)){
          RE.mat[k,se,1:29] <- RE.ind$data[RE.ind$names == inds[k],se,]
          RE.mat[k,se,30] <- RE.ind$ensmean[RE.ind$names == inds[k],se]
          corr.mat[k*2,se,1:29] <- acorr.ind$data[acorr.ind$names == inds[k],se,]
          corr.mat[k*2,se,30] <- acorr.ind$ensmean[acorr.ind$names == inds[k],se]
          corr.mat[k*2-1,se,1:29] <- ecorr.ind$data[ecorr.ind$names == inds[k],se,]
          corr.mat[k*2-1,se,30] <- ecorr.ind$ensmean[ecorr.ind$names == inds[k],se]
        }
      }
    } else {
      RE.mat <- array(NA, c(length(inds), 2, 31))
      corr.mat <- array(NA, c(length(inds)*2, 2, 31))
      ind.names <- list()
      for (se in 1:2){
        #        if (se == 2){
        #          inds <- c('NH.temp2', 'NEU.temp2', 'NEU.precip', 'HC', 'SJ', 'Z100', 'Z300', 'PWC')
        #        }
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
      #boxplot(t(RE.mat[,se,1:29]), at=seq(inds), col=ind.col[[se]],
      #    add=T, xaxt='n', yaxt='n', xlab='', ylab='', lwd=1, boxwex=0.5, range=0, lty=1)
      if (pseudoproxy){
        boxplot(t(RE.mat[,se,1:29]), at=seq(inds), col=hcl(rep(seq(indices)/length(indices)*360, each=1), c=45, l=75),add=T, xaxt='n', yaxt='n', xlab='', ylab='', lwd=1, boxwex=0.5, range=0, lty=1)
        points(seq(inds) - 0.35, RE.mat[,se,30], pch='>', lwd=3, cex=1.4)
      } else {
        boxplot(t(RE.mat[,se,1:30]), at=seq(inds), col=hcl(rep(seq(indices)/length(indices)*360, each=1), c=45, l=75),add=T, xaxt='n', yaxt='n', xlab='', ylab='', lwd=1, boxwex=0.5, range=0, lty=1)
        points(seq(inds) - 0.35, RE.mat[,se,31], pch='>', lwd=3, cex=1.4)
        
      }
      text(seq(inds), rep(-0.9,3), ind.names[[se]], cex=1.2, adj=c(0.5,1))
      text(0.5,1.1,paste(pnames[se], ' Skill (RE) in ', c('October to March', 'April to September')[se], sep='') , cex=1.4, adj=c(0,1))
    }
    
    #    indind <- c(1,2,4,5,7,8,10,11,13,14,16,17)
    indind <- c(1,2,4,5,7,8,10,11)
    for (se in 1:2){
      plot(0, type='n', xaxt='n', ylim=c(-1.1,1.1), xlim=c(0,max(indind)+1), xlab='', ylab='', yaxt=if (se == 2) 'n' else 's')
      polygon(c(-1,-1,max(indind)+c(2,2)), c(0,-2,-2,0), border=NA, col=grey(0.9))
      #boxplot(t(corr.mat[,se,1:29]), at=indind, col=ind.col2[[se]],
      #    add=T, xaxt='n', yaxt='n', xlab='', ylab='', lwd=1, range=0, lty=1)
      if (pseudoproxy){
        boxplot(t(corr.mat[,se,1:29]), at=indind, col=hcl(rep(seq(indices)/length(indices)*360, each=2), c=c(30,60), l=c(95,50)),
                add=T, xaxt='n', yaxt='n', xlab='', ylab='', lwd=1, range=0, lty=1)
        points(seq(1.5,max(indind),3) - 1.2, corr.mat[seq(1,nrow(corr.mat),2),se,30], pch='>', lwd=3, cex=1.4)
        points(seq(1.5,max(indind),3) + 1.2, corr.mat[seq(2,nrow(corr.mat),2),se,30], pch='<', lwd=3, cex=1.4)
      } else {
        boxplot(t(corr.mat[,se,1:30]), at=indind, col=hcl(rep(seq(indices)/length(indices)*360, each=2), c=c(30,60), l=c(95,50)),
                add=T, xaxt='n', yaxt='n', xlab='', ylab='', lwd=1, range=0, lty=1)
        points(seq(1.5,max(indind),3) - 1.2, corr.mat[seq(1,nrow(corr.mat),2),se,31], pch='>', lwd=3, cex=1.4)
        points(seq(1.5,max(indind),3) + 1.2, corr.mat[seq(2,nrow(corr.mat),2),se,31], pch='<', lwd=3, cex=1.4)
      }
      text(seq(1.5,max(indind),3), rep(-0.9,3), ind.names[[se]], cex=1.2, adj=c(0.5,1))
      text(0,1.1,paste(pnames[se + 2], ' Correlation in ', c('October to March', 'April to September')[se], sep=''), cex=1.4, adj=c(0,1))
    }
    dev.off()
    
    # old version
    
    # ind <- nrow(corr[[4]]$data) - rev(seq(indices)) + 1
    # ind.arr <- array(NA, c(length(indices)*2, 2, dim(echam$data)[3]))
    # ind.arr[seq(1,dim(ind.arr)[1],2),,] <- corr.ech$data[ind,,]
    # ind.arr[seq(2,dim(ind.arr)[1],2),,] <- corr[[4]]$data[ind,,]
    # xind <- sort(c(seq(1, length(indices)*3, 3), seq(2, length(indices)*3, 3)))
    #
    # par(mfrow=c(2,1), oma=c(6,0,1,0), mar=c(1,3,1,2), cex.axis=1.4, cex.lab=1.4, cex.main=1.4)
    # for (i in 1:2) {
    #   plot(0, type='n', xlim=range(xind), ylim=c(-1,1), xaxt='n')
    #   polygon(c(-1,rep(max(xind)+2, 2), -1), c(0,0,-2,-2), border=NA, col=grey(0.8))
    #   box()
    #   boxplot(t(ind.arr[,i,]), at=xind, add=T, xaxt='n', col=hcl(rep(seq(indices)/length(indices)*360, each=2), c=c(30,60), l=c(95,50)))
    #   text(1,1,c('Winter', 'Summer')[i], adj=c(0,1), font=1, cex=1.4)
    # #  points(seq(2.8, length(indices)*3, 3), corr.analogue[[2]]$data[ind,i], pch='<', cex=1.4)
    # }
    # axis(1, at=seq(1.5, length(indices)*3, 3), labels=indices, las=2, tick=F)
    
    #    \caption{Correlation of the indices in the ECHAM ensemble (light boxes) and the localized pseudo proxy analysis (dark boxes) with the respective validation time series.}
    
    
    
    ################################################################################
    # Tab x: RE for indices
    
    #<<label=tabREwinter, echo=FALSE, results=tex>>=
    # require(xtable)
    # RE.ind <- RE.tot$data[nrow(RE.tot$data) - length(indices):1 + 1,1,seq(1,8,2)]
    # 
    # rownames(RE.ind) <- indices
    # colnames(RE.ind) <- pnames[seq(along=analysis)]
    # 
    # print(xtable(RE.ind, caption=paste('RE for different indices in winter (October to March), the columns contain the results for', paste(pnames[seq(along=analysis)], names(analysis), collapse=', '), collapse=' '), label="tab:REwinter", digits=3, table.placement='htbp', caption.placement='top'))
    # 
    # #<<label=tabREsummer, echo=FALSE, results=tex>>=
    # require(xtable)
    # RE.ind <- RE.tot$data[nrow(RE.tot$data) - length(indices):1 + 1,1,seq(2,8,2)]
    # 
    # rownames(RE.ind) <- indices
    # colnames(RE.ind) <- pnames[seq(along=analysis)]
    # 
    # print(xtable(RE.ind, caption=paste('RE for different indices in winter (April to September), the columns contain the results for', paste(pnames[seq(along=analysis)], names(analysis), collapse=', '), collapse=' '), label="tab:REsummer", digits=3, table.placement='htbp', caption.placement='top'))
    
  } # end if(echam_ind)
  
  ################################################################################
  # # Fig. 9a: temperature reconstruction skill 
  # if (!monthly_out){
  # #<<label=Giorgitemp2, echo=FALSE, fig=TRUE, eps=FALSE, width=12, height=14>>=
  # if (pseudoproxy){
  #   giorgi.RE <- compute_avg_RE_pseudoproxy(H.giorgi, echam, analysis[[2]], validate)
  #   giorgi.corr <- compute_avg_corr_pseudoproxy(H.giorgi, echam, analysis[[i]], validate)
  # } else {
  #   #  giorgi.RE <- compute_avg_RE(H.giorgi, echam, analysis[[2]], validate)
  #   #  giorgi.corr <- compute_avg_corr(H.giorgi, echam, analysis[[2]], validate)
  #   giorgi.RE <- compute_avg_RE(H.giorgi, echam, analysis, validate)
  #   giorgi.corr <- compute_avg_corr(H.giorgi, echam, analysis, validate)
  # }
  # varname <- c('temp2', 'precip', 'slp')
  # varn <- 1
  # 
  # pdf(paste(figpath,'RE_giorgi_temp.pdf',sep='/'), width=9, height=6, paper='special')
  # #if (pseudoproxy){
  # par(mfrow=c(4,1), mar=c(1,3,1,1), oma=c(0,0,0,0), cex.axis=1.2, cex.lab=1.2)
  # #} else {
  # #  par(mfrow=c(2,1), mar=c(1,3,1,1), oma=c(0,0,0,0), cex.axis=1.2, cex.lab=1.2)
  # #}
  # ncorr <- length(giorgi.short)
  # RE.mat <- giorgi.RE[ncorr*(varn-1) + 1:ncorr,,]
  # indind <- 1:ncorr
  # for (se in 1:2){
  #   plot(0, type='n', xaxt='n', ylim=c(-1,1), xlim=c(1,max(indind)), xlab='', ylab='')
  #   polygon(c(-3, -3, max(indind)+5, max(indind+5)), c(0,-2,-2,0), border=NA, col=grey(0.8))
  #   if (pseudoproxy){
  #     boxplot(t(RE.mat[,se,1:29]), at=indind, col=hcl((varn-1)*120, c=40, l=50),
  #             add=T, xaxt='n', yaxt='n', xlab='', ylab='', lwd=1, boxwex=0.5)
  #     points(indind - 0.3, RE.mat[,se,30], pch='>', lwd=3, cex=1.2)
  #   } else {
  #     boxplot(t(RE.mat[,se,1:30]), at=indind, col=hcl((varn-1)*120, c=40, l=50),
  #             add=T, xaxt='n', yaxt='n', xlab='', ylab='', lwd=1, boxwex=0.5)
  #     points(indind - 0.3, RE.mat[,se,31], pch='>', lwd=3, cex=1.2)
  #   }
  #   text(indind, rep(-0.8,length(indind)), giorgi.short, cex=1, adj=c(0.5,1.3))
  #   text(0.5,1,c('Winter', 'Summer')[se], cex=1.2, adj=c(0,1))
  #   abline(v=seq(1.5, max(indind), 1), lty=3)
  # }
  # 
  # 
  # ncorr <- length(giorgi.short)*2
  # corr.mat <- giorgi.corr[ncorr*(varn-1) + 1:ncorr,,]
  # indind <- rbind(seq(1, by=4, length=ncorr/2), seq(2, by=4, length=ncorr/2))
  # 
  # for (se in 1:2){
  #   plot(0, type='n', xaxt='n', ylim=c(-1,1), xlim=c(2,max(indind)-1), xlab='', ylab='')
  #   polygon(c(-3, -3, max(indind)+5, max(indind+5)), c(0,-2,-2,0), border=NA, col=grey(0.8))
  #   if (pseudoproxy){
  #     boxplot(t(corr.mat[,se,1:29]), at=as.vector(indind), col=hcl((varn-1)*120, c=40, l=c(90,50)), add=T, xaxt='n', yaxt='n', xlab='', ylab='', lwd=1)
  #     points(apply(indind, 2, mean) - 1.4, corr.mat[seq(1,ncorr,2),se,30], pch='>', lwd=3, cex=1.2)
  #     points(apply(indind, 2, mean) + 1.4, corr.mat[seq(2,ncorr,2),se,30], pch='<', lwd=3, cex=1.2)
  #   } else {
  #     boxplot(t(corr.mat[,se,1:30]), at=as.vector(indind), col=hcl((varn-1)*120, c=40, l=c(90,50)), add=T, xaxt='n', yaxt='n', xlab='', ylab='', lwd=1)
  #     points(apply(indind, 2, mean) - 1.4, corr.mat[seq(1,ncorr,2),se,31], pch='>', lwd=3, cex=1.2)
  #     points(apply(indind, 2, mean) + 1.4, corr.mat[seq(2,ncorr,2),se,31], pch='<', lwd=3, cex=1.2) 
  #   }  
  #   text(apply(indind, 2, mean), rep(-0.8,ncol(indind)), giorgi.short, cex=1, adj=c(0.5,1.3))
  #   text(0.5,1,c('Winter', 'Summer')[se], cex=1.2, adj=c(0,1))
  #   abline(v=seq(3.5,max(indind), by=4), lty=3)
  # }
  # dev.off() 
  # #    \caption{Skill in reconstructing area-average temperature in different subcontinental regions as defined in \citet{Giorgi2000}}
  # }
  
  #  
  # ################################################################################
  # # Fig. 9b: precip reconstruction skill 
  # 
  # #<<label=Giorgiprecip, echo=FALSE, fig=TRUE, eps=FALSE, width=12, height=14>>=
  # giorgi.RE <- compute_avg_RE(H.giorgi, echam, analysis, validate)
  # varname <- c('temp2', 'precip', 'slp')
  # varn <- 2 # 2 for precip
  # 
  # ncorr <- length(giorgi.short)
  # RE.mat <- giorgi.RE[ncorr*(varn-1) + 1:ncorr,,]
  # indind <- 1:ncorr
  # 
  # pdf(paste(figpath,'RE_giorgi_precip.pdf',sep='/'), width=9, height=6, paper='special')
  # par(mfrow=c(4,1), mar=c(1,3,1,1), oma=c(0,0,0,0), cex.axis=1.2, cex.lab=1.2)
  # for (se in 1:2){
  #   plot(0, type='n', xaxt='n', ylim=c(-1,1), xlim=c(1,max(indind)), xlab='', ylab='')
  #   polygon(c(-3, -3, max(indind)+5, max(indind+5)), c(0,-2,-2,0), border=NA, col=grey(0.8))
  #   boxplot(t(RE.mat[,se,1:29]), at=indind, col=hcl((varn-1)*120, c=40, l=50),
  #           add=T, xaxt='n', yaxt='n', xlab='', ylab='', lwd=1, boxwex=0.5)
  #   points(indind - 0.3, RE.mat[,se,30], pch='>', lwd=3, cex=1.2)
  #   text(indind, rep(-0.8,length(indind)), giorgi.short, cex=1, adj=c(0.5,1.3))
  #   text(0.5,1,c('Winter', 'Summer')[se], cex=1.2, adj=c(0,1))
  #   abline(v=seq(1.5, max(indind), 1), lty=3)
  # }
  # ncorr <- length(giorgi.short)*2
  # corr.mat <- giorgi.corr[ncorr*(varn-1) + 1:ncorr,,]
  # indind <- rbind(seq(1, by=4, length=ncorr/2), seq(2, by=4, length=ncorr/2))
  # 
  # for (se in 1:2){
  #   plot(0, type='n', xaxt='n', ylim=c(-1,1), xlim=c(2,max(indind)-1), xlab='', ylab='')
  #   polygon(c(-3, -3, max(indind)+5, max(indind+5)), c(0,-2,-2,0), border=NA, col=grey(0.8))
  #   boxplot(t(corr.mat[,se,1:29]), at=as.vector(indind), col=hcl((varn-1)*120, c=40, l=c(90,50)),
  #           add=T, xaxt='n', yaxt='n', xlab='', ylab='', lwd=1)
  #   points(apply(indind, 2, mean) - 1.4, corr.mat[seq(1,ncorr,2),se,30], pch='>', lwd=3, cex=1.2)
  #   points(apply(indind, 2, mean) + 1.4, corr.mat[seq(2,ncorr,2),se,30], pch='<', lwd=3, cex=1.2)
  #   text(apply(indind, 2, mean), rep(-0.8,ncol(indind)), giorgi.short, cex=1, adj=c(0.5,1.3))
  #   text(0.5,1,c('Winter', 'Summer')[se], cex=1.2, adj=c(0,1))
  #   abline(v=seq(3.5,max(indind), by=4), lty=3)
  # }
  # #    \caption{Skill in reconstructing area-average precipitation in different subcontinental regions as defined in \citet{Giorgi2000}}
  # 
  # dev.off()
  # 
  
  
  
  
  ################################################################################
  # plot distance weights
  if (plot_dweights) {
    library(raster)
    pdf(paste(figpath,'distance_weights.pdf',sep='/'), width=4.5, height=4.5, 
        paper='special')
    image(d.weights_all)
    dev.off()
  }
  
  
  
  
  
  
  
  
  
  
  
  
  #   if ((sixmonstatevector) & (instrumental) & (cyr == syr2)) {
  #     echamatprox.arr <- array(NA,c(length(proxies$lon), length(proxies$time)))
  #     echamatprox.arr.allts <- array(NA,c(length(proxies.allts$lon), 
  #                                         length(proxies.allts$time)))
  #     for(i in 1:(length(proxies$lon)/6)){
  #       plon <- proxies$lon[i]
  #       plat <- proxies$lat[i]
  #       clon <- echam$lon
  #       clat <- echam$lat
  #       k=which(abs(clon-plon+0.001)==min(abs(clon-plon+0.001))) 
  #       # +0.001 to avoid find 2 locations with same distance
  #       l=which(abs(clat-plat)==min(abs(clat-plat)))
  #       if (max(match(k,l,nomatch=-99999))==-99999) {
  #         k=which(abs(clon-plon+2.001)==min(abs(clon-plon+2.001)))
  #         l=which(abs(clat-plat)==min(abs(clat-plat)))
  #       }
  #       if (max(match(k,l,nomatch=-99999))==-99999) {
  #         k=which(abs(clon-plon-2.001)==min(abs(clon-plon-2.001)))
  #         l=which(abs(clat-plat)==min(abs(clat-plat)))
  #       }
  #       if (max(match(k,l,nomatch=-99999))==-99999) {
  #         k=which(abs(clon-plon+0.001)==min(abs(clon-plon+0.001)))
  #         l=which(abs(clat-plat+2)==min(abs(clat-plat+2)))
  #       }
  #       if (max(match(k,l,nomatch=-99999))==-99999) {
  #         k=which(abs(clon-plon+0.001)==min(abs(clon-plon+0.001)))
  #         l=which(abs(clat-plat-2)==min(abs(clat-plat-2)))
  #       }
  #       if (max(match(k,l,nomatch=-99999))==-99999) {
  #         k=which(abs(clon-plon+2.001)==min(abs(clon-plon+2.001)))
  #         l=which(abs(clat-plat+2)==min(abs(clat-plat+2)))
  #       }
  #       if (max(match(k,l,nomatch=-99999))==-99999) {
  #         k=which(abs(clon-plon-2.001)==min(abs(clon-plon-2.001)))
  #         l=which(abs(clat-plat-2)==min(abs(clat-plat-2)))
  #       }
  #       if (max(match(k,l,nomatch=-99999))==-99999) {
  #         k=which(abs(clon-plon-2.001)==min(abs(clon-plon-2.001)))
  #         l=which(abs(clat-plat+2)==min(abs(clat-plat+2)))
  #       }
  #       if (max(match(k,l,nomatch=-99999))==-99999) {
  #         k=which(abs(clon-plon+2.001)==min(abs(clon-plon+2.001)))
  #         l=which(abs(clat-plat-2)==min(abs(clat-plat-2)))
  #       }
  #       m=k[which(match(k,l)>0)]
  #       if (length(m) > 0) {
  #         if (proxies$names[i]=="precip") { 
  #           m <- m+length(echam$lon) 
  #         } else if (proxies$names[i]=="slp") { 
  #           m <- m+(2*length(echam$lon)) 
  #         }
  #         echamatprox.arr[c(i,i+(length(proxies$lon)/6),i+2*(length(proxies$lon)/6),
  #                 i+3*(length(proxies$lon)/6),i+4*(length(proxies$lon)/6),
  #                 i+5*(length(proxies$lon)/6)),] <- echam$ensmean[c(m,
  #                 m+(dim(echam$ensmean)[1]/6),m+2*(dim(echam$ensmean)[1]/6),
  #                 m+3*(dim(echam$ensmean)[1]/6),m+4*(dim(echam$ensmean)[1]/6),
  #                 m+5*(dim(echam$ensmean)[1]/6)),]
  #         echamatprox.arr.allts[c(i,i+(length(proxies$lon)/6),i+2*(length(proxies$lon)
  #                 /6),i+3*(length(proxies$lon)/6),i+4*(length(proxies$lon)/6),
  #                 i+5*(length(proxies$lon)/6)),] <- echam.allts$ensmean[c(m,
  #                 m+(dim(echam$ensmean)[1]/6),m+2*(dim(echam$ensmean)[1]/6),
  #                 m+3*(dim(echam$ensmean)[1]/6),m+4*(dim(echam$ensmean)[1]/6),
  #                 m+5*(dim(echam$ensmean)[1]/6)),]
  #         #        echamatprox.arr[i,which(is.na(proxies$data[i,]))] <- NA
  #         #        echamatprox.arr.allts[i,which(is.na(proxies$data[i,]))] <- NA
  #       } else {
  #         echamatprox.arr[i,] <- NA
  #         echamatprox.arr.allts[i,] <- NA
  #         proxies$data[i,] <- NA
  #         proxies.allts$data[i,] <- NA
  #         #        proxies.anom$data[i,] <- NA
  #       }
  #     }
  #     pdims <- c(prod(dim(proxies.allts$data)[1],s),dim(proxies.allts$data)[2]/s)
  #     prox.echam.corr <- array(diag(cor(t(array(proxies.allts$data,pdims)), 
  #          t(array(echamatprox.arr.allts, pdims)), use='pairwise.complete.obs')), 
  #          c(dim(proxies$data)[1],s)) 
  #     prox.echam.bias <- array(apply(array(proxies.allts$data,pdims),1,mean, na.rm=T)-
  #          (apply(array(echamatprox.arr.allts,pdims),1,mean, na.rm=T)),
  #          c(dim(proxies$data)[1],s)) 
  #   }
  #   
  #   if ((!sixmonstatevector) & (!real_proxies) & (cyr == syr2)) { 
  #     # calculated at first time step for entire series
  #     echamatprox.arr <- array(NA,c(length(proxies$lon), length(proxies$time)))
  #     echamatprox.arr.allts <- array(NA,c(length(proxies.allts$lon), 
  #                                         length(proxies.allts$time)))
  #     for(i in 1:length(proxies$lon)){
  #       plon <- proxies$lon[i]
  #       plat <- proxies$lat[i]
  #       clon <- echam$lon
  #       clat <- echam$lat
  #       k=which(abs(clon-plon+0.001)==min(abs(clon-plon+0.001))) 
  #       # +0.001 to avoid find 2 locations with same distance
  #       l=which(abs(clat-plat)==min(abs(clat-plat))) 
  #       if (max(match(k,l,nomatch=-99999))==-99999) {
  #         k=which(abs(clon-plon+2.001)==min(abs(clon-plon+2.001)))
  #         l=which(abs(clat-plat)==min(abs(clat-plat)))
  #       }
  #       if (max(match(k,l,nomatch=-99999))==-99999) {
  #         k=which(abs(clon-plon-2.001)==min(abs(clon-plon-2.001)))
  #         l=which(abs(clat-plat)==min(abs(clat-plat)))
  #       }
  #       if (max(match(k,l,nomatch=-99999))==-99999) {
  #         k=which(abs(clon-plon+0.001)==min(abs(clon-plon+0.001)))
  #         l=which(abs(clat-plat+2)==min(abs(clat-plat+2)))
  #       }
  #       if (max(match(k,l,nomatch=-99999))==-99999) {
  #         k=which(abs(clon-plon+0.001)==min(abs(clon-plon+0.001)))
  #         l=which(abs(clat-plat-2)==min(abs(clat-plat-2)))
  #       }
  #       if (max(match(k,l,nomatch=-99999))==-99999) {
  #         k=which(abs(clon-plon+2.001)==min(abs(clon-plon+2.001)))
  #         l=which(abs(clat-plat+2)==min(abs(clat-plat+2)))
  #       }
  #       if (max(match(k,l,nomatch=-99999))==-99999) {
  #         k=which(abs(clon-plon-2.001)==min(abs(clon-plon-2.001)))
  #         l=which(abs(clat-plat-2)==min(abs(clat-plat-2)))
  #       }
  #       if (max(match(k,l,nomatch=-99999))==-99999) {
  #         k=which(abs(clon-plon-2.001)==min(abs(clon-plon-2.001)))
  #         l=which(abs(clat-plat+2)==min(abs(clat-plat+2)))
  #       }
  #       if (max(match(k,l,nomatch=-99999))==-99999) {
  #         k=which(abs(clon-plon+2.001)==min(abs(clon-plon+2.001)))
  #         l=which(abs(clat-plat-2)==min(abs(clat-plat-2)))
  #       }
  #       m=k[which(match(k,l)>0)]
  #       if (length(m) > 0) {
  #         if (proxies$names[i]=="temp") {
  #           echamatprox.arr[i,] <- echam$ensmean[m,]
  #           echamatprox.arr.allts[i,] <- echam.allts$ensmean[m,]
  #         } else if (proxies$names[i]=="precip") {  
  #           echamatprox.arr[i,] <- echam$ensmean[(m+length(echam$lon)),]
  #           echamatprox.arr.allts[i,] <- echam.allts$ensmean[(m+length(echam$lon)),]
  #         } else if (proxies$names[i]=="slp") {  
  #           echamatprox.arr[i,] <- echam$ensmean[(m+(2*length(echam$lon))),]
  #           echamatprox.arr.allts[i,] <- echam.allts$ensmean[(m+(2*length(echam$lon))),]
  #         }
  # #        echamatprox.arr[i,which(is.na(proxies$data[i,]))] <- NA
  # #        echamatprox.arr.allts[i,which(is.na(proxies$data[i,]))] <- NA
  #       } else {
  #         echamatprox.arr[i,] <- NA
  #         echamatprox.arr.allts[i,] <- NA
  #         proxies$data[i,] <- NA
  #         proxies.allts$data[i,] <- NA
  # #        proxies.anom$data[i,] <- NA
  #       }
  #     }
  # #    pdims <- c(prod(dim(pdata.anom)[1:2]), dim(pdata.anom)[3])
  #     pdims <- c(prod(dim(proxies.allts$data)[1],s),dim(proxies.allts$data)[2]/s)
  #     prox.echam.corr <- array(diag(cor(t(array(proxies.allts$data,pdims)), 
  #           t(array(echamatprox.arr.allts, pdims)), use='pairwise.complete.obs')), 
  #           c(dim(proxies$data)[1],s)) 
  #     prox.echam.bias <- array(apply(array(proxies.allts$data,pdims),1,mean, na.rm=T)-
  #           (apply(array(echamatprox.arr.allts,pdims),1,mean, na.rm=T)),
  #           c(dim(proxies$data)[1],s)) 
  # #   } else {}
  #   }
  
  # ################################################################################
  # # Fig. 10: Continuous Ranked Probability Score (CRPS)
  # pdf(paste(figpath,'CRPS_new.pdf',sep='/'), width=9, height=6, paper='special')
  
  
  ### plot CRPS score difference ech-ana (positive values are good) with plot_echam: Temp
  if(CRPS){
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
    
  }
  
}
