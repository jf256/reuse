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

figpath=paste0('../figures/',expname,'_',syr,'-',eyr,'/') #format(Sys.time(), "%Y%m%d_%H-%M_")
dir.create(figpath)

pnames <- paste(c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k'), ')', sep='')
data.dim <- dim(echam$data)[c(1,2,2,3)]
data.dim[2:3] <- c(s,data.dim[3]/s)
ens.dim <- c(nrow(echam$ensmean), s, ncol(echam$ensmean)/s)

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
stat.pos1 <- which(calibrate$lon>105&calibrate$lon<120&calibrate$lat>70&calibrate$lat<80) #sibiria

for (i in 1:length(stat.pos1)){
  if (monthly_out) {
    pdf(paste(figpath,'/example_timeseries_scandinavia_mon.pdf',sep=''), 
        width=4.5, height=6, paper='special')
  } else {
    pdf(paste(figpath,'/example_timeseries_sibiria',i,'.pdf',sep=''), 
        width=7, height=6, paper='special')
  }
  par(oma=c(0,0,0,0),mar=c(2,4,2,0.5),mfrow=c(3,1))
  stat.pos<-stat.pos1[i]  
  pos <- getgridboxnum(calibrate,echam)[stat.pos]
  
  
  par(oma=c(0,0,0,0),mar=c(2,4,2,0.5),mfrow=c(3,1))
  
  
  # stat.pos <- which(calibrate$lon>5&calibrate$lon<30&calibrate$lat>56&calibrate$lat<70) #scandinavia
  if (monthly_out) {
    #3yr monthly temp
    pos2 <- seq(1,36,1)
    period <- as.Date(paste(rep(seq(from=1950,to=1952),each=12),rep(seq(1,12),3),rep(15,36),sep='-'))
  } else {
    #30yr summer temp
    pos2 <- seq(2,ncol(validate$data),2)
    period <- validate$time[pos2]
  }
  plot(period,validate$data[pos,pos2],ylim=c(min(cbind(echam$data[pos,pos2,],validate$data[pos,pos2])),
                                             max(cbind(echam$data[pos,pos2,],validate$data[pos,pos2]))),ty='l',col=rgb(0,0,10,10,maxColorValue=10),
       main="Temperature",xlab='',ylab='ºC',xaxt='n')
  polygon(c(period,rev(period)),c(apply(echam$data[pos,pos2,],1,max),
                                  rev(apply(echam$data[pos,pos2,],1,min))),density=NA,
          col=rgb(1,1,1,3,maxColorValue=10))
  lines(period,echam$ensmean[pos,pos2],col=rgb(0,0,0,10,maxColorValue=10))
  polygon(c(period,rev(period)),c(apply(analysis$data[pos,pos2,],1,max),
                                  rev(apply(analysis$data[pos,pos2,],1,min))),density=NA,
          col=rgb(10,0,0,3,maxColorValue=10))
  lines(period,analysis$ensmean[pos,pos2],col=rgb(10,0,0,10,maxColorValue=10))
  lines(period,validate$data[pos,pos2],ylim=c(min(echam$data[pos,pos2,]),
                                              max(echam$data[pos,pos2,])),col=rgb(0,0,10,10,maxColorValue=10))
  lines(period,calibrate$data[stat.pos,pos2]-3,col="yellow")
  
  # plot(stat$time[seq(2,48,2)],stat$data[pos,][201,seq(2,48,2)],ylim=c(min(stat$data[201,pos2]),
  #   max(stat$data[201,pos2])),ty='l',col=rgb(0,0,10,10,maxColorValue=10),
  #   main="Temperature",xlab='',ylab='ºC',xaxt='n')
  # lines(stat$time[seq(2,48,2)],stat$data[pos,][216,seq(2,48,2)],col=rgb(10,0,0,10,maxColorValue=10))
  # lines(stat$time[seq(2,48,2)],stat$data[pos,][219,seq(2,48,2)],col=rgb(0,10,0,10,maxColorValue=10))
  # lines(analysis.anom$time[seq(2,48,2)],analysis.anom$ensmean[676,seq(2,48,2)],
  #       col=rgb(10,0,0,10,maxColorValue=10),ty='l')
#  legend("topleft", c('Instrumental CRU TS3', 'CCC400', "EFK400"),col=c("blue", "black", "red"), 
#         lty=c(1,1,1),lwd=c(1,1,1),pt.cex=1, pt.lwd=1,inset=0.005, bg='white', 
#         box.col='white', cex=1)
  #precip
  plot(period,validate$data[4608+pos,pos2],ylim=c(min(echam$data[4608+pos,pos2,]),
    max(echam$data[4608+pos,pos2,])),ty='l',col=rgb(0,0,10,10,maxColorValue=10),
    main="Precipitation",xlab='',ylab='mm',xaxt='n')
  polygon(c(period,rev(period)),c(apply(echam$data[4608+pos,pos2,],1,max),
    rev(apply(echam$data[4608+pos,pos2,],1,min))),density=NA,
    col=rgb(1,1,1,3,maxColorValue=10))
  lines(period,echam$ensmean[4608+pos,pos2],col=rgb(0,0,0,10,maxColorValue=10))
  polygon(c(period,rev(period)),c(apply(analysis$data[4608+pos,pos2,],1,max),
    rev(apply(analysis$data[4608+pos,pos2,],1,min))),density=NA,
    col=rgb(10,0,0,3,maxColorValue=10))
  lines(period,analysis$ensmean[4608+pos,pos2],col=rgb(10,0,0,10,maxColorValue=10))
  lines(period,validate$data[4608+pos,pos2],ylim=c(min(echam$data[pos,pos2,]),
    max(echam$data[4608+pos,pos2,])),col=rgb(0,0,10,10,maxColorValue=10))
  lines(period,calibrate$data[stat.pos,pos2]+74,col="yellow")
  # lines(period,apply(calibrate$data[stat.pos,pos2],2,mean)+74,col="yellow")
  legend("topleft", c('CRU TS3', 'CCC400', "EFK400"),col=c("blue", "black", "red"), 
         lty=c(1,1,1),lwd=c(1,1,1),pt.cex=1, pt.lwd=1,inset=0.005, bg='white', 
         box.col='white', cex=1)
  #slp
  plot(period,validate$data[2*4608+pos,pos2],ylim=c(min(echam$data[2*4608+pos,pos2,]),
    max(echam$data[2*4608+pos,pos2,])),ty='l',col=rgb(0,0,10,10,maxColorValue=10),
    main="SLP",xlab='',ylab='hPa')
  polygon(c(period,rev(period)),c(apply(echam$data[2*4608+pos,pos2,],1,max),
    rev(apply(echam$data[2*4608+pos,pos2,],1,min))),density=NA,
    col=rgb(1,1,1,3,maxColorValue=10))
  lines(period,echam$ensmean[2*4608+pos,pos2],col=rgb(0,0,0,10,maxColorValue=10))
  polygon(c(period,rev(period)),c(apply(analysis$data[2*4608+pos,pos2,],1,max),
    rev(apply(analysis$data[2*4608+pos,pos2,],1,min))),density=NA,
    col=rgb(10,0,0,3,maxColorValue=10))
  lines(period,analysis$ensmean[2*4608+pos,pos2],col=rgb(10,0,0,10,maxColorValue=10))
  lines(period,validate$data[2*4608+pos,pos2],ylim=c(min(echam$data[pos,pos2,]),
    max(echam$data[2*4608+pos,pos2,])),col=rgb(0,0,10,10,maxColorValue=10))
  # legend("bottomleft", c('Instrumental CRU TS3', 'CCC400', "EFK400"),col=c("blue", "black", "red"), 
  #   lty=c(1,1,1),lwd=c(1,1,1),pt.cex=1, pt.lwd=1,inset=0.005, bg='white', 
  #   box.col='white', cex=1)
dev.off()
}

# plot sample time series for lon/lat:
paste(validate$lon[278],validate$lat[278])
pdf(paste(figpath,'/example_timeseries_greenland.pdf',sep=''), width=4.5, height=6, paper='special')
par(oma=c(0,0,0,0),mar=c(2,4,2,0.5),mfrow=c(3,1))
#50yr summer temp
plot(validate$data[278,pos2],ylim=c(-15,-12),ty='l',col="blue",main="Temperature",
     xlab='',ylab='ºC',xaxt='n')
lines(echam$ensmean[278,pos2],col="black")
lines(analysis$ensmean[278,pos2],col="red")
#precip
plot(validate$data[4608+278,pos2],ylim=c(11,14),ty='l',col="blue",main="Precipitation",
     xlab='',ylab='mm',xaxt='n')
lines(echam$ensmean[4608+278,pos2],col="black")
lines(analysis$ensmean[4608+278,pos2],col="red")
#slp
plot(validate$time[pos2],validate$data[2*4608+278,pos2],
     ylim=c(1013,1021),ty='l',col="blue",main="SLP",xlab='Year',ylab='hPA')
lines(echam$time[pos2],echam$ensmean[2*4608+278,pos2],col="black")
lines(analysis$time[pos2],analysis$ensmean[2*4608+278,pos2],col="red")
legend("bottomleft", c('Instrumental CRU TS3', 'CCC400', "EFK400"),col=c("blue", "black", "red"),
       lty=c(1,1,1),lwd=c(1,1,1),pt.cex=1, pt.lwd=1,inset=0.005, bg='white',
       box.col='white', cex=1)
dev.off()

 
#chose timestep for sample year plots
t=2

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




################################################################################
# Fig. 0: station locations 
#
if (countseries) {
  # count number of assimilation data of each type and add to plot
  nmxd <- ntrw <- nprox <- ndoc <-ninstslp <-ninsttemp <- ninst <- rep(NA,length(1603:2004))
  i <- 1
  for (cyr in 1603:2004) {
    print(cyr)
    load(paste0(datadir,"EnSRF_analysis/prepplot_v3_seasonal/analysis_",cyr,".Rdata"))
    ninst[i] <- length(which(calibrate$sour=="inst"))/6
    ninsttemp[i] <- length(which(calibrate$sour=="inst"&calibrate$names=="temp2"))/6
    ninstslp[i] <- length(which(calibrate$sour=="inst"&calibrate$names=="slp"))/6
    ndoc[i] <- length(which(calibrate$sour=="doc"))/6
    nprox[i] <- length(which(calibrate$sour=="prox"))/6
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
  nrecords <- cbind((1603:2004),ninst,ninsttemp,ninstslp,ndoc,nprox,ntrw,nmxd,rep(0,length(1603:2004)))
  nrecords[41,] <- nrecords[40,] <- nrecords[39,] <- nrecords[38,] # check input data 1639-41!
  pdf(paste(figpath,'record_no.pdf',sep='/'), width=9, height=3.5, paper='special')
  par(oma=c(0,0,0,0),mar=c(4,4,0.4,0.4))
  plot(nrecords[,1],nrecords[,9],ty="l",col='white',ylim=c(0,350),xlab="year",ylab="No of records")
  polygon(c(nrecords[,1],rev(nrecords[,1])),c(nrecords[,7],rev(nrecords[,9])),col="brown")
  polygon(c(nrecords[,1],rev(nrecords[,1])),c((nrecords[,7]+nrecords[,8]),rev(nrecords[,7])),col="green")
  polygon(c(nrecords[,1],rev(nrecords[,1])),
          c((nrecords[,5]+nrecords[,7]+nrecords[,8]),rev((nrecords[,7]+nrecords[,8]))),col="red")
  polygon(c(nrecords[,1],rev(nrecords[,1])),
          c((nrecords[,3]+nrecords[,5]+nrecords[,7]+nrecords[,8]),
            rev((nrecords[,5]+nrecords[,7]+nrecords[,8]))),col="blue")
  polygon(c(nrecords[,1],rev(nrecords[,1])),
          c((nrecords[,4]+nrecords[,3]+nrecords[,5]+nrecords[,7]+nrecords[,8]),
            rev((nrecords[,3]+nrecords[,5]+nrecords[,7]+nrecords[,8]))),col="orange")
  legend("topleft", c('TRW', 'MXD', "Docum. temp.", "Instr. temp.", "Instr. SLP"), 
         pch=rep(15,5), col=c("brown", "green", "red", "blue", "orange"), pt.cex=1, pt.lwd=1, 
         inset=0.005, bg='white', box.col='white', cex=1)
  dev.off()
  
  statyr=1804
  #statsyr=1790
  #stateyr=1804
  load(paste0(datadir,"EnSRF_analysis/prepplot_v3_seasonal/analysis_",statsyr,".Rdata"))
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
  pdf(paste0(figpath,'stat_locations_',statyr,'.pdf'), width=9, height=3.5, paper='special')
  #layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
  layout(matrix(c(1,2), 1, 2, byrow = TRUE))
  par(oma=c(0,0,0,0))
  levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0.05,0.2,0.4,0.6,0.8,1)
  plot_echam(dat,type='data',cex.pt=1.5,names=pnames[1:dim(dat$data)[3]],lev=levs,colorbar=F,
             stations=calibrate,st.col=NULL,statpanel=c(1,2),add=T)
  dev.off()
}



###############################################################################
# Fig. xx: spread-error ratio analysis
plotdata=echam
#plotdata$ensmean <- cbind(ech_spr_err_ratio,ana_spr_err_ratio)
#plotdata$ensmean <- cbind(sprerr.win,sprerr.sum)
plotdata$data <- array(cbind(sprerr.win,sprerr.sum,sprerr.c.win,sprerr.c.sum), 
                   c(length(sprerr.c.win), 1, 4))
plotdata$names <- plotdata$names[echam$names=="temp2"]
plotdata$lon <- plotdata$lon[echam$names=="temp2"]
plotdata$lat <- plotdata$lat[echam$names=="temp2"]
# v3: corrected for uncertainties in instrumental data
# v4: only 1901-1980 because of validation data error after 1980
# v5: obs error not taken into account (top), taken in account (bottom)

pdf(paste0(figpath,'/spread_error_ratio_anom_tmean.pdf'), width=9, height=7, paper='special')
layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE), height=c(3,3,1))
#layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
par(oma=c(0,2,3,0))
#levs <- c(0,0.33,0.5,0.57,0.67,0.77,0.83,0.91,1.1,1.2,1.3,1.5,1.75,2,3,Inf)
levs <- c(0,0.33,0.5,0.67,0.83,1.2,1.5,2,3,Inf)
plot_echam(plotdata, varname='temp2', type='data', cex.pt=1.5, names=pnames[1:dim(plotdata$data)[3]], 
           lev=levs, st.col=NULL, stations=NULL, add=T, #ti=(1:2),
           colnames=c("Oct.-Apr.","May-Sep."),rownames=c("w/o obs. error","w/ obs. error"))
dev.off()


# Fig. xx: spread-error ratio echam
plotdata=echam
plotdata$data <- array(cbind(ech.sprerr.win,ech.sprerr.sum,ech.sprerr.c.win,ech.sprerr.c.sum), 
                       c(length(ech.sprerr.c.win), 1, 4))
plotdata$names <- plotdata$names[echam$names=="temp2"]
plotdata$lon <- plotdata$lon[echam$names=="temp2"]
plotdata$lat <- plotdata$lat[echam$names=="temp2"]
# v3: corrected for uncertainties in instrumental data
# v4: only 1901-1980 because of validation data error after 1980
# v5: obs error not taken into account (top), taken in account (bottom)
pdf(paste0(figpath,'/spread_error_ratio_echam_anom_tmean.pdf'), width=9, height=7, paper='special')
layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE), height=c(3,3,1))
#layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
par(oma=c(0,2,3,0))
#levs <- c(0,0.33,0.5,0.57,0.67,0.77,0.83,0.91,1.1,1.2,1.3,1.5,1.75,2,3,Inf)
levs <- c(0,0.33,0.5,0.67,0.83,1.2,1.5,2,3,Inf)
plot_echam(plotdata, varname='temp2', type='data', cex.pt=1.5, names=pnames[1:dim(plotdata$data)[3]], 
           lev=levs, st.col=NULL, stations=NULL, add=T, #ti=(1:2),
           colnames=c("Oct.-Apr.","May-Sep."),rownames=c("w/o obs. error","w/ obs. error"))
dev.off()

# plotdata$data <- array(mes_obserr,c(nrow(mes_obserr),1,2))
# pdf(paste(figpath,'/spread_error_ratio_obserr.pdf',sep=''), width=9, height=4.5, paper='special')
# layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
# par(oma=c(0,0,0,0))
# levs <- c(-Inf, seq(0,2,0.25), Inf)
# plot_echam(plotdata, varname='temp2', type='data', cex.pt=1.5, names=pnames[1:dim(plotdata$data)[3]], lev=levs, st.col=NULL, stations=NULL, add=T)
# dev.off()


################################################################################
# Fig. xx: Talagrant diagram
if (vali) {
  if (!recon_vali) {
    if (anomaly_assim) {
      ereliable <- ereliable.anom
      areliable <- ereliable.anom
      erel_obserr <- erel_obserr.anom
      arel_obserr <- erel_obserr.anom
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
espread <- echam
espread$data <- array(ech.spread, c(nrow(ech.spread), 1, ncol(ech.spread)))
aspread <- echam
if (pseudoproxy) { 
  ana.spread.bkp <- ana.spread
  ana.spread <- ana.spread.bkp[4] 
}
#if ((instrumental) || (inst_at_proxy) || (real_proxies)) {ana.spread <- ana.spread.bkp[2] }
#aspread$data <- array(unlist(ana.spread)/as.vector(ech.spread)*100, c(nrow(ech.spread), 1, ncol(ech.spread)*(length(ana.spread))))
aspread$data <- array(unlist(ana.spread)/as.vector(ech.spread)*100, c(nrow(ech.spread), 1, ncol(ech.spread)))

#pdf('figures/inst/spread.pdf', width=9, height=6, paper='special')
pdf(paste(figpath,'spread_temp.pdf',sep='/'), width=9, height=6, paper='special')
oldpar <- par(no.readonly=TRUE)
#layout(matrix(c(1:3,3,3+ seq(1,length(ana.spread)*2), rep(4,2) + length(ana.spread)*2),length(ana.spread)+3, 2, byrow=T), height=c(5,lcm(2), rep(5, length(ana.spread)),lcm(2)))
#layout(matrix(c(1:3,3,3+ seq(1,2), rep(4,2) + 2),3, 2, byrow=T), height=c(5,lcm(2), rep(5, 1),lcm(2)))
layout(matrix(c(1,2,3,3,4,5,6,6), 4, 2, byrow = TRUE), height=c(3,1,3,1))
par(oma=c(0,0,0,0))
plot_echam(espread, varname='temp2', names=pnames[1:dim(espread$data)[3]], symmetric=FALSE, cex.pt=1.3, add=TRUE)
plot_echam(aspread, varname='temp2', names=pnames[dim(espread$data)[3] + 1:dim(aspread$data)[3]], lev=seq(0,100,10), cex.pt=1.3, add=TRUE, st.col=NULL, stations=plstat) #calibrate)
#par(oldpar)
#plot_echam(espread, varname='temp2', names=pnames[1:dim(espread$data)[3]], symmetric=FALSE, cex.pt=1.3, add=TRUE)
#plot_echam(aspread, varname='temp2', names=pnames[dim(espread$data)[3] + 1:dim(aspread$data)[3]], lev=seq(0,100,10), cex.pt=1.3, add=TRUE, st.col=1, stations=calibrate) #calibrate)
#par(oldpar)
# \caption{Average temperature spread of the ECHAM ensemble in a) winter (October to March) and b) summer (April to September). Percentage of spread in the analysis ensembles with respect to the ECHAM ensemble for the EnSRF analysis with perfect proxies in c) and d), with perfect proxies and localization in e) and f), with pseudoproxies in g) and h), and with pseudoproxies and localization in i) and j). }
dev.off()

pdf(paste(figpath,'spread_temp_2.pdf',sep='/'), width=9, height=3, paper='special')
oldpar <- par(no.readonly=TRUE)
layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
par(oma=c(0,0,0,0))
plot_echam(aspread, varname='temp2', names=pnames[dim(espread$data)[3] + 1:dim(aspread$data)[3]], lev=seq(0,100,10), cex.pt=1.3, add=TRUE, st.col=NULL, stations=plstat)
dev.off()

################################################################################
# Fig. 2b: maps of precipitation spread
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


################################################################################
# Fig. 2c: maps of SLP spread
pdf(paste(figpath,'spread_slp.pdf',sep='/'), width=9, height=6, paper='special')
#<<label=spread_slp, echo=FALSE, fig=TRUE, width=8, height=11, results=hide, eps=FALSE>>=
oldpar <- par(no.readonly=TRUE)
layout(matrix(c(1,2,3,3,4,5,6,6), 4, 2, byrow = TRUE), height=c(3,1,3,1))
par(oma=c(0,0,0,0))
plot_echam(espread, varname='slp', names=pnames[1:dim(espread$data)[3]], symmetric=FALSE, cex.pt=1.3, add=TRUE)
plot_echam(aspread, varname='slp', names=pnames[dim(espread$data)[3] + 1:dim(aspread$data)[3]], lev=seq(0,100,10), cex.pt=1.3, add=TRUE, st.col=NULL, stations=plstat)
#oldpar <- par(no.readonly=TRUE)
#layout(matrix(c(1:3,3,3+ seq(1,length(ana.spread)*2), rep(4,2) + length(ana.spread)*2),
#              length(ana.spread)+3, 2, byrow=T), height=c(5,lcm(2), rep(5, length(ana.spread)),lcm(2)))
#par(oma=c(0,0,0,0))
#plot_echam(espread, varname='slp', names=pnames[1:dim(espread$data)[3]], symmetric=FALSE, cex.pt=1.3, add=TRUE)
#plot_echam(aspread, varname='slp', cols=rbfun(10), names=pnames[dim(espread$data)[3] + 1:dim(aspread$data)[3]], levs=seq(0,100,10), cex.pt=1.3, add=TRUE, st.col=1, stations=plstat) #calibrate)
#par(oldpar)
dev.off()
#    \caption{According to figure \ref{fig:spread} but for mean sea level pressure in winter (left) and summer (right). }



################################################################################
# Fig. 3a: average temperature update

#<<label=update, echo=FALSE, fig=TRUE, width=8, height=9, results=hide, eps=FALSE>>=
#pdf(paste(figpath,'avg_upd_temp.pdf',sep='/'), width=4, height=9, paper='special')
pdf(paste(figpath,'avg_upd_temp.pdf',sep='/'), width=4, height=4.5, paper='special')
#upd <- lapply(analysis, function(x) apply(array(x$data - echam$data, data.dim), 1:2, mean, na.rm=T))
upd <- apply(array(analysis$data - echam$data, data.dim), 1:2, mean, na.rm=T)

update <- echam
#update$data <- array(unlist(upd), c(nrow(echam$data), 1, ncol(upd[[1]])*length(upd)))
update$data <- array(upd,c(dim(upd)[1],1,dim(upd)[2]))

plot_echam(update, symmetric=T, names=pnames[1:dim(update$data)[3]], cex.pt=1.3, st.col=NULL, stations=plstat)
dev.off()
#    \caption{Average update of the echam ensemble members according to EnSRF with perfect proxies in a) and b), localized perfect proxies in c) and d), pseudoproxies in e) and f), and localized pseudoproxies in g) and h). Results for winter (October to March) are shown in the left column, results for summer (April to October) in the right column. In case of instrumental data there are no pseudoproxies: a) NO localization winter, b) NO localization summer, c) WITH localization winter, d WITH localization summer}



################################################################################
# Fig. 3b: average precipitation update

pdf(paste(figpath,'avg_upd_precip.pdf',sep='/'), width=4, height=4.5, paper='special')
#<<label=update_precip, echo=FALSE, fig=TRUE, width=8, height=9, results=hide, eps=FALSE>>=
plot_echam(update, varname='precip', symmetric=T, names=pnames[1:dim(update$data)[3]], cex.pt=1.3, st.col=NULL, stations=plstat)
dev.off()
#    \caption{According to figure \ref{fig:update} but for precipitation in winter (left) and summer (right)}



################################################################################
# Fig. 3c: average SLP update

#<<update_slp, echo=FALSE, fig=TRUE, width=8, height=9, results=hide, eps=FALSE>>=
pdf(paste(figpath,'avg_upd_slp.pdf',sep='/'), width=4, height=4.5, paper='special')
plot_echam(update, varname='slp', symmetric=T, names=pnames[1:dim(update$data)[3]], cex.pt=1.3, st.col=NULL, stations=plstat)
dev.off()
#    \caption{According to figure \ref{fig:update} but for mean sea level pressure in winter (left) and summer (right)}



################################################################################
# Fig. 4a: average temperature corr echam and analysis vs validation wrt ens. mean 
#<<label=rmse, echo=FALSE, fig=TRUE, width=8, height=9, results=hide, eps=FALSE>>=

corr.tot <- echam
corr.tot$data <- array(cbind(corr.ech$ensmean,corr$ensmean),c(dim(corr$ensmean)[1],1,dim(corr$ensmean)[2]*2)) #[,,1:2,drop=F]
#corr.ech <- corr               # analysis vs. validate
#corr.ech$Analysis <- corr.ech  # echam vs. validate
#names(corr.ech) <- c('echam', 'analysis_localized')
#corr.tot$data <- array(sapply(corr.ech, function(x) x$ensmean), c(nrow(corr.ech[[1]]$ensmean), 1, ncol(corr.ech[[1]]$ensmean)*length(corr.ech)))

#corr.tot$data <- array(corr.tot$data[,,5:8], c(nrow(corr[[1]]$ensmean), 1, 4))
#corr.tot$data <- array(corr.tot$data[,,], c(nrow(corr[[1]]$ensmean), 1, 4))

pdf(paste(figpath,'corr_echam_anal-cru_temp.pdf',sep='/'), width=9, height=6, paper='special')
#pdf(paste(figpath,'corr_echam_anal-cru_temp.pdf',sep='/'), width=9, height=12, paper='special')
layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE), height=c(3,3,1))
#layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,13), 7, 2, byrow = TRUE), height=c(3,3,3,3,3,3,1))
par(oma=c(0,0,0,0))
levs <- c(-1,seq(-0.9,0.9,0.2),1)
#levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0,0.1,0.2,0.3, 0.5, 0.7, 1)
plot_echam(corr.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(corr.tot$data)[3]], lev=levs, st.col=NULL, stations=plstat,add=TRUE)
#plot_echam(corr.tot, cex.pt=1.5, names=pnames[1:dim(corr.tot$data)[3]], lev=levs, st.col=1, stations=plstat,add=TRUE)
dev.off()
#    \caption{temp corr echam ens mean with cru validation (top) and analysis ens mean with cru validation (bottom). Winter left, summer right.}


################################################################################
# Fig. 4b: average precip corr echam and analysis vs validation wrt ens. mean 

pdf(paste(figpath,'corr_echam_anal-cru_precip.pdf',sep='/'), width=9, height=6, paper='special')
layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE), height=c(3,3,1))
par(oma=c(0,0,0,0))
levs <- c(-1,seq(-0.9,0.9,0.2),1)
plot_echam(corr.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(corr.tot$data)[3]], lev=levs, st.col=NULL, stations=plstat, add=TRUE)
dev.off()
#    \caption{precip corr echam ens mean with cru validation (top) and analysis ens mean with cru validation (bottom). Winter left, summer right.}


################################################################################
# Fig. 4b: average slp corr echam and analysis vs validation wrt ens. mean 

pdf(paste(figpath,'corr_echam_anal-cru_slp.pdf',sep='/'), width=9, height=6, paper='special')
layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE), height=c(3,3,1))
par(oma=c(0,0,0,0))
levs <- c(-1,seq(-0.9,0.9,0.2),1)
plot_echam(corr.tot, varname='slp', cex.pt=1.5, names=pnames[1:dim(corr.tot$data)[3]], lev=levs, st.col=NULL, stations=plstat, add=TRUE)
dev.off()
#    \caption{slp corr echam ens mean with cru validation (top) and analysis ens mean with cru validation (bottom). Winter left, summer right.}


################################################################################
# Fig. 5a: average temperature bias echam and analysis vs. validate wrt ens. mean 
#<<label=rmse, echo=FALSE, fig=TRUE, width=8, height=9, results=hide, eps=FALSE>>=
bias.tot <- echam
bias.tot$data <- array(cbind(bias.ech$ensmean,bias$ensmean),c(dim(bias$ensmean)[1],1,dim(bias$ensmean)[2]*2))
#bias.ech <- bias
#bias.ech$Analysis <- bias.ech
#names(bias.ech) <- c('echam', 'analysis_localized')
#bias.tot$data <- array(sapply(bias.ech, function(x) x$ensmean), c(nrow(bias.ech[[1]]$ensmean), 1, ncol(bias.ech[[1]]$ensmean)*length(bias.ech)))

#bias.tot$data <- array(bias.tot$data[,,5:8], c(nrow(bias[[1]]$ensmean), 1, 4))
#bias.tot$data <- array(bias.tot$data[,,], c(nrow(bias[[1]]$ensmean), 1, 4))

pdf(paste(figpath,'bias_echam_anal-cru_temp.pdf',sep='/'), width=9, height=6, paper='special')
#pdf(paste(figpath,'bias_echam_anal-cru_temp.pdf',sep='/'), width=9, height=12, paper='special')
layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE), height=c(3,3,1))
#layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,13), 7, 2, byrow = TRUE), height=c(3,3,3,3,3,3,1))
par(oma=c(0,0,0,0))
levs <- c(-Inf,-10,-5,-2,-1,0,1,2,5,10,Inf)
#if (check_on_anomaly){levs <- c(-Inf,-0.1,-0.05,-0.02,-0.01,0,0.01,0.02,0.05,0.1,Inf)}
#levs <- c(-Inf,seq(-12,12,3),Inf)
#plot_echam(bias.tot, cex.pt=1.5, names=pnames[1:dim(bias.tot$data)[3]], lev=levs, st.col=1, stations=plstat,add=TRUE)
plot_echam(bias.tot, cex.pt=1.5, names=pnames[1:dim(bias.tot$data)[3]], lev=levs, st.col=NULL, stations=plstat,add=TRUE)
dev.off()
#    \caption{}


################################################################################
# Fig. 5b: average precip bias echam and analysis vs. validate wrt ens. mean 

pdf(paste(figpath,'bias_echam_anal-cru_precip.pdf',sep='/'), width=9, height=6, paper='special')
layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE), height=c(3,3,1))
par(oma=c(0,0,0,0))
levs <- c(-Inf,-40,-20,-10,-5,0,5,10,20,40,Inf)
plot_echam(bias.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(bias.tot$data)[3]], lev=levs, st.col=NULL, stations=plstat,add=TRUE)
#plot_echam(bias.tot, cex.pt=1.5, names=pnames[1:dim(bias.tot$data)[3]], lev=levs, st.col=1, stations=plstat,add=TRUE)
dev.off()


################################################################################
# Fig. 5c: average precip bias echam and analysis vs. validate wrt ens. mean 

pdf(paste(figpath,'bias_echam_anal-cru_slp.pdf',sep='/'), width=9, height=6, paper='special')
layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE), height=c(3,3,1))
par(oma=c(0,0,0,0))
levs <- c(-Inf,-10,-5,-2,-1,0,1,2,5,10,Inf)
plot_echam(bias.tot, varname='slp', cex.pt=1.5, names=pnames[1:dim(bias.tot$data)[3]], lev=levs, st.col=NULL, stations=plstat,add=TRUE)
dev.off()




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
pdf(paste(figpath,'re_anom_echam_anal-cru_temp.pdf',sep='/'), width=9, height=3.5, paper='special')
#pdf(paste(figpath,'re_clim_echam-cru_temp.pdf',sep='/'), width=9, height=6, paper='special')
#pdf(paste(figpath,'re_echam_anal-cru_temp_12mon.pdf',sep='/'), width=9, height=12, paper='special')
layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
#layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE), height=c(3,3,1))
#layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,13), 7, 2, byrow = TRUE), height=c(3,3,3,3,3,3,1))
par(oma=c(0,0,0,0))
levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0.05,0.2,0.4,0.6,0.8,1)
#levs <- c(-Inf,seq(-1,1,0.2),Inf)
plot_echam(RE.tot, cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, st.col=NULL, 
           stations=calibrate,add=TRUE)
#plot_echam(RE.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, st.col=1, stations=plstat,add=TRUE)
dev.off()
#    \caption{According to figure \ref{fig:update} but for average reduction of error of the analysis ensemble mean with respect to the ECHAM5 ensemble mean. Positive values indicate that the analysis is closer to the validation measurement than the unconstrained ensemble mean in the respective season.}



################################################################################
# Fig. 6b: average precip RE wrt ens. mean 
#<<label=rmse, echo=FALSE, fig=TRUE, width=8, height=9, results=hide, eps=FALSE>>=

pdf(paste(figpath,'re_anom_echam_anal-cru_precip.pdf',sep='/'), width=9, height=3.5, paper='special')
layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
#layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE), height=c(3,3,1))
par(oma=c(0,0,0,0))
levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0,0.2,0.4,0.6,0.8,1)
plot_echam(RE.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, st.col=NULL, stations=plstat,add=TRUE)
dev.off()
#    \caption{According to figure \ref{fig:update} but for average reduction of error of the analysis ensemble mean with respect to the ECHAM5 ensemble mean. Positive values indicate that the analysis is closer to the validation measurement than the unconstrained ensemble mean in the respective season.}


################################################################################
# Fig. 6c: average slp RE wrt ens. mean 

pdf(paste(figpath,'re_anom_echam_anal-cru_slp.pdf',sep='/'), width=9, height=3.5, paper='special')
#layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE), height=c(3,3,1))
layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
par(oma=c(0,0,0,0))
levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0,0.2,0.4,0.6,0.8,1)
plot_echam(RE.tot, varname='slp', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, st.col=NULL, stations=plstat,add=TRUE)
dev.off()
#    \caption{According to figure \ref{fig:update} but for average reduction of error of the analysis ensemble mean with respect to the ECHAM5 ensemble mean. Positive values indicate that the analysis is closer to the validation measurement than the unconstrained ensemble mean in the respective season.}


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
pdf(paste(figpath,'re_echam_anal-cru_temp.pdf',sep='/'), width=9, height=3.5, paper='special')
#pdf(paste(figpath,'re_clim_echam-cru_temp.pdf',sep='/'), width=9, height=6, paper='special')
#pdf(paste(figpath,'re_echam_anal-cru_temp_12mon.pdf',sep='/'), width=9, height=12, paper='special')
layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
#layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE), height=c(3,3,1))
#layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,13), 7, 2, byrow = TRUE), height=c(3,3,3,3,3,3,1))
par(oma=c(0,0,0,0))
levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0.05,0.2,0.4,0.6,0.8,1)
#levs <- c(-Inf,seq(-1,1,0.2),Inf)
plot_echam(RE.tot, cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, st.col=NULL, 
           stations=calibrate,add=TRUE)
#plot_echam(RE.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, st.col=1, stations=plstat,add=TRUE)
dev.off()
#    \caption{According to figure \ref{fig:update} but for average reduction of error of the analysis ensemble mean with respect to the ECHAM5 ensemble mean. Positive values indicate that the analysis is closer to the validation measurement than the unconstrained ensemble mean in the respective season.}



################################################################################
# Fig. 6b: average precip RE wrt ens. mean 
#<<label=rmse, echo=FALSE, fig=TRUE, width=8, height=9, results=hide, eps=FALSE>>=

pdf(paste(figpath,'re_echam_anal-cru_precip.pdf',sep='/'), width=9, height=3.5, paper='special')
layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
#layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE), height=c(3,3,1))
par(oma=c(0,0,0,0))
levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0,0.2,0.4,0.6,0.8,1)
plot_echam(RE.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, st.col=NULL, stations=plstat,add=TRUE)
dev.off()
#    \caption{According to figure \ref{fig:update} but for average reduction of error of the analysis ensemble mean with respect to the ECHAM5 ensemble mean. Positive values indicate that the analysis is closer to the validation measurement than the unconstrained ensemble mean in the respective season.}


################################################################################
# Fig. 6c: average slp RE wrt ens. mean 

pdf(paste(figpath,'re_echam_anal-cru_slp.pdf',sep='/'), width=9, height=3.5, paper='special')
#layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE), height=c(3,3,1))
layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
par(oma=c(0,0,0,0))
levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0,0.2,0.4,0.6,0.8,1)
plot_echam(RE.tot, varname='slp', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, st.col=NULL, stations=plstat,add=TRUE)
dev.off()
#    \caption{According to figure \ref{fig:update} but for average reduction of error of the analysis ensemble mean with respect to the ECHAM5 ensemble mean. Positive values indicate that the analysis is closer to the validation measurement than the unconstrained ensemble mean in the respective season.}


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
pdf(paste(figpath,'re_clim_echam_anal-cru_temp.pdf',sep='/'), width=9, height=3.5, paper='special')
#pdf(paste(figpath,'re_clim_echam-cru_temp.pdf',sep='/'), width=9, height=6, paper='special')
#pdf(paste(figpath,'re_echam_anal-cru_temp_12mon.pdf',sep='/'), width=9, height=12, paper='special')
layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
#layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE), height=c(3,3,1))
#layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,13), 7, 2, byrow = TRUE), height=c(3,3,3,3,3,3,1))
par(oma=c(0,0,0,0))
levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0.05,0.2,0.4,0.6,0.8,1)
#levs <- c(-Inf,seq(-1,1,0.2),Inf)
plot_echam(RE.tot, cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, st.col=NULL, 
           stations=calibrate,add=TRUE)
#plot_echam(RE.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, st.col=1, stations=plstat,add=TRUE)
dev.off()
#    \caption{According to figure \ref{fig:update} but for average reduction of error of the analysis ensemble mean with respect to the ECHAM5 ensemble mean. Positive values indicate that the analysis is closer to the validation measurement than the unconstrained ensemble mean in the respective season.}



################################################################################
# Fig. 6b: average precip RE wrt ens. mean 
#<<label=rmse, echo=FALSE, fig=TRUE, width=8, height=9, results=hide, eps=FALSE>>=

pdf(paste(figpath,'re_clim_echam_anal-cru_precip.pdf',sep='/'), width=9, height=3.5, paper='special')
layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
#layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE), height=c(3,3,1))
par(oma=c(0,0,0,0))
levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0,0.2,0.4,0.6,0.8,1)
plot_echam(RE.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, st.col=NULL, stations=plstat,add=TRUE)
dev.off()
#    \caption{According to figure \ref{fig:update} but for average reduction of error of the analysis ensemble mean with respect to the ECHAM5 ensemble mean. Positive values indicate that the analysis is closer to the validation measurement than the unconstrained ensemble mean in the respective season.}


################################################################################
# Fig. 6c: average slp RE wrt ens. mean 

pdf(paste(figpath,'re_clim_echam_anal-cru_slp.pdf',sep='/'), width=9, height=3.5, paper='special')
#layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE), height=c(3,3,1))
layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
par(oma=c(0,0,0,0))
levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0,0.2,0.4,0.6,0.8,1)
plot_echam(RE.tot, varname='slp', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, st.col=NULL, stations=plstat,add=TRUE)
dev.off()
#    \caption{According to figure \ref{fig:update} but for average reduction of error of the analysis ensemble mean with respect to the ECHAM5 ensemble mean. Positive values indicate that the analysis is closer to the validation measurement than the unconstrained ensemble mean in the respective season.}


RE.tot <- echam
#RE.bkp <- RE
#RE <- RE.bkp
#anomalies and only 1 analysis with distance weighting
RE.tot$data <- array(RE.echam.clim.anom$ensmean,c(dim(RE$ensmean)[1], 1, dim(RE$ensmean)[2]))
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
pdf(paste(figpath,'re_clim_anom_echam_anal-cru_temp.pdf',sep='/'), width=9, height=3.5, paper='special')
#pdf(paste(figpath,'re_clim_echam-cru_temp.pdf',sep='/'), width=9, height=6, paper='special')
#pdf(paste(figpath,'re_echam_anal-cru_temp_12mon.pdf',sep='/'), width=9, height=12, paper='special')
layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
#layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE), height=c(3,3,1))
#layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,13), 7, 2, byrow = TRUE), height=c(3,3,3,3,3,3,1))
par(oma=c(0,0,0,0))
levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0.05,0.2,0.4,0.6,0.8,1)
#levs <- c(-Inf,seq(-1,1,0.2),Inf)
plot_echam(RE.tot, cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, st.col=NULL, 
           stations=calibrate,add=TRUE)
#plot_echam(RE.tot, varname='temp2', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, st.col=1, stations=plstat,add=TRUE)
dev.off()
#    \caption{According to figure \ref{fig:update} but for average reduction of error of the analysis ensemble mean with respect to the ECHAM5 ensemble mean. Positive values indicate that the analysis is closer to the validation measurement than the unconstrained ensemble mean in the respective season.}



################################################################################
# Fig. 6b: average precip RE wrt ens. mean 
#<<label=rmse, echo=FALSE, fig=TRUE, width=8, height=9, results=hide, eps=FALSE>>=

pdf(paste(figpath,'re_clim_anom_echam_anal-cru_precip.pdf',sep='/'), width=9, height=3.5, paper='special')
layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
#layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE), height=c(3,3,1))
par(oma=c(0,0,0,0))
levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0,0.2,0.4,0.6,0.8,1)
plot_echam(RE.tot, varname='precip', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, st.col=NULL, stations=plstat,add=TRUE)
dev.off()
#    \caption{According to figure \ref{fig:update} but for average reduction of error of the analysis ensemble mean with respect to the ECHAM5 ensemble mean. Positive values indicate that the analysis is closer to the validation measurement than the unconstrained ensemble mean in the respective season.}


################################################################################
# Fig. 6c: average slp RE wrt ens. mean 

pdf(paste(figpath,'re_clim_anom_echam_anal-cru_slp.pdf',sep='/'), width=9, height=3.5, paper='special')
#layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE), height=c(3,3,1))
layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
par(oma=c(0,0,0,0))
levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0,0.2,0.4,0.6,0.8,1)
plot_echam(RE.tot, varname='slp', cex.pt=1.5, names=pnames[1:dim(RE.tot$data)[3]], lev=levs, st.col=NULL, stations=plstat,add=TRUE)
dev.off()
#    \caption{According to figure \ref{fig:update} but for average reduction of error of the analysis ensemble mean with respect to the ECHAM5 ensemble mean. Positive values indicate that the analysis is closer to the validation measurement than the unconstrained ensemble mean in the respective season.}


# 
# 
# ################################################################################
# # Fig. 7a: validation time series 
# 
# #<<label=indices, echo=FALSE, fig=TRUE, width=8, height=6, results=hide, eps=FALSE>>=
# #
# pdf(paste(figpath,'vali_ind_1.pdf',sep='/'), width=9, height=6, paper='special')
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
# pdf(paste(figpath,'vali_ind_2.pdf',sep='/'), width=9, height=6, paper='special') 
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
  
pdf(paste(figpath,'vali_ind.pdf',sep='/'), width=9, height=6, paper='special') 
par(mfrow=c(4,2), cex.axis=1.4, cex.lab=1.4, mar=c(3,5,1,1), oma=c(0,0,0,0))

vind2=vind.allts
vind2$data=array(vind.allts$data,c(dim(vind.allts$data)[1],2,dim(vind.allts$data)[2]/2))
eind2=eind.allts
eind2$ensmean=array(eind.allts$ensmean,c(dim(eind.allts$ensmean)[1],2,dim(eind.allts$ensmean)[2]/2))
eind2$min=array(apply(eind.allts$data,1:2,min),c(dim(eind.allts$ensmean)[1],2,dim(eind.allts$ensmean)[2]/2))
eind2$max=array(apply(eind.allts$data,1:2,max),c(dim(eind.allts$ensmean)[1],2,dim(eind.allts$ensmean)[2]/2))
aind2=aind.allts
aind2$ensmean=array(aind.allts$ensmean,c(dim(aind.allts$ensmean)[1],2,dim(aind.allts$ensmean)[2]/2))
aind2$min=array(apply(aind.allts$data,1:2,min),c(dim(aind.allts$ensmean)[1],2,dim(aind.allts$ensmean)[2]/2))
aind2$max=array(apply(aind.allts$data,1:2,max),c(dim(aind.allts$ensmean)[1],2,dim(aind.allts$ensmean)[2]/2))
t<-vind.allts$time[seq(1,length(vind.allts$time),by=2)]

for (ind in c(31,9,19,29)) {
if (ind == 31) {mainname='NHt2m'}
if (ind == 9) {mainname='NEUt2m'}
if (ind == 19) {mainname='NEUpr'}
if (ind == 29) {mainname='NEUslp'}
 for (seas in c(1,2)) {
  if (seas == 1) {color='black'; color3='blue'; color2='cyan'}
  if (seas == 2) {color='black'; color3='red'; color2='orange'}
   
    ymin=min(eind2$min[ind,seas,])
    ymax=max(eind2$max[ind,seas,])

  plot(t,vind2$data[ind,seas,],ty='l',col=color,lty=1,ylim=c(ymin,ymax),main=mainname)
  lines(t,eind2$ensmean[ind,seas,],ty='l',col=color2,lwd=2,lty=2,main='')
  lines(t,eind2$min[ind,seas,],ty='l',col=color2,lty=2,lwd=1,main='')
  lines(t,eind2$max[ind,seas,],ty='l',col=color2,lty=2,lwd=1,main='')
  lines(t,aind2$ensmean[ind,seas,],ty='l',col=color3,lwd=2,lty=3,main='')
  lines(t,aind2$min[ind,seas,],ty='l',col=color3,lty=3,lwd=1,main='')
  lines(t,aind2$max[ind,seas,],ty='l',col=color3,lty=3,lwd=1,main='')
 }
}
dev.off()



################################################################################
# Fig. 8:  index correlation 
# new like in the paper
#<<label=indices_corr, echo=FALSE, fig=TRUE, width=12, height=7, results=hide, eps=FALSE>>=

inds <- c('NH.temp2', 'NEU.temp2', 'NEU.precip', 'NAM.temp2', 'AFR.temp2')
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
          RE.mat[k,se,1:29] <- RE.ind$data[RE.ind$name == inds[k],se,]
          RE.mat[k,se,30] <- RE.ind$ensmean[RE.ind$name == inds[k],se]
          corr.mat[k*2,se,1:29] <- acorr.ind$data[acorr.ind$name == inds[k],se,]
          corr.mat[k*2,se,30] <- acorr.ind$ensmean[acorr.ind$name == inds[k],se]
          corr.mat[k*2-1,se,1:29] <- ecorr.ind$data[ecorr.ind$name == inds[k],se,]
          corr.mat[k*2-1,se,30] <- ecorr.ind$ensmean[ecorr.ind$name == inds[k],se]
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
          RE.mat[k,se,1:30] <- RE.ind$data[RE.ind$name == inds[k],se,]
          RE.mat[k,se,31] <- RE.ind$ensmean[RE.ind$name == inds[k],se]
          corr.mat[k*2,se,1:30] <- acorr.ind$data[acorr.ind$name == inds[k],se,]
          corr.mat[k*2,se,31] <- acorr.ind$ensmean[acorr.ind$name == inds[k],se]
          corr.mat[k*2-1,se,1:30] <- ecorr.ind$data[ecorr.ind$name == inds[k],se,]
          corr.mat[k*2-1,se,31] <- ecorr.ind$ensmean[ecorr.ind$name == inds[k],se]
        }
    }
}
    ind.col <- list()
    ind.col2 <- list()
    ind.col[[1]] <- hcl(c(0,0,120,240,270,300), c=40, l=50)
    ind.col[[2]] <- hcl(c(0,0,120,240,270,200), c=40, l=50)
    ind.col2[[1]] <- hcl(rep(c(0,0,120,240,270,300), each=2), c=40, l=c(90,50))
    ind.col2[[2]] <- hcl(rep(c(0,0,120,240,270,200), each=2), c=40, l=c(90,50))

#pdf('figures/inst/indices.pdf', width=9, height=6, paper='special')
pdf(paste(figpath,'indices.pdf',sep='/'), width=9, height=6, paper='special')
    par(mfrow=c(2,2), mar=c(1,1,0,0), oma=c(0,2,1,0), cex.axis=1.4, cex.lab=1.4)
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
indind <- c(1,2,4,5,7,8,10,11,13,14)
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


# 
# ################################################################################
# # Tab x: RE for indices 
# 
# #<<label=tabREwinter, echo=FALSE, results=tex>>=
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
# Fig. 9a: temperature reconstruction skill 

#<<label=Giorgitemp2, echo=FALSE, fig=TRUE, eps=FALSE, width=12, height=14>>=
if (pseudoproxy){
  giorgi.RE <- compute_avg_RE_pseudoproxy(H.giorgi, echam, analysis[[2]], validate)
  giorgi.corr <- compute_avg_corr_pseudoproxy(H.giorgi, echam, analysis[[i]], validate)
} else {
#  giorgi.RE <- compute_avg_RE(H.giorgi, echam, analysis[[2]], validate)
#  giorgi.corr <- compute_avg_corr(H.giorgi, echam, analysis[[2]], validate)
  giorgi.RE <- compute_avg_RE(H.giorgi, echam, analysis, validate)
  giorgi.corr <- compute_avg_corr(H.giorgi, echam, analysis, validate)
}
varname <- c('temp2', 'precip', 'slp')
varn <- 1

pdf(paste(figpath,'RE_giorgi_temp.pdf',sep='/'), width=9, height=6, paper='special')
#if (pseudoproxy){
  par(mfrow=c(4,1), mar=c(1,3,1,1), oma=c(0,0,0,0), cex.axis=1.2, cex.lab=1.2)
#} else {
#  par(mfrow=c(2,1), mar=c(1,3,1,1), oma=c(0,0,0,0), cex.axis=1.2, cex.lab=1.2)
#}
  ncorr <- length(giorgi.short)
  RE.mat <- giorgi.RE[ncorr*(varn-1) + 1:ncorr,,]
  indind <- 1:ncorr
  for (se in 1:2){
    plot(0, type='n', xaxt='n', ylim=c(-1,1), xlim=c(1,max(indind)), xlab='', ylab='')
    polygon(c(-3, -3, max(indind)+5, max(indind+5)), c(0,-2,-2,0), border=NA, col=grey(0.8))
    if (pseudoproxy){
      boxplot(t(RE.mat[,se,1:29]), at=indind, col=hcl((varn-1)*120, c=40, l=50),
          add=T, xaxt='n', yaxt='n', xlab='', ylab='', lwd=1, boxwex=0.5)
      points(indind - 0.3, RE.mat[,se,30], pch='>', lwd=3, cex=1.2)
    } else {
      boxplot(t(RE.mat[,se,1:30]), at=indind, col=hcl((varn-1)*120, c=40, l=50),
            add=T, xaxt='n', yaxt='n', xlab='', ylab='', lwd=1, boxwex=0.5)
      points(indind - 0.3, RE.mat[,se,31], pch='>', lwd=3, cex=1.2)
    }
    text(indind, rep(-0.8,length(indind)), giorgi.short, cex=1, adj=c(0.5,1.3))
    text(0.5,1,c('Winter', 'Summer')[se], cex=1.2, adj=c(0,1))
    abline(v=seq(1.5, max(indind), 1), lty=3)
  }


ncorr <- length(giorgi.short)*2
corr.mat <- giorgi.corr[ncorr*(varn-1) + 1:ncorr,,]
indind <- rbind(seq(1, by=4, length=ncorr/2), seq(2, by=4, length=ncorr/2))

for (se in 1:2){
  plot(0, type='n', xaxt='n', ylim=c(-1,1), xlim=c(2,max(indind)-1), xlab='', ylab='')
  polygon(c(-3, -3, max(indind)+5, max(indind+5)), c(0,-2,-2,0), border=NA, col=grey(0.8))
  if (pseudoproxy){
    boxplot(t(corr.mat[,se,1:29]), at=as.vector(indind), col=hcl((varn-1)*120, c=40, l=c(90,50)), add=T, xaxt='n', yaxt='n', xlab='', ylab='', lwd=1)
  points(apply(indind, 2, mean) - 1.4, corr.mat[seq(1,ncorr,2),se,30], pch='>', lwd=3, cex=1.2)
  points(apply(indind, 2, mean) + 1.4, corr.mat[seq(2,ncorr,2),se,30], pch='<', lwd=3, cex=1.2)
  } else {
   boxplot(t(corr.mat[,se,1:30]), at=as.vector(indind), col=hcl((varn-1)*120, c=40, l=c(90,50)), add=T, xaxt='n', yaxt='n', xlab='', ylab='', lwd=1)
  points(apply(indind, 2, mean) - 1.4, corr.mat[seq(1,ncorr,2),se,31], pch='>', lwd=3, cex=1.2)
  points(apply(indind, 2, mean) + 1.4, corr.mat[seq(2,ncorr,2),se,31], pch='<', lwd=3, cex=1.2) 
  }  
  text(apply(indind, 2, mean), rep(-0.8,ncol(indind)), giorgi.short, cex=1, adj=c(0.5,1.3))
  text(0.5,1,c('Winter', 'Summer')[se], cex=1.2, adj=c(0,1))
  abline(v=seq(3.5,max(indind), by=4), lty=3)
}
dev.off() 
#    \caption{Skill in reconstructing area-average temperature in different subcontinental regions as defined in \citet{Giorgi2000}}
 
 
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
# layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE), height=c(3,3,1))
# par(mfrow=c(2, 2))
# for (s in 1:2) { ## if s==1 then crps is crps.ana which is the score of the analysis, if s==2 crps is the difference 
#   if (s==1) {    ## in score of echam and analysis. It should be positive. 
#     crps <- crps.ana[,seq(2,ncol(crps.ana),2)] ## the one that's updated (rem: proxies only half a year)
#     crps[which(is.na(crps))]<-0
#     title<-"analysis"
#   } else { 
#     crps <- crps.ech[,seq(2,ncol(crps.ech),2)]-crps.ana[,seq(2,ncol(crps.ana),2)] ## difference of score between updated and not, should be positive
#     if (PAGES){
#       crps <- crps.ech-crps.ana
#     }
#     crps[which(is.na(crps))]<-0
#     title<-"Score Diff: CCC400 - EKF400"
#   }
#   crps.tot<-data.frame(averaged=apply(crps,1,mean))
#   crps.tot$lon<-echam$lon
#   crps.tot$lat<-echam$lat
#   crps.tot$names<-echam$names
#   
#   
#   
#   
#   crps.plot<-crps.tot[which(crps.tot$names=="temp2"),]
#   if (s==1){
#     lev <- pretty(crps.plot$averaged,12) 
#     br <- length(lev)
#     colpal <- two.colors(n=br,start="white",middle="orange", end="red", alpha=1.0)
#   }else{
#     # lev <- c(-0.05,-0.04,-0.03,-0.02,-0.01,0.03,0.06,0.09,0.12,0.15)
#     lev <- pretty(crps.plot$averaged,12) # for PAGES and NTREND
#     br <- length(lev)
#     colpal <- two.colors(n=br,start="blue",middle="white", end="red", alpha=1.0)
#   }
#   
#   datcol <- colpal[as.numeric(cut(crps.plot$averaged,breaks=lev))] 
#   datcol[which(crps.plot$averaged>(-0.02)&crps.plot$averaged<0.02)]<-"#FFFFFF" #NTREND
#   
#   # png(filename=paste0(plotintdir,expname,'/map_DOC_month_',m,'.png'),width=1400,height = 800)
#   plot(crps.plot$lon, crps.plot$lat,
#        xlim=c(-180,180),ylim=c(-90,90),   
#        cex=1.2,                   # cex=point size~precip
#        pch=15,                    # point type: 15 is square with possible color filling
#        col=datcol[which(!is.na(crps.plot$averaged))],
#        xlab="Longitude",ylab="Latitude",main=paste0("Temperature, ",title))# point fill color
#   map("world",interior=F,add=T,ty='l',col='black',xlim=c(-180,180),ylim=c(-90,90))
#   legend("bottomleft", inset=0.01, as.character(lev),fill=colpal,bty="o",
#          bg="white",box.col="white",box.lwd=0,cex=0.7)
#   
#   crps.plot<-crps.tot[which(crps.tot$names=="precip"),]
#   if (s==1){
#     lev <- pretty(crps.plot$averaged,12) 
#     br <- length(lev)
#     colpal <- two.colors(n=br,start="white",middle="orange", end="red", alpha=1.0)
#   }else{
#     # lev <- c(-0.5,-0.4,-0.3,-0.2,-0.1,0.3,0.6,0.9,1.2,1.5) # for normal proxies
#     # lev <- c(-2.2,-2,-1.8,-1.6,-1.4,-1.2,-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1,1.2,1.6,1.8,2,2.2) # for PAGES
#     lev <- c(-1.2,-1,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1,1.2) #NTREND
#     br <- length(lev)
#     colpal <- two.colors(n=br,start="blue",middle="white", end="red", alpha=1.0)
#   }
#   datcol <- colpal[as.numeric(cut(crps.plot$averaged,breaks=lev))] 
#   datcol[which(crps.plot$averaged>(-0.2)&crps.plot$averaged<0.2)]<-"#FFFFFF" #NTREND
#   # png(filename=paste0(plotintdir,expname,'/map_DOC_month_',m,'.png'),width=1400,height = 800)
#   plot(crps.plot$lon, crps.plot$lat,
#        xlim=c(-180,180),ylim=c(-90,90),   
#        cex=1.2,                   # cex=point size~precip
#        pch=15,                    # point type: 15 is square with possible color filling
#        col=datcol[which(!is.na(crps.plot$averaged))],
#        xlab="Longitude",ylab="Latitude",main=paste0("Precip, ",title))# point fill color
#   map("world",interior=F,add=T,ty='l',col='black',xlim=c(-180,180),ylim=c(-90,90))
#   legend("bottomleft", inset=0.01, as.character(lev),fill=colpal,bty="o",
#          bg="white",box.col="white",box.lwd=0,cex=0.7)
#   # dev.off()
#   
# }
# dev.off()
# par(mfrow=c(1, 1))



### plot CRPS with plot_echam Temp


crps <- echam

crps$data<-array(cbind(crps.ana.winter,crps.ana.summer),c(dim(echam$ensmean)[1], 1, 2))

pdf(paste(figpath,'crps.ana.clim_temp.pdf',sep='/'), width=9, height=3.5, paper='special')

layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))

par(oma=c(0,0,0,0))
levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0.05,0.2,0.4,0.6,0.8,1)

plot_echam(crps, cex.pt=1.5, names=pnames[1:dim(crps$data)[3]], lev=levs, st.col=NULL, 
           stations=calibrate,add=TRUE)

dev.off()

### plot CRPS with plot_echam Precip

pdf(paste(figpath,'crps.ana.clim_precip.pdf',sep='/'), width=9, height=3.5, paper='special')

layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))

par(oma=c(0,0,0,0))
levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0.05,0.2,0.4,0.6,0.8,1)

plot_echam(crps, varname='precip', cex.pt=1.5, names=pnames[1:dim(crps$data)[3]], lev=levs, st.col=NULL, 
           stations=calibrate,add=TRUE)

dev.off()

### plot CRPS score difference ech-ana (positive values are good) with plot_echam: Temp

crps <- echam


crps_ech_ana.winter <- crps.ech.winter-crps.ana.winter
crps_ech_ana.summer <- crps.ech.summer-crps.ana.summer


crps$data<-array(cbind(crps_ech_ana.winter,crps_ech_ana.summer),c(dim(echam$ensmean)[1], 1, 2))

pdf(paste(figpath,'crps.clim_ech-ana_temp.pdf',sep='/'), width=9, height=3.5, paper='special')

layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))

par(oma=c(0,0,0,0))
levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0.01,0.02,0.04,0.08,0.1,0.12)

plot_echam(crps, cex.pt=1.5, names=pnames[1:dim(crps$data)[3]], lev=levs, st.col=NULL, 
           stations=calibrate,add=TRUE)

dev.off()

### plot CRPS score difference ech-ana (positive values are good) with plot_echam: Precip

pdf(paste(figpath,'crps.clim_ech-ana_precip.pdf',sep='/'), width=9, height=3.5, paper='special')

layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))

par(oma=c(0,0,0,0))
levs <- c(-Inf, -10,-3,-1,-0.3, -.1, 0.05,0.2,0.4,0.6,0.8,1)

plot_echam(crps,varname="precip", cex.pt=1.5, names=pnames[1:dim(crps$data)[3]], lev=levs, st.col=NULL, 
           stations=calibrate,add=TRUE)

dev.off()