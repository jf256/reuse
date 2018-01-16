

## This script contains 3 loops:

##  1. the yearly loop: 
##  it goes from syr to eyr. For each loop it loads the echam, analysis 
##  and calibrate. then during the first loop it only calculate the sigma and its derivates in 
##  order to estimate the error. The errors (and its derivates) are then saved in all1 for each 
##  year. After the loop has reached eyr it will calculate Rcalnew which is the final error averaged 
##  over all stations and all time. Watch that you take a corenr that is a multiple of nyr for obvious 
##  reasons.
##  Also do library(parallel) if the script doesnt seem to be working

##  2. the iterary loop:
##  the iterary loop then uses the calculated Rcalnew as the new error Rcal. then again the model 
##  is run for all the years. For each iterary loop the error (Rcalnew) is estimated a bit better 
##  meaning the value converges.

##  3. the errtime loop:
##  the errtime loop is used to calculate the converged error not for the whole period (nyr) but 
##  for increments of this time. For example if you want to know the evolution of the error over 100 
##  years you can't just run the model for 100 years and the do (e.g.) three iterations because then 
##  you will get the converged error estimation of those 100 years, but you want the converged error 
##  for every single year (or every second year). the errtime loop does that.
##  An example: initsyr=1940, initeyr=1949, niter=3, errtime=2
##  for (n in seq(from=1, to=nyr, by=errtime)) calculates the error with niter iterations for 2 year 
##  time blocks. You'll get five Rcalnew's each one a result of 3 iterations.


##  data
rm(list=ls())


workdir='/scratch4/lucaf/reuse/reuse_git/'
setwd(workdir)
dataextdir='/mnt/climstor/giub/EKF400/'
dataintdir='/scratch4/lucaf/reuse/data/analysis/obs.error/'
plotintdir='/scratch4/lucaf/reuse/figures/obs.error/'
ana_dir="/scratch4/joerg/projects/reuse/data/analysis/EKF400_v1.3_corr_echam_clim/"

## input
expname="obs_error_inflation_factor_1604-2003_iter_8_errtime_40"
initsyr=1604
initeyr=2003 # always watch the nr of year you pick! because if nyr for ex. = 5 then 
         # errtime of 2 makes little sense
niter=8 # 1660-1999, errtime=20,niter=4, nrcores=10 geht 11:40h
errtime=40 # error is calculated for errtime years
errtimeloop=TRUE #if F then you also need to comment lines 66-74
calc.delta=TRUE
nrcores=5 ## take a multiple of errtime
notop4=FALSE
cutoffoutliers=FALSE




## go
syr=initsyr
eyr=initeyr
yr <- syr:eyr
nyr=(eyr-syr)+1

dir.create(paste0("../data/analysis/obs.error"))
dir.create(paste0("../data/analysis/obs.error/",expname))
dir.create(paste0("../figures/obs.error"))
dir.create(paste0("../figures/obs.error/",expname))

ptm <- proc.time()

for (n in seq(from=1, to=nyr, by=errtime)){ # calculates the error of errtime year until eyr to
  syr=initsyr                               # come up how the error changes over time
  eyr=initeyr
  nyr=(eyr-syr)+1
syr=syr-1+n
eyr=syr-1+errtime
nyr=eyr-syr+1
yr <- syr:eyr
print(paste0("years ",syr,'-',eyr))

Rcalnew <- array(0,dim=c(niter,3))
Pnew <- array(0,dim=c(niter,3))
deltanew <- array(NA,dim=c(niter,3))


# here, we write the results to files
comb <- function(x, ...) {
  lapply(seq_along(x),
         function(r) c(x[[r]], lapply(list(...), function(y) y[[r]])))
}

for (iter in 1:niter) {
  print(paste0("iteration ",iter))
  ## {doParallel} depends on {foreach} and {iterators}
  # library(doParallel)
  
  ## reserve the number of cores
  registerDoParallel(cores=nrcores)
  
  ## the yearly foreach loop 
  all<- foreach(r=1:length(yr),.combine='comb', .multicombine=TRUE,
                .init=list(list(),list(), list(), list(), list(), list(),list(),list(),list(),list(),list(),list(), list(), list())) %dopar% {

                  source('EnSRF_switches.R')
                  source('EnSRF_functions.R')
                  
                  if (syr-1+r > 1659) {
                    instrumental=T        # all instrumental stations
                  } else {
                    instrumental=F
                  }
                  if (syr-1+r < 1960) {
                    real_proxies=T         # Proxy data experiment (regression NOT H operator)
                  } else {
                    real_proxies=F
                  }
                  if (syr-1+r > 1853) {
                    docum=F                 # read documentary based data
                  } else {
                    docum=T
                  }
                  if (import_luca==TRUE & syr-1+r>1740 & syr-1+r<1900) {
                    load(paste0(workdir,'../data/analysis/new_data_luca_1740-1899/analysis_',syr-1+r,'_2ndgrid.Rdata'))
                    docum=T
                  }else{
                    load(paste0(ana_dir,'analysis_',syr-1+r,'_2ndgrid.Rdata'))
                  }
        echam <- echam.anom
        analysis <- analysis.anom
        rm(echam.abs,echam.anom)
        
        if (notop4) {
          load(paste0(workdir,'../data/analysis/archive/j-loop-archive_1660-1999/top_4.Rdata'))
          notop4<- which(calibrate$lon==top_4[1,2]&top_4[1,3]|calibrate$lon==top_4[2,2]&top_4[2,3]
                |calibrate$lon==top_4[3,2]&top_4[3,3]|calibrate$lon==top_4[4,2]&top_4[4,3])
          calibrate$data <- calibrate$data[-c(notop4),]
          calibrate$names <- calibrate$names[-c(notop4)]
          calibrate$sour <- calibrate$sour[-c(notop4)]
          calibrate$lon <- calibrate$lon[-c(notop4)]
          calibrate$lat <- calibrate$lat[-c(notop4)]
        }
        
        
        
        if (tps_only) {
          tpspos <- sort(c(which(echam$names=='temp2'), which(echam$names=='precip'), 
                           which(echam$names=='slp'), which(echam$names=='bias')))
          echam$data <- echam$data[tpspos,,]
          echam$ensmean <- echam$ensmean[tpspos,]
          echam$lon <- echam$lon[tpspos]
          echam$lat <- echam$lat[tpspos]
          echam$names <- echam$names[tpspos]
          
          
          tpspos <- sort(c(which(analysis$names=='temp2'), which(analysis$names=='precip'), 
                           which(analysis$names=='slp'), which(analysis$names=='bias')))
          analysis$data <- analysis$data[tpspos,,]
          analysis$ensmean <- analysis$ensmean[tpspos,]
          analysis$lon <- analysis$lon[tpspos]
          analysis$lat <- analysis$lat[tpspos]
          analysis$names <- analysis$names[tpspos]
          
        }
        
        
        lvec <- rep(l_dist_ind, length(unique(echam$names)))
        names(lvec) <- unique(echam$names)
        lvec['slp'] <- l_dist_slp
        lvec['precip'] <- l_dist_precip
        lvec['temp2'] <- l_dist_temp2 
        lvec['gph500'] <- l_dist_gph500
        lvec['gph100'] <- l_dist_gph100
        lvec['u850'] <- l_dist_u850
        lvec['u200'] <- l_dist_u200
        lvec['v850'] <- l_dist_v850
        lvec['v200'] <- l_dist_v200
        lvec['omega500'] <- l_dist_omega500
        lvec['t850'] <- l_dist_t850
        
        
        
        
        inst=calibrate
        inst$data=calibrate$data[calibrate$sour=="inst",]
        inst$lon=calibrate$lon[calibrate$sour=="inst"]
        inst$lat=calibrate$lat[calibrate$sour=="inst"]
        inst$names=calibrate$names[calibrate$sour=="inst"]
        
        docall=calibrate
        docall$data=docall$data[docall$sour=="doc",]
        docall$lon=docall$lon[docall$sour=="doc"]
        docall$lat=docall$lat[docall$sour=="doc"]
        docall$names=docall$names[docall$sour=="doc"]

        realprox=calibrate
        realprox$data=calibrate$data[calibrate$names=="prox",]
        realprox$lon=calibrate$lon[calibrate$names=="prox"]
        realprox$lat=calibrate$lat[calibrate$names=="prox"]
        realprox$var_residu=calibrate$var_residu[which(!is.na(calibrate$var_residu))]
        realprox$mr=calibrate$mr[calibrate$names=="prox",]
        realprox$names=calibrate$names[calibrate$names=="prox"]
        
        etmp <- echam
        etmp$lon[is.na(etmp$lon)] <- 0
        etmp$lat[is.na(etmp$lat)] <- -90
        
        if (instrumental) {
          print(paste('Calculating Hcal1'))
          Hcal1 <- array(NA,dim=c(dim(inst$data)[1],2))
          Hcal1 <- compute_Hi_Hredux_sixmonstatevector(inst, etmp, threshold=700)
        }
        if (real_proxies) {
          print(paste('Calculating Hcal2'))
          Hcal2 <- array(NA,dim=c(dim(realprox$data)[1],14))
          Hcal2 <- compute_Hi_Hredux_proxy(realprox, etmp, realprox$mr, threshold=700)
        }
        if (docum) {
          print(paste('Calculating Hcal3'))
          Hcal3 <- array(NA,dim=c(dim(docall$data)[1],2))
          Hcal3 <- compute_Hi_Hredux_sixmonstatevector(docall, etmp, threshold=700)
          Hcal3[Hcal3==0] <- NA
        }
        rm(etmp)
        if (instrumental & !docum & sixmonstatevector & !real_proxies) {
          H.i <- array(NA,c(nrow(Hcal1),1))
        } else if (instrumental & docum & sixmonstatevector & !real_proxies) {
          H.i <- array(NA,c((nrow(Hcal1)+nrow(Hcal3)),1))
        } else if (!instrumental & docum & sixmonstatevector & real_proxies) {
          H.i <- array(NA,c((nrow(Hcal2)+nrow(Hcal3)),7))
        } else if (instrumental & !docum & sixmonstatevector & real_proxies) {
          H.i <- array(NA,c((nrow(Hcal1)+nrow(Hcal2)),7))
        } else { 
          H.i <- array(NA,c((nrow(Hcal1)+nrow(Hcal2)+nrow(Hcal3)),7))
        }
        Hredux <- H.i
        if (instrumental){
          H.i[1:nrow(Hcal1),1] <- Hcal1[,1] 
          Hredux[1:nrow(Hcal1),1] <- Hcal1[,2] 
        } 
        if (instrumental & real_proxies & docum) {
          H.i[(nrow(Hcal1)+1):(nrow(Hcal1)+nrow(Hcal2)),] <- Hcal2[,c(1,3,5,7,9,11,13)] 
          Hredux[(nrow(Hcal1)+1):(nrow(Hcal1)+nrow(Hcal2)),] <- Hcal2[,c(2,4,6,8,10,12,14)]   
          H.i[((nrow(Hcal1)+nrow(Hcal2)+1):nrow(H.i)),1] <- Hcal3[,1]
          Hredux[((nrow(Hcal1)+nrow(Hcal2)+1):nrow(H.i)),1] <- Hcal3[,2]
        }
        if (!instrumental & real_proxies & docum) {
          H.i[1:nrow(Hcal2),] <- Hcal2[,c(1,3,5,7,9,11,13)] 
          Hredux[1:nrow(Hcal2),] <- Hcal2[,c(2,4,6,8,10,12,14)]   
          H.i[((nrow(Hcal2)+1):nrow(H.i)),1] <- Hcal3[,1]
          Hredux[((nrow(Hcal2)+1):nrow(H.i)),1] <- Hcal3[,2]
        }
        if (instrumental & real_proxies & !docum) {
          H.i[(nrow(Hcal1)+1):(nrow(Hcal1)+nrow(Hcal2)),] <- Hcal2[,c(1,3,5,7,9,11,13)] 
          Hredux[(nrow(Hcal1)+1):(nrow(Hcal1)+nrow(Hcal2)),] <- Hcal2[,c(2,4,6,8,10,12,14)] 
        }
        if (instrumental & !real_proxies & docum) {
          H.i[(nrow(Hcal1)+1):(nrow(H.i)),] <- Hcal3[,1]
          Hredux[(nrow(Hcal1)+1):(nrow(H.i)),] <- Hcal3[,2]
        }
        H.i[H.i==0] <- NA
        Hredux[Hredux==0] <- NA
        
        # if ((real_proxies) & ((instrumental) | (docum))) { is it okay to leave this out because 
        #                                                    for 1961-1965 real_proxies is FALSE?!?
        if (iter==1){
          Rcal <- c(temp2=0.9, precip=50, slp=10)[calibrate$names]
          Rcal[calibrate$sour=="doc"] <- 0.25 # equals 0.5 std. dev.
          
          #        if (avg_prox_per_grid) {Rcal <- Rcal*(1/calibrate$numavg)}
        } else {
          Rcal <- c(temp2=Rcalnew[(iter-1),1], precip=50, slp=Rcalnew[(iter-1),2])[calibrate$names]
          Rcal[calibrate$sour=="doc"] <- Rcalnew[iter-1,3]
        }
        Rcal[calibrate$names=="prox"] <- realprox$var_residu/2 
        # previously used residuals/2 for 1. paper version to give proxies more weight
        # better delete "/2"
        # probably should have given instrumentals more error instead!
        # }
        
        
        ntim <- ncol(analysis$data)
        nens <- dim(analysis$data)[3]
        nprox <- nrow(calibrate$data) 
        
      
        
        sigma <-array(NA,dim=c(length(which(calibrate$sour!="prox")),ntim))
        sigmab <-array(NA,dim=c(length(which(calibrate$sour!="prox")),ntim))
        echam.rel <- array(NA,dim=c(length(which(calibrate$sour!="prox")),ntim))
        analysis.rel <- array(NA,dim=c(length(which(calibrate$sour!="prox")),ntim))
        deltaup <- array(NA,dim=dim(calibrate$data))
        deltadown <- array(NA,dim=dim(calibrate$data))
        delta<- array(NA, dim=c(1,3))
        

        
        
          analysis <- echam
          analysis$data <- echam$data - as.vector(echam$ensmean)
          for (j in 1:nprox){
            # for (j in 1:500){
            if (j %% 100 == 0) {
              print(paste('Assimilating observation', j))
            }
            hisna <- is.na(H.i[j,])
            h.i <- H.i[j,!hisna,drop=F]
            H <-  t(as.matrix(Hredux[j,!hisna]))
            
            if (!is.na(h.i[1])) {
              pos <- which(is.na(echam$lon))
              echam$lon[pos] <- 0
              echam$lat[pos] <- -90
              dist <- compute_dist_2d(echam$lon, echam$lat, echam$lon[h.i], echam$lat[h.i]) 
              # weights are a matrix of ndim x nh (non-null elements in H, here different months)
              # NO temporal correlation for use of monthly instr. data and 1 col H operator          
              wgt <- corr_function(dist,outer(lvec[echam$names], lvec[echam$names[h.i]], pmin))
              wgt[which(echam$names %in% c("DIMI","z100","z300","PWC","HC","SJ")),] <- 1
              # temporal correlation quickly drops to zero (~ 0.6 for 1-month lag, 
              # ~0.4 for two months, etc.)
            } else {
              dist <- NA
              wgt <- NA
            }
            if ((!is.na(dist[1])) & (!is.na(wgt[1]))){
              for (i in 1:ntim){
                if (!is.na(calibrate$data[j,i])) {
                  
                  x2tmp <- analysis$data[,i,] # entire state vector at time step i, all ensemble members
                  x2 <- x2tmp[h.i,,drop=F] # state vector at time step i and h.i, all ensemble members
                  PH <- (analysis$data[,i,] %*% t(x2) / (nens - 1) * wgt) %*% t(H)
                  if (iter!=1&calc.delta) {
                    if (calibrate$names[j]=="temp2"&calibrate$sour[j]=="inst") {
                    PH<- PH*deltanew[iter-1,1]
                    }else if (calibrate$names[j]=="slp") {
                      PH<- PH*deltanew[iter-1,2]
                    }else if (calibrate$sour[j]=="doc"){
                      PH<- PH*deltanew[iter-1,3]
                    }
                    }
                  HPHR <- as.vector(H %*% PH[h.i,] + Rcal[j])
                  K <- PH / HPHR
                  Ktilde <- K / (1 + sqrt(Rcal[j]/HPHR))
                  analysis$ensmean[,i] <- analysis$ensmean[,i] + K[,1] * (calibrate$data[j,i] -
                                                                            H %*% analysis$ensmean[h.i,i])
                  analysis$data[,i,] <- analysis$data[,i,] - Ktilde %*% H %*% analysis$data[h.i,i,]
                  
                  if (calibrate$sour[j]!="prox") {
                  
                  deltaup[j,i] <- ((calibrate$data[j,i] - echam$ensmean[h.i,i])*(calibrate$data[j,i] - echam$ensmean[h.i,i]))
                  deltadown[j,i] <- H %*% PH[h.i]
                  }
                  
                  if (length(which(!is.nan(analysis$data)))==0){
                    stop()
                    }
                  }
              }
            }
          }
          
          if (instrumental) {
          for (j in which(calibrate$sour=="inst")){
            # for (j in 1:500){
            if (j %% 100 == 0) {
              print(paste('Instrumental Errors', j))
            }
            hisna <- is.na(H.i[j,])
            h.i <- H.i[j,!hisna,drop=F]
            H <-  t(as.matrix(Hredux[j,!hisna]))
            
            for (i in 1:ntim){
              if (!is.na(calibrate$data[j,i])) {
                sigmab[j,i] <-  (analysis$ensmean[h.i,i] - echam$ensmean[h.i,i]) *
                  (calibrate$data[j,i] - echam$ensmean[h.i,i])
                sigma[j,i] <-  ((calibrate$data[j,i] - analysis$ensmean[h.i,i]) *
                                  (calibrate$data[j,i] - echam$ensmean[h.i,i]))
               
                analysis.rel[j,i] <- analysis$ensmean[h.i,i]
                echam.rel[j,i] <- echam$ensmean[h.i,i]
               
                }
              }
            }
          }
          
          
          nproxies <- length(which(calibrate$sour=="prox"))
          if (docum) {
          for (j in which(calibrate$sour=="doc")) {
            if (j %% 10 == 0) {
              print(paste('Document Errors', j))
            }
            hisna <- is.na(H.i[j,])
            h.i <- H.i[j,!hisna,drop=F]
            H <-  t(as.matrix(Hredux[j,!hisna]))
            
            for (i in 1:ntim){
              if (!is.na(calibrate$data[j,i])) {
                sigmab[j-nproxies,i] <-  (analysis$ensmean[h.i,i] - echam$ensmean[h.i,i]) *
                  (calibrate$data[j,i] - echam$ensmean[h.i,i])
                sigma[j-nproxies,i] <-  ((calibrate$data[j,i] - analysis$ensmean[h.i,i]) *
                                  (calibrate$data[j,i] - echam$ensmean[h.i,i]))

                analysis.rel[j-nproxies,i] <- analysis$ensmean[h.i,i]
                echam.rel[j-nproxies,i] <- echam$ensmean[h.i,i]
                
              }
            }
          }
          }

       
        
        sigma <- array(sigma,dim=c(dim(sigma)[1]*2))
        sigmab <- array(sigmab,dim=c(dim(sigmab)[1]*2))
        analysis.rel <- array(analysis.rel,dim=c(dim(analysis.rel)[1]*2))
        echam.rel <- array(echam.rel,dim=c(dim(echam.rel)[1]*2))
        calibrateio <- array(calibrate$data[which(calibrate$sour!="prox"),],dim=c(dim(calibrate$data[which(calibrate$sour!="prox"),])[1]*2))
        
        if (docum&!instrumental){
          doctemp <- which(calibrate$sour=="doc"&calibrate$names=="temp2")
          delta[,3] <- (sum(deltaup[doctemp,],na.rm=TRUE)-sum(cbind(Rcal[doctemp],Rcal[doctemp])[which(!is.na(deltaup[doctemp,]))])) / sum(deltadown[doctemp,],na.rm=TRUE)
          
          calibrate$lon <- append(calibrate$lon[calibrate$sour=="doc"],calibrate$lon[calibrate$sour=="doc"]) ## cut calibrate
          calibrate$lat <- append(calibrate$lat[calibrate$sour=="doc"],calibrate$lat[calibrate$sour=="doc"])
          calibrate$names <- append(calibrate$names[calibrate$sour=="doc"],calibrate$names[calibrate$sour=="doc"])
          calibrate$sour <- append(calibrate$sour[calibrate$sour=="doc"],calibrate$sour[calibrate$sour=="doc"])
          
          month <- rep(1:12,each=5)
          year <- rep(syr-1+r,length(sigma))
          station <- rep(1:5,12)
          obsperyear <- rep(1:length(sigma))
          
        }else if (docum & instrumental) { 
          
          tempinst<- which(calibrate$sour=="inst"&calibrate$names=="temp2")
          slpinst <- which(calibrate$sour=="inst"&calibrate$names=="slp")
          doctemp <- which(calibrate$sour=="doc"&calibrate$names=="temp2")
          delta[,1] <- (sum(deltaup[tempinst,],na.rm=TRUE)-sum(cbind(Rcal[tempinst],Rcal[tempinst])[which(!is.na(deltaup[tempinst,]))])) / sum(deltadown[tempinst,],na.rm=TRUE)
          delta[,2] <- (sum(deltaup[slpinst,],na.rm=TRUE)-sum(cbind(Rcal[slpinst],Rcal[slpinst])[which(!is.na(deltaup[slpinst,]))])) / sum(deltadown[slpinst,],na.rm=TRUE)
          delta[,3] <- (sum(deltaup[doctemp,],na.rm=TRUE)-sum(cbind(Rcal[doctemp],Rcal[doctemp])[which(!is.na(deltaup[doctemp,]))])) / sum(deltadown[doctemp,],na.rm=TRUE)
          
            
          calibrate$lon <- array(append(calibrate$lon[calibrate$sour=="inst"],calibrate$lon[calibrate$sour=="doc"]),dim=length(calibrateio)) ## cut calibrate
          calibrate$lat <- array(append(calibrate$lat[calibrate$sour=="inst"],calibrate$lat[calibrate$sour=="doc"]),dim=length(calibrateio))
          calibrate$names <- array(append(calibrate$names[calibrate$sour=="inst"],calibrate$names[calibrate$sour=="doc"]),dim=length(calibrateio))
          calibrate$sour <- array(append(calibrate$sour[calibrate$sour=="inst"],calibrate$sour[calibrate$sour=="doc"]),dim=length(calibrateio))
          
          month1 <- append(rep(1:6,each=(length(sigma)-60)/12),rep(1:6,each=5))
          month2 <- append(rep(7:12,each=(length(sigma)-60)/12),rep(7:12,each=5))
          month <- append(month1,month2)
          year <- rep(syr-1+r,length(sigma))
          station <- array(append(rep(1:((length(sigma)-60)/12),6),rep(1:5,6)),dim=length(calibrateio))
          obsperyear <- rep(1:length(sigma))
        
        }else if (!docum&instrumental){
          tempinst<- which(calibrate$sour=="inst"&calibrate$names=="temp2")
          slpinst <- which(calibrate$sour=="inst"&calibrate$names=="slp")
          delta[,1] <- (sum(deltaup[tempinst,],na.rm=TRUE)-sum(cbind(Rcal[tempinst],Rcal[tempinst])[which(!is.na(deltaup[tempinst,]))])) / sum(deltadown[tempinst,],na.rm=TRUE)
          delta[,2] <- (sum(deltaup[slpinst,],na.rm=TRUE)-sum(cbind(Rcal[slpinst],Rcal[slpinst])[which(!is.na(deltaup[slpinst,]))])) / sum(deltadown[slpinst,],na.rm=TRUE)
          
          calibrate$lon <- append(calibrate$lon[calibrate$sour=="inst"],calibrate$lon[calibrate$sour=="inst"]) ## cut calibrate
          calibrate$lat <- append(calibrate$lat[calibrate$sour=="inst"],calibrate$lat[calibrate$sour=="inst"])
          calibrate$names <- append(calibrate$names[calibrate$sour=="inst"],calibrate$names[calibrate$sour=="inst"])
          calibrate$sour <- append(calibrate$sour[calibrate$sour=="inst"],calibrate$sour[calibrate$sour=="inst"])
          
          month <- rep(1:12,each=(length(sigma)/12))
          year <- rep(syr-1+r,length(sigma))
          station <- rep(1:(length(sigma)/12),12)
          obsperyear <- rep(1:length(sigma))
        }
          
      
        
        if (cutoffoutliers) {
          q<-quantile(sigma,0.95,na.rm=TRUE) ## here we could also add a cutoff value instead of quant.

          analysis.rel[which(sigma>q)] <- NA
          echam.rel[which(sigma>q)] <- NA
          calibrateio[which(sigma>q)] <- NA
          calibrate$lon[which(sigma>q)] <- NA
          calibrate$lat[which(sigma>q)] <- NA
          calibrate$names[which(sigma>q)] <- NA
          calibrate$sour[which(sigma>q)] <- NA
          month[which(sigma>q)] <- NA
          year[which(sigma>q)] <- NA
          station[which(sigma>q)] <- NA
          obsperyear[which(sigma>q)] <- NA
          sigmab[which(sigma>q)] <- NA
          sigma[which(sigma>q)] <- NA
        }
        
       
        all<-list(station=station, calibrate=calibrateio, analysis=analysis.rel, echam=echam.rel, sigma=sigma, sigmab=sigmab, lon=calibrate$lon, lat=calibrate$lat, names=calibrate$names, sour=calibrate$sour, month=month, year=year, obsperyear = obsperyear,delta = delta)
      
        save(all, delta,file=paste0(dataintdir,'../all/',syr-1+r,'.Rdata'))
        
      } ##end of yearly loop
  for (r in 1:length(yr)){
        load(paste0(dataintdir,'../all/',syr-1+r,'.Rdata'))
    if (r==1){
      all1<-all
      delta1<-delta
    }else{
      all1<-mapply(c, all1, all, SIMPLIFY=FALSE)
      delta1<-rbind(delta1,delta)
    }
  }
  # new data is to be included here

  station <- matrix(sapply(all1[[1]],c),nrow=length(sapply(all1[[1]],c)))
  calibrate <- matrix(sapply(all1[[2]],c),nrow=length(sapply(all1[[2]],c)))
  analysis <- matrix(sapply(all1[[3]],c),nrow=length(sapply(all1[[3]],c)))
  echam <- matrix(sapply(all1[[4]],c),nrow=length(sapply(all1[[4]],c)))
  sigma <- matrix(sapply(all1[[5]],c),nrow=length(sapply(all1[[5]],c)))
  sigmab <- matrix(sapply(all1[[6]],c),nrow=length(sapply(all1[[6]],c)))
  lon <- matrix(sapply(all1[[7]],c),nrow=length(sapply(all1[[7]],c)))
  lat <- matrix(sapply(all1[[8]],c),nrow=length(sapply(all1[[8]],c)))
  names <- matrix(sapply(all1[[9]],c),nrow=length(sapply(all1[[9]],c)))
  sour <- matrix(sapply(all1[[10]],c),nrow=length(sapply(all1[[10]],c)))
  month <- matrix(sapply(all1[[11]],c),nrow=length(sapply(all1[[11]],c)))
  year <- matrix(sapply(all1[[12]],c),nrow=length(sapply(all1[[12]],c)))
  obsperyear <- matrix(sapply(all1[[13]],c),nrow=length(sapply(all1[[13]],c)))

  # make a dataframe
archive <- data.frame(station,calibrate,analysis,echam,sigma,sigmab,lon,lat,names,sour,month,year,obsperyear,stringsAsFactors = FALSE)

  ## Results

  sigma <- tapply(archive$sigma[which(archive$sour=="inst")],archive$station[which(archive$sour=="inst")],mean,na.rm=TRUE)
  sigma.doc <- tapply(archive$sigma[which(archive$sour=="doc")],archive$station[which(archive$sour=="doc")],mean,na.rm=TRUE)
 
  Rcalnew[iter,1] <- abs(mean(archive$sigma[which(archive$names=="temp2"&archive$sour=="inst")],na.rm=TRUE)) ## this is only possible if all the lists of names are equal
  Rcalnew[iter,2] <- abs(mean(archive$sigma[archive$names=="slp"&archive$sour=="inst"],na.rm=TRUE))
  Rcalnew[iter,3] <- abs(mean(archive$sigma[archive$sour=="doc"],na.rm=TRUE))
  Pnew[iter,1] <- abs(mean(archive$sigmab[archive$names=="temp2"&archive$sour=="inst"],na.rm=TRUE))
  Pnew[iter,2] <- abs(mean(archive$sigmab[archive$names=="slp"&archive$sour=="inst"],na.rm=TRUE))
  Pnew[iter,3] <- abs(mean(archive$sigmab[archive$sour=="doc"],na.rm=TRUE))
  deltanew[iter,] <- apply(delta1,2,mean,na.rm=TRUE)

} # end of iterary loop


print(paste0('Rcalnew for the years ',syr,'-',eyr,':'))
print(Rcalnew)
print(paste0('Pnew for the years ',syr,'-',eyr,':'))
print(Pnew)
if (calc.delta){
  print(paste0('deltanew for the years ',syr,'-',eyr,':'))
  print(deltanew)
 }

if (length(sigma)==0) {
  sigma<- array(NA, dim=c(278))
}
if (length(sigma.doc)==0) {
  sigma.doc <- array(NA,dim=c(5))
}
  if (errtimeloop){
  if (n==1){
      archive.alltime <- archive
      Rcalnew.time<- Rcalnew[niter,]
      sigma.time <- sigma
      sigma.doctime<- sigma.doc
      delta.time <- deltanew[niter,]
    } else {
      archive.alltime<- rbind(archive.alltime,archive) #contains the fourth iteration sigma of every observation!
      Rcalnew.time <- rbind(Rcalnew.time,Rcalnew[iter,])
      sigma.time <- cbind(sigma.time,sigma)
      sigma.doctime<- cbind(sigma.doctime,sigma.doc)
      delta.time <- rbind(delta.time,deltanew[niter,])
    }
  }
} # end of errtime-loop




print(paste("time for sigma calc:",(proc.time() -ptm)[3]))

##save and load ##
save(calc.delta,initsyr,initeyr,errtime,sigma.time,sigma.doctime,delta.time,archive.alltime, Rcalnew.time, file=paste0(dataintdir,expname,'/sigma.Rdata'))
load(paste0(dataintdir,expname,'/sigma.Rdata'))
##variable description ####

# variable                 dim               Error Variance             averaged over               T/SLP

#                                                                 years      months  stations     
# 
# sigma                   [278]           Observation              yes        yes       no           Both
# tsigmonths.all          [192,12]        Observation              yes        no        no            T
# ssigmonths.all          [86,12]         Observation              yes        no        no           SLP 
# tsigmonths              [12]            Observation              yes        no        yes           T
# ssigmonths              [12]            Observation              yes        no        yes          SLP
# Rcalnew.time            [nloops,2]      Observation              yes        yes       yes          Both



 
 
 
 
## plot and see for instruments#######################################################################################################


  ## prepplots
if (any(!is.na(Rcalnew.time[,1]))){
  sigma.T <- tapply(archive.alltime$sigma[archive.alltime$names=="temp2"&archive.alltime$sour=="inst"],archive.alltime$station[archive.alltime$names=="temp2"&archive.alltime$sour=="inst"],mean,na.rm=TRUE)
  
  tsigmonths <- tapply(archive.alltime$sigma[archive.alltime$names=="temp2"&archive.alltime$sour=="inst"],archive.alltime$month[archive.alltime$names=="temp2"&archive.alltime$sour=="inst"],mean,na.rm=TRUE)
  
  tsigmonths.all <-  tapply(archive.alltime$sigma[archive.alltime$names=="temp2"&archive.alltime$sour=="inst"],archive.alltime$obsperyear[archive.alltime$names=="temp2"&archive.alltime$sour=="inst"],mean,na.rm=TRUE)
  tsigmonths.all <- array(tsigmonths.all,dim=c(length(tsigmonths.all)/12,12))
  
} 
if (any(!is.na(Rcalnew.time[,2]))){
  sigma.SLP <- tapply(archive.alltime$sigma[archive.alltime$names=="slp"&archive.alltime$sour=="inst"],archive.alltime$station[archive.alltime$names=="slp"&archive.alltime$sour=="inst"],mean,na.rm=TRUE)
  
  ssigmonths <- tapply(archive.alltime$sigma[archive.alltime$names=="slp"&archive.alltime$sour=="inst"],archive.alltime$month[archive.alltime$names=="slp"&archive.alltime$sour=="inst"],mean,na.rm=TRUE)
  
  ssigmonths.all <-  tapply(archive.alltime$sigma[archive.alltime$names=="slp"],archive.alltime$obsperyear[archive.alltime$names=="slp"],mean,na.rm=TRUE)
  
  ssigmonths.all <- array(ssigmonths.all,dim=c(length(ssigmonths.all)/12,12))
}
if (any(!is.na(Rcalnew.time[,1:2]))){
  names <- archive.alltime$names[which(archive.alltime$sour=="inst")][1:278]
  lon <- archive.alltime$lon[which(archive.alltime$sour=="inst")][1:278]
  lat <- archive.alltime$lat[which(archive.alltime$sour=="inst")][1:278]
  sour <-archive.alltime$sour[which(archive.alltime$sour=="inst")][1:278]
  
  nrobs.temp<- array(NA, dim=c(initeyr-initsyr+1))
  nrobs.slp<- nrobs.temp
  for(r in initsyr:initeyr){
    nrobs.temp[r-initsyr+1] <- length(which(!is.na(tapply(archive.alltime$calibrate[which(archive.alltime$year==r&archive.alltime$sour=="inst")],archive.alltime$station[which(archive.alltime$year==r&archive.alltime$sour=="inst")],mean,na.rm=TRUE))[1:192]))
    nrobs.slp[r-initsyr+1] <- length(which(!is.na(tapply(archive.alltime$calibrate[which(archive.alltime$year==r&archive.alltime$sour=="inst")],archive.alltime$station[which(archive.alltime$year==r&archive.alltime$sour=="inst")],mean,na.rm=TRUE))[193:278]))
  }
}


  
  
#Evolution of Sigma² ###################################################################################################
#for T###################################################################################################
  if (any(!is.na(Rcalnew.time[,1]))){

 png(filename=paste0(plotintdir,expname,'/Evolution_T.png'),width=1400,height = 800)
par(mfrow=c(2, 1))
op <- par(mar = c(5,4,4,4) + 0.1)
x <- seq(from=initsyr, to=initeyr, by=errtime)+errtime/2## I assume that the first 192 obs are always T and the last 89 are always SLP
x2<- seq(initsyr,initeyr)
plot(x, Rcalnew.time[,1], main=expression(paste("Instrument's Evolution of ",sigma[T],"²"," Method: Desroziers et al.  Intervall: 20 yrs Iterations: 4")), sub="",xlab="Year" , ylab=expression(paste(sigma[T],"²")), xaxp  = c(1650, 2000, 35), cex=1, pch=21, bg="red", ylim=c(0,max(Rcalnew.time[,1],na.rm=TRUE)))
if (calc.delta){
par(new = TRUE)
plot(x, delta.time[,1], type = "l", axes = FALSE,  xlab = "", ylab = "")
axis(side=4, at = pretty(range(delta.time[,1])))
mtext("Inflation Factor", side=4, line=3)
}
par(op)
plot(x2,nrobs.temp,type="l", ylab="Number of observations", xlab="Year" , xaxp  = c(1650, 2000, 35))
dev.off()
}
#for SLP####################################################################################################################################

if (any(!is.na(Rcalnew.time[,2]))){

png(filename=paste0(plotintdir,expname,'/Evolution_SLP.png'),width=1400,height = 800)
par(mfrow=c(2, 1))
op <- par(mar = c(5,4,4,4) + 0.1)
plot(x, Rcalnew.time[,2], main=expression(paste("Evolution of ",sigma[SLP],"²","  Method: Desroziers et al.  Intervall: 20 yrs Iterations: 4")), sub="",xlab="Year" , ylab=expression(paste(sigma[SLP],"²")), xaxp  = c(1650, 2000, 35), cex=1, pch=21, bg="blue", ylim=c(0,max(Rcalnew.time[,2],na.rm=TRUE)))
if (calc.delta) {
par(new = TRUE)
plot(x, delta.time[,2], type = "l", axes = FALSE,  xlab = "", ylab = "")
axis(side=4, at = pretty(range(delta.time[,2])))
mtext("Inflation Factor", side=4, line=3)
}
par(op)
plot(x2,nrobs.slp,type="l", ylab="Number of observations", xlab="Year" , xaxp  = c(1650, 2000, 35))
dev.off()

par(mfrow=c(1,1))
}
# map T, for each time block###################################################################################################
  if (length(which(!is.na(Rcalnew.time[,1])))>1){
  x <- seq(from=initsyr, to=initeyr, by=errtime)+errtime/2
tsigma.time <- abs(sigma.time)[names=="temp2",]
lev <- pretty(tsigma.time[,],13) 
br <- length(lev)
colpal <- two.colors(n=br,start="yellow", end="purple", middle="red", alpha=1.0)
datcol <- colpal[as.numeric(cut(tsigma.time[,],breaks=br))] 
datcol <- array(datcol,dim=dim(tsigma.time))
for (m in 1:((initeyr-initsyr+1)/errtime)) { ##nloops
  png(filename=paste0(plotintdir,expname,'/map_T_',x[m]-(errtime/2),"-",x[m]+(errtime/2-1),'.png'),width=1400,height = 800)
  plot(lon[names=="temp2"][which(!is.na(tsigma.time[,m]))], lat[names=="temp2"][which(!is.na(tsigma.time[,m]))],
       xlim=c(-180,180),ylim=c(-90,90),
       cex=1.2,                   # cex=point size~precip
       pch=15,                    # point type: 15 is square with possible color filling
       col=datcol[,m][which(!is.na(tsigma.time[,m]))],
       xlab="Longitude",ylab="Latitude",main=paste0("Sigma^2 of Temperature: years ", x[m]-10,"-",x[m]+9))# point fill color
  map("world",interior=F,add=T,ty='l',col='black',xlim=c(-180,180),ylim=c(-90,90))
  legend("bottomleft", inset=0.01, as.character(lev),fill=colpal,bty="o",
         bg="white",box.col="white",box.lwd=0,cex=0.7)
  
  dev.off()
}
}
#map SLP for each time block###################################################################################################
  if (length(which(!is.na(Rcalnew.time[,2])))>1){
x <- seq(from=initsyr, to=initeyr, by=errtime)+errtime/2
ssigma.time <- abs(sigma.time)[names=="slp",]
lev <- pretty(ssigma.time,13) 
br <- length(lev)
colpal <- two.colors(n=br,start="yellow", end="purple", middle="red", alpha=1.0)
datcol <- colpal[as.numeric(cut(ssigma.time,breaks=br))] 
datcol <- array(datcol,dim=dim(ssigma.time))
for (m in 1:((initeyr-initsyr+1)/errtime)) { ##nloops
  png(filename=paste0(plotintdir,expname,'/map_SLP_',x[m]-(errtime/2),"-",x[m]+(errtime/2-1),'.png'),width=1400,height = 800)
  plot(lon[names=="slp"][which(!is.na(ssigma.time[,m]))], lat[names=="slp"][which(!is.na(ssigma.time[,m]))],
       xlim=c(-180,180),ylim=c(-90,90),
       cex=1.2,                   # cex=point size~precip
       pch=15,                    # point type: 15 is square with possible color filling
       col=datcol[,m][which(!is.na(ssigma.time[,m]))],
       xlab="Longitude",ylab="Latitude",main=paste0("Sigma^2 of SLP: years ", x[m]-10,"-",x[m]+9))# point fill color
  map("world",interior=F,add=T,ty='l',col='black',xlim=c(-180,180),ylim=c(-90,90))
  legend("bottomleft", inset=0.01, as.character(lev),fill=colpal,bty="o",
         bg="white",box.col="white",box.lwd=0,cex=0.7)
  
  dev.off()
}
}
## highlatlowlat etc###################################################################################################
  if (any(!is.na(Rcalnew.time[,1]))&any(!is.na(Rcalnew.time[,2]))){
lat.crit <- 45

highlat<- which(archive.alltime$names=="temp2"&archive.alltime$sour=="inst"&(archive.alltime$lat>lat.crit | archive.alltime$lat<(-lat.crit)))
lowlat <- which(archive.alltime$names=="temp2"&archive.alltime$sour=="inst"&(archive.alltime$lat<lat.crit & archive.alltime$lat>(-lat.crit)))
tsig.highlat <- tapply(archive.alltime$sigma[highlat],archive.alltime$month[highlat],mean,na.rm=TRUE)
tsig.lowlat <- tapply(archive.alltime$sigma[lowlat],archive.alltime$month[lowlat],mean,na.rm=TRUE)

highlat<- which(archive.alltime$names=="slp"&archive.alltime$sour=="inst"&(archive.alltime$lat>lat.crit | archive.alltime$lat<(-lat.crit)))
lowlat <- which(archive.alltime$names=="slp"&archive.alltime$sour=="inst"&(archive.alltime$lat<lat.crit & archive.alltime$lat>(-lat.crit)))
ssig.highlat <- tapply(archive.alltime$sigma[highlat],archive.alltime$month[highlat],mean,na.rm=TRUE)
ssig.lowlat <- tapply(archive.alltime$sigma[lowlat],archive.alltime$month[lowlat],mean,na.rm=TRUE)

png(filename=paste0(plotintdir,expname,'/monthly&lat.png'),width=1400,height = 800)

old.par <- par(mfrow=c(2, 2))

plot(tsig.highlat,pch=21, bg="pink",xlab="Month",ylab=expression(paste(sigma[T],"²")),ylim=c(0,range(tsig.highlat)[2]),cex.lab=1.5)
lines(tsig.highlat)
points(tsigmonths, pch=21, bg="black")
lines(tsigmonths)
points(tsig.lowlat,pch=21, bg="violet")
lines(tsig.lowlat)

legend(
  "topright", 
  pch = 21,
  pt.bg=c("pink", "black","violet"), 
  legend = c("Stations above 45 and below -45","All Stations", "Stations between 45 and  -45"))

plot(ssig.highlat,pch=21, bg="lightblue",xlab="Month",ylab=expression(paste(sigma[SLP],"²")),ylim=c(0,range(ssig.highlat)[2]),cex.lab=1.5)
lines(ssig.highlat)
points(ssigmonths, pch=21, bg="black")
lines(ssigmonths)
points(ssig.lowlat,pch=21, bg="darkblue")
lines(ssig.lowlat)
legend(
  "topright", 
  pch = 21,
  pt.bg=c("lightblue","black", "darkblue"), 
  legend = c("Stations above 45 and below -45","All Stations", "Stations between 45 and  -45"))

plot(archive.alltime$sigma[names=="temp2"],archive.alltime$lat[names=="temp2"],col="red",xlab=expression(paste(sigma,"²")),ylab="Latitude",xlim=range(archive.alltime$sigma[names=="temp2"],na.rm=TRUE),cex.lab=1.5)
points(archive.alltime$sigma[names=="slp"],archive.alltime$lat[names=="slp"],col="blue")
legend("bottomright",pch = 1,col=c("red", "blue"),legend = c("Temp", "Slp"))

par(old.par)
dev.off()

par(mfrow=c(1,1))
}
##map T sigma for all years###################################################################################################
  if (any(!is.na(Rcalnew.time[,1]))){
png(filename=paste0(plotintdir,expname,'/map_T.png'),width=1400,height = 800)

lev <- pretty(sigma.T,13)
br <- length(lev)
colpal <- two.colors(n=br,start="yellow", end="purple", middle="red", alpha=1.0)
datcol <- colpal[as.numeric(cut(sigma.T,breaks=br))][which(!is.na(sigma.T))]
plot(lon[names=="temp2"][which(!is.na(sigma.T))], lat[names=="temp2"][which(!is.na(sigma.T))],
     xlim=c(-180,180),ylim=c(-90,90),
     cex=1.2,                   # cex=point size~precip
     pch=15,                    # point type: 15 is square with possible color filling
     col=datcol,
     xlab="Longitude",ylab="Latitude",main="Sigma^2 of Temperature")                # point fill color
map("world",interior=F,add=T,ty='l',col='black',xlim=c(-180,180),ylim=c(-90,90))
legend("bottomleft", inset=0.01, as.character(lev),fill=colpal,bty="o",
       bg="white",box.col="white",box.lwd=0,cex=0.7)
dev.off()
}
## monthwise WATCH! -> only absolute values ###################################################################################################
  if (any(!is.na(Rcalnew.time[,1]))){
tsigmonths.all <- abs(tsigmonths.all)
lev <- pretty(tsigmonths.all[,],13) 
br <- length(lev)
colpal <- two.colors(n=br,start="yellow", end="purple", middle="red", alpha=1.0)
datcol <- colpal[as.numeric(cut(tsigmonths.all[,],breaks=br))] 
datcol <- array(datcol,dim=dim(tsigmonths.all))
for (m in 1:12) {
  png(filename=paste0(plotintdir,expname,'/map_T_month_',m,'.png'),width=1400,height = 800)
  
  plot(lon[names=="temp2"][which(!is.na(tsigmonths.all[,m]))], lat[names=="temp2"][which(!is.na(tsigmonths.all[,m]))],
       xlim=c(-180,180),ylim=c(-90,90),
       cex=1.2,                   # cex=point size~precip
       pch=15,                    # point type: 15 is square with possible color filling
       col=datcol[,m][which(!is.na(tsigmonths.all[,m]))],
       xlab="Longitude",ylab="Latitude",main=paste0("Sigma^2 of Temperature: Month ", m))# point fill color
  map("world",interior=F,add=T,ty='l',col='black',xlim=c(-180,180),ylim=c(-90,90))
  legend("bottomleft", inset=0.01, as.character(lev),fill=colpal,bty="o",
         bg="white",box.col="white",box.lwd=0,cex=0.7)
  
  dev.off()
}
}
##map SLP sigma for all years###################################################################################################
  if (any(!is.na(Rcalnew.time[,2]))){
png(filename=paste0(plotintdir,expname,'/map_SLP.png'),width=1400,height = 800)

lev <- pretty(sigma.SLP,13)
br <- length(lev)
colpal <- two.colors(n=br,start="yellow", end="purple", middle="red", alpha=1.0)
datcol <- colpal[as.numeric(cut(sigma.SLP,breaks=br))][which(!is.na(sigma.SLP))]
plot(lon[names=="slp"][which(!is.na(sigma.SLP))], lat[names=="slp"][which(!is.na(sigma.SLP))],
     xlim=c(-180,180),ylim=c(-90,90),
     cex=1.2,                   # cex=point size~precip
     pch=15,                    # point type: 15 is square with possible color filling
     col=datcol,
     xlab="Longitude",ylab="Latitude",main="Sigma² of SLP")                # point fill color
map("world",interior=F,add=T,ty='l',col='black',xlim=c(-180,180),ylim=c(-90,90))
legend("bottomleft", inset=0.01, as.character(lev),fill=colpal,bty="o",
       bg="white",box.col="white",box.lwd=0,cex=0.7)

dev.off()
}
## monthwise WATCH! -> only absolute values ###################################################################################################
  if (any(!is.na(Rcalnew.time[,2]))){
ssigmonths.all <- abs(ssigmonths.all)
lev <- pretty(ssigmonths.all[,],13)
br <- length(lev)
colpal <- two.colors(n=br,start="yellow", end="purple", middle="red", alpha=1.0)
datcol <- colpal[as.numeric(cut(ssigmonths.all[,],breaks=br))]
datcol <- array(datcol,dim=dim(ssigmonths.all))

for (m in 1:12) {
  png(filename=paste0(plotintdir,expname,'/map_SLP_month_',m,'.png'),width=1400,height = 800)
  plot(lon[names=="slp"][which(!is.na(ssigmonths.all[,m]))], lat[names=="slp"][which(!is.na(ssigmonths.all[,m]))],
       xlim=c(-180,180),ylim=c(-90,90),
       cex=1.2,                   # cex=point size~precip
       pch=15,                    # point type: 15 is square with possible color filling
       col=datcol[,m][which(!is.na(ssigmonths.all[,m]))],
       xlab="Longitude",ylab="Latitude",main=paste0("Sigma^2 of SLP: Month ", m))                # point fill color
  map("world",interior=F,add=T,ty='l',col='black',xlim=c(-180,180),ylim=c(-90,90))
  legend("bottomleft", inset=0.01, as.character(lev),fill=colpal,bty="o",
         bg="white",box.col="white",box.lwd=0,cex=0.7)
  
  dev.off()
}
}

## plot and see for documents#######################################################################################################


    ##prepplots
if (any(!is.na(Rcalnew.time[,3]))){
  sigma.DOC <- tapply(archive.alltime$sigma[archive.alltime$sour=="doc"],archive.alltime$station[archive.alltime$sour=="doc"],mean,na.rm=TRUE)
  dsigmonths <- tapply(archive.alltime$sigma[archive.alltime$sour=="doc"],archive.alltime$month[archive.alltime$sour=="doc"],mean,na.rm=TRUE)
  dsigmonths.all <-  tapply(archive.alltime$sigma[archive.alltime$sour=="doc"],archive.alltime$obsperyear[archive.alltime$sour=="doc"],mean,na.rm=TRUE)
  dsigmonths.all <- array(dsigmonths.all,dim=c(length(dsigmonths.all)/12,12))
  
  names <- archive.alltime$names[which(archive.alltime$sour=="doc")][1:60]
  lon <- archive.alltime$lon[which(archive.alltime$sour=="doc")][1:60]
  lat <- archive.alltime$lat[which(archive.alltime$sour=="doc")][1:60]
  sour <-archive.alltime$sour[which(archive.alltime$sour=="doc")][1:60]
  
  nrobs.doc <- array(NA, dim=c(initeyr-initsyr+1))
  for(r in initsyr:initeyr){
    nrobs.doc[r-initsyr+1] <- length(which(!is.na(tapply(archive.alltime$calibrate[which(archive.alltime$year==r&archive.alltime$sour=="doc")],archive.alltime$station[which(archive.alltime$year==r&archive.alltime$sour=="doc")],mean,na.rm=TRUE))))
  }

#Evolution of Sigma² ###################################################################################################

  png(filename=paste0(plotintdir,expname,'/Evolution_DOC.png'),width=1400,height = 800)
  par(mfrow=c(2, 1))
  op <- par(mar = c(5,4,4,4) + 0.1)
  x <- seq(from=initsyr, to=initeyr, by=errtime)+errtime/2## I assume that the first 192 obs are always T and the last 89 are always SLP
  x2<- seq(initsyr,initeyr)
  plot(x, Rcalnew.time[,3], main=expression(paste("Instrument's Evolution of ",sigma[doc],"²"," Method: Desroziers et al.  Intervall: 20 yrs Iterations: 4")), sub="",xlab="Year" , ylab=expression(paste(sigma[T],"²")), xaxp  = c(1650, 2000, 35), cex=1, pch=21, bg="red", ylim=c(0,max(Rcalnew.time[,3],na.rm=TRUE)))
  if (calc.delta) {
  par(new = TRUE)
  plot(x, delta.time[,3], type = "l", axes = FALSE,  xlab = "", ylab = "")
  axis(side=4, at = pretty(range(delta.time[,3])))
  mtext("Inflation Factor", side=4, line=3)
  }
  par(op)
  plot(x2,nrobs.doc,type="l", ylab="Number of observations", xlab="Year" , xaxp  = c(1650, 2000, 35))
  par(mfrow=c(1, 1))
  dev.off()
  
  par(mfrow=c(1,1))

# map doc, for each time block###################################################################################################
  if (length(which(!is.na(Rcalnew.time[,3])))>1){
  x <- seq(from=initsyr, to=initeyr, by=errtime)+errtime/2
  tsigma.doctime <- abs(sigma.doctime)
  lev <- pretty(tsigma.doctime[,],13) 
  br <- length(lev)
  colpal <- two.colors(n=br,start="yellow", end="purple", middle="red", alpha=1.0)
  datcol <- colpal[as.numeric(cut(tsigma.doctime[,],breaks=br))] 
  datcol <- array(datcol,dim=dim(tsigma.doctime))
  for (m in 1:((initeyr-initsyr+1)/errtime)) { ##nloops
    png(filename=paste0(plotintdir,expname,'/map_DOC_',x[m]-(errtime/2),"-",x[m]+(errtime/2-1),'.png'),width=1400,height = 800)
    plot(lon[names=="temp2"][which(!is.na(tsigma.doctime[,m]))], lat[names=="temp2"][which(!is.na(tsigma.doctime[,m]))],
         xlim=c(-180,180),ylim=c(-90,90),
         cex=1.2,                   # cex=point size~precip
         pch=15,                    # point type: 15 is square with possible color filling
         col=datcol[,m][which(!is.na(tsigma.doctime[,m]))],
         xlab="Longitude",ylab="Latitude",main=paste0("Sigma² of documents: years ", x[m]-10,"-",x[m]+9))# point fill color
    map("world",interior=F,add=T,ty='l',col='black',xlim=c(-180,180),ylim=c(-90,90))
    legend("bottomleft", inset=0.01, as.character(lev),fill=colpal,bty="o",
           bg="white",box.col="white",box.lwd=0,cex=0.7)
    
    dev.off()
  }
  }
#map doc sigma for all years###################################################################################################
  
  png(filename=paste0(plotintdir,expname,'/map_DOC.png'),width=1400,height = 800)
  
  lev <- pretty(sigma.DOC,13)
  br <- length(lev)
  colpal <- two.colors(n=br,start="yellow", end="purple", middle="red", alpha=1.0)
  datcol <- colpal[as.numeric(cut(sigma.DOC,breaks=br))][which(!is.na(sigma.DOC))]
  plot(lon[which(!is.na(sigma.DOC))], lat[which(!is.na(sigma.DOC))],
       xlim=c(-180,180),ylim=c(-90,90),
       cex=1.2,                   # cex=point size~precip
       pch=15,                    # point type: 15 is square with possible color filling
       col=datcol,
       xlab="Longitude",ylab="Latitude",main="Sigma² of documents")                # point fill color
  map("world",interior=F,add=T,ty='l',col='black',xlim=c(-180,180),ylim=c(-90,90))
  legend("bottomleft", inset=0.01, as.character(lev),fill=colpal,bty="o",
         bg="white",box.col="white",box.lwd=0,cex=0.7)
  dev.off()
  
# monthwise WATCH! -> only absolute values ###################################################################################################
  
  dsigmonths.all <- abs(dsigmonths.all)
  lev <- pretty(dsigmonths.all[,],13) 
  br <- length(lev)
  colpal <- two.colors(n=br,start="yellow", end="purple", middle="red", alpha=1.0)
  datcol <- colpal[as.numeric(cut(dsigmonths.all[,],breaks=br))] 
  datcol <- array(datcol,dim=dim(dsigmonths.all))
  for (m in 1:12) {
    png(filename=paste0(plotintdir,expname,'/map_DOC_month_',m,'.png'),width=1400,height = 800)
    plot(lon[which(!is.na(dsigmonths.all[,m]))], lat[which(!is.na(dsigmonths.all[,m]))],
         xlim=c(-180,180),ylim=c(-90,90),   
         cex=1.2,                   # cex=point size~precip
         pch=15,                    # point type: 15 is square with possible color filling
         col=datcol[,m][which(!is.na(dsigmonths.all[,m]))],
         xlab="Longitude",ylab="Latitude",main=paste0("Sigma² of documents: Month ", m))# point fill color
    map("world",interior=F,add=T,ty='l',col='black',xlim=c(-180,180),ylim=c(-90,90))
    legend("bottomleft", inset=0.01, as.character(lev),fill=colpal,bty="o",
           bg="white",box.col="white",box.lwd=0,cex=0.7)
    
    dev.off()
  }
  }
  
  
  




###################################################################################################################################
## extrapolation (Wartenburger) ###################################################################################################
###################################################################################################################################

cutoff=1500 # in km

## go

syr=initsyr
eyr=initeyr
yr <- syr:eyr
nyr=(eyr-syr)+1
source('EnSRF_switches.R')
source('EnSRF_functions.R')
expname2=paste0("extrapolation_",syr,"-",eyr,"_cutoff_",cutoff,"_errtime_",errtimeloop)
dir.create(paste0(dataintdir,expname,"/",expname2))
dir.create(paste0(plotintdir,expname,"/",expname2))

instrumental=TRUE
docum=FALSE ## no docum used because too complicated, we concentrate on inst. 
real_proxies=FALSE

ptm <- proc.time()

for (z in 1:2) {
  if (z==1){
    Temp=TRUE ## if false then SLP is taken, also (!!!): change colours and ylab from T <-> SLP in last plot
    name="temp2"
    savein="T"
    print("Calculating Temperature")
    }else { 
    Temp=FALSE
    name="slp"
    savein="SLP"
    print("Calculating SLP")
    }
nyr=(initeyr-initsyr)+1  

k=0
## a function to calculate the distance between to coordinate points


# dist <- function (lat,lon,lat1,lon1) {
#   R = 6371 # Radius of the earth in km
#   dLat = deg2rad(lat-lat1) #deg2rad below
#   dLon = deg2rad(lon-lon1)
#   a=sin(dLat/2) * sin(dLat/2) + cos(deg2rad(lat1)) * cos(deg2rad(lat)) * sin(dLon/2) * sin(dLon/2)
#   c = 2 * atan2(sqrt(a), sqrt(1-a))
#   d = R * c
# }
# deg2rad <- function(deg) {(deg * pi) / (180)}# commented because compute_dist already exist in functions -> (lon,lat,lon0,lat0)
nona <- function(matr){
  matr[which(matr==0)] <- NA
  matr<-which(!is.na(matr))
  matr<-length(matr)
}
smallerequal<- function(matr) {
  if (matr>cutoff){
    matr=NA
  } else if (matr<1){
    matr=NA
  } else {matr=matr}
  
}

#here we load any year and calculate the distances of each station to each station
load(paste0(ana_dir,'analysis_',1940,'_2ndgrid.Rdata'))
rm(analysis.abs,analysis.anom,echam.abs,echam.anom)

calibrate$lat <- calibrate$lat[which(calibrate$sour=="inst"&calibrate$names==name)] # take only the instrument data and only
lat <- calibrate$lat[1:(length(calibrate$lat)/6)] # one month because all stations are in one month

calibrate$lon <- calibrate$lon[which(calibrate$sour=="inst"&calibrate$names==name)]
lon <- calibrate$lon[1:(length(calibrate$lon)/6)]
len<-length(lon)## if Temp==TRUE this is 192 else 86
distances<-array(NA,dim=c(len,len))


for (s in 1:len){
  # im Umkreis von x kilometer:    (x - center_x)^2 + (y - center_y)^2 < radius^2
  lat0 <- lat[s]
  lon0 <- lon[s]
  
  distances[,s]<-mapply(compute_dist,lon,lat,lon0,lat0)
}
#here we only take the ones that are within the cutoff distance

wewant<- sapply(distances,smallerequal)
wewant<- array(wewant,dim=c(len,len))
howmany <- apply(wewant,2,nona) # this gives a list indicating how many stations are surrounding it

coef.alltime <- array(NA,dim=c(nyr/errtime,9))
nrobs <- array(NA,dim=c(nyr/errtime,3))
for (n in seq(from=1, to=nyr, by=errtime)){ # calculates the error of errtime year until eyr to
  syr=initsyr                               # come up how the error changes over time
  eyr=initeyr
  nyr=(eyr-syr)+1
  syr=syr-1+n
  eyr=syr-1+errtime
  nyr=eyr-syr+1
  yr <- syr:eyr
  print(paste0("years ",syr,'-',eyr))
  k=k+1  
  if (syr < 1660) {
    next      
  } 
  
  
  
  ## here we go through the years taking the stations that have data and taking 
  ## their surrounding stations with data
  
  
  for (r in 1:length(yr)) {
    print(paste0("loading year: ",syr-1+r))
    load(paste0(ana_dir,'analysis_',syr-1+r,'_2ndgrid.Rdata'))
    rm(analysis.abs,analysis.anom,echam.abs,echam.anom)
    calibrate$data <- calibrate$data[which(calibrate$sour=="inst" & calibrate$names==name),] #only instrumental data
    calibrate$data <- array(calibrate$data,dim=c((dim(calibrate$data)[1]*2))) # rearrange data
    
    month <- rep(1:12,each=(dim(calibrate$data)/12)) #define which rows belong to which month
    d <- array(NA,dim=c(len,12,len))
    diff <- array(NA,dim=c(len,12,len))
    
    
    for (m in 1:12){                                  # monthly loop
      data <- calibrate$data[which(month==m)]                # data is only one month at a time
      
      active <- array(NA,dim=c(len,len))
      
      
      for (s in 1:length(data)){                      # loop over all stations
        if (!is.na(data[s])) {                        # loop only executed if middle station is not NA
          active[,s][which(!is.na(wewant[,s]))] <- data[which(!is.na(wewant[,s]))] # here the data for all stations within perimeter are selected (remember distance to self is NA)
          active[s,s] <- data[s]                                                   # data of station itself is added
          d[,m,s][which(!is.na(active[,s]))] <- wewant[which(!is.na(active[,s])),s] # here the distance to the station is taken but only the ones that have data
          d[s,m,s]<-0                                                               # now we see: if in the diagonal ther is a 0 then the station has a value but if it's NA it has none
          diff[,m,s] <- active[,s]-active[s,s]                                      # here we substract the middle station's value from the surrounding stations's
          diff[s,m,s] <- NA                                                        # the diff would be zero but make it NA so we don't plot it at the end
        }
      }   
      # the diff has a different format: each thrid dimension is a station with the coloumns 
      # being the months, this way it is easier later on to calculate the variance of the diff
      
      
      ## in the resulting matrix active.tot the third dimension is for each month
      ## one coloumn represents the "middle" station and the rows all the stations
      ## surrounding it. In the diagonal (1,1 ; 2,2 ..) you have the value
      ## of the station itself. the other values in the row are of the stations
      ## within the perimeter.
      ## Now if for instance in the first coloumn the third row [3,1] has value, this 
      ## means the third row is within the perimeter AND has a value (is not NA)
      ## Therefore the value at [3,3] which is the value of the third station itself
      ## should be exactly the same, which it is :)
      ## d.tot is the same matrix but instead of the data it has the distance to
      ## the station of the coloumn. There the distance to the station itself, 
      ## the diagonal is zero, but I changed it to NA
      
    } # end of monthly loop
    
    # add calculation of deltasigma² here 
    # produce a file with 3 dimensions: first coloumn for deltasigma² of first station to all its 
    # surrounding stations in the perimeter. Second coloumn the corresponding distances. There should
    # 278 of these matrices in the third dimension, for each stations one. 
    
    if (r==1) {
      diff.tot <- diff
      dist.tot <- d
    }else{
      diff.tot <- abind(diff.tot,diff,along=2)
      dist.tot <- abind(dist.tot,d, along=2) 
    } # diff.tot has the months of the one station as coloumns
    # the rows are the differences to the surrounding stations 
    # and the 278 "middle" stations are in the third dimension
  } #end of yearly loop  
  for (r in 1:(dim(diff.tot)[2]/12))  {
    diff.winter<-diff.tot[,c((1+12*(r-1)):(1+12*(r-1)+5)),]
    diff.summer<-diff.tot[,c((7+12*(r-1)):(7+12*(r-1)+5)),]
  }
  sigmadiff.winter <- apply(diff.winter,c(1,3),var,na.rm=TRUE)
  sigmadiff.summer <- apply(diff.summer,c(1,3),var,na.rm=TRUE)
  sigmadiff.tot <- apply(diff.tot,c(1,3),var,na.rm=TRUE)
  
  ### Fits plotten und y-Achsenabschnitte speichern ####
  
  coef.summer<-array(NA,dim=c(len,3))
  coef.winter<-array(NA,dim=c(len,3))
  coef.tot<-array(NA,dim=c(len,3))
  
  for (s in 1:len){ #plot sigmadiff.summer and winter
    x <- wewant[,s][which(!is.na(sigmadiff.summer[,s]))]*wewant[,s][which(!is.na(sigmadiff.summer[,s]))] #if x is length=1 then none of the surrounding stations have data or it has no surrounding stations
    y <- sigmadiff.summer[,s][which(!is.na(sigmadiff.summer[,s]))] 
    x1 <- wewant[,s][which(!is.na(sigmadiff.winter[,s]))]*wewant[,s][which(!is.na(sigmadiff.winter[,s]))] #if x is length=1 then none of the surrounding stations have data or it has no surrounding stations
    y1 <- sigmadiff.winter[,s][which(!is.na(sigmadiff.winter[,s]))]  
    x2 <- wewant[,s][which(!is.na(sigmadiff.tot[,s]))]*wewant[,s][which(!is.na(sigmadiff.tot[,s]))]
    y2 <- sigmadiff.tot[,s][which(!is.na(sigmadiff.tot[,s]))]
    
    if ( length(x2)>4 & length(y2)>4){
      ylim <- c(0,max(sigmadiff.tot),na.rm=TRUE)
      fit2 <- lm(y2 ~ x2)
      coef.tot[s,]<-c("value"=coef(fit2)[1],"2.5"=confint(fit2)[1,1],"97.5"=confint(fit2)[1,2])
      if (errtimeloop==FALSE) {
        # png(filename=paste0(plotintdir,expname,'/',expname2,'/extrapolation','_',savein,'_','/station',s,'.png'),width=1400,height = 800)
        plot(x2,y2, main=paste0("Middle station: ",s),xlab="Distance²", ylab=expression(paste(sigma[Delta],"²")), cex=1, pch=21, bg="black", xlim=c(0,cutoff*cutoff),ylim=ylim)
        abline(fit2, col="red")
        legend(
          "topleft",
          pch = 21,
          pt.bg=c("black", "red"),
          legend = c("Surrounding Stations","Regression"))
        # dev.off()
      }
    }
    
    if ( length(x)>4 & length(y)>4 & length(x1)>4 & length(y1)>4){  
      ylim <- c(0,max(cbind(sigmadiff.summer,sigmadiff.winter),na.rm=TRUE))
      fit <- lm(y ~ x)
      fit1 <- lm(y1 ~ x1)
      coef.summer[s,]<-c("value"=coef(fit)[1],"2.5"=confint(fit)[1,1],"97.5"=confint(fit)[1,2])
      coef.winter[s,]<-c("value"=coef(fit1)[1],"2.5"=confint(fit1)[1,1],"97.5"=confint(fit1)[1,2])
      if (errtimeloop==FALSE) {
        # png(filename=paste0(plotintdir,expname,'/',expname2,'/',extrapolation,'_',savein,'_','/station',s,'.png'),width=1400,height = 800)
        plot(x,y, main=paste0("Summer - Middle station: ",s),xlab="Distance²", ylab=expression(paste(sigma[Delta],"²")), cex=1, pch=21, bg="green", xlim=c(0,cutoff*cutoff),ylim=ylim)
        points(x1, y1, cex=1, pch=21, bg="blue", xlim=c(0,cutoff*cutoff),ylim=ylim)
        abline(fit, col="green")
        abline(fit1, col="blue")
        legend(
          "topleft", 
          pch = 21,
          pt.bg=c("green", "blue"), 
          legend = c("Summer","Winter"))
        # dev.off()
      }
    }
  }
  # save(coef.summer,coef.winter, coef.tot, sigmadiff.winter, sigmadiff.summer,sigmadiff.tot, file=paste0(dataintdir,expname,'/',expname2,'/extrapolation_',savein,'_data_',syr,'-',eyr,'.Rdata'))
  
  
  
  
  
  if (errtimeloop==TRUE) {
    coef.alltime[k,1:3]<- cbind("Summer"=mean(coef.summer[,1],na.rm=TRUE),"2.5"=mean(coef.summer[,2],na.rm=TRUE),"97.5"=mean(coef.summer[,3],na.rm=TRUE))
    coef.alltime[k,4:6]<- cbind("Winter"=mean(coef.winter[,1],na.rm=TRUE),"2.5"=mean(coef.winter[,2],na.rm=TRUE),"97.5"=mean(coef.winter[,3],na.rm=TRUE))
    coef.alltime[k,7:9]<- cbind("Tot"=mean(coef.tot[,1],na.rm=TRUE),"2.5"=mean(coef.tot[,2],na.rm=TRUE),"97.5"=mean(coef.tot[,3],na.rm=TRUE))
    nrobs[k,1] <- length(which(!is.na(apply(sigmadiff.summer,1,mean,na.rm=TRUE))))
    nrobs[k,2] <- length(which(!is.na(apply(sigmadiff.winter,1,mean,na.rm=TRUE))))
    nrobs[k,3]<-length(which(!is.na(apply(sigmadiff.tot,1,mean,na.rm=TRUE))))
  } 
  
} # end of errtime loop

if (errtimeloop==TRUE&length(which(nrobs!=0))!=0&length(which(!is.na(coef.alltime)))!=0){
  x <- seq(from=initsyr, to=initeyr, by=errtime)+errtime/2
  # png(filename=paste0(plotintdir,expname,'/',expname2,'/extrapolation_',savein,'_Evolution.png'),width=1400,height = 800)
  op <- par(mar = c(5,4,4,4) + 0.1)
  plot(x, coef.alltime[,7], cex=1, pch=21, bg="red", main=expression(paste("Evolution of ", sigma,"²  Method: Wartenburger et al. Intervall: 20 yrs  Cutoff: 1500km")), xlab="Year", ylab=expression(paste(sigma[T],"²")), xaxp  = c(1650, 2000, 35), ylim=c(min(coef.alltime[,7:9],na.rm=TRUE),max(coef.alltime[,7:9],na.rm=TRUE)))
  arrows(x,coef.alltime[,8],x,coef.alltime[,9], code=3, angle=90, length=0.05)
  par(new = TRUE)
  plot(x, nrobs[,3], type = "l", col="black", axes = FALSE,  xlab = "", ylab = "")
  axis(side=4, at = pretty(range(nrobs[,3])))
  mtext("Number of midlle stations", side=4, line=3)
  par(op)
  # legend(
  #  "topleft", 
  #   pch = c(21,2),
  #   pt.bg=c("red","black"), 
  #   legend = c("Value","95% Confidence Intervall"))
  # dev.off()
  
  # data<-data.frame(Year=x,upper=coef.alltime[,9],intercept=coef.alltime[,7],lower=coef.alltime[,8])
  # data1<- melt(data, id = "Year")
  # g1<-ggplot(data1, aes(x=Year,y=value,color=variable))+ geom_point(na.rm=TRUE)  + 
  #   geom_ribbon(data=merge(data1,data),aes(x=Year, ymax=upper, ymin=lower), fill="pink", alpha=.5, colour=NA) + 
  #   xlab("Year") + ylab(expression(paste(sigma[SLP],"²")))+
  #   labs(title=expression(paste("Evolution of ",sigma[SLP],"²")),subtitle= "Intervall: 20 yrs ")+
  #   labs(caption="Method: Wartenburger et al.")+
  #   scale_x_continuous(breaks = round(seq(min(data$Year), max(data$Year), by = 20),1))+
  #   theme(legend.text=element_text(size=12))+
  #   theme(axis.text=element_text(size=10),axis.title=element_text(size=12))
  # data<-data.frame(Year=x,number.of.observations=nrobs[,3])
  # 
  # g2<-ggplot(data,aes(x=Year,y=number.of.observations))+geom_line()+
  #   scale_x_continuous(breaks = round(seq(min(data$Year), max(data$Year), by = 20),1))+
  #   theme(axis.text=element_text(size=10),axis.title=element_text(size=12))
  # # g2<-gtable_add_cols(ggplotGrob(g2),unit(1,"null"))
  # theme_set(theme_minimal())
  # # dev.new()
  # plot_grid(
  #   plot_grid(
  #     g1 + theme(legend.position = "none")
  #     ,g2+ theme(legend.position = "none")
  #     , ncol = 1
  #     , align = "hv")
  #   , plot_grid(
  #     get_legend(g1)
  #     , ggplot()
  #     , ncol =1)
  #   , rel_widths = c(7,3)
  # )
  # ggsave(filename=paste0(plotintdir,expname,'/',expname2,'/extrapolation_',savein,'_Evolution_gg.png'),plot=last_plot(),width=10,height = 7)

  
}
# save(coef.alltime, nrobs, file=paste0(dataintdir,expname,'/',expname2,'/extrapolation_',savein,'_alldata_',initsyr,'-',initeyr,'.Rdata'))

} # end of z loop
print(paste("time for extrapolation:",(proc.time() -ptm)[3]))

# load(paste0(workdir,'../data/analysis/archive/',expname,'/',expname2,'/',extrapolation,'_',savein,'_alldata_',initsyr,'-',initeyr,'.Rdata'))
# load(paste0(workdir,'../data/analysis/format'))
