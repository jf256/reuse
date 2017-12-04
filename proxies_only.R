# for now I only took H.i[1,] to calculate the Covariance, the proxy location of this is somewhere in eastern EU 
# it's marked black on the map plot

rm(list=ls())
workdir='/scratch4/lucaf/reuse/src/'
setwd(workdir)
dataextdir='/mnt/climstor/giub/EKF400/'
dataintdir='/scratch4/joerg/projects/reuse/src/../data/'
ana_dir="/scratch4/joerg/projects/reuse/data/analysis/EKF400_v1.3_corr_echam_clim/"

## input
initsyr=1902
initeyr=1921
nrcores=10

## go
syr=initsyr
eyr=initeyr
yr <- syr:eyr
nyr=(eyr-syr)+1
source('EnSRF_switches.R')
source('EnSRF_functions.R')
# dir.create(paste0("../data/analysis/all/Cov.ana"))
docum=FALSE ## no docum used because too complicated, we concentrate on inst. 
real_proxies=TRUE
instrumental=FALSE


ptm <- proc.time()

 
    # library(doParallel)

    # reserve the number of cores
    registerDoParallel(cores=nrcores)
    
    ## the yearly foreach loop 
    all<-foreach(r=1:length(yr)) %dopar% {
                    
                    
                    load(paste0(ana_dir,'analysis_',syr-1+r,'_2ndgrid.Rdata'))
                    echam <- echam.anom
                    analysis <- analysis.anom
                    rm(echam.abs,echam.anom)
                    
                   
                    
                    
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
                    realprox$var_residu=calibrate$var_residu[calibrate$names=="prox"]
                    realprox$mr=calibrate$mr[calibrate$names=="prox",]
                    realprox$names=calibrate$names[calibrate$names=="prox"]
                    
                    etmp <- echam
                    etmp$lon[is.na(etmp$lon)] <- 0
                    etmp$lat[is.na(etmp$lat)] <- -90
                    
                    print(paste('Calculating Hcal1'))
                    if (instrumental) {
                      Hcal1 <- array(NA,dim=c(dim(inst$data)[1],2))
                      Hcal1 <- compute_Hi_Hredux_sixmonstatevector(inst, etmp, threshold=700)
                    }
                    print(paste('Calculating Hcal2'))
                    if (real_proxies) {
                      Hcal2 <- array(NA,dim=c(dim(realprox$data)[1],14))
                      Hcal2 <- compute_Hi_Hredux_proxy(realprox, etmp, realprox$mr, threshold=700)
                    }
                    if (docum) {
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
                    } else if (!instrumental & !docum & real_proxies) { #new
                      H.i <- array(NA,c(nrow(Hcal2),7))
                      }else { 
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
                    if (!instrumental & !docum & real_proxies){ #new
                      H.i[1:nrow(Hcal2),]<-Hcal2[,c(1,3,5,7,9,11,13)] 
                      Hredux[1:nrow(Hcal2),]<-Hcal2[,c(1,3,5,7,9,11,13)] 
                    }
                    
                    H.i[H.i==0] <- NA
                    Hredux[Hredux==0] <- NA
                    
                 

                      Rcal <- c(temp2=0.9, precip=50, slp=10)[calibrate$names]
                      #        if (avg_prox_per_grid) {Rcal <- Rcal*(1/calibrate$numavg)}
                    

                    
                    Rcal[calibrate$names=="prox"] <- realprox$var_residu/2 
                    # previously used residuals/2 for 1. paper version to give proxies more weight
                    # better delete "/2"
                    # probably should have given instrumentals more error instead!
                    Rcal[calibrate$sour=="doc"] <- 0.25 # equals 0.5 std. dev.
                    
                    
                    
                    nens <- dim(analysis$data)[3]
                    ntim=2 
                    
                    
                    
                    calibrate$data <- calibrate$data[which(calibrate$sour=="prox"),]
                    calibrate$lon <- calibrate$lon[which(calibrate$sour=="prox")]
                    calibrate$lat <- calibrate$lat[which(calibrate$sour=="prox")]
                    calibrate$names <- calibrate$names[which(calibrate$sour=="prox")]
                    calibrate$sour <- calibrate$sour[which(calibrate$sour=="prox")]
                    
                    
                    
                    analysis <- echam
                    for (j in 1:length(which(calibrate$sour=="prox"))){
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
                            
                            HPHR <- as.vector(H %*% PH[h.i,] + Rcal[j])
                            K <- PH / HPHR
                            Ktilde <- K / (1 + sqrt(Rcal[j]/HPHR))
                            analysis$ensmean[,i] <- analysis$ensmean[,i] + K[,1] * as.vector((calibrate$data[j,i] -
                                                                                      H %*% analysis$ensmean[h.i,i]))
                            analysis$data[,i,] <- analysis$data[,i,] - Ktilde %*% H %*% analysis$data[h.i,i,]
                            
                            
                            
                            
                          }
                        }
                      }
                    }
                    
                    


                    
                   save(analysis,echam, file=paste0(workdir,'../data/analysis/all/Cov.ana/',syr-1+r,'.Rdata'))

                  } ##end of yearly loop
    
###################################################################################################################
################--------PREPPLOT-------######################################################################################
###################################################################################################################
    source('EnSRF_switches.R')
    source('EnSRF_functions.R')
    
     calc_prepplot=TRUE
     calc_vali_stat=TRUE
     load_prepplot=TRUE
     
     if (monthly_out) {
      if(calc_prepplot){
        dir.create(paste0("../data/prepplot/",expname))
        dir.create(paste0("../data/prepplot/",expname,'/prepplot_monthly/'))
      }
      prepplotdir=paste0("../data/prepplot/",expname,'/prepplot_monthly/')
    } else {
      if(calc_prepplot){
        dir.create(paste0("../data/prepplot/",expname))
        dir.create(paste0("../data/prepplot/",expname,'/prepplot_seasonal/'))
      }
      prepplotdir=paste0("../data/prepplot/",expname,'/prepplot_seasonal/') 
    }
    if (load_prepplot) {dir.create(paste0("../data/image/",expname))}
    
    if (calc_prepplot) {
      # read yearly analysis, calc indices, cut temp, precip, slp, merge timesteps
      for (cyr in syr:eyr) {
        if ((cyr < 1902)) { #1751)) { #} |  (cyr > 1979)) {
          vali=F                 # switch off prepplot if no vali data selected
        } else {
          vali=T
        }
        #   if ((cyr > 1900) &  (cyr < 1903)) {
        #     vali=F                 # switch off prepplot if no vali data selected
        #   }
        if ((cyr > 1901) & (cyr < 2005)) { #  & (syr > 1900) & (eyr < 2005)) {
          cru_vali=T             # monthly CRU TS3 temp, precip and HADSLP2 gridded instrumentals (1901-2004)
          #    ind_recon=T            # Stefan's reconstructed indices until 1948 and NCAR reanalysis later added to CRU and NCEP
        } else {
          cru_vali=F 
          #    ind_recon=F
        }
        #if ((cyr < 1901) & (cyr > 1749)) { # & (syr < 1901) & (eyr > 1749)) {
        #  recon_vali=T           # seasonal luterbacher, pauling, kuettel recons (1750-1999)
        #} else {
        recon_vali=F
        #}
        # t=syr2
        print(cyr)
        print(paste("recon_vali=",recon_vali))
        print(paste("vali=",vali))
        if (every2grid){
          load(file=paste0('../data/analysis/',expname,'/analysis_',cyr,'_2ndgrid.Rdata'))
        } else {  
          rm(validate,analysis,echam)
          load(file=paste0('../data/analysis/',expname,'/analysis_',cyr,'.Rdata'))
        }
        if (anomaly_assim){
          echam <- echam.abs
          analysis <- analysis.abs
          echam.abs <- NULL
          analysis.abs <- NULL   
        }
        if (write_coor) {
          # proxy data coordinates at each time step
          lonlat_year=list()
          for (t in 1:length(calibrate$time)) {
            tf=as.logical(!is.na(calibrate$data[,t]))
            outstr <- cbind(calibrate$lon[tf],calibrate$lat[tf],calibrate$names[tf])
            outprox <- outstr[outstr[,3]=="prox",]
            outprox2 <- unique(t(apply(outprox[,1:2],1,paste0)))
            outstat <- outstr[outstr[,3]!="prox",]
            outstat2 <- unique(t(apply(outstat[,1:2],1,paste0)))
            write.table(outstat2, file=paste0('../data/coor/stat_coor_',calibrate$time[t],'.csv'),
                        row.names=F,col.names=F)
            write.table(outprox2, file=paste0('../data/coor/prox_coor_',calibrate$time[t],'.csv'),
                        row.names=F,col.names=F)
          }
        }
        
        if (ind_ECHAM) {
          # extract the indices from output 
          # 1. Analysis indices
          ech_ind <- echam.anom
          ech_ind$data <- echam.anom$data[which(echam.anom$names=="DIMI" | echam.anom$names=="z100" | 
                                                  echam.anom$names=="z300" | echam.anom$names=="PWC" | echam.anom$names=="HC" | 
                                                  echam.anom$names=="SJ"),,]
          ech_ind$ensmean <- echam.anom$ensmean[which(echam.anom$names=="DIMI" | echam.anom$names=="z100" |  
                                                        echam.anom$names=="z300" | echam.anom$names=="PWC" | echam.anom$names=="HC" | 
                                                        echam.anom$names=="SJ"),]
          ech_ind$names <- echam.anom$names[which(echam.anom$names=="DIMI" | echam.anom$names=="z100" | 
                                                    echam.anom$names=="z300" | echam.anom$names=="PWC" | echam.anom$names=="HC" | 
                                                    echam.anom$names=="SJ")]
          ech_ind$lon <- echam.anom$lon[which(echam.anom$names=="DIMI" | echam.anom$names=="z100" | 
                                                echam.anom$names=="z300" | echam.anom$names=="PWC" | echam.anom$names=="HC" | 
                                                echam.anom$names=="SJ")]
          ech_ind$lat <- echam.anom$lat[which(echam.anom$names=="DIMI" | echam.anom$names=="z100" | 
                                                echam.anom$names=="z300" | echam.anom$names=="PWC" | echam.anom$names=="HC" | 
                                                echam.anom$names=="SJ")]
          echam$data <- echam$data[which(echam$names!="DIMI" & echam$names!="z100" & echam$names!="z300" &
                                           echam$names!="PWC" & echam$names!="HC" & echam$names!="SJ"),,]
          echam$ensmean <- echam$ensmean[which(echam$names!="DIMI" & echam$names!="z100" & echam$names!="z300" &
                                                 echam$names!="PWC" & echam$names!="HC" & echam$names!="SJ"),]
          echam$names <- echam$names[which(echam$names!="DIMI" & echam$names!="z100" & echam$names!="z300" &
                                             echam$names!="PWC" & echam$names!="HC" & echam$names!="SJ")]
          echam$lon <- echam$lon[which(echam$names!="DIMI" & echam$names!="z100" & echam$names!="z300" &
                                         echam$names!="PWC" & echam$names!="HC" & echam$names!="SJ")]
          echam$lat <- echam$lat[which(echam$names!="DIMI" & echam$names!="z100" & echam$names!="z300" &
                                         echam$names!="PWC" & echam$names!="HC" & echam$names!="SJ")]
          echam.anom$data <- echam.anom$data[which(echam.anom$names!="DIMI" & echam.anom$names!="z100" & 
                                                     echam.anom$names!="z300" & echam.anom$names!="PWC" & 
                                                     echam.anom$names!="HC" & echam.anom$names!="SJ"),,]
          echam.anom$ensmean <- echam.anom$ensmean[which(echam.anom$names!="DIMI" & echam.anom$names!="z100" &
                                                           echam.anom$names!="z300" & echam.anom$names!="PWC" & 
                                                           echam.anom$names!="HC" & echam.anom$names!="SJ"),]
          echam.anom$names <- echam.anom$names[which(echam.anom$names!="DIMI" & echam.anom$names!="z100" & 
                                                       echam.anom$names!="z300" & echam.anom$names!="PWC" & 
                                                       echam.anom$names!="HC" & echam.anom$names!="SJ")]
          echam.anom$lon <- echam.anom$lon[which(echam.anom$names!="DIMI" & echam.anom$names!="z100" & 
                                                   echam.anom$names!="z300" & echam.anom$names!="PWC" & 
                                                   echam.anom$names!="HC" & echam.anom$names!="SJ")]
          echam.anom$lat <- echam.anom$lat[which(echam.anom$names!="DIMI" & echam.anom$names!="z100" & 
                                                   echam.anom$names!="z300" & echam.anom$names!="PWC" & 
                                                   echam.anom$names!="HC" & echam.anom$names!="SJ")]
          ana_ind <- analysis.anom
          ana_ind$data <- analysis.anom$data[which(analysis.anom$names=="DIMI" | analysis.anom$names=="z100" |
                                                     analysis.anom$names=="z300" |
                                                     analysis.anom$names=="PWC" | analysis.anom$names=="HC" | analysis.anom$names=="SJ"),,]
          ana_ind$ensmean <- analysis.anom$ensmean[which(analysis.anom$names=="DIMI" | analysis.anom$names=="z100" | 
                                                           analysis.anom$names=="z300" | analysis.anom$names=="PWC" | analysis.anom$names=="HC" | 
                                                           analysis.anom$names=="SJ"),]
          ana_ind$names <- analysis.anom$names[which(analysis.anom$names=="DIMI" | analysis.anom$names=="z100" | 
                                                       analysis.anom$names=="z300" | analysis.anom$names=="PWC" | analysis.anom$names=="HC" | 
                                                       analysis.anom$names=="SJ")]
          ana_ind$lon <- analysis.anom$lon[which(analysis.anom$names=="DIMI" | analysis.anom$names=="z100" | analysis.anom$names=="z300" |
                                                   analysis.anom$names=="PWC" | analysis.anom$names=="HC" | analysis.anom$names=="SJ")]
          ana_ind$lat <- analysis.anom$lat[which(analysis.anom$names=="DIMI" | analysis.anom$names=="z100" | analysis.anom$names=="z300" |
                                                   analysis.anom$names=="PWC" | analysis.anom$names=="HC" | analysis.anom$names=="SJ")]
          analysis$data <- analysis$data[which(analysis$names!="DIMI" & analysis$names!="z100" & 
                                                 analysis$names!="z300" & analysis$names!="PWC" & analysis$names!="HC" & 
                                                 analysis$names!="SJ"),,]
          analysis$ensmean <- analysis$ensmean[which(analysis$names!="DIMI" & analysis$names!="z100" & 
                                                       analysis$names!="z300" & analysis$names!="PWC" & analysis$names!="HC" & 
                                                       analysis$names!="SJ"),]
          analysis$names <- analysis$names[which(analysis$names!="DIMI" & analysis$names!="z100" & 
                                                   analysis$names!="z300" & analysis$names!="PWC" & analysis$names!="HC" & 
                                                   analysis$names!="SJ")]
          analysis$lon <- analysis$lon[which(analysis$names!="DIMI" & analysis$names!="z100" & analysis$names!="z300" &
                                               analysis$names!="PWC" & analysis$names!="HC" & analysis$names!="SJ")]
          analysis$lat <- analysis$lat[which(analysis$names!="DIMI" & analysis$names!="z100" & analysis$names!="z300" &
                                               analysis$names!="PWC" & analysis$names!="HC" & analysis$names!="SJ")]
          analysis.anom$data <- analysis.anom$data[which(analysis.anom$names!="DIMI" & analysis.anom$names!="z100" & 
                                                           analysis.anom$names!="z300" & analysis.anom$names!="PWC" & analysis.anom$names!="HC" & 
                                                           analysis.anom$names!="SJ"),,]
          analysis.anom$ensmean <- analysis.anom$ensmean[which(analysis.anom$names!="DIMI" & analysis.anom$names!="z100" & 
                                                                 analysis.anom$names!="z300" & analysis.anom$names!="PWC" & analysis.anom$names!="HC" & 
                                                                 analysis.anom$names!="SJ"),]
          analysis.anom$names <- analysis.anom$names[which(analysis.anom$names!="DIMI" & analysis.anom$names!="z100" & 
                                                             analysis.anom$names!="z300" & analysis.anom$names!="PWC" & analysis.anom$names!="HC" & 
                                                             analysis.anom$names!="SJ")]
          analysis.anom$lon <- analysis.anom$lon[which(analysis.anom$names!="DIMI" & analysis.anom$names!="z100" & analysis.anom$names!="z300" &
                                                         analysis.anom$names!="PWC" & analysis.anom$names!="HC" & analysis.anom$names!="SJ")]
          analysis.anom$lat <- analysis.anom$lat[which(analysis.anom$names!="DIMI" & analysis.anom$names!="z100" & analysis.anom$names!="z300" &
                                                         analysis.anom$names!="PWC" & analysis.anom$names!="HC" & analysis.anom$names!="SJ")]
        }
        if (ind_recon) {
          if (vali) {  
            vali_ind <- validate
            vali_ind$data <- validate$data[which(validate$names=="ind_rec_dimi" | validate$names=="ind_rec_z100" |
                                                   validate$names=="ind_rec_z300" | validate$names=="ind_rec_pwc" | 
                                                   validate$names=="ind_rec_hc" | validate$names=="ind_rec_sj"),]
            vali_ind$names <- validate$names[which(validate$names=="ind_rec_dimi" | validate$names=="ind_rec_z100" |
                                                     validate$names=="ind_rec_z300" | validate$names=="ind_rec_pwc" | 
                                                     validate$names=="ind_rec_hc" | validate$names=="ind_rec_sj")]
            vali_ind$lon <- validate$lon[which(validate$names=="ind_rec_dimi" | validate$names=="ind_rec_z100" | 
                                                 validate$names=="ind_rec_z300" | validate$names=="ind_rec_pwc" | 
                                                 validate$names=="ind_rec_hc" | validate$names=="ind_rec_sj")]
            vali_ind$lat <- validate$lat[which(validate$names=="ind_rec_dimi" | validate$names=="ind_rec_z100" | 
                                                 validate$names=="ind_rec_z300" | validate$names=="ind_rec_pwc" | 
                                                 validate$names=="ind_rec_hc" | validate$names=="ind_rec_sj")]
            validate$data <- validate$data[which(validate$names=="temp2" | validate$names=="precip" |
                                                   validate$names=="slp"),]
            validate$ensmean <- validate$data # defined later
            validate$names <- validate$names[which(validate$names=="temp2" | validate$names=="precip" |
                                                     validate$names=="slp")]
            validate$lon <- validate$lon[which(validate$names=="temp2" | validate$names=="precip" |
                                                 validate$names=="slp")]
            validate$lat <- validate$lat[which(validate$names=="temp2" | validate$names=="precip" |
                                                 validate$names=="slp")]
          }
          # ACHTUNG does NOT work if variables are in different order as here with sixmonstatevector
          # therefore "sort" has been added
          #   tpspos <- sort(c(which(echam$names=='temp2'), which(echam$names=='precip'), 
          #               which(echam$names=='slp')))
          #   echam$data <- echam$data[tpspos,,]
          #   echam$ensmean <- echam$ensmean[tpspos,]
          #   echam$names <- echam$names[tpspos]
          #   echam$lon <- echam$lon[tpspos]
          #   echam$lat <- echam$lat[tpspos]
          #   tpspos1 <- sort(c(which(analysis$names=='temp2'), which(analysis$names=='precip'), 
          #                which(analysis$names=='slp')))
          #   analysis$data <- analysis$data[tpspos1,,]    
          #   analysis$ensmean <- analysis$ensmean[tpspos1,]
          #   analysis$names <- analysis$names[tpspos1]
          #   analysis$lon <- analysis$lon[tpspos1]
          #   analysis$lat <- analysis$lat[tpspos1]
          #   if (vali) {
          #     tpspos2 <- sort(c(which(validate$names=='temp2'), which(validate$names=='precip'), 
          #                  which(validate$names=='slp')))
          #     validate$data <- validate$data[tpspos2,]
          #     validate$names <- validate$names[tpspos2]
          #     validate$lon <- validate$lon[tpspos2]
          #     validate$lat <- validate$lat[tpspos2]
        }
        
        # convert data back to old format for plotting and analysis
        if (sixmonstatevector) {
          print('convert 6-months state vector')
          if (monthly_out) {
            print(paste('seasons =',s))
            s=12 # set back to 12 months for plotting
            echam$data <- array(echam$data,c((dim(echam$data)[1]/6),dim(echam$data)[2]*6,
                                             dim(echam$data)[3]))
            echam$ensmean <- array(echam$ensmean,c((dim(echam$ensmean)[1]/6),
                                                   dim(echam$ensmean)[2]*6))
            tmptime <- seq(syr,(eyr+1),by=1/12)
            echam$time <- tmptime[(season[2]+1):(length(tmptime)-4)]
            echam$names <- echam$names[1:dim(echam$data)[1]]
            echam$lon <- echam$lon[1:(dim(echam$data)[1])] #/length(unique(echam$names)))]
            echam$lat <- echam$lat[1:(dim(echam$data)[1])] #/length(unique(echam$names)))]
            echam.anom$data <- array(echam.anom$data,c((dim(echam.anom$data)[1]/6),dim(echam.anom$data)[2]*6,
                                                       dim(echam.anom$data)[3]))
            echam.anom$ensmean <- array(echam.anom$ensmean,c((dim(echam.anom$ensmean)[1]/6),
                                                             dim(echam.anom$ensmean)[2]*6))
            echam.anom$time <- tmptime[(season[2]+1):(length(tmptime)-4)]
            echam.anom$names <- echam.anom$names[1:dim(echam.anom$data)[1]]
            echam.anom$lon <- echam.anom$lon[1:(dim(echam.anom$data)[1])] #/length(unique(echam.anom$names)))]
            echam.anom$lat <- echam.anom$lat[1:(dim(echam.anom$data)[1])] #/length(unique(echam.anom$names)))]
            if (ind_ECHAM) {
              ech_ind$data <- array(ech_ind$data,c((dim(ech_ind$data)[1]/6),dim(ech_ind$data)[2]*6,
                                                   dim(ech_ind$data)[3]))
              ech_ind$ensmean <- array(ech_ind$ensmean,c((dim(ech_ind$ensmean)[1]/6),
                                                         dim(ech_ind$ensmean)[2]*6))    
              ech_ind$time <- tmptime[(season[2]+1):(length(tmptime)-4)]
              ech_ind$names <- ech_ind$names[1:dim(ech_ind$data)[1]]
              ech_ind$lon <- ech_ind$lon[1:(dim(ech_ind$data)[1])] #/length(unique(ech_ind$names)))]
              ech_ind$lat <- ech_ind$lat[1:(dim(ech_ind$data)[1])]
            }
            analysis$data <- array(analysis$data,c((dim(analysis$data)[1]/6),
                                                   dim(analysis$data)[2]*6,dim(analysis$data)[3]))
            analysis$ensmean <- array(analysis$ensmean,c((dim(analysis$ensmean)[1]/6),
                                                         dim(analysis$ensmean)[2]*6))
            analysis$time <- tmptime[(season[2]+1):(length(tmptime)-4)]
            analysis$names <- analysis$names[1:dim(analysis$data)[1]]
            analysis$lon <- analysis$lon[1:(dim(analysis$data)[1])] #/length(unique(analysis$names)))]
            analysis$lat <- analysis$lat[1:(dim(analysis$data)[1])] 
            analysis.anom$data <- array(analysis.anom$data,c((dim(analysis.anom$data)[1]/6),
                                                             dim(analysis.anom$data)[2]*6,dim(analysis.anom$data)[3]))
            analysis.anom$ensmean <- array(analysis.anom$ensmean,c((dim(analysis.anom$ensmean)[1]/6),
                                                                   dim(analysis.anom$ensmean)[2]*6))
            analysis.anom$time <- tmptime[(season[2]+1):(length(tmptime)-4)]
            analysis.anom$names <- analysis.anom$names[1:dim(analysis.anom$data)[1]]
            analysis.anom$lon <- analysis.anom$lon[1:(dim(analysis.anom$data)[1])] #/length(unique(analysis.anom$names)))]
            analysis.anom$lat <- analysis.anom$lat[1:(dim(analysis.anom$data)[1])] 
            if (ind_ECHAM) {
              ana_ind$data <- array(ana_ind$data,c((dim(ana_ind$data)[1]/6),dim(ana_ind$data)[2]*6,
                                                   dim(ana_ind$data)[3]))
              ana_ind$ensmean <- array(ana_ind$ensmean,c((dim(ana_ind$ensmean)[1]/6),
                                                         dim(ana_ind$ensmean)[2]*6))    
              ana_ind$time <- tmptime[(season[2]+1):(length(tmptime)-4)]
              ana_ind$names <- ana_ind$names[1:dim(ana_ind$data)[1]]
              ana_ind$lon <- ana_ind$lon[1:(dim(ana_ind$data)[1])] #/length(unique(ana_ind$names)))]
              ana_ind$lat <- ana_ind$lat[1:(dim(ana_ind$data)[1])]
            }
            if (vali) {
              if (!recon_vali) {
                #          if (cyr==1901) {
                #            validate$data <- cbind(rep(NA,length(validate$data)),validate$data)
                #            validate$ensmean <- cbind(rep(NA,length(validate$ensmean)),validate$ensmean)
                #          }
                validate$data <- array(validate$data,c((dim(validate$data)[1]/6),
                                                       dim(validate$data)[2]*6))
                validate$ensmean <- array(validate$ensmean,c((dim(validate$ensmean)[1]/6),
                                                             dim(validate$ensmean)[2]*6))
                validate$time <- tmptime[(season[2]+1):(length(tmptime)-4)]
                validate$names <- validate$names[1:dim(validate$data)[1]]
                validate$lon <- validate$lon[1:(dim(validate$data)[1])] #/length(unique(validate$names)))]
                validate$lat <- validate$lat[1:(dim(validate$data)[1])]
                if (ind_recon) {
                  #          vali_ind$data <- array(vali_ind$data,c((dim(vali_ind$data)[1]/6),dim(vali_ind$data)[2]*6,
                  #                                                 dim(vali_ind$data)[3]))
                  vali_ind$data <- array(vali_ind$data,c((dim(vali_ind$data)[1]/6),
                                                         dim(vali_ind$data)[2]*6))    
                  vali_ind$time <- tmptime[(season[2]+1):(length(tmptime)-4)]
                  vali_ind$names <- vali_ind$names[1:dim(vali_ind$data)[1]]
                  vali_ind$lon <- vali_ind$lon[1:(dim(vali_ind$data)[1])] #/length(unique(vali_ind$names)))]
                  vali_ind$lat <- vali_ind$lat[1:(dim(vali_ind$data)[1])]
                }
              }
            }
            #       if (!real_proxies) {
            #         calibrate$data <- array(calibrate$data,c((dim(calibrate$data)[1]/6),
            #                                                  dim(calibrate$data)[2]*6))
            #         calibrate$time <- tmptime[(season[2]+1):(length(tmptime)-4)]
            #         calibrate$names <- calibrate$names[1:dim(calibrate$data)[1]]
            #         calibrate$lon <- calibrate$lon[1:dim(calibrate$data)[1]]
            #         calibrate$lat <- calibrate$lat[1:dim(calibrate$data)[1]]
            #       }
          } else { # seasonal output, probably summer/winter averages
            s=length(season)
            print(paste('seasons =',s))
            etmp <- array(echam$data,c((dim(echam$data)[1]/6),6,dim(echam$data)[2],
                                       dim(echam$data)[3]))
            echam$data <- array(NA,c(dim(etmp)[1],dim(etmp)[3:4]))
            for (ensmem in 1:(dim(etmp)[4])) {
              print(paste('ECHAM member',ensmem))
              echam$data[,,ensmem] <- apply(etmp[,,,ensmem],c(1,3),mean)
            }
            #      echam$data <- apply(etmp,c(1,3,4),mean)
            etmp2 <- array(echam$ensmean,c((dim(echam$ensmean)[1]/6),6,dim(echam$ensmean)[2]))
            echam$ensmean <- apply(etmp2,c(1,3),mean)
            #    echam$names <- echam$names[1:dim(echam$data)[1]/length(unique(echam$names))]
            echam$names <- rep(unique(echam$names),each=dim(echam$data)[1]/
                                 length(unique(echam$names)))  
            echam$lon <- echam$lon[1:(dim(echam$data)[1])] #/length(unique(echam$names)))]
            echam$lat <- echam$lat[1:(dim(echam$data)[1])] #/length(unique(echam$names)))]
            if (ind_ECHAM) {
              etmp <- array(ech_ind$data,c((dim(ech_ind$data)[1]/6),6,dim(ech_ind$data)[2],
                                           dim(ech_ind$data)[3]))  
              ech_ind$data <- apply(etmp,c(1,3,4),mean)
              etmp2 <- array(ech_ind$ensmean,c((dim(ech_ind$ensmean)[1]/6),6,dim(ech_ind$ensmean)[2]))
              ech_ind$ensmean <- apply(etmp2,c(1,3),mean)
              #    ech_ind$names <- ech_ind$names[1:dim(ech_ind$data)[1]/length(unique(ech_ind$names))]
              ech_ind$names <- rep(unique(ech_ind$names),each=dim(ech_ind$data)[1]/
                                     length(unique(ech_ind$names)))  
              ech_ind$lon <- ech_ind$lon[1:(dim(ech_ind$data)[1])] #/length(unique(ech_ind$names)))]
              ech_ind$lat <- ech_ind$lat[1:(dim(ech_ind$data)[1])] #/length(unique(ech_ind$names)))] 
            }
            etmp <- NULL
            etmp2 <- NULL
            
            atmp <- array(analysis$data,c((dim(analysis$data)[1]/6),6,dim(analysis$data)[2],
                                          dim(analysis$data)[3]))  
            analysis$data <- array(NA,c(dim(atmp)[1],dim(atmp)[3:4]))
            for (ensmem in 1:(dim(atmp)[4])) {
              print(paste('Analysis member',ensmem))
              analysis$data[,,ensmem] <- apply(atmp[,,,ensmem],c(1,3),mean)
            }
            #      analysis$data <- apply(atmp,c(1,3,4),mean)
            #    analysis$names <- analysis$names[1:dim(analysis$data)[1]]
            analysis$names <- rep(unique(analysis$names),each=dim(analysis$data)[1]/
                                    length(unique(analysis$names)))
            analysis$lon <- analysis$lon[1:(dim(analysis$data)[1])] #/length(unique(analysis$names)))]
            analysis$lat <- analysis$lat[1:(dim(analysis$data)[1])] #/length(unique(analysis$names)))]
            atmp2 <- array(analysis$ensmean,c((dim(analysis$ensmean)[1]/6),6,
                                              dim(analysis$ensmean)[2]))
            analysis$ensmean <- apply(atmp2,c(1,3),mean)
            if (ind_ECHAM) {
              atmp <- array(ana_ind$data,c((dim(ana_ind$data)[1]/6),6,dim(ana_ind$data)[2],
                                           dim(ana_ind$data)[3]))  
              ana_ind$data <- apply(atmp,c(1,3,4),mean)
              #    ana_ind$names <- ana_ind$names[1:dim(ana_ind$data)[1]]
              ana_ind$names <- rep(unique(ana_ind$names),each=dim(ana_ind$data)[1]/
                                     length(unique(ana_ind$names)))
              ana_ind$lon <- ana_ind$lon[1:(dim(ana_ind$data)[1])] #/length(unique(ana_ind$names)))]
              ana_ind$lat <- ana_ind$lat[1:(dim(ana_ind$data)[1])] #/length(unique(ana_ind$names)))]
              atmp2 <- array(ana_ind$ensmean,c((dim(ana_ind$ensmean)[1]/6),6,
                                               dim(ana_ind$ensmean)[2]))
              ana_ind$ensmean <- apply(atmp2,c(1,3),mean)
            }
            atmp <- NULL
            atmp2 <- NULL
            
            # same for anomalies
            etmp <- array(echam.anom$data,c((dim(echam.anom$data)[1]/6),6,dim(echam.anom$data)[2],
                                            dim(echam.anom$data)[3]))
            echam.anom$data <- array(NA,c(dim(etmp)[1],dim(etmp)[3:4]))
            for (ensmem in 1:(dim(etmp)[4])) {
              print(paste('ECHAM anomaly member',ensmem))
              echam.anom$data[,,ensmem] <- apply(etmp[,,,ensmem],c(1,3),mean)
            }
            #    echam.anom$data <- apply(etmp,c(1,3,4),mean)
            etmp2 <- array(echam.anom$ensmean,c((dim(echam.anom$ensmean)[1]/6),6,dim(echam.anom$ensmean)[2]))
            echam.anom$ensmean <- apply(etmp2,c(1,3),mean)
            #    echam.anom$names <- echam.anom$names[1:dim(echam.anom$data)[1]/length(unique(echam.anom$names))]
            echam.anom$names <- rep(unique(echam.anom$names),each=dim(echam.anom$data)[1]/
                                      length(unique(echam.anom$names)))  
            echam.anom$lon <- echam.anom$lon[1:(dim(echam.anom$data)[1])] #/length(unique(echam.anom$names)))]
            echam.anom$lat <- echam.anom$lat[1:(dim(echam.anom$data)[1])] #/length(unique(echam.anom$names)))]
            if (ind_ECHAM) {
              etmp <- array(ech_ind$data,c((dim(ech_ind$data)[1]/6),6,dim(ech_ind$data)[2],
                                           dim(ech_ind$data)[3])) 
            }
            etmp <- NULL
            
            atmp <- array(analysis.anom$data,c((dim(analysis.anom$data)[1]/6),6,dim(analysis.anom$data)[2],
                                               dim(analysis.anom$data)[3])) 
            analysis.anom$data <- array(NA,c(dim(atmp)[1],dim(atmp)[3:4]))
            for (ensmem in 1:(dim(atmp)[4])) {
              print(paste('Analysis anomaly member',ensmem))
              analysis.anom$data[,,ensmem] <- apply(atmp[,,,ensmem],c(1,3),mean)
            }
            #      analysis.anom$data <- apply(atmp,c(1,3,4),mean)
            #    analysis.anom$names <- analysis.anom$names[1:dim(analysis.anom$data)[1]]
            analysis.anom$names <- rep(unique(analysis.anom$names),each=dim(analysis.anom$data)[1]/
                                         length(unique(analysis.anom$names)))
            analysis.anom$lon <- analysis.anom$lon[1:(dim(analysis.anom$data)[1])] #/length(unique(analysis.anom$names)))]
            analysis.anom$lat <- analysis.anom$lat[1:(dim(analysis.anom$data)[1])] #/length(unique(analysis.anom$names)))]
            atmp2 <- array(analysis.anom$ensmean,c((dim(analysis.anom$ensmean)[1]/6),6,
                                                   dim(analysis.anom$ensmean)[2]))
            analysis.anom$ensmean <- apply(atmp2,c(1,3),mean)
            atmp <- NULL
            
            if (vali) {    
              if (!recon_vali) {
                #           if (cyr==1901) {
                #             validate$data <- cbind(rep(NA,length(validate$data)),validate$data)
                #             validate$ensmean <- cbind(rep(NA,length(validate$ensmean)),validate$ensmean)
                #           }
                vtmp <- array(validate$data,c((dim(validate$data)[1]/6),6,dim(validate$data)[2]))  
                validate$data <- apply(vtmp,c(1,3),mean)
                validate$ensmean <- validate$data
                #      validate$names <- validate$names[1:dim(validate$data)[1]]
                validate$names <- rep(unique(validate$names),each=dim(validate$data)[1]/
                                        length(unique(validate$names)))
                validate$lon <- rep(validate$lon,length(unique(validate$names)))
                validate$lat <- rep(validate$lat,length(unique(validate$names)))
                #           validate$lon <- validate$lon[1:(dim(validate$data)[1])] #/length(unique(validate$names)))]
                #           validate$lat <- validate$lat[1:(dim(validate$data)[1])] #/length(unique(validate$names)))]
                # #          vtmp <- array(vali_ind$data,c((dim(vali_ind$data)[1]/6),6,dim(vali_ind$data)[2]))  
                # #          vali_ind$data <- apply(vtmp,c(1,3),mean)
                # #          vali_ind$ensmean <- vali_ind$data
                # #          #      vali_ind$names <- vali_ind$names[1:dim(vali_ind$data)[1]]
                # #          vali_ind$names <- rep(unique(vali_ind$names),each=dim(vali_ind$data)[1]/
                # #                                  length(unique(vali_ind$names)))
                # #          vali_ind$lon <- vali_ind$lon[1:(dim(vali_ind$data)[1])] #/length(unique(vali_ind$names)))]
                # #          vali_ind$lat <- vali_ind$lat[1:(dim(vali_ind$data)[1])] #/length(unique(vali_ind$names)))]
              }
            }
            # #       if (!real_proxies) {
            # #         # ERROR, NOT ADJUSTED FOR seasonal/annual DOCUMENTARY DATA, YET !!!
            # #         ctmp <- array(calibrate$data,c((dim(calibrate$data)[1]/6),6,
            # #                                        dim(calibrate$data)[2]))  
            # #         calibrate$data <- apply(ctmp,c(1,3),mean)
            # #         calibrate$names <- calibrate$names[1:dim(calibrate$data)[1]] 
            # #         calibrate$lon <- calibrate$lon[1:dim(calibrate$data)[1]]
            # #         calibrate$lat <- calibrate$lat[1:dim(calibrate$data)[1]]
            # #       }
            #    }
          }
          print("transformed 6-mon state vector")
          
          
          #  if (save_prepplot){
          if (vali) {
            print(paste("dim validate data:",paste(nrow(validate$data),ncol(validate$data))))
            save(analysis,analysis.anom,echam,echam.anom,validate,calibrate,
                 file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
            #save(analysis,analysis.anom,ana_ind,echam,echam.anom,ech_ind,validate,vali_ind,calibrate,
            #     file=paste0(prepplotdir,'/analysis_',cyr,'.Rdata'))
            #paste0('../data/prepplot/analysis_',cyr,'.Rdata'))
          } else {
            save(analysis,analysis.anom,echam,echam.anom,calibrate,
                 file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
            #save(analysis,analysis.anom,ana_ind,echam,echam.anom,ech_ind,calibrate,
            #     file=paste0(prepplotdir,'/analysis_',cyr,'.Rdata'))
            #paste0('../data/prepplot/analysis_',cyr,'.Rdata'))
          }
          #  } else {
          if (cyr == syr) {
            analysis.allts=analysis
            analysis.anom.allts=analysis.anom
            if (ind_ECHAM) {
              ana_ind.allts=ana_ind
            }
            if (vali) {
              validate.allts=validate
              if (ind_recon) {
                vali_ind.allts=vali_ind
              }
            }
            echam.allts=echam
            echam.anom.allts=echam.anom
            if (ind_ECHAM) {
              ech_ind.allts=ech_ind
            }
          } else {
            analysis.allts$data=abind(analysis.allts$data,analysis$data,along=2)
            analysis.allts$ensmean=cbind(analysis.allts$ensmean,analysis$ensmean)
            analysis.allts$time=c(analysis.allts$time,analysis$time)
            analysis.anom.allts$data=abind(analysis.anom.allts$data,analysis.anom$data,along=2)
            analysis.anom.allts$ensmean=cbind(analysis.anom.allts$ensmean,analysis.anom$ensmean)
            analysis.anom.allts$time=c(analysis.anom.allts$time,analysis.anom$time)
            if (ind_ECHAM) {
              ana_ind.allts$data=abind(ana_ind.allts$data,ana_ind$data,along=2)
              ana_ind.allts$ensmean=cbind(ana_ind.allts$ensmean,ana_ind$ensmean)
              ana_ind.allts$time=c(ana_ind.allts$time,ana_ind$time)
            }
            if (vali) {
              validate.allts$data=cbind(validate.allts$data,validate$data)
              validate.allts$ensmean=cbind(validate.allts$ensmean,validate$ensmean)
              validate.allts$time=c(validate.allts$time,validate$time)
              if (ind_recon) {
                vali_ind.allts$data=cbind(vali_ind.allts$data,vali_ind$data)
                vali_ind.allts$ensmean=cbind(vali_ind.allts$ensmean,vali_ind$ensmean)
                vali_ind.allts$time=c(vali_ind.allts$time,vali_ind$time)
              }
            }
            echam.allts$data=abind(echam.allts$data,echam$data,along=2)
            echam.allts$ensmean=cbind(echam.allts$ensmean,echam$ensmean)
            echam.allts$time=c(echam.allts$time,echam$time)
            echam.anom.allts$data=abind(echam.anom.allts$data,echam.anom$data,along=2)
            echam.anom.allts$ensmean=cbind(echam.anom.allts$ensmean,echam.anom$ensmean)
            echam.anom.allts$time=c(echam.anom.allts$time,echam.anom$time)
            if (ind_ECHAM) {
              ech_ind.allts$data=abind(ech_ind.allts$data,ech_ind$data,along=2)
              ech_ind.allts$ensmean=cbind(ech_ind.allts$ensmean,ech_ind$ensmean)
              ech_ind.allts$time=c(ech_ind.allts$time,ech_ind$time)
            }
          }  
        }
      } # loop over years
    } # end calc_prepplot=F loop
     
########################################################################################################
     
     if (load_prepplot){
       # load stat data network for specific year
       stat_yr=statyr #1850
       load(paste0(prepplotdir,'/analysis_',stat_yr,'.Rdata'))
       cali <- calibrate
       for (cyr in syr:(eyr-1)) {
         print(paste('year',cyr))
         if ((cyr < 1902)) { #1751)) { #} |  (cyr > 1979)) {
           vali=F                 # switch off prepplot if no vali data selected
         } else {
           vali=T
         }
         if ((cyr > 1901) & (cyr < 2005)) { #  & (syr > 1900) & (eyr < 2005)) {
           cru_vali=T             # monthly CRU TS3 temp, precip and HADSLP2 gridded instrumentals (1901-2004)
           #      if (ind_recon){
           #        ind_recon=T            # Stefan's reconstructed indices until 1948 and NCAR reanalysis later added to CRU and NCEP
           #      }
         } else {
           cru_vali=F 
           #      if (ind_recon){
           #        ind_recon=F
           #      }
         }
         #if ((cyr < 1901) & (cyr > 1749)) { # & (syr < 1901) & (eyr > 1749)) {
         #  recon_vali=T           # seasonal luterbacher, pauling, kuettel recons (1750-1999)
         #} else {
         recon_vali=F
         #}
         # t=syr2
         print(paste("recon_vali=",recon_vali))
         print(paste("vali=",vali))
         
         # load data and make normal calendar year Jan-Dec
         #    if (cyr==eyr) {stop("last year cannot be created because oct-dec data is missing")}
         load(file=paste0(prepplotdir,'analysis_',(cyr+1),'.Rdata')) 
         echam2 <- echam
         analysis2 <- analysis
         load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata')) 
         if (monthly_out) {
           echam$data <- abind(echam$data[,4:12,],echam2$data[,1:3,],along=2)
           echam$ensmean <- abind(echam$ensmean[,4:12],echam2$ensmean[,1:3],along=2)
           analysis$data <- abind(analysis$data[,4:12,],analysis2$data[,1:3,],along=2)
           analysis$ensmean <- abind(analysis$ensmean[,4:12],analysis2$ensmean[,1:3],along=2)
         }
         
         #    lenvar2 <- length(c(which(validate$names=="temp2"), which(validate$names=="precip"), 
         #                        which(validate$names=="slp")))
         lenvar2 <- length(c(which(echam$names=="temp2"), which(echam$names=="precip"), 
                             which(echam$names=="slp")))
         #    analysis_noindex=analysis
         analysis$data=analysis$data[1:lenvar2,,]
         analysis$ensmean=analysis$ensmean[1:lenvar2,]
         analysis$names=analysis$names[1:lenvar2]
         analysis$lon=analysis$lon[1:lenvar2]
         analysis$lat=analysis$lat[1:lenvar2]
         echam=echam
         echam$data=echam$data[1:lenvar2,,]
         echam$ensmean=echam$ensmean[1:lenvar2,]
         echam$names=echam$names[1:lenvar2]
         echam$lon=echam$lon[1:lenvar2]
         echam$lat=echam$lat[1:lenvar2]
         analysis.anom=analysis.anom
         analysis.anom$data=analysis.anom$data[1:lenvar2,,]
         analysis.anom$ensmean=analysis.anom$ensmean[1:lenvar2,]
         analysis.anom$names=analysis.anom$names[1:lenvar2]
         analysis.anom$lon=analysis.anom$lon[1:lenvar2]
         analysis.anom$lat=analysis.anom$lat[1:lenvar2]
         echam.anom=echam.anom
         echam.anom$data=echam.anom$data[1:lenvar2,,]
         echam.anom$ensmean=echam.anom$ensmean[1:lenvar2,]
         echam.anom$names=echam.anom$names[1:lenvar2]
         echam.anom$lon=echam.anom$lon[1:lenvar2]
         echam.anom$lat=echam.anom$lat[1:lenvar2]
         
         if (cyr == syr) {
           analysis.allts=analysis
           analysis.anom.allts=analysis.anom
           #      ana_ind.allts=ana_ind
           if (vali) {
             validate.allts=validate
             #        vali_ind.allts=vali_ind
           }
           echam.allts=echam
           echam.anom.allts=echam.anom
           #      ech_ind.allts=ech_ind
           # only instr. calibration data for period with fixed network
           # because proxies not in fixed grid format. Hence number of records 
           # and position/row in data/ensmean matrix changes over time
           calibrate.allts=calibrate
           pos <- which(calibrate$sour=="inst")
           if (length(pos)==0){
             pos <- which(cali$sour=="inst")
             calibrate.allts$data=array(NA,dim=dim(cali$data[pos,]))
             calibrate.allts$lon=as.matrix(cali$lon[pos])
             calibrate.allts$lat=as.matrix(cali$lat[pos])
             calibrate.allts$names=as.matrix(cali$names[pos])
             calibrate.allts$sour=as.matrix(cali$sour[pos])
           } else {
             calibrate.allts$data=calibrate$data[pos,]
             calibrate.allts$lon=as.matrix(calibrate$lon[pos])
             calibrate.allts$lat=as.matrix(calibrate$lat[pos])
             calibrate.allts$names=as.matrix(calibrate$names[pos])
             calibrate.allts$sour=as.matrix(calibrate$sour[pos])
             #      calibrate.anom.allts=calibrate.anom
             #      calibrate.anom.allts$data=calibrate.anom$data[pos,]
             #      calibrate.anom.allts$lon=as.matrix(calibrate.anom$lon[pos])
             #      calibrate.anom.allts$lat=as.matrix(calibrate.anom$lat[pos])
             #      calibrate.anom.allts$names=as.matrix(calibrate.anom$names[pos])
             #      calibrate.anom.allts$sour=as.matrix(calibrate.anom$sour[pos])
           }
         } else {
           if (cyr==1751 & vali) {
             validate.allts=validate 
           }
           analysis.allts$data=abind(analysis.allts$data,analysis$data,along=2)
           analysis.allts$ensmean=cbind(analysis.allts$ensmean,analysis$ensmean)
           analysis.allts$time=c(analysis.allts$time,analysis$time)
           analysis.anom.allts$data=abind(analysis.anom.allts$data,analysis.anom$data,along=2)
           analysis.anom.allts$ensmean=cbind(analysis.anom.allts$ensmean,analysis.anom$ensmean)
           analysis.anom.allts$time=c(analysis.anom.allts$time,analysis.anom$time)
           #      ana_ind.allts$data=abind(ana_ind.allts$data,ana_ind$data,along=2)
           #      ana_ind.allts$ensmean=cbind(ana_ind.allts$ensmean,ana_ind$ensmean)
           #      ana_ind.allts$time=c(ana_ind.allts$time,ana_ind$time)
           if (vali & cyr>1751) {
             validate.allts$data=cbind(validate.allts$data,validate$data)
             validate.allts$ensmean=cbind(validate.allts$ensmean,validate$ensmean)
             validate.allts$time=c(validate.allts$time,validate$time)
             #        vali_ind.allts$data=cbind(vali_ind.allts$data,vali_ind$data)
             #        vali_ind.allts$ensmean=cbind(vali_ind.allts$ensmean,vali_ind$ensmean)
             #        vali_ind.allts$time=c(vali_ind.allts$time,vali_ind$time)
           }
           echam.allts$data=abind(echam.allts$data,echam$data,along=2)
           echam.allts$ensmean=cbind(echam.allts$ensmean,echam$ensmean)
           echam.allts$time=c(echam.allts$time,echam$time)
           echam.anom.allts$data=abind(echam.anom.allts$data,echam.anom$data,along=2)
           echam.anom.allts$ensmean=cbind(echam.anom.allts$ensmean,echam.anom$ensmean)
           echam.anom.allts$time=c(echam.anom.allts$time,echam.anom$time)
           #      ech_ind.allts$data=abind(ech_ind.allts$data,ech_ind$data,along=2)
           #      ech_ind.allts$ensmean=cbind(ech_ind.allts$ensmean,ech_ind$ensmean)
           #      ech_ind.allts$time=c(ech_ind.allts$time,ech_ind$time)
           pos <- which(calibrate$sour=="inst")
           if (length(pos)==0){
             pos <- which(cali$sour=="inst")
             calibrate.allts$data=abind(calibrate.allts$data[pos,],
                                        array(NA,dim=dim(cali$data))[pos,],along=2)
             calibrate.allts$time=c(calibrate.allts$time,calibrate$time)
             calibrate.allts$lon=cbind(calibrate.allts$lon[pos,],cali$lon[pos])
             calibrate.allts$lat=cbind(calibrate.allts$lat[pos,],cali$lat[pos])
             calibrate.allts$names=cbind(calibrate.allts$names[pos,],cali$names[pos])
             calibrate.allts$sour=cbind(calibrate.allts$sour[pos,],cali$sour[pos])
           } else {
             calibrate.allts$data=abind(calibrate.allts$data[pos,],calibrate$data[pos,],along=2)
             calibrate.allts$time=c(calibrate.allts$time,calibrate$time)
             calibrate.allts$lon=cbind(calibrate.allts$lon[pos,],calibrate$lon[pos])
             calibrate.allts$lat=cbind(calibrate.allts$lat[pos,],calibrate$lat[pos])
             calibrate.allts$names=cbind(calibrate.allts$names[pos,],calibrate$names[pos])
             calibrate.allts$sour=cbind(calibrate.allts$sour[pos,],calibrate$sour[pos])     
             #      calibrate.anom.allts$data=abind(calibrate.anom.allts$data[pos,],calibrate.anom$data[pos,],along=2)
             #      calibrate.anom.allts$time=c(calibrate.anom.allts$time,calibrate.anom$time)
             #      calibrate.anom.allts$lon=cbind(calibrate.anom.allts$lon[pos,],calibrate.anom$lon[pos])
             #      calibrate.anom.allts$lat=cbind(calibrate.anom.allts$lat[pos,],calibrate.anom$lat[pos])
             #      calibrate.anom.allts$names=cbind(calibrate.anom.allts$names[pos,],calibrate.anom$names[pos])
             #      calibrate.anom.allts$sour=cbind(calibrate.anom.allts$sour[pos,],calibrate.anom$sour[pos])
           }
         }  
       }
       
       # if (cyr == syr) {
       #   analysis.allts=analysis
       #   analysis.anom.allts=analysis.anom
       #   ana_ind.allts=ana_ind
       #   if (vali) {
       #     validate.allts=validate
       #     vali_ind.allts=vali_ind
       #   }
       #   echam.allts=echam
       #   echam.anom.allts=echam.anom
       #   ech_ind.allts=ech_ind
       # } else {
       #   analysis.allts$data=abind(analysis.allts$data,analysis$data,along=2)
       #   analysis.allts$ensmean=cbind(analysis.allts$ensmean,analysis$ensmean)
       #   analysis.allts$time=c(analysis.allts$time,analysis$time)
       #   analysis.anom.allts$data=abind(analysis.anom.allts$data,analysis.anom$data,along=2)
       #   analysis.anom.allts$ensmean=cbind(analysis.anom.allts$ensmean,analysis.anom$ensmean)
       #   analysis.anom.allts$time=c(analysis.anom.allts$time,analysis.anom$time)
       #   ana_ind.allts$data=abind(ana_ind.allts$data,ana_ind$data,along=2)
       #   ana_ind.allts$ensmean=cbind(ana_ind.allts$ensmean,ana_ind$ensmean)
       #   ana_ind.allts$time=c(ana_ind.allts$time,ana_ind$time)
       #   if (vali) {
       #     validate.allts$data=cbind(validate.allts$data,validate$data)
       #     validate.allts$ensmean=cbind(validate.allts$ensmean,validate$ensmean)
       #     validate.allts$time=c(validate.allts$time,validate$time)
       #     vali_ind.allts$data=cbind(vali_ind.allts$data,vali_ind$data)
       #     vali_ind.allts$ensmean=cbind(vali_ind.allts$ensmean,vali_ind$ensmean)
       #     vali_ind.allts$time=c(vali_ind.allts$time,vali_ind$time)
       #   }
       #   echam.allts$data=abind(echam.allts$data,echam$data,along=2)
       #   echam.allts$ensmean=cbind(echam.allts$ensmean,echam$ensmean)
       #   echam.allts$time=c(echam.allts$time,echam$time)
       #   echam.anom.allts$data=abind(echam.anom.allts$data,echam.anom$data,along=2)
       #   echam.anom.allts$ensmean=cbind(echam.anom.allts$ensmean,echam.anom$ensmean)
       #   echam.anom.allts$time=c(echam.anom.allts$time,echam.anom$time)
       #   ech_ind.allts$data=abind(ech_ind.allts$data,ech_ind$data,along=2)
       #   ech_ind.allts$ensmean=cbind(ech_ind.allts$ensmean,ech_ind$ensmean)
       #   ech_ind.allts$time=c(ech_ind.allts$time,ech_ind$time)
       # }
       #}
       
       
       
       
       echam <- echam.allts
       echam.anom <- echam.anom.allts
       #ech_ind <- ech_ind.allts
       rm(echam.allts,echam.anom.allts,ech_ind.allts)
       analysis <- analysis.allts
       analysis.anom <- analysis.anom.allts
       #ana_ind <- ana_ind.allts
       rm(analysis.allts,analysis.anom.allts,ana_ind.allts)
       calibrate <- calibrate.allts
       #calibrate.anom <- calibrate.anom.allts
       if (vali) {
         validate <- validate.allts
         #  vali_ind <- vali_ind.allts
         rm(validate.allts,vali_ind.allts)
       }
       print("calc time for a year")

       
       # validate$data[validate$names=="temp2",]<- validate$data[validate$names=="temp2",]-273.15
       # validate$ensmean[validate$names=="temp2",]<- validate$ensmean[validate$names=="temp2",]-273.15
       # validate$data[validate$names=="slp",]<- validate$data[validate$names=="slp",]/100
       # validate$ensmean[validate$names=="slp",]<- validate$ensmean[validate$names=="slp",]/100
       
       # # convert data back to old format for plotting and analysis
       # if (sixmonstatevector) {
       #   if (monthly_out) {
       #     s=12 # set back to 12 months for plotting
       #     echam$data <- array(echam$data,c((dim(echam$data)[1]/6),dim(echam$data)[2]*6,
       #                                      dim(echam$data)[3]))
       #     echam$ensmean <- array(echam$ensmean,c((dim(echam$ensmean)[1]/6),
       #                                            dim(echam$ensmean)[2]*6))
       #     tmptime <- seq(syr,(eyr+1),by=1/12)
       #     echam$time <- tmptime[(season[2]+1):(length(tmptime)-4)]
       #     echam$names <- echam$names[1:dim(echam$data)[1]]
       #     ech_ind$data <- array(ech_ind$data,c((dim(ech_ind$data)[1]/6),dim(ech_ind$data)[2]*6,
       #                           dim(ech_ind$data)[3]))
       #     ech_ind$ensmean <- array(ech_ind$ensmean,c((dim(ech_ind$ensmean)[1]/6),
       #                              dim(ech_ind$ensmean)[2]*6))    
       #     ech_ind$time <- tmptime[(season[2]+1):(length(tmptime)-4)]
       #     ech_ind$names <- ech_ind$names[1:dim(ech_ind$data)[1]]
       #     analysis$data <- array(analysis$data,c((dim(analysis$data)[1]/6),
       #                                            dim(analysis$data)[2]*6,dim(analysis$data)[3]))
       #     analysis$ensmean <- array(analysis$ensmean,c((dim(analysis$ensmean)[1]/6),
       #                                                  dim(analysis$ensmean)[2]*6))
       #     analysis$time <- tmptime[(season[2]+1):(length(tmptime)-4)]
       #     analysis$names <- analysis$names[1:dim(analysis$data)[1]]
       #     ana_ind$data <- array(ana_ind$data,c((dim(ana_ind$data)[1]/6),dim(ana_ind$data)[2]*6,
       #                                          dim(ana_ind$data)[3]))
       #     ana_ind$ensmean <- array(ana_ind$ensmean,c((dim(ana_ind$ensmean)[1]/6),
       #                                                dim(ana_ind$ensmean)[2]*6))    
       #     ana_ind$time <- tmptime[(season[2]+1):(length(tmptime)-4)]
       #     ana_ind$names <- ana_ind$names[1:dim(ana_ind$data)[1]]
       #     if (vali) {
       #      if (!recon_vali) {
       #       validate$data <- array(validate$data,c((dim(validate$data)[1]/6),
       #                                              dim(validate$data)[2]*6))
       #       validate$ensmean <- array(validate$ensmean,c((dim(validate$ensmean)[1]/6),
       #                                                    dim(validate$ensmean)[2]*6))
       #       validate$time <- tmptime[(season[2]+1):(length(tmptime)-4)]
       #       validate$names <- validate$names[1:dim(validate$data)[1]]
       #       vali_ind$data <- array(vali_ind$data,c((dim(vali_ind$data)[1]/6),dim(vali_ind$data)[2]*6,
       #                                            dim(vali_ind$data)[3]))
       #       vali_ind$ensmean <- array(vali_ind$ensmean,c((dim(vali_ind$ensmean)[1]/6),
       #                                                  dim(vali_ind$ensmean)[2]*6))    
       #       vali_ind$time <- tmptime[(season[2]+1):(length(tmptime)-4)]
       #       vali_ind$names <- vali_ind$names[1:dim(vali_ind$data)[1]]
       #      }
       #     }
       #     if (!real_proxies) {
       #       calibrate$data <- array(calibrate$data,c((dim(calibrate$data)[1]/6),
       #                                                dim(calibrate$data)[2]*6))
       #       calibrate$time <- tmptime[(season[2]+1):(length(tmptime)-4)]
       #       calibrate$names <- calibrate$names[1:dim(calibrate$data)[1]]
       #       calibrate$lon <- calibrate$lon[1:dim(calibrate$data)[1]]
       #       calibrate$lat <- calibrate$lat[1:dim(calibrate$data)[1]]
       #     }
       #   } else { # seasonal output, probably summer/winter averages
       #     s=length(season) 
       #     etmp <- array(echam$data,c((dim(echam$data)[1]/6),6,dim(echam$data)[2],
       #                                dim(echam$data)[3]))  
       #     echam$data <- apply(etmp,c(1,3,4),mean)
       #     etmp2 <- array(echam$ensmean,c((dim(echam$ensmean)[1]/6),6,dim(echam$ensmean)[2]))
       #     echam$ensmean <- apply(etmp2,c(1,3),mean)
       # #    echam$names <- echam$names[1:dim(echam$data)[1]/length(unique(echam$names))]
       #     echam$names <- rep(unique(echam$names),each=dim(echam$data)[1]/
       #                          length(unique(echam$names)))  
       #     echam$lon <- echam$lon[1:(dim(echam$data)[1])] #/length(unique(echam$names)))]
       #     echam$lat <- echam$lat[1:(dim(echam$data)[1])] #/length(unique(echam$names)))]
       #     etmp <- array(ech_ind$data,c((dim(ech_ind$data)[1]/6),6,dim(ech_ind$data)[2],
       #                                dim(ech_ind$data)[3]))  
       #     ech_ind$data <- apply(etmp,c(1,3,4),mean)
       #     etmp2 <- array(ech_ind$ensmean,c((dim(ech_ind$ensmean)[1]/6),6,dim(ech_ind$ensmean)[2]))
       #     ech_ind$ensmean <- apply(etmp2,c(1,3),mean)
       #     #    ech_ind$names <- ech_ind$names[1:dim(ech_ind$data)[1]/length(unique(ech_ind$names))]
       #     ech_ind$names <- rep(unique(ech_ind$names),each=dim(ech_ind$data)[1]/
       #                          length(unique(ech_ind$names)))  
       #     ech_ind$lon <- ech_ind$lon[1:(dim(ech_ind$data)[1])] #/length(unique(ech_ind$names)))]
       #     ech_ind$lat <- ech_ind$lat[1:(dim(ech_ind$data)[1])] #/length(unique(ech_ind$names)))]    
       #     atmp <- array(analysis$data,c((dim(analysis$data)[1]/6),6,dim(analysis$data)[2],
       #                                   dim(analysis$data)[3]))  
       #     analysis$data <- apply(atmp,c(1,3,4),mean)
       # #    analysis$names <- analysis$names[1:dim(analysis$data)[1]]
       #     analysis$names <- rep(unique(analysis$names),each=dim(analysis$data)[1]/
       #                             length(unique(analysis$names)))
       #     analysis$lon <- analysis$lon[1:(dim(analysis$data)[1])] #/length(unique(analysis$names)))]
       #     analysis$lat <- analysis$lat[1:(dim(analysis$data)[1])] #/length(unique(analysis$names)))]
       #     atmp2 <- array(analysis$ensmean,c((dim(analysis$ensmean)[1]/6),6,
       #                                       dim(analysis$ensmean)[2]))
       #     analysis$ensmean <- apply(atmp2,c(1,3),mean)
       #     atmp <- array(ana_ind$data,c((dim(ana_ind$data)[1]/6),6,dim(ana_ind$data)[2],
       #                                   dim(ana_ind$data)[3]))  
       #     ana_ind$data <- apply(atmp,c(1,3,4),mean)
       #     #    ana_ind$names <- ana_ind$names[1:dim(ana_ind$data)[1]]
       #     ana_ind$names <- rep(unique(ana_ind$names),each=dim(ana_ind$data)[1]/
       #                             length(unique(ana_ind$names)))
       #     ana_ind$lon <- ana_ind$lon[1:(dim(ana_ind$data)[1])] #/length(unique(ana_ind$names)))]
       #     ana_ind$lat <- ana_ind$lat[1:(dim(ana_ind$data)[1])] #/length(unique(ana_ind$names)))]
       #     atmp2 <- array(ana_ind$ensmean,c((dim(ana_ind$ensmean)[1]/6),6,
       #                                       dim(ana_ind$ensmean)[2]))
       #     ana_ind$ensmean <- apply(atmp2,c(1,3),mean)
       #     if (vali) {    
       #      if (!recon_vali) {
       #       vtmp <- array(validate$data,c((dim(validate$data)[1]/6),6,dim(validate$data)[2]))  
       #       validate$data <- apply(vtmp,c(1,3),mean)
       #       validate$ensmean <- validate$data
       # #      validate$names <- validate$names[1:dim(validate$data)[1]]
       #       validate$names <- rep(unique(validate$names),each=dim(validate$data)[1]/
       #                               length(unique(validate$names)))
       #       validate$lon <- validate$lon[1:(dim(validate$data)[1])] #/length(unique(validate$names)))]
       #       validate$lat <- validate$lat[1:(dim(validate$data)[1])] #/length(unique(validate$names)))]
       #       vtmp <- array(vali_ind$data,c((dim(vali_ind$data)[1]/6),6,dim(vali_ind$data)[2]))  
       #       vali_ind$data <- apply(vtmp,c(1,3),mean)
       #       vali_ind$ensmean <- vali_ind$data
       #       #      vali_ind$names <- vali_ind$names[1:dim(vali_ind$data)[1]]
       #       vali_ind$names <- rep(unique(vali_ind$names),each=dim(vali_ind$data)[1]/
       #                               length(unique(vali_ind$names)))
       #       vali_ind$lon <- vali_ind$lon[1:(dim(vali_ind$data)[1])] #/length(unique(vali_ind$names)))]
       #       vali_ind$lat <- vali_ind$lat[1:(dim(vali_ind$data)[1])] #/length(unique(vali_ind$names)))]
       #      }
       #     }
       #     if (!real_proxies) {
       #       # ERROR, NOT ADJUSTED FOR seasonal/annual DOCUMENTARY DATA, YET !!!
       #       ctmp <- array(calibrate$data,c((dim(calibrate$data)[1]/6),6,
       #                                      dim(calibrate$data)[2]))  
       #       calibrate$data <- apply(ctmp,c(1,3),mean)
       #       calibrate$names <- calibrate$names[1:dim(calibrate$data)[1]] 
       #       calibrate$lon <- calibrate$lon[1:dim(calibrate$data)[1]]
       #       calibrate$lat <- calibrate$lat[1:dim(calibrate$data)[1]]
       #     }
       #   }
       # }
       # print("transformed 6-mon state vector")
       # print(proc.time() - ptm1)
       
       #lenvar <- length(c(which(validate$names=="temp2"), which(validate$names=="precip"), 
       #                   which(validate$names=="slp")))
       #validate$data=validate$data[1:lenvar,]
       #validate$ensmean=validate$data[1:lenvar,]
       
       # is building anom still necessary?
       # analysis.anom <- analysis
       # analysis.anom$data <- array(analysis$data - as.vector(apply(array(analysis$data, c(nrow(analysis$data), s, ncol(analysis$data)/s,30)), 1:2, mean)), c(nrow(analysis$data), s, ncol(analysis$data)/s, 30))
       # analysis.anom$data <- (array(analysis.anom$data,c(dim(analysis.anom$data)[1],dim(analysis.anom$data)[2]*dim(analysis.anom$data)[3],dim(analysis.anom$data)[4])))
       # analysis.anom$ensmean <- array(analysis$ensmean - as.vector(apply(array(analysis$ensmean, c(nrow(analysis$ensmean), s, ncol(analysis$ensmean)/s)), 1:2, mean)), c(nrow(analysis$ensmean), s, ncol(analysis$ensmean)/s))
       # analysis.anom$ensmean <- (array(analysis.anom$ensmean,c(dim(analysis.anom$ensmean)[1],dim(analysis.anom$ensmean)[2]*dim(analysis.anom$ensmean)[3])))  
       # 
       # echam.anom <- echam
       # echam.anom$data <- array(echam$data - as.vector(apply(array(echam$data, c(nrow(echam$data), s, ncol(echam$data)/s,30)), 1:2, mean)), c(nrow(echam$data), s, ncol(echam$data)/s, 30))
       # echam.anom$data <- (array(echam.anom$data,c(dim(echam.anom$data)[1],dim(echam.anom$data)[2]*dim(echam.anom$data)[3],dim(echam.anom$data)[4])))
       # echam.anom$ensmean <- array(echam$ensmean - as.vector(apply(array(echam$ensmean, c(nrow(echam$ensmean), s, ncol(echam$ensmean)/s)), 1:2, mean)), c(nrow(echam$ensmean), s, ncol(echam$ensmean)/s))
       # echam.anom$ensmean <- (array(echam.anom$ensmean,c(dim(echam.anom$ensmean)[1],dim(echam.anom$ensmean)[2]*dim(echam.anom$ensmean)[3])))  
       
       if (vali) {
         validate.clim <- validate
         validate.clim$data <- apply(array(validate$data, c(nrow(validate$data), s, 
                                                            ncol(validate$data)/s)), 1:2, mean, na.rm=T)
         validate.anom <- validate
         validate.anom$data <- array(validate$data - as.vector(apply(array(validate$data, 
                                                                           c(nrow(validate$data), s, ncol(validate$data)/s)), 1:2, 
                                                                     mean)), c(nrow(validate$data), ncol(validate$data)))
         validate.anom$ensmean <- array(validate$ensmean - as.vector(apply(array(validate$ensmean, 
                                                                                 c(nrow(validate$ensmean), s, ncol(validate$ensmean)/s)), 1:2, 
                                                                           mean)), c(nrow(validate$ensmean), ncol(validate$ensmean)))
       }
       
       calibrate.clim <- calibrate
       calibrate.clim$data <- apply(array(calibrate$data, c(nrow(calibrate$data), s, ncol(calibrate$data)/s)), 1:2, mean, na.rm=T)
       calibrate.anom <- calibrate
       calibrate.anom$data <- array(calibrate$data - as.vector(calibrate.clim$data), c(nrow(calibrate$data), ncol(calibrate$data)))
       
       print("calc anomalies")
       print(proc.time() - ptm1)
       
       if (every2grid) {
         if (monthly_out) {
           save.image(file=paste0("../data/image/",expname,"/prepplot_validation_image_",
                                  syr,"-",eyr,"_monthly_2ndgrid.Rdata"))  
         } else {
           save.image(file=paste0("../data/image/",expname,"/prepplot_validation_image_",
                                  syr,"-",eyr,"_seasonal_2ndgrid.Rdata"))
         }  
       } else {
         if (monthly_out) {
           save.image(file=paste("../data/image/",expname,"/prepplot_validation_image_",
                                 syr,"-",eyr,"_monthly.Rdata",sep=""))
         } else {
           save.image(file=paste("../data/image/",expname,"/prepplot_validation_image_",
                                 syr,"-",eyr,"_seasonal.Rdata",sep=""))
         }
       }
     } # end load_prepplot
########################################################################################################
     if (calc_vali_stat){
       # year 1902-1980 only because of missing validation data
       # if (eyr>1980) {
       #   print('ATTENTION: period shorted!')
       #   echam$data <- echam$data[,1:160,]
       #   echam.anom$data <- echam.anom$data[,1:160,]
       #   echam$ensmean <- echam$ensmean[,1:160]
       #   echam.anom$ensmean <- echam.anom$ensmean[,1:160]
       #   analysis$data <- analysis$data[,1:160,]
       #   analysis.anom$data <- analysis.anom$data[,1:160,]
       #   analysis$ensmean <- analysis$ensmean[,1:160]
       #   analysis.anom$ensmean <- analysis.anom$ensmean[,1:160]
       #   validate$data <- validate$data[,1:160]
       #   validate.anom$data <- validate.anom$data[,1:160]
       #   validate$ensmean <- validate$ensmean[,1:160]
       #   validate.anom$ensmean <- validate.anom$ensmean[,1:160]
       # }
       
       # if ((instrumental) & (!real_proxies)) {
       #   # find and cut overlapping period for cru validation
       #   #cruatprox.arr <- array(NA,c(length(proxies$lon), length(inst$time)))
       #   #cruatprox.arr.allts <- array(NA,c(length(proxies.allts$lon), length(inst.allts$time)))
       #   cruatprox.arr.allts <- array(NA,c(length(calibrate$lon), length(calibrate$time)))
       #   if ((cru_vali) | (ncep_vali) | (recon_vali)) {
       #     for(i in 1:length(calibrate$lon)){
       #       plon <- calibrate$lon[i]
       #       plat <- calibrate$lat[i]
       #       clon <- valiall$lon
       #       clat <- valiall$lat  
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
       #         #      cruatprox.arr[i,] <- cruall$data[m,]
       #         #      cruatprox.arr.allts[i,] <- cruall.allts$data[m,]
       #         if (calibrate$names[i]=='temp2') {
       #           cruatprox.arr.allts[i,] <- validate$data[m,]
       #         } else if (calibrate$names[i]=='precip') {
       #           cruatprox.arr.allts[i,] <- validate$data[(m+length(validate$lon)),]
       #         } else if (calibrate$names[i]=='slp') {
       #           cruatprox.arr.allts[i,] <- validate$data[(m+(2*length(validate$lon))),]
       #         } else {
       #           cruatprox.arr.allts[i,] <- NA
       #         }  
       #       }
       #     }  
       #   } else {
       #     #  cruatprox.arr[i,] <- NA
       #     cruatprox.arr.allts[i,] <- NA
       #   }
       #   
       #   pdims2 <- c(prod(dim(calibrate$data)[1],s), (dim(calibrate$data)[2])/s)
       #   prox.cru.corr <- array(diag(cor(t(array(calibrate$data,pdims2)), 
       #                                   t(array(cruatprox.arr.allts, pdims2)),use='pairwise.complete.obs')), 
       #                          c(dim(calibrate$data)[1],s)) # just check temp, 
       #   prox.cru.bias <- array(apply(array(calibrate$data,pdims2),1,mean, na.rm=T) - 
       #                            (apply(array(cruatprox.arr.allts,pdims2),1,mean, na.rm=T)), 
       #                          c(dim(calibrate$data)[1],s))
       #   prox.cru.corr.temp <- prox.cru.corr[calibrate$names=='temp2',]
       #   prox.cru.corr.precip <- prox.cru.corr[calibrate$names=='precip',]
       #   prox.cru.corr.slp <- prox.cru.corr[calibrate$names=='slp',]
       #   prox.cru.bias.temp <- prox.cru.bias[calibrate$names=='temp2',]
       #   prox.cru.bias.precip <- prox.cru.bias[calibrate$names=='precip',]
       #   prox.cru.bias.slp <- prox.cru.bias[calibrate$names=='slp',]
       # }
       # 
       # #if (!tps_only) {
       # lenvar2 <- length(c(which(validate$names=="temp2"), which(validate$names=="precip"), 
       #                      which(validate$names=="slp")))
       # #lenvar3 <- length(c(which(echam$names=="temp2"), which(echam$names=="precip"), 
       # #                    which(echam$names=="slp")))
       # # if (lenvar2!=lenvar3) {print("ACHTUNG: validate and echam have different dimensions")}  
       # analysis_noindex=analysis
       # analysis_noindex$data=analysis$data[1:lenvar2,,]
       # analysis_noindex$ensmean=analysis$ensmean[1:lenvar2,]
       # analysis_noindex$names=analysis$names[1:lenvar2]
       # analysis_noindex$lon=analysis$lon[1:lenvar2]
       # analysis_noindex$lat=analysis$lat[1:lenvar2]
       # echam_noindex=echam
       # echam_noindex$data=echam$data[1:lenvar2,,]
       # echam_noindex$ensmean=echam$ensmean[1:lenvar2,]
       # echam_noindex$names=echam$names[1:lenvar2]
       # echam_noindex$lon=echam$lon[1:lenvar2]
       # echam_noindex$lat=echam$lat[1:lenvar2]
       # analysis.anom_noindex=analysis.anom
       # analysis.anom_noindex$data=analysis.anom$data[1:lenvar2,,]
       # analysis.anom_noindex$ensmean=analysis.anom$ensmean[1:lenvar2,]
       # analysis.anom_noindex$names=analysis.anom$names[1:lenvar2]
       # analysis.anom_noindex$lon=analysis.anom$lon[1:lenvar2]
       # analysis.anom_noindex$lat=analysis.anom$lat[1:lenvar2]
       # echam.anom_noindex=echam.anom
       # echam.anom_noindex$data=echam.anom$data[1:lenvar2,,]
       # echam.anom_noindex$ensmean=echam.anom$ensmean[1:lenvar2,]
       # echam.anom_noindex$names=echam.anom$names[1:lenvar2]
       # echam.anom_noindex$lon=echam.anom$lon[1:lenvar2]
       # echam.anom_noindex$lat=echam.anom$lat[1:lenvar2]
       ##rmse <- lapply(analysis, rmse_fun, y=validate, seas=s)
       #rmse <- rmse_fun(analysis_noindex, y=validate, seas=s)
       #rmse.ech <- rmse_fun(echam_noindex, y=validate, seas=s)
       if (vali) {
         #  rmse <- rmse_fun(analysis_noindex, y=validate, seas=s)
         #  rmse.ech <- rmse_fun(echam_noindex, y=validate, seas=s)
         rmse <- rmse_fun(analysis, y=validate, seas=s)
         rmse.ech <- rmse_fun(echam, y=validate, seas=s)
         # vali.clim <- multiyear_seas_average(validate, seas=s)
         # vali.clim.data <- vali.clim$data
         # vali.clim.time <- vali.clim$time
         # for (i in 1:((ncol(validate$data)/s)-1)) {
         #   vali.clim.data <- cbind(vali.clim.data,vali.clim$data)
         #   vali.clim.time <- c(vali.clim.time,vali.clim$time)
         # }
         # vali.clim$data <- vali.clim.data
         # vali.clim$time <- vali.clim.time
         # rmse.clim <- rmse_fun(vali.clim, y=validate, seas=s)
         #  rmseobserr <- rmse_obserr_fun(analysis, y=validate, seas=s)
         #  rmseobserr.ech <- rmse_obserr_fun(echam, y=validate, seas=s)
       }
       ##if (cru_vali) { 
       # validation of bias makes little sense with regional reconstruction
       ##}  
       ##RE <- lapply(rmse, RE_fun, y=rmse.ech)
       if (vali) { 
         #  rmse.anom <- rmse_fun(analysis.anom_noindex, y=validate.anom, seas=s)
         #  rmse.ech.anom <- rmse_fun(echam.anom_noindex, y=validate.anom, seas=s)    
         rmse.anom <- rmse_fun(analysis.anom, y=validate.anom, seas=s)
         rmse.ech.anom <- rmse_fun(echam.anom, y=validate.anom, seas=s) 
         # vali.clim.anom <- multiyear_seas_average(validate.anom, seas=s)
         # vali.clim.anom.data <- vali.clim.anom$data
         # vali.clim.anom.time <- vali.clim.anom$time
         # for (i in 1:((ncol(validate.anom$data)/s)-1)) {
         #   vali.clim.anom.data <- cbind(vali.clim.anom.data,vali.clim.anom$data)
         #   vali.clim.anom.time <- c(vali.clim.anom.time,vali.clim.anom$time)
         # }
         # vali.clim.anom$data <- vali.clim.anom$ensmean <- vali.clim.anom.data
         # vali.clim.anom$time <- vali.clim.anom.time
         # rmse.clim.anom <- rmse_fun(vali.clim.anom, y=validate.anom, seas=s)
         echam.clim.anom <- multiyear_seas_average(echam.anom, seas=s)
         echam.clim.anom.data <- echam.clim.anom$data
         echam.clim.anom.time <- echam.clim.anom$time
         for (i in 1:((ncol(echam.anom$data)/s)-1)) {
           print(i)
           echam.clim.anom.data <- cbind(echam.clim.anom.data,echam.clim.anom$data)
           echam.clim.anom.time <- c(echam.clim.anom.time,echam.clim.anom$time)
         }
         echam.clim.anom$data <- echam.clim.anom$ensmean <- echam.clim.anom.data
         echam.clim.anom$time <- echam.clim.anom.time
         rmse.echam.clim.anom <- rmse_fun(echam.clim.anom, y=validate.anom, seas=s)
         echam.clim <- multiyear_seas_average(echam, seas=s)
         echam.clim.data <- echam.clim$data
         echam.clim.time <- echam.clim$time
         for (i in 1:((ncol(echam$data)/s)-1)) {
           print(i)
           echam.clim.data <- cbind(echam.clim.data,echam.clim$data)
           echam.clim.time <- c(echam.clim.time,echam.clim$time)
         }
         echam.clim$data <- echam.clim$ensmean <- echam.clim.data
         echam.clim$time <- echam.clim.time
         rmse.echam.clim <- rmse_fun(echam.clim, y=validate, seas=s)
         #  rmseobserr.anom <- rmse_obserr_fun(analysis.anom, y=validate.anom, seas=s)
         #  rmseobserr.ech.anom <- rmse_obserr_fun(echam.anom, y=validate.anom, seas=s)    
         RE <- RE_fun(rmse, y=rmse.ech)
         #  RE.clim <- RE_fun(rmse, y=rmse.clim)
         RE.anom <- RE_fun(rmse.anom, y=rmse.ech.anom)
         #  RE.clim.anom <- RE_fun(rmse.anom, y=rmse.clim.anom)
         RE.echam.clim.anom <- RE_fun(rmse.anom, y=rmse.echam.clim.anom)
         #  RE.ech.clim.anom <- RE_fun(rmse.ech.anom, y=rmse.clim.anom)
         RE.echam.clim <- RE_fun(rmse, y=rmse.echam.clim)
       }
       ##corr <- lapply(analysis, corr_fun, y=validate, seas=s)
       if (vali) {
         #  corr <-  corr_fun(analysis_noindex, y=validate, seas=s)
         #  corr.ech <- corr_fun(echam_noindex, y=validate, seas=s)
         corr <-  corr_fun(analysis, y=validate, seas=s)
         corr.ech <- corr_fun(echam, y=validate, seas=s)
       }
       ##bias <- lapply(analysis, bias_fun, y=validate, seas=s)
       if (vali) {
         #  bias <- bias_fun(analysis_noindex, y=validate, seas=s)
         #  bias.ech <- bias_fun(echam_noindex, y=validate, seas=s)
         bias <- bias_fun(analysis, y=validate, seas=s)
         bias.ech <- bias_fun(echam, y=validate, seas=s)
       }
       #if (check_dist) {
       #  rank <- lapply(analysis_noindex, rank_fun, y=validate, seas=s)
       #  rank.ech <- rank_fun(echam, y=validate, seas=s)
       #  rmse.dist <- lapply(ana.dist, rmse_fun, y=validate, seas=s)
       #  corr.dist <- lapply(ana.dist, corr_fun, y=validate, seas=s)
       #  bias.dist <- lapply(ana.dist, bias_fun, y=validate, seas=s)
       #  rank.dist <- lapply(ana.dist, rank_fun, y=validate, seas=s)
       #  RE.dist <- lapply(rmse.dist, RE_fun, y=rmse.ech)
       #}
       
       data.dim <- dim(echam$data)[c(1,2,2,3)]
       data.dim[2:3] <- c(s,data.dim[3]/s)
       ens.dim <- c(nrow(echam$ensmean), s, ncol(echam$ensmean)/s)
       
       ech.spread <- apply(sqrt(apply(array(echam[['data']] - as.vector(echam[['ensmean']]), 
                                            data.dim)**2, 1:3, mean,na.rm=T)), 1:2, mean, na.rm=T)
       #ana.spread <- lapply(analysis, function(x) apply(sqrt(apply(array(x[['data']] - 
       #                  as.vector(x[['ensmean']]), data.dim)**2, 1:3, mean, na.rm=T)), 
       #                  1:2, mean, na.rm=T))
       ana.spread <- apply(sqrt(apply(array(analysis[['data']] - as.vector(analysis[['ensmean']]), 
                                            data.dim)**2, 1:3, mean,na.rm=T)), 1:2, mean, na.rm=T)
       #analysis.pre1800 <-
       #analysis.post1800 <-
       #analysis.post1900 <-
       
       
       
       # spread to error ratio rechnen (also das zeitliche mittel der Ensemblevarianz dividiert 
       # durch den mittleren quadrierten Fehler der Beobachtungen vom ensemble mean). Diese sollte 
       # ca. 1 sein wenn die Fehlerabschätzung ok ist (kleiner 1 deutet dann auf over-fitting hin, 
       # grösser 1 auf zu strikte Lokalisierung oder Fehler in R, der Kovarianz der Beobachtungsfehler)
       # momentan nur 1902-2004.
       # error-spread ratio definition of Weigel (http://onlinelibrary.wiley.com/store/10.1002/9781119960003.ch8/asset/ch8.pdf?v=1&t=ir1stg2j&s=d8c8a7df156a195e5631d4d086c22040702cd267):
       # Here we define the spread-error rate as the square root of the ratio of mean ensemble variance to the mean squared error of the ensemble mean with the verifying observations. 
       #install.packages("easyVerification")
       #library("easyVerification")
       #tm <- toymodel()
       #ens: n x k matrix of n forecasts for k ensemble members
       #obs:	vector with n verifying observations
       #FairSprErr(tm$fcst, tm$obs)
       # let's do the same for each grid box of your analysis:
       #   n are the 100 forecasts for the years 1901-2000
       #   k are the 30 ensemble members
       if (every2grid) {
         load(file="../data/cru/cru4_ens_sd_2ndgrid.Rdata")
       } else {
         load(file="../data/cru/cru4_ens_sd.Rdata")  
       }
       land <- cru4_oct_apr$lsm.i
       obs.spread.win <- cru4_oct_apr$data[,1:80,1]
       obs.spread.sum <- cru4_may_sep$data[,1:80,1]
       obs.spread.win[obs.spread.win==0] <- NA
       obs.spread.sum[obs.spread.sum==0] <- NA
       sprerr.c.sum <- rep(NA, nrow(obs.spread.sum))
       sprerr.c.win <- rep(NA, nrow(obs.spread.win))
       sprerr.sum <- rep(NA, nrow(obs.spread.sum))
       sprerr.win <- rep(NA, nrow(obs.spread.win))
       ech.sprerr.c.sum <- rep(NA, nrow(obs.spread.sum))
       ech.sprerr.c.win <- rep(NA, nrow(obs.spread.win))
       ech.sprerr.sum <- rep(NA, nrow(obs.spread.sum))
       ech.sprerr.win <- rep(NA, nrow(obs.spread.win))
       e.tmp <- array(echam.anom$data[,,],dim=c(dim(echam.anom$data[,,])[1],2,
                                                (dim(echam.anom$data[,,])[2]/2), dim(echam.anom$data[,,])[3])) 
       a.tmp <- array(analysis.anom$data[,,],dim=c(dim(analysis.anom$data[,,])[1],2,
                                                   (dim(analysis.anom$data[,,])[2]/2), dim(analysis.anom$data[,,])[3])) 
       e.sum <- e.tmp[,2,,]
       e.win <- e.tmp[,1,,]
       rm(e.tmp)
       a.sum <- a.tmp[,2,,]
       a.win <- a.tmp[,1,,]
       rm(a.tmp)
       v.tmp <- array(validate.anom$data[,],dim=c(dim(validate.anom$data[,])[1],2,
                                                  (dim(validate.anom$data[,])[2]/2)))
       v.sum <- v.tmp[,2,]
       v.win <- v.tmp[,1,]
       rm(v.tmp)
       
       for (i in which(analysis.anom$names=="temp2")) {
         if (i %in% cru4_oct_apr$lsm.i) {
           sprerr.sum[i] <- mean(FairSprErr2(ens=a.sum[i,,],obs=v.sum[i,]),na.rm=T)
           sprerr.win[i] <- mean(FairSprErr2(ens=a.win[i,,],obs=v.win[i,]),na.rm=T)
         } else {  
           sprerr.sum[i] <- sprerr.win[i] <- NA
         }  
       }
       sprerr <- cbind(sprerr.win,sprerr.sum) # corrected for obs. uncertainties
       
       for (i in which(echam$names=="temp2")) {
         if (i %in% cru4_oct_apr$lsm.i) {
           ech.sprerr.sum[i] <- mean(FairSprErr2(ens=e.sum[i,,],obs=v.sum[i,]),na.rm=T)
           ech.sprerr.win[i] <- mean(FairSprErr2(ens=e.win[i,,],obs=v.win[i,]),na.rm=T)
         } else {  
           ech.sprerr.sum[i] <- ech.sprerr.win[i] <- NA
         }  
       }
       ech.sprerr <- cbind(ech.sprerr.win,ech.sprerr.sum) # corrected for obs. uncertainties
       
       # This does NOT account for observations errors. Hence we apply the correction suggested by Bowler (2007)
       # (http://research.metoffice.gov.uk/research/nwp/publications/papers/technical_reports/reports/506.pdf):
       # true rmse (RMSE_t), given an additive observation error can be estimated with the SD of the 
       # observations, which are estimate from the CRUTEM ensemble here: RMSE_t = sqrt(RMSE**2 - o_sd**2)
       # obs are the cru validation data for the specific grid box and the same 100 years
       # a.tmp <- array(analysis$data[,,],dim=c(dim(analysis$data[,,])[1],2,
       #            (dim(analysis$data[,,])[2]/2), dim(analysis$data[,,])[3])) 
       #            # added 1:160 because analysis matrix has 1 yr too much and 
       #            # validate$data has no data after 1980
       # a.sum <- a.tmp[,2,,]
       # a.win <- a.tmp[,1,,]
       # rm(a.tmp)
       # v.tmp <- array(validate$data[,1:160],dim=c(dim(validate$data[,1:160])[1],2,
       #            (dim(validate$data[,1:160])[2]/2)))
       # v.sum <- v.tmp[,2,]
       # v.win <- v.tmp[,1,]
       # rm(v.tmp)
       
       for (i in which(analysis.anom$names=="temp2")) {
         if (i %in% cru4_oct_apr$lsm.i) {
           sprerr.c.sum[i] <- mean(FairSprErr2(ens=a.sum[i,,],obs=v.sum[i,],obs.sd=obs.spread.sum[i,]),na.rm=T)
           sprerr.c.win[i] <- mean(FairSprErr2(ens=a.win[i,,],obs=v.win[i,],obs.sd=obs.spread.win[i,]),na.rm=T)
         } else {  
           sprerr.c.sum[i] <- sprerr.c.win[i] <- NA
         }
       }
       sprerr.corr <- cbind(sprerr.c.win,sprerr.c.sum) # corrected for obs. uncertainties
       
       for (i in which(echam$names=="temp2")) {
         if (i %in% cru4_oct_apr$lsm.i) {
           ech.sprerr.c.sum[i] <- mean(FairSprErr2(ens=e.sum[i,,],obs=v.sum[i,],
                                                   obs.sd=obs.spread.sum[i,]),na.rm=T)
           ech.sprerr.c.win[i] <- mean(FairSprErr2(ens=e.win[i,,],obs=v.win[i,],
                                                   obs.sd=obs.spread.win[i,]),na.rm=T)
         } else {  
           ech.sprerr.c.sum[i] <- ech.sprerr.c.win[i] <- NA
         }
       }
       ech.sprerr.corr <- cbind(sprerr.c.win,sprerr.c.sum) # corrected for obs. uncertainties
       
       
       
       
       
       
       
       
       if (vali) {
         #   ereliable <- tapply(apply(echam_noindex$data[1:(dim(validate$data)[1]),,] > 
         #                   as.vector(validate$data),1:2,mean), rep(echam_noindex$names,
         #                   length=length(validate$data)), table)
         #   areliable <- tapply(apply(analysis_noindex$data[1:(dim(validate$data)[1]),,] > 
         #                   as.vector(validate$data),1:2,mean), rep(analysis_noindex$names,
         #                   length=length(validate$data)), table)
         #  ereliable <- tapply(apply(echam$data[1:(dim(validate$data)[1]),1:160,] > 
         #                     as.vector(validate$data[,1:160]),1:2,mean), rep(echam$names,
         #                     length=length(validate$data[,1:160])), table)
         #  areliable <- tapply(apply(analysis$data[1:(dim(validate$data)[1]),1:160,] > 
         #                     as.vector(validate$data[,1:160]),1:2,mean), rep(analysis$names,
         #                     length=length(validate$data[,1:160])), table) 
         ereliable <- tapply(apply(echam$data[1:(dim(validate$data)[1]),,] > 
                                     as.vector(validate$data[,]),1:2,mean), rep(echam$names,
                                                                                length=length(validate$data[,])), table)
         areliable <- tapply(apply(analysis$data[1:(dim(validate$data)[1]),,] > 
                                     as.vector(validate$data[,]),1:2,mean), rep(analysis$names,
                                                                                length=length(validate$data[,])), table)   
         # this does not take model biases into account, hence better with anomalies
         ereliable.anom <- tapply(apply(echam.anom$data[1:(dim(validate.anom$data)[1]),,] > 
                                          as.vector(validate.anom$data[,]),1:2,mean), rep(echam.anom$names,
                                                                                          length=length(validate.anom$data[,])), table)
         areliable.anom <- tapply(apply(analysis.anom$data[1:(dim(validate.anom$data)[1]),,] > 
                                          as.vector(validate.anom$data[,]),1:2,mean), rep(analysis.anom$names,
                                                                                          length=length(validate.anom$data[,])), table) 
         #   ereliable.anom <- tapply(apply(echam.anom$data[1:(dim(validate.anom$data)[1]),1:160,] > 
         #                      as.vector(validate.anom$data[,1:160]),1:2,mean), rep(echam.anom$names,
         #                      length=length(validate.anom$data[,1:160])), table)
         #   areliable.anom <- tapply(apply(analysis.anom$data[1:(dim(validate.anom$data)[1]),1:160,] > 
         #                      as.vector(validate.anom$data[,1:160]),1:2,mean), rep(analysis.anom$names,
         #                      length=length(validate.anom$data[,1:160])), table) 
         # rank histogram that takes error in instrumental data into account
         #rh = apply(ens + rnorm(length(ens), mean=0, sd=o_sd) < obs, 1, sum)  
         ## unter der Annahme das ens eine n x nens matrix und obs ein entsprechender Vektor der Länge n ist
         #obsdat <- array(validate$data, c(nrow(validate$data), s, ncol(validate$data)/s, 
         #            length(validate$data)/prod(dim(validate$data)[1:2])))
         # ERROR: next line needs to be changed to contain SD of CRUTEM4 ensemble, 
         #        i.e. instr. obs. error instead of obs. SD!!!
         # obs_sd <- apply(obsdat,c(1,2,4), sd, na.rm=T)
         # NEW from CRU ensemble
         obs_sd <- cbind(apply(obs.spread.win,1,mean),apply(obs.spread.sum,1,mean))
         obs.sd.nona <- obs_sd
         obs.sd.nona[is.na(obs_sd)] <- 0
         #obs_sd_mean <- mean(obs_sd)
         rdat <- array(NA,dim(analysis.anom$data[which(analysis.anom$names=="temp2"),,]))
         #rdat <- array(NA,dim(analysis.anom$data[which(analysis.anom$names=="temp2"),1:160,]))
         #for (i in 1:dim(obs_sd)[1]) {
         #  if (i %% 100 == 0) {print(i)}
         for (j in 1:dim(analysis.anom$data[,,])[2]) {
           #for (j in 1:dim(analysis.anom$data[,1:160,])[2]) {
           print(j)
           if (is.even(j)) {k=2} else {k=1}
           #if (!is.na(obs_sd[i,k])) {
           #rdat[i,j,] <- rnorm(dim(analysis.anom$data)[3],mean=0,sd=obs_sd[i,k])
           rdat[,j,] <- rnorm(dim(obs_sd)[1]*dim(analysis.anom$data)[3],mean=0,sd=as.vector(obs.sd.nona[,k]))
           #} else {
           #    rdat[,j,] <- rnorm(dim(obs_sd)[1]*dim(analysis.anom$data)[3],mean=0,sd=obs_sd_mean)
           #  }
           #}
         }
         atala <- analysis.anom$data[which(analysis.anom$names=="temp2"),,] + rdat
         etala <- echam.anom$data[which(echam.anom$names=="temp2"),,] + rdat
         erel_obserr.anom <- tapply(apply(etala[1:(dim(obs_sd)[1]),,] > 
                                            as.vector(validate.anom$data[which(validate.anom$names=="temp2"),]),
                                          1:2,mean), rep(echam.anom$names,length=length(validate.anom$data[which(
                                            validate.anom$names=="temp2"),])), table)
         arel_obserr.anom <- tapply(apply(atala[1:(dim(obs_sd)[1]),,] > 
                                            as.vector(validate.anom$data[which(validate.anom$names=="temp2"),]),
                                          1:2,mean), rep(analysis.anom$names,length=length(validate.anom$data[which(
                                            validate.anom$names=="temp2"),])), table) 
         #   atala <- analysis.anom$data[which(analysis.anom$names=="temp2"),1:160,] + rdat
         #   etala <- echam.anom$data[which(echam.anom$names=="temp2"),1:160,] + rdat
         #   erel_obserr.anom <- tapply(apply(etala[1:(dim(obs_sd)[1]),1:160,] > 
         #                         as.vector(validate.anom$data[which(validate.anom$names=="temp2"),1:160]),
         #                         1:2,mean), rep(echam.anom$names,length=length(validate.anom$data[which(
         #                         validate.anom$names=="temp2"),1:160])), table)
         #   arel_obserr.anom <- tapply(apply(atala[1:(dim(obs_sd)[1]),1:160,] > 
         #                         as.vector(validate.anom$data[which(validate.anom$names=="temp2"),1:160]),
         #                         1:2,mean), rep(analysis.anom$names,length=length(validate.anom$data[which(
         #                         validate.anom$names=="temp2"),1:160])), table) 
       }
       
       # correlation of observation errors
       # 1. difference between obs and analysis ensemble mean at obs. locations 
       # run for 1901-1979 period to have non-changing instr obs network
       calibrate$lon <- as.vector(calibrate$lon[,1])
       calibrate$lat <- as.vector(calibrate$lat[,1])
       calibrate$sour <- as.vector(calibrate$sour[,1])
       calibrate$names <- as.vector(calibrate$names[,1])
       calibrate.anom$lon <- as.vector(calibrate.anom$lon[,1])
       calibrate.anom$lat <- as.vector(calibrate.anom$lat[,1])
       calibrate.anom$sour <- as.vector(calibrate.anom$sour[,1])
       calibrate.anom$names <- as.vector(calibrate.anom$names[,1])
       
       #obserr <- diff_fun(analysis.anom,calibrate.anom,seas=s)
       #pos <- !is.na(apply(obserr$ensmean,1,sum))
       # compute_H <- function(stations, echam, threshold=700){
       #   H <- array(0, c(nrow(stations$data), nrow(echam$data)))
       #   for (i in seq(stations$lon)){
       #     dist <- compute_dist(echam$lon, echam$lat, stations$lon[i], stations$lat[i])
       #     H[i, which.min(dist)] <- if (min(dist) < threshold) 1 else 0
       #   }
       #   H
       # }
       
       stat=calibrate.anom
       stat$data=stat$data[which(stat$sour=="inst"),]
       stat$lon=stat$lon[which(stat$sour=="inst")]
       stat$lat=stat$lat[which(stat$sour=="inst")]
       stat$names=stat$names[which(stat$sour=="inst")]
       stat$sour=stat$sour[which(stat$sour=="inst")]
       ctmp <- array(stat$data,c((dim(stat$data)[1]/6),6,
                                 dim(stat$data)[2]))  
       stat$data <- apply(ctmp,c(1,3),mean,na.rm=T)
       stat$names <- stat$names[1:dim(stat$data)[1]] 
       stat$lon <- stat$lon[1:dim(stat$data)[1]]
       stat$lat <- stat$lat[1:dim(stat$data)[1]]
       
       dlist=NA
       for(i in 1:length(stat$lon)){
         plon <- stat$lon[i]
         plat <- stat$lat[i]
         clon <- rep(validate.anom$lon[1:length(validate.anom$lsm.i)],3)
         clat <- rep(validate.anom$lat[1:length(validate.anom$lsm.i)],3)
         k=which(abs(clon-plon+0.001)==min(abs(clon-plon+0.001))) 
         # +0.001 to avoid find 2 locations with same distance
         l=which(abs(clat-plat)==min(abs(clat-plat))) 
         m=k[which(match(k,l)>0)]
         if (length(m) > 0) {
           if (stat$names[i]=="temp2") {
             dlist[i]=m[1]
           } else if (stat$names[i]=="slp") {
             dlist[i]=m[3]  
           }
         } else {
           dlist[i]=NA
         }
       }
       
       #gpos <- getgridboxnum(calibrate.anom,analysis.anom)
       obserr <- matrix(NA,nrow=length(dlist),ncol=ncol(validate.anom$ensmean))
       for (i in 1:length(dlist)) {
         m <- dlist[i]
         obserr[i,] <- validate.anom$ensmean[m,]-stat$data[i,] 
       }
       # pos=NA
       # for (i in 1:nrow(obserr)) {
       #   if (length(which(!is.na(obserr[i,])))>0) {pos=c(pos,i)}   # >30
       # }
       # pos=pos[-1]
       # obserr2=obserr[pos,]
       cmat <- cor(t(obserr),use="pairwise.complete.obs")   #obserr2
       # there are correlated observations
       # hist(cmat,breaks=seq(-1,1,0.2))
       # for (i in 1:nrow(cmat)) {
       #   print(i)
       #   print(which(cmat[i,]>0.8))
       # }
       # loc=c(194,195,199,232)
       # plot(obserr[loc[1],],ty='l')
       # lines(obserr[loc[2],],ty='l',col='red')
       # lines(obserr[loc[3],],ty='l',col='blue')
       # lines(obserr[loc[4],],ty='l',col='green')
       # plot(stat$lon[loc],stat$lat[loc])
       # cor(obserr[loc[1],],obserr[loc[2],],use="pairwise.complete.obs")
       # cbind(stat$lon[loc],stat$lat[loc]) 
       # dlist[loc]
       # cbind(echam$lon[dlist[loc]],echam$lat[dlist[loc]])
       # par(mfrow=c(4,1))
       # plot(stat$data[loc[1],seq(2,40,2)],ty='l')
       # lines(stat$data[loc[2],seq(2,40,2)],ty='l',col='red')
       # lines(stat$data[loc[3],seq(2,40,2)],ty='l',col='blue')
       # lines(stat$data[loc[4],seq(2,40,2)],ty='l',col='green')
       # plot(validate.anom$ensmean[dlist[loc[1]],seq(2,40,2)],ty='l')
       # lines(validate.anom$ensmean[dlist[loc[2]],seq(2,40,2)],ty='l',col='red')
       # lines(validate.anom$ensmean[dlist[loc[3]],seq(2,40,2)],ty='l',col='blue')
       # lines(validate.anom$ensmean[dlist[loc[4]],seq(2,40,2)],ty='l',col='green')
       # plot(analysis.anom$ensmean[dlist[loc[1]],seq(2,40,2)],ty='l')
       # lines(analysis.anom$ensmean[dlist[loc[2]],seq(2,40,2)],ty='l',col='red')
       # lines(analysis.anom$ensmean[dlist[loc[3]],seq(2,40,2)],ty='l',col='blue')
       # lines(analysis.anom$ensmean[dlist[loc[4]],seq(2,40,2)],ty='l',col='green')
       # plot(obserr[loc[1],seq(2,40,2)],ty='l')
       # lines(obserr[loc[2],seq(2,40,2)],ty='l',col='red')
       # lines(obserr[loc[3],seq(2,40,2)],ty='l',col='blue')
       # lines(obserr[loc[4],seq(2,40,2)],ty='l',col='green')
       
       
       
       print("calc validation statistics")
       print(proc.time() - ptm1)
       
       
       
       
       
       # extract the indices from output 
       # 1. Analysis indices
       H.giorgi <- compute_giorgi_H_sixmon(giorgi, echam) #, numvar=3) # 3 vars temp, precip, slp
       #  H.giorgi <- compute_giorgi_H_v2(giorgi, echam_noindex) 
       
       if (!tps_only & ind_ECHAM) {
         # add weights for indices # NOT WORKING because indices are not simple means 
         #tmpweights <- cos(echam$latstream/180*pi)
         #tmpweights <- tmpweights/sum(tmpweights)
         #tmpstream <- rep(0,dim(echam$ensmean)[1])
         #tmpstream[which(echam$names=='stream')] <- tmpweights
         #Htmp <- rbind(H.giorgi,tmpstream,tmpgph100,tmpgph300,tmpu200,tmpomega,tmpomega_2) 
         indices <- c('NH.temp2', 'NH.precip', 'NH.slp', 'NEU.temp2', 'NEU.precip', 'NEU.slp',
                      'DIMI', 'HC', 'SJ', 'Z100', 'Z300', 'PWC')
         Hind <- matrix(0,nrow=12,ncol=nrow(echam$data))
         Hind[1,which(echam$names=="temp2")] <- H.giorgi[which(giorgi.short == 'NH'),
                                                         which(echam$names=="temp2")]
         Hind[2,which(echam$names=="precip")] <- H.giorgi[which(giorgi.short == 'NH'),
                                                          which(echam$names=="precip")]
         Hind[3,which(echam$names=="slp")] <- H.giorgi[which(giorgi.short == 'NH'),
                                                       which(echam$names=="slp")]
         Hind[4,which(echam$names=="temp2")] <- H.giorgi[which(giorgi.short == 'NEU'),
                                                         which(echam$names=="temp2")]
         Hind[5,which(echam$names=="precip")] <- H.giorgi[which(giorgi.short == 'NEU'),
                                                          which(echam$names=="precip")]
         Hind[6,which(echam$names=="slp")] <- H.giorgi[which(giorgi.short == 'NEU'),
                                                       which(echam$names=="slp")] 
         
         #   Hind <- matrix(0,nrow=12,ncol=nrow(echam_noindex$data))
         #   Hind[1,which(echam_noindex$names=="temp2")] <- H.giorgi[which(giorgi.short == 'NH'),
         #                                                           which(echam_noindex$names=="temp2")]
         #   Hind[2,which(echam_noindex$names=="precip")] <- H.giorgi[which(giorgi.short == 'NH'),
         #                                                           which(echam_noindex$names=="precip")]
         #   Hind[3,which(echam_noindex$names=="slp")] <- H.giorgi[which(giorgi.short == 'NH'),
         #                                                            which(echam_noindex$names=="slp")]
         #   Hind[4,which(echam_noindex$names=="temp2")] <- H.giorgi[which(giorgi.short == 'NEU'),
         #                                                           which(echam_noindex$names=="temp2")]
         #   Hind[5,which(echam_noindex$names=="precip")] <- H.giorgi[which(giorgi.short == 'NEU'),
         #                                                            which(echam_noindex$names=="precip")]
         #   Hind[6,which(echam_noindex$names=="slp")] <- H.giorgi[which(giorgi.short == 'NEU'),
         #                                                         which(echam_noindex$names=="slp")] 
         ##  Hind <- rbind(H.giorgi[c(which(giorgi.short == 'NH')+ c(0,1,2)*length(giorgi.short),
         ##                           which(giorgi.short == 'NEU') + c(0,1,2)*
         ##                           length(giorgi.short)),], 0, 0, 0, 0, 0)
         #  # Hadley cell strength (HC)
         #  # do -s fldmax -sellonlatbox,-180,180,0,30 -zonmean -sellevel,50000 -selvar,stream in out
         #  # zonal mean
         #   if (!no_stream) {
         #     tmphc <- aggregate(echam$ensmean[which(echam$names=='stream'),],
         #                      list(echam$latstream),mean)
         #     # fieldmax
         #     echHC <- apply(tmphc[,2:ncol(tmphc)],2,max)
         #   }
         #   # midlatitude circulation (Z300)
         #   # cdo -s fldmean -sellonlatbox,-180,180,30,60 -sellevel,3000 -selvar,geopoth in out
         #   tmpweights <- cos(echam$latgph300/180*pi)
         #   tmpweights <- tmpweights/sum(tmpweights)
         #   echZ300 <- apply(echam$ensmean[which(echam$names=='gph300'),]*tmpweights,2,sum) 
         #   #z300 <- apply(echam$ensmean[which(echam$names=='gph300'),],2,mean)
         #   # subtropical jet (SJ)
         #   # cdo -s fldmax -sellonlatbox,-180,180,0,50 -zonmean -sellevel,20000 -selvar,u in out
         #   tmpu200 <- aggregate(echam$ensmean[which(echam$names=='u200'),],
         #                        list(echam$latu200),mean)
         #   # fieldmax
         #   echSJ <- apply(tmpu200[,2:ncol(tmpu200)],2,max)
         #   # stratospheric polar vortex (Z100)
         #   # cdo -s sub -fldmean -sellonlatbox,-180,180,75,90 ${STOR_DIR}/gph100.$f -fldmean -sellonlatbox,-180,180,40,55 in out
         #   tmpweights <- cos(echam$latgph100/180*pi)
         #   tmpweights <- tmpweights/sum(tmpweights)
         #   tmpz100 <- apply(echam$ensmean[which(echam$names=='gph100'),]*tmpweights,2,sum) 
         #   tmpweights <- cos(echam$latgph100_2/180*pi)
         #   tmpweights <- tmpweights/sum(tmpweights)
         #   tmpz100_2 <- apply(echam$ensmean[which(echam$names=='gph100_2'),]*tmpweights,2,sum) 
         #   echZ100 <- tmpz100 - tmpz100_2
         #   #echZ100 <- rep(0,length(tmpz100)) # set to zero until working
         #   # Pacific walker circulation (PWC)
         #   # cdo -s sellevel,50000 -selvar,omega in out
         #   # cdo -s sub -fldmean -sellonlatbox,-180,-100,-10,10 ${STOR_DIR}/omega.$f -fldmean -sellonlatbox,100,150,-10,10 in out
         #   tmpweights <- cos(echam$latomega/180*pi)
         #   tmpweights <- tmpweights/sum(tmpweights)
         #   tmpomega <- apply(echam$ensmean[which(echam$names=='omega'),]*tmpweights,2,sum) 
         #   tmpweights <- cos(echam$latomega_2/180*pi)
         #   tmpweights <- tmpweights/sum(tmpweights)
         #   tmpomega_2 <- apply(echam$ensmean[which(echam$names=='omega_2'),]*tmpweights,2,sum) 
         #   echPWC <- tmpomega - tmpomega_2
         #   
         #  #for (n in 5:length(indices)) Hind[n,which(echam$names == indices[n])] <- 1
         eind <- list(ensmean=Hind%*%echam$ensmean, names=indices)
         #  eind <- list(ensmean=Hind%*%echam_noindex$ensmean, names=indices)
         #  #eind$data <- array(Hind %*% array(echam$data, c(nrow(echam$data), length(echam$data)/nrow(echam$data))), c(nrow(Hind), dim(echam$data)[2:3]))
         #  eind$ensmean[which(indices=='HC'),] <- echHC
         #  eind$ensmean[which(indices=='SJ'),] <- echSJ
         #  eind$ensmean[which(indices=='Z100'),] <- echZ100
         #  eind$ensmean[which(indices=='Z300'),] <- echZ300
         #  eind$ensmean[which(indices=='PWC'),] <- echPWC
         eind$ensmean[which(indices=='DIMI'),] <- ech_ind$ensmean[which(ech_ind$names=='DIMI'),]
         eind$ensmean[which(indices=='HC'),] <- ech_ind$ensmean[which(ech_ind$names=='HC'),]
         eind$ensmean[which(indices=='SJ'),] <- ech_ind$ensmean[which(ech_ind$names=='SJ'),]
         eind$ensmean[which(indices=='Z100'),] <- ech_ind$ensmean[which(ech_ind$names=='z100'),]
         eind$ensmean[which(indices=='Z300'),] <- ech_ind$ensmean[which(ech_ind$names=='z300'),]
         eind$ensmean[which(indices=='PWC'),] <- ech_ind$ensmean[which(ech_ind$names=='PWC'),]
         # data for ensemble member
         eind$data <- array(Hind %*% array(echam$data, c(nrow(echam$data), 
                                                         length(echam$data)/nrow(echam$data))), c(nrow(Hind), 
                                                                                                  dim(echam$data)[2:3]))
         #  eind$data <- array(Hind %*% array(echam_noindex$data, c(nrow(echam_noindex$data), 
         #                  length(echam_noindex$data)/nrow(echam_noindex$data))), c(nrow(Hind), 
         #                  dim(echam_noindex$data)[2:3]))
         #  for (i in 1:nmem) {
         #     tmphc <- aggregate(echam$data[which(echam$names=='stream'),,i],
         #                        list(echam$latstream),mean)
         #     echHC <- apply(tmphc[,2:ncol(tmphc)],2,max)
         #     tmpweights <- cos(echam$latgph300/180*pi)
         #     tmpweights <- tmpweights/sum(tmpweights)
         #     echZ300 <- apply(echam$data[which(echam$names=='gph300'),,i]*tmpweights,2,sum) 
         #     tmpu200 <- aggregate(echam$data[which(echam$names=='u200'),,i],
         #                          list(echam$latu200),mean)
         #     echSJ <- apply(tmpu200[,2:ncol(tmpu200)],2,max)
         #     tmpweights <- cos(echam$latgph100/180*pi)
         #     tmpweights <- tmpweights/sum(tmpweights)
         #     tmpz100 <- apply(echam$data[which(echam$names=='gph100'),,i]*tmpweights,2,sum) 
         #     tmpweights <- cos(echam$latgph100_2/180*pi)
         #     tmpweights <- tmpweights/sum(tmpweights)
         #     tmpz100_2 <- apply(echam$data[which(echam$names=='gph100_2'),,i]*tmpweights,2,sum) 
         #     echz100 <- tmpz100 - tmpz100_2
         #     #echZ100 <- rep(0,length(tmpz100)) # set to zero until working
         #     tmpweights <- cos(echam$latomega/180*pi)
         #     tmpweights <- tmpweights/sum(tmpweights)
         #     tmpomega <- apply(echam$data[which(echam$names=='omega'),,i]*tmpweights,2,sum) 
         #     tmpweights <- cos(echam$latomega_2/180*pi)
         #     tmpweights <- tmpweights/sum(tmpweights)
         #     tmpomega_2 <- apply(echam$data[which(echam$names=='omega_2'),,i]*tmpweights,2,sum) 
         #     echPWC <- tmpomega - tmpomega_2
         #     eind$data[which(indices=='HC'),,i] <- echHC
         #     eind$data[which(indices=='SJ'),,i] <- echSJ
         #     eind$data[which(indices=='Z100'),,i] <- echZ100
         #     eind$data[which(indices=='Z300'),,i] <- echZ300
         #     eind$data[which(indices=='PWC'),,i] <- echPWC
         eind$data[which(indices=='DIMI'),,] <- ech_ind$data[which(ech_ind$names=='DIMI'),,]
         eind$data[which(indices=='HC'),,] <- ech_ind$data[which(ech_ind$names=='HC'),,]
         eind$data[which(indices=='SJ'),,] <- ech_ind$data[which(ech_ind$names=='SJ'),,]
         eind$data[which(indices=='Z100'),,] <- ech_ind$data[which(ech_ind$names=='z100'),,]
         eind$data[which(indices=='Z300'),,] <- ech_ind$data[which(ech_ind$names=='z300'),,]
         eind$data[which(indices=='PWC'),,] <- ech_ind$data[which(ech_ind$names=='PWC'),,]
         #   } 
         #   
         #   tmphc <- aggregate(analysis$ensmean[which(analysis$names=='stream'),],
         #                      list(analysis$latstream),mean)
         #   # fieldmax
         #   anaHC <- apply(tmphc[,2:ncol(tmphc)],2,max)
         #   # midlatitude circulation (Z300)
         #   tmpweights <- cos(analysis$latgph300/180*pi)
         #   tmpweights <- tmpweights/sum(tmpweights)
         #   anaZ300 <- apply(analysis$ensmean[which(analysis$names=='gph300'),]*tmpweights,2,sum) 
         #   # subtropical jet (SJ)
         #   tmpu200 <- aggregate(analysis$ensmean[which(analysis$names=='u200'),],
         #                        list(analysis$latu200),mean)
         #   # fieldmax
         #   anaSJ <- apply(tmpu200[,2:ncol(tmpu200)],2,max)
         #   # stratospheric polar vortex (Z100)
         #   tmpweights <- cos(analysis$latgph100/180*pi)
         #   tmpweights <- tmpweights/sum(tmpweights)
         #   tmpz100 <- apply(analysis$ensmean[which(analysis$names=='gph100'),]*tmpweights,2,sum) 
         #   tmpweights <- cos(analysis$latgph100_2/180*pi)
         #   tmpweights <- tmpweights/sum(tmpweights)
         #   tmpz100_2 <- apply(analysis$ensmean[which(analysis$names=='gph100_2'),]*tmpweights,
         #                      2,sum) 
         #   anaZ100 <- tmpz100 - tmpz100_2
         #   #anaZ100 <- rep(0,length(tmpz100)) # set to zero until working
         #   # Pacific walker circulation (PWC)
         #   tmpweights <- cos(analysis$latomega/180*pi)
         #   tmpweights <- tmpweights/sum(tmpweights)
         #   tmpweights <- cos(analysis$latomega_2/180*pi)
         #   tmpweights <- tmpweights/sum(tmpweights)
         #   tmpomega_2 <- apply(analysis$ensmean[which(analysis$names=='omega_2'),]*tmpweights,
         #                       2,sum) 
         #   anaPWC <- tmpomega - tmpomega_2
         aind <- list(ensmean=Hind%*%analysis$ensmean, names=indices)
         #   aind <- list(ensmean=Hind%*%analysis_noindex$ensmean, names=indices)
         #   aind$ensmean[which(indices=='HC'),] <- anaHC
         #   aind$ensmean[which(indices=='SJ'),] <- anaSJ
         #   aind$ensmean[which(indices=='Z100'),] <- anaZ100
         #   aind$ensmean[which(indices=='Z300'),] <- anaZ300
         #   aind$ensmean[which(indices=='PWC'),] <- anaPWC
         aind$ensmean[which(indices=='DIMI'),] <- ana_ind$ensmean[which(ana_ind$names=='DIMI'),]
         aind$ensmean[which(indices=='HC'),] <- ana_ind$ensmean[which(ana_ind$names=='HC'),]
         aind$ensmean[which(indices=='SJ'),] <- ana_ind$ensmean[which(ana_ind$names=='SJ'),]
         aind$ensmean[which(indices=='Z100'),] <- ana_ind$ensmean[which(ana_ind$names=='z100'),]
         aind$ensmean[which(indices=='Z300'),] <- ana_ind$ensmean[which(ana_ind$names=='z300'),]
         aind$ensmean[which(indices=='PWC'),] <- ana_ind$ensmean[which(ana_ind$names=='PWC'),]  
         #   #aind <- list(ensmean=Hind%*%analysis$ensmean, names=indices)
         #   #aind$data <- array(Hind %*% array(analysis$data, c(nrow(analysis$data), length(analysis$data)/nrow(analysis$data))), c(nrow(Hind), dim(analysis$data)[2:3]))
         #   
         #   # data for ensemble member
         aind$data <- array(Hind %*% array(analysis$data, c(nrow(analysis$data), 
                                                            length(analysis$data)/nrow(analysis$data))), c(nrow(Hind), 
                                                                                                           dim(analysis$data)[2:3]))
         #  aind$data <- array(Hind %*% array(analysis_noindex$data, c(nrow(analysis_noindex$data), 
         #                 length(analysis_noindex$data)/nrow(analysis_noindex$data))), c(nrow(Hind), 
         #                 dim(analysis_noindex$data)[2:3]))
         #   for (i in 1:nmem) {
         #     tmphc <- aggregate(analysis$data[which(analysis$names=='stream'),,i],
         #                        list(analysis$latstream),mean)
         #     anaHC <- apply(tmphc[,2:ncol(tmphc)],2,max)
         #     tmpweights <- cos(analysis$latgph300/180*pi)
         #     tmpweights <- tmpweights/sum(tmpweights)
         #     anaZ300 <- apply(analysis$data[which(analysis$names=='gph300'),,i]*tmpweights,2,sum) 
         #     tmpu200 <- aggregate(analysis$data[which(analysis$names=='u200'),,i],
         #                          list(analysis$latu200),mean)
         #     anaSJ <- apply(tmpu200[,2:ncol(tmpu200)],2,max)
         #     tmpweights <- cos(analysis$latgph100/180*pi)
         #     tmpweights <- tmpweights/sum(tmpweights)
         #     tmpz100 <- apply(analysis$data[which(analysis$names=='gph100'),,i]*tmpweights,2,sum) 
         #     tmpweights <- cos(analysis$latgph100_2/180*pi)
         #     tmpweights <- tmpweights/sum(tmpweights)
         #     tmpz100_2 <- apply(analysis$data[which(analysis$names=='gph100_2'),,i]*tmpweights,
         #                        2,sum) 
         #     anaz100 <- tmpz100 - tmpz100_2
         #     #anaZ100 <- rep(0,length(tmpz100)) # set to zero until working
         #     tmpweights <- cos(analysis$latomega/180*pi)
         #     tmpweights <- tmpweights/sum(tmpweights)
         #     tmpomega <- apply(analysis$data[which(analysis$names=='omega'),,i]*tmpweights,2,sum) 
         #     tmpweights <- cos(analysis$latomega_2/180*pi)
         #     tmpweights <- tmpweights/sum(tmpweights)
         #     tmpomega_2 <- apply(analysis$data[which(analysis$names=='omega_2'),,i]*tmpweights,
         #                         2,sum) 
         #     anaPWC <- tmpomega - tmpomega_2
         #     aind$data[which(indices=='HC'),,i] <- anaHC
         #     aind$data[which(indices=='SJ'),,i] <- anaSJ
         #     aind$data[which(indices=='Z100'),,i] <- anaZ100
         #     aind$data[which(indices=='Z300'),,i] <- anaZ300
         #     aind$data[which(indices=='PWC'),,i] <- anaPWC
         aind$data[which(indices=='DIMI'),,] <- ana_ind$data[which(ana_ind$names=='DIMI'),,]
         aind$data[which(indices=='HC'),,] <- ana_ind$data[which(ana_ind$names=='HC'),,]
         aind$data[which(indices=='SJ'),,] <- ana_ind$data[which(ana_ind$names=='SJ'),,]
         aind$data[which(indices=='Z100'),,] <- ana_ind$data[which(ana_ind$names=='z100'),,]
         aind$data[which(indices=='Z300'),,] <- ana_ind$data[which(ana_ind$names=='z300'),,]
         aind$data[which(indices=='PWC'),,] <- ana_ind$data[which(ana_ind$names=='PWC'),,]
         #   }  
         
         # 3. Validation data set indices
         #  H.giorgi2 <- compute_giorgi_H(giorgi, validate) #, numvar=3) # 3 vars temp, precip, slp
         #  Hind2 <- rbind(H.giorgi2[c(which(giorgi.short == 'NH')+ c(0,1,2)*length(giorgi.short),
         #                             which(giorgi.short == 'NEU') + c(0,1,2)*length(giorgi.short)),], 
         #                              0, 0, 0, 0, 0)          
         
         if (vali & ((cru_vali) | (ncep_vali) | (recon_vali))) { 
           # there is already stefans indices in the data, only calc NH and EU
           #  vind <- list(data=Hind2 %*% validate$data[1:length(
           #               which(cruall$names=="temp2"))*3,], names=indices)
           #    validatanona = validate$data
           #    validatanona[is.na(validate$data)]=9e36
           # land sea mask to take care of missing validata in the ocean
           vpos <- which(!is.na(validate$data[,1]))
           Hind2 <- Hind[,vpos]
           validatanona <- validate$data[vpos,]
           vind <- list(data=Hind2 %*% validatanona, names=indices)
           #    vind <- list(data=Hind2 %*% validatanona, names=indices)
           #  vind <- list(data=Hind2 %*% validate$data, names=indices)
           #     vind$data[7,] <- valiall.allts$data[which(valiall$names=="ind_rec_hc"),]
           #     vind$data[8,] <- valiall.allts$data[which(valiall$names=="ind_rec_sj"),]
           #     vind$data[9,] <- valiall.allts$data[which(valiall$names=="ind_rec_z100"),]
           #     vind$data[10,] <- valiall.allts$data[which(valiall$names=="ind_rec_z300"),]
           #     vind$data[11,] <- valiall.allts$data[which(valiall$names=="ind_rec_pwc"),]
           vind$data[7,] <- vali_ind$data[which(vali_ind$names=="ind_rec_hc"),]
           vind$data[8,] <- vali_ind$data[which(vali_ind$names=="ind_rec_sj"),]
           vind$data[9,] <- vali_ind$data[which(vali_ind$names=="ind_rec_z100"),]
           vind$data[10,] <- vali_ind$data[which(vali_ind$names=="ind_rec_z300"),]
           vind$data[11,] <- vali_ind$data[which(vali_ind$names=="ind_rec_pwc"),]
           
         } 
         #   if (ncep_vali) { # there is already stefans indices in the data, only calc NH and EU
         #     #  vind <- list(data=Hind2 %*% validate$data[1:length(which(ncepall$names=="temp2"))*3,], names=indices)
         #     validatanona = validate$data
         #     validatanona[is.na(validate$data)]=9e36
         #     vind <- list(data=Hind2 %*% validatanona, names=indices)
         #     #  vind <- list(data=Hind2 %*% validate$data, names=indices)
         #     vind$data[7,] <- ncepall.allts$data[which(ncepall$names=="ind_rec_hc"),]
         #     vind$data[8,] <- ncepall.allts$data[which(ncepall$names=="ind_rec_sj"),]
         #     vind$data[9,] <- ncepall.allts$data[which(ncepall$names=="ind_rec_z100"),]
         #     vind$data[10,] <- ncepall.allts$data[which(ncepall$names=="ind_rec_z300"),]
         #     vind$data[11,] <- ncepall.allts$data[which(ncepall$names=="ind_rec_pwc"),]
         #   } 
         #   if (recon_vali) { # no indices available, only calc NH and EU
         #     validate$data[which(is.na(validate$data))] <- 0
         #     vind <- list(data=Hind2 %*% validate$data, names=indices)
         #     vind$data[7:11,] <- 0
         # #  vind$data[1:3,] <- 0 
         # #  vind$data[5,] <- vind$data[5,]/3 # convert seasonal sum to month sum
       }
       
       # compute indices from cut-off distance runs
       if (check_dist) {
         adist.ind <- lapply(ana.dist, function(x) {
           out <- list(ensmean = Hind %*% x$ensmean, names=indices)
           out$data <- array(Hind %*% array(x$data, c(nrow(echam$data), 
                                                      length(echam$data)/nrow(echam$data))), c(nrow(Hind), 
                                                                                               dim(echam$data)[2:3]))
           return(out)})
       }
       
       #valnona = validate$data
       #posnona = which(is.na(valnona))
       #valnona[posnona] = 0
       #vind <- list(data=Hind2 %*% valnona, names=indices)
       
       #for (n in vind$names) {
       #  validate
       #}
       # this calculates the indices 'NHt2m', 'NEUt2m', 'NEUpr', 'NEUslp' as there are no regional averages in the original validation data (CRU TS)
       # 'HC', 'SJ', 'z100', 'DIMI' and simply taken from the validation data set and were calculated based on Stefan's recon pre 1948 and from NCEP reanalysis afterwards
       
       
       # stratospheric polar vortex
       #    cdo -s -sellevel,10000 -selvar,geopoth $f ${STOR_DIR}/gph100.$f
       #    cdo -s sub -fldmean -sellonlatbox,-180,180,75,90 ${STOR_DIR}/gph100.$f -fldmean -sellonlatbox,-180,180,40,55 ${STOR_DIR}/gph100.$f ${STOR_DIR}/z100.$f
       #    ncrename -h -v geopoth,z100 ${STOR_DIR}/z100.$f  
       # Hadley cell strength
       #    cdo -s fldmax -sellonlatbox,-180,180,0,30 -zonmean -sellevel,50000 -selvar,stream $f ${STOR_DIR}/HC.$f
       #    ncrename -h -v stream,HC ${STOR_DIR}/HC.$f 
       # subtropical jet
       #    cdo -s fldmax -sellonlatbox,-180,180,0,50 -zonmean -sellevel,20000 -selvar,u $f ${STOR_DIR}/SJ.$f
       #    ncrename -h -v u,SJ ${STOR_DIR}/SJ.$f  
       
       #      ncrename -h -v geopoth,z100 ${STOR_DIR}/z100.$f
       
       if (vali) {
         if (ind_ECHAM) {
           # compute validation statistics on indices
           ecorr.ind <- corr_fun(eind, vind, seas=2)
           acorr.ind <- corr_fun(aind, vind, seas=2)
           ermse.ind <- rmse_fun(eind, vind, seas=2)
           armse.ind <- rmse_fun(aind, vind, seas=2)
           ebias.ind <- bias_fun(eind, vind, seas=2)
           abias.ind <- bias_fun(aind, vind, seas=2)
           RE.ind <- RE_fun(armse.ind, y=ermse.ind)
         }
       }
       
       # compute validation stats on indices from distance runs
       if (check_dist) {
         if (vali) {
           erank.ind <- rank_fun(eind, vind, seas=2)
           arank.ind <- rank_fun(aind, vind, seas=2)
           acorr.dist.ind <- lapply(adist.ind, corr_fun, y=vind, seas=2)
           armse.dist.ind <- lapply(adist.ind, rmse_fun, y=vind, seas=2)
           abias.dist.ind <- lapply(adist.ind, bias_fun, y=vind, seas=2)
           RE.dist.ind <- lapply(armse.dist.ind, RE_fun, y=ermse.ind)
           arank.dist.ind <- lapply(adist.ind, rank_fun, y=vind, seas=2)
         }
       }
       #}
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       # ---------------------------------------------------------------------------------------------
       # Run EnSRF with varying ensemble size
       # ------------------------------------------------------------------------------
       if (ana.enssize){
         nensn <- seq(5,25,5)
         calibrate$data <- pdata + dist.arr
         R <- apply(array(dist.arr, c(nrow(Hcal), 2, ncol(pdata)/2))**2, c(1,2), sum) / (ncol(pdata)/2 - 1)
         
         enssize <- list()
         for (k in 1:5){
           nens <- nensn[k]
           outlist <- list()
           enssize[[k]] <- list()
           for (anai in 1:10){
             # select nens ensemble members randomly from simulations 1:29
             print(paste('Running the ',anai,'-th analysis for ', nensn[k], ' ensemble members', sep=''))
             nsim <- dim(echam$data)[3]
             simi <- sample(1:nsim, nens, replace=FALSE)
             echtmp <- echam
             echtmp$data <- echam$data[,,simi]
             echtmp$ensmean <- apply(echam$data, 1:2, mean, na.rm=T)
             anaens <- EnSRF(echtmp, calibrate, R=R, Hcal=Hcal, weights=d.weights_all)
             
             # compute the output statistics
             outlist[[anai]] <- list()
             outlist[[anai]]$simi <- simi
             outlist[[anai]]$echcor <- corr_fun(echtmp, y=validate, seas=2)
             outlist[[anai]]$anacor <- corr_fun(anaens, y=validate, seas=2)
             outlist[[anai]]$echbias <- bias_fun(echtmp, y=validate, seas=2)
             outlist[[anai]]$anabias <- bias_fun(anaens, y=validate, seas=2)
             outlist[[anai]]$echrmse <- rmse_fun(echtmp, y=validate, seas=2)
             outlist[[anai]]$anarmse <- rmse_fun(anaens, y=validate, seas=2)
             outlist[[anai]]$RE <- RE_fun(outlist[[anai]]$anarmse, outlist[[anai]]$echrmse)
             
             # compute/extract the derived indices for manuscript
             etmpind <- list(ensmean=Hind%*%echtmp$ensmean, names=indices)
             etmpind$data <- array(Hind %*% array(echtmp$data, c(nrow(echtmp$data), length(echtmp$data)/nrow(echtmp$data))), c(nrow(Hind), dim(echtmp$data)[2:3]))
             atmpind <- list(ensmean=Hind%*%anaens$ensmean, names=indices)
             atmpind$data <- array(Hind %*% array(anaens$data, c(nrow(anaens$data), length(anaens$data)/nrow(anaens$data))), c(nrow(Hind), dim(anaens$data)[2:3]))
             vind <- list(data=Hind %*% validate$data, names=indices)
             
             # write RE results to enssize
             ermse.tmp <- rmse_fun(etmpind, vind, seas=2)
             armse.tmp <- rmse_fun(atmpind, vind, seas=2)
             enssize[[k]][[anai]] <- RE_fun(armse.tmp, y=ermse.tmp)      
           }
           save(outlist, armse.tmp, ermse.tmp, file=paste('data/ensemble_size_', nens, '_', paste(timlim, collapse='-'),'.Rdata', sep=''))
         }
         rm(outlist)
       }
       
       
       # ------------------------------------------------------------------------------
       # run EnSRF with NCEP/NCAR and SOCOL
       # ------------------------------------------------------------------------------
       # set up new pseudo-proxies for NNR
       if (NCEP_SOCOL){
         Hnr <- compute_H(proxies.mn, nnr.mn, 1000)
         nr.prox <- nnr.mn
         nr.i <- apply(Hnr, 1, function(x) which(x == 1))
         nr.prox$data <- Hnr %*% nnr.mn$data[,,1]
         nr.prox$ensmean <- NULL
         nr.prox$lon <- nnr.mn$lon[nr.i]
         nr.prox$lat <- nnr.mn$lat[nr.i]
         Hnr <- compute_H(nr.prox, socol.mn)
         
         # alternatively use CRU data for pseudo-proxies
         Hcru <- compute_H(proxies.mn, cru, threshold=1500)
         cru.i <- apply(Hcru, 1, function(x) if (any(x == 1)) which(x == 1) else NA)
         Hnnr <- compute_H(proxies.mn, nnr, threshold=1500)
         nnr.i <- apply(Hnnr, 1, function(x) which(x == 1))
         cruprox <- cru.mn
         cruprox$lon <- cru.mn$lon[cru.i]
         cruprox$lat <- cru.mn$lat[cru.i]
         cruprox$data <- cru.mn$data[cru.i,,] - apply(cru.mn$data[cru.i,,], 1, mean, na.rm=T) + apply(nnr.mn$data[nnr.i,,], 1, mean, na.rm=T)
         Hcru <- compute_H(cruprox, socol.mn)
         
         # compute new localization for NNR
         dist.outer <- array(0, rep(ncol(Hnr), 2))
         for (i in 1:length(socol.mn$lon)) dist.outer[i,1:length(socol.mn$lon)] <- compute_dist(socol.mn$lon, socol.mn$lat, socol.mn$lon[i], socol.mn$lat[i])
         d.weights <- corr_function(dist.outer, L=5000)
         nvar <- nrow(dist.outer) %/% length(socol.mn$lon)
         d.weights2.nnr <- d.weights
         d.weights2.nnr[1:(length(socol.mn$lon)*nvar), 1:(length(socol.mn$lon)*nvar)] <-
           d.weights[rep(1:length(socol.mn$lon), nvar), rep(1:length(socol.mn$lon),nvar)]
         rm(dist.outer)
         
         # ------------------------------------------------------------------------------
         # run EnSRF with SOCOL and NCEP/NCAR reanalysis
         # ------------------------------------------------------------------------------
         ana.mn <- EnSRF(socol.mn, nr.prox, R=rep(0.1, nrow(Hnr)), Hcal=Hnr, weights=d.weights2.nnr)   
         # ------------------------------------------------------------------------------
         # compute rmse, corr and bias
         # ------------------------------------------------------------------------------
         ana.rmse <- rmse_fun(ana.mn, nnr.mn, seas=2)
         ana.bias <- bias_fun(ana.mn, nnr.mn, seas=2)
         ana.corr <- corr_fun(ana.mn, nnr.mn, seas=2)
         soc.rmse <- rmse_fun(socol.mn, nnr.mn, seas=2)
         soc.bias <- bias_fun(socol.mn, nnr.mn, seas=2)
         soc.corr <- corr_fun(socol.mn, nnr.mn, seas=2)
         ana.RE <- RE_fun(ana.rmse, soc.rmse)
         
         rm(echam, calibrate, validate, proxies)
         rm(socol, nnr)
         rm(d.weights, d.weights2.nnr, d.weights2, Hval)
         rm(analysis.mn)
         rm(ana.dist)
       }
       
     } # end calc_vali_stat
    
     print(paste("time for whole:",(proc.time() -ptm)[3]))
    
#############################################################################################################################
#############################################################################################################################
    
        # Cov.ana <- all1
    # rm(all1)
    # meanCov <- apply(Cov.ana,2,function(x){tapply(x,nrobs,mean,na.rm=TRUE)})
    # 
    
    ## Vgl der Proxykoordinaten im Modell und in den Beobachtungen.
    
    # load(paste0(ana_dir,'analysis_',1940,'_2ndgrid.Rdata'))
    # calibratelon<-calibrate$lon[which(calibrate$sour=="prox")]
    # calibratelat<-calibrate$lat[which(calibrate$sour=="prox")]
    # png(filename=paste0(workdir,'../data/analysis/test1.png'),width=800,height = 480)
    # plot(calibrate$lon[which(calibrate$sour=="prox")], calibrate$lat[which(calibrate$sour=="prox")], main='1940 station locations',pch=1,col='red')
    # points(analysis$lon[H.i[,1]],analysis$lat[H.i[,1]],pch=1,col='blue')
    # map("world",interior=F,add=T,ty='l',col='black',xlim=c(-180,180),ylim=c(-90,90))
    # map.axes()
    
    # 
    # legend(
    #   "bottomright", 
    #   pch = 21,
    #   pt.bg=c("red", "blue"), 
    #   legend = c("Calibrate Coord.","Modell Coord"))
    # dev.off()
    


    # library(fields)

    for (m in 1:ncol(h.i)) {
    lon<-analysis$lon[H.i[,m]]
    lat<-analysis$lat[H.i[,m]]
    Cov.ana1<-Cov.ana[H.i[,m],m]
    
    lev <- pretty(Cov.ana1,13)
    br <- length(lev)
    colpal <- two.colors(n=br,start="green", end="purple", middle="red", alpha=1.0)
    datcol <- colpal[as.numeric(cut(Cov.ana1,breaks=br))]
    datcol[1] <- "#000000"
    png(filename=paste0(workdir,'../data/analysis/cov.anal_',m,'.png'),width=800,height = 480)
    
    plot(lon, lat,
                            # cex=point size~precip
         pch=15,                    # point type: 15 is square with possible color filling
         col=datcol,
         xlab="Longitude",ylab="Latitude",main="Cov.ana")                # point fill color
    map("world",interior=F,add=T,ty='l',col='black',xlim=c(-180,180),ylim=c(-90,90))
    legend("bottomleft", inset=0.01, as.character(lev),fill=colpal,bty="o",
           bg="white",box.col="white",box.lwd=0,cex=0.7)
    
    dev.off()
    }

##plot difference: analysis-echam, to see whether there has been great difference
    
      lon<-analysis$lon[H.i[,1]]
      lat<-analysis$lat[H.i[,1]]
      diff1<-diff[H.i[,1],2]
      
      lev <- pretty(diff1,13)
      br <- length(lev)
      colpal <- two.colors(n=br,start="green", end="purple", middle="red", alpha=1.0)
      datcol <- colpal[as.numeric(cut(diff1,breaks=br))]
      png(filename=paste0(workdir,'../data/analysis/abs_diff_1900-1959.png'),width=800,height = 480)
      
      plot(lon, lat,
           # cex=point size~precip
           pch=15,                    # point type: 15 is square with possible color filling
           col=datcol,
           xlab="Longitude",ylab="Latitude",main="Abs(Analysis-Echam)")                # point fill color
      map("world",interior=F,add=T,ty='l',col='black',xlim=c(-180,180),ylim=c(-90,90))
      legend("bottomleft", inset=0.01, as.character(lev),fill=colpal,bty="o",
             bg="white",box.col="white",box.lwd=0,cex=0.7)
      
      dev.off()
    