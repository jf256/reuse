# code can be run from bash:
# Rscript EnSRF_data.R syr eyr
# analysis start from syr+1 if sixmonths statevector==TRUE

rm(list=ls())
# uncomment following lines to run on Joerg's Macbook
#machine="climcal" #"macbook" 
#if (machine=="macbook") {
#  syr=1603     # set to 1603 to process analysis; set to >=1901 or 1902? for cru validation
#  eyr=2004     # max 2004 ?
#  workdir='~/unibe/projects/EnSRF/src/'
#  dataextdir='/Volumes/DATA/climdata/EKF400/'
#  dataintdir=paste0(workdir,'../data/')
#}

# enter syr ane eyr manually


syr=1900
eyr=1990

# read syr and eyr from Rscript parameters entered in bash and 
# if existing overwrite manually entered years 
args <- commandArgs(TRUE)
if (length(args)>0) {
  syr = as.numeric(args[1])
  eyr = as.numeric(args[2])
}
print(paste('period',syr, 'to', eyr))

user <- system("echo $USER",intern=T)
print(paste('User:',user))
if (user=="veronika") {
  # workdir('/scratch/veronika/rerun/r_code')
  workdir ='/scratch3/veronika/reuse/reuse_git/' # where are the scripts from github
} else if (user=="lucaf") {
  workdir='/scratch3/lucaf/reuse/reuse_git/'
} else if (user=="joerg") {
  workdir='/scratch3/joerg/projects/reuse/reuse_git/'
} else if (user == "nevin"){
  workdir = '/scratch3/nevin/reuse_climcal/reuse_git/'
} else{
  stop("Unknown user!")
  
}
dataextdir='/mnt/climstor/giub/EKF400/'
dataintdir=paste0(workdir,'../data/')
setwd(workdir)

source('EnSRF_switches.R')
source('EnSRF_functions.R')

dir.create(paste0("../data/analysis/",expname))
# create logfile
logfn <- paste0(expname,'_',syr,'-',eyr,'_',format(Sys.time(),"%Y%m%d_%H%M"),'.log')
write(c(user),file=paste0('../log/',logfn),append=F)
if (loo) {dir.create(paste0("../data/loo/",expname))}

# it can be useful to save how the switches were set
con <- file(paste0(workdir,"../data/analysis/",expname,"/switches_",format(Sys.time(),"%Y%m%d_%H%M"),".log"))
sink(con, append=TRUE)
sink(con, append=TRUE, type="message")
# This will echo all input and not truncate 150+ character lines...
source('EnSRF_switches.R', echo = TRUE, max.deparse.length=10000)
# Restore output to console
sink() 
sink(type="message")
# And look at the log...
# cat(readLines("test.log"), sep="\n")

source('EnSRF_generate.R')


##########################################################################################
# 0. Loop over years to reduce size of state vector
##########################################################################################
if (sixmonstatevector) {syr2=syr+1} else {syr2=syr}
for (cyr in syr2:eyr) {
  print(expname)
  print(cyr)
  write(cyr,file=paste0('../log/',logfn),append=T)
  if (cyr > 1659&(yuri_temp|yuri_slp|ghcn_temp|isti_instead_ghcn|ghcn_prec)) {
    #Nevin June2018: I added this control sequence instead of the "proxies_only" in the experiment name. Like that: at least one of the instrumental
    #datasets has to be true in the switches such that instrumental = T here. Otherwhise there was an error when only proxies were used.
    instrumental=T        # all instrumental stations 
  } else {
    instrumental=F
  }

  if (TRW|MXD|SCHWEINGR|PAGES|NTREND|TRW_PETRA) {        
    real_proxies=T         # Proxy data experiment (regression NOT H operator) 
  } else {
    real_proxies=F
  }
  if (cyr <= 1853&import_luca) {
    docum=T                 # read documentary based data
  } else {
    docum=F
  }
  
  # if (substring(expname,1,12)=="proxies_only") { ### Luca: this is for experiments where only proxies are used otherwhise there is an error
  #   docum=F
  #   instrumental=F
  #   real_proxies=T
  # }
  # next line not included yet: 
  if (cyr > min(c(syr_cru,syr_twentycr,syr_recon)[c(vali_cru, vali_twentycr, vali_recon)]) & cyr<=max(c(eyr_cru,eyr_twentycr,eyr_recon)[c(vali_cru, vali_twentycr, vali_recon)])) {        # if we don't use reconvali, the eyr here should be changed (Error in valiall : object 'valiall' not found) -> but then instead of the eyr we should use cyr
    vali=T                 # switch off prepplot if no vali data selected
  } else {
    vali=F
  }
  if ((cyr > syr_cru) & (cyr <=eyr_cru) & vali_cru) {
    cru_vali=T             # monthly CRU TS3 temp, precip and HADSLP2 gridded instrumentals (1901-2004)
    #  ind_recon=T            # Stefan's reconstructed indices until 1948 and NCAR reanalysis later added to CRU
  } else {
    cru_vali=F 
    #  ind_recon=F
  }
  if ((cyr > syr_twentycr) & (cyr <=eyr_twentycr) & vali_twentycr) {
    twentycr_vali=T             
  } else {
    twentycr_vali=F 
  }
  
  if ((cyr > syr_recon) & (cyr <=eyr_recon) & vali_recon) {
    recon_vali=T           # seasonal luterbacher, pauling, kuettel recons (1750-1999)
  } else {
    recon_vali=F
  }

  print(paste0("instr:",instrumental, "; proxies:",real_proxies, "; documentary:",docum, 
               "; validation data:",vali,"; CRU:",cru_vali,"; 20cr:",twentycr_vali,"; Recon:",recon_vali))
  write(paste("instr:",instrumental, "; proxies:",real_proxies, "; documentary:",docum, 
              "; validation data:",vali,"; CRU:",cru_vali,"; 20cr:",twentycr_vali,"; Recon:",recon_vali),
        file=paste0('../log/',logfn),append=T)
  asyr <- cyr-35 # "a" for anomaly
  if (asyr < 1601) {asyr = 1601}
  aeyr <- cyr+35
  if (aeyr > 2005) {aeyr = 2005}
  ptm1 <- proc.time()
  
  ##########################################################################################
  # 1. Echam Part
  # 1.1 Loading echam, echam_anom, echam_clim, landcorrected_anom, landcorrected_clim
  # 1.1.2 Load/creat bigger ensemble size 
  # 1.2 Choose which variables want to use from the model
  # 1.3 Calc echam st. dev. for each grid point and month over ens memb. to scale docu data
  # 1.4 Calculate decorrelation distance
  # 1.5 Just leave data for one year in memory and convert to sixmonstatevector data format
  # 1.6 Set up the cutoff distance for localisation
  ##########################################################################################
  
  # 1.1 Loading echam, echam_anom, echam_clim, landcorrected_anom, landcorrected_clim
  if (every2grid) {
    load(paste(dataextdir,"echam/echam_",(cyr-1),"-",cyr,"_2ndgrid.Rdata",sep=""))
  } else {
    load(paste(dataextdir,"echam/echam_",(cyr-1),"-",cyr,".Rdata",sep=""))
  }
  if ((anomaly_assim) & (!no_forc_big_ens)) {  
    # anomalies calculated efficiently with cdo, slow calculation within R has been removed
    yr1 <- cyr-1
    yr2 <- cyr
    yr3 <- yr1
    yr4 <- yr2
    if (cyr < 1637) {yr3 <- 1636}
    if (cyr < 1637) {yr4 <- 1637}
    if (cyr > 1970) {yr3 <- 1969}
    if (cyr > 1970) {yr4 <- 1970}
    if (every2grid) {
      load(paste0(echanompath,'echam_anom_',yr1,'-',yr2,'_2ndgrid.Rdata'))
      load(paste0(echclimpath,'echam_clim_',yr3,'-',yr4,'_2ndgrid.Rdata'))
      if (landcorr) {
        load(paste0(echanompath,'echam_anom_103_',yr1,'-',yr2,'_2ndgrid.Rdata'))
        load(paste0(echclimpath,'echam_clim_103_',yr3,'-',yr4,'_2ndgrid.Rdata'))
      }
    } else {
      load(paste0(echanompath,'echam_anom_',yr1,'-',yr2,'.Rdata'))
      load(paste0(echclimpath,'echam_clim_',yr3,'-',yr4,'.Rdata'))
      if (landcorr) {
        load(paste0(echanompath,'echam_anom_103_',yr1,'-',yr2,'.Rdata'))
        load(paste0(echclimpath,'echam_clim_103_',yr3,'-',yr4,'.Rdata'))
      }
    }
  }
  
  # 1.1.2 Load/creat bigger ensemble size 
  # covarclim: estimate of the background covariance matrix by blending with the climatology
  # no_forc_big_ens: use in the whole assimilation
  # just use limited number of years (n_covar) to make calculation faster
  # state ="static": already in 6monstatevector format and units are correct
  if (covarclim>0) {
    echanomallts <- background_matrix(state ,n_covar, echanomallts)
  } else if (no_forc_big_ens) {
    echam_anom = background_matrix(state ,n_covar, echam_anom)
    load(paste0(dataextdir,"echam_400yr_ensmean_clim/400yr_monthly_ensmean.clim_2ndgrid.Rdata"))
    echam_clim = echam_clim_400yr
    echam_clim$ensmean = echam_clim_400yr$data
    echam_clim$data = array(echam_clim$data, c(dim(echam_clim$data),1))
    print(state)
    print(echam_anom$data[1,1,])
  }
  
  
  # 1.2 Choose which variables want to use from the model
  # just leave temp precip slp in state vector
  if (tps_only) {
    tpspos <- c(which(echam$names=='temp2'), which(echam$names=='precip'),
                which(echam$names=='slp'))
    echam$data <- echam$data[tpspos,,]
    echam$ensmean <- echam$ensmean[tpspos,]
    echam$names <- echam$names[tpspos]
    if (anomaly_assim){
      # tpspos <- c(which(echam_anom$names=='temp2'), which(echam_anom$names=='precip'),
      #             which(echam_anom$names=='slp'))
      tpspos <- which(echam_anom$names=='temp2'|echam_anom$names=='precip'|echam_anom$names=='slp')
      echam_anom$data <- echam_anom$data[tpspos,,]
      echam_anom$ensmean <- echam_anom$ensmean[tpspos,]
      echam_anom$names <- echam_anom$names[tpspos]
      if(state=="static" & ((covarclim>0 | no_forc_big_ens))){
        if(no_forc_big_ens){
          echam_anom$lat<-echam_anom$lat[tpspos]
          echam_anom$lon<-echam_anom$lon[tpspos]
        }else if(covarclim>0){
          echanomallts$lon<-echanomallts$lon[tpspos]
          echanomallts$lat<-echanomallts$lat[tpspos]
        }
      }
      tpspos<-c(which(echam_clim$names=='temp2'), which(echam_clim$names=='precip'),
                which(echam_clim$names=='slp'))
      echam_clim$data <- echam_clim$data[tpspos,,]
      echam_clim$ensmean <- echam_clim$ensmean[tpspos,]
      echam_clim$names <- echam_clim$names[tpspos]
    }
  }
  if (no_stream) {
    # ACHTUNG stream var has ERROR because the 5/9 levels before/after 1880 have a lat dimension
    tpspos <- c(which(echam$names!='stream'))
    echam$data <- echam$data[tpspos,,]
    echam$ensmean <- echam$ensmean[tpspos,]
    echam$names <- echam$names[tpspos]
    if (anomaly_assim){
      tpspos <- c(which(echam_anom$names!='stream'))
      echam_anom$data <- echam_anom$data[tpspos,,]
      echam_anom$ensmean <- echam_anom$ensmean[tpspos,]
      echam_anom$names <- echam_anom$names[tpspos]
      if (covarclim>0 | no_forc_big_ens) {
        if (state == "changing") {
          tpspos <- c(which(echanomallts$names!='stream'))
          echanomallts$data <- echanomallts$data[tpspos,,]
          echanomallts$ensmean <- echanomallts$ensmean[tpspos,]
          echanomallts$names <- echanomallts$names[tpspos]
        }
      if (state == "static" & no_forc_big_ens ) {
        tpspos <- c(which(echam_clim$names!='stream'))
      }
    }
    echam_clim$data <- echam_clim$data[tpspos,,]
    echam_clim$ensmean <- echam_clim$ensmean[tpspos,]
    echam_clim$names <- echam_clim$names[tpspos]
      if (landcorr) {
        landcorrected_anom$data <- landcorrected_anom$data[tpspos,,]
        landcorrected_anom$names <- landcorrected_anom$names[tpspos]
        landcorrected_clim$data <- landcorrected_clim$data[tpspos,,]
        landcorrected_clim$names <- landcorrected_clim$names[tpspos]
      }
    }
  }
  if (fasttest) {
    mulc <- 4 # choose every 4th grid box
    loi <- seq(1:length(echam_anom$lon))
    lai <- seq(1:length(echam_anom$lat))
    di <- seq(1:dim(echam_anom$data)[1])
    ni <- seq(1:length(echam_anom$names))
    loi <- loi[seq(ceiling(mulc/2),length(loi),mulc)]
    lai <- lai[seq(ceiling(mulc/2), length(lai),mulc)]
    di <- di[seq(ceiling(mulc/2), length(di),mulc)]
    ni <- ni[seq(ceiling(mulc/2), length(ni),mulc)]
    echam$lon <- echam$lon[loi]
    echam$lat <- echam$lat[lai]
    echam$data <- echam$data[di,,]
    echam$ensmean <- echam$ensmean[di,]
    echam$names <- echam$names[ni]    
    if (anomaly_assim){
      echam_anom$lon <- echam_anom$lon[loi]
      echam_anom$lat <- echam_anom$lat[lai]
      echam_anom$data <- echam_anom$data[di,,]
      echam_anom$ensmean <- echam_anom$ensmean[di,]
      echam_anom$names <- echam_anom$names[ni]
      echam_clim$lon <- echam_clim$lon[loi]
      echam_clim$lat <- echam_clim$lat[lai]
      echam_clim$data <- echam_clim$data[di,,]
      echam_clim$ensmean <- echam_clim$ensmean[di,]
      echam_clim$names <- echam_clim$names[ni]
    }
  }
  if (anomaly_assim){
    # new echam_clim made by veronika have correct units
    # error where echam_anom and echam_clim were generated initially: 
    # no unit correction happened, thus here
    if (!(state =="static" & no_forc_big_ens)) {
        echam_anom$data[echam_anom$names=='precip',,] <- echam_anom$data[echam_anom$names=='precip',,] * 3600 * 24 * 30
        echam_anom$ensmean[echam_anom$names=='precip',] <- echam_anom$ensmean[echam_anom$names=='precip',] * 3600 * 24 * 30
        echam_anom$data[echam_anom$names=='slp',,] <- echam_anom$data[echam_anom$names=='slp',,]/100
        echam_anom$ensmean[echam_anom$names=='slp',] <- echam_anom$ensmean[echam_anom$names=='slp',]/100
        echam.anom <- echam_anom
        if (state=="changing" & covarclim > 0) {
          echanomallts$data[echanomallts$names=='precip',,] <- echanomallts$data[echanomallts$names=='precip',,] * 3600 * 24 * 30
          echanomallts$ensmean[echanomallts$names=='precip',] <- echanomallts$ensmean[echanomallts$names=='precip',] * 3600 * 24 * 30
          echanomallts$data[echanomallts$names=='slp',,] <- echanomallts$data[echanomallts$names=='slp',,]/100
          echanomallts$ensmean[echanomallts$names=='slp',] <- echanomallts$ensmean[echanomallts$names=='slp',]/100
        }
    }else{
      echam.anom<-echam_anom
    }
    if (no_forc_big_ens) {
      echam_clim$data = array(echam_clim$data, c(dim(echam_clim$data),1))  
      echam_clim$data[echam_clim$names=='temp2',,] <- echam_clim$data[echam_clim$names=='temp2',,]-273.15
      echam_clim$ensmean[echam_clim$names=='temp2',] <- echam_clim$ensmean[echam_clim$names=='temp2',]-273.15
      echam_clim$data[echam_clim$names=='precip',,] <- echam_clim$data[echam_clim$names=='precip',,]*3600 * 24 * 30
      echam_clim$ensmean[echam_clim$names=='precip',] <- echam_clim$ensmean[echam_clim$names=='precip',]*3600 * 24 * 30
      echam_clim$data[echam_clim$names=='slp',,] <- echam_clim$data[echam_clim$names=='slp',,]/100
      echam_clim$ensmean[echam_clim$names=='slp',] <- echam_clim$ensmean[echam_clim$names=='slp',]/100
    }
    echam.clim <- echam_clim
    # echam_clim_mon_ensmean (next line) stays 12 months version for instr screening
    if (!no_forc_big_ens) {
      echam_clim_mon_ensmean <- echam_clim$ensmean[,10:21] 
    } else {
      echam_clim_mon_ensmean <- echam_clim$ensmean[,c(10,11,12,1,2,3,4,5,6,7,8,9)]
    }
    if (state=="changing") { # Roni: why I delete it only if state = "changing ?? -> I think because if it is "static" I need it in the background_matrix function
      rm (echam_anom)
    }
    rm(echam_clim)
    if (landcorr) {
      # They (temp2, precip, slp) are corrected for the climatology
      # They (precip, slp) are NOT corrected for the anomaly
      landcorrected_anom$data[landcorrected_anom$names=='precip',] <-
        landcorrected_anom$data[landcorrected_anom$names=='precip',] * 3600 * 24 * 30
      landcorrected_anom$data[landcorrected_anom$names=='slp',] <-
        landcorrected_anom$data[landcorrected_anom$names=='slp',] / 100
      landcorrected.anom <- landcorrected_anom
      landcorrected.clim <- landcorrected_clim
      rm(landcorrected_anom,landcorrected_clim)
    }
  }
  print(proc.time() - ptm1)
  
  
  # 1.3 Calc echam st. dev. for each grid point and month over ens memb. to scale docu data
  echam.sd <- echam # to have non-sixmon echam in docum section -> I would move it to 1.3
  if (sixmonstatevector) {
    echam.sd$data <- apply(echam$data[,10:21,],1:2,sd)
    echam.sd$time <- echam$time[10:21]
  } else {
    echam.sd$data <- apply(echam$data[,13:24,],1:2,sd)
    echam.sd$time <- echam$time[13:24]
  }
  echam.sd$ensmean <- NULL
  print('calc time for standard deviations')
  print(proc.time() - ptm1)
  
  
  # 1.4 Calculate decorrelation distance
  if (calc_decorr_dist) { 
    if (covarclim == 0) {
      d <- compute_dist_2d(echam$lon,echam$lat,echam$lon,echam$lat, region) 
    } else if (covarclim > 0) {
      d <- compute_dist_2d(echanomallts$lon[1:4608],echanomallts$lat[1:4608],echanomallts$lon[1:4608],echanomallts$lat[1:4608], region)
    }
    if (covarclim == 0) {
      ech = echam.anom
    } else if (covarclim > 0) {
      ech = echanomallts
    }
    for (i in unique(ech$names)) {
      print(i)
      if (covarclim == 0) {
        tmp <- corr_over_region(echam.anom,-90,90,"global",cor_length_period)
      } else if (covarclim > 0) {
        if (region == "global") {
          tmp = corr_over_region(echanomallts,-90,90,"global",cor_length_period)
        } else if (region == "ENH") {
          tmp = corr_over_region(echanomallts,20,90,"ENH",cor_length_period)
        } else if (region == "ESH") {
          tmp = corr_over_region(echanomallts,-90,-20,"ESH",cor_length_period)
        } else if (region == "tropics") {
          tmp = corr_over_region(echanomallts,-20,20,"tropics",cor_length_period)
        }
        if (region == "lat_band" ) {
          for (k in unique(echanomallts$lat)) {
            if (k == unique(echanomallts$lat)[1]) {
              tmp1 = corr_over_region(echanomallts,k-0.5,k+0.5,"lat_band",cor_length_period)
              tmp1 = array(tmp1, c(dim(tmp1),1))
              tmp = tmp1
            } else {
              tmp1 = corr_over_region(echanomallts,k-0.5,k+0.5,"lat_band",cor_length_period)
              tmp1 = array(tmp1, c(dim(tmp1),1))
              tmp = abind(tmp,tmp1, along=3)
            }
          }
        } else if (region == "lon_band") {
          for (k in unique(echanomallts$lon)) {
            if (k == unique(echanomallts$lon)[1]) {
              tmp1 = corr_over_region(echanomallts,k-0.5,k+0.5,"lon_band",cor_length_period)
              tmp1 = array(tmp1, c(dim(tmp1),1))
              tmp = tmp1
            } else {
              tmp1 = corr_over_region(echanomallts,k-0.5,k+0.5,"lon_band", cor_length_period)
              tmp1 = array(tmp1, c(dim(tmp1),1))
              tmp = abind(tmp,tmp1, along=3)
            }
          }
        }
      }
      # if (region == "lat_band" | region == "lon_band") {
      #   for (k in 1:dim(tmp)[3]) {
      #     if (region == "lat_band") {
      #       png(paste0('../figures/decorr_lat_',cor_length_period,'_',i,'_lat:',unique(echanomallts$lat)[k],'.png'), width = 1024, height = 768)
      #     } else {
      #       png(paste0('../figures/decorr_lon_',cor_length_period,'_',i,'_lon:',unique(echanomallts$lon)[k],'.png'), width = 1024, height = 768)
      #     }
      #     d_dist = d[,,k]
      #     tmp_cor = tmp[,,k]
      #     corens <- cor(t(tmp_cor[,]))
      #     plot(as.vector(d_dist),as.vector(corens),col='blue', xlim=c(0,5000),ylim=c(0,1), pch = 16)
      #     lines(exp(-1/2 * (seq(1:5000)/get(paste0('l_dist_',i)))**2),col='red')
      #     dev.off()
      #   }
      # } else {
      #   corens <- cor(t(tmp[,]))
      #   png(paste0('../figures/decorr_',cor_length_period,'_',region,'_',i,'.png'), width = 1024, height = 768)
      #   plot(as.vector(d),as.vector(corens),col='#4c8bff01',
      #        xlim=c(0,5000),ylim=c(0,1))
      #   lines(exp(-1/2 * (seq(1:5000)/get(paste0('l_dist_',i)))**2),col='red')
      #   dev.off()
      # }
      # make a plot as in Ingebly 2001
      # select unique lats from matrix d and corens
      ind0 = seq(1,ncol(d),96)
      d96 =d[,ind0]
      corens <- cor(t(tmp[,]))
      cor96 =corens[,ind0]
      
      l = 1
      coruj = matrix(NA,nrow=49,ncol=48)
      for (g in ind0){
      vmi = matrix(NA,nrow=49,ncol=96)
      k=1
      for (vv in g:(g+95)) {
        if (vv+48 > 4608) {
          dif = (vv+48) - 4608
          vmi[,k] = corens[c(vv:((vv+48)-dif),1:dif),vv]
        } else {
          vmi[,k] = corens[vv:(vv+48),vv]
        }
          k=k+1
      }
      coruj[,l] = rowMeans(vmi)
      l=l+1
    }
      
      # select for the growing longitudes (from 0 till 180) the data
      # first distance should be 0, first cor = 1
      duj = matrix(NA,nrow=49,ncol=48)
      #coruj = matrix(NA,nrow=49,ncol=48)
      for(v in 1:48){
        duj[,v] = d96[ind0[v]:(ind0[v]+48),v]
        #coruj[,v] = cor96[ind0[v]:(ind0[v]+48),v]
      }
      
      # need the distance along a latitude and not on a sphere that in in d
      for (h in 2:48) {
      duj[h+1,] = duj[2,]*h
      }
      
      # interpolating the cor values to fix distances
      fixdist <- seq(100,3000,100) # define the fix distances
      fixcor <- matrix(NA,nrow=48,ncol=length(fixdist)+1)
      fixcor[,1] <- 1
      for (ilat in 1:48) {
        for (j in 1:length(fixdist)) {
          if (max(duj[,ilat] >= fixdist[j])) {
            # find the minimum positive distance
            difference <- fixdist[j] - duj[,ilat]
            difference[difference<0] <- 10e+10
            k <- which.min(difference)
            # interpolate between the closest value to j and the one after that
            x <- duj[k:(k+1),ilat]
            y <- coruj[k:(k+1),ilat]
            lmodel <- lm(y~x)
            fixcor[ilat,j+1] <- lmodel$coefficients[1] + lmodel$coefficients[2]*fixdist[j]
          }
        }
      }
      # end interpolation
      lats = unique(echanomallts$lat)[48:1]
      corplot = t(fixcor)[2:31,48:1]
      colfunc <- colorRampPalette(c("dodgerblue4","dodgerblue1","white","brown1","brown4"))
      # png(paste0('../figures/contourplot_',i,'.png'), width = 768, height = 1024)
      pdf(paste0('../figures/contourplot_',i,'.pdf'),paper = "a4") # saving as pdf displayed with white interior grid
      # use only every second lat
      # filled.contour(x=fixdist,y=lats[seq(2,length(lats),2)],z=corplot[,seq(2,ncol(corplot),2)], col= colfunc(10),levels=seq(-1,1,0.2),
      #                plot.axes = { contour(x=fixdist,y=lats[seq(2,length(lats),2)],z=corplot[,seq(2,ncol(corplot),2)], nlevels = 5, 
      #                                      drawlabels = TRUE, axes = FALSE, frame.plot = FALSE, add = TRUE);
      #                                      axis(1); axis(2,at=seq(-90,90,10)) },
      #                key.axes = axis(4, seq(-1,1,0.2)) )
      # use all the lat
      filled.contour(x=fixdist,y=lats,z=corplot, col= colfunc(10),levels=seq(-1,1,0.2),
                     plot.axes = { contour(x=fixdist,y=lats,z=corplot, nlevels = 5, labcex = 1, drawlabels = TRUE, axes = FALSE, frame.plot = FALSE, add = TRUE);
                                   axis(1); axis(2,at=seq(-90,90,10)) },
                     key.axes = axis(4, seq(-1,1,0.2)),
                     plot.title = title(main =paste(i), xlab = "Distance", ylab = "Latitude")
                     )
      
      dev.off()
    }
  }
  
  
  # 1.5 Just leave data for one year in memory and convert to sixmonstatevector data format
  if (sixmonstatevector) {
    # change array to have 6 months in state vector for winter and summer
    # first winter starts in oct of syr
    # 6 mon stat vectors for oct-mar and apr and sep
    
    echam <-convert_to_2_seasons(echam,source="echam")
    if (!no_forc_big_ens & covarclim== 0) {
      echam.anom <- convert_to_2_seasons(echam.anom,source="echam")
    } else if (state != "static" & no_forc_big_ens) {
      echam.anom <- convert_to_2_seasons(echam.anom,source="echam")
    } else if (covarclim>0) {
      echam.anom <- convert_to_2_seasons(echam.anom,source="echam")
    }
    
    if (state=="changing" & (covarclim > 0)) {
      echanomallts = convert_to_2_seasons(echanomallts,source="echam")
    }
    
    
    if (landcorr) {
      land41 <- array(landcorrected.anom$data,c(dim(landcorrected.anom$data)[1]*
                                                  dim(landcorrected.anom$data)[2], 1))
      land42 <- land41[((9*dim(landcorrected.anom$data)[1]+1):(dim(land41)[1]-
                                                                 (3*dim(landcorrected.anom$data)[1]))),]
      land42 = array(land42, c(length(land42),1))
      landcorrected.anom$data <- array(land42,c(dim(land42)[1]/
                                                  (((dim(landcorrected.anom$data)[2]/12)-1)*2),
                                                (((dim(landcorrected.anom$data)[2]/12)-1)*2),dim(land42)[2]))
      landcorrected.anom$time <- c(cyr,cyr+0.5)
    }
    
    
    if (no_forc_big_ens) {
      # Special because climatology is only 1 year long -> reorder months
      tmp61 <- array(echam.clim$data,c(dim(echam.clim$data)[1]*dim(echam.clim$data)[2],
                                       dim(echam.clim$data)[3]))
      tmp62 <- tmp61[((9*dim(echam.clim$data)[1]+1):(dim(tmp61)[1])),] # cut from oct till dec
      tmp63 <- tmp61[(1:(9*dim(echam.clim$data)[1])),] # from jan till sept
      tmp64 <- c(tmp62,tmp63)
      echam.clim$data <- array(tmp64,c(length(tmp64)/(((dim(echam.clim$data)[2]/12))*2),
                                       (((dim(echam.clim$data)[2]/12))*2),dim(echam.clim$data)[3])) 
      tmp71 <- array(echam.clim$ensmean,c(dim(echam.clim$ensmean)[1]*dim(echam.clim$ensmean)[2]))
      tmp72 <- tmp71[((9*dim(echam.clim$ensmean)[1]+1):(dim(tmp71)[1]))] # cut oct syr to sep eyr
      tmp73 <- tmp71[(1:(9*dim(echam.clim$ensmean)[1]))]
      tmp74 <- c(tmp72,tmp73)
      echam.clim$ensmean <- array(tmp74,c(length(tmp74)/(((dim(echam.clim$ensmean)[2]/12))*2),
                                          (((dim(echam.clim$ensmean)[2]/12))*2)))  
      echam.clim$time <-  c(cyr,cyr+0.5) # what shall be the year???
      rm(tmp71);rm(tmp72);rm(tmp73);rm(tmp74)
      if(state=="static"){
        echam.anom$time<-c(cyr,cyr+0.5)
      }
    } else {
      echam.clim<-convert_to_2_seasons(echam.clim,source="echam")
    }    
    
    echam.clim$names <- rep(echam.clim$names,6)
    
    if (landcorr) {
      land61 <- array(landcorrected.clim$data,c(dim(landcorrected.clim$data)[1]*dim(landcorrected.clim$data)[2],1))
      land62 <- land61[((9*dim(landcorrected.clim$data)[1]+1):(dim(land61)[1]-(3*dim(landcorrected.clim$data)[1]))),]
      land62 = array(land62, c(length(land62),1))
      landcorrected.clim$data <- array(land62,c(dim(land62)[1]/(((dim(landcorrected.clim$data)[2]/12)-1)*2),
                                                (((dim(landcorrected.clim$data)[2]/12)-1)*2),dim(land62)[2]))
      landcorrected.clim$time <- c(cyr,cyr+0.5)
      landcorrected.clim$names <- rep(landcorrected.clim$names,6)
      rm(land41);rm(land42);rm(land61);rm(land62)
    }
  }
  
  # rename echam.amon to echam if anomaly_assim==T
  if (anomaly_assim){
    echam=echam.anom
    rm(echam.anom)
  } 
  
  # convert lon, lat, names sixmonstatevector data format
  # ACHTUNG: only 11 vars withOUT stream function included so far
  if (sixmonstatevector) {
    if(!(state=="static" & no_forc_big_ens)){
      one_var_dim<-length(echam$lon)
      numvar <- length(unique(echam$names)) 
      echam$lon <-  c(rep(echam$lon, (numvar*6)))
      echam$lat <- c(rep(echam$lat, (numvar*6)))
      echam$names <- c(rep(echam$names, 6))
      echam <- echam[c('data', 'ensmean', 'lon', 'lat', 'height', 'lsm.i', 'time', 'names')]
    }else{
      one_var_dim<-length(echam$lon)/6
      numvar <- length(unique(echam$names)) 
    }
    if (state =="changing" & (covarclim > 0)) {
      echanomallts$lon <-  c(rep(echanomallts$lon, (numvar*6)))
      echanomallts$lat <- c(rep(echanomallts$lat, (numvar*6)))
      echanomallts$names <- c(rep(echanomallts$names, 6))
      echanomallts <- echanomallts[c('data', 'ensmean', 'lon', 'lat', 'height', 'lsm.i', 'time', 'names')]
    }
    if (landcorr) {
      landcorrected.anom$lon <-  c(rep(landcorrected.anom$lon, (numvar*6)))
      landcorrected.anom$lat <- c(rep(landcorrected.anom$lat, (numvar*6)))
      landcorrected.anom$names <- c(rep(landcorrected.anom$names, 6))
      landcorrected.anom <- landcorrected.anom[c('data', 'ensmean', 'lon', 'lat', 'height', 'lsm.i', 'time', 'names')]
    }
  }
  
  
  # 1.6 Set up the cutoff distance for localisation
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
  # lvec['t850'] <- l_dist_t850 # replaced by t500
  lvec['t500'] <- l_dist_t500
  
  print('calc time for loading data')
  print(proc.time() - ptm1)
  
  if (covarclim > 0) {
    if (PHclim_loc) {
      PHclim_lvec = lvec*PHclim_lvec_factor
    }
  }
  
  
  ##########################################################################################
  # 2. Validate Part
  # 2.1 Loading the validation data set
  # 2.2 Choose which variables want to use from the data set
  # 2.3 Cut out the current year
  # 2.4 Convert it to 2 season per year
  # 2.5 Warning message
  # 2.6 Set validate$ensmean equal to validate$data
  ##########################################################################################
  
  # 2.1 Loading the validation data set
  if (vali) {
    validate = list()
    l=0
    ## this part enables experiments with multiple vali data sets. valiname is a variable with all the vali data sets  
    ## set to true. The for loop adds the validata set to a list. 
    valiname = c("cru_vali","recon_vali","twentycr_vali")[c(cru_vali,recon_vali,twentycr_vali)] 
    for (v in valiname) {
      l=l+1
      print(v)
      cru_vali=F
      recon_vali=F
      twentycr_vali=F
      if (v=="cru_vali"){
        cru_vali=T
        print("cru is true")
      } else if (v=="recon_vali"){
        recon_vali=T
        print("recon is true")
      } else if (v=="twentycr_vali"){
        twentycr_vali=T
        print("20cr is true")
      }
      
      if (every2grid) {
        if (recon_vali) {load(paste(dataextdir,"vali_data/recon/recon_allvar_",syr_recon,"-",eyr_recon,"_2ndgrid.Rdata",sep=""))
        } else if (cru_vali) {load(paste(dataextdir,"vali_data/cru/cru_allvar_",syr_cru,"-",eyr_cru,"_2ndgrid.Rdata",sep=""))
        } else if (twentycr_vali){load(paste0(twentycrpath,"twentycr_allvar_",syr_twentycr,"-",eyr_twentycr,"_2ndgrid.Rdata"))}
      } else {
        if (recon_vali) {load(paste(dataextdir,"vali_data/recon/recon_allvar_",syr_recon,"-",eyr_recon,".Rdata",sep=""))
        } else if (cru_vali) {load(paste(dataextdir,"vali_data/cru/cru_allvar_",syr_cru,"-",eyr_cru,".Rdata",sep=""))} 
      }
      # if (ind_recon) {
      #   load(file=paste("../data/indices/indices_recon_",syr,"-",eyr,".Rdata",sep=""))
      # }
      if (cru_vali) {
        valiall <- cruall
        valiall$data <- valiall$data[,,1]
      } else if (recon_vali) {
        valiall <- reconall
        valiall$data <- valiall$data[,,1]
      } else if (twentycr_vali){
        valiall <- twentycr.all
        valiall$data <- valiall$data[,,1]
      } else { vali = F }
    
    # 2.2 Choose which variables want to use from the data set
    if (tps_only) {

      tpspos2 <- c(which(valiall$names=='temp2'), which(valiall$names=='precip'), 
                   which(valiall$names=='slp'))
      valiall$data <- valiall$data[tpspos2,]
      valiall$names <- valiall$names[tpspos2]
    } 
    if (fasttest) {
      mulc <- 4 # choose every 4th grid box
      loi <- seq(1:length(valiall$lon))
      lai <- seq(1:length(valiall$lat))
      di <- seq(1:dim(valiall$data)[1])
      ni <- seq(1:length(valiall$names))
      loi <- loi[seq(ceiling(mulc/2),length(loi),mulc)]
      lai <- lai[seq(ceiling(mulc/2), length(lai),mulc)]
      di <- di[seq(ceiling(mulc/2), length(di),mulc)]
      ni <- ni[seq(ceiling(mulc/2), length(ni),mulc)]
      valiall$lon <- valiall$lon[loi]
      valiall$lat <- valiall$lat[lai]
      valiall$data <- valiall$data[di,,]
      valiall$names <- valiall$names[ni]
    }
    
    # 2.3 Cut out the 24 months around current year 
    valiall.allts=valiall

      ti=which(floor(valiall$time)==(cyr-1) | floor(valiall$time)==cyr) 
      sts=ti[1]
      ets=ti[length(ti)]
      valiall$data=valiall$data[,sts:ets]
      valiall$time=valiall$time[sts:ets]
      
      
      # 2.4 Convert it to 2 season per year
      if (!recon_vali) {
        tmp21 <- array(valiall$data,c(dim(valiall$data)[1]*dim(valiall$data)[2]))
        tmp22 <- tmp21[((9*dim(valiall$data)[1]+1):(dim(tmp21)[1]-(3*dim(valiall$data)[1])))] 
        valiall$data <- array(tmp22,c(dim(tmp22)[1]/(((dim(valiall$data)[2]/12)-1)*2),
                                      (((dim(valiall$data)[2]/12)-1)*2))) # reconvert to 2 seasons per year
        valiall$time <- seq((floor(valiall$time[1])+1),
                            (floor(valiall$time[length(valiall$time)])+0.5),0.5) 
        valiall$names <- rep(valiall$names,6) 
        valiall$lon <- rep(valiall$lon,6) 
        valiall$lat <- rep(valiall$lat,6) 
        rm(tmp21);rm(tmp22)
      }
      if (recon_vali) {
        # only keep winter and summer season from luterbacher and co and remove spring and autumn
        pos <- sort(c(agrep('.042',as.character(valiall$time)), agrep('.542',as.character(valiall$time))))
        valiall$time <- round(valiall$time[pos],digits=1)
        valiall$data <- valiall$data[,pos]      
      }
      
      # 2.5 Warning message
      # if (sum(c(cru_vali,recon_vali))>1) {
      #   print("WARNING: more than 1 validation data set selected!")
      #   write("WARNING: more than 1 validation data set selected!",
      #         file=paste0('../log/',logfn),append=T)
      # } else {
      #   write(paste("cru_vali:",cru_vali,"; recon_vali:",recon_vali),
      #         file=paste0('../log/',logfn),append=T)
      # }
      # 
      # 2.6 Set validate$ensmean equal to validate$data
      validate[[paste(v)]] <- valiall
      validate[[paste(v)]]$ensmean <- validate[[paste(v)]]$data
    }
    if ("cru_vali"%in%valiname) cru_vali=T
    if ("recon_vali"%in%valiname) recon_vali=T
    if ("twentycr_vali"%in%valiname) twentycr_vali=T
  }
  
  
  
  
  ##########################################################################################
  # 3. Real Proxy Part
  # 3.1 Loading proxy data
  # 3.2 Screen the proxy data: more than 5 std. dev. from mean
  # 3.3 Convert to 2 season per year
  # 3.4 Calculate the anomalies
  # 3.5 Create a list named proxies ("combine assimilation data")
  # 3.6 Set real_proxies to FALSE if there is no data
  ##########################################################################################
  
  # 3.1 Loading proxy data
  if (real_proxies){
    if (!generate_PROXIESnew) {
      load(paste0(dataextdir,"assimil_data/rdata_files/real_proxies_",fsyr,"-",feyr,".Rdata")) 
    } else {
      load(paste0("../data/proxies/real_proxies_",fsyr,"-",feyr,".Rdata"))
    }
    
    # 3.2 Screen the proxy data
    if (check_assimdata) {
      # correlation screening already where multiple regression coefficients are calculated
      # thus, screen if value at current time step is more than 5 std. dev. from mean
      # in this case treated as outlier and set to NA
      realprox <- screenstd(realprox,cyr=cyr,source="proxy")
      realprox$data <- t(realprox$data)
      
      # 3.3 Convert to 2 season per year
      # no scaling because regresion takes care of it
      
      
      # the function takes realprox and output is 2 season conversion and realprox.allts(which is needed below)
      listoftwo <-convert_to_2_seasons(realprox,source="proxy")
      realprox.allts <- listoftwo$x.allts
      realprox <- listoftwo$x
      
      # 3.4 Calculate the anomalies
      if (anomaly_assim){
        ti=which(floor(realprox.allts$time)>=(cyr-35) &
                   floor(realprox.allts$time)<=(cyr+35))
        sts=ti[1]
        ets=ti[length(ti)]
        realprox.tmp=realprox.allts
        realprox.tmp$data=realprox.allts$data[,sts:ets]
        realprox.tmp$time=realprox.allts$time[sts:ets]
        realprox.anom <- realprox
        realprox.anom$data <- (realprox$data - (apply(array(realprox.tmp$data,
                                                            c(nrow(realprox.tmp$data), 2, ncol(realprox.tmp$data)/2)), 1:2, mean,na.rm=T)))
        realprox <- realprox.anom
      }
      
      
      # 3.5 Create a list named proxies
      # real trw proxy multiple regression approach
      if (!instrumental) {
        realprox$sour <- rep('prox',length(realprox$lon))
        realprox.allts$sour <- rep('prox',length(realprox$lon))
        if (reduced_proxies) {
          every <- 12
          redpos <-seq(1,length(realprox$lon),every)
          realprox$lon <- realprox$lon[redpos]
          realprox$lat <- realprox$lat[redpos]
          realprox$data <- realprox$data[redpos,]
          realprox$names <- realprox$names[redpos]
          realprox$sour <- realprox$sour[redpos]
          realprox.allts$lon <- realprox.allts$lon[redpos]
          realprox.allts$lat <- realprox.allts$lat[redpos]
          realprox.allts$data <- realprox.allts$data[redpos,]
          realprox.allts$names <- realprox.allts$names[redpos]
          realprox.allts$sour <- realprox.allts$sour[redpos]
        }
        proxies<-list(data=realprox$data, lon=realprox$lon,
                      lat=realprox$lat, names=realprox$names, sour=realprox$sour,
                      height=realprox$elevation, time=realprox$time,
                      mr=realprox$mr, var_residu=realprox$var_residu,
                      numavg=rep(1,length(realprox$lon)))
      }
      
      
      # 3.6 Set real_proxies to FALSE if there is no data
      if (dim(realprox$data)[1]==0) { real_proxies=F }
    }
    
  }
  
  
  
  
  ##########################################################################################
  # 4. Documentary Part
  # 4.1 Loading documentary data
  # 4.2 Calculate the 71 anomaly from cyr-1 October to cyr September
  # 4.3 Scale the documentary data
  # 4.4 Reconvert to 2 seasons per year
  # 4.5 Combine assimilation data into variable named "proxies"
  # 4.6 Combine assimilation data into variable named "proxies"
  ##########################################################################################
  
  if (docum) {
    # 4.1 Loading documentary data
    if (import_luca) {
      load('../assimil_data/t_docu/t_docu_monthly_luca.Rdata')
    }else{
      load(paste0(dataextdir,"assimil_data/rdata_files/t_docu_monthly.Rdata"))
    }
    if (!any(!is.na(t$data))) { docu=F }
    
    if (sixmonstatevector) {
      
      # 4.2 Calculate the 71 anomaly from cyr-1 October to cyr September
      
      if (anomaly_assim){
        # 4.2.1 Calculate the climatology
        
        t.clim<-calculate_climatology(t,cyr,36,35,source="inst")
        t.clim$data<-t(t.clim$data)
        
        t.clim$data = apply(array(t.clim$data, c(nrow(t.clim$data), 12, ncol(t.clim$data)/12)), 1:2, mean,na.rm=T)
        
        # 4.2.3 Cut out the cyr-1 and cyr time window
        t.cyr<-calculate_climatology(t,cyr,1,0,source="inst")
        t.cyr$data<-t(t.cyr$data)
        
        # 4.2.5 Calculate the anomaly for cyr
        t.anom <-t.cyr
        t.anom$data <- t(t.anom$data - t.clim$data)
        doc_t_mon <- t.anom
      } else {
        stop("docum section currently only coded for anomaly_assim!")
      } # end of anomaly assim
      
      # 4.3 Scale the documentary data
      if (scaleprox){
        gpos <- getgridboxnum(doc_t_mon,echam.sd) 
        for (i in 1:length(gpos)) {
          m <- gpos[i]
          if (length(m) > 0) { 
            m = m[echam.sd$names[m]=="temp2"]
            scalefac <- 1/echam.sd$data[m,]
          } else {
            scalefac <- array(NA, c(length(doc_t_mon$lon), length(doc_t_mon$time)))
          }
          if (i == 1) { 
            echam.scale = scalefac
          } else {
            echam.scale = rbind(echam.scale, scalefac)
          }
        }
        # scale only by deviding, already centered to 71yr anom
        doc_t_mon$data = doc_t_mon$data / t(echam.scale)
      }
      
      # 4.4 Reconvert to 2 seasons per year
      
      listoftwo<-convert_to_2_seasons(doc_t_mon,source="doc")
      doc_t_mon <- listoftwo$x
      
    } else {
      stop("docum section currently only coded for sixmonstatevector!")
    } # end of 6monstatevector
    
    
    # 4.5 Combine assimilation data into variable named "proxies"
    if ((!instrumental) & (!real_proxies)) {
      proxies <- doc_t_mon
    } else {
      docall <- doc_t_mon
    }
    
    # 4.6 Combine assimilation data into variable named "proxies"
    if (!instrumental & real_proxies ) {
      docall$sour <- rep('doc',length(docall$lon))
      realprox$sour <- rep('prox',length(realprox$lon))
      tmpnum2 <- rep(1,length(realprox$lon)) #realprox.numavg
      tmpmr <- matrix(NA,nrow=length(docall$lon),ncol=ncol(proxies$mr))
      tmpres <- rep(NA,length(proxies$var_residu))
      tmpelev <- rep(NA,length(proxies$height))
      tmpnum3 <- rep(1,length(docall$lon))
      proxies<-list(data=rbind(realprox$data,docall$data), lon=c(realprox$lon,docall$lon),
                    lat=c(realprox$lat,docall$lat), names=c(realprox$names,docall$names), 
                    sour=c(realprox$sour,docall$sour), 
                    height=c(realprox$height,tmpelev), time=realprox$time,
                    mr=rbind(realprox$mr,tmpmr), var_residu=c(realprox$var_residu,tmpres),
                    numavg=c(tmpnum2,tmpnum3))
    }
  }  # end of docu
  
  
  
  ####################################################################
  # 5. Instrumental Part
  # 5.1 Loading the files
  # 5.2 Qualtiy check of the data
  # 5.3 Calculating the anomalies and transforming them to Oct-Sept time period
  # 5.4 Treat multiple data in the same grid box
  # 5.5 Combaining all type of instrumental data
  # 5.6 Mask other data if there is instrumental data in same grid box at that time
  # 5.7 Use every ??th (see code below) proxy record
  # 5.8 Add data source
  # 5.9 Convert to 2 season per year
  # 5.10 Combine assimilation data into variable named "proxies"
  ####################################################################
  
  if (instrumental){
    # 5.1 Loading the files
    if (yuri_slp) {
      load(paste0(dataextdir,'assimil_data/rdata_files/slp.Rdata')) # monthly slp collection from yuri, histalp included
      inst_slp <- slp
    }
    if (yuri_temp) {
      load(paste0(dataextdir,'assimil_data/rdata_files/t.Rdata')) # monthly temp collection from yuri, histalp included
      inst_t <- t
    }
    if (ghcn_temp) {
      if (isti_instead_ghcn) {
        load(paste0(dataextdir,"assimil_data/isti/isti_1600-2005.Rdata"))
        print("ACHTUNG: ISTI instr. temp. data will be assimilated NOT GHCN")
      } else {
        load(paste0(dataextdir,"assimil_data/ghcn/ghcn_temp_1600-2005.Rdata"))
      }
    }
    if (ghcn_prec) {
      load(paste0(dataextdir,"assimil_data/ghcn/ghcn_precip_1600-2005.Rdata"))
    }
    
    # 5.2 Qualtiy check of the data
    if (check_assimdata) {
      if (ghcn_prec) {
        varlist <- c("inst_t","inst_slp","ghcn","ghcn_precip")
      } else {
        varlist <- c("inst_t","inst_slp","ghcn")
      }
      for (varname in varlist) {
        
        var<-screendistance(echam.sd,varname)
        if (varname=="inst_t") {inst_t$data<- var$data}
        if (varname=="inst_slp") {inst_slp$data<- var$data}
        if (varname=="ghcn") {ghcn$data <- var$data}
        if (ghcn_prec){
          if (varname=="ghcn_precip") {ghcn_precip$data <- var$data}
        }
        
        
        var<-screenstd(get(varname),cyr,source="inst")
        if (varname=="inst_t") {inst_t$data <- var$data}
        if (varname=="inst_slp") {inst_slp$data <- var$data}
        if (varname=="ghcn") {ghcn$data<- var$data}
        if (ghcn_prec){
          if (varname=="ghcn_precip") {ghcn_precip$data <- var$data}
        }
        
        
      }
    }
    
    # 5.3 Calculating the anomalies and transforming them to Oct-Sept time period
    if ((ghcn_temp) & (dim(ghcn$data)[2]>0)) {
      if (anomaly_assim & sixmonstatevector){
        # calculate the climatology
        ghcn.clim<-calculate_climatology(ghcn,cyr,36,35,source="inst")
        ghcn.clim$data<-t(ghcn.clim$data)
        
        ghcn.clim$data = apply(array(ghcn.clim$data, c(nrow(ghcn.clim$data), 12, ncol(ghcn.clim$data)/12)), 1:2, mean,na.rm=T)
        
        # cut out the cyr-1 and cyr time window
        
        ghcn.cyr<-calculate_climatology(ghcn,cyr,1,0,source="inst")
        ghcn.cyr$data<-t(ghcn.cyr$data)
        
        # calculate the anomaly for cyr
        ghcn.anom <-ghcn.cyr
        ghcn.anom$data <- t(ghcn.anom$data - ghcn.clim$data)
        ghcn <- ghcn.anom
      } else {
        stop("instr section currently only coded for anom_assim and sixmonstatevector!")
      }
    }
    if (ghcn_prec) {
      if( dim(ghcn_precip$data)[2]>0) {
        if (anomaly_assim & sixmonstatevector){
          # calculate the climatology
          ghcn_prec.clim=ghcn_prec
          ti=which(floor(ghcn_prec$time)>=(cyr-36) & floor(ghcn_prec$time)<=(cyr+35))
          sts=ti[1]
          ets=ti[length(ti)]
          ghcn_prec.clim$data=t(ghcn_prec$data[sts:ets,])
          ghcn_prec.clim$time=ghcn_prec$time[sts:ets]
          # transform it to Oct-Sept  years
          if (cyr < 1636) {
            y_start = agrep(paste0(1600,'.792'),as.character(ghcn_prec.clim$time))
          } else {
            y_start = agrep(paste0(cyr-36,'.792'),as.character(ghcn_prec.clim$time))
          }
          if (cyr > 1970) {
            y_end =  grep(paste0(2005,'.708'),as.character(ghcn_prec.clim$time))
          } else {
            y_end = grep(paste0(cyr+35,'.708'),as.character(ghcn_prec.clim$time))
          }
          ghcn_prec.clim$data = ghcn_prec.clim$data[,y_start:y_end]
          ghcn_prec.clim$time = ghcn_prec.clim$time[y_start:y_end]
          ghcn_prec.clim$data = apply(array(ghcn_prec.clim$data, c(nrow(ghcn_prec.clim$data), 12, ncol(ghcn_prec.clim$data)/12)), 1:2, mean,na.rm=T)
          
          # cut out the cyr-1 and cyr time window
          ghcn_prec.cyr = ghcn_prec
          ti=which(floor(ghcn_prec$time)>=cyr-1 & floor(ghcn_prec$time)<cyr+1)
          sts=ti[1]
          ets=ti[length(ti)]
          ghcn_prec.cyr$data=t(ghcn_prec$data[sts:ets, ])
          ghcn_prec.cyr$time=ghcn_prec$time[sts:ets]
          # transform it to Oct-Sept  years
          y_start = agrep(paste0(cyr-1,'.792'),as.character(ghcn_prec.cyr$time))
          y_end = grep(paste0(cyr,'.708'),as.character(ghcn_prec.cyr$time))
          ghcn_prec.cyr$data = ghcn._prec.cyr$data[,y_start:y_end]
          ghcn_prec.cyr$time = ghcn_prec.cyr$time[y_start:y_end]
          # calculate the anomaly for cyr
          ghcn_prec.anom <-ghcn_prec.cyr
          ghcn_prec.anom$data <- t(ghcn_prec.anom$data - ghcn_prec.clim$data)
          ghcn <- ghcn_prec.anom
        } else {
          stop("instr section currently only coded for anom_assim and sixmonstatevector!")
        }
      }
    }
    if ((yuri_slp) & (dim(inst_slp$data)[2]>0)) {
      if (anomaly_assim & sixmonstatevector) {
        # calculate the climatology
        
        inst_slp.clim<-calculate_climatology(inst_slp,cyr,36,35,source="inst")
        inst_slp.clim$data<-t(inst_slp.clim$data)
        
        inst_slp.clim$data = apply(array(inst_slp.clim$data, c(nrow(inst_slp.clim$data), 12, ncol(inst_slp.clim$data)/12)), 1:2, mean,na.rm=T)
        
        # cut out the cyr-1 and cyr time window
        
        inst_slp.cyr<-calculate_climatology(inst_slp,cyr,1,0,source="inst")
        inst_slp.cyr$data<-t(inst_slp.cyr$data)
        
        # calculate the anomaly for cyr
        inst_slp.anom <-inst_slp.cyr
        inst_slp.anom$data <- t(inst_slp.anom$data - inst_slp.clim$data)
        inst_slp <-inst_slp.anom
      } else {
        stop("instr section currently only coded for anom_assim and sixmonstatevector!")
      }
    }
    if ((yuri_temp) & (dim(inst_t$data)[2]>0)) {
      if (anomaly_assim & sixmonstatevector){
        
        # calculate the climatology
        
        inst_t.clim<-calculate_climatology(inst_t,cyr,36,35,source="inst")
        inst_t.clim$data<-t(inst_t.clim$data)
        
        inst_t.clim$data = apply(array(inst_t.clim$data, c(nrow(inst_t.clim$data), 12, ncol(inst_t.clim$data)/12)), 1:2, mean,na.rm=T)
        
        # cut out the cyr-1 and cyr time window
        
        inst_t.cyr<-calculate_climatology(inst_t,cyr,1,0,source="inst")
        inst_t.cyr$data<-t(inst_t.cyr$data)
        
        # calculate the anomaly for cyr
        inst_t.anom <- inst_t.cyr
        inst_t.anom$data <- t(inst_t.anom$data - inst_t.clim$data)
        inst_t <-inst_t.anom
      } else {
        stop("instr section currently only coded for anom_assim and sixmonstatevector!")
      }
    }
    if (ghcn_prec) {
      if ((dim(inst_slp$data)[2]==0) & (dim(inst_t$data)[2]==0) & (dim(ghcn$data)[2]==0) &
          (dim(ghcn_precip$data)[2]==0)) {instrumental=F}
    } else {
      if ((dim(inst_slp$data)[2]==0) & (dim(inst_t$data)[2]==0) & (dim(ghcn$data)[2]==0)) {
        instrumental=F}
    }
    
    
    # 5.4 Treat multiple data in the same grid box
    # at the moment only in case of instrumental data
    # for proxy or docum. data assim. without instr., all series get serially assimilated
    # 5.4.1 Combining temp data from different sources
    if ((ghcn_temp) & (yuri_temp)) {
      inst_t<-list(data=cbind(ghcn$data,inst_t$data), lon=c(ghcn$lon,inst_t$lon),
                   lat=c(ghcn$lat,inst_t$lat), names=c(ghcn$names,inst_t$names),
                   height=c(ghcn$height,inst_t$height), time=ghcn$time)
    }
    if ((ghcn_temp) & (!yuri_temp)) {
      inst_t <- ghcn
    }
    if (ghcn_prec) {
      inst_p <- ghcn_precip
    }
    
    # 5.4.2 First proxy per echam grid box
    if (first_prox_per_grid) {
      res <- firstproxres
      newgrid <- echam
      newgrid$data <- NULL
      newgrid$ensmean <- NULL
      newgrid$lon <- echam$lon[seq(1,length(echam$lon),res)]
      newgrid$lat <- echam$lat[seq(1,length(echam$lat),res)]
      
      tmp_t <- inst_t
      tmp_t$data <- t(inst_t$data)
      Hproxy <- compute_H(tmp_t, newgrid, threshold=(700*res))  
      p.i <- as.logical(apply(Hproxy, 1, sum)) & !duplicated(Hproxy, margin=1)
      inst_t$lon <- inst_t$lon[p.i]
      inst_t$lat <- inst_t$lat[p.i]
      inst_t$data <- inst_t$data[,p.i]
      inst_t$names <- inst_t$names[p.i] 
      
      if (ghcn_prec) {
        tmp_precip <- inst_p
        tmp_precip$data <- t(inst_p$data)
        Hproxy <- compute_H(tmp_precip, newgrid, threshold=(700*res))  
        p.i <- as.logical(apply(Hproxy, 1, sum)) & !duplicated(Hproxy, margin=1)
        inst_p$lon <- inst_p$lon[p.i]
        inst_p$lat <- inst_p$lat[p.i]
        inst_p$data <- inst_p$data[,p.i]
        inst_p$names <- inst_p$names[p.i] 
      }
      
      if (yuri_slp) {  
        tmp_slp <- inst_slp
        tmp_slp$data <- t(inst_slp$data)
        Hproxy <- compute_H(tmp_slp, newgrid, threshold=(700*res))  
        p.i <- as.logical(apply(Hproxy, 1, sum)) & !duplicated(Hproxy, margin=1)
        inst_slp$lon <- inst_slp$lon[p.i]
        inst_slp$lat <- inst_slp$lat[p.i]
        inst_slp$data <- inst_slp$data[,p.i]
        inst_slp$names <- inst_slp$names[p.i] 
      }
      
      if (real_proxies) {
        print("DON'T USE FIRST_PROX_PER_GRID IF YOU WANT TO INCLUDE REAL PROXY DATA AND 
              INSTRUMENTALS AT THE SAME TIME")
      }  
      }
    
    # 5.4.3 Average more than one proxy per echam grid box 
    if (avg_prox_per_grid) {
      # average data in same grid box and count number of avg. series as error estimate
      # separate temp, precip, slp
      # makes no sense for realprox: proxies would need to be calibrated before building 
      # regession model
      if (ghcn_prec) {
        varlist <- c("inst_t","inst_p","inst_slp")
      } else {
        varlist <- c("inst_t","inst_slp") 
      }
      for (varname in varlist) {
        stat <- get(varname)
        dlist=NA
        for(i in 1:length(stat$lon)){
          plon <- stat$lon[i]
          plat <- stat$lat[i]
          clon <- echam$lon[!is.na(echam$lon)]
          clat <- echam$lat[!is.na(echam$lat)]
          k=which(abs(clon-plon+0.001)==min(abs(clon-plon+0.001))) 
          # +0.001 to avoid find 2 locations with same distance
          l=which(abs(clat-plat)==min(abs(clat-plat))) 
          m=k[which(match(k,l)>0)]
          if (length(m) > 0) {
            dlist[i]=m[1]
            # it gives warning message because now the echam is in 6monstatevector format the length(m)=66 (6 month * 11 variables)
          } else {
            dlist[i]=NA
          }
        }
        stat.avg=rep(NA,12) 
        stat.numavg=NA
        stat.lon=NA
        stat.lat=NA
        stat.names=NA
        
        for(i in unique(dlist)[(!is.na(unique(dlist)))]){
          no=which(dlist==i)
          mask=as.logical(dlist==i)
          mask[is.na(mask)]=FALSE
          if (length(stat$data[1,mask])>1) {
            stat.avg=rbind(stat.avg,apply(stat$data[,mask],1,mean,na.rm=T))
            stat.numavg=rbind(stat.numavg,(dim(stat$data[,mask])[2]))
            stat.lon=rbind(stat.lon,mean(stat$lon[mask],na.rm=T))
            stat.lat=rbind(stat.lat,mean(stat$lat[mask],na.rm=T))
            stat.names=rbind(stat.names,stat$names[mask][1])
          } else {
            stat.avg=rbind(stat.avg,stat$data[,mask])
            stat.numavg=rbind(stat.numavg,1)
            stat.lon=rbind(stat.lon,stat$lon[mask])
            stat.lat=rbind(stat.lat,stat$lat[mask])
            stat.names=rbind(stat.names,stat$names[mask])
          }
          stat.avg[stat.avg=='NaN']=NA
        }
        stat$lon <- stat.lon[2:length(stat.lon)]
        stat$lat <- stat.lat[2:length(stat.lat)]
        stat$data <- t(stat.avg[2:dim(stat.avg)[1],])
        stat$numavg <- stat.numavg[2:length(stat.numavg)]
        stat$names <- stat.names[2:length(stat.names)]
        stat$ensmean <- NULL
        if (varname == "inst_t") {
          inst_t = stat
        }
        if (varname == "inst_p") {
          inst_p = stat
        }    
        if (varname == "inst_slp") {
          inst_slp = stat
        }
      }
    } 
    
    # 5.5 Combaining all type of instrumental data
    if (ghcn_prec) {
      inst<-list(data=t(cbind(inst_t$data,inst_p$data,inst_slp$data)),
                 lon=c(inst_t$lon,inst_p$lon,inst_slp$lon),
                 lat=c(inst_t$lat,inst_p$lat,inst_slp$lat),
                 numavg=c(inst_t$numavg,inst_p$numavg,inst_slp$numavg),
                 names=c(inst_t$names,inst_p$names,inst_slp$names),
                 height=c(inst_t$height,inst_p$height,inst_slp$height), time=inst_t$time)
    } else { # only t und slp
      inst<-list(data=t(cbind(inst_t$data,inst_slp$data)), lon=c(inst_t$lon,inst_slp$lon),
                 lat=c(inst_t$lat,inst_slp$lat), names=c(inst_t$names,inst_slp$names),
                 numavg=c(inst_t$numavg,inst_slp$numavg),
                 height=c(inst_t$height,inst_slp$height), time=inst_t$time)
    }
    
    # 5.6 Mask other data if there is instrumental data in same grid box at that time
    # 5.6.1 Mask proxy data
    if (instmaskprox) {
      if ((real_proxies) & (instrumental)) {
        # remove inst data that is all NA
        pos <- apply(!is.na(inst$data),1,any)
        inst$data = inst$data[pos,]
        if (length(which(pos))==1) {
          inst$data = t(matrix(inst$data))
        }
        inst$lon <- inst$lon[pos]
        inst$lat <- inst$lat[pos]
        inst$names <- inst$names[pos] 
        for(i in 1:length(realprox$lon)){
          plon <- realprox$lon[i]
          plat <- realprox$lat[i]
          clon <- inst$lon
          clat <- inst$lat
          dist <- compute_dist(plon, plat, clon, clat)
          if (min(dist)<700){
            realprox$data[i,] = NA
            realprox$lon[i] <- NA
            realprox$lat[i] <- NA
            realprox$mr[i,] <- NA
            realprox$var_residu[i] <- NA
          }
        }
        pos <- apply(!is.na(realprox$data),1,any)
        realprox$data = realprox$data[pos,]
        realprox$lon <- realprox$lon[pos]
        realprox$lat <- realprox$lat[pos]
        realprox$mr <- realprox$mr[pos,]
        realprox$var_residu <- realprox$var_residu[pos]
        realprox$names <- realprox$names[pos] 
      }
      # 5.6.2 Mask docu data
      if ((docum) & (instrumental)) {
        ti <- which(floor(inst$time)==cyr)
        sts <- ti[1]
        ets <- ti[length(ti)]
        if (any(!is.na(inst$data[,sts:ets]))) {
          itmp <- apply(inst$data[,sts:ets],1,mean,na.rm=T)
          itmp[inst$names!="temp2"] <- NaN
          clon <- inst$lon[!is.nan(itmp)]
          clat <- inst$lat[!is.nan(itmp)]
          for(i in 1:length(docall$lon)){
            plon <- docall$lon[i]
            plat <- docall$lat[i]
            if ((!is.na(plon)) & (!is.na(plat))) {    
              dist <- compute_dist(plon, plat, clon, clat)
              if (min(dist)<700){
                docall$data[i,] = NA
                docall$lon[i] <- NA
                docall$lat[i] <- NA
              }
            }
          }
        }       
      }
    } # end instmaskprox
    
    # 5.7 Use every ??th (see code below) proxy record
    if (reduced_proxies) {
      every <- 12 
      if (real_proxies) {
        redpos <-seq(1,length(realprox$lon),every)
        realprox$lon <- realprox$lon[redpos]
        realprox$lat <- realprox$lat[redpos]
        realprox$data <- realprox$data[redpos,]
        realprox$names <- realprox$names[redpos]
        realprox.allts$lon <- realprox.allts$lon[redpos]
        realprox.allts$lat <- realprox.allts$lat[redpos]
        realprox.allts$data <- realprox.allts$data[redpos,]
        realprox.allts$names <- realprox.allts$names[redpos]
      }
      if (instrumental) {
        redpos <-seq(1,length(inst$lon),every)
        inst$lon <- inst$lon[redpos]
        inst$lat <- inst$lat[redpos]
        inst$data <- inst$data[redpos,]
        inst$names <- inst$names[redpos]
      }
    }
    
    # 5.8 Add data source/type information to list
    if (instrumental) {inst$sour <- rep('inst',length(inst$lon))}
    if (docum) {
      docall$sour <- rep('doc',length(docall$lon))
      # docall.allts$sour <- rep('doc',length(docall$lon)) # doesn't exists anymore   
    }
    if (real_proxies) {
      realprox$sour <- rep('prox',length(realprox$lon))
      realprox.allts$sour <- rep('prox',length(realprox$lon))
    }
    
    # 5.9 Convert to 2 season per year
    if (sixmonstatevector) { 
      # change instr array that only 2 ts instead of 12 monthly ts but 6 times as 
      # many variables, one for each month
      tmp1 <- array(inst$data,c(dim(inst$data)[1]*dim(inst$data)[2]))
      inst$data <- array(tmp1,c(dim(tmp1)[1]/(dim(inst$data)[2]/6),2))
      rm(tmp1)
      # data already in oct to sep format, hence comment next line
      #tmp2 <- tmp1[((9*dim(inst$data)[1]+1):(length(tmp1)[1]-(3*dim(inst$data)[1])))] 
      #inst$data <- array(tmp2,c(dim(tmp2)[1]/(((dim(inst$data)[2]
      #               /12)-1)*2),(((dim(inst$data)[2]/12)-1)*2))) 
      # reconvert to 2 seasons per year
      inst$time <- c(cyr,cyr+0.5) 
      inst$names <- rep(inst$names,6)
      inst$sour <- rep(inst$sour,6)
      inst$lon <- rep(inst$lon,6)
      inst$lat <- rep(inst$lat,6)
      if (avg_prox_per_grid) {
        inst$numavg <- rep(inst$numavg,6)
      }
    }
    
    # 5.10 Combine assimilation data into variable named "proxies"
    if (!real_proxies) {        
      proxies <- inst
    }
    
    if (real_proxies) {
      tmpmr <- matrix(NA,nrow=length(inst$lon),ncol=ncol(realprox$mr))
      tmpres <- rep(NA,length(realprox$var_residu))
      if (avg_prox_per_grid) {
        tmpnum1 <- inst$numavg
        tmpnum2 <- rep(1,length(realprox$lon))
      } else {
        tmpnum1 <- rep(1,length(inst$lon))
        tmpnum2 <- rep(1,length(realprox$lon)) #realprox.numavg
      }
      if (real_proxies & instrumental) {
        proxies<-list(data=rbind(inst$data,realprox$data), lon=c(inst$lon,realprox$lon), 
                      lat=c(inst$lat,realprox$lat), names=c(inst$names,realprox$names), 
                      sour=c(inst$sour,realprox$sour), 
                      height=c(inst$height,realprox$elevation), time=inst$time,
                      mr=rbind(tmpmr,realprox$mr), var_residu=c(tmpres,realprox$var_residu),
                      numavg=c(tmpnum1,tmpnum2))
      }
    }
    
    if (docum) {
      tmpmr <- matrix(NA,nrow=length(docall$lon),ncol=ncol(proxies$mr))
      tmpres <- rep(NA,length(proxies$var_residu))
      tmpelev <- rep(NA,length(proxies$height))
      tmpnum3 <- rep(1,length(docall$lon))
      proxies<-list(data=rbind(proxies$data,docall$data), lon=c(proxies$lon,docall$lon), 
                    lat=c(proxies$lat,docall$lat), names=c(proxies$names,docall$names), 
                    sour=c(proxies$sour,docall$sour), 
                    height=c(proxies$height,tmpelev), time=proxies$time,
                    mr=rbind(proxies$mr,tmpmr), var_residu=c(realprox$var_residu,tmpres),
                    numavg=c(proxies$numavg,tmpnum3))
    }
    } # end "if (instrumental)"
  
  
  
  
  #########################################################################################
  # 6. All assimilation data
  #########################################################################################
  
  calibrate <- proxies
  print(paste('number of proxies/observations:',dim(calibrate$data)[1]))
  print("calc time preparing proxies")
  print(proc.time() - ptm1)
  
  
  
  
  #########################################################################################
  # 7. Compute H (forward operator)
  #########################################################################################
  if ((real_proxies) & (!instrumental) & (!docum) & (!sixmonstatevector)) {
    Hcal <- Matrix(compute_H_proxy(realprox, echam, realprox$mr, threshold=700), 
                   sparse=T)
  }
  if (sixmonstatevector) {  
    # next 3 lines: solve problem of distance calc with NA in compute_dist function
    etmp <- echam
    etmp$lon[is.na(etmp$lon)] <- 0
    etmp$lat[is.na(etmp$lat)] <- -90
    if (instrumental) {
      Hcal1 <- compute_Hi_Hredux_sixmonstatevector(inst, etmp, threshold=700)
    }
    if (real_proxies) {
      Hcal2 <- compute_Hi_Hredux_proxy(realprox, etmp, realprox$mr, threshold=700)
    }
    if (docum) {
      Hcal3 <- compute_Hi_Hredux_sixmonstatevector(docall, etmp, threshold=700)
      Hcal3[Hcal3==0] <- NA
    }
    rm(etmp)
    if (instrumental & !docum & sixmonstatevector & !real_proxies) {
      H.i <- array(NA,c(nrow(Hcal1),1))
    } else if (instrumental & docum & sixmonstatevector & !real_proxies) {
      H.i <- array(NA,c((nrow(Hcal1)+nrow(Hcal3)),1))
    } else if (!instrumental & docum & sixmonstatevector & real_proxies) {
      H.i <- array(NA,c((nrow(Hcal2)+nrow(Hcal3)),(ncol(Hcal2))/2))
    } else if (instrumental & !docum & sixmonstatevector & real_proxies) {
      H.i <- array(NA,c((nrow(Hcal1)+nrow(Hcal2)),(ncol(Hcal2))/2))
    } else if (!instrumental & !docum & real_proxies) { #new
      H.i <- array(NA,c(nrow(Hcal2),(ncol(Hcal2))/2))
    } else { 
      H.i <- array(NA,c((nrow(Hcal1)+nrow(Hcal2)+nrow(Hcal3)),(ncol(Hcal2))/2))
    }
    
    Hredux <- H.i
    if (instrumental){
      H.i[1:nrow(Hcal1),1] <- Hcal1[,1] 
      Hredux[1:nrow(Hcal1),1] <- Hcal1[,2] 
    } 
    if (instrumental & real_proxies & docum) {
      H.i[(nrow(Hcal1)+1):(nrow(Hcal1)+nrow(Hcal2)),] <- Hcal2[,seq(1,ncol(Hcal2),2)] 
      Hredux[(nrow(Hcal1)+1):(nrow(Hcal1)+nrow(Hcal2)),] <- Hcal2[,seq(2,ncol(Hcal2),2)]   
      H.i[((nrow(Hcal1)+nrow(Hcal2)+1):nrow(H.i)),1] <- Hcal3[,1]
      Hredux[((nrow(Hcal1)+nrow(Hcal2)+1):nrow(H.i)),1] <- Hcal3[,2]
    }
    if (!instrumental & real_proxies & docum) {
      H.i[1:nrow(Hcal2),] <- Hcal2[,seq(1,ncol(Hcal2),2)] 
      Hredux[1:nrow(Hcal2),] <- Hcal2[,seq(2,ncol(Hcal2),2)]   
      H.i[((nrow(Hcal2)+1):nrow(H.i)),1] <- Hcal3[,1]
      Hredux[((nrow(Hcal2)+1):nrow(H.i)),1] <- Hcal3[,2]
    }
    if (instrumental & real_proxies & !docum) {
      H.i[(nrow(Hcal1)+1):(nrow(Hcal1)+nrow(Hcal2)),] <- Hcal2[,seq(1,ncol(Hcal2),2)] 
      Hredux[(nrow(Hcal1)+1):(nrow(Hcal1)+nrow(Hcal2)),] <- Hcal2[,seq(2,ncol(Hcal2),2)] 
    }
    if (!instrumental & !docum & real_proxies){ #new
      H.i[1:nrow(Hcal2),]<-Hcal2[,seq(1,ncol(Hcal2),2)] 
      Hredux[1:nrow(Hcal2),]<-Hcal2[,seq(2,ncol(Hcal2),2)] 
    }
    H.i[H.i==0] <- NA
    Hredux[Hredux==0] <- NA
  }
  
  print("calc time H operator")
  print(proc.time() - ptm1)
  
  
  
  
  #########################################################################################
  # 8. Run analysis
  #########################################################################################
  
  if ((instrumental) || (real_proxies) || (docum)) {
    # R for perfect data would be:
    # R <- rep(0, nrow(Hcal))
    # R is simply the squared random error assumed for instr. data because we assume
    # the error would be spatially uncorrelated
    # set squared error R to 1 for 1degC measurement error
    if (sixmonstatevector) {
      if ((real_proxies) & ((instrumental) | (docum))) { 
        Rcal <- c(temp2=0.9, precip=50, slp=10)[calibrate$names]
        #        if (avg_prox_per_grid) {Rcal <- Rcal*(1/calibrate$numavg)}
        Rcal[calibrate$names=="prox"] <- realprox$var_residu
        # previously used residuals/2 for 1. paper version to give proxies more weight
        # better delete "/2"
        # probably should have given instrumentals more error instead!
        Rcal[calibrate$sour=="doc"] <- 0.25 # equals 0.5 std. dev.
      } else if ((real_proxies) & (!instrumental) & (!docum)) { 
        Rcal <- realprox$var_residu
        Rcal[which(is.na(Rcal))] <- 0
      } else if (((instrumental) | (docum)) & (!real_proxies)) { 
        Rcal <- c(temp2=0.9, precip=50, slp=10)[calibrate$names]
        Rcal[calibrate$sour=="doc"] <- 0.25 # equals 0.5 std. dev.
        #        if (avg_prox_per_grid) {Rcal <- Rcal*(1/calibrate$numavg)}
        #      # set squared error R for precip measurements to 25% of data value
        #      # R <- abs(calibrate$data*0.25)
      }
      
      analysis <- echam
      # take anomalies
      analysis$data <- echam$data - as.vector(echam$ensmean)
      if (cov_inflate){
        analysis$data= analysis$data * inflate_fac
      }
      if (covarclim>0 & covarclim<=100) { 
        ananomallts = echanomallts
        ananomallts$data <- echanomallts$data - as.vector(echanomallts$ensmean)
      }
      nmonths <- 6
      ndim <- nrow(analysis$data)
      ntim <- ncol(analysis$data)
      nens <- dim(analysis$data)[3]
      nprox <- nrow(calibrate$data) 
      ndimold <- length(echam$lon) / nmonths
      itime <- rep(1:nmonths, each=ndimold)
      if (landcorr) {
        corland_analysis <- landcorrected.anom
      }
      i=1 
      # compute loop over obs first to minimize repetition of computation (e.g. H and weights)
      for (j in 1:nprox){
        #if (j %% 100 == 0) {
        #  print(paste('Assimilating observation', j))
        #}
        ## assume proxy location in state vector is known
        ## otherwise specify / compute h.i and reduced-size H here
        hisna <- is.na(H.i[j,])
        h.i <- H.i[j,!hisna,drop=F]
        H <-  t(as.matrix(Hredux[j,!hisna]))
        if (!is.na(h.i[1])) {
          pos <- which(is.na(echam$lon))
          echam$lon[pos] <- 0
          echam$lat[pos] <- -90
          if (shape_wgt == "circle") {
            dist <- compute_dist_2d(echam$lon, echam$lat, echam$lon[h.i], echam$lat[h.i],region) 
            # weights are a matrix of ndim x nh (non-null elements in H, here different months)
            # NO temporal correlation for use of monthly instr. data and 1 col H operator          
            wgt <- corr_function(dist,outer(lvec[echam$names], lvec[echam$names[h.i]], pmin))
            wgt[which(echam$names %in% c("DIMI","z100","z300","PWC","HC","SJ")),] <- 1
            # temporal correlation quickly drops to zero (~ 0.6 for 1-month lag, 
            # ~0.4 for two months, etc.)
          } else if (shape_wgt == "ellipse") {
            wgt = compute_dist_2d_ellipse (echam$lon, echam$lat, echam$lon[h.i], echam$lat[h.i],echam$names,lvec,calibrate$sour[j])
          }
          if (covarclim>0 & covarclim<=100) { 
            if (PHclim_loc) {
              dist <- compute_dist_2d(echam$lon, echam$lat, echam$lon[h.i], echam$lat[h.i],region) 
              wgt_PHclim <- corr_function(dist,outer(PHclim_lvec[echam$names], PHclim_lvec[echam$names[h.i]], pmin))
              wgt_PHclim[which(echam$names %in% c("DIMI","z100","z300","PWC","HC","SJ")),] <- 1
            }
          }
        } else {
          dist <- NA
          wgt <- NA
        }
        if (!landcorr) {
          if (!is.na(wgt[1])){ # for shape_wgt = ellipse there is no dist
            for (i in 1:ntim){
              if (!is.na(calibrate$data[j,i])) {
                # if (cov_inflate){
                #   analysis$data[,i,] = analysis$data[,i,] * inflate_fac
                # }
                x2tmp <- analysis$data[,i,] # entire state vector at time step i, all members
                x2 <- x2tmp[h.i,,drop=F] # state vector at time step i and h.i, all members
                if(mixed_loc){ # changed for the B exps: first combining and than localizing
                  PH <- (analysis$data[,i,] %*% t(x2) / (nens - 1)) # %*% t(H)
                } else { # original
                  PH <- (analysis$data[,i,] %*% t(x2) / (nens - 1) * wgt) %*% t(H)
                }
                if (covarclim>0 & covarclim<100) { 
                  x2climtmp <- ananomallts$data[,i,] 
                  x2clim <- x2climtmp[h.i,,drop=F] 
                  if (mixed_loc) {
                    PHclim <- (ananomallts$data[,i,] %*% t(x2clim) / ((dim(ananomallts$data)[3]) - 1) ) #  %*% t(H)
                    PH <- (PH*(1-(covarclim/100))) + (PHclim*(covarclim/100))
                    PH <- (PH * wgt_PHclim) %*% t(H)
                  } else {
                    if (PHclim_loc) {
                      PHclim <- (ananomallts$data[,i,] %*% t(x2clim) /
                                   ((dim(ananomallts$data)[3]) - 1) * wgt_PHclim) %*% t(H)
                    } else {
                      PHclim <- (ananomallts$data[,i,] %*% t(x2clim) /
                                   ((dim(ananomallts$data)[3]) - 1) ) %*% t(H)
                    }
                    PH = (PH*(1-(covarclim/100))) + (PHclim*(covarclim/100))
                  }
                } else if (covarclim==100) {
                  x2climtmp <- ananomallts$data[,i,] 
                  x2clim <- x2climtmp[h.i,,drop=F] 
                  if (PHclim_loc) {
                    PH <- (ananomallts$data[,i,] %*% t(x2clim) / 
                             ((dim(ananomallts$data)[3]) - 1) * wgt_PHclim) %*% t(H) 
                  } else {
                    PH <- (ananomallts$data[,i,] %*% t(x2clim) / 
                             ((dim(ananomallts$data)[3]) - 1) ) %*% t(H) 
                  }
                } 
                if (ins_tim_loc) {
                  if (calibrate$sour[j] == "inst") {
                    # time localization on PH for instrumental data
                    # let influence only the variables of the current month of the observations
                    month_start = seq(1,dim(PH)[1],dim(PH)[1]/6)
                    month_end = month_start + (dim(PH)[1]/6 - 1)
                    time_loc = array(NA,dim(PH))
                    for (k in 1:6) {
                      if (month_start[k]<=h.i & month_end[k]>=h.i) {
                        time_loc[month_start[k]:month_end[k],] = 1
                      } else {
                        time_loc[month_start[k]:month_end[k],] = 0
                      }
                    }
                    PH <- time_loc * PH
                  }
                }
                HPHR <- as.vector(H %*% PH[h.i,] + Rcal[j])
                K <- PH / HPHR
                Ktilde <- K / (1 + sqrt(Rcal[j]/HPHR))
                analysis$ensmean[,i] <- analysis$ensmean[,i] + K[,1] * (calibrate$data[j,i] -
                                                                          H %*% analysis$ensmean[h.i,i])
                analysis$data[,i,] <- analysis$data[,i,] - Ktilde %*% H %*% analysis$data[h.i,i,]
                if (update_PHclim){
                  ananomallts$ensmean[,i] <- ananomallts$ensmean[,i] + K[,1] * (calibrate$data[j,i] -
                                                                                  H %*% ananomallts$ensmean[h.i,i])
                  ananomallts$data[,i,] <- ananomallts$data[,i,] - Ktilde %*% H %*% ananomallts$data[h.i,i,]
                }
              }
            }
          }  
        } else { # landcorr is TRUE
          if ((!is.na(dist[1])) & (!is.na(wgt[1]))){
            for (i in 1:ntim){
              if (!is.na(calibrate$data[j,i])) {
                x2tmp <- analysis$data[,i,] # entire state vector at time step i, all members
                x2 <- x2tmp[h.i,,drop=F] # state vector at time step i and h.i, all members
                PH <- (analysis$data[,i,] %*% t(x2) / (nens - 1) * wgt) %*% t(H)
                HPHR <- as.vector(H %*% PH[h.i,] + Rcal[j])
                K <- PH / HPHR
                corland_analysis$data[,i,] <- corland_analysis$data[,i,] + K[,1] %*%
                  (calibrate$data[j,i] - H %*% corland_analysis$data[h.i,i,])
                # get warrning although didn't get previously
                # In K[, 1] * (calibrate$data[j, i] - H %*% corland_analysis$data[h.i,  ... :
                # Recycling array of length 1 in vector-array arithmetic is deprecated.
                # Use c() or as.vector() instead.
                # I think it is because dim(K) = 304128 1, and we multiply it only with a matrix [1 1] -> though still giving the good result
                # using the inner product gives the same result, without warning
              }
            }
          } 
        } 
      }
      ## add ensemble mean analysis back 
      if (!landcorr) {
        analysis$data <- analysis$data + as.vector(analysis$ensmean)
        if (save_ananomallts == T) {
          ananomallts$data <- ananomallts$data + as.vector(ananomallts$ensmean) 
        }
        if (anomaly_assim){
          analysis.anom <- analysis
          analysis.abs <- analysis
          analysis.abs$data <- analysis$data + as.vector(echam.clim$data)
          analysis.abs$ensmean <- analysis$ensmean + as.vector(echam.clim$ensmean)
          echam.anom <- echam
          echam.abs <- echam
          echam.abs$data <- echam$data + as.vector(echam.clim$data)
          echam.abs$ensmean <- echam$ensmean + as.vector(echam.clim$ensmean)
        }
      } else { # landcorr part
        if (anomaly_assim){
          analysis.anom <- corland_analysis
          analysis.abs <- corland_analysis
          analysis.abs$data <- corland_analysis$data + as.vector(landcorrected.clim$data)
          echam.anom <- echam
          echam.abs <- echam
          echam.abs$data <- echam$data + as.vector(echam.clim$data)
          echam.abs$ensmean <- echam$ensmean + as.vector(echam.clim$ensmean)
        }
      }
      if (!landcorr) {
        if (vali){
          if (every2grid){
            if(save_ananomallts == T) {
              save(analysis.anom,analysis.abs,echam.anom,echam.abs,validate,calibrate, ananomallts, 
                   file=paste0(dataintdir,'analysis/',expname,'/analysis_',cyr,'_2ndgrid.Rdata'))
            } else {
              save(analysis.anom,analysis.abs,echam.anom,echam.abs,validate,calibrate, 
                   file=paste0(dataintdir,'analysis/',expname,'/analysis_',cyr,'_2ndgrid.Rdata'))
            }
          } else {
            save(analysis.anom,analysis.abs,echam.anom,echam.abs,validate,calibrate,
                 file=paste0(dataintdir,'analysis/',expname,'/analysis_',cyr,'.Rdata'))
          }
        } else {
          if (every2grid){    
            save(analysis.anom,analysis.abs,echam.anom,echam.abs,calibrate,
                 file=paste0(dataintdir,'analysis/',expname,'/analysis_',cyr,'_2ndgrid.Rdata'))
          } else {
            save(analysis.anom,analysis.abs,echam.anom,echam.abs,calibrate,
                 file=paste0(dataintdir,'analysis/',expname,'/analysis_',cyr,'.Rdata'))
          }
        }
      } else {
        if (vali) {
          save(analysis.anom,analysis.abs,echam.anom,echam.abs,landcorrected.anom,
               landcorrected.clim,validate,calibrate,
               file=paste0(dataintdir,'analysis/',expname,'/analysis_',cyr,'_2ndgrid.Rdata'))
        } else {
          save(analysis.anom,analysis.abs,echam.anom,echam.abs,landcorrected.anom,
               landcorrected.clim,calibrate,
               file=paste0(dataintdir,'analysis/',expname,'/analysis_',cyr,'_2ndgrid.Rdata'))
        }   
      }
      if (loo) { # leave one out validation
        ptm2 <- proc.time()
        proxyloc <- getgridboxnum(calibrate, echam)
        proxyloc <- proxyloc[1:length(which(calibrate$sour=="inst"))]
        np <- length(proxyloc)/6
        ne <- dim(echam$data)[1]/6
        proxyloc <- proxyloc+rep(ne*seq(0,5),each=np)
        ll <- np*seq(0,5)
        ltmp <- 0
        for (l in 1:np) {
          print(paste(l,"of",np,"loo runs"))
          print(proc.time() - ptm2)      
          lloc <- l+ll
          if (!all(is.na(calibrate$data[lloc,]))) {
            ltmp <- ltmp+1
            print(ltmp)
            loocalibrate <- calibrate
            loocalibrate$data <- calibrate$data[-lloc,]
            loocalibrate$lon <- calibrate$lon[-lloc]
            loocalibrate$lat <- calibrate$lat[-lloc]
            looH.i <- H.i[-lloc,]
            looHredux <- Hredux[-lloc,]
            looRcal <- Rcal[-lloc]
            looanalysis <- echam
            # take anomalies
            looanalysis$data <- echam$data - as.vector(echam$ensmean)
            nmonths <- 6
            ndim <- nrow(looanalysis$data)
            ntim <- ncol(looanalysis$data)
            nens <- dim(looanalysis$data)[3]
            nprox <- nrow(loocalibrate$data) 
            ndimold <- length(echam$lon) / nmonths
            itime <- rep(1:nmonths, each=ndimold)
            i=1
            # compute loop over obs first to minimize repetition of computation (e.g. H and weights)
            for (j in 1:nprox){
              if (j %% 100 == 0) {
                print(paste('Assimilating observation', j))
              }
              ## assume proxy location in state vector is known
              ## otherwise specify / compute h.i and reduced-size H here
              hisna <- is.na(looH.i[j,])
              h.i <- looH.i[j,!hisna,drop=F]
              H <-  t(as.matrix(looHredux[j,!hisna]))
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
                  if (!is.na(loocalibrate$data[j,i])) {
                    x2tmp <- looanalysis$data[,i,] # entire state vector at time step i, all members
                    x2 <- x2tmp[h.i,,drop=F] # state vector at time step i and h.i, all members
                    PH <- (looanalysis$data[,i,] %*% t(x2) / (nens - 1) * wgt) %*% t(H)
                    HPHR <- as.vector(H %*% PH[h.i,] + looRcal[j])
                    K <- PH / HPHR
                    Ktilde <- K / (1 + sqrt(looRcal[j]/HPHR))
                    looanalysis$ensmean[,i] <- looanalysis$ensmean[,i] + K[,1] * (loocalibrate$data[j,i] -
                                                                                    H %*% looanalysis$ensmean[h.i,i])
                    looanalysis$data[,i,] <- looanalysis$data[,i,] - Ktilde %*% H %*% 
                      looanalysis$data[h.i,i,]
                  }
                } # end time step loop
              } # if dist and wgt not NA 
            } # end serial assimilation
          } # end if cali is not NA
          # save ID (proxy number), lon/lat and difference between grid box temp/slp/precip 
          # of full analysis and analysis minus station at this grid box
          # every time step add columns for 6 new months to matrix (nrow=#ofproxies)
          if (ltmp==1) {
            looensmeanres <- c(calibrate$lon[l],calibrate$lat[l],
                               as.vector(looanalysis$ensmean[proxyloc[lloc],])-
                                 as.vector(calibrate$data[lloc,]),
                               as.vector(echam$ensmean[proxyloc[lloc],]),
                               as.vector(analysis$ensmean[proxyloc[lloc],]),
                               as.vector(looanalysis$ensmean[proxyloc[lloc],]),
                               as.vector(calibrate$data[lloc,]))
            loodatares <- c(calibrate$lon[l],calibrate$lat[l],
                            as.vector(looanalysis$data[proxyloc[lloc],,])-
                              rep(as.vector(calibrate$data[lloc,]),30),
                            as.vector(echam$data[proxyloc[lloc],,]),
                            as.vector(analysis$data[proxyloc[lloc],,]),
                            as.vector(looanalysis$data[proxyloc[lloc],,]))
            looensmeanres <- t(looensmeanres)
            loodatares <- t(loodatares)
            colnames(looensmeanres) <- c("lon","lat","diff apr","diff may","diff jun",
                                         "diff jul","diff aug","diff sep","diff oct","diff nov","diff dec",
                                         "diff jan","diff feb","diff mar","cali apr","cali may","cali jun",
                                         "cali jul","cali aug","cali sep","cali oct","cali nov","cali dec",
                                         "cali jan","cali feb","cali mar","ech apr","ech may","ech jun",
                                         "ech jul","ech aug","ech sep","ech oct","ech nov","ech dec",
                                         "ech jan","ech feb","ech mar","ana apr","ana may","ana jun",
                                         "ana jul","ana aug","ana sep","ana oct","ana nov","ana dec",
                                         "ana jan","ana feb","ana mar","loo apr","loo may","loo jun",
                                         "loo jul","loo aug","loo sep","loo oct","loo nov","loo dec",
                                         "loo jan","loo feb","loo mar")
            colnames(loodatares) <- c("lon","lat",paste(rep(seq(1,30),each=12),rep(c(
              "diff apr","diff may","diff jun","diff jul","diff aug","diff sep",
              "diff oct","diff nov","diff dec","diff jan","diff feb","diff mar"),30)),
              paste(rep(seq(1,30),each=12),rep(c("ech apr","ech may","ech jun",
                                                 "ech jul","ech aug","ech sep","ech oct","ech nov","ech dec",
                                                 "ech jan","ech feb","ech mar"),30)),
              paste(rep(seq(1,30),each=12),rep(c("ana apr","ana may","ana jun",
                                                 "ana jul","ana aug","ana sep","ana oct","ana nov","ana dec",
                                                 "ana jan","ana feb","ana mar"),30)),
              paste(rep(seq(1,30),each=12),rep(c("loo apr","loo may","loo jun",
                                                 "loo jul","loo aug","loo sep","loo oct","loo nov","loo dec",
                                                 "loo jan","loo feb","loo mar"),30)))
          } else if (ltmp>1) {
            looensmeanres <- rbind(looensmeanres,c(calibrate$lon[l],calibrate$lat[l],
                                                   as.vector(looanalysis$ensmean[proxyloc[lloc],])-
                                                     as.vector(calibrate$data[lloc,]),
                                                   as.vector(echam$ensmean[proxyloc[lloc],]),
                                                   as.vector(analysis$ensmean[proxyloc[lloc],]),
                                                   as.vector(looanalysis$ensmean[proxyloc[lloc],]),
                                                   as.vector(calibrate$data[lloc,])))
            loodatares <- rbind(loodatares,c(calibrate$lon[l],calibrate$lat[l],
                                             as.vector(looanalysis$data[proxyloc[lloc],,])-
                                               rep(as.vector(calibrate$data[lloc,]),30),
                                             as.vector(echam$data[proxyloc[lloc],,]),
                                             as.vector(analysis$data[proxyloc[lloc],,]),
                                             as.vector(looanalysis$data[proxyloc[lloc],,])))
          }
        } # end of loo loop
        save(looensmeanres,loodatares,
             file=paste0(dataintdir,'loo/',expname,'/loo_results_',cyr,'.Rdata'))
        # 
        #         # quick plots to check loo
        #         #cyr=1811
        #         load(paste0(dataintdir,'loo/loo_results_',cyr,'.Rdata'))
        #         #lev <- pretty(looensmeanres[,3:14],11)
        #         lev <- c(-Inf,-2,-1,-0.5,0.5,1,2,Inf) #pretty(looensmeanres[,3:14],11)
        #         br <- length(lev)
        #         colpal <- two.colors(n=br,start="darkblue", end="darkred", alpha=0.5)
        #         colpal[4:5] <- "#FFFFFFFF"
        #         pdf(file=paste0("../figures/loo/test_",cyr,".pdf"),width=5,height=36)
        #           set.panel(12,1)
        #           par(oma=c(.1,.1,.1,.1))
        #           for (i in 3:14){
        #             datcol <- colpal[as.numeric(cut(looensmeanres[,i],breaks=br))]
        #             #plot(NULL,xlim=c(-170,40),ylim=c(20,80),xlab="lon",ylab="lon")
        #             plot(NULL,xlim=c(-180,180),ylim=c(-60,80),xlab="lon",ylab="lon")
        #             points(looensmeanres[,1], looensmeanres[,2],cex=0.5,pch=15,col=datcol,
        #               xlab="Longitude",ylab="Latitude")                # point fill color 
        #             map("world",interior=F,add=T,ty='l',col='black',
        #               xlim=c(-180,180),ylim=c(-60,80))
        #               #xlim=c(-170,40),ylim=c(20,80))
        #             legend("bottomleft", inset=0.01, as.character(lev),fill=colpal,bty="o",
        #               bg="white",box.col="white",box.lwd=0,cex=0.7)
        #           }
        #         dev.off()
        print('calc time for leave one out')
        print(proc.time() - ptm2)
      } # end loo if
      
    } else {      
      R = array(1,c(nrow(Hcal), s))
      R[calibrate$names=='precip',]<-try(abs(calibrate$data[calibrate$names=='precip',]
                                             *0.5))  
      R[calibrate$names=='slp',]<-try(abs(calibrate$data[calibrate$names=='slp',]
                                          *0.1))  
      R[which(is.na(R))] <- 0
      if (fortran_ensrf) {calibrate$data[is.na(calibrate$data)]=9e36}  
      analysis <- EnSRF_new(echam, calibrate,  R=R, Hcal=Hcal, weights=d.weights_all)
      #      analysis <- EnSRF4(echam, calibrate,  R=R, Hcal=Hcal, weights=d.weights_all)
      if (vali) {
        if (every2grid){
          save(analysis,echam,validate,calibrate, #calibrate.allts,
               file=paste0('../data/analysis/',expname,'/analysis_',cyr,'_2ndgrid.Rdata'))
        } else {
          save(analysis,echam,validate,calibrate, #calibrate.allts,
               file=paste0('../data/analysis/',expname,'/analysis_',cyr,'.Rdata'))
        }
      } else {
        if (every2grid){
          save(analysis,echam,calibrate, #calibrate.allts,
               file=paste0('../data/analysis/',expname,'/analysis_',cyr,'_2ndgrid.Rdata'))
        } else {  
          save(analysis,echam,calibrate, #calibrate.allts,
               file=paste0('../data/analysis/',expname,'/analysis_',cyr,'.Rdata'))
        }
      }
    }
  }
  print("calc time for a year")
  print(proc.time() - ptm1)
} # end of yearly analysis loop

warnings()
#quit(save='no')


