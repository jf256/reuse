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



syr=1950
eyr=1951


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
#} else if (user=="lucaf") {
#  workdir='/scratch3/lucaf/reuse/reuse_git/'
} else if (user=="joerg") {
  workdir='/scratch3/joerg/projects/reuse/git/'
#} else if (user == "nevin"){
#  workdir = '/scratch3/nevin/reuse_climcal/reuse_git/'
} else{
  stop("Unknown user!")
  
}
dataextdir='/mnt/climstor/giub/EKF400/'
dataintdir=paste0(workdir,'../data/')
setwd(workdir)

source('EnSRF_switches.R')
source('EnSRF_functions.R')

dir.create(paste0("../data/analysis/EKF400_",version,'_',expname))
# create logfile
logfn <- paste0("EKF400_",version,'_',expname,'_',syr,'-',eyr,'_',format(Sys.time(),"%Y%m%d_%H%M"),'.log')
write(c(user),file=paste0('../log/',logfn),append=F)
if (loo) {dir.create(paste0("../data/loo/",expname))}

# it can be useful to save how the switches were set
con <- file(paste0(workdir,"../data/analysis/EKF400_",version,'_',expname,"/switches_",format(Sys.time(),"%Y%m%d_%H%M"),".log"))
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
  if (cyr > 1659 & (yuri_temp|yuri_slp|ghcn_temp|isti_instead_ghcn|ghcn_prec)) {
    #Nevin June2018: I added this control sequence instead of the "proxies_only" in the experiment name. Like that: at least one of the instrumental
    #datasets has to be true in the switches such that instrumental = T here. Otherwhise there was an error when only proxies were used.
    instrumental=T        # all instrumental stations 
  } else {
    instrumental=F
  }

  if (TRW|MXD|SCHWEINGR|PAGES|NTREND|TRW_PETRA|pseudo_prox) {        
    real_proxies=T         # Proxy data experiment (regression NOT H operator) 
  } else {
    real_proxies=F
  }
  if (cyr <= 1900 & pseudo_prox) {
    real_proxies=F
    stop('DAPS pseudo proxy experiment is currently only working from 1901 onwards')
  }
  if (assim_docu) { # cyr <= 1853 & assim_docu) { #import_luca) {
    docum=T                 # read documentary based data
  } else {
    docum=F
  }
  
  # next line not included yet: 
  if (cyr > min(c(syr_cru,syr_twentycr,syr_recon)[c(vali_cru, vali_twentycr, vali_recon)]) & 
      cyr <= max(c(eyr_cru,eyr_twentycr,eyr_recon)[c(vali_cru, vali_twentycr, vali_recon)])) {        # if we don't use reconvali, the eyr here should be changed (Error in valiall : object 'valiall' not found) -> but then instead of the eyr we should use cyr
    vali=T                 # switch off prepplot if no vali data selected
  } else {
    vali=F
  }
  if ((cyr > syr_cru) & (cyr <=eyr_cru) & vali_cru) {
    cru_vali=T             # monthly CRU TS3 temp, precip and HADSLP2 gridded instrumentals (1901-2004)
    #  ind_recon=T         # Stefan's reconstructed indices until 1948 and NCAR reanalysis later added to CRU
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
  # 1.1 Loading (echam commented), echam_anom, echam_clim, landcorrected_anom, landcorrected_clim
  # 1.1.2 Load/create bigger ensemble size 
  # 1.2 Choose which variables want to use from the model
  # 1.3 Calc echam st. dev. for each grid point and month over ens memb. to scale docu data
  # 1.4 Calculate decorrelation distance
  # 1.5 Just leave data for one year in memory and convert to sixmonstatevector data format
  # 1.6 Set up the cutoff distance for localisation
  ##########################################################################################
  
  # 1.1 Loading echam, echam_anom, echam_clim, landcorrected_anom, landcorrected_clim
  # JF commented 09/2019
  #if (every2grid) {
  #  load(paste(echpath,"echam_",(cyr-1),"-",cyr,"_2ndgrid.Rdata",sep=""))
  #  #load(paste(dataextdir,"echam/echam_",(cyr-1),"-",cyr,"_2ndgrid.Rdata",sep=""))
  #  #load(paste0('/scratch3/veronika/reuse/data/echam/new_statvec/echam/echam_',(cyr-1),'-',cyr,'_2ndgrid.Rdata'))
  #} else {
  #  load(paste0(echpath,"echam_",(cyr-1),"-",cyr,".Rdata"))
  #  # load(paste0(dataextdir,"echam/echam_",(cyr-1),"-",cyr,".Rdata"))
  #}
  if ((anomaly_assim) & (!no_forc_big_ens)) {  
    # anomalies calculated efficiently with cdo, slow calculation within R has been removed
    yr1 <- cyr-1
    yr2 <- cyr
    # since the anomalies in the ncdf files are calculated by using the running mean 
    #   (first year only 36 yrs used) by using this we are not consistent
    # yr3 <- yr1
    # yr4 <- yr2
    # if (cyr < 1637) {yr3 <- 1636}
    # if (cyr < 1637) {yr4 <- 1637}
    # if (cyr > 1970) {yr3 <- 1969}
    # if (cyr > 1970) {yr4 <- 1970}
    if (every2grid) {
      load(paste0(echanompath,'echam_anom_',yr1,'-',yr2,'_2ndgrid.Rdata'))
      ## load(paste0(echclimpath,'echam_clim_',yr3,'-',yr4,'_2ndgrid.Rdata')) # old version
      #load(paste0('/scratch3/veronika/reuse/data/echam/new_statvec/echam_anom/echam_anom_',yr1,'-',yr2,'_2ndgrid.Rdata'))
      load(paste0(echclimpath,'echam_clim_',yr1,'-',yr2,'_2ndgrid.Rdata'))
      #load(paste0('/scratch3/veronika/reuse/data/echam/new_statvec/echam_clim/echam_clim_',yr1,'-',yr2,'_2ndgrid.Rdata'))
      if (landcorr) {
        load(paste0(echanompath,'echam_anom_103_',yr1,'-',yr2,'_2ndgrid.Rdata'))
        load(paste0(echclimpath,'echam_clim_103_',yr1,'-',yr2,'_2ndgrid.Rdata'))
      }
    } else {
      load(paste0(echanompath,'echam_anom_',yr1,'-',yr2,'.Rdata'))
      load(paste0(echclimpath,'echam_clim_',yr1,'-',yr2,'.Rdata'))
      if (landcorr) {
        load(paste0(echanompath,'echam_anom_103_',yr1,'-',yr2,'.Rdata'))
        load(paste0(echclimpath,'echam_clim_103_',yr1,'-',yr2,'.Rdata'))
      }
    }
    # correct units if old state vector with error
    if (old_statvec) {
      echam_anom$data[echam_anom$names=='precip',,] <- echam_anom$data[echam_anom$names=='precip',,] * 3600 * 24 * 30
      echam_anom$ensmean[echam_anom$names=='precip',] <- echam_anom$ensmean[echam_anom$names=='precip',] * 3600 * 24 * 30
      echam_anom$data[echam_anom$names=='slp',,] <- echam_anom$data[echam_anom$names=='slp',,]/100
      echam_anom$ensmean[echam_anom$names=='slp',] <- echam_anom$ensmean[echam_anom$names=='slp',]/100
    }
  }
  
  # 1.1.2 Load/create bigger ensemble size 
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
  }
  
  
  # 1.2 Choose which variables want to use from the model
  # just leave temp precip slp in state vector
  if (tps_only) {
    # JF commented 09/2019
    #tpspos <- c(which(echam$names=='temp2'), which(echam$names=='precip'),
    #            which(echam$names=='slp'))
    #echam$data <- echam$data[tpspos,,]
    #echam$ensmean <- echam$ensmean[tpspos,]
    #echam$names <- echam$names[tpspos]
    if (anomaly_assim){
      tpspos <- c(which(echam_anom$names=='temp2'), which(echam_anom$names=='precip'),
                  which(echam_anom$names=='slp'))
      echam_anom$data <- echam_anom$data[tpspos,,]
      echam_anom$ensmean <- echam_anom$ensmean[tpspos,]
      echam_anom$names <- echam_anom$names[tpspos]
      if (state == "changing" & covarclim>0) {
        tpspos <- c(which(echanomallts$names=='temp2'), which(echanomallts$names=='precip'),
                    which(echanomallts$names=='slp'))
        echanomallts$data <- echanomallts$data[tpspos,,]
        echanomallts$ensmean <- echanomallts$ensmean[tpspos,]
        echanomallts$names <- echanomallts$names[tpspos]
      }
      tpspos <- c(which(echam_clim$names=='temp2'), which(echam_clim$names=='precip'),
                  which(echam_clim$names=='slp'))
      echam_clim$data <- echam_clim$data[tpspos,,]
      echam_clim$ensmean <- echam_clim$ensmean[tpspos,]
      echam_clim$names <- echam_clim$names[tpspos]
    }
  }

  if (tpsw_only) {
    # JF commented 09/2019
    #tpswpos <- c(which(echam$names=='temp2'), which(echam$names=='precip'), which(echam$names=='slp'),
    #            which(echam$names=='wdays'), which(echam$names=='bias'))
    #echam$data <- echam$data[tpswpos,,]
    #echam$ensmean <- echam$ensmean[tpswpos,]
    #echam$names <- echam$names[tpswpos]
    if (anomaly_assim){
      tpswpos <- c(which(echam_anom$names=='temp2'), which(echam_anom$names=='precip'), which(echam_anom$names=='slp'),
                  which(echam_anom$names=='wdays'), which(echam_anom$names=='bias'))
      echam_anom$data <- echam_anom$data[tpswpos,,]
      echam_anom$ensmean <- echam_anom$ensmean[tpswpos,]
      echam_anom$names <- echam_anom$names[tpswpos]
      if (state == "changing" & covarclim>0) {
        tpswpos <- c(which(echanomallts$names=='temp2'), which(echanomallts$names=='precip'),
                    which(echanomallts$names=='slp'), which(echam$names=='wdays'), which(echanomallts$names=='bias'))
        echanomallts$data <- echanomallts$data[tpswpos,,]
        echanomallts$ensmean <- echanomallts$ensmean[tpswpos,]
        echanomallts$names <- echanomallts$names[tpswpos]
      }
      # Blending was not tested
      # if (state == "static" & no_forc_big_ens ) {
      #   tpspos <- c(which(echam_clim$names=='temp2'), which(echam_clim$names=='precip'),
      #               which(echam_clim$names=='slp'), which(echam_clim$names=='bias'))
      # }
      echam_clim$data <- echam_clim$data[tpswpos,,]
      echam_clim$ensmean <- echam_clim$ensmean[tpswpos,]
      echam_clim$names <- echam_clim$names[tpswpos]
    }
  }

  if (no_stream) {
    # ACHTUNG stream var has ERROR because the 5/9 levels before/after 1880 have a lat dimension
    # JF commented 09/2019
    #tpspos <- c(which(echam$names!='stream'))
    #echam$data <- echam$data[tpspos,,]
    #echam$ensmean <- echam$ensmean[tpspos,]
    #echam$names <- echam$names[tpspos]
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
      }
      tpspos <- c(which(echam_clim$names!='stream'))
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

  if (anomaly_assim){
    # new echam_clim made by veronika have correct units
    # error where echam_anom and echam_clim were generated initially: 
    # no unit correction happened, thus here
    if (!(state =="static" & no_forc_big_ens)) {
       # if (old_statvec) {
       #   echam_anom$data[echam_anom$names=='precip',,] <- echam_anom$data[echam_anom$names=='precip',,] * 3600 * 24 * 30
       #   echam_anom$ensmean[echam_anom$names=='precip',] <- echam_anom$ensmean[echam_anom$names=='precip',] * 3600 * 24 * 30
       #   echam_anom$data[echam_anom$names=='slp',,] <- echam_anom$data[echam_anom$names=='slp',,]/100
       #   echam_anom$ensmean[echam_anom$names=='slp',] <- echam_anom$ensmean[echam_anom$names=='slp',]/100
       # }
       echam.anom <- echam_anom
       if (state=="changing" & covarclim > 0) {
         if(old_statvec) {
           echanomallts$data[echanomallts$names=='precip',,] <- echanomallts$data[echanomallts$names=='precip',,] * 3600 * 24 * 30
           echanomallts$ensmean[echanomallts$names=='precip',] <- echanomallts$ensmean[echanomallts$names=='precip',] * 3600 * 24 * 30
           echanomallts$data[echanomallts$names=='slp',,] <- echanomallts$data[echanomallts$names=='slp',,]/100
           echanomallts$ensmean[echanomallts$names=='slp',] <- echanomallts$ensmean[echanomallts$names=='slp',]/100
         }
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
    if (state=="changing") { # Roni: why I delete it only if state = "changing ?? 
                             # -> I think because if it is "static" I need it in 
                             #    the background_matrix function
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
  print('calc time reading echam and ev. calc anom.')
  print(proc.time() - ptm1)
  
  
  # 1.3 Calc echam st. dev. for each grid point and month over ens memb. to scale docu data
  # JF added .anom 09/2019
  echam.sd <- echam.anom # to have non-sixmon echam in docum section -> I would move it to 1.3
  if (sixmonstatevector) {
    echam.sd$data <- apply(echam.anom$data[,10:21,],1:2,sd)
    echam.sd$time <- echam.anom$time[10:21]
  } else {
    echam.sd$data <- apply(echam.anom$data[,13:24,],1:2,sd)
    echam.sd$time <- echam.anom$time[13:24]
  }
  echam.sd$ensmean <- NULL
  print('calc time for standard deviations')
  print(proc.time() - ptm1)
  
  
  # 1.4 Calculate decorrelation distance
  if (calc_decorr_dist) { 
    if (covarclim == 0) {
      # JF added .anom 09/2019
      d <- compute_dist_2d(echam.anom$lon,echam.anom$lat,echam.anom$lon,echam.anom$lat, region) 
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
  #     or annual mean for DAPS pseudoproxy experiment
  if (pseudo_prox) {
    # JF commented 09/2019
    #echam$data <- array(apply(echam$data[,13:24,],c(1,3),mean),dim=c(dim(echam$data)[1],1,dim(echam$data)[3]))  
    #echam$ensmean <- array(apply(echam$ensmean[,13:24],1,mean),dim=c(dim(echam$data)[1],1))    
    #echam$time <- floor(echam$time[13])
    if (anomaly_assim) {
      echam.anom$data <- array(apply(echam.anom$data[,13:24,],c(1,3),mean),
                               dim=c(dim(echam.anom$data)[1],1,dim(echam.anom$data)[3]))  
      echam.anom$ensmean <- array(apply(echam.anom$ensmean[,13:24],1,mean),dim=c(dim(echam.anom$data)[1],1))  
      echam.anom$time <- floor(echam.anom$time[13]) 
      if (no_forc_big_ens) {
        echam.clim$data <- array(apply(echam.clim$data[,1:12,],1,mean),
                                 dim=c(dim(echam.clim$data)[1],1,dim(echam.clim$data)[3])) 
        #echam.clim.tmp=echam.clim$data
        #for (i in 1:29) {
        #  echam.clim.tmp=abind(echam.clim.tmp,echam.clim$data,along=3)
        #}
        #echam.clim$data <- echam.clim.tmp
        #rm(echam.clim.tmp)
        echam.clim$ensmean <- array(apply(echam.clim$ensmean[,1:12],1,mean),dim=c(dim(echam.clim$data)[1],1))  
        echam.clim$time <- cyr 
      } else {
        echam.clim$data <- array(apply(echam.clim$data[,13:24,],c(1,3),mean),
                                 dim=c(dim(echam.clim$data)[1],1,dim(echam.clim$data)[3]))  
        echam.clim$ensmean <- array(apply(echam.clim$ensmean[,13:24],1,mean),dim=c(dim(echam.clim$data)[1],1))  
        echam.clim$time <- floor(echam.clim$time[13]) 
      }
    }
    if (state=="changing" & (covarclim > 0)) {
      echanomallts$data <- array(apply(echanomallts$data[,13:24,],c(1,3),mean),
                                 dim=c(dim(echanomallts$data)[1],1,dim(echanomallts$data)[3]))    
      echanomallts$ensmean <- array(apply(echanomallts$ensmean[,13:24],1,mean),dim=c(dim(echanomallts$data)[1],1))      
      echanomallts$time <- floor(echanomallts$time[13]) 
    }
  }
  
  if (sixmonstatevector) {
    # change array to have 6 months in state vector for winter and summer
    # first winter starts in oct of syr
    # 6 mon stat vectors for oct-mar and apr and sep
    
    # JF commented 09/2019
    #echam <-convert_to_2_seasons(echam,source="echam")
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
  if (pseudo_prox) { 
    one_var_dim<-length(echam$lon)
    numvar <- length(unique(echam$names)) 
    echam$lon <-  c(rep(echam$lon, (numvar)))
    echam$lat <- c(rep(echam$lat, (numvar)))
  }
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
  lvec['omega500'] <- l_dist_omega500
  if (old_statvec) {
    # lvec['t850'] <- l_dist_t850 # replaced by t500
    lvec['t500'] <- l_dist_t500
    lvec['v200'] <- l_dist_v200
  }
  if (new_statvec) {
    lvec['wdays'] <- l_dist_wdays
    lvec['blocks'] <- l_dist_blocks
    lvec['cycfreq'] <- l_dist_cycfreq
  }
  
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
        #print("cru is true")
      } else if (v=="recon_vali"){
        recon_vali=T
        #print("recon is true")
      } else if (v=="twentycr_vali"){
        twentycr_vali=T
        #print("20cr is true")
      }
      
      if (every2grid) {
        if (recon_vali) {load(paste(dataextdir,"vali_data/recon/recon_allvar_",syr_recon,"-",eyr_recon,"_2ndgrid.Rdata",sep=""))
        } else if (cru_vali) {load(paste(dataextdir,"vali_data/cru/cru_allvar_",syr_cru,"-",eyr_cru,"_2ndgrid.Rdata",sep=""))
        # } else if (cru_vali) {load(paste(dataextdir,"vali_data/cru/cru_allvar_",syr_cru,"-",eyr_cru,"_2ndgrid_v2017.Rdata",sep=""))
        } else if (twentycr_vali){load(paste0(twentycrpath,"twentycr_allvar_",syr_twentycr,"-",eyr_twentycr,"_2ndgrid.Rdata"))}
      } else {
        if (recon_vali) {load(paste(dataextdir,"vali_data/recon/recon_allvar_",syr_recon,"-",eyr_recon,".Rdata",sep=""))
        } else if (cru_vali) {load(paste(dataextdir,"vali_data/cru/cru_allvar_",syr_cru,"-",eyr_cru,".Rdata",sep="")) 
        } else if (twentycr_vali){load(paste0(twentycrpath,"twentycr_allvar_",syr_twentycr,"-",eyr_twentycr,".Rdata"))}
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
    
      # 2.3 Cut out the 24 months around current year 
      valiall.allts=valiall
      ti=which(floor(valiall$time)==(cyr-1) | floor(valiall$time)==cyr) 
      sts=ti[1]
      ets=ti[length(ti)]
      valiall$data=valiall$data[,sts:ets]
      valiall$time=valiall$time[sts:ets]
      
      # 2.4 Convert it to 2 season per year or annual in case of DAPS pseudoproxy experiment
      if (pseudo_prox) {
        valiall$data <- array(apply(valiall$data[,13:24],1,mean,na.rm=T),dim=c(dim(valiall$data)[1],1))
        valiall$time <- floor(valiall$time[13])
      }
      # MISSING part for not sixmonstatevector or DAPS part of monthly analysis with 1-mon statevector
      if (sixmonstatevector) {
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
  # 3.2.1 Only use best proxies per grid cell, i.e. smallest residuals
  # 3.3 Convert to 2 season per year
  # 3.4 Calculate the anomalies
  # 3.5 Create a list named proxies ("combine assimilation data")
  # 3.6 Set real_proxies to FALSE if there is no data
  ##########################################################################################
  
  # 3.1 Loading proxy data
  if (real_proxies){
    if (pseudo_prox) {
      load(file=paste0("../data/proxies/DAPS_pseudoproxies_",fsyr,"-",feyr,".Rdata"))
    } else {
      if (!generate_PROXIES & (length(regression_months)==12) & (length(type)==2) & NTREND & TRW_PETRA & PAGES) {
        load(paste0(dataextdir,"assimil_data/rdata_files_v2/real_proxies_NTREND_PETRA_PAGES_tree_coral_1602-2004.Rdata"))
      } else {  
        load(paste0("../data/proxies/real_proxies_",expname,"_",fsyr,"-",feyr,".Rdata"))
      }
    }
    
    # 3.2 Screen the proxy data
    if (check_assimdata) {
      # correlation screening already where multiple regression coefficients are calculated
      # thus, screen if value at current time step is more than 5 std. dev. from mean
      # in this case treated as outlier and set to NA
      realprox <- screenstd(realprox,cyr=cyr,source="proxy")
      realprox$data <- t(realprox$data)
    }
    
    # 3.2.1 Only use best proxies per grid cell, i.e. smallest residuals
    if (best_tree_per_grid) {
      res=besttreeres
      newgrid <- echam.sd
      newgrid$ensmean <- NULL
      #newgrid$lon <- echam.sd$lon[seq(1,length(echam.sd$lon),res)]
      #newgrid$lat <- echam.sd$lat[seq(1,length(echam.sd$lat),res)]
      lons <- seq((-180+(res/20)),(180-(res/20)),res)
      lats <- seq((-90+(res/20)),(90-(res/20)),res)
      newgrid$lon <- rep(lons, length(lats))
      newgrid$lat <- rep(lats, each=length(lons))
      newgrid$data <- as.matrix(rep(NA,length(newgrid$lon)))
      
      tpos=which(realprox$archivetype=='tree')
      cpos=which(realprox$archivetype=='coral')
      if (length(cpos) > 0) {
        # save corals
        coralprox <- realprox
        coralprox$data <- realprox$data[,cpos]
        coralprox$lon <- realprox$lon[cpos]
        coralprox$lat <- realprox$lat[cpos]
        coralprox$archivetype <- realprox$archivetype[cpos]
        coralprox$elevation <- realprox$elevation[cpos]
        coralprox$mr <- realprox$mr[cpos,]
        coralprox$var_residu <- realprox$var_residu[cpos,]
        coralprox$datasource <- realprox$datasource[cpos]
      }
      if (length(tpos) > 0) {
        # run duplicate removal on trees only
        treeprox <- realprox
        treeprox$data <- realprox$data[,tpos]
        treeprox$lon <- realprox$lon[tpos]
        treeprox$lat <- realprox$lat[tpos]
        treeprox$archivetype <- realprox$archivetype[tpos]
        treeprox$elevation <- realprox$elevation[tpos]
        treeprox$mr <- realprox$mr[tpos,]
        treeprox$var_residu <- realprox$var_residu[tpos,]
        treeprox$datasource <- realprox$datasource[tpos]
        tmp_prox <- treeprox
        tmp_vec <- vector()
        for (i in 1:ncol(tmp_prox$data)) {
          if (i %% 1000 == 0) {
            print(paste('Remove duplicate proxies, i.e. best_prox_per_grid. Number:', i))
          }
          pos <- which(newgrid$lon>(tmp_prox$lon[i]-(2*res))&newgrid$lon<(tmp_prox$lon[i]+(2*res))&
                         newgrid$lat>(tmp_prox$lat[i]-(2*res))&newgrid$lat<(tmp_prox$lat[i]+(2*res))) 
          dist <- compute_dist(newgrid$lon[pos], newgrid$lat[pos], tmp_prox$lon[i], tmp_prox$lat[i])
          tmp_vec[i] <- pos[which.min(dist)]
        }
        #n_prox_per_grid <- apply(tmp_array, 2, sum)
        #Hproxy <- compute_H(tmp_prox, newgrid, threshold=(700*res))  
        #n_prox_per_grid <- apply(Hproxy, 2, sum)
        removepos=NULL
        tiv=which(floor(treeprox$time)==cyr)
        for (i in unique(tmp_vec)) {
          pos <- which(tmp_vec==i)
          if (length(pos) > 1) {
            #for (i in 1:ncol(Hproxy)){
            #  if (n_prox_per_grid[i] > 1) {
            #    pos <- which(Hproxy[,i]>0) # find prox in same grid box
            posnotna <- !is.na(treeprox$data[tiv,pos])
            if (treeprox$lat[pos][1] > 0) {
              tcol=2
            } else {
              tcol=1
            }
            pos2 <- pos[posnotna]
            pos_t <- rep(NA,length=length(pos2))
            pos_n <- rep(NA,length=length(pos2))
            # split into rather t and p sensitive
            for (j in 1:length(pos2)) {
              pos_t[j] <- (length(which(!is.na(treeprox$mr[pos2[j],2:7])))>=
                             length(which(!is.na(treeprox$mr[pos2[j],8:13]))))
              pos_n[j] <- (length(which(!is.na(treeprox$mr[pos2[j],2:7])))<
                             length(which(!is.na(treeprox$mr[pos2[j],8:13]))))
            }
            if (all(is.na(treeprox$var_residu[pos2,]))) {
              removepos <- c(removepos,pos)
            } else if (any(pos_t) & any(pos_n)) {
              removepos <- c(removepos,pos[-c(which(treeprox$var_residu[pos2,tcol]==
                                                      min(treeprox$var_residu[pos2[pos_t],tcol],na.rm=T)),
                                              which(treeprox$var_residu[pos,tcol]==
                                                      min(treeprox$var_residu[pos[pos_n],tcol],na.rm=T)))])
            } else if (any(pos_t) & !any(pos_n)) {
              removepos <- c(removepos,pos[-which(treeprox$var_residu[pos,tcol]==
                                                    min(treeprox$var_residu[pos2[pos_t],tcol],na.rm=T))])
            } else if (!any(pos_n) & any(pos_n)) {
              removepos <- c(removepos,pos[-which(treeprox$var_residu[pos,tcol]==
                                                    min(treeprox$var_residu[pos2[pos_n],tcol],na.rm=T))])
            } 
          }
        }
        # remove stations, will not be shown on station map 
        # treeprox$lon <- treeprox$lon[-removepos] 
        # treeprox$lat <- treeprox$lat[-removepos] 
        # treeprox$data <- treeprox$data[,-removepos] 
        # treeprox$names <- treeprox$names[-removepos]
        # treeprox$mr <- treeprox$mr[-removepos,]
        # treeprox$var_residu <- treeprox$var_residu[-removepos] 
        # set NA to keep station map symbol instead of removing
        treeprox$data[tiv,removepos] <- NA
        treeprox$mr[removepos,] <- NA
        treeprox$var_residu[removepos,] <- NA
      } # end if any tree data exists
      if ((length(cpos) > 0) & (length(tpos) > 0)) {
        # merge trees and corals to realprox again
        realprox$data <- cbind(treeprox$data,coralprox$data)
        realprox$lon <- c(treeprox$lon,coralprox$lon)
        realprox$lat <- c(treeprox$lat,coralprox$lat)
        realprox$archivetype <- c(treeprox$archivetype,coralprox$archivetype)
        realprox$elevation <- c(treeprox$elevation,coralprox$elevation)
        realprox$mr <- rbind(treeprox$mr,coralprox$mr)
        realprox$var_residu <- rbind(treeprox$var_residu,coralprox$var_residu)
        realprox$datasource <- c(treeprox$datasource,coralprox$datasource)
      } else if (length(tpos) > 0) {
        realprox <- treeprox
      } else if (length(cpos) > 0) {
        realprox <- coralprox
      }
    }
    
    # 3.3 Convert to 2 season per year or leave annual average for DAPS pseudoproxy exp.
    # no scaling because regresion takes care of it
      
    if (pseudo_prox) { # i.e. just only value representing the annual mean and not the growing season
      realprox.allts <- realprox
      realprox.allts$data <- t(realprox.allts$data) 
      ti=which(floor(realprox.allts$time)==cyr)
      realprox$time <- realprox.allts$time[ti]
      realprox$data <- realprox.allts$data[,ti,drop=F]
    } else {
      # treat corals with seasonal and trees with annual values separately
      #coral_pos <- which(realprox$archivetype=='coral')
      # should we have them completely separate?
      # how does pages data files look with combined corals and trees, trees T; and P, corals only T?
      
      # the function takes realprox and output is 2 season conversion and realprox.allts(which is needed below)
      listoftwo <-convert_to_2_seasons(realprox,source="proxy")
      realprox.allts <- listoftwo$x.allts
      realprox <- listoftwo$x
    }
    
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
      if (pseudo_prox) {
        realprox.anom$data <- (realprox$data - apply(realprox.tmp$data,1,mean,na.rm=T))
      } else {
        realprox.anom$data <- (realprox$data - (apply(array(realprox.tmp$data,
                               c(nrow(realprox.tmp$data), 2, ncol(realprox.tmp$data)/2)), 
                               1:2, mean,na.rm=T)))
      }
      realprox <- realprox.anom
    }
    
    # 3.4.5 Sort proxies by residuals 
    if (best2worst) {
      pos <- order(rowMeans(realprox$var_residu,na.rm=T))
    } else {
      pos <- rev(order(rowMeans(realprox$var_residu,na.rm=T)))
    }
    #if (pseudo_prox) {
    realprox$data <- realprox$data[pos,,drop=F]
    #} else {
    #  realprox$data <- realprox$data[pos,]
    #}
    realprox$lon <- realprox$lon[pos]
    realprox$lat <- realprox$lat[pos]
    realprox$archivetype <- NULL
    realprox$elevation <- NULL
    realprox$parameter <- NULL
    realprox$mr <- realprox$mr[pos,]
    realprox$names <- realprox$names[pos]
    realprox$var_residu <- realprox$var_residu[pos,]
    realprox$error <- sqrt(realprox$var_residu)
    realprox$sour <- rep('prox',length(realprox$lon))
    realprox.allts$sour <- rep('prox',length(realprox$lon))
    
    # 3.5 Create a list named proxies
    # real trw proxy multiple regression approach
    if (!instrumental) {
      proxies<-list(data=realprox$data, lon=realprox$lon,
                    lat=realprox$lat, names=realprox$names, sour=realprox$sour,
                    height=realprox$elevation, time=realprox$time,
                    mr=realprox$mr, error=realprox$error,
                    numavg=rep(1,length(realprox$lon)))
    }
    
    # 3.6 Set real_proxies to FALSE if there is no data  
    if (!pseudo_prox) {  
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
  ##########################################################################################
  
  if (docum) {
    #print('prepare documentary data')
    # 4.1 Loading documentary data: angie 2019 contains data from v1 + new monthly series
    load(paste0(dataextdir,'assimil_data/rdata_files_v2/temp2_docu_monthly_angie2017_2019.Rdata'))
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
      doc_t_mon$sour <- rep('doc',length(doc_t_mon$lon))
      doc_t_mon$error <- matrix(docu_err,nrow=length(doc_t_mon$lon),ncol=2)
      docall <- doc_t_mon
    } else {
      stop("docum section currently only coded for sixmonstatevector!")
    } # end of 6monstatevector
    
    # 4.5 Combine assimilation data into variable named "proxies"
    if ((!instrumental) & (!real_proxies)) {
      proxies <- doc_t_mon
    } else {
      docall <- doc_t_mon
    }
    
    if (!instrumental & real_proxies ) {
      tmpnum2 <- rep(1,length(realprox$lon)) #realprox.numavg
      tmpmr <- matrix(NA,nrow=length(docall$lon),ncol=ncol(proxies$mr))
      #tmpres <- rep(NA,length(proxies$var_residu))
      tmpelev <- rep(NA,length(proxies$height))
      tmpnum3 <- rep(1,length(docall$lon))
      proxies<-list(data=rbind(realprox$data,docall$data), lon=c(realprox$lon,docall$lon),
                    lat=c(realprox$lat,docall$lat), names=c(realprox$names,docall$names), 
                    sour=c(realprox$sour,docall$sour), 
                    height=c(realprox$height,tmpelev), time=realprox$time,
                    mr=rbind(realprox$mr,tmpmr), error=rbind(realprox$error,docall$error),
                    #var_residu=rbind(realprox$var_residu,tmpres),
                    numavg=c(tmpnum2,tmpnum3))
    }
  }  # end of docu
  
  
  ####################################################################
  # 5. Instrumental Part
  # 5.1 Loading the files
  # 5.2 Qualtiy check of the data
  # 5.3 Calculating the anomalies and transforming them to Oct-Sept time period
  # 5.4 Treat multiple data in the same grid box
  # 5.5 Combining all type of instrumental data
  # 5.6 Mask other data if there is instrumental data in same grid box at that time
  # 5.8 Add data source
  # 5.9 Convert to 2 season per year
  # 5.10 Combine assimilation data into variable named "proxies"
  ####################################################################
  
  if (instrumental){
    # 5.1 Loading the files
    if (yuri_temp) {
      # monthly temp collection from yuri, histalp included
      load(paste0(dataextdir,'assimil_data/rdata_files_v2/temp2_instr_monthly_yuri2017.Rdata')) 
      t2017 <- t
      load(paste0(dataextdir,'assimil_data/rdata_files_v2/temp2_instr_monthly_yuri2019.Rdata')) 
      t2019 <- temp2
      tpos <- t2019$time %in% t2017$time
      inst_t <- list()
      inst_t$data <- cbind(t2017$data,t2019$data[tpos,])
      inst_t$lon <- c(t2017$lon,t2019$lon)
      inst_t$lat <- c(t2017$lat,t2019$lat)
      inst_t$names <- c(t2017$names,t2019$names)
      inst_t$height <- c(t2017$height,t2019$height)
      inst_t$time <- t2017$time
    }
    if (yuri_slp) {
      # monthly slp collection from yuri, histalp included
      load(paste0(dataextdir,'assimil_data/rdata_files_v2/slp_instr_monthly_yuri2017.Rdata')) 
      slp2017 <- slp
      load(paste0(dataextdir,'assimil_data/rdata_files_v2/slp_instr_monthly_yuri2019.Rdata')) 
      slp2019 <- slp
      slppos <- slp2019$time %in% slp2017$time
      inst_slp <- list()
      inst_slp$data <- cbind(slp2017$data,slp2019$data[slppos,])
      inst_slp$lon <- c(slp2017$lon,slp2019$lon)
      inst_slp$lat <- c(slp2017$lat,slp2019$lat)
      inst_slp$names <- c(slp2017$names,slp2019$names)
      inst_slp$height <- c(slp2017$height,slp2019$height)
      inst_slp$time <- slp2017$time
    }
    if (yuri_prec) {
      # monthly prec collection from yuri, histalp included
      load(paste0(dataextdir,'assimil_data/rdata_files_v2/prec_instr_monthly_yuri2017.Rdata')) 
      prec2017 <- r
      load(paste0(dataextdir,'assimil_data/rdata_files_v2/prec_instr_monthly_yuri2019.Rdata')) 
      prec2019 <- prec
      precpos <- prec2019$time %in% prec2017$time
      inst_prec <- list()
      inst_prec$data <- cbind(prec2017$data,prec2019$data[precpos,])
      inst_prec$lon <- c(prec2017$lon,prec2019$lon)
      inst_prec$lat <- c(prec2017$lat,prec2019$lat)
      inst_prec$names <- c(prec2017$name,prec2019$names)
      inst_prec$height <- c(prec2017$height,prec2019$height)
      inst_prec$time <- prec2017$time
    }
    if (ghcn_temp) {
      if (isti_instead_ghcn) {
        load(paste0(dataextdir,"assimil_data/rdata_files_v2/isti_freeze1880_1600-2005.RData"))
        print("ACHTUNG: ISTI instr. temp. data will be assimilated NOT GHCN")
      } else {
        load(paste0(dataextdir,"assimil_data/rdata_files_v2/ghcn_temp_freeze1880_1602-2004.Rdata"))
      }
    }
    if (ghcn_prec) {
      load(paste0(dataextdir,"assimil_data/rdata_files_v2/ghcn_precip_freeze1880_1602-2004.Rdata"))
      #load(paste0("/scratch3/veronika/reuse/assimil_data/ghcn/ghcn_precip_1600-2005.Rdata"))
      #load("/scratch3/veronika/reuse/assimil_data/ghcn_d/ghcn_precip.RData")
      #load("/scratch3/veronika/reuse/assimil_data/ghnc_d_better_coverage/ghcn_precip.RData")
    }
    if (ghcn_wday) {
      # load("/scratch3/brugnara/ekf400/ghcn_wetdays.RData")
      # load("/scratch3/veronika/reuse/assimil_data/ghcn_d/ghcn_wetdays.RData")
      load("/scratch3/veronika/reuse/assimil_data/ghnc_d_better_coverage/ghcn_wetdays.RData")
    }
    
    # 5.2 Qualtiy check of the data
    if (check_assimdata) {
      varlist <- c("inst_t","inst_slp","inst_prec","ghcn","ghcn_precip","ghcn_wetdays")[c(yuri_temp,yuri_slp,yuri_prec,ghcn_temp,ghcn_prec,ghcn_wday)]
      # if (ghcn_prec) { 
      #   if (ghcn_temp & yuri_slp) {
      #     varlist <- c("inst_t","inst_slp","ghcn","ghcn_precip")
      #   } else {
      #     varlist <- c("ghcn_precip") # assimilating only precip
      #   }
      # } else if (ghcn_wday) {
      #   if (ghcn_temp & yuri_slp) {
      #     varlist <- c("inst_t","inst_slp","ghcn","ghcn_wetdays")
      #   } else {
      #     varlist <- c("ghcn_wetdays") # assimilating only wetdays
      #   }
      # } else if (ghcn_temp & yuri_slp) {
      #   varlist <- c("inst_t","inst_slp","ghcn")
      # }
      for (varname in varlist) {
        var<-screendistance(echam.sd,varname)
        if (varname=="inst_t") {inst_t$data<- var$data}
        if (varname=="inst_slp") {inst_slp$data<- var$data}
        if (varname=="inst_prec") {inst_prec$data<- var$data}
        if (varname=="ghcn") {ghcn$data <- var$data}
        if (varname=="ghcn_precip") {ghcn_precip$data <- var$data}
        if (varname=="ghcn_wetdays") {ghcn_wetdays$data <- var$data}
        var<-screenstd(get(varname),cyr,source="inst",sdlim=5) # screen 5SD for all temp, slp, precip???
        if (varname=="inst_t") {inst_t$data <- var$data}
        if (varname=="inst_slp") {inst_slp$data <- var$data}
        if (varname=="inst_prec") {inst_prec$data<- var$data}
        if (varname=="ghcn") {ghcn$data<- var$data}
        if (varname=="ghcn_precip") {ghcn_precip$data <- var$data}
        if (varname=="ghcn_wetdays") {ghcn_wetdays$data <- var$data}
      }
    }
    
    # 5.3 Calculating the anomalies and transforming them to Oct-Sept time period
    if (ghcn_temp) {
      if (dim(ghcn$data)[2]>0) {
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
    }
    if (ghcn_prec) {
      if(dim(ghcn_precip$data)[2]>0) {
        if (anomaly_assim & sixmonstatevector){
          if (gauss_ana) {
            # select the first obs per grid
            if (first_prox_per_grid) {
              res <- firstproxres
              newgrid <- echam
              newgrid$data <- NULL
              newgrid$ensmean <- NULL
              newgrid$lon <- echam$lon[seq(1,length(echam$lon),res)]
              newgrid$lat <- echam$lat[seq(1,length(echam$lat),res)]
              tmp_precip <-ghcn_precip
              tmp_precip$data <- t(ghcn_precip$data)
              Hproxy <- compute_H(tmp_precip, newgrid, threshold=(700*res))  
              p.i <- as.logical(apply(Hproxy, 1, sum)) & !duplicated(Hproxy, margin=1)
              ghcn_precip$lon <-ghcn_precip$lon[p.i]
              ghcn_precip$lat <- ghcn_precip$lat[p.i]
              ghcn_precip$data <- ghcn_precip$data[,p.i]
              ghcn_precip$names <- ghcn_precip$names[p.i] 
            }
            
            # Gaussian transformation is done on the whole period for each month
            for(GAyr in 1601:2005) {
              # calculate the climatology
              ghcn_prec.clim=ghcn_precip
              ti=which(floor(ghcn_precip$time)>=(GAyr-36) & floor(ghcn_precip$time)<=(GAyr+35))
              sts=ti[1]
              ets=ti[length(ti)]
              ghcn_prec.clim$data=t(ghcn_precip$data[sts:ets,])
              ghcn_prec.clim$time=ghcn_precip$time[sts:ets]
              # transform it to Oct-Sept  years
              if (GAyr < 1636) {
                y_start = agrep(paste0(1600,'.792'),as.character(ghcn_prec.clim$time))
              } else {
                y_start = agrep(paste0(GAyr-36,'.792'),as.character(ghcn_prec.clim$time))
              }
              if (GAyr > 1970) {
                y_end =  grep(paste0(2005,'.708'),as.character(ghcn_prec.clim$time))
              } else {
                y_end = grep(paste0(GAyr+35,'.708'),as.character(ghcn_prec.clim$time))
              }
              ghcn_prec.clim$data = ghcn_prec.clim$data[,y_start:y_end]
              ghcn_prec.clim$time = ghcn_prec.clim$time[y_start:y_end]
              ghcn_prec.clim$data = apply(array(ghcn_prec.clim$data, c(nrow(ghcn_prec.clim$data), 12, ncol(ghcn_prec.clim$data)/12)), 1:2, mean,na.rm=T)
              
              # cut out the GAyr-1 and GAyr time window
              ghcn_prec.cyr = ghcn_precip
              ti=which(floor(ghcn_precip$time)>=GAyr-1 & floor(ghcn_precip$time)<GAyr+1)
              sts=ti[1]
              ets=ti[length(ti)]
              ghcn_prec.cyr$data=t(ghcn_precip$data[sts:ets, ])
              ghcn_prec.cyr$time=ghcn_precip$time[sts:ets]
              # transform it to Oct-Sept  years
              y_start = agrep(paste0(GAyr-1,'.792'),as.character(ghcn_prec.cyr$time))
              y_end = grep(paste0(GAyr,'.708'),as.character(ghcn_prec.cyr$time))
              ghcn_prec.cyr$data = ghcn_prec.cyr$data[,y_start:y_end]
              ghcn_prec.cyr$time = ghcn_prec.cyr$time[y_start:y_end]
              # calculate the anomaly for GAyr and combine them
              if (GAyr == 1601)  {
                ghcn_prec.anom <-ghcn_prec.cyr
                if (precip_ratio) {
                  ghcn_prec.anom$data <- t(ghcn_prec.cyr$data / ghcn_prec.clim$data)
                } else {
                  ghcn_prec.anom$data <- t(ghcn_prec.cyr$data - ghcn_prec.clim$data)
                }
                ghcn_p <- ghcn_prec.anom
              } else {
                ghcn_prec.anom <-ghcn_prec.cyr
                if (precip_ratio) {
                  ghcn_prec.anom$data <- t(ghcn_prec.cyr$data / ghcn_prec.clim$data)
                } else {
                  ghcn_prec.anom$data <- t(ghcn_prec.cyr$data - ghcn_prec.clim$data)
                }
                ghcn_pp <- ghcn_prec.anom
                ghcn_p$data = rbind(ghcn_p$data, ghcn_pp$data)
                ghcn_p$time = c(ghcn_p$time, ghcn_pp$time)
              }
            }
            # Gaussian anamorphosis
            ghcn_p$GA_data = ghcn_p$data
            for (p.stat in 1:dim(ghcn_p$data)[2]) {
              for (mm in 1:12) {
                m = seq(mm,dim(ghcn_p$data)[1],12)
                prec <- ghcn_p$data[m,p.stat]
                prec[is.nan(prec)] = NA
                if (!all(is.na(prec))) {
                  pdb <- db.create(prec,ndim=1,autoname=F)
                  pdb.gaus <- anam.fit(pdb,name="prec",type="gaus",draw=F)
                  pdb.trans <- anam.z2y(pdb,'prec',anam=pdb.gaus)
                  prec.trans = pdb.trans@items$Gaussian.prec
                  anam.write(pdb.gaus, filename=paste0("/scratch3/veronika/reuse/assimil_data/gauss_ana/GA_obs_pdiff/obs_", p.stat,"_",mm,".ascii"), 
                             verbose = 0, flag.calcul = 1)
                  if (check_norm) {
                    x=shapiro.test(prec.trans)
                    if (x$p.value >= 0.05) {
                      ghcn_p$GA_data[m,p.stat] = prec.trans
                    } else {
                      ghcn_p$GA_data[m,p.stat] = NA
                    }
                  } else {
                    ghcn_p$GA_data[m,p.stat] = prec.trans
                  }
                  # print(prec.trans,flag.stats=TRUE,names="Gaussian.prec")
                  # hist(prec.trans)
                  # qqnorm(prec.trans)
                  # qqline(prec.trans)
                }
              }
            }
            # keep the current year that it is assimilated
            y_start = agrep(paste0(cyr-1,'.792'),as.character(ghcn_p$time))
            y_end = grep(paste0(cyr,'.708'),as.character(ghcn_p$time))
            ghcn_p$data = ghcn_p$GA_data[y_start:y_end,]
            ghcn_p$err = ghcn_p$data * 0.1
            ghcn_p$time = ghcn_p$time[y_start:y_end]
          } else {
            # calculate the climatology
            ghcn_prec.clim=ghcn_precip
            ti=which(floor(ghcn_precip$time)>=(cyr-36) & floor(ghcn_precip$time)<=(cyr+35))
            sts=ti[1]
            ets=ti[length(ti)]
            ghcn_prec.clim$data=t(ghcn_precip$data[sts:ets,])
            ghcn_prec.clim$time=ghcn_precip$time[sts:ets]
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
            ghcn_prec.clim$data = apply(array(ghcn_prec.clim$data, c(nrow(ghcn_prec.clim$data), 
                                    12, ncol(ghcn_prec.clim$data)/12)), 1:2, mean,na.rm=T)
            # cut out the cyr-1 and cyr time window
            ghcn_prec.cyr = ghcn_precip
            ti=which(floor(ghcn_precip$time)>=cyr-1 & floor(ghcn_precip$time)<cyr+1)
            sts=ti[1]
            ets=ti[length(ti)]
            ghcn_prec.cyr$data=t(ghcn_precip$data[sts:ets, ])
            ghcn_prec.cyr$time=ghcn_precip$time[sts:ets]
            # transform it to Oct-Sept  years
            y_start = agrep(paste0(cyr-1,'.792'),as.character(ghcn_prec.cyr$time))
            y_end = grep(paste0(cyr,'.708'),as.character(ghcn_prec.cyr$time))
            ghcn_prec.cyr$data = ghcn_prec.cyr$data[,y_start:y_end]
            ghcn_prec.cyr$time = ghcn_prec.cyr$time[y_start:y_end]
            # define precip error as %
            ghcn_prec.cyr$err = ghcn_prec.cyr$data * ghcn_p_err # error is 30% or min 10 mm
            ghcn_prec.cyr$err[ghcn_prec.cyr$err < ghcn_p_min] = ghcn_p_min 
            # calculate the anomaly for cyr
            if (precip_ratio) {
              ghcn_prec.anom <- ghcn_prec.cyr
              ghcn_prec.anom$data <- t(ghcn_prec.cyr$data / ghcn_prec.clim$data)
              ghcn_p <- ghcn_prec.anom
            } else {
              ghcn_prec.anom <- ghcn_prec.cyr
              ghcn_prec.anom$data <- t(ghcn_prec.anom$data - ghcn_prec.clim$data)
              ghcn_prec.anom$err <- t(ghcn_prec.cyr$err)
              ghcn_p <- ghcn_prec.anom
            }
          }
        } else {
          stop("instr section currently only coded for anom_assim and sixmonstatevector!")
        }
      }
    }
    if (ghcn_wday) {
      if (dim(ghcn_wetdays$data)[2]>0) {
        if (anomaly_assim & sixmonstatevector){
          # calculate the climatology
          ghcn_wetdays.clim<-calculate_climatology(ghcn_wetdays,cyr,36,35,source="inst")
          ghcn_wetdays.clim$data<-t(ghcn_wetdays.clim$data)
          ghcn_wetdays.clim$data = apply(array(ghcn_wetdays.clim$data, c(nrow(ghcn_wetdays.clim$data), 12, ncol(ghcn_wetdays.clim$data)/12)), 
                                         1:2, mean,na.rm=T)
          
          # cut out the cyr-1 and cyr time window
          ghcn_wetdays.cyr<-calculate_climatology(ghcn_wetdays,cyr,1,0,source="inst")
          ghcn_wetdays.cyr$data<-t(ghcn_wetdays.cyr$data)
          
          # define wetdays error as 10%
          ghcn_wetdays.cyr$err = ghcn_wetdays.cyr$data * 0.1
          
          # calculate the anomaly for cyr
          ghcn_wetdays.anom <-ghcn_wetdays.cyr
          ghcn_wetdays.anom$data <- t(ghcn_wetdays.anom$data - ghcn_wetdays.clim$data)
          ghcn_wetdays.anom$err <- t(ghcn_wetdays.cyr$err)
          ghcn_w <- ghcn_wetdays.anom
        } else {
          stop("instr section currently only coded for anom_assim and sixmonstatevector!")
        }
      }
    }
    if (yuri_slp) {
      if (dim(inst_slp$data)[2]>0) {
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
    }
    if (yuri_prec) {
      if (dim(inst_prec$data)[2]>0) {
        if (anomaly_assim & sixmonstatevector) {
          # calculate the climatology
          inst_prec.clim<-calculate_climatology(inst_prec,cyr,36,35,source="inst")
          inst_prec.clim$data<-t(inst_prec.clim$data)
          inst_prec.clim$data = apply(array(inst_prec.clim$data, c(nrow(inst_prec.clim$data), 
                                  12, ncol(inst_prec.clim$data)/12)), 1:2, mean,na.rm=T)
          # cut out the cyr-1 and cyr time window
          inst_prec.cyr<-calculate_climatology(inst_prec,cyr,1,0,source="inst")
          inst_prec.cyr$data<-t(inst_prec.cyr$data)
          # define precip error as %
          inst_prec.cyr$err = inst_prec.cyr$data * ghcn_p_err # error is 30% or min 10 mm
          inst_prec.cyr$err[inst_prec.cyr$err < ghcn_p_min] = ghcn_p_min 
          # calculate the anomaly for cyr
          if (precip_ratio) {
            inst_prec.anom <- inst_prec.cyr
            inst_prec.anom$data <- t(inst_prec.cyr$data / inst_prec.clim$data)
            inst_p <- inst_prec.anom
          } else {
            inst_prec.anom <- inst_prec.cyr
            inst_prec.anom$data <- t(inst_prec.anom$data - inst_prec.clim$data)
            inst_prec.anom$err <- t(inst_prec.cyr$err)
            inst_p <- inst_prec.anom
          }
        } else {
          stop("instr section currently only coded for anom_assim and sixmonstatevector!")
        }
      }
    }
    if (yuri_temp) {
      if (dim(inst_t$data)[2]>0) {
        if (anomaly_assim & sixmonstatevector){
          # calculate the climatology
          inst_t.clim<-calculate_climatology(inst_t,cyr,36,35,source="inst")
          inst_t.clim$data<-t(inst_t.clim$data)
          inst_t.clim$data = apply(array(inst_t.clim$data, c(nrow(inst_t.clim$data), 
                               12, ncol(inst_t.clim$data)/12)), 1:2, mean,na.rm=T)
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
    }
    # if (ghcn_prec) { 
    #   if (ghcn_temp & yuri_slp) {
    #     if ((dim(inst_slp$data)[2]==0) & (dim(inst_t$data)[2]==0) & (dim(ghcn$data)[2]==0) &
    #         (dim(ghcn_precip$data)[2]==0)) {instrumental=F}
    #   } else {
    #     if ((dim(ghcn_precip$data)[2]==0)) {instrumental=F} # assimilating only precip
    #   }
    # } else if (ghcn_wday) {
    #   if (ghcn_temp & yuri_slp) {
    #     if ((dim(inst_slp$data)[2]==0) & (dim(inst_t$data)[2]==0) & (dim(ghcn$data)[2]==0) &
    #         (dim(ghcn_w$data)[2]==0)) {instrumental=F}
    #   } else {
    #     if (dim(ghcn_w$data)[2]==0) {instrumental=F}
    #   }
    # } else if (ghcn_temp & yuri_slp) {
    #   if ((dim(inst_slp$data)[2]==0) & (dim(inst_t$data)[2]==0) & (dim(ghcn$data)[2]==0)) {
    #     instrumental=F}
    # }
  
    
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
    if ((ghcn_prec) & (yuri_prec)) {
      inst_p<-list(data=cbind(ghcn_p$data,inst_p$data), lon=c(ghcn_p$lon,inst_p$lon),
                   lat=c(ghcn_p$lat,inst_p$lat), names=c(ghcn_p$names,inst_p$names),
                   height=c(ghcn_p$height,inst_p$height), time=ghcn_p$time,
                   err=cbind(ghcn_p$err,inst_p$err))
    }
    if ((ghcn_prec) & (!yuri_prec)) {
      inst_p <- ghcn_p
    }
    if (ghcn_wday) {
      inst_w = ghcn_w
    }
    
    # 5.4.2 First proxy per echam grid box
    if (first_inst_per_grid) {
      res <- firstinstres
      newgrid <- echam
      newgrid$data <- NULL
      newgrid$ensmean <- NULL
      newgrid$lon <- echam$lon[seq(1,length(echam$lon),res)]
      newgrid$lat <- echam$lat[seq(1,length(echam$lat),res)]
      
      if (ghcn_temp) {
        tmp_t <- inst_t
        tmp_t$data <- t(inst_t$data)
        Hproxy <- compute_H(tmp_t, newgrid, threshold=(700*res))  
        p.i <- as.logical(apply(Hproxy, 1, sum)) & !duplicated(Hproxy, margin=1)
        inst_t$lon <- inst_t$lon[p.i]
        inst_t$lat <- inst_t$lat[p.i]
        inst_t$data <- inst_t$data[,p.i]
        inst_t$names <- inst_t$names[p.i]
      }
      
      if (!gauss_ana) {  # if GA is T, this selection happens earlier 
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
        print("DON'T USE FIRST_INST_PER_GRID IF YOU WANT TO INCLUDE REAL PROXY DATA AND 
              INSTRUMENTALS AT THE SAME TIME")
      }  
    }
    
    # 5.4.3 Average more than one proxy per echam grid box 
    if (avg_obs_per_grid) {
      # average data in same grid box and count number of avg. series as error estimate
      # separate temp, precip, slp
      # makes no sense for realprox: proxies would need to be calibrated before building 
      # regession model
      if (ghcn_prec) { 
        if (ghcn_temp & yuri_slp) {
          varlist <- c("inst_t","inst_p","inst_slp")
        } else {
          varlist <- c("inst_p") # assimilating only precip
        }
      } else if (ghcn_wday) {
        if (ghcn_temp & yuri_slp) {
          varlist <- c("inst_t","inst_slp","inst_w") 
        } else {
          varlist = "inst_w"
        }
      } else if (ghcn_temp & yuri_slp) {
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
        stat.err=rep(NA,12) 
        
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
            if (varname == "inst_t") {stat.err = rbind(stat.err,rep(inst_t_err,12))} # define error here
            if (varname == "inst_slp") {stat.err = rbind(stat.err,rep(inst_slp_err,12))} # define error here
            if (varname == "inst_p") {stat.err = rbind(stat.err,apply(stat$err[,mask],1,mean,na.rm=T))} 
            if (varname == "inst_w") {stat.err = rbind(stat.err,rep(ghcn_w_err,12))} # define error here
          } else {
            stat.avg=rbind(stat.avg,stat$data[,mask])
            stat.numavg=rbind(stat.numavg,1)
            stat.lon=rbind(stat.lon,stat$lon[mask])
            stat.lat=rbind(stat.lat,stat$lat[mask])
            stat.names=rbind(stat.names,stat$names[mask])
            if (varname == "inst_t") {stat.err = rbind(stat.err,rep(inst_t_err,12))} # define error here
            if (varname == "inst_slp") {stat.err = rbind(stat.err,rep(inst_slp_err,12))} # define error here
            if (varname == "inst_p") {stat.err = rbind(stat.err,stat$err[,mask])} 
            if (varname == "inst_w") {stat.err = rbind(stat.err,rep(ghcn_w_err,12))} # define error here
          }
          stat.avg[stat.avg=='NaN']=NA
          stat.err[is.na(stat.avg)]=NA
        }
        stat$lon <- stat.lon[2:length(stat.lon)]
        stat$lat <- stat.lat[2:length(stat.lat)]
        stat$data <- t(stat.avg[2:dim(stat.avg)[1],])
        stat$numavg <- stat.numavg[2:length(stat.numavg)]
        stat$names <- stat.names[2:length(stat.names)]
        stat$err <- t(stat.err[2:dim(stat.err)[1],])
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
        if (varname == "inst_w") {
          inst_w = stat
        }
      }
    } 
    
    # 5.5 Combining all type of instrumental data
    if (ghcn_prec & !ghcn_wday) {
      if (ghcn_temp & yuri_slp) { # since the order matters: combine them as: temp, slp, precip
        inst<-list(data=t(cbind(inst_t$data,inst_slp$data,inst_p$data)),
                   error=t(cbind(inst_t$err,inst_slp$err,inst_p$err)), # new define error
                   lon=c(inst_t$lon,inst_slp$lon,inst_p$lon),
                   lat=c(inst_t$lat,inst_slp$lat,inst_p$lat),
                   numavg=c(inst_t$numavg,inst_slp$numavg,inst_p$numavg),
                   names=c(inst_t$names,inst_slp$names,inst_p$names),
                   height=c(inst_t$height,inst_slp$height,inst_p$height), time=inst_t$time)
      } else { # assimilating only precip
      inst<-list(data=t(cbind(inst_p$data)),error=t(cbind(inst_p$err)), 
                             lon=c(inst_p$lon),
                             lat=c(inst_p$lat),
                             numavg=c(inst_p$numavg),
                             names=c(inst_p$names),
                             height=c(inst_p$height), time=inst_p$time)
      }
    } else if (ghcn_wday & !ghcn_prec) {
      if (ghcn_temp & yuri_slp) { # since the order matters: combine them as: temp, slp, wdays
        inst<-list(data=t(cbind(inst_t$data,inst_slp$data,inst_w$data)),
                   error=t(cbind(inst_t$err,inst_slp$err,inst_w$err)), # new define error
                   lon=c(inst_t$lon,inst_slp$lon,inst_w$lon),
                   lat=c(inst_t$lat,inst_slp$lat,inst_w$lat),
                   numavg=c(inst_t$numavg,inst_slp$numavg,inst_w$numavg),
                   names=c(inst_t$names,inst_slp$names,inst_w$names),
                   height=c(inst_t$height,inst_slp$height,inst_w$height), time=inst_t$time)
      } else { # assimilating only wetdays
        inst<-list(data=t(cbind(inst_w$data)), error=t(cbind(inst_w$err)), # new define error
                   lon=c(inst_w$lon),
                   lat=c(inst_w$lat),
                   numavg=c(inst_w$numavg),
                   names=c(inst_w$names),
                   height=c(inst_w$height), time=inst_w$time)
      }
    } else if (ghcn_temp & yuri_slp & !ghcn_prec) {
      inst<-list(data=t(cbind(inst_t$data,inst_slp$data)), lon=c(inst_t$lon,inst_slp$lon),
                 lat=c(inst_t$lat,inst_slp$lat), names=c(inst_t$names,inst_slp$names),
                 numavg=c(inst_t$numavg,inst_slp$numavg),
                 height=c(inst_t$height,inst_slp$height), time=inst_t$time)
    } else if (ghcn_temp & yuri_slp & ghcn_prec & ghcn_wday) {
      inst<-list(data=t(cbind(inst_t$data,inst_slp$data,inst_w$data,inst_p$data)),
                 error=t(cbind(inst_t$err,inst_slp$err,inst_w$err,inst_p$err)), # new define error
                 lon=c(inst_t$lon,inst_slp$lon,inst_w$lon,inst_p$lon),
                 lat=c(inst_t$lat,inst_slp$lat,inst_w$lat,inst_p$lat),
                 numavg=c(inst_t$numavg,inst_slp$numavg,inst_w$numavg,inst_p$numavg),
                 names=c(inst_t$names,inst_slp$names,inst_w$names,inst_p$names),
                 height=c(inst_t$height,inst_slp$height,inst_w$height,inst_p$height), time=inst_t$time)
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
            realprox$var_residu[i,] <- NA
            realprox$error[i,] <- NA
          }
        }
        pos <- apply(!is.na(realprox$data),1,any)
        realprox$data = realprox$data[pos,]
        realprox$lon <- realprox$lon[pos]
        realprox$lat <- realprox$lat[pos]
        realprox$mr <- realprox$mr[pos,]
        realprox$var_residu <- realprox$var_residu[pos,]
        realprox$error <- realprox$error[pos,]        
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
    
    # 5.8 Add data source/type information to list
    if (instrumental) {inst$sour <- rep('inst',length(inst$lon))}
    # if (docum) {
    #   docall$sour <- rep('doc',length(docall$lon))
    #   # docall.allts$sour <- rep('doc',length(docall$lon)) # doesn't exists anymore   
    # }
    # if (real_proxies) {
    #   realprox$sour <- rep('prox',length(realprox$lon))
    #   realprox.allts$sour <- rep('prox',length(realprox$lon))
    # }
    
    # 5.9 Convert to 2 season per year
    if (sixmonstatevector) { 
      # change instr array that only 2 ts instead of 12 monthly ts but 6 times as 
      # many variables, one for each month
      tmp1 <- array(inst$data,c(dim(inst$data)[1]*dim(inst$data)[2]))
      inst$data <- array(tmp1,c(dim(tmp1)[1]/(dim(inst$data)[2]/6),2))
      rm(tmp1)
      # transform the errors as well to 6month-state vector format
      tmp2 <- array(inst$error,c(dim(inst$error)[1]*dim(inst$error)[2]))
      inst$error <- array(tmp2,c(dim(tmp2)[1]/(dim(inst$error)[2]/6),2))
      rm(tmp2)
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
      if (avg_obs_per_grid) {
        inst$numavg <- rep(inst$numavg,6)
      }
    }
    
    # 5.10 Combine assimilation data into variable named "proxies"
    if (!real_proxies) {        
      proxies <- inst
    }
    
    if (real_proxies) {
      tmpmr <- matrix(NA,nrow=length(inst$lon),ncol=ncol(realprox$mr))
      tmpres <- matrix(NA,nrow=nrow(realprox$var_residu),ncol=ncol(realprox$var_residu))
      if (avg_obs_per_grid) {
        tmpnum1 <- inst$numavg
        tmpnum2 <- rep(1,length(realprox$lon))
      } else {
        tmpnum1 <- rep(1,length(inst$lon))
        tmpnum2 <- rep(1,length(realprox$lon)) #realprox.numavg
      }
      if (real_proxies & instrumental) {
        proxies<-list(data=rbind(inst$data,realprox$data), error=rbind(inst$err,realprox$error),
                      lon=c(inst$lon,realprox$lon), lat=c(inst$lat,realprox$lat), 
                      names=c(inst$names,realprox$names), 
                      sour=c(inst$sour,realprox$sour), 
                      height=c(inst$height,realprox$elevation), time=inst$time,
                      mr=rbind(tmpmr,realprox$mr), #var_residu=c(tmpres,realprox$var_residu),
                      numavg=c(tmpnum1,tmpnum2))
      }
    }
    
    if (docum) {
      tmpmr <- matrix(NA,nrow=length(docall$lon),ncol=ncol(proxies$mr))
      # tmpres <- rep(NA,length(proxies$var_residu))
      tmpelev <- rep(NA,length(proxies$height))
      tmpnum3 <- rep(1,length(docall$lon))
      proxies<-list(data=rbind(proxies$data,docall$data), error=rbind(proxies$error,docall$error),
                    lon=c(proxies$lon,docall$lon), lat=c(proxies$lat,docall$lat), 
                    names=c(proxies$names,docall$names), sour=c(proxies$sour,docall$sour), 
                    height=c(proxies$height,tmpelev), time=proxies$time,
                    mr=rbind(proxies$mr,tmpmr), #var_residu=c(realprox$var_residu,tmpres),
                    numavg=c(proxies$numavg,tmpnum3))
    }
  } # end "if (instrumental)"
  
  
  
  
  #########################################################################################
  # 6. All assimilation data
  #########################################################################################
  
  calibrate <- proxies
  if (is.null(dim(proxies$data))) {
    calibrate$data <- t(proxies$data)
    print(paste('number of proxies/observations:',length(calibrate$data)))
  } else {
    print(paste('number of proxies/observations:',dim(calibrate$data)[1]))  
  }
  print("calc time preparing proxies")
  print(proc.time() - ptm1)
  
  
  
  
  #########################################################################################
  # 7. Compute H (forward operator)
  #########################################################################################
  if ((real_proxies) & (!instrumental) & (!docum) & (!sixmonstatevector)) {
    if (pseudo_prox) {
      Hcal <- compute_Hi_Hredux_pseudoproxy(realprox, echam, realprox$mr, threshold=700)
      H.i <- array(NA,c(nrow(Hcal),(ncol(Hcal))/2))
      Hredux <- H.i
      H.i[1:nrow(Hcal),]<-Hcal[,seq(1,ncol(Hcal),2)] 
      Hredux[1:nrow(Hcal),]<-Hcal[,seq(2,ncol(Hcal),2)] 
      H.i[H.i==0] <- NA
      Hredux[Hredux==0] <- NA
    } else {
      # does this work with Hi and Hredox already? Guess not
      Hcal <- Matrix(compute_H_proxy(realprox, echam, realprox$mr, threshold=700), 
                     sparse=T)
    }
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
    if (instrumental & !docum & !real_proxies) {
      H.i <- array(NA,c(nrow(Hcal1),1))
    } else if (instrumental & docum & !real_proxies) {
      H.i <- array(NA,c((nrow(Hcal1)+nrow(Hcal3)),1))
    } else if (!instrumental & docum & real_proxies) {
      H.i <- array(NA,c((nrow(Hcal2)+nrow(Hcal3)),(ncol(Hcal2))/2))
    } else if (instrumental & !docum & real_proxies) {
      H.i <- array(NA,c((nrow(Hcal1)+nrow(Hcal2)),(ncol(Hcal2))/2))
    } else if (!instrumental & !docum & real_proxies) { 
      H.i <- array(NA,c(nrow(Hcal2),(ncol(Hcal2))/2))
    } else if (!instrumental & docum & !real_proxies) { 
      H.i <- array(NA,c(nrow(Hcal3),(ncol(Hcal3))/2))  
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
    if (!instrumental & !docum & real_proxies){ 
      H.i[1:nrow(Hcal2),]<-Hcal2[,seq(1,ncol(Hcal2),2)] 
      Hredux[1:nrow(Hcal2),]<-Hcal2[,seq(2,ncol(Hcal2),2)] 
    }
    if (!instrumental & docum & !real_proxies){ 
      H.i[1:nrow(Hcal3),]<-Hcal3[,seq(1,ncol(Hcal3),2)] 
      Hredux[1:nrow(Hcal3),]<-Hcal3[,seq(2,ncol(Hcal3),2)] 
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
    if (sixmonstatevector | pseudo_prox) {
      #if ((real_proxies) & ((instrumental) | (docum))) { 
      Rcal <- calibrate$error^2
       # Rcal <- c(temp2=0.9, precip=50, slp=10)[calibrate$names]
       #  #        if (avg_prox_per_grid) {Rcal <- Rcal*(1/calibrate$numavg)}
       # Rcal[calibrate$names=="prox"] <- realprox$var_residu 
       #  # previously used residuals/2 for 1. paper version to give proxies more weight
       #  # better delete "/2"
       #  # probably should have given instrumentals more error instead!
       #  Rcal[calibrate$sour=="doc"] <- 0.25 # equals 0.5 std. dev.
      #} else if ((real_proxies) & (!instrumental) & (!docum)) { 
      #  Rcal <- cbind(rep(NA,length(realprox$var_residu)),realprox$var_residu)
        # *2 to increase the possibly underestimated 
        # error and prohibit overcorrection if using thousands of proxy records
        #/ 10 # devision to make proxy error smaller and update larger
      Rcal[which(is.na(Rcal))] <- 0
      #} else if (((instrumental) | (docum)) & (!real_proxies)) { 

        # original: Rcal <- c(temp2=0.9, precip=50, slp=10)[calibrate$names]
        # Rcal[calibrate$sour=="doc"] <- 0.25 # equals 0.5 std. dev
      #  Rcal <- calibrate$error^2

        #      # set squared error R for precip measurements to 25% of data value
        #      # R <- abs(calibrate$data*0.25)
      #}
      
      if (gauss_ana) {
        GAecham = echam
        for (i in 1:2) {
          if (i == 1) {
            for (mm in c(10,11,12,1,2,3)) {
              for (p.stat in 4609:9216){
                pdb.gaus = anam.read(filename=paste0("/scratch3/veronika/reuse/assimil_data/gauss_ana/clim_70yr_pdiff/anam_",p.stat,"_",mm,".ascii"))
                if (mm>=10) {
                  k = mm-10
                } else {
                  k=mm+2
                }
                pos = k * (one_var_dim*numvar) + p.stat
                echam_mm = GAecham$data[pos,i,]
                pdb_mm <- db.create(echam_mm,ndim=1,autoname=F)
                prec.mm = anam.z2y(pdb_mm,'echam_mm',anam=pdb.gaus)
                GAecham$data[pos,i,] = prec.mm@items$Gaussian.echam_mm
              }
            }
            GAecham$ensmean[,i] = rowMeans(GAecham$data[,i,])
          }
          if (i == 2) {
            for (mm in 4:9) {
              for (p.stat in 4609:9216){
                pdb.gaus = anam.read(filename=paste0("/scratch3/veronika/reuse/assimil_data/gauss_ana/clim_70yr_pdiff/anam_",p.stat,"_",mm,".ascii"))
                k = mm-4
                pos = k * (one_var_dim*numvar) + p.stat
                echam_mm =  GAecham$data[pos,i,]
                pdb_mm <- db.create(echam_mm,ndim=1,autoname=F)
                prec.mm = anam.z2y(pdb_mm,'echam_mm',anam=pdb.gaus)
                GAecham$data[pos,i,] = prec.mm@items$Gaussian.echam_mm
              }
            } 
            GAecham$ensmean[,i] = rowMeans(GAecham$data[,i,])
          }
        }
        analysis <- GAecham
        # translate the mean of the rv N(0,1) to the mean of the original variable
        # analysis$data[analysis$names=="precip",1,] = GAecham$data[GAecham$names == "precip",1,] + echam$ensmean[echam$names == "precip", 1]
        # analysis$ensmean[analysis$names == "precip",1] = rowMeans(analysis$data[analysis$names=="precip",1,])
        # analysis$data[analysis$names=="precip",2,] = GAecham$data[GAecham$names == "precip",2,] + echam$ensmean[echam$names == "precip", 2]
        # analysis$ensmean[analysis$names == "precip",2] = rowMeans(analysis$data[analysis$names=="precip",2,])
      } else { #original
        analysis <- echam
      }
      # take anomalies
      if (gauss_ana) {
        analysis$data <- analysis$data - as.vector(analysis$ensmean)
      } else {
        analysis$data <- echam$data - as.vector(echam$ensmean)
      }
      if (cov_inflate){
        analysis$data= analysis$data * inflate_fac
      }
      if (covarclim>0 & covarclim<=100) { 
        ananomallts = echanomallts
        ananomallts$data <- echanomallts$data - as.vector(echanomallts$ensmean)
      }
      #nmonths <- 6
      ndim <- nrow(analysis$data)
      ntim <- ncol(analysis$data)
      nens <- dim(analysis$data)[3]
      #if (pseudo_prox) {
      #  nprox <- length(calibrate$data)  
      #} else {
      nprox <- nrow(calibrate$data)    
      #}
      #ndimold <- length(echam$lon) / nmonths
      #itime <- rep(1:nmonths, each=ndimold)
      if (landcorr) {
        corland_analysis <- landcorrected.anom
      }
      i=1 
      # compute loop over obs first to minimize repetition of computation (e.g. H and weights)
      for (j in 1:nprox) {
        if (j %% 100 == 0) {
          print(paste('Assimilating observation', j))
        }
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
                if (mixed_loc) { # changed for the B exps: first combining and than localizing
                  PH <- (analysis$data[,i,] %*% t(x2) / (nens - 1)) # %*% t(H)
                } else { 
                  if (gauss_ana) {
                    P <- analysis$data[,i,] %*% t(x2) / (nens - 1)
                    GA_P = P
                    # converting the P back to physical space
                    for (p.stat in 4609:(4608+4608)) {
                      if(i == 1) {
                        for (mm in c(10,11,12,1,2,3)) {
                          fname = paste0("/scratch3/veronika/reuse/assimil_data/gauss_ana/clim_70yr_pdiff/anam_",p.stat,"_",mm,".ascii")
                          if (file.exists(fname)) {
                            pdb.gaus = anam.read(fname)
                            if (mm>=10) {
                              k = mm-10
                            } else {
                              k=mm+2
                            }
                            pos = k * (one_var_dim*numvar) + p.stat
                            anam_mm = P[pos,1]
                            prec.anam = anam.eval(anam=pdb.gaus, val=anam_mm)
                            GA_P[pos,1] = prec.anam
                          }
                        }
                      }
                      if (i == 2) {
                        for (mm in 4:9) {
                          fname = paste0("/scratch3/veronika/reuse/assimil_data/gauss_ana/clim_70yr_pdiff/anam_",p.stat,"_",mm,".ascii")
                          if (file.exists(fname)) {
                            pdb.gaus = anam.read(fname)
                            k = mm-4
                            pos = k * (one_var_dim*numvar) + p.stat
                            anam_mm = P[pos,1]
                            prec.anam = anam.eval(anam=pdb.gaus, val=anam_mm)
                            GA_P[pos,1] = prec.anam
                          }
                        } 
                      }
                    }
                    # should multiply with H, but H is 1
                    # multiplying GA_P with the GA_obs, only the month in which the obs was taken, because we use time localization for instrumental data
                    if (i == 1) {
                      if (j<=length(calibrate$lat)/6) {
                        mm = 10
                        pobs = j
                      } else if (j>length(calibrate$lat)/6 & j<=length(calibrate$lat)/6*2) {
                        mm = 11
                        pobs = j - length(calibrate$lat)/6
                      } else if (j>length(calibrate$lat)/6*2 & j<=length(calibrate$lat)/6*3) {
                        mm = 12
                        pobs = j - length(calibrate$lat)/6*2
                      } else if (j>length(calibrate$lat)/6*3 & j<=length(calibrate$lat)/6*4) {
                        mm = 1
                        pobs = j - length(calibrate$lat)/6*3
                      } else if (j>length(calibrate$lat)/6*4 & j<=length(calibrate$lat)/6*5) {
                        mm = 2
                        pobs = j - length(calibrate$lat)/6*4
                      } else if (j>length(calibrate$lat)/6*5 & j<=length(calibrate$lat)/6*6) {
                        mm = 3
                        pobs = j - length(calibrate$lat)/6*5
                      }
                        fname = paste0("/scratch3/veronika/reuse/assimil_data/gauss_ana/GA_obs_pdiff/obs_",pobs,"_",mm,".ascii")
                        if (file.exists(fname)) {
                          pdb.gaus = anam.read(fname)
                          if (mm>=10) {
                            k = mm-10
                          } else {
                            k=mm+2
                          }
                          pos = k * (one_var_dim*numvar) + seq(4609,9216,1)
                          echam_mm = GA_P[pos,i]
                          pdb_mm <- db.create(echam_mm,ndim=1,autoname=F)
                          prec.mm = anam.z2y(pdb_mm,'echam_mm',anam=pdb.gaus)
                          GA_P[pos,i] = prec.mm@items$Gaussian.echam_mm
                        }
                      # after the transformation -> the cov matrix is localized
                      PH = GA_P * wgt
                    }
                    if (i == 2) {
                      if (j<=length(calibrate$lat)/6) {
                        mm = 4
                        pobs = j
                      } else if (j>length(calibrate$lat)/6 & j<=length(calibrate$lat)/6*2) {
                        mm = 5
                        pobs = j - length(calibrate$lat)/6
                      } else if (j>length(calibrate$lat)/6*2 & j<=length(calibrate$lat)/6*3) {
                        mm = 6
                        pobs = j - length(calibrate$lat)/6*2
                      } else if (j>length(calibrate$lat)/6*3 & j<=length(calibrate$lat)/6*4) {
                        mm = 7
                        pobs = j - length(calibrate$lat)/6*3
                      } else if (j>length(calibrate$lat)/6*4 & j<=length(calibrate$lat)/6*5) {
                        mm = 8
                        pobs = j - length(calibrate$lat)/6*4
                      } else if (j>length(calibrate$lat)/6*5 & j<=length(calibrate$lat)/6*6) {
                        mm = 9
                        pobs = j - length(calibrate$lat)/6*5
                      }
                        fname = paste0("/scratch3/veronika/reuse/assimil_data/gauss_ana/GA_obs_pdiff/obs_",pobs,"_",mm,".ascii")
                        if (file.exists(fname)) {
                          pdb.gaus = anam.read(fname)
                          k = mm-4
                          pos = k * (one_var_dim*numvar) + seq(4609,9216,1)
                          echam_mm = GA_P[pos,1]
                          pdb_mm <- db.create(echam_mm,ndim=1,autoname=F)
                          prec.mm = anam.z2y(pdb_mm,'echam_mm',anam=pdb.gaus)
                          GA_P[pos,1] = prec.mm@items$Gaussian.echam_mm
                        }
                       
                      # after the transformation -> the cov matrix is localized
                      PH = GA_P * wgt
                    }
                  } else { #original
                    PH <- (analysis$data[,i,] %*% t(x2) / (nens - 1) * wgt) %*% t(H)
                  }
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
                  if (calibrate$sour[j] == "inst" | calibrate$sour[j] == "doc") {
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
                if (gauss_ana) { # calculating HPHR, analysis$mean, analysis$data with "GA_H"
                  if (i == 1) {
                    pos = floor(h.i/(one_var_dim*numvar))
                    if (pos==0 | pos==1 | pos==2 ){
                      mm = pos + 10
                    } 
                    if (pos==3 | pos==4 | pos==5 ){
                      mm = pos - 2 
                    }
                  } 
                  if (i == 2) {
                    pos = floor(h.i/(one_var_dim*numvar))
                    mm = pos + 4
                  }
                  p.stat = h.i - pos * (one_var_dim*numvar)
                  fname = paste0("/scratch3/veronika/reuse/assimil_data/gauss_ana/clim_70yr_pdiff/anam_",p.stat,"_",mm,".ascii")
                  pdb.gaus = anam.read(fname)
                  #### transforming PH with GA_model
                  anam_mm = PH[h.i,]
                  prec.anam = anam.eval(anam=pdb.gaus, val=anam_mm)
                  GA_PH = prec.anam
                  # transforming analysis$ensmean with GA_model
                  anam_mm = analysis$ensmean[h.i,i]
                  prec.anam = anam.eval(anam=pdb.gaus, val=anam_mm)
                  GA_ana_ensmean = prec.anam
                  # transforming analysis$data with GA_model
                  anam_mm = analysis$data[h.i,i,]
                  prec.anam = anam.eval(anam=pdb.gaus, val=anam_mm)
                  GA_ana_data = prec.anam
                  #### transforming GA_PH with GA_obs
                  if (i == 1) {
                    if (j<=length(calibrate$lat)/6) {
                      mm = 10
                      pobs = j
                    } else if (j>length(calibrate$lat)/6 & j<=length(calibrate$lat)/6*2) {
                      mm = 11
                      pobs = j - length(calibrate$lat)/6
                    } else if (j>length(calibrate$lat)/6*2 & j<=length(calibrate$lat)/6*3) {
                      mm = 12
                      pobs = j - length(calibrate$lat)/6*2
                    } else if (j>length(calibrate$lat)/6*3 & j<=length(calibrate$lat)/6*4) {
                      mm = 1
                      pobs = j - length(calibrate$lat)/6*3
                    } else if (j>length(calibrate$lat)/6*4 & j<=length(calibrate$lat)/6*5) {
                      mm = 2
                      pobs = j - length(calibrate$lat)/6*4
                    } else if (j>length(calibrate$lat)/6*5 & j<=length(calibrate$lat)/6*6) {
                      mm = 3
                      pobs = j - length(calibrate$lat)/6*5
                    }
                  }
                  if (i == 2) {
                    if (j<=length(calibrate$lat)/6) {
                      mm = 4
                      pobs = j
                    } else if (j>length(calibrate$lat)/6 & j<=length(calibrate$lat)/6*2) {
                      mm = 5
                      pobs = j - length(calibrate$lat)/6
                    } else if (j>length(calibrate$lat)/6*2 & j<=length(calibrate$lat)/6*3) {
                      mm = 6
                      pobs = j - length(calibrate$lat)/6*2
                    } else if (j>length(calibrate$lat)/6*3 & j<=length(calibrate$lat)/6*4) {
                      mm = 7
                      pobs = j - length(calibrate$lat)/6*3
                    } else if (j>length(calibrate$lat)/6*4 & j<=length(calibrate$lat)/6*5) {
                      mm = 8
                      pobs = j - length(calibrate$lat)/6*4
                    } else if (j>length(calibrate$lat)/6*5 & j<=length(calibrate$lat)/6*6) {
                      mm = 9
                      pobs = j - length(calibrate$lat)/6*5
                    }
                  }
                  fname = paste0("/scratch3/veronika/reuse/assimil_data/gauss_ana/GA_obs_pdiff/obs_",pobs,"_",mm,".ascii")
                  pdb.gaus = anam.read(fname)
                  echam_mm = GA_PH
                  pdb_mm <- db.create(echam_mm,ndim=1,autoname=F)
                  prec.mm = anam.z2y(pdb_mm,'echam_mm',anam=pdb.gaus)
                  HPH = prec.mm@items$Gaussian.echam_mm
                  HPHR = as.vector(HPH + Rcal[j,i])
                  # transforming GA_ana_ensmean with GA_obs
                  echam_mm =  GA_ana_ensmean
                  pdb_mm <- db.create(echam_mm,ndim=1,autoname=F)
                  prec.mm = anam.z2y(pdb_mm,'echam_mm',anam=pdb.gaus)
                  GA_ana_ensmean = prec.mm@items$Gaussian.echam_mm
                  # transforming GA_ana_data with GA_obs
                  echam_mm =  GA_ana_data
                  pdb_mm <- db.create(echam_mm,ndim=1,autoname=F)
                  prec.mm = anam.z2y(pdb_mm,'echam_mm',anam=pdb.gaus)
                  GA_ana_data= prec.mm@items$Gaussian.echam_mm
                } else { # original
                  HPHR <- as.vector(H %*% PH[h.i,] + Rcal[j,i]) # original was only R[j] replaced by Rcal[j,i]   
                }
                if ( abs(HPHR) < 0.001 ) { # it happened when assimilating wetdays that both PH and Rcal were equal to 0 (eg. year=1986,i=2,j=1205)
                  K = matrix(0,dim(PH))
                  Ktilde = matrix(0,dim(PH))
                  write(paste("obs was excluded beacuse HPHR is smaller < 0.001 R value:",Rcal[j,i],"; lon:",calibrate$lon[j], 
                              "; lat:",calibrate$lat[j]), file=paste0('../log/',logfn),append=T)
                } else {
                  K <- PH / HPHR
                  Ktilde <- K / (1 + sqrt(Rcal[j,i]/HPHR)) # original was only R[j] replaced by Rcal[j,i]
                }
                if (gauss_ana) {
                  analysis$ensmean[,i] <- analysis$ensmean[,i] + K[,1] * (calibrate$data[j,i] - GA_ana_ensmean)
                  analysis$data[,i,] <- analysis$data[,i,] - Ktilde %*% GA_ana_data
                } else { # original
                  analysis$ensmean[,i] <- analysis$ensmean[,i] + K[,1] * (calibrate$data[j,i] -
                                                                            H %*% analysis$ensmean[h.i,i])
                  analysis$data[,i,] <- analysis$data[,i,] - Ktilde %*% H %*% analysis$data[h.i,i,]
                }
                if (update_PHclim){
                  ananomallts$ensmean[,i] <- ananomallts$ensmean[,i] + K[,1] * (calibrate$data[j,i] -
                                                                                  H %*% ananomallts$ensmean[h.i,i])
                  ananomallts$data[,i,] <- ananomallts$data[,i,] - Ktilde %*% H %*% ananomallts$data[h.i,i,]
                }
              } # end of !is.na(calibrate$data[j,i])
            } # end of ntim
          } # end if (!is.na(wgt[1]))
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
        if (gauss_ana) {
          # Transform an array of values from Gaussian into Raw values
          # the "full" ensembles have to be transformed back (dev + mean)
          analysis$data <- analysis$data + as.vector(analysis$ensmean)
          for (i in 1:2) {
            if (i == 1) {
              for (mm in c(10,11,12,1,2,3)) {
                for (p.stat in 4609:9216){
                  pdb.gaus = anam.read(filename=paste0("/scratch3/veronika/reuse/assimil_data/gauss_ana/clim_70yr_pdiff/anam_",p.stat,"_",mm,".ascii"))
                  if (mm>=10) {
                    k = mm-10
                  } else {
                    k=mm+2
                  }
                  pos = k * (one_var_dim*numvar) + p.stat
                  # analysis data
                  anam_mm = analysis$data[pos,i,]
                  prec.anam = anam.eval(anam=pdb.gaus ,val=anam_mm)
                  analysis$data[pos,i,] = prec.anam
                }
              }
              analysis$ensmean[,i] = rowMeans(analysis$data[,i,])
            }
            if (i == 2) {
              for (mm in 4:9) {
                for (p.stat in 4609:9216){
                  pdb.gaus = anam.read(filename=paste0("/scratch3/veronika/reuse/assimil_data/gauss_ana/clim_70yr_pdiff/anam_",p.stat,"_",mm,".ascii"))
                  k = mm-4
                  pos = k * (one_var_dim*numvar) + p.stat
                  # analysis data
                  anam_mm = analysis$data[pos,i,]
                  prec.anam = anam.eval(anam=pdb.gaus ,val=anam_mm)
                  analysis$data[pos,i,] = prec.anam
                }
              }  
              analysis$ensmean[,i] = rowMeans(analysis$data[,i,])
            }
          }
        } else { # original
          analysis$data <- analysis$data + as.vector(analysis$ensmean)
        }
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
          if (precip_ratio) {
            analysis.abs$data[analysis.abs$names=='precip',,] <- 
              analysis.anom$data[analysis.anom$names=='precip',,] * echam.clim$data[echam.clim$names=='precip',,]
            analysis.abs$ensmean[analysis.abs$names=='precip',] <- 
              analysis.anom$ensmean[analysis.anom$names=='precip',] * echam.clim$ensmean[echam.clim$names=='precip',]
            echam.abs$data[echam.abs$names=='precip',,] <- 
              echam.anom$data[echam.anom$names=='precip',,] * echam.clim$data[echam.clim$names=='precip',,]
            echam.abs$ensmean[echam.abs$names=='precip',] <- 
              echam.anom$ensmean[echam.anom$names=='precip',] * echam.clim$ensmean[echam.clim$names=='precip',]
          }
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
                   file=paste0(dataintdir,'analysis/EKF400_',version,'_',expname,'/analysis_',cyr,'_2ndgrid.Rdata'),version=2)
            } else {
              save(analysis.anom,analysis.abs,echam.anom,echam.abs,validate,calibrate, 
                   file=paste0(dataintdir,'analysis/EKF400_',version,'_',expname,'/analysis_',cyr,'_2ndgrid.Rdata'),version=2)
            }
          } else {
            save(analysis.anom,analysis.abs,echam.anom,echam.abs,validate,calibrate,
                 file=paste0(dataintdir,'analysis/EKF400_',version,'_',expname,'/analysis_',cyr,'.Rdata'),version=2)
          }
        } else {
          if (every2grid){    
            save(analysis.anom,analysis.abs,echam.anom,echam.abs,calibrate,
                 file=paste0(dataintdir,'analysis/EKF400_',version,'_',expname,'/analysis_',cyr,'_2ndgrid.Rdata'),version=2)
          } else {
            save(analysis.anom,analysis.abs,echam.anom,echam.abs,calibrate,
                 file=paste0(dataintdir,'analysis/EKF400_',version,'_',expname,'/analysis_',cyr,'.Rdata'),version=2)
          }
        }
      } else {
        if (vali) {
          save(analysis.anom,analysis.abs,echam.anom,echam.abs,landcorrected.anom,
               landcorrected.clim,validate,calibrate,
               file=paste0(dataintdir,'analysis/EKF400_',version,'_',expname,'/analysis_',cyr,'_2ndgrid.Rdata'),version=2)
        } else {
          save(analysis.anom,analysis.abs,echam.anom,echam.abs,landcorrected.anom,
               landcorrected.clim,calibrate,
               file=paste0(dataintdir,'analysis/EKF400_',version,'_',expname,'/analysis_',cyr,'_2ndgrid.Rdata'),version=2)
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
            #nmonths <- 6
            ndim <- nrow(looanalysis$data)
            ntim <- ncol(looanalysis$data)
            nens <- dim(looanalysis$data)[3]
            nprox <- nrow(loocalibrate$data) 
            #ndimold <- length(echam$lon) / nmonths
            #itime <- rep(1:nmonths, each=ndimold)
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
                    looanalysis$ensmean[,i] <- looanalysis$ensmean[,i] + K[,1] * as.vector(loocalibrate$data[j,i] -
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
             file=paste0(dataintdir,'loo/EKF400_',version,'_',expname,'/loo_results_',cyr,'.Rdata'))
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
      stop('old fortran function not tested for a long time')
      # R = array(1,c(nrow(Hcal), s))
      # R[calibrate$names=='precip',]<-try(abs(calibrate$data[calibrate$names=='precip',]
      #                                        *0.5))  
      # R[calibrate$names=='slp',]<-try(abs(calibrate$data[calibrate$names=='slp',]
      #                                     *0.1))  
      # R[which(is.na(R))] <- 0
      # if (fortran_ensrf) {calibrate$data[is.na(calibrate$data)]=9e36}  
      # analysis <- EnSRF_new(echam, calibrate,  R=R, Hcal=Hcal, weights=d.weights_all)
      # #      analysis <- EnSRF4(echam, calibrate,  R=R, Hcal=Hcal, weights=d.weights_all)
      if (vali) {
        if (every2grid){
          save(analysis,echam,validate,calibrate, #calibrate.allts,
               file=paste0('../data/analysis/EKF400_',version,'_',expname,'/analysis_',cyr,'_2ndgrid.Rdata'))
        } else {
          save(analysis,echam,validate,calibrate, #calibrate.allts,
               file=paste0('../data/analysis/EKF400_',version,'_',expname,'/analysis_',cyr,'.Rdata'))
        }
      } else {
        if (every2grid){
          save(analysis,echam,calibrate, #calibrate.allts,
               file=paste0('../data/analysis/EKF400_',version,'_',expname,'/analysis_',cyr,'_2ndgrid.Rdata'))
        } else {  
          save(analysis,echam,calibrate, #calibrate.allts,
               file=paste0('../data/analysis/EKF400_',version,'_',expname,'/analysis_',cyr,'.Rdata'))
        }
      }
    }
  }
  print("calc time for a year")
  print(proc.time() - ptm1)
} # end of yearly analysis loop

warnings()
#quit(save='no')


