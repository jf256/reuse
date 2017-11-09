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
syr=1942
eyr=1943
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
  workdir('/scratch/veronika/rerun/r_code')
} else if (user=="joerg") {
  workdir='/scratch3/joerg/projects/reuse/reuse_git/'
} else {
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

##########################################################################################
# start loading data  
##########################################################################################  
if (generate_ECHAM){
  print("generate_ECHAM")
  # create ECHAM .RData from .nc
  # ACHTUNG: members 001 to 030 with land surface bug, 103 and 2xx without
  read_echam4('EnSRF.ccc400_0', timlim=c(1601,2005), small=every2grid, landonly=land_only)
} 

if (generate_ECHAM_anom){
  # read echam 71yr anom, clim and sd calculated with cdo from .nc files to .RData
  print("generate_ECHAM_anom")
  read_echam4('ano', path=echanompath, timlim=c(1601,2005), small=every2grid, 
              landonly=land_only, anom=T)
  read_echam4('EnSRF', path=echclimpath, timlim=c(1635,1970), small=every2grid, 
              landonly=land_only, clim=T)
  read_echam4('EnSRF', path=echsdpath, timlim=c(1601,2005), small=every2grid, 
              landonly=land_only, std=T)
}

if (generate_ECHAM_1901_70){
  print("generate_ECHAM_1901_70")
  # ECHAM data for bias calculation with real proxy data
  echam1901_70 <- read_echam_ensmean('EnSRF', timlim=c(1901,1970),small=F)
  save(echam1901_70, file="../data/echam_1911-70.Rdata")
} 

if (generate_ECHAM_103){
  print("generate_ECHAM ens. mem. 103")
    # ens. member WITHOUT land surface bug
    read_echam4(filehead='EnSRF', path=paste0(ext,'echam/echam103'), timlim=c(syr,eyr), 
                small=every2grid, landonly=land_only)
} 

if (generate_ECHAM_covar){
  print("generate_ECHAM all time step array long-term covariance")
  echam_covar(syr=1603,eyr=2004)
}

# if (generate_ind_recon){
#   # read Broennimann et al. 2009 atm. indices from .txt to .Rdata for comparison
#   #if (syr<1901){syr_ind=1901} else {syr_ind=syr}
#   #if (eyr>2004){eyr_ind=2004} else {eyr_ind=eyr}
#   #syr_ind=1901
#   #eyr_ind=2004
#   ind=read.table(file=paste(dataintdir,'/indices/stefan/stefan_monthly_indices.txt'
#                             ,sep=''),header=T)   
#   ind_rec_dimi = window(ts(c(rep(NA,length(ind[,colnames(ind) == 'Z100'])),rep(NA,132)),
#                            start=ind[1,colnames(ind) == 'yr'],freq=12),syr_ind,c(eyr_ind,12))
#   ind_rec_z100 = window(ts(ind[,colnames(ind) == 'Z100'],start=ind[1, 
#                            colnames(ind) == 'yr'],freq=12),syr_ind,freq=12,c(eyr_ind,12))
#   ind_rec_z300 = window(ts(ind[,colnames(ind) == 'Z300'],start=ind[1,
#                            colnames(ind) == 'yr'],freq=12),syr_ind,c(eyr_ind,12))
#   ind_rec_pwc = window(ts(ind[,colnames(ind) == 'PWC'],start=ind[1,
#                            colnames(ind) == 'yr'],freq=12),syr_ind,c(eyr_ind,12))
#   ind_rec_hc = window(ts(ind[,colnames(ind) == 'HCL'],start=ind[1,
#                            colnames(ind) == 'yr'],freq=12),syr_ind,c(eyr_ind,12))
#   ind_rec_sj = window(ts(ind[,colnames(ind) == 'SJ'],start=ind[1,
#                            colnames(ind) == 'yr'],freq=12),syr_ind,c(eyr_ind,12))
#   indall = t(cbind(ind_rec_dimi, ind_rec_z100, ind_rec_z300, ind_rec_pwc, ind_rec_hc, 
#                   ind_rec_sj))
#   save(indall, file=paste("../data/indices/indices_recon_",syr_ind,"-",eyr_ind,".Rdata",sep=""))
# }

if (generate_NCEP) {
  print("generate_NCEP")
  #  see script in EnSRF/script/merge_ncep.sh for regridding and co of orig. ncep data set
  ncepall <- read_echam1('ncep_allvar_1948-2009',timlim=c(syr_ncep,eyr_ncep),
                          path=nceppath,small=every2grid)
  if (every2grid) {
    save(ncepall, file=paste0("../data/ncep/ncep_allvar_",syr_ncep,"-",eyr_ncep,"_2ndgrid.Rdata"))
  } else {
    save(ncepall, file=paste0("../data/ncep/ncep_allvar_",syr_ncep,"-",eyr_ncep,".Rdata"))
  }
}

if (generate_CRUALLVAR) {
  print("generate_CRUALLVAR")
  #  see script in EnSRF/script/merge_cru.sh for regridding and co of orig. cru data set
  cruall <- read_echam1('cru_allvar_abs_1901-2004.nc',timlim=c(syr_cru,eyr_cru),
                         path=crupath,small=every2grid,landonly=land_only)
  if (every2grid) {
    save(cruall, file=paste0(dataintdir,"cru/cru_allvar_",syr_cru,"-",eyr_cru,"_2ndgrid.Rdata"))
  } else {
    save(cruall, file=paste0(dataintdir,"cru/cru_allvar_",syr_cru,"-",eyr_cru,".Rdata"))  
  }
}

if (generate_HadCRU4){ 
  print("generate_HadCRU4 standard deviations")
  # from instr. CRU ensemble, precalculated with cdo
  cru4_may_sep <- read_echam1('HadCRUT4_ens_sd_may-sep_yrmean.nc',
                              timlim=c(1901,2002),path=crupath,small=every2grid,landonly=F)
  cru4_oct_apr <- read_echam1('HadCRUT4_ens_sd_oct-apr_yrmean.nc',
                              timlim=c(1901,2002),path=crupath,small=every2grid,landonly=F)
  cru4_may_sep$ensmean <- cru4_may_sep$data
  cru4_oct_apr$ensmean <- cru4_oct_apr$data
  if (every2grid) {
    save(cru4_may_sep,cru4_oct_apr, file="../data/cru/cru4_ens_sd_2ndgrid.Rdata")
  } else {
    save(cru4_may_sep,cru4_oct_apr, file="../data/cru/cru4_ens_sd.Rdata")  
  }
}

if (generate_LUTPAULKUT){ 
  print("generate_LUTPAULKUT")
  # ATTENTION: seasonal resolution
  #  see script in EnSRF/script/merge_recon.sh for regridding and co of orig. recon data set
  reconall <- read_echam1('recon_allvar_1750-1999',xlim=c(-180,180), ylim=c(-90,90), 
                          timlim=c(syr_recon,eyr_recon),path=reconpath,small=every2grid,
                          landonly=land_only)
  reconall$data[reconall$names=="precip",,]<-reconall$data[reconall$names=="precip",,]/3
  if (every2grid) {
    save(reconall, file=paste0(dataintdir,"recon/recon_allvar_",syr_recon,"-",eyr_recon,"_2ndgrid.Rdata"))
  } else {
    save(reconall, file=paste0(dataintdir,"recon/recon_allvar_",syr_recon,"-",eyr_recon,".Rdata"))  
  }
}

if (generate_GHCN){
  print("generate_GHCN")
  # created only once for the full period and then cut after loading
  ghcn <- read_ghcn_refyr(1600,2005,1600,1869)
  ghcn$names <-rep('temp2',length(ghcn$names))
  save(ghcn, file=paste0("../assim_data/ghcn/ghcn_temp",fsyr,"-",feyr,".Rdata"))
}

if (generate_GHCN_precip){
  print("generate_GHCN_precip")
  # created only once for the full period and then cut after loading
  ghcn_precip <- read_ghcn_refyr_precip(1600,2005,1600,1869)
  ghcn_precip$data <- ghcn_precip$data / 10 # to make units echam conform
  ghcn_precip$names <-rep('precip',length(ghcn_precip$names))
  save(ghcn_precip, file=paste0("../assim_data/ghcn/ghcn_precip_",fsyr,"-",feyr,".Rdata")) 
}

if (generate_t_yuri){
  print("generate_t_yuri")
  source("../assim_data/data_yuri/t_assimil/read_all.R")
}

if (generate_slp_yuri){
  print("generate_slp_yuri")
  source("../assim_data/data_yuri/slp_assimil/read_all.R")
}

if (generate_DOCUM){
  print("generate_DOCUM")
  source("../assim_data/data_yuri/t_docu/read_seas.R")
  source("../assim_data/data_yuri/t_docu/read_monthly.R")
  source("../assim_data/data_yuri/t_docu/read_JFMA.R")
  source("../assim_data/data_yuri/t_docu/read_AMJJA.R")
}

if (generate_PROXIES){
  print("generate_PROXIES")
  # real trw proxy multiple regression approach
  # only with monthly state vector of 6 months
  if (trw_only) {
     realprox <- read_proxy2(fsyr,feyr)
  } else if (schweingr_only) {
    realprox <- read_proxy_schweingr(fsyr,feyr)
  } else if (mxd_only) {
    realprox <- read_proxy_mxd(fsyr,feyr)
  } else {
    schprox <- read_proxy_schweingr(fsyr,feyr)
    mxdprox <- read_proxy_mxd(fsyr,feyr)
    trwprox <- read_proxy2(fsyr,feyr)
# ORIG VERSION JOERG
    ## add NA for TRW data that ends 1970
    #trwprox$data <- rbind(trwprox$data,matrix(NA,nrow=length(seq(1971,2004)),
    #                  ncol=dim(trwprox$data)[2]))
# NEW VERSION VERONIKA: PLEASE CHECK IF YOUR VERSION IS CORRECT
    ############# # because the dim is not equal, the year for trwprox stops in 1970, other two in 2005
    trwprox$data <- rbind(trwprox$data,matrix(data=NA, nrow=dim(mxdprox$data)[1]-
                      dim(trwprox$data)[1], ncol=dim(trwprox$data)[2]))
    #############
    trwprox$time <- c(trwprox$time,seq(1971,2004))
    realprox <- list()
    realprox$data <- cbind(mxdprox$data, schprox$data, trwprox$data)
    realprox$lon <- c(mxdprox$lon, schprox$lon, trwprox$lon)
    realprox$lat <- c(mxdprox$lat, schprox$lat, trwprox$lat)
    realprox$time <- mxdprox$time
    realprox$mr <- rbind(mxdprox$mr, schprox$mr, trwprox$mr)
    realprox$var_residu <- c(mxdprox$var_residu, schprox$var_residu, trwprox$var_residu)
  }
  save(realprox, file=paste0("../data/proxies/real_proxies_",fsyr,"-",feyr,".Rdata"))
} 






                                     

##########################################################################################
# 0. Loop over years to reduce size of state vector
##########################################################################################
if (sixmonstatevector) {syr2=syr+1} else {syr2=syr}
for (cyr in syr2:eyr) {
  print(cyr)
  write(cyr,file=paste0('../log/',logfn),append=T)
  if (cyr > 1659) {
    instrumental=T        # all instrumental stations
  } else {
    instrumental=F
  }
  if (cyr < 1960) {
    real_proxies=T         # Proxy data experiment (regression NOT H operator)
  } else {
    real_proxies=F
  }
  if (cyr > 1853) {
    docum=F                 # read documentary based data
  } else {
    docum=T
  }
  # next line not included yet: 
  if (eyr < 1750) {
    vali=F                 # switch off prepplot if no vali data selected
  } else {
    vali=T
  }
  if ((cyr > 1901) & (eyr < 2006)) {
    cru_vali=T             # monthly CRU TS3 temp, precip and HADSLP2 gridded instrumentals (1901-2004)
    #  ind_recon=T            # Stefan's reconstructed indices until 1948 and NCAR reanalysis later added to CRU and NCEP
  } else {
    cru_vali=F 
    #  ind_recon=F
  }
  #if ((syr < 1901) & (eyr > 1749)) {
  #  recon_vali=T           # seasonal luterbacher, pauling, kuettel recons (1750-1999)
  #} else {
    recon_vali=F
  #}
  print(paste("instr:",instrumental, "; proxies:",real_proxies, "; documentary:",docum, 
              "; validation data:",vali,"; CRU:",cru_vali,"; Recon:",recon_vali))
  write(paste("instr:",instrumental, "; proxies:",real_proxies, "; documentary:",docum, 
              "; validation data:",vali,"; CRU:",cru_vali,"; Recon:",recon_vali),
        file=paste0('../log/',logfn),append=T)
  asyr <- cyr-35 # "a" for anomaly
  if (asyr < 1601) {asyr = 1601}
  aeyr <- cyr+35
  if (aeyr > 2005) {aeyr = 2005}
  ptm1 <- proc.time()
  
  ##########################################################################################
  # 1. Echam Part
  # 1.1 Loading echam, echam_anom, echam_clim, landcorrected_anom, landcorrected_clim
  # 1.2 Choose which variables want to use from the model
  # 1.3 Calc echam st. dev. for each grid point and month over ens memb. to scale docu data
  # 1.4 Calculate decorrelation distance
  # 1.5 Just leave data for one year in memory and convert to sixmonstatevector data format
  # 1.6 Set up the cutoff distance for localisation
  ##########################################################################################
  
  # 1.1 Loading echam, echam_anom, echam_clim, landcorrected_anom, landcorrected_clim
  if ((cyr==syr2) & (covarclim>0)) {
    load(file="../data/echam/echallts_for_covar.Rdata")
    #just use limited number of years (n_covar) to make calculation faster
    echanomallts$data <- echanomallts$data[,,sample(seq(1,dim(echanomallts$data)[3]),n_covar)]
  }
  if (every2grid) {
    load(paste(dataextdir,"echam/echam_",(cyr-1),"-",cyr,"_2ndgrid.Rdata",sep=""))
    #load(paste(dataintdir,"echam/echam_",(cyr-1),"-",cyr,"_2ndgrid.Rdata",sep=""))
  } else {
    load(paste(dataextdir,"echam/echam_",(cyr-1),"-",cyr,".Rdata",sep=""))
    #load(paste(dataintdir,"echam/echam_",(cyr-1),"-",cyr,".Rdata",sep=""))
  }
  echam.sd <- echam # to have non-sixmon echam in docum section
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
  
  if (no_forc_big_ens) {
    yrs <- floor(runif(n_no_forc,1602,2004))
    m <- floor(runif(n_no_forc,1,30))
    for (n in 1:n_no_forc) {
      yr1 <- yrs[n]
      yr2 <- yr1+1
      if (every2grid) {
        load(paste0(echanompath,'echam_anom_',yr1,'-',yr2,'_2ndgrid.Rdata'))
      } else {
        load(paste0(echanompath,'echam_anom_',yr1,'-',yr2,'.Rdata'))
      }
      if (n==1){
        echam_anom_data <- echam_anom$data[,,m[n]]
      } else {
        echam_anom_data <- abind(echam_anom_data,echam_anom$data[,,m[n]],along=3)
      }
    }
    echam_anom$data <- echam_anom_data
    echam_anom$ensmean <- apply(echam_anom_data,1:2,mean)
  } # end no_forc_big_ens
    
  # 1.2 Choose which variables want to use from the model
  # just leave temp precip slp in state vector
  if (tps_only) {
    tpspos <- c(which(echam$names=='temp2'), which(echam$names=='precip'),
                which(echam$names=='slp'), which(echam$names=='bias'))
    echam$data <- echam$data[tpspos,,]
    echam$ensmean <- echam$ensmean[tpspos,]
    echam$names <- echam$names[tpspos]
    if (anomaly_assim){
      tpspos <- c(which(echam_anom$names=='temp2'), which(echam_anom$names=='precip'),
                  which(echam_anom$names=='slp'), which(echam_anom$names=='bias'))
      echam_anom$data <- echam_anom$data[tpspos,,]
      echam_anom$ensmean <- echam_anom$ensmean[tpspos,]
      echam_anom$names <- echam_anom$names[tpspos]
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
    # new echam_clim made by veronika have correct units, 
    # hence comment following corrections
    # error where echam_anom and echam_clim were generated initially: 
    # no unit correction happened, thus here
    echam_anom$data[echam_anom$names=='precip',,] <- 
      echam_anom$data[echam_anom$names=='precip',,] * 3600 * 24 * 30
    echam_anom$ensmean[echam_anom$names=='precip',] <- 
      echam_anom$ensmean[echam_anom$names=='precip',] * 3600 * 24 * 30
    echam_anom$data[echam_anom$names=='slp',,] <- echam_anom$data[echam_anom$names=='slp',,]/100
    echam_anom$ensmean[echam_anom$names=='slp',] <- echam_anom$ensmean[echam_anom$names=='slp',]/100
    echam.anom <- echam_anom
    # echam_clim$data[echam_clim$names=='temp2',,] <- echam_clim$data[echam_clim$names=='temp2',,]-273.15
    # echam_clim$ensmean[echam_clim$names=='temp2',] <- echam_clim$ensmean[echam_clim$names=='temp2',]-273.15
    # echam_clim$data[echam_clim$names=='precip',,] <- echam_clim$data[echam_clim$names=='precip',,]*3600 * 24 * 30
    # echam_clim$ensmean[echam_clim$names=='precip',] <- echam_clim$ensmean[echam_clim$names=='precip',]*3600 * 24 * 30
    # echam_clim$data[echam_clim$names=='slp',,] <- echam_clim$data[echam_clim$names=='slp',,]/100
    # echam_clim$ensmean[echam_clim$names=='slp',] <- echam_clim$ensmean[echam_clim$names=='slp',]/100
    echam.clim <- echam_clim
    # echam_clim_mon_ensmean (next line) stays 12 months version for instr screening
    echam_clim_mon_ensmean <- echam_clim$ensmean[,10:21]
    rm(echam_anom,echam_clim)
    if (landcorr) {
      # Veronika thinks they (temp2, precip, slp) are corrected for the climatology
      # Veronika thinks they (precip, slp) are NOT corrected for the anomaly
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

# VERONIKA: PLEASE CHECK IF THE NEXT LINES ARE NEEDED
  #tmp1 <- array(echam$data,c(dim(echam$data)[1]*dim(echam$data)[2],dim(echam$data)[3]))
  # cut oct cyr-1 to sep cyr
  #tmp2 <- tmp1[((9*dim(echam$data)[1]+1):(dim(tmp1)[1]-(3*dim(echam$data)[1]))),] 
  #if (docum) { # probabaly later this if will be deleted, now I just would like to test for the docu data
  #    # if I want to use the echam.sd in the 6monstatevector format then maybe I could merge this point to 1.5
  #    # or maybe simply enough to do the first if in 1.5
  #    echam$data <- array(tmp2,c(dim(tmp2)[1]/(((dim(echam$data)[2]/12)-1)*2), (((dim(echam$data)[2]/12)-1)*2),dim(tmp2)[2]))
  #} else {
  #    echam$data <- array(tmp2,c(dim(tmp2)[1]/12,12,dim(tmp2)[2]))
  #}
  #rm(tmp1);rm(tmp2)
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
    d <- compute_dist_2d(echam$lon,echam$lat,echam$lon,echam$lat)
    for (i in unique(echam$names)) {
      tmp <- echam.anom$ensmean[echam$names==i,]
      corens <- cor(t(tmp[,]))
      png(paste0('../figures/decorr_',i,'.png'), width = 1024, height = 768)
      plot(as.vector(d),as.vector(corens),col='#4c8bff01',
           xlim=c(0,5000),ylim=c(0,1))
      lines(exp(-1/2 * (seq(1:5000)/get(paste0('l_dist_',i)))**2),col='red')
      dev.off()
    }
  }
  
  
  # 1.5 Just leave data for one year in memory and convert to sixmonstatevector data format
  if (sixmonstatevector) {
    # change array to have 6 months in state vector for winter and summer
    # first winter starts in oct of syr
    # 6 mon stat vectors for oct-mar and apr and sep
    tmp21 <- array(echam$data,c(dim(echam$data)[1]*dim(echam$data)[2],
               dim(echam$data)[3]))
    tmp22 <- tmp21[((9*dim(echam$data)[1]+1):(dim(tmp21)[1]-(3*dim(echam$data)[1]))),]
    echam$data <- array(tmp22,c(dim(tmp22)[1]/(((dim(echam$data)[2]/12)-1)*2),
                    (((dim(echam$data)[2]/12)-1)*2),dim(tmp22)[2]))
    tmp31 <- array(echam$ensmean,c(dim(echam$ensmean)[1]*dim(echam$ensmean)[2]))
    tmp32 <- tmp31[((9*dim(echam$ensmean)[1]+1):(dim(tmp31)[1]-
               (3*dim(echam$ensmean)[1])))]
    echam$ensmean <- array(tmp32,c(dim(tmp32)[1]/(((dim(echam$ensmean)[2]/12)-1)*2),
                       (((dim(echam$ensmean)[2]/12)-1)*2)))
    echam$time <- c(cyr,cyr+0.5)
    tmp41 <- array(echam.anom$data,c(dim(echam.anom$data)[1]*dim(echam.anom$data)[2],
               dim(echam.anom$data)[3]))
    tmp42 <- tmp41[((9*dim(echam.anom$data)[1]+1):(dim(tmp41)[1]-(3*dim(echam.anom$data)[1]))),]
    echam.anom$data <- array(tmp42,c(dim(tmp42)[1]/(((dim(echam.anom$data)[2]/12)-1)*2),
                         (((dim(echam.anom$data)[2]/12)-1)*2),dim(tmp42)[2]))
    tmp51 <- array(echam.anom$ensmean,c(dim(echam.anom$ensmean)[1]*dim(echam.anom$ensmean)[2]))
    tmp52 <- tmp51[((9*dim(echam.anom$ensmean)[1]+1):(dim(tmp51)[1]-
               (3*dim(echam.anom$ensmean)[1])))]
    echam.anom$ensmean <- array(tmp52,c(dim(tmp52)[1]/(((dim(echam.anom$ensmean)[2]/12)-1)*2),
                            (((dim(echam.anom$ensmean)[2]/12)-1)*2)))
    echam.anom$time <- c(cyr,cyr+0.5)
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
    tmp61 <- array(echam.clim$data,c(dim(echam.clim$data)[1]*dim(echam.clim$data)[2],
               dim(echam.clim$data)[3]))
    tmp62 <- tmp61[((9*dim(echam.clim$data)[1]+1):(dim(tmp61)[1]-(3*dim(echam.clim$data)[1]))),]
    echam.clim$data <- array(tmp62,c(dim(tmp62)[1]/(((dim(echam.clim$data)[2]/12)-1)*2),
                         (((dim(echam.clim$data)[2]/12)-1)*2),dim(tmp62)[2]))
    tmp71 <- array(echam.clim$ensmean,c(dim(echam.clim$ensmean)[1]*dim(echam.clim$ensmean)[2]))
    tmp72 <- tmp71[((9*dim(echam.clim$ensmean)[1]+1):(dim(tmp71)[1]-
               (3*dim(echam.clim$ensmean)[1])))]
    echam.clim$ensmean <- array(tmp72,c(dim(tmp72)[1]/(((dim(echam.clim$ensmean)[2]/12)-1)*2),
                            (((dim(echam.clim$ensmean)[2]/12)-1)*2)))
    echam.clim$time <- c(cyr,cyr+0.5)
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
    rm(tmp41);rm(tmp42);rm(tmp51);rm(tmp52);rm(tmp61);rm(tmp62);rm(tmp71);rm(tmp72)
  }
  
  # rename echam.amon to echam if anomaly_assim==T
  if (anomaly_assim){
    echam=echam.anom
    rm(echam.anom)
  } 
  
  # convert lon, lat, names sixmonstatevector data format
  # ACHTUNG: only 11 vars withOUT stream function included so far
  if (sixmonstatevector) {
    numvar <- length(unique(echam$names)) 
    echam$lon <-  c(rep(echam$lon, (numvar*6)))
    echam$lat <- c(rep(echam$lat, (numvar*6)))
    echam$names <- c(rep(echam$names, 6))
    echam <- echam[c('data', 'ensmean', 'lon', 'lat', 'height', 'lsm.i', 'time', 'names')]
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
  lvec['t850'] <- l_dist_t850

  print('calc time for loading data')
  print(proc.time() - ptm1)
  
 
  
   
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
    if (every2grid) {
      if (ncep_vali) {load(paste("../data/ncep/ncep_allvar_",syr_ncep,"-",eyr_ncep,"_2ndgrid.Rdata",sep=""))
      } else if (recon_vali) {load(paste("../data/recon/recon_allvar_",syr_recon,"-",eyr_recon,"_2ndgrid.Rdata",sep=""))
      } else if (cru_vali) {load(paste("../data/cru/cru_allvar_",syr_cru,"-",eyr_cru,"_2ndgrid.Rdata",sep=""))} 
    } else {
      if (ncep_vali) {load(paste("../data/ncep/ncep_allvar_",syr_ncep,"-",eyr_ncep,".Rdata",sep=""))
      } else if (recon_vali) {load(paste("../data/recon/recon_allvar_",syr_recon,"-",eyr_recon,".Rdata",sep=""))
      } else if (cru_vali) {load(paste("../data/cru/cru_allvar_",syr_cru,"-",eyr_cru,".Rdata",sep=""))} 
    }
    # if (ind_recon) {
    #   load(file=paste("../data/indices/indices_recon_",syr,"-",eyr,".Rdata",sep=""))
    # }
    if (cru_vali) {
      valiall <- cruall
      valiall$data <- valiall$data[,,1]
    } else if (ncep_vali) {
      valiall <- ncepall
    } else if (recon_vali) {
      valiall <- reconall
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
    if (cru_vali) {
      ti=which(floor(valiall$time)==(cyr-1) | floor(valiall$time)==cyr) 
      sts=ti[1]
      ets=ti[length(ti)]
      valiall$data=valiall$data[,sts:ets]
      valiall$time=valiall$time[sts:ets]
    }
    if (ncep_vali) {
      ti=which(floor(valiall$time)==(cyr-1) | floor(valiall$time)==cyr) 
      sts=ti[1]
      ets=ti[length(ti)]
      valiall$data=valiall$data[,sts:ets]
      valiall$time=valiall$time[sts:ets]
    }
    if (recon_vali) {
      ti=which(floor(valiall$time)==(cyr-1) | floor(valiall$time)==cyr) 
      sts=ti[1]
      ets=ti[length(ti)]
      valiall$data=valiall$data[,sts:ets]
      valiall$time=valiall$time[sts:ets]
    }
  
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
    if (sum(c(ncep_vali,cru_vali,recon_vali))>1) {
      print("WARNING: more than 1 validation data set selected!")
      write("WARNING: more than 1 validation data set selected!",
            file=paste0('../log/',logfn),append=T)
    } else {
      write(paste("ncep_vali:",ncep_vali,"; cru_vali:",cru_vali,"; recon_vali:",recon_vali),
            file=paste0('../log/',logfn),append=T)
    }

    # 2.6 Set validate$ensmean equal to validate$data
    validate=valiall
    validate$ensmean=validate$data
  }
  
  
  
  
  ##########################################################################################
  # 3. Real Proxy Part
  # 3.1 Loading documentary data
  # 3.2 Screen the proxy data: more than 5 std. dev. from mean
  # 3.3 Convert to 2 season per year
  # 3.4 Calculate the anomalies
  # 3.5 Create a list named proxies ("combine assimilation data")
  # 3.6 Set real_proxies to FALSE if there is no data
  ##########################################################################################

  # 3.1 Loading documentary data
  if (real_proxies){
    load(paste0("../data/proxies/real_proxies_",fsyr,"-",feyr,".Rdata"))  
  

    # 3.2 Screen the proxy data
    if (check_assimdata) {
      # correlation screening already where multiple regression coefficients are calculated
      # thus, screen if value at current time step is more than 5 std. dev. from mean
      # in this case treated as outlier and set to NA
      for (i in 1:length(realprox$lon)) {
        tiv=which(floor(realprox$time)==cyr)
        ti=which(floor(realprox$time)>=(cyr-35) & floor(realprox$time)<=(cyr+35))
        sts=ti[1]
        ets=ti[length(ti)]
        rpmean <- mean(realprox$data[sts:ets,i],na.rm=T)
        rpsd   <- sd(realprox$data[sts:ets,i],na.rm=T)
        #rpmean <- mean(realprox$data[,i],na.rm=T)
        #rpsd   <- sd(realprox$data[,i],na.rm=T)
        if ((!is.na(rpmean)) & (!is.na(rpsd))) {
          if ((!is.na(realprox$data[tiv,i])) & ((realprox$data[tiv,i] < rpmean-5*rpsd)
                                                | (realprox$data[tiv,i] > rpmean+5*rpsd))) {
            realprox$data[,i] <- NA
            print(paste('proxy data', i, 'out of range'))
            write(paste('proxy data', i, 'out of range'),file=paste0('../log/',logfn),append=T)
          }
        }
      }
    }
  
    # 3.3 Convert to 2 season per year
    # no scaling because regresion takes care of it
    realprox.allts <- realprox
    tmp1=t(realprox$data)
    tmp2=array(NA,c(dim(tmp1)[1],2,dim(tmp1)[2]))
    tmp2[,2,]=tmp1
    realprox.allts$data=array(tmp2,c(dim(tmp2)[1],dim(tmp2)[2]*dim(tmp2)[3]))
    realprox.allts$data=realprox.allts$data[,3:dim(realprox.allts$data)[2]] # why is first year removed?
    realprox.allts$time=seq(fsyr+1,feyr+0.5,0.5) 
    ti=which(floor(realprox.allts$time)==cyr)
    sts=ti[1]
    ets=ti[length(ti)]      
    realprox$data=realprox.allts$data[,sts:ets]
    realprox$time=realprox.allts$time[sts:ets]
    realprox$names=rep("prox",dim(realprox.allts$data)[1])
    
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
                    lat=realprox$lat, names=realprox$names,
                    height=realprox$elevation, time=realprox$time,
                    mr=realprox$mr, var_residu=realprox$var_residu,
                    numavg=rep(1,length(realprox$lon)))
    }
  
  
    # 3.6 Set real_proxies to FALSE if there is no data
    if (dim(realprox$data)[1]==0) { real_proxies=F }
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
    load('../assimil_data/data_yuri/t_docu_monthly.Rdata')
    if (!any(!is.na(t$data))) { docu=F }
    
    if (sixmonstatevector) {
      # 4.2 Calculate the 71 anomaly from cyr-1 October to cyr September
      if (anomaly_assim){
        # 4.2.1 Calculate the climatology
        t.clim = t
        ti=which(floor(t$time)>=(cyr-36) & floor(t$time)<=(cyr+35))
        sts=ti[1]
        ets=ti[length(ti)]
        t.clim$data=t(t$data[sts:ets,])
        t.clim$time=t$time[sts:ets]
        # 4.2.2 Transform it to Oct-Sept years
        if (cyr < 1636) {
          y_start = agrep(paste0(1600,'.792'),as.character(t.clim$time))
        } else {
          y_start = agrep(paste0(cyr-36,'.792'),as.character(t.clim$time))
        }
        if (cyr > 1970) {
          y_end =  grep(paste0(2005,'.708'),as.character(t.clim$time))
        } else {
          y_end = grep(paste0(cyr+35,'.708'),as.character(t.clim$time))
        }
        t.clim$data = t.clim$data[,y_start:y_end]
        t.clim$time = t.clim$time[y_start:y_end]
        t.clim$data = apply(array(t.clim$data, c(nrow(t.clim$data), 12, ncol(t.clim$data)/12)), 1:2, mean,na.rm=T)
        # 4.2.3 Cut out the cyr-1 and cyr time window
        t.cyr = t
        ti=which(floor(t$time)>=cyr-1 & floor(t$time)<cyr+1)
        sts=ti[1]
        ets=ti[length(ti)]
        t.cyr$data=t(t$data[sts:ets,])
        t.cyr$time=t$time[sts:ets]
        # 4.2.4 Transform it to Oct-Sept  years
        y_start = agrep(paste0(cyr-1,'.792'),as.character(t.cyr$time))
        y_end = grep(paste0(cyr,'.708'),as.character(t.cyr$time))
        t.cyr$data = t.cyr$data[,y_start:y_end]
        t.cyr$time = t.cyr$time[y_start:y_end]
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
      tmp1 <- array(t(doc_t_mon$data),c(dim(doc_t_mon$data)[1] *  dim(doc_t_mon$data)[2]))
      doc_t_mon$data <- array(tmp1,c(dim(tmp1)[1]/(dim(doc_t_mon$data)[1]/6),2))
      rm(tmp1)
      doc_t_mon$time <- c(cyr,cyr+0.5)
      doc_t_mon$names <- rep(doc_t_mon$names,6)
      doc_t_mon$lon <- rep(doc_t_mon$lon,6)
      doc_t_mon$lat <- rep(doc_t_mon$lat,6)
      doc_t_mon$des=rep('mon',nrow(doc_t_mon$data))
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
      load('../assimil_data/data_yuri/slp.Rdata') # monthly slp collection from yuri, histalp included
      inst_slp <- slp
    }
    if (yuri_temp) {
      load('../assimil_data/data_yuri/t.Rdata') # monthly temp collection from yuri, histalp included
      inst_t <- t
    }
    if (ghcn_temp) {
      load("../assimil_data/ghcn/ghcn_temp_1600-2005.Rdata")
    }
    if (ghcn_prec) {
      load("../assimil_data/ghcn/ghcn_precip_1600-2005.Rdata")
    }
  
    # 5.2 Qualtiy check of the data
    if (check_assimdata) {
      if (ghcn_prec) {
        varlist <- c("inst_t","inst_slp","ghcn","ghcn_precip")
      } else {
        varlist <- c("inst_t","inst_slp","ghcn")
      }
      for (varname in varlist) {
        var <- get(varname)
        gpos <- getgridboxnum(var,echam.sd)
        for (i in 1:length(gpos)) { # order of this and next if statement changed 2017-07-25
          if (!is.na(gpos[i])) {
            m <- gpos[i]
            d <- compute_dist(var$lon[i],var$lat[i],echam.sd$lon[m],echam.sd$lat[m])
            if ((!is.na(d)) & (d > 600)) { # check distance of assim data to next model grid box
              m=NA
              print(paste('inst data', varname, i, '>600km from echam grid box; set to NA'))
              write(paste('inst data', varname, i, '>600km from echam grid box; set to NA'),
                    file=paste0('../log/',logfn),append=T)
              if (varname=="inst_t") {inst_t$data[,i] <- NA}
              if (varname=="inst_slp") {inst_slp$data[,i] <- NA}
              if (varname=="ghcn") {ghcn$data[,i] <- NA}
              if (ghcn_prec){
                if (varname=="ghcn_precip") {ghcn_precip$data[,i] <- NA}
              }
            }
            if (!is.na(m)) {
              # year from Oct to Sept for the climatology
              ti=which(floor(var$time)>=(cyr-36) & floor(var$time)<=(cyr+35))
              sts=ti[1]
              ets=ti[length(ti)]
              var.clim = var
              var.clim$data=t(var$data[sts:ets,])
              var.clim$time=var$time[sts:ets]
              # transform it to Oct-Septyears
              if (cyr < 1636) {
                  y_start = agrep(paste0(1600,'.792'),as.character(var.clim$time))
              } else {
                  y_start = agrep(paste0(cyr-36,'.792'),as.character(var.clim$time))
              }
              if (cyr > 1970) {
                  y_end =  grep(paste0(2005,'.708'),as.character(var.clim$time)) # 2012 should be changed to 2005, everywhere
              } else {
                  y_end = grep(paste0(cyr+35,'.708'),as.character(var.clim$time))
              }
              var.clim$data = var.clim$data[,y_start:y_end]
              var.clim$time = var.clim$time[y_start:y_end]
              var.clim$data = t(var.clim$data)
              vtmp <- array(var.clim$data,c(12, nrow(var.clim$data)/12, dim(var.clim$data)[2]))
              stsv = round(dim(vtmp)[2]/2)

              for (j in 1:12) {
                # if ( !is.na(vtmp[j,((stsv-1)/12+1),i]) ) { # I added this, cause it can happen that the station in a certain month of cyr year has no data -> not good yet
                #if (all(is.na(vtmp[j,,i])) ) { # good
                #  print(paste0("There is no data for month:", j, " at station:",i))
                #} else {
# VERONIKA PLEASE CHECK:  I spent quite a bit of time and think I fixed everything now!
#                  echam.clim_2 = echam.clim
#                  tmp.echam.clim = array(echam.clim_2$ensmean, dim(echam.clim_2$ensmean)[1]*dim(echam.clim_2$ensmean)[2])
#                  echam.clim_2$ensmean = array(tmp.echam.clim, c(dim(tmp.echam.clim)/12, 12))
#                  # Probably at this point Veronika throws out a lot of data that Jörg still uses to assimilate
#                  # Therfore Veronika commented following if
#                  # if ( sum(is.na(vtmp[j,,i])) > 41) { # I added this if, but maybe one if would be enough (all - 41)
#                  #  print("Inst data will not be assimilated at this time step, because vtmp has more than 41 NA values" )
#                  #  biasm = NA
#                  #} else {
#                  biasm <- echam.clim_2$ensmean[m,j] - mean(vtmp[j,,i], na.rm = T)
                
                # if bias corrected proxy/inst is outside echam ens range +- 5SD,
                # data point will not be assimilated at this time step
                if (!all(is.na(vtmp[j,,i])) ) { 
                  biasm <- echam_clim_mon_ensmean[m,j] - mean(vtmp[j,,i],na.rm=T) 
                  if (!is.na(biasm) & !is.na(vtmp[j,stsv,i]) ) {  # Veronika added the second term of the if
                  #  if (((vtmp[j,((stsv-1)/12+1),i]+ biasm) < echam$ensmean[m,(j+12)]-5*echam.sd$data[m,j]) |
                  #     ((vtmp[j,((stsv-1)/12+1),i]+ biasm) > echam$ensmean[m,(j+12)]+5*echam.sd$data[m,j])) {
                  if (((vtmp[j,stsv,i]+biasm) < echam_clim_mon_ensmean[m,j]-5*echam.sd$data[m,j]) |
                       ((vtmp[j,stsv,i]+ biasm) > echam_clim_mon_ensmean[m,j]+5*echam.sd$data[m,j])) {
                      print(paste('inst data',varname,'#:',i,'mon:',j,'out of range'))
                      write(paste('inst data', varname,'#:',i,'mon:',j,'out of range'),
                            file=paste0('../log/',logfn),append=T)
                      write(paste('data lon/lat',var$lon[i],var$lat[i]),
                            file=paste0('../log/',logfn),append=T)
                      write(paste('echam lon/lat', echam.sd$lon[m],echam.sd$lat[m]),
                            file=paste0('../log/',logfn),append=T)
                      write(paste('bias corr. data', vtmp[j,stsv,i]+biasm),
                            file=paste0('../log/',logfn),append=T)
                      write(paste('echam clim', echam_clim_mon_ensmean[m,j]),
                            file=paste0('../log/',logfn),append=T)
                      write(paste('echam sd', echam.sd$data[m,j]),
                            file=paste0('../log/',logfn),append=T)
                      if (varname=="inst_t") {inst_t$data[,i] <- NA}
                      if (varname=="inst_slp") {inst_slp$data[,i] <- NA}
                      if (varname=="ghcn") {ghcn$data[,i] <- NA}
                      if (ghcn_prec){
                        if (varname=="ghcn_precip") {ghcn_precip$data[,i] <- NA}
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    # 5.3 Calculating the anomalies and transforming them to Oct-Sept time period
    if ((ghcn_temp) & (dim(ghcn$data)[2]>0)) {
      if (anomaly_assim & sixmonstatevector){
        # calculate the climatology
        ghcn.clim=ghcn
        ti=which(floor(ghcn$time)>=(cyr-36) & floor(ghcn$time)<=(cyr+35))
        sts=ti[1]
        ets=ti[length(ti)]
        ghcn.clim$data=t(ghcn$data[sts:ets,])
        ghcn.clim$time=ghcn$time[sts:ets]
        # transform it to Oct-Sept  years
        if (cyr < 1636) {
            y_start = agrep(paste0(1600,'.792'),as.character(ghcn.clim$time))
        } else {
            y_start = agrep(paste0(cyr-36,'.792'),as.character(ghcn.clim$time))
        }
        if (cyr > 1970) {
            y_end =  grep(paste0(2005,'.708'),as.character(ghcn.clim$time))
        } else {
            y_end = grep(paste0(cyr+35,'.708'),as.character(ghcn.clim$time))
        }
        ghcn.clim$data = ghcn.clim$data[,y_start:y_end]
        ghcn.clim$time = ghcn.clim$time[y_start:y_end]
        ghcn.clim$data = apply(array(ghcn.clim$data, c(nrow(ghcn.clim$data), 12, ncol(ghcn.clim$data)/12)), 1:2, mean,na.rm=T)
        
        # cut out the cyr-1 and cyr time window
        ghcn.cyr = ghcn
        ti=which(floor(ghcn$time)>=cyr-1 & floor(ghcn$time)<cyr+1)
        sts=ti[1]
        ets=ti[length(ti)]
        ghcn.cyr$data=t(ghcn$data[sts:ets, ])
        ghcn.cyr$time=ghcn$time[sts:ets]
        # transform it to Oct-Sept  years
        y_start = agrep(paste0(cyr-1,'.792'),as.character(ghcn.cyr$time))
        y_end = grep(paste0(cyr,'.708'),as.character(ghcn.cyr$time))
        ghcn.cyr$data = ghcn.cyr$data[,y_start:y_end]
        ghcn.cyr$time = ghcn.cyr$time[y_start:y_end]
        
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
        inst_slp.clim=inst_slp
        ti=which(floor(inst_slp$time)>=(cyr-36) & floor(inst_slp$time)<=(cyr+35))
        sts=ti[1]
        ets=ti[length(ti)]
        inst_slp.clim$data=t(inst_slp$data[sts:ets,])
        inst_slp.clim$time=inst_slp$time[sts:ets]
        # transform it to Oct-Sept  years
        if (cyr < 1636) {
          y_start = agrep(paste0(1600,'.792'),as.character(inst_slp.clim$time))
        } else {
          y_start = agrep(paste0(cyr-36,'.792'),as.character(inst_slp.clim$time))
        }
        if (cyr > 1970) {
          y_end =  grep(paste0(2005,'.708'),as.character(inst_slp.clim$time))
        } else {
          y_end = grep(paste0(cyr+35,'.708'),as.character(inst_slp.clim$time))
        }
        inst_slp.clim$data = inst_slp.clim$data[,y_start:y_end]
        inst_slp.clim$time = inst_slp.clim$time[y_start:y_end]
        inst_slp.clim$data = apply(array(inst_slp.clim$data, c(nrow(inst_slp.clim$data), 12, ncol(inst_slp.clim$data)/12)), 1:2, mean,na.rm=T)
        # cut out the cyr-1 and cyr time window
        inst_slp.cyr = inst_slp
        ti=which(floor(inst_slp$time)>=cyr-1 & floor(inst_slp$time)<cyr+1)
        sts=ti[1]
        ets=ti[length(ti)]
        inst_slp.cyr$data=t(inst_slp$data[sts:ets, ])
        inst_slp.cyr$time=inst_slp$time[sts:ets]
        # transform it to Oct-Sept  years
        y_start = agrep(paste0(cyr-1,'.792'),as.character(inst_slp.cyr$time))
        y_end = grep(paste0(cyr,'.708'),as.character(inst_slp.cyr$time))
        inst_slp.cyr$data = inst_slp.cyr$data[,y_start:y_end]
        inst_slp.cyr$time = inst_slp.cyr$time[y_start:y_end]
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
        inst_t.clim=inst_t
        ti=which(floor(inst_t$time)>=(cyr-36) & floor(inst_t$time)<=(cyr+35))
        sts=ti[1]
        ets=ti[length(ti)]
        inst_t.clim$data=t(inst_t$data[sts:ets,])
        inst_t.clim$time=inst_t$time[sts:ets]
        # transform it to Oct-Sept  years
        if (cyr < 1636) {
          y_start = agrep(paste0(1600,'.792'),as.character(inst_t.clim$time))
        } else {
          y_start = agrep(paste0(cyr-36,'.792'),as.character(inst_t.clim$time))
        }
        if (cyr > 1970) {
          y_end =  grep(paste0(2005,'.708'),as.character(inst_t.clim$time))
        } else {
          y_end = grep(paste0(cyr+35,'.708'),as.character(inst_t.clim$time))
        }
        inst_t.clim$data = inst_t.clim$data[,y_start:y_end]
        inst_t.clim$time = inst_t.clim$time[y_start:y_end]
        inst_t.clim$data = apply(array(inst_t.clim$data, c(nrow(inst_t.clim$data), 12, ncol(inst_t.clim$data)/12)), 1:2, mean,na.rm=T)
        # cut out the cyr-1 and cyr time window
        inst_t.cyr = inst_t
        ti=which(floor(inst_t$time)>=cyr-1 & floor(inst_t$time)<cyr+1)
        sts=ti[1]
        ets=ti[length(ti)]
        inst_t.cyr$data=t(inst_t$data[sts:ets, ])
        inst_t.cyr$time=inst_t$time[sts:ets]
        # transform it to Oct-Sept  years
        y_start = agrep(paste0(cyr-1,'.792'),as.character(inst_t.cyr$time))
        y_end = grep(paste0(cyr,'.708'),as.character(inst_t.cyr$time))
        inst_t.cyr$data = inst_t.cyr$data[,y_start:y_end]
        inst_t.cyr$time = inst_t.cyr$time[y_start:y_end]
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
            dlist[i]=m   #[1]
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
      docall.allts$sour <- rep('doc',length(docall$lon))    
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
      Hcal1 <- array(NA,dim=c(dim(inst$data)[1],2))
      Hcal1 <- compute_Hi_Hredux_sixmonstatevector(inst, etmp, threshold=700)
    }
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
        if (!landcorr) {
          if ((!is.na(dist[1])) & (!is.na(wgt[1]))){
            for (i in 1:ntim){
              if (!is.na(calibrate$data[j,i])) {
                x2tmp <- analysis$data[,i,] # entire state vector at time step i, all members
                x2 <- x2tmp[h.i,,drop=F] # state vector at time step i and h.i, all members
                PH <- (analysis$data[,i,] %*% t(x2) / (nens - 1) * wgt) %*% t(H)
                if (covarclim>0 & covarclim<100) { 
                  x2climtmp <- echanomallts$data[,i,] 
                  x2clim <- x2climtmp[h.i,,drop=F] 
                  PHclim <- (echanomallts$data[,i,] %*% t(x2clim) / 
                              ((dim(echanomallts$data)[3]) - 1) * wgt) %*% t(H)
                  PH <- (PH*(1-(covarclim/100))) + (PHclim*(covarclim/100))
                } else if (covarclim==100) {
                  x2climtmp <- echanomallts$data[,i,] 
                  x2clim <- x2climtmp[h.i,,drop=F] 
                  PH <- (echanomallts$data[,i,] %*% t(x2clim) / 
                          ((dim(echanomallts$data)[3]) - 1) * wgt) %*% t(H) 
                } 
                HPHR <- as.vector(H %*% PH[h.i,] + Rcal[j])
                K <- PH / HPHR
                Ktilde <- K / (1 + sqrt(Rcal[j]/HPHR))
                analysis$ensmean[,i] <- analysis$ensmean[,i] + K[,1] * (calibrate$data[j,i] -
                                          H %*% analysis$ensmean[h.i,i])
                analysis$data[,i,] <- analysis$data[,i,] - Ktilde %*% H %*% analysis$data[h.i,i,]
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
                corland_analysis$data[,i,] <- corland_analysis$data[,i,] + K[,1] * 
                  (calibrate$data[j,i] - H %*% corland_analysis$data[h.i,i,])
              }
            }
          } 
        } 
      }
      ## add ensemble mean analysis back 
      if (!landcorr) {
        analysis$data <- analysis$data + as.vector(analysis$ensmean)
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
            save(analysis.anom,analysis.abs,echam.anom,echam.abs,validate,calibrate,
              file=paste0(dataintdir,'analysis/',expname,'/analysis_',cyr,'_2ndgrid.Rdata'))
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
               file=paste0('analysis/analysis_',cyr,'.Rdata'))
        } else {
          save(analysis.anom,analysis.abs,echam.anom,echam.abs,landcorrected.anom,
               landcorrected.clim,calibrate,
               file=paste0('analysis/analysis_',cyr,'.Rdata'))
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
