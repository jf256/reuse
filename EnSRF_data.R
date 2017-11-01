rm(list=ls())
machine="climcal" #"macbook" #"climcal3" # "climcal"
if (machine=="macbook") {
  # analysis start from syr+1 if sixmonths statevector
  # use 1659 as eyr because first instr. data available
  # use 1749 as eyr becasue Luterbacher vali data avail.
  # use 1854 as syr when no docu are available anymore
  # use 1901 as syr because of cru validation data start
  # use 1961 as syr when no proxies are assimilated anymore
  syr=1603     # set to 1603 to process analysis; set to >=1901 or 1902? for cru validation
  eyr=2004     # max 2004 ?
  workdir='~/unibe/projects/EnSRF/src/'
  dataextdir='/Volumes/DATA/climdata/EKF400/'
  dataintdir=paste0(workdir,'../data/')
} else {
  args <- commandArgs(TRUE)
  syr = as.numeric(args[1])
  eyr = as.numeric(args[2])
  print(paste('period',syr, 'to', eyr))
  workdir='/scratch3/joerg/projects/reuse/src/'
  dataextdir='/mnt/climstor/giub/EKF400/'
  dataintdir=paste0(workdir,'../data/')
}
setwd(workdir)

source('EnSRF_switches.R')
source('EnSRF_functions.R')

dir.create(paste0("../data/analysis/",expname))
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
  # read echam 70yr anom, clim and sd calculated with cdo from .nc files to .RData
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
  if ((!eyr<1948) & (!syr>2009)) {
    #if (syr<1948){syr_ncep=1948} else {syr_ncep=syr}
    #if (eyr>2009){eyr_ncep=2009} else {eyr_ncep=eyr}
    #syr_ncep=1948
    #eyr_ncep=2009
    #  see script in EnSRF/script/merge_ncep.sh for regridding and co of orig. ncep data set
    ncepall <- read_echam1('ncep_allvar_1948-2009',timlim=c(syr_ncep,eyr_ncep),
                            path=nceppath,small=every2grid)
    if (every2grid) {
      save(ncepall, file=paste0("../data/ncep/ncep_allvar_",syr_ncep,"-",eyr_ncep,"_2ndgrid.Rdata"))
    } else {
      save(ncepall, file=paste0("../data/ncep/ncep_allvar_",syr_ncep,"-",eyr_ncep,".Rdata"))
    }
  } 
}

if (generate_CRUALLVAR) {
  print("generate_CRUALLVAR")
#  if ((!eyr<1901) & (!syr>2004)) {
#    if (syr<1901){syr_cru=1901} else {syr_cru=syr}
#    syr_cru=1901
#    if (eyr>2004){eyr_cru=2004} else {eyr_cru=eyr}
#    eyr_cru=2004
#  see script in EnSRF/script/merge_cru.sh for regridding and co of orig. cru data set
    cruall <- read_echam1('cru_allvar_abs_1901-2004.nc',timlim=c(syr_cru,eyr_cru),
                           path=crupath,small=every2grid,landonly=land_only)
    if (every2grid) {
      save(cruall, file=paste0(dataintdir,"cru/cru_allvar_",syr_cru,"-",eyr_cru,"_2ndgrid.Rdata"))
    } else {
      save(cruall, file=paste0(dataintdir,"cru/cru_allvar_",syr_cru,"-",eyr_cru,".Rdata"))  
    }
#  }
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
  #if (syr<1750){syr_recon=1750} else {syr_recon=syr}
  #if (eyr>1999){eyr_recon=1999} else {eyr_recon=eyr}
  #syr_recon=1750
  #eyr_recon=1900
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
    # add NA for TRW data that ends 1970
    trwprox$data <- rbind(trwprox$data,
                      matrix(NA,nrow=length(seq(1971,2004)),ncol=dim(trwprox$data)[2]))
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
# loading echam timeslice (loop over years to reduce size of state vector)
##########################################################################################
if (sixmonstatevector) {syr2=syr+1} else {syr2=syr}
for (cyr in syr2:eyr) {
  print(cyr)
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
  print(paste("instr:",instrumental, "proxies:",real_proxies, "documentary:",docum))
  asyr <- cyr-35 # "a" for anomaly
  if (asyr < 1601) {asyr = 1601}
  aeyr <- cyr+35
  if (aeyr > 2005) {aeyr = 2005}
  ptm1 <- proc.time()
  if ((cyr==syr2) & (covarclim>0)) {
    load(file="../data/echam/echallts_for_covar.Rdata")
    #just use every 3. data point to make calculation faster
    echanomallts$data <- echanomallts$data[,,seq(1,dim(echanomallts$data)[3],n_covar)] 
  }
  if (every2grid) {
    load(paste(dataextdir,"echam/echam_",(cyr-1),"-",cyr,"_2ndgrid.Rdata",sep=""))
    #load(paste(dataintdir,"echam/echam_",(cyr-1),"-",cyr,"_2ndgrid.Rdata",sep=""))
  } else {
    load(paste(dataextdir,"echam/echam_",(cyr-1),"-",cyr,".Rdata",sep=""))
    #load(paste(dataintdir,"echam/echam_",(cyr-1),"-",cyr,".Rdata",sep=""))
  }
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
  if (instrumental){
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
  }
  
  if (real_proxies){
    load(paste0("../data/proxies/real_proxies_",fsyr,"-",feyr,".Rdata"))  
  }

  if (docum){
    load('../assimil_data/data_yuri/t_docu_monthly.Rdata') 
    doc_t_mon <- t
    if (!any(!is.na(doc_t_mon$data))) { docum=F }
# next section for seasonally resolved ducumentary data is commented 
# because we only read monthly docum data at the moment    
#     load('../data_yuri/t_docu_seas.Rdata') 
#     doc_t_seas <- t
#     # cut time period as not done yet
#     ti=which(floor(doc_t_seas$time)>=syr & floor(doc_t_seas$time)<=eyr)
#     sts=ti[1]
#     ets=ti[length(ti)]
#     doc_t_seas$data=doc_t_seas$data[sts:ets,]
#     doc_t_seas$time=doc_t_seas$time[sts:ets]
#     load('../data_yuri/t_docu_JFMA.Rdata') 
#     doc_t_JFMA <- t
#     ti=which(floor(doc_t_JFMA$time)>=syr & floor(doc_t_JFMA$time)<=eyr)
#     sts=ti[1]
#     ets=ti[length(ti)]
#     doc_t_JFMA$data=doc_t_JFMA$data[sts:ets,]      
#     doc_t_JFMA$time=doc_t_JFMA$time[sts:ets]
#     load('../data_yuri/t_docu_AMJJA.Rdata') 
#     doc_t_AMJJA <- t
#     ti=which(floor(doc_t_AMJJA$time)>=syr & floor(doc_t_AMJJA$time)<=eyr)
#     sts=ti[1]
#     ets=ti[length(ti)]
#     doc_t_AMJJA$data=doc_t_AMJJA$data[sts:ets,]       
#     doc_t_AMJJA$time=doc_t_AMJJA$time[sts:ets]
  }

  print('calc time for loading data')
  print(proc.time() - ptm1)




  
  
  
  
  
# just leave temp precip slp in state vector  
  if (tps_only) {
    tpspos <- c(which(echam$names=='temp2'), which(echam$names=='precip'), 
      which(echam$names=='slp'), which(echam$names=='bias'))
    echam$data <- echam$data[tpspos,,]
    echam$ensmean <- echam$ensmean[tpspos,]
    echam$names <- echam$names[tpspos]
    if ((cyr==syr2) & (covarclim>0)) {
      tpspos <- c(which(echanomallts$names=='temp2'), which(echanomallts$names=='precip'), 
                  which(echanomallts$names=='slp'), which(echanomallts$names=='bias'))
      echanomallts$data <- echanomallts$data[tpspos,,]
      echanomallts$ensmean <- echanomallts$ensmean[tpspos,]
      echanomallts$names <- echanomallts$names[tpspos]
    }
    if (vali) {
      tpspos2 <- c(which(valiall$names=='temp2'), which(valiall$names=='precip'), 
                   which(valiall$names=='slp'))
      valiall$data <- valiall$data[tpspos2,]
      valiall$names <- valiall$names[tpspos2]
    }
  } else if (no_stream) {
    # ACHTUNG stream var has ERROR because the 5/9 levels before/after 1880 have a lat dimension
    tpspos <- c(which(echam$names!='stream'))
    echam$data <- echam$data[tpspos,,]
    echam$ensmean <- echam$ensmean[tpspos,]
    echam$names <- echam$names[tpspos]
  }
  if (fasttest) {
    mulc <- 4 # choose every 4th grid box
    loi <- seq(1:length(echam$lon))
    lai <- seq(1:length(echam$lat))
    di <- seq(1:dim(echam$data)[1])
    ni <- seq(1:length(echam$names))
    loi <- loi[seq(ceiling(mulc/2),length(loi),mulc)]
    lai <- lai[seq(ceiling(mulc/2), length(lai),mulc)]
    di <- di[seq(ceiling(mulc/2), length(di),mulc)]
    ni <- ni[seq(ceiling(mulc/2), length(ni),mulc)]
    echam$lon <- echam$lon[loi]
    echam$lat <- echam$lat[lai]
    echam$data <- echam$data[di,,]
    echam$ensmean <- echam$ensmean[di,]
    echam$names <- echam$names[ni]
    if (vali) {
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
  }

  if (anomaly_assim){
    if (load_71yr_anom) { # anomalies calculated efficiently with cdo 
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
      } else {
        load(paste0(echanompath,'echam_anom_',yr1,'-',yr2,'.Rdata'))
        load(paste0(echclimpath,'echam_clim_',yr3,'-',yr4,'.Rdata'))
      }
      if (no_forc_big_ens) {
        yrs <- floor(runif(n_no_forc,1602,2004))
        m <- floor(runif(n_no_forc,1,30))
        for (n in 1:n_no_forc) {
          yr1 <- yrs[n]
          yr2 <- yr1+1
          if (every2grid) {
            load(paste0(echanompath,'echam_anom_',yr1,'-',yr2,'_2ndgrid.Rdata'))
            #  load(paste0(echclimpath,'echam_clim_',yr3,'-',yr4,'_2ndgrid.Rdata'))
          } else {
            load(paste0(echanompath,'echam_anom_',yr1,'-',yr2,'.Rdata'))
            #  load(paste0(echclimpath,'echam_clim_',yr3,'-',yr4,'.Rdata'))
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
      
      # just leave temp precip slp in state vector  
      if (tps_only) {
        tpspos <- c(which(echam_anom$names=='temp2'), which(echam_anom$names=='precip'), 
                    which(echam_anom$names=='slp'), which(echam_anom$names=='bias'))
        echam_anom$data <- echam_anom$data[tpspos,,]
        echam_anom$ensmean <- echam_anom$ensmean[tpspos,]
        echam_anom$names <- echam_anom$names[tpspos]
        echam_clim$data <- echam_clim$data[tpspos,,]
        echam_clim$ensmean <- echam_clim$ensmean[tpspos,]
        echam_clim$names <- echam_clim$names[tpspos]
      } else if (no_stream) {
        # ACHTUNG stream var has ERROR because the 5/9 levels before/after 1880 have a lat dimension
        tpspos <- c(which(echam_anom$names!='stream'))
        echam_anom$data <- echam_anom$data[tpspos,,]
        echam_anom$ensmean <- echam_anom$ensmean[tpspos,]
        echam_anom$names <- echam_anom$names[tpspos]
        echam_clim$data <- echam_clim$data[tpspos,,]
        echam_clim$ensmean <- echam_clim$ensmean[tpspos,]
        echam_clim$names <- echam_clim$names[tpspos]
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
      #new echam_clim made by veronika have correct units, hence comment following corrections
      #error where echam_anom and echam_clim were generated: no unit correction happened, thus here
      echam_anom$data[echam_anom$names=='precip',,] <- echam_anom$data[echam_anom$names=='precip',,]*
                                                         3600 * 24 * 30
      echam_anom$ensmean[echam_anom$names=='precip',] <- echam_anom$ensmean[echam_anom$names=='precip',]*                                                            3600 * 24 * 30
      echam_anom$data[echam_anom$names=='slp',,] <- echam_anom$data[echam_anom$names=='slp',,]/100
      echam_anom$ensmean[echam_anom$names=='slp',] <- echam_anom$ensmean[echam_anom$names=='slp',]/100
      echam.anom <- echam_anom
      # echam_clim$data[echam_clim$names=='temp2',,] <- echam_clim$data[echam_clim$names=='temp2',,]- 
      #                                                   273.15
      # echam_clim$ensmean[echam_clim$names=='temp2',] <- echam_clim$ensmean[echam_clim$names=='temp2',]-
      #                                                     273.15
      # echam_clim$data[echam_clim$names=='precip',,] <- echam_clim$data[echam_clim$names=='precip',,]* 
      #                                                    3600 * 24 * 30
      # echam_clim$ensmean[echam_clim$names=='precip',] <- echam_clim$ensmean[echam_clim$names=='precip',]*
      #                                                      3600 * 24 * 30
      # echam_clim$data[echam_clim$names=='slp',,] <- echam_clim$data[echam_clim$names=='slp',,]/100 
      # echam_clim$ensmean[echam_clim$names=='slp',] <- echam_clim$ensmean[echam_clim$names=='slp',]/100 
      echam.clim <- echam_clim
      rm(echam_anom,echam_clim)
    } else if (anom_reload) { # anomalies calulated with following slow r code
      load(file=paste0("../data/anom/EnSRF_anom_",cyr,".Rdata")) 
    } else { # calulate anomalies with r: ATTENTION super slow!!!
    # calc running 71yr anomalies and climatology
    # 2. dimension contains e.g. 71 years of 12 monthly data 
      echam.backup <- echam
      if (machine == "climcal3") {
        etd <- array(0,dim=c(dim(echam$data)[1],((aeyr-asyr+1)*nseas),
                          dim(echam$data)[3]))
      } else {
        # not enough memory, thus use harddrive -> SLOW
        etd <- ff(0,dim=c(dim(echam$data)[1],((aeyr-asyr+1)*nseas),
                               dim(echam$data)[3]))
      }
      ete <- array(NA,dim=c(dim(echam$ensmean)[1],((aeyr-asyr+1)*nseas)))
      for (i in seq(0,(aeyr-asyr-1),by=1)) {
        yr1 <- asyr+i
        yr2 <- asyr+i+1
        if (every2grid) {
          load(paste0(dataextdir,'echam/echam_',yr1,'-',yr2,'_2ndgrid.Rdata'))
        } else {
          load(paste0(dataextdir,'echam/echam_',yr1,'-',yr2,'.Rdata'))
        }
        if (tps_only) {
          tpspos <- c(which(echam$names=='temp2'), which(echam$names=='precip'), 
                      which(echam$names=='slp'), which(echam$names=='bias'))
          echam$data <- echam$data[tpspos,,]
          echam$ensmean <- echam$ensmean[tpspos,]
          echam$names <- echam$names[tpspos]
        } else if (no_stream) {
          tpspos <- c(which(echam$names!='stream'))
          echam$data <- echam$data[tpspos,,]
          echam$ensmean <- echam$ensmean[tpspos,]
          echam$names <- echam$names[tpspos]
        } 
        if (fasttest) {
          loi <- seq(1:length(echam$lon))
          lai <- seq(1:length(echam$lat))
          di <- seq(1:dim(echam$data)[1])
          ni <- seq(1:length(echam$names))
          loi <- loi[seq(ceiling(mulc/2),length(loi),mulc)]
          lai <- lai[seq(ceiling(mulc/2), length(lai),mulc)]
          di <- di[seq(ceiling(mulc/2), length(di),mulc)]
          ni <- ni[seq(ceiling(mulc/2), length(ni),mulc)]
          echam$lon <- echam$lon[loi]
          echam$lat <- echam$lat[lai]
          echam$data <- echam$data[di,,]
          echam$ensmean <- echam$ensmean[di,]
          echam$names <- echam$names[ni]
        }
        if (i == (aeyr-asyr-1)) {
          etd[,((i*nseas+1):((i*nseas)+(2*nseas))),] <- echam$data[,1:(2*nseas),]
          ete[,((i*nseas+1):((i*nseas)+(2*nseas)))] <- echam$ensmean[,1:(2*nseas)]
        } else {
          etd[,((i*nseas+1):(i*nseas+nseas)),] <- echam$data[,1:nseas,]
          ete[,((i*nseas+1):(i*nseas+nseas))] <- echam$ensmean[,1:nseas]
        }
      } 
      nclim <- dim(etd)[2]/nseas # climatology is the 71 year mean of echam data
      nclime <- dim(ete)[2]/nseas # climatology is the 71 year mean of echam ensmean
    # calc climatology
      emn <- array(0, dim(etd[,1:nseas,]))
      emne <- array(0, dim(ete[,1:nseas]))
      for (i in 1:nclim){
        emn <- emn + etd[,(1:nseas + (i - 1)*nseas),]/nclim
        emne <- emne + ete[,(1:nseas + (i - 1)*nseas)]/nclime
      }
      eanom <- array(0, dim(etd))
      eanome <- array(0, dim(ete))
      for (i in 1:nclim){
        eanom[,1:nseas + (i - 1)*nseas,] <- etd[,1:nseas + (i - 1)*nseas,] - emn
        eanome[,1:nseas + (i - 1)*nseas] <- ete[,1:nseas + (i - 1)*nseas] - emne
      }
      echam.anom <- echam
      echam.anom$data <- eanom
      echam.anom$ensmean <- eanome
      echam.anom$time <- seq(asyr,aeyr+11/12,by=1/12)
      echam.clim <- echam
      echam.clim$data <- emn
      echam.clim$ensmean <- emne
      echam.clim$time <- echam.anom$time
      echam <- echam.backup
      rm(echam.backup) 
      rm(etd);rm(emn);rm(emne);rm(ete);rm(eanom);rm(eanome) 
      if (anom_save) {
        save(echam.anom,echam.clim,file=paste0("../data/anom/EnSRF_anom_",cyr,".Rdata"))
      }
    }
  }
  print('calc time for echam anomalies')
  print(proc.time() - ptm1)

  
  
  
# calc echam st. dev. for each grid point and month over ens memb. to scale docu data
  echam.sd <- apply(echam$data[,13:24,],1:2,sd)
  print('calc time for standard deviations')
  print(proc.time() - ptm1)

  
  
  
  

#########################################################################################
# screen proxy/instr. assimilation data 
#########################################################################################
  if ((instrumental) & (check_assimdata)) {
    if (ghcn_prec) {
      varlist <- c("inst_t","inst_slp","ghcn","ghcn_precip")
    } else {
      varlist <- c("inst_t","inst_slp","ghcn")
    }
    for (varname in varlist) {
      var <- get(varname)
      gpos <- getgridboxnum(var,echam)
      for (i in 1:length(gpos)) { # order of this and next if statement changed 2017-07-25
        if (!is.na(gpos[i])) {
          m <- gpos[i]
          d <- compute_dist(var$lon[i],var$lat[i],echam$lon[m],echam$lat[m])
          if ((!is.na(d)) & (d > 600)) { # check distance of assim data to next model grid box
            m=NA
            print(paste('inst data', varname, i, '>600km from echam grid box; set to NA'))
            if (varname=="inst_t") {inst_t$data[,i] <- NA}
            if (varname=="inst_slp") {inst_slp$data[,i] <- NA}
            if (varname=="ghcn") {ghcn$data[,i] <- NA}
            if (ghcn_prec){
              if (varname=="ghcn_precip") {ghcn_precip$data[,i] <- NA}
            }
          }
          if (!is.na(m)) {
            tiv=which(floor(var$time)==cyr)
            stsv=tiv[1]
            vtmp <- array(var$data,c(12, nrow(var$data)/12, dim(var$data)[2]))
            if (cyr<(floor(var$time[1])+35)) {
              vper=(1:(((stsv-1)/12+1)+35))
            } else if (cyr>(floor(var$time[length(var$time)])-35)){
              vper=((((stsv-1)/12+1)-35):dim(vtmp)[2])
            } else{
              vper=((((stsv-1)/12+1)-35):(((stsv-1)/12+1)+35))
            }  
            for (j in 1:12) {
              # if bias corrected proxy/inst is outside echam ens range +- 5SD, 
              # data point will not be assimilated at this time step
              # vtmp should be 71 yr period and not just one year
              biasm <- echam.clim$ensmean[m,j] - 
                mean(vtmp[j,vper,i])  
              #biasm <- echam.clim$ensmean[m,j] - mean(vtmp[j,((stsv-1)/12+1),i])            
              if (!is.na(biasm)) {
                #print("biasm not NA")
                if (((vtmp[j,((stsv-1)/12+1),i]+ biasm) < echam$ensmean[m,(j+12)]-5*echam.sd[m,j]) | 
                   ((vtmp[j,((stsv-1)/12+1),i]+ biasm) > echam$ensmean[m,(j+12)]+5*echam.sd[m,j])) {
                  print(paste('inst data', varname, i, j, 'out of range'))
                  if (cyr == syr2) {
                    write(paste('inst data', varname, i, j, 'out of range'),
                      file=paste0('../log/screening_instr_',cyr,'.log'),append=F)
                  } else {
                    write(paste('inst data', varname, i, j, 'out of range'),
                      file=paste0('../log/screening_instr_',cyr,'.log'),append=T)
                  }
                  write(paste('data lon/lat',var$lon[i],var$lat[i]),
                      file=paste0('../log/screening_instr_',cyr,'.log'),append=T)
                  write(paste('echam lon/lat', echam$lon[m],echam$lat[m]),
                       file=paste0('../log/screening_instr_',cyr,'.log'),append=T)
                  write(paste('bias corr. data', vtmp[j,((stsv-1)/12+1),i]+ biasm),
                       file=paste0('../log/screening_instr_',cyr,'.log'),append=T)
                  write(paste('echam mean', echam$ensmean[m,(j+12)]),
                       file=paste0('../log/screening_instr_',cyr,'.log'),append=T)
                  write(paste('echam sd', echam.sd[m,j]),
                       file=paste0('../log/screening_instr_',cyr,'.log'),append=T)
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


  # no check for docu data at the moment as good quality already manually checked
  
  if ((real_proxies) & (check_assimdata)) {
  # correlation screening already where multiple regression coefficients are calculated
  # thus, screen if value at current time step is more than 5 std. dev. from mean
  # in this case treated as outlier and set to NA
    for (i in 1:length(realprox$lon)) {
      tiv=which(floor(realprox$time)==cyr)
      rpmean <- mean(realprox$data[,i],na.rm=T)
      rpsd   <- sd(realprox$data[,i],na.rm=T)
      if ((!is.na(rpmean)) & (!is.na(rpsd))) {
        if ((!is.na(realprox$data[tiv,i])) & ((realprox$data[tiv,i] < rpmean-5*rpsd) 
                                           | (realprox$data[tiv,i] > rpmean+5*rpsd))) {
          realprox$data[,i] <- NA
          print(paste('proxy data', i, 'out of range'))
          if (cyr == syr2) {
            write(paste('proxy data', i, 'out of range'),file='../log/screening_proxies.log',append=F)
          } else {
            write(paste('proxy data', i, 'out of range'),file='../log/screening_proxies.log',append=T)
          }
        }  
      }
    }
  }

  print('calc time for screening proxies')
  print(proc.time() - ptm1)

  
  


# just leave data for one year (max 12 months) and correct time resolution in memory
  if (sixmonstatevector) {
    # change array to have 6 months in state vector for winter and summer
    # first winter starts in oct of syr
    # 6 mon stat vectors for oct-mar and apr and sep
    if (!anomaly_assim) {
      tmp1 <- array(echam$data,c(dim(echam$data)[1]*dim(echam$data)[2],dim(echam$data)[3]))
      # cut oct syr to sep eyr
      tmp2 <- tmp1[((9*dim(echam$data)[1]+1):(dim(tmp1)[1]-(3*dim(echam$data)[1]))),] 
      echam$data <- array(tmp2,c(dim(tmp2)[1]/(((dim(echam$data)[2]/12)-1)*2),
                                 (((dim(echam$data)[2]/12)-1)*2),dim(tmp2)[2])) 
      # reconvert to 2 seasons per year
      tmp11 <- array(echam$ensmean,c(dim(echam$ensmean)[1]*dim(echam$ensmean)[2]))
      # cut oct syr to sep eyr
      tmp12 <- tmp11[((9*dim(echam$ensmean)[1]+1):(dim(tmp11)[1]-(3*dim(echam$ensmean)[1])))] 
      # reconvert to 2 seasons per year
      echam$ensmean <- array(tmp12,c(dim(tmp12)[1]/(((dim(echam$ensmean)[2]/12)-1)*2),
                          (((dim(echam$ensmean)[2]/12)-1)*2))) 
      echam$time <- seq(cyr,cyr+1.5,0.5) 
      rm(tmp1);rm(tmp2);rm(tmp11);rm(tmp12)
    }
    if (anomaly_assim) {
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
      if (load_71yr_anom) { 
        echam.anom$time <- c(cyr,cyr+0.5) 
      } else {  
        echam.anom$time <- seq(asyr+1,aeyr+0.5,0.5) 
      }
      if (load_71yr_anom) { 
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
      } else { 
        stop("ERROR: load_71yr_anom set to FALSE")
        # tmp61 <- array(echam.clim$data,c(dim(echam.clim$data)[1]*dim(echam.clim$data)[2],
        #                                dim(echam.clim$data)[3]))
        # # Special because climatology is only 1 year lang -> reorder months
        # # cut oct syr to sep eyr
        # tmp62 <- tmp61[((9*dim(echam.clim$data)[1]+1):(dim(tmp61)[1])),] 
        # tmp63 <- tmp61[(1:(9*dim(echam.clim$data)[1])),]
        # tmp64 <- rbind(tmp62,tmp63)
        # echam.clim$data <- array(tmp64,c(dim(tmp64)[1]/(((dim(echam.clim$data)[2]/12))*2),
        #                                (((dim(echam.clim$data)[2]/12))*2),dim(tmp64)[2])) 
        # tmp71 <- array(echam.clim$ensmean,c(dim(echam.clim$ensmean)[1]*dim(echam.clim$ensmean)[2]))
        # tmp72 <- tmp71[((9*dim(echam.clim$ensmean)[1]+1):(dim(tmp71)[1]))] # cut oct syr to sep eyr
        # tmp73 <- tmp71[(1:(9*dim(echam.clim$ensmean)[1]))]
        # tmp74 <- c(tmp72,tmp73)
        # echam.clim$ensmean <- array(tmp74,c(length(tmp74)/(((dim(echam.clim$ensmean)[2]/12))*2),
        #                                   (((dim(echam.clim$ensmean)[2]/12))*2)))  
        # echam.clim$time <- seq(syr+1,eyr+0.5,0.5) 
      }
      echam.clim$names <- rep(echam.clim$names,6)
      rm(tmp41);rm(tmp42);rm(tmp51);rm(tmp52);rm(tmp61);rm(tmp62);rm(tmp71);rm(tmp72)
      if (!load_71yr_anom) {rm(tmp63);rm(tmp64);rm(tmp73);rm(tmp74)}
    }
    if (!recon_vali & vali) {
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
    } else if (vali) {
      # only keep winter and summer season from luterbacher and co and remove spring and autumn
      pos <- sort(c(agrep('.042',as.character(valiall$time)), agrep('.542',as.character(valiall$time))))
      valiall$time <- round(valiall$time[pos],digits=1)
      valiall$data <- valiall$data[,pos]      
    }
  }
  
  # cut one year time slice  
  if (anomaly_assim){
    if (load_71yr_anom) {
      echam=echam.anom
      rm(echam.anom)
    } else {  
      ti=which(floor(echam.anom$time)==cyr)
      sts=ti[1]
      ets=ti[length(ti)]
      echam$data=echam.anom$data[,sts:ets,]
      echam$ensmean=echam.anom$ensmean[,sts:ets]
      echam$time=echam.anom$time[sts:ets]
      rm(echam.anom)
    }
  } else {  
    ti=which(floor(echam$time)==cyr)
    sts=ti[1]
    ets=ti[length(ti)]
    echam$data=echam$data[,sts:ets,]
    echam$ensmean=echam$ensmean[,sts:ets]
    echam$time=echam$time[sts:ets]    
  }   
  
  if (vali) {
    valiall.allts=valiall
    if (cru_vali) {
      if (cyr>1901 && cyr<2005){
        ti=which(floor(valiall$time)==cyr)
        sts=ti[1]
        ets=ti[length(ti)]
        valiall$data=valiall$data[,sts:ets]
        valiall$time=valiall$time[sts:ets]
      } else {
        ti=which(floor(valiall$time)==cyr)
        sts=ti[1]
        ets=ti[length(ti)]
        valiall$data=array(NA,c(dim(valiall$data)[1],ets-sts+1))
        valiall$time=echam$time
      }
    }
    if (ncep_vali) {
      if (cyr>1947 && cyr<2010){
        ti=which(floor(valiall$time)==cyr)
        sts=ti[1]
        ets=ti[length(ti)]
        valiall$data=valiall$data[,sts:ets]
        valiall$time=valiall$time[sts:ets]
      } else {
        ti=which(floor(valiall$time)==cyr)
        sts=ti[1]
        ets=ti[length(ti)]
        valiall$data=array(NA,c(dim(valiall$data)[1],ets-sts+1))
        valiall$time=echam$time
      }   
    }
    if (recon_vali) {
      if (cyr>1750 && cyr<1901){
        ti=which(floor(valiall$time)==cyr)
        sts=ti[1]
        ets=ti[length(ti)]
        valiall$data=valiall$data[,sts:ets]
        valiall$time=valiall$time[sts:ets]
      } else {
        ti=which(floor(valiall$time)==cyr)
        sts=ti[1]
        ets=ti[length(ti)]
        valiall$data=array(NA,c(dim(valiall$data)[1],ets-sts+1))
        valiall$time=echam$time
      }  
    }
  }

  if (instrumental) {
    if ((ghcn_temp) & (dim(ghcn$data)[2]>0)) {
      if (anomaly_assim){
        ti=which(floor(ghcn$time)>=(cyr-35) &
                   floor(ghcn$time)<=(cyr+35))
        sts=ti[1]
        ets=ti[length(ti)]
        ghcn.tmp=ghcn
        ghcn.tmp$data=t(ghcn$data[sts:ets,])
        ghcn.tmp$time=ghcn$time[sts:ets]
        ghcn.anom <- ghcn
        ghcn.anom$data <- (t(ghcn$data) - matrix(rep(apply(array(ghcn.tmp$data,
           c(nrow(ghcn.tmp$data), 12, ncol(ghcn.tmp$data)/12)), 1:2, mean,na.rm=T),
           (length(ghcn$time)/12)), nrow=nrow(t(ghcn$data))))
        ghcn <- ghcn.anom
        ghcn$data <- t(ghcn.anom$data)
      }
      # cut 2-year time period as not done yet
      #ti=which(floor(ghcn$time)>=syr & floor(ghcn$time)<=eyr) # old version
      ti=which(floor(ghcn$time)>=(cyr-1) & floor(ghcn$time)<=cyr)
      sts=ti[1]
      ets=ti[length(ti)]
      ghcn$data=ghcn$data[sts:ets,]       
      ghcn$time=ghcn$time[sts:ets]
    }
    if (ghcn_prec) {
      if (dim(ghcn_precip$data)[2]>0) {  
        if (anomaly_assim){
          ti=which(floor(ghcn_precip$time)>=(cyr-35) &
                     floor(ghcn_precip$time)<=(cyr+35))
          sts=ti[1]
          ets=ti[length(ti)]
          ghcn_precip.tmp=ghcn_precip
          ghcn_precip.tmp$data=t(ghcn_precip$data[sts:ets,])
          ghcn_precip.tmp$time=ghcn_precip$time[sts:ets]
          ghcn_precip.anom <- ghcn_precip
          ghcn_precip.anom$data <- (t(ghcn_precip$data) - matrix(rep(apply(array(ghcn_precip.tmp$data,
            c(nrow(ghcn_precip.tmp$data), 12, ncol(ghcn_precip.tmp$data)/12)), 1:2, mean,na.rm=T),
            (length(ghcn_precip$time)/12)), nrow=nrow(t(ghcn_precip$data)))) 
          ghcn_precip <- ghcn_precip.anom
          ghcn_precip$data <- t(ghcn_precip.anom$data)
        }
        # cut time period as not done yet
        ti=which(floor(ghcn_precip$time)>=(cyr-1) & floor(ghcn_precip$time)<=cyr)
        sts=ti[1]
        ets=ti[length(ti)]
        ghcn_precip$data=ghcn_precip$data[sts:ets,]       
        ghcn_precip$time=ghcn_precip$time[sts:ets]
      }  
    }
    if ((yuri_slp) & (dim(inst_slp$data)[2]>0)) {  
      if (anomaly_assim){
        ti=which(floor(inst_slp$time)>=(cyr-35) &
                   floor(inst_slp$time)<=(cyr+35))
        sts=ti[1]
        ets=ti[length(ti)]
        inst_slp.tmp=inst_slp
        inst_slp.tmp$data=t(inst_slp$data[sts:ets,])
        inst_slp.tmp$time=inst_slp$time[sts:ets]
        inst_slp.anom <- inst_slp
        inst_slp.anom$data <- (t(inst_slp$data) - matrix(rep(apply(array(inst_slp.tmp$data,
          c(nrow(inst_slp.tmp$data), 12, ncol(inst_slp.tmp$data)/12)), 1:2, mean,na.rm=T),
          (length(inst_slp$time)/12)), nrow=nrow(t(inst_slp$data)))) 
        inst_slp <- inst_slp.anom
        inst_slp$data <- t(inst_slp.anom$data)
      }
      # cut time period as not done yet
      ti=which(floor(inst_slp$time)>=(cyr-1) & floor(inst_slp$time)<=cyr)
      sts=ti[1]
      ets=ti[length(ti)]
      inst_slp$data=inst_slp$data[sts:ets,]       
      inst_slp$time=inst_slp$time[sts:ets]
    }
    if ((yuri_temp) & (dim(inst_t$data)[2]>0)) {
      if (anomaly_assim){
        ti=which(floor(inst_t$time)>=(cyr-35) &
                   floor(inst_t$time)<=(cyr+35))
        sts=ti[1]
        ets=ti[length(ti)]
        inst_t.tmp=inst_t
        inst_t.tmp$data=t(inst_t$data[sts:ets,])
        inst_t.tmp$time=inst_t$time[sts:ets]
        inst_t.anom <- inst_t
        inst_t.anom$data <- (t(inst_t$data) - matrix(rep(apply(array(inst_t.tmp$data,
          c(nrow(inst_t.tmp$data), 12, ncol(inst_t.tmp$data)/12)), 1:2, mean,na.rm=T),
          (length(inst_t$time)/12)), nrow=nrow(t(inst_t$data))))
        inst_t <- inst_t.anom
        inst_t$data <- t(inst_t.anom$data)
      }
      # cut time period as not done yet
      ti=which(floor(inst_t$time)>=(cyr-1) & floor(inst_t$time)<=cyr)
      sts=ti[1]
      ets=ti[length(ti)]
      inst_t$data=inst_t$data[sts:ets,]       
      inst_t$time=inst_t$time[sts:ets]
    }
    if (ghcn_prec) {  
      if ((dim(inst_slp$data)[2]==0) & (dim(inst_t$data)[2]==0) & (dim(ghcn$data)[2]==0) & 
          (dim(ghcn_precip$data)[2]==0)) {instrumental=F}
    } else {
      if ((dim(inst_slp$data)[2]==0) & (dim(inst_t$data)[2]==0) & (dim(ghcn$data)[2]==0)) {
        instrumental=F}
    }  
  }


  if (docum) {
    if (anomaly_assim){
      ti=which(floor(doc_t_mon$time)>=(cyr-35) &
                 floor(doc_t_mon$time)<=(cyr+35))
      sts=ti[1]
      ets=ti[length(ti)]
      doc_t.tmp=doc_t_mon
      doc_t.tmp$data=t(doc_t_mon$data[sts:ets,])
      doc_t.tmp$time=doc_t_mon$time[sts:ets]
      doc_t_mon.anom <- doc_t_mon
      doc_t_mon.anom$data <- (t(doc_t_mon$data) - matrix(rep(apply(array(doc_t.tmp$data,
        c(nrow(doc_t.tmp$data), 12, ncol(doc_t.tmp$data)/12)), 1:2, mean,na.rm=T),
        (length(doc_t_mon$time)/12)), nrow=nrow(t(doc_t_mon$data))))
      doc_t_mon <- doc_t_mon.anom
      doc_t_mon$data <- t(doc_t_mon.anom$data)
      doc_t_mon$season <- seq(1,12)
      if (scaleprox){
        # calc echam variability at data location
        # short way (not integrated yet!):
        #    m=getgridboxnum(doc_t_mon,echam)
        # long way
        echamatprox.arr <- array(NA,c(length(doc_t_mon$lon), length(doc_t_mon$time)))
        for(i in 1:(length(doc_t_mon$lon))){
          plon <- doc_t_mon$lon[i]
          plat <- doc_t_mon$lat[i]
          clon <- echam$lon[!is.na(echam$lon)]
          clat <- echam$lat[!is.na(echam$lat)]
          k=which(abs(clon-plon+0.001)==min(abs(clon-plon+0.001))) 
          l=which(abs(clat-plat)==min(abs(clat-plat)))
          # search at surrounding grid boxes if no grid box matches (coastal problem)
          if (max(match(k,l,nomatch=-99999))==-99999) {
            k=which(abs(clon-plon+2.001)==min(abs(clon-plon+2.001)))
            l=which(abs(clat-plat)==min(abs(clat-plat)))
          }
          if (max(match(k,l,nomatch=-99999))==-99999) {
            k=which(abs(clon-plon-2.001)==min(abs(clon-plon-2.001)))
            l=which(abs(clat-plat)==min(abs(clat-plat)))
          }
          if (max(match(k,l,nomatch=-99999))==-99999) {
            k=which(abs(clon-plon+0.001)==min(abs(clon-plon+0.001)))
            l=which(abs(clat-plat+2)==min(abs(clat-plat+2)))
          }
          if (max(match(k,l,nomatch=-99999))==-99999) {
            k=which(abs(clon-plon+0.001)==min(abs(clon-plon+0.001)))
            l=which(abs(clat-plat-2)==min(abs(clat-plat-2)))
          }
          if (max(match(k,l,nomatch=-99999))==-99999) {
            k=which(abs(clon-plon+2.001)==min(abs(clon-plon+2.001)))
            l=which(abs(clat-plat+2)==min(abs(clat-plat+2)))
          }
          if (max(match(k,l,nomatch=-99999))==-99999) {
            k=which(abs(clon-plon-2.001)==min(abs(clon-plon-2.001)))
            l=which(abs(clat-plat-2)==min(abs(clat-plat-2)))
          }
          if (max(match(k,l,nomatch=-99999))==-99999) {
            k=which(abs(clon-plon-2.001)==min(abs(clon-plon-2.001)))
            l=which(abs(clat-plat+2)==min(abs(clat-plat+2)))
          }
          if (max(match(k,l,nomatch=-99999))==-99999) {
            k=which(abs(clon-plon+2.001)==min(abs(clon-plon+2.001)))
            l=which(abs(clat-plat-2)==min(abs(clat-plat-2)))
          }
          m=k[which(match(k,l)>0)]
          if (length(m) > 0) {
            if (doc_t_mon$names[i]=="precip") { 
              m <- m+length(echam$lon) 
            } else if (doc_t_mon$names[i]=="slp") { 
              m <- m+(2*length(echam$lon)) 
            }
            scalefac <- 1/echam.sd[m,]
          } else {
            scalefac <- NA
          }
        }
        doc_t_mon$data <- t(scale(t(doc_t_mon$data),center=F,scale=rep(scalefac,
                             (dim(doc_t_mon$data)[1]/12))))
      }
# next lines commented because NO seasonal documentary data is used at the moment      
#       ti=which(floor(doc_t_seas$time)>=(cyr-35) &
#                  floor(doc_t_seas$time)<=(cyr+35))
#       sts=ti[1]
#       ets=ti[length(ti)]
#       doc_t.tmp=doc_t_seas
#       doc_t.tmp$data=t(doc_t_seas$data[sts:ets,])
#       doc_t.tmp$time=doc_t_seas$time[sts:ets]
#       doc_t_seas.anom <- doc_t_seas
#       doc_t_seas.anom$data <- (t(doc_t_seas$data) - matrix(rep(apply(array(doc_t.tmp$data,
#         c(nrow(doc_t.tmp$data), 4, ncol(doc_t.tmp$data)/4)), 1:2, mean,na.rm=T),
#         (length(doc_t.tmp$time)/4)), nrow=nrow(t(doc_t_seas$data))))
#       doc_t_seas <- doc_t_seas.anom
#       doc_t_seas$data <- t(doc_t_seas.anom$data)
#       doc_t_seas$season <- c('win','spr','sum','aut')
#       if (scaleprox){
#         # calc echam variability at data location
#         echamatprox.arr <- array(NA,c(length(doc_t_seas$lon), length(doc_t_seas$time)))
#         for(i in 1:(length(doc_t_seas$lon))){
#           plon <- doc_t_seas$lon[i]
#           plat <- doc_t_seas$lat[i]
#           clon <- echam$lon[!is.na(echam$lon)]
#           clat <- echam$lat[!is.na(echam$lat)]
#           k=which(abs(clon-plon+0.001)==min(abs(clon-plon+0.001))) 
#           l=which(abs(clat-plat)==min(abs(clat-plat)))
#           if (max(match(k,l,nomatch=-99999))==-99999) {
#             k=which(abs(clon-plon+2.001)==min(abs(clon-plon+2.001)))
#             l=which(abs(clat-plat)==min(abs(clat-plat)))
#           }
#           if (max(match(k,l,nomatch=-99999))==-99999) {
#             k=which(abs(clon-plon-2.001)==min(abs(clon-plon-2.001)))
#             l=which(abs(clat-plat)==min(abs(clat-plat)))
#           }
#           if (max(match(k,l,nomatch=-99999))==-99999) {
#             k=which(abs(clon-plon+0.001)==min(abs(clon-plon+0.001)))
#             l=which(abs(clat-plat+2)==min(abs(clat-plat+2)))
#           }
#           if (max(match(k,l,nomatch=-99999))==-99999) {
#             k=which(abs(clon-plon+0.001)==min(abs(clon-plon+0.001)))
#             l=which(abs(clat-plat-2)==min(abs(clat-plat-2)))
#           }
#           if (max(match(k,l,nomatch=-99999))==-99999) {
#             k=which(abs(clon-plon+2.001)==min(abs(clon-plon+2.001)))
#             l=which(abs(clat-plat+2)==min(abs(clat-plat+2)))
#           }
#           if (max(match(k,l,nomatch=-99999))==-99999) {
#             k=which(abs(clon-plon-2.001)==min(abs(clon-plon-2.001)))
#             l=which(abs(clat-plat-2)==min(abs(clat-plat-2)))
#           }
#           if (max(match(k,l,nomatch=-99999))==-99999) {
#             k=which(abs(clon-plon-2.001)==min(abs(clon-plon-2.001)))
#             l=which(abs(clat-plat+2)==min(abs(clat-plat+2)))
#           }
#           if (max(match(k,l,nomatch=-99999))==-99999) {
#             k=which(abs(clon-plon+2.001)==min(abs(clon-plon+2.001)))
#             l=which(abs(clat-plat-2)==min(abs(clat-plat-2)))
#           }
#           m=k[which(match(k,l)>0)]
#           if (length(m) > 0) {
#             if (doc_t_seas$names[i]=="precip") { 
#               m <- m+length(echam$lon) 
#             } else if (doc_t_seas$names[i]=="slp") { 
#               m <- m+(2*length(echam$lon)) 
#             }
#             scalefac <- 1/c(mean(echam.sd[m,c(12,1,2)]),mean(echam.sd[m,c(3,4,5)]),
#                             mean(echam.sd[m,c(6,7,8)]),mean(echam.sd[m,c(9,10,11)]))
#           } else {
#             scalefac <- NA
#           }
#         }
#         doc_t_seas$data <- t(scale(t(doc_t_seas$data),center=F,scale=rep(scalefac,
#                                  (dim(doc_t_seas$data)[1]/4))))
#       }
# 
# # just scale if one value per year!      
#       ti=which(floor(doc_t_JFMA$time)>=(cyr-35) &
#                  floor(doc_t_JFMA$time)<=(cyr+35))
#       sts=ti[1]
#       ets=ti[length(ti)]
#       doc_t.mean=mean(doc_t_JFMA$data[sts:ets])
#       doc_t_JFMA$data <- doc_t_JFMA$data - doc_t.mean
#       doc_t_JFMA$season <- 'JFMA'
#       if (scaleprox){
#         # calc echam variability at data location
#         echamatprox.arr <- array(NA,c(length(doc_t_JFMA$lon), length(doc_t_JFMA$time)))
#         for(i in 1:(length(doc_t_JFMA$lon))){
#           plon <- doc_t_JFMA$lon[i]
#           plat <- doc_t_JFMA$lat[i]
#           clon <- echam$lon[!is.na(echam$lon)]
#           clat <- echam$lat[!is.na(echam$lat)]
#           k=which(abs(clon-plon+0.001)==min(abs(clon-plon+0.001))) 
#           l=which(abs(clat-plat)==min(abs(clat-plat)))
#           if (max(match(k,l,nomatch=-99999))==-99999) {
#             k=which(abs(clon-plon+2.001)==min(abs(clon-plon+2.001)))
#             l=which(abs(clat-plat)==min(abs(clat-plat)))
#           }
#           if (max(match(k,l,nomatch=-99999))==-99999) {
#             k=which(abs(clon-plon-2.001)==min(abs(clon-plon-2.001)))
#             l=which(abs(clat-plat)==min(abs(clat-plat)))
#           }
#           if (max(match(k,l,nomatch=-99999))==-99999) {
#             k=which(abs(clon-plon+0.001)==min(abs(clon-plon+0.001)))
#             l=which(abs(clat-plat+2)==min(abs(clat-plat+2)))
#           }
#           if (max(match(k,l,nomatch=-99999))==-99999) {
#             k=which(abs(clon-plon+0.001)==min(abs(clon-plon+0.001)))
#             l=which(abs(clat-plat-2)==min(abs(clat-plat-2)))
#           }
#           if (max(match(k,l,nomatch=-99999))==-99999) {
#             k=which(abs(clon-plon+2.001)==min(abs(clon-plon+2.001)))
#             l=which(abs(clat-plat+2)==min(abs(clat-plat+2)))
#           }
#           if (max(match(k,l,nomatch=-99999))==-99999) {
#             k=which(abs(clon-plon-2.001)==min(abs(clon-plon-2.001)))
#             l=which(abs(clat-plat-2)==min(abs(clat-plat-2)))
#           }
#           if (max(match(k,l,nomatch=-99999))==-99999) {
#             k=which(abs(clon-plon-2.001)==min(abs(clon-plon-2.001)))
#             l=which(abs(clat-plat+2)==min(abs(clat-plat+2)))
#           }
#           if (max(match(k,l,nomatch=-99999))==-99999) {
#             k=which(abs(clon-plon+2.001)==min(abs(clon-plon+2.001)))
#             l=which(abs(clat-plat-2)==min(abs(clat-plat-2)))
#           }
#           m=k[which(match(k,l)>0)]
#           if (length(m) > 0) {
#             if (doc_t_JFMA$names[i]=="precip") { 
#               m <- m+length(echam$lon) 
#             } else if (doc_t_JFMA$names[i]=="slp") { 
#               m <- m+(2*length(echam$lon)) 
#             }
#             scalefac <- 1/mean(echam.sd[m,1:4])
#           } else {
#             scalefac <- NA
#           }
#         }
#         if (!is.na(scalefac)) {
#           doc_t_JFMA$data <- as.vector(scale(doc_t_JFMA$data,center=F,scale=scalefac))
#         }
#       }
#     } 
# 
#     ti=which(floor(doc_t_AMJJA$time)>=(cyr-35) &
#            floor(doc_t_AMJJA$time)<=(cyr+35))
#     sts=ti[1]
#     ets=ti[length(ti)]
#     doc_t.mean=mean(doc_t_AMJJA$data[sts:ets])
#     doc_t_AMJJA$data <- doc_t_AMJJA$data - doc_t.mean
#     doc_t_AMJJA$season <- 'AMJJA'
#     if (scaleprox){
#       # calc echam variability at data location
#       echamatprox.arr <- array(NA,c(length(doc_t_AMJJA$lon), length(doc_t_AMJJA$time)))
#       for(i in 1:(length(doc_t_AMJJA$lon))){
#         plon <- doc_t_AMJJA$lon[i]
#         plat <- doc_t_AMJJA$lat[i]
#         clon <- echam$lon[!is.na(echam$lon)]
#         clat <- echam$lat[!is.na(echam$lat)]
#         k=which(abs(clon-plon+0.001)==min(abs(clon-plon+0.001))) 
#         l=which(abs(clat-plat)==min(abs(clat-plat)))
#         if (max(match(k,l,nomatch=-99999))==-99999) {
#           k=which(abs(clon-plon+2.001)==min(abs(clon-plon+2.001)))
#           l=which(abs(clat-plat)==min(abs(clat-plat)))
#         }
#         if (max(match(k,l,nomatch=-99999))==-99999) {
#           k=which(abs(clon-plon-2.001)==min(abs(clon-plon-2.001)))
#           l=which(abs(clat-plat)==min(abs(clat-plat)))
#         }
#         if (max(match(k,l,nomatch=-99999))==-99999) {
#           k=which(abs(clon-plon+0.001)==min(abs(clon-plon+0.001)))
#           l=which(abs(clat-plat+2)==min(abs(clat-plat+2)))
#         }
#         if (max(match(k,l,nomatch=-99999))==-99999) {
#           k=which(abs(clon-plon+0.001)==min(abs(clon-plon+0.001)))
#           l=which(abs(clat-plat-2)==min(abs(clat-plat-2)))
#         }
#         if (max(match(k,l,nomatch=-99999))==-99999) {
#           k=which(abs(clon-plon+2.001)==min(abs(clon-plon+2.001)))
#           l=which(abs(clat-plat+2)==min(abs(clat-plat+2)))
#         }
#         if (max(match(k,l,nomatch=-99999))==-99999) {
#           k=which(abs(clon-plon-2.001)==min(abs(clon-plon-2.001)))
#           l=which(abs(clat-plat-2)==min(abs(clat-plat-2)))
#         }
#         if (max(match(k,l,nomatch=-99999))==-99999) {
#           k=which(abs(clon-plon-2.001)==min(abs(clon-plon-2.001)))
#           l=which(abs(clat-plat+2)==min(abs(clat-plat+2)))
#         }
#         if (max(match(k,l,nomatch=-99999))==-99999) {
#           k=which(abs(clon-plon+2.001)==min(abs(clon-plon+2.001)))
#           l=which(abs(clat-plat-2)==min(abs(clat-plat-2)))
#         }
#         m=k[which(match(k,l)>0)]
#         if (length(m) > 0) {
#           if (doc_t_AMJJA$names[i]=="precip") { 
#             m <- m+length(echam$lon) 
#           } else if (doc_t_AMJJA$names[i]=="slp") { 
#             m <- m+(2*length(echam$lon)) 
#           }
#           scalefac <- 1/mean(echam.sd[m,4:8])
#         } else {
#           scalefac <- NA
#         }
#       }
#       if (!is.na(scalefac)) {
#         doc_t_AMJJA$data <- as.vector(scale(doc_t_AMJJA$data,center=F,scale=scalefac))
#       }
    }
  }





  if (real_proxies) {
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
  }  


#########################################################################################
# prepare proxy/instr. assimilation data
#########################################################################################
  
# real trw proxy multiple regression approach
  if ((real_proxies) & (!instrumental)) {  
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
#    proxies.allts<-list(data=realprox.allts$data, 
#                    lon=realprox.allts$lon, 
#                    lat=realprox.allts$lat, 
#                    names=realprox.allts$names, 
#                    time=realprox.allts$time)  
  } 
  
  if (docum) {
    if (sixmonstatevector) { 
      doc_t_mon$data <- t(doc_t_mon$data)
      tmp1 <- array(doc_t_mon$data,c(dim(doc_t_mon$data)[1] *  dim(doc_t_mon$data)[2]))
      # cut oct syr to sep eyr
      tmp2 <- tmp1[((9*dim(doc_t_mon$data)[1]+1):(length(tmp1)[1] - 
                                               (3*dim(doc_t_mon$data)[1])))] 
      doc_t_mon$data <- array(tmp2,c(dim(tmp2)[1]/(((dim(doc_t_mon$data)[2]
                                                /12)-1)*2),(((dim(doc_t_mon$data)[2]/12)-1)*2))) 
      # reconvert to 2 seasons per year
      # first all doc_t_mon for one month, than all prox for next month
      doc_t_mon$time <- seq((floor(doc_t_mon$time[1])+1),
                          floor(doc_t_mon$time[length(doc_t_mon$time)])+0.5,0.5) 
      doc_t_mon$names <- rep(doc_t_mon$names,6)
      doc_t_mon$lon <- rep(doc_t_mon$lon,6)
      doc_t_mon$lat <- rep(doc_t_mon$lat,6)
    }
    ti=which(floor(doc_t_mon$time)==cyr)
    sts=ti[1]
    ets=ti[length(ti)]
    doc_t_mon.allts=doc_t_mon
    doc_t_mon$data=doc_t_mon$data[,sts:ets]
    doc_t_mon$time=doc_t_mon$time[sts:ets]
    doc_t_mon$des=rep('mon',nrow(doc_t_mon$data))
  
#    # seasonal documentary data: attribiute JFM to winter season and later create fitting H operator
#     if (sixmonstatevector) { #allows for merging differnt temp. resolutions
#       doc_t_seas$data <- t(doc_t_seas$data)
#       tmp1 <- doc_t_seas$data[,((doc_t_seas$time-trunc(doc_t_seas$time)==0.125) | 
#                                   (doc_t_seas$time-trunc(doc_t_seas$time)==0.625))]
#       tmp2 <- array(tmp1,c(dim(tmp1)[1] *  dim(tmp1)[2]))
#       tmp3 <- tmp2[(2*dim(tmp1)[1]+1):(length(tmp2)[1])] # cut oct syr to sep eyr
#       doc_t_seas$data <- array(tmp3,c(dim(tmp3)[1]/(((dim(tmp1)[2]/2)-1)*2),(((dim(tmp1)[2]/2)-1)*2)))
#       doc_t_seas$time <- seq(syr+1,eyr+0.5,0.5) 
# #      doc_t_seas$names <- rep(doc_t_seas$names,2)
# #      doc_t_seas$lon <- rep(doc_t_seas$lon,2)
# #      doc_t_seas$lat <- rep(doc_t_seas$lat,2)
#     }
#     ti=which(floor(doc_t_seas$time)==cyr)
#     sts=ti[1]
#     ets=ti[length(ti)]
#     doc_t_seas.allts=doc_t_seas
#     doc_t_seas$data=doc_t_seas$data[,sts:ets]
#     doc_t_seas$time=doc_t_seas$time[sts:ets]
#     doc_t_seas$des=rep('seas',nrow(doc_t_seas$data))
#     
#     # attribute JFMA documentary data (1 value per year) to winter season and later create fitting H operator
#     ti=which(floor(doc_t_JFMA$time)==cyr)
#     doc_t_JFMA.allts=doc_t_JFMA
#     if (sixmonstatevector) {
#       doc_t_JFMA$data=t(matrix(c(doc_t_JFMA$data[ti],NA)))
#       doc_t_JFMA.allts$data=t(matrix(as.vector(rbind(doc_t_JFMA.allts$data[2:length(doc_t_JFMA.allts$data)],rep(NA,(length(doc_t_JFMA.allts$data)-1))))))
#     } else {
#       doc_t_JFMA$data=doc_t_JFMA$data[ti]
#       print("ACHTUNG: documentary proxies only included with option sixmonstatevector so far!!!")
#     }
#     doc_t_JFMA$time=doc_t_JFMA$time[ti]-0.5
#     doc_t_JFMA$des=rep('JFMA',nrow(doc_t_JFMA$data))
#     
#     # attribute AMJJA documentary data (1 value per year) to summer season and later create fitting H operator
#     ti=which(floor(doc_t_AMJJA$time)==cyr)
#     doc_t_AMJJA.allts=doc_t_AMJJA
#     if (sixmonstatevector) {
#       doc_t_AMJJA$data=t(matrix(c(NA,doc_t_AMJJA$data[ti])))
#       doc_t_AMJJA.allts$data=t(matrix(as.vector(rbind(doc_t_AMJJA.allts$data[2:length(doc_t_AMJJA.allts$data)],rep(NA,(length(doc_t_AMJJA.allts$data)-1))))))
#     } else {
#       doc_t_JFMA$data=doc_t_AMJJA$data[ti]
#       print("ACHTUNG: documentary proxies only included with option sixmonstatevector so far!!!")
#     }
#     doc_t_AMJJA$time=doc_t_AMJJA$time[ti]
#     doc_t_AMJJA$des=rep('AMJJA',nrow(doc_t_AMJJA$data))

    if ((!instrumental) & (!real_proxies)) {
      proxies <- doc_t_mon
#       proxies <- list(lon=c(doc_t_mon$lon,doc_t_seas$lon,doc_t_JFMA$lon,doc_t_AMJJA$lon), 
#                         lat=c(doc_t_mon$lat,doc_t_seas$lat,doc_t_JFMA$lat,doc_t_AMJJA$lat),
#                         data=rbind(doc_t_mon$data,doc_t_seas$data,doc_t_JFMA$data,doc_t_AMJJA$data), 
#                         names=c(doc_t_mon$names,doc_t_seas$names,doc_t_JFMA$names,doc_t_AMJJA$names),
#                         des=c(doc_t_mon$des,doc_t_seas$des,doc_t_JFMA$des,doc_t_AMJJA$des),
#                         time=doc_t_mon$time)
#      proxies.allts <- doc_t_mon.allts
#       proxies.allts <- list(lon=c(doc_t_mon.allts$lon,doc_t_seas.allts$lon,doc_t_JFMA.allts$lon,
#                          doc_t_AMJJA.allts$lon),lat=c(doc_t_mon.allts$lat,doc_t_seas.allts$lat,
#                          doc_t_JFMA.allts$lat,doc_t_AMJJA.allts$lat),
#                          data=rbind(doc_t_mon.allts$data,doc_t_seas.allts$data,
#                          doc_t_JFMA.allts$data,doc_t_AMJJA.allts$data), 
#                          names=c(doc_t_mon.allts$names,doc_t_seas.allts$names,
#                          doc_t_JFMA.allts$names,doc_t_AMJJA.allts$names),
#                          des=c(doc_t_mon.allts$des,doc_t_seas.allts$des,
#                          doc_t_JFMA.allts$des,doc_t_AMJJA.allts$des),
#                          time=doc_t_mon.allts$time)        
    } else {
      docall <- doc_t_mon
#       docall <- list(lon=c(doc_t_mon$lon,doc_t_seas$lon,doc_t_JFMA$lon,doc_t_AMJJA$lon), 
#                       lat=c(doc_t_mon$lat,doc_t_seas$lat,doc_t_JFMA$lat,doc_t_AMJJA$lat),
#                       data=rbind(doc_t_mon$data,doc_t_seas$data,doc_t_JFMA$data,doc_t_AMJJA$data), 
#                       names=c(doc_t_mon$names,doc_t_seas$names,doc_t_JFMA$names,doc_t_AMJJA$names),
#                       time=doc_t_mon$time) 
      docall.allts <- doc_t_mon.allts
    }
  }

  
  
  
    
# treat multiple data in the same grid box:
#   at the moment only in case of instrumental data
#   for proxy or docum. data assim. without instr., all series get serially assimilated
  if (instrumental) {
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

    # mask proxy data if there is instrumental data in same grid box at that time
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
#          print(paste('set proxy',i,'to NA'))
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
#              print(paste('set doc',i,'to NA'))
                docall$data[i,] = NA
                docall$lon[i] <- NA
                docall$lat[i] <- NA
              }
            }
#           pos <- apply(!is.na(docall$data),1,any)
#           docall$data = docall$data[pos,]
#           docall$lon <- docall$lon[pos]
#           docall$lat <- docall$lat[pos]
#           docall$height <- docall$height[pos]
#           docall$season <- docall$season[pos]
#           docall$des <- docall$des[pos]
#           docall$names <- docall$names[pos]
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
    
    # add data source/type information to list
    if (instrumental) {inst$sour <- rep('inst',length(inst$lon))}
    if (docum) {
      docall$sour <- rep('doc',length(docall$lon))
      docall.allts$sour <- rep('doc',length(docall$lon))    
    }
    if (real_proxies) {
      realprox$sour <- rep('prox',length(realprox$lon))
      realprox.allts$sour <- rep('prox',length(realprox$lon))
    }

    if (sixmonstatevector) { 
      # change proxy array that only 2 ts instead of 12 monthly ts but 6 times as 
      # many variables, one for each month
      tmp1 <- array(inst$data,c(dim(inst$data)[1] *  dim(inst$data)[2]))
      tmp2 <- tmp1[((9*dim(inst$data)[1]+1):(length(tmp1)[1] - 
                    (3*dim(inst$data)[1])))] # cut oct syr to sep eyr
      inst$data <- array(tmp2,c(dim(tmp2)[1]/(((dim(inst$data)[2]
                                 /12)-1)*2),(((dim(inst$data)[2]/12)-1)*2))) 
      # reconvert to 2 seasons per year
      # first all inst for one month, than all prox for next month
      #inst$time <- seq(syr+1,eyr+0.5,0.5) # changes 2017-07-26
      inst$time <- c(cyr,cyr+0.5) 
      inst$names <- rep(inst$names,6)
      inst$sour <- rep(inst$sour,6)
      inst$lon <- rep(inst$lon,6)
      inst$lat <- rep(inst$lat,6)
      if (avg_prox_per_grid) {
        inst$numavg <- rep(inst$numavg,6)
      }
    }
    
    # ti=which(floor(inst$time)==cyr) # commented 2017-07-26
    # sts=ti[1]
    # ets=ti[length(ti)]
    # inst.allts=inst
    # inst$data=inst$data[,sts:ets]
    # inst$time=inst$time[sts:ets]
    
    # combine assimilation data into variable named "proxies")
    if (!real_proxies) {        
      proxies <- inst
#      proxies.allts <- inst.allts
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
#        proxies.allts<-list(data=rbind(inst.allts$data,realprox.allts$data), 
#                          lon=c(inst.allts$lon,realprox.allts$lon), 
#                          lat=c(inst.allts$lat,realprox.allts$lat), 
#                          names=c(inst.allts$names,realprox.allts$names), 
#                          sour=c(inst.allts$sour,realprox.allts$sour), 
#                          time=inst.allts$time)  
                          # WHAT TO DO WITH $MR AND $VAR_RESIDU and INST.NUMAVG???      
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
#      proxies.allts<-list(data=rbind(proxies.allts$data,docall.allts$data), 
#                    lon=c(proxies.allts$lon,docall.allts$lon), 
#                    lat=c(proxies.allts$lat,docall.allts$lat), 
#                    names=c(proxies.allts$names,docall.allts$names), 
#                    sour=c(proxies.allts$sour,docall.allts$sour), 
#                    time=proxies.allts$time)  
    }
  } # end "if (instrumental)"

  if (!instrumental & real_proxies & docum) {
    docall$sour <- rep('doc',length(docall$lon))
    docall.allts$sour <- rep('doc',length(docall$lon))    
    realprox$sour <- rep('prox',length(realprox$lon))
    realprox.allts$sour <- rep('prox',length(realprox$lon)) 
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
#    proxies.allts<-list(data=rbind(realprox.allts$data,docall.allts$data), 
#                        lon=c(realprox.allts$lon,docall.allts$lon), 
#                        lat=c(realprox.allts$lat,docall.allts$lat), 
#                        names=c(realprox.allts$names,docall.allts$names), 
#                        sour=c(realprox.allts$sour,docall.allts$sour), 
#                        time=realprox.allts$time)  
  }

  if (real_proxies) {
    if (dim(realprox$data)[1]==0) { real_proxies=F }
  }
  
  calibrate <- proxies
#  calibrate.allts <- proxies.allts
  if (sum(c(ncep_vali,cru_vali,recon_vali))>1) {
    print("WARNING: more than 1 validation data set selected!")
  }
  if (vali) {
    validate=valiall
    validate$ensmean=validate$data
  }
  
  print(paste('number of proxies/observations:',dim(calibrate$data)[1]))
  print("calc time preparing proxies")
  print(proc.time() - ptm1)




## convert to new data format
# ACHTUNG only 11 vars withOUT stream function included so far
  if (sixmonstatevector) {
    numvar <- length(unique(echam$names)) 
    echam$lon <-  c(rep(echam$lon, (numvar*6)))
    echam$lat <- c(rep(echam$lat, (numvar*6)))
    echam$names <- c(rep(echam$names, 6))
    echam <- echam[c('data', 'ensmean', 'lon', 'lat', 'height', 'lsm.i', 'time', 'names')]
  }
  
  ## set up the cutoff distance for localisation
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

  
# ------------------------------------------------------------------------------
# Compute H (forward operator)
# ------------------------------------------------------------------------------
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

  
  
  

# ------------------------------------------------------------------------------
# run analysis
# ------------------------------------------------------------------------------
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
        if ((!is.na(dist[1])) & (!is.na(wgt[1]))){
          for (i in 1:ntim){
            if (!is.na(calibrate$data[j,i])) {
            #            print(i)
#              ptm10=proc.time()
              x2tmp <- analysis$data[,i,] # entire state vector at time step i, all members
              x2 <- x2tmp[h.i,,drop=F] # state vector at time step i and h.i, all members
              PH <- (analysis$data[,i,] %*% t(x2) / (nens - 1) * wgt) %*% t(H)
#              print(proc.time() - ptm10)
              if (covarclim>0 & covarclim<100) { # (covarclim>50) {
                x2climtmp <- echanomallts$data[,i,] 
                x2clim <- x2climtmp[h.i,,drop=F] 
                PHclim <- (echanomallts$data[,i,] %*% t(x2clim) / 
                             ((dim(echanomallts$data)[3]) - 1) * wgt) %*% t(H)
                #PH <- (PH + PHclim) / 2
                PH <- (PH*(1-(covarclim/100))) + (PHclim*(covarclim/100))
              } else if (covarclim==100) {
                x2climtmp <- echanomallts$data[,i,] 
                x2clim <- x2climtmp[h.i,,drop=F] 
                PH <- (echanomallts$data[,i,] %*% t(x2clim) / 
                         ((dim(echanomallts$data)[3]) - 1) * wgt) %*% t(H) 
 #               print(proc.time() - ptm10)
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
      }    
      ## add ensemble mean analysis back in
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
                    #            print(i)
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
            # add echam.sd[lloc],] and analysis SD (not yet calc.) to output!
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
              colnames(looensmeanres) <- c("lon","lat",
                                           "diff apr","diff may","diff jun",
                                           "diff jul","diff aug","diff sep",
                                           "diff oct","diff nov","diff dec",
                                           "diff jan","diff feb","diff mar",
                                           "cali apr","cali may","cali jun",
                                           "cali jul","cali aug","cali sep",
                                           "cali oct","cali nov","cali dec",
                                           "cali jan","cali feb","cali mar",
                                           "ech apr","ech may","ech jun",
                                           "ech jul","ech aug","ech sep",
                                           "ech oct","ech nov","ech dec",
                                           "ech jan","ech feb","ech mar",
                                           "ana apr","ana may","ana jun",
                                           "ana jul","ana aug","ana sep",
                                           "ana oct","ana nov","ana dec",
                                           "ana jan","ana feb","ana mar",
                                           "loo apr","loo may","loo jun",
                                           "loo jul","loo aug","loo sep",
                                           "loo oct","loo nov","loo dec",
                                           "loo jan","loo feb","loo mar")
              colnames(loodatares) <- c("lon","lat",paste(rep(seq(1,30),each=12),rep(c(
                                        "diff apr","diff may","diff jun",
                                        "diff jul","diff aug","diff sep",
                                        "diff oct","diff nov","diff dec",
                                        "diff jan","diff feb","diff mar"),30)),
                                        paste(rep(seq(1,30),each=12),rep(c(
                                          "ech apr","ech may","ech jun",
                                          "ech jul","ech aug","ech sep",
                                          "ech oct","ech nov","ech dec",
                                          "ech jan","ech feb","ech mar"),30)),
                                        paste(rep(seq(1,30),each=12),rep(c(
                                          "ana apr","ana may","ana jun",
                                          "ana jul","ana aug","ana sep",
                                          "ana oct","ana nov","ana dec",
                                          "ana jan","ana feb","ana mar"),30)),
                                        paste(rep(seq(1,30),each=12),rep(c(
                                          "loo apr","loo may","loo jun",
                                          "loo jul","loo aug","loo sep",
                                          "loo oct","loo nov","loo dec",
                                          "loo jan","loo feb","loo mar"),30)))
#            }
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
