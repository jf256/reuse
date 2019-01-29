##########################################################################################
# generate Rdata files from orig sources
########################################################################################## 

if (generate_ECHAM){
  print("generate_ECHAM")
  # create ECHAM .RData from .nc
  # ACHTUNG: members 001 to 030 with land surface bug, 103 and 2xx without
  # for 60 ensm: timlim=c(1941,1970)
  read_echam4('EnSRF.ccc400_0', timlim=c(1601,2005), small=every2grid, landonly=land_only)
} 

if (generate_ECHAM_anom){
  # read echam 71yr anom, clim and sd calculated with cdo from .nc files to .RData
  # for 60 ensm: timlim=c(1941,1970)
  print("generate_ECHAM_anom") 
  read_echam4('ano', path=echanompath, timlim=c(1601,2005), small=every2grid, 
              landonly=land_only, anom=T)
  read_echam4('EnSRF', path=echclimpath, timlim=c(1635,1970), small=every2grid, 
              landonly=land_only, clim=T)
  read_echam4('EnSRF', path=echsdpath, timlim=c(1601,2005), small=every2grid, 
              landonly=land_only, std=T)
}

# if (generate_ECHAM_1901_70){
#   print("generate_ECHAM_1901_70")
#   # ECHAM data for bias calculation with real proxy data
#   echam1901_70 <- read_echam_ensmean('EnSRF', timlim=c(1901,1970),small=F)
#   save(echam1901_70, file="../data/echam_1911-70.Rdata")
# } 

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

if (generate_ind_recon){
  # read Broennimann et al. 2009 atm. indices from .txt to .Rdata for comparison
  #if (syr<1901){syr_ind=1901} else {syr_ind=syr}
  #if (eyr>2004){eyr_ind=2004} else {eyr_ind=eyr}
  syr_ind=1900
  eyr_ind=2000
  dataextdir='/mnt/climstor/giub/EKF400/'
  indicespath <- paste0(dataextdir,'vali_data/indices/')
  # ind=read.table(file=paste0(indicespath,'stefan_monthly_indices.txt'
  #                           ,sep=''),header=T)
  # ind_rec_dimi = window(ts(c(rep(NA,length(ind[,colnames(ind) == 'Z100'])),rep(NA,132)),
  #                          start=ind[1,colnames(ind) == 'yr'],freq=12),syr_ind,c(eyr_ind,12))
  # ind_rec_z100 = window(ts(ind[,colnames(ind)=='Z100'],frequency = 12, start = ind[1,colnames(ind)=='yr']), start=syr_ind, end =c(eyr_ind, 12))
  # 
  # ind_rec_z300 = window(ts(ind[,colnames(ind) == 'Z300'],start=ind[1,
  #                          colnames(ind) == 'yr'],freq=12),syr_ind,c(eyr_ind,12))
  # ind_rec_pwc = window(ts(ind[,colnames(ind) == 'PWC'],start=ind[1,
  #                          colnames(ind) == 'yr'],freq=12),syr_ind,c(eyr_ind,12))
  # ind_rec_hc = window(ts(ind[,colnames(ind) == 'HCL'],start=ind[1,
  #                          colnames(ind) == 'yr'],freq=12),syr_ind,c(eyr_ind,12))
  # ind_rec_sj = window(ts(ind[,colnames(ind) == 'SJ'],start=ind[1,
  #                          colnames(ind) == 'yr'],freq=12),syr_ind,c(eyr_ind,12))
  # indall = t(cbind(ind_rec_dimi, ind_rec_z100, ind_rec_z300, ind_rec_pwc, ind_rec_hc,
  #                 ind_rec_sj))
  # save(indall, file=paste0(dataintdir,'/indices/indices_recon_',syr_ind,'-',eyr_ind,'_monthly.Rdata',sep=''))
  
  #generate seasonal indices
  
  # syr_ind=1900
  # eyr_ind=2000
  # 
  # ind=read.table(file=paste0(indicespath,'stefan_seasonal_indices.txt'
  #                            ,sep=''),header=T)
  # ind_rec_dimi = window(ts(c(rep(NA,length(ind[,colnames(ind) == 'Z100'])),rep(NA,132)),
  #                          start=ind[1,colnames(ind) == 'yr'],freq=12),syr_ind,c(eyr_ind,12))
  # ind_rec_z100 = window(ts(ind[,colnames(ind)=='Z100'],frequency = 12, start = ind[1,colnames(ind)=='yr']), start=syr_ind, end =c(eyr_ind, 12))
  # 
  # ind_rec_z300 = window(ts(ind[,colnames(ind) == 'Z300'],start=ind[1,
  #                                                                  colnames(ind) == 'yr'],freq=12),syr_ind,c(eyr_ind,12))
  # ind_rec_pwc = window(ts(ind[,colnames(ind) == 'PWC'],start=ind[1,
  #                                                                colnames(ind) == 'yr'],freq=12),syr_ind,c(eyr_ind,12))
  # ind_rec_hc = window(ts(ind[,colnames(ind) == 'HCL'],start=ind[1,
  #                                                               colnames(ind) == 'yr'],freq=12),syr_ind,c(eyr_ind,12))
  # ind_rec_sj = window(ts(ind[,colnames(ind) == 'SJ'],start=ind[1,
  #                                                              colnames(ind) == 'yr'],freq=12),syr_ind,c(eyr_ind,12))
  # indall = t(cbind(ind_rec_dimi, ind_rec_z100, ind_rec_z300, ind_rec_pwc, ind_rec_hc,
  #                  ind_rec_sj))
  # save(indall, file=paste0(dataintdir,'/indices/indices_recon_',syr_ind,'-',eyr_ind,'_seasonal.Rdata',sep=''))
  # 
  # 
  
  
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
  source(paste0(dataextdir,"assim_data/data_yuri/t_assimil/read_all.R"))
}

if (generate_slp_yuri){
  print("generate_slp_yuri")
  source(paste0(dataextdir,"assim_data/data_yuri/slp_assimil/read_all.R"))
}

if (generate_DOCUM){
  print("generate_DOCUM")
  source(paste0(dataextdir,"assim_data/data_yuri/t_docu/read_seas.R"))
  source(paste0(dataextdir,"assim_data/data_yuri/t_docu/read_monthly.R"))
  source(paste0(dataextdir,"assim_data/data_yuri/t_docu/read_JFMA.R"))
  source(paste0(dataextdir,"assim_data/data_yuri/t_docu/read_AMJJA.R"))
}

# generate PAGES v2 from Raphi's 2018/01 export (version 1.6.1 from Raphi's data base):
#   /tank/exports/giub/EKF400/assimil_data/proxies/PAGES/read_pages_2018.R # Skript from Veronika

# generate NTREND
#   /tank/exports/giub/EKF400/assimil_data/proxies/NTREND/N-Trend.R        # Skript from Veronika

# generate petra, schweingr, mxd, etc. missing

if (generate_PROXIES){
  print("generate_PROXIES")
  read.these <- c("trw","mxd","schweingr","pages","ntrend","trw_petra")[c(TRW,MXD,SCHWEINGR,PAGES,NTREND,TRW_PETRA)]
  if(exists("realprox")){rm(realprox)}
  for (varname in read.these){
    if (varname=="trw") {
      print("reading Petra's 35 best TRW records used in first EKF400 version")
      trwprox <- read_proxy2(fsyr,feyr)
      realprox<-trwprox
      # NEW VERSION VERONIKA: PLEASE CHECK IF YOUR VERSION IS CORRECT
      ############# # because the dim is not equal, the year for trwprox stops in 1970, other two in 2005
      realprox$data <- rbind(trwprox$data,matrix(data=NA, nrow=length(seq(1971,2004)), ncol=dim(trwprox$data)[2])) # 34 not very elegant
      #############
      realprox$time <- c(realprox$time,seq(1971,2004))
    }
    
    if (varname=="mxd") {
      print("reading mxd")
      mxdprox <- read_proxy_mxd(fsyr,feyr)
      
      if (exists("realprox")){
        realprox$data <- cbind(realprox$data, mxdprox$data)
        realprox$lon <- c(realprox$lon, mxdprox$lon)
        realprox$lat <- c(realprox$lat, mxdprox$lat)
        realprox$mr <- rbind(realprox$mr, mxdprox$mr)
        realprox$var_residu <- c(realprox$var_residu, mxdprox$var_residu)
        
      } else { realprox<-mxdprox}
    }
    
    if (varname=="schweingr") {
      
      print("reading schweingr")
      schprox <- read_proxy_schweingr(fsyr,feyr)
      
      if (exists("realprox")){
        
        realprox$data <- cbind(realprox$data, schprox$data)
        realprox$lon <- c(realprox$lon, schprox$lon)
        realprox$lat <- c(realprox$lat, schprox$lat)
        realprox$mr <- rbind(realprox$mr, schprox$mr)
        realprox$var_residu <- c(realprox$var_residu, schprox$var_residu)
        
      } else { realprox<-schprox}
    }
    
    if (varname=="pages") {
      print("reading pages")
      pagesprox <- setup_read_pages(type) 
      
      if (exists("realprox")){
        
        realprox$data <- cbind(realprox$data, pagesprox$data)
        realprox$lon <- c(realprox$lon, pagesprox$lon)
        realprox$lat <- c(realprox$lat, pagesprox$lat)
        realprox$mr <- rbind(realprox$mr, pagesprox$mr)
        realprox$var_residu <- c(realprox$var_residu, pagesprox$var_residu)
        
      } else { realprox<-pagesprox}
     
    }
    
    if (varname=="ntrend") {
      print("reading ntrend")
      ntrend = read_ntrend(fsyr,feyr, validate=lm_fit_data)
      
      if (exists("realprox")){
        
        realprox$data <- cbind(realprox$data, ntrend$data)
        realprox$lon <- c(realprox$lon, ntrend$lon)
        realprox$lat <- c(realprox$lat, ntrend$lat)
        realprox$mr <- rbind(realprox$mr, ntrend$mr)
        realprox$var_residu <- c(realprox$var_residu, ntrend$var_residu)
        
      } else { realprox<-ntrend}
    }
    if (varname=="trw_petra") {
      print("reading trw_petra")
      trw_petra <- read_trw_petra(fsyr,feyr, validate=lm_fit_data) 
      
      if (exists("realprox")){
        
        realprox$data <- cbind(realprox$data, trw_petra$data)
        realprox$lon <- c(realprox$lon, trw_petra$lon)
        realprox$lat <- c(realprox$lat, trw_petra$lat)
        realprox$mr <- rbind(realprox$mr, trw_petra$mr)
        realprox$var_residu <- c(realprox$var_residu, trw_petra$var_residu)
        
      } else { realprox<-trw_petra}
    }
  }
  save(realprox, file=paste0("../data/proxies/real_proxies_",fsyr,"-",feyr,".Rdata"))
}

if (pseudo_prox) {
  if (generate_PSEUDO){
    print("generate_PSEUDO")
    pseudoprox<-read_pseudo()
    realprox<-pseudoprox
    save(realprox, file=paste0("../data/proxies/DAPS_pseudoproxies_",fsyr,"-",feyr,".Rdata"))
  } else {
    # load DAPS pseudo proxies in object called 'realprox'
    load(file=paste0("../data/proxies/DAPS_pseudoproxies_",fsyr,"-",feyr,".Rdata"))
  }  
}

