##########################################################################################
# generate Rdata files from orig sources
########################################################################################## 

if (generate_ECHAM){
  print("generate_ECHAM")
  print("ATTENTION: manually copy generated .Rdata file to climstor!")
  # create ECHAM .RData from .nc
  # ACHTUNG: members 001 to 030 with land surface bug, 103 and 2xx without
  # for 60 ensm: timlim=c(1941,1970)
  # read_echam4('EnSRF.ccc400_0', timlim=c(1601,2005), small=every2grid, landonly=land_only) # original
  read_echam4('EnSRF', timlim=c(1601,2004), small=every2grid, landonly=land_only)
} 

if (generate_ECHAM_anom){
  # Version 2 (newstatevector switch)
  # CHECK read_newCompo() function for Roni's code to generate newstatevector 
  # where vars can be choosen
  # probably not fully working and some var names need to be un/commented
  
  # Version 1 (oldstatevector switch)
  # read echam 71yr anom, clim and sd calculated with cdo from .nc files to .RData
  # for 60 ensm: timlim=c(1941,1970)
  # print("generate_ECHAM_anom") 
  # print("ATTENTION: manually copy generated .Rdata file to climstor!")
  # read_echam4('EnSRF', path=echanompath, timlim=c(1601,2004), small=every2grid, # originally was 'ano', path=echanompath, timlim=c(1601,2005)
  #             landonly=land_only, anom=T)
  # read_echam4('EnSRF', path=echclimpath, timlim=c(1601,2004), small=every2grid, # originally was timlim=c(1635,1970)
  #             landonly=land_only, clim=T)
  # read_echam4('EnSRF', path=echsdpath, timlim=c(1601,2005), small=every2grid, 
  #             landonly=land_only, std=T)
}

# if (generate_ECHAM_1901_70){
#   print("generate_ECHAM_1901_70")
#   # ECHAM data for bias calculation with real proxy data
#   echam1901_70 <- read_echam_ensmean('EnSRF', timlim=c(1901,1970),small=F)
#   save(echam1901_70, file="../data/echam_1911-70.Rdata")
# } 

if (generate_ECHAM_103){
  print("generate_ECHAM ens. mem. 103")
  print("ATTENTION: manually copy generated .Rdata file to climstor!")
  # ens. member WITHOUT land surface bug
  read_echam4(filehead='EnSRF', path=paste0(ext,'echam/echam103'), timlim=c(syr,eyr), 
              small=every2grid, landonly=land_only)
} 

if (generate_ECHAM_covar){
  print("generate_ECHAM all time step array long-term covariance")
  print("ATTENTION: manually copy generated .Rdata file to climstor!")
  echam_covar(syr=1603,eyr=2004)
}

if (generate_CCSM_last_mill_ens){
  print("generate CCSM last millennium ensemble")
  print("ATTENTION: manually copy generated .Rdata file to climstor!")
  # following filehead including "001" just reads first member
  # only works with landonly=F
  read_last_mill_ens(filehead="ncar_last_mill_001")
}

if (generate_ind_recon){
  print("generate Broennimann et al. 2009 atm. indices")
  print("ATTENTION: manually copy generated .Rdata file to climstor!")
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
  print("ATTENTION: manually copy generated .Rdata file to climstor!")
  #  see script in EnSRF/script/merge_cru.sh for regridding and co of orig. cru data set
  # cruall <- read_echam1('cru_allvar_abs_1901-2004_v2019.nc',timlim=c(syr_cru,eyr_cru),
  #                       path=crupath,small=every2grid,landonly=land_only)
  # load 1º BEST temp instead for 5º CRU 
  cruall <- read_echam1('cru_allvar_abs_1901-2004_berkeley2019.nc',timlim=c(syr_cru,eyr_cru),
                        path=bestpath,small=every2grid,landonly=land_only)
  if (every2grid) {
    # save(cruall, file=paste0(dataintdir,"cru/cru_allvar_",syr_cru,"-",eyr_cru,"_2ndgrid.Rdata"))
    # load 1º BEST temp instead for 5º CRU 
    save(cruall, file=paste0(dataintdir,"data/cru_allvar_",syr_cru,"-",eyr_cru,"_berkeley2019_2ndgrid.Rdata"))
    print("copy validation file manually to climstor!")
  } else {
    # save(cruall, file=paste0(dataintdir,"cru/cru_allvar_",syr_cru,"-",eyr_cru,".Rdata"))  
    # load 1º BEST temp instead for 5º CRU 
    save(cruall, file=paste0(dataintdir,"data/cru_allvar_",syr_cru,"-",eyr_cru,"_berkeley2019.Rdata"))  
    print("copy validation file manually to climstor!")
  }
}

if (generate_HadCRU4SD){ 
  print("generate_HadCRU4 standard deviations")
  print("ATTENTION: manually copy generated .Rdata file to climstor!")
  # from instr. CRU ensemble, precalculated with cdo
  cru4_may_sep <- read_echam1('HadCRUT4_ens_sd_may-sep_yrmean.nc',
                              timlim=c(1901,2002),path=crupath,small=every2grid,landonly=F)
  cru4_oct_apr <- read_echam1('HadCRUT4_ens_sd_oct-apr_yrmean.nc',
                              timlim=c(1901,2002),path=crupath,small=every2grid,landonly=F)
  cru4_may_sep$ensmean <- cru4_may_sep$data
  cru4_oct_apr$ensmean <- cru4_oct_apr$data
  if (every2grid) {
    save(cru4_may_sep,cru4_oct_apr, file="cru4_ens_sd_2ndgrid.Rdata")
  } else {
    save(cru4_may_sep,cru4_oct_apr, file="cru4_ens_sd.Rdata")  
  }
}

# old version 2
# if (generate_20CR) {
#   twentycr.all <- read_20cr(filehead="twentycr_allvar_1901-2004.nc",path=twentycrv2path,xlim=c(-180,180), 
#                             ylim=c(-90,90), timlim=c(1901,2004), small=every2grid, landonly=F, 
#                             calc_ensmean=T)
#   save(twentycr.all,file='twentycrv2_allvar_1901-2004.Rdata')
# }

# new version 3
if (generate_20CR) {
  twentycr.all <- read_20cr(filehead="20CRv3_1836-2015.nc",path=twentycrpath,xlim=c(-180,180), 
                            ylim=c(-90,90), timlim=c(1836,2004), small=every2grid, landonly=F, 
                            calc_ensmean=F)
  if (every2grid) {
    save(twentycr.all,file='twentycrv3_allvar_1836-2004_2ndgrid.Rdata')
  } else {
    save(twentycr.all,file='twentycrv3_allvar_1836-2004.Rdata')
  }
}

# if (generate_LUTPAULKUT){ 
#   print("generate_LUTPAULKUT")
#   # ATTENTION: seasonal resolution
#   #  see script in EnSRF/script/merge_recon.sh for regridding and co of orig. recon data set
#   reconall <- read_echam1('recon_allvar_1750-1999',xlim=c(-180,180), ylim=c(-90,90), 
#                           timlim=c(syr_recon,eyr_recon),path=reconpath,small=every2grid,
#                           landonly=land_only)
#   reconall$data[reconall$names=="precip",,]<-reconall$data[reconall$names=="precip",,]/3
#   if (every2grid) {
#     save(reconall, file=paste0(dataintdir,"recon/recon_allvar_",syr_recon,"-",eyr_recon,"_2ndgrid.Rdata"))
#   } else {
#     save(reconall, file=paste0(dataintdir,"recon/recon_allvar_",syr_recon,"-",eyr_recon,".Rdata"))  
#   }
# }

if (generate_GHCN){
  print("generate_GHCN")
  print("ATTENTION: manually copy generated .Rdata file to climstor!")
  # created only once for the full period and then cut after loading
  ghcn <- read_ghcn_refyr(1600,2005,1600,statyr)
  # create table of assimilated data including first and last year with data
  fyr <- lyr <- vector()
  for (i in 1:ncol(ghcn$data)) {
    pos <- which(!is.na(ghcn$data[,i]))
    fyr[i] <- floor(ghcn$time[pos[1]])
    lyr[i] <- floor(ghcn$time[pos[length(pos)]])
  }
  ghcn_tab <- cbind(ghcn$id,ghcn$names,ghcn$lon,ghcn$lat,fyr,lyr)
  write.table(ghcn_tab,file=paste0('EKF400_v2_assim_GHCN_freeze',statyr,'.txt'))
  
  ghcn$names <-rep('temp2',length(ghcn$names))
  save(ghcn, file=paste0("ghcn_temp_freeze",statyr,"_",fsyr,"-",feyr,".Rdata"))
}

if (generate_GHCN_precip){
  print("generate_GHCN_precip")
  print("ATTENTION: manually copy generated .Rdata file to climstor!")
  # created only once for the full period and then cut after loading
  ghcn_precip <- read_ghcn_refyr_precip(1600,2005,1600,statyr)
  ghcn_precip$data <- ghcn_precip$data / 10 # to make units echam conform
  ghcn_precip$names <-rep('precip',length(ghcn_precip$names))
  # old: save(ghcn_precip, file=paste0("../assim_data/ghcn/ghcn_precip_",fsyr,"-",feyr,".Rdata")) 
  save(ghcn_precip, file=paste0("ghcn_precip_freeze",statyr,"_",fsyr,"-",feyr,".Rdata"))
}

if (generate_ISTI){
  print("generate_ISTI")
  print("ATTENTION: manually copy generated .Rdata file to climstor!")
  source("convert_ISTI_to_EKF400.R")
}

if (generate_instr_yuri){
  print("generate_t_yuri")
  print("ATTENTION: manually copy generated .Rdata file to climstor!")
  source(paste0(dataextdir,"assim_data/v1/t/read_all.R"))
  source('read_instr_monthly_v2_2019.R')
  source(paste0(dataextdir,"assim_data/v1/slp/read_all.R"))
  source(paste0(dataextdir,"assim_data/v1/prec/read_all.R"))
}

if (generate_DOCUM){
  print("generate_DOCUM")
  # new version 2 with angie's collection
  source('read_docum_monthly_v2_2019.R')
  print("ATTENTION: manually copy generated .Rdata file to climstor!")
  # old version 1
  #source(paste0(dataextdir,"assim_data/data_yuri/t_docu/read_seas.R"))
  #source(paste0(dataextdir,"assim_data/data_yuri/t_docu/read_monthly.R"))
  #source(paste0(dataextdir,"assim_data/data_yuri/t_docu/read_JFMA.R"))
  #source(paste0(dataextdir,"assim_data/data_yuri/t_docu/read_AMJJA.R"))
}

# generate PAGES v2 from Raphi's 2018/01 export (version 1.6.1 from Raphi's data base):
#   /tank/exports/giub/EKF400/assimil_data/proxies/PAGES/read_pages_2018.R # Skript from Veronika

# generate NTREND
#   /tank/exports/giub/EKF400/assimil_data/proxies/NTREND/N-Trend.R        # Skript from Veronika

# generate petra, schweingr, mxd, etc. missing

if (generate_PROXIES){
  print("generate_PROXIES")
  print("ATTENTION: manually copy generated .Rdata file to climstor!")
  #read.these <- c("trw","mxd","schweingr","ntrend","pages","trw_petra")[c(TRW,MXD,SCHWEINGR,NTREND,PAGES,TRW_PETRA)]
  read.these <- c("ntrend","pages","trw_petra","mxd","schweingr","trw")[c(NTREND,PAGES,TRW_PETRA,MXD,SCHWEINGR,TRW)]
  if(exists("realprox")){rm(realprox)}
  for (varname in read.these){
    if (varname=="trw") {
      print("reading Petra's 35 best TRW records used in first EKF400 version")
      trwprox <- read_proxy2(fsyr,feyr)
      trwprox$archivetype <- rep('tree',length(trwproxd$lon))
      trwprox$datasource <- rep('mxdprox',length(trwprox$lon))
      realprox<-trwprox
      # NEW VERSION VERONIKA: PLEASE CHECK IF YOUR VERSION IS CORRECT
      ############# # because the dim is not equal, the year for trwprox stops in 1970, other two in 2005
      realprox$data <- rbind(trwprox$data,matrix(data=NA, nrow=length(seq(1971,2004)), ncol=dim(trwprox$data)[2])) # 34 not very elegant
      #############
      realprox$time <- c(realprox$time,seq(1971,2004))
    }
    
    if (varname=="mxd") {
      print("reading mxd")
      mxdprox <- read_proxy_mxd(fsyr,feyr,validate=lm_fit_data)
      mxdprox$archivetype <- rep('tree',length(mxdprox$lon))
      mxdprox$datasource <- rep('mxdprox',length(mxdprox$lon))
      
      if (exists("realprox")){
        realprox$data <- cbind(realprox$data, mxdprox$data)
        realprox$lon <- c(realprox$lon, mxdprox$lon)
        realprox$lat <- c(realprox$lat, mxdprox$lat)
        realprox$mr <- rbind(realprox$mr, mxdprox$mr)
        realprox$var_residu <- rbind(realprox$var_residu, mxdprox$var_residu)
        realprox$archivetype <- c(realprox$archivetype,mxdprox$archivetype)
        realprox$datasource <- c(realprox$datasource,mxdprox$datasource)
        
      } else { realprox<-mxdprox}
    }
    
    if (varname=="schweingr") {
      
      print("reading schweingr")
      schprox <- read_proxy_schweingr(fsyr,feyr,validate=lm_fit_data)
      schprox$archivetype <- rep('tree',length(schprox$lon))
      schprox$datasource <- rep('schweingruber_mxd',length(schprox$lon))
      
      if (exists("realprox")){
        
        realprox$data <- cbind(realprox$data, schprox$data)
        realprox$lon <- c(realprox$lon, schprox$lon)
        realprox$lat <- c(realprox$lat, schprox$lat)
        realprox$mr <- rbind(realprox$mr, schprox$mr)
        realprox$var_residu <- rbind(realprox$var_residu, schprox$var_residu)
        realprox$archivetype <- c(realprox$archivetype,schprox$archivetype)
        realprox$datasource <- c(realprox$datasource,schprox$datasource)
        
      } else { realprox<-schprox}
    }
    
    if (varname=="ntrend") {
      print("reading ntrend")
      ntrend = read_ntrend(fsyr,feyr, validate=lm_fit_data)
      ntrend$archivetype <- rep('tree',length(ntrend$lon))
      ntrend$datasource <- rep('ntrend',length(ntrend$lon))
      
      if (exists("realprox")){
        realprox$data <- cbind(realprox$data, ntrend$data)
        realprox$lon <- c(realprox$lon, ntrend$lon)
        realprox$lat <- c(realprox$lat, ntrend$lat)
        realprox$mr <- rbind(realprox$mr, ntrend$mr)
        realprox$var_residu <- rbind(realprox$var_residu, ntrend$var_residu)
        realprox$archivetype <- c(realprox$archivetype,ntrend$archivetype)
        realprox$datasource <- c(realprox$datasource,ntrend$datasource)
        
      } else { realprox<-ntrend}
    }
    
    if (varname=="pages") {
      print("reading pages")
      pagesprox <- setup_read_pages(type)
      pagesprox$datasource <- rep('pages_db',length(pagesprox$lon))
      
      if (exists("realprox")){
        
        realprox$data <- cbind(realprox$data, pagesprox$data)
        realprox$lon <- c(realprox$lon, pagesprox$lon)
        realprox$lat <- c(realprox$lat, pagesprox$lat)
        realprox$mr <- rbind(realprox$mr, pagesprox$mr)
        realprox$var_residu <- rbind(realprox$var_residu, pagesprox$var_residu)
        realprox$archivetype <- c(realprox$archivetype,pagesprox$archivetype)
        realprox$datasource <- c(realprox$datasource,pagesprox$datasource)
        
      } else { realprox<-pagesprox}
     
    }
    
    if (varname=="trw_petra") {
      print("reading trw_petra")
      trw_petra <- read_trw_petra(fsyr,feyr, validate=lm_fit_data) 
      trw_petra$archivetype <- rep('tree',length(trw_petra$lon))
      trw_petra$datasource <- rep('petra',length(trw_petra$lon))
      
      if (exists("realprox")){
        
        realprox$data <- cbind(realprox$data, trw_petra$data)
        realprox$lon <- c(realprox$lon, trw_petra$lon)
        realprox$lat <- c(realprox$lat, trw_petra$lat)
        realprox$mr <- rbind(realprox$mr, trw_petra$mr)
        realprox$var_residu <- rbind(realprox$var_residu, trw_petra$var_residu)
        realprox$archivetype <- c(realprox$archivetype,trw_petra$archivetype)
        realprox$datasource <- c(realprox$datasource,trw_petra$datasource)
        
      } else { realprox<-trw_petra}
    }
  }
  save(realprox, file=paste0("../data/proxies/real_proxies_",expname,"_",fsyr,"-",feyr,".Rdata"))
}

if (pseudo_prox) {
  if (generate_PSEUDO){
    print("generate_PSEUDO")
    print("ATTENTION: manually copy generated .Rdata file to climstor!")
    pseudoprox<-read_pseudo()
    realprox<-pseudoprox
    save(realprox, file=paste0("../data/proxies/DAPS_pseudoproxies_",fsyr,"-",feyr,".Rdata"))
  # } else {
  #   # load DAPS pseudo proxies in object called 'realprox'
  #   load(file=paste0("../data/proxies/DAPS_pseudoproxies_",fsyr,"-",feyr,".Rdata"))
  }  
}

