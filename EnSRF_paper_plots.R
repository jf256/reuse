rm(list=ls()) 

machine="climcal" #"macbook" #"climcal3" # "climcal"
if (machine=="macbook") {
  workdir='~/unibe/projects/EnSRF/src/'
  dataextdir='/Volumes/DATA/climdata/EKF400/'
  dataintdir=paste0(workdir,'../data/')
} else {
  workdir='/scratch3/joerg/projects/reuse/reuse_git/' #src/'
  dataextdir='/mnt/climstor/giub/EKF400/'
  dataintdir=paste0(workdir,'../data/')
}
setwd(workdir)

syr=1603; syr2=syr+1  # startyear e.g. 1602
eyr=2003

source('EnSRF_switches.R')
source('EnSRF_functions.R')
suppressMessages(library(lubridate)) # decinmal year to date conversion
suppressMessages(library(birk)) # which.closest function


save_prepplot=F
# load already transformed data to seasonal resolution
  load_precalc=F 
mergetime=F
load_prepplot=T
# monthly or seasonal analysis
if (monthly_out) {
  prepplotdir=paste0("../data/prepplot/",expname,'/prepplot_monthly/')
} else {
  prepplotdir=paste0("../data/prepplot/",expname,'/prepplot_seasonal/') 
}
prepplotdirseas=paste0("../data/prepplot/",expname,'/prepplot_seasonal/') 
prepplotdirmon=paste0("../data/prepplot/",expname,'/prepplot_monthly/')
dir.create(paste0('../data/indices/',expname))
dir.create(paste0('../figures/',expname))

timeseriesplots=F
  ind_anom=T # indices calculated based on 70yr anomalies
  land_only=F # overwrite switch
drought_1790=F
plots1816=F
multiann=F            # check multiannual drought, not finished but see code below
glacier_advances=F    # check period 1830-50 of great glacier_advances
warmcold_decades=F    # warm, cold periods global
warmcold_us=F         # warm, cold periods US
warmcold_eu=F         # warm, cold periods EU
warmcold_nh=F         # warm, cold periods NH
dryhumid_us=F         # dry, humid periods US
dryhumid_eu=F         # dry, humid periods EU
warmcold_1790_1830=F  # warm, cold periods NH
landusebugbias=F      # check bias due to land use bug, see also Marco's analysis
arctic_slp_trend=F    # check for trend in arctic as seen by Vavrus 2013 (http://onlinelibrary.wiley.com/doi/10.1002/2013GL058161/epdf)
updowntrend_decades=F # plot periods with steepest decadal trends
laki=F                # analyze Laki eruption 1783 in Island 
#cold_winters=F        # analyze cold winters in 19th century eruption
iceland=T             # look at eruptions in iceland only


if (drought_1790) {
  syr <- 1790
  eyr <- 1799  
} else if (plots1816) {
  syr <- 1816
  eyr <- 1817  
#} else {
#  syr=1603     
#  eyr=2004
}

if (ind_anom) {
  if (land_only) {
    if (monthly) {
      filenameext <- paste0('_anom_landonly_mon_')  
    } else {
      filenameext <- paste0('_anom_landonly_seas_')  
    }
  } else {
    if (monthly) {
      filenameext <- paste0('_anom_landocean_mon_')  
    } else {
      filenameext <- paste0('_anom_landocean_seas_')  
    }
  } 
} else {
  if (land_only) {
    if (monthly) {
      filenameext <- paste0('_abs_landonly_mon_')  
    } else {
      filenameext <- paste0('_abs_landonly_seas_')  
    }
  } else {
    if (monthly) {
      filenameext <- paste0('_abs_landocean_mon_')  
    } else {
      filenameext <- paste0('_abs_landocean_seas_')  
    }
  }
}
print(filenameext)








##########################################################################################
# start timeseriesplots
##########################################################################################
if (timeseriesplots) { 
  
  if (save_prepplot) {
    ptm1 <- proc.time()
    for (cyr in syr:eyr) {
#      if (cyr < 1902) {
      vali=F                 # switch off prepplot if no vali data selected
#      } else {
#        vali=T
#      }
      recon_vali=F
#       if ((cyr < 1751) | (cyr > 1979)) {
#         vali=F                 # switch off prepplot if no vali data selected
#       } else {
#         vali=T
#       }
#       if ((cyr < 1901) & (cyr > 1750)) {
#         recon_vali=T           # seasonal luterbacher, pauling, kuettel recons (1750-1999)
#       } else {
#         recon_vali=F
#       }
      
      print(cyr)
      if (load_precalc) { # load already calculated seasonal data
        load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
        echam.abs <- echam
        analysis.abs <- analysis
        bemlist <- read.table(file='../data/bem/bem.txt',header=F)
        bem <- bemlist[which(bemlist[,1]==cyr),2]
        echam.abs$bem <- echam.abs$data[,,bem]
        echam.anom$bem <- echam.anom$data[,,bem]
        analysis.abs$bem <- analysis.abs$data[,,bem]
        analysis.anom$bem <- analysis.anom$data[,,bem]
#        ech_ind$bem <- ech_ind$data[,,bem]
#        ana_ind$bem <- ana_ind$data[,,bem]
      } else {   
        load(file=paste0('../data/analysis/',expname,'analysis_',cyr,'.Rdata'))
        bemlist <- read.table(file='../data/bem/bem.txt',header=F)
        bem <- bemlist[which(bemlist[,1]==cyr),2]
        echam.abs$bem <- echam.abs$data[,,bem]
        echam.anom$bem <- echam.anom$data[,,bem]
        analysis.abs$bem <- analysis.abs$data[,,bem]
        analysis.anom$bem <- analysis.anom$data[,,bem]
        #  if (anomaly_assim){
        #    echam <- echam.abs
        #    analysis <- analysis.abs
        #    echam.abs <- NULL
        #    analysis.abs <- NULL   
        #   }
        
        # if (ind_ECHAM) {
        #   # extract the indices from output 
        #   # 1. Analysis indices
        #   ech_ind <- echam.anom
        #   ech_ind$data <- echam.anom$data[which(echam.anom$names=="DIMI" | echam.anom$names=="z100" | 
        #     echam.anom$names=="z300" | echam.anom$names=="PWC" | echam.anom$names=="HC" | 
        #     echam.anom$names=="SJ"),,]
        #   ech_ind$ensmean <- echam.anom$ensmean[which(echam.anom$names=="DIMI" | 
        #     echam.anom$names=="z100" | echam.anom$names=="z300" | echam.anom$names=="PWC" | 
        #     echam.anom$names=="HC" |  echam.anom$names=="SJ"),]
        #   ech_ind$bem <- echam.anom$bem[which(echam.anom$names=="DIMI" | echam.anom$names=="z100" | 
        #     echam.anom$names=="z300" | echam.anom$names=="PWC" | echam.anom$names=="HC" | 
        #     echam.anom$names=="SJ"),]
        #   ech_ind$names <- echam.anom$names[which(echam.anom$names=="DIMI" | echam.anom$names=="z100" | 
        #     echam.anom$names=="z300" | echam.anom$names=="PWC" | echam.anom$names=="HC" | 
        #     echam.anom$names=="SJ")]
        #   ech_ind$lon <- echam.anom$lon[which(echam.anom$names=="DIMI" | echam.anom$names=="z100" | 
        #     echam.anom$names=="z300" | echam.anom$names=="PWC" | echam.anom$names=="HC" | 
        #     echam.anom$names=="SJ")]
        #   ech_ind$lat <- echam.anom$lat[which(echam.anom$names=="DIMI" | echam.anom$names=="z100" | 
        #     echam.anom$names=="z300" | echam.anom$names=="PWC" | echam.anom$names=="HC" | 
        #     echam.anom$names=="SJ")]
        #   ana_ind <- analysis.anom
        #   ana_ind$data <- analysis.anom$data[which(analysis.anom$names=="DIMI" | 
        #     analysis.anom$names=="z100" | analysis.anom$names=="z300" | 
        #     analysis.anom$names=="PWC" | analysis.anom$names=="HC" | 
        #     analysis.anom$names=="SJ"),,]
        #   ana_ind$ensmean <- analysis.anom$ensmean[which(analysis.anom$names=="DIMI" | 
        #     analysis.anom$names=="z100" | analysis.anom$names=="z300" | 
        #     analysis.anom$names=="PWC" | analysis.anom$names=="HC" | 
        #     analysis.anom$names=="SJ"),]
        #   ana_ind$bem <- analysis.anom$bem[which(analysis.anom$names=="DIMI" | 
        #     analysis.anom$names=="z100" | analysis.anom$names=="z300" | 
        #     analysis.anom$names=="PWC" | analysis.anom$names=="HC" | 
        #     analysis.anom$names=="SJ"),]
        #   ana_ind$names <- analysis.anom$names[which(analysis.anom$names=="DIMI" | 
        #     analysis.anom$names=="z100" | analysis.anom$names=="z300" | 
        #     analysis.anom$names=="PWC" | analysis.anom$names=="HC" | 
        #     analysis.anom$names=="SJ")]
        #   ana_ind$lon <- analysis.anom$lon[which(analysis.anom$names=="DIMI" | 
        #     analysis.anom$names=="z100" | analysis.anom$names=="z300" | 
        #     analysis.anom$names=="PWC" | analysis.anom$names=="HC" | 
        #     analysis.anom$names=="SJ")]
        #   ana_ind$lat <- analysis.anom$lat[which(analysis.anom$names=="DIMI" | 
        #     analysis.anom$names=="z100" | analysis.anom$names=="z300" | 
        #     analysis.anom$names=="PWC" | analysis.anom$names=="HC" | 
        #     analysis.anom$names=="SJ")]
        # }
        
        
        # remove precalc indices from echam
        #  echam <- echam.anom
        echam.anom$data <- echam.anom$data[which(echam.anom$names!="DIMI" & 
          echam.anom$names!="z100" & echam.anom$names!="z300" & echam.anom$names!="PWC" & 
            echam.anom$names!="HC" & echam.anom$names!="SJ"),,]
        echam.anom$ensmean <- echam.anom$ensmean[which(echam.anom$names!="DIMI" & 
          echam.anom$names!="z100" & echam.anom$names!="z300" & echam.anom$names!="PWC" & 
          echam.anom$names!="HC" & echam.anom$names!="SJ"),]
        echam.anom$bem <- echam.anom$bem[which(echam.anom$names!="DIMI" & 
          echam.anom$names!="z100" & echam.anom$names!="z300" & echam.anom$names!="PWC" & 
          echam.anom$names!="HC" & echam.anom$names!="SJ"),]
        echam.anom$names <- echam.anom$names[which(echam.anom$names!="DIMI" & 
          echam.anom$names!="z100" & echam.anom$names!="z300" & echam.anom$names!="PWC" & 
          echam.anom$names!="HC" & echam.anom$names!="SJ")]
        echam.anom$lon <- echam.anom$lon[which(echam.anom$names!="DIMI" & 
          echam.anom$names!="z100" & echam.anom$names!="z300" & echam.anom$names!="PWC" & 
            echam.anom$names!="HC" & echam.anom$names!="SJ")]
        echam.anom$lat <- echam.anom$lat[which(echam.anom$names!="DIMI" & 
           echam.anom$names!="z100" & echam.anom$names!="z300" & echam.anom$names!="PWC" & 
             echam.anom$names!="HC" & echam.anom$names!="SJ")]
        #  analysis <- analysis.anom
        analysis.anom$data <- analysis.anom$data[which(analysis.anom$names!="DIMI" & 
          analysis.anom$names!="z100" & analysis.anom$names!="z300" & 
          analysis.anom$names!="PWC" & analysis.anom$names!="HC" & analysis.anom$names!="SJ"),,]
        analysis.anom$ensmean <- analysis.anom$ensmean[which(analysis.anom$names!="DIMI" & 
          analysis.anom$names!="z100" & analysis.anom$names!="z300" & 
          analysis.anom$names!="PWC" & analysis.anom$names!="HC" & analysis.anom$names!="SJ"),]
        analysis.anom$bem <- analysis.anom$bem[which(analysis.anom$names!="DIMI" & 
          analysis.anom$names!="z100" & analysis.anom$names!="z300" & 
          analysis.anom$names!="PWC" & analysis.anom$names!="HC" & analysis.anom$names!="SJ"),]
        analysis.anom$names <- analysis.anom$names[which(analysis.anom$names!="DIMI" & 
          analysis.anom$names!="z100" & analysis.anom$names!="z300" & 
          analysis.anom$names!="PWC" & analysis.anom$names!="HC" & analysis.anom$names!="SJ")]
        analysis.anom$lon <- analysis.anom$lon[which(analysis.anom$names!="DIMI" & 
          analysis.anom$names!="z100" & analysis.anom$names!="z300" & 
          analysis.anom$names!="PWC" & analysis.anom$names!="HC" & analysis.anom$names!="SJ")]
        analysis.anom$lat <- analysis.anom$lat[which(analysis.anom$names!="DIMI" & 
          analysis.anom$names!="z100" & analysis.anom$names!="z300" & 
          analysis.anom$names!="PWC" & analysis.anom$names!="HC" & analysis.anom$names!="SJ")]
        echam.abs$data <- echam.abs$data[which(echam.abs$names!="DIMI" & 
          echam.abs$names!="z100" & echam.abs$names!="z300" & echam.abs$names!="PWC" & 
          echam.abs$names!="HC" & echam.abs$names!="SJ"),,]
        echam.abs$ensmean <- echam.abs$ensmean[which(echam.abs$names!="DIMI" & 
          echam.abs$names!="z100" & echam.abs$names!="z300" & echam.abs$names!="PWC" & 
          echam.abs$names!="HC" & echam.abs$names!="SJ"),]
        echam.abs$bem <- echam.abs$bem[which(echam.abs$names!="DIMI" & echam.abs$names!="z100" & 
          echam.abs$names!="z300" & echam.abs$names!="PWC" & echam.abs$names!="HC" & 
          echam.abs$names!="SJ"),]
        echam.abs$names <- echam.abs$names[which(echam.abs$names!="DIMI" & 
          echam.abs$names!="z100" & echam.abs$names!="z300" & echam.abs$names!="PWC" & 
          echam.abs$names!="HC" & echam.abs$names!="SJ")]
        echam.abs$lon <- echam.abs$lon[which(echam.abs$names!="DIMI" & echam.abs$names!="z100" & 
          echam.abs$names!="z300" & echam.abs$names!="PWC" & echam.abs$names!="HC" & 
          echam.abs$names!="SJ")]
        echam.abs$lat <- echam.abs$lat[which(echam.abs$names!="DIMI" & echam.abs$names!="z100" & 
          echam.abs$names!="z300" & echam.abs$names!="PWC" & echam.abs$names!="HC" & 
          echam.abs$names!="SJ")]
        analysis.abs$data <- analysis.abs$data[which(analysis.abs$names!="DIMI" & 
          analysis.abs$names!="z100" & analysis.abs$names!="z300" & 
          analysis.abs$names!="PWC" & analysis.abs$names!="HC" & analysis.abs$names!="SJ"),,]
        analysis.abs$ensmean <- analysis.abs$ensmean[which(analysis.abs$names!="DIMI" & 
          analysis.abs$names!="z100" & analysis.abs$names!="z300" & 
          analysis.abs$names!="PWC" & analysis.abs$names!="HC" & analysis.abs$names!="SJ"),]
        analysis.abs$bem <- analysis.abs$bem[which(analysis.abs$names!="DIMI" & 
          analysis.abs$names!="z100" & analysis.abs$names!="z300" & 
          analysis.abs$names!="PWC" & analysis.abs$names!="HC" & analysis.abs$names!="SJ"),]
        analysis.abs$names <- analysis.abs$names[which(analysis.abs$names!="DIMI" & 
          analysis.abs$names!="z100" & analysis.abs$names!="z300" & 
          analysis.abs$names!="PWC" & analysis.abs$names!="HC" & analysis.abs$names!="SJ")]
        analysis.abs$lon <- analysis.abs$lon[which(analysis.abs$names!="DIMI" & 
          analysis.abs$names!="z100" & analysis.abs$names!="z300" & 
          analysis.abs$names!="PWC" & analysis.abs$names!="HC" & analysis.abs$names!="SJ")]
        analysis.abs$lat <- analysis.abs$lat[which(analysis.abs$names!="DIMI" & 
          analysis.abs$names!="z100" & analysis.abs$names!="z300" & 
          analysis.abs$names!="PWC" & analysis.abs$names!="HC" & analysis.abs$names!="SJ")]
        if (vali) {  
          vali_ind <- validate
          vali_ind$data <- validate$data[which(validate$names=="ind_rec_dimi" | 
            validate$names=="ind_rec_z100" | validate$names=="ind_rec_z300" | 
            validate$names=="ind_rec_pwc" | validate$names=="ind_rec_hc" | 
            validate$names=="ind_rec_sj"),]
          vali_ind$names <- validate$names[which(validate$names=="ind_rec_dimi" | 
            validate$names=="ind_rec_z100" | validate$names=="ind_rec_z300" | 
            validate$names=="ind_rec_pwc" | validate$names=="ind_rec_hc" | 
            validate$names=="ind_rec_sj")]
          vali_ind$lon <- validate$lon[which(validate$names=="ind_rec_dimi" | 
            validate$names=="ind_rec_z100" | validate$names=="ind_rec_z300" | 
            validate$names=="ind_rec_pwc" | validate$names=="ind_rec_hc" | 
            validate$names=="ind_rec_sj")]
          vali_ind$lat <- validate$lat[which(validate$names=="ind_rec_dimi" | 
            validate$names=="ind_rec_z100" | validate$names=="ind_rec_z300" | 
            validate$names=="ind_rec_pwc" | validate$names=="ind_rec_hc" | 
            validate$names=="ind_rec_sj")]
          validate$data <- validate$data[which(validate$names=="temp2" | 
            validate$names=="precip" | validate$names=="slp"),]
          validate$names <- validate$names[which(validate$names=="temp2" | 
            validate$names=="precip" | validate$names=="slp")]
          validate$lon <- validate$lon[which(validate$names=="temp2" | 
            validate$names=="precip" | validate$names=="slp")]
          validate$lat <- validate$lat[which(validate$names=="temp2" | 
            validate$names=="precip" | validate$names=="slp")]
        }
      }
      
      if (land_only) {
        xlim=c(-180,180)
        ylim=c(-90,90)
        nc <- nc_open(paste(echmaskpath, 'landseamask.nc', sep='/'))
        lon <- nc$dim$lon$vals
        lat <- nc$dim$lat$vals
        lon[lon > 180] <- lon[lon > 180] - 360
        loi <- which(lon >= xlim[1] & lon <= xlim[2])
        lai <- which(lat >= ylim[1] & lat <= ylim[2])
        # mulc for reading each 2rd grid cell to avoid memory problems
        mulc <- floor(length(loi)/96)
        loi <- loi[seq(ceiling(mulc/2),length(loi),mulc)]
        lai <- lai[seq(ceiling(mulc/2), length(lai),mulc)]
        lsm <- ncvar_get(nc)[loi, lai]
        nc_close(nc)
        lonlist=round(rep(lon[loi], length(lai))[lsm > 0.5],digits=2)
        latlist=round(rep(lat[lai], each=length(loi))[lsm > 0.5],digits=2)
        landll <- paste(lonlist,latlist)
        
        e.abs.ll <- paste(round(echam.abs$lon,digits=2),round(echam.abs$lat,digits=2))
        landpos <- match(landll,e.abs.ll)
        landpos_init <- landpos
        if (!load_precalc) {  
          ngrid <- dim(echam.abs$data)[1]/(nseas/s)/length(unique(echam.abs$names))
          nrep <-  (nseas/s)*length(unique(echam.abs$names))
        } else {
          # NEW if prepplot data are loaded
          ngrid <- dim(echam.abs$data)[1]/length(unique(echam.abs$names))
          nrep <-  length(unique(echam.abs$names))
        }
        for (i in 1:(nrep-1)) {
          landpos <- c(landpos,(landpos_init+ngrid*i))
        }
        echam.anom$data <- echam.anom$data[landpos,,] 
        echam.anom$ensmean <- echam.anom$ensmean[landpos,] 
        echam.anom$bem <- echam.anom$bem[landpos,] 
        echam.anom$lon <- echam.anom$lon[landpos] 
        echam.anom$lat <- echam.anom$lat[landpos] 
        echam.anom$names <- echam.anom$names[landpos] 
        echam.abs$data <- echam.abs$data[landpos,,] 
        echam.abs$ensmean <- echam.abs$ensmean[landpos,] 
        echam.abs$bem <- echam.abs$bem[landpos,] 
        echam.abs$lon <- echam.abs$lon[landpos] 
        echam.abs$lat <- echam.abs$lat[landpos] 
        echam.abs$names <- echam.abs$names[landpos] 
        analysis.anom$data <- analysis.anom$data[landpos,,] 
        analysis.anom$ensmean <- analysis.anom$ensmean[landpos,] 
        analysis.anom$bem <- analysis.anom$bem[landpos,] 
        analysis.anom$lon <- analysis.anom$lon[landpos] 
        analysis.anom$lat <- analysis.anom$lat[landpos] 
        analysis.anom$names <- analysis.anom$names[landpos] 
        analysis.abs$data <- analysis.abs$data[landpos,,] 
        analysis.abs$ensmean <- analysis.abs$ensmean[landpos,] 
        analysis.abs$bem <- analysis.abs$bem[landpos,] 
        analysis.abs$lon <- analysis.abs$lon[landpos] 
        analysis.abs$lat <- analysis.abs$lat[landpos] 
        analysis.abs$names <- analysis.abs$names[landpos] 
      }

      
      if (!load_precalc) {    
        # change to seasonal resolution
        if (sixmonstatevector)
          #  echam.anom.mon <- echam.anom        
          #  analysis.anom.mon <- analysis.anom
          #  echam.abs.mon <- echam.abs        
          #  analysis.abs.mon <- analysis.abs
          #  validate.mon <- validate
          #  if (landonly) {
          #    for (p in 1:length(which(!is.na(echam$lat))))
          #    vpos <- which(!is.na(validate$data[,1]))
          #    Hind2 <- Hind[,vpos]
          #  }
          
        etmp <- array(echam.abs$data,c((dim(echam.abs$data)[1]/6),6,dim(echam.abs$data)[2],
                                       dim(echam.abs$data)[3]))
        echam.abs$data <- array(NA,c(dim(etmp)[1],dim(etmp)[3:4]))
        for (ensmem in 1:(dim(etmp)[4])) {
          print(paste('ECHAM member',ensmem))
          echam.abs$data[,,ensmem] <- apply(etmp[,,,ensmem],c(1,3),mean)
        }
        #      echam.abs$data <- apply(etmp,c(1,3,4),mean)
        etmp2 <- array(echam.abs$ensmean,c((dim(echam.abs$ensmean)[1]/6),6,dim(echam.abs$ensmean)[2]))
        echam.abs$ensmean <- apply(etmp2,c(1,3),mean)
        etmp3 <- array(echam.abs$bem,c((dim(echam.abs$bem)[1]/6),6,dim(echam.abs$bem)[2]))
        echam.abs$bem <- apply(etmp3,c(1,3),mean)
        #    echam.abs$names <- echam.abs$names[1:dim(echam.abs$data)[1]/length(unique(echam.abs$names))]
        echam.abs$names <- rep(unique(echam.abs$names),each=dim(echam.abs$data)[1]/
                                 length(unique(echam.abs$names)))  
        echam.abs$lon <- echam.abs$lon[1:(dim(echam.abs$data)[1])] #/length(unique(echam.abs$names)))]
        echam.abs$lat <- echam.abs$lat[1:(dim(echam.abs$data)[1])] #/length(unique(echam$names)))]
#         etmp <- array(ech_ind$data,c((dim(ech_ind$data)[1]/6),6,dim(ech_ind$data)[2],
#                                      dim(ech_ind$data)[3]))  
#         ech_ind$data <- apply(etmp,c(1,3,4),mean)
#         etmp2 <- array(ech_ind$ensmean,c((dim(ech_ind$ensmean)[1]/6),6,dim(ech_ind$ensmean)[2]))
#         ech_ind$ensmean <- apply(etmp2,c(1,3),mean)
#         etmp3 <- array(ech_ind$bem,c((dim(ech_ind$bem)[1]/6),6,dim(ech_ind$bem)[2]))
#         ech_ind$bem <- apply(etmp3,c(1,3),mean)
#         #    ech_ind$names <- ech_ind$names[1:dim(ech_ind$data)[1]/length(unique(ech_ind$names))]
#         ech_ind$names <- rep(unique(ech_ind$names),each=dim(ech_ind$data)[1]/
#                                length(unique(ech_ind$names)))  
#         ech_ind$lon <- ech_ind$lon[1:(dim(ech_ind$data)[1])] #/length(unique(ech_ind$names)))]
#         ech_ind$lat <- ech_ind$lat[1:(dim(ech_ind$data)[1])] #/length(unique(ech_ind$names)))] 
#         etmp <- NULL
#         etmp2 <- NULL
        
        atmp <- array(analysis.abs$data,c((dim(analysis.abs$data)[1]/6),6,dim(analysis.abs$data)[2],
                                          dim(analysis.abs$data)[3]))  
        analysis.abs$data <- array(NA,c(dim(atmp)[1],dim(atmp)[3:4]))
        for (ensmem in 1:(dim(atmp)[4])) {
          print(paste('Analysis member',ensmem))
          analysis.abs$data[,,ensmem] <- apply(atmp[,,,ensmem],c(1,3),mean)
        }
        #      analysis.abs$data <- apply(atmp,c(1,3,4),mean)
        #    analysis.abs$names <- analysis.abs$names[1:dim(analysis.abs$data)[1]]
        analysis.abs$names <- rep(unique(analysis.abs$names),each=dim(analysis.abs$data)[1]/
                                    length(unique(analysis.abs$names)))
        analysis.abs$lon <- analysis.abs$lon[1:(dim(analysis.abs$data)[1])] #/length(unique(analysis.abs$names)))]
        analysis.abs$lat <- analysis.abs$lat[1:(dim(analysis.abs$data)[1])] #/length(unique(analysis.abs$names)))]
        atmp2 <- array(analysis.abs$ensmean,c((dim(analysis.abs$ensmean)[1]/6),6,
                                              dim(analysis.abs$ensmean)[2]))
        analysis.abs$ensmean <- apply(atmp2,c(1,3),mean)
        atmp3 <- array(analysis.abs$bem,c((dim(analysis.abs$bem)[1]/6),6,
                                          dim(analysis.abs$bem)[2]))
        analysis.abs$bem <- apply(atmp3,c(1,3),mean)
#         atmp <- array(ana_ind$data,c((dim(ana_ind$data)[1]/6),6,dim(ana_ind$data)[2],
#                                      dim(ana_ind$data)[3]))  
#         ana_ind$data <- apply(atmp,c(1,3,4),mean)
#         #    ana_ind$names <- ana_ind$names[1:dim(ana_ind$data)[1]]
#         ana_ind$names <- rep(unique(ana_ind$names),each=dim(ana_ind$data)[1]/
#                                length(unique(ana_ind$names)))
#         ana_ind$lon <- ana_ind$lon[1:(dim(ana_ind$data)[1])] #/length(unique(ana_ind$names)))]
#         ana_ind$lat <- ana_ind$lat[1:(dim(ana_ind$data)[1])] #/length(unique(ana_ind$names)))]
#         atmp2 <- array(ana_ind$ensmean,c((dim(ana_ind$ensmean)[1]/6),6,
#                                          dim(ana_ind$ensmean)[2]))
#         ana_ind$ensmean <- apply(atmp2,c(1,3),mean)
#         atmp3 <- array(ana_ind$bem,c((dim(ana_ind$bem)[1]/6),6,
#                                      dim(ana_ind$bem)[2]))
#         ana_ind$bem <- apply(atmp3,c(1,3),mean)
#         atmp <- NULL
#         atmp2 <- NULL
        
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
        etmp3 <- array(echam.anom$bem,c((dim(echam.anom$bem)[1]/6),6,dim(echam.anom$bem)[2]))
        echam.anom$bem <- apply(etmp3,c(1,3),mean)
        #    echam.anom$names <- echam.anom$names[1:dim(echam.anom$data)[1]/length(unique(echam.anom$names))]
        echam.anom$names <- rep(unique(echam.anom$names),each=dim(echam.anom$data)[1]/
                                  length(unique(echam.anom$names)))  
        echam.anom$lon <- echam.anom$lon[1:(dim(echam.anom$data)[1])] #/length(unique(echam.anom$names)))]
        echam.anom$lat <- echam.anom$lat[1:(dim(echam.anom$data)[1])] #/length(unique(echam.anom$names)))]
#        etmp <- array(ech_ind$data,c((dim(ech_ind$data)[1]/6),6,dim(ech_ind$data)[2],
#                                     dim(ech_ind$data)[3]))  
#        etmp <- NULL
        
#       atmp <- array(analysis.anom$data,c((dim(analysis.anom$data)[1]/6),6,dim(analysis.anom$data)[2],
#                                           dim(analysis.anom$data)[3])) 
#        analysis.anom$data <- array(NA,c(dim(atmp)[1],dim(atmp)[3:4]))
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
        atmp3 <- array(analysis.anom$bem,c((dim(analysis.anom$bem)[1]/6),6,
                                           dim(analysis.anom$bem)[2]))
        analysis.anom$bem <- apply(atmp3,c(1,3),mean)
        atmp <- NULL
        
        if (vali) {    
          if (!recon_vali) {
            vtmp <- array(validate$data,c((dim(validate$data)[1]/6),6,dim(validate$data)[2]))  
            validate$data <- apply(vtmp,c(1,3),mean)
            validate$ensmean <- validate$data
            #      validate$names <- validate$names[1:dim(validate$data)[1]]
            validate$names <- rep(unique(validate$names),each=dim(validate$data)[1]/
                                    length(unique(validate$names)))
            validate$lon <- validate$lon[1:(dim(validate$data)[1])] #/length(unique(validate$names)))]
            validate$lat <- validate$lat[1:(dim(validate$data)[1])] #/length(unique(validate$names)))]
            vtmp <- array(vali_ind$data,c((dim(vali_ind$data)[1]/6),6,dim(vali_ind$data)[2]))  
            #      vali_ind=validate
            vali_ind$data <- apply(vtmp,c(1,3),mean)
            vali_ind$ensmean <- vali_ind$data
            #      vali_ind$names <- vali_ind$names[1:dim(vali_ind$data)[1]]
            vali_ind$names <- rep(unique(vali_ind$names),each=dim(vali_ind$data)[1]/
                                    length(unique(vali_ind$names)))
            vali_ind$lon <- vali_ind$lon[1:(dim(vali_ind$data)[1])] #/length(unique(vali_ind$names)))]
            vali_ind$lat <- vali_ind$lat[1:(dim(vali_ind$data)[1])] #/length(unique(vali_ind$names)))]
          }
        }
        if (!real_proxies) {
          # ERROR, NOT ADJUSTED FOR seasonal/annual DOCUMENTARY DATA, YET !!!
          ctmp <- array(calibrate$data,c((dim(calibrate$data)[1]/6),6,
                                         dim(calibrate$data)[2]))  
          calibrate$data <- apply(ctmp,c(1,3),mean)
          calibrate$names <- calibrate$names[1:dim(calibrate$data)[1]] 
          calibrate$lon <- calibrate$lon[1:dim(calibrate$data)[1]]
          calibrate$lat <- calibrate$lat[1:dim(calibrate$data)[1]]
        }
      }
      
      
      
      
      
      
      if (ind_anom) {
        echam <- echam.anom
        analysis <- analysis.anom
      } else {
        echam <- echam.abs
        analysis <- analysis.abs
      }    
      
      
      
      
      
      
      
      # calculate the indices from output 
      H.giorgi <- compute_giorgi_H_v2(giorgi, echam) #, numvar=3) # 3 vars temp, precip, slp
      #    H.giorgi.anom <- compute_giorgi_H_sixmon(giorgi, echam.anom) #, numvar=3) # 3 vars temp, precip, slp
      #    indices <- c('NH.temp2', 'NH.precip', 'NH.slp', 'NEU.temp2', 'NEU.precip', 'NEU.slp',
      #                 'DIMI', 'HC', 'SJ', 'Z100', 'Z300', 'PWC')
      #    Hind <- matrix(0,nrow=12,ncol=nrow(echam$data))
      indices_tmp <- c('ENH.temp2','NAM.temp2','SAM.temp2','AFR.temp2',
                       'ASI.temp2','AUS.temp2','ARC.temp2','ANT.temp2',
                       'NEU.temp2','MED.temp2',
                       'GLO.temp2','NAM.precip','SAM.precip','AFR.precip',
                       'ASI.precip','AUS.precip','ARC.precip','ANT.precip',
                       'NEU.precip','MED.precip',
                       'SH.temp2','NAM.slp','SAM.slp','AFR.slp',
                       'ASI.slp','AUS.slp','ARC.slp','ANT.slp',
                       'NEU.slp','MED.slp',
                       'NH.temp2', 'EU.temp2', 'EU.precip', 'EU.slp', 
                       'PV1.calc', 'PV2.calc', 'PWC1.calc', 'PWC2.calc', 
                       'DIMI1.calc', 'DIMI2.calc', 'NAO1.calc', 'NAO2.calc',  
                       'PNA1.calc', 'PNA2.calc', 'PNA3.calc', 'PNA4.calc')
      indices <- c('ENH.temp2','NAM.temp2','SAM.temp2','AFR.temp2',
                   'ASI.temp2','AUS.temp2','ARC.temp2','ANT.temp2',
                   'NEU.temp2','MED.temp2',
                   'GLO.temp2','NAM.precip','SAM.precip','AFR.precip',
                   'ASI.precip','AUS.precip','ARC.precip','ANT.precip',
                   'NEU.precip','MED.precip',
                   'SH.temp2','NAM.slp','SAM.slp','AFR.slp',
                   'ASI.slp','AUS.slp','ARC.slp','ANT.slp',
                   'NEU.slp','MED.slp',
                   'NH.temp2', 'EU.temp2', 'EU.precip', 'EU.slp', 
                   'HC.calc', 'ITCZ.calc', 'SJ_u200.calc', 'SJ_slp.calc', 'PV.calc', 
                   'PWC.calc', 'DIMI.calc', 'NAO.calc', 'PNA.calc') #,
#                   'DIMI', 'HC', 'SJ', 'Z100', 'Z300', 'PWC')
      # "MC" kein z300 in output! 'HCL.calc',
      # SO WIRD INDEX NUR MITTELWERT ¨UBER REGION UND ZEIT! 
      # Deshalb Indicies die min max etc basiert sind extra rechnen
      Hind <- matrix(0,nrow=length(indices_tmp),ncol=nrow(echam$data))
      Hind[1,which(echam$names=="temp2")] <- H.giorgi[which(giorgi.short == 'ENH'),
                                                          which(echam$names=="temp2")]
      Hind[2,which(echam$names=="temp2")] <- H.giorgi[which(giorgi.short == 'NAM'),
                                                          which(echam$names=="temp2")]
      Hind[3,which(echam$names=="temp2")] <- H.giorgi[which(giorgi.short == 'SAM'),
                                                          which(echam$names=="temp2")]
      Hind[4,which(echam$names=="temp2")] <- H.giorgi[which(giorgi.short == 'AFR'),
                                                          which(echam$names=="temp2")]
      Hind[5,which(echam$names=="temp2")] <- H.giorgi[which(giorgi.short == 'ASI'),
                                                          which(echam$names=="temp2")]
      Hind[6,which(echam$names=="temp2")] <- H.giorgi[which(giorgi.short == 'AUS'),
                                                          which(echam$names=="temp2")]
      Hind[7,which(echam$names=="temp2")] <- H.giorgi[which(giorgi.short == 'ARC'),
                                                          which(echam$names=="temp2")]
      Hind[8,which(echam$names=="temp2")] <- H.giorgi[which(giorgi.short == 'ANT'),
                                                          which(echam$names=="temp2")]
      Hind[9,which(echam$names=="temp2")] <- H.giorgi[which(giorgi.short == 'NEU'),
                                                          which(echam$names=="temp2")]
      Hind[10,which(echam$names=="temp2")] <- H.giorgi[which(giorgi.short == 'MED'),
                                                           which(echam$names=="temp2")]
      Hind[11,which(echam$names=="temp2")] <- H.giorgi[which(giorgi.short == 'GLO'),
                                                           which(echam$names=="temp2")]
      Hind[12,which(echam$names=="precip")] <- H.giorgi[which(giorgi.short == 'NAM'),
                                                            which(echam$names=="precip")]
      Hind[13,which(echam$names=="precip")] <- H.giorgi[which(giorgi.short == 'SAM'),
                                                            which(echam$names=="precip")]
      Hind[14,which(echam$names=="precip")] <- H.giorgi[which(giorgi.short == 'AFR'),
                                                            which(echam$names=="precip")]
      Hind[15,which(echam$names=="precip")] <- H.giorgi[which(giorgi.short == 'ASI'),
                                                            which(echam$names=="precip")]
      Hind[16,which(echam$names=="precip")] <- H.giorgi[which(giorgi.short == 'AUS'),
                                                            which(echam$names=="precip")]
      Hind[17,which(echam$names=="precip")] <- H.giorgi[which(giorgi.short == 'ARC'),
                                                            which(echam$names=="precip")]
      Hind[18,which(echam$names=="precip")] <- H.giorgi[which(giorgi.short == 'ANT'),
                                                            which(echam$names=="precip")]
      Hind[19,which(echam$names=="precip")] <- H.giorgi[which(giorgi.short == 'NEU'),
                                                            which(echam$names=="precip")]
      Hind[20,which(echam$names=="precip")] <- H.giorgi[which(giorgi.short == 'MED'),
                                                            which(echam$names=="precip")]
      Hind[21,which(echam$names=="temp2")] <- H.giorgi[which(giorgi.short == 'SH'),
                                                           which(echam$names=="temp2")]
      Hind[22,which(echam$names=="slp")] <- H.giorgi[which(giorgi.short == 'NAM'),
                                                         which(echam$names=="slp")]
      Hind[23,which(echam$names=="slp")] <- H.giorgi[which(giorgi.short == 'SAM'),
                                                         which(echam$names=="slp")]
      Hind[24,which(echam$names=="slp")] <- H.giorgi[which(giorgi.short == 'AFR'),
                                                         which(echam$names=="slp")]
      Hind[25,which(echam$names=="slp")] <- H.giorgi[which(giorgi.short == 'ASI'),
                                                         which(echam$names=="slp")]
      Hind[26,which(echam$names=="slp")] <- H.giorgi[which(giorgi.short == 'AUS'),
                                                         which(echam$names=="slp")]
      Hind[27,which(echam$names=="slp")] <- H.giorgi[which(giorgi.short == 'ARC'),
                                                         which(echam$names=="slp")]
      Hind[28,which(echam$names=="slp")] <- H.giorgi[which(giorgi.short == 'ANT'),
                                                         which(echam$names=="slp")]
      Hind[29,which(echam$names=="slp")] <- H.giorgi[which(giorgi.short == 'NEU'),
                                                         which(echam$names=="slp")]
      Hind[30,which(echam$names=="slp")] <- H.giorgi[which(giorgi.short == 'MED'),
                                                         which(echam$names=="slp")]  
      Hind[31,which(echam$names=="temp2")] <- H.giorgi[which(giorgi.short == 'NH'),
                                                           which(echam$names=="temp2")]
      Hind[32,which(echam$names=="temp2")] <- H.giorgi[which(giorgi.short == 'EU'),
                                                           which(echam$names=="temp2")]
      Hind[33,which(echam$names=="precip")] <- H.giorgi[which(giorgi.short == 'EU'),
                                                            which(echam$names=="precip")]
      Hind[34,which(echam$names=="slp")] <- H.giorgi[which(giorgi.short == 'EU'),
                                                         which(echam$names=="slp")] 
      # calc indices from model fields
      # midlatitude circulation
      #    Hind[5,which(echam$names=="gph300")] <- H.giorgi[which(giorgi.short == 'MC'),
      #                                              which(echam$names=="gph300")] 
      # cdo -s fldmean -sellonlatbox,0,360,30,60 -sellevel,3000 -selvar,geopoth $f ${STOR_DIR}/z300.$f
      
      # stratospheric polar vortex
      Hind[35,which(echam$names=="gph100")] <- H.giorgi[which(giorgi.short == 'PV1'),
                                                            which(echam$names=="gph100")] 
      Hind[36,which(echam$names=="gph100")] <- H.giorgi[which(giorgi.short == 'PV2'),
                                                            which(echam$names=="gph100")] 
      # cdo -s -sellevel,10000 -selvar,geopoth $f ${STOR_DIR}/gph100.$f
      # cdo -s sub -fldmean -sellonlatbox,0,360,75,90 ${STOR_DIR}/gph100.$f -fldmean -sellonlatbox,0,360,40,55 ${STOR_DIR}/gph100.$f ${STOR_DIR}/z100.$f
      
      # Pacific walker circulation
      Hind[37,which(echam$names=="omega500")] <- H.giorgi[which(giorgi.short == 'PWC1'),
                                                              which(echam$names=="omega500")] 
      # cdo -s sellevel,50000 -selvar,omega $f ${STOR_DIR}/omega.$f
      # cdo -s sub -fldmean -sellonlatbox,180,260,-10,10 ${STOR_DIR}/omega.$f -fldmean -sellonlatbox,100,150,-10,10 ${STOR_DIR}/omega.$f ${STOR_DIR}/PWC.$f
      
      # Pacific walker circulation
      Hind[38,which(echam$names=="omega500")] <- H.giorgi[which(giorgi.short == 'PWC2'),
                                                              which(echam$names=="omega500")] 
      # cdo -s sellevel,50000 -selvar,omega $f ${STOR_DIR}/omega.$f
      # cdo -s sub -fldmean -sellonlatbox,180,260,-10,10 ${STOR_DIR}/omega.$f -fldmean -sellonlatbox,100,150,-10,10 ${STOR_DIR}/omega.$f ${STOR_DIR}/PWC.$f
      
      # Dynamic Indian Monsoon index
      Hind[39,which(echam$names=="u850")] <- H.giorgi[which(giorgi.short == 'DIMI1'),
                                                          which(echam$names=="u850")] 
      # cdo -s sellevel,85000 -selvar,u $f ${STOR_DIR}/u850.$f
      # cdo -s sub -fldmean -sellonlatbox,40,80,5,15 ${STOR_DIR}/u850.$f -fldmean -sellonlatbox,70,90,20,30 ${STOR_DIR}/u850.$f ${STOR_DIR}/DIMI.$f
      
      # Dynamic Indian Monsoon index
      Hind[40,which(echam$names=="u850")] <- H.giorgi[which(giorgi.short == 'DIMI2'),
                                                          which(echam$names=="u850")] 
      # cdo -s sellevel,85000 -selvar,u $f ${STOR_DIR}/u850.$f
      # cdo -s sub -fldmean -sellonlatbox,40,80,5,15 ${STOR_DIR}/u850.$f -fldmean -sellonlatbox,70,90,20,30 ${STOR_DIR}/u850.$f ${STOR_DIR}/DIMI.$f
      
      # NAO index
      Hind[41,which(echam.anom$names=="slp")] <- H.giorgi[which(giorgi.short == 'NAO1'),
                                                          which(echam.anom$names=="slp")] 
      Hind[42,which(echam.anom$names=="slp")] <- H.giorgi[which(giorgi.short == 'NAO2'),
                                                          which(echam.anom$names=="slp")] 
      Hind[43,which(echam.anom$names=="gph500")] <- H.giorgi[which(giorgi.short == 'PNA1'),
                                                             which(echam.anom$names=="gph500")] 
      Hind[44,which(echam.anom$names=="gph500")] <- H.giorgi[which(giorgi.short == 'PNA2'),
                                                             which(echam.anom$names=="gph500")] 
      Hind[45,which(echam.anom$names=="gph500")] <- H.giorgi[which(giorgi.short == 'PNA3'),
                                                             which(echam.anom$names=="gph500")] 
      Hind[46,which(echam.anom$names=="gph500")] <- H.giorgi[which(giorgi.short == 'PNA4'),
                                                             which(echam.anom$names=="gph500")] 
      eind_tmp <- list(ensmean=Hind%*%echam$ensmean, names=indices_tmp)
      
      # Hadley cell strength
      #    Hind[5,which(echam$names=="stream")] <- H.giorgi[which(giorgi.short == 'HC'),
      #                                              which(echam$names=="stream")] 
      # Hind[5,which(echam$names=="omega500")] <- H.giorgi[which(giorgi.short == 'HC'),
      #                                                   which(echam$names=="omega500")] 
      c=1
      zmean <- matrix(0,nrow=length(unique(echam$lat)),ncol=2)
      for (l in unique(echam$lat)) {
        pos1 <- which(trunc(echam$lat,digits=8)==trunc(l,digits=8))
        pos2 <- which(echam$names=='omega500')
        pos <- intersect(pos1,pos2)
        if (length(echam$ensmean[pos,])==2) { zmean[c,] <- echam$ensmean[pos,] }
        if (length(echam$ensmean[pos,])>2) { zmean[c,] <- apply(echam$ensmean[pos,],2,mean) }
        c=c+1
      }
      latlist <- unique(echam$lat)[which(which(unique(echam$lat) > -30) %in% which(unique(echam$lat) < 30))]
      pos3 <- match(round(latlist,digits=8),round(unique(echam$lat),digits=8)) 
      zmean <- zmean[pos3,]
      hc <- apply(zmean,2,min) # max upward wind is neg. omega value -> min function
      # cdo -s fldmax -sellonlatbox,0,360,0,30 -zonmean -sellevel,50000 -selvar,stream $f ${STOR_DIR}/HC.$f
      
      # Hadley cell poleward extend
      # Hind[6,which(echam$names=="u200")] <- H.giorgi[which(giorgi.short == 'HCL'),
      #                                               which(echam$names=="u200")] 
      
      # ITCZ location
      #Hind[7,which(echam$names=="omega500")] <- H.giorgi[which(giorgi.short == 'ITCZ'),
      #                                                   which(echam$names=="omega500")] 
      c=1
      zmean <- matrix(0,nrow=length(unique(echam$lat)),ncol=2)
      for (l in unique(echam$lat)) {
        pos1 <- which(trunc(echam$lat,digits=8)==trunc(l,digits=8))
        pos2 <- which(echam$names=='omega500')
        pos <- intersect(pos1,pos2)
        if (length(echam$ensmean[pos,])==2) { zmean[c,] <- echam$ensmean[pos,] }
        if (length(echam$ensmean[pos,])>2) { zmean[c,] <- apply(echam$ensmean[pos,],2,mean) }
        c=c+1
      }
      latlist <- unique(echam$lat)[which(which(unique(echam$lat) > -30) %in% 
                    which(unique(echam$lat) < 30))]
      pos3 <- match(round(latlist,digits=8),round(unique(echam$lat),digits=8)) 
      zmean <- zmean[pos3,]
      itcz <- apply(zmean[2:(nrow(zmean)-1),],2,min) 
      dlat <-  mean(latlist[1:(length(latlist)-1)]-latlist[2:(length(latlist))])    
      gridminpos1 <- which(round(zmean[,1],digits=8)==round(itcz[1],digits=8)) 
      gridmin1 <- latlist[gridminpos1]
      gridminpos2 <- which(round(zmean[,2],digits=8)==round(itcz[2],digits=8)) 
      gridmin2 <- latlist[gridminpos2]
      latmin1 <- gridmin1 - dlat * (((zmean[(gridminpos1-1),1])-(zmean[(gridminpos1+1),1]))/
        (2*(abs(zmean[gridminpos1]-max(c(zmean[(gridminpos1-1),1],zmean[(gridminpos1+1),1]))))))
      latmin2 <- gridmin2 - dlat * (((zmean[(gridminpos2-1),2])-(zmean[(gridminpos2+1),2]))/
        (2*(abs(zmean[gridminpos2]-max(c(zmean[(gridminpos2-1),2],zmean[(gridminpos2+1),2]))))))
      if (length(latmin1) > 1) {latmin1 <- mean(latmin1)}
      if (length(latmin2) > 1) {latmin2 <- mean(latmin2)}
      itcz <- c(latmin1,latmin2)
      
      # subtropical jet from subtrop. jet (u200) max
      #Hind[8,which(echam$names=="u200")] <- H.giorgi[which(giorgi.short == 'SJ'),
      #                                               which(echam$names=="u200")] 
      c=1
      # zonal means:
      zmean <- matrix(0,nrow=length(unique(echam$lat)),ncol=2)
      for (l in unique(echam$lat)) {
        pos1 <- which(trunc(echam$lat,digits=8)==trunc(l,digits=8))
        pos2 <- which(echam$names=='u200')
        pos <- intersect(pos1,pos2)
        if (length(echam$ensmean[pos,])==2) { zmean[c,] <- echam$ensmean[pos,] }
        if (length(echam$ensmean[pos,])>2) { zmean[c,] <- apply(echam$ensmean[pos,],2,mean) }
        c=c+1
      }
      latlist <- unique(echam$lat)[which(which(unique(echam$lat) > 0) %in% which(unique(echam$lat) < 57))]
      pos3 <- match(round(latlist,digits=8),round(unique(echam$lat),digits=8)) 
      zmean <- zmean[pos3,]
      sj <- apply(zmean[2:(nrow(zmean)-1),],2,max) # max uwind
      dlat <-  mean(latlist[1:(length(latlist)-1)]-latlist[2:(length(latlist))])    
      # pos1 for season 1; would need adjustment for monthly data
      gridminpos1 <- which(round(zmean[,1],digits=8)==round(sj[1],digits=8)) 
      gridmin1 <- latlist[gridminpos1]
      # pos2 for season 2
      gridminpos2 <- which(round(zmean[,2],digits=8)==round(sj[2],digits=8)) 
      gridmin2 <- latlist[gridminpos2]
      zmean <- zmean * -1 # because next equation is to find minimum and here we search for max
      sj <- sj * -1
      latmin1 <- gridmin1 - dlat * (((zmean[(gridminpos1-1),1])-(zmean[(gridminpos1+1),1]))/
                                      (2*(abs(sj[1]-max(c(zmean[(gridminpos1-1),1],zmean[(gridminpos1+1),1]))))))
      latmin2 <- gridmin2 - dlat * (((zmean[(gridminpos2-1),2])-(zmean[(gridminpos2+1),2]))/
                                      (2*(abs(sj[2]-max(c(zmean[(gridminpos2-1),2],zmean[(gridminpos2+1),2]))))))
      if (length(latmin1) > 1) {latmin1 <- mean(latmin1)}
      if (length(latmin2) > 1) {latmin2 <- mean(latmin2)}
      sj_u200 <- c(latmin1,latmin2)
      # cdo -s fldmax -sellonlatbox,0,360,0,50 -zonmean -sellevel,20000 -selvar,u $f ${STOR_DIR}/SJ.$f
      
      # subtropical jet from zonal SLP max
      c=1
      zmean <- matrix(0,nrow=length(unique(echam$lat)),ncol=2)
      for (l in unique(echam$lat)) {
        pos1 <- which(trunc(echam$lat,digits=8)==trunc(l,digits=8))
        pos2 <- which(echam$names=='slp')
        pos <- intersect(pos1,pos2)
        if (length(echam$ensmean[pos,])==2) { zmean[c,] <- echam$ensmean[pos,] }
        if (length(echam$ensmean[pos,])>2) { zmean[c,] <- apply(echam$ensmean[pos,],2,mean) }
        c=c+1
      }
      latlist <- unique(echam$lat)[which(which(unique(echam$lat) > 0) %in% 
                    which(unique(echam$lat) < 57))]
      pos3 <- match(round(latlist,digits=8),round(unique(echam$lat),digits=8)) 
      zmean <- zmean[pos3,]
      sj <- apply(zmean[2:(nrow(zmean)-1),],2,max) # max slp
      dlat <-  mean(latlist[1:(length(latlist)-1)]-latlist[2:(length(latlist))])    
      gridminpos1 <- which(round(zmean[,1],digits=8)==round(sj[1],digits=8)) 
      gridmin1 <- latlist[gridminpos1]
      gridminpos2 <- which(round(zmean[,2],digits=8)==round(sj[2],digits=8)) 
      gridmin2 <- latlist[gridminpos2]
      zmean <- zmean * -1 # because next equation is to find minimum and here we search for max
      sj <- sj * -1
      latmin1 <- gridmin1 - dlat * (((zmean[(gridminpos1-1),1])-(zmean[(gridminpos1+1),1]))/
        (2*(abs(sj[1]-max(c(zmean[(gridminpos1-1),1],zmean[(gridminpos1+1),1]))))))
      latmin2 <- gridmin2 - dlat * (((zmean[(gridminpos2-1),2])-(zmean[(gridminpos2+1),2]))/
        (2*(abs(sj[2]-max(c(zmean[(gridminpos2-1),2],zmean[(gridminpos2+1),2]))))))
      if (length(latmin1) > 1) {latmin1 <- mean(latmin1)}
      if (length(latmin2) > 1) {latmin2 <- mean(latmin2)}
      sj_slp <- c(latmin1,latmin2)
      # cdo -s fldmax -sellonlatbox,0,360,0,50 -zonmean -sellevel,20000 -selvar,u $f ${STOR_DIR}/SJ.$f
      
      
      # stratospheric polar vortex
      pv_reg1 <- eind_tmp$ensmean[which(eind_tmp$names=="PV1.calc"),]
      pv_reg2 <- eind_tmp$ensmean[which(eind_tmp$names=="PV2.calc"),]
      pv <- pv_reg1 - pv_reg2 
      #cdo -s -sellevel,10000 -selvar,geopoth $f ${STOR_DIR}/gph100.$f
      #cdo -s sub -fldmean -sellonlatbox,0,360,75,90 ${STOR_DIR}/gph100.$f -fldmean -sellonlatbox,0,360,40,55 ${STOR_DIR}/gph100.$f ${STOR_DIR}/z100.$f
      
      # Pacific walker circulation
      #Hind[10,which(echam$names=="omega500")] <- H.giorgi[which(giorgi.short == 'PWC1'),
      #                                                    which(echam$names=="omega500")] 
      pwc_reg1 <- eind_tmp$ensmean[which(eind_tmp$names=="PWC1.calc"),]
      pwc_reg2 <- eind_tmp$ensmean[which(eind_tmp$names=="PWC2.calc"),]
      pwc <- pwc_reg1 - pwc_reg2 
      # cdo -s sellevel,50000 -selvar,omega $f ${STOR_DIR}/omega.$f
      # cdo -s sub -fldmean -sellonlatbox,180,260,-10,10 ${STOR_DIR}/omega.$f -fldmean -sellonlatbox,100,150,-10,10 ${STOR_DIR}/omega.$f ${STOR_DIR}/PWC.$f
      
      # Dynamic Indian Monsoon index
      dimi_reg1 <- eind_tmp$ensmean[which(eind_tmp$names=="DIMI1.calc"),]
      dimi_reg2 <- eind_tmp$ensmean[which(eind_tmp$names=="DIMI2.calc"),]
      dimi <- dimi_reg1 - dimi_reg2 
      #cdo -s sellevel,85000 -selvar,u $f ${STOR_DIR}/u850.$f
      #cdo -s sub -fldmean -sellonlatbox,40,80,5,15 ${STOR_DIR}/u850.$f -fldmean -sellonlatbox,70,90,20,30 ${STOR_DIR}/u850.$f ${STOR_DIR}/DIMI.$f
      
      
      # NAO based on anomalies
      nao_reg1 <- eind_tmp$ensmean[which(eind_tmp$names=="NAO1.calc"),]
      nao_reg2 <- eind_tmp$ensmean[which(eind_tmp$names=="NAO2.calc"),]
      nao <- nao_reg1 - nao_reg2 
      
      # PNA based on anomalies
      pna_reg1 <- eind_tmp$ensmean[which(eind_tmp$names=="PNA1.calc"),]
      pna_reg2 <- eind_tmp$ensmean[which(eind_tmp$names=="PNA2.calc"),]
      pna_reg3 <- eind_tmp$ensmean[which(eind_tmp$names=="PNA3.calc"),]
      pna_reg4 <- eind_tmp$ensmean[which(eind_tmp$names=="PNA4.calc"),]
      pna <- 0.25*(pna_reg1 - pna_reg2 + pna_reg3 -pna_reg4) 
      # PNA = 0.25*(Z[20◦ N,160◦ W] – Z[45◦ N,165◦ W] + Z[55◦ N,115◦ W] – Z[30◦ N,85◦ W]) 
      
      eind <- list(ensmean=matrix(0,nrow=length(indices),ncol=2), names=indices)  
      eind$ensmean[which(indices=='ENH.temp2'),] <- eind_tmp$ensmean[which(eind$names=='ENH.temp2'),]
      eind$ensmean[which(indices=='NAM.temp2'),] <- eind_tmp$ensmean[which(eind$names=='NAM.temp2'),]
      eind$ensmean[which(indices=='SAM.temp2'),] <- eind_tmp$ensmean[which(eind$names=='SAM.temp2'),]
      eind$ensmean[which(indices=='AFR.temp2'),] <- eind_tmp$ensmean[which(eind$names=='AFR.temp2'),]
      eind$ensmean[which(indices=='ASI.temp2'),] <- eind_tmp$ensmean[which(eind$names=='ASI.temp2'),]
      eind$ensmean[which(indices=='AUS.temp2'),] <- eind_tmp$ensmean[which(eind$names=='AUS.temp2'),]
      eind$ensmean[which(indices=='ARC.temp2'),] <- eind_tmp$ensmean[which(eind$names=='ARC.temp2'),]
      eind$ensmean[which(indices=='ANT.temp2'),] <- eind_tmp$ensmean[which(eind$names=='ANT.temp2'),]
      eind$ensmean[which(indices=='NEU.temp2'),] <- eind_tmp$ensmean[which(eind$names=='NEU.temp2'),]
      eind$ensmean[which(indices=='MED.temp2'),] <- eind_tmp$ensmean[which(eind$names=='MED.temp2'),]
      eind$ensmean[which(indices=='GLO.temp2'),] <- eind_tmp$ensmean[which(eind$names=='GLO.temp2'),]
      eind$ensmean[which(indices=='NAM.precip'),] <- eind_tmp$ensmean[which(eind$names=='NAM.precip'),]
      eind$ensmean[which(indices=='SAM.precip'),] <- eind_tmp$ensmean[which(eind$names=='SAM.precip'),]
      eind$ensmean[which(indices=='AFR.precip'),] <- eind_tmp$ensmean[which(eind$names=='AFR.precip'),]
      eind$ensmean[which(indices=='ASI.precip'),] <- eind_tmp$ensmean[which(eind$names=='ASI.precip'),]
      eind$ensmean[which(indices=='AUS.precip'),] <- eind_tmp$ensmean[which(eind$names=='AUS.precip'),]
      eind$ensmean[which(indices=='ARC.precip'),] <- eind_tmp$ensmean[which(eind$names=='ARC.precip'),]
      eind$ensmean[which(indices=='ANT.precip'),] <- eind_tmp$ensmean[which(eind$names=='ANT.precip'),]
      eind$ensmean[which(indices=='NEU.precip'),] <- eind_tmp$ensmean[which(eind$names=='NEU.precip'),]
      eind$ensmean[which(indices=='MED.precip'),] <- eind_tmp$ensmean[which(eind$names=='MED.precip'),]
      eind$ensmean[which(indices=='SH.temp2'),] <- eind_tmp$ensmean[which(eind$names=='SH.temp2'),]
      eind$ensmean[which(indices=='NAM.slp'),] <- eind_tmp$ensmean[which(eind$names=='NAM.slp'),]
      eind$ensmean[which(indices=='SAM.slp'),] <- eind_tmp$ensmean[which(eind$names=='SAM.slp'),]
      eind$ensmean[which(indices=='AFR.slp'),] <- eind_tmp$ensmean[which(eind$names=='AFR.slp'),]
      eind$ensmean[which(indices=='ASI.slp'),] <- eind_tmp$ensmean[which(eind$names=='ASI.slp'),]
      eind$ensmean[which(indices=='AUS.slp'),] <- eind_tmp$ensmean[which(eind$names=='AUS.slp'),]
      eind$ensmean[which(indices=='ARC.slp'),] <- eind_tmp$ensmean[which(eind$names=='ARC.slp'),]
      eind$ensmean[which(indices=='ANT.slp'),] <- eind_tmp$ensmean[which(eind$names=='ANT.slp'),]
      eind$ensmean[which(indices=='NEU.slp'),] <- eind_tmp$ensmean[which(eind$names=='NEU.slp'),]
      eind$ensmean[which(indices=='MED.slp'),] <- eind_tmp$ensmean[which(eind$names=='MED.slp'),]
      eind$ensmean[which(indices=='NH.temp2'),] <- eind_tmp$ensmean[which(eind$names=='NH.temp2'),]
      eind$ensmean[which(indices=='EU.temp2'),] <- eind_tmp$ensmean[which(eind$names=='EU.temp2'),]
      eind$ensmean[which(indices=='EU.precip'),] <- eind_tmp$ensmean[which(eind$names=='EU.precip'),]
      eind$ensmean[which(indices=='EU.slp'),] <- eind_tmp$ensmean[which(eind$names=='EU.slp'),]
      eind$ensmean[which(indices=='HC.calc'),] <- hc #eind_tmp$ensmean[which(eind$names=='HC.calc'),]
      eind$ensmean[which(indices=='ITCZ.calc'),] <- itcz #eind_tmp$ensmean[which(eind$names=='ITCZ.calc'),]
      eind$ensmean[which(indices=='SJ_u200.calc'),] <- sj_u200 #eind_tmp$ensmean[which(eind$names=='SJ.calc'),]
      eind$ensmean[which(indices=='SJ_slp.calc'),] <- sj_slp
      eind$ensmean[which(indices=='PV.calc'),] <- pv #eind_tmp$ensmean[which(eind$names=='PV.calc'),]
      eind$ensmean[which(indices=='PWC.calc'),] <- pwc #eind_tmp$ensmean[which(eind$names=='PWC.calc'),]
      eind$ensmean[which(indices=='DIMI.calc'),] <- dimi #eind_tmp$ensmean[which(eind$names=='DIMI.calc'),]
      eind$ensmean[which(indices=='NAO.calc'),] <- nao #eind_tmp$ensmean[which(eind$names=='DIMI.calc'),] #c(0,0)
      eind$ensmean[which(indices=='PNA.calc'),] <- pna #eind_tmp$ensmean[which(eind$names=='DIMI.calc'),] #c(0,0)
#      eind$ensmean[which(indices=='DIMI'),] <- ech_ind$ensmean[which(ech_ind$names=='DIMI'),]
#      eind$ensmean[which(indices=='HC'),] <- ech_ind$ensmean[which(ech_ind$names=='HC'),]
#      eind$ensmean[which(indices=='SJ'),] <- ech_ind$ensmean[which(ech_ind$names=='SJ'),]
#      eind$ensmean[which(indices=='Z100'),] <- ech_ind$ensmean[which(ech_ind$names=='z100'),]
#      eind$ensmean[which(indices=='Z300'),] <- ech_ind$ensmean[which(ech_ind$names=='z300'),]
#      eind$ensmean[which(indices=='PWC'),] <- ech_ind$ensmean[which(ech_ind$names=='PWC'),]
      
      
      
      # BEM
      eind_tmp <- list(bem=Hind%*%echam$bem, names=indices_tmp)
      
      # Hadley cell strength
      #    Hind[5,which(echam$names=="stream")] <- H.giorgi[which(giorgi.short == 'HC'),
      #                                              which(echam$names=="stream")] 
      # Hind[5,which(echam$names=="omega500")] <- H.giorgi[which(giorgi.short == 'HC'),
      #                                                   which(echam$names=="omega500")] 
      c=1
      zmean <- matrix(0,nrow=length(unique(echam$lat)),ncol=2)
      for (l in unique(echam$lat)) {
        pos1 <- which(trunc(echam$lat,digits=8)==trunc(l,digits=8))
        pos2 <- which(echam$names=='omega500')
        pos <- intersect(pos1,pos2)
        if (length(echam$bem[pos,])==2) { zmean[c,] <- echam$bem[pos,] }
        if (length(echam$bem[pos,])>2) { zmean[c,] <- apply(echam$bem[pos,],2,mean) }
        c=c+1
      }
      latlist <- unique(echam$lat)[which(which(unique(echam$lat) > -30) %in% which(unique(echam$lat) < 30))]
      pos3 <- match(round(latlist,digits=8),round(unique(echam$lat),digits=8)) 
      zmean <- zmean[pos3,]
      hc <- apply(zmean,2,min) # max upward wind is neg. omega value -> min function
      # cdo -s fldmax -sellonlatbox,0,360,0,30 -zonmean -sellevel,50000 -selvar,stream $f ${STOR_DIR}/HC.$f
      
      # Hadley cell poleward extend
      # Hind[6,which(echam$names=="u200")] <- H.giorgi[which(giorgi.short == 'HCL'),
      #                                               which(echam$names=="u200")] 
      
      # ITCZ location
      #Hind[7,which(echam$names=="omega500")] <- H.giorgi[which(giorgi.short == 'ITCZ'),
      #                                                   which(echam$names=="omega500")] 
      c=1
      zmean <- matrix(0,nrow=length(unique(echam$lat)),ncol=2)
      for (l in unique(echam$lat)) {
        pos1 <- which(trunc(echam$lat,digits=8)==trunc(l,digits=8))
        pos2 <- which(echam$names=='omega500')
        pos <- intersect(pos1,pos2)
        if (length(echam$bem[pos,])==2) { zmean[c,] <- echam$bem[pos,] }
        if (length(echam$bem[pos,])>2) { zmean[c,] <- apply(echam$bem[pos,],2,mean) }
        c=c+1
      }
      latlist <- unique(echam$lat)[which(which(unique(echam$lat) > -30) %in% which(unique(echam$lat) < 30))]
      pos3 <- match(round(latlist,digits=8),round(unique(echam$lat),digits=8)) 
      zmean <- zmean[pos3,]
      itcz <- apply(zmean[2:(nrow(zmean)-1),],2,min) 
      dlat <-  mean(latlist[1:(length(latlist)-1)]-latlist[2:(length(latlist))])    
      gridminpos1 <- which(round(zmean[,1],digits=8)==round(itcz[1],digits=8)) 
      gridmin1 <- latlist[gridminpos1]
      gridminpos2 <- which(round(zmean[,2],digits=8)==round(itcz[2],digits=8)) 
      gridmin2 <- latlist[gridminpos2]
      latmin1 <- gridmin1 - dlat * (((zmean[(gridminpos1-1),1])-(zmean[(gridminpos1+1),1]))/
                                      (2*(abs(itcz[1]-max(c(zmean[(gridminpos1-1),1],zmean[(gridminpos1+1),1]))))))
      latmin2 <- gridmin2 - dlat * (((zmean[(gridminpos2-1),2])-(zmean[(gridminpos2+1),2]))/
                                      (2*(abs(itcz[2]-max(c(zmean[(gridminpos2-1),2],zmean[(gridminpos2+1),2]))))))
      if (length(latmin1) > 1) {latmin1 <- mean(latmin1)}
      if (length(latmin2) > 1) {latmin2 <- mean(latmin2)}
      itcz <- c(latmin1,latmin2)
      
      # subtropical jet (u200)
      #Hind[8,which(echam$names=="u200")] <- H.giorgi[which(giorgi.short == 'SJ'),
      #                                               which(echam$names=="u200")] 
      c=1
      zmean <- matrix(0,nrow=length(unique(echam$lat)),ncol=2)
      for (l in unique(echam$lat)) {
        pos1 <- which(trunc(echam$lat,digits=8)==trunc(l,digits=8))
        pos2 <- which(echam$names=='u200')
        pos <- intersect(pos1,pos2)
        if (length(echam$bem[pos,])==2) { zmean[c,] <- echam$bem[pos,] }
        if (length(echam$bem[pos,])>2) { zmean[c,] <- apply(echam$bem[pos,],2,mean) }
        c=c+1
      }
      latlist <- unique(echam$lat)[which(which(unique(echam$lat) > 0) %in% which(unique(echam$lat) < 57))]
      pos3 <- match(round(latlist,digits=8),round(unique(echam$lat),digits=8)) 
      zmean <- zmean[pos3,]
      sj <- apply(zmean[2:(nrow(zmean)-1),],2,max) # max uwind
      dlat <-  mean(latlist[1:(length(latlist)-1)]-latlist[2:(length(latlist))])    
      gridminpos1 <- which(round(zmean[,1],digits=8)==round(sj[1],digits=8)) 
      gridmin1 <- latlist[gridminpos1]
      gridminpos2 <- which(round(zmean[,2],digits=8)==round(sj[2],digits=8)) 
      gridmin2 <- latlist[gridminpos2]
      zmean <- zmean * -1 # because next equation is to find minimum and here we search for max
      sj <- sj * -1
      latmin1 <- gridmin1 - dlat * (((zmean[(gridminpos1-1),1])-(zmean[(gridminpos1+1),1]))/
                                      (2*(abs(sj[1]-max(c(zmean[(gridminpos1-1),1],zmean[(gridminpos1+1),1]))))))
      latmin2 <- gridmin2 - dlat * (((zmean[(gridminpos2-1),2])-(zmean[(gridminpos2+1),2]))/
                                      (2*(abs(sj[2]-max(c(zmean[(gridminpos2-1),2],zmean[(gridminpos2+1),2]))))))
      if (length(latmin1) > 1) {latmin1 <- mean(latmin1)}
      if (length(latmin2) > 1) {latmin2 <- mean(latmin2)}
      sj_u200 <- c(latmin1,latmin2)
      # cdo -s fldmax -sellonlatbox,0,360,0,50 -zonmean -sellevel,20000 -selvar,u $f ${STOR_DIR}/SJ.$f
      
      # subtropical jet (slp)
      c=1
      zmean <- matrix(0,nrow=length(unique(echam$lat)),ncol=2)
      for (l in unique(echam$lat)) {
        pos1 <- which(trunc(echam$lat,digits=8)==trunc(l,digits=8))
        pos2 <- which(echam$names=='slp')
        pos <- intersect(pos1,pos2)
        if (length(echam$bem[pos,])==2) { zmean[c,] <- echam$bem[pos,] }
        if (length(echam$bem[pos,])>2) { zmean[c,] <- apply(echam$bem[pos,],2,mean) }
        c=c+1
      }
      latlist <- unique(echam$lat)[which(which(unique(echam$lat) > 0) %in% which(unique(echam$lat) < 57))]
      pos3 <- match(round(latlist,digits=8),round(unique(echam$lat),digits=8)) 
      zmean <- zmean[pos3,]
      sj <- apply(zmean[2:(nrow(zmean)-1),],2,max) # max uwind
      dlat <-  mean(latlist[1:(length(latlist)-1)]-latlist[2:(length(latlist))])    
      gridminpos1 <- which(round(zmean[,1],digits=8)==round(sj[1],digits=8)) 
      gridmin1 <- latlist[gridminpos1]
      gridminpos2 <- which(round(zmean[,2],digits=8)==round(sj[2],digits=8)) 
      gridmin2 <- latlist[gridminpos2]
      zmean <- zmean * -1 # because next equation is to find minimum and here we search for max
      sj <- sj * -1
      latmin1 <- gridmin1 - dlat * (((zmean[(gridminpos1-1),1])-(zmean[(gridminpos1+1),1]))/
                                      (2*(abs(sj[1]-max(c(zmean[(gridminpos1-1),1],zmean[(gridminpos1+1),1]))))))
      latmin2 <- gridmin2 - dlat * (((zmean[(gridminpos2-1),2])-(zmean[(gridminpos2+1),2]))/
                                      (2*(abs(sj[2]-max(c(zmean[(gridminpos2-1),2],zmean[(gridminpos2+1),2]))))))
      if (length(latmin1) > 1) {latmin1 <- mean(latmin1)}
      if (length(latmin2) > 1) {latmin2 <- mean(latmin2)}
      sj_slp <- c(latmin1,latmin2)
      
      # stratospheric polar vortex
      pv_reg1 <- eind_tmp$bem[which(eind_tmp$names=="PV1.calc"),]
      pv_reg2 <- eind_tmp$bem[which(eind_tmp$names=="PV2.calc"),]
      pv <- pv_reg1 - pv_reg2 
      #cdo -s -sellevel,10000 -selvar,geopoth $f ${STOR_DIR}/gph100.$f
      #cdo -s sub -fldmean -sellonlatbox,0,360,75,90 ${STOR_DIR}/gph100.$f -fldmean -sellonlatbox,0,360,40,55 ${STOR_DIR}/gph100.$f ${STOR_DIR}/z100.$f
      
      # Pacific walker circulation
      #Hind[10,which(echam$names=="omega500")] <- H.giorgi[which(giorgi.short == 'PWC1'),
      #                                                    which(echam$names=="omega500")] 
      pwc_reg1 <- eind_tmp$bem[which(eind_tmp$names=="PWC1.calc"),]
      pwc_reg2 <- eind_tmp$bem[which(eind_tmp$names=="PWC2.calc"),]
      pwc <- pwc_reg1 - pwc_reg2 
      # cdo -s sellevel,50000 -selvar,omega $f ${STOR_DIR}/omega.$f
      # cdo -s sub -fldmean -sellonlatbox,180,260,-10,10 ${STOR_DIR}/omega.$f -fldmean -sellonlatbox,100,150,-10,10 ${STOR_DIR}/omega.$f ${STOR_DIR}/PWC.$f
      
      # Dynamic Indian Monsoon index
      dimi_reg1 <- eind_tmp$bem[which(eind_tmp$names=="DIMI1.calc"),]
      dimi_reg2 <- eind_tmp$bem[which(eind_tmp$names=="DIMI2.calc"),]
      dimi <- dimi_reg1 - dimi_reg2 
      #cdo -s sellevel,85000 -selvar,u $f ${STOR_DIR}/u850.$f
      #cdo -s sub -fldmean -sellonlatbox,40,80,5,15 ${STOR_DIR}/u850.$f -fldmean -sellonlatbox,70,90,20,30 ${STOR_DIR}/u850.$f ${STOR_DIR}/DIMI.$f
      
      
      # NAO based on anomalies
      nao_reg1 <- eind_tmp$bem[which(eind_tmp$names=="NAO1.calc"),]
      nao_reg2 <- eind_tmp$bem[which(eind_tmp$names=="NAO2.calc"),]
      nao <- nao_reg1 - nao_reg2 
      
      # PNA based on anomalies
      pna_reg1 <- eind_tmp$bem[which(eind_tmp$names=="PNA1.calc"),]
      pna_reg2 <- eind_tmp$bem[which(eind_tmp$names=="PNA2.calc"),]
      pna_reg3 <- eind_tmp$bem[which(eind_tmp$names=="PNA3.calc"),]
      pna_reg4 <- eind_tmp$bem[which(eind_tmp$names=="PNA4.calc"),]
      pna <- 0.25*(pna_reg1 - pna_reg2 + pna_reg3 -pna_reg4) 
      # PNA = 0.25*(Z[20◦ N,160◦ W] – Z[45◦ N,165◦ W] + Z[55◦ N,115◦ W] – Z[30◦ N,85◦ W]) 
      
      eind$bem <- matrix(0,nrow=length(indices),ncol=2)
      #eind <- list(bem=matrix(0,nrow=length(indices),ncol=2), names=indices)  
      eind$bem[which(indices=='ENH.temp2'),] <- eind_tmp$bem[which(eind$names=='ENH.temp2'),]
      eind$bem[which(indices=='NAM.temp2'),] <- eind_tmp$bem[which(eind$names=='NAM.temp2'),]
      eind$bem[which(indices=='SAM.temp2'),] <- eind_tmp$bem[which(eind$names=='SAM.temp2'),]
      eind$bem[which(indices=='AFR.temp2'),] <- eind_tmp$bem[which(eind$names=='AFR.temp2'),]
      eind$bem[which(indices=='ASI.temp2'),] <- eind_tmp$bem[which(eind$names=='ASI.temp2'),]
      eind$bem[which(indices=='AUS.temp2'),] <- eind_tmp$bem[which(eind$names=='AUS.temp2'),]
      eind$bem[which(indices=='ARC.temp2'),] <- eind_tmp$bem[which(eind$names=='ARC.temp2'),]
      eind$bem[which(indices=='ANT.temp2'),] <- eind_tmp$bem[which(eind$names=='ANT.temp2'),]
      eind$bem[which(indices=='NEU.temp2'),] <- eind_tmp$bem[which(eind$names=='NEU.temp2'),]
      eind$bem[which(indices=='MED.temp2'),] <- eind_tmp$bem[which(eind$names=='MED.temp2'),]
      eind$bem[which(indices=='GLO.temp2'),] <- eind_tmp$bem[which(eind$names=='GLO.temp2'),]
      eind$bem[which(indices=='NAM.precip'),] <- eind_tmp$bem[which(eind$names=='NAM.precip'),]
      eind$bem[which(indices=='SAM.precip'),] <- eind_tmp$bem[which(eind$names=='SAM.precip'),]
      eind$bem[which(indices=='AFR.precip'),] <- eind_tmp$bem[which(eind$names=='AFR.precip'),]
      eind$bem[which(indices=='ASI.precip'),] <- eind_tmp$bem[which(eind$names=='ASI.precip'),]
      eind$bem[which(indices=='AUS.precip'),] <- eind_tmp$bem[which(eind$names=='AUS.precip'),]
      eind$bem[which(indices=='ARC.precip'),] <- eind_tmp$bem[which(eind$names=='ARC.precip'),]
      eind$bem[which(indices=='ANT.precip'),] <- eind_tmp$bem[which(eind$names=='ANT.precip'),]
      eind$bem[which(indices=='NEU.precip'),] <- eind_tmp$bem[which(eind$names=='NEU.precip'),]
      eind$bem[which(indices=='MED.precip'),] <- eind_tmp$bem[which(eind$names=='MED.precip'),]
      eind$bem[which(indices=='SH.temp2'),] <- eind_tmp$bem[which(eind$names=='SH.temp2'),]
      eind$bem[which(indices=='NAM.slp'),] <- eind_tmp$bem[which(eind$names=='NAM.slp'),]
      eind$bem[which(indices=='SAM.slp'),] <- eind_tmp$bem[which(eind$names=='SAM.slp'),]
      eind$bem[which(indices=='AFR.slp'),] <- eind_tmp$bem[which(eind$names=='AFR.slp'),]
      eind$bem[which(indices=='ASI.slp'),] <- eind_tmp$bem[which(eind$names=='ASI.slp'),]
      eind$bem[which(indices=='AUS.slp'),] <- eind_tmp$bem[which(eind$names=='AUS.slp'),]
      eind$bem[which(indices=='ARC.slp'),] <- eind_tmp$bem[which(eind$names=='ARC.slp'),]
      eind$bem[which(indices=='ANT.slp'),] <- eind_tmp$bem[which(eind$names=='ANT.slp'),]
      eind$bem[which(indices=='NEU.slp'),] <- eind_tmp$bem[which(eind$names=='NEU.slp'),]
      eind$bem[which(indices=='MED.slp'),] <- eind_tmp$bem[which(eind$names=='MED.slp'),]
      eind$bem[which(indices=='NH.temp2'),] <- eind_tmp$bem[which(eind$names=='NH.temp2'),]
      eind$bem[which(indices=='EU.temp2'),] <- eind_tmp$bem[which(eind$names=='EU.temp2'),]
      eind$bem[which(indices=='EU.precip'),] <- eind_tmp$bem[which(eind$names=='EU.precip'),]
      eind$bem[which(indices=='EU.slp'),] <- eind_tmp$bem[which(eind$names=='EU.slp'),]
      eind$bem[which(indices=='HC.calc'),] <- hc #eind_tmp$bem[which(eind$names=='HC.calc'),]
      eind$bem[which(indices=='ITCZ.calc'),] <- itcz #eind_tmp$bem[which(eind$names=='ITCZ.calc'),]
      eind$bem[which(indices=='SJ_u200.calc'),] <- sj_u200 #eind_tmp$bem[which(eind$names=='SJ.calc'),]
      eind$bem[which(indices=='SJ_slp.calc'),] <- sj_slp
      eind$bem[which(indices=='PV.calc'),] <- pv #eind_tmp$bem[which(eind$names=='PV.calc'),]
      eind$bem[which(indices=='PWC.calc'),] <- pwc #eind_tmp$bem[which(eind$names=='PWC.calc'),]
      eind$bem[which(indices=='DIMI.calc'),] <- dimi #eind_tmp$bem[which(eind$names=='DIMI.calc'),]
      eind$bem[which(indices=='NAO.calc'),] <- nao #eind_tmp$bem[which(eind$names=='DIMI.calc'),] #c(0,0)
      eind$bem[which(indices=='PNA.calc'),] <- pna #eind_tmp$bem[which(eind$names=='DIMI.calc'),] #c(0,0)
#      eind$bem[which(indices=='DIMI'),] <- ech_ind$bem[which(ech_ind$names=='DIMI'),]
#      eind$bem[which(indices=='HC'),] <- ech_ind$bem[which(ech_ind$names=='HC'),]
#      eind$bem[which(indices=='SJ'),] <- ech_ind$bem[which(ech_ind$names=='SJ'),]
#      eind$bem[which(indices=='Z100'),] <- ech_ind$bem[which(ech_ind$names=='z100'),]
#      eind$bem[which(indices=='Z300'),] <- ech_ind$bem[which(ech_ind$names=='z300'),]
#      eind$bem[which(indices=='PWC'),] <- ech_ind$bem[which(ech_ind$names=='PWC'),]
      
      
      
      
      
      
      
      # data for ensemble member
      #    eind$data <- array(Hind %*% array(echam$data, c(nrow(echam$data), 
      #                       length(echam$data)/nrow(echam$data))), c(nrow(Hind), 
      #                       dim(echam$data)[2:3]))
      eind_tmp <- list(data=array(Hind %*% array(echam$data, c(nrow(echam$data), 
                                                                   length(echam$data)/nrow(echam$data))), c(nrow(Hind), 
                                                                                                                    dim(echam$data)[2:3])), names=indices_tmp)
      # Hadley cell strength
      hc <- matrix(0,nrow=nmem,ncol=2)
      for (m in 1:nmem) {
        c=1
        zmean <- matrix(0,nrow=length(unique(echam$lat)),ncol=2)
        for (l in unique(echam$lat)) {
          pos1 <- which(trunc(echam$lat,digits=8)==trunc(l,digits=8))
          pos2 <- which(echam$names=='omega500')
          pos <- intersect(pos1,pos2)
          if (length(echam$data[pos,,m])==2) { zmean[c,] <- echam$data[pos,,m] }
          if (length(echam$data[pos,,m])>2) { zmean[c,] <- apply(echam$data[pos,,m],2,mean) }
          #        zmean[c,] <- apply(echam$data[pos,,m],2,mean)
          c=c+1
        }
        latlist <- unique(echam$lat)[which(which(unique(echam$lat) > -30) %in% 
                                                 which(unique(echam$lat) < 30))]
        pos3 <- match(round(latlist,digits=8),round(unique(echam$lat),digits=8)) 
        zmean <- zmean[pos3,]
        hc[m,] <- apply(zmean,2,min) # max upward wind is neg. omega value -> min function
      }
      # Hadley cell poleward extend
      
      # ITCZ location
      itcz <- matrix(0,nrow=nmem,ncol=2)
      for (m in 1:nmem) {
        c=1
        zmean <- matrix(0,nrow=length(unique(echam$lat)),ncol=2)
        for (l in unique(echam$lat)) {
          pos1 <- which(trunc(echam$lat,digits=8)==trunc(l,digits=8))
          pos2 <- which(echam$names=='omega500')
          pos <- intersect(pos1,pos2)
          if (length(echam$data[pos,,m])==2) { zmean[c,] <- echam$data[pos,,m] }
          if (length(echam$data[pos,,m])>2) { zmean[c,] <- apply(echam$data[pos,,m],2,mean) }
          c=c+1
        }
        latlist <- unique(echam$lat)[which(which(unique(echam$lat) > -30) %in% 
                                                 which(unique(echam$lat) < 30))]
        pos3 <- match(round(latlist,digits=8),round(unique(echam$lat),digits=8)) 
        zmean <- zmean[pos3,]
        itcz_tmp <- apply(zmean[2:(nrow(zmean)-1),],2,min) 
        dlat <-  mean(latlist[1:(length(latlist)-1)]-latlist[2:(length(latlist))])    
        gridminpos1 <- which(round(zmean[,1],digits=8)==round(itcz_tmp[1],digits=8)) 
        gridmin1 <- latlist[gridminpos1]
        gridminpos2 <- which(round(zmean[,2],digits=8)==round(itcz_tmp[2],digits=8)) 
        gridmin2 <- latlist[gridminpos2]
        latmin1 <- gridmin1 - dlat * (((zmean[(gridminpos1-1),1])-(zmean[(gridminpos1+1),1]))/
          (2*(abs(itcz_tmp[1]-max(c(zmean[(gridminpos1-1),1],zmean[(gridminpos1+1),1]))))))
        latmin2 <- gridmin2 - dlat * (((zmean[(gridminpos2-1),2])-(zmean[(gridminpos2+1),2]))/
          (2*(abs(itcz_tmp[2]-max(c(zmean[(gridminpos2-1),2],zmean[(gridminpos2+1),2]))))))
        if (length(latmin1) > 1) {latmin1 <- mean(latmin1)}
        if (length(latmin2) > 1) {latmin2 <- mean(latmin2)}
        itcz[m,] <- c(latmin1,latmin2)
      }
      
      # subtropical jet (u200)
      sj_u200 <- matrix(0,nrow=nmem,ncol=2)
      for (m in 1:nmem) {
        c=1
        zmean <- matrix(0,nrow=length(unique(echam$lat)),ncol=2)
        for (l in unique(echam$lat)) {
          pos1 <- which(trunc(echam$lat,digits=8)==trunc(l,digits=8))
          pos2 <- which(echam$names=='u200')
          pos <- intersect(pos1,pos2)
          if (length(echam$data[pos,,m])==2) { zmean[c,] <- echam$data[pos,,m] }
          if (length(echam$data[pos,,m])>2) { zmean[c,] <- apply(echam$data[pos,,m],2,mean) }
          c=c+1
        }
        latlist <- unique(echam$lat)[which(which(unique(echam$lat) > 0) %in% 
                                                 which(unique(echam$lat) < 57))]
        pos3 <- match(round(latlist,digits=8),round(unique(echam$lat),digits=8)) 
        zmean <- zmean[pos3,]
        sj_tmp <- apply(zmean[2:(nrow(zmean)-1),],2,max) # max uwind
        dlat <-  mean(latlist[1:(length(latlist)-1)]-latlist[2:(length(latlist))])    
        gridminpos1 <- which(round(zmean[,1],digits=8)==round(sj_tmp[1],digits=8)) 
        gridmin1 <- latlist[gridminpos1]
        gridminpos2 <- which(round(zmean[,2],digits=8)==round(sj_tmp[2],digits=8)) 
        gridmin2 <- latlist[gridminpos2]
        zmean <- zmean * -1 # because next equation is to find minimum and here we search for max
        sj_tmp <- sj_tmp * -1
        latmin1 <- gridmin1 - dlat * (((zmean[(gridminpos1-1),1])-(zmean[(gridminpos1+1),1]))/
                                        (2*(abs(sj_tmp[1]-max(c(zmean[(gridminpos1-1),1],zmean[(gridminpos1+1),1]))))))
        latmin2 <- gridmin2 - dlat * (((zmean[(gridminpos2-1),2])-(zmean[(gridminpos2+1),2]))/
                                        (2*(abs(sj_tmp[2]-max(c(zmean[(gridminpos2-1),2],zmean[(gridminpos2+1),2]))))))
        if (length(latmin1) > 1) {latmin1 <- mean(latmin1)}
        if (length(latmin2) > 1) {latmin2 <- mean(latmin2)}
        sj_u200[m,] <- c(latmin1,latmin2)
      }
      
      # subtropical jet (slp)
      sj_slp <- matrix(0,nrow=nmem,ncol=2)
      for (m in 1:nmem) {
        c=1
        zmean <- matrix(0,nrow=length(unique(echam$lat)),ncol=2)
        for (l in unique(echam$lat)) {
          pos1 <- which(trunc(echam$lat,digits=8)==trunc(l,digits=8))
          pos2 <- which(echam$names=='slp')
          pos <- intersect(pos1,pos2)
          if (length(echam$data[pos,,m])==2) { zmean[c,] <- echam$data[pos,,m] }
          if (length(echam$data[pos,,m])>2) { zmean[c,] <- apply(echam$data[pos,,m],2,mean) }
          c=c+1
        }
        latlist <- unique(echam$lat)[which(which(unique(echam$lat) > 0) %in% 
                                                 which(unique(echam$lat) < 57))]
        pos3 <- match(round(latlist,digits=8),round(unique(echam$lat),digits=8)) 
        zmean <- zmean[pos3,]
        sj_tmp <- apply(zmean[2:(nrow(zmean)-1),],2,max) # max slp
        dlat <-  mean(latlist[1:(length(latlist)-1)]-latlist[2:(length(latlist))])    
        gridminpos1 <- which(round(zmean[,1],digits=8)==round(sj_tmp[1],digits=8)) 
        gridmin1 <- latlist[gridminpos1]
        gridminpos2 <- which(round(zmean[,2],digits=8)==round(sj_tmp[2],digits=8)) 
        gridmin2 <- latlist[gridminpos2]
        zmean <- zmean * -1 # because next equation is to find minimum and here we search for max
        sj_tmp <- sj_tmp * -1
        latmin1 <- gridmin1 - dlat * (((zmean[(gridminpos1-1),1])-(zmean[(gridminpos1+1),1]))/
                                        (2*(abs(sj_tmp[1]-max(c(zmean[(gridminpos1-1),1],zmean[(gridminpos1+1),1]))))))
        latmin2 <- gridmin2 - dlat * (((zmean[(gridminpos2-1),2])-(zmean[(gridminpos2+1),2]))/
                                        (2*(abs(sj_tmp[2]-max(c(zmean[(gridminpos2-1),2],zmean[(gridminpos2+1),2]))))))
        if (length(latmin1) > 1) {latmin1 <- mean(latmin1)}
        if (length(latmin2) > 1) {latmin2 <- mean(latmin2)}
        sj_slp[m,] <- c(latmin1,latmin2)
      }    
      
      # stratospheric polar vortex
      pv <- matrix(0,nrow=nmem,ncol=2)
      for (m in 1:nmem) {
        pv_reg1 <- eind_tmp$data[which(eind_tmp$names=="PV1.calc"),,m]
        pv_reg2 <- eind_tmp$data[which(eind_tmp$names=="PV2.calc"),,m]
        pv[m,] <- pv_reg1 - pv_reg2 
      }
      
      # Pacific walker circulation
      pwc <- matrix(0,nrow=nmem,ncol=2)
      for (m in 1:nmem) {
        pwc_reg1 <- eind_tmp$data[which(eind_tmp$names=="PWC1.calc"),,m]
        pwc_reg2 <- eind_tmp$data[which(eind_tmp$names=="PWC2.calc"),,m]
        pwc[m,] <- pwc_reg1 - pwc_reg2 
      }
      
      # Dynamic Indian Monsoon index
      dimi <- matrix(0,nrow=nmem,ncol=2)
      for (m in 1:nmem) {
        dimi_reg1 <- eind_tmp$data[which(eind_tmp$names=="DIMI1.calc"),,m]
        dimi_reg2 <- eind_tmp$data[which(eind_tmp$names=="DIMI2.calc"),,m]
        dimi[m,] <- dimi_reg1 - dimi_reg2 
      }
      
      # NAO based on anomalies
      nao <- matrix(0,nrow=nmem,ncol=2)
      for (m in 1:nmem) {
        nao_reg1 <- eind_tmp$data[which(eind_tmp$names=="NAO1.calc"),,m]
        nao_reg2 <- eind_tmp$data[which(eind_tmp$names=="NAO2.calc"),,m]
        nao[m,] <- nao_reg1 - nao_reg2 
      }
      
      # PNA based on anomalies
      pna <- matrix(0,nrow=nmem,ncol=2)
      for (m in 1:nmem) {
        pna_reg1 <- eind_tmp$data[which(eind_tmp$names=="PNA1.calc"),,m]
        pna_reg2 <- eind_tmp$data[which(eind_tmp$names=="PNA2.calc"),,m]
        pna_reg3 <- eind_tmp$data[which(eind_tmp$names=="PNA3.calc"),,m]
        pna_reg4 <- eind_tmp$data[which(eind_tmp$names=="PNA4.calc"),,m]
        pna[m,] <- 0.25*(pna_reg1 - pna_reg2 + pna_reg3 -pna_reg4) 
      }
      
      eind$data <- array(0,c(length(indices),2,nmem))
      for (m in 1:nmem) {
        eind$data[which(indices=='ENH.temp2'),,m] <- eind_tmp$data[which(eind$names=='ENH.temp2'),,m]
        eind$data[which(indices=='NAM.temp2'),,m] <- eind_tmp$data[which(eind$names=='NAM.temp2'),,m]
        eind$data[which(indices=='SAM.temp2'),,m] <- eind_tmp$data[which(eind$names=='SAM.temp2'),,m]
        eind$data[which(indices=='AFR.temp2'),,m] <- eind_tmp$data[which(eind$names=='AFR.temp2'),,m]
        eind$data[which(indices=='ASI.temp2'),,m] <- eind_tmp$data[which(eind$names=='ASI.temp2'),,m]
        eind$data[which(indices=='AUS.temp2'),,m] <- eind_tmp$data[which(eind$names=='AUS.temp2'),,m]
        eind$data[which(indices=='ARC.temp2'),,m] <- eind_tmp$data[which(eind$names=='ARC.temp2'),,m]
        eind$data[which(indices=='ANT.temp2'),,m] <- eind_tmp$data[which(eind$names=='ANT.temp2'),,m]
        eind$data[which(indices=='NEU.temp2'),,m] <- eind_tmp$data[which(eind$names=='NEU.temp2'),,m]
        eind$data[which(indices=='MED.temp2'),,m] <- eind_tmp$data[which(eind$names=='MED.temp2'),,m]
        eind$data[which(indices=='GLO.temp2'),,m] <- eind_tmp$data[which(eind$names=='GLO.temp2'),,m]
        eind$data[which(indices=='NAM.precip'),,m] <- eind_tmp$data[which(eind$names=='NAM.precip'),,m]
        eind$data[which(indices=='SAM.precip'),,m] <- eind_tmp$data[which(eind$names=='SAM.precip'),,m]
        eind$data[which(indices=='AFR.precip'),,m] <- eind_tmp$data[which(eind$names=='AFR.precip'),,m]
        eind$data[which(indices=='ASI.precip'),,m] <- eind_tmp$data[which(eind$names=='ASI.precip'),,m]
        eind$data[which(indices=='AUS.precip'),,m] <- eind_tmp$data[which(eind$names=='AUS.precip'),,m]
        eind$data[which(indices=='ARC.precip'),,m] <- eind_tmp$data[which(eind$names=='ARC.precip'),,m]
        eind$data[which(indices=='ANT.precip'),,m] <- eind_tmp$data[which(eind$names=='ANT.precip'),,m]
        eind$data[which(indices=='NEU.precip'),,m] <- eind_tmp$data[which(eind$names=='NEU.precip'),,m]
        eind$data[which(indices=='MED.precip'),,m] <- eind_tmp$data[which(eind$names=='MED.precip'),,m]
        eind$data[which(indices=='SH.temp2'),,m] <- eind_tmp$data[which(eind$names=='SH.temp2'),,m]
        eind$data[which(indices=='NAM.slp'),,m] <- eind_tmp$data[which(eind$names=='NAM.slp'),,m]
        eind$data[which(indices=='SAM.slp'),,m] <- eind_tmp$data[which(eind$names=='SAM.slp'),,m]
        eind$data[which(indices=='AFR.slp'),,m] <- eind_tmp$data[which(eind$names=='AFR.slp'),,m]
        eind$data[which(indices=='ASI.slp'),,m] <- eind_tmp$data[which(eind$names=='ASI.slp'),,m]
        eind$data[which(indices=='AUS.slp'),,m] <- eind_tmp$data[which(eind$names=='AUS.slp'),,m]
        eind$data[which(indices=='ARC.slp'),,m] <- eind_tmp$data[which(eind$names=='ARC.slp'),,m]
        eind$data[which(indices=='ANT.slp'),,m] <- eind_tmp$data[which(eind$names=='ANT.slp'),,m]
        eind$data[which(indices=='NEU.slp'),,m] <- eind_tmp$data[which(eind$names=='NEU.slp'),,m]
        eind$data[which(indices=='MED.slp'),,m] <- eind_tmp$data[which(eind$names=='MED.slp'),,m]
        eind$data[which(indices=='NH.temp2'),,m] <- eind_tmp$data[which(eind$names=='NH.temp2'),,m]
        eind$data[which(indices=='EU.temp2'),,m] <- eind_tmp$data[which(eind$names=='EU.temp2'),,m]
        eind$data[which(indices=='EU.precip'),,m] <- eind_tmp$data[which(eind$names=='EU.precip'),,m]
        eind$data[which(indices=='EU.slp'),,m] <- eind_tmp$data[which(eind$names=='EU.slp'),,m]
        eind$data[which(indices=='HC.calc'),,m] <- hc[m,] #eind_tmp$data[which(eind$names=='HC.calc'),,m]
        eind$data[which(indices=='ITCZ.calc'),,m] <- itcz[m,] #eind_tmp$data[which(eind$names=='ITCZ.calc'),,m]
        eind$data[which(indices=='SJ_u200.calc'),,m] <- sj_u200[m,] #eind_tmp$data[which(eind$names=='SJ.calc'),,m]
        eind$data[which(indices=='SJ_slp.calc'),,m] <- sj_slp[m,]
        eind$data[which(indices=='PV.calc'),,m] <- pv[m,] #eind_tmp$data[which(eind$names=='PV.calc'),,m]
        eind$data[which(indices=='PWC.calc'),,m] <- pwc[m,] #eind_tmp$data[which(eind$names=='PWC.calc'),,m]
        eind$data[which(indices=='DIMI.calc'),,m] <- dimi[m,] #eind_tmp$data[which(eind$names=='DIMI.calc'),,m]
        eind$data[which(indices=='NAO.calc'),,m] <- nao[m,] #eind_tmp$data[which(eind$names=='DIMI.calc'),,m] #c(0,0)
        eind$data[which(indices=='PNA.calc'),,m] <- pna[m,] #eind_tmp$data[which(eind$names=='DIMI.calc'),,m] #c(0,0)
#        eind$data[which(indices=='DIMI'),,m] <- ech_ind$data[which(ech_ind$names=='DIMI'),,m]
#        eind$data[which(indices=='HC'),,m] <- ech_ind$data[which(ech_ind$names=='HC'),,m]
#        eind$data[which(indices=='SJ'),,m] <- ech_ind$data[which(ech_ind$names=='SJ'),,m]
#        eind$data[which(indices=='Z100'),,m] <- ech_ind$data[which(ech_ind$names=='z100'),,m]
#        eind$data[which(indices=='Z300'),,m] <- ech_ind$data[which(ech_ind$names=='z300'),,m]
#        eind$data[which(indices=='PWC'),,m] <- ech_ind$data[which(ech_ind$names=='PWC'),,m]
      }
      
 #     eind$data[which(indices=='DIMI'),,] <- ech_ind$data[which(ech_ind$names=='DIMI'),,]
#      eind$data[which(indices=='HC'),,] <- ech_ind$data[which(ech_ind$names=='HC'),,]
#      eind$data[which(indices=='SJ'),,] <- ech_ind$data[which(ech_ind$names=='SJ'),,]
#      eind$data[which(indices=='Z100'),,] <- ech_ind$data[which(ech_ind$names=='z100'),,]
#      eind$data[which(indices=='Z300'),,] <- ech_ind$data[which(ech_ind$names=='z300'),,]
#      eind$data[which(indices=='PWC'),,] <- ech_ind$data[which(ech_ind$names=='PWC'),,]
      
      
      
      
      # SAME FOR ANALYSIS
      H.giorgi <- compute_giorgi_H_v2(giorgi, analysis) #, numvar=3) # 3 vars temp, precip, slp
      #    H.giorgi <- compute_giorgi_H_sixmon(giorgi, analysis) #, numvar=3) # 3 vars temp, precip, slp
      #    H.giorgi.anom <- compute_giorgi_H_sixmon(giorgi, analysis.anom) #, numvar=3) # 3 vars temp, precip, slp
      #    indices <- c('NH.temp2', 'NH.precip', 'NH.slp', 'NEU.temp2', 'NEU.precip', 'NEU.slp',
      #                 'DIMI', 'HC', 'SJ', 'Z100', 'Z300', 'PWC')
      #    Hind <- matrix(0,nrow=12,ncol=nrow(analysis$data))
      indices_tmp <- c('ENH.temp2','NAM.temp2','SAM.temp2','AFR.temp2',
                       'ASI.temp2','AUS.temp2','ARC.temp2','ANT.temp2',
                       'NEU.temp2','MED.temp2',
                       'GLO.temp2','NAM.precip','SAM.precip','AFR.precip',
                       'ASI.precip','AUS.precip','ARC.precip','ANT.precip',
                       'NEU.precip','MED.precip',
                       'SH.temp2','NAM.slp','SAM.slp','AFR.slp',
                       'ASI.slp','AUS.slp','ARC.slp','ANT.slp',
                       'NEU.slp','MED.slp',
                       'NH.temp2', 'EU.temp2', 'EU.precip', 'EU.slp', 
                       'PV1.calc', 'PV2.calc', 'PWC1.calc', 'PWC2.calc', 
                       'DIMI1.calc', 'DIMI2.calc', 'NAO1.calc', 'NAO2.calc',  
                       'PNA1.calc', 'PNA2.calc', 'PNA3.calc', 'PNA4.calc')
      indices <- c('ENH.temp2','NAM.temp2','SAM.temp2','AFR.temp2',
                   'ASI.temp2','AUS.temp2','ARC.temp2','ANT.temp2',
                   'NEU.temp2','MED.temp2',
                   'GLO.temp2','NAM.precip','SAM.precip','AFR.precip',
                   'ASI.precip','AUS.precip','ARC.precip','ANT.precip',
                   'NEU.precip','MED.precip',
                   'SH.temp2','NAM.slp','SAM.slp','AFR.slp',
                   'ASI.slp','AUS.slp','ARC.slp','ANT.slp',
                   'NEU.slp','MED.slp',
                   'NH.temp2', 'EU.temp2', 'EU.precip', 'EU.slp', 
                   'HC.calc', 'ITCZ.calc', 'SJ_u200.calc', 'SJ_slp.calc', 'PV.calc', 
                   'PWC.calc', 'DIMI.calc', 'NAO.calc', 'PNA.calc') #,
 #                  'DIMI', 'HC', 'SJ', 'Z100', 'Z300', 'PWC')
      # "MC" kein z300 in output! 'HCL.calc',
      # SO WIRD INDEX NUR MITTELWERT ¨UBER REGION UND ZEIT! 
      # Deshalb Indicies die min max etc basiert sind extra rechnen
      Hind <- matrix(0,nrow=length(indices_tmp),ncol=nrow(analysis$data))
      Hind[1,which(analysis$names=="temp2")] <- H.giorgi[which(giorgi.short == 'ENH'),
                                                             which(analysis$names=="temp2")]
      Hind[2,which(analysis$names=="temp2")] <- H.giorgi[which(giorgi.short == 'NAM'),
                                                             which(analysis$names=="temp2")]
      Hind[3,which(analysis$names=="temp2")] <- H.giorgi[which(giorgi.short == 'SAM'),
                                                             which(analysis$names=="temp2")]
      Hind[4,which(analysis$names=="temp2")] <- H.giorgi[which(giorgi.short == 'AFR'),
                                                             which(analysis$names=="temp2")]
      Hind[5,which(analysis$names=="temp2")] <- H.giorgi[which(giorgi.short == 'ASI'),
                                                             which(analysis$names=="temp2")]
      Hind[6,which(analysis$names=="temp2")] <- H.giorgi[which(giorgi.short == 'AUS'),
                                                             which(analysis$names=="temp2")]
      Hind[7,which(analysis$names=="temp2")] <- H.giorgi[which(giorgi.short == 'ARC'),
                                                             which(analysis$names=="temp2")]
      Hind[8,which(analysis$names=="temp2")] <- H.giorgi[which(giorgi.short == 'ANT'),
                                                             which(analysis$names=="temp2")]
      Hind[9,which(analysis$names=="temp2")] <- H.giorgi[which(giorgi.short == 'NEU'),
                                                             which(analysis$names=="temp2")]
      Hind[10,which(analysis$names=="temp2")] <- H.giorgi[which(giorgi.short == 'MED'),
                                                              which(analysis$names=="temp2")]
      Hind[11,which(analysis$names=="temp2")] <- H.giorgi[which(giorgi.short == 'GLO'),
                                                              which(analysis$names=="temp2")]
      Hind[12,which(analysis$names=="precip")] <- H.giorgi[which(giorgi.short == 'NAM'),
                                                               which(analysis$names=="precip")]
      Hind[13,which(analysis$names=="precip")] <- H.giorgi[which(giorgi.short == 'SAM'),
                                                               which(analysis$names=="precip")]
      Hind[14,which(analysis$names=="precip")] <- H.giorgi[which(giorgi.short == 'AFR'),
                                                               which(analysis$names=="precip")]
      Hind[15,which(analysis$names=="precip")] <- H.giorgi[which(giorgi.short == 'ASI'),
                                                               which(analysis$names=="precip")]
      Hind[16,which(analysis$names=="precip")] <- H.giorgi[which(giorgi.short == 'AUS'),
                                                               which(analysis$names=="precip")]
      Hind[17,which(analysis$names=="precip")] <- H.giorgi[which(giorgi.short == 'ARC'),
                                                               which(analysis$names=="precip")]
      Hind[18,which(analysis$names=="precip")] <- H.giorgi[which(giorgi.short == 'ANT'),
                                                               which(analysis$names=="precip")]
      Hind[19,which(analysis$names=="precip")] <- H.giorgi[which(giorgi.short == 'NEU'),
                                                               which(analysis$names=="precip")]
      Hind[20,which(analysis$names=="precip")] <- H.giorgi[which(giorgi.short == 'MED'),
                                                               which(analysis$names=="precip")]
      Hind[21,which(analysis$names=="temp2")] <- H.giorgi[which(giorgi.short == 'SH'),
                                                              which(analysis$names=="temp2")]
      Hind[22,which(analysis$names=="slp")] <- H.giorgi[which(giorgi.short == 'NAM'),
                                                            which(analysis$names=="slp")]
      Hind[23,which(analysis$names=="slp")] <- H.giorgi[which(giorgi.short == 'SAM'),
                                                            which(analysis$names=="slp")]
      Hind[24,which(analysis$names=="slp")] <- H.giorgi[which(giorgi.short == 'AFR'),
                                                            which(analysis$names=="slp")]
      Hind[25,which(analysis$names=="slp")] <- H.giorgi[which(giorgi.short == 'ASI'),
                                                            which(analysis$names=="slp")]
      Hind[26,which(analysis$names=="slp")] <- H.giorgi[which(giorgi.short == 'AUS'),
                                                            which(analysis$names=="slp")]
      Hind[27,which(analysis$names=="slp")] <- H.giorgi[which(giorgi.short == 'ARC'),
                                                            which(analysis$names=="slp")]
      Hind[28,which(analysis$names=="slp")] <- H.giorgi[which(giorgi.short == 'ANT'),
                                                            which(analysis$names=="slp")]
      Hind[29,which(analysis$names=="slp")] <- H.giorgi[which(giorgi.short == 'NEU'),
                                                            which(analysis$names=="slp")]
      Hind[30,which(analysis$names=="slp")] <- H.giorgi[which(giorgi.short == 'MED'),
                                                            which(analysis$names=="slp")]  
      
      Hind[31,which(analysis$names=="temp2")] <- H.giorgi[which(giorgi.short == 'NH'),
                                                              which(analysis$names=="temp2")]
      Hind[32,which(analysis$names=="temp2")] <- H.giorgi[which(giorgi.short == 'EU'),
                                                              which(analysis$names=="temp2")]
      Hind[33,which(analysis$names=="precip")] <- H.giorgi[which(giorgi.short == 'EU'),
                                                               which(analysis$names=="precip")]
      Hind[34,which(analysis$names=="slp")] <- H.giorgi[which(giorgi.short == 'EU'),
                                                            which(analysis$names=="slp")] 
      # calc indices from model fields
      # midlatitude circulation
      #    Hind[5,which(analysis$names=="gph300")] <- H.giorgi[which(giorgi.short == 'MC'),
      #                                              which(analysis$names=="gph300")] 
      # cdo -s fldmean -sellonlatbox,0,360,30,60 -sellevel,3000 -selvar,geopoth $f ${STOR_DIR}/z300.$f
      
      # stratospheric polar vortex
      Hind[35,which(analysis$names=="gph100")] <- H.giorgi[which(giorgi.short == 'PV1'),
                                                               which(analysis$names=="gph100")] 
      Hind[36,which(analysis$names=="gph100")] <- H.giorgi[which(giorgi.short == 'PV2'),
                                                               which(analysis$names=="gph100")] 
      # cdo -s -sellevel,10000 -selvar,geopoth $f ${STOR_DIR}/gph100.$f
      # cdo -s sub -fldmean -sellonlatbox,0,360,75,90 ${STOR_DIR}/gph100.$f -fldmean -sellonlatbox,0,360,40,55 ${STOR_DIR}/gph100.$f ${STOR_DIR}/z100.$f
      
      # Pacific walker circulation
      Hind[37,which(analysis$names=="omega500")] <- H.giorgi[which(giorgi.short == 'PWC1'),
                                                                 which(analysis$names=="omega500")] 
      # cdo -s sellevel,50000 -selvar,omega $f ${STOR_DIR}/omega.$f
      # cdo -s sub -fldmean -sellonlatbox,180,260,-10,10 ${STOR_DIR}/omega.$f -fldmean -sellonlatbox,100,150,-10,10 ${STOR_DIR}/omega.$f ${STOR_DIR}/PWC.$f
      
      # Pacific walker circulation
      Hind[38,which(analysis$names=="omega500")] <- H.giorgi[which(giorgi.short == 'PWC2'),
                                                                 which(analysis$names=="omega500")] 
      # cdo -s sellevel,50000 -selvar,omega $f ${STOR_DIR}/omega.$f
      # cdo -s sub -fldmean -sellonlatbox,180,260,-10,10 ${STOR_DIR}/omega.$f -fldmean -sellonlatbox,100,150,-10,10 ${STOR_DIR}/omega.$f ${STOR_DIR}/PWC.$f
      
      # Dynamic Indian Monsoon index
      Hind[39,which(analysis$names=="u850")] <- H.giorgi[which(giorgi.short == 'DIMI1'),
                                                             which(analysis$names=="u850")] 
      # cdo -s sellevel,85000 -selvar,u $f ${STOR_DIR}/u850.$f
      # cdo -s sub -fldmean -sellonlatbox,40,80,5,15 ${STOR_DIR}/u850.$f -fldmean -sellonlatbox,70,90,20,30 ${STOR_DIR}/u850.$f ${STOR_DIR}/DIMI.$f
      
      # Dynamic Indian Monsoon index
      Hind[40,which(analysis$names=="u850")] <- H.giorgi[which(giorgi.short == 'DIMI2'),
                                                             which(analysis$names=="u850")] 
      # cdo -s sellevel,85000 -selvar,u $f ${STOR_DIR}/u850.$f
      # cdo -s sub -fldmean -sellonlatbox,40,80,5,15 ${STOR_DIR}/u850.$f -fldmean -sellonlatbox,70,90,20,30 ${STOR_DIR}/u850.$f ${STOR_DIR}/DIMI.$f
      
      # NAO index
      Hind[41,which(analysis.anom$names=="slp")] <- H.giorgi[which(giorgi.short == 'NAO1'),
                                                             which(analysis.anom$names=="slp")] 
      Hind[42,which(analysis.anom$names=="slp")] <- H.giorgi[which(giorgi.short == 'NAO2'),
                                                             which(analysis.anom$names=="slp")] 
      Hind[43,which(analysis.anom$names=="gph500")] <- H.giorgi[which(giorgi.short == 'PNA1'),
                                                                which(analysis.anom$names=="gph500")] 
      Hind[44,which(analysis.anom$names=="gph500")] <- H.giorgi[which(giorgi.short == 'PNA2'),
                                                                which(analysis.anom$names=="gph500")] 
      Hind[45,which(analysis.anom$names=="gph500")] <- H.giorgi[which(giorgi.short == 'PNA3'),
                                                                which(analysis.anom$names=="gph500")] 
      Hind[46,which(analysis.anom$names=="gph500")] <- H.giorgi[which(giorgi.short == 'PNA4'),
                                                                which(analysis.anom$names=="gph500")] 
      aind_tmp <- list(ensmean=Hind%*%analysis$ensmean, names=indices_tmp)
      
      # Hadley cell strength
      #    Hind[5,which(analysis$names=="stream")] <- H.giorgi[which(giorgi.short == 'HC'),
      #                                              which(analysis$names=="stream")] 
      # Hind[5,which(analysis$names=="omega500")] <- H.giorgi[which(giorgi.short == 'HC'),
      #                                                   which(analysis$names=="omega500")] 
      c=1
      zmean <- matrix(0,nrow=length(unique(analysis$lat)),ncol=2)
      for (l in unique(analysis$lat)) {
        pos1 <- which(trunc(analysis$lat,digits=8)==trunc(l,digits=8))
        pos2 <- which(analysis$names=='omega500')
        pos <- intersect(pos1,pos2)
        if (length(analysis$ensmean[pos,])==2) { zmean[c,] <- analysis$ensmean[pos,] }
        if (length(analysis$ensmean[pos,])>2) { zmean[c,] <- apply(analysis$ensmean[pos,],2,mean) }
        #      zmean[c,] <- apply(analysis$ensmean[pos,],2,mean)
        c=c+1
      }
      latlist <- unique(analysis$lat)[which(which(unique(analysis$lat) > -30) %in% 
                                                  which(unique(analysis$lat) < 30))]
      pos3 <- match(round(latlist,digits=8),round(unique(analysis$lat),digits=8)) 
      zmean <- zmean[pos3,]
      hc <- apply(zmean,2,min) # max upward wind is neg. omega value -> min function
      # cdo -s fldmax -sellonlatbox,0,360,0,30 -zonmean -sellevel,50000 -selvar,stream $f ${STOR_DIR}/HC.$f
      
      # Hadley cell poleward extend
      # Hind[6,which(analysis$names=="u200")] <- H.giorgi[which(giorgi.short == 'HCL'),
      #                                               which(analysis$names=="u200")] 
      
      # ITCZ location
      #Hind[7,which(analysis$names=="omega500")] <- H.giorgi[which(giorgi.short == 'ITCZ'),
      #                                                   which(analysis$names=="omega500")] 
      c=1
      zmean <- matrix(0,nrow=length(unique(analysis$lat)),ncol=2)
      for (l in unique(analysis$lat)) {
        pos1 <- which(trunc(analysis$lat,digits=8)==trunc(l,digits=8))
        pos2 <- which(analysis$names=='omega500')
        pos <- intersect(pos1,pos2)
        if (length(analysis$ensmean[pos,])==2) { zmean[c,] <- analysis$ensmean[pos,] }
        if (length(analysis$ensmean[pos,])>2) { zmean[c,] <- apply(analysis$ensmean[pos,],2,mean) }
        c=c+1
      }
      latlist <- unique(analysis$lat)[which(which(unique(analysis$lat) > -30) %in% 
                                                  which(unique(analysis$lat) < 30))]
      pos3 <- match(round(latlist,digits=8),round(unique(analysis$lat),digits=8)) 
      zmean <- zmean[pos3,]
      itcz <- apply(zmean[2:(nrow(zmean)-1),],2,min) 
      dlat <-  mean(latlist[1:(length(latlist)-1)]-latlist[2:(length(latlist))])    
      gridminpos1 <- which(round(zmean[,1],digits=8)==round(itcz[1],digits=8)) 
      gridmin1 <- latlist[gridminpos1]
      gridminpos2 <- which(round(zmean[,2],digits=8)==round(itcz[2],digits=8)) 
      gridmin2 <- latlist[gridminpos2]
      latmin1 <- gridmin1 - dlat * (((zmean[(gridminpos1-1),1])-(zmean[(gridminpos1+1),1]))/
                                      (2*(abs(itcz[1]-min(c(zmean[(gridminpos1-1),1],zmean[(gridminpos1+1),1]))))))
      latmin2 <- gridmin2 - dlat * (((zmean[(gridminpos2-1),2])-(zmean[(gridminpos2+1),2]))/
                                      (2*(abs(itcz[2]-min(c(zmean[(gridminpos2-1),2],zmean[(gridminpos2+1),2]))))))
      if (length(latmin1) > 1) {latmin1 <- mean(latmin1)}
      if (length(latmin2) > 1) {latmin2 <- mean(latmin2)}
      itcz <- c(latmin1,latmin2)
      
      # subtropical jet (u200)
      #Hind[8,which(analysis$names=="u200")] <- H.giorgi[which(giorgi.short == 'SJ'),
      #                                               which(analysis$names=="u200")] 
      c=1
      zmean <- matrix(0,nrow=length(unique(analysis$lat)),ncol=2)
      for (l in unique(analysis$lat)) {
        pos1 <- which(trunc(analysis$lat,digits=8)==trunc(l,digits=8))
        pos2 <- which(analysis$names=='u200')
        pos <- intersect(pos1,pos2)
        if (length(analysis$ensmean[pos,])==2) { zmean[c,] <- analysis$ensmean[pos,] }
        if (length(analysis$ensmean[pos,])>2) { zmean[c,] <- apply(analysis$ensmean[pos,],2,mean) }
        c=c+1
      }
      latlist <- unique(analysis$lat)[which(which(unique(analysis$lat) > 0) %in% which(unique(analysis$lat) < 57))]
      pos3 <- match(round(latlist,digits=8),round(unique(analysis$lat),digits=8)) 
      zmean <- zmean[pos3,]
      sj <- apply(zmean[2:(nrow(zmean)-1),],2,max) # max uwind
      dlat <-  mean(latlist[1:(length(latlist)-1)]-latlist[2:(length(latlist))])    
      gridminpos1 <- which(round(zmean[,1],digits=8)==round(sj[1],digits=8)) 
      gridmin1 <- latlist[gridminpos1]
      gridminpos2 <- which(round(zmean[,2],digits=8)==round(sj[2],digits=8)) 
      gridmin2 <- latlist[gridminpos2]
      zmean <- zmean * -1 # because next equation is to find minimum and here we search for max
      sj <- sj * -1
      latmin1 <- gridmin1 - dlat * (((zmean[(gridminpos1-1),1])-(zmean[(gridminpos1+1),1]))/
                                      (2*(abs(sj[1]-max(c(zmean[(gridminpos1-1),1],zmean[(gridminpos1+1),1]))))))
      latmin2 <- gridmin2 - dlat * (((zmean[(gridminpos2-1),2])-(zmean[(gridminpos2+1),2]))/
                                      (2*(abs(sj[2]-max(c(zmean[(gridminpos2-1),2],zmean[(gridminpos2+1),2]))))))
      if (length(latmin1) > 1) {latmin1 <- mean(latmin1)}
      if (length(latmin2) > 1) {latmin2 <- mean(latmin2)}
      sj_u200 <- c(latmin1,latmin2)
      # cdo -s fldmax -sellonlatbox,0,360,0,50 -zonmean -sellevel,20000 -selvar,u $f ${STOR_DIR}/SJ.$f
      
      # subtropical jet (slp)
      c=1
      zmean <- matrix(0,nrow=length(unique(analysis$lat)),ncol=2)
      for (l in unique(analysis$lat)) {
        pos1 <- which(trunc(analysis$lat,digits=8)==trunc(l,digits=8))
        pos2 <- which(analysis$names=='slp')
        pos <- intersect(pos1,pos2)
        if (length(analysis$ensmean[pos,])==2) { zmean[c,] <- analysis$ensmean[pos,] }
        if (length(analysis$ensmean[pos,])>2) { zmean[c,] <- apply(analysis$ensmean[pos,],2,mean) }
        c=c+1
      }
      latlist <- unique(analysis$lat)[which(which(unique(analysis$lat) > 0) %in% which(unique(analysis$lat) < 57))]
      pos3 <- match(round(latlist,digits=8),round(unique(analysis$lat),digits=8)) 
      zmean <- zmean[pos3,]
      sj <- apply(zmean[2:(nrow(zmean)-1),],2,max) # max slp
      dlat <-  mean(latlist[1:(length(latlist)-1)]-latlist[2:(length(latlist))])    
      gridminpos1 <- which(round(zmean[,1],digits=8)==round(sj[1],digits=8)) 
      gridmin1 <- latlist[gridminpos1]
      gridminpos2 <- which(round(zmean[,2],digits=8)==round(sj[2],digits=8)) 
      gridmin2 <- latlist[gridminpos2]
      zmean <- zmean * -1 # because next equation is to find minimum and here we search for max
      sj <- sj * -1
      latmin1 <- gridmin1 - dlat * (((zmean[(gridminpos1-1),1])-(zmean[(gridminpos1+1),1]))/
                                      (2*(abs(sj[1]-max(c(zmean[(gridminpos1-1),1],zmean[(gridminpos1+1),1]))))))
      latmin2 <- gridmin2 - dlat * (((zmean[(gridminpos2-1),2])-(zmean[(gridminpos2+1),2]))/
                                      (2*(abs(sj[2]-max(c(zmean[(gridminpos2-1),2],zmean[(gridminpos2+1),2]))))))
      if (length(latmin1) > 1) {latmin1 <- mean(latmin1)}
      if (length(latmin2) > 1) {latmin2 <- mean(latmin2)}
      sj_slp <- c(latmin1,latmin2)    
      
      # stratospheric polar vortex
      pv_reg1 <- aind_tmp$ensmean[which(aind_tmp$names=="PV1.calc"),]
      pv_reg2 <- aind_tmp$ensmean[which(aind_tmp$names=="PV2.calc"),]
      pv <- pv_reg1 - pv_reg2 
      #cdo -s -sellevel,10000 -selvar,geopoth $f ${STOR_DIR}/gph100.$f
      #cdo -s sub -fldmean -sellonlatbox,0,360,75,90 ${STOR_DIR}/gph100.$f -fldmean -sellonlatbox,0,360,40,55 ${STOR_DIR}/gph100.$f ${STOR_DIR}/z100.$f
      
      # Pacific walker circulation
      #Hind[10,which(analysis$names=="omega500")] <- H.giorgi[which(giorgi.short == 'PWC1'),
      #                                                    which(analysis$names=="omega500")] 
      pwc_reg1 <- aind_tmp$ensmean[which(aind_tmp$names=="PWC1.calc"),]
      pwc_reg2 <- aind_tmp$ensmean[which(aind_tmp$names=="PWC2.calc"),]
      pwc <- pwc_reg1 - pwc_reg2 
      # cdo -s sellevel,50000 -selvar,omega $f ${STOR_DIR}/omega.$f
      # cdo -s sub -fldmean -sellonlatbox,180,260,-10,10 ${STOR_DIR}/omega.$f -fldmean -sellonlatbox,100,150,-10,10 ${STOR_DIR}/omega.$f ${STOR_DIR}/PWC.$f
      
      # Dynamic Indian Monsoon index
      dimi_reg1 <- aind_tmp$ensmean[which(aind_tmp$names=="DIMI1.calc"),]
      dimi_reg2 <- aind_tmp$ensmean[which(aind_tmp$names=="DIMI2.calc"),]
      dimi <- dimi_reg1 - dimi_reg2 
      #cdo -s sellevel,85000 -selvar,u $f ${STOR_DIR}/u850.$f
      #cdo -s sub -fldmean -sellonlatbox,40,80,5,15 ${STOR_DIR}/u850.$f -fldmean -sellonlatbox,70,90,20,30 ${STOR_DIR}/u850.$f ${STOR_DIR}/DIMI.$f
      
      
      # NAO based on anomalies
      nao_reg1 <- aind_tmp$ensmean[which(aind_tmp$names=="NAO1.calc"),]
      nao_reg2 <- aind_tmp$ensmean[which(aind_tmp$names=="NAO2.calc"),]
      nao <- nao_reg1 - nao_reg2 
      
      # PNA based on anomalies
      pna_reg1 <- aind_tmp$ensmean[which(aind_tmp$names=="PNA1.calc"),]
      pna_reg2 <- aind_tmp$ensmean[which(aind_tmp$names=="PNA2.calc"),]
      pna_reg3 <- aind_tmp$ensmean[which(aind_tmp$names=="PNA3.calc"),]
      pna_reg4 <- aind_tmp$ensmean[which(aind_tmp$names=="PNA4.calc"),]
      pna <- 0.25*(pna_reg1 - pna_reg2 + pna_reg3 -pna_reg4) 
      # PNA = 0.25*(Z[20◦ N,160◦ W] – Z[45◦ N,165◦ W] + Z[55◦ N,115◦ W] – Z[30◦ N,85◦ W]) 
      
      aind <- list(ensmean=matrix(0,nrow=length(indices),ncol=2), names=indices)
      aind$ensmean[which(indices=='ENH.temp2'),] <- aind_tmp$ensmean[which(aind$names=='ENH.temp2'),]
      aind$ensmean[which(indices=='NAM.temp2'),] <- aind_tmp$ensmean[which(aind$names=='NAM.temp2'),]
      aind$ensmean[which(indices=='SAM.temp2'),] <- aind_tmp$ensmean[which(aind$names=='SAM.temp2'),]
      aind$ensmean[which(indices=='AFR.temp2'),] <- aind_tmp$ensmean[which(aind$names=='AFR.temp2'),]
      aind$ensmean[which(indices=='ASI.temp2'),] <- aind_tmp$ensmean[which(aind$names=='ASI.temp2'),]
      aind$ensmean[which(indices=='AUS.temp2'),] <- aind_tmp$ensmean[which(aind$names=='AUS.temp2'),]
      aind$ensmean[which(indices=='ARC.temp2'),] <- aind_tmp$ensmean[which(aind$names=='ARC.temp2'),]
      aind$ensmean[which(indices=='ANT.temp2'),] <- aind_tmp$ensmean[which(aind$names=='ANT.temp2'),]
      aind$ensmean[which(indices=='NEU.temp2'),] <- aind_tmp$ensmean[which(aind$names=='NEU.temp2'),]
      aind$ensmean[which(indices=='MED.temp2'),] <- aind_tmp$ensmean[which(aind$names=='MED.temp2'),]
      aind$ensmean[which(indices=='GLO.temp2'),] <- aind_tmp$ensmean[which(aind$names=='GLO.temp2'),]
      aind$ensmean[which(indices=='NAM.precip'),] <- aind_tmp$ensmean[which(aind$names=='NAM.precip'),]
      aind$ensmean[which(indices=='SAM.precip'),] <- aind_tmp$ensmean[which(aind$names=='SAM.precip'),]
      aind$ensmean[which(indices=='AFR.precip'),] <- aind_tmp$ensmean[which(aind$names=='AFR.precip'),]
      aind$ensmean[which(indices=='ASI.precip'),] <- aind_tmp$ensmean[which(aind$names=='ASI.precip'),]
      aind$ensmean[which(indices=='AUS.precip'),] <- aind_tmp$ensmean[which(aind$names=='AUS.precip'),]
      aind$ensmean[which(indices=='ARC.precip'),] <- aind_tmp$ensmean[which(aind$names=='ARC.precip'),]
      aind$ensmean[which(indices=='ANT.precip'),] <- aind_tmp$ensmean[which(aind$names=='ANT.precip'),]
      aind$ensmean[which(indices=='NEU.precip'),] <- aind_tmp$ensmean[which(aind$names=='NEU.precip'),]
      aind$ensmean[which(indices=='MED.precip'),] <- aind_tmp$ensmean[which(aind$names=='MED.precip'),]
      aind$ensmean[which(indices=='SH.temp2'),] <- aind_tmp$ensmean[which(aind$names=='SH.temp2'),]
      aind$ensmean[which(indices=='NAM.slp'),] <- aind_tmp$ensmean[which(aind$names=='NAM.slp'),]
      aind$ensmean[which(indices=='SAM.slp'),] <- aind_tmp$ensmean[which(aind$names=='SAM.slp'),]
      aind$ensmean[which(indices=='AFR.slp'),] <- aind_tmp$ensmean[which(aind$names=='AFR.slp'),]
      aind$ensmean[which(indices=='ASI.slp'),] <- aind_tmp$ensmean[which(aind$names=='ASI.slp'),]
      aind$ensmean[which(indices=='AUS.slp'),] <- aind_tmp$ensmean[which(aind$names=='AUS.slp'),]
      aind$ensmean[which(indices=='ARC.slp'),] <- aind_tmp$ensmean[which(aind$names=='ARC.slp'),]
      aind$ensmean[which(indices=='ANT.slp'),] <- aind_tmp$ensmean[which(aind$names=='ANT.slp'),]
      aind$ensmean[which(indices=='NEU.slp'),] <- aind_tmp$ensmean[which(aind$names=='NEU.slp'),]
      aind$ensmean[which(indices=='MED.slp'),] <- aind_tmp$ensmean[which(aind$names=='MED.slp'),]
      aind$ensmean[which(indices=='NH.temp2'),] <- aind_tmp$ensmean[which(aind$names=='NH.temp2'),]
      aind$ensmean[which(indices=='EU.temp2'),] <- aind_tmp$ensmean[which(aind$names=='EU.temp2'),]
      aind$ensmean[which(indices=='EU.precip'),] <- aind_tmp$ensmean[which(aind$names=='EU.precip'),]
      aind$ensmean[which(indices=='EU.slp'),] <- aind_tmp$ensmean[which(aind$names=='EU.slp'),]
      aind$ensmean[which(indices=='HC.calc'),] <- hc #aind_tmp$ensmean[which(aind$names=='HC.calc'),]
      aind$ensmean[which(indices=='ITCZ.calc'),] <- itcz #aind_tmp$ensmean[which(aind$names=='ITCZ.calc'),]
      aind$ensmean[which(indices=='SJ_u200.calc'),] <- sj_u200 #aind_tmp$ensmean[which(aind$names=='SJ.calc'),]
      aind$ensmean[which(indices=='SJ_slp.calc'),] <- sj_slp
      aind$ensmean[which(indices=='PV.calc'),] <- pv #aind_tmp$ensmean[which(aind$names=='PV.calc'),]
      aind$ensmean[which(indices=='PWC.calc'),] <- pwc #aind_tmp$ensmean[which(aind$names=='PWC.calc'),]
      aind$ensmean[which(indices=='DIMI.calc'),] <- dimi #aind_tmp$ensmean[which(aind$names=='DIMI.calc'),]
      aind$ensmean[which(indices=='NAO.calc'),] <- nao #aind_tmp$ensmean[which(aind$names=='DIMI.calc'),] #c(0,0)
      aind$ensmean[which(indices=='PNA.calc'),] <- pna #aind_tmp$ensmean[which(aind$names=='DIMI.calc'),] #c(0,0)
#      aind$ensmean[which(indices=='DIMI'),] <- ana_ind$ensmean[which(ana_ind$names=='DIMI'),]
#      aind$ensmean[which(indices=='HC'),] <- ana_ind$ensmean[which(ana_ind$names=='HC'),]
#      aind$ensmean[which(indices=='SJ'),] <- ana_ind$ensmean[which(ana_ind$names=='SJ'),]
#      aind$ensmean[which(indices=='Z100'),] <- ana_ind$ensmean[which(ana_ind$names=='z100'),]
#      aind$ensmean[which(indices=='Z300'),] <- ana_ind$ensmean[which(ana_ind$names=='z300'),]
#      aind$ensmean[which(indices=='PWC'),] <- ana_ind$ensmean[which(ana_ind$names=='PWC'),]
      
      
      
      # BEM
      aind_tmp <- list(bem=Hind%*%analysis$bem, names=indices_tmp)
      
      # Hadley cell strength
      #    Hind[5,which(analysis$names=="stream")] <- H.giorgi[which(giorgi.short == 'HC'),
      #                                              which(analysis$names=="stream")] 
      # Hind[5,which(analysis$names=="omega500")] <- H.giorgi[which(giorgi.short == 'HC'),
      #                                                   which(analysis$names=="omega500")] 
      c=1
      zmean <- matrix(0,nrow=length(unique(analysis$lat)),ncol=2)
      for (l in unique(analysis$lat)) {
        pos1 <- which(trunc(analysis$lat,digits=8)==trunc(l,digits=8))
        pos2 <- which(analysis$names=='omega500')
        pos <- intersect(pos1,pos2)
        if (length(analysis$bem[pos,])==2) { zmean[c,] <- analysis$bem[pos,] }
        if (length(analysis$bem[pos,])>2) { zmean[c,] <- apply(analysis$bem[pos,],2,mean) }
        #      zmean[c,] <- apply(analysis$bem[pos,],2,mean)
        c=c+1
      }
      latlist <- unique(analysis$lat)[which(which(unique(analysis$lat) > -30) %in% 
                                                  which(unique(analysis$lat) < 30))]
      pos3 <- match(round(latlist,digits=8),round(unique(analysis$lat),digits=8)) 
      zmean <- zmean[pos3,]
      hc <- apply(zmean,2,min) # max upward wind is neg. omega value -> min function
      # cdo -s fldmax -sellonlatbox,0,360,0,30 -zonmean -sellevel,50000 -selvar,stream $f ${STOR_DIR}/HC.$f
      
      # Hadley cell poleward extend
      # Hind[6,which(analysis$names=="u200")] <- H.giorgi[which(giorgi.short == 'HCL'),
      #                                               which(analysis$names=="u200")] 
      
      # ITCZ location
      #Hind[7,which(analysis$names=="omega500")] <- H.giorgi[which(giorgi.short == 'ITCZ'),
      #                                                   which(analysis$names=="omega500")] 
      c=1
      zmean <- matrix(0,nrow=length(unique(analysis$lat)),ncol=2)
      for (l in unique(analysis$lat)) {
        pos1 <- which(trunc(analysis$lat,digits=8)==trunc(l,digits=8))
        pos2 <- which(analysis$names=='omega500')
        pos <- intersect(pos1,pos2)
        if (length(analysis$bem[pos,])==2) { zmean[c,] <- analysis$bem[pos,] }
        if (length(analysis$bem[pos,])>2) { zmean[c,] <- apply(analysis$bem[pos,],2,mean) }
        c=c+1
      }
      latlist <- unique(analysis$lat)[which(which(unique(analysis$lat) > -30) %in% 
                                                  which(unique(analysis$lat) < 30))]
      pos3 <- match(round(latlist,digits=8),round(unique(analysis$lat),digits=8)) 
      zmean <- zmean[pos3,]
      itcz <- apply(zmean[2:(nrow(zmean)-1),],2,min) 
      dlat <-  mean(latlist[1:(length(latlist)-1)]-latlist[2:(length(latlist))])    
      gridminpos1 <- which(round(zmean[,1],digits=8)==round(itcz[1],digits=8)) 
      gridmin1 <- latlist[gridminpos1]
      gridminpos2 <- which(round(zmean[,2],digits=8)==round(itcz[2],digits=8)) 
      gridmin2 <- latlist[gridminpos2]
      latmin1 <- gridmin1 - dlat * (((zmean[(gridminpos1-1),1])-(zmean[(gridminpos1+1),1]))/
                                      (2*(abs(itcz[1]-min(c(zmean[(gridminpos1-1),1],zmean[(gridminpos1+1),1]))))))
      latmin2 <- gridmin2 - dlat * (((zmean[(gridminpos2-1),2])-(zmean[(gridminpos2+1),2]))/
                                      (2*(abs(itcz[2]-min(c(zmean[(gridminpos2-1),2],zmean[(gridminpos2+1),2]))))))
      if (length(latmin1) > 1) {latmin1 <- mean(latmin1)}
      if (length(latmin2) > 1) {latmin2 <- mean(latmin2)}
      itcz <- c(latmin1,latmin2)
      
      # subtropical jet (u200)
      #Hind[8,which(analysis$names=="u200")] <- H.giorgi[which(giorgi.short == 'SJ'),
      #                                               which(analysis$names=="u200")] 
      c=1
      zmean <- matrix(0,nrow=length(unique(analysis$lat)),ncol=2)
      for (l in unique(analysis$lat)) {
        pos1 <- which(trunc(analysis$lat,digits=8)==trunc(l,digits=8))
        pos2 <- which(analysis$names=='u200')
        pos <- intersect(pos1,pos2)
        if (length(analysis$bem[pos,])==2) { zmean[c,] <- analysis$bem[pos,] }
        if (length(analysis$bem[pos,])>2) { zmean[c,] <- apply(analysis$bem[pos,],2,mean) }
        c=c+1
      }
      latlist <- unique(analysis$lat)[which(which(unique(analysis$lat) > 0) %in% which(unique(analysis$lat) < 57))]
      pos3 <- match(round(latlist,digits=8),round(unique(analysis$lat),digits=8)) 
      zmean <- zmean[pos3,]
      sj <- apply(zmean[2:(nrow(zmean)-1),],2,max) # max uwind
      dlat <-  mean(latlist[1:(length(latlist)-1)]-latlist[2:(length(latlist))])    
      gridminpos1 <- which(round(zmean[,1],digits=8)==round(sj[1],digits=8)) 
      gridmin1 <- latlist[gridminpos1]
      gridminpos2 <- which(round(zmean[,2],digits=8)==round(sj[2],digits=8)) 
      gridmin2 <- latlist[gridminpos2]
      zmean <- zmean * -1
      sj <- sj * -1
      latmin1 <- gridmin1 - dlat * (((zmean[(gridminpos1-1),1])-(zmean[(gridminpos1+1),1]))/
                                      (2*(abs(sj[1]-max(c(zmean[(gridminpos1-1),1],zmean[(gridminpos1+1),1]))))))
      latmin2 <- gridmin2 - dlat * (((zmean[(gridminpos2-1),2])-(zmean[(gridminpos2+1),2]))/
                                      (2*(abs(sj[2]-max(c(zmean[(gridminpos2-1),2],zmean[(gridminpos2+1),2]))))))
      if (length(latmin1) > 1) {latmin1 <- mean(latmin1)}
      if (length(latmin2) > 1) {latmin2 <- mean(latmin2)}
      sj_u200 <- c(latmin1,latmin2)
      # cdo -s fldmax -sellonlatbox,0,360,0,50 -zonmean -sellevel,20000 -selvar,u $f ${STOR_DIR}/SJ.$f
      
      # subtropical jet (slp)
      c=1
      zmean <- matrix(0,nrow=length(unique(analysis$lat)),ncol=2)
      for (l in unique(analysis$lat)) {
        pos1 <- which(trunc(analysis$lat,digits=8)==trunc(l,digits=8))
        pos2 <- which(analysis$names=='slp')
        pos <- intersect(pos1,pos2)
        if (length(analysis$bem[pos,])==2) { zmean[c,] <- analysis$bem[pos,] }
        if (length(analysis$bem[pos,])>2) { zmean[c,] <- apply(analysis$bem[pos,],2,mean) }
        c=c+1
      }
      latlist <- unique(analysis$lat)[which(which(unique(analysis$lat) > 0) %in% which(unique(analysis$lat) < 57))]
      pos3 <- match(round(latlist,digits=8),round(unique(analysis$lat),digits=8)) 
      zmean <- zmean[pos3,]
      sj <- apply(zmean[2:(nrow(zmean)-1),],2,max) # max slp
      dlat <-  mean(latlist[1:(length(latlist)-1)]-latlist[2:(length(latlist))])    
      gridminpos1 <- which(round(zmean[,1],digits=8)==round(sj[1],digits=8)) 
      gridmin1 <- latlist[gridminpos1]
      gridminpos2 <- which(round(zmean[,2],digits=8)==round(sj[2],digits=8)) 
      gridmin2 <- latlist[gridminpos2]
      zmean <- zmean * -1
      sj <- sj * -1
      latmin1 <- gridmin1 - dlat * (((zmean[(gridminpos1-1),1])-(zmean[(gridminpos1+1),1]))/
                                      (2*(abs(sj[1]-max(c(zmean[(gridminpos1-1),1],zmean[(gridminpos1+1),1]))))))
      latmin2 <- gridmin2 - dlat * (((zmean[(gridminpos2-1),2])-(zmean[(gridminpos2+1),2]))/
                                      (2*(abs(sj[2]-max(c(zmean[(gridminpos2-1),2],zmean[(gridminpos2+1),2]))))))
      if (length(latmin1) > 1) {latmin1 <- mean(latmin1)}
      if (length(latmin2) > 1) {latmin2 <- mean(latmin2)}
      sj_slp <- c(latmin1,latmin2)
      
      # stratospheric polar vortex
      pv_reg1 <- aind_tmp$bem[which(aind_tmp$names=="PV1.calc"),]
      pv_reg2 <- aind_tmp$bem[which(aind_tmp$names=="PV2.calc"),]
      pv <- pv_reg1 - pv_reg2 
      #cdo -s -sellevel,10000 -selvar,geopoth $f ${STOR_DIR}/gph100.$f
      #cdo -s sub -fldmean -sellonlatbox,0,360,75,90 ${STOR_DIR}/gph100.$f -fldmean -sellonlatbox,0,360,40,55 ${STOR_DIR}/gph100.$f ${STOR_DIR}/z100.$f
      
      # Pacific walker circulation
      #Hind[10,which(analysis$names=="omega500")] <- H.giorgi[which(giorgi.short == 'PWC1'),
      #                                                    which(analysis$names=="omega500")] 
      pwc_reg1 <- aind_tmp$bem[which(aind_tmp$names=="PWC1.calc"),]
      pwc_reg2 <- aind_tmp$bem[which(aind_tmp$names=="PWC2.calc"),]
      pwc <- pwc_reg1 - pwc_reg2 
      # cdo -s sellevel,50000 -selvar,omega $f ${STOR_DIR}/omega.$f
      # cdo -s sub -fldmean -sellonlatbox,180,260,-10,10 ${STOR_DIR}/omega.$f -fldmean -sellonlatbox,100,150,-10,10 ${STOR_DIR}/omega.$f ${STOR_DIR}/PWC.$f
      
      # Dynamic Indian Monsoon index
      dimi_reg1 <- aind_tmp$bem[which(aind_tmp$names=="DIMI1.calc"),]
      dimi_reg2 <- aind_tmp$bem[which(aind_tmp$names=="DIMI2.calc"),]
      dimi <- dimi_reg1 - dimi_reg2 
      #cdo -s sellevel,85000 -selvar,u $f ${STOR_DIR}/u850.$f
      #cdo -s sub -fldmean -sellonlatbox,40,80,5,15 ${STOR_DIR}/u850.$f -fldmean -sellonlatbox,70,90,20,30 ${STOR_DIR}/u850.$f ${STOR_DIR}/DIMI.$f
      
      
      # NAO based on anomalies
      nao_reg1 <- aind_tmp$bem[which(aind_tmp$names=="NAO1.calc"),]
      nao_reg2 <- aind_tmp$bem[which(aind_tmp$names=="NAO2.calc"),]
      nao <- nao_reg1 - nao_reg2 
      
      # PNA based on anomalies
      pna_reg1 <- aind_tmp$bem[which(aind_tmp$names=="PNA1.calc"),]
      pna_reg2 <- aind_tmp$bem[which(aind_tmp$names=="PNA2.calc"),]
      pna_reg3 <- aind_tmp$bem[which(aind_tmp$names=="PNA3.calc"),]
      pna_reg4 <- aind_tmp$bem[which(aind_tmp$names=="PNA4.calc"),]
      pna <- 0.25*(pna_reg1 - pna_reg2 + pna_reg3 -pna_reg4) 
      # PNA = 0.25*(Z[20◦ N,160◦ W] – Z[45◦ N,165◦ W] + Z[55◦ N,115◦ W] – Z[30◦ N,85◦ W]) 
      
      aind$bem <- matrix(0,nrow=length(indices),ncol=2)
      #aind <- list(bem=matrix(0,nrow=length(indices),ncol=2), names=indices)
      aind$bem[which(indices=='ENH.temp2'),] <- aind_tmp$bem[which(aind$names=='ENH.temp2'),]
      aind$bem[which(indices=='NAM.temp2'),] <- aind_tmp$bem[which(aind$names=='NAM.temp2'),]
      aind$bem[which(indices=='SAM.temp2'),] <- aind_tmp$bem[which(aind$names=='SAM.temp2'),]
      aind$bem[which(indices=='AFR.temp2'),] <- aind_tmp$bem[which(aind$names=='AFR.temp2'),]
      aind$bem[which(indices=='ASI.temp2'),] <- aind_tmp$bem[which(aind$names=='ASI.temp2'),]
      aind$bem[which(indices=='AUS.temp2'),] <- aind_tmp$bem[which(aind$names=='AUS.temp2'),]
      aind$bem[which(indices=='ARC.temp2'),] <- aind_tmp$bem[which(aind$names=='ARC.temp2'),]
      aind$bem[which(indices=='ANT.temp2'),] <- aind_tmp$bem[which(aind$names=='ANT.temp2'),]
      aind$bem[which(indices=='NEU.temp2'),] <- aind_tmp$bem[which(aind$names=='NEU.temp2'),]
      aind$bem[which(indices=='MED.temp2'),] <- aind_tmp$bem[which(aind$names=='MED.temp2'),]
      aind$bem[which(indices=='GLO.temp2'),] <- aind_tmp$bem[which(aind$names=='GLO.temp2'),]
      aind$bem[which(indices=='NAM.precip'),] <- aind_tmp$bem[which(aind$names=='NAM.precip'),]
      aind$bem[which(indices=='SAM.precip'),] <- aind_tmp$bem[which(aind$names=='SAM.precip'),]
      aind$bem[which(indices=='AFR.precip'),] <- aind_tmp$bem[which(aind$names=='AFR.precip'),]
      aind$bem[which(indices=='ASI.precip'),] <- aind_tmp$bem[which(aind$names=='ASI.precip'),]
      aind$bem[which(indices=='AUS.precip'),] <- aind_tmp$bem[which(aind$names=='AUS.precip'),]
      aind$bem[which(indices=='ARC.precip'),] <- aind_tmp$bem[which(aind$names=='ARC.precip'),]
      aind$bem[which(indices=='ANT.precip'),] <- aind_tmp$bem[which(aind$names=='ANT.precip'),]
      aind$bem[which(indices=='NEU.precip'),] <- aind_tmp$bem[which(aind$names=='NEU.precip'),]
      aind$bem[which(indices=='MED.precip'),] <- aind_tmp$bem[which(aind$names=='MED.precip'),]
      aind$bem[which(indices=='SH.temp2'),] <- aind_tmp$bem[which(aind$names=='SH.temp2'),]
      aind$bem[which(indices=='NAM.slp'),] <- aind_tmp$bem[which(aind$names=='NAM.slp'),]
      aind$bem[which(indices=='SAM.slp'),] <- aind_tmp$bem[which(aind$names=='SAM.slp'),]
      aind$bem[which(indices=='AFR.slp'),] <- aind_tmp$bem[which(aind$names=='AFR.slp'),]
      aind$bem[which(indices=='ASI.slp'),] <- aind_tmp$bem[which(aind$names=='ASI.slp'),]
      aind$bem[which(indices=='AUS.slp'),] <- aind_tmp$bem[which(aind$names=='AUS.slp'),]
      aind$bem[which(indices=='ARC.slp'),] <- aind_tmp$bem[which(aind$names=='ARC.slp'),]
      aind$bem[which(indices=='ANT.slp'),] <- aind_tmp$bem[which(aind$names=='ANT.slp'),]
      aind$bem[which(indices=='NEU.slp'),] <- aind_tmp$bem[which(aind$names=='NEU.slp'),]
      aind$bem[which(indices=='MED.slp'),] <- aind_tmp$bem[which(aind$names=='MED.slp'),]
      aind$bem[which(indices=='NH.temp2'),] <- aind_tmp$bem[which(aind$names=='NH.temp2'),]
      aind$bem[which(indices=='EU.temp2'),] <- aind_tmp$bem[which(aind$names=='EU.temp2'),]
      aind$bem[which(indices=='EU.precip'),] <- aind_tmp$bem[which(aind$names=='EU.precip'),]
      aind$bem[which(indices=='EU.slp'),] <- aind_tmp$bem[which(aind$names=='EU.slp'),]
      aind$bem[which(indices=='HC.calc'),] <- hc #aind_tmp$bem[which(aind$names=='HC.calc'),]
      aind$bem[which(indices=='ITCZ.calc'),] <- itcz #aind_tmp$bem[which(aind$names=='ITCZ.calc'),]
      aind$bem[which(indices=='SJ_u200.calc'),] <- sj_u200 #aind_tmp$bem[which(aind$names=='SJ.calc'),]
      aind$bem[which(indices=='SJ_slp.calc'),] <- sj_slp
      aind$bem[which(indices=='PV.calc'),] <- pv #aind_tmp$bem[which(aind$names=='PV.calc'),]
      aind$bem[which(indices=='PWC.calc'),] <- pwc #aind_tmp$bem[which(aind$names=='PWC.calc'),]
      aind$bem[which(indices=='DIMI.calc'),] <- dimi #aind_tmp$bem[which(aind$names=='DIMI.calc'),]
      aind$bem[which(indices=='NAO.calc'),] <- nao #aind_tmp$bem[which(aind$names=='DIMI.calc'),] #c(0,0)
      aind$bem[which(indices=='PNA.calc'),] <- pna #aind_tmp$bem[which(aind$names=='DIMI.calc'),] #c(0,0)
#      aind$bem[which(indices=='DIMI'),] <- ana_ind$bem[which(ana_ind$names=='DIMI'),]
#      aind$bem[which(indices=='HC'),] <- ana_ind$bem[which(ana_ind$names=='HC'),]
#      aind$bem[which(indices=='SJ'),] <- ana_ind$bem[which(ana_ind$names=='SJ'),]
#      aind$bem[which(indices=='Z100'),] <- ana_ind$bem[which(ana_ind$names=='z100'),]
#      aind$bem[which(indices=='Z300'),] <- ana_ind$bem[which(ana_ind$names=='z300'),]
#      aind$bem[which(indices=='PWC'),] <- ana_ind$bem[which(ana_ind$names=='PWC'),]
      
      
      
      
      
      # data for ensemble member
      #    aind$data <- array(Hind %*% array(analysis$data, c(nrow(analysis$data), 
      #                       length(analysis$data)/nrow(analysis$data))), c(nrow(Hind), 
      #                       dim(analysis$data)[2:3]))
      aind_tmp <- list(data=array(Hind %*% array(analysis$data, c(nrow(analysis$data), 
                                                                      length(analysis$data)/nrow(analysis$data))), c(nrow(Hind), 
                                                                                                                             dim(analysis$data)[2:3])), names=indices_tmp)
      # Hadley cell strength
      hc <- matrix(0,nrow=nmem,ncol=2)
      for (m in 1:nmem) {
        c=1
        zmean <- matrix(0,nrow=length(unique(analysis$lat)),ncol=2)
        for (l in unique(analysis$lat)) {
          pos1 <- which(trunc(analysis$lat,digits=8)==trunc(l,digits=8))
          pos2 <- which(analysis$names=='omega500')
          pos <- intersect(pos1,pos2)
          if (length(analysis$data[pos,,m])==2) { zmean[c,] <- analysis$data[pos,,m] }
          if (length(analysis$data[pos,,m])>2) { zmean[c,] <- apply(analysis$data[pos,,m],2,mean) }
          c=c+1
        }
        latlist <- unique(analysis$lat)[which(which(unique(analysis$lat) > -30) %in% 
                                                    which(unique(analysis$lat) < 30))]
        pos3 <- match(round(latlist,digits=8),round(unique(analysis$lat),digits=8)) 
        zmean <- zmean[pos3,]
        hc[m,] <- apply(zmean,2,min) # max upward wind is neg. omega value -> min function
      }
      # Hadley cell poleward extend
      
      # ITCZ location
      itcz <- matrix(0,nrow=nmem,ncol=2)
      for (m in 1:nmem) {
        c=1
        zmean <- matrix(0,nrow=length(unique(analysis$lat)),ncol=2)
        for (l in unique(analysis$lat)) {
          pos1 <- which(trunc(analysis$lat,digits=8)==trunc(l,digits=8))
          pos2 <- which(analysis$names=='omega500')
          pos <- intersect(pos1,pos2)
          if (length(analysis$data[pos,,m])==2) { zmean[c,] <- analysis$data[pos,,m] }
          if (length(analysis$data[pos,,m])>2) { zmean[c,] <- apply(analysis$data[pos,,m],2,mean) }
          c=c+1
        }
        latlist <- unique(analysis$lat)[which(which(unique(analysis$lat) > -30) %in% 
                                                    which(unique(analysis$lat) < 30))]
        pos3 <- match(round(latlist,digits=8),round(unique(analysis$lat),digits=8)) 
        zmean <- zmean[pos3,]
        itcz_tmp <- apply(zmean[2:(nrow(zmean)-1),],2,min) 
        dlat <-  mean(latlist[1:(length(latlist)-1)]-latlist[2:(length(latlist))])    
        gridminpos1 <- which(round(zmean[,1],digits=8)==round(itcz_tmp[1],digits=8)) 
        gridmin1 <- latlist[gridminpos1]
        gridminpos2 <- which(round(zmean[,2],digits=8)==round(itcz_tmp[2],digits=8)) 
        gridmin2 <- latlist[gridminpos2]
        latmin1 <- gridmin1 - dlat * (((zmean[(gridminpos1-1),1])-(zmean[(gridminpos1+1),1]))/
                                        (2*(abs(itcz_tmp[1]-max(c(zmean[(gridminpos1-1),1],zmean[(gridminpos1+1),1]))))))
        latmin2 <- gridmin2 - dlat * (((zmean[(gridminpos2-1),2])-(zmean[(gridminpos2+1),2]))/
                                        (2*(abs(itcz_tmp[2]-max(c(zmean[(gridminpos2-1),2],zmean[(gridminpos2+1),2]))))))
        if (length(latmin1) > 1) {latmin1 <- mean(latmin1)}
        if (length(latmin2) > 1) {latmin2 <- mean(latmin2)}
        itcz[m,] <- c(latmin1,latmin2)
      }
      
      # subtropical jet (u200)
      sj_u200 <- matrix(0,nrow=nmem,ncol=2)
      for (m in 1:nmem) {
        c=1
        zmean <- matrix(0,nrow=length(unique(analysis$lat)),ncol=2)
        for (l in unique(analysis$lat)) {
          pos1 <- which(trunc(analysis$lat,digits=8)==trunc(l,digits=8))
          pos2 <- which(analysis$names=='u200')
          pos <- intersect(pos1,pos2)
          if (length(analysis$data[pos,,m])==2) { zmean[c,] <- analysis$data[pos,,m] }
          if (length(analysis$data[pos,,m])>2) { zmean[c,] <- apply(analysis$data[pos,,m],2,mean) }
          c=c+1
        }
        latlist <- unique(analysis$lat)[which(which(unique(analysis$lat) > 0) %in% 
                                                    which(unique(analysis$lat) < 57))]
        pos3 <- match(round(latlist,digits=8),round(unique(analysis$lat),digits=8)) 
        zmean <- zmean[pos3,]
        sj_tmp <- apply(zmean[2:(nrow(zmean)-1),],2,max) # max uwind
        dlat <-  mean(latlist[1:(length(latlist)-1)]-latlist[2:(length(latlist))])    
        gridminpos1 <- which(round(zmean[,1],digits=8)==round(sj_tmp[1],digits=8)) 
        gridmin1 <- latlist[gridminpos1]
        gridminpos2 <- which(round(zmean[,2],digits=8)==round(sj_tmp[2],digits=8)) 
        gridmin2 <- latlist[gridminpos2]
        zmean <- zmean * -1 # because next equation is to find minimum and here we search for max
        sj_tmp <- sj_tmp * -1 
        latmin1 <- gridmin1 - dlat * (((zmean[(gridminpos1-1),1])-(zmean[(gridminpos1+1),1]))/
                                        (2*(abs(sj_tmp[1]-max(c(zmean[(gridminpos1-1),1],zmean[(gridminpos1+1),1]))))))
        latmin2 <- gridmin2 - dlat * (((zmean[(gridminpos2-1),2])-(zmean[(gridminpos2+1),2]))/
                                        (2*(abs(sj_tmp[2]-max(c(zmean[(gridminpos2-1),2],zmean[(gridminpos2+1),2]))))))
        if (length(latmin1) > 1) {latmin1 <- mean(latmin1)}
        if (length(latmin2) > 1) {latmin2 <- mean(latmin2)}
        sj_u200[m,] <- c(latmin1,latmin2)
      }
      
      # subtropical jet (slp)
      sj_slp <- matrix(0,nrow=nmem,ncol=2)
      for (m in 1:nmem) {
        c=1
        zmean <- matrix(0,nrow=length(unique(analysis$lat)),ncol=2)
        for (l in unique(analysis$lat)) {
          pos1 <- which(trunc(analysis$lat,digits=8)==trunc(l,digits=8))
          pos2 <- which(analysis$names=='slp')
          pos <- intersect(pos1,pos2)
          if (length(analysis$data[pos,,m])==2) { zmean[c,] <- analysis$data[pos,,m] }
          if (length(analysis$data[pos,,m])>2) { zmean[c,] <- apply(analysis$data[pos,,m],2,mean) }
          c=c+1
        }
        latlist <- unique(analysis$lat)[which(which(unique(analysis$lat) > 0) %in% 
                                                    which(unique(analysis$lat) < 57))]
        pos3 <- match(round(latlist,digits=8),round(unique(analysis$lat),digits=8)) 
        zmean <- zmean[pos3,]
        sj_tmp <- apply(zmean[2:(nrow(zmean)-1),],2,max) # max uwind
        dlat <-  mean(latlist[1:(length(latlist)-1)]-latlist[2:(length(latlist))])    
        gridminpos1 <- which(round(zmean[,1],digits=8)==round(sj_tmp[1],digits=8)) 
        gridmin1 <- latlist[gridminpos1]
        gridminpos2 <- which(round(zmean[,2],digits=8)==round(sj_tmp[2],digits=8)) 
        gridmin2 <- latlist[gridminpos2]
        zmean <- zmean * -1 # because next equation is to find minimum and here we search for max
        sj_tmp <- sj_tmp * -1 
        latmin1 <- gridmin1 - dlat * (((zmean[(gridminpos1-1),1])-(zmean[(gridminpos1+1),1]))/
                                        (2*(abs(sj_tmp[1]-max(c(zmean[(gridminpos1-1),1],zmean[(gridminpos1+1),1]))))))
        latmin2 <- gridmin2 - dlat * (((zmean[(gridminpos2-1),2])-(zmean[(gridminpos2+1),2]))/
                                        (2*(abs(sj_tmp[2]-max(c(zmean[(gridminpos2-1),2],zmean[(gridminpos2+1),2]))))))
        if (length(latmin1) > 1) {latmin1 <- mean(latmin1)}
        if (length(latmin2) > 1) {latmin2 <- mean(latmin2)}
        sj_slp[m,] <- c(latmin1,latmin2)
      }
      
      # stratospheric polar vortex
      pv <- matrix(0,nrow=nmem,ncol=2)
      for (m in 1:nmem) {
        pv_reg1 <- aind_tmp$data[which(aind_tmp$names=="PV1.calc"),,m]
        pv_reg2 <- aind_tmp$data[which(aind_tmp$names=="PV2.calc"),,m]
        pv[m,] <- pv_reg1 - pv_reg2 
      }
      
      # Pacific walker circulation
      pwc <- matrix(0,nrow=nmem,ncol=2)
      for (m in 1:nmem) {
        pwc_reg1 <- aind_tmp$data[which(aind_tmp$names=="PWC1.calc"),,m]
        pwc_reg2 <- aind_tmp$data[which(aind_tmp$names=="PWC2.calc"),,m]
        pwc[m,] <- pwc_reg1 - pwc_reg2 
      }
      
      # Dynamic Indian Monsoon index
      dimi <- matrix(0,nrow=nmem,ncol=2)
      for (m in 1:nmem) {
        dimi_reg1 <- aind_tmp$data[which(aind_tmp$names=="DIMI1.calc"),,m]
        dimi_reg2 <- aind_tmp$data[which(aind_tmp$names=="DIMI2.calc"),,m]
        dimi[m,] <- dimi_reg1 - dimi_reg2 
      }
      
      # NAO based on anomalies
      nao <- matrix(0,nrow=nmem,ncol=2)
      for (m in 1:nmem) {
        nao_reg1 <- aind_tmp$data[which(aind_tmp$names=="NAO1.calc"),,m]
        nao_reg2 <- aind_tmp$data[which(aind_tmp$names=="NAO2.calc"),,m]
        nao[m,] <- nao_reg1 - nao_reg2 
      }
      
      # PNA based on anomalies
      pna <- matrix(0,nrow=nmem,ncol=2)
      for (m in 1:nmem) {
        pna_reg1 <- aind_tmp$data[which(aind_tmp$names=="PNA1.calc"),,m]
        pna_reg2 <- aind_tmp$data[which(aind_tmp$names=="PNA2.calc"),,m]
        pna_reg3 <- aind_tmp$data[which(aind_tmp$names=="PNA3.calc"),,m]
        pna_reg4 <- aind_tmp$data[which(aind_tmp$names=="PNA4.calc"),,m]
        pna[m,] <- 0.25*(pna_reg1 - pna_reg2 + pna_reg3 -pna_reg4) 
      }
      
      aind$data <- array(0,c(length(indices),2,nmem))
      for (m in 1:nmem) {
        aind$data[which(indices=='ENH.temp2'),,m] <- aind_tmp$data[which(aind$names=='ENH.temp2'),,m]
        aind$data[which(indices=='NAM.temp2'),,m] <- aind_tmp$data[which(aind$names=='NAM.temp2'),,m]
        aind$data[which(indices=='SAM.temp2'),,m] <- aind_tmp$data[which(aind$names=='SAM.temp2'),,m]
        aind$data[which(indices=='AFR.temp2'),,m] <- aind_tmp$data[which(aind$names=='AFR.temp2'),,m]
        aind$data[which(indices=='ASI.temp2'),,m] <- aind_tmp$data[which(aind$names=='ASI.temp2'),,m]
        aind$data[which(indices=='AUS.temp2'),,m] <- aind_tmp$data[which(aind$names=='AUS.temp2'),,m]
        aind$data[which(indices=='ARC.temp2'),,m] <- aind_tmp$data[which(aind$names=='ARC.temp2'),,m]
        aind$data[which(indices=='ANT.temp2'),,m] <- aind_tmp$data[which(aind$names=='ANT.temp2'),,m]
        aind$data[which(indices=='NEU.temp2'),,m] <- aind_tmp$data[which(aind$names=='NEU.temp2'),,m]
        aind$data[which(indices=='MED.temp2'),,m] <- aind_tmp$data[which(aind$names=='MED.temp2'),,m]
        aind$data[which(indices=='GLO.temp2'),,m] <- aind_tmp$data[which(aind$names=='GLO.temp2'),,m]
        aind$data[which(indices=='NAM.precip'),,m] <- aind_tmp$data[which(aind$names=='NAM.precip'),,m]
        aind$data[which(indices=='SAM.precip'),,m] <- aind_tmp$data[which(aind$names=='SAM.precip'),,m]
        aind$data[which(indices=='AFR.precip'),,m] <- aind_tmp$data[which(aind$names=='AFR.precip'),,m]
        aind$data[which(indices=='ASI.precip'),,m] <- aind_tmp$data[which(aind$names=='ASI.precip'),,m]
        aind$data[which(indices=='AUS.precip'),,m] <- aind_tmp$data[which(aind$names=='AUS.precip'),,m]
        aind$data[which(indices=='ARC.precip'),,m] <- aind_tmp$data[which(aind$names=='ARC.precip'),,m]
        aind$data[which(indices=='ANT.precip'),,m] <- aind_tmp$data[which(aind$names=='ANT.precip'),,m]
        aind$data[which(indices=='NEU.precip'),,m] <- aind_tmp$data[which(aind$names=='NEU.precip'),,m]
        aind$data[which(indices=='MED.precip'),,m] <- aind_tmp$data[which(aind$names=='MED.precip'),,m]
        aind$data[which(indices=='SH.temp2'),,m] <- aind_tmp$data[which(aind$names=='SH.temp2'),,m]
        aind$data[which(indices=='NAM.slp'),,m] <- aind_tmp$data[which(aind$names=='NAM.slp'),,m]
        aind$data[which(indices=='SAM.slp'),,m] <- aind_tmp$data[which(aind$names=='SAM.slp'),,m]
        aind$data[which(indices=='AFR.slp'),,m] <- aind_tmp$data[which(aind$names=='AFR.slp'),,m]
        aind$data[which(indices=='ASI.slp'),,m] <- aind_tmp$data[which(aind$names=='ASI.slp'),,m]
        aind$data[which(indices=='AUS.slp'),,m] <- aind_tmp$data[which(aind$names=='AUS.slp'),,m]
        aind$data[which(indices=='ARC.slp'),,m] <- aind_tmp$data[which(aind$names=='ARC.slp'),,m]
        aind$data[which(indices=='ANT.slp'),,m] <- aind_tmp$data[which(aind$names=='ANT.slp'),,m]
        aind$data[which(indices=='NEU.slp'),,m] <- aind_tmp$data[which(aind$names=='NEU.slp'),,m]
        aind$data[which(indices=='MED.slp'),,m] <- aind_tmp$data[which(aind$names=='MED.slp'),,m]    
        aind$data[which(indices=='NH.temp2'),,m] <- aind_tmp$data[which(aind$names=='NH.temp2'),,m]
        aind$data[which(indices=='EU.temp2'),,m] <- aind_tmp$data[which(aind$names=='EU.temp2'),,m]
        aind$data[which(indices=='EU.precip'),,m] <- aind_tmp$data[which(aind$names=='EU.precip'),,m]
        aind$data[which(indices=='EU.slp'),,m] <- aind_tmp$data[which(aind$names=='EU.slp'),,m]
        aind$data[which(indices=='HC.calc'),,m] <- hc[m,] #aind_tmp$data[which(aind$names=='HC.calc'),,m]
        aind$data[which(indices=='ITCZ.calc'),,m] <- itcz[m,] #aind_tmp$data[which(aind$names=='ITCZ.calc'),,m]
        aind$data[which(indices=='SJ_u200.calc'),,m] <- sj_u200[m,] #aind_tmp$data[which(aind$names=='SJ.calc'),,m]
        aind$data[which(indices=='SJ_slp.calc'),,m] <- sj_slp[m,]
        aind$data[which(indices=='PV.calc'),,m] <- pv[m,] #aind_tmp$data[which(aind$names=='PV.calc'),,m]
        aind$data[which(indices=='PWC.calc'),,m] <- pwc[m,] #aind_tmp$data[which(aind$names=='PWC.calc'),,m]
        aind$data[which(indices=='DIMI.calc'),,m] <- dimi[m,] #aind_tmp$data[which(aind$names=='DIMI.calc'),,m]
        aind$data[which(indices=='NAO.calc'),,m] <- nao[m,] #aind_tmp$data[which(aind$names=='DIMI.calc'),,m] #c(0,0)
        aind$data[which(indices=='PNA.calc'),,m] <- pna[m,] #aind_tmp$data[which(aind$names=='DIMI.calc'),,m] #c(0,0)
#        aind$data[which(indices=='DIMI'),,m] <- ana_ind$data[which(ana_ind$names=='DIMI'),,m]
#        aind$data[which(indices=='HC'),,m] <- ana_ind$data[which(ana_ind$names=='HC'),,m]
#        aind$data[which(indices=='SJ'),,m] <- ana_ind$data[which(ana_ind$names=='SJ'),,m]
#        aind$data[which(indices=='Z100'),,m] <- ana_ind$data[which(ana_ind$names=='z100'),,m]
#        aind$data[which(indices=='Z300'),,m] <- ana_ind$data[which(ana_ind$names=='z300'),,m]
#        aind$data[which(indices=='PWC'),,m] <- ana_ind$data[which(ana_ind$names=='PWC'),,m]
      }
      
#      aind$data[which(indices=='DIMI'),,] <- ana_ind$data[which(ana_ind$names=='DIMI'),,]
#      aind$data[which(indices=='HC'),,] <- ana_ind$data[which(ana_ind$names=='HC'),,]
#      aind$data[which(indices=='SJ'),,] <- ana_ind$data[which(ana_ind$names=='SJ'),,]
#      aind$data[which(indices=='Z100'),,] <- ana_ind$data[which(ana_ind$names=='z100'),,]
#      aind$data[which(indices=='Z300'),,] <- ana_ind$data[which(ana_ind$names=='z300'),,]
#      aind$data[which(indices=='PWC'),,] <- ana_ind$data[which(ana_ind$names=='PWC'),,]
      
      #     aind <- list(ensmean=Hind%*%analysis$ensmean, names=indices)
      #     aind$ensmean[which(indices=='DIMI'),] <- ana_ind$ensmean[which(ana_ind$names=='DIMI'),]
      #     aind$ensmean[which(indices=='HC'),] <- ana_ind$ensmean[which(ana_ind$names=='HC'),]
      #     aind$ensmean[which(indices=='SJ'),] <- ana_ind$ensmean[which(ana_ind$names=='SJ'),]
      #     aind$ensmean[which(indices=='Z100'),] <- ana_ind$ensmean[which(ana_ind$names=='z100'),]
      #     aind$ensmean[which(indices=='Z300'),] <- ana_ind$ensmean[which(ana_ind$names=='z300'),]
      #     aind$ensmean[which(indices=='PWC'),] <- ana_ind$ensmean[which(ana_ind$names=='PWC'),]  
      # # data for ensemble member
      #     aind$data <- array(Hind %*% array(analysis$data, c(nrow(analysis$data), 
      #                        length(analysis$data)/nrow(analysis$data))), c(nrow(Hind), 
      #                        dim(analysis$data)[2:3]))
      #     aind$data[which(indices=='DIMI'),,] <- ana_ind$data[which(ana_ind$names=='DIMI'),,]
      #     aind$data[which(indices=='HC'),,] <- ana_ind$data[which(ana_ind$names=='HC'),,]
      #     aind$data[which(indices=='SJ'),,] <- ana_ind$data[which(ana_ind$names=='SJ'),,]
      #     aind$data[which(indices=='Z100'),,] <- ana_ind$data[which(ana_ind$names=='z100'),,]
      #     aind$data[which(indices=='Z300'),,] <- ana_ind$data[which(ana_ind$names=='z300'),,]
      #     aind$data[which(indices=='PWC'),,] <- ana_ind$data[which(ana_ind$names=='PWC'),,]
      
      
      
      
      
      
      
      if (vali) { 
        H.giorgi <- compute_giorgi_H_v2(giorgi, validate) #, numvar=3) # 3 vars temp, precip, slp
        indices <- c('ENH.temp2','NAM.temp2','SAM.temp2','AFR.temp2',
                     'ASI.temp2','AUS.temp2','ARC.temp2','ANT.temp2',
                     'NEU.temp2','MED.temp2',
                     'GLO.temp','NAM.precip','SAM.precip','AFR.precip',
                     'ASI.precip','AUS.precip','ARC.precip','ANT.precip',
                     'NEU.precip','MED.precip',
                     'SH.temp2','NAM.slp','SAM.slp','AFR.slp',
                     'ASI.slp','AUS.slp','ARC.slp','ANT.slp',
                     'NEU.slp','MED.slp',
                     'NH.temp2', 'EU.temp2', 'EU.precip', 'EU.slp') #, 
                     #                   'HC.calc', 'ITCZ.calc', 'SJ.calc', 'PV.calc', 
                     #                   'PWC.calc', 'DIMI.calc', 'NAO.calc', 'PNA.calc',
#                     'HC', 'SJ', 'Z100', 'Z300', 'PWC')
        Hind <- matrix(0,nrow=length(indices),ncol=nrow(validate$data))
        Hind[1,which(validate$names=="temp2")] <- H.giorgi[which(giorgi.short == 'ENH'),
                                                           which(validate$names=="temp2")]
        Hind[2,which(validate$names=="temp2")] <- H.giorgi[which(giorgi.short == 'NAM'),
                                                           which(validate$names=="temp2")]
        Hind[3,which(validate$names=="temp2")] <- H.giorgi[which(giorgi.short == 'SAM'),
                                                           which(validate$names=="temp2")]
        Hind[4,which(validate$names=="temp2")] <- H.giorgi[which(giorgi.short == 'AFR'),
                                                           which(validate$names=="temp2")]
        Hind[5,which(validate$names=="temp2")] <- H.giorgi[which(giorgi.short == 'ASI'),
                                                           which(validate$names=="temp2")]
        Hind[6,which(validate$names=="temp2")] <- H.giorgi[which(giorgi.short == 'AUS'),
                                                           which(validate$names=="temp2")]
        Hind[7,which(validate$names=="temp2")] <- H.giorgi[which(giorgi.short == 'ARC'),
                                                           which(validate$names=="temp2")]
        Hind[8,which(validate$names=="temp2")] <- H.giorgi[which(giorgi.short == 'ANT'),
                                                           which(validate$names=="temp2")]
        Hind[9,which(validate$names=="temp2")] <- H.giorgi[which(giorgi.short == 'NEU'),
                                                           which(validate$names=="temp2")]
        Hind[10,which(validate$names=="temp2")] <- H.giorgi[which(giorgi.short == 'MED'),
                                                            which(validate$names=="temp2")]
        Hind[11,which(validate$names=="temp2")] <- H.giorgi[which(giorgi.short == 'GLO'),
                                                            which(validate$names=="temp2")]
        Hind[12,which(validate$names=="precip")] <- H.giorgi[which(giorgi.short == 'NAM'),
                                                             which(validate$names=="precip")]
        Hind[13,which(validate$names=="precip")] <- H.giorgi[which(giorgi.short == 'SAM'),
                                                             which(validate$names=="precip")]
        Hind[14,which(validate$names=="precip")] <- H.giorgi[which(giorgi.short == 'AFR'),
                                                             which(validate$names=="precip")]
        Hind[15,which(validate$names=="precip")] <- H.giorgi[which(giorgi.short == 'ASI'),
                                                             which(validate$names=="precip")]
        Hind[16,which(validate$names=="precip")] <- H.giorgi[which(giorgi.short == 'AUS'),
                                                             which(validate$names=="precip")]
        Hind[17,which(validate$names=="precip")] <- H.giorgi[which(giorgi.short == 'ARC'),
                                                             which(validate$names=="precip")]
        Hind[18,which(validate$names=="precip")] <- H.giorgi[which(giorgi.short == 'ANT'),
                                                             which(validate$names=="precip")]
        Hind[19,which(validate$names=="precip")] <- H.giorgi[which(giorgi.short == 'NEU'),
                                                             which(validate$names=="precip")]
        Hind[20,which(validate$names=="precip")] <- H.giorgi[which(giorgi.short == 'MED'),
                                                             which(validate$names=="precip")]
        Hind[21,which(validate$names=="temp2")] <- H.giorgi[which(giorgi.short == 'SH'),
                                                            which(validate$names=="temp2")]
        Hind[22,which(validate$names=="slp")] <- H.giorgi[which(giorgi.short == 'NAM'),
                                                          which(validate$names=="slp")]
        Hind[23,which(validate$names=="slp")] <- H.giorgi[which(giorgi.short == 'SAM'),
                                                          which(validate$names=="slp")]
        Hind[24,which(validate$names=="slp")] <- H.giorgi[which(giorgi.short == 'AFR'),
                                                          which(validate$names=="slp")]
        Hind[25,which(validate$names=="slp")] <- H.giorgi[which(giorgi.short == 'ASI'),
                                                          which(validate$names=="slp")]
        Hind[26,which(validate$names=="slp")] <- H.giorgi[which(giorgi.short == 'AUS'),
                                                          which(validate$names=="slp")]
        Hind[27,which(validate$names=="slp")] <- H.giorgi[which(giorgi.short == 'ARC'),
                                                          which(validate$names=="slp")]
        Hind[28,which(validate$names=="slp")] <- H.giorgi[which(giorgi.short == 'ANT'),
                                                          which(validate$names=="slp")]
        Hind[29,which(validate$names=="slp")] <- H.giorgi[which(giorgi.short == 'NEU'),
                                                          which(validate$names=="slp")]
        Hind[30,which(validate$names=="slp")] <- H.giorgi[which(giorgi.short == 'MED'),
                                                          which(validate$names=="slp")]       
        Hind[31,which(validate$names=="temp2")] <- H.giorgi[which(giorgi.short == 'NH'),
                                                            which(validate$names=="temp2")]
        Hind[32,which(validate$names=="temp2")] <- H.giorgi[which(giorgi.short == 'EU'),
                                                            which(validate$names=="temp2")]
        Hind[33,which(validate$names=="precip")] <- H.giorgi[which(giorgi.short == 'EU'),
                                                             which(validate$names=="precip")]
        Hind[34,which(validate$names=="slp")] <- H.giorgi[which(giorgi.short == 'EU'),
                                                          which(validate$names=="slp")] 
        #      # land sea mask to take care of missing validata in the ocean      
        vpos <- which(!is.na(validate$data[,1]))
        Hind2 <- Hind[,vpos]
        validatanona <- validate$data[vpos,]
        vind <- list(data=Hind2 %*% validatanona, names=indices)
        
        try(vind$data[35,] <- vali_ind$data[which(vali_ind$names=="ind_rec_hc"),])
        try(vind$data[36,] <- vali_ind$data[which(vali_ind$names=="ind_rec_sj"),])
        try(vind$data[37,] <- vali_ind$data[which(vali_ind$names=="ind_rec_z100"),])
        try(vind$data[38,] <- vali_ind$data[which(vali_ind$names=="ind_rec_z300"),])
        try(vind$data[39,] <- vali_ind$data[which(vali_ind$names=="ind_rec_pwc"),])
      }
      #  }  # end sixmonstatevector loop
      
      # save annual file
      if (vali) {
            save(aind,eind,vind, file=paste0('../data/indices/',expname,'/indices',filenameext,cyr,'.Rdata'))
      } else {
            save(aind,eind, file=paste0('../data/indices/',expname,'/indices',filenameext,cyr,'.Rdata'))
      }
    }    # end year loop
  }    # end save prepplot







  if (mergetime) {
    for (cyr in syr:eyr){
      print(cyr)
      ptm1 <- proc.time()
#      if ((cyr < 1751) | (cyr > 1979)) {
        vali=F                 # switch off prepplot if no vali data selected
#      } else {
#        vali=T
#      }
      
      load(file=paste0('../data/indices/',expname,'/indices',filenameext,cyr,'.Rdata'))
      # merge timesteps 
      
      if (vali & cyr == 1751) {
        vind.allts=vind
        vind.allts$time=cyr
      }
      if (cyr == syr) {
        aind.allts=aind
        eind.allts=eind
        aind.allts$time=cyr
        eind.allts$time=cyr
      } else {
        aind.allts$data=abind(aind.allts$data,aind$data,along=2)
        aind.allts$ensmean=cbind(aind.allts$ensmean,aind$ensmean)
        aind.allts$bem=cbind(aind.allts$bem,aind$bem)
        aind.allts$time=c(aind.allts$time,cyr)
        if (vali) {
          vind.allts$data=cbind(vind.allts$data,vind$data)
          vind.allts$ensmean=cbind(vind.allts$ensmean,vind$ensmean)
          if (cyr != 1751) {
            vind.allts$time=c(vind.allts$time,cyr)
          }
        }
        eind.allts$data=abind(eind.allts$data,eind$data,along=2)
        eind.allts$ensmean=cbind(eind.allts$ensmean,eind$ensmean)
        eind.allts$bem=cbind(eind.allts$bem,eind$bem)
        eind.allts$time=c(eind.allts$time,cyr)
      } 
    }
    
    # save file
    if (vali) {
      save(aind.allts,eind.allts,vind.allts, 
           file=paste0('../data/indices/',expname,'/indices',filenameext,'allts.Rdata')) 
    } else {
      save(aind.allts,eind.allts, 
           file=paste0('../data/indices/',expname,'/indices',filenameext,'allts.Rdata'))
    }
  } #end mergetime
  
  








if (load_prepplot) {
  load(paste0('../data/indices/',expname,'/indices',filenameext,'allts.Rdata'))

  # calc global mean temp spread/sd pre post 1800
  e.glomean <- array(eind.allts$data[1,,],c(2,(dim(eind.allts$data)[2]/2),
                                            dim(eind.allts$data)[3]))
  a.glomean <- array(aind.allts$data[1,,],c(2,(dim(aind.allts$data)[2]/2),
    dim(aind.allts$data)[3]))
  if (every2grid) {
    load(file="../data/cru4_ens_sd_2ndgrid.Rdata")
  } else {
    load(file="../data/cru4_ens_sd.Rdata")  
  }
  obs.spread.temp <- cbind(as.vector(cru4_oct_apr$data),as.vector(cru4_may_sep$data))
  apply(obs.spread.temp,2,mean,na.rm=T)
  # winter
  mean(apply(e.glomean[1,,],1,sd)[1:200])
  mean(apply(a.glomean[1,,],1,sd)[1:200])
  mean(apply(a.glomean[1,,],1,sd)[201:400])
  #sommer
  mean(apply(e.glomean[2,,],1,sd)[1:200])
  mean(apply(a.glomean[2,,],1,sd)[1:200])
  mean(apply(a.glomean[2,,],1,sd)[201:400])
  #pdf(file='../figures/nat_data_paper/glo_mean_ana_ens_spread.pdf',width=9,height=3.5)
  pdf(file='../figures/',expname,'/glo_mean_ana_ens_spread.pdf',width=9,height=3.5)
    par(mfrow=c(1,2))
   #winter
    plot(aind.allts$time,apply(e.glomean[1,,],1,sd),col=rgb(0,0,0,10,maxColorValue=10),ty='l',
         ylim=c(0.05,0.5),ylab="Ens. std. dev.",xlab='Year',main='Oct. - Apr.')
    lines(aind.allts$time,apply(a.glomean[1,,],1,sd),ty='l',col=rgb(10,0,0,7,maxColorValue=10))
  # summer  
    plot(aind.allts$time,apply(e.glomean[1,,],1,sd),col=rgb(0,0,0,10,maxColorValue=10),ty='l',
         ylim=c(0.05,0.5),ylab="Ens. std. dev.",xlab='Year',main='May - Sep.')
    lines(aind.allts$time,apply(a.glomean[1,,],1,sd),ty='l',col=rgb(10,0,0,7,maxColorValue=10))
  dev.off()
  
  # plot global mean temp
#  pdf(file='../figures/nat_data_paper/eNH_summer_temp_recon4.pdf',width=10,height=4)
  pdf(file='../figures/',expname,'/eNH_summer_temp_recon4.pdf',width=10,height=4)
  wlen <- 1
  linew=1
  linew1=1
  linew2=2
  par(mfrow=c(1,1))
  par(mai=c(0.5,0.8,0.4,0.2))
  fs <- 1.0
  scal=F
  plbem=F
  plmem=F
  attributes(giorgi)
  anomper=(1901:1980)

  noaa <- read.table('../comparison_data/NOAA_PCN_1600_onward_stefan.txt',skip=3,header=F,sep='\t',
                     stringsAsFactors=FALSE)
  colnames(noaa) <- read.table('../comparison_data/NOAA_PCN_1600_onward_stefan.txt',nrow=1,sep='\t',
                               stringsAsFactors=FALSE)
  data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='ENH.temp2'),seas='sum',anomper=(1901:1980),
                     lw=linew,sca=scal,xa='s',plotbem=F,plotmem=T,plotech=F,
                     title='Extratropical Northern Hemisphere Summer Temperature Anomalies')
  ts_e <- ts(data$allind[,"eindmean"],start=data$period[1])
  ts_a <- ts(data$allind[,"aindmean"],start=data$period[1])

  # northern hemisphere extra-tropics (>20N)
  nc=open.nc('../comparison_data/CRUTEM4_A-Smean_ENHmean.nc', write=F)
  print.nc(nc)
  tmp1 <- var.get.nc(nc, "TEMPERATURE_ANOMALY") # for CRU temp
  tlist <- var.get.nc(nc, "TIME")
  unitstr <- "days since 1850-01-15 00:00:00"
  yrv <- utcal.nc(unitstr, tlist,type="n")[,1]
  pos <- yrv %in% data[['period']]
  tmp1 <- tmp1[pos]
  yrv <- yrv[pos]
  anopos <- yrv %in% anomper
  centerfac <- mean(tmp1[anopos])
  scalefaccru <- sd(tmp1[anopos])
  ts_i <- ts(tmp1,yrv[1])
  if (wlen > 1) {
    tmp1 <- runmean(scale(tmp1,center=centerfac,scale=F),wlen)
  } else {
    tmp1 <- scale(tmp1,center=centerfac,scale=F)
  }
  lines(yrv[30:length(yrv)],tmp1[30:length(tmp1)],ty='l',lwd=linew,lty=1,
        col=rgb(0,0,10,7,maxColorValue=10), ylab="",xlab="")
  
  cr <- read.table('../comparison_data/crowley_NH_temp_recon-2014.txt',sep='\t',header=T)
  tmp2 <- cr[,2]
  yrv <- cr[,1]
  pos <- yrv %in% data[['period']]
  tmp2 <- tmp2[pos]
  yrv <- yrv[pos]
  pos <- data[['period']] %in% yrv 
  cordata <- data[['allind']][pos,]
  data[['allind']][,1]
  anopos <- yrv %in% anomper
  centerfac <- mean(tmp2[anopos])
  ts_crow <- ts(tmp2,yrv[1])
  if (wlen > 1) {
    tmp2 <- runmean(scale(tmp2,center=centerfac,scale=F),wlen)
  } else {
    tmp2 <- scale(tmp2,center=centerfac,scale=F)
  }
  lines(yrv,tmp2,ty='l',lwd=linew,lty=1,col=rgb(0,5,5,10,maxColorValue=10), ylab="",xlab="")

  tmp4 <- noaa[,which(colnames(noaa)=="darrigo2006b")]
  yrv <- noaa[,1]
  pos <- yrv %in% data[['period']]
  tmp4 <- tmp4[pos]
  yrv <- yrv[pos]
  pos <- data[['period']] %in% yrv 
  cordata4 <- data[['allind']][pos,]
  anopos <- yrv %in% anomper
  centerfac <- mean(tmp4[anopos])
  ts_darr <- ts(tmp4,yrv[1])
  if (wlen > 1) {
    tmp4 <- runmean(scale(tmp4,center=centerfac,scale=F),wlen)
  } else {
    tmp4 <- scale(tmp4,center=centerfac,scale=F)
  }
  lines(yrv,tmp4,ty='l',lwd=linew,lty=1,col=rgb(0,10,10,5,maxColorValue=10), ylab="",xlab="")
  legend('topleft', c(
         #'ECHAM ens. mean','ECHAM spread',
         'Analysis ens. mean','Analysis spread',
         'CRUTEM4','Crowley et al. 2014','D Arrigo et al. 2006'),
         lwd=c(
          # 2,12,
          2,12,2,2,2),lty=c(1,1,1,1,1,1,1), col=c(
          #'black','lightgrey',
         'red','indianred1','blue','darkcyan','cyan'),cex=fs, bty='o', bg='white',
         box.col='white')
  
  # mark warmest/coldest decade before 1800
#  polygon(c(1790,1800,1800,1790),c(-2,-2,-1.9,-1.9),density=NA,col=rgb(5,0,0,3,maxColorValue=10))
#  polygon(c(1830,1840,1840,1830),c(-2,-2,-1.9,-1.9),density=NA,col=rgb(0,0,5,3,maxColorValue=10))
#  polygon(c(1640,1650,1650,1640),c(-2,-2,-1.9,-1.9),density=NA,col=rgb(0,0,5,3,maxColorValue=10))
# calc corr, RE between series
  ts_all <- ts.union(ts_e,ts_a,ts_i,ts_crow,ts_darr)
  cor(ts_all,use="pairwise.complete.obs",method='spearman')
  cor(ts_all,use="pairwise.complete.obs",method='pearson')
#  ts_all_detr <- detrend(na.omit(ts_all))) #,bp=150)
  # look at high freq. changes only:
  ts_all_detr <- ts_all-(filter(ts_all,rep((1/11),11),method="convolution",sides=2))
  cor(ts_all_detr,use="pairwise.complete.obs",method='spearman')
  cor(ts_all_detr,use="pairwise.complete.obs",method='pearson')
#  1- sum((ts_all[,2]-ts_all[,1])**2,na.rm=T)  / sum((ts_all[,1]-mean(ts_all[,1]))**2,na.rm=T)
  library("dcv")
  nona_i <- window(ts_all,1850,2004)
  nona_p <- window(ts_all,1782,1984)
  nona_pi <- window(ts_all,1850,1984)
  # D'Arrigo 1608-1990
#  test.RE(nona_i[,2],nona_i[,1])
#  test.RE(nona_p[,3],nona_p[,1])
#  test.RE(nona_pi[,2],nona_pi[,3])
  1-((sum(nona_i[,2]-nona_i[,3])^2)/(sum(nona_i[,1]-nona_i[,3])^2)) # RE CRU
  1-((sum(nona_p[,2]-nona_p[,4])^2)/(sum(nona_p[,1]-nona_p[,4])^2)) # RE Crowley
  1-((sum(nona_pi[,2]-nona_pi[,4])^2)/(sum(nona_pi[,1]-nona_pi[,4])^2)) # RE Crowley
  1-((sum(nona_p[,2]-nona_p[,5])^2)/(sum(nona_p[,1]-nona_p[,4])^2)) # RE D'Arrigo
  1-((sum(nona_pi[,2]-nona_pi[,5])^2)/(sum(nona_pi[,1]-nona_pi[,4])^2)) # RE D'Arrigo
dev.off() # end Fig. 1 current version with 1 panel
    
     
  
  
  
  
  
  
  
  
  
# plot difference between seasons, ENH, NH, Glob
  reg3 <- which(eind.allts$names=="GLO.temp2")
  reg2 <- which(eind.allts$names=="NH.temp2")
  reg <- which(eind.allts$names=="ENH.temp2")
  period <- eind.allts$time
  anomper=(1901:1980)
# annual mean
  years <- rep(eind.allts$time,each=2)  
  eindmeanann <- aggregate(eind.allts$ensmean[reg,],list(years),mean,na.rm=T)[,2]
  aindmeanann <- aggregate(aind.allts$ensmean[reg,],list(years),mean,na.rm=T)[,2]
  eindmeanannNH <- aggregate(eind.allts$ensmean[reg2,],list(years),mean,na.rm=T)[,2]
  aindmeanannNH <- aggregate(aind.allts$ensmean[reg2,],list(years),mean,na.rm=T)[,2]
  eindmeanannGL <- aggregate(eind.allts$ensmean[reg3,],list(years),mean,na.rm=T)[,2]
  aindmeanannGL <- aggregate(aind.allts$ensmean[reg3,],list(years),mean,na.rm=T)[,2]
  # summer
  somcol <- is.even(seq(1,ncol(aind.allts$ensmean)))  
  somcolv <- is.even(seq(1,ncol(vind.allts$data)))    
  eindmeansum <- eind.allts$ensmean[reg,somcol]  
  aindmeansum <- aind.allts$ensmean[reg,somcol]  
# winter
  wincol <- is.odd(seq(1,ncol(aind.allts$ensmean)))
  wincolv <- is.odd(seq(1,ncol(vind.allts$data)))
  eindmeanwin <- eind.allts$ensmean[reg,wincol]  
  aindmeanwin <- aind.allts$ensmean[reg,wincol]  
# join all indices
  allind <- cbind(eindmeanannGL,aindmeanannGL,eindmeanannNH,aindmeanannNH,
              eindmeanann,aindmeanann,eindmeansum,aindmeansum,eindmeanwin,aindmeanwin)
  anopos <- period %in% anomper
  centerfac <- apply(allind[anopos,],2,mean)
  centerfac[2] <- centerfac[1] 
  centerfac[3] <- centerfac[4]
  centerfac[5] <- centerfac[6]
  allind <- scale(allind,center=centerfac,scale=F)  
#  pdf(file=paste0('../figures/indices/summer-winter_temp_comp_EKF_ENH',filenameext,'.pdf'),
  pdf(file=paste0('../figures/',expname,'/summer-winter_temp_comp_EKF_ENH',filenameext,'.pdf'),      
      width=10,height=10)
    plot(period[],allind[,2],ty='l',col=rgb(0,0,0,0,maxColorValue=10), ylab="Temperature",
       ylim=c(min(allind[,2],na.rm=T)-0.5,max(allind[,2],na.rm=T)+0.5),
       main='',xaxt='s',bty='n')
    lines(period[],allind[,6],ty='l',lwd=1,lty=1,col=rgb(0,0,0,10,maxColorValue=10), ylab="",xlab="")
    lines(period[],allind[,8],ty='l',lwd=1,lty=1,col=rgb(10,0,0,8,maxColorValue=10), ylab="",xlab="")
    lines(period[],allind[,10],ty='l',lwd=1,lty=1,col=rgb(0,0,10,8,maxColorValue=10), ylab="",xlab="")
    cor(allind[,c(6,8,10)])
    mean(allind[1:100,7]-allind[1:100,9]) # sum vs win diff
    mean(allind[1:100,5]-allind[1:100,7]) # ann vs sum
#  mean(allind[1:100,8]-allind[1:100,10]) # sum vs win diff
#  mean(allind[1:100,6]-allind[1:100,8]) # ann vs sum
    legend('topleft', c('EKF400 ens. mean annual mean','EKF400 ens. mean Winter',
      'EKF400 ens. mean Summer'),lwd=c(2,2,2),lty=c(1,1,1), 
      col=c('black','blue','red'),cex=fs, bty='o', bg='white',box.col='white')
  dev.off()
    
# find and mark coldest/warmest decade
    d <- aind.allts$ensmean[(which(eind.allts$names=="ENH.temp2")),
                     (is.even(seq(1,ncol(aind.allts$ensmean))))]
    dec <- rep(NA,250)
    dectime <- rep(NA,250)
    for (ti in 1:length(dec)) {
      if (ti==1) {
        dec[ti] <- mean(d[1:8])
        dectime[ti] <- 1603
        i=2
      } else {
        dec[ti] <- mean(d[i:(i+7)])
        dectime[ti] <- dectime[ti-1] + 1 
        i=i+1
      }
    } 
#  plot(dectime,dec,ty='l',col='red')
    tmp=cbind(dec,dectime)
    tmp[order(dec),]
  
# mark warmest/coldest decade before 1800
# depends on definition, based on abs. values or anomalies?    
#    polygon(c(1790,1800,1800,1790),c(0.9,0.9,1,1),density=NA,col=rgb(5,0,0,3,maxColorValue=10))
#    polygon(c(1640,1650,1650,1640),c(0.9,0.9,1,1),density=NA,col=rgb(0,0,5,3,maxColorValue=10))

#   runm <- runmean(aind.allts$ensmean[(which(eind.allts$names=="ENH.temp2")),
#                      (is.even(seq(1,ncol(aind.allts$ensmean))))],11)
#   sortrunm <- sort(runm)
#   aind.allts$time[which(runm==sortrunm[1])]
#   aind.allts$time[which(runm==sortrunm[2])]
#   aind.allts$time[which(runm==sortrunm[3])]
 
#   aind.allts$time[which(runm==sortrunm[length(sortrunm)])]
#   aind.allts$time[which(runm==sortrunm[length(sortrunm)-1])]
#   aind.allts$time[which(runm==sortrunm[length(sortrunm)-2])]
  dev.off()

  
  
#   
#   # ADD NAO plots + reconstructions + skill (corr and RE, ...?)
#   load('../data/prepplot_season/indices_allts.Rdata') # NAO and PNA only in seasonal! why?
#   # NAO
#   data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='NAO.calc'),seas='win',anomper=(1901:1980),
#                      lw=linew,sca=T,xa='s',yl='',plotbem=F,plotmem=F,pag=F,plotech=T,title='NAO')
#   # some ocean based indices only NOT in landonly version
#   stefind <- read.table("/Users/joerg/Documents/climdata/indices/stefan/stefan_seasonal_indices.txt",
#                         header=T,skip=1,na.string=-9999)[,c(1,seq(2,5),seq(18,21),seq(34,37),seq(50,53),
#                         seq(66,69),seq(82,85),seq(98,101),seq(114,117),seq(130,133),seq(146,149))]
#   head <- read.table("/Users/joerg/Documents/climdata/indices/stefan/stefan_seasonal_indices.txt",
#                      nrow=2)[,c(1,seq(2,5),seq(18,21),seq(34,37),seq(50,53),seq(66,69),seq(82,85),
#                                 seq(98,101),seq(114,117),seq(130,133),seq(146,149))]
#   newhead <- rep(NA,ncol(head))
#   for (i in 1:ncol(head)) {
#     newhead[i] <- paste(head[1,i],head[2,i],sep='_')
#   }
#   colnames(stefind) <- newhead
#   vdata <- cbind(stefind['NAO_20C'],stefind['NAO_REC'],stefind['NAO_NCEP'],stefind['NAO_ERA-40'])
#   yrv <- stefind[,1]
#   centerfac <- apply(vdata,2,mean,na.rm=T)
#   scalefac <- apply(vdata,2,sd,na.rm=T)
#   if (wlen > 1) {
#     vdata <- runmean(scale(vdata,center=centerfac,scale=scalefac),wlen)
#   } else {
#     vdata <- scale(vdata,center=centerfac,scale=scalefac)
#   }
#   lines(yrv,vdata[,1],ty='l',lwd=linew,lty=1,col=rgb(0,0,10,8,maxColorValue=10), ylab="",xlab="")
# #  lines(yrv,vdata[,2],ty='l',lwd=linew,lty=1,col=rgb(10,0,10,10,maxColorValue=10), ylab="",xlab="")
#   lines(yrv,vdata[,3],ty='l',lwd=linew,lty=1,col=rgb(0,0,5,8,maxColorValue=10), ylab="",xlab="")
# #  lines(yrv,vdata[,4],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,10,maxColorValue=10), ylab="",xlab="")
#   ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],
#                         use="pairwise.complete.obs"),digits=2) 
#   acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],
#                         use="pairwise.complete.obs"),digits=2) 
#   ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],
#                        use="pairwise.complete.obs"),digits=2) 
#   acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],
#                        use="pairwise.complete.obs"),digits=2)
#   legend('bottomleft', paste('ECHAM(1,2)/EKF(3,4) ens. mean - 20CR/NNR:',
#     ecormean[1],ecormean[3],acormean[1],acormean[3]),
# #    paste('ECHAM/EKF BEM - 20CR:',ecorbem[1],acorbem[1])),
#     cex=fs, bty='o', bg='white', box.col='white')
#   
# #  legend('bottomleft', c(paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',
# #    ecormean[1],ecormean[2],ecormean[3],ecormean[4]),paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',
# #    acormean[1],acormean[2],acormean[3],acormean[4]),paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',
# #    ecorbem[1],ecorbem[2],ecorbem[3],ecorbem[4]),paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',
# #    acorbem[1],acorbem[2],acorbem[3],acorbem[4])),cex=fs, bty='o', bg='white', box.col='white')
#   
#   vf <- c(1600,1650,1660,1783,1809,1815,1835,1883,1886,1902,1912,1980,1991) 
#   nao_trou <- read.table("/Users/joerg/Documents/climdata/indices/www/nao-trouet2009.txt",skip=93,header=T)
#   nao_cook <- read.table("/Users/joerg/Documents/climdata/indices/www/nao_cook2002.txt",skip=50,header=F)
#   nao_lut <- read.table("/Users/joerg/Documents/climdata/indices/www/nao_sea_luterbacher_2002.txt",
#                         skip=19,header=T)
#   nao_lut2 <- read.table("/Users/joerg/Documents/climdata/indices/www/nao_mon_luterbacher_2002.txt",
#                         skip=27,header=T)
#   pna_trou <- read.table("/Users/joerg/Documents/climdata/indices/www/pna-trouet2009.txt",skip=75,header=F)
#   
#   # data has to be loaded just before "plot_pages"! why?
#   naol <- nao_lut[nao_lut[,2]=="wi",]
#   naol <- ts(naol[naol[,1]>=1600,3],start=1600)
#   agglist <- c(rep(1659,3),rep(seq(1660,2001),each=3))
#   naol2 <- nao_lut2[(nao_lut2[,2]=="Dec" | nao_lut2[,2]=="Jan" | nao_lut2[,2]=="Feb"),]
# #  naol2 <- nao_lut2[(nao_lut2[,2]=="Oct" | nao_lut2[,2]=="Nov" | nao_lut2[,2]=="Dec" | 
# #                     nao_lut2[,2]=="Jan" | nao_lut2[,2]=="Feb" | nao_lut2[,2]=="Mar"),]
# #  agglist <- c(rep(1659,4),rep(seq(1660,2001),each=6))
#   naol2agg <- aggregate(naol2,list(agglist),mean)
#   naol2 <- ts(naol2agg[,4],start=1659,freq=1)
#   naol3 <- ts(c(naol,naol2),start=1600)
#   naot <- ts(nao_trou[nao_trou[,1]>=1600,2],start=1600)
#   naoc <- ts(nao_cook[nao_cook[,1]>=1600,2],start=1600)
#   vdata <- ts.union(naol3,naot,naoc)
#   yrv <- seq(1600,1600+nrow(vdata)-1)
#   centerfac <- apply(vdata,2,mean,na.rm=T)
#   scalefac <- apply(vdata,2,sd,na.rm=T)
#   if (wlen > 1) {
#     vdata <- runmean(scale(vdata,center=centerfac,scale=scalefac),wlen)
#   } else {
#     vdata <- scale(vdata,center=centerfac,scale=scalefac)
#   }
#   lines(yrv,vdata[,1],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,8,maxColorValue=10), ylab="",xlab="")
# #  lines(yrv,vdata[,2],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,10,maxColorValue=10), ylab="",xlab="")
# #  lines(yrv,vdata[,3],ty='l',lwd=linew,lty=1,col=rgb(10,10,0,10,maxColorValue=10), ylab="",xlab="")
# #  lines(yrv,vdata[,4],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,10,maxColorValue=10), ylab="",xlab="")
# #  ecormean <- round(cor(data[['allind']][,1],vdata[,],use="pairwise.complete.obs"),digits=2) 
# #  acormean <- round(cor(data[['allind']][,2],vdata[,],use="pairwise.complete.obs"),digits=2) 
# #  ecorbem <- round(cor(data[['allind']][,3],vdata[,],use="pairwise.complete.obs"),digits=2) 
# #  acorbem <- round(cor(data[['allind']][,4],vdata[,],use="pairwise.complete.obs"),digits=2)
#   round(cor(vdata,use="pairwise.complete.obs"),digits=2)
#   round(cor(data[['allind']][,c(1,2)],vdata[,],use="pairwise.complete.obs"),digits=2)
# #  legend('topleft', c(paste(''Analysis ens. mean - Luterb., Trouet, Cook:',
# #    ecormean[1],ecormean[2],ecormean[3],ecormean[4]),paste('Analysis ens.mean - Lut.,Trouet,Cook:',
# #    acormean[1],acormean[2],acormean[3],acormean[4]),paste('ECHAM BEM - Lut.,Trouet,Cook:',
# #    ecorbem[1],ecorbem[2],ecorbem[3],ecorbem[4]),paste('Analysis ens. bem - Lut.,Trouet,Cook:',
# #    acorbem[1],acorbem[2],acorbem[3],acorbem[4])),cex=fs, bty='o', bg='white', box.col='white')
#   
#   
#   # PNA
#   data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='PNA.calc'),seas='win',anomper=(1901:1980),lw=linew,sca=T,xa='s',yl='')
#   #abline(v=vf,col='black',lty=2)
#   vdata <- cbind(stefind['PNA_20C'],stefind['PNA_REC'],stefind['PNA_NCEP'],stefind['PNA_ERA-40'])
#   yrv <- stefind[,1]
#   centerfac <- apply(vdata,2,mean,na.rm=T)
#   scalefac <- apply(vdata,2,sd,na.rm=T)
#   if (wlen > 1) {
#     vdata <- runmean(scale(vdata,center=centerfac,scale=scalefac),wlen)
#   } else {
#     vdata <- scale(vdata,center=centerfac,scale=scalefac)
#   }
#   lines(yrv,vdata[,1],ty='l',lwd=linew,lty=1,col=rgb(0,0,10,10,maxColorValue=10), ylab="",xlab="")
#   lines(yrv,vdata[,2],ty='l',lwd=linew,lty=1,col=rgb(10,0,10,10,maxColorValue=10), ylab="",xlab="")
#   lines(yrv,vdata[,3],ty='l',lwd=linew,lty=1,col=rgb(10,10,0,10,maxColorValue=10), ylab="",xlab="")
#   lines(yrv,vdata[,4],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,10,maxColorValue=10), ylab="",xlab="")
#   ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
#   acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
#   ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
#   acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
#   legend('bottomleft', c(paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',ecormean[1],ecormean[2],ecormean[3],ecormean[4]),paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',acormean[1],acormean[2],acormean[3],acormean[4]),paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],ecorbem[2],ecorbem[3],ecorbem[4]),paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],acorbem[2],acorbem[3],acorbem[4])),cex=fs, bty='o', bg='white', box.col='white')
#   
#   dev.off()
  
  
  
  

  
  
  
  
  
#  } else {
#    legend('bottomleft', c('ECHAM ens. mean','ECHAM spread','Analysis ens. mean','Analysis spread'),
#           lwd=c(2,12,2,12),lty=c(1,1,1,1), col=c('black','lightgrey',
#           'red','indianred1'),cex=fs, bty='o',bg='white', box.col='white')  
#  }
  #data <- plot_pages(wl=wlen,reg=1,seas='sum')
  #data <- plot_pages(wl=wlen,reg=1,seas='win')
  data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='ENH.temp2'),seas='sum',anomper=(1901:1980),
                     lw=linew,sca=scal,xa='s',plotbem=F,plotmem=T,plotech=F,yl='')  
#  par(new=T)
#  data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='ENH.temp2'),seas='sum',anomper=(1901:1980),
#                   lw=linew,sca=scal,xa='n',plotbem=F,plotmem=F,plotech=F,yl='')  
#  par(new=T)  
#  data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='ENH.temp2'),seas='win',anomper=(1901:1980),
#                   lw=linew,sca=scal,xa='n',plotbem=F,plotmem=F,plotech=F,yl='')  
  nc=open.nc('../comparison_data/cru3_tmp_25_annmean_enhmean.nc', write=F)
  print.nc(nc)
  tmp1 <- var.get.nc(nc, "TMP") # for CRU temp
  tlist <- var.get.nc(nc, "TIME")
  unitstr <- "months since 1901-01-15 00:00:00"
  yrv <- utcal.nc(unitstr, tlist,type="n")[,1]
  pos <- yrv %in% data[['period']]
  tmp1 <- tmp1[pos]
  yrv <- yrv[pos]
  anopos <- yrv %in% anomper
  centerfac <- mean(tmp1[anopos])
  scalefaccru <- sd(tmp1[anopos])
  if (wlen > 1) {
    tmp1 <- runmean(scale(tmp1,center=centerfac,scale=F),wlen)
  } else {
    tmp1 <- scale(tmp1,center=centerfac,scale=F)
  }
  lines(yrv,tmp1,ty='l',lwd=linew,lty=1,col=rgb(0,0,10,7,maxColorValue=10), ylab="",xlab="")
  
  cr <- read.table('../comparison_data/crowley_NH_temp_recon-2014.txt',sep='\t',header=T)
  tmp2 <- cr[,2]
  yrv <- cr[,1]
  pos <- yrv %in% data[['period']]
  tmp2 <- tmp2[pos]
  yrv <- yrv[pos]
  pos <- data[['period']] %in% yrv 
  cordata <- data[['allind']][pos,]
  data[['allind']][,1]
  anopos <- yrv %in% anomper
  centerfac <- mean(tmp2[anopos])
  if (wlen > 1) {
    tmp2 <- runmean(scale(tmp2,center=centerfac,scale=F),wlen)
  } else {
    tmp2 <- scale(tmp2,center=centerfac,scale=F)
  }
  lines(yrv,tmp2,ty='l',lwd=linew,lty=1,col=rgb(0,10,10,8,maxColorValue=10), ylab="",xlab="")
  
#   tmp3 <- noaa[,which(colnames(noaa)=="esper2002")]
#   yrv <- noaa[,1]
#   pos <- yrv %in% data[['period']]
#   tmp3 <- tmp3[pos]
#   yrv <- yrv[pos]
#   pos <- data[['period']] %in% yrv 
#   cordata3 <- data[['allind']][pos,]
#   anopos <- yrv %in% anomper
#   centerfac <- mean(tmp3[anopos])
#   if (wlen > 1) {
#     tmp3 <- runmean(scale(tmp3,center=centerfac,scale=F),wlen)
#   } else {
#     tmp3 <- scale(tmp3,center=centerfac,scale=(scalefaccru*2))
#   }
#   lines(yrv,tmp3,ty='l',lwd=linew,lty=1,col=rgb(10,0,10,5,maxColorValue=10), ylab="",xlab="")
  
  tmp4 <- noaa[,which(colnames(noaa)=="darrigo2006b")]
  yrv <- noaa[,1]
  pos <- yrv %in% data[['period']]
  tmp4 <- tmp4[pos]
  yrv <- yrv[pos]
  pos <- data[['period']] %in% yrv 
  cordata4 <- data[['allind']][pos,]
  anopos <- yrv %in% anomper
  centerfac <- mean(tmp4[anopos])
  if (wlen > 1) {
    tmp4 <- runmean(scale(tmp4,center=centerfac,scale=F),wlen)
  } else {
    tmp4 <- scale(tmp4,center=centerfac,scale=F)
  }
  lines(yrv,tmp4,ty='l',lwd=linew,lty=1,col=rgb(0,2,0,5,maxColorValue=10), ylab="",xlab="")
  
  # tmp5 <- noaa[,which(colnames(noaa)=="wilson2007")]
  # yrv <- noaa[,1]
  # pos <- yrv %in% data[['period']]
  # tmp5 <- tmp5[pos]
  # yrv <- yrv[pos]
  # anopos <- yrv %in% anomper
  # centerfac <- mean(tmp5[anopos])
  # if (wlen > 1) {
  #   tmp5 <- runmean(scale(tmp5,center=centerfac,scale=F),wlen)
  # } else {
  #   tmp5 <- scale(tmp5,center=centerfac,scale=F)
  # }
  # lines(yrv,tmp5,ty='l',lwd=linew,lty=1,col=rgb(0,10,0,5,maxColorValue=10), ylab="",xlab="")
  
  #ecormean <- round(cor(data[['allind']][,1],tmp2,use="pairwise.complete.obs"),digits=2) 
  #acormean <- round(cor(data[['allind']][,2],tmp2,use="pairwise.complete.obs"),digits=2) 
  ecormean <- round(cor(cordata[,1],tmp2,use="pairwise.complete.obs"),digits=2) 
  acormean <- round(cor(cordata[,2],tmp2,use="pairwise.complete.obs"),digits=2) 
  #ecormean3 <- round(cor(cordata3[,1],tmp3,use="pairwise.complete.obs"),digits=2) 
  #acormean3 <- round(cor(cordata3[,2],tmp3,use="pairwise.complete.obs"),digits=2) 
  ecormean4 <- round(cor(cordata4[,1],tmp4,use="pairwise.complete.obs"),digits=2) 
  acormean4 <- round(cor(cordata4[,2],tmp4,use="pairwise.complete.obs"),digits=2) 
  ecorcru <- round(cor(data[['allind']][299:length(data[['period']]),1],tmp1,use="pairwise.complete.obs"),digits=2) 
  acorcru <- round(cor(data[['allind']][299:length(data[['period']]),2],tmp1,use="pairwise.complete.obs"),digits=2) 
  legend('topleft', c('Correlation:',
         paste('ECHAM/Analysis ens. mean - Crowley 2014:',ecormean,'/',acormean), 
#         paste('ECHAM/Analysis ens. mean - Esper 2002:',ecormean3,'/',acormean3), 
         paste('ECHAM/Analysis ens. mean - DArrigo 2007:',ecormean4,'/',acormean4), 
         paste('ECHAM/Analysis ens. mean - CRU TS3:',ecorcru,'/',acorcru)), 
         cex=fs, bty='o', bg='white', box.col='white')                    
  legend('bottomright', c('Analysis ens. mean','Analysis spread','CRU TS3','Crowley 2014',
         'D Arrigo 2006'),lwd=c(1,12,2,2,2,2),lty=c(1,1,1,1,1,1), 
         col=c('red','indianred1','blue','cyan','darkgreen'),cex=fs, bty='o', bg='white', box.col='white')  
#legend('bottomright', c('Analysis ens. mean','Analysis spread','CRU TS3','Crowley 2014',
#                        'Esper 2002','D Arrigo 2006'),lwd=c(1,12,2,2,2,2),lty=c(1,1,1,1,1,1), 
#       col=c('red','indianred1','blue','cyan','violet','orange'),cex=fs, bty='o', bg='white',
#       box.col='white')  


  
  # Europe
  data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='EU.temp2'),lw=linew,sca=scal, yl='Temperature',xa='s',
                     plotbem=plbem,plotmem=plmem,plotech=F,seas='sum')
  reg <- which(vind.allts$names=='EU.temp2')
  years <- rep(vind.allts$time,each=2) 
  years <- years[(length(years)-155):length(years)]
  vdat <- vind.allts$data[reg,(ncol(vind.allts$data)-155):ncol(vind.allts$data)]
  vindmean <- aggregate(vdat,list(years),mean)[,2]
  mi <- apply(cbind(data[['pages']][,3],data[['pages']][,4]),1,min)
  ma <- apply(cbind(data[['pages']][,3],data[['pages']][,4]),1,max)
  polygon(c(data[['periodv']][-nrow(data[['pages']])],rev(data[['periodv']][-nrow(data[['pages']])])),
          c(mi[-nrow(data[['pages']])],rev(ma[-nrow(data[['pages']])])),
          density=NA,col=rgb(0,5,5,3,maxColorValue=10))
  #lines(data[['periodv']],data[['pages']][,1],ty='l',lwd=linew,lty=1,col=rgb(0,5,5,10,maxColorValue=10), 
  #      ylab="",xlab="")
  lines(data[['periodv']],data[['pages']][,2],ty='l',lwd=linew2,lty=1,col=rgb(0,5,5,10,maxColorValue=10), 
        ylab="",xlab="")
  lines(seq(1902,1979),scale(vindmean,center=T,scale=scal),ty='l', 
        lwd=linew2, lty=1, col=rgb(0,0,10,7,maxColorValue=10), ylab="",xlab="")
  lines(data[['periodv']],data[['allind']][,2],ty='l',lwd=linew,lty=1,col=rgb(10,0,0,8,maxColorValue=10),
        ylab="",xlab="")
  ecormean <- round(cor(data[['allind']][,1],data[['pages']][,2],use="pairwise.complete.obs"),digits=2) 
  acormean <- round(cor(data[['allind']][,2],data[['pages']][,2],use="pairwise.complete.obs"),digits=2) 
  ecorbem <- round(cor(data[['allind']][,3],data[['pages']][,2],use="pairwise.complete.obs"),digits=2) 
  acorbem <- round(cor(data[['allind']][,4],data[['pages']][,2],use="pairwise.complete.obs"),digits=2)
  ecorval <- round(cor(data[['allind']][300:377,1],vindmean,use="pairwise.complete.obs"),digits=2) 
  acorval <- round(cor(data[['allind']][300:377,2],vindmean,use="pairwise.complete.obs"),digits=2) 
  ecorvalbem <- round(cor(data[['allind']][300:377,3],vindmean,use="pairwise.complete.obs"),digits=2) 
  acorvalbem <- round(cor(data[['allind']][300:377,4],vindmean,use="pairwise.complete.obs"),digits=2) 
  legend('topleft', c('Correlation:',
         paste('ECHAM/Analysis ens. mean - PAGES recon.:',ecormean,'/',acormean), 
         paste('ECHAM/Analysis ens. mean - CRU TS3:',ecorcru,'/',acorcru)), 
         cex=fs, bty='o', bg='white', box.col='white')                    
#  legend('bottomright', c('Analysis ens. mean','PAGES','PAGES spread','CRU TS3'),
#         lwd=c(2,2,12,2,2),lty=c(1,1,1,1,1), 
#         col=c('red','lightcyan','darkcyan','blue'),cex=fs, bty='o', bg='white', box.col='white')  
#
#    legend('topleft', c(paste('ECHAM ens. mean - PAGES recon, CRU TS3:',ecormean,ecorval),
#                        paste('Analysis ens. mean - PAGES recon, CRU TS3:',acormean,acorval)),
#           cex=fs, bty='o', bg='white', box.col='white')   
  
  # North America
  data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='NAM.temp2'),lw=linew,sca=scal, yl='',xa='s',
                     plotbem=plbem,plotmem=plmem,plotech=F,seas='sum')
  reg <- which(vind.allts$names=='NAM.temp2')
  years <- rep(vind.allts$time,each=2) 
  years <- years[(length(years)-155):length(years)]
  vdat <- vind.allts$data[reg,(ncol(vind.allts$data)-155):ncol(vind.allts$data)]
  vindmean <- aggregate(vdat,list(years),mean)[,2]
  mi <- apply(cbind(data[['pages']][,3],data[['pages']][,4]),1,min)
  ma <- apply(cbind(data[['pages']][,3],data[['pages']][,4]),1,max)
  polygon(c(data[['periodv']][-nrow(data[['pages']])],rev(data[['periodv']][-nrow(data[['pages']])])),
          c(mi[-nrow(data[['pages']])],rev(ma[-nrow(data[['pages']])])),
          density=NA,col=rgb(0,5,5,3,maxColorValue=10))
  #lines(data[['periodv']],data[['pages']][,1],ty='l',lwd=linew,lty=1,col=rgb(0,5,5,10,maxColorValue=10), 
  #      ylab="",xlab="")
  lines(data[['periodv']],data[['pages']][,2],ty='l',lwd=linew2,lty=1,col=rgb(0,5,5,10,maxColorValue=10), 
        ylab="",xlab="")
  lines(seq(1902,1979),scale(vindmean,center=T,scale=scal),ty='l', 
        lwd=linew2, lty=1, col=rgb(0,0,10,7,maxColorValue=10), ylab="",xlab="")
  lines(data[['periodv']],data[['allind']][,2],ty='l',lwd=linew,lty=1,col=rgb(10,0,0,8,maxColorValue=10), 
        ylab="",xlab="")
  ecormean <- round(cor(data[['allind']][,1],data[['pages']][,2],use="pairwise.complete.obs"),digits=2) 
  acormean <- round(cor(data[['allind']][,2],data[['pages']][,2],use="pairwise.complete.obs"),digits=2) 
  ecorbem <- round(cor(data[['allind']][,3],data[['pages']][,2],use="pairwise.complete.obs"),digits=2) 
  acorbem <- round(cor(data[['allind']][,4],data[['pages']][,2],use="pairwise.complete.obs"),digits=2)
  ecorval <- round(cor(data[['allind']][300:377,1],vindmean,use="pairwise.complete.obs"),digits=2) 
  acorval <- round(cor(data[['allind']][300:377,2],vindmean,use="pairwise.complete.obs"),digits=2) 
  ecorvalbem <- round(cor(data[['allind']][300:377,3],vindmean,use="pairwise.complete.obs"),digits=2) 
  acorvalbem <- round(cor(data[['allind']][300:377,4],vindmean,use="pairwise.complete.obs"),digits=2) 
  legend('topleft', c('Correlation:',
         paste('ECHAM/Analysis ens. mean - PAGES recon.:',ecormean,'/',acormean), 
         paste('ECHAM/Analysis ens. mean - CRU TS3:',ecorcru,'/',acorcru)), 
         cex=fs, bty='o', bg='white', box.col='white')                    
  legend('bottomright', c('Analysis ens. mean','PAGES','PAGES spread','CRU TS3'),
         lwd=c(2,2,12,2,2),lty=c(1,1,1,1,1), 
         col=c('red','lightcyan3','darkcyan','blue'),cex=fs, bty='o', bg='white', box.col='white') 
  dev.off()
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
# continental temp comparison
pdf(file=paste0('../figures/',expname,'/pages_cont_temp_comp',filenameext,'.pdf'),width=10,height=10)
wlen <- 1
linew=1
linew1=1
linew2=2
par(mfrow=c(4,2))
par(mai=c(0.3,0.6,0.2,0))
fs <- 1.0
scal=F
plbem=F
plmem=F
attributes(giorgi)



# Arctic
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='ARC.temp2'),lw=linew,sca=scal, yl='',xa='n',plotbem=plbem)
reg <- which(vind.allts$names=='ARC.temp2')
years <- rep(vind.allts$time,each=2) 
years <- years[(length(years)-155):length(years)]
vdat <- vind.allts$data[reg,(ncol(vind.allts$data)-155):ncol(vind.allts$data)]
vindmean <- aggregate(vdat,list(years),mean)[,2]
lines(seq(1902,1979),scale(vindmean,center=T,scale=scal),ty='l', 
      lwd=linew, lty=1, col=rgb(0,0,10,5,maxColorValue=10), ylab="",xlab="")
lines(data[['periodv']],data[['pages']][,22],ty='l',lwd=linew2,lty=2,col=rgb(10,0,10,5,maxColorValue=10), 
      ylab="",xlab="")
lines(data[['periodv']],data[['pages']][,23],ty='l',lwd=linew2,lty=2,col=rgb(0,10,10,5,maxColorValue=10), 
      ylab="",xlab="")
ecormean <- round(cor(data[['allind']][,1],data[['pages']][,23],use="pairwise.complete.obs"),digits=2)
acormean <- round(cor(data[['allind']][,2],data[['pages']][,23],use="pairwise.complete.obs"),digits=2) 
ecorbem <- round(cor(data[['allind']][,3],data[['pages']][,23],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][,4],data[['pages']][,23],use="pairwise.complete.obs"),digits=2)
ecorval <- round(cor(data[['allind']][300:377,1],vindmean,use="pairwise.complete.obs"),digits=2) 
acorval <- round(cor(data[['allind']][300:377,2],vindmean,use="pairwise.complete.obs"),digits=2) 
ecorvalbem <- round(cor(data[['allind']][300:377,3],vindmean,use="pairwise.complete.obs"),digits=2) 
acorvalbem <- round(cor(data[['allind']][300:377,4],vindmean,use="pairwise.complete.obs"),digits=2) 
if (plbem) {
  legend('topleft', c(paste('ECHAM ens. mean - PAGES recon, CRU TS3:',ecormean,ecorval),
                    paste('Analysis ens. mean - PAGES recon, CRU TS3:',acormean,acorval),
                    paste('ECHAM BEM - PAGES recon, CRU TS3:',ecorbem,ecorvalbem),
                    paste('Analysis BEM - PAGES recon, CRU TS3:',acorbem,acorvalbem)), 
       cex=fs, bty='o', bg='white', box.col='white') 
} else {
  legend('topleft', c(paste('ECHAM ens. mean - PAGES recon, CRU TS3:',ecormean,ecorval),
                      paste('Analysis ens. mean - PAGES recon, CRU TS3:',acormean,acorval)),
         cex=fs, bty='o', bg='white', box.col='white') 
}

# Europe
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='EU.temp2'),lw=linew,sca=scal, yl='',xa='n',
                   plotbem=plbem,plotmem=plmem,plotech=F,seas='yrmean')
reg <- which(vind.allts$names=='EU.temp2')
years <- rep(vind.allts$time,each=2) 
years <- years[(length(years)-155):length(years)]
vdat <- vind.allts$data[reg,(ncol(vind.allts$data)-155):ncol(vind.allts$data)]
vindmean <- aggregate(vdat,list(years),mean)[,2]
mi <- apply(cbind(data[['pages']][,3],data[['pages']][,4]),1,min)
ma <- apply(cbind(data[['pages']][,3],data[['pages']][,4]),1,max)
polygon(c(data[['periodv']][-nrow(data[['pages']])],rev(data[['periodv']][-nrow(data[['pages']])])),
        c(mi[-nrow(data[['pages']])],rev(ma[-nrow(data[['pages']])])),
        density=NA,col=rgb(3,3,3,2,maxColorValue=10))
#lines(data[['periodv']],data[['pages']][,1],ty='l',lwd=linew,lty=1,col=rgb(0,5,5,10,maxColorValue=10), 
#      ylab="",xlab="")
lines(data[['periodv']],data[['pages']][,2],ty='l',lwd=linew2,lty=1,col=rgb(3,3,3,7,maxColorValue=10), 
      ylab="",xlab="")
lines(seq(1902,1979),scale(vindmean,center=T,scale=scal),ty='l', 
      lwd=linew2, lty=1, col=rgb(0,0,10,7,maxColorValue=10), ylab="",xlab="")
lines(period,data[['allind']][,2],ty='l',lwd=linew,lty=1,col=rgb(10,0,0,5,maxColorValue=10), ylab="",xlab="")
ecormean <- round(cor(data[['allind']][,1],data[['pages']][,2],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][,2],data[['pages']][,2],use="pairwise.complete.obs"),digits=2) 
ecorbem <- round(cor(data[['allind']][,3],data[['pages']][,2],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][,4],data[['pages']][,2],use="pairwise.complete.obs"),digits=2)
ecorval <- round(cor(data[['allind']][300:377,1],vindmean,use="pairwise.complete.obs"),digits=2) 
acorval <- round(cor(data[['allind']][300:377,2],vindmean,use="pairwise.complete.obs"),digits=2) 
ecorvalbem <- round(cor(data[['allind']][300:377,3],vindmean,use="pairwise.complete.obs"),digits=2) 
acorvalbem <- round(cor(data[['allind']][300:377,4],vindmean,use="pairwise.complete.obs"),digits=2) 
if (plbem) {
  legend('topleft', c(paste('ECHAM ens. mean - PAGES recon, CRU TS3:',ecormean,ecorval),
                    paste('Analysis ens. mean - PAGES recon, CRU TS3:',acormean,acorval),
                    paste('ECHAM BEM - PAGES recon, CRU TS3:',ecorbem,ecorvalbem),
                    paste('Analysis BEM - PAGES recon, CRU TS3:',acorbem,acorvalbem)), 
       cex=fs, bty='o', bg='white', box.col='white') 
} else {
  legend('topleft', c(paste('ECHAM ens. mean - PAGES recon, CRU TS3:',ecormean,ecorval),
                      paste('Analysis ens. mean - PAGES recon, CRU TS3:',acormean,acorval)),
         cex=fs, bty='o', bg='white', box.col='white') 
}

# North America
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='NAM.temp2'),lw=linew,sca=scal, yl='',xa='n',plotbem=plbem)
reg <- which(vind.allts$names=='NAM.temp2')
years <- rep(vind.allts$time,each=2) 
years <- years[(length(years)-155):length(years)]
vdat <- vind.allts$data[reg,(ncol(vind.allts$data)-155):ncol(vind.allts$data)]
vindmean <- aggregate(vdat,list(years),mean)[,2]
lines(seq(1902,1979),scale(vindmean,center=T,scale=scal),ty='l', 
      lwd=linew, lty=1, col=rgb(0,0,10,5,maxColorValue=10), ylab="",xlab="")
lines(data[['periodv']],data[['pages']][,6],ty='l',lwd=linew2,lty=2,col=rgb(10,0,10,5,maxColorValue=10), 
      ylab="",xlab="")
lines(data[['periodv']],data[['pages']][,7],ty='l',lwd=linew2,lty=2,col=rgb(0,10,10,5,maxColorValue=10), 
      ylab="",xlab="")
ecormean <- round(cor(data[['allind']][,1],data[['pages']][,7],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][,2],data[['pages']][,7],use="pairwise.complete.obs"),digits=2) 
ecorbem <- round(cor(data[['allind']][,3],data[['pages']][,7],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][,4],data[['pages']][,7],use="pairwise.complete.obs"),digits=2)
ecorval <- round(cor(data[['allind']][300:377,1],vindmean,use="pairwise.complete.obs"),digits=2) 
acorval <- round(cor(data[['allind']][300:377,2],vindmean,use="pairwise.complete.obs"),digits=2) 
ecorvalbem <- round(cor(data[['allind']][300:377,3],vindmean,use="pairwise.complete.obs"),digits=2) 
acorvalbem <- round(cor(data[['allind']][300:377,4],vindmean,use="pairwise.complete.obs"),digits=2) 
if (plbem) {
  legend('topleft', c(paste('ECHAM ens. mean - PAGES recon, CRU TS3:',ecormean,ecorval),
                    paste('Analysis ens. mean - PAGES recon, CRU TS3:',acormean,acorval),
                    paste('ECHAM BEM - PAGES recon, CRU TS3:',ecorbem,ecorvalbem),
                    paste('Analysis BEM - PAGES recon, CRU TS3:',acorbem,acorvalbem)), 
       cex=fs, bty='o', bg='white', box.col='white') 
} else {
  legend('topleft', c(paste('ECHAM ens. mean - PAGES recon, CRU TS3:',ecormean,ecorval),
                      paste('Analysis ens. mean - PAGES recon, CRU TS3:',acormean,acorval)),
         cex=fs, bty='o', bg='white', box.col='white') 
}

# Asia
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='ASI.temp2'),lw=linew,sca=scal, yl='',xa='n',plotbem=plbem)
reg <- which(vind.allts$names=='ASI.temp2')
years <- rep(vind.allts$time,each=2) 
years <- years[(length(years)-155):length(years)]
vdat <- vind.allts$data[reg,(ncol(vind.allts$data)-155):ncol(vind.allts$data)]
vindmean <- aggregate(vdat,list(years),mean)[,2]
lines(seq(1902,1979),scale(vindmean,center=T,scale=scal),ty='l', 
      lwd=linew, lty=1, col=rgb(0,0,10,5,maxColorValue=10), ylab="",xlab="")
lines(data[['periodv']],data[['pages']][,26],ty='l',lwd=linew2,lty=2,col=rgb(10,0,10,5,maxColorValue=10), 
      ylab="",xlab="")
lines(data[['periodv']],data[['pages']][,27],ty='l',lwd=linew2,lty=2,col=rgb(0,10,10,5,maxColorValue=10), 
      ylab="",xlab="")
ecormean <- round(cor(data[['allind']][,1],data[['pages']][,27],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][,2],data[['pages']][,27],use="pairwise.complete.obs"),digits=2) 
ecorbem <- round(cor(data[['allind']][,3],data[['pages']][,27],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][,4],data[['pages']][,27],use="pairwise.complete.obs"),digits=2)
ecorval <- round(cor(data[['allind']][300:377,1],vindmean,use="pairwise.complete.obs"),digits=2) 
acorval <- round(cor(data[['allind']][300:377,2],vindmean,use="pairwise.complete.obs"),digits=2) 
ecorvalbem <- round(cor(data[['allind']][300:377,3],vindmean,use="pairwise.complete.obs"),digits=2) 
acorvalbem <- round(cor(data[['allind']][300:377,4],vindmean,use="pairwise.complete.obs"),digits=2) 
if (plbem) {
  legend('topleft', c(paste('ECHAM ens. mean - PAGES recon, CRU TS3:',ecormean,ecorval),
                    paste('Analysis ens. mean - PAGES recon, CRU TS3:',acormean,acorval),
                    paste('ECHAM BEM - PAGES recon, CRU TS3:',ecorbem,ecorvalbem),
                    paste('Analysis BEM - PAGES recon, CRU TS3:',acorbem,acorvalbem)), 
       cex=fs, bty='o', bg='white', box.col='white') 
} else {
  legend('topleft', c(paste('ECHAM ens. mean - PAGES recon, CRU TS3:',ecormean,ecorval),
                      paste('Analysis ens. mean - PAGES recon, CRU TS3:',acormean,acorval)),
         cex=fs, bty='o', bg='white', box.col='white') 
}
  
# Africa
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='AFR.temp2'),lw=linew,sca=scal, yl='',xa='n',plotbem=plbem)
reg <- which(vind.allts$names=='AFR.temp2')
years <- rep(vind.allts$time,each=2) 
years <- years[(length(years)-155):length(years)]
vdat <- vind.allts$data[reg,(ncol(vind.allts$data)-155):ncol(vind.allts$data)]
vindmean <- aggregate(vdat,list(years),mean)[,2]
lines(seq(1902,1979),scale(vindmean,center=T,scale=scal),ty='l', 
      lwd=linew, lty=1, col=rgb(0,0,10,5,maxColorValue=10), ylab="",xlab="")
ecorval <- round(cor(data[['allind']][300:377,1],vindmean,use="pairwise.complete.obs"),digits=2) 
acorval <- round(cor(data[['allind']][300:377,2],vindmean,use="pairwise.complete.obs"),digits=2) 
ecorvalbem <- round(cor(data[['allind']][300:377,3],vindmean,use="pairwise.complete.obs"),digits=2) 
acorvalbem <- round(cor(data[['allind']][300:377,4],vindmean,use="pairwise.complete.obs"),digits=2) 
if (plbem) {
  legend('topleft', c(paste('ECHAM ens. mean - CRU TS3:',ecorval),
                    paste('Analysis ens. mean - CRU TS3:',acorval),
                    paste('ECHAM BEM - CRU TS3:',ecorvalbem),
                    paste('Analysis BEM - CRU TS3:',acorvalbem)), 
       cex=fs, bty='o', bg='white', box.col='white') 
} else {
  legend('topleft', c(paste('ECHAM ens. mean - PAGES recon, CRU TS3:',ecormean,ecorval),
                      paste('Analysis ens. mean - PAGES recon, CRU TS3:',acormean,acorval)),
         cex=fs, bty='o', bg='white', box.col='white') 
}

# South America
data <- plot_pages(wl=wlen,reg=which(attributes(giorgi)$names=='SAM'),lw=linew,sca=scal, yl='',xa='n',plotbem=plbem)
reg <- which(vind.allts$names=='SAM.temp2')
years <- rep(vind.allts$time,each=2) 
years <- years[(length(years)-155):length(years)]
vdat <- vind.allts$data[reg,(ncol(vind.allts$data)-155):ncol(vind.allts$data)]
vindmean <- aggregate(vdat,list(years),mean)[,2]
lines(seq(1902,1979),scale(vindmean,center=T,scale=scal),ty='l', 
      lwd=linew, lty=1, col=rgb(0,0,10,5,maxColorValue=10), ylab="",xlab="")
lines(data[['periodv']],data[['pages']][,18],ty='l',lwd=linew2,lty=2,col=rgb(10,0,10,5,maxColorValue=10), 
      ylab="",xlab="")
lines(data[['periodv']],data[['pages']][,19],ty='l',lwd=linew2,lty=2,col=rgb(0,10,10,5,maxColorValue=10), 
      ylab="",xlab="")
ecormean <- round(cor(data[['allind']][,1],data[['pages']][,19],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][,2],data[['pages']][,19],use="pairwise.complete.obs"),digits=2) 
ecorbem <- round(cor(data[['allind']][,3],data[['pages']][,19],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][,4],data[['pages']][,19],use="pairwise.complete.obs"),digits=2)
ecorval <- round(cor(data[['allind']][300:377,1],vindmean,use="pairwise.complete.obs"),digits=2) 
acorval <- round(cor(data[['allind']][300:377,2],vindmean,use="pairwise.complete.obs"),digits=2) 
ecorvalbem <- round(cor(data[['allind']][300:377,3],vindmean,use="pairwise.complete.obs"),digits=2) 
acorvalbem <- round(cor(data[['allind']][300:377,4],vindmean,use="pairwise.complete.obs"),digits=2) 
if (plbem) {
  legend('topleft', c(paste('ECHAM ens. mean - PAGES recon, CRU TS3:',ecormean,ecorval),
                    paste('Analysis ens. mean - PAGES recon, CRU TS3:',acormean,acorval),
                    paste('ECHAM BEM - PAGES recon, CRU TS3:',ecorbem,ecorvalbem),
                    paste('Analysis BEM - PAGES recon, CRU TS3:',acorbem,acorvalbem)), 
       cex=fs, bty='o', bg='white', box.col='white') 
} else {
  legend('topleft', c(paste('ECHAM ens. mean - PAGES recon, CRU TS3:',ecormean,ecorval),
                      paste('Analysis ens. mean - PAGES recon, CRU TS3:',acormean,acorval)),
         cex=fs, bty='o', bg='white', box.col='white') 
}

# Australia
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='AUS.temp2'),lw=linew,sca=scal, yl='',xa='s',plotbem=plbem)
reg <- which(vind.allts$names=='AUS.temp2')
years <- rep(vind.allts$time,each=2) 
years <- years[(length(years)-155):length(years)]
vdat <- vind.allts$data[reg,(ncol(vind.allts$data)-155):ncol(vind.allts$data)]
vindmean <- aggregate(vdat,list(years),mean)[,2]
lines(seq(1902,1979),scale(vindmean,center=T,scale=scal),ty='l', 
      lwd=linew, lty=1, col=rgb(0,0,10,5,maxColorValue=10), ylab="",xlab="")
lines(data[['periodv']],data[['pages']][,10],ty='l',lwd=linew2,lty=2,col=rgb(10,0,10,5,maxColorValue=10), 
      ylab="",xlab="")
lines(data[['periodv']],data[['pages']][,12],ty='l',lwd=linew2,lty=2,col=rgb(0,10,10,5,maxColorValue=10), 
      ylab="",xlab="")
ecormean <- round(cor(data[['allind']][,1],data[['pages']][,12],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][,2],data[['pages']][,12],use="pairwise.complete.obs"),digits=2) 
ecorbem <- round(cor(data[['allind']][,3],data[['pages']][,12],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][,4],data[['pages']][,12],use="pairwise.complete.obs"),digits=2)
ecorval <- round(cor(data[['allind']][300:377,1],vindmean,use="pairwise.complete.obs"),digits=2) 
acorval <- round(cor(data[['allind']][300:377,2],vindmean,use="pairwise.complete.obs"),digits=2) 
ecorvalbem <- round(cor(data[['allind']][300:377,3],vindmean,use="pairwise.complete.obs"),digits=2) 
acorvalbem <- round(cor(data[['allind']][300:377,4],vindmean,use="pairwise.complete.obs"),digits=2) 
if (plbem) {
  legend('topleft', c(paste('ECHAM ens. mean - PAGES recon, CRU TS3:',ecormean,ecorval),
                    paste('Analysis ens. mean - PAGES recon, CRU TS3:',acormean,acorval),
                    paste('ECHAM BEM - PAGES recon, CRU TS3:',ecorbem,ecorvalbem),
                    paste('Analysis BEM - PAGES recon, CRU TS3:',acorbem,acorvalbem)), 
       cex=fs, bty='o', bg='white', box.col='white') 
} else {
  legend('topleft', c(paste('ECHAM ens. mean - PAGES recon, CRU TS3:',ecormean,ecorval),
                      paste('Analysis ens. mean - PAGES recon, CRU TS3:',acormean,acorval)),
         cex=fs, bty='o', bg='white', box.col='white') 
}

# Antarctica
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='ANT.temp2'),box='n',lw=linew,sca=scal, yl='',xa='s',plotbem=plbem)
#reg <- which(vind.allts$names=='ANT.temp2')
#years <- rep(vind.allts$time,each=2) 
#years <- years[(length(years)-155):length(years)]
#vdat <- vind.allts$data[reg,(ncol(vind.allts$data)-155):ncol(vind.allts$data)]
#vindmean <- aggregate(vdat,list(years),mean)[,2]
#lines(seq(1902,1979),scale(vindmean,center=T,scale=scal),ty='l', 
#lwd=linew, lty=1, col=rgb(0,10,10,5,maxColorValue=10), ylab="",xlab="")
lines(data[['periodv']],data[['pages']][,15],ty='l',lwd=linew2,lty=2,col=rgb(0,10,10,5,maxColorValue=10), 
      ylab="",xlab="")
ecormean <- round(cor(data[['allind']][,1],data[['pages']][,15],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][,2],data[['pages']][,15],use="pairwise.complete.obs"),digits=2) 
ecorbem <- round(cor(data[['allind']][,3],data[['pages']][,15],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][,4],data[['pages']][,15],use="pairwise.complete.obs"),digits=2)
#ecorval <- round(cor(data[['allind']][300:377,1],vindmean,use="pairwise.complete.obs"),digits=2) 
#acorval <- round(cor(data[['allind']][300:377,2],vindmean,use="pairwise.complete.obs"),digits=2) 
#ecorvalbem <- round(cor(data[['allind']][300:377,3],vindmean,use="pairwise.complete.obs"),digits=2) 
#acorvalbem <- round(cor(data[['allind']][300:377,4],vindmean,use="pairwise.complete.obs"),digits=2) 
if (plbem) {
  legend('topleft', c(paste('ECHAM ens. mean - PAGES recon:',ecormean),
                    paste('Analysis ens. mean - PAGES recon:',acormean),
                    paste('ECHAM BEM - PAGES recon:',ecorbem),
                    paste('Analysis BEM - PAGES recon:',acorbem)), 
       cex=fs, bty='o', bg='white', box.col='white') 
} else {
  legend('topleft', c(paste('ECHAM ens. mean - PAGES recon, CRU TS3:',ecormean,ecorval),
                      paste('Analysis ens. mean - PAGES recon, CRU TS3:',acormean,acorval)),
         cex=fs, bty='o', bg='white', box.col='white') 
}

legend('bottomleft', c('ECHAM ens. mean','ECHAM spread','Analysis ens. mean','Analysis spread',
#'Best ECHAM ensemble member','Best Analysis ensemble member',
'CRU TS3','PAGES reconstruction','Other reconstruction'),lwd=c(2,12,2,12,2,2,2),lty=c(1,1,1,1,1,2,2), col=c('black','lightgrey','red','indianred1',
#'grey','darkorange',                                                                                                                                                                                                                                                                              
'blue','cyan','violet'),cex=fs, bty='o', bg='white', box.col='white')

dev.off()












# plot global and hemisperic averages
pdf(file=paste0('../figures/',expname,'/glo_hemis_temp_comp',filenameext,'.pdf'),width=5,height=10)
wlen <- 1
linew <- 1
par(mfrow=c(4,1))
par(mai=c(0.3,0.6,0.2,0))
fs <- 1.0
anomper=(1901:1980)
scal=F
plbem=F
#eind.allts$names
noaa <- read.table('../comparison_data/NOAA_PCN_1600_onward_stefan.txt',skip=3,header=F,sep='\t',stringsAsFactors=FALSE)
colnames(noaa) <- read.table('../comparison_data/NOAA_PCN_1600_onward_stefan.txt',nrow=1,sep='\t',stringsAsFactors=FALSE)

# global mean land only
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='GLO.temp2'),seas='yrmean',anomper=(1901:1980),lw=linew,sca=scal,xa='n',plotbem=plbem)

nc=open.nc('../comparison_data/cru3_tmp_25_annmean_glomean.nc', write=F)
print.nc(nc)
tmp1 <- var.get.nc(nc, "TMP") # for CRU temp
tlist <- var.get.nc(nc, "TIME")
unitstr <- "months since 1901-01-15 00:00:00"
yrv <- utcal.nc(unitstr, tlist,type="n")[,1]
pos <- yrv %in% data[['period']]
tmp1 <- tmp1[pos]
yrv <- yrv[pos]
anopos <- yrv %in% anomper
centerfac <- mean(tmp1[anopos])
if (wlen > 1) {
  tmp1 <- runmean(scale(tmp1,center=centerfac,scale=F),wlen)
} else {
  tmp1 <- scale(tmp1,center=centerfac,scale=F)
}
lines(yrv,tmp1,ty='l',lwd=linew,lty=1,col=rgb(0,0,10,5,maxColorValue=10), ylab="",xlab="")

lec <- read.table('../comparison_data/leclercq2012global.txt',header=T,skip=83,stringsAsFactors=FALSE)
tmp2 <- lec[,2]
yrv <- lec[,1]
pos <- yrv %in% data[['period']]
tmp2 <- tmp2[pos]
yrv <- yrv[pos]
anopos <- yrv %in% anomper
centerfac <- mean(tmp2[anopos])
if (wlen > 1) {
  tmp2 <- runmean(scale(tmp2,center=centerfac,scale=F),wlen)
} else {
  tmp2 <- scale(tmp2,center=centerfac,scale=F)
}
lines(yrv,tmp2,ty='l',lwd=linew,lty=1,col=rgb(0,10,10,5,maxColorValue=10), ylab="",xlab="")

tmp3 <- noaa[,which(colnames(noaa)=="mann2008e")]
yrv <- noaa[,1]
pos <- yrv %in% data[['period']]
tmp3 <- tmp3[pos]
yrv <- yrv[pos]
anopos <- yrv %in% anomper
centerfac <- mean(tmp3[anopos])
if (wlen > 1) {
  tmp3 <- runmean(scale(tmp3,center=centerfac,scale=F),wlen)
} else {
  tmp3 <- scale(tmp3,center=centerfac,scale=F)
}
lines(yrv,tmp3,ty='l',lwd=linew,lty=1,col=rgb(10,0,10,5,maxColorValue=10), ylab="",xlab="")

tmp4 <- noaa[,which(colnames(noaa)=="oerlemans2005")]
yrv <- noaa[,1]
pos <- yrv %in% data[['period']]
tmp4 <- tmp4[pos]
yrv <- yrv[pos]
anopos <- yrv %in% anomper
centerfac <- mean(tmp4[anopos])
if (wlen > 1) {
  tmp4 <- runmean(scale(tmp4,center=centerfac,scale=F),wlen)
} else {
  tmp4 <- scale(tmp4,center=centerfac,scale=F)
}
lines(yrv,tmp4,ty='l',lwd=linew,lty=1,col=rgb(10,10,0,5,maxColorValue=10), ylab="",xlab="")

ecormean <- round(cor(data[['allind']][,1],tmp3,use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][,2],tmp3,use="pairwise.complete.obs"),digits=2) 
ecorcru <- round(cor(data[['allind']][299:length(data[['period']]),1],tmp1,use="pairwise.complete.obs"),digits=2) 
acorcru <- round(cor(data[['allind']][299:length(data[['period']]),2],tmp1,use="pairwise.complete.obs"),digits=2) 
legend('topleft', c(paste('ECHAM ens. mean - Mann 2008, CRU TS3:',ecormean,ecorcru),
                    paste('Analysis ens. mean - Mann 2008, CRU TS3:',acormean,acorcru)), 
       cex=fs, bty='o', bg='white', box.col='white')

legend('bottomleft', c('Leclercq 2012','Mann 2008','Oerlemans 2005'),lwd=c(2,2,2),lty=c(1,1,1), col=c('cyan','violet','yellow'),cex=fs, bty='o', bg='white', box.col='white')





# northern hemisphere extra-tropics (>20N)
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='ENH.temp2'),seas='yrmean',anomper=(1901:1980),lw=linew,sca=scal,xa='n',plotbem=plbem)
#data <- plot_pages(wl=wlen,reg=1,seas='sum')
#data <- plot_pages(wl=wlen,reg=1,seas='win')
nc=open.nc('../comparison_data/cru3_tmp_25_annmean_enhmean.nc', write=F)
print.nc(nc)
tmp1 <- var.get.nc(nc, "TMP") # for CRU temp
tlist <- var.get.nc(nc, "TIME")
unitstr <- "months since 1901-01-15 00:00:00"
yrv <- utcal.nc(unitstr, tlist,type="n")[,1]
pos <- yrv %in% data[['period']]
tmp1 <- tmp1[pos]
yrv <- yrv[pos]
anopos <- yrv %in% anomper
centerfac <- mean(tmp1[anopos])
if (wlen > 1) {
  tmp1 <- runmean(scale(tmp1,center=centerfac,scale=F),wlen)
} else {
  tmp1 <- scale(tmp1,center=centerfac,scale=F)
}
lines(yrv,tmp1,ty='l',lwd=linew,lty=1,col=rgb(0,0,10,5,maxColorValue=10), ylab="",xlab="")

cr <- read.table('../comparison_data/crowley_NH_temp_recon-2014.txt',sep='\t',header=T)
tmp2 <- cr[,2]
yrv <- cr[,1]
pos <- yrv %in% data[['period']]
tmp2 <- tmp2[pos]
yrv <- yrv[pos]
pos <- data[['period']] %in%yrv 
cordata <- data[['allind']][pos,]
data[['allind']][,1]
anopos <- yrv %in% anomper
centerfac <- mean(tmp2[anopos])
if (wlen > 1) {
  tmp2 <- runmean(scale(tmp2,center=centerfac,scale=F),wlen)
} else {
  tmp2 <- scale(tmp2,center=centerfac,scale=F)
}
lines(yrv,tmp2,ty='l',lwd=linew,lty=1,col=rgb(0,10,10,5,maxColorValue=10), ylab="",xlab="")

tmp3 <- noaa[,which(colnames(noaa)=="esper2002")]
yrv <- noaa[,1]
pos <- yrv %in% data[['period']]
tmp3 <- tmp3[pos]
yrv <- yrv[pos]
cordata3 <- data[['allind']][pos,]
anopos <- yrv %in% anomper
centerfac <- mean(tmp3[anopos])
if (wlen > 1) {
  tmp3 <- runmean(scale(tmp3,center=centerfac,scale=F),wlen)
} else {
  tmp3 <- scale(tmp3,center=centerfac,scale=F)
}
lines(yrv,tmp3,ty='l',lwd=linew,lty=1,col=rgb(10,0,10,5,maxColorValue=10), ylab="",xlab="")

tmp4 <- noaa[,which(colnames(noaa)=="darrigo2006b")]
yrv <- noaa[,1]
pos <- yrv %in% data[['period']]
tmp4 <- tmp4[pos]
yrv <- yrv[pos]
cordata4 <- data[['allind']][pos,]
anopos <- yrv %in% anomper
centerfac <- mean(tmp4[anopos])
if (wlen > 1) {
  tmp4 <- runmean(scale(tmp4,center=centerfac,scale=F),wlen)
} else {
  tmp4 <- scale(tmp4,center=centerfac,scale=F)
}
lines(yrv,tmp4,ty='l',lwd=linew,lty=1,col=rgb(10,10,0,5,maxColorValue=10), ylab="",xlab="")

# tmp5 <- noaa[,which(colnames(noaa)=="wilson2007")]
# yrv <- noaa[,1]
# pos <- yrv %in% data[['period']]
# tmp5 <- tmp5[pos]
# yrv <- yrv[pos]
# anopos <- yrv %in% anomper
# centerfac <- mean(tmp5[anopos])
# if (wlen > 1) {
#   tmp5 <- runmean(scale(tmp5,center=centerfac,scale=F),wlen)
# } else {
#   tmp5 <- scale(tmp5,center=centerfac,scale=F)
# }
# lines(yrv,tmp5,ty='l',lwd=linew,lty=1,col=rgb(0,10,0,5,maxColorValue=10), ylab="",xlab="")

#ecormean <- round(cor(data[['allind']][,1],tmp2,use="pairwise.complete.obs"),digits=2) 
#acormean <- round(cor(data[['allind']][,2],tmp2,use="pairwise.complete.obs"),digits=2) 
ecormean <- round(cor(cordata[,1],tmp2,use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(cordata[,2],tmp2,use="pairwise.complete.obs"),digits=2) 
ecorcru <- round(cor(data[['allind']][299:length(data[['period']]),1],tmp1,use="pairwise.complete.obs"),digits=2) 
acorcru <- round(cor(data[['allind']][299:length(data[['period']]),2],tmp1,use="pairwise.complete.obs"),digits=2) 
legend('topleft', c(paste('ECHAM ens. mean - Crowley 2014, CRU TS3:',ecormean,ecorcru),
                    paste('Analysis ens. mean - Crowley 2014, CRU TS3:',acormean,acorcru)), 
       cex=fs, bty='o', bg='white', box.col='white')

legend('bottomleft', c('Crowley 2014','Esper 2002','D Arrigo 2006'),lwd=c(2,2,2),lty=c(1,1,1), col=c('cyan','violet','yellow'),cex=fs, bty='o', bg='white', box.col='white')




# northern hemisphere complete land
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='NH.temp2'),seas='yrmean',anomper=(1901:1980),lw=linew,sca=scal,xa='n',plotbem=plbem,use18=F,pag=F)

nc=open.nc('../comparison_data/cru3_tmp_25_annmean_nhmean.nc', write=F)
print.nc(nc)
tmp1 <- var.get.nc(nc, "TMP") # for CRU temp
tlist <- var.get.nc(nc, "TIME")
unitstr <- "months since 1901-01-15 00:00:00"
yrv <- utcal.nc(unitstr, tlist,type="n")[,1]
pos <- yrv %in% data[['period']]
tmp1 <- tmp1[pos]
yrv <- yrv[pos]
anopos <- yrv %in% anomper
centerfac <- mean(tmp1[anopos])
if (wlen > 1) {
  tmp1 <- runmean(scale(tmp1,center=centerfac,scale=F),wlen)
} else {
  tmp1 <- scale(tmp1,center=centerfac,scale=F)
}
lines(yrv,tmp1,ty='l',lwd=linew,lty=1,col=rgb(0,0,10,5,maxColorValue=10), ylab="",xlab="")

lec <- read.table('../comparison_data/leclercq2012nh.txt',header=T,skip=83,stringsAsFactors=FALSE)
tmp2 <- lec[,2]
yrv <- lec[,1]
pos <- yrv %in% data[['period']]
tmp2 <- tmp2[pos]
yrv <- yrv[pos]
anopos <- yrv %in% anomper
centerfac <- mean(tmp2[anopos])
if (wlen > 1) {
  tmp2 <- runmean(scale(tmp2,center=centerfac,scale=F),wlen)
} else {
  tmp2 <- scale(tmp2,center=centerfac,scale=F)
}
lines(yrv,tmp2,ty='l',lwd=linew,lty=1,col=rgb(0,10,10,5,maxColorValue=10), ylab="",xlab="")

tmp3 <- noaa[,which(colnames(noaa)=="mann2008g")]
yrv <- noaa[,1]
pos <- yrv %in% data[['period']]
tmp3 <- tmp3[pos]
yrv <- yrv[pos]
anopos <- yrv %in% anomper
centerfac <- mean(tmp3[anopos])
if (wlen > 1) {
  tmp3 <- runmean(scale(tmp3,center=centerfac,scale=F),wlen)
} else {
  tmp3 <- scale(tmp3,center=centerfac,scale=F)
}
lines(yrv,tmp3,ty='l',lwd=linew,lty=1,col=rgb(10,0,10,5,maxColorValue=10), ylab="",xlab="")

tmp4 <- noaa[,which(colnames(noaa)=="ammann2007")]
yrv <- noaa[,1]
pos <- yrv %in% data[['period']]
tmp4 <- tmp4[pos]
yrv <- yrv[pos]
anopos <- yrv %in% anomper
centerfac <- mean(tmp4[anopos])
if (wlen > 1) {
  tmp4 <- runmean(scale(tmp4,center=centerfac,scale=F),wlen)
} else {
  tmp4 <- scale(tmp4,center=centerfac,scale=F)
}
lines(yrv,tmp4,ty='l',lwd=linew,lty=1,col=rgb(10,10,0,5,maxColorValue=10), ylab="",xlab="")

tmp5 <- noaa[,which(colnames(noaa)=="moberg2005")]
yrv <- noaa[,1]
pos <- yrv %in% data[['period']]
tmp5 <- tmp5[pos]
yrv <- yrv[pos]
anopos <- yrv %in% anomper
centerfac <- mean(tmp5[anopos],na.rm=T)
if (wlen > 1) {
  tmp5 <- runmean(scale(tmp5,center=centerfac,scale=F),wlen)
} else {
  tmp5 <- scale(tmp5,center=centerfac,scale=F)
}
lines(yrv,tmp5,ty='l',lwd=linew,lty=1,col=rgb(0,10,0,5,maxColorValue=10), ylab="",xlab="")

ecormean <- round(cor(data[['allind']][,1],tmp4,use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][,2],tmp4,use="pairwise.complete.obs"),digits=2)
ecorcru <- round(cor(data[['allind']][299:length(data[['period']]),1],tmp1,use="pairwise.complete.obs"),digits=2) 
acorcru <- round(cor(data[['allind']][299:length(data[['period']]),2],tmp1,use="pairwise.complete.obs"),digits=2) 
legend('topleft', c(paste('ECHAM ens. mean - Ammann 2007, CRU TS3:',ecormean,ecorcru),
                    paste('Analysis ens. mean - Ammann 2007, CRU TS3:',acormean,acorcru)), 
       cex=fs, bty='o', bg='white', box.col='white')

legend('bottomleft', c('Leclercq 2012','Mann 2008', 'Ammann 2007', 'Moberg 2005'),lwd=c(2,2,2,2),lty=c(1,1,1,1), col=c('cyan','violet','yellow','green'),cex=fs, bty='o', bg='white', box.col='white')




# southern hemisphere extra-tropical land ( < -20N )
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='SH.temp2'),seas='yrmean',anomper=(1901:1980),lw=linew,sca=scal,xa='s',plotbem=plbem, pag=F, use18=F, plotmem=T)

nc=open.nc('../comparison_data/cru3_tmp_25_annmean_shmean.nc', write=F)
print.nc(nc)
tmp1 <- var.get.nc(nc, "TMP") # for CRU temp
tlist <- var.get.nc(nc, "TIME")
unitstr <- "months since 1901-01-15 00:00:00"
yrv <- utcal.nc(unitstr, tlist,type="n")[,1]
pos <- yrv %in% data[['period']]
tmp1 <- tmp1[pos]
yrv <- yrv[pos]
anopos <- yrv %in% anomper
centerfac <- mean(tmp1[anopos])
if (wlen > 1) {
  tmp1 <- runmean(scale(tmp1,center=centerfac,scale=F),wlen)
} else {
  tmp1 <- scale(tmp1,center=centerfac,scale=F)
}
lines(yrv,tmp1,ty='l',lwd=linew,lty=1,col=rgb(0,0,10,5,maxColorValue=10), ylab="",xlab="")

lec <- read.table('../comparison_data/leclercq2012sh.txt',header=T,skip=83,stringsAsFactors=FALSE)
tmp2 <- lec[,2]
yrv <- lec[,1]
pos <- yrv %in% data[['period']]
tmp2 <- tmp2[pos]
yrv <- yrv[pos]
anopos <- yrv %in% anomper
centerfac <- mean(tmp2[anopos])
if (wlen > 1) {
  tmp2 <- runmean(scale(tmp2,center=centerfac,scale=F),wlen)
} else {
  tmp2 <- scale(tmp2,center=centerfac,scale=F)
}
lines(yrv,tmp2,ty='l',lwd=linew,lty=1,col=rgb(0,10,10,5,maxColorValue=10), ylab="",xlab="")

tmp3 <- noaa[,which(colnames(noaa)=="mann2008i")]
yrv <- noaa[,1]
pos <- yrv %in% data[['period']]
tmp3 <- tmp3[pos]
yrv <- yrv[pos]
anopos <- yrv %in% anomper
centerfac <- mean(tmp3[anopos])
if (wlen > 1) {
  tmp3 <- runmean(scale(tmp3,center=centerfac,scale=F),wlen)
} else {
  tmp3 <- scale(tmp3,center=centerfac,scale=F)
}
lines(yrv,tmp3,ty='l',lwd=linew,lty=1,col=rgb(10,0,10,5,maxColorValue=10), ylab="",xlab="")

nk <- read.table('../comparison_data/SH_temp_recon_neukom2014.txt',header=T,skip=106,stringsAsFactors=FALSE)
tmp4 <- nk[,2]
yrv <- nk[,1]
pos <- yrv %in% data[['period']]
tmp4 <- tmp4[pos]
yrv <- yrv[pos]
anopos <- yrv %in% anomper
centerfac <- mean(tmp4[anopos])
if (wlen > 1) {
  tmp4 <- runmean(scale(tmp4,center=centerfac,scale=F),wlen)
} else {
  tmp4 <- scale(tmp4,center=centerfac,scale=F)
}
lines(yrv,tmp4,ty='l',lwd=linew,lty=1,col=rgb(10,10,0,5,maxColorValue=10), ylab="",xlab="")

ecormean <- round(cor(data[['allind']][,1],tmp4,use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][,2],tmp4,use="pairwise.complete.obs"),digits=2) 
ecorcru <- round(cor(data[['allind']][299:length(data[['period']]),1],tmp1,use="pairwise.complete.obs"),digits=2) 
acorcru <- round(cor(data[['allind']][299:length(data[['period']]),2],tmp1,use="pairwise.complete.obs"),digits=2) 
legend('topleft', c(paste('ECHAM ens. mean - Neukom 2014, CRU TS3:',ecormean,ecorcru),
                    paste('Analysis ens. mean - Neukom 2014, CRU TS3:',acormean,acorcru)), 
       cex=fs, bty='o', bg='white', box.col='white')
if (plbem) {
  legend('bottomleft', c('ECHAM ens. mean','ECHAM spread','Analysis ens. mean','Analysis spread','Best ECHAM ensemble member','Best Analysis ensemble member','CRU TS3','Leclercq 2012','Mann 2008','Neukom 2014'),lwd=c(2,12,2,12,2,2,2,2,2,2),lty=c(1,1,1,1,2,2,1,1,1,1), col=c('black','lightgrey','red','indianred1','grey','darkorange','blue','cyan','violet','yellow'),cex=fs, bty='o', bg='white', box.col='white')
} else {
  legend('bottomleft', c('ECHAM ens. mean','ECHAM spread','Analysis ens. mean','Analysis spread','CRU TS3','Leclercq 2012','Mann 2008','Neukom 2014'),lwd=c(2,12,2,12,2,2,2,2),lty=c(1,1,1,1,1,1,1,1), col=c('black','lightgrey','red','indianred1','blue','cyan','violet','yellow'),cex=fs, bty='o', bg='white', box.col='white')
}  

# load instr. and recon data for comparison
#system("cdo -fldmean -sellonlatbox,-180,180,20,90 /Users/joerg/Documents/climdata/cru/HadCRUT.4.2.0.0.median.nc ../comparison_data/crutem4_enh_mean.nc")
 
dev.off()







# PLOT atm indices
load(paste0(dataintdir,'indices',expname,'/indices',filenameext,'allts.Rdata'))  
#load(paste0(prepplotdir,'indices',filenameext,'allts.Rdata'))  
#load(paste0(prepplotdir,'indices_allts.Rdata')) 
# some ocean based indices only NOT in landonly version
#load('../data/prepplot_season/indices_allts_landonly.Rdata')
#stefindmon <- read.table("/Users/joerg/Documents/climdata/indices/stefan/stefan_monthly_indices.txt",header=T)
#stefindyrmean <- aggregate(stefindmon[,c(6,9,12,15,18,21)],list(stefindmon[,1]),mean,na.rm=T)
stefind <- read.table("/Users/joerg/Documents/climdata/indices/stefan/stefan_seasonal_indices.txt",header=T,skip=1,na.string=-9999)[,c(1,seq(2,5),seq(18,21),seq(34,37),seq(50,53),seq(66,69),seq(82,85),seq(98,101),seq(114,117),seq(130,133),seq(146,149))]
head <- read.table("/Users/joerg/Documents/climdata/indices/stefan/stefan_seasonal_indices.txt",nrow=2)[,c(1,seq(2,5),seq(18,21),seq(34,37),seq(50,53),seq(66,69),seq(82,85),seq(98,101),seq(114,117),seq(130,133),seq(146,149))]
# (1,2,3,18,19,34,35,50,51,66,67,82,83,98,99,114,115,130,131,146,147)]
newhead <- rep(NA,ncol(head))
for (i in 1:ncol(head)) {
  newhead[i] <- paste(head[1,i],head[2,i],sep='_')
}
colnames(stefind) <- newhead
#load('../data/indices_recon.Rdata')
vf <- c(1600,1650,1660,1783,1809,1815,1835,1883,1886,1902,1912,1980,1991) #read.table('/Users/joerg/Documents/climdata/forcings/forc-total-crowley_2001.txt',skip=51,header=T)[,1:2]

#pdf(file=paste0('../figures/indices/atm_index_comp',filenameext,'.pdf'),width=20,height=22)
pdf(file=paste0('../figures/',expname,'/atm_index_comp',filenameext,'.pdf'),width=20,height=22)
wlen <- 1
linew <- 1
par(mfrow=c(9,2))
par(mai=c(0.3,0.6,0.2,0))
fs <- 1.0
anomper=(1901:1980)
#eind.allts$names

# HC
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='HC.calc'),seas='win',anomper=(1901:1980),lw=linew,sca=T,yl='')
vdata1 <- stefind['HC_20C'] # [,'SJ_20C'] 
vdata2 <- stefind['HC_REC'] # [,'SJ_REC'] 
vdata3 <- stefind['HC_NCEP']
vdata4 <- stefind['HC_ERA-40']
vdata <- cbind(vdata1,vdata2,vdata3,vdata4)
yrv <- stefind[,1]
centerfac <- apply(vdata,2,mean,na.rm=T)
scalefac <- apply(vdata,2,sd,na.rm=T)
#centerfac2 <- mean(vdata2,na.rm=T)
#scalefac2 <- sd(vdata2,na.rm=T)
if (wlen > 1) {
  vdata <- runmean(scale(vdata,center=centerfac,scale=scalefac),wlen)
#  vdata1 <- runmean(scale(vdata1,center=centerfac1,scale=scalefac1),wlen)
#  vdata2 <- runmean(scale(vdata2,center=centerfac2,scale=scalefac2),wlen)
} else {
  vdata <- scale(vdata,center=centerfac,scale=scalefac)
#  vdata1 <- scale(vdata1,center=centerfac1,scale=scalefac1)
#  vdata2 <- scale(vdata2,center=centerfac2,scale=scalefac2)
}
lines(yrv,vdata[,1],ty='l',lwd=linew,lty=1,col=rgb(0,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,2],ty='l',lwd=linew,lty=1,col=rgb(10,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,3],ty='l',lwd=linew,lty=1,col=rgb(10,10,0,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,4],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,10,maxColorValue=10), ylab="",xlab="")
ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
legend('bottomleft', c(paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',ecormean[1],ecormean[2],ecormean[3],ecormean[4]),paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',acormean[1],acormean[2],acormean[3],acormean[4]),paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],ecorbem[2],ecorbem[3],ecorbem[4]),paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],acorbem[2],acorbem[3],acorbem[4])),cex=fs, bty='o', bg='white', box.col='white')
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='HC'),seas='win',anomper=(1901:1980),lw=linew,sca=T,yl='')
ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
legend('bottomleft', c(paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',ecormean[1],ecormean[2],ecormean[3],ecormean[4]),paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',acormean[1],acormean[2],acormean[3],acormean[4]),paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],ecorbem[2],ecorbem[3],ecorbem[4]),paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],acorbem[2],acorbem[3],acorbem[4])),cex=fs, bty='o', bg='white', box.col='white')


# SJ_u200
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='SJ_u200.calc'),seas='win',anom=T,anomper=(1901:1980),lw=linew,sca=F,yl='',use18=T)
vdata <- cbind(stefind['SJ_20C'],stefind['SJ_REC'],stefind['SJ_NCEP'],stefind['SJ_ERA-40'])
yrv <- stefind[,1]
centerfac <- apply(vdata,2,mean,na.rm=T)
scalefac <- apply(vdata,2,sd,na.rm=T)
if (wlen > 1) {
  vdata <- runmean(scale(vdata,center=centerfac,scale=scalefac),wlen)
} else {
  vdata <- scale(vdata,center=centerfac,scale=scalefac)
}
lines(yrv,vdata[,1],ty='l',lwd=linew,lty=1,col=rgb(0,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,2],ty='l',lwd=linew,lty=1,col=rgb(10,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,3],ty='l',lwd=linew,lty=1,col=rgb(10,10,0,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,4],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,10,maxColorValue=10), ylab="",xlab="")
ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
legend('bottomleft', c(paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',ecormean[1],ecormean[2],ecormean[3],ecormean[4]),paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',acormean[1],acormean[2],acormean[3],acormean[4]),paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],ecorbem[2],ecorbem[3],ecorbem[4]),paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],acorbem[2],acorbem[3],acorbem[4])),cex=fs, bty='o', bg='white', box.col='white')
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='SJ'),seas='win',anomper=(1901:1980),lw=linew,sca=T,yl='')
ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
legend('bottomleft', c(paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',ecormean[1],ecormean[2],ecormean[3],ecormean[4]),paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',acormean[1],acormean[2],acormean[3],acormean[4]),paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],ecorbem[2],ecorbem[3],ecorbem[4]),paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],acorbem[2],acorbem[3],acorbem[4])),cex=fs, bty='o', bg='white', box.col='white')

# SJ_slp
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='SJ_slp.calc'),seas='win',anomper=(1901:1980),lw=linew,sca=T,yl='',use18=T)
vdata <- cbind(stefind['SJ_20C'],stefind['SJ_REC'],stefind['SJ_NCEP'],stefind['SJ_ERA-40'])
yrv <- stefind[,1]
centerfac <- apply(vdata,2,mean,na.rm=T)
scalefac <- apply(vdata,2,sd,na.rm=T)
if (wlen > 1) {
  vdata <- runmean(scale(vdata,center=centerfac,scale=scalefac),wlen)
} else {
  vdata <- scale(vdata,center=centerfac,scale=scalefac)
}
lines(yrv,vdata[,1],ty='l',lwd=linew,lty=1,col=rgb(0,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,2],ty='l',lwd=linew,lty=1,col=rgb(10,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,3],ty='l',lwd=linew,lty=1,col=rgb(10,10,0,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,4],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,10,maxColorValue=10), ylab="",xlab="")
ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
legend('bottomleft', c(paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',ecormean[1],ecormean[2],ecormean[3],ecormean[4]),paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',acormean[1],acormean[2],acormean[3],acormean[4]),paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],ecorbem[2],ecorbem[3],ecorbem[4]),paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],acorbem[2],acorbem[3],acorbem[4])),cex=fs, bty='o', bg='white', box.col='white')
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='SJ'),seas='win',anomper=(1901:1980),lw=linew,sca=T,yl='')
ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
legend('bottomleft', c(paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',ecormean[1],ecormean[2],ecormean[3],ecormean[4]),paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',acormean[1],acormean[2],acormean[3],acormean[4]),paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],ecorbem[2],ecorbem[3],ecorbem[4]),paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],acorbem[2],acorbem[3],acorbem[4])),cex=fs, bty='o', bg='white', box.col='white')

# PWC
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='PWC.calc'),seas='win',anomper=(1901:1980),lw=linew,sca=T,yl='')
vdata <- cbind(stefind['PWC_20C'],stefind['PWC_REC'],stefind['PWC_NCEP'],stefind['PWC_ERA-40'])
yrv <- stefind[,1]
centerfac <- apply(vdata,2,mean,na.rm=T)
scalefac <- apply(vdata,2,sd,na.rm=T)
if (wlen > 1) {
  vdata <- runmean(scale(vdata,center=centerfac,scale=scalefac),wlen)
} else {
  vdata <- scale(vdata,center=centerfac,scale=scalefac)
}
lines(yrv,vdata[,1],ty='l',lwd=linew,lty=1,col=rgb(0,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,2],ty='l',lwd=linew,lty=1,col=rgb(10,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,3],ty='l',lwd=linew,lty=1,col=rgb(10,10,0,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,4],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,10,maxColorValue=10), ylab="",xlab="")
ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
legend('bottomleft', c(paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',ecormean[1],ecormean[2],ecormean[3],ecormean[4]),paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',acormean[1],acormean[2],acormean[3],acormean[4]),paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],ecorbem[2],ecorbem[3],ecorbem[4]),paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],acorbem[2],acorbem[3],acorbem[4])),cex=fs, bty='o', bg='white', box.col='white')
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='PWC'),seas='win',anomper=(1901:1980),lw=linew,sca=T,yl='')
ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
legend('bottomleft', c(paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',ecormean[1],ecormean[2],ecormean[3],ecormean[4]),paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',acormean[1],acormean[2],acormean[3],acormean[4]),paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],ecorbem[2],ecorbem[3],ecorbem[4]),paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],acorbem[2],acorbem[3],acorbem[4])),cex=fs, bty='o', bg='white', box.col='white')

# DIMI
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='DIMI.calc'),seas='sum',anomper=(1901:1980),lw=linew,sca=T,yl='')
vdata <- cbind(stefind['DIMI_20C'],stefind['DIMI_REC'],stefind['DIMI_NCEP'],stefind['DIMI_ERA-40'])
yrv <- stefind[,1]
centerfac <- apply(vdata,2,mean,na.rm=T)
scalefac <- apply(vdata,2,sd,na.rm=T)
if (wlen > 1) {
  vdata <- runmean(scale(vdata,center=centerfac,scale=scalefac),wlen)
} else {
  vdata <- scale(vdata,center=centerfac,scale=scalefac)
}
lines(yrv,vdata[,1],ty='l',lwd=linew,lty=1,col=rgb(0,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,2],ty='l',lwd=linew,lty=1,col=rgb(10,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,3],ty='l',lwd=linew,lty=1,col=rgb(10,10,0,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,4],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,10,maxColorValue=10), ylab="",xlab="")
ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
legend('bottomleft', c(paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',ecormean[1],ecormean[2],ecormean[3],ecormean[4]),paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',acormean[1],acormean[2],acormean[3],acormean[4]),paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],ecorbem[2],ecorbem[3],ecorbem[4]),paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],acorbem[2],acorbem[3],acorbem[4])),cex=fs, bty='o', bg='white', box.col='white')
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='DIMI'),seas='sum',anomper=(1901:1980),lw=linew,sca=T,yl='')
ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
legend('bottomleft', c(paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',ecormean[1],ecormean[2],ecormean[3],ecormean[4]),paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',acormean[1],acormean[2],acormean[3],acormean[4]),paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],ecorbem[2],ecorbem[3],ecorbem[4]),paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],acorbem[2],acorbem[3],acorbem[4])),cex=fs, bty='o', bg='white', box.col='white')


# PV
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='PV.calc'),seas='win',anomper=(1901:1980),lw=linew,sca=T,yl='')
#abline(v=vf,col='black',lty=2)
vdata <- cbind(stefind['Z100_20C'],stefind['Z100_REC'],stefind['Z100_NCEP'],stefind['Z100_ERA-40'])
yrv <- stefind[,1]
centerfac <- apply(vdata,2,mean,na.rm=T)
scalefac <- apply(vdata,2,sd,na.rm=T)
if (wlen > 1) {
  vdata <- runmean(scale(vdata,center=centerfac,scale=scalefac),wlen)
} else {
  vdata <- scale(vdata,center=centerfac,scale=scalefac)
}
lines(yrv,vdata[,1],ty='l',lwd=linew,lty=1,col=rgb(0,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,2],ty='l',lwd=linew,lty=1,col=rgb(10,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,3],ty='l',lwd=linew,lty=1,col=rgb(10,10,0,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,4],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,10,maxColorValue=10), ylab="",xlab="")
ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
legend('bottomleft', c(paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',ecormean[1],ecormean[2],ecormean[3],ecormean[4]),paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',acormean[1],acormean[2],acormean[3],acormean[4]),paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],ecorbem[2],ecorbem[3],ecorbem[4]),paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],acorbem[2],acorbem[3],acorbem[4])),cex=fs, bty='o', bg='white', box.col='white')
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='Z100'),seas='win',anomper=(1901:1980),lw=linew,sca=T,yl='')
#abline(v=vf,col='black',lty=2)
ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
legend('bottomleft', c(paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',ecormean[1],ecormean[2],ecormean[3],ecormean[4]),paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',acormean[1],acormean[2],acormean[3],acormean[4]),paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],ecorbem[2],ecorbem[3],ecorbem[4]),paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],acorbem[2],acorbem[3],acorbem[4])),cex=fs, bty='o', bg='white', box.col='white')


# z300
#data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='z300.calc'),seas='win',anomper=(1901:1980),lw=linew,sca=T)
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='Z300'),seas='win',anomper=(1901:1980),lw=linew,sca=T,yl='')
vdata <- cbind(stefind['Z300_20C'],stefind['Z300_REC'],stefind['Z300_NCEP'],stefind['Z300_ERA-40'])
yrv <- stefind[,1]
centerfac <- apply(vdata,2,mean,na.rm=T)
scalefac <- apply(vdata,2,sd,na.rm=T)
if (wlen > 1) {
  vdata <- runmean(scale(vdata,center=centerfac,scale=scalefac),wlen)
} else {
  vdata <- scale(vdata,center=centerfac,scale=scalefac)
}
lines(yrv,vdata[,1],ty='l',lwd=linew,lty=1,col=rgb(0,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,2],ty='l',lwd=linew,lty=1,col=rgb(10,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,3],ty='l',lwd=linew,lty=1,col=rgb(10,10,0,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,4],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,10,maxColorValue=10), ylab="",xlab="")
ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
legend('bottomleft', c(paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',ecormean[1],ecormean[2],ecormean[3],ecormean[4]),paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',acormean[1],acormean[2],acormean[3],acormean[4]),paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],ecorbem[2],ecorbem[3],ecorbem[4]),paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],acorbem[2],acorbem[3],acorbem[4])),cex=fs, bty='o', bg='white', box.col='white')


# ITCZ
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='ITCZ.calc'),seas='win',anomper=(1901:1980),lw=linew,sca=T,yl='')


# NAO
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='NAO.calc'),seas='win',anomper=(1901:1980),lw=linew,sca=T,xa='s',yl='')
#abline(v=vf,col='black',lty=2)
vdata <- cbind(stefind['NAO_20C'],stefind['NAO_REC'],stefind['NAO_NCEP'],stefind['NAO_ERA-40'])
yrv <- stefind[,1]
centerfac <- apply(vdata,2,mean,na.rm=T)
scalefac <- apply(vdata,2,sd,na.rm=T)
if (wlen > 1) {
  vdata <- runmean(scale(vdata,center=centerfac,scale=scalefac),wlen)
} else {
  vdata <- scale(vdata,center=centerfac,scale=scalefac)
}
lines(yrv,vdata[,1],ty='l',lwd=linew,lty=1,col=rgb(0,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,2],ty='l',lwd=linew,lty=1,col=rgb(10,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,3],ty='l',lwd=linew,lty=1,col=rgb(10,10,0,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,4],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,10,maxColorValue=10), ylab="",xlab="")
ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
legend('bottomleft', c(paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',ecormean[1],ecormean[2],ecormean[3],ecormean[4]),paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',acormean[1],acormean[2],acormean[3],acormean[4]),paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],ecorbem[2],ecorbem[3],ecorbem[4]),paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],acorbem[2],acorbem[3],acorbem[4])),cex=fs, bty='o', bg='white', box.col='white')


# PNA
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='PNA.calc'),seas='win',anomper=(1901:1980),lw=linew,sca=T,xa='s',yl='')
#abline(v=vf,col='black',lty=2)
vdata <- cbind(stefind['PNA_20C'],stefind['PNA_REC'],stefind['PNA_NCEP'],stefind['PNA_ERA-40'])
yrv <- stefind[,1]
centerfac <- apply(vdata,2,mean,na.rm=T)
scalefac <- apply(vdata,2,sd,na.rm=T)
if (wlen > 1) {
  vdata <- runmean(scale(vdata,center=centerfac,scale=scalefac),wlen)
} else {
  vdata <- scale(vdata,center=centerfac,scale=scalefac)
}
lines(yrv,vdata[,1],ty='l',lwd=linew,lty=1,col=rgb(0,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,2],ty='l',lwd=linew,lty=1,col=rgb(10,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,3],ty='l',lwd=linew,lty=1,col=rgb(10,10,0,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,4],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,10,maxColorValue=10), ylab="",xlab="")
ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
legend('bottomleft', c(paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',ecormean[1],ecormean[2],ecormean[3],ecormean[4]),paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',acormean[1],acormean[2],acormean[3],acormean[4]),paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],ecorbem[2],ecorbem[3],ecorbem[4]),paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],acorbem[2],acorbem[3],acorbem[4])),cex=fs, bty='o', bg='white', box.col='white')

dev.off()










#rm(list=ls())
#setwd('~/unibe/projects/EnSRF/r/src')
#source('EnSRF_functions.R')
# PLOT atm indices
load(paste0('../data/indices/',expname,'/indices',filenameext,'allts.Rdata'))
#load(paste0(prepplotdir,'indices',filenameext,'allts.Rdata'))
#load(paste0(prepplotdir,'indices_allts.Rdata')) # some ocean based indices only NOT in landonly version
#load('../data/prepplot_season/indices_allts_landonly.Rdata')
#stefindmon <- read.table("/Users/joerg/Documents/climdata/indices/stefan/stefan_monthly_indices.txt",header=T)
#stefindyrmean <- aggregate(stefindmon[,c(6,9,12,15,18,21)],list(stefindmon[,1]),mean,na.rm=T)
stefind <- read.table("/Users/joerg/Documents/climdata/indices/stefan/stefan_seasonal_indices.txt",header=T,skip=1,na.string=-9999)[,c(1,seq(2,5),seq(18,21),seq(34,37),seq(50,53),seq(66,69),seq(82,85),seq(98,101),seq(114,117),seq(130,133),seq(146,149))]
head <- read.table("/Users/joerg/Documents/climdata/indices/stefan/stefan_seasonal_indices.txt",nrow=2)[,c(1,seq(2,5),seq(18,21),seq(34,37),seq(50,53),seq(66,69),seq(82,85),seq(98,101),seq(114,117),seq(130,133),seq(146,149))]
# (1,2,3,18,19,34,35,50,51,66,67,82,83,98,99,114,115,130,131,146,147)]
newhead <- rep(NA,ncol(head))
for (i in 1:ncol(head)) {
  newhead[i] <- paste(head[1,i],head[2,i],sep='_')
}
colnames(stefind) <- newhead
#load('../data/indices_recon.Rdata')
vf <- c(1600,1650,1660,1783,1809,1815,1835,1883,1886,1902,1912,1980,1991) #read.table('/Users/joerg/Documents/climdata/forcings/forc-total-crowley_2001.txt',skip=51,header=T)[,1:2]



#pdf(file=paste0('../figures/indices/atm_index_comp_short',filenameext,'.pdf'),width=10,height=10)
pdf(file=paste0('../figures/',expname,'/atm_index_comp_short',filenameext,'.pdf'),width=10,height=10)
wlen <- 1
linew <- 1
par(mfrow=c(3,1))
par(mai=c(0.3,0.6,0.2,0))
fs <- 1.0
anomper=(1901:1980)
#eind.allts$names

# PV
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='PV.calc'),seas='win',anomper=(1901:1980),
                   lw=linew,sca=T,yl='',plotbem=F,plotmem=T,plotech=T,pag=F)
abline(v=vf,col='black',lty=2)
abline(h=0,col='black',lty=2)
vdata <- cbind(stefind['Z100_20C'],stefind['Z100_REC'],stefind['Z100_NCEP'],stefind['Z100_ERA-40'])
yrv <- stefind[,1]
centerfac <- apply(vdata,2,mean,na.rm=T)
scalefac <- apply(vdata,2,sd,na.rm=T)
if (wlen > 1) {
  vdata <- runmean(scale(vdata,center=centerfac,scale=scalefac),wlen)
} else {
  vdata <- scale(vdata,center=centerfac,scale=scalefac)
}
lines(yrv,vdata[,1],ty='l',lwd=linew,lty=1,col=rgb(0,0,10,8,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,2],ty='l',lwd=linew,lty=1,col=rgb(10,0,10,8,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,3],ty='l',lwd=linew,lty=1,col=rgb(10,10,0,8,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,4],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,8,maxColorValue=10), ylab="",xlab="")
ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
legend('bottomleft', c(
  paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',ecormean[1],',',ecormean[2],',',ecormean[3],',',ecormean[4]),
  paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',acormean[1],',',acormean[2],',',acormean[3],',',acormean[4])),
#  paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],ecorbem[2],ecorbem[3],ecorbem[4]),
#  paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],acorbem[2],acorbem[3],acorbem[4])
  cex=fs, bty='o', bg='white', box.col='white')
legend('topleft', c('ECHAM ens. mean','ECHAM spread','Analysis ens. mean','Analysis spread',
       '20th cent. reanal.','Brönnimann 2009, NCEP reanal.','ERA40'),
       lwd=c(2,12,2,12,2,2,2,2),lty=c(1,1,1,1,1,1,1,1), 
       col=c('black','lightgrey','red','indianred1','blue','magenta','yellow','cyan'),
       cex=fs, bty='o', bg='white', box.col='white') 

# data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='Z100'),seas='win',anomper=(1901:1980),lw=linew,sca=T,yl='')
# abline(v=vf,col='black',lty=2)
# abline(h=0,col='black',lty=2)
# ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
# acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
# ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
# acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
# legend('bottomleft', c(paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',ecormean[1],ecormean[2],ecormean[3],ecormean[4]),paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',acormean[1],acormean[2],acormean[3],acormean[4]),paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],ecorbem[2],ecorbem[3],ecorbem[4]),paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],acorbem[2],acorbem[3],acorbem[4])),cex=fs, bty='o', bg='white', box.col='white')


# NAO
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='NAO.calc'),seas='win',anomper=(1901:1980),lw=linew,sca=T,yl='',plotbem=F,plotmem=T,plotech=T,pag=F)
abline(v=vf,col='black',lty=2)
abline(h=0,col='black',lty=2)
vdata <- cbind(stefind['NAO_20C'],stefind['NAO_REC'],stefind['NAO_NCEP'],stefind['NAO_ERA-40'])
yrv <- stefind[,1]
centerfac <- apply(vdata,2,mean,na.rm=T)
scalefac <- apply(vdata,2,sd,na.rm=T)
if (wlen > 1) {
  vdata <- runmean(scale(vdata,center=centerfac,scale=scalefac),wlen)
} else {
  vdata <- scale(vdata,center=centerfac,scale=scalefac)
}
lines(yrv,vdata[,1],ty='l',lwd=linew,lty=1,col=rgb(0,0,10,8,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,2],ty='l',lwd=linew,lty=1,col=rgb(10,0,10,8,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,3],ty='l',lwd=linew,lty=1,col=rgb(10,10,0,8,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,4],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,8,maxColorValue=10), ylab="",xlab="")
ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
legend('bottomleft', c(
  paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',ecormean[1],',',ecormean[2],',',ecormean[3],',',ecormean[4]),
  paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',acormean[1],',',acormean[2],',',acormean[3],',',acormean[4])),
#  paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],',',ecorbem[2],',',ecorbem[3],',',ecorbem[4]),
#  paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],',',acorbem[2],',',acorbem[3],',',acorbem[4])),
  cex=fs, bty='o', bg='white', box.col='white')


# PNA
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='PNA.calc'),seas='win',anomper=(1901:1980),lw=linew,sca=T,xa='s',yl='',plotbem=F,plotmem=T,plotech=T,pag=F)
abline(v=vf,col='black',lty=2)
abline(h=0,col='black',lty=2)
vdata <- cbind(stefind['PNA_20C'],stefind['PNA_REC'],stefind['PNA_NCEP'],stefind['PNA_ERA-40'])
yrv <- stefind[,1]
centerfac <- apply(vdata,2,mean,na.rm=T)
scalefac <- apply(vdata,2,sd,na.rm=T)
if (wlen > 1) {
  vdata <- runmean(scale(vdata,center=centerfac,scale=scalefac),wlen)
} else {
  vdata <- scale(vdata,center=centerfac,scale=scalefac)
}
lines(yrv,vdata[,1],ty='l',lwd=linew,lty=1,col=rgb(0,0,10,8,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,2],ty='l',lwd=linew,lty=1,col=rgb(10,0,10,8,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,3],ty='l',lwd=linew,lty=1,col=rgb(10,10,0,8,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,4],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,8,maxColorValue=10), ylab="",xlab="")
ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
legend('bottomleft', c(
  paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',ecormean[1],',',ecormean[2],',',ecormean[3],',',ecormean[4]),
  paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',acormean[1],',',acormean[2],',',acormean[3],',',acormean[4])),
#  paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],',',ecorbem[2],',',ecorbem[3],',',ecorbem[4]),
#  paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],',',acorbem[2],',',acorbem[3],',',acorbem[4])),
  cex=fs, bty='o', bg='white', box.col='white')

dev.off()

# Why is recon corr with ECHAM so high (NH, SH, GLO mean temp.)?

# DO Jonas plots for general skill!


} # end load_prepplot











 
# PLOT atm indices
#if (ind_anom) {
  # some ocean based indices only NOT in landonly version  
  load(paste0('../data/indices/',expname,'/indices',filenameext,'allts.Rdata'))
# load(paste0(prepplotdir,'indices'_anom_',filenameext,'allts.Rdata')) 
#load('../data/prepplot_season/indices_allts_landonly.Rdata')
#stefindmon <- read.table("/Users/joerg/Documents/climdata/indices/stefan/stefan_monthly_indices.txt",header=T)
#stefindyrmean <- aggregate(stefindmon[,c(6,9,12,15,18,21)],list(stefindmon[,1]),mean,na.rm=T)
stefind <- read.table(paste0(dataextdir,"indices/stefan/stefan_seasonal_indices.txt"),header=T,
                      skip=1,na.string=-9999)[,c(1,seq(2,5),seq(18,21),seq(34,37),seq(50,53),
                      seq(66,69),seq(82,85),seq(98,101),seq(114,117),seq(130,133),seq(146,149))]
header <- read.table(paste0(dataextdir,"indices/stefan/stefan_seasonal_indices.txt"),
                   nrow=2)[,c(1,seq(2,5),seq(18,21),seq(34,37),seq(50,53),seq(66,69),seq(82,85),
                   seq(98,101),seq(114,117),seq(130,133),seq(146,149))]
                   # (1,2,3,18,19,34,35,50,51,66,67,82,83,98,99,114,115,130,131,146,147)]
newhead <- rep(NA,ncol(header))
for (i in 1:ncol(header)) {
 newhead[i] <- paste(header[1,i],header[2,i],sep='_')
}
colnames(stefind) <- newhead
# #load('../data/indices_recon.Rdata')
# vf <- c(1600,1650,1660,1783,1809,1815,1835,1883,1886,1902,1912,1980,1991) #read.table('/Users/joerg/Documents/climdata/forcings/forc-total-crowley_2001.txt',skip=51,header=T)[,1:2]
#
nao_trou <- read.table(paste0(dataextdir,"indices/www/nao-trouet2009.txt"),skip=93,header=T)
nao_cook <- read.table(paste0(dataextdir,"indices/www/nao_cook2002.txt"),skip=50,header=F)
nao_lut_seas <- read.table(paste0(dataextdir,"indices/www/nao_sea_luterbacher_2002.txt"),skip=19,header=T)
nao_lut_mon <- read.table(paste0(dataextdir,"indices/www/nao_mon_luterbacher_2002.txt"),skip=27,header=T)
pna_trou <- read.table(paste0(dataextdir,"indices/www/pna-trouet2009.txt"),skip=75,header=F)
pna_ewen_mon <- read.table(paste0(dataextdir,"indices/www/pna_ewen_nnr.txt"),skip=0,header=T)
pna_ewen_yr <- aggregate(pna_ewen_mon[,3],by=list(pna_ewen_mon[,1]),mean)
# read new file from stefan because other version seem to have 20CRv1 and trouble with ewen data
nao_stefan_mon <- read.table(paste0(dataintdir,"../comparison_data/nao_stefan.txt"),skip=0,header=T)
nao_stefan_yr <- aggregate(nao_stefan_mon[,3:10],by=list(nao_stefan_mon[,1]),mean)
pna_stefan_mon <- read.table(paste0(dataintdir,"../comparison_data/pna_stefan_2017.txt"),skip=0,header=T)
pna_stefan_yr <- aggregate(pna_stefan_mon[,3:10],by=list(pna_stefan_mon[,1]),mean)

#pdf(file=paste0('../figures/nat_data_paper/nao_pna_index_v3',filenameext,'.pdf'),width=7,height=10) 
pdf(file=paste0('../figures/',expname,'/nao_pna_index_v3',filenameext,'.pdf'),width=7,height=10) 
wlen <- 1
linew <- 1
par(mfrow=c(4,1))
par(mai=c(0.3,0.6,0.2,0))
fs <- 1.0
anomwrt=(1961:1980)
#eind.allts$names

# NAO
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='NAO.calc'),seas='win',anomper=(1901:1980),lw=2,sca=T,xa='s',yl='',plotbem=F,plotmem=F,plotech=F,plotts=c(1,300),title='NAO')
#abline(v=vf,col='black',lty=2)
naols <- nao_lut_seas[nao_lut_seas[,2]=="wi",] # seas. averages until 1659
naols <- ts(naols[naols[,1]>=1600,3],start=1600)
pos <- nao_lut_mon[,2]=="Dec" | nao_lut_mon[,2]=="Jan" | 
       nao_lut_mon[,2]=="Feb" | nao_lut_mon[,2]=="Mar"
naolm <- nao_lut_mon[pos,]
pos2 <- which(naolm[,2]=="Dec")
naolm[pos2,1] <- (naolm[pos2,1]+1) # move Dec to next year to make winter season annual avg
naolm <- aggregate(naolm[,3],by=list(naolm[,1]),mean) # monthly data since 1659
naolm <- ts(naolm[,2],start=naolm[1,1],freq=1)
pos <- nao_stefan_mon[,2]==12 | nao_stefan_mon[,2]==2 | 
  nao_stefan_mon[,2]==2 | nao_stefan_mon[,2]==3
naosm <- nao_stefan_mon[pos,]
pos2 <- which(naosm[,2]==12)
naosm[pos2,1] <- (naosm[pos2,1]+1) # move Dec to next year to make winter season annual avg
naosm <- aggregate(naosm[,3:13],by=list(naosm[,1]),mean) # monthly data since 1659
nao20cr <- ts(naosm[1:104,4],start=naosm[1,1],freq=1)
naol <- ts(c(naols,naolm),start=tsp(naols)[1]) # annual winter NAO
naot <- ts(nao_trou[nao_trou[,1]>=1600,2],start=1600)
naoc <- ts(nao_cook[nao_cook[,1]>=1600,2],start=1600)
vdata=NULL
#vdata <- ts(cbind(stefind['NAO_20C'],stefind['NAO_REC'],stefind['NAO_NCEP'],stefind['NAO_ERA-40']),
#            start=stefind[1,1],freq=1)
vdata <- ts(cbind(stefind['NAO_REC'],stefind['NAO_NCEP'],stefind['NAO_ERA-40']),
           start=stefind[1,1],freq=1)
vdata <- ts.union(naol,naot,naoc,nao20cr,vdata)
yrv <- seq(1600,1600+nrow(vdata)-1)
centerfac <- apply(vdata,2,mean,na.rm=T)
scalefac <- apply(vdata,2,sd,na.rm=T)
vdata <- scale(vdata,center=centerfac,scale=scalefac)
if (wlen > 1) { vdata <- runmean(vdata,wlen) }
lines(vdata[,1],ty='l',lwd=linew,lty=1,col=rgb(0,0,10,7,maxColorValue=10), ylab="",xlab="")
lines(vdata[,2],ty='l',lwd=linew,lty=1,col=rgb(3,3,3,7,maxColorValue=10), ylab="",xlab="")
lines(vdata[,3],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,7,maxColorValue=10), ylab="",xlab="")
pos <- which(data[['period']]>1602 & data[['period']]<1901)
ecormean <- round(cor(data[['allind']][pos,1],window(vdata,1603,1900),
                      use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][pos,2],window(vdata,1603,1900),
                      use="pairwise.complete.obs"),digits=2) 

#ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
#acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
legend('topleft', c('EKF400','Luterbacher','Trouet','Cook'),
       lty=rep(1,5),lwd=c(2,1,1,1,1),
       col=c(rgb(10,0,0,10,maxColorValue=10),
             rgb(0,0,10,10,maxColorValue=10),rgb(3,3,3,10,maxColorValue=10),
             rgb(0,10,10,10,maxColorValue=10)),
       cex=fs, bty='o', bg='white', box.col='white')
legend('bottomleft',c('Correlation with Lut., Trou., Cook',
                      paste('CCC400 ens. mean:', ecormean[1],ecormean[2],ecormean[3]),
                      paste('EKF400 ens. mean:',acormean[1],acormean[2],acormean[3])),
       cex=fs, bty='o', bg='white', box.col='white')
#,paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],ecorbem[2],ecorbem[3],ecorbem[4]),paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],acorbem[2],acorbem[3],acorbem[4]))

data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='NAO.calc'),seas='win',anomper=anomwrt,
                   lw=2,sca=T,xa='s',yl='',plotbem=F,plotmem=F,plotech=F,
                   plotts=c(299,402),title='NAO 20th century')
lines(vdata[,4],ty='l',lwd=linew,lty=1,col=rgb(0,0,10,7,maxColorValue=10), ylab="",xlab="")
#lines(vdata[,5],ty='l',lwd=linew,lty=1,col=rgb(10,0,10,7,maxColorValue=10), ylab="",xlab="")
#lines(vdata[,6],ty='l',lwd=linew,lty=1,col=rgb(10,5,0,7,maxColorValue=10), ylab="",xlab="")
lines(vdata[,7],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,7,maxColorValue=10), ylab="",xlab="")
pos <- which(data[['period']]>1900 & data[['period']]<2001)
ecormean <- round(cor(data[['allind']][pos,1],window(vdata,1901,2000),
                      use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][pos,2],window(vdata,1901,2000),
                      use="pairwise.complete.obs"),digits=2) 
legend('topleft', c('EKF400','20CR','ERA40'), #'REC'),
       lty=rep(1,3),lwd=c(2,1,1),
       col=c(rgb(10,0,0,10,maxColorValue=10),
             rgb(0,0,10,10,maxColorValue=10),
             rgb(0,10,10,10,maxColorValue=10)),
             #rgb(10,0,10,10,maxColorValue=10))),
       cex=fs, bty='o', bg='white', box.col='white')
legend('bottomleft', c('Correlation with 20CR, ERA40', #REC',
    paste('CCC400 ens. mean:',ecormean[4],ecormean[7]), #ecormean[5]),
    paste('EKF400 ens. mean:',acormean[4],acormean[7])), #acormean[5])),
  cex=fs, bty='o', bg='white', box.col='white')
#  paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[4],ecorbem[5],ecorbem[3],ecorbem[6]),
#  paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[4],acorbem[5],acorbem[6],acorbem[7])


# PNA
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='PNA.calc'),seas='win',anomper=anomwrt,lw=2,sca=T,xa='s',yl='',plotbem=F,plotmem=F,plotech=F,title='PNA')
#abline(v=vf,col='black',lty=2)
#vdata <- cbind(stefind['PNA_20C'],stefind['PNA_REC'],stefind['PNA_NCEP'],stefind['PNA_ERA-40'])
pos <- pna_stefan_mon[,2]==12 | pna_stefan_mon[,2]==2 | 
  pna_stefan_mon[,2]==2 | pna_stefan_mon[,2]==3
pnasm <- pna_stefan_mon[pos,]
pos2 <- which(pnasm[,2]==12)
pnasm[pos2,1] <- (pnasm[pos2,1]+1) # move Dec to next year to make winter season annual avg
pnasm <- aggregate(pnasm[,3:10],by=list(pnasm[,1]),mean) # monthly data since 1659
pna20cr <- ts(pnasm[,3],start=pnasm[1,1],freq=1)
pnaewen <- ts(pnasm[,6],start=pnasm[1,1],freq=1)
pnabroe <- ts(pnasm[,7],start=pnasm[1,1],freq=1)
pnaera40 <- ts(pnasm[,9],start=pnasm[1,1],freq=1)
vdata <- cbind(pna20cr,pnaewen,pnabroe,pnaera40)
yrv <- pnasm[,1]
vdata <- ts(vdata,start=yrv[1],freq=1)
pnat <-ts(pna_trou[,2],start=pna_trou[1,1],freq=1)
#pnae <-ts(pna_ewen_yr[,2],start=pna_ewen_yr[1,1],freq=1)
vdata <- ts.union(vdata,pnat) #,pnae)
centerfac <- apply(vdata,2,mean,na.rm=T)
scalefac <- apply(vdata,2,sd,na.rm=T)
if (wlen > 1) {
  vdata <- runmean(scale(vdata,center=centerfac,scale=scalefac),wlen)
} else {
  vdata <- scale(vdata,center=centerfac,scale=scalefac)
}
#lines(vdata[,6],ty='l',lwd=linew,lty=1,col=rgb(3,3,3,7,maxColorValue=10), ylab="",xlab="")
lines(vdata[,5],ty='l',lwd=linew,lty=1,col=rgb(3,3,3,7,maxColorValue=10), ylab="",xlab="")
#lines(vdata[,1],ty='l',lwd=linew,lty=1,col=rgb(0,0,10,7,maxColorValue=10), ylab="",xlab="")
#lines(vdata[,2],ty='l',lwd=linew,lty=1,col=rgb(10,0,10,7,maxColorValue=10), ylab="",xlab="")
#lines(vdata[,3],ty='l',lwd=linew,lty=1,col=rgb(10,10,0,7,maxColorValue=10), ylab="",xlab="")
#lines(vdata[,4],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,7,maxColorValue=10), ylab="",xlab="")
pos <- which(data[['period']]>1900 & data[['period']]<2001)
ecormean <- round(cor(data[['allind']][pos,1],window(vdata[,1:5],1901,2000),
                      use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][pos,2],window(vdata[,1:5],1901,2000),
                      use="pairwise.complete.obs"),digits=2) 
pos2 <- which(data[['period']]>1724 & data[['period']]<2001)
ecormean[5] <- round(cor(data[['allind']][pos2,1],window(vdata[,5],1725,2000),
                      use="pairwise.complete.obs"),digits=2) 
acormean[5] <- round(cor(data[['allind']][pos2,2],window(vdata[,5],1725,2000),
                      use="pairwise.complete.obs"),digits=2) 
#legend('bottomleft', c('Correlation with 20CR, REC, NCEP, ERA40, Trouet, Ewen',
legend('bottomleft', c('Correlation with Trouet',paste('CCC400 ens. mean:',ecormean[5]),
  paste('EKF400 ens. mean:',acormean[5])),
  cex=fs, bty='o', bg='white', box.col='white')
#legend('topleft', c('EKF400','20CR','REC','NCEP','ERA40','Trouet','Ewen'),
legend('topleft', c('EKF400','Trouet'),
  lty=rep(1,6),lwd=c(2,1),col=c(rgb(10,0,0,10,maxColorValue=10),
  #rgb(0,0,10,10,maxColorValue=10),rgb(10,0,10,10,maxColorValue=10),
  #rgb(10,5,0,10,maxColorValue=10),rgb(0,10,10,10,maxColorValue=10),
  rgb(3,3,3,10,maxColorValue=10)), #,rgb(3,3,3,10,maxColorValue=10)),
  cex=fs, bty='o', bg='white', box.col='white')

data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='PNA.calc'),seas='win',anomper=anomwrt,
                   lw=2,sca=T,xa='s',yl='',plotbem=F,plotmem=F,plotech=F,
                   plotts=c(299,402),title='PNA 20th century')
#lines(vdata[,6],ty='l',lwd=linew,lty=1,col=rgb(3,3,3,7,maxColorValue=10), ylab="",xlab="")
lines(vdata[,1],ty='l',lwd=linew,lty=1,col=rgb(0,0,10,7,maxColorValue=10), ylab="",xlab="")
#lines(vdata[,2],ty='l',lwd=linew,lty=1,col=rgb(10,5,0,7,maxColorValue=10), ylab="",xlab="")
#lines(vdata[,3],ty='l',lwd=linew,lty=1,col=rgb(10,10,0,7,maxColorValue=10), ylab="",xlab="")
lines(vdata[,4],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,7,maxColorValue=10), ylab="",xlab="")
legend('bottomleft', c('Correlation with 20CR, ERA40', #Ewen',
  paste('CCC400 ens. mean:',ecormean[1],ecormean[4]), #ecormean[2]),
  paste('EKF400 ens. mean:',acormean[1],acormean[4])), #acormean[2])),
  cex=fs, bty='o', bg='white', box.col='white')
legend('topleft', c('EKF400','20CR','ERA40'), #'Ewen'),
  lty=rep(1,3),lwd=c(2,1,1),
  col=c(rgb(10,0,0,10,maxColorValue=10),rgb(0,0,10,10,maxColorValue=10),
        rgb(0,10,10,10,maxColorValue=10)), #rgb(10,5,0,10,maxColorValue=10)),
  cex=fs, bty='o', bg='white', box.col='white')
dev.off()








# PLOT atm indices new short version
  # some ocean based indices only NOT in landonly version
load(paste0('../data/indices/',expname,'/indices',filenameext,'allts.Rdata'))
#load(paste0(prepplotdir,'indices',filenameext,'allts.Rdata')) 
stefind <- read.table("/Users/joerg/Documents/climdata/indices/stefan/stefan_seasonal_indices.txt"
                      ,header=T,skip=1,na.string=-9999)[,c(1,seq(2,5),seq(18,21),seq(34,37),
                      seq(50,53),seq(66,69),seq(82,85),seq(98,101),seq(114,117),seq(130,133),
                      seq(146,149))]
head <- read.table("/Users/joerg/Documents/climdata/indices/stefan/stefan_seasonal_indices.txt",
                   nrow=2)[,c(1,seq(2,5),seq(18,21),seq(34,37),seq(50,53),seq(66,69),
                   seq(82,85),seq(98,101),seq(114,117),seq(130,133),seq(146,149))]
newhead <- rep(NA,ncol(head))
for (i in 1:ncol(head)) {
  newhead[i] <- paste(head[1,i],head[2,i],sep='_')
}
colnames(stefind) <- newhead

#pdf(file=paste0('../figures/indices/nao_pna_index_comp4',filenameext,'.pdf'),width=10,height=10) 
pdf(file=paste0('../figures/',expname,'/nao_pna_index_comp4',filenameext,'.pdf'),width=10,height=10) 
wlen <- 1
linew <- 2
par(mfrow=c(2,1))
par(mai=c(0.5,0.5,0.1,0.1))
fs <- 1.0
anomwrt=(1961:1990)

# NAO
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='NAO.calc'),seas='win',anomper=anomwrt,
                   lw=linew,sca=T,xa='s',yl='',plotbem=F,plotmem=F,plotts=c(299,402),title=' ')
vdata <- cbind(stefind['NAO_20C'],stefind['NAO_REC'],stefind['NAO_NCEP'],stefind['NAO_ERA-40'])
yrv <- stefind[,1]
centerfac <- apply(vdata,2,mean,na.rm=T)
scalefac <- apply(vdata,2,sd,na.rm=T)
if (wlen > 1) {
  vdata <- runmean(scale(vdata,center=centerfac,scale=scalefac),wlen)
} else {
  vdata <- scale(vdata,center=centerfac,scale=scalefac)
}
lines(yrv,vdata[,1],ty='l',lwd=linew,lty=1,col=rgb(0,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,2],ty='l',lwd=linew,lty=1,col=rgb(10,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,3],ty='l',lwd=linew,lty=1,col=rgb(10,5,0,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,4],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,10,maxColorValue=10), ylab="",xlab="")
ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],
                      use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],
                      use="pairwise.complete.obs"),digits=2) 
legend('bottomleft', c(paste('CCC400 ens. mean - 20CR, REC, NCEP, ERA40:',
                             ecormean[1],ecormean[2],ecormean[3],ecormean[4]),
                       paste('EKF400 ens. mean - 20CR, REC, NCEP, ERA40:',
                             acormean[1],acormean[2],acormean[3],acormean[4])),
                       cex=fs, bty='o', bg='white', box.col='white')
legend('topleft', c('CCC400 NAO','EKF400 NAO','20CR NAO','REC NAO','NCEP NAO','ERA40 NAO'),
       lty=rep(1,6),lwd=c(2,2,1,1,1,1),
       col=c(rgb(0,0,0,10,maxColorValue=10),rgb(10,0,0,10,maxColorValue=10),
             rgb(0,0,10,10,maxColorValue=10),rgb(10,0,10,10,maxColorValue=10),
             rgb(10,5,0,10,maxColorValue=10),rgb(0,10,10,10,maxColorValue=10)),
       cex=fs, bty='o', bg='white', box.col='white')

# PNA
data <- plot_pages(wl=wlen,reg=which(eind.allts$names=='PNA.calc'),seas='win',anomper=anomwrt,
                   lw=linew,sca=T,xa='s',yl='',plotbem=F,plotmem=F,plotts=c(299,402),title=' ')
vdata <- cbind(stefind['PNA_20C'],stefind['PNA_REC'],stefind['PNA_NCEP'],stefind['PNA_ERA-40'])
yrv <- stefind[,1]
centerfac <- apply(vdata,2,mean,na.rm=T)
scalefac <- apply(vdata,2,sd,na.rm=T)
if (wlen > 1) {
  vdata <- runmean(scale(vdata,center=centerfac,scale=scalefac),wlen)
} else {
  vdata <- scale(vdata,center=centerfac,scale=scalefac)
}
lines(yrv,vdata[,1],ty='l',lwd=linew,lty=1,col=rgb(0,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,2],ty='l',lwd=linew,lty=1,col=rgb(10,0,10,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,3],ty='l',lwd=linew,lty=1,col=rgb(10,5,0,10,maxColorValue=10), ylab="",xlab="")
lines(yrv,vdata[,4],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,10,maxColorValue=10), ylab="",xlab="")
ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],
                      use="pairwise.complete.obs"),digits=2) 
acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],
                      use="pairwise.complete.obs"),digits=2) 
legend('bottomleft', c(paste('CCC400 ens. mean - 20CR, REC, NCEP, ERA40:',
                             ecormean[1],ecormean[2],ecormean[3],ecormean[4]),
                       paste('EKF400 ens. mean - 20CR, REC, NCEP, ERA40:',
                             acormean[1],acormean[2],acormean[3],acormean[4])),
                       cex=fs, bty='o', bg='white', box.col='white')
legend('topleft', c('CCC400 PNA','EKF400 PNA','20CR PNA','REC PNA','NCEP PNA','ERA40 PNA'),
              lty=rep(1,6),lwd=c(2,2,1,1,1,1),
              col=c(rgb(0,0,0,10,maxColorValue=10),rgb(10,0,0,10,maxColorValue=10),
                    rgb(0,0,10,10,maxColorValue=10),rgb(10,0,10,10,maxColorValue=10),
                    rgb(10,5,0,10,maxColorValue=10),rgb(0,10,10,10,maxColorValue=10)),
              cex=fs, bty='o', bg='white', box.col='white')
dev.off()



} # end timeseriesplots









if (drought_1790) {
  ptm1 <- proc.time()
  # calc validate clim (70 yr period) to substract for anomalies
  for (cyr in seq(1760,1830)) {
    load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
    if (cyr==1760){
      vtmp <- validate$data
      c=1
    } else {
      vtmp <- vtmp+validate$data
      c <- c+1
    }
  }
  valclim <- vtmp/c
    
  # Asia PDSI comparison years 1790, 1792-96
#  for (cyr in c(1790, 1792, 1793, 1794, 1795, 1796)) {
# define drought period
  syr <- 1790
  eyr <- 1799

  for (cyr in seq(syr,eyr)) {
    load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
    echam.abs <- echam
    analysis.abs <- analysis
    bemlist <- read.table(file='../data/bem/bem.txt',header=F)
    bem <- bemlist[which(bemlist[,1]==cyr),2]
    echam.abs$bem <- echam.abs$data[,,bem]
    echam.anom$bem <- echam.anom$data[,,bem]
    analysis.abs$bem <- analysis.abs$data[,,bem]
    analysis.anom$bem <- analysis.anom$data[,,bem]
#    ech_ind$bem <- ech_ind$data[,,bem]
#    ana_ind$bem <- ana_ind$data[,,bem]
    
    if (cyr==syr){
      atmp <- analysis.anom$ensmean
      etmp <- echam.anom$ensmean
      a2tmp <- analysis.abs$ensmean
      e2tmp <- echam.abs$ensmean
      a3tmp <- analysis.anom$bem
      e3tmp <- echam.anom$bem
      vtmp <- validate$data
      c=1
    } else {
      atmp <- atmp+analysis.anom$ensmean
      etmp <- etmp+echam.anom$ensmean
      a2tmp <- a2tmp+analysis.abs$ensmean
      e2tmp <- e2tmp+echam.abs$ensmean
      a3tmp <- a3tmp+analysis.anom$bem
      e3tmp <- e3tmp+echam.anom$bem
      vtmp <- vtmp+validate$data
      c=c+1
    }
  }
  anameananom <- atmp/c
  echmeananom <- etmp/c
  eadiff <- anameananom-echmeananom
  anameanabs <- a2tmp/c
  echmeanabs <- e2tmp/c
  eadiffabs <- anameanabs-echmeanabs
  anameananombem <- a3tmp/c
  echmeananombem <- e3tmp/c
  eadiffbem <- anameananombem-echmeananombem
  valmean <- vtmp/c
  valmeananom <- valmean-valclim

  pnames <- paste(c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k'), ')', sep='')
  data.dim <- dim(echam$data)[c(1,2,2,3)]
  data.dim[2:3] <- c(s,data.dim[3]/s)
  ens.dim <- c(nrow(echam$ensmean), s, ncol(echam$ensmean)/s)



  # open cook pdsi recons
  # first interpolate to echam global grid:
  # cdo -m 'NaNf' -setmissval,9e33 MADApdsi.nc madapdsi2.nc
  # cdo -r remapbil,t63grid madapdsi2.nc mada_pdsi_echamgrid.nc
  # cdo -m 'NaNf' -setmissval,9e33 NADAv2-2008.nc nadapdsi2.nc
  # cdo -r remapbil,t63grid nadapdsi2.nc nada_pdsi_echamgrid.nc
  mada <- read_pdsi('mada_pdsi_echamgrid.nc',timlim=c(syr,eyr), 
                      path="/Users/joerg/Documents/climdata/recon/cook/pdsi_asia/",
                      small=T,landonly=F)
  nada <- read_pdsi('nada_pdsi_echamgrid.nc',timlim=c(syr,eyr), 
                  path="/Users/joerg/Documents/climdata/recon/cook/pdsi_usa/",
                  small=T,landonly=F)
  nada$data[is.na(nada$data)] <- 0
  mada$data[is.na(mada$data)] <- 0
  lut <- valmeananom[(which(validate$names=='precip')),2]
  lut[is.na(lut)] <- 0
#  nada$data <- scale(nada$data,scale=sd(lut))
#  mada$data <- scale(mada$data,scale=sd(lut))
  pdsi <- nada
  # scale factor to make precip and pdsi comparable
  f=2
  tmp <- nada$data*f + mada$data*f + lut
  tmp[tmp==0] <- NA
  pdsi$data <- apply(tmp[,,1],1,mean,na.rm=T)
  # erst 4608 temp und dann erst precip
  pdsi$data <- rep(pdsi$data,11)
  pdsi$names <- rep('precip',length(pdsi$names))


  
# Figure 2
  pdata <- echam
  ti=2 # summer season as drought recon is JJA
  pdf(paste0('../figures/',expname,'/drought_JJA_',syr,'-',eyr,'_',filenameext,'.pdf'), width=6, height=7, paper='special')
#  layout(matrix(c(1,2,3,3,4,5,6,6,7,8,9,9), 6, 2, byrow = TRUE), height=c(3,1,3,1,3,1))
  layout(matrix(c(1,2,3,3,4,4,5,6,7,7), 5, 2, byrow = TRUE), height=c(3,1,1,3,1))
#  layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
  par(oma=c(0.5,4,4,0))
  levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
  contlevs <- seq(-4,4,2)
  pdata$data <- array(c(anameananom[,ti],eadiff[,ti]), c(nrow(echam.anom$ensmean),1,2))
  plot_echam(pdata, varname='temp2', type='data', cex.pt=1.5, latlim=c(-10,90),wcol='darkgrey',
             lev=levs, st.col=NULL, stations=calibrate, names=c('a)', 'b)'),
             colnames=c("EKF400 anomaly","EKF400 update"),
             statpanel=2, add=T, rownames='temperature and GHP500', main='1790-1799',
             addcontours=T, contvarname='gph500',conttype='data',contcol='black',
             contlev=contlevs)
  plot(NA,xlim=c(0,1),ylim=c(0,1),axes=F,xlab='',ylab='')
  levs <- c(-Inf, seq(-2,2,0.5), Inf)
  pdata$data <- array(c(anameananom[,ti],pdsi$data), c(nrow(echam.anom$ensmean),1,2))
  plot_echam(pdata, varname='precip', type='data', cex.pt=1.5, latlim=c(-10,90),
             names=c('c)', 'd)'), lev=levs, st.col=NULL, stations=NULL, add=T,
             addcontours=F,wcol='darkgrey',rownames='precipitation/drought',
             colnames=c("EKF400 anomaly","reconstructions"))
  
  dev.off()
#   
#   # read mann et al temp recon and plot 1816
#   mann <- read.table("/Users/joerg/Documents/climdata/mann_temp_recon/allproxyfieldrecon",
#                      na.string="NaN")
#   ll <- read.table("/Users/joerg/Documents/climdata/mann_temp_recon/longlat.txt",
#                    na.string="NaN")
#   ts <- which(mann[,1]==1816)
#   plotmann <- plotdata
#   mdata <- as.vector(t(mann[ts,-1]))
#   #plotmann$data <- array(rep(mdata,2),c(length(mdata),1,2))
#   plotmann$data <- array(mdata,c(length(mdata),1))
#   #plotmann$data <- array(mdata,,c(length(mdata),1,1))
#   #plotmann$ensmean <- array(mdata,c(length(mdata),1,1))
#   #plotmann$time <- rep(mann[ts,1],2)
#   plotmann$time <- mann[ts,1]
#   plotmann$names <- rep('temp2',dim(plotmann$data)[1])
#   plotmann$lon <- ll[,1] 
#   plotmann$lat <- ll[,2]
#   plotmann$height <- ll[,3]   
#   plotmann$lsm.i <- ll[,3]  
#   layout(matrix(c(1,2,3,3), 2, 1, byrow = TRUE), height=c(3,1))
#   #layout(matrix(c(1,2,2), 2, 2, byrow = TRUE), height=c(3,1))
#   par(oma=c(0,0,0,0))
#   levs <- c(-Inf, seq(-2,2,0.5), Inf)
#   plot_echam(plotmann, varname='temp2', type='data', cex.pt=1.5, names='',
#               lev=levs, st.col=NULL, stations=NULL, add=T,units='ºC')

  
  
  
  
  
  
  
  
  
  # make plots
  ti=2 # summer season as drought recon is JJA
  plotdata=echam
  #  plotdata$data <- array(c(echmeananom[,ti],anameananom[,ti],eadiff[,ti],valmeananom[,ti]), c(nrow(echmeananom),1,4))
  plotdata$data <- array(c(anameananom[,ti],valmeananom[,ti]), c(nrow(echmeananom),1,2))  
  pdf(paste0('../figures/',expname,'/drought_',syr,'-',eyr,filenameext,'.pdf'), width=9, height=9, paper='special')
  layout(matrix(c(1,2,3,3,4,5,6,6), 4, 2, byrow = TRUE), height=c(3,1,3,1))
  par(oma=c(0,0,0,0))
  levs <- c(-Inf, seq(-8,8,2), Inf)
  plot_echam(plotdata, varname='precip', type='data', cex.pt=1.5, names=pnames[1:dim(plotdata$data)[3]],
             lev=levs, st.col=NULL, stations=NULL, add=T) #, addcontours=T, contvarname='gph500')
  
#  add_contour(plotdata, varname = "temp2", ti=1)
  plotdata$data <- array(c(eadiff[,ti],pdsi$data), c(nrow(echmeananom),1,2))
  levs <- c(-Inf, seq(-2,2,0.5), Inf)
  plot_echam(plotdata, varname='precip', type='data', cex.pt=1.5, names=pnames[3:4],
             lev=levs, st.col=NULL, stations=NULL, add=T)
  dev.off()


# search explanation for droughts: do HC index, omega500 and gph500 anomaly indicate subtropical subsidence?
  ti=2 # summer season as drought recon is JJA
  plotdata$data <- array(c(echmeananom[,ti],anameananom[,ti]), c(nrow(echam.anom$ensmean),1,2))
  pdf(paste0('../figures/',expname,'/gph500_anom_',syr,'-',eyr,filenameext,'.pdf'), width=9, height=4.5, paper='special')
  layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
  par(oma=c(0,0,0,0))
  levs <- c(-Inf, seq(-8,8,2), Inf)
  plot_echam(plotdata, varname='gph500', type='data', cex.pt=1.5, 
             names=pnames[1:dim(plotdata$data)[3]], lev=levs, st.col=NULL, stations=NULL, add=T) #,
#             addcontours=T, contvarname='gph500',conttype='data',contcol='white')
  dev.off()

  plotdata$data <- array(c(echmeananombem[,ti],anameananombem[,ti]), c(nrow(echam.anom$ensmean),1,2))
  pdf(paste0('../figures/',expname,'/gph500_anom_bem_',syr,'-',eyr,filenameext,'.pdf'), width=9, height=4.5, paper='special')
  layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
  par(oma=c(0,0,0,0))
  levs <- c(-Inf, seq(-12,12,2), Inf)
  plot_echam(plotdata, varname='gph500', type='data', cex.pt=1.5, 
           names=pnames[1:dim(plotdata$data)[3]], lev=levs, st.col=NULL, 
           stations=NULL, add=T)
  dev.off()

  plotdata$data <- array(c(echmeananom[,ti],anameananom[,ti]), c(nrow(echam.anom$ensmean),1,2))
  pdf(paste0('../figures/',expname,'/temp_anom_',syr,'-',eyr,filenameext,'.pdf'), width=9, height=4.5, paper='special')
  layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
  par(oma=c(0,0,0,0))
  levs <- c(-Inf, seq(-1,1,0.1), Inf)
  plot_echam(plotdata, varname='temp2', type='data', cex.pt=1.5, 
           names=pnames[1:dim(plotdata$data)[3]], lev=levs, st.col=NULL, 
           stations=NULL, add=T) #, addcontours=T, contvarname='gph500')
  dev.off()

  plotdata$data <- array(c(echmeananom[,ti],anameananom[,ti]), c(nrow(echam.anom$ensmean),1,2))
  pdf(paste0('../figures/',expname,'/slp_anom_',syr,'-',eyr,filenameext,'.pdf'), width=9, height=4.5, paper='special')
  layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
  par(oma=c(0,0,0,0))
  levs <- c(-Inf, seq(-1,1,0.2), Inf)
  plot_echam(plotdata, varname='slp', type='data', cex.pt=1.5, 
           names=pnames[1:dim(plotdata$data)[3]], lev=levs, st.col=NULL, 
           stations=NULL, add=T)
  dev.off()

  plotdata$data <- array(c(eadiff[,ti],anameananom[,ti]), c(nrow(echam.anom$ensmean),1,2))
  pdf(paste0('../figures/',expname,'/omega500_anom_',syr,'-',eyr,filenameext,'.pdf'), width=9, height=4.5, paper='special')
  layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
  par(oma=c(0,0,0,0))
  levs <- c(-Inf, seq(-0.001,0.001,0.0001), Inf)
  plot_echam(plotdata, varname='omega500', type='data', cex.pt=1.5, 
           names=pnames[1:dim(plotdata$data)[3]], lev=levs, st.col=NULL, 
           stations=NULL, add=T)
  dev.off()

} # end 1790 droughts





#if (cold_1640) {}
#if (cold_1690) {} 

if (warmcold_decades) {
  ptm1 <- proc.time()
#  syr <- 1640
#  eyr <- 1649
#  syr <- 1641
#  eyr <- 1645
  syr <- 1641
  eyr <- 1648
  
  for (cyr in seq(syr,eyr)) {
    load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
    if (cyr==syr){
      atmp <- analysis.anom$ensmean
      etmp <- echam.anom$ensmean
      c=1
    } else {
      atmp <- atmp+analysis.anom$ensmean
      etmp <- etmp+echam.anom$ensmean
      c=c+1
    }
  }
#  anameananom1640 <- atmp/c
#  echmeananom1640 <- etmp/c
  anameananom1641 <- atmp/c
  echmeananom1641 <- etmp/c

#  syr <- 1658
#  eyr <- 1662
  syr <- 1810
  eyr <- 1817

  for (cyr in seq(syr,eyr)) {
    load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
    if (cyr==syr){
      atmp <- analysis.anom$ensmean
      etmp <- echam.anom$ensmean
      c=1
    } else {
      atmp <- atmp+analysis.anom$ensmean
      etmp <- etmp+echam.anom$ensmean
      c=c+1
    }
  }
#  anameananom1658 <- atmp/c
#  echmeananom1658 <- etmp/c
  anameananom1810 <- atmp/c
  echmeananom1810 <- etmp/c

  
#  syr <- 1660
#  eyr <- 1669
  syr <- 1695
#  eyr <- 1699
  eyr <- 1702
  
  for (cyr in seq(syr,eyr)) {
    load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
    if (cyr==syr){
      atmp <- analysis.anom$ensmean
      etmp <- echam.anom$ensmean
      c=1
    } else {
      atmp <- atmp+analysis.anom$ensmean
      etmp <- etmp+echam.anom$ensmean
      c=c+1
    }
  }
#  anameananom1660 <- atmp/c
#  echmeananom1660 <- etmp/c
  anameananom1695 <- atmp/c
  echmeananom1695 <- etmp/c
  
#  syr <- 1690
#  eyr <- 1699
#  syr <- 1676
#  eyr <- 1680
  syr <- 1831
  eyr <- 1838
  
  for (cyr in seq(syr,eyr)) {
    load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
    if (cyr==syr){
      atmp <- analysis.anom$ensmean
      etmp <- echam.anom$ensmean
      c=1
    } else {
      atmp <- atmp+analysis.anom$ensmean
      etmp <- etmp+echam.anom$ensmean
      c=c+1
    }
  }
#  anameananom1690 <- atmp/c
#  echmeananom1690 <- etmp/c
#  anameananom1676 <- atmp/c
#  echmeananom1676 <- etmp/c
  anameananom1831 <- atmp/c
  echmeananom1838 <- etmp/c
  
#  syr <- 1790
#  eyr <- 1799
#  syr <- 1791
#  eyr <- 1795
  syr <- 1655
  eyr <- 1662
  
  for (cyr in seq(syr,eyr)) {
    load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
    if (cyr==syr){
      atmp <- analysis.anom$ensmean
      etmp <- echam.anom$ensmean
      c=1
    } else {
      atmp <- atmp+analysis.anom$ensmean
      etmp <- etmp+echam.anom$ensmean
      c=c+1
    }
  }
#  anameananom1790 <- atmp/c
#  echmeananom1790 <- etmp/c
#  anameananom1791 <- atmp/c
#  echmeananom1791 <- etmp/c
  anameananom1655 <- atmp/c
  echmeananom1655 <- etmp/c
  
#  syr <- 1800
#  eyr <- 1809
#  syr <- 1801
#  eyr <- 1805
  syr <- 1798
  eyr <- 1805
  
  for (cyr in seq(syr,eyr)) {
    load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
    if (cyr==syr){
      atmp <- analysis.anom$ensmean
      etmp <- echam.anom$ensmean
      c=1
    } else {
      atmp <- atmp+analysis.anom$ensmean
      etmp <- etmp+echam.anom$ensmean
      c=c+1
    }
  }
#  anameananom1800 <- atmp/c
#  echmeananom1800 <- etmp/c
#  anameananom1801 <- atmp/c
#  echmeananom1801 <- etmp/c
  anameananom1798 <- atmp/c
  echmeananom1798 <- etmp/c
  
#  syr <- 1810
#  eyr <- 1819
#  syr <- 1815
#  eyr <- 1819
  syr <- 1774
  eyr <- 1781
  
  for (cyr in seq(syr,eyr)) {
    load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
    if (cyr==syr){
      atmp <- analysis.anom$ensmean
      etmp <- echam.anom$ensmean
      c=1
    } else {
      atmp <- atmp+analysis.anom$ensmean
      etmp <- etmp+echam.anom$ensmean
      c=c+1
    }
  }
#  anameananom1810 <- atmp/c
#  echmeananom1810 <- etmp/c
#  anameananom1815 <- atmp/c
#  echmeananom1815 <- etmp/c
  anameananom1774 <- atmp/c
  echmeananom1774 <- etmp/c
  
#  syr <- 1820
#  eyr <- 1829
#  syr <- 1826
#  eyr <- 1830
  syr <- 1822
  eyr <- 1829
  
  for (cyr in seq(syr,eyr)) {
    load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
    if (cyr==syr){
      atmp <- analysis.anom$ensmean
      etmp <- echam.anom$ensmean
      c=1
    } else {
      atmp <- atmp+analysis.anom$ensmean
      etmp <- etmp+echam.anom$ensmean
      c=c+1
    }
  }
#  anameananom1820 <- atmp/c
#  echmeananom1820 <- etmp/c
#  anameananom1826 <- atmp/c
#  echmeananom1826 <- etmp/c
  anameananom1822 <- atmp/c
  echmeananom1822 <- etmp/c
  
#  syr <- 1830
#  eyr <- 1839
#  syr <- 1832
#  eyr <- 1836
  syr <- 1847
  eyr <- 1854
  
  for (cyr in seq(syr,eyr)) {
    load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
    if (cyr==syr){
      atmp <- analysis.anom$ensmean
      etmp <- echam.anom$ensmean
      c=1
    } else {
      atmp <- atmp+analysis.anom$ensmean
      etmp <- etmp+echam.anom$ensmean
      c=c+1
    }
  }
#  anameananom1830 <- atmp/c
#  echmeananom1830 <- etmp/c
#  anameananom1832 <- atmp/c
#  echmeananom1832 <- etmp/c
  anameananom1847 <- atmp/c
  echmeananom1847 <- etmp/c
  
  syr <- 1686
  eyr <- 1693
  
  for (cyr in seq(syr,eyr)) {
    load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
    if (cyr==syr){
      atmp <- analysis.anom$ensmean
      etmp <- echam.anom$ensmean
      c=1
    } else {
      atmp <- atmp+analysis.anom$ensmean
      etmp <- etmp+echam.anom$ensmean
      c=c+1
    }
  }
  anameananom1686 <- atmp/c
  echmeananom1686 <- etmp/c
  
  pnames <- paste(c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k'), ')', sep='')
  data.dim <- dim(echam$data)[c(1,2,2,3)]
  data.dim[2:3] <- c(s,data.dim[3]/s)
  ens.dim <- c(nrow(echam$ensmean), s, ncol(echam$ensmean)/s)
  
# # convert monthly_out from 6-mon state vector because of now fixed bug in EnSRF_prepplot
#  echmeananom1791 <- array(echmeananom1791,c((dim(echmeananom1791)[1]/6),
#                       dim(echmeananom1791)[2]*6))
#  anameananom1791 <- array(anameananom1791,c((dim(anameananom1791)[1]/6),
#                                             dim(anameananom1791)[2]*6))


  # 8-yr cold period sommer
  pdata <- echam
  ti=2 # 
  pdf('../figures/',expname,'/decadel_cold_period_temp_anom_som.pdf', 
      width=12, height=7, paper='special')
    layout(matrix(c(1,2,3,4,5,5,5,5,6,7,8,9,10,10,10,10), 4, 4, byrow = TRUE), height=c(3,1,3,1))
    par(oma=c(0.5,4,4,0))
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    pdata$data <- array(c(anameananom1641[,ti],anameananom1695[,ti],anameananom1810[,ti],
                        anameananom1831[,ti]), c(nrow(echam.anom$ensmean),1,4))
    plot_echam(pdata, varname='temp2', type='data', cex.pt=1.5, latlim=c(-90,90),wcol='darkgrey',
             lev=levs, st.col=NULL, stations=calibrate, names=c('a)', 'b)', 'c)', 'd)'),
             colnames=c("EKF400 1641-48","EKF400 1695-1702","EKF400 1810-17","EKF400 1831-38"),
             statpanel=NULL, add=T, rownames='temperature and GHP500', main='Cold 8-yr periods',
             addcontours=T, contvarname='gph500',conttype='data',contcol='black',
             contlev=contlevs)
    levs <- c(-Inf, seq(-5,5,1), Inf)
    plot_echam(pdata, varname='precip', type='data', cex.pt=1.5, latlim=c(-90,90),
             names=c('e)', 'f)','g)', 'h)'), lev=levs, st.col=NULL, stations=NULL, add=T,
             addcontours=F,wcol='darkgrey',rownames='precipitation',colnames=rep('',4))
  
  dev.off()
  
  # 8-yr cold period winter
  pdata <- echam
  ti=1 # see volc high lat winter warming
  pdf('../figures/',expname,'/decadel_cold_period_temp_anom_win.pdf', 
      width=12, height=7, paper='special')
    layout(matrix(c(1,2,3,4,5,5,5,5,6,7,8,9,10,10,10,10), 4, 4, byrow = TRUE), height=c(3,1,3,1))
    par(oma=c(0.5,4,4,0))
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    pdata$data <- array(c(anameananom1641[,ti],anameananom1695[,ti],anameananom1810[,ti],
                          anameananom1831[,ti]), c(nrow(echam.anom$ensmean),1,4))
    plot_echam(pdata, varname='temp2', type='data', cex.pt=1.5, latlim=c(-90,90),wcol='darkgrey',
               lev=levs, st.col=NULL, stations=calibrate, names=c('a)', 'b)', 'c)', 'd)'),
               colnames=c("EKF400 1641-48","EKF400 1695-1702","EKF400 1810-17","EKF400 1831-38"),
               statpanel=NULL, add=T, rownames='temperature and GHP500', main='Cold 8-yr periods',
               addcontours=T, contvarname='gph500',conttype='data',contcol='black',
               contlev=contlevs)
    levs <- c(-Inf, seq(-5,5,1), Inf)
    plot_echam(pdata, varname='precip', type='data', cex.pt=1.5, latlim=c(-90,90),
             names=c('e)', 'f)','g)', 'h)'), lev=levs, st.col=NULL, stations=NULL, add=T,
             addcontours=F,wcol='darkgrey',rownames='precipitation',colnames=rep('',4))
  dev.off()
  
  
  # 8-yr warm period summer
  pdata <- echam
  ti=2 # summer season as drought recon is JJA
  pdf('../figures/',expname,'/decadel_warm_period_temp_anom_som.pdf', 
      width=12, height=7, paper='special')
    layout(matrix(c(1,2,3,4,5,5,5,5,6,7,8,9,10,10,10,10), 4, 4, byrow = TRUE), height=c(3,1,3,1))
    par(oma=c(0.5,4,4,0))
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    pdata$data <- array(c(anameananom1655[,ti],anameananom1686[,ti],anameananom1774[,ti],
                          anameananom1798[,ti]), c(nrow(echam.anom$ensmean),1,4))
    plot_echam(pdata, varname='temp2', type='data', cex.pt=1.5, latlim=c(-90,90),wcol='darkgrey',
             lev=levs, st.col=NULL, stations=calibrate, names=c('a)', 'b)', 'c)', 'd)'),
             colnames=c("EKF400 1655-62","EKF400 1686-93","EKF400 1774-81","EKF400 1798-1805"),
             statpanel=NULL, add=T, rownames='temperature and GHP500', main='Warm 5-yr periods',
             addcontours=T, contvarname='gph500',conttype='data',contcol='black',
             contlev=contlevs)
    levs <- c(-Inf, seq(-5,5,1), Inf)
    plot_echam(pdata, varname='precip', type='data', cex.pt=1.5, latlim=c(-90,90),
             names=c('e)', 'f)','g)', 'h)'), lev=levs, st.col=NULL, stations=NULL, add=T,
             addcontours=F,wcol='darkgrey',rownames='precipitation',colnames=rep('',4))
  
  dev.off()
  
  # 8-yr warm period winter
  pdata <- echam
  ti=1 
  pdf('../figures/',expname,'/decadel_warm_period_temp_anom_win.pdf', 
      width=12, height=7, paper='special')
    layout(matrix(c(1,2,3,4,5,5,5,5,6,7,8,9,10,10,10,10), 4, 4, byrow = TRUE), height=c(3,1,3,1))
    par(oma=c(0.5,4,4,0))
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    pdata$data <- array(c(anameananom1655[,ti],anameananom1686[,ti],anameananom1774[,ti],
                          anameananom1798[,ti]), c(nrow(echam.anom$ensmean),1,4))
    plot_echam(pdata, varname='temp2', type='data', cex.pt=1.5, latlim=c(-90,90),wcol='darkgrey',
               lev=levs, st.col=NULL, stations=calibrate, names=c('a)', 'b)', 'c)', 'd)'),
               colnames=c("EKF400 1655-62","EKF400 1686-93","EKF400 1774-81","EKF400 1798-1805"),
               statpanel=NULL, add=T, rownames='temperature and GHP500', main='Warm 5-yr periods',
               addcontours=T, contvarname='gph500',conttype='data',contcol='black',
               contlev=contlevs)
    levs <- c(-Inf, seq(-5,5,1), Inf)
    plot_echam(pdata, varname='precip', type='data', cex.pt=1.5, latlim=c(-90,90),
             names=c('e)', 'f)','g)', 'h)'), lev=levs, st.col=NULL, stations=NULL, add=T,
             addcontours=F,wcol='darkgrey',rownames='precipitation',colnames=rep('',4))
  
  dev.off()

#     
#   # Figure 2e SON, DJF, MAM, JJA to see if better Z500 signal in winter
#   pdata <- echam
#   ti=2 # look at winters of 5-yr warm periods
#   #  pdf('../figures/figure_2c_v3.pdf', width=12, height=7, paper='special')
#   pdf('../figures/figure_2c_v4.pdf', width=12, height=7, paper='special')
#   #  layout(matrix(c(1,2,3,3,4,5,6,6,7,8,9,9), 6, 2, byrow = TRUE), height=c(3,1,3,1,3,1))
#   # layout(matrix(c(1,2,3,3,4,4,5,6,7,7), 5, 2, byrow = TRUE), height=c(3,1,1,3,1))
#   layout(matrix(c(1,2,3,4,5,5,5,5,6,7,8,9,10,10,10,10), 4, 4, byrow = TRUE), height=c(3,1,3,1))
#   #  layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
#   par(oma=c(0.5,4,4,0))
#   levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
#   contlevs <- seq(-4,4,2)
#   #  pdata$data <- array(c(anameananom1660[,ti],anameananom1790[,ti],anameananom1800[,ti],anameananom1820[,ti]), c(nrow(echam.anom$ensmean),1,4))
#   pdata$data <- array(c(anameananom1658[,ti],anameananom1791[,ti],anameananom1801[,ti],anameananom1826[,ti]), c(nrow(echam.anom$ensmean),1,4))
#   
#   plot_echam(pdata, varname='temp2', type='data', cex.pt=1.5, latlim=c(-90,90),wcol='darkgrey',
#              lev=levs, st.col=NULL, stations=calibrate, names=c('a)', 'b)', 'c)', 'd)'),
#              colnames=c("EKF400 1658-62","EKF400 1791-95","EKF400 1801-05","EKF400 1826-30"),
#              statpanel=NULL, add=T, rownames='temperature and GHP500', main='Warm 5-yr periods',
#              addcontours=T, contvarname='gph500',conttype='data',contcol='black',
#              contlev=contlevs)
#   #  plot(NA,xlim=c(0,1),ylim=c(0,1),axes=F,xlab='',ylab='')
#   levs <- c(-Inf, seq(-5,5,1), Inf)
#   #  pdata$data <- array(c(anameananom1640[,ti],anameananom1690[,ti]), c(nrow(echam.anom$ensmean),1,2))
#   plot_echam(pdata, varname='precip', type='data', cex.pt=1.5, latlim=c(-90,90),
#              names=c('e)', 'f)','g)', 'h)'), lev=levs, st.col=NULL, stations=NULL, add=T,
#              addcontours=F,wcol='darkgrey',rownames='precipitation',colnames=rep('',4))
#   
#   dev.off()  
#}
} # end warmcold_decades









if (plots1816) {
  # check yuri's 1816 data compilation
  #load("../comparison_data/yuri_1816_all_station_and_travel_pressure_data.Rdata")
  syr=1816     
  eyr=1817
  
  ptm1 <- proc.time()
  # calc validate clim (70 yr period) to substract for anomalies
  for (cyr in seq(1781,1851)) {
    load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
    if (cyr==1781){
      vtmp <- validate$data
      c=1
    } else {
      vtmp <- vtmp+validate$data
      c <- c+1
    }
  }
  valclim <- vtmp/c

  # European T, P and SLP comparison year 1816
  cyr = 1816
  load(file=paste0(prepplotdirseas,'analysis_',cyr,'.Rdata'))
  #load(file=paste0(prepplotdirmon,'analysis_',cyr,'.Rdata'))
  validate.anom <- echam
  validate.anom$data <- echam$ensmean
  validate.anom$data[(1:nrow(validate$data)),] <- validate$data-valclim
  validate.anom$data[(nrow(validate$data)+1):nrow(echam$data),] <- NA
  
  pnames <- paste(c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k'), ')', sep='')
  data.dim <- dim(echam$data)[c(1,2,2,3)]
  data.dim[2:3] <- c(s,data.dim[3]/s)
  ens.dim <- c(nrow(echam$ensmean), s, ncol(echam$ensmean)/s)
  
  # load station location where data was assimilated
  pcoor <- read.table("../data/coor/prox_coor_1816.5.csv")
  scoor <- read.table("../data/coor/stat_coor_1816.5.csv")
  station <- list(lon=c(scoor[,1],pcoor[,1]),lat=c(scoor[,2],pcoor[,2]),
                  data=matrix(c(rep(1,nrow(pcoor)),rep(2,nrow(scoor))),ncol=1))
  # 
  # # make plots
  # #echam.anom$data <- array(echam.anom$data,c((dim(echam.anom$data)[1]/6),dim(echam.anom$data)[2]*6,
  # #                       dim(echam.anom$data)[3]))
  # echam.anom$ensmean <- array(echam.anom$ensmean,c((dim(echam.anom$ensmean)[1]/6),
  #                                                  dim(echam.anom$ensmean)[2]*6))
  # tmptime <- seq(1815,1817,by=1/12)
  # echam.anom$time <- tmptime[(season[2]+1):(length(tmptime)-4)]
  # echam.anom$names <- echam.anom$names[1:dim(echam.anom$data)[1]]
  # echam.anom$lon <- echam.anom$lon[1:(dim(echam.anom$data)[1])] #/length(unique(echam.anom$names)))]
  # echam.anom$lat <- echam.anom$lat[1:(dim(echam.anom$data)[1])] #/length(unique(echam.anom$names)))]
  # #analysis.anom$data <- array(analysis.anom$data,c((dim(analysis.anom$data)[1]/6),
  # #                         dim(analysis.anom$data)[2]*6,dim(analysis.anom$data)[3]))
  # analysis.anom$ensmean <- array(analysis.anom$ensmean,c((dim(analysis.anom$ensmean)[1]/6),
  #                                                        dim(analysis.anom$ensmean)[2]*6))
  # analysis.anom$time <- tmptime[(season[2]+1):(length(tmptime)-4)]
  # analysis.anom$names <- analysis.anom$names[1:dim(analysis.anom$data)[1]]
  # analysis.anom$lon <- analysis.anom$lon[1:(dim(analysis.anom$data)[1])] 
  # analysis.anom$lat <- analysis.anom$lat[1:(dim(analysis.anom$data)[1])]
  tim=2 # summer season as we look at year without a summer
#  for (tim in seq(10,12)) {
    pdata=echam
    #pdata$lon <- pdata$lon[1:length(which(validate.anom$names=="temp2"))]
    #pdata$lat <- pdata$lat[1:length(which(validate.anom$names=="temp2"))]
    #pdata$names <- pdata$names[1:length(which(validate.anom$names=="temp2"))]
#    pdf(paste0('../figures/nat_data_paper/eu_1816_',tim,filenameext,'.pdf'), width=9, height=9, paper='special')
    pdf(paste0('../figures/',expname,'/eu_1816_',tim,filenameext,'.pdf'), width=9, height=9, paper='special')
    layout(matrix(c(1,2,3,4,4,4,5,6,7,8,8,8), 4, 3, byrow = TRUE), height=c(3,1,3,1))
    #  layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
    par(oma=c(2,4,4,0))
    levs <- c(-Inf, seq(-1.75,1.75,0.5), Inf)
    contlevs <- seq(-2,2,0.5)
    pdata$data <- array(c(echam.anom$ensmean[,tim],analysis.anom$ensmean[,tim],
                      validate.anom$data[,tim]), 
                      c(nrow(validate.anom$data),1,3))
    plot_echam3(pdata, varname='temp2', type='data', cex.pt=5.2, lonlim=c(-20,40),
             latlim=c(30,75), wcol='darkgrey',
             lev=levs, st.col=NULL, stations=calibrate, names=c('a)', 'b)','c)'),
             colnames=c("CCC400 anomaly","EKF400 anomaly","Reconstruction anomaly"),
             statpanel=2, add=T, rownames='Temperature and SLP', main='1816', #no ghp in vali data
             addcontours=T, contvarname='slp',conttype='data',contcol='black',
             contlev=contlevs)
    scfac <- max(echam.anom$ensmean[echam$names=='u200'])
    levs <- c(-Inf, seq(-9,9,2), Inf)
    plot_echam3(pdata, varname='precip', type='data', cex.pt=5.2, lonlim=c(-20,40),
             latlim=c(30,75), wcol='darkgrey',
             lev=levs, st.col=NULL, stations=NULL, names=c('d)', 'e)','f)'),
             colnames='',
             statpanel=NULL, add=T, rownames='Precipitation and 850hPa wind', main='',
             addvectors=T, vecnames=c('u850','v850'), veccol='black', 
             veclen=scfac*0.01, vecscale=scfac*0.9, vecwd=0.95, every_x_vec=1,
             colorbar=T)
    dev.off()
#  }
  
  
  #  plotdata$data <- array(c(echmeananom[,ti],anameananom[,ti],eadiff[,ti],valmeananom[,ti]), c(nrow(echmeananom),1,4))
  #plotdata$data <- array(c(anameananom[,ti],valmeananom[,ti]), c(nrow(echmeananom),1,2))
#  plotdata$data <- array(c(anameananom[,ti],anameananombem[,ti],valmeananom[,ti]), c(nrow(echmeananom),1,3))  
# 
#   pdf('../figures/1816/1816.pdf', width=13.5, height=9, paper='special')
# #  layout(matrix(c(1,2,3,3,4,5,6,6,7,8,9,9), 6, 2, byrow = TRUE), height=c(3,1,3,1,3,1))
#   layout(matrix(c(1,2,3,4,4,4,5,6,7,8,8,8,9,10,11,12,12,12), 6, 3, byrow = TRUE), height=c(3,1,3,1,3,1))
#   par(oma=c(0,0,0,0))
#   levs <- c(-Inf, seq(-2,2,0.5), Inf)
#   plot_echam(plotdata, varname='temp2', type='data', cex.pt=1.5, names=pnames[1:2],
#              lev=levs, st.col=NULL, stations=NULL, add=T)
#   levs <- c(-Inf, seq(-8,8,2), Inf)
#   plot_echam(plotdata, varname='precip', type='data', cex.pt=1.5, names=pnames[3:4],
#              lev=levs, st.col=NULL, stations=NULL, add=T)
#   levs <- c(-Inf, seq(-4,4,1), Inf)
#   plot_echam(plotdata, varname='slp', type='data', cex.pt=1.5, names=pnames[5:6],
#              lev=levs, st.col=NULL, stations=NULL, add=T)
#   dev.off()
  
  
  
  
  
  
# Figure 3: 1816 alternative  
  ti=2 # for summer 1816
  pdata<-echam
#   pdata$data <- array(c(echmeananom[,ti],anameananom[,ti],anameananombem[,ti]), c(nrow(echmeananom),1,3))  
#   pdf('../figures/figure_3_1816_alltemp.pdf', width=10, height=6, paper='special')
#   layout(matrix(c(1,2,3,4,4,4,5,6,7,8,8,8), 4, 3, byrow = TRUE), height=c(3,0.5,3,1))
#   par(oma=c(2,2,4,2))
#   levs <- c(-Inf, seq(-2,2,0.5), Inf)
#   plot_echam(pdata, varname='temp2', type='data', cex.pt=2.0,names=c('a)', 'b)', 'c)'),
#               lev=levs, st.col=NULL, stations=calibrate, statpanel=2,
#               add=T ,units='ºC',seas=1,colnames=c("CCC","EKF","BEM"),
#               rownames=c("this study","reconstructions"),
#               latlim=c(0,90),main='1816',colorbar=F)
  pdata$data <- array(c(echmeananom[,ti],anameananom[,ti]), c(nrow(echmeananom),1,2))  
  pdf(paste0('../figures/',expname,'/1816_alltemp_v3',filenameext,'.pdf'), width=10, height=6, paper='special')
  layout(matrix(c(1,2,3,4,4,4,5,6,7,8,8,8), 4, 3, byrow = TRUE), height=c(3,0.5,3,1))
  par(oma=c(2,2,4,2))
  levs <- c(-Inf, seq(-2,2,0.5), Inf)
  plot_echam(pdata, varname='temp2', type='data', cex.pt=2.0,names=c('a)', 'b)'),
                lev=levs, st.col=NULL, stations=calibrate, statpanel=2,
                add=T ,units='ºC',seas=1,colnames=c("CCC","EKF"),
                rownames="this study",latlim=c(0,90),main='1816',colorbar=F)  
  pdata$data <- array(anameananombem[,ti], c(nrow(echmeananom),1,1))  
  coor <- read.csv('/Users/joerg/Documents/unibe/projects/2013_small/pages_paper/ccc400_coor_all.csv')
  stat35 <- calibrate
  stat35$lon <- coor$"echam_lon"
  stat35$lon[stat35$lon>180] <- stat35$lon[stat35$lon>180]-360
  stat35$lat <- coor$"echam_lat"
  stat35$names <- rep('prox',length(stat35$lon))
  plot_echam(pdata, varname='temp2', type='data', cex.pt=2.0,names='c)',
               lev=levs, st.col=NULL, stations=stat35, statpanel=1,
               add=T ,units='ºC',seas=1,colnames="BEM",
               latlim=c(0,90),main='',colorbar=F) #,centercol='lightgrey')    
    
  plot(NA,axes=F,bty='n',ylim=c(1,2))
  
# Luterbacher recon  
  pdata$data <- array(valmeananom[,ti], c(nrow(valmeananom),1,1)) 
  pdata$lon <- validate$lon  
  pdata$lat <- validate$lat
  pdata$names <- validate$names
#  levs <- c(-Inf, seq(-2,2,0.5), Inf)
  plot_echam(pdata, varname='temp2', type='data', cex.pt=1.5,latlim=c(0,90),
              seas=1,colorbar=F,rownames='reconstructions',colnames='Luterbacher',
              lev=levs, st.col=NULL, stations=NULL, add=T ,units='ºC',
              names='d)') #,centercol='lightgrey')
  
# read mann et al temp recon and plot 1816
  mann <- read.table("/Users/joerg/Documents/climdata/mann_temp_recon/allproxyfieldrecon",
                     na.string="NaN")
  ll <- read.table("/Users/joerg/Documents/climdata/mann_temp_recon/longlat.txt",
                   na.string="NaN")
  ts <- which(mann[,1]==1816)
  plotmann <- pdata
  mdata <- as.vector(t(mann[ts,-1]))
  mdata2 <- apply(mann[(ts-35):(ts+35),-1],2,mean)
  plotmann$data <- array(mdata-mdata2,c(length(mdata),1))
  plotmann$time <- mann[ts,1]
  plotmann$names <- rep('temp2',dim(plotmann$data)[1])
  plotmann$lon <- ll[,1] 
  plotmann$lat <- ll[,2]
  plotmann$height <- ll[,3]   
  plotmann$lsm.i <- ll[,3]  
#  levs <- c(-Inf, seq(-2,2,0.5), Inf)
  plot_echam(plotmann, varname='temp2', type='data', cex.pt=2.,names='e)',
              latlim=c(0,90),seas=1,colorbar=F,wcol="black",colnames='Mann',
              lev=levs, st.col=NULL, stations=NULL, add=T,units='ºC') #,centercol='lightgrey')
  
  schw <- read.table("/Users/joerg/Documents/climdata/recon/schweingr_mxd_grid/schweingruber_mxdabd_grid.dat",
                     na.string="-9.990")
  lon <- read.table("~/Documents/climdata/recon/schweingr_mxd_grid/schweingruber_lon.dat")
  lat <- read.table("~/Documents/climdata/recon/schweingr_mxd_grid/schweingruber_lat.dat")
  yrs <- seq(1400,1994)
  ts <- which(yrs==1816)
  plotschw <- echam
  sdata <- as.vector(t(schw[ts,]))
  sdata2 <- apply(schw[(ts-35):(ts+35),],2,mean)
  plotschw$data <- array(scale(sdata-sdata2,center=F),c(length(sdata),1))
  plotschw$ensmean <- array(sdata,c(length(sdata),1))
  plotschw$time <- schw[ts,1]
  plotschw$names <- rep('temp2',dim(plotschw$data)[1])
  plotschw$lon <- as.vector(t(lon))
  plotschw$lat <- as.vector(t(lat))
  plotschw$height <- rep(NA,length(lon))
  plotschw$lsm.i <- rep(NA,length(lon))
#  levs <- c(-Inf, seq(-3,3,0.5), Inf)
  plot_echam(plotschw, varname='temp2', type='data', cex.pt=2.,names='f)',
             latlim=c(0,90),seas=1,colorbar=T,wcol="black",colnames='Briffa',
             lev=levs, st.col=NULL, stations=NULL, add=T,units='') #,centercol='lightgrey')
  
dev.off()
  
 
  
  
  
  


  
  
plotdata$data <- array(echmeananom[,ti], c(nrow(echmeananom),1))  
pdf(paste0('../figures/',expname,'/1816_CCC400',filenameext,'.pdf'), width=4, height=10, paper='special')
layout(matrix(c(1,2,3,4,5,6), 6, 1, byrow = TRUE), height=c(3,1,3,1,3,1))
par(oma=c(0,0,4,0))
levs <- c(-Inf, seq(-2,2,0.5), Inf)
plot_echam2(plotdata, varname='temp2', type='data', cex.pt=1.5, names='',
            lev=levs, st.col=NULL, stations=NULL, add=T ,units='ºC',main='CCC400')
levs <- c(-Inf, seq(-4,4,1), Inf)
plot_echam2(plotdata, varname='slp', type='data', cex.pt=1.5, names='',
            lev=levs, st.col=NULL, stations=NULL, add=T ,units='hPa')
levs <- c(-Inf, seq(-15,15,3), Inf)
plot_echam2(plotdata, varname='precip', type='data', cex.pt=1.5, names='',
            lev=levs, st.col=NULL, stations=NULL, add=T, ,units='mm')
dev.off()

plotdata$data <- array(anameananom[,ti], c(nrow(echmeananom),1))  
pdf(paste0('../figures/',expname,'/1816_EKF',filenameext,'.pdf'), width=4, height=10, paper='special')
layout(matrix(c(1,2,3,4,5,6), 6, 1, byrow = TRUE), height=c(3,1,3,1,3,1))
par(oma=c(0,0,4,0))
levs <- c(-Inf, seq(-2,2,0.5), Inf)
plot_echam2(plotdata, varname='temp2', type='data', cex.pt=1.5, names='',
           lev=levs, st.col=NULL, stations=NULL, add=T ,units='ºC',main='EKF')
levs <- c(-Inf, seq(-4,4,1), Inf)
plot_echam2(plotdata, varname='slp', type='data', cex.pt=1.5, names='',
           lev=levs, st.col=NULL, stations=NULL, add=T ,units='hPa')
levs <- c(-Inf, seq(-15,15,3), Inf)
plot_echam2(plotdata, varname='precip', type='data', cex.pt=1.5, names='',
           lev=levs, st.col=NULL, stations=NULL, add=T, ,units='mm')
dev.off()

plotdata$data <- array(eadiff[,ti], c(nrow(echmeananom),1))  
pdf(paste0('../figures/',expname,'/1816_EKF-CCC400_all',filenameext,'.pdf'), width=8, height=6, paper='special')
layout(matrix(c(1,2,3,4,5,6,7,8), 4, 2, byrow = F), height=c(3,1,3,1))
par(oma=c(0,0,4,0))
levs <- c(-Inf, seq(-20,20,5), Inf)
plot_echam2(plotdata, varname='gph500', type='data', cex.pt=1.5, names='',
            lev=levs, st.col="cyan", stations=station, add=T ,units='hPa')
levs <- c(-Inf, seq(-1,1,0.25), Inf)
plot_echam2(plotdata, varname='u850', type='data', cex.pt=1.5, names='',
            lev=levs, st.col="cyan", stations=station, add=T ,units='m/s')
plot_echam2(plotdata, varname='v850', type='data', cex.pt=1.5, names='',
            lev=levs, st.col="cyan", stations=station, add=T ,units='m/s')
levs <- c(-Inf, seq(-0.005,0.005,0.001), Inf)
plot_echam2(plotdata, varname='omega500', type='data', cex.pt=1.5, names='',
            lev=levs, st.col="cyan", stations=NULL, add=T ,units='m/s')
dev.off()

plotdata$data <- array(anameananombem[,ti], c(nrow(echmeananom),1))  
pdf(paste0('../figures/',expname,'/1816_BEM',filenameext,'.pdf'), width=4, height=10, paper='special')
layout(matrix(c(1,2,3,4,5,6), 6, 1, byrow = TRUE), height=c(3,1,3,1,3,1))
par(oma=c(0,0,4,0))
levs <- c(-Inf, seq(-2,2,0.5), Inf)
plot_echam2(plotdata, varname='temp2', type='data', cex.pt=1.5, names='',
            lev=levs, st.col=NULL, stations=NULL, add=T ,units='ºC',main='BEM')
levs <- c(-Inf, seq(-4,4,1), Inf)
plot_echam2(plotdata, varname='slp', type='data', cex.pt=1.5, names='',
            lev=levs, st.col=NULL, stations=NULL, add=T ,units='hPa')
levs <- c(-Inf, seq(-15,15,3), Inf)
plot_echam2(plotdata, varname='precip', type='data', cex.pt=1.5, names='',
            lev=levs, st.col=NULL, stations=NULL, add=T, ,units='mm')
dev.off()

plotdata$data <- array(valmeananom[,ti], c(nrow(valmeananom),1))  
pdf(paste0('../figures/',expname,'/1816_LUT',filenameext,'.pdf'), width=4, height=10, paper='special')
layout(matrix(c(1,2,3,4,5,6), 6, 1, byrow = TRUE), height=c(3,1,3,1,3,1))
par(oma=c(0,0,4,0))
levs <- c(-Inf, seq(-2,2,0.5), Inf)
plot_echam2(plotdata, varname='temp2', type='data', cex.pt=1.5, names='',
            lev=levs, st.col=NULL, stations=NULL, add=T ,units='ºC',main='Reconstructions')
levs <- c(-Inf, seq(-4,4,1), Inf)
plot_echam2(plotdata, varname='slp', type='data', cex.pt=1.5, names='',
            lev=levs, st.col=NULL, stations=NULL, add=T ,units='hPa')
levs <- c(-Inf, seq(-15,15,3), Inf)
plot_echam2(plotdata, varname='precip', type='data', cex.pt=1.5, names='',
            lev=levs, st.col=NULL, stations=NULL, add=T, ,units='mm')
dev.off()


  # read mann et al temp recon and plot 1816
  mann <- read.table("/Users/joerg/Documents/climdata/mann_temp_recon/allproxyfieldrecon",
                   na.string="NaN")
  ll <- read.table("/Users/joerg/Documents/climdata/mann_temp_recon/longlat.txt",
                 na.string="NaN")
  ts <- which(mann[,1]==1816)
  plotmann <- plotdata
  mdata <- as.vector(t(mann[ts,-1]))
  #plotmann$data <- array(rep(mdata,2),c(length(mdata),1,2))
  plotmann$data <- array(mdata,c(length(mdata),1))
  #plotmann$data <- array(mdata,,c(length(mdata),1,1))
  #plotmann$ensmean <- array(mdata,c(length(mdata),1,1))
  #plotmann$time <- rep(mann[ts,1],2)
  plotmann$time <- mann[ts,1]
  plotmann$names <- rep('temp2',dim(plotmann$data)[1])
  plotmann$lon <- ll[,1] 
  plotmann$lat <- ll[,2]
  plotmann$height <- ll[,3]   
  plotmann$lsm.i <- ll[,3]  
  pdf(paste0('../figures/',expname,'/1816mann',filenameext,'.pdf'), width=8, height=6, paper='special')
    layout(matrix(c(1,2,3,3), 2, 1, byrow = TRUE), height=c(3,1))
  #layout(matrix(c(1,2,2), 2, 2, byrow = TRUE), height=c(3,1))
    par(oma=c(0,0,0,0))
    levs <- c(-Inf, seq(-2,2,0.5), Inf)
    plot_echam2(plotmann, varname='temp2', type='data', cex.pt=1.5, names='',
             lev=levs, st.col=NULL, stations=NULL, add=T,units='ºC')
  dev.off()
  
  
# schweingruber mxd grid by rutherford and briffa
#  The (text or ASCII format file) has 115 columns and 595 rows. 
#  Each column represents one grid-box time series, 
#  while each row represents one year, starting in 1400 and ending in 1994
  schw <- read.table("/Users/joerg/Documents/climdata/recon/schweingr_mxd_grid/schweingruber_mxdabd_grid.dat",na.string="-9.990")
  lon <- read.table("~/Documents/climdata/recon/schweingr_mxd_grid/schweingruber_lon.dat")
  lat <- read.table("~/Documents/climdata/recon/schweingr_mxd_grid/schweingruber_lat.dat")
  yrs <- seq(1400,1994)
  ts <- which(yrs==1816)
  plotschw <- plotdata
  sdata <- as.vector(t(schw[ts,]))
# plotschw$data <- array(rep(sdata,2),c(length(sdata),1,2))
# plotschw$ensmean <- array(rep(sdata,2),c(length(sdata),1,2))
# plotschw$time <- rep(mann[ts,1],2)
  plotschw$data <- array(sdata,c(length(sdata),1))
  plotschw$ensmean <- array(sdata,c(length(sdata),1))
  plotschw$time <- schw[ts,1]
  plotschw$names <- rep('temp2',dim(plotschw$data)[1])
  plotschw$lon <- as.vector(t(lon))
  plotschw$lat <- as.vector(t(lat))
  plotschw$height <- rep(NA,length(lon))
  plotschw$lsm.i <- rep(NA,length(lon))
  pdf(paste0('../figures/',expname,'/1816schw',filenameext,'.pdf'), width=8, height=4, paper='special')
    layout(matrix(c(1,2), 2, 1, byrow = TRUE), height=c(2,1))
    par(oma=c(0,0,0,0))
    levs <- c(-Inf, seq(-3,3,0.5), Inf)
    plot_echam2(plotschw, varname='temp2', type='data', cex.pt=3.0, names='',
             lev=levs, st.col=NULL, stations=NULL, add=T, units='TRW')
  dev.off()  


# pdsi 1816
  mada <- read_pdsi('mada_pdsi_echamgrid.nc',timlim=c(syr,eyr), 
                  path="/Users/joerg/Documents/climdata/recon/cook/pdsi_asia/",
                  small=T,landonly=F)
  nada <- read_pdsi('nada_pdsi_echamgrid.nc',timlim=c(syr,eyr), 
                  path="/Users/joerg/Documents/climdata/recon/cook/pdsi_usa/",
                  small=T,landonly=F)
  nada$data[is.na(nada$data)] <- 0
  mada$data[is.na(mada$data)] <- 0
  pdsi <- nada
  tmp <- nada$data + mada$data
  tmp[tmp==0] <- NA
  pdsi$data <- apply(tmp[,,1],1,mean,na.rm=T)
  # erst 4608 temp und dann erst precip
  #pdsi$data <- rep(pdsi$data,2)
  pdsi$names <- rep('precip',length(pdsi$names))

  # make plots
  ti=1 # summer season as drought recon is JJA
  plotdata=echam
  #  plotdata$data <- array(c(echmeananom[,ti],anameananom[,ti],eadiff[,ti],valmeananom[,ti]), c(nrow(echmeananom),1,4))
#  plotdata$data <- array(c(anameananom[,ti],valmeananom[,ti]), c(nrow(echmeananom),1,2))

  pdf(paste0('../figures/',expname,'/1816_pdsi',filenameext,'.pdf'), width=8, height=6, paper='special')
    layout(matrix(c(1,2), 2, 1, byrow = TRUE), height=c(3,1))
    par(oma=c(0,0,0,0))
    plotdata$data <- array(pdsi$data, c(nrow(echmeananom),1))
    levs <- c(-Inf, seq(-2,2,0.5), Inf)
    plot_echam2(plotdata, varname='precip', type='data', cex.pt=1.5, names='',
           lev=levs, st.col=NULL, stations=NULL, add=T, units='PDSI')
  dev.off()



# 20CR 1816 plot, anomalies to 20cr 2c 1851-80, because no good data for earlier period
#cdo merge t_AMJJAS_1816_monmean_ano.nc p_AMJJAS_1816_monmean_ano.nc slp_AMJJAS_1816_monmean_ano.nc tpslp_AMJJAS_1816_monmean_ano.nc
scout <- read_20cr('tpslp_AMJJAS_1816_summermean_ano.nc', path=twentycrpath, 
            xlim=c(-180,180), ylim=c(-90,90), timlim=c(syr, eyr), small=F, landonly=F)

# make plots
ti=1 # summer season as drought recon is JJA
plotdata=scout
#  plotdata$data <- array(c(echmeananom[,ti],anameananom[,ti],eadiff[,ti],valmeananom[,ti]), c(nrow(echmeananom),1,4))
#plotdata$data <- array(c(scout$data,scout$data), c(length(scout$data),1,2))
#plotdata$data <- array(scout$data, c(length(scout$data),1,1))
plotdata$data <- array(scout$data, c(length(scout$data),1))

pdf(paste0('../figures/',expname,'/20cr_1816',filenameext,'.pdf'), width=4, height=10, paper='special')
#  layout(matrix(c(1,2,3,3,4,5,6,6,7,8,9,9), 6, 2, byrow = TRUE), height=c(3,1,3,1,3,1))
  layout(matrix(c(1,2,3,4,5,6), 6, 1, byrow = TRUE), height=c(3,1,3,1,3,1))
  par(oma=c(0,0,4,0))
  levs <- c(-Inf, seq(-2,2,0.5), Inf)
  plot_echam2(plotdata, varname='temp2', type='data', cex.pt=1.5, seas=1,
        names='', lev=levs, st.col=NULL, stations=NULL, add=T, 
        main='20CR', units='ºC')
  levs <- c(-Inf, seq(-4,4,1), Inf)
  plot_echam2(plotdata, varname='slp', type='data', cex.pt=1.5, 
        names='', lev=levs, st.col=NULL, stations=NULL, add=T,
        units='hPa')
levs <- c(-Inf, seq(-15,15,3), Inf)
plot_echam2(plotdata, varname='precip', type='data', cex.pt=1.5, 
            names='', lev=levs, st.col=NULL, stations=NULL, add=T,
            units='mm')
dev.off()


# 
# # all 1816 plots together
# plotdata$data <- array(c(anameananom[,2],anameananombem[,2],valmeananom[,2]),scout$data, c(nrow(echmeananom),1,4))
# plotdata$lon <- c(rep(echam$lon,3),scout$lon)
# plotdata$lat <- c(rep(echam$lat,3),scout$lat)
# plotdata$names <- c(rep(echam$names,3),scout$names)
# 
# pdf('../figures/1816/1816.pdf', width=13.5, height=9, paper='special')
# #  layout(matrix(c(1,2,3,3,4,5,6,6,7,8,9,9), 6, 2, byrow = TRUE), height=c(3,1,3,1,3,1))
# layout(matrix(c(1,2,3,4,4,4,5,6,7,8,8,8,9,10,11,12,12,12), 6, 3, byrow = TRUE), height=c(3,1,3,1,3,1))
# par(oma=c(0,0,0,0))
# levs <- c(-Inf, seq(-2,2,0.5), Inf)
# plot_echam(plotdata, varname='temp2', type='data', cex.pt=1.5, names=pnames[1:2],
#            lev=levs, st.col=NULL, stations=NULL, add=T)
# levs <- c(-Inf, seq(-8,8,2), Inf)
# plot_echam(plotdata, varname='precip', type='data', cex.pt=1.5, names=pnames[3:4],
#            lev=levs, st.col=NULL, stations=NULL, add=T)
# levs <- c(-Inf, seq(-4,4,1), Inf)
# plot_echam(plotdata, varname='slp', type='data', cex.pt=1.5, names=pnames[5:6],
#            lev=levs, st.col=NULL, stations=NULL, add=T)
# dev.off()

  
  
} # end 1816 plots 









# find dryest, warmest, ... years in NH and specific regions
if (multiann) {
  # SW USA droughts
  # no index defined yet for that region. First check somewhere else if it works!
  # Not ideal that simulation finish 2005 and we do not have this drought in your analysis
  # northern hemisphere extra-tropics (>20N)
  par(mfrow=c(1,1))
  par(mai=c(1,1,0.5,0.5))
  data <- plot_pages(wl=1,reg=which(aind.allts$names=='MED.precip'),seas='sum',
                     anomper=anomwrt,lw=2,sca=F,xa='s',plotbem=F,plotmem=F,
                     plotech=T,pag=F,use18=T,yl="Precip.")
  data <- plot_pages(wl=1,reg=which(aind.allts$names=='MED.slp'),seas='sum',
                     anomper=anomwrt,lw=2,sca=F,xa='s',plotbem=F,plotmem=T,
                     plotech=T,pag=F,use18=T,yl="SLP")
  data <- plot_pages(wl=1,reg=which(aind.allts$names=='ENH.temp2'),seas='sum',
                     anomper=anomwrt,lw=2,sca=F,xa='s',plotbem=F,plotmem=T,
                     plotech=T,pag=F,use18=T,yl="Temp.")
  plot(runmean(data$allindann[,"aindmean"],31),ty='l')
  datdetr<-data$allindann[,"aindmean"]-runmean(data$allindann[,"aindmean"],31)
  plot(datdetr,ty='l')
  data$period[which(datdetr < quantile(datdetr, probs = 0.1))]
} # multiann, see more code below
  
  
    






if (glacier_advances) {  
  # 1830-1850 glacier advances
  # NAM, NEU and MED winter precip and summer temp
  wlen=1
  pdf(file=paste0('../figures/',expname,'/1830-50_winprecip_sumtemp_',filenameext,
                  wlen,'yr_runmean.pdf'),width=10,height=6)
    par(mfrow=c(2,3))
    par(mai=c(0.3,0.6,0.2,0))
    pdata <- aind.allts$ensmean[which(aind.allts$names=='NEU.precip'),
                          !is.even(seq(1,ncol(aind.allts$ensmean)))]
    plot(seq(1750,1930),runmean(pdata[which((aind.allts$time>1749)&(aind.allts$time<1931))],wlen),
       ty='l',col='red',main='NEU Winter precip.',ylab='')
    pdata <- eind.allts$ensmean[which(aind.allts$names=='NEU.precip'),
                              !is.even(seq(1,ncol(aind.allts$ensmean)))]
    lines(seq(1750,1930),runmean(pdata[which((aind.allts$time>1749)&(aind.allts$time<1931))],wlen),
       ty='l',col='blue',main='',ylab='')
    abline(v=c(1830,1850),lty=2)
    pdata <- aind.allts$ensmean[which(aind.allts$names=='MED.precip'),
                              !is.even(seq(1,ncol(aind.allts$ensmean)))]
    plot(seq(1750,1930),runmean(pdata[which((aind.allts$time>1749)&(aind.allts$time<1931))],wlen),
       ty='l',col='red',main='MED Winter precip.',ylab='')
    pdata <- eind.allts$ensmean[which(aind.allts$names=='MED.precip'),
                              !is.even(seq(1,ncol(aind.allts$ensmean)))]
    lines(seq(1750,1930),runmean(pdata[which((aind.allts$time>1749)&(aind.allts$time<1931))],wlen),
        ty='l',col='blue',main='',ylab='')
    abline(v=c(1830,1850),lty=2)
    pdata <- aind.allts$ensmean[which(aind.allts$names=='NAM.precip'),
                              !is.even(seq(1,ncol(aind.allts$ensmean)))]
    plot(seq(1750,1930),runmean(pdata[which((aind.allts$time>1749)&(aind.allts$time<1931))],wlen),
       ty='l',col='red',main='NAM Winter precip.',ylab='')
    pdata <- eind.allts$ensmean[which(aind.allts$names=='NAM.precip'),
                              !is.even(seq(1,ncol(aind.allts$ensmean)))]
    lines(seq(1750,1930),runmean(pdata[which((aind.allts$time>1749)&(aind.allts$time<1931))],wlen),
        ty='l',col='blue',main='',ylab='')
    abline(v=c(1830,1850),lty=2)
  
    pdata <- aind.allts$ensmean[which(aind.allts$names=='NEU.temp2'),
                              is.even(seq(1,ncol(aind.allts$ensmean)))]
    plot(seq(1750,1930),runmean(pdata[which((aind.allts$time>1749)&(aind.allts$time<1931))],wlen),
       ty='l',col='red',main='NEU Summer Temp.',ylab='')
    pdata <- eind.allts$ensmean[which(aind.allts$names=='NEU.temp2'),
                              is.even(seq(1,ncol(aind.allts$ensmean)))]
    lines(seq(1750,1930),runmean(pdata[which((aind.allts$time>1749)&(aind.allts$time<1931))],wlen),
        ty='l',col='blue',main='',ylab='')
    abline(v=c(1830,1850),lty=2)
    pdata <- aind.allts$ensmean[which(aind.allts$names=='MED.temp2'),
                              is.even(seq(1,ncol(aind.allts$ensmean)))]
    plot(seq(1750,1930),runmean(pdata[which((aind.allts$time>1749)&(aind.allts$time<1931))],wlen),
       ty='l',col='red',main='MED Summer Temp.',ylab='')
    pdata <- eind.allts$ensmean[which(aind.allts$names=='MED.temp2'),
                              is.even(seq(1,ncol(aind.allts$ensmean)))]
    lines(seq(1750,1930),runmean(pdata[which((aind.allts$time>1749)&(aind.allts$time<1931))],wlen),
        ty='l',col='blue',main='',ylab='')
    abline(v=c(1830,1850),lty=2)
    pdata <- aind.allts$ensmean[which(aind.allts$names=='NAM.temp2'),
                              is.even(seq(1,ncol(aind.allts$ensmean)))]
    plot(seq(1750,1930),runmean(pdata[which((aind.allts$time>1749)&(aind.allts$time<1931))],wlen),
       ty='l',col='red',main='NAM Summer Temp.',ylab='')
    pdata <- eind.allts$ensmean[which(aind.allts$names=='NAM.temp2'),
                              is.even(seq(1,ncol(aind.allts$ensmean)))]
    lines(seq(1750,1930),runmean(pdata[which((aind.allts$time>1749)&(aind.allts$time<1931))],wlen),
        ty='l',col='blue',main='',ylab='')
    abline(v=c(1830,1850),lty=2) 
    legend('topleft', c('EKF400','CCC400'),lwd=c(2,2), col=c('red','blue'), bty='o', 
         bg='white', box.col='white')  
  dev.off()
  

  
  
  ptm1 <- proc.time()
  # calc validate clim (70 yr period) to substract for anomalies
  for (cyr in seq(1800,1880)) {
    load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
    if (cyr==1800){
      vtmp <- validate$data
      c=1
    } else {
      vtmp <- vtmp+validate$data
      c <- c+1
    }
  }
  valclim <- vtmp/c
  
  # 1830-50 glacier advances
  #  for (cyr in c(1790, 1792, 1793, 1794, 1795, 1796)) {
  # define period
  syr <- 1830
  eyr <- 1850
  
  for (cyr in seq(syr,eyr)) {
    load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
    echam.abs <- echam
    analysis.abs <- analysis
    bemlist <- read.table(file='../data/bem/bem.txt',header=F)
    bem <- bemlist[which(bemlist[,1]==cyr),2]
    echam.abs$bem <- echam.abs$data[,,bem]
    echam.anom$bem <- echam.anom$data[,,bem]
    analysis.abs$bem <- analysis.abs$data[,,bem]
    analysis.anom$bem <- analysis.anom$data[,,bem]
#    ech_ind$bem <- ech_ind$data[,,bem]
#    ana_ind$bem <- ana_ind$data[,,bem]
    
    if (cyr==syr){
      atmp <- analysis.anom$ensmean
      etmp <- echam.anom$ensmean
      a2tmp <- analysis.abs$ensmean
      e2tmp <- echam.abs$ensmean
      a3tmp <- analysis.anom$bem
      e3tmp <- echam.anom$bem
      vtmp <- validate$data
      c=1
    } else {
      atmp <- atmp+analysis.anom$ensmean
      etmp <- etmp+echam.anom$ensmean
      a2tmp <- a2tmp+analysis.abs$ensmean
      e2tmp <- e2tmp+echam.abs$ensmean
      a3tmp <- a3tmp+analysis.anom$bem
      e3tmp <- e3tmp+echam.anom$bem
      vtmp <- vtmp+validate$data
      c=c+1
    }
  }
  anameananom <- atmp/c
  echmeananom <- etmp/c
  eadiff <- anameananom-echmeananom
  anameanabs <- a2tmp/c
  echmeanabs <- e2tmp/c
  eadiffabs <- anameanabs-echmeanabs
  anameananombem <- a3tmp/c
  echmeananombem <- e3tmp/c
  eadiffbem <- anameananombem-echmeananombem
  valmean <- vtmp/c
  valmeananom <- valmean-valclim
  
  pnames <- paste(c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k'), ')', sep='')
  data.dim <- dim(echam$data)[c(1,2,2,3)]
  data.dim[2:3] <- c(s,data.dim[3]/s)
  ens.dim <- c(nrow(echam$ensmean), s, ncol(echam$ensmean)/s)
  
  
  
  # open cook pdsi recons
  # first interpolate to echam global grid:
  # cdo -m 'NaNf' -setmissval,9e33 MADApdsi.nc madapdsi2.nc
  # cdo -r remapbil,t63grid madapdsi2.nc mada_pdsi_echamgrid.nc
  # cdo -m 'NaNf' -setmissval,9e33 NADAv2-2008.nc nadapdsi2.nc
  # cdo -r remapbil,t63grid nadapdsi2.nc nada_pdsi_echamgrid.nc
  mada <- read_pdsi('mada_pdsi_echamgrid.nc',timlim=c(syr,eyr), 
                    path="/Users/joerg/Documents/climdata/recon/cook/pdsi_asia/",
                    small=T,landonly=F)
  nada <- read_pdsi('nada_pdsi_echamgrid.nc',timlim=c(syr,eyr), 
                    path="/Users/joerg/Documents/climdata/recon/cook/pdsi_usa/",
                    small=T,landonly=F)
  nada$data[is.na(nada$data)] <- 0
  mada$data[is.na(mada$data)] <- 0
  pdsi <- nada
  tmp <- nada$data + mada$data
  tmp[tmp==0] <- NA
  pdsi$data <- apply(tmp[,,1],1,mean,na.rm=T)
  # erst 4608 temp und dann erst precip
  pdsi$data <- rep(pdsi$data,2)
  pdsi$names <- rep('precip',length(pdsi$names))
  
  # make plots
  ti=2 # summer season as drought recon is JJA
  plotdata=echam
  #  plotdata$data <- array(c(echmeananom[,ti],anameananom[,ti],eadiff[,ti],valmeananom[,ti]), c(nrow(echmeananom),1,4))
  plotdata$data <- array(c(anameananom[,ti],valmeananom[,ti]), c(nrow(echmeananom),1,2))
  
  pdf(paste0('../figures/',expname,'/drought_',syr,'-',eyr,filenameext,'.pdf'), 
      width=9, height=9, paper='special')
    layout(matrix(c(1,2,3,3,4,5,6,6), 4, 2, byrow = TRUE), height=c(3,1,3,1))
    par(oma=c(0,0,0,0))
    levs <- c(-Inf, seq(-8,8,2), Inf)
    plot_echam(plotdata, varname='precip', type='data', cex.pt=1.5, 
               names=pnames[1:dim(plotdata$data)[3]],
               lev=levs, st.col=NULL, stations=NULL, add=T)
  
    plotdata$data <- array(c(eadiff[,ti],pdsi$data), c(nrow(echmeananom),1,2))
    levs <- c(-Inf, seq(-2,2,0.5), Inf)
    plot_echam(plotdata, varname='precip', type='data', cex.pt=1.5, names=pnames[3:4],
             lev=levs, st.col=NULL, stations=NULL, add=T)
  dev.off()
  
  
  # search explanation for droughts: do HC index, omega500 and gph500 anomaly 
  # indicate subtropical subsidence?
  ti=2 # summer season as drought recon is JJA
  plotdata$data <- array(c(echmeananom[,ti],anameananom[,ti]), c(nrow(echam.anom$ensmean),1,2))
  pdf(paste0('../figures/',expname,'/gph500_anom_',syr,'-',eyr,filenameext,'.pdf'), 
      width=9, height=4.5, paper='special')
    layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
    par(oma=c(0,0,0,0))
    levs <- c(-Inf, seq(-8,8,2), Inf)
    plot_echam(plotdata, varname='gph500', type='data', cex.pt=1.5, 
             names=pnames[1:dim(plotdata$data)[3]], lev=levs, st.col=NULL, 
             stations=NULL, add=T)
  dev.off()
  
  plotdata$data <- array(c(echmeananombem[,ti],anameananombem[,ti]), c(nrow(echam.anom$ensmean),1,2))
  pdf(paste0('../figures/',expname,'/gph500_anom_bem_',syr,'-',eyr,filenameext,'.pdf'), 
      width=9, height=4.5, paper='special')
    layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
    par(oma=c(0,0,0,0))
    levs <- c(-Inf, seq(-12,12,2), Inf)
    plot_echam(plotdata, varname='gph500', type='data', cex.pt=1.5, 
             names=pnames[1:dim(plotdata$data)[3]], lev=levs, st.col=NULL, 
             stations=NULL, add=T)
  dev.off()
  
  plotdata$data <- array(c(echmeananom[,ti],anameananom[,ti]), c(nrow(echam.anom$ensmean),1,2))
  pdf(paste0('../figures/',expname,'/temp_anom_',syr,'-',eyr,filenameext,'.pdf'), 
      width=9, height=4.5, paper='special')
    layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
    par(oma=c(0,0,0,0))
    levs <- c(-Inf, seq(-1,1,0.1), Inf)
    plot_echam(plotdata, varname='temp2', type='data', cex.pt=1.5, 
             names=pnames[1:dim(plotdata$data)[3]], lev=levs, st.col=NULL, 
             stations=NULL, add=T)
  dev.off()
  
  plotdata$data <- array(c(echmeananom[,ti],anameananom[,ti]), c(nrow(echam.anom$ensmean),1,2))
  pdf(paste0('../figures/',expname,'/slp_anom_',syr,'-',eyr,filenameext,'.pdf'),
      width=9, height=4.5, paper='special')
    layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
    par(oma=c(0,0,0,0))
    levs <- c(-Inf, seq(-1,1,0.2), Inf)
    plot_echam(plotdata, varname='slp', type='data', cex.pt=1.5, 
             names=pnames[1:dim(plotdata$data)[3]], lev=levs, st.col=NULL, 
             stations=NULL, add=T)
  dev.off()
  
  plotdata$data <- array(c(eadiff[,ti],anameananom[,ti]), c(nrow(echam.anom$ensmean),1,2))
  pdf(paste0('../figures/',expname,'/omega500_anom_',syr,'-',eyr,filenameext,'.pdf'), 
      width=9, height=4.5, paper='special')
    layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
    par(oma=c(0,0,0,0))
    levs <- c(-Inf, seq(-0.001,0.001,0.0001), Inf)
    plot_echam(plotdata, varname='omega500', type='data', cex.pt=1.5, 
             names=pnames[1:dim(plotdata$data)[3]], lev=levs, st.col=NULL, 
             stations=NULL, add=T)
  dev.off()
  
  
  
  
  
  # 1830-1840 eruptions/cold decade (or 1816 for anniversary version?)
  pdata <- echam
  pdf(paste0('../figures/',expname,'/cold_decade_',syr,'-',eyr,filenameext,'.pdf'), width=6, height=6, paper='special')
    layout(matrix(c(1,2,3,3,4,5,6,6), 4, 2, byrow = TRUE), height=c(3,1,3,1))
    par(oma=c(0,0,0,0))
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    pdata$data <- array(anameananom[,1:2], c(nrow(echam.anom$ensmean),1,2))
    plot_echam(pdata, varname='temp2', type='data', cex.pt=1.5, latlim=c(-10,90),wcol='darkgrey',
             names=pnames[1:dim(pdata$data)[3]], lev=levs, st.col=NULL, stations=calibrate, 
             statpanel=c(1,2), add=T,main="Temp. & GPH500 anom.",
             addcontours=T, contvarname='gph500',conttype='data',contcol='black',contlev=contlevs)
    levs <- c(-Inf, seq(-2,2,0.5), Inf)
    pdata$data <- array(anameananom[,1:2], c(nrow(echam.anom$ensmean),1,2))
    plot_echam(pdata, varname='precip', type='data', cex.pt=1.5, latlim=c(-10,90),
             names=pnames[3:4], lev=levs, st.col=NULL, stations=NULL, add=T,
             addcontours=F,wcol='darkgrey')
  
  dev.off()
  

} # end glacier_advances









if (warmcold_us) {
  load(paste0('../data/indices/',expname,'/indices',filenameext,'allts.Rdata'))
   # find and mark coldest/warmest 8yr period in US
    d <- aind.allts$ensmean[(which(eind.allts$names=="NAM.temp2")),
                            (is.even(seq(1,ncol(aind.allts$ensmean))))]
    dec <- rep(NA,250)
    dectime <- rep(NA,250)
    for (ti in 1:length(dec)) {
      if (ti==1) {
        dec[ti] <- mean(d[1:8])
        dectime[ti] <- 1603
        i=2
      } else {
        dec[ti] <- mean(d[i:(i+7)])
        dectime[ti] <- dectime[ti-1] + 1
        i=i+1
      }
    }
    #  plot(dectime,dec,ty='l',col='red')
    tmp=cbind(dec,dectime)
    tmp[order(dec),]
  #coldper <- c(1641, 1697, 1835, 1667, 1738, 1836)
  coldper <- c(1641, 1835, 1809, 1697,1738)
  
  for (syr in coldper) {
    eyr <- syr+7
    # load nada pdsi recon
    nada <- read_pdsi('nada_pdsi_echamgrid.nc',timlim=c(syr,eyr), 
                      path="/Users/joerg/Documents/climdata/recon/cook/pdsi_usa/",
                      small=T,landonly=F)
    nada$data[is.na(nada$data)] <- 0
    # change to precip in % 
    # load monthly 70yr climatology data
    if (syr < 1636) {
      load(file=paste0(echclimpath,'echam_clim_1636-1637_2ndgrid.Rdata'))
    } else {
      load(file=paste0(echclimpath,'echam_clim_',syr+3,'-',(syr+4),'_2ndgrid.Rdata'))  
    }
    pos <- which(echam$names=="precip")[1:4608]
    # seasonal mean
    echmeanclim <- cbind(rowMeans(echam_clim$ensmean[,10:15]),
                         rowMeans(echam_clim$ensmean[,16:21]))
    #      echam$ensmean <- echmeanclim
    for (cyr in seq(syr,eyr)) {
      load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
      # change to precip in % 
      echam.abs <- echam
      echam.anom$ensmean[pos,] <- echam.abs$ensmean[pos,]/(echmeanclim[pos,]*3600*24*30)*100
      analysis.abs <- analysis
      analysis.anom$ensmean[pos,] <- analysis.abs$ensmean[pos,]/
                                        (echmeanclim[pos,]*3600*24*30)*100
      if (cyr==syr){
        atmp <- analysis.anom$ensmean
        etmp <- echam.anom$ensmean
        c=1
      } else {
        atmp <- atmp+analysis.anom$ensmean
        etmp <- etmp+echam.anom$ensmean
        c=c+1
      }
    }
    if (syr==coldper[1]){
      anameananom <- atmp/c
      echmeananom <- etmp/c
      pdsimeananom <- apply(nada$data,1,mean)
    } else {
      anameananom <- cbind(anameananom,(atmp/c))
      echmeananom <- cbind(echmeananom,(etmp/c))
      pdsimeananom <- cbind(pdsimeananom,(apply(nada$data,1,mean)))
    }
  } # end cold periods loop  
  colnames(anameananom) <- rep(coldper,each=2)
  colnames(echmeananom) <- rep(coldper,each=2)
    #
    #  plotdata <- echam
    #  plotdata$ensmean <- echmeananom # 1641
    #  scfac <- max(echam.anom$ensmean[echam$names=='u200'])
    #  plot_echam3(plotdata, varname='precip', type='ensmean', cex.pt=6, ti=2,
    #              levs=c(50,90,95,100,105,110,200),
    #              latlim=c(20,60),lonlim=c(-140,-40),
    #              main='Precip. & 200hPa wind anomalies', units='%',
    #              addvectors=T, vecnames=c('u200','v200'), veccol='black', 
    #              veclen=scfac*0.01, vecscale=scfac*0.9, vecwd=0.95, every_x_vec=1,
    #              wcol='darkgrey')
    #  
    #   # change to precip in % 
    #   load(file=paste0('../data/prepplot_season/analysis_',cyr,'.Rdata'))
    #   echam.abs <- echam
    #   load(file=paste0('../data/echam_clim/echam_clim_',cyr,'-',(cyr+1),'_2ndgrid.Rdata'))
    #   # seasonal mean
    #   echmeanclim <- cbind(rowMeans(echam_clim$ensmean[,10:15]),rowMeans(echam_clim$ensmean[,16:21]))
    #   echam$ensmean <- echmeanclim
    #   plot_echam(echam, varname='precip', type='ensmean', cex.pt=2, ti=1)
    #   pos <- which(echam$names=="precip")[1:4608]
    #   echam$ensmean[pos,] <- echam.abs$ensmean[pos,]/(echmeanclim[pos,]*3600*24*30)*100
    ##   quartz()
    #   plot_echam(echam,levs=c(50,85,90,95,100,105,110,115,200), varname='precip', 
    #              type='ensmean', cex.pt=2, ti=1)
    
    # 8-yr cold period winter
  pdata <- echam
  for (ti in 1:2){
    if (ti==1){s="win"} else {s="som"}
    #   ti=2 # summer
    pdf(paste0('../figures/',expname,'/cold_us_',s,filenameext,'.pdf'), width=12, 
          height=14, paper='special')
#      layout(matrix(c(1,2,3,4,5,6,7,8), 4, 2, byrow = F), height=c(3,1,3,1))
      layout(matrix(c(1,2,3,4,5,5,5,5,6,7,8,9,10,10,10,10,11,12,13,14,15,15,15,15,
             16,17,18,19,20,20,20,20),8,4,byrow=TRUE),height=c(3,1,3,1,3,1,3,1))
      par(oma=c(0.5,4,4,0))
      levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
      contlevs <- seq(-4,4,2)
      scfac <- max(echam.anom$ensmean[echam$names=='u200'])
#      pdata$data <- array(c(echmeananom[,ti],anameananom[,ti]), c(nrow(echam.anom$ensmean),1,2))
      pdata$data <- array(c(anameananom[,ti],anameananom[,ti+2],anameananom[,ti+4],
                            anameananom[,ti+6]), c(nrow(echam.anom$ensmean),1,4))
#      pdata$data <- array(c(echmeananom[,ti],echmeananom[,ti+2],echmeananom[,ti+4],
#                            echmeananom[,ti+6]), c(nrow(echam.anom$ensmean),1,4))      
#      pdata$dataechmeanclim[pos,]/echam$ensmean[pos,]*100
      plot_echam3(pdata, varname='temp2', type='data', cex.pt=3.7, 
                  latlim=c(20,60),lonlim=c(-160,-40),
                  wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                  names=c('a)', 'b)', 'c)', 'd)'),colnames=coldper[1:4],
                  statpanel=NULL, add=T, rownames='Analysis T2m, Z500, UV200', 
                  main='EKF400 5-yr cold periods in North America', units='K',
                  addcontours=T, contvarname='gph500', conttype='data',
                  contcol='black', contlev=contlevs,
                  addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                  veclen=scfac*0.02, vecscale=scfac*0.98, vecwd=0.95, every_x_vec=2)
      # diff echam analysis
      levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
      contlevs <- seq(-4,4,2)
      pdata$data <- array(c((anameananom[,ti]-echmeananom[,ti]),
                            (anameananom[,ti+2]-echmeananom[,ti+2]),
                            (anameananom[,ti+4]-echmeananom[,ti+4]),
                            (anameananom[,ti+6]-echmeananom[,ti+6])),
                            c(nrow(echam.anom$ensmean),1,4))
      plot_echam3(pdata, varname='temp2', type='data', cex.pt=3.7, 
                  latlim=c(20,60),lonlim=c(-160,-40),
                  wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                  names=c('e)', 'f)', 'g)', 'h)'),colnames=rep(" ",4),
                  statpanel=NULL, add=T, rownames='Update (Anal.-ECHAM)', 
                  units='K',
                  addcontours=T, contvarname='gph500', conttype='data',
                  contcol='black', contlev=contlevs,
                  addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                  veclen=scfac*0.02, vecscale=scfac*0.98, vecwd=0.95, every_x_vec=2)
      pdata$data <- array(c(anameananom[,ti],anameananom[,ti+2],anameananom[,ti+4],
                            anameananom[,ti+6]), c(nrow(echam.anom$ensmean),1,4))
#      #    levs <- c(-Inf, seq(-5,5,1), Inf)
      levs=c(50,80,85,90,95,100,105,110,115,120,200)
      plot_echam3(pdata, varname='precip', type='data', cex.pt=3.7, 
                  latlim=c(20,60),lonlim=c(-160,-40),units="%",
                  names=c('i)', 'j)','k)', 'l)'), lev=levs, 
                  st.col=NULL, stations=NULL, add=T,addcontours=F,
                  wcol='darkgrey',rownames='Precipitation',colnames=rep('',4))
      pdata$data <- array(c(pdsimeananom[,1],pdsimeananom[,2],pdsimeananom[,3],
                            pdsimeananom[,4]), c(nrow(echam.anom$ensmean),1,4))
      levs <- seq(-4,4,1)
      plot_echam3(pdata, varname='precip', type='data', cex.pt=3.7, 
                  latlim=c(20,60),lonlim=c(-160,-40),units="PDSI index",
                  names=c('m)', 'n)','o)', 'p)'), lev=levs, 
                  st.col=NULL, stations=NULL, add=T,addcontours=F,
                  wcol='darkgrey',rownames='PDSI',colnames=rep('',4))
    dev.off()  
  } # end summer winter loop
#  } # end cold periods loop 
  
  # warmper <- c(1825, 1790, 1842, 1800, 1772, 1733)
  warmper <- c(1659, 1825, 1733, 1850, 1618,1798)
  for (syr in warmper) {
    eyr <- syr+7
    nada <- read_pdsi('nada_pdsi_echamgrid.nc',timlim=c(syr,eyr), 
                      path="/Users/joerg/Documents/climdata/recon/cook/pdsi_usa/",
                      small=T,landonly=F)
    nada$data[is.na(nada$data)] <- 0
    # change to precip in % 
    # load monthly 70yr climatology data
    if (syr < 1636) {
      load(file=paste0(echclimpath,'echam_clim_1636-1637_2ndgrid.Rdata'))
    } else {
      load(file=paste0(echclimpath,'echam_clim_',syr+3,'-',(syr+4),'_2ndgrid.Rdata'))  
    }
    pos <- which(echam$names=="precip")[1:4608]
    # seasonal mean
    echmeanclim <- cbind(rowMeans(echam_clim$ensmean[,10:15]),
                         rowMeans(echam_clim$ensmean[,16:21]))
    for (cyr in seq(syr,eyr)) {
      load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
      # change to precip in % 
      echam.abs <- echam
      echam.anom$ensmean[pos,] <- echam.abs$ensmean[pos,]/(echmeanclim[pos,]*3600*24*30)*100
      analysis.abs <- analysis
      analysis.anom$ensmean[pos,] <- analysis.abs$ensmean[pos,]/
                                        (echmeanclim[pos,]*3600*24*30)*100
      if (cyr==syr){
        atmp <- analysis.anom$ensmean
        etmp <- echam.anom$ensmean
        c=1
      } else {
        atmp <- atmp+analysis.anom$ensmean
        etmp <- etmp+echam.anom$ensmean
        c=c+1
      }
    }
    if (syr==warmper[1]){
      anameananom <- atmp/c
      echmeananom <- etmp/c
      pdsimeananom <- apply(nada$data,1,mean)
    } else {
      anameananom <- cbind(anameananom,(atmp/c))
      echmeananom <- cbind(echmeananom,(etmp/c))
      pdsimeananom <- cbind(pdsimeananom,(apply(nada$data,1,mean)))
    }
  } # end warm periods loop  
  colnames(anameananom) <- rep(warmper,each=2)
  colnames(echmeananom) <- rep(warmper,each=2)
  
  # 8-yr warm period summer and winter
  pdata <- echam
  for (ti in 1:2){
    if (ti==1){s="win"} else {s="som"}
    pdf(paste0('../figures/',expname,'/warm_us_',s,filenameext,'.pdf'), width=12, 
        height=10.5, paper='special')
    #      layout(matrix(c(1,2,3,4,5,6,7,8), 4, 2, byrow = F), height=c(3,1,3,1))
    layout(matrix(c(1,2,3,4,5,5,5,5,6,7,8,9,10,10,10,10,11,12,13,14,15,15,15,15,
                    16,17,18,19,20,20,20,20),8,4,byrow=TRUE),height=c(3,1,3,1,3,1,3,1))
    par(oma=c(0.5,4,4,0))
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    scfac <- max(echam.anom$ensmean[echam$names=='u200'])
    pdata$data <- array(c(anameananom[,ti],anameananom[,ti+2],anameananom[,ti+4],
                          anameananom[,ti+6]), c(nrow(echam.anom$ensmean),1,4))
    plot_echam3(pdata, varname='temp2', type='data', cex.pt=3.7, 
                latlim=c(20,60),lonlim=c(-160,-40),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('a)', 'b)', 'c)', 'd)'),colnames=warmper[1:4],
                statpanel=NULL, add=T, rownames='Analysis T2m, Z500, UV200', 
                main='EKF400 5-yr warm periods in North America', units='K',
                addcontours=T, contvarname='gph500', conttype='data',
                contcol='black', contlev=contlevs,
                addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                veclen=scfac*0.02, vecscale=scfac*0.98, vecwd=0.95, every_x_vec=2)
    # diff echam analysis
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    pdata$data <- array(c((anameananom[,ti]-echmeananom[,ti]),
                          (anameananom[,ti+2]-echmeananom[,ti+2]),
                          (anameananom[,ti+4]-echmeananom[,ti+4]),
                          (anameananom[,ti+6]-echmeananom[,ti+6])),
                        c(nrow(echam.anom$ensmean),1,4))
    plot_echam3(pdata, varname='temp2', type='data', cex.pt=3.7, 
                latlim=c(20,60),lonlim=c(-160,-40),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('e)', 'f)', 'g)', 'h)'),colnames=rep(" ",4),
                statpanel=NULL, add=T, rownames='Update (Anal.-ECHAM)', 
                units='K',
                addcontours=T, contvarname='gph500', conttype='data',
                contcol='black', contlev=contlevs,
                addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                veclen=scfac*0.02, vecscale=scfac*0.98, vecwd=0.95, every_x_vec=2)
    pdata$data <- array(c(anameananom[,ti],anameananom[,ti+2],anameananom[,ti+4],
                          anameananom[,ti+6]), c(nrow(echam.anom$ensmean),1,4))
    levs=c(50,80,85,90,95,100,105,110,115,120,200)
    plot_echam3(pdata, varname='precip', type='data', cex.pt=3.7, 
                latlim=c(20,60),lonlim=c(-160,-40),units="%",
                names=c('e)', 'f)','g)', 'h)'), lev=levs, 
                st.col=NULL, stations=NULL, add=T,addcontours=F,
                wcol='darkgrey',rownames='Precipitation',colnames=rep('',4))
    pdata$data <- array(c(pdsimeananom[,1],pdsimeananom[,2],pdsimeananom[,3],
                          pdsimeananom[,4]), c(nrow(echam.anom$ensmean),1,4))
    levs <- seq(-4,4,1)
    plot_echam3(pdata, varname='precip', type='data', cex.pt=3.7, 
                latlim=c(20,60),lonlim=c(-160,-40),units="PDSI index",
                names=c('m)', 'n)','o)', 'p)'), lev=levs, 
                st.col=NULL, stations=NULL, add=T,addcontours=F,
                wcol='darkgrey',rownames='PDSI',colnames=rep('',4))
    dev.off()  
  } # end summer winter loop
} # end warm, cold periods US









if (warmcold_eu) {
  load(paste0('../data/indices/',expname,'/indices',filenameext,'allts.Rdata'))
 # find and mark coldest/warmest 5yr period in EU
  d <- aind.allts$ensmean[(which(eind.allts$names=="EU.temp2")),
                          (is.even(seq(1,ncol(aind.allts$ensmean))))]
  dec <- rep(NA,250)
  dectime <- rep(NA,250)
  for (ti in 1:length(dec)) {
    if (ti==1) {
      dec[ti] <- mean(d[1:8])
      dectime[ti] <- 1603
      i=2
    } else {
      dec[ti] <- mean(d[i:(i+7)])
      dectime[ti] <- dectime[ti-1] + 1
      i=i+1
    }
  }
  #  plot(dectime,dec,ty='l',col='red')
  tmp=cbind(dec,dectime)
  tmp[order(dec),]

  #coldper <- c(1695, 1641, 1674, 1739, 1713, 1701)
  coldper <- c(1810, 1695, 1640, 1832, 1740)
  for (syr in coldper) {
#    # calc validate clim (70 yr period) to substract for anomalies
#     if (syr>1645) {asyr<-syr-35} else {asyr=1610}
#     for (cyr in seq(asyr,(syr+39))) {  
#       load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
#       if (cyr==asyr){
#         vtmp <- validate$data
#         c=1
#       } else {
#         vtmp <- vtmp+validate$data
#         c <- c+1
#       }
#     }
#     valclim <- vtmp/c  
    eyr <- syr+7
    # change to precip in % 
    # load monthly 70yr climatology data
    if (syr < 1636) {
      load(file=paste0(echclimpath,'echam_clim_1636-1637_2ndgrid.Rdata'))
    } else {
      load(file=paste0(echclimpath,'echam_clim_',syr+3,'-',(syr+4),'_2ndgrid.Rdata'))  
    }
    pos <- which(echam_clim$names=="precip")[1:4608]
    # seasonal mean
    echmeanclim <- cbind(rowMeans(echam_clim$ensmean[,10:15]),
                         rowMeans(echam_clim$ensmean[,16:21]))
    #      echam$ensmean <- echmeanclim
    for (cyr in seq(syr,eyr)) {
      load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
      # change to precip in % 
      echam.abs <- echam
      echam.anom$ensmean[pos,] <- echam.abs$ensmean[pos,]/(echmeanclim[pos,]*3600*24*30)*100
      analysis.abs <- analysis
      analysis.anom$ensmean[pos,] <- analysis.abs$ensmean[pos,]/
                                        (echmeanclim[pos,]*3600*24*30)*100
      if (cyr==syr){
        atmp <- analysis.anom$ensmean
        etmp <- echam.anom$ensmean
#        vtmp <- validate$data
        c=1
      } else {
        atmp <- atmp+analysis.anom$ensmean
        etmp <- etmp+echam.anom$ensmean
#        vtmp ...        
        c=c+1
      }
    }
    if (syr==coldper[1]){
      anameananom <- atmp/c
      echmeananom <- etmp/c
#      valmean <- vtmp/c
#      valmeananom <- valmean-valclim
    } else {
      anameananom <- cbind(anameananom,(atmp/c))
      echmeananom <- cbind(echmeananom,(etmp/c))
#      valmean <- vtmp/c
#      valmeananom <- cbind(valmeananom,(valmean-valclim))
    }
  } # end cold periods loop  
  colnames(anameananom) <- rep(coldper,each=2)
  colnames(echmeananom) <- rep(coldper,each=2)

  # 8-yr cold period winter
  pdata <- echam
  for (ti in 1:2){
    if (ti==1){s="win"} else {s="som"}
    #   ti=2 # summer
    pdf(paste0('../figures/',expname,'/cold_eu_',s,filenameext,'.pdf'), width=12, 
        height=10.5, paper='special')
    #      layout(matrix(c(1,2,3,4,5,6,7,8), 4, 2, byrow = F), height=c(3,1,3,1))
    layout(matrix(c(1,2,3,4,5,5,5,5,6,7,8,9,10,10,10,10,11,12,13,14,15,15,15,15),
                  6,4,byrow=TRUE),height=c(3,1,3,1,3,1))
    par(oma=c(0.5,4,4,0))
#     pdata$data <- array(c(valmeananom[,1],valmeananom[,2],valmeananom[,3],
#                           valmeananom[,4]), c(nrow(echam.anom$ensmean),1,4))
#     levs <- seq(-4,4,1)
#     plot_echam3(pdata, varname='temp2', type='data', cex.pt=3.7, 
#                 latlim=c(30,75),lonlim=c(-20,40),units="K",
#                 names=c('m)', 'n)','o)', 'p)'), lev=levs, 
#                 st.col=NULL, stations=NULL, add=T,addcontours=F,
#                 wcol='darkgrey',rownames='Temp. reconstruction',colnames=rep('',4))
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    scfac <- max(echam.anom$ensmean[echam$names=='u200'])
    pdata$data <- array(c(anameananom[,ti],anameananom[,ti+2],anameananom[,ti+4],
                          anameananom[,ti+6]), c(nrow(echam.anom$ensmean),1,4))
    plot_echam3(pdata, varname='temp2', type='data', cex.pt=3.7, 
                latlim=c(30,75),lonlim=c(-20,40),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('a)', 'b)', 'c)', 'd)'),colnames=coldper[1:4],
                statpanel=NULL, add=T, rownames='Analysis T2m, Z500, UV200', 
                main='EKF400 5-yr cold periods in Europe', units='K',
                addcontours=T, contvarname='gph500', conttype='data',
                contcol='black', contlev=contlevs,
                addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                veclen=scfac*0.02, vecscale=scfac*0.98, vecwd=0.95, every_x_vec=2)
    # diff echam analysis
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    pdata$data <- array(c((anameananom[,ti]-echmeananom[,ti]),
                          (anameananom[,ti+2]-echmeananom[,ti+2]),
                          (anameananom[,ti+4]-echmeananom[,ti+4]),
                          (anameananom[,ti+6]-echmeananom[,ti+6])),
                        c(nrow(echam.anom$ensmean),1,4))
    plot_echam3(pdata, varname='temp2', type='data', cex.pt=3.7, 
                latlim=c(30,75),lonlim=c(-20,40),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('e)', 'f)', 'g)', 'h)'),colnames=rep(" ",4),
                statpanel=NULL, add=T, rownames='Update (Anal.-ECHAM)', 
                units='K',
                addcontours=T, contvarname='gph500', conttype='data',
                contcol='black', contlev=contlevs,
                addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                veclen=scfac*0.02, vecscale=scfac*0.98, vecwd=0.95, every_x_vec=2)
    pdata$data <- array(c(anameananom[,ti],anameananom[,ti+2],anameananom[,ti+4],
                          anameananom[,ti+6]), c(nrow(echam.anom$ensmean),1,4))
    #      #    levs <- c(-Inf, seq(-5,5,1), Inf)
    levs=c(50,80,85,90,95,100,105,110,115,120,200)
    plot_echam3(pdata, varname='precip', type='data', cex.pt=3.7, 
                latlim=c(30,75),lonlim=c(-20,40),units="%",
                names=c('i)', 'j)','k)', 'l)'), lev=levs, 
                st.col=NULL, stations=NULL, add=T,addcontours=F,
                wcol='darkgrey',rownames='Prec. analysis',colnames=rep('',4))
#     pdata$data <- array(c(valmeananom[,1],valmeananom[,2],valmeananom[,3],
#                           valmeananom[,4]), c(nrow(echam.anom$ensmean),1,4))
#     levs <- seq(-4,4,1)
#     plot_echam3(pdata, varname='precip', type='data', cex.pt=3.7, 
#                 latlim=c(30,75),lonlim=c(-20,40),units="mm",
#                 names=c('m)', 'n)','o)', 'p)'), lev=levs, 
#                 st.col=NULL, stations=NULL, add=T,addcontours=F,
#                 wcol='darkgrey',rownames='Prec. reconstruction',colnames=rep('',4))
    dev.off()  
  } # end summer winter loop

  
  #warmper <- c(1794, 1800, 1846, 1826, 1819, 1779)
  warmper <- c(1654, 1794, 1726, 1774, 1846, 1800)
  for (syr in warmper) {
    eyr <- syr+7
    # change to precip in % 
    # load monthly 70yr climatology data
    if (syr < 1636) {
      load(file=paste0(echclimpath,'echam_clim_1636-1637_2ndgrid.Rdata'))
    } else {
      load(file=paste0(echclimpath,'echam_clim_',syr+3,'-',(syr+4),'_2ndgrid.Rdata'))  
    }
    pos <- which(echam_clim$names=="precip")[1:4608]
    # seasonal mean
    echmeanclim <- cbind(rowMeans(echam_clim$ensmean[,10:15]),
                         rowMeans(echam_clim$ensmean[,16:21]))
    for (cyr in seq(syr,eyr)) {
      load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
      # change to precip in % 
      echam.abs <- echam
      echam.anom$ensmean[pos,] <- echam.abs$ensmean[pos,]/(echmeanclim[pos,]*3600*24*30)*100
      analysis.abs <- analysis
      analysis.anom$ensmean[pos,] <- analysis.abs$ensmean[pos,]/
                                        (echmeanclim[pos,]*3600*24*30)*100
      if (cyr==syr){
        atmp <- analysis.anom$ensmean
        etmp <- echam.anom$ensmean
        c=1
      } else {
        atmp <- atmp+analysis.anom$ensmean
        etmp <- etmp+echam.anom$ensmean
        c=c+1
      }
    }
    if (syr==warmper[1]){
      anameananom <- atmp/c
      echmeananom <- etmp/c
    } else {
      anameananom <- cbind(anameananom,(atmp/c))
      echmeananom <- cbind(echmeananom,(etmp/c))
    }
  } # end warm periods loop  
  colnames(anameananom) <- rep(warmper,each=2)
  colnames(echmeananom) <- rep(warmper,each=2)
  
  # 8-yr warm period winter
  pdata <- echam
  for (ti in 1:2){
    if (ti==1){s="win"} else {s="som"}
    pdf(paste0('../figures/',expname,'/warm_eu_',s,filenameext,'.pdf'), width=12, 
        height=10.5, paper='special')
    layout(matrix(c(1,2,3,4,5,5,5,5,6,7,8,9,10,10,10,10,11,12,13,14,15,15,15,15),
                  6,4,byrow=TRUE),height=c(3,1,3,1,3,1))
    par(oma=c(0.5,4,4,0))
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    scfac <- max(echam.anom$ensmean[echam$names=='u200'])
    pdata$data <- array(c(anameananom[,ti],anameananom[,ti+2],anameananom[,ti+4],
                          anameananom[,ti+6]), c(nrow(echam.anom$ensmean),1,4))
    plot_echam3(pdata, varname='temp2', type='data', cex.pt=3.7, 
                latlim=c(30,75),lonlim=c(-20,40),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('a)', 'b)', 'c)', 'd)'),colnames=warmper[1:4],
                statpanel=NULL, add=T, rownames='Analysis T2m, Z500, UV200', 
                main='EKF400 5-yr warm periods in Europe', units='K',
                addcontours=T, contvarname='gph500', conttype='data',
                contcol='black', contlev=contlevs,
                addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                veclen=scfac*0.02, vecscale=scfac*0.98, vecwd=0.95, every_x_vec=2)
    # diff echam analysis
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    pdata$data <- array(c((anameananom[,ti]-echmeananom[,ti]),
                          (anameananom[,ti+2]-echmeananom[,ti+2]),
                          (anameananom[,ti+4]-echmeananom[,ti+4]),
                          (anameananom[,ti+6]-echmeananom[,ti+6])),
                        c(nrow(echam.anom$ensmean),1,4))
    plot_echam3(pdata, varname='temp2', type='data', cex.pt=3.7, 
                latlim=c(30,75),lonlim=c(-20,40),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('e)', 'f)', 'g)', 'h)'),colnames=rep(" ",4),
                statpanel=NULL, add=T, rownames='Update (Anal.-ECHAM)', 
                units='K',
                addcontours=T, contvarname='gph500', conttype='data',
                contcol='black', contlev=contlevs,
                addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                veclen=scfac*0.02, vecscale=scfac*0.98, vecwd=0.95, every_x_vec=2)
    pdata$data <- array(c(anameananom[,ti],anameananom[,ti+2],anameananom[,ti+4],
                          anameananom[,ti+6]), c(nrow(echam.anom$ensmean),1,4))
    levs=c(50,80,85,90,95,100,105,110,115,120,200)
    plot_echam3(pdata, varname='precip', type='data', cex.pt=3.7, 
                latlim=c(30,75),lonlim=c(-20,40),units="%",
                names=c('i)', 'j)','k)', 'l)'), lev=levs, 
                st.col=NULL, stations=NULL, add=T,addcontours=F,
                wcol='darkgrey',rownames='Prec. analysis',colnames=rep('',4))
    dev.off()  
  } # end summer winter loop
} # end warm, cold periods EU











if (warmcold_nh) {
  # oct 2016 10 periods 1600-2000 + composite
  load(paste0('../data/indices/',expname,'/indices',filenameext,'allts.Rdata'))
   # find and mark warm/coldest 8yr periods in NH based on 70yr anomalies (land and ocean)
  d <- aind.allts$ensmean[(which(eind.allts$names=="NH.temp2")),
                          (is.even(seq(1,ncol(aind.allts$ensmean))))]
  dec <- rep(NA,390) #250)
  dectime <- rep(NA,390) #250)
  for (ti in 1:length(dec)) {
    if (ti==1) {
      dec[ti] <- mean(d[1:8])
      dectime[ti] <- 1603
      i=2
    } else {
      dec[ti] <- mean(d[i:(i+7)])
      dectime[ti] <- dectime[ti-1] + 1
      i=i+1
    }
  }
  #  plot(dectime,dec,ty='l',col='red')
  tmp=cbind(dec,dectime)
  cold <- tmp
  warm <- tmp
  for (i in 1:10) { # for 10 coldest/warmest periods
    pos_cold <- which(cold[,1]==min(cold[,1],na.rm=T))
    if (i==1) { 
      coldper <- cold[pos_cold,] 
    } else {
      coldper <- rbind(coldper,cold[pos_cold,]) 
    }
    cold[(pos_cold-8):(pos_cold+8),1] <- NA
    pos_warm <- which(warm[,1]==max(warm[,1],na.rm=T))
    if (i==1) { 
      warmper <- warm[pos_warm,] 
    } else {
      warmper <- rbind(warmper,warm[pos_warm,]) 
    }
    if (pos_warm>382) {
      warm[(pos_warm-8):390] <- NA
    } else {
      warm[(pos_warm-8):(pos_warm+8),1] <- NA  
    }
  }
  
  #load absolute indices for plotting  
  load(paste0('../data/indices/',expname,'/indices','_abs_landocean_seas_','allts.Rdata'))
  #load(paste0('../data/indices/indices','_abs_landonly_seas_','allts.Rdata'))
  
  pos_warm <- which(rep(eind.allts$time,each=2) %in%
                      (sort(warmper[,2]+rep(seq(0,7),each=nrow(warmper)))))
  pos_cold <- which(rep(eind.allts$time,each=2) %in%
                      (sort(coldper[,2]+rep(seq(0,7),each=nrow(coldper)))))
  #warm_ind_ann <- array(aind.allts$ensmean[,pos_warm],dim=c(nrow(aind.allts$data),
  #                  dim(aind.allts$data)[2])) #,1,mean)
  warm_ind_ann <- array(aind.allts$data[,pos_warm,],dim=c(nrow(aind.allts$data),
                    dim(aind.allts$data)[2]*dim(aind.allts$data)[3])) #,1,mean)
  warm_ind_sum <- array(aind.allts$data[,pos_warm[is.even(pos_warm)],],dim=c(nrow(aind.allts$data),
                    dim(aind.allts$data)[2]*dim(aind.allts$data)[3])) #,1,mean)
  warm_ind_win <- array(aind.allts$data[,pos_warm[!is.even(pos_warm)],],dim=c(nrow(aind.allts$data),
                    dim(aind.allts$data)[2]*dim(aind.allts$data)[3])) #,1,mean)
  #cold_ind_ann <- array(aind.allts$ensmean[,pos_cold],dim=c(nrow(aind.allts$data),
  #                  dim(aind.allts$data)[2]))# ,1,mean)
  
  cold_ind_ann <- array(aind.allts$data[,pos_cold,],dim=c(nrow(aind.allts$data),
                    dim(aind.allts$data)[2]*dim(aind.allts$data)[3]))# ,1,mean)
  cold_ind_sum <- array(aind.allts$data[,pos_cold[is.even(pos_cold)],],dim=c(nrow(aind.allts$data),
                    dim(aind.allts$data)[2]*dim(aind.allts$data)[3]))# ,1,mean)
  cold_ind_win <- array(aind.allts$data[,pos_cold[!is.even(pos_cold)],],dim=c(nrow(aind.allts$data),
                    dim(aind.allts$data)[2]*dim(aind.allts$data)[3]))#,1,mean)
  
  #quartz()
  pdf(paste0('../figures/',expname,'/warm_cold_indices_abs_landocean_seas.pdf'), width=10,
  #pdf(paste0('../figures/decadal_warm_cold/warm_cold_indices_abs_landocean_seas_summer.pdf'), width=10,
  #pdf(paste0('../figures/decadal_warm_cold/warm_cold_indices_abs_landocean_seas_winter.pdf'), width=10,
  #pdf(paste0('../figures/decadal_warm_cold/warm_cold_indices_abs_landonly_seas.pdf'), width=10,     
      height=7, paper='special')
    par(mfrow=c(2,5))
    for (i in seq(35,43)) {
    #for (i in c(35,36,37,38,39,41)) {
      boxplot(cbind(warm_ind_ann[i,],cold_ind_ann[i,]),main=eind.allts$names[i],
      #boxplot(cbind(warm_ind_sum[i,],cold_ind_sum[i,]),main=eind.allts$names[i],
      #boxplot(cbind(warm_ind_win[i,],cold_ind_win[i,]),main=eind.allts$names[i],
      #paste(eind.allts$names[i],"p-value:",round(as.numeric(t.test(warm_ind_ann[i,],
      #cold_ind_ann[i,])[3]),digits=2)),
      las=2,col=c('red','blue'),outline=F)
    }
  dev.off()
  #beanplot(t(data.frame(warm_ind_ann[c(1:4),1:100])))

  
#  coldper <- c(1641, 1695, 1676, 1815, 1667, 1701, 1738)
#  coldper <- c(1641, 1695, 1815, 1832) #, 1676, 1667)
#  coldper <- c(1812, 1641, 1831, 1695, 1739) 
  coldper <- (as.vector(coldper[,2]))
  for (syr in coldper) {
    eyr <- syr+7
    # change to precip in % 
    # load monthly 70yr climatology data
    if (syr < 1636) {
      load(file=paste0(echclimpath,'echam_clim_1636-1637_2ndgrid.Rdata'))
    } else {
      load(file=paste0(echclimpath,'echam_clim_',syr+3,'-',(syr+4),'_2ndgrid.Rdata'))  
    }
    pos <- which(echam_clim$names=="precip")[1:4608]
    # seasonal mean
    echmeanclim <- cbind(rowMeans(echam_clim$ensmean[,10:15]),
                         rowMeans(echam_clim$ensmean[,16:21]))
    for (cyr in seq(syr,eyr)) {
      load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
      # change to precip in % 
      echam.abs <- echam
      echam.anom$ensmean[pos,] <- echam.abs$ensmean[pos,]/(echmeanclim[pos,]*3600*24*30)*100
      analysis.abs <- analysis
      analysis.anom$ensmean[pos,] <- analysis.abs$ensmean[pos,]/
                                        (echmeanclim[pos,]*3600*24*30)*100
      if (cyr==syr){
        atmp <- analysis.anom$ensmean
        etmp <- echam.anom$ensmean
        c=1
      } else {
        atmp <- atmp+analysis.anom$ensmean
        etmp <- etmp+echam.anom$ensmean
        c=c+1
      }
    }
    # read mann et al temp recon
    mann <- read.table("/Users/joerg/Documents/climdata/mann_temp_recon/allproxyfieldrecon",
                       na.string="NaN")
    ll <- read.table("/Users/joerg/Documents/climdata/mann_temp_recon/longlat.txt",
                     na.string="NaN")
    ts <- which(mann[,1]==syr)
    plotmann <- list() #echam
    mclim <- as.vector(t(colMeans(mann[((ts-35):(ts+39)),-1],na.rm=F)))
    mdata <- as.vector(t(colMeans(mann[ts:(ts+4),-1],na.rm=F)))
    manom <- mdata-mclim
    plotmann$data <- array(manom,c(length(manom),1))
    plotmann$time <- mann[ts,1]
    plotmann$names <- rep('temp2',dim(plotmann$data)[1])
    plotmann$lon <- ll[,1] 
    plotmann$lat <- ll[,2]
    plotmann$height <- ll[,3]   
    plotmann$lsm.i <- ll[,3]  
    if (syr==coldper[1]){
      anameananom <- atmp/c
      echmeananom <- etmp/c
      mannmeananom <- plotmann$data
    } else {
      anameananom <- cbind(anameananom,(atmp/c))
      echmeananom <- cbind(echmeananom,(etmp/c))
      mannmeananom <- cbind(mannmeananom,plotmann$data)
    }
  } # end cold periods loop  
  colnames(anameananom) <- rep(coldper,each=2)
  colnames(echmeananom) <- rep(coldper,each=2)
  colnames(mannmeananom) <- coldper
  
  # 8-yr cold period winter
  pdata <- echam
  for (ti in 1:2){
    if (ti==1){s="win"} else {s="som"}
    pdf(paste0('../figures/',expname,'/10_cold_nh_',s,filenameext,'.pdf'), 
        width=40, height=12, paper='special')
#    layout(matrix(c(1,2,3,4,5,5,5,5,6,7,8,9,10,10,10,10,11,12,13,14,15,15,15,15,
#                    16,17,18,19,20,20,20,20),8,4,byrow=TRUE),height=c(3,1,3,1,3,1,3,1))
#    layout(matrix(rep(seq(1,44,1),rep(c(rep(1,10),10),4)),8,10,byrow=TRUE),
#           height=rep(c(3,1),4))
    layout(matrix(rep(seq(1,84,1),rep(c(rep(1,20),20),4)),8,20,byrow=TRUE),
           height=rep(c(3,1),4),width=rep(c(4,1),10))
    par(oma=c(0.5,2,2,0))
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    #plotmann$data <- array(c(mannmeananom[,1],mannmeananom[,2],mannmeananom[,3],
    #                         mannmeananom[,4]), c(nrow(plotmann$data),1,4))
    #plotmann$data <- array(mannmeananom[,],c(nrow(mannmeananom),1,ncol(mannmeananom)))
    plotmann$data <- array(c(mannmeananom[,1],mannmeananom[,2],mannmeananom[,3],
                             mannmeananom[,4],mannmeananom[,5],mannmeananom[,6],
                             mannmeananom[,7],mannmeananom[,8],mannmeananom[,9],
                             rowMeans(mannmeananom[,1:9])), c(nrow(plotmann$data),1,10))
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    plot_echam4(plotmann, varname='temp2', type='data', cex.pt=2.5, 
                latlim=c(0,90),lonlim=c(-180,180),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('a)', 'b)', 'c)', 'd)'),colnames=c(coldper[1:9],"composite"),
                statpanel=NULL, add=T, rownames='Mann recon. temp.', 
                main='', units='K',
                zonalmean=T,zmvarname='temp2')
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    scfac <- max(echam.anom$ensmean[echam$names=='u200'])
    pdata$data <- array(c(anameananom[,ti],anameananom[,ti+2],anameananom[,ti+4],
                        anameananom[,ti+6],anameananom[,ti+8],anameananom[,ti+10],
                        anameananom[,ti+12],anameananom[,ti+14],anameananom[,ti+16],
                        rowMeans(anameananom[,c(ti,ti+2,ti+4,ti+6,ti+8,ti+10,ti+12,
                        ti+14,ti+16)])),
                        c(nrow(echam.anom$ensmean),1,10))
    plot_echam4(pdata, varname='temp2', type='data', cex.pt=1.5, 
                latlim=c(0,90),lonlim=c(-180,180),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('e)', 'f)', 'g)', 'h)'),colnames=rep("",4),
                statpanel=NULL, add=T, rownames='Analysis temp, Z500, UV200', 
                main='EKF400 5-yr cold periods', units='K',
                addcontours=T, contvarname='gph500', conttype='data',
                contcol='black', contlev=contlevs,
                addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                veclen=scfac*0.01, vecscale=scfac*2.0, vecwd=0.95, every_x_vec=4,
                zonalmean=T,zmvarname='gph500')
    # diff echam analysis
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    pdata$data <- array(c((anameananom[,ti]-echmeananom[,ti]),
                          (anameananom[,ti+2]-echmeananom[,ti+2]),
                          (anameananom[,ti+4]-echmeananom[,ti+4]),
                          (anameananom[,ti+6]-echmeananom[,ti+6]),
                          (anameananom[,ti+8]-echmeananom[,ti+8]),
                          (anameananom[,ti+10]-echmeananom[,ti+10]),
                          (anameananom[,ti+12]-echmeananom[,ti+12]),
                          (anameananom[,ti+14]-echmeananom[,ti+14]),
                          (anameananom[,ti+16]-echmeananom[,ti+16]),
                          (rowMeans(anameananom[,c(ti,ti+2,ti+4,ti+6,
                          ti+8,ti+10,ti+12,ti+14,ti+16)])-
                          rowMeans(echmeananom[,c(ti,ti+2,ti+4,ti+6,
                          ti+8,ti+10,ti+12,ti+14,ti+16)]))),
                          c(nrow(echam.anom$ensmean),1,10))
    plot_echam4(pdata, varname='temp2', type='data', cex.pt=1.5, 
                latlim=c(0,90),lonlim=c(-180,180),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('i)', 'j)', 'k)', 'l)'),colnames=rep(" ",4),
                statpanel=NULL, add=T, rownames='Update (Anal.-ECHAM)', 
                units='K',
                addcontours=T, contvarname='gph500', conttype='data',
                contcol='black', contlev=contlevs,
                addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                veclen=scfac*0.01, vecscale=scfac*2.0, vecwd=0.95, every_x_vec=4,
                zonalmean=T,zmvarname='gph500')
    pdata$data <- array(c(anameananom[,ti],anameananom[,ti+2],anameananom[,ti+4],
                          anameananom[,ti+6],anameananom[,ti+8],anameananom[,ti+10],
                          anameananom[,ti+12],anameananom[,ti+14],anameananom[,ti+16],
                          rowMeans(anameananom[,c(ti,ti+2,ti+4,ti+6,ti+8,ti+10,ti+12,
                          ti+14,ti+16)])),c(nrow(echam.anom$ensmean),1,10))
    levs=c(0,80,85,90,95,100,105,110,115,120,Inf)
    plot_echam4(pdata, varname='precip', type='data', cex.pt=1.5, 
                latlim=c(0,90),lonlim=c(-180,180), units="%",
                names=c('m)', 'n)','o)', 'p)'), lev=levs, 
                st.col=NULL, stations=NULL, add=T,addcontours=F,
                wcol='darkgrey',rownames='Prec. analysis',colnames=rep('',4),
                zonalmean=T,zmvarname='precip')
    dev.off()  
  } # end summer winter loop
  

  
  #warmper <- c(1801, 1826, 1791, 1778, 1844, 1779)
  #warmper <- c(1656, 1798, 1686, 1722, 1849)
  warmper <- (as.vector(warmper[,2]))
  for (syr in warmper) {
    eyr <- syr+7
    # load monthly 70yr climatology data
    if (syr < 1636) {
      load(file=paste0(echclimpath,'echam_clim_1636-1637_2ndgrid.Rdata'))
    } else {
      load(file=paste0(echclimpath,'echam_clim_',syr+3,'-',(syr+4),'_2ndgrid.Rdata'))  
    }
    pos <- which(echam_clim$names=="precip")[1:4608]
    # seasonal mean
    echmeanclim <- cbind(rowMeans(echam_clim$ensmean[,10:15]),
                         rowMeans(echam_clim$ensmean[,16:21]))
    for (cyr in seq(syr,eyr)) {
      load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
      # change to precip in % 
      echam.abs <- echam
      echam.anom$ensmean[pos,] <- echam.abs$ensmean[pos,]/(echmeanclim[pos,]*3600*24*30)*100
      analysis.abs <- analysis
      analysis.anom$ensmean[pos,] <- analysis.abs$ensmean[pos,]/
                                        (echmeanclim[pos,]*3600*24*30)*100
      if (cyr==syr){
        atmp <- analysis.anom$ensmean
        etmp <- echam.anom$ensmean
        c=1
      } else {
        atmp <- atmp+analysis.anom$ensmean
        etmp <- etmp+echam.anom$ensmean
        c=c+1
      }
    }
    # read mann et al temp recon
    mann <- read.table("/Users/joerg/Documents/climdata/mann_temp_recon/allproxyfieldrecon",
                       na.string="NaN")
    ll <- read.table("/Users/joerg/Documents/climdata/mann_temp_recon/longlat.txt",
                     na.string="NaN")
    ts <- which(mann[,1]==syr)
    plotmann <- list() #echam
    mclim <- as.vector(t(colMeans(mann[((ts-35):(ts+39)),-1],na.rm=F)))
    mdata <- as.vector(t(colMeans(mann[ts:(ts+4),-1],na.rm=F)))
    manom <- mdata-mclim
    plotmann$data <- array(manom,c(length(manom),1))
    plotmann$time <- mann[ts,1]
    plotmann$names <- rep('temp2',dim(plotmann$data)[1])
    plotmann$lon <- ll[,1] 
    plotmann$lat <- ll[,2]
    plotmann$height <- ll[,3]   
    plotmann$lsm.i <- ll[,3]  
    if (syr==warmper[1]){
      anameananom <- atmp/c
      echmeananom <- etmp/c
      mannmeananom <- plotmann$data
    } else {
      anameananom <- cbind(anameananom,(atmp/c))
      echmeananom <- cbind(echmeananom,(etmp/c))
      mannmeananom <- cbind(mannmeananom,plotmann$data)
    }
  } # end warm periods loop  
  colnames(anameananom) <- rep(warmper,each=2)
  colnames(echmeananom) <- rep(warmper,each=2)
  colnames(mannmeananom) <- warmper

  # 8-yr warm period winter
  pdata <- echam
  for (ti in 1:2){
    if (ti==1){s="win"} else {s="som"}
    pdf(paste0('../figures/',expname,'/10_warm_nh_',s,filenameext,'.pdf'), width=40, 
        height=12, paper='special')
    #layout(matrix(c(1,2,3,4,5,5,5,5,6,7,8,9,10,10,10,10,11,12,13,14,15,15,15,15,
    #                16,17,18,19,20,20,20,20),8,4,byrow=TRUE),height=c(3,1,3,1,3,1,3,1))
    #layout(matrix(rep(seq(1,44,1),rep(c(rep(1,10),10),4)),8,10,byrow=TRUE),height=rep(c(3,1),4))
    #par(oma=c(0.5,4,4,0))
    layout(matrix(rep(seq(1,84,1),rep(c(rep(1,20),20),4)),8,20,byrow=TRUE),
           height=rep(c(3,1),4),width=rep(c(4,1),10))
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    plotmann$data <- array(c(mannmeananom[,1],mannmeananom[,2],mannmeananom[,3],mannmeananom[,4],
                             mannmeananom[,5],mannmeananom[,6],mannmeananom[,7],mannmeananom[,8],
                             mannmeananom[,9],rowMeans(mannmeananom[,1:9])),
                             c(nrow(plotmann$data),1,10))
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
#    plot_echam3(plotmann, varname='temp2', type='data', cex.pt=1.5, names='',
#                lev=levs, st.col=NULL, stations=NULL, add=T,units='ºC')
    plot_echam4(plotmann, varname='temp2', type='data', cex.pt=2.5, 
                latlim=c(0,90),lonlim=c(-180,180),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('a)', 'b)', 'c)', 'd)'),colnames=c(warmper[1:9],"composite"),
                statpanel=NULL, add=T, rownames='Mann recon. temp.', 
                main='', units='K',
                zonalmean=T,zmvarname='temp2')
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    scfac <- max(echam.anom$ensmean[echam$names=='u200'])
    pdata$data <- array(c(anameananom[,ti],anameananom[,ti+2],anameananom[,ti+4],
                          anameananom[,ti+6],anameananom[,ti+8],anameananom[,ti+10],
                          anameananom[,ti+12],anameananom[,ti+14],anameananom[,ti+16],
                          rowMeans(anameananom[,c(ti,ti+2,ti+4,ti+6,ti+8,ti+10,ti+12,
                          ti+14,ti+16)])),c(nrow(echam.anom$ensmean),1,10))
    plot_echam4(pdata, varname='temp2', type='data', cex.pt=1.5, 
                latlim=c(0,90),lonlim=c(-180,180),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('e)', 'f)', 'g)', 'h)'),colnames=rep("",4),
                statpanel=NULL, add=T, rownames='Analysis temp, Z500, UV200', 
                main='EKF400 5-yr warm periods', units='K',
                addcontours=T, contvarname='gph500', conttype='data',
                contcol='black', contlev=contlevs,
                addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                veclen=scfac*0.01, vecscale=scfac*2.0, vecwd=0.95, every_x_vec=4,
                zonalmean=T,zmvarname='gph500')
    # diff echam analysis
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    pdata$data <- array(c((anameananom[,ti]-echmeananom[,ti]),
                          (anameananom[,ti+2]-echmeananom[,ti+2]),
                          (anameananom[,ti+4]-echmeananom[,ti+4]),
                          (anameananom[,ti+6]-echmeananom[,ti+6]),
                          (anameananom[,ti+8]-echmeananom[,ti+8]),
                          (anameananom[,ti+10]-echmeananom[,ti+10]),
                          (anameananom[,ti+12]-echmeananom[,ti+12]),
                          (anameananom[,ti+14]-echmeananom[,ti+14]),
                          (anameananom[,ti+16]-echmeananom[,ti+16]),
                          (rowMeans(anameananom[,c(ti,ti+2,ti+4,ti+6,
                           ti+8,ti+10,ti+12,ti+14,ti+16)])-
                           rowMeans(echmeananom[,c(ti,ti+2,ti+4,ti+6,
                           ti+8,ti+10,ti+12,ti+14,ti+16)]))),
                        c(nrow(echam.anom$ensmean),1,10))
    plot_echam4(pdata, varname='temp2', type='data', cex.pt=1.5, 
                latlim=c(0,90),lonlim=c(-180,180),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('i)', 'j)', 'k)', 'l)'),colnames=rep(" ",4),
                statpanel=NULL, add=T, rownames='Update (Anal.-ECHAM)', 
                units='K',
                addcontours=T, contvarname='gph500', conttype='data',
                contcol='black', contlev=contlevs,
                addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                veclen=scfac*0.01, vecscale=scfac*2.0, vecwd=0.95, every_x_vec=4,
                zonalmean=T,zmvarname='gph500')
    pdata$data <- array(c(anameananom[,ti],anameananom[,ti+2],anameananom[,ti+4],
                          anameananom[,ti+6],anameananom[,ti+8],anameananom[,ti+10],
                          anameananom[,ti+12],anameananom[,ti+14],anameananom[,ti+16],
                          rowMeans(anameananom[,c(ti,ti+2,ti+4,ti+6,ti+8,ti+10,ti+12,
                          ti+14,ti+16)])),c(nrow(echam.anom$ensmean),1,10))
    levs=c(0,80,85,90,95,100,105,110,115,120,Inf)
    plot_echam4(pdata, varname='precip', type='data', cex.pt=1.5, 
                latlim=c(0,90),lonlim=c(-180,180),units="%",
                names=c('m)', 'n)','o)', 'p)'), lev=levs, 
                st.col=NULL, stations=NULL, add=T,addcontours=F,
                wcol='darkgrey',rownames='Prec. analysis',colnames=rep('',4),
                zonalmean=T,zmvarname='precip')
    dev.off()  
  } # end summer winter loop
} # end warm, cold periods NH


  
  

  




if (warmcold_1790_1830) {
    load(paste0('../data/indices/',expname,'/indices',filenameext,'allts.Rdata'))
    # find and mark warm/coldest 8yr periods in NH based on 70yr summer anomalies (land and ocean)
    d <- aind.allts$ensmean[(which(eind.allts$names=="NH.temp2")),
                            (is.even(seq(1,ncol(aind.allts$ensmean))))]
    dec <- rep(NA,250)
    dectime <- rep(NA,250)
    for (ti in 1:length(dec)) {
      if (ti==1) {
        dec[ti] <- mean(d[1:8])
        dectime[ti] <- 1603
        i=2
      } else {
        dec[ti] <- mean(d[i:(i+7)])
        dectime[ti] <- dectime[ti-1] + 1
        i=i+1
      }
    }
    #  plot(dectime,dec,ty='l',col='red')
    tmp=cbind(dec,dectime)
    tmp[order(dec),]
    
    # find and mark warm/coldest 8yr periods in NH based on 70yr ANNUAL anomalies (land and ocean)
    d <- aggregate(aind.allts$ensmean[(which(eind.allts$names=="NH.temp2")),],
                   list(rep(aind.allts$time,each=2)),mean)[,2]
    dec <- rep(NA,250)
    dectime <- rep(NA,250)
    for (ti in 1:length(dec)) {
      if (ti==1) {
        dec[ti] <- mean(d[1:8])
        dectime[ti] <- 1603
        i=2
      } else {
        dec[ti] <- mean(d[i:(i+7)])
        dectime[ti] <- dectime[ti-1] + 1
        i=i+1
      }
    }
    #  plot(dectime,dec,ty='l',col='red')
    tmp=cbind(dec,dectime)
    tmp[order(dec),]
    
    warmper <- c(1775, 1791, 1798, 1824) 
    for (syr in warmper) {
      eyr <- syr+7
      # change to precip in % 
      # load monthly 70yr climatology data
      if (syr < 1636) {
        load(file=paste0(echclimpath,'echam_clim_1636-1637_2ndgrid.Rdata'))
      } else {
        load(file=paste0(echclimpath,'echam_clim_',syr+3,'-',(syr+4),'_2ndgrid.Rdata'))  
      }
      pos <- which(echam_clim$names=="precip")[1:4608]
      # seasonal mean
      echmeanclim <- cbind(rowMeans(echam_clim$ensmean[,10:15]),
                           rowMeans(echam_clim$ensmean[,16:21]))
      for (cyr in seq(syr,eyr)) {
        load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
        # change to precip in % 
        echam.abs <- echam
        echam.anom$ensmean[pos,] <- echam.abs$ensmean[pos,]/(echmeanclim[pos,]*3600*24*30)*100
        analysis.abs <- analysis
        analysis.anom$ensmean[pos,] <- analysis.abs$ensmean[pos,]/
          (echmeanclim[pos,]*3600*24*30)*100
        if (cyr==syr){
          atmp <- analysis.anom$ensmean
          etmp <- echam.anom$ensmean
          c=1
        } else {
          atmp <- atmp+analysis.anom$ensmean
          etmp <- etmp+echam.anom$ensmean
          c=c+1
        }
      }
      # read mann et al temp recon
      mann <- read.table(paste0(datadir,"mann_temp_recon/allproxyfieldrecon"),na.string="NaN")
      ll <- read.table(paste0(datadir,"mann_temp_recon/longlat.txt"),na.string="NaN")
      ts <- which(mann[,1]==syr)
      plotmann <- list() #echam
      mclim <- as.vector(t(colMeans(mann[((ts-35):(ts+39)),-1],na.rm=F)))
      mdata <- as.vector(t(colMeans(mann[ts:(ts+4),-1],na.rm=F)))
      manom <- mdata-mclim
      plotmann$data <- array(manom,c(length(manom),1))
      plotmann$time <- mann[ts,1]
      plotmann$names <- rep('temp2',dim(plotmann$data)[1])
      plotmann$lon <- ll[,1] 
      plotmann$lat <- ll[,2]
      plotmann$height <- ll[,3]   
      plotmann$lsm.i <- ll[,3]  
      if (syr==warmper[1]){
        anameananom <- atmp/c
        echmeananom <- etmp/c
        mannmeananom <- plotmann$data
      } else {
        anameananom <- cbind(anameananom,(atmp/c))
        echmeananom <- cbind(echmeananom,(etmp/c))
        mannmeananom <- cbind(mannmeananom,plotmann$data)
      }
    } # end periods loop  
    colnames(anameananom) <- rep(warmper,each=2)
    colnames(echmeananom) <- rep(warmper,each=2)
    colnames(mannmeananom) <- warmper
    
    # 8-yr period winter
    pdata <- echam
    for (ti in 1:3){
      if (ti==1){s="winter"} else if (ti==2){s="summer"} else {s="annual"}
      pdf(paste0('../figures/',expname,'/1740-1840_warm_',s,filenameext,'.pdf'), width=16, 
          height=12, paper='special')
      layout(matrix(c(1,2,3,4,5,5,5,5,6,7,8,9,10,10,10,10,11,12,13,14,15,15,15,15,
                      16,17,18,19,20,20,20,20),8,4,byrow=TRUE),height=c(3,1,3,1,3,1,3,1))
      par(oma=c(0.5,4,4,0))
      levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
      plotmann$data <- array(c(mannmeananom[,1],mannmeananom[,2],mannmeananom[,3],
                               mannmeananom[,4]), c(nrow(plotmann$data),1,4))
      levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
      plot_echam3(plotmann, varname='temp2', type='data', cex.pt=2.5, 
                  latlim=c(0,90),lonlim=c(-180,180),
                  wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                  names=c('a)', 'b)', 'c)', 'd)'),colnames=warmper[1:4],
                  statpanel=NULL, add=T, rownames='Mann recon. temp.', 
                  main='', units='K')
      levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
      contlevs <- seq(-4,4,2)
      scfac <- max(echam.anom$ensmean[echam$names=='u200'])
      if (ti==3){
        pdata$data <- array(c(rowMeans(anameananom[,1:2]),rowMeans(anameananom[,3:4]),
                              rowMeans(anameananom[,5:6]),rowMeans(anameananom[,7:8])),
                              c(nrow(echam.anom$ensmean),1,4))
      } else {
        pdata$data <- array(c(anameananom[,ti],anameananom[,ti+2],anameananom[,ti+4],
                            anameananom[,ti+6]), c(nrow(echam.anom$ensmean),1,4))
      }
      plot_echam3(pdata, varname='temp2', type='data', cex.pt=1.5, 
                  latlim=c(0,90),lonlim=c(-180,180),
                  wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                  names=c('e)', 'f)', 'g)', 'h)'),colnames=rep("",4),
                  statpanel=NULL, add=T, rownames='Analysis temp, Z500, UV200', 
                  main=paste('EKF400 8-yr',s,'periods'), units='K',
                  addcontours=T, contvarname='gph500', conttype='data',
                  contcol='black', contlev=contlevs,
                  addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                  veclen=scfac*0.01, vecscale=scfac*2.0, vecwd=0.95, every_x_vec=4)
      # diff echam analysis
      levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
      contlevs <- seq(-4,4,2)
      if (ti==3){
        pdata$data <- array(c((rowMeans(anameananom[,1:2])-rowMeans(echmeananom[,1:2])),
                              (rowMeans(anameananom[,3:4])-rowMeans(echmeananom[,3:4])),
                              (rowMeans(anameananom[,5:6])-rowMeans(echmeananom[,5:6])),
                              (rowMeans(anameananom[,7:8])-rowMeans(echmeananom[,7:8]))),
                              c(nrow(echam.anom$ensmean),1,4))
      } else {
        pdata$data <- array(c((anameananom[,ti]-echmeananom[,ti]),
                            (anameananom[,ti+2]-echmeananom[,ti+2]),
                            (anameananom[,ti+4]-echmeananom[,ti+4]),
                            (anameananom[,ti+6]-echmeananom[,ti+6])),
                            c(nrow(echam.anom$ensmean),1,4))
      }
      plot_echam3(pdata, varname='temp2', type='data', cex.pt=1.5, 
                  latlim=c(0,90),lonlim=c(-180,180),
                  wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                  names=c('i)', 'j)', 'k)', 'l)'),colnames=rep(" ",4),
                  statpanel=NULL, add=T, rownames='Update (Anal.-ECHAM)', 
                  units='K',
                  addcontours=T, contvarname='gph500', conttype='data',
                  contcol='black', contlev=contlevs,
                  addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                  veclen=scfac*0.01, vecscale=scfac*2.0, vecwd=0.95, every_x_vec=4)
      if (ti==3){
        pdata$data <- array(c(rowMeans(anameananom[,1:2]),rowMeans(anameananom[,3:4]),
                              rowMeans(anameananom[,5:6]),rowMeans(anameananom[,7:8])),
                            c(nrow(echam.anom$ensmean),1,4))
      } else {
        pdata$data <- array(c(anameananom[,ti],anameananom[,ti+2],anameananom[,ti+4],
                              anameananom[,ti+6]), c(nrow(echam.anom$ensmean),1,4))
      }
      levs=c(0,80,85,90,95,100,105,110,115,120,Inf)
      plot_echam3(pdata, varname='precip', type='data', cex.pt=1.5, 
                  latlim=c(0,90),lonlim=c(-180,180),units="%",
                  names=c('m)', 'n)','o)', 'p)'), lev=levs, 
                  st.col=NULL, stations=NULL, add=T,addcontours=F,
                  wcol='darkgrey',rownames='Prec. analysis',colnames=rep('',4))
      dev.off()  
    } # end summer winter loop
    
    
    
    coldper <- c(1739, 1809, 1815, 1831) 
    for (syr in coldper) {
      eyr <- syr+7
      # change to precip in % 
      # load monthly 70yr climatology data
      if (syr < 1636) {
        load(file=paste0(echclimpath,'echam_clim_1636-1637_2ndgrid.Rdata'))
      } else {
        load(file=paste0(echclimpath,'echam_clim_',syr+3,'-',(syr+4),'_2ndgrid.Rdata'))  
      }
      pos <- which(echam_clim$names=="precip")[1:4608]
      # seasonal mean
      echmeanclim <- cbind(rowMeans(echam_clim$ensmean[,10:15]),
                           rowMeans(echam_clim$ensmean[,16:21]))
      for (cyr in seq(syr,eyr)) {
        load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
        # change to precip in % 
        echam.abs <- echam
        echam.anom$ensmean[pos,] <- echam.abs$ensmean[pos,]/(echmeanclim[pos,]*3600*24*30)*100
        analysis.abs <- analysis
        analysis.anom$ensmean[pos,] <- analysis.abs$ensmean[pos,]/
          (echmeanclim[pos,]*3600*24*30)*100
        if (cyr==syr){
          atmp <- analysis.anom$ensmean
          etmp <- echam.anom$ensmean
          c=1
        } else {
          atmp <- atmp+analysis.anom$ensmean
          etmp <- etmp+echam.anom$ensmean
          c=c+1
        }
      }
      # read mann et al temp recon
      mann <- read.table(paste0(datadir,"mann_temp_recon/allproxyfieldrecon"),na.string="NaN")
      ll <- read.table(paste0(datadir,"mann_temp_recon/longlat.txt"),na.string="NaN")
      ts <- which(mann[,1]==syr)
      plotmann <- list() #echam
      mclim <- as.vector(t(colMeans(mann[((ts-35):(ts+39)),-1],na.rm=F)))
      mdata <- as.vector(t(colMeans(mann[ts:(ts+4),-1],na.rm=F)))
      manom <- mdata-mclim
      plotmann$data <- array(manom,c(length(manom),1))
      plotmann$time <- mann[ts,1]
      plotmann$names <- rep('temp2',dim(plotmann$data)[1])
      plotmann$lon <- ll[,1] 
      plotmann$lat <- ll[,2]
      plotmann$height <- ll[,3]   
      plotmann$lsm.i <- ll[,3]  
      if (syr==coldper[1]){
        anameananom <- atmp/c
        echmeananom <- etmp/c
        mannmeananom <- plotmann$data
      } else {
        anameananom <- cbind(anameananom,(atmp/c))
        echmeananom <- cbind(echmeananom,(etmp/c))
        mannmeananom <- cbind(mannmeananom,plotmann$data)
      }
    } # end periods loop  
    colnames(anameananom) <- rep(coldper,each=2)
    colnames(echmeananom) <- rep(coldper,each=2)
    colnames(mannmeananom) <- coldper
    
    # 8-yr period winter
    pdata <- echam
    for (ti in 1:3){
      if (ti==1){s="winter"} else if (ti==2){s="summer"} else {s="annual"}
      pdf(paste0('../figures/',expname,'/1740-1840_cold_',s,filenameext,'.pdf'), width=16, 
          height=12, paper='special')
      layout(matrix(c(1,2,3,4,5,5,5,5,6,7,8,9,10,10,10,10,11,12,13,14,15,15,15,15,
                      16,17,18,19,20,20,20,20),8,4,byrow=TRUE),height=c(3,1,3,1,3,1,3,1))
      par(oma=c(0.5,4,4,0))
      levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
      plotmann$data <- array(c(mannmeananom[,1],mannmeananom[,2],mannmeananom[,3],
                               mannmeananom[,4]), c(nrow(plotmann$data),1,4))
      levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
      plot_echam3(plotmann, varname='temp2', type='data', cex.pt=2.5, 
                  latlim=c(0,90),lonlim=c(-180,180),
                  wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                  names=c('a)', 'b)', 'c)', 'd)'),colnames=coldper[1:4],
                  statpanel=NULL, add=T, rownames='Mann recon. temp.', 
                  main='', units='K')
      levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
      contlevs <- seq(-4,4,2)
      scfac <- max(echam.anom$ensmean[echam$names=='u200'])
      if (ti==3){
        pdata$data <- array(c(rowMeans(anameananom[,1:2]),rowMeans(anameananom[,3:4]),
                              rowMeans(anameananom[,5:6]),rowMeans(anameananom[,7:8])),
                            c(nrow(echam.anom$ensmean),1,4))
      } else {
        pdata$data <- array(c(anameananom[,ti],anameananom[,ti+2],anameananom[,ti+4],
                              anameananom[,ti+6]), c(nrow(echam.anom$ensmean),1,4))
      }
      plot_echam3(pdata, varname='temp2', type='data', cex.pt=1.5, 
                  latlim=c(0,90),lonlim=c(-180,180),
                  wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                  names=c('e)', 'f)', 'g)', 'h)'),colnames=rep("",4),
                  statpanel=NULL, add=T, rownames='Analysis temp, Z500, UV200', 
                  main=paste('EKF400 8-yr',s,'periods'), units='K',
                  addcontours=T, contvarname='gph500', conttype='data',
                  contcol='black', contlev=contlevs,
                  addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                  veclen=scfac*0.01, vecscale=scfac*2.0, vecwd=0.95, every_x_vec=4)
      # diff echam analysis
      levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
      contlevs <- seq(-4,4,2)
      if (ti==3){
        pdata$data <- array(c((rowMeans(anameananom[,1:2])-rowMeans(echmeananom[,1:2])),
                              (rowMeans(anameananom[,3:4])-rowMeans(echmeananom[,3:4])),
                              (rowMeans(anameananom[,5:6])-rowMeans(echmeananom[,5:6])),
                              (rowMeans(anameananom[,7:8])-rowMeans(echmeananom[,7:8]))),
                            c(nrow(echam.anom$ensmean),1,4))
      } else {
        pdata$data <- array(c((anameananom[,ti]-echmeananom[,ti]),
                              (anameananom[,ti+2]-echmeananom[,ti+2]),
                              (anameananom[,ti+4]-echmeananom[,ti+4]),
                              (anameananom[,ti+6]-echmeananom[,ti+6])),
                            c(nrow(echam.anom$ensmean),1,4))
      }
      plot_echam3(pdata, varname='temp2', type='data', cex.pt=1.5, 
                  latlim=c(0,90),lonlim=c(-180,180),
                  wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                  names=c('i)', 'j)', 'k)', 'l)'),colnames=rep(" ",4),
                  statpanel=NULL, add=T, rownames='Update (Anal.-ECHAM)', 
                  units='K',
                  addcontours=T, contvarname='gph500', conttype='data',
                  contcol='black', contlev=contlevs,
                  addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                  veclen=scfac*0.01, vecscale=scfac*2.0, vecwd=0.95, every_x_vec=4)
      if (ti==3){
        pdata$data <- array(c(rowMeans(anameananom[,1:2]),rowMeans(anameananom[,3:4]),
                              rowMeans(anameananom[,5:6]),rowMeans(anameananom[,7:8])),
                            c(nrow(echam.anom$ensmean),1,4))
      } else {
        pdata$data <- array(c(anameananom[,ti],anameananom[,ti+2],anameananom[,ti+4],
                              anameananom[,ti+6]), c(nrow(echam.anom$ensmean),1,4))
      }
      levs=c(0,80,85,90,95,100,105,110,115,120,Inf)
      plot_echam3(pdata, varname='precip', type='data', cex.pt=1.5, 
                  latlim=c(0,90),lonlim=c(-180,180),units="%",
                  names=c('m)', 'n)','o)', 'p)'), lev=levs, 
                  st.col=NULL, stations=NULL, add=T,addcontours=F,
                  wcol='darkgrey',rownames='Prec. analysis',colnames=rep('',4))
      dev.off()  
    } # end summer winter loop
} # end 1790th 1830th warm, cold periods NH
  
  
  
  
  
    
  
  
  


if (dryhumid_us) {
  load(paste0('../data/indices/',expname,'/indices',filenameext,'allts.Rdata'))
 # find and mark driest/wettest 8yr period in US
  d <- aind.allts$ensmean[(which(eind.allts$names=="NAM.precip")),
                          (is.even(seq(1,ncol(aind.allts$ensmean))))]
  dec <- rep(NA,250)
  dectime <- rep(NA,250)
  for (ti in 1:length(dec)) {
    if (ti==1) {
      dec[ti] <- mean(d[1:8])
      dectime[ti] <- 1603
      i=2
    } else {
      dec[ti] <- mean(d[i:(i+7)])
      dectime[ti] <- dectime[ti-1] + 1
      i=i+1
    }
  }
  #  plot(dectime,dec,ty='l',col='red')
  tmp=cbind(dec,dectime)
  tmp[order(dec),]


  # dryper <- c(1819, 1843, 1635, 1671, 1760, 1743)
  dryper <- c(1816, 1743, 1607,1760,1641)
  for (syr in dryper) {
    eyr <- syr+7
    # load nada pdsi recon
    nada <- read_pdsi('nada_pdsi_echamgrid.nc',timlim=c(syr,eyr), 
                      path="/Users/joerg/Documents/climdata/recon/cook/pdsi_usa/",
                      small=T,landonly=F)
    nada$data[is.na(nada$data)] <- 0
    # change to precip in % 
    # load monthly 70yr climatology data
    if (syr < 1636) {
      load(file=paste0(echclimpath,'echam_clim_1636-1637_2ndgrid.Rdata'))
    } else {
      load(file=paste0(echclimpath,'echam_clim_',syr+3,'-',(syr+4),'_2ndgrid.Rdata'))  
    }
    pos <- which(echam$names=="precip")[1:4608]
    # seasonal mean
    echmeanclim <- cbind(rowMeans(echam_clim$ensmean[,10:15]),
                         rowMeans(echam_clim$ensmean[,16:21]))
    #      echam$ensmean <- echmeanclim
    for (cyr in seq(syr,eyr)) {
      load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
      # change to precip in % 
      echam.abs <- echam
      echam.anom$ensmean[pos,] <- echam.abs$ensmean[pos,]/(echmeanclim[pos,]*3600*24*30)*100
      analysis.abs <- analysis
      analysis.anom$ensmean[pos,] <- analysis.abs$ensmean[pos,]/
                                        (echmeanclim[pos,]*3600*24*30)*100
      if (cyr==syr){
        atmp <- analysis.anom$ensmean
        etmp <- echam.anom$ensmean
        c=1
      } else {
        atmp <- atmp+analysis.anom$ensmean
        etmp <- etmp+echam.anom$ensmean
        c=c+1
      }
    }
    if (syr==dryper[1]){
      anameananom <- atmp/c
      echmeananom <- etmp/c
      pdsimeananom <- apply(nada$data,1,mean)
    } else {
      anameananom <- cbind(anameananom,(atmp/c))
      echmeananom <- cbind(echmeananom,(etmp/c))
      pdsimeananom <- cbind(pdsimeananom,(apply(nada$data,1,mean)))
    }
  } # end dry periods loop  
  colnames(anameananom) <- rep(dryper,each=2)
  colnames(echmeananom) <- rep(dryper,each=2)

  pdata <- echam
  for (ti in 1:2){
    if (ti==1){s="win"} else {s="som"}
    pdf(paste0('../figures/',expname,'/dry_us_',s,filenameext,'.pdf'), width=12, 
        height=14, paper='special')
    #      layout(matrix(c(1,2,3,4,5,6,7,8), 4, 2, byrow = F), height=c(3,1,3,1))
    layout(matrix(c(1,2,3,4,5,5,5,5,6,7,8,9,10,10,10,10,11,12,13,14,15,15,15,15,
                    16,17,18,19,20,20,20,20),8,4,byrow=TRUE),height=c(3,1,3,1,3,1,3,1))
    par(oma=c(0.5,4,4,0))
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    scfac <- max(echam.anom$ensmean[echam$names=='u200'])
    pdata$data <- array(c(anameananom[,ti],anameananom[,ti+2],anameananom[,ti+4],
                          anameananom[,ti+6]), c(nrow(echam.anom$ensmean),1,4))
    plot_echam3(pdata, varname='temp2', type='data', cex.pt=3.7, 
                latlim=c(20,60),lonlim=c(-160,-40),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('a)', 'b)', 'c)', 'd)'),colnames=dryper[1:4],
                statpanel=NULL, add=T, rownames='Analysis T2m, Z500, UV200', 
                main='EKF400 5-yr dry periods in North America', units='K',
                addcontours=T, contvarname='gph500', conttype='data',
                contcol='black', contlev=contlevs,
                addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                veclen=scfac*0.02, vecscale=scfac*0.98, vecwd=0.95, every_x_vec=2)
    # diff echam analysis
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    pdata$data <- array(c((anameananom[,ti]-echmeananom[,ti]),
                          (anameananom[,ti+2]-echmeananom[,ti+2]),
                          (anameananom[,ti+4]-echmeananom[,ti+4]),
                          (anameananom[,ti+6]-echmeananom[,ti+6])),
                        c(nrow(echam.anom$ensmean),1,4))
    plot_echam3(pdata, varname='temp2', type='data', cex.pt=3.7, 
                latlim=c(20,60),lonlim=c(-160,-40),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('e)', 'f)', 'g)', 'h)'),colnames=rep(" ",4),
                statpanel=NULL, add=T, rownames='Update (Anal.-ECHAM)', 
                units='K',
                addcontours=T, contvarname='gph500', conttype='data',
                contcol='black', contlev=contlevs,
                addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                veclen=scfac*0.02, vecscale=scfac*0.98, vecwd=0.95, every_x_vec=2)
    pdata$data <- array(c(anameananom[,ti],anameananom[,ti+2],anameananom[,ti+4],
                          anameananom[,ti+6]), c(nrow(echam.anom$ensmean),1,4))
    #      #    levs <- c(-Inf, seq(-5,5,1), Inf)
    levs=c(50,80,85,90,95,100,105,110,115,120,200)
    plot_echam3(pdata, varname='precip', type='data', cex.pt=3.7, 
                latlim=c(20,60),lonlim=c(-160,-40),units="%",
                names=c('i)', 'j)','k)', 'l)'), lev=levs, 
                st.col=NULL, stations=NULL, add=T,addcontours=F,
                wcol='darkgrey',rownames='Precipitation',colnames=rep('',4))
    pdata$data <- array(c(pdsimeananom[,1],pdsimeananom[,2],pdsimeananom[,3],
                          pdsimeananom[,4]), c(nrow(echam.anom$ensmean),1,4))
    levs <- seq(-4,4,1)
    plot_echam3(pdata, varname='precip', type='data', cex.pt=3.7, 
                latlim=c(20,60),lonlim=c(-160,-40),units="PDSI index",
                names=c('m)', 'n)','o)', 'p)'), lev=levs, 
                st.col=NULL, stations=NULL, add=T,addcontours=F,
                wcol='darkgrey',rownames='PDSI',colnames=rep('',4))
    dev.off()  
  } # end summer winter loop

  
  #humidper <- c(1688, 1729, 1656, 1802, 1698) #, 1603)
  humidper <- c(1655, 1725, 1777, 1799, 1686)
  for (syr in humidper) {
    eyr <- syr+7
    nada <- read_pdsi('nada_pdsi_echamgrid.nc',timlim=c(syr,eyr), 
                      path="/Users/joerg/Documents/climdata/recon/cook/pdsi_usa/",
                      small=T,landonly=F)
    nada$data[is.na(nada$data)] <- 0
    # change to precip in % 
    # load monthly 70yr climatology data
    if (syr < 1636) {
      load(file=paste0(echclimpath,'echam_clim_1636-1637_2ndgrid.Rdata'))
    } else {
      load(file=paste0(echclimpath,'echam_clim_',syr+3,'-',(syr+4),'_2ndgrid.Rdata'))  
    }
    pos <- which(echam$names=="precip")[1:4608]
    # seasonal mean
    echmeanclim <- cbind(rowMeans(echam_clim$ensmean[,10:15]),
                         rowMeans(echam_clim$ensmean[,16:21]))
    for (cyr in seq(syr,eyr)) {
      load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
      # change to precip in % 
      echam.abs <- echam
      echam.anom$ensmean[pos,] <- echam.abs$ensmean[pos,]/(echmeanclim[pos,]*3600*24*30)*100
      analysis.abs <- analysis
      analysis.anom$ensmean[pos,] <- analysis.abs$ensmean[pos,]/
                                        (echmeanclim[pos,]*3600*24*30)*100
      if (cyr==syr){
        atmp <- analysis.anom$ensmean
        etmp <- echam.anom$ensmean
        c=1
      } else {
        atmp <- atmp+analysis.anom$ensmean
        etmp <- etmp+echam.anom$ensmean
        c=c+1
      }
    }
    if (syr==humidper[1]){
      anameananom <- atmp/c
      echmeananom <- etmp/c
      pdsimeananom <- apply(nada$data,1,mean)
    } else {
      anameananom <- cbind(anameananom,(atmp/c))
      echmeananom <- cbind(echmeananom,(etmp/c))
      pdsimeananom <- cbind(pdsimeananom,(apply(nada$data,1,mean)))
    }
  } # end humid periods loop  
  colnames(anameananom) <- rep(humidper,each=2)
  colnames(echmeananom) <- rep(humidper,each=2)
  
  pdata <- echam
  for (ti in 1:2){
    if (ti==1){s="win"} else {s="som"}
    pdf(paste0('../figures/',expname,'/humid_us_',s,filenameext,'.pdf'), width=12, 
        height=10.5, paper='special')
    #      layout(matrix(c(1,2,3,4,5,6,7,8), 4, 2, byrow = F), height=c(3,1,3,1))
    layout(matrix(c(1,2,3,4,5,5,5,5,6,7,8,9,10,10,10,10,11,12,13,14,15,15,15,15,
                    16,17,18,19,20,20,20,20),8,4,byrow=TRUE),height=c(3,1,3,1,3,1,3,1))
    par(oma=c(0.5,4,4,0))
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    scfac <- max(echam.anom$ensmean[echam$names=='u200'])
    pdata$data <- array(c(anameananom[,ti],anameananom[,ti+2],anameananom[,ti+4],
                          anameananom[,ti+6]), c(nrow(echam.anom$ensmean),1,4))
    plot_echam3(pdata, varname='temp2', type='data', cex.pt=3.7, 
                latlim=c(20,60),lonlim=c(-160,-40),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('a)', 'b)', 'c)', 'd)'),colnames=humidper[1:4],
                statpanel=NULL, add=T, rownames='Analysis T2m, Z500, UV200', 
                main='EKF400 5-yr humid periods in North America', units='K',
                addcontours=T, contvarname='gph500', conttype='data',
                contcol='black', contlev=contlevs,
                addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                veclen=scfac*0.02, vecscale=scfac*0.98, vecwd=0.95, every_x_vec=2)
    # diff echam analysis
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    pdata$data <- array(c((anameananom[,ti]-echmeananom[,ti]),
                          (anameananom[,ti+2]-echmeananom[,ti+2]),
                          (anameananom[,ti+4]-echmeananom[,ti+4]),
                          (anameananom[,ti+6]-echmeananom[,ti+6])),
                        c(nrow(echam.anom$ensmean),1,4))
    plot_echam3(pdata, varname='temp2', type='data', cex.pt=3.7, 
                latlim=c(20,60),lonlim=c(-160,-40),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('e)', 'f)', 'g)', 'h)'),colnames=rep(" ",4),
                statpanel=NULL, add=T, rownames='Update (Anal.-ECHAM)', 
                units='K',
                addcontours=T, contvarname='gph500', conttype='data',
                contcol='black', contlev=contlevs,
                addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                veclen=scfac*0.02, vecscale=scfac*0.98, vecwd=0.95, every_x_vec=2)
    pdata$data <- array(c(anameananom[,ti],anameananom[,ti+2],anameananom[,ti+4],
                          anameananom[,ti+6]), c(nrow(echam.anom$ensmean),1,4))
    levs=c(50,80,85,90,95,100,105,110,115,120,200)
    plot_echam3(pdata, varname='precip', type='data', cex.pt=3.7, 
                latlim=c(20,60),lonlim=c(-160,-40),units="%",
                names=c('e)', 'f)','g)', 'h)'), lev=levs, 
                st.col=NULL, stations=NULL, add=T,addcontours=F,
                wcol='darkgrey',rownames='Precipitation',colnames=rep('',4))
    pdata$data <- array(c(pdsimeananom[,1],pdsimeananom[,2],pdsimeananom[,3],
                          pdsimeananom[,4]), c(nrow(echam.anom$ensmean),1,4))
    levs <- seq(-4,4,1)
    plot_echam3(pdata, varname='precip', type='data', cex.pt=3.7, 
                latlim=c(20,60),lonlim=c(-160,-40),units="PDSI index",
                names=c('m)', 'n)','o)', 'p)'), lev=levs, 
                st.col=NULL, stations=NULL, add=T,addcontours=F,
                wcol='darkgrey',rownames='PDSI',colnames=rep('',4))
    dev.off()  
  } # end summer winter loop  
} # end dry, humid periods US










if (dryhumid_eu) {
  load(paste0('../data/indices/',expname,'/indices',filenameext,'allts.Rdata'))
 # find and mark driest/wettest 8yr period in EU
  d <- aind.allts$ensmean[(which(eind.allts$names=="EU.precip")),
                          (is.even(seq(1,ncol(aind.allts$ensmean))))]
  dec <- rep(NA,250)
  dectime <- rep(NA,250)
  for (ti in 1:length(dec)) {
    if (ti==1) {
      dec[ti] <- mean(d[1:8])
      dectime[ti] <- 1603
      i=2
    } else {
      dec[ti] <- mean(d[i:(i+7)])
      dectime[ti] <- dectime[ti-1] + 1
      i=i+1
    }
  }
  #  plot(dectime,dec,ty='l',col='red')
  tmp=cbind(dec,dectime)
  tmp[order(dec),]

  
  #dryper <- c(1846, 1823, 1808, 1837, 1816, 1794)
  dryper <- c(1843, 1807, 1743, 1822, 1630)
  for (syr in dryper) {
    #    # calc validate clim (70 yr period) to substract for anomalies
    #     if (syr>1645) {asyr<-syr-35} else {asyr=1610}
    #     for (cyr in seq(asyr,(syr+39))) {  
    #       load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
    #       if (cyr==asyr){
    #         vtmp <- validate$data
    #         c=1
    #       } else {
    #         vtmp <- vtmp+validate$data
    #         c <- c+1
    #       }
    #     }
    #     valclim <- vtmp/c  
    eyr <- syr+7
    # change to precip in % 
    # load monthly 70yr climatology data
    if (syr < 1636) {
      load(file=paste0(echclimpath,'echam_clim_1636-1637_2ndgrid.Rdata'))
    } else {
      load(file=paste0(echclimpath,'echam_clim_',syr+3,'-',(syr+4),'_2ndgrid.Rdata'))  
    }
    pos <- which(echam_clim$names=="precip")[1:4608]
    # seasonal mean
    echmeanclim <- cbind(rowMeans(echam_clim$ensmean[,10:15]),
                         rowMeans(echam_clim$ensmean[,16:21]))
    #      echam$ensmean <- echmeanclim
    for (cyr in seq(syr,eyr)) {
      load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
      # change to precip in % 
      echam.abs <- echam
      echam.anom$ensmean[pos,] <- echam.abs$ensmean[pos,]/(echmeanclim[pos,]*3600*24*30)*100
      analysis.abs <- analysis
      analysis.anom$ensmean[pos,] <- analysis.abs$ensmean[pos,]/
                                        (echmeanclim[pos,]*3600*24*30)*100
      if (cyr==syr){
        atmp <- analysis.anom$ensmean
        etmp <- echam.anom$ensmean
        #        vtmp <- validate$data
        c=1
      } else {
        atmp <- atmp+analysis.anom$ensmean
        etmp <- etmp+echam.anom$ensmean
        #        vtmp ...        
        c=c+1
      }
    }
    if (syr==dryper[1]){
      anameananom <- atmp/c
      echmeananom <- etmp/c
      #      valmean <- vtmp/c
      #      valmeananom <- valmean-valclim
    } else {
      anameananom <- cbind(anameananom,(atmp/c))
      echmeananom <- cbind(echmeananom,(etmp/c))
      #      valmean <- vtmp/c
      #      valmeananom <- cbind(valmeananom,(valmean-valclim))
    }
  } # end dry periods loop  
  colnames(anameananom) <- rep(dryper,each=2)
  colnames(echmeananom) <- rep(dryper,each=2)
  
  pdata <- echam
  for (ti in 1:2){
    if (ti==1){s="win"} else {s="som"}
    #   ti=2 # summer
    pdf(paste0('../figures/',expname,'/dry_eu_',s,filenameext,'.pdf'), width=12, 
        height=10.5, paper='special')
    #      layout(matrix(c(1,2,3,4,5,6,7,8), 4, 2, byrow = F), height=c(3,1,3,1))
    layout(matrix(c(1,2,3,4,5,5,5,5,6,7,8,9,10,10,10,10,11,12,13,14,15,15,15,15),
                  6,4,byrow=TRUE),height=c(3,1,3,1,3,1))
    par(oma=c(0.5,4,4,0))
    #     pdata$data <- array(c(valmeananom[,1],valmeananom[,2],valmeananom[,3],
    #                           valmeananom[,4]), c(nrow(echam.anom$ensmean),1,4))
    #     levs <- seq(-4,4,1)
    #     plot_echam3(pdata, varname='temp2', type='data', cex.pt=3.7, 
    #                 latlim=c(30,75),lonlim=c(-20,40),units="K",
    #                 names=c('m)', 'n)','o)', 'p)'), lev=levs, 
    #                 st.col=NULL, stations=NULL, add=T,addcontours=F,
    #                 wcol='darkgrey',rownames='Temp. reconstruction',colnames=rep('',4))
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    scfac <- max(echam.anom$ensmean[echam$names=='u200'])
    pdata$data <- array(c(anameananom[,ti],anameananom[,ti+2],anameananom[,ti+4],
                          anameananom[,ti+6]), c(nrow(echam.anom$ensmean),1,4))
    plot_echam3(pdata, varname='temp2', type='data', cex.pt=3.7, 
                latlim=c(30,75),lonlim=c(-20,40),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('a)', 'b)', 'c)', 'd)'),colnames=dryper[1:4],
                statpanel=NULL, add=T, rownames='Analysis T2m, Z500, UV200', 
                main='EKF400 5-yr dry periods in Europe', units='K',
                addcontours=T, contvarname='gph500', conttype='data',
                contcol='black', contlev=contlevs,
                addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                veclen=scfac*0.02, vecscale=scfac*0.98, vecwd=0.95, every_x_vec=2)
    # diff echam analysis
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    pdata$data <- array(c((anameananom[,ti]-echmeananom[,ti]),
                          (anameananom[,ti+2]-echmeananom[,ti+2]),
                          (anameananom[,ti+4]-echmeananom[,ti+4]),
                          (anameananom[,ti+6]-echmeananom[,ti+6])),
                        c(nrow(echam.anom$ensmean),1,4))
    plot_echam3(pdata, varname='temp2', type='data', cex.pt=3.7, 
                latlim=c(30,75),lonlim=c(-20,40),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('e)', 'f)', 'g)', 'h)'),colnames=rep(" ",4),
                statpanel=NULL, add=T, rownames='Update (Anal.-ECHAM)', 
                units='K',
                addcontours=T, contvarname='gph500', conttype='data',
                contcol='black', contlev=contlevs,
                addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                veclen=scfac*0.02, vecscale=scfac*0.98, vecwd=0.95, every_x_vec=2)
    pdata$data <- array(c(anameananom[,ti],anameananom[,ti+2],anameananom[,ti+4],
                          anameananom[,ti+6]), c(nrow(echam.anom$ensmean),1,4))
    #      #    levs <- c(-Inf, seq(-5,5,1), Inf)
    levs=c(50,80,85,90,95,100,105,110,115,120,200)
    plot_echam3(pdata, varname='precip', type='data', cex.pt=3.7, 
                latlim=c(30,75),lonlim=c(-20,40),units="%",
                names=c('i)', 'j)','k)', 'l)'), lev=levs, 
                st.col=NULL, stations=NULL, add=T,addcontours=F,
                wcol='darkgrey',rownames='Prec. analysis',colnames=rep('',4))
    #     pdata$data <- array(c(valmeananom[,1],valmeananom[,2],valmeananom[,3],
    #                           valmeananom[,4]), c(nrow(echam.anom$ensmean),1,4))
    #     levs <- seq(-4,4,1)
    #     plot_echam3(pdata, varname='precip', type='data', cex.pt=3.7, 
    #                 latlim=c(30,75),lonlim=c(-20,40),units="mm",
    #                 names=c('m)', 'n)','o)', 'p)'), lev=levs, 
    #                 st.col=NULL, stations=NULL, add=T,addcontours=F,
    #                 wcol='darkgrey',rownames='Prec. reconstruction',colnames=rep('',4))
    dev.off()  
  } # end summer winter loop
  
  
  
  
  
  #humidper <- c(1617, 1664, 1638, 1701, 1713, 1625)
  humidper <- c(1661, 1735, 1688, 1767, 1708, 1826)
  for (syr in humidper) {
    eyr <- syr+7
    # change to precip in % 
    # load monthly 70yr climatology data
    if (syr < 1636) {
      load(file=paste0(echclimpath,'echam_clim_1636-1637_2ndgrid.Rdata'))
    } else {
      load(file=paste0(echclimpath,'echam_clim_',syr+3,'-',(syr+4),'_2ndgrid.Rdata'))  
    }
    pos <- which(echam_clim$names=="precip")[1:4608]
    # seasonal mean
    echmeanclim <- cbind(rowMeans(echam_clim$ensmean[,10:15]),
                         rowMeans(echam_clim$ensmean[,16:21]))
    for (cyr in seq(syr,eyr)) {
      load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
      # change to precip in % 
      echam.abs <- echam
      echam.anom$ensmean[pos,] <- echam.abs$ensmean[pos,]/(echmeanclim[pos,]*3600*24*30)*100
      analysis.abs <- analysis
      analysis.anom$ensmean[pos,] <- analysis.abs$ensmean[pos,]/
                                        (echmeanclim[pos,]*3600*24*30)*100
      if (cyr==syr){
        atmp <- analysis.anom$ensmean
        etmp <- echam.anom$ensmean
        c=1
      } else {
        atmp <- atmp+analysis.anom$ensmean
        etmp <- etmp+echam.anom$ensmean
        c=c+1
      }
    }
    if (syr==humidper[1]){
      anameananom <- atmp/c
      echmeananom <- etmp/c
    } else {
      anameananom <- cbind(anameananom,(atmp/c))
      echmeananom <- cbind(echmeananom,(etmp/c))
    }
  } # end humid periods loop  
  colnames(anameananom) <- rep(humidper,each=2)
  colnames(echmeananom) <- rep(humidper,each=2)
  
    pdata <- echam
  for (ti in 1:2){
    if (ti==1){s="win"} else {s="som"}
    pdf(paste0('../figures/',expname,'/humid_eu_',s,filenameext,'.pdf'), width=12, 
        height=10.5, paper='special')
    layout(matrix(c(1,2,3,4,5,5,5,5,6,7,8,9,10,10,10,10,11,12,13,14,15,15,15,15),
                  6,4,byrow=TRUE),height=c(3,1,3,1,3,1))
    par(oma=c(0.5,4,4,0))
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    scfac <- max(echam.anom$ensmean[echam$names=='u200'])
    pdata$data <- array(c(anameananom[,ti],anameananom[,ti+2],anameananom[,ti+4],
                          anameananom[,ti+6]), c(nrow(echam.anom$ensmean),1,4))
    plot_echam3(pdata, varname='temp2', type='data', cex.pt=3.7, 
                latlim=c(30,75),lonlim=c(-20,40),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('a)', 'b)', 'c)', 'd)'),colnames=humidper[1:4],
                statpanel=NULL, add=T, rownames='Analysis T2m, Z500, UV200', 
                main='EKF400 5-yr humid periods in Europe', units='K',
                addcontours=T, contvarname='gph500', conttype='data',
                contcol='black', contlev=contlevs,
                addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                veclen=scfac*0.02, vecscale=scfac*0.98, vecwd=0.95, every_x_vec=2)
    # diff echam analysis
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    pdata$data <- array(c((anameananom[,ti]-echmeananom[,ti]),
                          (anameananom[,ti+2]-echmeananom[,ti+2]),
                          (anameananom[,ti+4]-echmeananom[,ti+4]),
                          (anameananom[,ti+6]-echmeananom[,ti+6])),
                        c(nrow(echam.anom$ensmean),1,4))
    plot_echam3(pdata, varname='temp2', type='data', cex.pt=3.7, 
                latlim=c(30,75),lonlim=c(-20,40),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('e)', 'f)', 'g)', 'h)'),colnames=rep(" ",4),
                statpanel=NULL, add=T, rownames='Update (Anal.-ECHAM)', 
                units='K',
                addcontours=T, contvarname='gph500', conttype='data',
                contcol='black', contlev=contlevs,
                addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                veclen=scfac*0.02, vecscale=scfac*0.98, vecwd=0.95, every_x_vec=2)
    pdata$data <- array(c(anameananom[,ti],anameananom[,ti+2],anameananom[,ti+4],
                          anameananom[,ti+6]), c(nrow(echam.anom$ensmean),1,4))
    levs=c(50,80,85,90,95,100,105,110,115,120,200)
    plot_echam3(pdata, varname='precip', type='data', cex.pt=3.7, 
                latlim=c(30,75),lonlim=c(-20,40),units="%",
                names=c('i)', 'j)','k)', 'l)'), lev=levs, 
                st.col=NULL, stations=NULL, add=T,addcontours=F,
                wcol='darkgrey',rownames='Prec. analysis',colnames=rep('',4))
    dev.off()  
  } # end summer winter loop
} # end dry, humid periods EU




  
  
  
  

if (jetshift){
  load(paste0('../data/indices/',expname,'/indices',filenameext,'allts.Rdata'))
  # Stefans ideas: 
  # decadal shifts of subtropical belt. 1945 northern position, 
  # until 1980 southward movement, since 1980 northward trend
  wl=11 # 31yr low-pass filter
  ts=1:402
  anomper <- (1901:1980)
  seas1 <- 'annual mean'
  seas2 <- 'winter'
  seas3 <- 'summer'
  reg1 <- 37 # for SJ_u200.calc
  reg2 <- 38 # for SJ_slp.calc
  reg3 <- 46 # for SJ index
  period <- eind.allts$time #seq(syr,eyr)
  pdf(file=paste0('../figures/',expname,'/indices/N_hadley_cell_extend',filenameext,wl,'yr_runmean.pdf'),
                    width=10,height=14)
  par(oma=c(1.5,1,1,0))
  par(mfrow=c(3,1))  
#  stefind <- read.table("/Users/joerg/Documents/climdata/indices/stefan/stefan_seasonal_indices.txt",
#               header=T,skip=1,na.string=-9999)[,c(1,seq(2,5),seq(18,21),seq(34,37),seq(50,53),
#               seq(66,69),seq(82,85),seq(98,101),seq(114,117),seq(130,133),seq(146,149))]
#  head <- read.table("/Users/joerg/Documents/climdata/indices/stefan/stefan_seasonal_indices.txt",
#            nrow=2)[,c(1,seq(2,5),seq(18,21),seq(34,37),seq(50,53),seq(66,69),seq(82,85),
#            seq(98,101),seq(114,117),seq(130,133),seq(146,149))]
#  # (1,2,3,18,19,34,35,50,51,66,67,82,83,98,99,114,115,130,131,146,147)]
#  newhead <- rep(NA,ncol(head))
#  for (i in 1:ncol(head)) {
#    newhead[i] <- paste(head[1,i],head[2,i],sep='_')
#  }
#  colnames(stefind) <- newhead
#  #load('../data/indices_recon.Rdata')
#  vf <- c(1600,1650,1660,1783,1809,1815,1835,1883,1886,1902,1912,1980,1991) 
#  #read.table('/Users/joerg/Documents/climdata/forcings/forc-total-crowley_2001.txt',
#  #  skip=51,header=T)[,1:2]
  
  for (seas in c(seas1,seas2,seas3)) {
    for (reg in c(reg1,reg2)) { #,reg3)) {
      if (seas=='annual mean') {
        years <- rep(eind.allts$time,each=2)  
      } else {
        wincol <- is.odd(seq(1,ncol(aind.allts$ensmean)))
        somcol <- is.even(seq(1,ncol(aind.allts$ensmean)))  
        wincolv <- is.odd(seq(1,ncol(vind.allts$data)))
        somcolv <- is.even(seq(1,ncol(vind.allts$data)))    
      }
      if (seas=='annual mean') {
        print("annual mean")
        eindmean <- aggregate(apply(eind.allts$data[reg,,],1,mean),list(years),mean,na.rm=T)[,2]
        aindmean <- aggregate(apply(aind.allts$data[reg,,],1,mean),list(years),mean,na.rm=T)[,2]
#        eindmean <- aggregate(eind.allts$ensmean[reg,],list(years),mean,na.rm=T)[,2]
#        aindmean <- aggregate(aind.allts$ensmean[reg,],list(years),mean,na.rm=T)[,2]
        eindbem <- aggregate(eind.allts$bem[reg,],list(years),mean,na.rm=T)[,2]
        aindbem <- aggregate(aind.allts$bem[reg,],list(years),mean,na.rm=T)[,2]
        eindmin <- apply(aggregate(eind.allts$data[reg,,],list(years),mean,na.rm=T)[2:31],1,min,na.rm=T)
        eindmax <- apply(aggregate(eind.allts$data[reg,,],list(years),mean,na.rm=T)[2:31],1,max,na.rm=T)
        aindmin <- apply(aggregate(aind.allts$data[reg,,],list(years),mean,na.rm=T)[2:31],1,min,na.rm=T)
        aindmax <- apply(aggregate(aind.allts$data[reg,,],list(years),mean,na.rm=T)[2:31],1,max,na.rm=T)
      } else if (seas=='summer'){  
        print("summer season")
        eindmean <- apply(eind.allts$data[reg,somcol,],1,mean)  
        aindmean <- apply(aind.allts$data[reg,somcol,],1,mean)  
#        eindmean <- eind.allts$ensmean[reg,somcol]  
#        aindmean <- aind.allts$ensmean[reg,somcol]  
        eindbem <- eind.allts$bem[reg,somcol]
        aindbem <- aind.allts$bem[reg,somcol]
        eindmin <- apply(eind.allts$data[reg,somcol,],1,min,na.rm=T)
        eindmax <- apply(eind.allts$data[reg,somcol,],1,max,na.rm=T)
        aindmin <- apply(aind.allts$data[reg,somcol,],1,min,na.rm=T)
        aindmax <- apply(aind.allts$data[reg,somcol,],1,max,na.rm=T)
      } else if (seas=='winter'){
        print("winter season")
        eindmean <- apply(eind.allts$data[reg,wincol,],1,mean)  
        aindmean <- apply(aind.allts$data[reg,wincol,],1,mean)  
#        eindmean <- eind.allts$ensmean[reg,wincol]  
#        aindmean <- aind.allts$ensmean[reg,wincol]  
        eindbem <- eind.allts$bem[reg,wincol]
        aindbem <- aind.allts$bem[reg,wincol]
        eindmin <- apply(eind.allts$data[reg,wincol,],1,min,na.rm=T)
        eindmax <- apply(eind.allts$data[reg,wincol,],1,max,na.rm=T)
        aindmin <- apply(aind.allts$data[reg,wincol,],1,min,na.rm=T)
        aindmax <- apply(aind.allts$data[reg,wincol,],1,max,na.rm=T)
      }
      print("join indices to allind")
      allind <- cbind(eindmean,aindmean,eindbem,aindbem,eindmin,eindmax,aindmin,aindmax)
      anopos <- period %in% anomper
      centerfac <- apply(allind[anopos,],2,mean)
      centerfac[5] <- centerfac[6] <- centerfac[3] <- centerfac[1]
      centerfac[7] <- centerfac[8] <- centerfac[4] <- centerfac[2]
      #  if (sca) {
      #    scalefac <- apply(allind[anopos,],2,sd)
      #    scalefac[5] <- scalefac[6] <- scalefac[3] <- scalefac[1]
      #    scalefac[7] <- scalefac[8] <- scalefac[4] <- scalefac[2]
      #  } else { scalefac <- rep(1,8) }  
      allindann <- allind
      if (wl > 1) {
        for (i in 1:ncol(allind)){
          allind[,i] <- runmean(allind[,i],wl) #scalefac),wl)
        }
        #  allind <- runmean(allind,wl) #scalefac),wl)
      }
      #print("scale allind")
      #if (wl > 1) {
      #  allind <- runmean(scale(allind,center=centerfac,scale=F),wl) #scalefac),wl)
      #} else {
      #  allind <- scale(allind,center=centerfac,scale=F) #scalefac)
      #}
      if (reg==reg1){allind1 <- allind}
      if (reg==reg2){allind2 <- allind}
      if (reg==reg3){allind3 <- allind}
    } # end reg loop
    # plot all
    plot(period[ts],allind1[ts,1],ty='l',lwd=1,lty=1,col=rgb(0,0,0,10,maxColorValue=10), 
           ylab="lat",xlab="",axes=F,ylim=c(min(allind1[ts,1])-0.5,max(allind1[ts,1])+1.5),
           main=paste0('Shifts in the northern Hadley cell extend (',seas,')'))
    axis(1);axis(2)
    lines(period[ts],allind1[ts,2],ty='l',lwd=1,lty=1,col=rgb(10,0,0,10,maxColorValue=10),
            ylab="",xlab="",axes=F)
    par(new=T)
    plot(period[ts],allind2[ts,1],ty='l',lwd=1,lty=2,col=rgb(0,0,10,10,maxColorValue=10),
           ylab="",xlab="",axes=F,ylim=c(min(allind2[ts,1])-1.5,max(allind2[ts,1])+0.5))
    axis(4,col='blue',col.axis='blue')
    lines(period[ts],allind2[ts,2],ty='l',lwd=1,lty=2,col=rgb(10,0,0,10,maxColorValue=10),
          ylab="",xlab="",axes=F)       
    if (seas==seas1){
      legend('topleft', c('echam zonal mean u200 max','anal. u200 max','echam slp max','anal. slp max'),
             bg='white',box.col='white',lty=c(1,1,2,2),col=c(1,2,4,2))
    }
#     vdata <- cbind(stefind['SJ_20C'],stefind['SJ_REC'],stefind['SJ_NCEP'],stefind['SJ_ERA-40'])
#     yrv <- stefind[,1]
#     centerfac <- apply(vdata,2,mean,na.rm=T)
#     scalefac <- apply(vdata,2,sd,na.rm=T)
#     #if (wlen > 1) {
#     vdata <- runmean(scale(vdata,center=F,scale=F),wl) #centerfac,scale=scalefac),wlen)
#     #} else {
#     #  vdata <- scale(vdata,center=centerfac,scale=scalefac)
#     #}
#     lines(yrv,vdata[,1],ty='l',lwd=linew,lty=1,col=rgb(0,0,10,10,maxColorValue=10), ylab="",xlab="")
#     lines(yrv,vdata[,2],ty='l',lwd=linew,lty=1,col=rgb(10,0,10,10,maxColorValue=10), ylab="",xlab="")
#     lines(yrv,vdata[,3],ty='l',lwd=linew,lty=1,col=rgb(10,10,0,10,maxColorValue=10), ylab="",xlab="")
#     lines(yrv,vdata[,4],ty='l',lwd=linew,lty=1,col=rgb(0,10,10,10,maxColorValue=10), ylab="",xlab="")
    
    # ecormean <- round(cor(data[['allind']][seq(298,402),1],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
    # acormean <- round(cor(data[['allind']][seq(298,402),2],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
    # ecorbem <- round(cor(data[['allind']][seq(298,402),3],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2) 
    # acorbem <- round(cor(data[['allind']][seq(298,402),4],vdata[seq(1,105),],use="pairwise.complete.obs"),digits=2)
    # legend('bottomleft', c(paste('ECHAM ens. mean - 20CR, REC, NCEP, ERA40:',ecormean[1],ecormean[2],ecormean[3],
    #   ecormean[4]),paste('Analysis ens. mean - 20CR, REC, NCEP, ERA40:',acormean[1],acormean[2],acormean[3],
    #   acormean[4]),paste('ECHAM ens. bem - 20CR, REC, NCEP, ERA40:',ecorbem[1],ecorbem[2],ecorbem[3],
    #   ecorbem[4]),paste('Analysis ens. bem - 20CR, REC, NCEP, ERA40:',acorbem[1],acorbem[2],acorbem[3],
    #   acorbem[4])),cex=fs, bty='o', bg='white', box.col='white')
  }
  dev.off()
} # end jet shift









if (amo) {
# plot AMO index in warm/cold periods
#load(paste0(prepplotdir,'.Rdata'))
if (ind_anom) {
  load(paste0(prepplotdir,'indices_anom_allts_landonly.Rdata'))
} else {
  load(paste0(prepplotdir,'indices_abs_allts_landonly.Rdata')) 
}
v <- read.table('/Volumes/DATA/climdata/forcings/forc-total-crowley_2001.txt',skip=51,header=T)[,1:2]
# 95% quantile of largest eruptions
vf <- v[600:999,1][which(v[600:999,2]<quantile(v[600:999,2],0.05))]
nc=open.nc('../data/pdo_amo/joerg/amo_yr.nc', write=F)
print.nc(nc)
tmp <- as.vector(var.get.nc(nc, "sst"))
tlist <- var.get.nc(nc, "time")
unitstr <- "months since 1500-01-15 00:00:00"
tim <- utcal.nc(unitstr, tlist)[,1]*100+utcal.nc(unitstr, tlist)[,2]
close.nc(nc)
amo_j <- cbind(trunc(tim/100),tmp)

nc=open.nc('../data/pdo_amo/joerg_v1/pdo_yr.nc', write=F)
print.nc(nc)
tmp <- as.vector(var.get.nc(nc, "sst"))
tlist <- var.get.nc(nc, "time")
unitstr <- "months since 1500-01-15 00:00:00"
tim <- utcal.nc(unitstr, tlist)[,1]*100+utcal.nc(unitstr, tlist)[,2]
close.nc(nc)
pdo_j <- cbind(trunc(tim/100),tmp)

nc=open.nc('../data/pdo_amo/joerg_v2/eofcoeff_a.nc', write=F)
print.nc(nc)
tmp <- as.vector(var.get.nc(nc, "sst_res"))
tlist <- var.get.nc(nc, "time")
unitstr <- "months since 1500-01-15 00:00:00"
tim <- trunc(tlist/100)
close.nc(nc)
pdo_j2 <- cbind(tim,tmp)

tmp <- read.table('../data/pdo_amo/sina_pdo/PDO.txt',header=T)
tmp2 <- aggregate(tmp[,3],list(tmp[,1]),mean)
pdo_s <- tmp2[which(tmp2[,1]>1599),]

pdo_m <- read.table('../data/pdo_amo/clim_expl/ipdo_mann_1600_2000.dat',header=F)

tmp <- as.matrix(read.table('../data/pdo_amo/clim_expl/ipdo_hadsst3.dat'))
pos <- which(tmp < -999)
tmp[pos] <- NA
pdo_h <- cbind(tmp[,1],apply(tmp[,2:13],1,mean))

nc=open.nc('../data/pdo_amo/sina_amo/amo.nc', write=F)
print.nc(nc)
tmp <- as.vector(var.get.nc(nc, "sst"))
tlist <- var.get.nc(nc, "time")
unitstr <- "months since 1500-01-15 00:00:00"
tim <- utcal.nc(unitstr, tlist)[,1]*100+utcal.nc(unitstr, tlist)[,2]
tmp2 <- cbind(trunc(tim/100),tmp)
amo_s <- tmp2[tmp2[,1]>1599,]

amo_m <- read.table('../data/pdo_amo/clim_expl/iamo_mann_1600_2000.dat',header=F)

tmp <- as.matrix(read.table('../data/pdo_amo/clim_expl/iamo_hadsst.dat'))
pos <- which(tmp < -999)
tmp[pos] <- NA
amo_h <- cbind(tmp[,1],apply(tmp[,2:13],1,mean))

pdf(file='../figures/',expname,'/AMO_PDO_JF_Sina_Mann_HadSST.pdf',width=8,height=4)
par(oma=c(1.5,1,1,0.1),mar=c(1.5,1,1,0.1))
par(mfrow=c(2,1))  
plot(pdo_j[,1],scale(runmean(pdo_j[,2],1),center=T,scale=T),ty='l',lwd=1,lty=1,
     ylim=c(-3,3),xlim=c(1600,2000),cex=0.5,cex.axis=0.5,main='PDO',xaxt='n',
     col=rgb(0,10,10,7,maxColorValue=10), ylab="",xlab="")
lines(pdo_j2[,1],scale(runmean(-pdo_j2[,2],1),center=T,scale=T),ty='l',lwd=1,lty=1,
      col=rgb(10,0,10,7,maxColorValue=10), ylab="",xlab="")
#lines(pdo_s[,1],scale(runmean(-pdo_s[,2],1),center=T,scale=T),ty='l',lwd=1,lty=1,
#      col=rgb(10,0,10,7,maxColorValue=10), ylab="",xlab="")
lines(pdo_m[,1],scale(runmean(pdo_m[,2],1),center=T,scale=T),ty='l',lwd=1,lty=1,
      col=rgb(0,0,10,7,maxColorValue=10), ylab="",xlab="")
lines(pdo_h[,1],scale(runmean(pdo_h[,2],1),center=T,scale=T),ty='l',lwd=1,lty=1,
      col=rgb(0,0,0,7,maxColorValue=10), ylab="",xlab="")
cold <- rep(c(1641,1646,1695,1700,1815,1820,1832,1837),each=2) #*100
warm <- rep(c(1778,1783,1791,1796,1801,1806,1826,1831),each=2) #*100
minmax <- c(-3,3,3,-3)
for (i in c(1,5,9,13)) {
  polygon(cold[i:(i+3)],minmax,col=rgb(0,0,10,3,maxColorValue=10),border=NA)
}
for (i in c(1,5,9,13)) {
  polygon(warm[i:(i+3)],minmax,col=rgb(10,0,0,3,maxColorValue=10),border=NA)
}
abline(v=vf,col='black',lty=2)
legend('topleft',c('PDO Joerg','PDO Sina','PDO Mann','PDO HadSST'),lty=1,cex=0.5,
       col=c(rgb(0,10,10,7,maxColorValue=10),rgb(10,0,10,7,maxColorValue=10),
             rgb(0,0,10,7,maxColorValue=10),rgb(0,0,0,7,maxColorValue=10)),
       bg='white', box.col='white')

plot(amo_j[,1],scale(runmean(amo_j[,2],1),center=T,scale=T),ty='l',lwd=1,lty=1,
     ylim=c(-3,3),xlim=c(1600,2000),cex=0.5,cex.axis=0.5,main='AMO',
     col=rgb(10,0,0,7,maxColorValue=10), ylab="",xlab="")
lines(amo_s[,1],scale(runmean(amo_s[,2],1),center=T,scale=T),ty='l',lwd=1,lty=1,
      col=rgb(10,7,7,7,maxColorValue=10), ylab="",xlab="")
lines(amo_m[,1],scale(runmean(amo_m[,2],1),center=T,scale=T),ty='l',lwd=1,lty=1,
      col=rgb(10,10,0,10,maxColorValue=10), ylab="",xlab="")
lines(amo_h[,1],scale(runmean(amo_h[,2],1),center=T,scale=T),ty='l',lwd=1,lty=1,
      col=rgb(10,0,10,10,maxColorValue=10), ylab="",xlab="")
#years <- rep(eind.allts$time,each=2)  
##reg=which(aind.allts$names=="PNA.calc")
#reg=which(aind.allts$names=="ENH.temp2")
#eindyrmean <- aggregate(eind.allts$ensmean[reg,],list(years),mean,na.rm=T)[,2]
#aindyrmean <- aggregate(aind.allts$ensmean[reg,],list(years),mean,na.rm=T)[,2]
##lines(eind.allts$time*100,scale(runmean(eindyrmean,1),center=T,scale=2),
##      lwd=2,lty=1,col=rgb(0,0,0,5,maxColorValue=10))
##lines(aind.allts$time*100,scale(runmean(aindyrmean,1),center=T,scale=2),
##      lwd=2,lty=1,col=rgb(10,0,0,5,maxColorValue=10))
#minmax <- c(-3,3,3,-3)
for (i in c(1,5,9,13)) {
  polygon(cold[i:(i+3)],minmax,col=rgb(0,0,10,3,maxColorValue=10),border=NA)
}
for (i in c(1,5,9,13)) {
  polygon(warm[i:(i+3)],minmax,col=rgb(10,0,0,3,maxColorValue=10),border=NA)
}
#abline(v=vf*100,col='black',lty=2)
abline(v=vf,col='black',lty=2)
#legend("topleft",c("AMO","PDO"),col=c(rgb(0,10,10,10,maxColorValue=10),
#       rgb(10,0,10,10,maxColorValue=10)),lwd=c(2,2),lty=c(1,1), 
#       cex=1, bty='o', bg='white', box.col='white')
legend('topleft',c('AMO Joerg','AMO Sina','AMO Mann','AMO HadSST'),lty=1,cex=0.5,
       col=c(rgb(10,0,0,7,maxColorValue=10),rgb(10,7,7,7,maxColorValue=10),
             rgb(10,10,0,7,maxColorValue=10),rgb(10,0,10,7,maxColorValue=10)),
       bg='white', box.col='white')
dev.off()

} #end amo









if (forcings) {

###################################################################
# plot ECHAM global mean temp. vs. scaled forcings ################
###################################################################
# load("/Volumes/DATA/climdata/PAGES_DB_Raphi/ProxyDB_1_3_3.RData")
forc <- read.table("../data/forcings/ccc400_forcings.txt",header=T)
pdf(file='../figures/',expname,'/figure2_glo_mean_scal_forc2.pdf',width=5,height=3.5)
par(mfrow=c(1,1),mar=c(4,4,1,1))
plot(forc[1:405,1],forc[1:405,7],ty='l',col='red',lwd=2,xlab="Year",
     ylab="Temperature anomaly [ºC]",ylim=c(-0.9,1.1))
lines(forc[1:405,1],forc[1:405,6],col='blue',lwd=2)
#lines(forc[250:405,1],forc[250:405,8],col='black',lwd=2)
colnames(forc)
reg <- lm(forc[250:405,7]~forc[250:405,2]+forc[250:405,3]+forc[250:405,4]+forc[250:405,5])
summary(reg)
ti <- seq(1601,2005)
fit <- reg$coefficients[1]+reg$coefficients[2]*forc[1:405,2]+reg$coefficients[3]*forc[1:405,3]+
       reg$coefficients[4]*forc[1:405,4]+reg$coefficients[5]*forc[1:405,5]
lines(ti,fit,col='black',lwd=2)
legend('topleft', c('CRUTEM4','CCC400','scaled forcings'),
       lwd=c(2,1,1),lty=c(1,1,1), col=c('red','blue','black'),
       cex=1, bty='o', bg='white',box.col='white')
dev.off()

} # en









if (landusebugbias){
  figpath="../figures/land_bug/"
  # 1605-1635 and 1971-2000 temp2 bias
  for (yr in c(1605,1971)){
    syr <- yr
    eyr <- yr+29
    print(c(syr,eyr))
    # 1. load echam mem 103 with corrected land use
    i=0
    emean103 <- array(NA,dim=c(4608,30))
    for (cyr in syr:eyr) {
      i <- i+1
      load(paste0(datadir,"EnSRF_analysis/echam_103/echam_",(cyr-1),"-",(cyr),"_2ndgrid.Rdata"))
      emean103[,i] <- apply(echam$data[which(echam$names=="temp2"),13:24,1],1,mean)
    } 
    emean103_period <- apply(emean103,1,mean)
    # 2. load analysis (oct syr-1 - sep eyr)
    amean <- array(NA,dim=c(4608,30))
    emean <- array(NA,dim=c(4608,30))
    i=0
    for (cyr in syr:eyr) {
      i=i+1
      load(paste0(datadir,"EnSRF_analysis/prepplot_v3_seasonal/analysis_",cyr,".Rdata"))
      amean[,i] <- apply(analysis$ensmean[which(echam$names=="temp2"),],1,mean)
      emean[,i] <- apply(echam$ensmean[which(echam$names=="temp2"),],1,mean)
    } 
    amean_period <- apply(amean,1,mean)
    emean_period <- apply(emean,1,mean)
    # calc bias
    bias_a_e103_period <- amean_period-(emean103_period-273.15)
    bias_e_e103_period <- emean_period-(emean103_period-273.15)
    # calc NH average
    load("/Volumes/DATA/unibe/r/data/analysis/analysis_1951.Rdata") # for landseamask
    land <- which(!is.na(validate$ensmean[1:4608,1]))
    pos <- land[land < max(which(echam$lat[1:length(which(echam$names=="temp2"))]>20))]
    #pos <- which(echam$lat[1:length(which(echam$names=="temp2"))]>0)
    lats <- echam$lat[pos]
    weights <- cos(lats/180*pi)
    weights <- weights/sum(weights)
    nhmeanbias <- sum(bias_e_e103_period[pos]*weights)
    nhbiasrange <- range(bias_e_e103_period[pos])
    
    # plot bias
    plotbias <- echam
    plotbias$ensmean <- cbind(bias_e_e103_period,bias_a_e103_period) #,dim=c(4608,2,1))
    plotbias$lon <- echam$lon[which(echam$names=="temp2")]
    plotbias$lat <- echam$lat[which(echam$names=="temp2")]
    plotbias$names <- echam$names[which(echam$names=="temp2")]
    pdf(paste(figpath,"bias_anal-ech103_",syr,"-",eyr,".pdf",sep=''), 
        width=9, height=4.5, paper='special')
    layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), height=c(3,1))
    par(oma=c(0,0,4,0))
    levs <- c(-Inf, -3,-2,-1,-0.5, -0.2, 0.2, 0.5, 1, 2, 3, Inf)
    plot_echam(plotbias, varname='temp2', type='ensmean', lty=3, lev=levs, 
               ti=1:2, colnames=c("echam mean - 103", "analysis mean - 103"), 
               main="Land use bias", add=T)
    dev.off()
  }
  # result: bias does not improve with assimilation because I only work with anomalies 
  # and add the climatology back on
}









if (arctic_slp_trend) {
#  install.packages("openair")
#  library(openair)
#  install.packages("trend")
  library(trend)
#  load("/Users/joerg/unibe/projects/EnSRF/data/indices/indices_abs_landocean_seas_allts.Rdata")  
  load("/Users/joerg/unibe/projects/EnSRF/data/indices/indices_abs_landonly_seas_allts.Rdata")  
  win <- is.odd(1:(length(aind.allts$time)*2))
  som <- is.even(1:(length(aind.allts$time)*2))
  vwin <- is.odd(1:(length(vind.allts$time)*2))
  vsom <- is.even(1:(length(vind.allts$time)*2))
  winlen <- 11 # window length for moving average
#  pdf("../figures/2nd_analysis_paper/arctic_slp_landocean.pdf",width=6,height=5)
  pdf("../figures/2nd_analysis_paper/arctic_slp_landonly.pdf",width=6,height=5)
    par(mfrow=c(2,1),mar=c(2,2,1,0.2))
    plot(filter(aind.allts$ensmean[aind.allts$names=="NEU.slp",win],filter=rep(1/winlen,winlen),
                 method="convolution",sides=2),ty="l",col="red",ylim=c(1010,1018),main='Winter')
    abline(h=mean(aind.allts$ensmean[aind.allts$names=="NEU.slp",win]),lty=3,col="red")
    lines(filter(aind.allts$ensmean[aind.allts$names=="ARC.slp",win],filter=rep(1/winlen,winlen),
                 method="convolution",sides=2),ty="l",col="blue")
    abline(h=mean(aind.allts$ensmean[aind.allts$names=="ARC.slp",win]),lty=3,col="blue")
    lines(filter(eind.allts$ensmean[aind.allts$names=="ARC.slp",win],filter=rep(1/winlen,winlen),
                 method="convolution",sides=2),ty="l",col="black")
    plot(filter(aind.allts$ensmean[aind.allts$names=="NEU.slp",som],filter=rep(1/winlen,winlen),
                 method="convolution",sides=2),ty="l",col="red",ylim=c(1012,1016),main='Summer')
    abline(h=mean(aind.allts$ensmean[aind.allts$names=="NEU.slp",som]),lty=3,col="red")
    lines(filter(aind.allts$ensmean[aind.allts$names=="ARC.slp",som],filter=rep(1/winlen,winlen),
                 method="convolution",sides=2),ty="l",col="blue")
    abline(h=mean(aind.allts$ensmean[aind.allts$names=="ARC.slp",som]),lty=3,col="blue")
    legend('bottomleft',c('Arctic','N-Europe'),col=c('blue','red'),lty=c(1,1),bty='n') 
  dev.off()

#  tmpdf <- cbind(aind.allts$time,aind.allts$ensmean[aind.allts$names=="ARC.slp",win]) 
#  colnames(tmpdf) <- c("date","var")
#  TheilSen(tmpdf,pollutant="var")
  tmpts <- ts(aind.allts$ensmean[aind.allts$names=="ARC.slp",win],
              start=aind.allts$time[1],freq=1)
  sens.slope(tmpts) # all 400 yrs
  sens.slope(window(tmpts,1603,1700)) 
  sens.slope(window(tmpts,1701,1800))
  sens.slope(window(tmpts,1801,1900))
  sens.slope(window(tmpts,1901,2000))
  sens.slope(window(tmpts,1850,2000))
  
  # Z100
  load("/Users/joerg/unibe/projects/EnSRF/data/indices/indices_abs_landonly_seas_allts.Rdata")
  pdf("../figures/2nd_analysis_paper/z100_landonly.pdf",width=6,height=5)
  plot(filter(aind.allts$ensmean[aind.allts$names=="Z100",win],filter=rep(1/winlen,winlen),
              method="convolution",sides=2),ty="l",col="blue")
  lines(filter(eind.allts$ensmean[aind.allts$names=="Z100",win],filter=rep(1/winlen,winlen),
              method="convolution",sides=2),ty="l",col="black")
#  vindts <- 
#  lines(filter(vind.allts$data[vind.allts$names=="Z100",vwin],filter=rep(1/winlen,winlen),
#               method="convolution",sides=2),ty="l",col="blue")
  abline(h=mean(aind.allts$ensmean[aind.allts$names=="Z100",win]),lty=3,col="blue")
  
  plot(filter(aind.allts$ensmean[aind.allts$names=="Z100",win],filter=rep(1/winlen,winlen),
              method="convolution",sides=2),ty="l",col="blue")
  lines(filter(eind.allts$ensmean[aind.allts$names=="Z100",win],filter=rep(1/winlen,winlen),
               method="convolution",sides=2),ty="l",col="black")
  abline(h=mean(aind.allts$ensmean[aind.allts$names=="Z100",win]),lty=3,col="blue")
  
  # PNA; only in landocean indices
  load("/Users/joerg/unibe/projects/EnSRF/data/indices/indices_abs_landocean_seas_allts.Rdata")
  winlen=1
  plot(filter(aind.allts$ensmean[aind.allts$names=="PNA.calc",win],filter=rep(1/winlen,winlen),
              method="convolution",sides=2),ty="l",col="blue")
  lines(filter(eind.allts$ensmean[aind.allts$names=="PNA.calc",win],filter=rep(1/winlen,winlen),
              method="convolution",sides=2),ty="l",col="black")
  dev.off()
}









if (updowntrend_decades) {
  load(paste0(dataintdir,'indices/indices_abs_landonly_seas_allts.Rdata'))  
  
  # find and mark steepest warming/cooling decade 
  d_sum <- aind.allts$ensmean[(which(eind.allts$names=="ENH.temp2")),
                          (is.even(seq(1,ncol(aind.allts$ensmean))))]
  d_win <- aind.allts$ensmean[(which(eind.allts$names=="ENH.temp2")),
                              (!is.even(seq(1,ncol(aind.allts$ensmean))))]
  d_ann <- rowMeans(cbind(d_sum,d_win))
  year <- aind.allts$time
  # calc trend over decades of smoothed to not overweight single extreme years 
  plot(d_ann,ty='l')
  lines(predict(sm.spline(year, d_ann), year, 0),ty='l',col='red') # library(pspline) ,1) for 1st derivative
  d_ann_sm <- predict(sm.spline(year, d_ann), year, 0)
  dec_temp <- rep(NA,300)
  dec_temptr <- rep(NA,300)
  dectime <- rep(NA,300)
  for (ti in 1:length(dec_temp)) {
    if (ti==1) {
      dec_temp[ti] <- mean(d_ann[1:10])
      dec_temptr[ti] <- sens.slope(ts(d_ann_sm[1:10]))[[1]]
      dectime[ti] <- 1603
      i=2
    } else {
      dec_temp[ti] <- mean(d_ann[i:(i+9)])
      dec_temptr[ti] <- sens.slope(ts(d_ann_sm[i:(i+9)]))[[1]]
      dectime[ti] <- dectime[ti-1] + 1 
      i=i+1
    }
  } 
  #  plot(dectime,dec,ty='l',col='red')
  tmp=cbind(dec_temp,dectime)
  tmp[order(dec_temp),]
  # warmest/coldest decade before 1900
  # warm: 1894, 1798, 1874, 1786, 1773, 1822
  # cold: 1641, 1695, 1614, 1674, 1832, 1809
  tmp2=cbind(dec_temptr,dectime)
  tmp2[order(dec_temptr),]
  # steepest warming/cooling decade before 1900
  # warm: 1646, 1838, 1816, 1768, 1890, 1860
  # cool: 1803, 1634, 1899, 1826, 1689, 1877
  
  
  warmper <- c(1646, 1838, 1816, 1768) #, 1890, 1860) 
  for (syr in warmper) {
    eyr <- syr+9
    # change to precip in % 
    # load monthly 70yr climatology data
    if (syr < 1636) {
      load(file=paste0(echclimpath,'echam_clim_1636-1637_2ndgrid.Rdata'))
    } else {
      load(file=paste0(echclimpath,'echam_clim_',syr+3,'-',(syr+4),'_2ndgrid.Rdata'))  
    }
    pos <- which(echam_clim$names=="precip")[1:4608]
    # seasonal mean
    echmeanclim <- cbind(rowMeans(echam_clim$ensmean[,10:15]),
                         rowMeans(echam_clim$ensmean[,16:21]))
    for (cyr in seq(syr,eyr)) {
      load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
      # change to precip in % 
      echam.abs <- echam
      echam.anom$ensmean[pos,] <- echam.abs$ensmean[pos,]/(echmeanclim[pos,]*3600*24*30)*100
      analysis.abs <- analysis
      analysis.anom$ensmean[pos,] <- analysis.abs$ensmean[pos,]/
        (echmeanclim[pos,]*3600*24*30)*100
      if (cyr==syr){
        atmp <- analysis.anom$ensmean
        etmp <- echam.anom$ensmean
        c=1
      } else {
        atmp <- atmp+analysis.anom$ensmean
        etmp <- etmp+echam.anom$ensmean
        c=c+1
      }
    }
    # read mann et al temp recon
    mann <- read.table(paste0(dataextdir,"mann_temp_recon/allproxyfieldrecon"),na.string="NaN")
    ll <- read.table(paste0(dataextdir,"mann_temp_recon/longlat.txt"),na.string="NaN")
    ts <- which(mann[,1]==syr)
    plotmann <- list() #echam
    mclim <- as.vector(t(colMeans(mann[((ts-35):(ts+39)),-1],na.rm=F)))
    mdata <- as.vector(t(colMeans(mann[ts:(ts+4),-1],na.rm=F)))
    manom <- mdata-mclim
    plotmann$data <- array(manom,c(length(manom),1))
    plotmann$time <- mann[ts,1]
    plotmann$names <- rep('temp2',dim(plotmann$data)[1])
    plotmann$lon <- ll[,1] 
    plotmann$lat <- ll[,2]
    plotmann$height <- ll[,3]   
    plotmann$lsm.i <- ll[,3]  
    if (syr==warmper[1]){
      anameananom <- atmp/c
      echmeananom <- etmp/c
      mannmeananom <- plotmann$data
    } else {
      anameananom <- cbind(anameananom,(atmp/c))
      echmeananom <- cbind(echmeananom,(etmp/c))
      mannmeananom <- cbind(mannmeananom,plotmann$data)
    }
  } # end periods loop  
  colnames(anameananom) <- rep(warmper,each=2)
  colnames(echmeananom) <- rep(warmper,each=2)
  colnames(mannmeananom) <- warmper
  
  # 10-yr period winter
  pdata <- echam
  for (ti in 1:3){
    if (ti==1){s="winter"} else if (ti==2){s="summer"} else {s="annual"}
    pdf(paste0('../figures/',expname,'/warming_trend_',s,filenameext,'.pdf'), width=16, 
        height=12, paper='special')
    layout(matrix(c(1,2,3,4,5,5,5,5,6,7,8,9,10,10,10,10,11,12,13,14,15,15,15,15,
                    16,17,18,19,20,20,20,20),8,4,byrow=TRUE),height=c(3,1,3,1,3,1,3,1))
    par(oma=c(0.5,4,4,0))
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    plotmann$data <- array(c(mannmeananom[,1],mannmeananom[,2],mannmeananom[,3],
                             mannmeananom[,4]), c(nrow(plotmann$data),1,4))
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    plot_echam3(plotmann, varname='temp2', type='data', cex.pt=2.5, 
                latlim=c(0,90),lonlim=c(-180,180),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('a)', 'b)', 'c)', 'd)'),colnames=warmper[1:4],
                statpanel=NULL, add=T, rownames='Mann recon. temp.', 
                main='', units='K')
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    scfac <- max(echam.anom$ensmean[echam$names=='u200'])
    if (ti==3){
      pdata$data <- array(c(rowMeans(anameananom[,1:2]),rowMeans(anameananom[,3:4]),
                            rowMeans(anameananom[,5:6]),rowMeans(anameananom[,7:8])),
                          c(nrow(echam.anom$ensmean),1,4))
    } else {
      pdata$data <- array(c(anameananom[,ti],anameananom[,ti+2],anameananom[,ti+4],
                            anameananom[,ti+6]), c(nrow(echam.anom$ensmean),1,4))
    }
    plot_echam3(pdata, varname='temp2', type='data', cex.pt=1.5, 
                latlim=c(0,90),lonlim=c(-180,180),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('e)', 'f)', 'g)', 'h)'),colnames=rep("",4),
                statpanel=NULL, add=T, rownames='Analysis temp, Z500, UV200', 
                main=paste('EKF400 8-yr',s,'periods'), units='K',
                addcontours=T, contvarname='gph500', conttype='data',
                contcol='black', contlev=contlevs,
                addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                veclen=scfac*0.01, vecscale=scfac*2.0, vecwd=0.95, every_x_vec=4)
    # diff echam analysis
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    if (ti==3){
      pdata$data <- array(c((rowMeans(anameananom[,1:2])-rowMeans(echmeananom[,1:2])),
                            (rowMeans(anameananom[,3:4])-rowMeans(echmeananom[,3:4])),
                            (rowMeans(anameananom[,5:6])-rowMeans(echmeananom[,5:6])),
                            (rowMeans(anameananom[,7:8])-rowMeans(echmeananom[,7:8]))),
                          c(nrow(echam.anom$ensmean),1,4))
    } else {
      pdata$data <- array(c((anameananom[,ti]-echmeananom[,ti]),
                            (anameananom[,ti+2]-echmeananom[,ti+2]),
                            (anameananom[,ti+4]-echmeananom[,ti+4]),
                            (anameananom[,ti+6]-echmeananom[,ti+6])),
                          c(nrow(echam.anom$ensmean),1,4))
    }
    plot_echam3(pdata, varname='temp2', type='data', cex.pt=1.5, 
                latlim=c(0,90),lonlim=c(-180,180),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('i)', 'j)', 'k)', 'l)'),colnames=rep(" ",4),
                statpanel=NULL, add=T, rownames='Update (Anal.-ECHAM)', 
                units='K',
                addcontours=T, contvarname='gph500', conttype='data',
                contcol='black', contlev=contlevs,
                addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                veclen=scfac*0.01, vecscale=scfac*2.0, vecwd=0.95, every_x_vec=4)
    if (ti==3){
      pdata$data <- array(c(rowMeans(anameananom[,1:2]),rowMeans(anameananom[,3:4]),
                            rowMeans(anameananom[,5:6]),rowMeans(anameananom[,7:8])),
                          c(nrow(echam.anom$ensmean),1,4))
    } else {
      pdata$data <- array(c(anameananom[,ti],anameananom[,ti+2],anameananom[,ti+4],
                            anameananom[,ti+6]), c(nrow(echam.anom$ensmean),1,4))
    }
    levs=c(0,80,85,90,95,100,105,110,115,120,Inf)
    plot_echam3(pdata, varname='precip', type='data', cex.pt=1.5, 
                latlim=c(0,90),lonlim=c(-180,180),units="%",
                names=c('m)', 'n)','o)', 'p)'), lev=levs, 
                st.col=NULL, stations=NULL, add=T,addcontours=F,
                wcol='darkgrey',rownames='Prec. analysis',colnames=rep('',4))
    dev.off()  
  } # end summer winter loop
  
  
  
  coldper <- c(1803, 1634, 1899, 1826) #, 1689, 1877) 
  for (syr in coldper) {
    eyr <- syr+9
    # change to precip in % 
    # load monthly 70yr climatology data
    if (syr < 1636) {
      load(file=paste0(echclimpath,'echam_clim_1636-1637_2ndgrid.Rdata'))
    } else {
      load(file=paste0(echclimpath,'echam_clim_',syr+3,'-',(syr+4),'_2ndgrid.Rdata'))  
    }
    pos <- which(echam_clim$names=="precip")[1:4608]
    # seasonal mean
    echmeanclim <- cbind(rowMeans(echam_clim$ensmean[,10:15]),
                         rowMeans(echam_clim$ensmean[,16:21]))
    for (cyr in seq(syr,eyr)) {
      load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
      # change to precip in % 
      echam.abs <- echam
      echam.anom$ensmean[pos,] <- echam.abs$ensmean[pos,]/(echmeanclim[pos,]*3600*24*30)*100
      analysis.abs <- analysis
      analysis.anom$ensmean[pos,] <- analysis.abs$ensmean[pos,]/
        (echmeanclim[pos,]*3600*24*30)*100
      if (cyr==syr){
        atmp <- analysis.anom$ensmean
        etmp <- echam.anom$ensmean
        c=1
      } else {
        atmp <- atmp+analysis.anom$ensmean
        etmp <- etmp+echam.anom$ensmean
        c=c+1
      }
    }
    # read mann et al temp recon
    mann <- read.table(paste0(dataextdir,"mann_temp_recon/allproxyfieldrecon"),na.string="NaN")
    ll <- read.table(paste0(dataextdir,"mann_temp_recon/longlat.txt"),na.string="NaN")
    ts <- which(mann[,1]==syr)
    plotmann <- list() #echam
    mclim <- as.vector(t(colMeans(mann[((ts-35):(ts+39)),-1],na.rm=F)))
    mdata <- as.vector(t(colMeans(mann[ts:(ts+4),-1],na.rm=F)))
    manom <- mdata-mclim
    plotmann$data <- array(manom,c(length(manom),1))
    plotmann$time <- mann[ts,1]
    plotmann$names <- rep('temp2',dim(plotmann$data)[1])
    plotmann$lon <- ll[,1] 
    plotmann$lat <- ll[,2]
    plotmann$height <- ll[,3]   
    plotmann$lsm.i <- ll[,3]  
    if (syr==coldper[1]){
      anameananom <- atmp/c
      echmeananom <- etmp/c
      mannmeananom <- plotmann$data
    } else {
      anameananom <- cbind(anameananom,(atmp/c))
      echmeananom <- cbind(echmeananom,(etmp/c))
      mannmeananom <- cbind(mannmeananom,plotmann$data)
    }
  } # end periods loop  
  colnames(anameananom) <- rep(coldper,each=2)
  colnames(echmeananom) <- rep(coldper,each=2)
  colnames(mannmeananom) <- coldper
  
  # 10-yr period winter
  pdata <- echam
  for (ti in 1:3){
    if (ti==1){s="winter"} else if (ti==2){s="summer"} else {s="annual"}
    pdf(paste0('../figures/',expname,'/cooling_trend_',s,filenameext,'.pdf'), width=16, 
        height=12, paper='special')
    layout(matrix(c(1,2,3,4,5,5,5,5,6,7,8,9,10,10,10,10,11,12,13,14,15,15,15,15,
                    16,17,18,19,20,20,20,20),8,4,byrow=TRUE),height=c(3,1,3,1,3,1,3,1))
    par(oma=c(0.5,4,4,0))
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    plotmann$data <- array(c(mannmeananom[,1],mannmeananom[,2],mannmeananom[,3],
                             mannmeananom[,4]), c(nrow(plotmann$data),1,4))
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    plot_echam3(plotmann, varname='temp2', type='data', cex.pt=2.5, 
                latlim=c(0,90),lonlim=c(-180,180),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('a)', 'b)', 'c)', 'd)'),colnames=coldper[1:4],
                statpanel=NULL, add=T, rownames='Mann recon. temp.', 
                main='', units='K')
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    scfac <- max(echam.anom$ensmean[echam$names=='u200'])
    if (ti==3){
      pdata$data <- array(c(rowMeans(anameananom[,1:2]),rowMeans(anameananom[,3:4]),
                            rowMeans(anameananom[,5:6]),rowMeans(anameananom[,7:8])),
                          c(nrow(echam.anom$ensmean),1,4))
    } else {
      pdata$data <- array(c(anameananom[,ti],anameananom[,ti+2],anameananom[,ti+4],
                            anameananom[,ti+6]), c(nrow(echam.anom$ensmean),1,4))
    }
    plot_echam3(pdata, varname='temp2', type='data', cex.pt=1.5, 
                latlim=c(0,90),lonlim=c(-180,180),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('e)', 'f)', 'g)', 'h)'),colnames=rep("",4),
                statpanel=NULL, add=T, rownames='Analysis temp, Z500, UV200', 
                main=paste('EKF400 8-yr',s,'periods'), units='K',
                addcontours=T, contvarname='gph500', conttype='data',
                contcol='black', contlev=contlevs,
                addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                veclen=scfac*0.01, vecscale=scfac*2.0, vecwd=0.95, every_x_vec=4)
    # diff echam analysis
    levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
    contlevs <- seq(-4,4,2)
    if (ti==3){
      pdata$data <- array(c((rowMeans(anameananom[,1:2])-rowMeans(echmeananom[,1:2])),
                            (rowMeans(anameananom[,3:4])-rowMeans(echmeananom[,3:4])),
                            (rowMeans(anameananom[,5:6])-rowMeans(echmeananom[,5:6])),
                            (rowMeans(anameananom[,7:8])-rowMeans(echmeananom[,7:8]))),
                          c(nrow(echam.anom$ensmean),1,4))
    } else {
      pdata$data <- array(c((anameananom[,ti]-echmeananom[,ti]),
                            (anameananom[,ti+2]-echmeananom[,ti+2]),
                            (anameananom[,ti+4]-echmeananom[,ti+4]),
                            (anameananom[,ti+6]-echmeananom[,ti+6])),
                          c(nrow(echam.anom$ensmean),1,4))
    }
    plot_echam3(pdata, varname='temp2', type='data', cex.pt=1.5, 
                latlim=c(0,90),lonlim=c(-180,180),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('i)', 'j)', 'k)', 'l)'),colnames=rep(" ",4),
                statpanel=NULL, add=T, rownames='Update (Anal.-ECHAM)', 
                units='K',
                addcontours=T, contvarname='gph500', conttype='data',
                contcol='black', contlev=contlevs,
                addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                veclen=scfac*0.01, vecscale=scfac*2.0, vecwd=0.95, every_x_vec=4)
    if (ti==3){
      pdata$data <- array(c(rowMeans(anameananom[,1:2]),rowMeans(anameananom[,3:4]),
                            rowMeans(anameananom[,5:6]),rowMeans(anameananom[,7:8])),
                          c(nrow(echam.anom$ensmean),1,4))
    } else {
      pdata$data <- array(c(anameananom[,ti],anameananom[,ti+2],anameananom[,ti+4],
                            anameananom[,ti+6]), c(nrow(echam.anom$ensmean),1,4))
    }
    levs=c(0,80,85,90,95,100,105,110,115,120,Inf)
    plot_echam3(pdata, varname='precip', type='data', cex.pt=1.5, 
                latlim=c(0,90),lonlim=c(-180,180),units="%",
                names=c('m)', 'n)','o)', 'p)'), lev=levs, 
                st.col=NULL, stations=NULL, add=T,addcontours=F,
                wcol='darkgrey',rownames='Prec. analysis',colnames=rep('',4))
    dev.off()  
  } # end summer winter loop
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # boxplots of indices
  load(paste0(dataintdir,'indices/indices_abs_landocean_seas_allts.Rdata'))  
  
  seas1 <- 'annual mean'
  seas2 <- 'winter'
  seas3 <- 'summer'
  period <- eind.allts$time #seq(syr,eyr)
  for (seas in c(seas1,seas2,seas3)) {
    print(seas)
    pdf(file=paste0('../figures/',expname,'/indices_boxplots',filenameext,seas,'.pdf'),
      width=20,height=15)
      par(oma=c(.5,1,1,0),mar=c(2,2,2,2))
      par(mfrow=c(5,3))  
    
    # warmest/coldest decade before 1900
      warmest <- c(1894, 1798, 1874, 1786, 1773, 1822)
      coldest <- c(1641, 1695, 1614, 1674, 1832, 1809)
    # steepest warming/cooling decade before 1900
      warm_trend <- c(1646, 1838, 1816, 1768, 1890, 1860)
      cool_trend <- c(1803, 1634, 1899, 1826, 1689, 1877)
    
      for (reg in c(35:49)) { 
        # "HC.calc""ITCZ.calc""SJ_u200.calc""SJ_slp.calc""PV.calc""PWC.calc"
        # "DIMI.calc""NAO.calc""PNA.calc""DIMI""HC""SJ""Z100""Z300""PWC"
        years <- rep(eind.allts$time,each=2)
        
        if (seas=='annual mean') { 
          winlen <- 39 
          sampsize <- 3600
          b1900 <- 1:(2*which(aind.allts$time==1900))
          p1900 <- (2*which(aind.allts$time==1900)):(2*which(aind.allts$time==2004))
        } else {
          winlen <- 19 
          sampsize <- 1800
        }
        if (seas=='summer') { 
          b1900 <- (1:(2*which(aind.allts$time==1900)))[is.even(1:(2*which(aind.allts$time==1900)))]
          p1900 <- ((2*which(aind.allts$time==1900)):(2*which(aind.allts$time==2004)))[
            is.even((2*which(aind.allts$time==1900)):(2*which(aind.allts$time==2004)))]
        }
        if (seas=='winter') { 
          b1900 <- (1:(2*which(aind.allts$time==1900)))[is.odd(1:(2*which(aind.allts$time==1900)))]
          p1900 <- ((2*which(aind.allts$time==1900)):(2*which(aind.allts$time==2004)))[
            is.odd((2*which(aind.allts$time==1900)):(2*which(aind.allts$time==2004)))]
        }
         
        ltmp2 <- which(years%in%warmest)
        if (seas=='annual mean') { 
          ltmp <- ltmp2[is.odd(ltmp2)]
        } else if (seas=='summer') { 
          ltmp <- ltmp2[is.even(ltmp2)] 
        } else if (seas=='winter') { 
          ltmp <- ltmp2[is.odd(ltmp2)] 
        }
        lw=NULL
        for (pos in ltmp){
          if (is.null(lw)){
            lw=seq(pos,pos+winlen,2)
          } else {
            lw <- c(lw,seq(pos,pos+winlen,2))
          }
        }
        
        ltmp2 <- which(years%in%coldest)
        if (seas=='annual mean') { 
          ltmp <- ltmp2[is.odd(ltmp2)]
        } else if (seas=='summer') { 
          ltmp <- ltmp2[is.even(ltmp2)] 
        } else if (seas=='winter') { 
          ltmp <- ltmp2[is.odd(ltmp2)] 
        }
        lc=NULL
        for (pos in ltmp){
          if (is.null(lc)){
            lc=seq(pos,pos+winlen,2)
          } else {
            lc <- c(lc,seq(pos,pos+winlen,2))
          }
        }
        
        ltmp2 <- which(years%in%warm_trend)
        if (seas=='annual mean') { 
          ltmp <- ltmp2[is.odd(ltmp2)]
        } else if (seas=='summer') { 
          ltmp <- ltmp2[is.even(ltmp2)] 
        } else if (seas=='winter') { 
          ltmp <- ltmp2[is.odd(ltmp2)] 
        }
        lwt=NULL
        for (pos in ltmp){
          if (is.null(lwt)){
            lwt=seq(pos,pos+winlen,2)
          } else {
            lwt <- c(lwt,seq(pos,pos+winlen,2))
          }
        }
        
        ltmp2 <- which(years%in%cool_trend)
        if (seas=='annual mean') { 
          ltmp <- ltmp2[is.odd(ltmp2)]
        } else if (seas=='summer') { 
          ltmp <- ltmp2[is.even(ltmp2)] 
        } else if (seas=='winter') { 
          ltmp <- ltmp2[is.odd(ltmp2)] 
        }
        lct=NULL
        for (pos in ltmp){
          if (is.null(lct)){
            lct=seq(pos,pos+winlen,2)
          } else {
            lct <- c(lct,seq(pos,pos+winlen,2))
          }
        }
# ANNUAL MEAN STILL HAS BUG: SOM AND WIN COLS NEED TO BE AVERAGED PROPERLY
        datmat <- t(rbind(sample(aind.allts$data[reg,b1900,],sampsize),
                     sample(eind.allts$data[reg,b1900,],sampsize),
                     rep(NA,sampsize),
                     sample(as.vector(aind.allts$data[reg,p1900,]),sampsize),   
                     sample(as.vector(eind.allts$data[reg,p1900,]),sampsize), 
                     rep(NA,sampsize),
                     as.vector(aind.allts$data[reg,lw,]),as.vector(eind.allts$data[reg,lw,]),
                     rep(NA,sampsize),
                     as.vector(aind.allts$data[reg,lc,]),as.vector(eind.allts$data[reg,lc,]),
                     rep(NA,sampsize),
                     as.vector(aind.allts$data[reg,lwt,]),as.vector(eind.allts$data[reg,lwt,]),
                     rep(NA,sampsize),
                     as.vector(aind.allts$data[reg,lct,]),as.vector(eind.allts$data[reg,lct,])))
        if (!all(is.nan(datmat))){
        #datmat <- apply(aind.allts$data[reg,l,],1L,c)
          boxplot(datmat,col=c("green","darkgreen","white","cyan","darkcyan","white",
                    "red","darkred","white","blue","darkblue","white",
                    "orange","darkorange","white","magenta","darkmagenta"),
                   main=eind.allts$names[reg],xaxt='n',notch=T,outline=F)
          mtext(#1:12,rep(min(aind.allts$data[reg,,]),12),
                c("Clim. <1900",""," ","Clim. >1900","",'',"Warm Dec.",""," ","Cold Dec.",
                  ""," ","Warm. Trend",""," ","Cool. Trend",""),
                  side=1,at=1:17,cex=0.7) #srt=90,pos=1)
        }
      }
      dev.off()
    } # end seasonal loop
  
  # mark warmest/coldest decade before 1800
  # depends on definition, based on abs. values or anomalies?    
  #    polygon(c(1790,1800,1800,1790),c(0.9,0.9,1,1),density=NA,col=rgb(5,0,0,3,maxColorValue=10))
  #    polygon(c(1640,1650,1650,1640),c(0.9,0.9,1,1),density=NA,col=rgb(0,0,5,3,maxColorValue=10))
} # end warming/cooling trend









if (laki) {
  # load pre calc time series
  filenameext = '_abs_landonly_seas_'
  load(paste0('../data/indices/indices',filenameext,'allts.Rdata'))
  
  # define plot parameters
  wlen=1
  linew=1
  scal=F
  plbem=F
  syr=1750
  eyr=1850
  # seas='sum'
  # reg=which(aind.allts$names=='GLO.temp2')
  pdf('../figures/laki/t2m_ts_1750-1850.pdf',width=8,height=6)
  par(mfrow=c(2,1),mar=c(2,2,2,0.2))
  for (seas in c('sum','win')) {
    c=1
    for (reg in c(which(aind.allts$names=='ENH.temp2'),which(aind.allts$names=='NAM.temp2'),
                  which(aind.allts$names=='EU.temp2'))) {
      print(reg)
      if (seas=='yrmean') {
        years <- rep(eind.allts$time,each=2)  
      } else {
        wincol <- is.odd(seq(1,ncol(aind.allts$ensmean)))
        somcol <- is.even(seq(1,ncol(aind.allts$ensmean)))  
      }
      if (seas=='yrmean') {
        print("annual mean")
        yl=c(0,15)
        eindmean <- aggregate(eind.allts$ensmean[reg,],list(years),mean,na.rm=T)[,2]
        aindmean <- aggregate(aind.allts$ensmean[reg,],list(years),mean,na.rm=T)[,2]
      } else if (seas=='sum'){ 
        yl=c(13,16.5)
        print("summer season")
        eindmean <- eind.allts$ensmean[reg,somcol]  
        aindmean <- aind.allts$ensmean[reg,somcol]  
      } else if (seas=='win'){
        print("winter season")
        yl=c(-5,4)
        eindmean <- eind.allts$ensmean[reg,wincol]  
        aindmean <- aind.allts$ensmean[reg,wincol]  
      }
      # northern hemisphere extra-tropics (>20N) seas="yrmean" or "sum" or "win"
      yrpos=which(aind.allts$time>syr&aind.allts$time<eyr)
      yrs=aind.allts$time[yrpos]
      if (c==1) {
        plot(yrs,aindmean[yrpos],ty="l",ylim=yl,main=seas)
        abline(v=c(1783,1809,1815,1831,1835),col='grey')
        abline(h=mean(aindmean[yrpos]),col=c,lty=2)
      } else {
        lines(yrs,aindmean[yrpos],ty="l",col=c)
        abline(h=mean(aindmean[yrpos]),col=c,lty=2)
      }
      legend("topleft",c('ENH','NAM','EU'),col=c(1,2,4),lty=1,bg='white') #,box.col='white')
      c=c+1
      if (c==3) {c=c+1}
    } # end reg loop for various regions
  } # end sum win season loop
  dev.off()
  
 
  
  
  
  
  # spatial maps
  coldper <- c(1783,1809,1815,1831,1835)
  for (syr in coldper) {
    eyr <- syr+3
    # change to precip in % 
    # load monthly 70yr climatology data
    if (syr < 1636) {
      load(file=paste0(echclimpath,'echam_clim_1636-1637_2ndgrid.Rdata'))
    } else {
      load(file=paste0(echclimpath,'echam_clim_',syr+3,'-',(syr+4),'_2ndgrid.Rdata'))  
    }
    pos <- which(echam_clim$names=="precip")[1:4608]
    # seasonal mean
    echmeanclim <- cbind(rowMeans(echam_clim$ensmean[,10:15]),
                         rowMeans(echam_clim$ensmean[,16:21]))
    for (cyr in seq(syr,eyr)) {
      load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
      # change to precip in % 
      echam.abs <- echam
      echam.anom$ensmean[pos,] <- echam.abs$ensmean[pos,]/(echmeanclim[pos,]*3600*24*30)*100
      analysis.abs <- analysis
      analysis.anom$ensmean[pos,] <- analysis.abs$ensmean[pos,]/
        (echmeanclim[pos,]*3600*24*30)*100
      if (cyr==syr) {
        anameananom <- analysis.anom$ensmean
      } else {
        anameananom <- cbind(anameananom,analysis.anom$ensmean)
      }
    }
 
    # plot post volc years
    pdata <- echam
    ti=1
    pdf(paste0('../figures/laki/post_erup_anom_',syr,'.pdf'), width=16, 
        height=12, paper='special')
      layout(matrix(c(1,2,3,4,5,5,5,5,6,7,8,9,10,10,10,10,11,12,13,14,15,15,15,15,
                      16,17,18,19,20,20,20,20),8,4,byrow=TRUE),height=c(3,1,3,1,3,1,3,1))
      par(oma=c(0.5,4,4,0))
      levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
      contlevs <- seq(-4,4,2)
      scfac <- max(echam.anom$ensmean[echam$names=='u200'])
      pdata$data <- array(c(anameananom[,ti],anameananom[,ti+2],anameananom[,ti+4],
                            anameananom[,ti+6]), c(nrow(echam.anom$ensmean),1,4))
      plot_echam3(pdata, varname='temp2', type='data', cex.pt=1.5, 
                  latlim=c(0,90),lonlim=c(-180,180),
                  wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                  names=c('e)', 'f)', 'g)', 'h)'),colnames=seq(syr,eyr),
                  statpanel=NULL, add=T, rownames='T2m, Z500, UV200 (Oct-Mar)', 
                  main='EKF400', units='K',
                  addcontours=T, contvarname='gph500', conttype='data',
                  contcol='black', contlev=contlevs,
                  addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                  veclen=scfac*0.01, vecscale=scfac*2.0, vecwd=0.95, every_x_vec=4)
      pdata$data <- array(c(anameananom[,ti],anameananom[,ti+2],anameananom[,ti+4],
                              anameananom[,ti+6]), c(nrow(echam.anom$ensmean),1,4))
      levs=c(0,80,85,90,95,100,105,110,115,120,Inf)
      plot_echam3(pdata, varname='precip', type='data', cex.pt=1.5, 
                  latlim=c(0,90),lonlim=c(-180,180),units="%",
                  names=c('m)', 'n)','o)', 'p)'), lev=levs, 
                  st.col=NULL, stations=NULL, add=T,addcontours=F,
                  wcol='darkgrey',rownames='Prec. (Oct-Mar)',colnames=rep('',4))
      
      ti=ti+1 # Apr-Sep season
      levs <- c(-Inf, seq(-0.5,0.5,0.1), Inf)
      pdata$data <- array(c(anameananom[,ti],anameananom[,ti+2],anameananom[,ti+4],
                            anameananom[,ti+6]), c(nrow(echam.anom$ensmean),1,4))
      plot_echam3(pdata, varname='temp2', type='data', cex.pt=1.5, 
                  latlim=c(0,90),lonlim=c(-180,180),
                  wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                  names=c('e)', 'f)', 'g)', 'h)'),colnames=rep(' ',4),
                  statpanel=NULL, add=T, rownames='T2m, Z500, UV200  (Apr-Sep)', 
                  main=' ', units='K',
                  addcontours=T, contvarname='gph500', conttype='data',
                  contcol='black', contlev=contlevs,
                  addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                  veclen=scfac*0.01, vecscale=scfac*2.0, vecwd=0.95, every_x_vec=4)
      pdata$data <- array(c(anameananom[,ti],anameananom[,ti+2],anameananom[,ti+4],
                            anameananom[,ti+6]), c(nrow(echam.anom$ensmean),1,4))
      levs=c(0,80,85,90,95,100,105,110,115,120,Inf)
      plot_echam3(pdata, varname='precip', type='data', cex.pt=1.5, 
                  latlim=c(0,90),lonlim=c(-180,180),units="%",
                  names=c('m)', 'n)','o)', 'p)'), lev=levs, 
                  st.col=NULL, stations=NULL, add=T,addcontours=F,
                  wcol='darkgrey',rownames='Prec. (Apr-Sep)',colnames=rep('',4))
    dev.off()  
  } # end periods loop
  
  
  
  
  
  
  # mik's weathertype recon
  wt <- read.table('../comparison_data/CAP7_1763-2009.txt',stringsAsFactors=F,head=T)
  wtnt <- read.table('../comparison_data/CAP7_1763-2009_NoTemp.txt',stringsAsFactors=F,head=T)
  wtleg=c('NE','WSW','W','E','HP','N','WC')
  # pos <- which(as.numeric(substr((wt[,1]),1,4))==1783)
  # h = hist(as.numeric(wt[pos,2])-0.2,breaks=seq(0.66,7.33,0.33),plot=F)
  # h$density = h$counts/sum(h$counts)*100
  # plot(h,col='blue',xlim=c(0.5,7.5),ylim=c(0,25),main='CAP7 1783',freq=F,xaxt='n')
  # axis(side=1, at=seq(1,7), labels=wtleg)
  # h2 = hist(as.numeric(wtnt[pos,2])+0.2,breaks=seq(0.66,7.33,0.33),plot=F)
  # h2$density = h2$counts/sum(h2$counts)*100
  # par(new=T)
  # plot(h2,col='cyan',xlim=c(0.5,7.5),ylim=c(0,25),main='',freq=F,xaxt='n')
  
  
  coldper <- c(1783,1809,1815,1831,1835)
  for (syr in coldper) {
    eyr <- syr+3
    # change to precip in % 
    # load monthly 70yr climatology data
    if (syr < 1636) {
      load(file=paste0(echclimpath,'echam_clim_1636-1637_2ndgrid.Rdata'))
    } else {
      load(file=paste0(echclimpath,'echam_clim_',syr+3,'-',(syr+4),'_2ndgrid.Rdata'))  
    }
    pos <- which(echam_clim$names=="precip")[1:4608]
    # seasonal mean
    echmeanclim <- cbind(rowMeans(echam_clim$ensmean[,10:15]),
                         rowMeans(echam_clim$ensmean[,16:21]))
    for (cyr in seq(syr,eyr)) {
      load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
      # change to precip in % 
      echam.abs <- echam
      echam.anom$ensmean[pos,] <- echam.abs$ensmean[pos,]/(echmeanclim[pos,]*3600*24*30)*100
      analysis.abs <- analysis
      analysis.anom$ensmean[pos,] <- analysis.abs$ensmean[pos,]/
        (echmeanclim[pos,]*3600*24*30)*100
      if (cyr==syr) {
        anameananom <- analysis.anom$ensmean
      } else {
        anameananom <- cbind(anameananom,analysis.anom$ensmean)
      }
    }
    
    # plot post volc years
    pdata <- echam
    ti=1
    pdf(paste0('../figures/laki/post_erup_anom_hist_',syr,'.pdf'), width=16, 
        height=12, paper='special')
    #layout(matrix(c(1,2,3,4,5,5,5,5,6,7,8,9,10,10,10,10,11,12,13,14,15,15,15,15,
    #                16,17,18,19,20,20,20,20),8,4,byrow=TRUE),height=c(3,1,3,1,3,1,3,1))
    layout(matrix(c(1,2,3,4,5,5,5,5,6,7,8,9,10,11,12,13,14,14,14,14,15,
                    16,17,18),6,4,byrow=TRUE),height=c(3,1,3,3,1,3))
    par(oma=c(0.5,4,4,0))
    levs <- c(-Inf, seq(-2,2,0.5), Inf)
    contlevs <- seq(-4,4,2)
    scfac <- max(echam.anom$ensmean[echam$names=='u200'])
    pdata$data <- array(c(anameananom[,ti],anameananom[,ti+2],anameananom[,ti+4],
                          anameananom[,ti+6]), c(nrow(echam.anom$ensmean),1,4))
    plot_echam3(pdata, varname='temp2', type='data', cex.pt=1.7, 
                latlim=c(0,90),lonlim=c(-180,180),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('e)', 'f)', 'g)', 'h)'),colnames=seq(syr,eyr),
                statpanel=NULL, add=T, rownames='T2m, Z500, UV200 (Oct-Mar)', 
                main='EKF400', units='K',
                addcontours=T, contvarname='gph500', conttype='data',
                contcol='black', contlev=contlevs,
                addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                veclen=scfac*0.01, vecscale=scfac*2.0, vecwd=0.95, every_x_vec=4)
    
    # histograms
    for (cyr in seq(syr,eyr)) {
      wpos <- which(substr((wt[,1]),1,7)==paste((cyr-1),'10',sep='-')|
                     substr((wt[,1]),1,7)==paste((cyr-1),'11',sep='-')|
                     substr((wt[,1]),1,7)==paste((cyr-1),'12',sep='-')|
                     substr((wt[,1]),1,7)==paste(cyr,'01',sep='-')|
                     substr((wt[,1]),1,7)==paste(cyr,'02',sep='-')|
                     substr((wt[,1]),1,7)==paste(cyr,'03',sep='-'))
      wpos2 <- which(as.numeric(substr((wt[,1]),1,4))>=(cyr-35)&
                     as.numeric(substr((wt[,1]),1,4))<=(cyr+35)&
                     (as.numeric(substr((wt[,1]),6,7))>=10|
                     as.numeric(substr((wt[,1]),6,7))<=3))
      par(mar=c(2,2,1,1))
      h1 = hist(as.numeric(wt[wpos,2])-0.2,breaks=seq(0.66,7.33,0.33),plot=F)
      h1$density = h1$counts/sum(h1$counts)*100
      #plot(h,col='blue',xlim=c(0.5,7.5),ylim=c(0,25),main='CAP7',freq=F,xaxt='n')
      #axis(side=1, at=seq(1,7), labels=wtleg)
      h2 = hist(as.numeric(wt[wpos2,2])-0.2,breaks=seq(0.66,7.33,0.33),plot=F)
      h2$density = h2$counts/sum(h2$counts)*100
      h <- h1
      h$density <- h1$density-h2$density
      #par(new=T)
      plot(h,col='blue',xlim=c(0.5,7.5),ylim=c(-10,10),main='',freq=F,xaxt='n')
      axis(side=1, at=seq(1,7), labels=wtleg)
      h3 = hist(as.numeric(wtnt[wpos,2])+0.2,breaks=seq(0.66,7.33,0.33),plot=F)
      h3$density = h3$counts/sum(h3$counts)*100
      #par(new=T)
      #plot(h3,col='red',xlim=c(0.5,7.5),ylim=c(0,25),main=' ',freq=F,xaxt='n')
      h4 = hist(as.numeric(wtnt[wpos2,2])+0.2,breaks=seq(0.66,7.33,0.33),plot=F)
      h4$density = h4$counts/sum(h4$counts)*100
      h5 <- h3
      h5$density <- h3$density-h4$density
      par(new=T)
      plot(h5,col='cyan',xlim=c(0.5,7.5),ylim=c(-10,10),main='',freq=F,xaxt='n')
    }
    
    ti=ti+1 # Apr-Sep season
    levs <- c(-Inf, seq(-2,2,0.5), Inf)
    pdata$data <- array(c(anameananom[,ti],anameananom[,ti+2],anameananom[,ti+4],
                          anameananom[,ti+6]), c(nrow(echam.anom$ensmean),1,4))
    plot_echam3(pdata, varname='temp2', type='data', cex.pt=1.7, 
                latlim=c(0,90),lonlim=c(-180,180),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('e)', 'f)', 'g)', 'h)'),colnames=rep(' ',4),
                statpanel=NULL, add=T, rownames='T2m, Z500, UV200  (Apr-Sep)', 
                main=' ', units='K',
                addcontours=T, contvarname='gph500', conttype='data',
                contcol='black', contlev=contlevs,
                addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                veclen=scfac*0.01, vecscale=scfac*2.0, vecwd=0.95, every_x_vec=4)
    
    for (cyr in seq(syr,eyr)) {
      spos <- which(substr((wt[,1]),1,7)==paste(cyr,'04',sep='-')|
                      substr((wt[,1]),1,7)==paste(cyr,'05',sep='-')|
                      substr((wt[,1]),1,7)==paste(cyr,'06',sep='-')|
                      substr((wt[,1]),1,7)==paste(cyr,'07',sep='-')|
                      substr((wt[,1]),1,7)==paste(cyr,'08',sep='-')|
                      substr((wt[,1]),1,7)==paste(cyr,'09',sep='-'))
      spos2 <- which(as.numeric(substr((wt[,1]),1,4))>=(cyr-35)&
                       as.numeric(substr((wt[,1]),1,4))<=(cyr+35)&
                       (as.numeric(substr((wt[,1]),6,7))>=4&
                          as.numeric(substr((wt[,1]),6,7))<=9))
      par(mar=c(2,2,1,1))
      h1 = hist(as.numeric(wt[spos,2])-0.2,breaks=seq(0.66,7.33,0.33),plot=F)
      h1$density = h1$counts/sum(h1$counts)*100
      #plot(h,col='blue',xlim=c(0.5,7.5),ylim=c(0,25),main='CAP7',freq=F,xaxt='n')
      #axis(side=1, at=seq(1,7), labels=wtleg)
      h2 = hist(as.numeric(wt[spos2,2])-0.2,breaks=seq(0.66,7.33,0.33),plot=F)
      h2$density = h2$counts/sum(h2$counts)*100
      h <- h1
      h$density <- h1$density-h2$density
      #par(new=T)
      plot(h,col='blue',xlim=c(0.5,7.5),ylim=c(-10,10),main='',freq=F,xaxt='n')
      axis(side=1, at=seq(1,7), labels=wtleg)
      h3 = hist(as.numeric(wtnt[spos,2])+0.2,breaks=seq(0.66,7.33,0.33),plot=F)
      h3$density = h3$counts/sum(h3$counts)*100
      #par(new=T)
      #plot(h3,col='red',xlim=c(0.5,7.5),ylim=c(0,25),main=' ',freq=F,xaxt='n')
      h4 = hist(as.numeric(wtnt[spos2,2])+0.2,breaks=seq(0.66,7.33,0.33),plot=F)
      h4$density = h4$counts/sum(h4$counts)*100
      h5 <- h3
      h5$density <- h3$density-h4$density
      par(new=T)
      plot(h5,col='cyan',xlim=c(0.5,7.5),ylim=c(-10,10),main='',freq=F,xaxt='n')
    }
    dev.off()  
  } # end periods loop
  
  
  
  
  # monthly analysis Jun 1783-Mai
  # 12x3 matrix: 12 months x 3 temp/gph/wind, prec/slp, w-type hist
  syr=1783
  eyr=1784
  # mik's weathertype recon
  wt <- read.table('../comparison_data/CAP7_1763-2009.txt',stringsAsFactors=F,head=T)
  wtnt <- read.table('../comparison_data/CAP7_1763-2009_NoTemp.txt',stringsAsFactors=F,head=T)
  wtleg=c('NE','WSW','W','E','HP','N','WC')
  # read clim for anom
  if (syr < 1636) {
    load(file=paste0(echclimpath,'echam_clim_1636-1637_2ndgrid.Rdata'))
  } else {
    load(file=paste0(echclimpath,'echam_clim_',syr+0,'-',(syr+1),'_2ndgrid.Rdata'))  
  }
  posp <- which(echam_clim$names=="precip")[1:4608]
  post <- which(echam_clim$names=="temp2")[1:4608]
  # select Oct to Sep, NO Jun to May
  echmeanclim <- echam_clim$ensmean[,10:21]
  for (cyr in seq(syr,eyr)) {
    load(file=paste0(prepplotdirmon,'analysis_',cyr,'.Rdata'))
    # strangely monthly data is still in sixmonstatevector???
    #analysis.anom$ensmean <- array(analysis.anom$ensmean,c((dim(analysis.anom$ensmean)[1]/6),
    #                                                       dim(analysis.anom$ensmean)[2]*6))
    # change to precip in % 
    analysis.abs <- analysis
    echam.abs <- echam
    analysis.abs$ensmean[posp,] <- analysis.abs$ensmean[posp,]/
                                     (echmeanclim[posp,]*3600*24*30)*100
    echam.abs$ensmean[posp,] <- echam.abs$ensmean[posp,]/
                                     (echmeanclim[posp,]*3600*24*30)*100
    # change to abs temp to anom
    analysis.abs$ensmean[post,] <- analysis.abs$ensmean[post,]-(echmeanclim[post,]-273.15)
    echam.abs$ensmean[post,] <- echam.abs$ensmean[post,]-(echmeanclim[post,]-273.15)
    if (cyr==syr) {
      anameananom <- analysis.abs$ensmean
      echmeananom <- echam.abs$ensmean
    } else {
      anameananom <- cbind(anameananom,analysis.abs$ensmean)
      echmeananom <- cbind(echmeananom,echam.abs$ensmean)
    }
  }
  anameananom <- anameananom[,6:17] #10:21]
  echmeananom <- echmeananom[,6:17]
  # plot post volc years
  pdata <- echam
  pdf(paste0('../figures/laki/laki_monthly_CCC400.pdf'), width=36, 
      height=8, paper='special')
  #layout(matrix(c(1,2,3,4,5,5,5,5,6,7,8,9,10,10,10,10,11,12,13,14,15,15,15,15,
  #                16,17,18,19,20,20,20,20),8,4,byrow=TRUE),height=c(3,1,3,1,3,1,3,1))
  layout(matrix(c(c(seq(1,12),rep(13,12)),c(seq(14,25),rep(26,12)),seq(27,38)),5,12,byrow=T),
         height=c(3,1,3,1,3))
  par(oma=c(0.5,4,4,0))
  levs <- c(-Inf, seq(-2,2,0.5), Inf)
  contlevs <- seq(4500,5900,100)       #-4,4,2)
  scfac <- max(echam.anom$ensmean[echam$names=='u200'])
  #pdata$data <- array(c(anameananom[,ti],anameananom[,ti+2],anameananom[,ti+4],
  #                      anameananom[,ti+6]), c(nrow(echam.anom$ensmean),1,4))
  pdata$data <- array(anameananom, c(nrow(anameananom),1,12))
  plot_echam3(pdata, varname='temp2', type='data', cex.pt=1.7, 
              latlim=c(0,90),lonlim=c(-180,180),
              wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
              names='', #c('e)', 'f)', 'g)', 'h)'),
              colnames=c('Jun','Jul','Aug','Sep','Oct','Nov','Dec',
                         'Jan','Feb','Mar','Apr','May'), #seq(syr,eyr),
              statpanel=NULL, add=T, rownames='T2m anom,Z500 abs,UV200 abs', 
              main='EKF400 1783-1784', units='K',
              addcontours=T, contvarname='gph500', conttype='data',
              contcol='black', contlev=contlevs,
              addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
              veclen=scfac*0.01, vecscale=scfac*0.2, 
              vecwd=0.95, every_x_vec=4)
  # levs=c(0,80,85,90,95,100,105,110,115,120,Inf)
  # contlevs <- seq(980,1030,5) #-4,4,2)
  # plot_echam3(pdata, varname='precip', type='data', cex.pt=1.7, 
  #             latlim=c(0,90),lonlim=c(-180,180),
  #             wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
  #             names='', #c('e)', 'f)', 'g)', 'h)'),
  #             colnames=rep(' ',12), #seq(syr,eyr),
  #             statpanel=NULL, add=T, rownames='Prec anom, SLP abs', 
  #             main='', units='%',
  #             addcontours=T, contvarname='slp', conttype='data',
  #             contcol='black', contlev=contlevs) #,
  #             #addvectors=F, vecnames=c('u200','v200'), veccol='cyan', 
  #             #veclen=scfac*0.01, vecscale=scfac*2.0, 
  #             #vecwd=0.95, every_x_vec=4)
  pdata$data <- array(echmeananom, c(nrow(echmeananom),1,12))
  plot_echam3(pdata, varname='temp2', type='data', cex.pt=1.7, 
              latlim=c(0,90),lonlim=c(-180,180),
              wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
              names='', #c('e)', 'f)', 'g)', 'h)'),
              colnames=rep(' ',12), #c('Jun','Jul','Aug','Sep','Oct','Nov','Dec',
              #           'Jan','Feb','Mar','Apr','May'), #seq(syr,eyr),
              statpanel=NULL, add=T, rownames='same for CCC400', 
              main='', units='K',
              addcontours=T, contvarname='gph500', conttype='data',
              contcol='black', contlev=contlevs,
              addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
              veclen=scfac*0.01, vecscale=scfac*0.2, 
              vecwd=0.95, every_x_vec=4)
  
  # histograms
  for (mon in c(seq(6,12),seq(1,5))) { #seq(syr,eyr)) {
    if (mon > 5) {yr=syr} else {yr=syr+1}
    if (mon < 10) {mon=paste0('0',mon)}
    #print(c(yr,mon))
    pos1 <- which(substr((wt[,1]),1,7)==paste(yr,mon,sep='-'))
    pos2 <- which(as.numeric(substr((wt[,1]),1,4))>=(cyr-35)&
                     as.numeric(substr((wt[,1]),1,4))<=(cyr+35)&
                     (as.numeric(substr((wt[,1]),6,7))==as.numeric(mon)))
    par(mar=c(2,2,1,1))
    h1 = hist(as.numeric(wt[pos1,2])-0.2,breaks=seq(0.66,7.33,0.33),plot=F)
    h1$density = h1$counts/sum(h1$counts)*100
    h2 = hist(as.numeric(wt[pos2,2])-0.2,breaks=seq(0.66,7.33,0.33),plot=F)
    h2$density = h2$counts/sum(h2$counts)*100
    h <- h1
    h$density <- h1$density-h2$density
    plot(h,col='blue',xlim=c(0.5,7.5),ylim=c(-35,35),main='',freq=F,xaxt='n')
    axis(side=1, at=seq(1,7), labels=wtleg)
    h3 = hist(as.numeric(wtnt[pos1,2])+0.2,breaks=seq(0.66,7.33,0.33),plot=F)
    h3$density = h3$counts/sum(h3$counts)*100
    h4 = hist(as.numeric(wtnt[pos2,2])+0.2,breaks=seq(0.66,7.33,0.33),plot=F)
    h4$density = h4$counts/sum(h4$counts)*100
    h5 <- h3
    h5$density <- h3$density-h4$density
    par(new=T)
    plot(h5,col='cyan',xlim=c(0.5,7.5),ylim=c(-35,35),main='',freq=F,xaxt='n')
  }
  dev.off()
  
  
  
  # plot ccc400 crowley forcing and GISS aerosols
  vol <- read.table("/Volumes/data/climdata/forcings/echam5_input/crowley_volcanic_forcing.txt",
                    header=T)
  vol[35389:35424,]
  volENH_ts_nh <- ts(vol[,2],start=vol[1,1],freq=36)
  volENH_ts_sh <- ts(vol[,8],start=vol[1,1],freq=36)
  # troposperic aerosols do not have Laki aerosols but just natural monthly climatology
  # plus antropogenic aerosols scaled to population density before 1875
  # for (yr in seq(1783,1786)) {
  # nc=nc_open(paste0(
  #   '/Volumes/DATA/climdata/forcings/echam_ccc400_forcings/GISS_Aerosols/all_sulfate_',
  #   yr,'.nc'), write=F)
  #   # print(nc)
  #   stmp <- ncvar_get(nc, "sulfate")
  #   if (yr==1783){
  #     sulf=stmp
  #   } else {
  #     sulf=c(sulf,stmp)    
  #   }
  # nc_close(nc)
  pdf('../figures/laki/aod_forc.pdf',width=15,height=4)
    par(mfrow=c(1,5),mar=c(4,4,4,1)) #,oma=c(2,2,2,0.2))
    plot(window(volENH_ts_nh,1783,1786),col='red',xlab='Year',ylab='AOD',
         main='Laki',ylim=c(0,0.45))
    lines(window(volENH_ts_sh,1783,1786),col='blue')
    legend('topright',lty=1,col=c('red','blue'),c('AOD 30-90N','AOD 30-90S'))
    plot(window(volENH_ts_nh,1808.8,1812),col='red',xlab='Year',ylab='AOD',
         main='Unknown',ylim=c(0,0.45))
    lines(window(volENH_ts_sh,1808.8,1812),col='blue')
    plot(window(volENH_ts_nh,1815,1818),col='red',xlab='Year',ylab='AOD',
         main='Tambora',ylim=c(0,0.45))
    lines(window(volENH_ts_sh,1815,1818),col='blue')
    plot(window(volENH_ts_nh,1831,1834),col='red',xlab='Year',ylab='AOD',
         main='Babuyan',ylim=c(0,0.45))
    lines(window(volENH_ts_sh,1831,1834),col='blue')
    plot(window(volENH_ts_nh,1835,1838),col='red',xlab='Year',ylab='AOD',
         main='Cosigüina',ylim=c(0,0.45))
    lines(window(volENH_ts_sh,1835,1838),col='blue')
  dev.off()
  #pdf('../figures/laki/tambora_aod_forc.pdf',width=5,height=4)
  #  par(mfrow=c(1,1),mar=c(4,4,4,1)) #,oma=c(2,2,2,0.2))
  #  plot(window(volENH_ts,1815,1817),col='red',xlab='Year',ylab='AOD',
  #     main='Crowley volcanic forcing 30-90N')
  #dev.off()
  
  
  
  
  
  
  # calc min temp and gph500 at 50, 60, 70ºN summer/winter landocean ens mem and ens mean
  #filenameext <- paste0('_anom_landocean_seas_')
  for (cyr in syr:eyr) {
    print(cyr)
    load(file=paste0(prepplotdirseas,'analysis_',cyr,'.Rdata'))
    post <- which(floor(analysis.anom$lat)==62 & analysis.anom$names=='temp2') 
    posg <- which(floor(analysis.anom$lat)==62 & analysis.anom$names=='gph500') 
    # other possible lat==51.294377    #lat==69.946081
    if (cyr==syr) {
      tminmean <- apply(analysis.anom$ensmean[post,],2,min)
      tminmem <- apply(analysis.anom$data[post,,],2:3,min)
      gminmean <- apply(analysis.anom$ensmean[posg,],2,min)
      gminmem <- apply(analysis.anom$data[posg,,],2:3,min)
    } else {
      tminmean <- rbind(tminmean,apply(analysis.anom$ensmean[post,],2,min))
      tminmem <- abind(tminmem,apply(analysis.anom$data[post,,],2:3,min),along=3)
      gminmean <- rbind(gminmean,apply(analysis.anom$ensmean[posg,],2,min))
      gminmem <- abind(gminmem,apply(analysis.anom$data[posg,,],2:3,min),along=3)
    }
  }  
  # plot minima
  vf <- c(1600,1650,1660,1783,1809,1815,1835,1883,1886,1902,1912,1980,1991)
  year <- syr:eyr

  pdf('../figures/laki/t2m_z500_min_ts_seas.pdf',width=8,height=12)
    par(mfrow=c(4,1),mar=c(4,4,4,1)) #,oma=c(2,2,2,0.2))
    plot(year,tminmem[1,1,],ty='l',col=rgb(10,5,0,1,maxColorValue=10),
      xlab='Year',ylab='Temp. [ºC]', main='62ºN Zonal Minima of Winter Temperature',
      ylim=c(min(tminmem[1,,]),max(tminmem[1,,])))
    for (i in 2:30){
      lines(year,tminmem[1,i,],col=rgb(10,5,0,1,maxColorValue=10))
    }
    lines(year,tminmean[,1],col=rgb(10,0,0,10,maxColorValue=10))
    abline(v=vf,col=rgb(7,7,7,5,maxColorValue=10))
    legend('bottomright',lty=1,col=c('red','orange'),c('Ensemble mean','Ensemble members'), 
           bty='o', bg='white',box.col='white')

    plot(year,tminmem[2,1,],ty='l',col=rgb(10,5,0,1,maxColorValue=10),
         xlab='Year',ylab='Temperature [ºC]', main='62ºN Zonal Minima of Summer Temperature',
         ylim=c(min(tminmem[2,,]),max(tminmem[2,,])))
    for (i in 2:30){
      lines(year,tminmem[2,i,],col=rgb(10,5,0,1,maxColorValue=10))
    }
    lines(year,tminmean[,2],col=rgb(10,0,0,10,maxColorValue=10))
    abline(v=vf,col=rgb(7,7,7,5,maxColorValue=10))    
  
    plot(year,gminmem[1,1,],ty='l',col=rgb(0,10,10,1,maxColorValue=10),
         xlab='Year',ylab='500 hPa Geopotential Height [m]', main='62ºN Zonal Minima of Winter GPH500',
         ylim=c(min(gminmem[1,,]),max(gminmem[1,,])))
    for (i in 2:30){
      lines(year,gminmem[1,i,],col=rgb(0,10,10,1,maxColorValue=10))
    }
    lines(year,gminmean[,1],col=rgb(0,0,10,10,maxColorValue=10))
    abline(v=vf,col=rgb(7,7,7,5,maxColorValue=10))
    legend('bottomright',lty=1,col=c('blue','cyan'),c('Ensemble mean','Ensemble members'), 
           bty='o', bg='white',box.col='white')
    
    plot(year,gminmem[2,1,],ty='l',col=rgb(0,10,10,1,maxColorValue=10),
         xlab='Year',ylab='500 hPa Geopotential Height [m]', main='62ºN Zonal Minima of Summer GPH500',
         ylim=c(min(gminmem[2,,]),max(gminmem[2,,])))
    for (i in 2:30){
      lines(year,gminmem[2,i,],col=rgb(0,10,10,1,maxColorValue=10))
    }
    lines(year,gminmean[,2],col=rgb(0,0,10,10,maxColorValue=10))
    abline(v=vf,col=rgb(7,7,7,5,maxColorValue=10))   
  dev.off()
  
  
  # same but based on one monthly minimum instead of seasonal averages 
  # but still summer and winter separated
  for (cyr in syr:eyr) {
    print(cyr)
    load(file=paste0(prepplotdirmon,'analysis_',cyr,'.Rdata'))
    post <- which(floor(analysis.anom$lat)==62 & analysis.anom$names=='temp2') 
    posg <- which(floor(analysis.anom$lat)==62 & analysis.anom$names=='gph500') 
    # other possible lat==51.294377    #lat==69.946081
    if (cyr==syr) {
      tminmean <- apply(analysis.anom$ensmean[post,],2,min)
      tminmem <- apply(analysis.anom$data[post,,],2:3,min)
      gminmean <- apply(analysis.anom$ensmean[posg,],2,min)
      gminmem <- apply(analysis.anom$data[posg,,],2:3,min)
    } else {
      tminmean <- rbind(tminmean,apply(analysis.anom$ensmean[post,],2,min))
      tminmem <- abind(tminmem,apply(analysis.anom$data[post,,],2:3,min),along=3)
      gminmean <- rbind(gminmean,apply(analysis.anom$ensmean[posg,],2,min))
      gminmem <- abind(gminmem,apply(analysis.anom$data[posg,,],2:3,min),along=3)
    }
  }  
  # plot minima
  vf <- c(1600,1650,1660,1783,1809,1815,1835,1883,1886,1902,1912,1980,1991)
  year <- syr:eyr
  
  pdf('../figures/laki/t2m_z500_min_ts_mon.pdf',width=8,height=12)
  par(mfrow=c(4,1),mar=c(4,4,4,1)) #,oma=c(2,2,2,0.2))
  plot(year,tminmem[1,1,],ty='l',col=rgb(10,5,0,1,maxColorValue=10),
       xlab='Year',ylab='Temperature [ºC]', main='62ºN Zonal Minima of Winter Temperature',
       ylim=c(min(tminmem[1,,]),max(tminmem[1,,])))
  for (i in 2:30){
    lines(year,tminmem[1,i,],col=rgb(10,5,0,1,maxColorValue=10))
  }
  lines(year,tminmean[,1],col=rgb(10,0,0,10,maxColorValue=10))
  abline(v=vf,col=rgb(7,7,7,5,maxColorValue=10))
  legend('bottomright',lty=1,col=c('red','orange'),c('Ensemble mean','Ensemble members'), 
         bty='o', bg='white',box.col='white')
  
  plot(year,tminmem[2,1,],ty='l',col=rgb(10,5,0,1,maxColorValue=10),
       xlab='Year',ylab='Temperature [ºC]', main='62ºN Zonal Minima of Summer Temperature',
       ylim=c(min(tminmem[2,,]),max(tminmem[2,,])))
  for (i in 2:30){
    lines(year,tminmem[2,i,],col=rgb(10,5,0,1,maxColorValue=10))
  }
  lines(year,tminmean[,2],col=rgb(10,0,0,10,maxColorValue=10))
  abline(v=vf,col=rgb(7,7,7,5,maxColorValue=10))    
  
  plot(year,gminmem[1,1,],ty='l',col=rgb(0,10,10,1,maxColorValue=10),
       xlab='Year',ylab='500 hPa Geopotential Height [m]', main='62ºN Zonal Minima of Winter GPH500',
       ylim=c(min(gminmem[1,,]),max(gminmem[1,,])))
  for (i in 2:30){
    lines(year,gminmem[1,i,],col=rgb(0,10,10,1,maxColorValue=10))
  }
  lines(year,gminmean[,1],col=rgb(0,0,10,10,maxColorValue=10))
  abline(v=vf,col=rgb(7,7,7,5,maxColorValue=10))
  legend('bottomright',lty=1,col=c('blue','cyan'),c('Ensemble mean','Ensemble members'), 
         bty='o', bg='white',box.col='white')
  
  plot(year,gminmem[2,1,],ty='l',col=rgb(0,10,10,1,maxColorValue=10),
       xlab='Year',ylab='500 hPa Geopotential Height [m]', main='62ºN Zonal Minima of Summer GPH500',
       ylim=c(min(gminmem[2,,]),max(gminmem[2,,])))
  for (i in 2:30){
    lines(year,gminmem[2,i,],col=rgb(0,10,10,1,maxColorValue=10))
  }
  lines(year,gminmean[,2],col=rgb(0,0,10,10,maxColorValue=10))
  abline(v=vf,col=rgb(7,7,7,5,maxColorValue=10))   
  dev.off()
  
} end laki  
  
  







if (iceland) {
  # load pre calc time series
  expname="EKF400_v1.3_corr_echam_clim" #"EKF400_v1.1_correct_71yr_anom"
  filenameext = '_anom_landonly_seas_'
  load(paste0('../data/indices/',expname,'/indices',filenameext,'allts.Rdata'))
  
  # define plot parameters
  wlen=1
  linew=1
  scal=F
  plbem=F
  syr=1604
  eyr=2003
  erup_yrs=c(1612,1625,1660,1693,1721,1727,1755,1766,1783,1845,1860,
             1873,1875,1902,1918,1947) #,2010,2011)
  # erup mon >=2 erup yr + 1  
  erup_win=c(1613,1626,1661,1694,1722,1728,1756,1767,1784,1846,1861,
             1874,1875,1903,1919,1948)
  erup_win_nolaki=c(1613,1626,1661,1694,1722,1728,1756,1767,1846,1861,
             1874,1875,1903,1919,1948)
  # erup mon >=8 erup yr + 1  
  erup_sum=c(1613,1625,1661,1693,1721,1728,1756,1766,1783,1845,1861,
             1874,1876,1903,1919,1947)
  erup_sum_nolaki=c(1613,1625,1661,1693,1721,1728,1756,1766,1783,1845,1861,
             1874,1876,1903,1919,1947)
  # strong eruption (5)
  erup_win2=c(1626,1722,1756,1874)
  # erup mon >=8 erup yr + 1  
  erup_sum2=c(1625,1721,1756,1874)
#   Katla             4	1612	10	12
#   Katla	            5	1625	9	  2
#   Katla	            4	1660	11	3
#   Hekla	            4	1693	2	  13
#   Katla	            5	1721	5	  11
#   Oraefajokull	    4	1727	8 	3
#   Katla	            5	1755	10	17
#   Hekla	            4	1766	4	  5
#   Grimsvotn	        4	1783	5  	0
#   Hekla	            4	1845	9	  2
#   Katla	            4	1860	5	  8
#   Grimsvotn	        4	1873	1	  8
#   Askja	            5	1875	1  	1
#   Grimsvotn        	4	1902	12 	0
#   Katla	            4	1918	10	12
#   Hekla	            4	1947	3	  29
#   Eyjafjallajokull	4	2010	3	  20
#   Grimsvotn	        4	2011	5	  21

  # seas='sum'
  # reg=which(aind.allts$names=='GLO.temp2')
  pdf('../figures/iceland/t2m_ts_1600-2000.pdf',width=8,height=6)
  par(mfrow=c(2,1),mar=c(2,2,2,0.2))
  for (seas in c('sum','win')) {
    i=1
    for (reg in c(which(aind.allts$names=='ENH.temp2'),which(aind.allts$names=='NAM.temp2'),
                  which(aind.allts$names=='EU.temp2'))) {
      print(reg)
      if (seas=='yrmean') {
        years <- rep(eind.allts$time,each=2)  
      } else {
        wincol <- is.odd(seq(1,ncol(aind.allts$ensmean)))
        somcol <- is.even(seq(1,ncol(aind.allts$ensmean)))  
      }
      if (seas=='yrmean') {
        print("annual mean")
        #yl=c(0,15)
        yl=c(-1.5,1.5)
        eindmean <- aggregate(eind.allts$ensmean[reg,],list(years),mean,na.rm=T)[,2]
        aindmean <- aggregate(aind.allts$ensmean[reg,],list(years),mean,na.rm=T)[,2]
      } else if (seas=='sum'){ 
        #yl=c(13,16.5)
        yl=c(-1.5,1.5)
        print("summer season")
        eindmean <- eind.allts$ensmean[reg,somcol]  
        aindmean <- aind.allts$ensmean[reg,somcol]  
      } else if (seas=='win'){
        print("winter season")
        #yl=c(-5,4)
        yl=c(-2,2)
        eindmean <- eind.allts$ensmean[reg,wincol]  
        aindmean <- aind.allts$ensmean[reg,wincol]  
      }
      # northern hemisphere extra-tropics (>20N) seas="yrmean" or "sum" or "win"
      yrpos=which(aind.allts$time>syr&aind.allts$time<eyr)
      yrs=aind.allts$time[yrpos]
      if (i==1) {
        plot(yrs,aindmean[yrpos],ty="l",ylim=yl,main=seas)
        abline(v=erup_yrs,col='grey')
        abline(h=mean(aindmean[yrpos]),col=i,lty=2)
      } else {
        lines(yrs,aindmean[yrpos],ty="l",col=i)
        abline(h=mean(aindmean[yrpos]),col=i,lty=2)
      }
      legend("topleft",c('ENH','NAM','EU'),col=c(1,2,4),lty=1,bg='white') #,box.col='white')
      i=i+1
      if (i==3) {i=i+1}
    } # end reg loop for various regions
  } # end sum win season loop
  dev.off()
  
  
  
  
  
  
  # spatial maps
  # composite all eruptions
  # separated summer and winter to really have season after eruption
  for (s in c('som','win','som_nolaki','win_nolaki')) {
    if (s=='som') {yrs <- erup_sum}
    if (s=='win') {yrs <- erup_win}
    if (s=='som_nolaki') {yrs <- erup_sum_nolaki}
    if (s=='win_nolaki') {yrs <- erup_win_nolaki}
    i=1
    for (cyr in yrs) {
    # change to precip in % 
    # load monthly 70yr climatology data
      if (cyr < 1636) {
        load(file=paste0(echclimpath,'echam_clim_1636-1637_2ndgrid.Rdata'))
      } else {
        load(file=paste0(echclimpath,'echam_clim_',cyr,'-',(cyr+1),'_2ndgrid.Rdata'))  
      }
      pos <- which(echam_clim$names=="precip")[1:4608]
      # seasonal mean
      echmeanclim <- cbind(rowMeans(echam_clim$ensmean[,10:15]),
                           rowMeans(echam_clim$ensmean[,16:21]))
      load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
      # change to precip in % 
      echam.abs <- echam
      echam.anom$ensmean[pos,] <- (echam.abs$ensmean[pos,]/echmeanclim[pos,])-1 
        #*3600*24*30)*100
      analysis.abs <- analysis
      analysis.anom$ensmean[pos,] <- (analysis.abs$ensmean[pos,]/echmeanclim[pos,])-1 
        #*3600*24*30)*100
      if (i==1) {
        anameananom <- analysis.anom$ensmean
      } else {
        anameananom <- cbind(anameananom,analysis.anom$ensmean)
      }
      i=i+1
    }
    if (s=='som') {pos2 <- seq(2,dim(anameananom)[2],2)}
    if (s=='win') {pos2 <- seq(1,dim(anameananom)[2],2)}
    if (s=='som_nolaki') {pos2 <- seq(2,dim(anameananom)[2],2)}
    if (s=='win_nolaki') {pos2 <- seq(1,dim(anameananom)[2],2)}
    for (v in c('temp2','precip','slp')) {
      varmem <- anameananom[which(analysis.anom$names==v),pos2]
      varmean <- apply(anameananom[which(analysis.anom$names==v),pos2],1,mean)
      # test significance
      pmat <- rep(NA,nrow(varmem))
      for (i in 1:nrow(varmem)){
        pmat[i] <- unlist(t.test(varmem[i,],alternative="two.sided",mu=0,
                                 paired=F,var.equal=F,conf.level=0.95)[3])
      }
      pmat <- replace(pmat,pmat>=0.05,NA)
      pmat <- replace(pmat,pmat<0.05,1)
      if (v=='temp2') {
        lev <- seq(-2.2,2.2,0.4) 
        colorTable<- designer.colors((length(lev)-1), c( "darkblue","white", "darkred"))
      } else if (v=='slp') {  
        colorTable<- designer.colors(11, c( "blue","white", "darkorange"))
        lev <- c(-1.5,-1.0,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1.0,1.5) 
      } else if (v=='precip') {  
        colorTable<- designer.colors(11, c( "brown","white", "darkgreen"))
        #lev <- c(-1.5,-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9,1.5) 
        lev <- c(-0.5,-0.2,-0.1,-0.075,-0.05,-0.025,0.025,0.05,0.075,0.1,0.2,0.5) 
      }
      lon<-unique(analysis.anom$lon[which(analysis.anom$names==v)])
      lat<-unique(analysis.anom$lat[which(analysis.anom$names==v)])
      varmeanmat <- matrix(varmean,nrow=length(lon),ncol=length(lat))
      varmemmat <- array(varmem,dim=c(length(lon),length(lat),dim(varmem)[2]))
      varmeanmat <- varmeanmat[c(which(lon<0),which(lon>=0)),]
      varmemmat <- varmemmat[c(which(lon<0),which(lon>=0)),,]
      lon <- c(lon[which(lon<0)],lon[which(lon>=0)])
      if (lat[1]>lat[2]) {
        lat <- rev(lat) # if latitude are in reverse order
        varmeanmat <- varmeanmat[,(dim(varmeanmat)[2]:1)]
        varmemmat <- varmemmat[,(dim(varmeanmat)[2]:1),]
      }
      pdf(paste0('../figures/iceland/post_erup_',s,'_',v,'.pdf'),width=8,height=6)
        image.plot(lon,lat,varmeanmat,col=colorTable,breaks=lev, #nlevel=11,zlim=c(-zr,zr),
                     main=paste('Post eruption',s,'mean',v)) #tim.colors())
        points(analysis.anom$lon[which(analysis.anom$names==v)],
                analysis.anom$lat[which(analysis.anom$names==v)],
                cex=pmat/30,pch=20,col="black")
        map("world",interior=F,col='lightgrey',fill=F,add=T)
        legend('bottomright','p<0.05',pch=20,col='black')
      dev.off()
    } # end variable loop
  } # end sum/win/nolaki loop



      
      
      
      
      
      
      
      
       
#  anameananomwin <- apply(anameananom[,seq(1,length(erup_win),2)],1,mean)

  i=1
  for (cyr in erup_sum) {
    # change to precip in % 
    # load monthly 70yr climatology data
    if (cyr < 1636) {
      load(file=paste0(echclimpath,'echam_clim_1636-1637_2ndgrid.Rdata'))
    } else {
      load(file=paste0(echclimpath,'echam_clim_',cyr,'-',(cyr+1),'_2ndgrid.Rdata'))  
    }
    pos <- which(echam_clim$names=="precip")[1:4608]
    # seasonal mean
    echmeanclim <- cbind(rowMeans(echam_clim$ensmean[,10:15]),
                         rowMeans(echam_clim$ensmean[,16:21]))
    load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
    # change to precip in % 
    echam.abs <- echam
    echam.anom$ensmean[pos,] <- echam.abs$ensmean[pos,]/(echmeanclim[pos,]*3600*24*30)*100
    analysis.abs <- analysis
    analysis.anom$ensmean[pos,] <- analysis.abs$ensmean[pos,]/
      (echmeanclim[pos,]*3600*24*30)*100
    if (i==1) {
      anameananom <- analysis.anom$ensmean
    } else {
      anameananom <- cbind(anameananom,analysis.anom$ensmean)
    }
    i=i+1
  }
  anameananomsum <- apply(anameananom[,seq(2,length(erup_sum),2)],1,mean)

  # mik's weathertype recon
  wt <- read.table('../comparison_data/CAP7_1763-2009.txt',stringsAsFactors=F,head=T)
  wtnt <- read.table('../comparison_data/CAP7_1763-2009_NoTemp.txt',stringsAsFactors=F,head=T)
  wtleg=c('NE','WSW','W','E','HP','N','WC')
  
  i=1
  for (cyr in erup_win) {
    wpos <- which(substr((wt[,1]),1,7)==paste((cyr-1),'10',sep='-')|
                  substr((wt[,1]),1,7)==paste((cyr-1),'11',sep='-')|
                  substr((wt[,1]),1,7)==paste((cyr-1),'12',sep='-')|
                  substr((wt[,1]),1,7)==paste(cyr,'01',sep='-')|
                  substr((wt[,1]),1,7)==paste(cyr,'02',sep='-')|
                  substr((wt[,1]),1,7)==paste(cyr,'03',sep='-'))
    wpos2 <- which(as.numeric(substr((wt[,1]),1,4))>=(cyr-35)&
                   as.numeric(substr((wt[,1]),1,4))<=(cyr+35)&
                   (as.numeric(substr((wt[,1]),6,7))>=10|
                   as.numeric(substr((wt[,1]),6,7))<=3))
    if (i==1) {
      hwin <- as.numeric(wt[wpos,2])
      hwinclim <- as.numeric(wt[wpos2,2])
    } else {
      hwin <- c(hwin,as.numeric(wt[wpos,2]))
      hwinclim <- c(hwinclim,as.numeric(wt[wpos2,2]))
    }
    i=i+1
  }
  i=1
  for (cyr in erup_sum) {
    wpos <- which(substr((wt[,1]),1,7)==paste((cyr),'04',sep='-')|
                  substr((wt[,1]),1,7)==paste((cyr),'05',sep='-')|
                  substr((wt[,1]),1,7)==paste((cyr),'06',sep='-')|
                  substr((wt[,1]),1,7)==paste(cyr,'07',sep='-')|
                  substr((wt[,1]),1,7)==paste(cyr,'08',sep='-')|
                  substr((wt[,1]),1,7)==paste(cyr,'09',sep='-'))
    wpos2 <- which(as.numeric(substr((wt[,1]),1,4))>=(cyr-35)&
                   as.numeric(substr((wt[,1]),1,4))<=(cyr+35)&
                   (as.numeric(substr((wt[,1]),6,7))>=10|
                      as.numeric(substr((wt[,1]),6,7))<=3))
    if (i==1) {
      hsum <- as.numeric(wt[wpos,2])
      hsumclim <- as.numeric(wt[wpos2,2])
    } else {
      hsum <- c(hsum,as.numeric(wt[wpos,2]))
      hsumclim <- c(hsumclim,as.numeric(wt[wpos2,2]))
    }
    i=i+1
  } 
  hw = hist(hwin-0.2,breaks=seq(0.66,7.33,0.33),plot=F)
  hw$density = hw$counts/sum(hw$counts)*100
  hwc = hist(hwinclim-0.2,breaks=seq(0.66,7.33,0.33),plot=F)
  hwc$density = hwc$counts/sum(hwc$counts)*100
  h <- hw
  h$density <- hw$density-hwc$density
  hs = hist(hsum+0.2,breaks=seq(0.66,7.33,0.33),plot=F)
  hs$density = hs$counts/sum(hs$counts)*100
  hsc = hist(hsumclim+0.2,breaks=seq(0.66,7.33,0.33),plot=F)
  hsc$density = hsc$counts/sum(hsc$counts)*100
  h2 <- hs
  h2$density <- hs$density-hsc$density

    
  # plot composites
  pdata <- echam
  pdf('../figures/iceland/post_iceland_erup_anom_composite.pdf', width=8, 
        height=8, paper='special')
    layout(matrix(c(1,2,3,3,4,5,6,6,7,7),5,2,byrow=TRUE),height=c(3,1,3,1,3))
    par(oma=c(0.5,4,4,0))
    levs <- c(-Inf, seq(-1,1,0.25), Inf)
    contlevs <- seq(-4,4,2)
    scfac <- max(echam.anom$ensmean[echam$names=='u200'])
    pdata$data <- array(c(anameananomwin,anameananomsum), c(nrow(echam.anom$ensmean),1,2))
    plot_echam3(pdata, varname='temp2', type='data', cex.pt=1.5, 
                latlim=c(0,90),lonlim=c(-180,180),
                wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
                names=c('a)', 'b)'),colnames=c('Oct-Mar','Apr-Sep'),
                statpanel=NULL, add=T, rownames='T2m, Z500, UV200', 
                main='EKF400', units='K',
                addcontours=T, contvarname='gph500', conttype='data',
                contcol='black', contlev=contlevs,
                addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
                veclen=scfac*0.01, vecscale=scfac*2.0, vecwd=0.95, every_x_vec=4)
    levs=c(0,80,85,90,95,100,105,110,115,120,Inf)
    plot_echam3(pdata, varname='precip', type='data', cex.pt=1.5, 
                latlim=c(0,90),lonlim=c(-180,180),units="%",
                names=c('c)', 'd)'), lev=levs, 
                st.col=NULL, stations=NULL, add=T,addcontours=F,
                wcol='darkgrey',rownames='Precipitation',colnames=rep('',4))
    # histograms
    par(mar=c(2,2,1,1))
    plot(h,col='blue',xlim=c(0.5,7.5),ylim=c(-20,20),main='',freq=F,xaxt='n')
    axis(side=1, at=seq(1,7), labels=wtleg)
    par(new=T)
    plot(h2,col='cyan',xlim=c(0.5,7.5),ylim=c(-20,20),main='',freq=F,xaxt='n')
    legend('topright',c('Winter','Summer'),lty=c(1,1),lwd=c(5,5),col=c('blue','cyan')) 
  dev.off()  
#  } # end periods loop
  





# WITHOUT Laki
# load pre calc time series
expname="EKF400_v1.1_correct_71yr_anom"
filenameext = '_anom_landonly_seas_'
load(paste0('../data/indices/',expname,'/indices',filenameext,'allts.Rdata'))

# define plot parameters
wlen=1
linew=1
scal=F
plbem=F
syr=1604
eyr=2003
erup_yrs=c(1612,1625,1660,1693,1721,1727,1755,1766,1845,1860,
           1873,1875,1902,1918,1947) #,2010,2011)
# erup mon >=2 erup yr + 1  
erup_win=c(1613,1626,1661,1694,1722,1728,1756,1767,1846,1861,
           1874,1875,1903,1919,1948)
# erup mon >=8 erup yr + 1  
erup_sum=c(1613,1625,1661,1693,1721,1728,1756,1766,1845,1861,
           1874,1876,1903,1919,1947)
# strong eruption (5)
erup_win2=c(1626,1722,1756,1874)
# erup mon >=8 erup yr + 1  
erup_sum2=c(1625,1721,1756,1874)
#   Katla             4  1612	10	12
#   Katla	            5	1625	9	  2
#   Katla	            4	1660	11	3
#   Hekla	            4	1693	2	  13
#   Katla	            5	1721	5	  11
#   Oraefajokull	    4	1727	8 	3
#   Katla	            5	1755	10	17
#   Hekla	            4	1766	4	  5
#   Grimsvotn	        4	1783	5  	0
#   Hekla	            4	1845	9	  2
#   Katla	            4	1860	5	  8
#   Grimsvotn	        4	1873	1	  8
#   Askja	            5	1875	1  	1
#   Grimsvotn        	4	1902	12 	0
#   Katla	            4	1918	10	12
#   Hekla	            4	1947	3	  29
#   Eyjafjallajokull	4	2010	3	  20
#   Grimsvotn	        4	2011	5	  21

# spatial maps
# composite all eruptions
# separated summer and winter to really have season after eruption
i=1
for (cyr in erup_win) {
  # change to precip in % 
  # load monthly 70yr climatology data
  if (cyr < 1636) {
    load(file=paste0(echclimpath,'echam_clim_1636-1637_2ndgrid.Rdata'))
  } else {
    load(file=paste0(echclimpath,'echam_clim_',cyr,'-',(cyr+1),'_2ndgrid.Rdata'))  
  }
  pos <- which(echam_clim$names=="precip")[1:4608]
  # seasonal mean
  echmeanclim <- cbind(rowMeans(echam_clim$ensmean[,10:15]),
                       rowMeans(echam_clim$ensmean[,16:21]))
  load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
  # change to precip in % 
  echam.abs <- echam
  echam.anom$ensmean[pos,] <- echam.abs$ensmean[pos,]/(echmeanclim[pos,]*3600*24*30)*100
  analysis.abs <- analysis
  analysis.anom$ensmean[pos,] <- analysis.abs$ensmean[pos,]/
    (echmeanclim[pos,]*3600*24*30)*100
  if (i==1) {
    anameananom <- analysis.anom$ensmean
  } else {
    anameananom <- cbind(anameananom,analysis.anom$ensmean)
  }
  i=i+1
}
anameananomwin <- apply(anameananom[,seq(1,length(erup_win),2)],1,mean)

i=1
for (cyr in erup_sum) {
  # change to precip in % 
  # load monthly 70yr climatology data
  if (cyr < 1636) {
    load(file=paste0(echclimpath,'echam_clim_1636-1637_2ndgrid.Rdata'))
  } else {
    load(file=paste0(echclimpath,'echam_clim_',cyr,'-',(cyr+1),'_2ndgrid.Rdata'))  
  }
  pos <- which(echam_clim$names=="precip")[1:4608]
  # seasonal mean
  echmeanclim <- cbind(rowMeans(echam_clim$ensmean[,10:15]),
                       rowMeans(echam_clim$ensmean[,16:21]))
  load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
  # change to precip in % 
  echam.abs <- echam
  echam.anom$ensmean[pos,] <- echam.abs$ensmean[pos,]/(echmeanclim[pos,]*3600*24*30)*100
  analysis.abs <- analysis
  analysis.anom$ensmean[pos,] <- analysis.abs$ensmean[pos,]/
    (echmeanclim[pos,]*3600*24*30)*100
  if (i==1) {
    anameananom <- analysis.anom$ensmean
  } else {
    anameananom <- cbind(anameananom,analysis.anom$ensmean)
  }
  i=i+1
}
anameananomsum <- apply(anameananom[,seq(2,length(erup_sum),2)],1,mean)

# mik's weathertype recon
wt <- read.table('../comparison_data/CAP7_1763-2009.txt',stringsAsFactors=F,head=T)
wtnt <- read.table('../comparison_data/CAP7_1763-2009_NoTemp.txt',stringsAsFactors=F,head=T)
wtleg=c('NE','WSW','W','E','HP','N','WC')

i=1
for (cyr in erup_win) {
  wpos <- which(substr((wt[,1]),1,7)==paste((cyr-1),'10',sep='-')|
                  substr((wt[,1]),1,7)==paste((cyr-1),'11',sep='-')|
                  substr((wt[,1]),1,7)==paste((cyr-1),'12',sep='-')|
                  substr((wt[,1]),1,7)==paste(cyr,'01',sep='-')|
                  substr((wt[,1]),1,7)==paste(cyr,'02',sep='-')|
                  substr((wt[,1]),1,7)==paste(cyr,'03',sep='-'))
  wpos2 <- which(as.numeric(substr((wt[,1]),1,4))>=(cyr-35)&
                   as.numeric(substr((wt[,1]),1,4))<=(cyr+35)&
                   (as.numeric(substr((wt[,1]),6,7))>=10|
                      as.numeric(substr((wt[,1]),6,7))<=3))
  if (i==1) {
    hwin <- as.numeric(wt[wpos,2])
    hwinclim <- as.numeric(wt[wpos2,2])
  } else {
    hwin <- c(hwin,as.numeric(wt[wpos,2]))
    hwinclim <- c(hwinclim,as.numeric(wt[wpos2,2]))
  }
  i=i+1
}
i=1
for (cyr in erup_sum) {
  wpos <- which(substr((wt[,1]),1,7)==paste((cyr),'04',sep='-')|
                  substr((wt[,1]),1,7)==paste((cyr),'05',sep='-')|
                  substr((wt[,1]),1,7)==paste((cyr),'06',sep='-')|
                  substr((wt[,1]),1,7)==paste(cyr,'07',sep='-')|
                  substr((wt[,1]),1,7)==paste(cyr,'08',sep='-')|
                  substr((wt[,1]),1,7)==paste(cyr,'09',sep='-'))
  wpos2 <- which(as.numeric(substr((wt[,1]),1,4))>=(cyr-35)&
                   as.numeric(substr((wt[,1]),1,4))<=(cyr+35)&
                   (as.numeric(substr((wt[,1]),6,7))>=4|
                      as.numeric(substr((wt[,1]),6,7))<=9))
  if (i==1) {
    hsum <- as.numeric(wt[wpos,2])
    hsumclim <- as.numeric(wt[wpos2,2])
  } else {
    hsum <- c(hsum,as.numeric(wt[wpos,2]))
    hsumclim <- c(hsumclim,as.numeric(wt[wpos2,2]))
  }
  i=i+1
} 
hw = hist(hwin-0.2,breaks=seq(0.66,7.33,0.33),plot=F)
hw$density = hw$counts/sum(hw$counts)*100
hwc = hist(hwinclim-0.2,breaks=seq(0.66,7.33,0.33),plot=F)
hwc$density = hwc$counts/sum(hwc$counts)*100
h <- hw
h$density <- hw$density-hwc$density
hs = hist(hsum+0.2,breaks=seq(0.66,7.33,0.33),plot=F)
hs$density = hs$counts/sum(hs$counts)*100
hsc = hist(hsumclim+0.2,breaks=seq(0.66,7.33,0.33),plot=F)
hsc$density = hsc$counts/sum(hsc$counts)*100
h2 <- hs
h2$density <- hs$density-hsc$density


# plot composites
pdata <- echam
pdf('../figures/iceland/post_iceland_erup_anom_composite_nolaki.pdf', width=8, 
    height=8, paper='special')
layout(matrix(c(1,2,3,3,4,5,6,6,7,7),5,2,byrow=TRUE),height=c(3,1,3,1,3))
par(oma=c(0.5,4,4,0))
levs <- c(-Inf, seq(-1,1,0.25), Inf)
contlevs <- seq(-4,4,2)
scfac <- max(echam.anom$ensmean[echam$names=='u200'])
pdata$data <- array(c(anameananomwin,anameananomsum), c(nrow(echam.anom$ensmean),1,2))
plot_echam3(pdata, varname='temp2', type='data', cex.pt=1.5, 
            latlim=c(0,90),lonlim=c(-180,180),
            wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
            names=c('a)', 'b)'),colnames=c('Oct-Mar','Apr-Sep'),
            statpanel=NULL, add=T, rownames='T2m, Z500, UV200', 
            main='EKF400', units='K',
            addcontours=T, contvarname='gph500', conttype='data',
            contcol='black', contlev=contlevs,
            addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
            veclen=scfac*0.01, vecscale=scfac*2.0, vecwd=0.95, every_x_vec=4)
levs=c(0,80,85,90,95,100,105,110,115,120,Inf)
plot_echam3(pdata, varname='precip', type='data', cex.pt=1.5, 
            latlim=c(0,90),lonlim=c(-180,180),units="%",
            names=c('c)', 'd)'), lev=levs, 
            st.col=NULL, stations=NULL, add=T,addcontours=F,
            wcol='darkgrey',rownames='Precipitation',colnames=rep('',4))
# histograms
par(mar=c(2,2,1,1))
plot(h,col='blue',xlim=c(0.5,7.5),ylim=c(-10,10),main='',freq=F,xaxt='n')
axis(side=1, at=seq(1,7), labels=wtleg)
par(new=T)
plot(h2,col='cyan',xlim=c(0.5,7.5),ylim=c(-10,10),main='',freq=F,xaxt='n')
legend('topright',c('Winter','Summer'),lty=c(1,1),lwd=c(5,5),col=c('blue','cyan')) 
dev.off()  
#  } # end periods loop




# ONLY eruptions that also have forcing in crowley
# load pre calc time series
expname="EKF400_v1.1_correct_71yr_anom"
filenameext = '_anom_landonly_seas_'
load(paste0('../data/indices/',expname,'/indices',filenameext,'allts.Rdata'))

# define plot parameters
wlen=1
linew=1
scal=F
plbem=F
syr=1604
eyr=2003
erup_yrs=c(1721,1783,1845,1860,
           1873,1875,1902,1918) #,2010,2011)
# erup mon >=2 erup yr + 1  
erup_win=c(1722,1784,1846,1861,
           1874,1875,1903,1919)
# erup mon >=8 erup yr + 1  
erup_sum=c(1721,1783,1845,1861,
           1874,1876,1903,1919)
#   Katla             4  1612  10	12
#   Katla	            5	1625	9	  2
#   Katla	            4	1660	11	3
#   Hekla	            4	1693	2	  13
#   Katla	            5	1721	5	  11
#   Oraefajokull	    4	1727	8 	3
#   Katla	            5	1755	10	17
#   Hekla	            4	1766	4	  5
#   Grimsvotn	        4	1783	5  	0
#   Hekla	            4	1845	9	  2
#   Katla	            4	1860	5	  8
#   Grimsvotn	        4	1873	1	  8
#   Askja	            5	1875	1  	1
#   Grimsvotn        	4	1902	12 	0
#   Katla	            4	1918	10	12
#   Hekla	            4	1947	3	  29
#   Eyjafjallajokull	4	2010	3	  20
#   Grimsvotn	        4	2011	5	  21

# spatial maps
# composite all eruptions
# separated summer and winter to really have season after eruption
i=1
for (cyr in erup_win) {
  # change to precip in % 
  # load monthly 70yr climatology data
  if (cyr < 1636) {
    load(file=paste0(echclimpath,'echam_clim_1636-1637_2ndgrid.Rdata'))
  } else {
    load(file=paste0(echclimpath,'echam_clim_',cyr,'-',(cyr+1),'_2ndgrid.Rdata'))  
  }
  pos <- which(echam_clim$names=="precip")[1:4608]
  # seasonal mean
  echmeanclim <- cbind(rowMeans(echam_clim$ensmean[,10:15]),
                       rowMeans(echam_clim$ensmean[,16:21]))
  load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
  # change to precip in % 
  echam.abs <- echam
  echam.anom$ensmean[pos,] <- echam.abs$ensmean[pos,]/(echmeanclim[pos,]*3600*24*30)*100
  analysis.abs <- analysis
  analysis.anom$ensmean[pos,] <- analysis.abs$ensmean[pos,]/
    (echmeanclim[pos,]*3600*24*30)*100
  if (i==1) {
    anameananom <- analysis.anom$ensmean
  } else {
    anameananom <- cbind(anameananom,analysis.anom$ensmean)
  }
  i=i+1
}
anameananomwin <- apply(anameananom[,seq(1,length(erup_win),2)],1,mean)

i=1
for (cyr in erup_sum) {
  # change to precip in % 
  # load monthly 70yr climatology data
  if (cyr < 1636) {
    load(file=paste0(echclimpath,'echam_clim_1636-1637_2ndgrid.Rdata'))
  } else {
    load(file=paste0(echclimpath,'echam_clim_',cyr,'-',(cyr+1),'_2ndgrid.Rdata'))  
  }
  pos <- which(echam_clim$names=="precip")[1:4608]
  # seasonal mean
  echmeanclim <- cbind(rowMeans(echam_clim$ensmean[,10:15]),
                       rowMeans(echam_clim$ensmean[,16:21]))
  load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
  # change to precip in % 
  echam.abs <- echam
  echam.anom$ensmean[pos,] <- echam.abs$ensmean[pos,]/(echmeanclim[pos,]*3600*24*30)*100
  analysis.abs <- analysis
  analysis.anom$ensmean[pos,] <- analysis.abs$ensmean[pos,]/
    (echmeanclim[pos,]*3600*24*30)*100
  if (i==1) {
    anameananom <- analysis.anom$ensmean
  } else {
    anameananom <- cbind(anameananom,analysis.anom$ensmean)
  }
  i=i+1
}
anameananomsum <- apply(anameananom[,seq(2,length(erup_sum),2)],1,mean)

# mik's weathertype recon
wt <- read.table('../comparison_data/CAP7_1763-2009.txt',stringsAsFactors=F,head=T)
wtnt <- read.table('../comparison_data/CAP7_1763-2009_NoTemp.txt',stringsAsFactors=F,head=T)
wtleg=c('NE','WSW','W','E','HP','N','WC')

i=1
for (cyr in erup_win) {
  wpos <- which(substr((wt[,1]),1,7)==paste((cyr-1),'10',sep='-')|
                  substr((wt[,1]),1,7)==paste((cyr-1),'11',sep='-')|
                  substr((wt[,1]),1,7)==paste((cyr-1),'12',sep='-')|
                  substr((wt[,1]),1,7)==paste(cyr,'01',sep='-')|
                  substr((wt[,1]),1,7)==paste(cyr,'02',sep='-')|
                  substr((wt[,1]),1,7)==paste(cyr,'03',sep='-'))
  wpos2 <- which(as.numeric(substr((wt[,1]),1,4))>=(cyr-35)&
                   as.numeric(substr((wt[,1]),1,4))<=(cyr+35)&
                   (as.numeric(substr((wt[,1]),6,7))>=10|
                      as.numeric(substr((wt[,1]),6,7))<=3))
  if (i==1) {
    hwin <- as.numeric(wt[wpos,2])
    hwinclim <- as.numeric(wt[wpos2,2])
  } else {
    hwin <- c(hwin,as.numeric(wt[wpos,2]))
    hwinclim <- c(hwinclim,as.numeric(wt[wpos2,2]))
  }
  i=i+1
}
i=1
for (cyr in erup_sum) {
  wpos <- which(substr((wt[,1]),1,7)==paste((cyr),'04',sep='-')|
                  substr((wt[,1]),1,7)==paste((cyr),'05',sep='-')|
                  substr((wt[,1]),1,7)==paste((cyr),'06',sep='-')|
                  substr((wt[,1]),1,7)==paste(cyr,'07',sep='-')|
                  substr((wt[,1]),1,7)==paste(cyr,'08',sep='-')|
                  substr((wt[,1]),1,7)==paste(cyr,'09',sep='-'))
  wpos2 <- which(as.numeric(substr((wt[,1]),1,4))>=(cyr-35)&
                   as.numeric(substr((wt[,1]),1,4))<=(cyr+35)&
                   (as.numeric(substr((wt[,1]),6,7))>=4|
                      as.numeric(substr((wt[,1]),6,7))<=9))
  if (i==1) {
    hsum <- as.numeric(wt[wpos,2])
    hsumclim <- as.numeric(wt[wpos2,2])
  } else {
    hsum <- c(hsum,as.numeric(wt[wpos,2]))
    hsumclim <- c(hsumclim,as.numeric(wt[wpos2,2]))
  }
  i=i+1
} 
hw = hist(hwin-0.2,breaks=seq(0.66,7.33,0.33),plot=F)
hw$density = hw$counts/sum(hw$counts)*100
hwc = hist(hwinclim-0.2,breaks=seq(0.66,7.33,0.33),plot=F)
hwc$density = hwc$counts/sum(hwc$counts)*100
h <- hw
h$density <- hw$density-hwc$density
hs = hist(hsum+0.2,breaks=seq(0.66,7.33,0.33),plot=F)
hs$density = hs$counts/sum(hs$counts)*100
hsc = hist(hsumclim+0.2,breaks=seq(0.66,7.33,0.33),plot=F)
hsc$density = hsc$counts/sum(hsc$counts)*100
h2 <- hs
h2$density <- hs$density-hsc$density


# plot composites
pdata <- echam
pdf('../figures/iceland/post_iceland_erup_anom_composite_withforc.pdf', width=8, 
    height=8, paper='special')
layout(matrix(c(1,2,3,3,4,5,6,6,7,7),5,2,byrow=TRUE),height=c(3,1,3,1,3))
par(oma=c(0.5,4,4,0))
levs <- c(-Inf, seq(-1,1,0.25), Inf)
contlevs <- seq(-4,4,2)
scfac <- max(echam.anom$ensmean[echam$names=='u200'])
pdata$data <- array(c(anameananomwin,anameananomsum), c(nrow(echam.anom$ensmean),1,2))
plot_echam3(pdata, varname='temp2', type='data', cex.pt=1.5, 
            latlim=c(0,90),lonlim=c(-180,180),
            wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
            names=c('a)', 'b)'),colnames=c('Oct-Mar','Apr-Sep'),
            statpanel=NULL, add=T, rownames='T2m, Z500, UV200', 
            main='EKF400', units='K',
            addcontours=T, contvarname='gph500', conttype='data',
            contcol='black', contlev=contlevs,
            addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
            veclen=scfac*0.01, vecscale=scfac*2.0, vecwd=0.95, every_x_vec=4)
levs=c(0,80,85,90,95,100,105,110,115,120,Inf)
plot_echam3(pdata, varname='precip', type='data', cex.pt=1.5, 
            latlim=c(0,90),lonlim=c(-180,180),units="%",
            names=c('c)', 'd)'), lev=levs, 
            st.col=NULL, stations=NULL, add=T,addcontours=F,
            wcol='darkgrey',rownames='Precipitation',colnames=rep('',4))
# histograms
par(mar=c(2,2,1,1))
plot(h,col='blue',xlim=c(0.5,7.5),ylim=c(-10,10),main='',freq=F,xaxt='n')
axis(side=1, at=seq(1,7), labels=wtleg)
par(new=T)
plot(h2,col='cyan',xlim=c(0.5,7.5),ylim=c(-10,10),main='',freq=F,xaxt='n')
legend('topright',c('Winter','Summer'),lty=c(1,1),lwd=c(5,5),col=c('blue','cyan')) 
dev.off()  
#  } # end periods loop








# spatial maps
# composite all eruptions
# separated summer and winter to really have season after eruption
# 3-yr avg
i=1
for (cyr in erup_win) {
  print(cyr)
  # change to precip in % 
  # load monthly 70yr climatology data
  if (cyr < 1636) {
    load(file=paste0(echclimpath,'echam_clim_1636-1637_2ndgrid.Rdata'))
  } else {
    load(file=paste0(echclimpath,'echam_clim_',cyr,'-',(cyr+1),'_2ndgrid.Rdata'))  
  }
  pos <- which(echam_clim$names=="precip")[1:4608]
  # seasonal mean
  echmeanclim <- cbind(rowMeans(echam_clim$ensmean[,10:15]),
                       rowMeans(echam_clim$ensmean[,16:21]))
  for (j in 0:2) {
    print(j)
    load(file=paste0(prepplotdir,'analysis_',(cyr+j),'.Rdata'))
    # change to precip in % 
    echam.abs <- echam
    echam.anom$ensmean[pos,] <- echam.abs$ensmean[pos,]/(echmeanclim[pos,]*3600*24*30)*100
    analysis.abs <- analysis
    analysis.anom$ensmean[pos,] <- analysis.abs$ensmean[pos,]/
      (echmeanclim[pos,]*3600*24*30)*100
    if (i==1 & j==0) {
      anameananom <- analysis.anom$ensmean
    } else {
      anameananom <- cbind(anameananom,analysis.anom$ensmean)
    }
  }
  i=i+1
}
anameananomwin <- apply(anameananom[,seq(1,length(erup_win),2)],1,mean)

i=1
for (cyr in erup_sum) {
  # change to precip in % 
  # load monthly 70yr climatology data
  if (cyr < 1636) {
    load(file=paste0(echclimpath,'echam_clim_1636-1637_2ndgrid.Rdata'))
  } else {
    load(file=paste0(echclimpath,'echam_clim_',cyr,'-',(cyr+1),'_2ndgrid.Rdata'))  
  }
  pos <- which(echam_clim$names=="precip")[1:4608]
  # seasonal mean
  echmeanclim <- cbind(rowMeans(echam_clim$ensmean[,10:15]),
                       rowMeans(echam_clim$ensmean[,16:21]))
  for (j in 0:2) {
    print(j)
    load(file=paste0(prepplotdir,'analysis_',(cyr+j),'.Rdata'))
    # change to precip in % 
    echam.abs <- echam
    echam.anom$ensmean[pos,] <- echam.abs$ensmean[pos,]/(echmeanclim[pos,]*3600*24*30)*100
    analysis.abs <- analysis
    analysis.anom$ensmean[pos,] <- analysis.abs$ensmean[pos,]/
      (echmeanclim[pos,]*3600*24*30)*100
    if (i==1 & j==0) {
      anameananom <- analysis.anom$ensmean
    } else {
      anameananom <- cbind(anameananom,analysis.anom$ensmean)
    }
  }
  i=i+1
}
anameananomsum <- apply(anameananom[,seq(2,length(erup_sum),2)],1,mean)

# mik's weathertype recon
wt <- read.table('../comparison_data/CAP7_1763-2009.txt',stringsAsFactors=F,head=T)
wtnt <- read.table('../comparison_data/CAP7_1763-2009_NoTemp.txt',stringsAsFactors=F,head=T)
wtleg=c('NE','WSW','W','E','HP','N','WC')

i=1
for (cyr in erup_win) {
  wpos01 <- which(substr((wt[,1]),1,7)==paste((cyr-1),'10',sep='-')|
                  substr((wt[,1]),1,7)==paste((cyr-1),'11',sep='-')|
                  substr((wt[,1]),1,7)==paste((cyr-1),'12',sep='-')|
                  substr((wt[,1]),1,7)==paste(cyr,'01',sep='-')|
                  substr((wt[,1]),1,7)==paste(cyr,'02',sep='-')|
                  substr((wt[,1]),1,7)==paste(cyr,'03',sep='-'))
  wpos02 <- which(substr((wt[,1]),1,7)==paste((cyr),'10',sep='-')|
                  substr((wt[,1]),1,7)==paste((cyr),'11',sep='-')|
                  substr((wt[,1]),1,7)==paste((cyr),'12',sep='-')|
                  substr((wt[,1]),1,7)==paste((cyr+1),'01',sep='-')|
                  substr((wt[,1]),1,7)==paste((cyr+1),'02',sep='-')|
                  substr((wt[,1]),1,7)==paste((cyr+1),'03',sep='-'))
  wpos03 <- which(substr((wt[,1]),1,7)==paste((cyr+1),'10',sep='-')|
                  substr((wt[,1]),1,7)==paste((cyr+1),'11',sep='-')|
                  substr((wt[,1]),1,7)==paste((cyr+1),'12',sep='-')|
                  substr((wt[,1]),1,7)==paste((cyr+2),'01',sep='-')|
                  substr((wt[,1]),1,7)==paste((cyr+2),'02',sep='-')|
                  substr((wt[,1]),1,7)==paste((cyr+2),'03',sep='-'))
  wpos <- c(wpos01,wpos02,wpos03)
  wpos2 <- which(as.numeric(substr((wt[,1]),1,4))>=(cyr-35)&
                   as.numeric(substr((wt[,1]),1,4))<=(cyr+35)&
                   (as.numeric(substr((wt[,1]),6,7))>=10|
                      as.numeric(substr((wt[,1]),6,7))<=3))
  if (i==1) {
    hwin <- as.numeric(wt[wpos,2])
    hwinclim <- as.numeric(wt[wpos2,2])
  } else {
    hwin <- c(hwin,as.numeric(wt[wpos,2]))
    hwinclim <- c(hwinclim,as.numeric(wt[wpos2,2]))
  }
  i=i+1
}
i=1
for (cyr in erup_sum) {
  wpos01 <- which(substr((wt[,1]),1,7)==paste((cyr),'04',sep='-')|
                    substr((wt[,1]),1,7)==paste((cyr),'05',sep='-')|
                    substr((wt[,1]),1,7)==paste((cyr),'06',sep='-')|
                    substr((wt[,1]),1,7)==paste(cyr,'07',sep='-')|
                    substr((wt[,1]),1,7)==paste(cyr,'08',sep='-')|
                    substr((wt[,1]),1,7)==paste(cyr,'09',sep='-'))
  wpos02 <- which(substr((wt[,1]),1,7)==paste((cyr+1),'04',sep='-')|
                    substr((wt[,1]),1,7)==paste((cyr+1),'05',sep='-')|
                    substr((wt[,1]),1,7)==paste((cyr+1),'06',sep='-')|
                    substr((wt[,1]),1,7)==paste((cyr+1),'07',sep='-')|
                    substr((wt[,1]),1,7)==paste((cyr+1),'08',sep='-')|
                    substr((wt[,1]),1,7)==paste((cyr+1),'09',sep='-'))
  wpos03 <- which(substr((wt[,1]),1,7)==paste((cyr+2),'04',sep='-')|
                    substr((wt[,1]),1,7)==paste((cyr+2),'05',sep='-')|
                    substr((wt[,1]),1,7)==paste((cyr+2),'06',sep='-')|
                    substr((wt[,1]),1,7)==paste((cyr+2),'07',sep='-')|
                    substr((wt[,1]),1,7)==paste((cyr+2),'08',sep='-')|
                    substr((wt[,1]),1,7)==paste((cyr+2),'09',sep='-'))
  wpos <- c(wpos01,wpos02,wpos03)
  wpos2 <- which(as.numeric(substr((wt[,1]),1,4))>=(cyr-35)&
                   as.numeric(substr((wt[,1]),1,4))<=(cyr+35)&
                   (as.numeric(substr((wt[,1]),6,7))>=4|
                      as.numeric(substr((wt[,1]),6,7))<=9))
  if (i==1) {
    hsum <- as.numeric(wt[wpos,2])
    hsumclim <- as.numeric(wt[wpos2,2])
  } else {
    hsum <- c(hsum,as.numeric(wt[wpos,2]))
    hsumclim <- c(hsumclim,as.numeric(wt[wpos2,2]))
  }
  i=i+1
} 
hw = hist(hwin-0.2,breaks=seq(0.66,7.33,0.33),plot=F)
hw$density = hw$counts/sum(hw$counts)*100
hwc = hist(hwinclim-0.2,breaks=seq(0.66,7.33,0.33),plot=F)
hwc$density = hwc$counts/sum(hwc$counts)*100
h <- hw
h$density <- hw$density-hwc$density
hs = hist(hsum+0.2,breaks=seq(0.66,7.33,0.33),plot=F)
hs$density = hs$counts/sum(hs$counts)*100
hsc = hist(hsumclim+0.2,breaks=seq(0.66,7.33,0.33),plot=F)
hsc$density = hsc$counts/sum(hsc$counts)*100
h2 <- hs
h2$density <- hs$density-hsc$density


# plot composites
pdata <- echam
pdf('../figures/iceland/post_iceland_erup_anom_composite_3yr.pdf', width=8, 
    height=8, paper='special')
layout(matrix(c(1,2,3,3,4,5,6,6,7,7),5,2,byrow=TRUE),height=c(3,1,3,1,3))
par(oma=c(0.5,4,4,0))
levs <- c(-Inf, seq(-1,1,0.25), Inf)
contlevs <- seq(-4,4,2)
scfac <- max(echam.anom$ensmean[echam$names=='u200'])
pdata$data <- array(c(anameananomwin,anameananomsum), c(nrow(echam.anom$ensmean),1,2))
plot_echam3(pdata, varname='temp2', type='data', cex.pt=1.5, 
            latlim=c(0,90),lonlim=c(-180,180),
            wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
            names=c('a)', 'b)'),colnames=c('Oct-Mar','Apr-Sep'),
            statpanel=NULL, add=T, rownames='T2m, Z500, UV200', 
            main='EKF400', units='K',
            addcontours=T, contvarname='gph500', conttype='data',
            contcol='black', contlev=contlevs,
            addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
            veclen=scfac*0.01, vecscale=scfac*2.0, vecwd=0.95, every_x_vec=4)
levs=c(0,80,85,90,95,100,105,110,115,120,Inf)
plot_echam3(pdata, varname='precip', type='data', cex.pt=1.5, 
            latlim=c(0,90),lonlim=c(-180,180),units="%",
            names=c('c)', 'd)'), lev=levs, 
            st.col=NULL, stations=NULL, add=T,addcontours=F,
            wcol='darkgrey',rownames='Precipitation',colnames=rep('',4))
# histograms
par(mar=c(2,2,1,1))
plot(h,col='blue',xlim=c(0.5,7.5),ylim=c(-10,10),main='',freq=F,xaxt='n')
axis(side=1, at=seq(1,7), labels=wtleg)
par(new=T)
plot(h2,col='cyan',xlim=c(0.5,7.5),ylim=c(-10,10),main='',freq=F,xaxt='n')
legend('topright',c('Winter','Summer'),lty=c(1,1),lwd=c(5,5),col=c('blue','cyan')) 
dev.off()  
#  } # end periods loop







  
  
  
# just plot anomaly composite for 1-month after eruption
# 11x3 matrix: 1 months x 3 temp/gph/wind, prec/slp, w-type hist
  erup_yrs_orig <- c(1612,1625,1660,1693,1721,1727,1755,1766,1783,1845,1860,
                     1873,1875,1902,1918,1947) #,2010,2011)
  erup_yrs <- c(1612,1625,1660,1693,1721,1727,1755,1766,1783,1845,1860,
                1873,1875,1903,1918,1947) #,2010,2011)
  erup_yrs_nolaki <- c(1612,1625,1660,1693,1721,1727,1755,1766,1845,1860,
              1873,1875,1903,1918,1947) #,2010,2011)
  erup_mon_orig <- c(10,9,11,2,5,8,10,4,5,9,5,1,1,12,10,3) #,4,6)
  erup_mon <- c(11,10,12,3,6,9,11,5,6,10,6,2,2,1,11,4) #,4,6)
  erup_mon_nolaki <- c(11,10,12,3,6,9,11,5,10,6,2,2,1,11,4) #,4,6)
  erup_day <- c(12,2,3,13,11,3,17,5,15,2,8,8,1,15,12,29) #,20,21)
# day=15 for unknown days
#   Katla             4 1612	10	12
#   Katla	            5	1625	9	  2
#   Katla	            4	1660	11	3
#   Hekla	            4	1693	2	  13
#   Katla	            5	1721	5	  11
#   Oraefajokull	    4	1727	8 	3
#   Katla	            5	1755	10	17
#   Hekla	            4	1766	4	  5
#   Grimsvotn	        4	1783	5  	0
#   Hekla	            4	1845	9	  2
#   Katla	            4	1860	5	  8
#   Grimsvotn	        4	1873	1	  8
#   Askja	            5	1875	1  	1
#   Grimsvotn        	4	1902	12 	0
#   Katla	            4	1918	10	12
#   Hekla	            4	1947	3	  29
#   Eyjafjallajokull	4	2010	3	  20
#   Grimsvotn	        4	2011	5	  21

# read clim for anom
  i=1
  for (cyr in erup_yrs) {
    print(c(cyr,i))
    if (cyr < 1636) {
      load(file=paste0(echclimpath,'echam_clim_1636-1637_2ndgrid.Rdata'))
    } else {
      load(file=paste0(echclimpath,'echam_clim_',cyr,'-',(cyr+1),'_2ndgrid.Rdata'))  
    }
    posp <- which(echam_clim$names=="precip")[1:4608]
#    post <- which(echam_clim$names=="temp2")[1:4608]
#    poss <- which(echam_clim$names=="slp")[1:4608]
    # select month
    echmeanclim <- echam_clim$ensmean[,erup_mon[i]]
# for (cyr in seq(syr,eyr)) {
    if (erup_mon[i]>9){
      load(file=paste0(prepplotdirmon,'analysis_',cyr+1,'.Rdata'))
    } else {
      load(file=paste0(prepplotdirmon,'analysis_',cyr,'.Rdata'))
    }
#    # strangely monthly data is still in sixmonstatevector???
#   #analysis.anom$ensmean <- array(analysis.anom$ensmean,c((dim(analysis.anom$ensmean)[1]/6),
#   #                                                       dim(analysis.anom$ensmean)[2]*6))
    analysis.abs <- analysis
    echam.abs <- echam
    # select month
    if (erup_mon[i]>9){
      analysis.abs$ensmean <- analysis.abs$ensmean[,(erup_mon[i]-9)]
      echam.abs$ensmean <- echam.abs$ensmean[,(erup_mon[i]-9)]
      analysis.anom$ensmean <- analysis.anom$ensmean[,(erup_mon[i]-9)]
      echam.anom$ensmean <- echam.anom$ensmean[,(erup_mon[i]-9)]
    } else {
      analysis.abs$ensmean <- analysis.abs$ensmean[,(erup_mon[i]+3)]
      echam.abs$ensmean <- echam.abs$ensmean[,(erup_mon[i]+3)]
      analysis.anom$ensmean <- analysis.anom$ensmean[,(erup_mon[i]+3)]
      echam.anom$ensmean <- echam.anom$ensmean[,(erup_mon[i]+3)]
    }
    ## change from abs TEMP to anom
    #analysis.abs$ensmean[post] <- analysis.abs$ensmean[post]-(echmeanclim[post]-273.15)
    #echam.abs$ensmean[post] <- echam.abs$ensmean[post]-(echmeanclim[post]-273.15)
    #analysis.abs$ensmean[poss] <- analysis.abs$ensmean[post]-(echmeanclim[poss]/100)
    #echam.abs$ensmean[poss] <- echam.abs$ensmean[post]-(echmeanclim[poss]/100)
    # change to precip in % 
    analysis.anom$ensmean[posp] <- analysis.abs$ensmean[posp]/
      (echmeanclim[posp]*3600*24*30)*100
    echam.anom$ensmean[posp] <- echam.abs$ensmean[posp]/
      (echmeanclim[posp]*3600*24*30)*100
    # change other variables from abs to anom
    #analysis.abs$ensmean[-c(post,posp)] <- analysis.abs$ensmean[-c(post,posp,poss)]-
    #                                         echmeanclim[-c(post,posp,poss)]
    #echam.abs$ensmean[-c(post,posp)] <- echam.abs$ensmean[-c(post,posp,poss)]-
    #                                      echmeanclim[-c(post,posp,poss)]
    
    if (i==1) {
      anameananom <- analysis.anom$ensmean
      echmeananom <- echam.anom$ensmean
    } else {
      anameananom <- cbind(anameananom,analysis.anom$ensmean)
      echmeananom <- cbind(echmeananom,echam.anom$ensmean)
    }
  i=i+1
  }
  anameananom <- apply(anameananom,1,mean)
  echmeananom <- apply(echmeananom,1,mean)

  # mik's weathertype recon
  wt <- read.table('../comparison_data/CAP7_1763-2009.txt',stringsAsFactors=F,head=T)
  wtnt <- read.table('../comparison_data/CAP7_1763-2009_NoTemp.txt',stringsAsFactors=F,head=T)
  wtleg=c('NE','WSW','W','E','HP','N','WC')

  i=1
  for (cyr in erup_yrs) {
    #print(c(i,cyr))
    if (erup_mon[i]<10){
      wpos <- which(substr((wt[,1]),1,7)==paste(cyr,erup_mon[i],sep='-0'))
    } else {
      wpos <- which(substr((wt[,1]),1,7)==paste(cyr,erup_mon[i],sep='-'))
    }
    wpos2 <- which(as.numeric(substr((wt[,1]),1,4))>=(cyr-35)&
                   as.numeric(substr((wt[,1]),1,4))<=(cyr+35)&
                   as.numeric(substr((wt[,1]),6,7))==erup_mon[i])
    if (i==1) {
      h1 <- as.numeric(wt[wpos,2])
      hclim <- as.numeric(wt[wpos2,2])
    } else {
      h1 <- c(h1,as.numeric(wt[wpos,2]))
      hclim <- c(hclim,as.numeric(wt[wpos2,2]))
    }
    #print(h1)
    i=i+1
  }
  h2 = hist(h1-0.2,breaks=seq(0.66,7.33,0.33),plot=F)
  h2$density = h2$counts/sum(h2$counts)*100
  h2c = hist(hclim-0.2,breaks=seq(0.66,7.33,0.33),plot=F)
  h2c$density = h2c$counts/sum(h2c$counts)*100
  h <- h2
  h$density <- h2$density-h2c$density


  # plot post volc months composite
  pdata <- echam
  pdf(paste0('../figures/iceland/post_iceland_erup_monthly_composite_v3.pdf'), width=5, 
     height=10, paper='special')
    layout(matrix(c(1,2,3,4,5),5,1,byrow=TRUE),height=c(3,1,3,1,3))
    par(oma=c(0.5,2,2,0))
    levs <- c(-Inf, seq(-1,1,0.25), Inf)
    contlevs <- seq(-4,4,2) #4500,5900,100)       #-4,4,2)
    scfac <- max(echam.anom$ensmean[echam$names=='u200'])
    pdata$data <- array(anameananom, c(nrow(echam.anom$ensmean),1,1))
    plot_echam3(pdata, varname='temp2', type='data', cex.pt=1.7, 
              latlim=c(0,90),lonlim=c(-180,180),
              wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
              names='', #c('e)', 'f)', 'g)', 'h)'),
              colnames=rep(' ',12), #seq(syr,eyr),
              statpanel=NULL, add=T, rownames='T2m,Z500 anom',#,UV200 anom', 
              main='EKF400 1783-1784', units='K',
              addcontours=T, contvarname='gph500', conttype='data',
              contcol='black', contlev=contlevs,
              addvectors=F, #vecnames=c('u200','v200'), veccol='cyan', 
              #veclen=scfac*0.01, vecscale=scfac*0.2, 
              #vecwd=0.95, every_x_vec=4,
              ti=1)
    levs=c(0,80,85,90,95,100,105,110,115,120,Inf)
    contlevs <- seq(-4,4,2) #980,1030,5) #-4,4,2)
    plot_echam3(pdata, varname='precip', type='data', cex.pt=1.7, 
            latlim=c(0,90),lonlim=c(-180,180),
            wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
            names='', #c('e)', 'f)', 'g)', 'h)'),
            colnames=rep(' ',12), #seq(syr,eyr),
            statpanel=NULL, add=T, rownames='Prec,SLP anom', 
            main='', units='%',
            addcontours=T, contvarname='slp', conttype='data',
            contcol='black', contlev=contlevs,ti=1) #,
            #addvectors=F, vecnames=c('u200','v200'), veccol='cyan', 
            #veclen=scfac*0.01, vecscale=scfac*2.0, 
            #vecwd=0.95, every_x_vec=4)

  # histograms
    par(mar=c(2,2,1,1))
    plot(h,col='blue',xlim=c(0.5,7.5),ylim=c(-10,10),main='',freq=F,xaxt='n')
    axis(side=1, at=seq(1,7), labels=wtleg, xlab='')
  dev.off()  






# LOOK AT ECHAM VOLC FORCING OF ERUPTIONS
  if (machine=="macbook"){
    vol <- read.table(
      "/Volumes/data/climdata/forcings/echam5_input/crowley_volcanic_forcing.txt",
      header=T)  
  } else {
    vol <- read.table(
      "/mnt/climstor/ccc400/echam5_input/crowley_volcanic_forcing.txt",
      header=T)
  }
  erup_date <- as.Date(paste(erup_yrs_orig,erup_mon_orig,erup_day,sep='-'),"%Y-%m-%d")
  forc_date <- as.Date(format(date_decimal(vol[,1]), "%Y-%m-%d"))
  pdf('../figures/iceland/iceland_forcing.pdf',width=12,height=12,paper='special')
    par(mfrow=c(4,4),mar=c(2,2,2,1))
    for (i in 1:length(erup_yrs)) {
      pos <- which.closest(forc_date,erup_date[i])
      posper <- (pos-6):(pos+36)
      plot(forc_date[posper],vol[posper,2],ty='l',
           col='red',main='', #erup_yrs_orig[i],
           ylab='AOD',xlab='')
    }
  dev.off()


# LOOK AT NAO
  load(file=paste0('../data/indices/',expname,'/indices',filenameext,'allts.Rdata'))
  pos <- which(aind.allts$names=="NAO.calc")
  # winter NAO
  naow <-aind.allts$ensmean[pos,seq(1,ncol(aind.allts$ensmean),2)]
  pos2 <- rep(NA,length(erup_yrs))
  for (i in 1:length(erup_yrs)) {
    pos2[i] <- which(erup_yrs[i]==aind.allts$time)
  }
  pdf("../figures/iceland/NAO_winter_boxplots.pdf",width=5,height=4)
    boxplot(naow,
          naow[(pos2-rep(seq(3,1,-1),length(pos2)))],
          naow[pos2+rep(seq(0,2,1),length(pos2))],
          names=c('all','3yrs bef erup','3yrs aft erup'),
          col=c(2,4,5))
  dev.off()
  t.test(naow,naow[pos2+rep(seq(0,2,1),length(pos2))])
  t.test(naow[(pos2-rep(seq(3,1,-1),length(pos2)))],
       naow[pos2+rep(seq(0,2,1),length(pos2))])
# result: NAO not changed after iceland eruptions 

# summer NAO
  naos <-aind.allts$ensmean[pos,seq(2,ncol(aind.allts$ensmean),2)]
  pos2 <- rep(NA,length(erup_yrs))
  for (i in 1:length(erup_yrs)) {
    pos2[i] <- which(erup_yrs[i]==aind.allts$time)
  }
  pdf("../figures/iceland/NAO_summer_boxplots.pdf",width=5,height=4)
    boxplot(naos,
        naos[(pos2-rep(seq(3,1,-1),length(pos2)))],
        naos[pos2+rep(seq(0,2,1),length(pos2))],
        names=c('all','3yrs bef erup','3yrs aft erup'),
        col=c(2,4,5))
  dev.off()
  t.test(naos,naos[pos2+rep(seq(0,2,1),length(pos2))])
  t.test(naos[(pos2-rep(seq(3,1,-1),length(pos2)))],
       naos[pos2+rep(seq(0,2,1),length(pos2))])
# result: NAO not changed after iceland eruptions 






# LOOK AT REANALYSIS FOR 21TH AND 20TH CENTURY ICELAND ERUPTIONS
# done externally with netcdf files  
  
  
  
  
  
#   # single eruption plots
#   coldper <- erup_yrs
#   for (syr in coldper) {
#     eyr <- syr+3
#     # change to precip in % 
#     # load monthly 70yr climatology data
#     if (syr < 1636) {
#       load(file=paste0(echclimpath,'echam_clim_1636-1637_2ndgrid.Rdata'))
#     } else {
#       load(file=paste0(echclimpath,'echam_clim_',syr+3,'-',(syr+4),'_2ndgrid.Rdata'))  
#     }
#     pos <- which(echam_clim$names=="precip")[1:4608]
#     # seasonal mean
#     echmeanclim <- cbind(rowMeans(echam_clim$ensmean[,10:15]),
#                          rowMeans(echam_clim$ensmean[,16:21]))
#     for (cyr in seq(syr,eyr)) {
#       load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
#       # change to precip in % 
#       echam.abs <- echam
#       echam.anom$ensmean[pos,] <- echam.abs$ensmean[pos,]/(echmeanclim[pos,]*3600*24*30)*100
#       analysis.abs <- analysis
#       analysis.anom$ensmean[pos,] <- analysis.abs$ensmean[pos,]/
#         (echmeanclim[pos,]*3600*24*30)*100
#       if (cyr==syr) {
#         anameananom <- analysis.anom$ensmean
#       } else {
#         anameananom <- cbind(anameananom,analysis.anom$ensmean)
#       }
#     }
#     
# 
#     
#     # plot post volc years
#     pdata <- echam
#     ti=1
#     pdf(paste0('../figures/iceland/post_erup_anom_',syr,'.pdf'), width=16, 
#         height=12, paper='special')
#     layout(matrix(c(1,2,3,4,5,5,5,5,6,7,8,9,10,10,10,10,11,12,13,14,15,15,15,15,
#                     16,17,18,19,20,20,20,20),8,4,byrow=TRUE),height=c(3,1,3,1,3,1,3,1))
#     par(oma=c(0.5,4,4,0))
#     levs <- c(-Inf, seq(-1,1,0.25), Inf)
#     contlevs <- seq(-4,4,2)
#     scfac <- max(echam.anom$ensmean[echam$names=='u200'])
#     pdata$data <- array(c(anameananom[,ti],anameananom[,ti+2],anameananom[,ti+4],
#                           anameananom[,ti+6]), c(nrow(echam.anom$ensmean),1,4))
#     plot_echam3(pdata, varname='temp2', type='data', cex.pt=1.5, 
#                 latlim=c(0,90),lonlim=c(-180,180),
#                 wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
#                 names=c('e)', 'f)', 'g)', 'h)'),colnames=seq(syr,eyr),
#                 statpanel=NULL, add=T, rownames='T2m, Z500, UV200 (Oct-Mar)', 
#                 main='EKF400', units='K',
#                 addcontours=T, contvarname='gph500', conttype='data',
#                 contcol='black', contlev=contlevs,
#                 addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
#                 veclen=scfac*0.01, vecscale=scfac*2.0, vecwd=0.95, every_x_vec=4)
#     pdata$data <- array(c(anameananom[,ti],anameananom[,ti+2],anameananom[,ti+4],
#                           anameananom[,ti+6]), c(nrow(echam.anom$ensmean),1,4))
#     levs=c(0,80,85,90,95,100,105,110,115,120,Inf)
#     plot_echam3(pdata, varname='precip', type='data', cex.pt=1.5, 
#                 latlim=c(0,90),lonlim=c(-180,180),units="%",
#                 names=c('m)', 'n)','o)', 'p)'), lev=levs, 
#                 st.col=NULL, stations=NULL, add=T,addcontours=F,
#                 wcol='darkgrey',rownames='Prec. (Oct-Mar)',colnames=rep('',4))
#     
#     ti=ti+1 # Apr-Sep season
#     levs <- c(-Inf, seq(-1,1,0.25), Inf)
#     pdata$data <- array(c(anameananom[,ti],anameananom[,ti+2],anameananom[,ti+4],
#                           anameananom[,ti+6]), c(nrow(echam.anom$ensmean),1,4))
#     plot_echam3(pdata, varname='temp2', type='data', cex.pt=1.5, 
#                 latlim=c(0,90),lonlim=c(-180,180),
#                 wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
#                 names=c('e)', 'f)', 'g)', 'h)'),colnames=rep(' ',4),
#                 statpanel=NULL, add=T, rownames='T2m, Z500, UV200  (Apr-Sep)', 
#                 main=' ', units='K',
#                 addcontours=T, contvarname='gph500', conttype='data',
#                 contcol='black', contlev=contlevs,
#                 addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
#                 veclen=scfac*0.01, vecscale=scfac*2.0, vecwd=0.95, every_x_vec=4)
#     pdata$data <- array(c(anameananom[,ti],anameananom[,ti+2],anameananom[,ti+4],
#                           anameananom[,ti+6]), c(nrow(echam.anom$ensmean),1,4))
#     levs=c(0,80,85,90,95,100,105,110,115,120,Inf)
#     plot_echam3(pdata, varname='precip', type='data', cex.pt=1.5, 
#                 latlim=c(0,90),lonlim=c(-180,180),units="%",
#                 names=c('m)', 'n)','o)', 'p)'), lev=levs, 
#                 st.col=NULL, stations=NULL, add=T,addcontours=F,
#                 wcol='darkgrey',rownames='Prec. (Apr-Sep)',colnames=rep('',4))
#     dev.off()  
#   } # end periods loop
#   
#   
#   
#   
#   
#   
#   # mik's weathertype recon
#   wt <- read.table('../comparison_data/CAP7_1763-2009.txt',stringsAsFactors=F,head=T)
#   wtnt <- read.table('../comparison_data/CAP7_1763-2009_NoTemp.txt',stringsAsFactors=F,head=T)
#   wtleg=c('NE','WSW','W','E','HP','N','WC')
#   # pos <- which(as.numeric(substr((wt[,1]),1,4))==1783)
#   # h = hist(as.numeric(wt[pos,2])-0.2,breaks=seq(0.66,7.33,0.33),plot=F)
#   # h$density = h$counts/sum(h$counts)*100
#   # plot(h,col='blue',xlim=c(0.5,7.5),ylim=c(0,25),main='CAP7 1783',freq=F,xaxt='n')
#   # axis(side=1, at=seq(1,7), labels=wtleg)
#   # h2 = hist(as.numeric(wtnt[pos,2])+0.2,breaks=seq(0.66,7.33,0.33),plot=F)
#   # h2$density = h2$counts/sum(h2$counts)*100
#   # par(new=T)
#   # plot(h2,col='cyan',xlim=c(0.5,7.5),ylim=c(0,25),main='',freq=F,xaxt='n')
#   
#   
#   coldper <- erup_yrs
#   for (syr in coldper) {
#     eyr <- syr+3
#     # change to precip in % 
#     # load monthly 70yr climatology data
#     if (syr < 1636) {
#       load(file=paste0(echclimpath,'echam_clim_1636-1637_2ndgrid.Rdata'))
#     } else {
#       load(file=paste0(echclimpath,'echam_clim_',syr+3,'-',(syr+4),'_2ndgrid.Rdata'))  
#     }
#     pos <- which(echam_clim$names=="precip")[1:4608]
#     # seasonal mean
#     echmeanclim <- cbind(rowMeans(echam_clim$ensmean[,10:15]),
#                          rowMeans(echam_clim$ensmean[,16:21]))
#     for (cyr in seq(syr,eyr)) {
#       load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
#       # change to precip in % 
#       echam.abs <- echam
#       echam.anom$ensmean[pos,] <- echam.abs$ensmean[pos,]/(echmeanclim[pos,]*3600*24*30)*100
#       analysis.abs <- analysis
#       analysis.anom$ensmean[pos,] <- analysis.abs$ensmean[pos,]/
#         (echmeanclim[pos,]*3600*24*30)*100
#       if (cyr==syr) {
#         anameananom <- analysis.anom$ensmean
#       } else {
#         anameananom <- cbind(anameananom,analysis.anom$ensmean)
#       }
#     }
#     
#     # plot post volc years
#     pdata <- echam
#     ti=1
#     pdf(paste0('../figures/iceland/post_erup_anom_hist_',syr,'.pdf'), width=16, 
#         height=12, paper='special')
#     #layout(matrix(c(1,2,3,4,5,5,5,5,6,7,8,9,10,10,10,10,11,12,13,14,15,15,15,15,
#     #                16,17,18,19,20,20,20,20),8,4,byrow=TRUE),height=c(3,1,3,1,3,1,3,1))
#     layout(matrix(c(1,2,3,4,5,5,5,5,6,7,8,9,10,11,12,13,14,14,14,14,15,
#                     16,17,18),6,4,byrow=TRUE),height=c(3,1,3,3,1,3))
#     par(oma=c(0.5,4,4,0))
#     levs <- c(-Inf, seq(-2,2,0.5), Inf)
#     contlevs <- seq(-4,4,2)
#     scfac <- max(echam.anom$ensmean[echam$names=='u200'])
#     pdata$data <- array(c(anameananom[,ti],anameananom[,ti+2],anameananom[,ti+4],
#                           anameananom[,ti+6]), c(nrow(echam.anom$ensmean),1,4))
#     plot_echam3(pdata, varname='temp2', type='data', cex.pt=1.7, 
#                 latlim=c(0,90),lonlim=c(-180,180),
#                 wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
#                 names=c('e)', 'f)', 'g)', 'h)'),colnames=seq(syr,eyr),
#                 statpanel=NULL, add=T, rownames='T2m, Z500, UV200 (Oct-Mar)', 
#                 main='EKF400', units='K',
#                 addcontours=T, contvarname='gph500', conttype='data',
#                 contcol='black', contlev=contlevs,
#                 addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
#                 veclen=scfac*0.01, vecscale=scfac*2.0, vecwd=0.95, every_x_vec=4)
#     
#     # histograms
#     for (cyr in seq(syr,eyr)) {
#       wpos <- which(substr((wt[,1]),1,7)==paste((cyr-1),'10',sep='-')|
#                       substr((wt[,1]),1,7)==paste((cyr-1),'11',sep='-')|
#                       substr((wt[,1]),1,7)==paste((cyr-1),'12',sep='-')|
#                       substr((wt[,1]),1,7)==paste(cyr,'01',sep='-')|
#                       substr((wt[,1]),1,7)==paste(cyr,'02',sep='-')|
#                       substr((wt[,1]),1,7)==paste(cyr,'03',sep='-'))
#       wpos2 <- which(as.numeric(substr((wt[,1]),1,4))>=(cyr-35)&
#                        as.numeric(substr((wt[,1]),1,4))<=(cyr+35)&
#                        (as.numeric(substr((wt[,1]),6,7))>=10|
#                           as.numeric(substr((wt[,1]),6,7))<=3))
#       par(mar=c(2,2,1,1))
#       h1 = hist(as.numeric(wt[wpos,2])-0.2,breaks=seq(0.66,7.33,0.33),plot=F)
#       h1$density = h1$counts/sum(h1$counts)*100
#       #plot(h,col='blue',xlim=c(0.5,7.5),ylim=c(0,25),main='CAP7',freq=F,xaxt='n')
#       #axis(side=1, at=seq(1,7), labels=wtleg)
#       h2 = hist(as.numeric(wt[wpos2,2])-0.2,breaks=seq(0.66,7.33,0.33),plot=F)
#       h2$density = h2$counts/sum(h2$counts)*100
#       h <- h1
#       h$density <- h1$density-h2$density
#       #par(new=T)
#       plot(h,col='blue',xlim=c(0.5,7.5),ylim=c(-10,10),main='',freq=F,xaxt='n')
#       axis(side=1, at=seq(1,7), labels=wtleg)
#       h3 = hist(as.numeric(wtnt[wpos,2])+0.2,breaks=seq(0.66,7.33,0.33),plot=F)
#       h3$density = h3$counts/sum(h3$counts)*100
#       #par(new=T)
#       #plot(h3,col='red',xlim=c(0.5,7.5),ylim=c(0,25),main=' ',freq=F,xaxt='n')
#       h4 = hist(as.numeric(wtnt[wpos2,2])+0.2,breaks=seq(0.66,7.33,0.33),plot=F)
#       h4$density = h4$counts/sum(h4$counts)*100
#       h5 <- h3
#       h5$density <- h3$density-h4$density
#       par(new=T)
#       plot(h5,col='cyan',xlim=c(0.5,7.5),ylim=c(-10,10),main='',freq=F,xaxt='n')
#     }
#     
#     ti=ti+1 # Apr-Sep season
#     levs <- c(-Inf, seq(-2,2,0.5), Inf)
#     pdata$data <- array(c(anameananom[,ti],anameananom[,ti+2],anameananom[,ti+4],
#                           anameananom[,ti+6]), c(nrow(echam.anom$ensmean),1,4))
#     plot_echam3(pdata, varname='temp2', type='data', cex.pt=1.7, 
#                 latlim=c(0,90),lonlim=c(-180,180),
#                 wcol='darkgrey',lev=levs, st.col=NULL, stations=calibrate, 
#                 names=c('e)', 'f)', 'g)', 'h)'),colnames=rep(' ',4),
#                 statpanel=NULL, add=T, rownames='T2m, Z500, UV200  (Apr-Sep)', 
#                 main=' ', units='K',
#                 addcontours=T, contvarname='gph500', conttype='data',
#                 contcol='black', contlev=contlevs,
#                 addvectors=T, vecnames=c('u200','v200'), veccol='cyan', 
#                 veclen=scfac*0.01, vecscale=scfac*2.0, vecwd=0.95, every_x_vec=4)
#     
#     for (cyr in seq(syr,eyr)) {
#       spos <- which(substr((wt[,1]),1,7)==paste(cyr,'04',sep='-')|
#                       substr((wt[,1]),1,7)==paste(cyr,'05',sep='-')|
#                       substr((wt[,1]),1,7)==paste(cyr,'06',sep='-')|
#                       substr((wt[,1]),1,7)==paste(cyr,'07',sep='-')|
#                       substr((wt[,1]),1,7)==paste(cyr,'08',sep='-')|
#                       substr((wt[,1]),1,7)==paste(cyr,'09',sep='-'))
#       spos2 <- which(as.numeric(substr((wt[,1]),1,4))>=(cyr-35)&
#                        as.numeric(substr((wt[,1]),1,4))<=(cyr+35)&
#                        (as.numeric(substr((wt[,1]),6,7))>=4&
#                           as.numeric(substr((wt[,1]),6,7))<=9))
#       par(mar=c(2,2,1,1))
#       h1 = hist(as.numeric(wt[spos,2])-0.2,breaks=seq(0.66,7.33,0.33),plot=F)
#       h1$density = h1$counts/sum(h1$counts)*100
#       #plot(h,col='blue',xlim=c(0.5,7.5),ylim=c(0,25),main='CAP7',freq=F,xaxt='n')
#       #axis(side=1, at=seq(1,7), labels=wtleg)
#       h2 = hist(as.numeric(wt[spos2,2])-0.2,breaks=seq(0.66,7.33,0.33),plot=F)
#       h2$density = h2$counts/sum(h2$counts)*100
#       h <- h1
#       h$density <- h1$density-h2$density
#       #par(new=T)
#       plot(h,col='blue',xlim=c(0.5,7.5),ylim=c(-10,10),main='',freq=F,xaxt='n')
#       axis(side=1, at=seq(1,7), labels=wtleg)
#       h3 = hist(as.numeric(wtnt[spos,2])+0.2,breaks=seq(0.66,7.33,0.33),plot=F)
#       h3$density = h3$counts/sum(h3$counts)*100
#       #par(new=T)
#       #plot(h3,col='red',xlim=c(0.5,7.5),ylim=c(0,25),main=' ',freq=F,xaxt='n')
#       h4 = hist(as.numeric(wtnt[spos2,2])+0.2,breaks=seq(0.66,7.33,0.33),plot=F)
#       h4$density = h4$counts/sum(h4$counts)*100
#       h5 <- h3
#       h5$density <- h3$density-h4$density
#       par(new=T)
#       plot(h5,col='cyan',xlim=c(0.5,7.5),ylim=c(-10,10),main='',freq=F,xaxt='n')
#     }
#     dev.off()  
#   } # end periods loop
#   
#   
# 
#   
#  
#   
#   
#   # plot ccc400 crowley forcing and GISS aerosols
#   vol <- read.table("/Volumes/data/climdata/forcings/echam5_input/crowley_volcanic_forcing.txt",
#                     header=T)
#   vol[35389:35424,]
#   volENH_ts_nh <- ts(vol[,2],start=vol[1,1],freq=36)
#   volENH_ts_sh <- ts(vol[,8],start=vol[1,1],freq=36)
#   # troposperic aerosols do not have Laki aerosols but just natural monthly climatology
#   # plus antropogenic aerosols scaled to population density before 1875
#   # for (yr in seq(1783,1786)) {
#   # nc=nc_open(paste0(
#   #   '/Volumes/DATA/climdata/forcings/echam_ccc400_forcings/GISS_Aerosols/all_sulfate_',
#   #   yr,'.nc'), write=F)
#   #   # print(nc)
#   #   stmp <- ncvar_get(nc, "sulfate")
#   #   if (yr==1783){
#   #     sulf=stmp
#   #   } else {
#   #     sulf=c(sulf,stmp)    
#   #   }
#   # nc_close(nc)
#   pdf('../figures/iceland/aod_forc.pdf',width=15,height=4)
#   par(mfrow=c(1,5),mar=c(4,4,4,1)) #,oma=c(2,2,2,0.2))
#   plot(window(volENH_ts_nh,1783,1786),col='red',xlab='Year',ylab='AOD',
#        main='Laki',ylim=c(0,0.45))
#   lines(window(volENH_ts_sh,1783,1786),col='blue')
#   legend('topright',lty=1,col=c('red','blue'),c('AOD 30-90N','AOD 30-90S'))
#   plot(window(volENH_ts_nh,1808.8,1812),col='red',xlab='Year',ylab='AOD',
#        main='Unknown',ylim=c(0,0.45))
#   lines(window(volENH_ts_sh,1808.8,1812),col='blue')
#   plot(window(volENH_ts_nh,1815,1818),col='red',xlab='Year',ylab='AOD',
#        main='Tambora',ylim=c(0,0.45))
#   lines(window(volENH_ts_sh,1815,1818),col='blue')
#   plot(window(volENH_ts_nh,1831,1834),col='red',xlab='Year',ylab='AOD',
#        main='Babuyan',ylim=c(0,0.45))
#   lines(window(volENH_ts_sh,1831,1834),col='blue')
#   plot(window(volENH_ts_nh,1835,1838),col='red',xlab='Year',ylab='AOD',
#        main='Cosigüina',ylim=c(0,0.45))
#   lines(window(volENH_ts_sh,1835,1838),col='blue')
#   dev.off()
#   #pdf('../figures/laki/tambora_aod_forc.pdf',width=5,height=4)
#   #  par(mfrow=c(1,1),mar=c(4,4,4,1)) #,oma=c(2,2,2,0.2))
#   #  plot(window(volENH_ts,1815,1817),col='red',xlab='Year',ylab='AOD',
#   #     main='Crowley volcanic forcing 30-90N')
#   #dev.off()
#   
#   
#   
#   
#   
#   
#   # calc min temp and gph500 at 50, 60, 70ºN summer/winter landocean ens mem and ens mean
#   #filenameext <- paste0('_anom_landocean_seas_')
#   for (cyr in syr:eyr) {
#     print(cyr)
#     load(file=paste0(prepplotdirseas,'analysis_',cyr,'.Rdata'))
#     post <- which(floor(analysis.anom$lat)==62 & analysis.anom$names=='temp2') 
#     posg <- which(floor(analysis.anom$lat)==62 & analysis.anom$names=='gph500') 
#     # other possible lat==51.294377    #lat==69.946081
#     if (cyr==syr) {
#       tminmean <- apply(analysis.anom$ensmean[post,],2,min)
#       tminmem <- apply(analysis.anom$data[post,,],2:3,min)
#       gminmean <- apply(analysis.anom$ensmean[posg,],2,min)
#       gminmem <- apply(analysis.anom$data[posg,,],2:3,min)
#     } else {
#       tminmean <- rbind(tminmean,apply(analysis.anom$ensmean[post,],2,min))
#       tminmem <- abind(tminmem,apply(analysis.anom$data[post,,],2:3,min),along=3)
#       gminmean <- rbind(gminmean,apply(analysis.anom$ensmean[posg,],2,min))
#       gminmem <- abind(gminmem,apply(analysis.anom$data[posg,,],2:3,min),along=3)
#     }
#   }  
#   # plot minima
#   vf <- c(1600,1650,1660,1783,1809,1815,1835,1883,1886,1902,1912,1980,1991)
#   year <- syr:eyr
#   
#   pdf('../figures/laki/t2m_z500_min_ts_seas.pdf',width=8,height=12)
#   par(mfrow=c(4,1),mar=c(4,4,4,1)) #,oma=c(2,2,2,0.2))
#   plot(year,tminmem[1,1,],ty='l',col=rgb(10,5,0,1,maxColorValue=10),
#        xlab='Year',ylab='Temp. [ºC]', main='62ºN Zonal Minima of Winter Temperature',
#        ylim=c(min(tminmem[1,,]),max(tminmem[1,,])))
#   for (i in 2:30){
#     lines(year,tminmem[1,i,],col=rgb(10,5,0,1,maxColorValue=10))
#   }
#   lines(year,tminmean[,1],col=rgb(10,0,0,10,maxColorValue=10))
#   abline(v=vf,col=rgb(7,7,7,5,maxColorValue=10))
#   legend('bottomright',lty=1,col=c('red','orange'),c('Ensemble mean','Ensemble members'), 
#          bty='o', bg='white',box.col='white')
#   
#   plot(year,tminmem[2,1,],ty='l',col=rgb(10,5,0,1,maxColorValue=10),
#        xlab='Year',ylab='Temperature [ºC]', main='62ºN Zonal Minima of Summer Temperature',
#        ylim=c(min(tminmem[2,,]),max(tminmem[2,,])))
#   for (i in 2:30){
#     lines(year,tminmem[2,i,],col=rgb(10,5,0,1,maxColorValue=10))
#   }
#   lines(year,tminmean[,2],col=rgb(10,0,0,10,maxColorValue=10))
#   abline(v=vf,col=rgb(7,7,7,5,maxColorValue=10))    
#   
#   plot(year,gminmem[1,1,],ty='l',col=rgb(0,10,10,1,maxColorValue=10),
#        xlab='Year',ylab='500 hPa Geopotential Height [m]', main='62ºN Zonal Minima of Winter GPH500',
#        ylim=c(min(gminmem[1,,]),max(gminmem[1,,])))
#   for (i in 2:30){
#     lines(year,gminmem[1,i,],col=rgb(0,10,10,1,maxColorValue=10))
#   }
#   lines(year,gminmean[,1],col=rgb(0,0,10,10,maxColorValue=10))
#   abline(v=vf,col=rgb(7,7,7,5,maxColorValue=10))
#   legend('bottomright',lty=1,col=c('blue','cyan'),c('Ensemble mean','Ensemble members'), 
#          bty='o', bg='white',box.col='white')
#   
#   plot(year,gminmem[2,1,],ty='l',col=rgb(0,10,10,1,maxColorValue=10),
#        xlab='Year',ylab='500 hPa Geopotential Height [m]', main='62ºN Zonal Minima of Summer GPH500',
#        ylim=c(min(gminmem[2,,]),max(gminmem[2,,])))
#   for (i in 2:30){
#     lines(year,gminmem[2,i,],col=rgb(0,10,10,1,maxColorValue=10))
#   }
#   lines(year,gminmean[,2],col=rgb(0,0,10,10,maxColorValue=10))
#   abline(v=vf,col=rgb(7,7,7,5,maxColorValue=10))   
#   dev.off()
#   
#   
#   # same but based on one monthly minimum instead of seasonal averages 
#   # but still summer and winter separated
#   for (cyr in syr:eyr) {
#     print(cyr)
#     load(file=paste0(prepplotdirmon,'analysis_',cyr,'.Rdata'))
#     post <- which(floor(analysis.anom$lat)==62 & analysis.anom$names=='temp2') 
#     posg <- which(floor(analysis.anom$lat)==62 & analysis.anom$names=='gph500') 
#     # other possible lat==51.294377    #lat==69.946081
#     if (cyr==syr) {
#       tminmean <- apply(analysis.anom$ensmean[post,],2,min)
#       tminmem <- apply(analysis.anom$data[post,,],2:3,min)
#       gminmean <- apply(analysis.anom$ensmean[posg,],2,min)
#       gminmem <- apply(analysis.anom$data[posg,,],2:3,min)
#     } else {
#       tminmean <- rbind(tminmean,apply(analysis.anom$ensmean[post,],2,min))
#       tminmem <- abind(tminmem,apply(analysis.anom$data[post,,],2:3,min),along=3)
#       gminmean <- rbind(gminmean,apply(analysis.anom$ensmean[posg,],2,min))
#       gminmem <- abind(gminmem,apply(analysis.anom$data[posg,,],2:3,min),along=3)
#     }
#   }  
#   # plot minima
#   vf <- c(1600,1650,1660,1783,1809,1815,1835,1883,1886,1902,1912,1980,1991)
#   year <- syr:eyr
#   
#   pdf('../figures/laki/t2m_z500_min_ts_mon.pdf',width=8,height=12)
#   par(mfrow=c(4,1),mar=c(4,4,4,1)) #,oma=c(2,2,2,0.2))
#   plot(year,tminmem[1,1,],ty='l',col=rgb(10,5,0,1,maxColorValue=10),
#        xlab='Year',ylab='Temperature [ºC]', main='62ºN Zonal Minima of Winter Temperature',
#        ylim=c(min(tminmem[1,,]),max(tminmem[1,,])))
#   for (i in 2:30){
#     lines(year,tminmem[1,i,],col=rgb(10,5,0,1,maxColorValue=10))
#   }
#   lines(year,tminmean[,1],col=rgb(10,0,0,10,maxColorValue=10))
#   abline(v=vf,col=rgb(7,7,7,5,maxColorValue=10))
#   legend('bottomright',lty=1,col=c('red','orange'),c('Ensemble mean','Ensemble members'), 
#          bty='o', bg='white',box.col='white')
#   
#   plot(year,tminmem[2,1,],ty='l',col=rgb(10,5,0,1,maxColorValue=10),
#        xlab='Year',ylab='Temperature [ºC]', main='62ºN Zonal Minima of Summer Temperature',
#        ylim=c(min(tminmem[2,,]),max(tminmem[2,,])))
#   for (i in 2:30){
#     lines(year,tminmem[2,i,],col=rgb(10,5,0,1,maxColorValue=10))
#   }
#   lines(year,tminmean[,2],col=rgb(10,0,0,10,maxColorValue=10))
#   abline(v=vf,col=rgb(7,7,7,5,maxColorValue=10))    
#   
#   plot(year,gminmem[1,1,],ty='l',col=rgb(0,10,10,1,maxColorValue=10),
#        xlab='Year',ylab='500 hPa Geopotential Height [m]', main='62ºN Zonal Minima of Winter GPH500',
#        ylim=c(min(gminmem[1,,]),max(gminmem[1,,])))
#   for (i in 2:30){
#     lines(year,gminmem[1,i,],col=rgb(0,10,10,1,maxColorValue=10))
#   }
#   lines(year,gminmean[,1],col=rgb(0,0,10,10,maxColorValue=10))
#   abline(v=vf,col=rgb(7,7,7,5,maxColorValue=10))
#   legend('bottomright',lty=1,col=c('blue','cyan'),c('Ensemble mean','Ensemble members'), 
#          bty='o', bg='white',box.col='white')
#   
#   plot(year,gminmem[2,1,],ty='l',col=rgb(0,10,10,1,maxColorValue=10),
#        xlab='Year',ylab='500 hPa Geopotential Height [m]', main='62ºN Zonal Minima of Summer GPH500',
#        ylim=c(min(gminmem[2,,]),max(gminmem[2,,])))
#   for (i in 2:30){
#     lines(year,gminmem[2,i,],col=rgb(0,10,10,1,maxColorValue=10))
#   }
#   lines(year,gminmean[,2],col=rgb(0,0,10,10,maxColorValue=10))
#   abline(v=vf,col=rgb(7,7,7,5,maxColorValue=10))   
#   dev.off()
#   
} #end iceland


