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

syr=1951 #1902 #1941
eyr=2004 #2000 #2003 #1970


# syrtot and eyrtot are only used for the total 400 yr indices time series 
#syrtot=1902 #set to the same syr and eyr of the prepplots script (default 1602)
#eyrtot=1905 #(default 2000) 



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
} else if (user=="joerg") {
  workdir='/scratch3/joerg/projects/reuse/git/'
  #} else if (user=="lucaf") {
  #  workdir='/scratch3/lucaf/reuse/reuse_git/'
  #} else if (user == "nevin"){
  #  workdir = '/scratch3/nevin/reuse_climcal/reuse_git/'
} else {
  stop("Unknown user!")
}
dataextdir='/mnt/climstor/giub/EKF400/'
dataintdir=paste0(workdir,'../data/')
setwd(workdir)

source('EnSRF_switches.R')
source('EnSRF_functions.R')

if (pseudo_prox) {
  if (temporal_postproc) {
    dir.create(paste0("../data/prepplot/EKF400_",version,"_",expname))
    dir.create(paste0("../data/prepplot/EKF400_",version,"_",expname,'/prepplot_annual/'))
  }
  prepplotdir=paste0("../data/prepplot/EKF400_",version,"_",expname,'/prepplot_annual/')
  nseas=1
} else if (monthly_out) {
  if (temporal_postproc) {
    dir.create(paste0("../data/prepplot/EKF400_",version,"_",expname))
    dir.create(paste0("../data/prepplot/EKF400_",version,"_",expname,'/prepplot_monthly/'))
  }
  prepplotdir=paste0("../data/prepplot/EKF400_",version,"_",expname,'/prepplot_monthly/')
  nseas=12
} else {
  if (temporal_postproc) {
    dir.create(paste0("../data/prepplot/EKF400_",version,"_",expname))
    dir.create(paste0("../data/prepplot/EKF400_",version,"_",expname,'/prepplot_seasonal/'))
    #dir.create("/mnt/climstor/REUSE/prepplot")
    #dir.create(paste0("/mnt/climstor/REUSE/prepplot/EKF400_",version,"_",expname))
    #dir.create(paste0("/mnt/climstor/REUSE/prepplot/EKF400_",version,"_",expname,"/prepplot_seasonal/"))
  }
  prepplotdir=paste0("../data/prepplot/EKF400_",version,"_",expname,'/prepplot_seasonal/') 
  #prepplotdir=paste0("/mnt/climstor/REUSE/prepplot/EKF400_",version,"_",expname,'/prepplot_seasonal/')
  nseas=2
}
# VV commented it 2019 Oct
#if (mergetime_fields) {dir.create(paste0("../data/image/EKF400_",version,"_",expname))}
#dir.create(paste0('../data/indices/EKF400_",version,"_',expname))

# if (ind_anom) {
#   if (land_only) {
#     if (monthly_out) {
#       filenameext <- paste0('_anom_landonly_mon_') 
#     } else {
#       filenameext <- paste0('_anom_landonly_seas_') 
#     }
#   } else {
#     if (monthly_out) {
#       filenameext <- paste0('_anom_landocean_mon_') 
#     } else {
#       filenameext <- paste0('_anom_landocean_seas_') 
#     }
#   }
# } else {
#   if (land_only) {
#     if (monthly_out) {
#       filenameext <- paste0('_abs_landonly_mon_') 
#     } else {
#       filenameext <- paste0('_abs_landonly_seas_') 
#     }
#   } else {
#     if (monthly_out) {
#       filenameext <- paste0('_abs_landocean_mon_') 
#     } else {
#       filenameext <- paste0('_abs_landocean_seas_') 
#     }
#   }
# }
# print(filenameext)



#####################################################################################

ptm1 <- proc.time()

if (temporal_postproc) {
  print('temporal_postproc')
  # Reads yearly data (output of data script) and changes them to either monthly (abs, oct-sept) or seasonal (means, oct-mar, apr-sept).  Furthermore the indices (compare clac_indices in functions) are calculated,
  # tps_only can be set to T, even when it was set to F before in the data script: like that the complete data is shortened to tps. 
  # In the end the processed data is again saved for each year.
  for (cyr in syr:eyr) {
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
    if ((cyr > syr_twentycr) & (cyr <=eyr_twentycr)& vali_twentycr) {
      twentycr_vali=T             
    } else {
      twentycr_vali=F 
    }
    
    if ((cyr > syr_recon) & (cyr <=eyr_recon) & vali_recon) {
      recon_vali=T           # seasonal luterbacher, pauling, kuettel recons (1750-1999)
    } else {
      recon_vali=F
    }
    if (cyr %% 10 == 0) {
      print(paste('calc year',cyr))
    }
    #print(cyr)
    #print(paste("recon_vali=",recon_vali))
    #print(paste("vali=",vali))
    if (every2grid){
      load(file=paste0('../data/analysis/EKF400_',version,'_',expname,'/analysis_',cyr,'_2ndgrid.Rdata'))
      #load(file=paste0('/mnt/climstor/REUSE/EKF400_',version,'_',expname,'/analysis_',cyr,'_2ndgrid.Rdata'))
      } else {  
      rm(validate,analysis,echam)
      load(file=paste0('../data/analysis/EKF400_',version,'_',expname,'/analysis_',cyr,'.Rdata'))
    }
    if (tps_only&length(unique(echam.abs$names))!=3) {
      echam.abs<-convert_to_tps_only(echam.abs)
      echam.anom<-convert_to_tps_only(echam.anom)
      analysis.abs<-convert_to_tps_only(analysis.abs)
      analysis.anom<-convert_to_tps_only(analysis.anom)
      if (vali){
        
        valiname = names(validate)
        validate_init <- validate
        validate_all <- list()
        
        l=0
        for (v in valiname){  ## for multiple vali data sets
          l=l+1
          #print(v)
          
          validate<-validate_init[[v]]
          validate<-convert_to_tps_only(validate)
          validate_all[[l]] <-validate
        }
        names(validate_all)<-valiname
        validate <- validate_all
      }
    }
    
    
    if (anomaly_assim){
      if (!landcorr) {
        echam <- echam.abs
        analysis <- analysis.abs
        echam.abs <- NULL
        analysis.abs <- NULL  
      } else {
        landcorrected.abs = landcorrected.anom
        landcorrected.abs$data = landcorrected.anom$data + landcorrected.clim$data
        landcorrected.abs$ensmean = landcorrected.anom$data + landcorrected.clim$data
        echam <- landcorrected.abs
        analysis <- analysis.abs
        analysis$ensmean = analysis.abs$data # this is not a nice way cause I make R calculate the same things twice
        echam.abs <- NULL
        analysis.abs <- NULL  
        echam.anom = landcorrected.anom
        echam.anom$ensmean = landcorrected.anom$data
        analysis.anom$ensmean = analysis.anom$data
      }
    }
    if (write_coor) {
      #print('write_coor')
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
    
    
    
    # convert data back to old format for plotting and analysis
    if (sixmonstatevector) {
      #print('convert 6-months state vector')
      if (monthly_out) {
        #print(paste('seasons =',s))
        s=12 # set back to 12 months for plotting
        
        
        tmptime <- seq(cyr-1,(cyr+1),by=1/12)
        echam<-convert_to_monthly(echam)
        echam.anom<-convert_to_monthly(echam.anom)
        analysis<-convert_to_monthly(analysis)
        analysis.anom<-convert_to_monthly(analysis.anom)
        
        if (vali) {
          if (!recon_vali) {
            
            valiname = names(validate)
            validate_init <- validate
            validate_all <- list()
            l=0
            for (v in valiname){  ## for multiple vali data sets
              l=l+1
              #print(v)
              
              validate<-validate_init[[v]]
              validate<-convert_to_monthly(validate)
              validate_all[[l]] <-validate
            }
            names(validate_all)<-valiname
            validate <- validate_all
          }
        }
      } else { # seasonal output, probably summer/winter averages
        s=length(season)
        #print(paste('seasons =',s))
        etmp <- array(echam$data,c((dim(echam$data)[1]/6),6,dim(echam$data)[2],
                                   dim(echam$data)[3]))
        echam$data <- array(NA,c(dim(etmp)[1],dim(etmp)[3:4]))
        for (ensmem in 1:(dim(etmp)[4])) {
          #print(paste('ECHAM member',ensmem))
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
        etmp <- NULL
        etmp2 <- NULL
        
        atmp <- array(analysis$data,c((dim(analysis$data)[1]/6),6,dim(analysis$data)[2],
                                      dim(analysis$data)[3]))  
        analysis$data <- array(NA,c(dim(atmp)[1],dim(atmp)[3:4]))
        for (ensmem in 1:(dim(atmp)[4])) {
          #print(paste('Analysis member',ensmem))
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
        atmp <- NULL
        atmp2 <- NULL
        
        # same for anomalies
        etmp <- array(echam.anom$data,c((dim(echam.anom$data)[1]/6),6,dim(echam.anom$data)[2],
                                        dim(echam.anom$data)[3]))
        echam.anom$data <- array(NA,c(dim(etmp)[1],dim(etmp)[3:4]))
        for (ensmem in 1:(dim(etmp)[4])) {
          #print(paste('ECHAM anomaly member',ensmem))
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
        
        etmp <- NULL
        
        atmp <- array(analysis.anom$data,c((dim(analysis.anom$data)[1]/6),6,dim(analysis.anom$data)[2],
                                           dim(analysis.anom$data)[3])) 
        analysis.anom$data <- array(NA,c(dim(atmp)[1],dim(atmp)[3:4]))
        for (ensmem in 1:(dim(atmp)[4])) {
          #print(paste('Analysis anomaly member',ensmem))
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
            valiname = names(validate)
            validate_init <- validate
            validate_all <- list()
            l=0
            for (v in valiname){  ## for multiple vali data sets
              l=l+1
              #print(v)
              
              validate<-validate_init[[v]]
              
              vtmp <- array(validate$data,c((dim(validate$data)[1]/6),6,dim(validate$data)[2]))  
              validate$data <- apply(vtmp,c(1,3),mean)
              validate$ensmean <- validate$data
              #      validate$names <- validate$names[1:dim(validate$data)[1]]
              validate$names <- rep(unique(validate$names),each=dim(validate$data)[1]/
                                      length(unique(validate$names)))
              validate$lon <- rep(validate$lon,length(unique(validate$names)))
              validate$lat <- rep(validate$lat,length(unique(validate$names)))
              validate$lat <- validate$lat[1:(dim(validate$data)[1])]
              validate$lon <- validate$lon[1:(dim(validate$data)[1])]
              
              validate_all[[l]] <-validate
            }
            names(validate_all)<-valiname
            validate <- validate_all
          }
        }
      }
    }
    #print("transformed 6-mon state vector")
    
    
    
    ########### Indices Calculation ###################
    
    # probably a switch from here for indices calculation would be enough
    echam2.abs <- echam
    analysis2.abs <- analysis
    echam2.anom<-echam.anom
    analysis2.anom<-analysis.anom
    if(vali){
      validate2<-validate
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

      e.abs.ll <- paste(round(echam2.abs$lon,digits=2),round(echam2.abs$lat,digits=2))
      landpos <- match(landll,e.abs.ll)
      landpos_init <- landpos
      landposv<-landpos
      ngrid <- dim(echam2.abs$data)[1]/length(unique(echam2.abs$names))
      nrep <-  length(unique(echam2.abs$names))

      for (i in 1:(nrep-1)) {
        landpos <- c(landpos,(landpos_init+ngrid*i))
      }


      echam2.anom$data <- echam.anom$data[landpos,,]
      echam2.anom$ensmean <- echam.anom$ensmean[landpos,]
      echam2.anom$lon <- echam.anom$lon[landpos]
      echam2.anom$lat <- echam.anom$lat[landpos]
      echam2.anom$names <- echam.anom$names[landpos]
      echam2.abs$data <- echam$data[landpos,,]
      echam2.abs$ensmean <- echam$ensmean[landpos,]
      echam2.abs$lon <- echam$lon[landpos]
      echam2.abs$lat <- echam$lat[landpos]
      echam2.abs$names <- echam$names[landpos]
      analysis2.anom$data <- analysis.anom$data[landpos,,]
      analysis2.anom$ensmean <- analysis.anom$ensmean[landpos,]
      analysis2.anom$lon <- analysis.anom$lon[landpos]
      analysis2.anom$lat <- analysis.anom$lat[landpos]
      analysis2.anom$names <- analysis.anom$names[landpos]
      analysis2.abs$data <- analysis$data[landpos,,]
      analysis2.abs$ensmean <- analysis$ensmean[landpos,]
      analysis2.abs$lon <- analysis$lon[landpos]
      analysis2.abs$lat <- analysis$lat[landpos]
      analysis2.abs$names <- analysis$names[landpos]
      if (vali) {
        valiname = names(validate2)
        validate2_init <- validate2
        validate2_all <- list()
        landpos2<-landposv
        l=0
        for (v in valiname){  ## for multiple vali data sets
          l=l+1
          #print(v)

          validate2<-validate2_init[[v]]
          nrepv <- length(unique(validate2$names))
          for (i in 1:(nrepv-1)) {
            landposv<-c(landposv,(landpos_init+ngrid*i))
          }
          validate2$data <- validate2$data[landposv,]
          validate2$ensmean<- validate2$ensmean[landposv,]
          validate2$lon<- validate2$lon[landposv]
          validate2$lat<- validate2$lat[landposv]
          validate2$names<- validate2$names[landposv]
          validate2_all[[l]] <-validate2
          landposv<-landpos2
        }
        names(validate2_all)<-valiname
        validate2 <- validate2_all
      }
    }



    if (indices) {
          echam.anom2 <- echam2.anom
    analysis.anom2 <- analysis2.anom
    echam.abs2 <- echam2.abs
    analysis.abs2 <- analysis2.abs
      
    eind<-calc_indices(echam.abs2,"echam2")
    aind<-calc_indices(analysis.abs2,"analysis2")
    eind.anom<-calc_indices(echam.anom2,"echam2")
    aind.anom<-calc_indices(analysis.anom2,"analysis2")
    if (vali) {
      valiname = names(validate2)
      vind_all<-list()
      validate2_init <- validate2
      l=0
      for (v in valiname){  ## for multiple vali data sets
        l=l+1
        validate2<-validate2_init[[v]]
        vind_all[[l]]<-calc_indices(validate2,v)
      }
      names(vind_all)<-valiname
      vind<-vind_all
    }
    } # end indices
    
    # save annual file
    # if (vali) {
    #   save(aind,eind,vind, file=paste0('../data/indices/EKF400_',version,'_,expname,'/indices',filenameext,cyr,'.Rdata'))
    # } else {
    #   save(aind,eind, file=paste0('../data/indices/EKF400_',version,'_',expname,'/indices',filenameext,cyr,'.Rdata'))
    # }
    
    if (vali) {
      if (every2grid) {
        valiname = names(validate)
        validate_init <- validate
        validate_all <- list()
        l=0
        for (v in valiname){  ## for multiple vali data sets
          l=l+1
          #print(v)
          validate<-validate_init[[v]]
          #print(paste("dim validate data:",paste(nrow(validate$data),ncol(validate$data))))
          validate_all[[l]] <-validate
        }
        names(validate_all)<-valiname
        validate <- validate_all
        if (indices) {
          save(analysis,analysis.anom,echam,echam.anom,validate,calibrate,vind,eind,aind,aind.anom,eind.anom,
               file=paste0(prepplotdir,'analysis_',cyr,'_2ndgrid.Rdata'))
        } else {
          save(analysis,analysis.anom,echam,echam.anom,validate,calibrate,file=paste0(prepplotdir,'analysis_',cyr,'_2ndgrid.Rdata'))
        }
      } else {
        #print(paste("dim validate data:",paste(nrow(validate$data),ncol(validate$data))))
        if (indices) {
          save(analysis,analysis.anom,echam,echam.anom,validate,calibrate,vind,eind,aind,aind.anom,eind.anom,
               file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
        } else {
          save(analysis,analysis.anom,echam,echam.anom,validate,calibrate,file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
        }
        
        #save(analysis,analysis.anom,ana_ind,echam,echam.anom,ech_ind,validate,vali_ind,calibrate,
        #     file=paste0(prepplotdir,'/analysis_',cyr,'.Rdata'))
        #paste0('../data/prepplot/analysis_',cyr,'.Rdata'))
      }
    } else {
      if (every2grid) {
        if (indices) {
          save(analysis,analysis.anom,echam,echam.anom,calibrate,aind,eind,aind.anom,eind.anom,
               file=paste0(prepplotdir,'analysis_',cyr,'_2ndgrid.Rdata'))
        } else {
          save(analysis,analysis.anom,echam,echam.anom,calibrate,file=paste0(prepplotdir,'analysis_',cyr,'_2ndgrid.Rdata'))
        }
        
      } else {
        if (indices) {
          save(analysis,analysis.anom,echam,echam.anom,calibrate,aind,eind,aind.anom,eind.anom,
               file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
        } else {
          save(analysis,analysis.anom,echam,echam.anom,calibrate,file=paste0(prepplotdir,'analysis_',cyr,'.Rdata'))
        }
        
        #save(analysis,analysis.anom,ana_ind,echam,echam.anom,ech_ind,calibrate,
        #     file=paste0(prepplotdir,'/analysis_',cyr,'.Rdata'))
        #paste0('../data/prepplot/analysis_',cyr,'.Rdata'))
      }
    }
  } # loop over years
  # tps_only is set to F here because if it was T before it's not needed anymore from now on
  # set tps_only = T manually if in mergetime_fields the not-tps-data should be discarded
  #tps_only=F
    rm(echam.abs2,echam2.abs,echam2.anom,analysis.anom2,echam.anom2,analysis.abs2,analysis2.abs,
     analysis2.anom,validate2,validate2_all,validate2_init,vind_all,validate_all,validate_init)

} #end temporal_postproc















if (write_netcdf) {
  print("write_netcdf")
  dir.create(paste0(dataintdir,"netcdf/",version,"/",expname)) 
  dir.create(paste0(dataintdir,"netcdf/",version,"/",expname,"/CCC400_ensmean")) 
  dir.create(paste0(dataintdir,"netcdf/",version,"/",expname,"/CCC400_ens_mem"))
  dir.create(paste0(dataintdir,"netcdf/",version,"/",expname,"/EKF400_ensmean"))
  dir.create(paste0(dataintdir,"netcdf/",version,"/",expname,"/EKF400_ens_mem"))
  # load stat data network for specific year
  # stat_yr=1904
  # if (every2grid){
  #   load(paste0(prepplotdir,'/analysis_',stat_yr,'_2ndgrid.Rdata'))
  # }else{
  #   load(paste0(prepplotdir,'/analysis_',stat_yr,'.Rdata'))
  # }
  # cali <- calibrate
  for (cyr in syr:eyr) {
    if (cyr %% 10 == 0) {
      print(paste('write netcdf year',cyr))
    }
    
    # load data and make normal calendar year Jan-Dec
    if (cyr==feyr) {stop("last year cannot be created because oct-dec data is missing")}
    if (!pseudo_prox) {
      if (every2grid) {
        load(file=paste0(prepplotdir,'/analysis_',(cyr+1),'_2ndgrid.Rdata')) 
      } else {
        load(file=paste0(prepplotdir,'/analysis_',(cyr+1),'.Rdata')) 
      }
      echam2 <- echam
      analysis2 <- analysis
    }
    if (every2grid) {
      load(file=paste0(prepplotdir,'/analysis_',cyr,'_2ndgrid.Rdata')) 
    } else {
      load(file=paste0(prepplotdir,'/analysis_',cyr,'.Rdata')) 
    }
    #    if (monthly_out) {
    if (!pseudo_prox) {
      echam$data <- abind(echam$data[,4:12,],echam2$data[,1:3,],along=2)
      if (length(dim(echam$data)) == 2) {
        echam$data = array(echam$data,c(dim(echam$data),1))
      }
      echam$ensmean <- abind(echam$ensmean[,4:12],echam2$ensmean[,1:3],along=2)
      analysis$data <- abind(analysis$data[,4:12,],analysis2$data[,1:3,],along=2)
      if (length(dim(analysis$data)) == 2) {
        analysis$data = array(analysis$data,c(dim(analysis$data),1))
      }
      analysis$ensmean <- abind(analysis$ensmean[,4:12],analysis2$ensmean[,1:3],along=2)
    }    
    # write entire analysis to netcdf format
    if (write_netcdf) {
      for (m in 1:31) { 
        #print(m)
        for (v in unique(echam$names)) { # to export all vars and not only temp2
          # v="temp2"
          #print(v)
          #define dimensions once, all variables are on the same grid
          if (v=="temp2") { 
            unilon <- unique(analysis$lon)
            unilon[unilon<0] <- unilon[unilon<0]+360
            unilat <- unique(analysis$lat)
            # define dimensions  once
            lon <- ncdim_def( "longitude", "degrees_east",vals=unilon)
            lat <- ncdim_def( "latitude", "degrees_north", vals=unilat)
            if (pseudo_prox) {
              time <- ncdim_def("time", "days since 1600-01-01 00:00:00", calendar="360_day",
                                vals=(cyr-1600)*360+180) 
            } else {
              time <- ncdim_def("time", "days since 1600-01-01 00:00:00", calendar="360_day",
                                vals=seq(from=((cyr-1600)*360+15),to=((cyr-1600)*360+345),by=30))
            }
            lev_temp <- ncdim_def("level_temp", units="m and Pa", vals=c(2,500))
            lev_wind <- ncdim_def("pressure_level_wind", units="Pa", vals=c(850,200))
            lev_gph  <- ncdim_def("pressure_level_gph", units="Pa", vals=c(500,100))
          }
          outpos <- which(echam$names==v)
          if (m==dim(echam$data)[3]+1) {
            out_ech <- echam$ensmean[outpos,,drop=F]
          } else {
            out_ech <- echam$data[outpos,,m,drop=F]
          }
          out_ech[is.na(out_ech)] <- -99999 #set missing value to -999
          if (m==dim(analysis$data)[3]+1) {
            out_ana <- analysis$ensmean[outpos,,drop=F]
          } else {
            out_ana <- analysis$data[outpos,,m,drop=F]
          }
          out_ana[is.na(out_ana)]<- -99999 #set missing value to -999
          # create full 3D field with potential NA over ocean
          if (length(unilon)*length(unilat)!=nrow(out_ech)) {
            fullout_ana <- fullout_ech <- array(NA,dim=c(length(unilon),length(unilat),ncol(out_ana)))
            for (t in 1:ncol(out_ech)) {
              for (i in 1:length(unilon)) {
                for (j in 1:length(unilat)) { #print(c(i,j,t))
                  pos <- NA
                  pos <- intersect(which(analysis$lon[outpos]==unilon[i]),
                                   which(analysis$lat[outpos]==unilat[j]))
                  if (!is.na(pos)){
                    fullout_ech[i,j,t] <- out_ech[pos,t]
                    fullout_ana[i,j,t] <- out_ana[pos,t]
                  }
                }
              }
            }
          } else {
            fullout_ech <- array(out_ech,dim=c(length(unilon),length(unilat),ncol(out_ech)))  
            fullout_ana <- array(out_ana,dim=c(length(unilon),length(unilat),ncol(out_ana)))  
          }
          
          if (v=="temp2") {
            fullout_ech=fullout_ech+273.15 # conversion to std units [K]
            fullout_ana=fullout_ana+273.15 
            temp <- ncvar_def(name="air_temperature",longname="air temperature at 2m and 500hPa", 
                              units="K", dim=list(lon,lat,lev_temp,time))
            temp_data_ech <- array(fullout_ech,c(dim(fullout_ech)[1:2],1,dim(fullout_ech)[3]))
            temp_data_ana <- array(fullout_ana,c(dim(fullout_ana)[1:2],1,dim(fullout_ana)[3]))
            fullout_ech <- fullout_ana <- NULL
          }
          if (v=="t500") {
            #if (v=="t850") {
            temp_data_ech <- abind(temp_data_ech,array(fullout_ech,c(dim(fullout_ech)[1:2],
                                                                     1,dim(fullout_ech)[3])),along=3)
            temp_data_ana <- abind(temp_data_ana,array(fullout_ana,c(dim(fullout_ana)[1:2],
                                                                     1,dim(fullout_ana)[3])),along=3)
            fullout_ech <- fullout_ana <- NULL
          }
          if (v=="precip") {
            precip <- ncvar_def(name="precipitation_amount",longname="precipitation amount (rain and snow)", 
                                units="kg m-2", dim=list(lon,lat,time))
            precip_data_ech <- fullout_ech
            precip_data_ana <- fullout_ana
            fullout_ech <- fullout_ana <- NULL
          }
          if (v=="slp") {
            fullout_ech=fullout_ech*100 # conversion to std units [Pa]
            fullout_ana=fullout_ana*100
            slp <- ncvar_def(name="air_pressure_at_sea_level",longname="air pressure at sea level", 
                             units="Pa", dim=list(lon,lat,time))
            slp_data_ech <- fullout_ech
            slp_data_ana <- fullout_ana
            fullout_ech <- fullout_ana <- NULL
          }
          if (v=="gph500") {
            gph <- ncvar_def(name="geopotential_height",longname="geopotential height", 
                             units="m", dim=list(lon,lat,lev_gph,time))
            gph_data_ech <- array(fullout_ech,c(dim(fullout_ech)[1:2],1,dim(fullout_ech)[3]))
            gph_data_ana <- array(fullout_ana,c(dim(fullout_ana)[1:2],1,dim(fullout_ana)[3]))
            fullout_ech <- fullout_ana <- NULL
          }
          if (v=="gph100") {
            gph_data_ech <- abind(gph_data_ech,array(fullout_ech,c(dim(fullout_ech)[1:2],
                                                                   1,dim(fullout_ech)[3])),along=3)
            gph_data_ana <- abind(gph_data_ana,array(fullout_ana,c(dim(fullout_ana)[1:2],
                                                                   1,dim(fullout_ana)[3])),along=3)
            fullout_ech <- fullout_ana <- NULL
          }
          if (v=="u850") {
            uw <- ncvar_def(name="eastward_wind",longname="eastward wind", 
                            units="m s-1", dim=list(lon,lat,lev_wind,time))
            uw_data_ech <- array(fullout_ech,c(dim(fullout_ech)[1:2],1,dim(fullout_ech)[3]))
            uw_data_ana <- array(fullout_ana,c(dim(fullout_ana)[1:2],1,dim(fullout_ana)[3]))
            fullout_ech <- fullout_ana <- NULL
          }
          if (v=="u200") {
            uw_data_ech <- abind(uw_data_ech,array(fullout_ech,c(dim(fullout_ech)[1:2],
                                                                 1,dim(fullout_ech)[3])),along=3)
            uw_data_ana <- abind(uw_data_ana,array(fullout_ana,c(dim(fullout_ana)[1:2],
                                                                 1,dim(fullout_ana)[3])),along=3)
            fullout_ech <- fullout_ana <- NULL
          }
          if (v=="v850") {
            vw <- ncvar_def(name="northward_wind",longname="northward wind", 
                            units="m s-1", dim=list(lon,lat,lev_wind,time))
            vw_data_ech <- array(fullout_ech,c(dim(fullout_ech)[1:2],1,dim(fullout_ech)[3]))
            vw_data_ana <- array(fullout_ana,c(dim(fullout_ana)[1:2],1,dim(fullout_ana)[3]))
            fullout_ech <- fullout_ana <- NULL
          }
          if (v=="v200") {
            vw_data_ech <- abind(vw_data_ech,array(fullout_ech,c(dim(fullout_ech)[1:2],
                                                                 1,dim(fullout_ech)[3])),along=3)
            vw_data_ana <- abind(vw_data_ana,array(fullout_ana,c(dim(fullout_ana)[1:2],
                                                                 1,dim(fullout_ana)[3])),along=3)
            fullout_ech <- fullout_ana <- NULL
          }
          if (v=="omega500") {
            omega <- ncvar_def(name="lagrangian_tendency_of_air_pressure",
                               longname="vertical velocity alias omega at 500 hPa", 
                               units="Pa s-1", dim=list(lon,lat,time))
            omega_data_ech <- fullout_ech
            omega_data_ana <- fullout_ana
            fullout_ech <- fullout_ana <- NULL
          }
        } # end variables loop  
        
        if (m==dim(echam$data)[3]+1) {
          outfile_ech <- nc_create(filename=paste0(dataintdir,"netcdf/",version,"/",expname,"/CCC400_ensmean/CCC400_ensmean_",
                                                   cyr,"_",version,".nc"), vars=list(temp,precip,slp,gph,uw,vw,omega))
          outfile_ana <- nc_create(filename=paste0(dataintdir,"netcdf/",version,"/",expname,"/EKF400_ensmean/EKF400_ensmean_",
                                                   cyr,"_",version,".nc"), vars=list(temp,precip,slp,gph,uw,vw,omega))
        } else {
          outfile_ech <- nc_create(filename=paste0(dataintdir,"netcdf/",version,"/",expname,"/CCC400_ens_mem/CCC400_ens_mem_",
                                                   m,"_",cyr,"_",version,".nc"), vars=list(temp,precip,slp,gph,uw,vw,omega))
          outfile_ana <- nc_create(filename=paste0(dataintdir,"netcdf/",version,"/",expname,"/EKF400_ens_mem/EKF400_ens_mem_",
                                                   m,"_",cyr,"_",version,".nc"), vars=list(temp,precip,slp,gph,uw,vw,omega))
        }
        
        # Write data to the NetCDF file
        ncvar_put(outfile_ech, temp, temp_data_ech) 
        ncvar_put(outfile_ana, temp, temp_data_ana) 
        ncvar_put(outfile_ech, precip, precip_data_ech)
        ncvar_put(outfile_ana, precip, precip_data_ana)
        ncvar_put(outfile_ech, slp, slp_data_ech)
        ncvar_put(outfile_ana, slp, slp_data_ana)
        ncvar_put(outfile_ech, gph, gph_data_ech)
        ncvar_put(outfile_ana, gph, gph_data_ana)
        ncvar_put(outfile_ech, uw, uw_data_ech)
        ncvar_put(outfile_ana, uw, uw_data_ana)
        ncvar_put(outfile_ech, vw, vw_data_ech)
        ncvar_put(outfile_ana, vw, vw_data_ana)
        ncvar_put(outfile_ech, omega, omega_data_ech)
        ncvar_put(outfile_ana, omega, omega_data_ana)
        
        ## add global attributes
        ncatt_put(outfile_ech,0,"title","CCC400")
        ncatt_put(outfile_ech,0,"institution","University of Bern")
        ncatt_put(outfile_ech,0,"source","ECHAM5.4")
        ncatt_put(outfile_ech,0,"references","Franke et al., 2017")
        ncatt_put(outfile_ech,0,"Conventions","CF-1.6")
        ncatt_put(outfile_ana,0,"title","EKF400")
        ncatt_put(outfile_ana,0,"institution","University of Bern")
        ncatt_put(outfile_ana,0,"source","ECHAM5.4")
        ncatt_put(outfile_ana,0,"references","Franke et al., 2017")
        ncatt_put(outfile_ana,0,"Conventions","CF-1.6")
        
        ncatt_put(outfile_ech,temp,"_FillValue",-99999,prec="single") 
        ncatt_put(outfile_ech,precip,"_FillValue",-99999,prec="single")
        ncatt_put(outfile_ech,slp,"_FillValue",-99999,prec="single")
        ncatt_put(outfile_ech,gph,"_FillValue",-99999,prec="single")
        ncatt_put(outfile_ech,uw,"_FillValue",-99999,prec="single")
        ncatt_put(outfile_ech,vw,"_FillValue",-99999,prec="single")
        ncatt_put(outfile_ech,omega,"_FillValue",-99999,prec="single")
        
        ncatt_put(outfile_ana,temp,"_FillValue",-99999,prec="single")  
        ncatt_put(outfile_ana,precip,"_FillValue",-99999,prec="single")
        ncatt_put(outfile_ana,slp,"_FillValue",-99999,prec="single")
        ncatt_put(outfile_ana,gph,"_FillValue",-99999,prec="single")
        ncatt_put(outfile_ana,uw,"_FillValue",-99999,prec="single")
        ncatt_put(outfile_ana,vw,"_FillValue",-99999,prec="single")
        ncatt_put(outfile_ana,omega,"_FillValue",-99999,prec="single")
        
        ncatt_put(outfile_ech,temp,"cell_methods","time: mean within months")  
        ncatt_put(outfile_ech,precip,"cell_methods","time: sum within months")
        ncatt_put(outfile_ech,slp,"cell_methods","time: mean within months")
        ncatt_put(outfile_ech,gph,"cell_methods","time: mean within months")
        ncatt_put(outfile_ech,uw,"cell_methods","time: mean within months")
        ncatt_put(outfile_ech,vw,"cell_methods","time: mean within months")
        ncatt_put(outfile_ech,omega,"cell_methods","time: mean within months")
        
        ncatt_put(outfile_ana,temp,"cell_methods","time: mean within months")
        ncatt_put(outfile_ana,precip,"cell_methods","time: sum within months")
        ncatt_put(outfile_ana,slp,"cell_methods","time: mean within months")
        ncatt_put(outfile_ana,gph,"cell_methods","time: mean within months")
        ncatt_put(outfile_ana,uw,"cell_methods","time: mean within months")
        ncatt_put(outfile_ana,vw,"cell_methods","time: mean within months")
        ncatt_put(outfile_ana,omega,"cell_methods","time: mean within months")
        
        ncatt_put(outfile_ech,"pressure_level_wind","positive","down")
        ncatt_put(outfile_ana,"pressure_level_wind","positive","down")
        ncatt_put(outfile_ech,"pressure_level_gph","positive","down")
        ncatt_put(outfile_ana,"pressure_level_gph","positive","down")
        ncatt_put(outfile_ech,"longitude","axis","X")
        ncatt_put(outfile_ana,"longitude","axis","X")
        ncatt_put(outfile_ech,"latitude","axis","Y")
        ncatt_put(outfile_ana,"latitude","axis","Y")
        
        # Close your new file to finish writing
        nc_close(outfile_ech)
        nc_close(outfile_ana)
      } # end ens mem loop 
    } # end monthly_out 
  } # end syr:eyr
} # end write_netcdf













# 
# # loads yearly saved indices and merges timesteps (for whole period e.g. 1600-2004 should not generate too huge files)
# if (mergetime_indices){
#   print('mergetime_indices')
#   allvindts<-list()
#   
#   for (cyr in syrtot:(eyrtot)) {
#     #print(paste('year',cyr))
#     if (cyr >= min(c(syr_cru,syr_twentycr,syr_recon)[c(vali_cru, vali_twentycr, vali_recon)]) & cyr<=max(c(eyr_cru,eyr_twentycr,eyr_recon)[c(vali_cru, vali_twentycr, vali_recon)])) {        # if we don't use reconvali, the eyr here should be changed (Error in valiall : object 'valiall' not found) -> but then instead of the eyr we should use cyr
#       vali=T                 # switch off prepplot if no vali data selected
#     } else {
#       vali=F
#     }
#     if ((cyr > syr_cru) & (cyr <=eyr_cru) & vali_cru) {
#       cru_vali=T             # monthly CRU TS3 temp, precip and HADSLP2 gridded instrumentals (1901-2004)
#       #  ind_recon=T            # Stefan's reconstructed indices until 1948 and NCAR reanalysis later added to CRU
#     } else {
#       cru_vali=F 
#       #  ind_recon=F
#     }
#     if ((cyr > syr_twentycr) & (cyr <=eyr_twentycr)& vali_twentycr) {
#       twentycr_vali=T             
#     } else {
#       twentycr_vali=F 
#     }
#     
#     if ((cyr > syr_recon) & (cyr <=eyr_recon) & vali_recon) {
#       recon_vali=T           # seasonal luterbacher, pauling, kuettel recons (1750-1999)
#     } else {
#       recon_vali=F
#     }
#     #print(paste("recon_vali=",recon_vali))
#     #print(paste("vali=",vali))
#     
#     if (every2grid) {
#       load(file=paste0(prepplotdir,'analysis_',cyr,'_2ndgrid.Rdata')) 
#     } else {
#       load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata')) 
#     }
#     rm(analysis,analysis.anom,calibrate,echam,echam.anom,validate)
#     if (tps_only) {
#       eind<-convert_to_tps_only(eind)
#       aind<-convert_to_tps_only(aind)
#       eind.anom<-convert_to_tps_only(eind.anom)
#       aind.anom<-convert_to_tps_only(aind.anom)
#       if (vali){
#         valiname = names(vind)
#         vind_init <- vind
#         vind_all <- list()
#         l=0
#         for (v in valiname){  ## for multiple vali data sets
#           l=l+1
#           #print(v)
#           vind<-vind_init[[v]]
#           vind<-convert_to_tps_only(vind)
#           vind_all[[l]]<-vind
#         }
#         names(vind_all)<-valiname
#         vind<-vind_all
#       }
#     }
#     
#     # merges timesteps into allts variables
#     if (cyr == syrtot) {
#       aind.allts=aind
#       eind.allts=eind
#       c.allts<-eind.anom
#       aind.anom.allts<-aind.anom
#       eind.anom.allts<-eind.anom
#       if (monthly_out) {
#         aind.allts$time=seq(cyr-1,cyr+1,by=(1/nseas))[10:21]
#         eind.allts$time=seq(cyr-1,cyr+1,by=(1/nseas))[10:21]
#         aind.anom.allts$time=seq(cyr-1,cyr+1,by=(1/nseas))[10:21]
#         eind.anom.allts$time=seq(cyr-1,cyr+1,by=(1/nseas))[10:21]
#       } else {
#         aind.allts$time=seq(cyr,cyr+1,by=(1/nseas))[-(nseas+1)]
#         eind.allts$time=seq(cyr,cyr+1,by=(1/nseas))[-(nseas+1)]
#         aind.anom.allts$time=seq(cyr,cyr+1,by=(1/nseas))[-(nseas+1)]
#         eind.anom.allts$time=seq(cyr,cyr+1,by=(1/nseas))[-(nseas+1)]
#       }
#     } else {
#       aind.allts$data=abind(aind.allts$data,aind$data,along=2)
#       aind.allts$ensmean=cbind(aind.allts$ensmean,aind$ensmean)
#       eind.allts$data=abind(eind.allts$data,eind$data,along=2)
#       eind.allts$ensmean=cbind(eind.allts$ensmean,eind$ensmean)
#       aind.anom.allts$data=abind(aind.anom.allts$data,aind.anom$data,along=2)
#       aind.anom.allts$ensmean=cbind(aind.anom.allts$ensmean,aind.anom$ensmean)
#       eind.anom.allts$data=abind(eind.anom.allts$data,eind.anom$data,along=2)
#       eind.anom.allts$ensmean=cbind(eind.anom.allts$ensmean,eind.anom$ensmean)
#       
#       if (monthly_out) {
#         aind.allts$time=c(aind.allts$time, seq(cyr-1,cyr+1,by=(1/nseas))[10:21])
#         eind.allts$time=c(eind.allts$time, seq(cyr-1,cyr+1,by=(1/nseas))[10:21])
#         aind.anom.allts$time=c(aind.anom.allts$time, seq(cyr-1,cyr+1,by=(1/nseas))[10:21])
#         eind.anom.allts$time=c(eind.anom.allts$time, seq(cyr-1,cyr+1,by=(1/nseas))[10:21])
#       } else {
#         aind.allts$time=c(aind.allts$time,seq(cyr,cyr+1,by=(1/nseas))[-(nseas+1)])
#         eind.allts$time=c(eind.allts$time,seq(cyr,cyr+1,by=(1/nseas))[-(nseas+1)])
#         aind.anom.allts$time=c(aind.anom.allts$time,seq(cyr,cyr+1,by=(1/nseas))[-(nseas+1)])
#         eind.anom.allts$time=c(eind.anom.allts$time,seq(cyr,cyr+1,by=(1/nseas))[-(nseas+1)])
#       }
#     }
#     
#     # adapted by Nevin: 
#     # looks complicated but only combines the yearly data into allts variables
#     # control flows are complicated because of possible different time periods of validation sets which are combined together into one list (called validate and vind)
#     # only tested for CRU and 20CR combined
#     
#     # asks whether the cyr is a starting year of a validation dataset
#     if (cyr %in% c(syr_cru,syr_twentycr,syr_recon)[c(vali_cru, vali_twentycr, vali_recon)]|((syrtot==cyr)&syrtot>= min(c(syr_cru,syr_twentycr,syr_recon)[c(vali_cru, vali_twentycr, vali_recon)]))){
#       valiname = names(vind)
#       vind_init<-vind
#       # it still can jump into this part even when one validation set already is running. Here it checks if it's the first time a validation set is started.
#       if (exists("vind.allts")){
#         vind.allts_init<-vind.allts
#       } else {
#         vind.allts_init<-vind
#       }
#       vind_all<-list()
#       vind.allts_all<- list()
#       l=0
#       for (v in valiname){  ## for multiple vali data sets
#         l=l+1
#         #print(v)
#         vind<-vind_init[[v]]
#         # if one vali set already runs and a new one starts now it checks here if v is the new one or not.
#         if (v %in% names(allvindts)) {
#           vind.allts<-vind.allts_init[[v]]
#           vind.allts$data=cbind(vind.allts$data,vind$data)
#           vind.allts$ensmean=cbind(vind.allts$ensmean,vind$ensmean)
#           if (monthly_out) {
#             vind.allts$time=c(vind.allts$time,seq(cyr-1,cyr+1,by=(1/nseas))[10:21])
#           } else {
#             vind.allts$time=c(vind.allts$time,seq(cyr,cyr+1,by=(1/nseas))[-(nseas+1)])
#           }
#         } else {
#           vind.allts<-vind_init[[v]]
#           if (monthly_out) {
#             vind.allts$time=seq(cyr-1,cyr+1,by=(1/nseas))[10:21]
#           } else {
#             vind.allts$time=seq(cyr,cyr+1,by=(1/nseas))[-(nseas+1)]
#           }
#         }
#         vind_all[[l]] <- vind
#         vind.allts_all[[l]] <- vind.allts
#       }
#       names(vind_all) <- valiname
#       vind <- vind_all
#       names(vind.allts_all) <- valiname
#       vind.allts <- vind.allts_all
#       
#     } else if (vali) {
#       valiname <- names(vind)
#       vind_init <- vind
#       vind.allts_init <- vind.allts
#       vind_all <- list()
#       vind.allts_all <- list()
#       l=0
#       for (v in valiname) {  ## for multiple vali data sets
#         l=l+1
#         #print(v)
#         vind <- vind_init[[v]]
#         vind.allts <- vind.allts_init[[v]]
#         vind.allts$data <- cbind(vind.allts$data,vind$data)
#         vind.allts$ensmean <- cbind(vind.allts$ensmean,vind$ensmean)
#         if (monthly_out) {
#           vind.allts$time=c(vind.allts$time,seq(cyr-1,cyr+1,by=(1/nseas))[10:21])
#         } else {
#           vind.allts$time=c(vind.allts$time,seq(cyr,cyr+1,by=(1/nseas))[-(nseas+1)])
#         }
#         # if yearly output of indices is wanted this can be uncommented
#         # if(yearly_out){
#         #   vind.allts<-convert_to_yearly(vind.allts,s)
#         # }
#         vind_all[[l]] <- vind
#         vind.allts_all[[l]] <- vind.allts
#       }
#       names(vind_all) <- valiname
#       vind <- vind_all
#       names(vind.allts_all) <- valiname
#       vind.allts <- vind.allts_all
#     }
#     # stores validate.allts into allvalits for the next year to check then whether the validation set already was running in the previous year
#     if(vali){
#       allvindts <- vind.allts
#     }
#   } # end loop over years
#   
#   # JF 01/2019: calc vind.anom
#   vind.anom <- vind.clim <- vind.allts 
#   l=0
#   for (v in valiname) {  ## for multiple vali data sets
#     l=l+1
#     winpos <- seq(1,ncol(vind.allts[[l]]$data),2)
#     sumpos <- seq(2,ncol(vind.allts[[l]]$data),2)
#     vind.clim[[l]]$ensmean[,winpos] <- vind.clim[[l]]$data[,winpos] <- apply(vind.allts[[l]]$data[,winpos],1,mean,na.rm=T)
#     vind.clim[[l]]$ensmean[,sumpos] <- vind.clim[[l]]$data[,sumpos] <- apply(vind.allts[[l]]$data[,sumpos],1,mean,na.rm=T)
#     vind.anom[[l]]$ensmean <- vind.anom[[l]]$data <- vind.allts[[l]]$data - vind.clim[[l]]$data
#   }
#   vind.tot <- vind.allts
#   vind.anom.tot <- vind.anom
#   aind.tot <- aind.allts
#   eind.tot <- eind.allts
#   eind.anom.tot <- eind.anom.allts
#   aind.anom.tot <- aind.anom.allts
#   if (monthly_out) {
#     stop("insert spos and epos for monthly plots! Not coded yet")
#   } else if (pseudo_prox) {
#     spos <- which(vind.tot[[1]]$time == syrtot)
#     epos <- which(vind.tot[[1]]$time == eyrtot)
#   } else {
#     spos <- which(vind.tot[[1]]$time == syrtot)
#     epos <- which(vind.tot[[1]]$time == paste0(eyrtot,'.5'))
#   }
#   vind <- vind.tot
#   vind.anom <- vind.anom.tot
#   l=0
#   for (v in valiname) {  ## for multiple vali data sets
#     l=l+1
#     vind[[l]]$ensmean <- vind[[l]]$data <- vind.tot[[l]]$data[,spos:epos]
#     vind.anom[[l]]$ensmean <- vind.anom[[l]]$data <- vind.anom.tot[[l]]$data[,spos:epos]
#   }
#   vind.anom$time <- vind$time <- vind.tot$time[spos:epos]
#   
#   eind$ensmean <- eind.allts$ensmean[,spos:epos]
#   eind$data <- eind.tot$data[,spos:epos,]
#   eind$time <- eind.allts$time[spos:epos]
#   eind.anom$ensmean <- eind.anom.allts$ensmean[,spos:epos]
#   eind.anom$data <- eind.anom.allts$data[,spos:epos,]
#   eind.anom$time <- eind.allts$time[spos:epos]
#   
#   aind$ensmean <- aind.allts$ensmean[,spos:epos]
#   aind$data <- aind.tot$data[,spos:epos,]
#   aind$time <- aind.allts$time[spos:epos]
#   aind.anom$ensmean <- aind.anom.allts$ensmean[,spos:epos]
#   aind.anom$data <- aind.anom.allts$data[,spos:epos,]
#   aind.anom$time <- aind.allts$time[spos:epos] 
#   
#   if (every2grid) {
#     if (monthly_out) {
#       save(aind,vind,eind,aind.anom,eind.anom,vind.anom,
#            aind.tot,vind.tot,eind.tot,aind.anom.tot,eind.anom.tot,vind.anom.tot,
#            file=paste0("../data/image/EKF400_",version,"_",expname,"/indices_tot_",syrtot,"-",eyrtot,"_monthly_2ndgrid.Rdata"))
#     } else if (pseudo_prox) {
#       save(aind,vind,eind,aind.anom,eind.anom,vind.anom, 
#            aind.tot,vind.tot,eind.tot,aind.anom.tot,eind.anom.tot,vind.anom.tot,
#            file=paste0("../data/image/EKF400_",version,"_",expname,"/indices_tot_",syrtot,"-",eyrtot,"_annual_2ndgrid.Rdata"))
#     } else {
#       save(aind,vind,eind,aind.anom,eind.anom,vind.anom, 
#            aind.tot,vind.tot,eind.tot,aind.anom.tot,eind.anom.tot,vind.anom.tot,
#            file=paste0("../data/image/EKF400_",version,"_",expname,"/indices_tot_",syrtot,"-",eyrtot,"_seasonal_2ndgrid.Rdata"))
#     }  
#   } else {
#     if (monthly_out) {
#       save(aind,vind,eind,aind.anom,eind.anom,vind.anom,
#            aind.tot,vind.tot,eind.tot,aind.anom.tot,eind.anom.tot,vind.anom.tot,
#            file=paste0("../data/image/EKF400_",version,"_",expname,"/indices_tot_",syrtot,"-",eyrtot,"_monthly.Rdata"))
#     } else if (pseudo_prox) {
#       save(aind,vind,eind,aind.anom,eind.anom,vind.anom,
#            aind.tot,vind.tot,eind.tot,aind.anom.tot,eind.anom.tot,vind.anom.tot,
#            file=paste0("../data/image/EKF400_",version,"_",expname,"/indices_tot_",syrtot,"-",eyrtot,"_annual.Rdata"))  
#     } else {
#       save(aind,vind,eind,aind.anom,eind.anom,vind.anom,
#            aind.tot,vind.tot,eind.tot,aind.anom.tot,eind.anom.tot,vind.anom.tot,
#            file=paste0("../data/image/EKF400_",version,"_",expname,"/indices_tot_",syrtot,"-",eyrtot,"_seasonal.Rdata"))
#     }
#   }
#   rm(allvindts,vind.allts,eind.allts,aind.allts)
# } # end mergetime_indices











# loads saved data (yearly) from temporal_postproc part and merges the timesteps for the validation period, 
# such that the statistics can be calculated. 
# Again tps_only can be set to T here even when it was set to F before: like that the complete data is shortened to tps.
# If you have run temporal_postproc successfully a subperiod can be chosen and run from mergetime_fields. 
# After this part validate.anom are calculated and all the data is stored in an image file.
# calibrate.anom is NOT needed because we always use anomaly_assim=T, i.e. calibrate are already anomalies
# (If someone wants to change months to jan-dec look at write_netcdf-part)
# this part serves as a preparation for the validation statistics 
if (mergetime_fields){ 
  print('mergetime_fields')
  # allvalits variable predefined here (beneath) to enable different periods of validation datasets: 
  #   e.g. 20cr starts in 1850 and cru starts
  # in 1901: if syr = 1899 for the first two years validate only contains 20cr and from 1901 it also contains as well cru.
  # (works for vind.allts as well)
  allvalits<-list()
  if(syr<=min(c(syr_cru,syr_twentycr,syr_recon)[c(vali_cru, vali_twentycr, vali_recon)])){
    syr<-min(c(syr_cru,syr_twentycr,syr_recon)[c(vali_cru, vali_twentycr, vali_recon)])+1
  }
  for (cyr in syr:eyr) {
    if (cyr %% 10 == 0) {
      print(paste('load year',cyr))
    }
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
    
    if ((cyr > syr_twentycr) & (cyr <=eyr_twentycr)& vali_twentycr) {
      twentycr_vali=T             
    } else {
      twentycr_vali=F 
    }   
    if ((cyr > syr_recon) & (cyr <=eyr_recon) & vali_recon) {
      recon_vali=T           # seasonal luterbacher, pauling, kuettel recons (1750-1999)
    } else {
      recon_vali=F
    }
    #print(paste("recon_vali=",recon_vali))
    #print(paste("vali=",vali))
    
    if (every2grid) {
      load(file=paste0(prepplotdir,'analysis_',cyr,'_2ndgrid.Rdata')) 
      # VV 2019 Oct: to test old files where no cru_vali or twentyce_vali existed
      # load(paste0("/scratch3/joerg/projects/reuse/data/prepplot/EKF400_v1.3_corr_echam_clim/prepplot_seasonal/analysis_",cyr,".Rdata"))
      # validate$cru_vali = validate
      # validate$data =NULL
      # validate$ensmean=NULL
      # validate$lon = NULL
      # validate$lat =NULL
      # validate$height =NULL
      # validate$lsm.i = NULL
      # validate$names = NULL
      # validate$time = NULL
    } else {
      load(file=paste0(prepplotdir,'analysis_',cyr,'.Rdata')) 
    }
    
    # remove next "if" again: just for 1. pseudoprox stat. offl. exp with transposed calibrate$data bug
    if (pseudo_prox & expname=="DAPS_pseudoprox_stat_offl") {
      calibrate$data=t(calibrate$data[,,drop=F])
    }
    
    cali<-calibrate # 
    if (tps_only) {
      echam<-convert_to_tps_only(echam)
      echam.anom<-convert_to_tps_only(echam.anom)
      analysis<-convert_to_tps_only(analysis)
      analysis.anom<-convert_to_tps_only(analysis.anom)
      # commented JF 01/2018
      if (indices) {
       eind<-convert_to_tps_only(eind)
       aind<-convert_to_tps_only(aind)
       eind.anom<-convert_to_tps_only(eind.anom)
       aind.anom<-convert_to_tps_only(aind.anom)
      }
      if (vali) {
        valiname <- names(validate)
        validate_init <- validate
        #validate.anom_init <- validate.anom
        validate_all <- list()
        #validate.anom_all <- list()
        # commented JF 01/2018
        if (indices) {
         vind_init <- vind
        #vind.anom_init <- vind.anom # anom is commented
         vind_all <- list()
        #vind.anom_all <- list()
        }
        l=0
        for (v in valiname) {  ## for multiple vali data sets
          l=l+1
          #print(v)
          validate <- validate_init[[v]]
          #validate.anom <- validate.anom_init[[v]]
          validate <- convert_to_tps_only(validate)
          #validate.anom <-convert_to_tps_only(validate.anom)
          validate_all[[l]] <- validate
          #validate.anom_all <- validate.anom 
          #JF 08/2019
          #validate.anom<-convert_to_tps_only(validate.anom)
          # commented JF 01/2018
          if (indices) {
           vind<-vind_init[[v]]
          #vind.anom<-vind.anom_init[[v]] # anom is commented
           vind<-convert_to_tps_only(vind)
          #vind.anom<-convert_to_tps_only(vind.anom) # anom is commented
           vind_all[[l]]<-vind
          #vind.anom_all[[l]]<-vind.anom
          }
        }
        names(validate_all)<-valiname
        # commented JF 01/2018
        if (indices) {
         names(vind_all)<-valiname
         vind<-vind_all
        }
        validate <- validate_all
      }
    }
    
    #provisorily commented this because not sure if going back to Jan-Dec is suitable for calc_vali_stat
    # if (monthly_out) {
    #  echam$data <- abind(echam$data[,4:12,],echam2$data[,1:3,],along=2)
    #  echam$ensmean <- abind(echam$ensmean[,4:12],echam2$ensmean[,1:3],along=2)
    #  analysis$data <- abind(analysis$data[,4:12,],analysis2$data[,1:3,],along=2)
    #  analysis$ensmean <- abind(analysis$ensmean[,4:12],analysis2$ensmean[,1:3],along=2)
    #  if (landcorr) {
    #    echam$data = array (echam$data, c(dim(echam$data),1))
    #    analysis$data = array (analysis$data, c(dim(analysis$data),1))
    #  }
    # }
    
    #    lenvar2 <- length(c(which(validate$names=="temp2"), which(validate$names=="precip"), 
    #                        which(validate$names=="slp")))
    
    
    ## if validate=20cr the state vector needs to include more variables than just tps.
    ## it is assumed later on that tps is included in twentycr
    if (vali) {
      if ("twentycr_vali" %in% names(validate)) {
        ##need to change names of variables because the varnames were not right in the import.
        validate[["twentycr_vali"]]$names[which(validate[["twentycr_vali"]]$names=="omega")] <- "omega500"
        validate[["twentycr_vali"]]$names[which(validate[["twentycr_vali"]]$names=="geopoth")] <- "gph500"
        twentycr.var <- unique(validate[["twentycr_vali"]]$names)
        var.tmp<-which(echam$names%in%twentycr.var)
        analysis$data=analysis$data[var.tmp,,]
        analysis$ensmean=analysis$ensmean[var.tmp,]
        analysis$names=analysis$names[var.tmp]
        analysis$lon=analysis$lon[var.tmp]
        analysis$lat=analysis$lat[var.tmp]
        echam$data=echam$data[var.tmp,,]
        echam$ensmean=echam$ensmean[var.tmp,]
        echam$names=echam$names[var.tmp]
        echam$lon=echam$lon[var.tmp]
        echam$lat=echam$lat[var.tmp]
        analysis.anom=analysis.anom
        analysis.anom$data=analysis.anom$data[var.tmp,,]
        analysis.anom$ensmean=analysis.anom$ensmean[var.tmp,]
        analysis.anom$names=analysis.anom$names[var.tmp]
        analysis.anom$lon=analysis.anom$lon[var.tmp]
        analysis.anom$lat=analysis.anom$lat[var.tmp]
        echam.anom=echam.anom
        echam.anom$data=echam.anom$data[var.tmp,,]
        echam.anom$ensmean=echam.anom$ensmean[var.tmp,]
        echam.anom$names=echam.anom$names[var.tmp]
        echam.anom$lon=echam.anom$lon[var.tmp]
        echam.anom$lat=echam.anom$lat[var.tmp]
      } else {
        lenvar2 <- length(c(which(echam$names=="temp2"), which(echam$names=="precip"), 
                            which(echam$names=="slp")))
        #    analysis_noindex=analysis
        analysis$data=analysis$data[1:lenvar2,,,drop=F]
        analysis$ensmean=analysis$ensmean[1:lenvar2,,drop=F]
        analysis$names=analysis$names[1:lenvar2]
        analysis$lon=analysis$lon[1:lenvar2]
        analysis$lat=analysis$lat[1:lenvar2]
        echam=echam
        echam$data=echam$data[1:lenvar2,,,drop=F]
        echam$ensmean=echam$ensmean[1:lenvar2,,drop=F]
        echam$names=echam$names[1:lenvar2]
        echam$lon=echam$lon[1:lenvar2]
        echam$lat=echam$lat[1:lenvar2]
        analysis.anom=analysis.anom
        analysis.anom$data=analysis.anom$data[1:lenvar2,,,drop=F]
        analysis.anom$ensmean=analysis.anom$ensmean[1:lenvar2,,drop=F]
        analysis.anom$names=analysis.anom$names[1:lenvar2]
        analysis.anom$lon=analysis.anom$lon[1:lenvar2]
        analysis.anom$lat=analysis.anom$lat[1:lenvar2]
        echam.anom=echam.anom
        echam.anom$data=echam.anom$data[1:lenvar2,,,drop=F]
        echam.anom$ensmean=echam.anom$ensmean[1:lenvar2,,drop=F]
        echam.anom$names=echam.anom$names[1:lenvar2]
        echam.anom$lon=echam.anom$lon[1:lenvar2]
        echam.anom$lat=echam.anom$lat[1:lenvar2]
      }
    } # end of vali
    
    # merges timesteps into allts variables
    if (cyr == syr) {
      analysis.allts <- analysis
      analysis.anom.allts <- analysis.anom
      # commented JF 01/2018
      if (indices) {
       aind.allts=aind
       aind.anom.allts<-aind.anom
      }
      echam.allts <- echam
      echam.anom.allts <- echam.anom
      # commented JF 01/2018
      if (indices) {
        eind.allts=eind
        eind.anom.allts<-eind.anom
        if(monthly_out){
          aind.allts$time=seq(cyr-1,cyr+1,by=(1/nseas))[10:21]
          eind.allts$time=seq(cyr-1,cyr+1,by=(1/nseas))[10:21]
          aind.anom.allts$time=seq(cyr-1,cyr+1,by=(1/nseas))[10:21]
          eind.anom.allts$time=seq(cyr-1,cyr+1,by=(1/nseas))[10:21]
        }else{
          aind.allts$time=seq(cyr,cyr+1,by=(1/nseas))[-(nseas+1)]
          eind.allts$time=seq(cyr,cyr+1,by=(1/nseas))[-(nseas+1)]
          aind.anom.allts$time=seq(cyr,cyr+1,by=(1/nseas))[-(nseas+1)]
          eind.anom.allts$time=seq(cyr,cyr+1,by=(1/nseas))[-(nseas+1)]
        }
      }
      
      # only instr. calibration data for period with fixed network
      # because proxies not in fixed grid format. Hence number of records
      # and position/row in data/ensmean matrix changes over time
      calibrate.allts <- calibrate
      # if (substring(expname,1,12)=="proxies_only") {
      if ( length(unique(calibrate$sour)) == 1 & unique(calibrate$sour) == "prox") {
        pos <- which(calibrate$sour=="prox")
      } else {
        pos <- which(calibrate$sour=="inst")
      }
      if (length(pos)==0) {
        pos <- which(cali$sour=="inst")
        calibrate.allts$data=array(NA,dim=dim(cali$data[pos,]))
        calibrate.allts$lon=as.matrix(cali$lon[pos])
        calibrate.allts$lat=as.matrix(cali$lat[pos])
        calibrate.allts$names=as.matrix(cali$names[pos])
        calibrate.allts$sour=as.matrix(cali$sour[pos])
      } else {
        calibrate.allts$data=calibrate$data[pos,,drop=F]  
        calibrate.allts$lon=as.matrix(calibrate$lon[pos])
        calibrate.allts$lat=as.matrix(calibrate$lat[pos])
        if (!pseudo_prox) {
          calibrate.allts$names=as.matrix(calibrate$names[pos])
        }
        calibrate.allts$sour=as.matrix(calibrate$sour[pos])
        #      calibrate.anom.allts=calibrate.anom
        #      calibrate.anom.allts$data=calibrate.anom$data[pos,]
        #      calibrate.anom.allts$lon=as.matrix(calibrate.anom$lon[pos])
        #      calibrate.anom.allts$lat=as.matrix(calibrate.anom$lat[pos])
        #      calibrate.anom.allts$names=as.matrix(calibrate.anom$names[pos])
        #      calibrate.anom.allts$sour=as.matrix(calibrate.anom$sour[pos])
      }
    } else {
      analysis.allts$data=abind(analysis.allts$data,analysis$data,along=2)
      analysis.allts$ensmean=cbind(analysis.allts$ensmean,analysis$ensmean)
      analysis.allts$time=c(analysis.allts$time,analysis$time)
      analysis.anom.allts$data=abind(analysis.anom.allts$data,analysis.anom$data,along=2)
      analysis.anom.allts$ensmean=cbind(analysis.anom.allts$ensmean,analysis.anom$ensmean)
      analysis.anom.allts$time=c(analysis.anom.allts$time,analysis.anom$time)
      echam.allts$data=abind(echam.allts$data,echam$data,along=2)
      echam.allts$ensmean=cbind(echam.allts$ensmean,echam$ensmean)
      echam.allts$time=c(echam.allts$time,echam$time)
      echam.anom.allts$data=abind(echam.anom.allts$data,echam.anom$data,along=2)
      echam.anom.allts$ensmean=cbind(echam.anom.allts$ensmean,echam.anom$ensmean)
      echam.anom.allts$time=c(echam.anom.allts$time,echam.anom$time)
      # commented JF 01/2018
      if (indices) {
        eind.allts$data=abind(eind.allts$data,eind$data,along=2)
        eind.allts$ensmean=cbind(eind.allts$ensmean,eind$ensmean)
        eind.anom.allts$data=abind(eind.anom.allts$data,eind.anom$data,along=2)
        eind.anom.allts$ensmean=cbind(eind.anom.allts$ensmean,eind.anom$ensmean)
        aind.allts$data=abind(aind.allts$data,aind$data,along=2)
        aind.allts$ensmean=cbind(aind.allts$ensmean,aind$ensmean)
        aind.anom.allts$data=abind(aind.anom.allts$data,aind.anom$data,along=2)
        aind.anom.allts$ensmean=cbind(aind.anom.allts$ensmean,aind.anom$ensmean)
        if(monthly_out){
          aind.allts$time=c(aind.allts$time, seq(cyr-1,cyr+1,by=(1/nseas))[10:21])
          eind.allts$time=c(eind.allts$time, seq(cyr-1,cyr+1,by=(1/nseas))[10:21])
          aind.anom.allts$time=c(aind.anom.allts$time, seq(cyr-1,cyr+1,by=(1/nseas))[10:21])
          eind.anom.allts$time=c(eind.anom.allts$time, seq(cyr-1,cyr+1,by=(1/nseas))[10:21])
        } else {
          aind.allts$time=c(aind.allts$time,seq(cyr,cyr+1,by=(1/nseas))[-(nseas+1)])
          eind.allts$time=c(eind.allts$time,seq(cyr,cyr+1,by=(1/nseas))[-(nseas+1)])
          aind.anom.allts$time=c(aind.anom.allts$time,seq(cyr,cyr+1,by=(1/nseas))[-(nseas+1)])
          eind.anom.allts$time=c(eind.anom.allts$time,seq(cyr,cyr+1,by=(1/nseas))[-(nseas+1)])
        }
      }
      
      if (all(calibrate$sour=="prox")) {
        pos <- which(calibrate$sour=="prox")
      } else {
        pos <- which(calibrate$sour=="inst")
      }
      if (length(pos)==0){
        pos <- which(cali$sour=="inst")
        calibrate.allts$data=abind(calibrate.allts$data[pos,],
                                   array(NA,dim=dim(cali$data))[pos,],along=2)
        calibrate.allts$time=c(calibrate.allts$time,calibrate$time)
        calibrate.allts$lon=cbind(calibrate.allts$lon[pos,],cali$lon[pos])
        calibrate.allts$lat=cbind(calibrate.allts$lat[pos,],cali$lat[pos])
        if (!pseudo_prox){
          calibrate.allts$names=cbind(calibrate.allts$names[pos,],cali$names[pos])
        }
        calibrate.allts$sour=cbind(calibrate.allts$sour[pos,],cali$sour[pos])
      } else {
        calibrate.allts$data=abind(calibrate.allts$data[pos,],calibrate$data[pos,],along=2)
        calibrate.allts$time=c(calibrate.allts$time,calibrate$time)
        calibrate.allts$lon=cbind(calibrate.allts$lon[pos,],calibrate$lon[pos])
        calibrate.allts$lat=cbind(calibrate.allts$lat[pos,],calibrate$lat[pos])
        if (!pseudo_prox){
          calibrate.allts$names=cbind(calibrate.allts$names[pos,],calibrate$names[pos])
        }
        calibrate.allts$sour=cbind(calibrate.allts$sour[pos,],calibrate$sour[pos])
        #      calibrate.anom.allts$data=abind(calibrate.anom.allts$data[pos,],calibrate.anom$data[pos,],along=2)
        #      calibrate.anom.allts$time=c(calibrate.anom.allts$time,calibrate.anom$time)
        #      calibrate.anom.allts$lon=cbind(calibrate.anom.allts$lon[pos,],calibrate.anom$lon[pos])
        #      calibrate.anom.allts$lat=cbind(calibrate.anom.allts$lat[pos,],calibrate.anom$lat[pos])
        #      calibrate.anom.allts$names=cbind(calibrate.anom.allts$names[pos,],calibrate.anom$names[pos])
        #      calibrate.anom.allts$sour=cbind(calibrate.anom.allts$sour[pos,],calibrate.anom$sour[pos])
      }
    } # "end of" cyr=syr 
    
    # adapted by Nevin: 
    # looks complicated but only combines the yearly data into allts variables
    # control flows are complicated because of possible different time periods of validation 
    # sets which are combined together into one list (called validate and vind)
    # only tested for CRU and 20CR combined
    
    
    # asks whether the cyr is a starting year of a validation dataset
    if (cyr %in% c(syr_cru,syr_twentycr,syr_recon)[c(vali_cru, vali_twentycr, vali_recon)]|
        ((syr==cyr)&syr>= min(c(syr_cru,syr_twentycr,syr_recon)[c(vali_cru, vali_twentycr, vali_recon)]))){
      valiname = names(validate)
      validate_init <- validate
      # commented JF 01/2018
      if (indices) {
      vind_init<-vind
      }
      
      # it still can jump into this part even when one validation set already is running. Here it checks if it's the first time a validation set is started.
      if (exists("validate.allts")){
        validate.allts_init <- validate.allts
        # commented JF 01/2018
        if (indices) {
          vind.allts_init<-vind.allts
        }
      } else {
        validate.allts_init <- validate
        # commented JF 01/2018
        if (indices) {
         vind.allts_init<-vind
        }
      }
      validate_all <- list()
      validate.allts_all <- list()
      # commented JF 01/2018
      if (indices) {
       vind_all<-list()
       vind.allts_all<- list()
      }
      l=0
      for (v in valiname){  ## for multiple vali data sets
        l=l+1
        #print(v)
        validate<-validate_init[[v]]
        # commented JF 01/2018
        if (indices) {
         vind<-vind_init[[v]]
        }
        # if one vali set already runs and a new one starts now it checks here if v is the new one or not.
        if (v %in% names(allvalits)){
          validate.allts <- validate.allts_init[[v]]
          # commented JF 01/2018
          if (indices) {
          vind.allts<-vind.allts_init[[v]]
          validate.allts$data=cbind(validate.allts$data,validate$data)
          validate.allts$ensmean=cbind(validate.allts$ensmean,validate$ensmean)
          validate.allts$time=c(validate.allts$time,validate$time)
          # commented JF 01/2018
          vind.allts$data=cbind(vind.allts$data,vind$data)
          vind.allts$ensmean=cbind(vind.allts$ensmean,vind$ensmean)
          if(monthly_out){
            vind.allts$time=c(vind.allts$time,seq(cyr-1,cyr+1,by=(1/nseas))[10:21])
          }else{
            vind.allts$time=c(vind.allts$time,seq(cyr,cyr+1,by=(1/nseas))[-(nseas+1)])
          }
          }
        } else {
          validate.allts <- validate_init[[v]]
          # commented JF 01/2018
          if (indices) {
          vind.allts<-vind_init[[v]]
          if(monthly_out){
            vind.allts$time=seq(cyr-1,cyr+1,by=(1/nseas))[10:21]
          }else{
            vind.allts$time=seq(cyr,cyr+1,by=(1/nseas))[-(nseas+1)]
          }
          }
        }
        validate_all[[l]] <-validate
        validate.allts_all[[l]] <- validate.allts
        # commented JF 01/2018
        if (indices) {
         vind_all[[l]]<-vind
         vind.allts_all[[l]]<-vind.allts
        }
      }
      names(validate_all)<-valiname
      validate <- validate_all
      names(validate.allts_all)<- valiname
      validate.allts <- validate.allts_all
      # commented JF 01/2018
      if (indices) {
      names(vind_all)<-valiname
      vind<-vind_all
      names(vind.allts_all)<-valiname
      vind.allts<-vind.allts_all
    }
    } else if (vali){
      valiname = names(validate)
      validate_init <- validate
      validate.allts_init <- validate.allts
      validate_all <- list()
      validate.allts_all <- list()
      # commented JF 01/2018
      if (indices) {
       vind_init<-vind
       vind.allts_init<-vind.allts
       vind_all<-list()
       vind.allts_all<- list()
      }
      l=0
      for (v in valiname){  ## for multiple vali data sets
        l=l+1
        #print(v)
        validate<-validate_init[[v]]
        validate.allts <- validate.allts_init[[v]]
        validate.allts$data=cbind(validate.allts$data,validate$data)
        validate.allts$ensmean=cbind(validate.allts$ensmean,validate$ensmean)
        validate.allts$time=c(validate.allts$time,validate$time)
        # commented JF 01/2018
        if (indices) {
          vind<-vind_init[[v]]
          vind.allts<-vind.allts_init[[v]]
          vind.allts$data=cbind(vind.allts$data,vind$data)
          vind.allts$ensmean=cbind(vind.allts$ensmean,vind$ensmean)
          if(monthly_out){
            vind.allts$time=c(vind.allts$time,seq(cyr-1,cyr+1,by=(1/nseas))[10:21])
          }else{
            vind.allts$time=c(vind.allts$time,seq(cyr,cyr+1,by=(1/nseas))[-(nseas+1)])
          }
          vind_all[[l]]<-vind
          vind.allts_all[[l]]<-vind.allts
        }
        validate_all[[l]] <-validate
        validate.allts_all[[l]] <- validate.allts
      }
      names(validate_all)<-valiname
      validate <- validate_all
      names(validate.allts_all)<- valiname
      validate.allts <- validate.allts_all
      # commented JF 01/2018
      if (indices) {
       names(vind_all)<-valiname
       vind<-vind_all
       names(vind.allts_all)<-valiname
       vind.allts<-vind.allts_all
      }
    }
    # stores validate.allts into allvalits for the next year to check then whether the validation set already was running in the previous year
    allvalits<-validate.allts
  } # end of year loop
  
  
  # following part until save.image belongs to mergetime_fields
  echam <- echam.allts
  echam.anom <- echam.anom.allts

  #ech_ind <- ech_ind.allts
  # rm(ech_ind.allts)
  rm(echam.allts,echam.anom.allts) 
  analysis <- analysis.allts
  analysis.anom <- analysis.anom.allts
  # commented JF 01/2018
  if (indices) {
    eind<-eind.allts
    eind.anom<-eind.anom.allts
    aind<-aind.allts
    aind.anom<-aind.anom.allts
    ana_ind <- ana_ind.allts
    rm(ana_ind.allts, eind.allts,aind.anom.allts,eind.anom.allts)
  }
  
  rm(analysis.allts,analysis.anom.allts,aind.allts)
  
  calibrate <- calibrate.allts
  if (vali) {
    validate <- validate.allts
     rm(validate.allts) 
    # commented JF 01/2018
     if (indices) {
       vind<-vind.allts
       #  vali_ind <- vali_ind.allts
       rm(vali_ind.allts, vind.allts)
     }
   
  }
  #print("calc time for a year")
  
  # VV: added 2019 September
  if (monthly_out){ 
    s.plot <- 12 
  } else if (pseudo_prox) {
    s.plot <- 1
  } else {
    s.plot <- 2
  }
  
  # calc validate.anom, validate.clim, calibrate.clim, calibrate.anom, vind.anom
  if (vali) {
    valiname = names(validate)
    validate_init <- validate
    validate_all <- list()
    validate.clim_all <- list()
    validate.anom_all <- list()
    l=0
    for (v in valiname){  ## for multiple vali data sets
      l=l+1
      # print(s): VV: when I printed s in a monthly_out exp it was 2 instead of 12, therefore above s.plot was introduced (which originally was used later in the script)
      validate<-validate_init[[v]]
      validate.clim <- validate
      validate.clim$data <- apply(array(validate$data, c(nrow(validate$data), s.plot, ncol(validate$data)/s.plot)), 1:2, mean, na.rm=T) 
      validate.clim$ensmean<-apply(array(validate$ensmean, c(nrow(validate$ensmean), s.plot, ncol(validate$ensmean)/s.plot)), 1:2, mean, na.rm=T) 
      validate.anom <- validate
      validate.anom$data <- array(validate$data - as.vector(validate.clim$data), c(nrow(validate$data), ncol(validate$data)))
      validate.anom$ensmean <- array(validate$ensmean - as.vector(validate.clim$ensmean), c(nrow(validate$ensmean), ncol(validate$ensmean)))
      validate_all[[l]] <-validate
      validate.clim_all[[l]] <- validate.clim
      validate.anom_all[[l]] <- validate.anom
    }
    names(validate_all)<-valiname
    validate <- validate_all
    names(validate.clim_all)<- valiname
    validate.clim <- validate.clim_all
    names(validate.anom_all)<- valiname
    validate.anom <- validate.anom_all
  }
  
  if (!anomaly_assim) {
    calibrate.clim <- calibrate
    calibrate.clim$data <- apply(array(calibrate$data, c(nrow(calibrate$data), 2, ncol(calibrate$data)/2)), 1:2, mean, na.rm=T)
    calibrate.anom <- calibrate
    calibrate.anom$data <- array(calibrate$data - as.vector(calibrate.clim$data), c(nrow(calibrate$data), ncol(calibrate$data)))
  }
  valiname = names(validate.anom)
  validate.anom_init <- validate.anom
  # commented JF 01/2018
  if (indices) {
   vind.anom_all <- list()
   vind.anom<-list()
  # vind.anom.tot<-list() # vind.anom is commented
  }
  vanom.yrly<-list()
  l=1
  for (v in valiname){
    validate.anom<-validate.anom_init[[v]]
    validate.anom$data<-array(validate.anom$data,c(dim(validate.anom$data)[1],s.plot,dim(validate.anom$data)[2]/s.plot)) 
    validate.anom$ensmean<-array(validate.anom$ensmean,c(dim(validate.anom$ensmean)[1],s.plot,dim(validate.anom$ensmean)[2]/s.plot))
    validate.anom$time<-array(validate.anom$time,c(s.plot,length(validate.anom$time)/s.plot)) 
    # commented JF 01/2018
    valitmp<-validate.anom
    for (i in 1:dim(validate.anom$time)[2]){
      if (dim(validate.anom$data)[2]==1) {
        valitmp$data<- t(t(validate.anom$data[,,i]))
        valitmp$ensmean<-t(t(validate.anom$ensmean[,,i]))
      } else {
        valitmp$data<-validate.anom$data[,,i]
        valitmp$ensmean<-validate.anom$ensmean[,,i]
      }
      valitmp$lon<-validate.anom$lon
      valitmp$lat<-validate.anom$lat
      valitmp$time<-validate.anom$time[,i]
      valitmp$names<-validate.anom$names
      vanom.yrly[[i]]<-valitmp
    }
    # vind.anom_all<-lapply(vanom.yrly,calc_indices,setname=v)
    # vindtmp<-vind.anom_all[[1]]
    # vindtmp$time<-validate.anom$time[,1]
    # for(i in 2:length(vind.anom_all)){
    #   vindtmp$data<-cbind(vindtmp$data,vind.anom_all[[i]]$data)
    #   vindtmp$ensmean<-cbind(vindtmp$ensmean,vind.anom_all[[i]]$ensmean)
    #   vindtmp$time<-c(vindtmp$time,validate.anom$time[,i])
    # }
    # if (!tps_only) {
    #   vind.clim<-apply(vindtmp$data[35:40,],1,mean)
    #   vindtmp$data[35:40]<-vindtmp$data[35:40,]-vind.clim
    #   vindtmp$ensmean[35:40]<-vindtmp$ensmean[35:40,]-vind.clim
    # }
    # vind.anom[[l]]<-vindtmp
    #vind.anom.tot[[l]]<-vindtmp
    l=l+1
  }
  validate.anom<-validate.anom_init
  # commented JF 01/2018
  if (indices) {
  names(vind.anom)<-valiname
  names(vind.anom.tot)<-valiname
  }
  #print("calc anomalies")
  rm(v,valiname,validate_all,validate_init,validate.allts_all,validate.allts_init,validate.anom_all,
     validate.clim_all,validate.anom_init,vanom.yrly)
  
  if (every2grid) {
    if (monthly_out) {
      save.image(file=paste0("../data/image/EKF400_",version,"_",expname,"/prepplot_validation_image_",
                             syr,"-",eyr,"_monthly_2ndgrid.Rdata"))  
    } else if (pseudo_prox) {
      save.image(file=paste0("../data/image/EKF400_",version,"_",expname,"/prepplot_validation_image_",
                             syr,"-",eyr,"_annual_2ndgrid.Rdata"))
    } else  {
      save.image(file=paste0("../data/image/EKF400_",version,"_",expname,"/prepplot_validation_image_",
                             syr,"-",eyr,"_seasonal_2ndgrid.Rdata"))
      #save.image(file=paste0("/mnt/climstor/REUSE/image/EKF400_",version,"_",expname,"/prepplot_validation_image_",
      #                       syr,"-",eyr,"_seasonal_2ndgrid.Rdata"))
    }  
  } else {
    if (monthly_out) {
      save.image(file=paste("../data/image/EKF400_",version,"_",expname,"/prepplot_validation_image_",
                            syr,"-",eyr,"_monthly.Rdata",sep=""))
    } else if (pseudo_prox) {
      save.image(file=paste("../data/image/EKF400_",version,"_",expname,"/prepplot_validation_image_",
                            syr,"-",eyr,"_annual.Rdata",sep=""))  
    } else {
      save.image(file=paste("../data/image/EKF400_",version,"_",expname,"/prepplot_validation_image_",
                            syr,"-",eyr,"_seasonal.Rdata",sep=""))
    }
  }
  print(proc.time() - ptm1)
}# end mergetime_fields









# ------------------------------------------------------------------------------
# Compute validation statistics
#
# the validate which is loaded at (load_image) can consist of mulitple data sets
# (e.g cru and twentycr). They may be of different first dimensions because 
# one has more variables. if one validation set is 20cr then Echam and Analysis (and .anom) 
# have as many variables as 20cr. If twentycr is not one of the sets then Echam and Analysis (and .anom)
# only have tps -> this is done in (mergetime_fields)
# There is a loop in (calc_vali_stat) which goes for each validation set present in validate
# echam and analysis are shortened if needed (such as for cru). Then the statistic image alongside the 
# complete validate (all data sets) are saved. 
# ------------------------------------------------------------------------------
if (load_images){
  print('load_image of timemerged field and indices')
  if(syr<=min(c(syr_cru,syr_twentycr,syr_recon)[c(vali_cru, vali_twentycr, vali_recon)])){
    syr<-min(c(syr_cru,syr_twentycr,syr_recon)[c(vali_cru, vali_twentycr, vali_recon)])+1
  }
  
  if (every2grid) {
    if (monthly_out) {
#      load(file=paste0("../data/image/EKF400_",version,"_",expname,"/indices_tot_",syrtot,"-",eyrtot,"_monthly_2ndgrid.Rdata"))
      load(file=paste("../data/image/EKF400_",version,"_",expname,"/prepplot_validation_image_",syr,"-",eyr,
                      "_monthly_2ndgrid.Rdata",sep=""))  
    } else if (pseudo_prox) {
#      load(file=paste0("../data/image/EKF400_",version,"_",expname,"/indices_tot_",syrtot,"-",eyrtot,"_annual_2ndgrid.Rdata"))
      load(file=paste("../data/image/EKF400_",version,"_",expname,"/prepplot_validation_image_",syr,"-",eyr,
                      "_annual_2ndgrid.Rdata",sep=""))
    } else {
#      load(file=paste0("../data/image/EKF400_",version,"_",expname,"/indices_tot_",syrtot,"-",eyrtot,"_seasonal_2ndgrid.Rdata"))
      load(file=paste("../data/image/EKF400_",version,"_",expname,"/prepplot_validation_image_",syr,"-",eyr,
                      "_seasonal_2ndgrid.Rdata",sep=""))
    }  
  } else {
    if (monthly_out) {
#      load(file=paste0("../data/image/EKF400_",version,"_",expname,"/indices_tot_",syrtot,"-",eyrtot,"_monthly.Rdata"))
      load(file=paste("../data/image/EKF400_",version,"_",expname,"/prepplot_validation_image_",syr,"-",eyr,
                      "_monthly.Rdata",sep=""))
    } else if (pseudo_prox) {
#      load(file=paste0("../data/image/EKF400_",version,"_",expname,"/indices_tot_",syrtot,"-",eyrtot,"_annual.Rdata"))
      load(file=paste("../data/image/EKF400_",version,"_",expname,"/prepplot_validation_image_",syr,"-",eyr,
                      "_annual.Rdata",sep=""))
    } else {
#      load(file=paste0("../data/image/EKF400_",version,"_",expname,"/indices_tot_",syr,"-",eyr,"_seasonal.Rdata"))
      load(file=paste("../data/image/EKF400_",version,"_",expname,"/prepplot_validation_image_",syr,"-",eyr,
                      "_seasonal.Rdata",sep=""))
    }
  }
} # end load images # only needed if not calculated right before  













if (calc_vali_stat) {  
  print('calc_vali_stat')  
  if (!"cru_vali" %in% names(validate) & !"twentycr_vali" %in% names(validate)) {
    validate_new<-list()
    validate_new[[1]]<-validate 
    names(validate_new)<-"cru_vali"
    validate<-validate_new
  }
  if (!"cru_vali" %in% names(validate.anom) & !"twentycr_vali" %in% names(validate.anom)) {
    validate_new<-list()
    validate_new[[1]]<-validate.anom 
    names(validate_new)<-"cru_vali"
    validate.anom <- validate_new
  }
  if (!"cru_vali" %in% names(validate.clim) & !"twentycr_vali" %in% names(validate.clim)) {
    validate_new<-list()
    validate_new[[1]]<-validate.clim
    names(validate_new)<-"cru_vali"
    validate.clim<-validate_new
  }
  if (load_images) {
    print("ATTENTION: EnSRF_switches loaded again because they may differ in loaded image")
    source('EnSRF_switches.R') # set switches again because they may differ in loaded image  
  }
  rm(s,v,valiname,validate_all,validate_init,validate.allts_all,validate.allts_init,validate.anom_all,validate.clim_all)
  
  if (monthly_out){ 
    s.plot <- 12 
  } else if (pseudo_prox) {
    s.plot <- 1
  } else {
    s.plot <- 2
  }
  
  if(yearly_out){
    echam<-convert_to_yearly(echam,s.plot)
    analysis<-convert_to_yearly(analysis,s.plot)
    calibrate<-convert_to_yearly(calibrate,s.plot)
    echam.anom<-convert_to_yearly(echam.anom,s.plot)
    analysis.anom<-convert_to_yearly(analysis.anom,s.plot)
    if (!anomaly_assim) {
      calibrate.anom<-convert_to_yearly(calibrate.anom,s.plot)
    }
    if (indices) {
    aind<-convert_to_yearly(aind,s.plot)
    eind<-convert_to_yearly(eind,s.plot)
    vind.allts_all<- list()
    vind_init<-vind
    }
    validate.allts_all <- list()
    validate.anom.allts_all<-list()
    validate_init<-validate
    validate.anom_init<-validate.anom
    valiname<-names(validate)
    l=0
    for (v in valiname){  ## for multiple vali data sets
      l=l+1
      #print(v)
      validate<-validate_init[[v]]
      validate.anom<-validate.anom_init[[v]]
      validate<-convert_to_yearly(validate,s.plot)
      validate.anom<-convert_to_yearly(validate.anom,s.plot)
      validate.allts_all[[l]] <-validate
      validate.anom.allts_all[[l]] <- validate.anom
      if (indices) {
        vind<-vind_init[[v]]
        vind<-convert_to_yearly(vind,s.plot)
        vind.allts_all[[l]] <- vind
      }
    }
    names(validate.allts_all)<-valiname
    validate <- validate.allts_all
    names(validate.anom.allts_all)<- valiname
    validate.anom <- validate.anom.allts_all
    if (indices) {
    names(vind.allts_all)<- valiname
    vind <- vind.allts_all
    }
    s.plot=1
  }
  
  ## need to make init file so it can load from this one for every iteration coming below
  echam.init <- echam
  echam.anom.init <- echam.anom
  analysis.init <- analysis
  analysis.anom.init <- analysis.anom
  if (exists("ananomallts")== T){
    ananomallts.init = ananomallts
  }
  calibrate.init <- calibrate
  if (!anomaly_assim) {
    calibrate.anom.init <- calibrate.anom
    calibrate.clim.init <- calibrate.clim
  }
  valiname = names(validate)
  validate_init <- validate
  validate.anom_init <- validate.anom
  if (indices) {
  vind_init<-vind #.anom.tot
  vind.anom_init<-vind.anom
  }
  
  for (v in valiname){  ## for multiple vali data sets ## v -> to vname because image is imported where v might be different
    #print(v)
    ## the variables below may change throughout the loop so i saved it as .init to load them again here
    echam <- echam.init
    echam.anom <- echam.anom.init
    analysis <- analysis.init
    analysis.anom <- analysis.anom.init
    calibrate <- calibrate.init
    if (!anomaly_assim) {
      calibrate.anom <- calibrate.anom.init
      calibrate.clim <- calibrate.clim.init
    }
    validate<-validate_init[[v]]
    validate.anom <- validate.anom_init[[v]]
    if (indices){
      vind <- vind_init[[v]]
      vind.anom <- vind.anom_init[[v]]
      # delete next 3 lines again JF 08/2019
      if (tps_only) {
        vind.anom <- convert_to_tps_only(vind.anom)  
      }
    }
    
    if (nrow(echam$data)!=nrow(validate$data)) { # when 20cr_vali then nrow should be different when v="cru_vali"
      # probably could use here the shorten_func as well
      var.tmp<-which(echam$names%in%validate$names)
      analysis$data=analysis$data[var.tmp,,]
      analysis$ensmean=analysis$ensmean[var.tmp,]
      analysis$names=analysis$names[var.tmp]
      analysis$lon=analysis$lon[var.tmp]
      analysis$lat=analysis$lat[var.tmp]
      echam$data=echam$data[var.tmp,,]
      echam$ensmean=echam$ensmean[var.tmp,]
      echam$names=echam$names[var.tmp]
      echam$lon=echam$lon[var.tmp]
      echam$lat=echam$lat[var.tmp]
      analysis.anom=analysis.anom
      analysis.anom$data=analysis.anom$data[var.tmp,,]
      analysis.anom$ensmean=analysis.anom$ensmean[var.tmp,]
      analysis.anom$names=analysis.anom$names[var.tmp]
      analysis.anom$lon=analysis.anom$lon[var.tmp]
      analysis.anom$lat=analysis.anom$lat[var.tmp]
      echam.anom=echam.anom
      echam.anom$data=echam.anom$data[var.tmp,,]
      echam.anom$ensmean=echam.anom$ensmean[var.tmp,]
      echam.anom$names=echam.anom$names[var.tmp]
      echam.anom$lon=echam.anom$lon[var.tmp]
      echam.anom$lat=echam.anom$lat[var.tmp]
      if (exists("ananomallts") == T) {
        ananomallts=ananomallts
        ananomallts$data=ananomallts$data[var.tmp,,]
        ananomallts$ensmean=ananomallts$ensmean[var.tmp,]
        aananomallts$names=ananomalltsm$names[var.tmp]
        ananomallts$lon=aananomallts$lon[var.tmp]
        ananomallts$lat=ananomallts$lat[var.tmp]
      }
    }

    if (CRPS) {
      ## first option: make a loop
      crps.ana <- array(NA, dim=dim(validate.anom$data))
      crps.ech <- crps.ana
      k=0
      for (i in 1:nrow(analysis.anom$data)) {
        if (i %% 1000 == 0) {
          print(paste('analysis',i,'of',nrow(analysis.anom$data)))
        }
        dressed.ana <- DressEnsemble(analysis.anom$data[i,,])
        dressed.ech <- DressEnsemble(echam.anom$data[i,,])
        
        if (k==0&!any(is.na(dressed.ana$ens))){
          PlotDressedEns(dressed.ana)  ## as an illustration
          k=1
        }
        
        crps.ana[i,] <- DressCrps(dressed.ana,validate.anom$data[i,])
        crps.ech[i,] <- DressCrps(dressed.ech,validate.anom$data[i,])
        
        
      }
      crps.ana <- array(crps.ana,dim=c(nrow(crps.ana),s.plot,ncol(crps.ana)/s.plot))
      crps.ech <- array(crps.ech,dim=c(nrow(crps.ech),s.plot,ncol(crps.ech)/s.plot))
      crps.ana <- apply(crps.ana,c(1,2),mean,na.rm=T)
      crps.ech <- apply(crps.ech,c(1,2),mean,na.rm=T)
    }
    
    rmse <- rmse_fun(analysis, y=validate, seas=s.plot)
    rmse.ech <- rmse_fun(echam, y=validate, seas=s.plot)
    rmse.anom <- rmse_fun(analysis.anom, y=validate.anom, seas=s.plot)
    rmse.ech.anom <- rmse_fun(echam.anom, y=validate.anom, seas=s.plot) 
    echam.clim.anom <- multiyear_seas_average(echam.anom, seas=s.plot)
    echam.clim.anom.data <- echam.clim.anom$data
    echam.clim.anom.time <- echam.clim.anom$time
    for (i in 1:((ncol(echam.anom$data)/s.plot)-1)) {
      echam.clim.anom.data <- cbind(echam.clim.anom.data,echam.clim.anom$data)
      echam.clim.anom.time <- c(echam.clim.anom.time,echam.clim.anom$time)
    }
    echam.clim.anom$data <- echam.clim.anom$ensmean <- echam.clim.anom.data
    echam.clim.anom$time <- echam.clim.anom.time
    rmse.echam.clim.anom <- rmse_fun(echam.clim.anom, y=validate.anom, seas=s.plot)
    echam.clim <- multiyear_seas_average(echam, seas=s.plot)
    echam.clim.data <- echam.clim$data
    echam.clim.time <- echam.clim$time
    for (i in 1:((ncol(echam$data)/s.plot)-1)) {
      echam.clim.data <- cbind(echam.clim.data,echam.clim$data)
      echam.clim.time <- c(echam.clim.time,echam.clim$time)
    }
    echam.clim$data <- echam.clim$ensmean <- echam.clim.data
    echam.clim$time <- echam.clim.time
    rmse.echam.clim <- rmse_fun(echam.clim, y=validate, seas=s.plot)
    RE <- RE_fun(rmse, y=rmse.ech)
    RE.anom <- RE_fun(rmse.anom, y=rmse.ech.anom)
    #  RE.clim.anom <- RE_fun(rmse.anom, y=rmse.clim.anom)
    RE.echam.clim.anom <- RE_fun(rmse.anom, y=rmse.echam.clim.anom)
    #  RE.ech.clim.anom <- RE_fun(rmse.ech.anom, y=rmse.clim.anom)
    RE.echam.clim <- RE_fun(rmse, y=rmse.echam.clim)
    corr <-  corr_fun(analysis, y=validate, seas=s.plot)
    corr.ech <- corr_fun(echam, y=validate, seas=s.plot)
    bias <- bias_fun(analysis, y=validate, seas=s.plot)
    bias.ech <- bias_fun(echam, y=validate, seas=s.plot)
    
    # compute validation statistics on indices
    if (indices) {
    ecorr.ind <- corr_fun(eind.anom, vind.anom, seas=s.plot)
    acorr.ind <- corr_fun(aind.anom, vind.anom, seas=s.plot)
    ermse.ind <- rmse_fun(eind.anom, vind.anom, seas=s.plot)
    armse.ind <- rmse_fun(aind.anom, vind.anom, seas=s.plot)
    ebias.ind <- bias_fun(eind.anom, vind.anom, seas=s.plot)
    abias.ind <- bias_fun(aind.anom, vind.anom, seas=s.plot)
    RE.ind <- RE_fun(armse.ind, y=ermse.ind)
    }
    
    if (!landcorr) {
      data.dim <- dim(echam$data)[c(1,2,2,3)]
      data.dim[2:3] <- c(s.plot,data.dim[3]/s.plot)
      ens.dim <- c(nrow(echam$ensmean), s.plot, ncol(echam$ensmean)/s.plot)
      
      ech.spread <- apply(sqrt(apply(array(echam[['data']] - as.vector(echam[['ensmean']]), 
                                           data.dim)**2, 1:3, mean,na.rm=T)), 1:2, mean, na.rm=T)
      ana.spread <- apply(sqrt(apply(array(analysis[['data']] - as.vector(analysis[['ensmean']]), 
                                           data.dim)**2, 1:3, mean,na.rm=T)), 1:2, mean, na.rm=T)
    }
    
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
    if (!landcorr) {
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
    } # end if landcorr
    
    if (!landcorr) {
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
      # Compute the rank histogram only for summer
      # Added by added by Roni (2018.03)
      summer = seq(2,dim(validate.anom$data)[2],2)
      ereliable.anom.summer <- tapply(apply(echam.anom$data[1:(dim(validate.anom$data)[1]),summer,] > 
                                              as.vector(validate.anom$data[,summer]),1:2,mean), rep(echam.anom$names,
                                                                                                    length=length(validate.anom$data[,summer])), table)
      areliable.anom.summer <- tapply(apply(analysis.anom$data[1:(dim(validate.anom$data)[1]),summer,] > 
                                              as.vector(validate.anom$data[,summer]),1:2,mean), rep(analysis.anom$names,
                                                                                                    length=length(validate.anom$data[,summer])), table) 
      # rank histogram that takes error in instrumental data into account
      #rh = apply(ens + rnorm(length(ens), mean=0, sd=o_sd) < obs, 1, sum)  
      # unter der Annahme das ens eine n x nens matrix und obs ein entsprechender Vektor der Länge n ist
      # from CRU ensemble
      obs_sd <- cbind(apply(obs.spread.win,1,mean),apply(obs.spread.sum,1,mean))
      obs.sd.nona <- obs_sd
      obs.sd.nona[is.na(obs_sd)] <- 0
      rdat <- array(NA,dim(analysis.anom$data[which(analysis.anom$names=="temp2"),,]))
      for (j in 1:dim(analysis.anom$data[,,])[2]) {
        if (is.even(j)) {k=2} else {k=1}
        rdat[,j,] <- rnorm(dim(obs_sd)[1]*dim(analysis.anom$data)[3],mean=0,sd=as.vector(obs.sd.nona[,k]))
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
    } # end if !landcorr
    
    validate <- validate_init
    validate.anom <- validate.anom_init
    if (indices) {
    vind<-vind_init
    vind.anom<-vind.anom_init
    }
    rm(allvalits,cali,cru4_may_sep,cru4_oct_apr,v.sum,v.win,vind.allts_all,
       a.sum, a.win,atala,atmp,atmp2,e.abs.ll,etala,etmp,etmp2,land,landll,
       landpos,landpos2,landposv,rdat,vtmp)
    print(paste("saving",v,"..."))
    if (every2grid) {
      if (monthly_out) {
        save.image(file=paste0("../data/image/EKF400_",version,"_",expname,"/prepplot_",v,"_calc_vali_stat_image_",
                               syr,"-",eyr,"_monthly_2ndgrid.Rdata"))  
        file.remove(paste0("../data/image/EKF400_",version,"_",expname,"/prepplot_validation_image_",
                           syr,"-",eyr,"_monthly_2ndgrid.Rdata"))
      } else if (pseudo_prox) {
        save.image(file=paste0("../data/image/EKF400_",version,"_",expname,"/prepplot_",v,"_calc_vali_stat_image_",
                               syr,"-",eyr,"_annual_2ndgrid.Rdata"))
        file.remove(paste0("../data/image/EKF400_",version,"_",expname,"/prepplot_validation_image_",
                           syr,"-",eyr,"_annual_2ndgrid.Rdata"))
      } else {
         save.image(file=paste0("../data/image/EKF400_",version,"_",expname,"/prepplot_",v,"_calc_vali_stat_image_",
                                syr,"-",eyr,"_seasonal_2ndgrid.Rdata"))
         file.remove(paste0("../data/image/EKF400_",version,"_",expname,"/prepplot_validation_image_",
                            syr,"-",eyr,"_seasonal_2ndgrid.Rdata"))
        #save.image(file=paste0("/mnt/climstor/REUSE/image/EKF400_",version,"_",expname,"/prepplot_",v,"_calc_vali_stat_image_",
        #                      syr,"-",eyr,"_seasonal_2ndgrid.Rdata"))
      }  
    } else {
      if (monthly_out) {
        save.image(file=paste("../data/image/EKF400_",version,"_",expname,"/prepplot_",v,"_calc_vali_stat_image_",
                              syr,"-",eyr,"_monthly.Rdata",sep=""))
        file.remove(paste("../data/image/EKF400_",version,"_",expname,"/prepplot_validation_image_",
                          syr,"-",eyr,"_monthly.Rdata",sep=""))
      } else if (pseudo_prox) {
        save.image(file=paste("../data/image/EKF400_",version,"_",expname,"/prepplot_",v,"_calc_vali_stat_image_",
                              syr,"-",eyr,"_annual.Rdata",sep=""))
        file.remove(paste("../data/image/EKF400_",version,"_",expname,"/prepplot_validation_image_",
                          syr,"-",eyr,"_annual.Rdata",sep=""))
      } else {
        save.image(file=paste("../data/image/EKF400_",version,"_",expname,"/prepplot_",v,"_calc_vali_stat_image_",
                              syr,"-",eyr,"_seasonal.Rdata",sep=""))
        file.remove(paste("../data/image/EKF400_",version,"_",expname,"/prepplot_validation_image_",
                          syr,"-",eyr,"_seasonal.Rdata",sep=""))
      }
    } # end saving
  } # end multiple validation data sets loop
    rm(validate.anom_init,validate_init,analysis.anom.init,analysis.init,
       calibrate.anom.init,calibrate.clim.init,calibrate.init,echam.anom.init,echam.init,
       vind.allts_init,landpos_init) # ,vind.tot.init
} # end calc_vali_stat



  

if (vali_plots) {
  source('EnSRF_plots.R')
  print('ATTENTION: make sure that time periods in EnSRF_plot are set properly by hand!')
}


warnings()
#quit(save='no')

