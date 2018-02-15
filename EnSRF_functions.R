# set direcotries
#echmaskpath <- paste0(dataextdir,'echam/')
echmaskpath <- paste0(dataintdir,'echam/')
echpath <- paste0(dataextdir,'echam/1600-2005/')
echallvarpath <- paste0(dataextdir,'echam_nc_allvar7/')
echanompath <- paste0(dataextdir,'echam_anom/')
echclimpath <- paste0(dataextdir,'echam_clim/')
echsdpath <- paste0(dataextdir,'echam_sd/')
crupath <- paste0(dataintdir,'cru/')
gisspath = '/scratch/veronika/PAGES/climdata/giss/' # can be copied to climstore
#cru4path <- paste0(dataintdir,'cru4/')
reconpath <- paste0(dataintdir,'recon/')
erapath <- paste0(dataintdir,'era/')
ghcntemppath <- paste0(workdir,'../instr/ghcn/temp_v3/')
ghcnprecippath <- paste0(workdir,'../instr/ghcn/precip_v2/')
histalppath <- paste0(dataintdir,'instr/histalp/')
proxypath <- paste0(dataextdir,'assimil_data/proxies/petra/')
mxdpath <- paste0(dataextdir,'assimil_data/proxies/mxd/')
pagespath = '/scratch3/veronika/reuse/'
ntrendpath = '/scratch3/veronika/reuse/'
schweingrpath <- paste0(dataextdir,'assimil_data/proxies/schweingr/')
nceppath <- paste0(dataintdir,'reanalysis/ncep/')
twentycrpath <- paste0(workdir,'../comparison_data/20cr/')
indicespath <- paste0(dataextdir,'vali_data/indices/')


# install.packages("akima")
# install.packages("maps")
# install.packages("mapdata")
# install.packages("TTR")
# install.packages("ncdf4")
# install.packages("RNetCDF")
# install.packages("zoo")
# install.packages("abind")
# install.packages("Matrix")
# install.packages("ff")
# install.packages("caTools")
# install.packages("pracma")
# install.packages("easyVerification")
# install.packages("doParallel")
# install.packages("reshape2")
# install.packages("ggplot2")
# install.packages("grid")
# install.packages("cowplot")
      
suppressMessages(library(akima))         # for interpolation
suppressMessages(library(maps))
suppressMessages(library(mapdata))
suppressMessages(library(TTR))
suppressMessages(library(ncdf4))
suppressMessages(library(RNetCDF))
suppressMessages(library(zoo))
suppressMessages(library(abind))         # to bind multidimensional arrays
suppressMessages(library(Matrix))        # for sparse matrices (not saving value 0)
suppressMessages(library(ff))            # saves large matrices on disk and not im memory
suppressMessages(library(caTools))
suppressMessages(library(pracma))        # for wind vectors
suppressMessages(library(easyVerification)) # spread-error ratio
suppressMessages(library(fields))        # for designer colors and raster maps
suppressMessages(library(doParallel))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))
suppressMessages(library(grid))
suppressMessages(library(cowplot))

#library(pspline)
# first install ncdfUtils of Jonas with all his functions 
# install.packages('.../jonas/ncdfUtils_0.4-12.tar.gz', repos=NULL)
# or alternatively from bash: R CMD INSTALL .../jonas/ncdfUtils*tar.gz
# or load manually:
source('ncdfUtils/R/colfun.R')
source('ncdfUtils/R/frac_in_polygon.R')
source('ncdfUtils/R/geo2rot.R')
source('ncdfUtils/R/ncdf_times.R')
source('ncdfUtils/R/open.ncdf.R')
source('ncdfUtils/R/plot_colourbar.R')
source('ncdfUtils/R/plotmap.R')
source('ncdfUtils/R/plotmap_new.R')
source('ncdfUtils/R/plotmap_rot.R')
source('ncdfUtils/R/set_position.R')
source('ncdfUtils/R/shaded_contour.R')



# NEW read in year by year more efficiently, therefore put ensemble together later
read_echam4 <- function(filehead, path=echallvarpath, xlim=c(-180,180), ylim=c(-90,90), 
                        timlim=c(1600, 2005), small=F, landonly=F, anom=F, clim=F, std=F){
  # read in the land-sea mask of echam
  # mask out sea grid boxes (no variability)
#  if (landonly){
    nc <- nc_open(paste(echmaskpath, 'landseamask.nc', sep='/'))
    nc2 <-nc_open(paste(echmaskpath, 'orography.nc', sep="/"))
    lon <- nc$dim$lon$vals
    lat <- nc$dim$lat$vals
    lon[lon > 180] <- lon[lon > 180] - 360
    loi <- which(lon >= xlim[1] & lon <= xlim[2])
    lai <- which(lat >= ylim[1] & lat <= ylim[2])
    if (small==T) {
      # mulc for reading each 3rd grid cell to avoid memory problems
      mulc <- floor(length(loi)/96)
      loi <- loi[seq(ceiling(mulc/2),length(loi),mulc)]
      lai <- lai[seq(ceiling(mulc/2), length(lai),mulc)]
    }
    lsm <- ncvar_get(nc)[loi, lai]
    alt <- ncvar_get(nc2)[loi, lai]
    nc_close(nc)
    nc_close(nc2)
    lon2 <- lon
    lat2 <- lat
#  }
  
  files <- list.files(path, pattern=paste('^', filehead, sep=''), full.names=T)
  tmp <- list()
  ensmean <- 0
  num <- 0
  for (f in files){
    num <- num + 1
    print(f)
    nc <- nc_open(f)
    lon <- nc$dim$lon$vals
    lat <- nc$dim$lat$vals
    lon[lon > 180] <- lon[lon > 180] - 360
    loi <- which(lon >= xlim[1] & lon <= xlim[2])
    lai <- which(lat >= ylim[1] & lat <= ylim[2])
    latstream <- nc$dim$lat_2$vals
    laistream <- which(latstream >= ylim[1] & latstream <= ylim[2])   
    if (small==T) {
      mulc <- floor(length(loi)/96)
      loi <- loi[seq(ceiling(mulc/2),length(loi),mulc)]
      lai <- lai[seq(ceiling(mulc/2), length(lai),mulc)]
      laistream <- laistream[seq(ceiling(mulc/2), length(laistream),mulc)]
    }
    addc <-0
    year <- as.numeric(format(ncdf_times(nc) + addc, "%Y"))
    month <- as.numeric(format(ncdf_times(nc) + addc, "%m"))
    tim <- year + (month-0.5)/12
    for (t in timlim[1]:(timlim[2]-1)) {
      print(t)
      ti <- which(tim >= t & tim < (t[1]+2))
      #    }
      #    ti <- which(tim >= timlim[1] & tim < (timlim[2]+1))
      outdata <- NULL
      names <- NULL
      # T2m (land only), SLP, Precip, Stream (should be zonal mean), GPH500, GPH100, 
      # u200, omega500, u850, v850, T850
      for (varname in c('temp2', 'precip', 'slp', 'geopoth', 'u', 'v', 'omega',
                        'st', 'stream')){
        #      'geopoth' 1=50000 2=10000
        #      'u', 'v' 1=85000 2=20000
        #       'omega' 1=50000
        #      'st' 1=85000
        #      'stream' 5 lev but just 1 lon (zonal mean) 100000, 85000, 50000, 30000, 20000
        print(varname)
        if (varname %in% names(nc$var)){
          if ((varname=='temp2') || (varname=='precip') || (varname=='slp')) {
            if (length(nc$var[[varname]]$dim) == 3){
              data <- ncvar_get(nc, varname, start=c(1,1,min(ti)),
                                   count=c(-1,-1,length(ti)))[loi, lai,]
            } else if (length(nc$var[[varname]]$dim) == 4){
              data <- ncvar_get(nc, varname, start=c(1,1,1,min(ti)),
                                   count=c(-1,-1,1,length(ti)))[loi, lai,]
              length(nc$var[[varname]]$dim[[3]]$vals)
            }
            if (landonly){
              data <- array(data, c(length(lsm), dim(data)[3]))[lsm > 0.5, ]
            } else {
              data <- array(data, c(dim(data)[1]*dim(data)[2], dim(data)[3]))
            }
            if ((path == echpath) || (path == echallvarpath) || (path == echclimpath) || 
                (path == echsdpath)) {
              if (varname == 'temp2') data <- data - 273.15
              if (varname == 'precip') data <- data * 3600 * 24 * 30
              if (varname == 'slp') data <- data / 100
            }
            print(dim(data))
            outdata <- rbind(outdata, data)
            names <- c(names, rep(varname, nrow(data)))
            #          print(names)
          }
          
          if (varname=='geopoth') {     
            data <- ncvar_get(nc, varname, start=c(1,1,1,min(ti)),
                                 count=c(-1,-1,length(nc$var[[varname]]$dim[[3]]$vals),
                                         length(ti)))[loi, lai,,]
            data <- array(data,c(dim(data)[1]*dim(data)[2], dim(data)[3], dim(data)[4]))
            print(dim(data))
            outdata <- rbind(outdata, data[,1,])
            outdata <- rbind(outdata, data[,2,])
            names <- c(names, c(rep('gph500',nrow(data)),rep('gph100',nrow(data))))
            #          print(names)
          }
          if (varname=='u') {
            data <- ncvar_get(nc, varname, start=c(1,1,1,min(ti)),
                                 count=c(-1,-1,length(nc$var[[varname]]$dim[[3]]$vals),
                                         length(ti)))[loi, lai,,]
            data <- array(data,c(dim(data)[1]*dim(data)[2], dim(data)[3], dim(data)[4]))
            print(dim(data))
            outdata <- rbind(outdata, data[,1,])
            outdata <- rbind(outdata, data[,2,])
            names <- c(names, c(rep('u850',nrow(data)),rep('u200',nrow(data))))
            #         print(names)
          }
          if (varname=='v') {
            data <- ncvar_get(nc, varname, start=c(1,1,1,min(ti)),
                                 count=c(-1,-1,length(nc$var[[varname]]$dim[[3]]$vals),
                                         length(ti)))[loi, lai,,]
            data <- array(data,c(dim(data)[1]*dim(data)[2], dim(data)[3], dim(data)[4]))
            print(dim(data))
            outdata <- rbind(outdata, data[,1,])
            outdata <- rbind(outdata, data[,2,])
            names <- c(names, c(rep('v850',nrow(data)),rep('v200',nrow(data))))
            #          print(names)
          }
          if (varname=='omega') {     
            data <- ncvar_get(nc, varname, start=c(1,1,1,min(ti)),
                                 count=c(-1,-1,1,length(ti)))[loi, lai,]
            data <- array(data, c(dim(data)[1]*dim(data)[2], dim(data)[3]))
            print(dim(data))
            outdata <- rbind(outdata, data)
            names <- c(names, c(rep('omega500',nrow(data))))
            #          print(names)
          }
          if (varname=='st') {     
            data <- ncvar_get(nc, varname, start=c(1,1,1,min(ti)),
                                 count=c(-1,-1,1,length(ti)))[loi, lai,]
            data <- array(data, c(dim(data)[1]*dim(data)[2], dim(data)[3]))
            print(dim(data))
            outdata <- rbind(outdata, data)
            names <- c(names, c(rep('t850',nrow(data))))
            #          print(names)
          }
          if (varname=='stream') {     
            data <- ncvar_get(nc, varname, start=c(1,1,1,min(ti)),
                                 count=c(1,-1,1,length(ti)))[laistream,]
            #          data <- array(data, c(dim(data)[1]*dim(data)[2], dim(data)[3]))
            print(dim(data))
            outdata <- rbind(outdata, data)
            names <- c(names, rep(varname, nrow(data)))
            #          print(names)
          }
        }
      }  
#       if (small) {
#         save(outdata, file=paste("../data/echam/echam_",num,'_',t,"-",(t+1),
#                                  "_3rdgrid.Rdata",sep=""))
#       } else {
#         save(outdata, file=paste("../data/echam/echam_",num,'_',t,"-",(t+1),".Rdata",sep=""))
#       }
      if (anom) {
        if (small) {
          save(outdata, file=paste("../data/echam/echam_anom/echam_anom_",num,"_",t,"-",(t+1),
                                   "_2ndgrid.Rdata",sep=""))
        } else {
          save(outdata, file=paste("../data/echam/echam_anom/echam_anom_",num,"_",t,"-",(t+1),
                                   ".Rdata",sep=""))
        } 
      } else if (clim) {
        if (small) {
          save(outdata, file=paste("../data/echam/echam_clim/echam_clim_",num,"_",t,"-",(t+1),
                                   "_2ndgrid.Rdata",sep=""))
        } else {
          save(outdata, file=paste("../data/echam/echam_clim/echam_clim_",num,"_",t,"-",(t+1),".Rdata",sep=""))
        }
      } else if (std) {
        if (small) {
          save(outdata, file=paste("../data/echam/echam_sd/echam_sd_",num,"_",t,"-",(t+1),
                              "_2ndgrid.Rdata",sep=""))
        } else {
          save(outdata, file=paste("../data/echam/echam_sd/echam_sd_",num,"_",t,"-",(t+1),
                                   ".Rdata",sep=""))
        }
      } else {
        if (small) {
          save(outdata, file=paste("../data/echam/echam_",num,"_",t,"-",(t+1),
                                   "_2ndgrid.Rdata",sep=""))
        } else {
          save(outdata, file=paste("../data/echam/echam_",num,"_",t,"-",(t+1),".Rdata",sep=""))
        }
      }
    } # end time loop
  } # end files loop
  # put 30 members in 1 file and calc ens mean for each time step
# ACHTUNG now starts at timlim[1]+1 = 1601 because some members miss 1600, probably 018  
  for (t in (timlim[1]+1):(timlim[2]-1)) {
    ti <- which(tim >= t & tim < (t[1]+2))
    time <- tim[ti]
    ensmean <- 0
    for (i in 1:length(files)) {
#       if (small) {
#         load(file=paste("../data/echam/echam_",i,'_',t,"-",(t+1),
#                         "_3rdgrid.Rdata",sep=""))
#         tmp[[i]] <- outdata
#       } else {
#         load(file=paste("../data/echam/echam_",i,'_',t,"-",(t+1),
#                         ".Rdata",sep=""))
#         tmp[[i]] <- outdata
#       }
      if (anom) {
        if (small) {
          load(file=paste("../data/echam/echam_anom/echam_anom_",i,'_',t,"-",(t+1),
                                      "_2ndgrid.Rdata",sep=""))
          tmp[[i]] <- outdata
        } else {
          load(file=paste("../data/echam/echam_anom/echam_anom_",i,'_',t,"-",(t+1),".Rdata",sep=""))
          tmp[[i]] <- outdata
        }
      } else if (clim) {
        if (small) {
          load(file=paste("../data/echam/echam_clim/echam_clim_",i,'_',t,"-",(t+1),
                                      "_2ndgrid.Rdata",sep=""))
          tmp[[i]] <- outdata
        } else {
          load(file=paste("../data/echam/echam_clim/echam_clim_",i,'_',t,"-",(t+1),".Rdata",sep=""))
          tmp[[i]] <- outdata
        }
      } else if (std) {
        if (small) {
          load(file=paste("../data/echam/echam_sd/echam_sd_",i,'_',t,"-",(t+1),
                                    "_2ndgrid.Rdata",sep=""))
          tmp[[i]] <- outdata
        } else {
          load(file=paste("../data/echam/echam_sd/echam_sd_",i,'_',t,"-",(t+1),".Rdata",sep=""))
          tmp[[i]] <- outdata
        }
      } else {
        if (small) {
          load(file=paste("../data/echam/echam_",i,'_',t,"-",(t+1),
                                 "_2ndgrid.Rdata",sep=""))
          tmp[[i]] <- outdata
        } else {
          load(file=paste("../data/echam/echam_",i,'_',t,"-",(t+1),".Rdata",sep=""))
          tmp[[i]] <- outdata
        }
      }
      if (i != length(files)) ensmean <- ensmean + outdata  
    }
    ensmean <- ensmean/(length(tmp)) #-1) why -1 ???
    if (landonly) {
      echam <- list(data=array(unlist(tmp), c(dim(tmp[[1]]), length(files))),
                    ensmean=ensmean, lon=rep(lon[loi], length(lai))[lsm > 0.5],
                    lat=rep(lat[lai], each=length(loi))[lsm > 0.5],
                    lonstream=0,
                    latstream=rep(latstream[laistream], each=length(loistream)),
                    height=alt[lsm > 0.5], lsm.i=which(lsm > 0.5),
                    time=time, names=names)
    } else {
      echam <- list(data=array(unlist(tmp), c(dim(tmp[[1]]), length(files))),
                    ensmean=ensmean, lon=rep(lon[loi], length(lai)),
                    lat=rep(lat[lai], each=length(loi)),
                    lonstream=0, latstream=rep(latstream[laistream], each=length(laistream)),
                    height=alt[lsm > 0.5], lsm.i=which(lsm > 0.5),
                    #height=rep(-999,(length(loi)*length(lai))),
                    #lsm.i=rep(-999,(length(loi)*length(lai))), 
                    time=time, names=names)
    }
    if (anom) {
      echam_anom <- echam
      if (small) {
        save(echam_anom, file=paste("../data/echam/echam_anom/echam_anom_",t,"-",(t+1),
                               "_2ndgrid.Rdata",sep=""))
      } else {
        save(echam_anom, file=paste("../data/echam/echam_anom/echam_anom_",t,"-",(t+1),".Rdata",sep=""))
      }
    } else if (clim) {
      echam_clim <- echam
      if (small) {
        save(echam_clim, file=paste("../data/echam/echam_clim/echam_clim_",t,"-",(t+1),
                                    "_2ndgrid.Rdata",sep=""))
      } else {
        save(echam_clim, file=paste("../data/echam/echam_clim/echam_clim_",t,"-",(t+1),".Rdata",sep=""))
      }
    } else if (std) {
      echam_sd <- echam
      if (small) {
        save(echam_sd, file=paste("../data/echam/echam_sd/echam_sd_",t,"-",(t+1),
                                    "_2ndgrid.Rdata",sep=""))
      } else {
        save(echam_sd, file=paste("../data/echam/echam_sd/echam_sd_",t,"-",(t+1),".Rdata",sep=""))
      }
    } else {
      if (small) {
        save(echam, file=paste("../data/echam/echam_",t,"-",(t+1),
                               "_2ndgrid.Rdata",sep=""))
      } else {
        save(echam, file=paste("../data/echam/echam_",t,"-",(t+1),".Rdata",sep=""))
      }
    }
  }
}




# read in the ensemble of echam CCC400 simulations (T P SLP ensmean only)
read_echam_ensmean <- function(filehead, path=echallvarpath, xlim=c(-180,180), ylim=c(-90,90), timlim=c(1901, 1970),small=F){
  
  # read in the land-sea mask of echam
  # mask out sea grid boxes (no variability)
  nc <- nc_open(paste(echmaskpath, 'landseamask.nc', sep='/'))
  nc2 <-nc_open(paste(echmaskpath, 'orography.nc', sep="/"))
  lon <- nc$dim$lon$vals
  lat <- nc$dim$lat$vals
  lon[lon > 180] <- lon[lon > 180] - 360
  loi <- which(lon >= xlim[1] & lon <= xlim[2])
  lai <- which(lat >= ylim[1] & lat <= ylim[2])
  if (small==T) {
    # mulc for reading each 3rd grid cell to avoid memory problems
    mulc <- floor(length(loi)/96)
    loi <- loi[seq(ceiling(mulc/2),length(loi),mulc)]
    lai <- lai[seq(ceiling(mulc/2), length(lai),mulc)]
  }
  lsm <- ncvar_get(nc)[loi, lai]
  alt <- ncvar_get(nc2)[loi,lai]
  nc_close(nc)
  nc_close(nc2)
  lon2 <- lon
  lat2 <- lat
  
  files <- list.files(path, pattern=paste('^', filehead, sep=''), full.names=T)
  tmp <- list()
  ensmean <- 0
  for (f in files){
    print(f)
    nc <- nc_open(f)
    lon <- nc$dim$lon$vals
    lat <- nc$dim$lat$vals
    lon[lon > 180] <- lon[lon > 180] - 360
    loi <- which(lon >= xlim[1] & lon <= xlim[2])
    lai <- which(lat >= ylim[1] & lat <= ylim[2])
    if (small==T) {
      # mulc for reading each 3rd grid cell to avoid memory problems
      loi <- loi[seq(ceiling(mulc/2),length(loi),mulc)]
      lai <- lai[seq(ceiling(mulc/2), length(lai),mulc)]
    }
    addc <-0
    year <- as.numeric(format(ncdf_times(nc) + addc, "%Y"))
    month <- as.numeric(format(ncdf_times(nc) + addc, "%m"))
    tim <- year + (month-0.5)/12
    ti <- which(tim >= timlim[1] & tim < (timlim[2]+1))
    outdata <- NULL
    names <- NULL
    for (varname in c('temp2', 'precip', 'slp')){
      print(varname)
      if (varname %in% names(nc$var)){
        if ((varname=='temp2') || (varname=='precip') || (varname=='slp')) { 
          if (length(nc$var[[varname]]$dim) == 3){
            data <- ncvar_get(nc, varname, start=c(1,1,min(ti)), 
                                 count=c(-1,-1,length(ti)))[loi, lai,]
          } else if (length(nc$var[[varname]]$dim) == 4){
            data <- ncvar_get(nc, varname, start=c(1,1,1,min(ti)), 
                                 count=c(-1,-1,1,length(ti)))[loi, lai,]
          }
          data <- array(data, c(length(lsm), dim(data)[3]))[lsm > 0.5, ]
        }
        if ((path == echpath) || (path == echallvarpath)) {
          if (varname == 'temp2') data <- data - 273.15
          if (varname == 'precip') data <- data * 3600 * 24 * 30
          if (varname == 'slp') data <- data / 100
        }
        outdata <- rbind(outdata, data)
        names <- c(names, rep(varname, nrow(data)))
      }
    }

    time <- tim[ti]
    tmp[[which(files == f)]] <- outdata
    print(dim(outdata))
    if (f != files[length(files)]) ensmean <- ensmean + outdata
  }
  ensmean <- ensmean/(length(tmp)) #-1) why -1 ???
  rm(tmp)
  data <- list(ensmean=ensmean, lon=rep(lon[loi], length(lai))[lsm > 0.5],
               lat=rep(lat[lai], each=length(loi))[lsm > 0.5],
               height=alt[lsm > 0.5], lsm.i=which(lsm > 0.5),
               time=time, names=names)
  invisible(data) 
}
# end read echem temp precip slp ensmean only





# OLD (T, P, SLP, indices version) read in the ensemble of echam simulations
read_echam1 <- function(filehead, path=echpath, xlim=c(-180,180), ylim=c(-90,90), timlim=c(1600, 1630), small=F, landonly=F){

  # read in the land-sea mask of echam
  # mask out sea grid boxes (no variability)
#  if (landonly){
    nc <- nc_open(paste(echmaskpath,'landseamask.nc', sep='/'))
    nc2 <-nc_open(paste(echmaskpath,'orography.nc', sep="/"))
    lon <- nc$dim$lon$vals
    lat <- nc$dim$lat$vals
    lon[lon > 180] <- lon[lon > 180] - 360
    loi <- which(lon >= xlim[1] & lon <= xlim[2])
    lai <- which(lat >= ylim[1] & lat <= ylim[2])
    if (small==T) {
    # mulc for reading each 3rd grid cell to avoid memory problems
      mulc <- floor(length(loi)/96)
      loi <- loi[seq(ceiling(mulc/2),length(loi),mulc)]
      lai <- lai[seq(ceiling(mulc/2), length(lai),mulc)]
    }  
    lsm <- ncvar_get(nc)[loi, lai]
    alt <- ncvar_get(nc2)[loi,lai]
    nc_close(nc)
    nc_close(nc2)
    lon2 <- lon
    lat2 <- lat
#  } 
  files <- list.files(path, pattern=paste('^', filehead, sep=''), full.names=T)
  tmp <- list()
  ensmean <- 0
  for (f in files){
    print(f)
    nc <- nc_open(f)
    lon <- nc$dim$lon$vals
    lat <- nc$dim$lat$vals
    lon[lon > 180] <- lon[lon > 180] - 360
    loi <- which(lon >= xlim[1] & lon <= xlim[2])
    lai <- which(lat >= ylim[1] & lat <= ylim[2])
    if (small==T) {
      mulc <- floor(length(loi)/96)
      loi <- loi[seq(ceiling(mulc/2),length(loi),mulc)]
      lai <- lai[seq(ceiling(mulc/2), length(lai),mulc)]
    }
    # change processing of the time coordinate
    # tim <- seq(1600 + 1/24, by=1/12, length=nc$dim$time$len)
    # ti  <- which(tim >= 1600 & tim < 1631)
    # glitches in time representation
#    if (timlim[1] < 1900){
#      addc <- -15
#    } else {
#      addc <- 3
#    }
    addc <-0
    year <- as.numeric(format(ncdf_times(nc) + addc, "%Y"))
    month <- as.numeric(format(ncdf_times(nc) + addc, "%m"))
    tim <- year + (month-0.5)/12
    ti <- which(tim >= timlim[1] & tim < (timlim[2]+1))
    outdata <- NULL
    names <- NULL
    for (varname in c('temp2', 'precip', 'slp')){
      if (varname %in% names(nc$var)){
        if (length(nc$var[[varname]]$dim) == 3){
          if (length(ti)==1) { # added 2016 July 11 to read cru4
            data <- ncvar_get(nc, varname, start=c(1,1,min(ti)), 
                              count=c(-1,-1,length(ti)))[loi, lai]  
          } else {
            data <- ncvar_get(nc, varname, start=c(1,1,min(ti)), 
                    count=c(-1,-1,length(ti)))[loi, lai,]
          }  
        } else if (length(nc$var[[varname]]$dim) == 4){
          data <- ncvar_get(nc, varname, start=c(1,1,1,min(ti)), 
                    count=c(-1,-1,1,length(ti)))[loi, lai,]
        }
        if (landonly){
          data <- array(data, c(length(lsm), dim(data)[3]))[lsm > 0.5, ]
        } else {
          if (length(ti)==1) {
            data <- array(data, c(dim(data)[1]*dim(data)[2]))
          } else {
            data <- array(data, c(dim(data)[1]*dim(data)[2], dim(data)[3]))
          }
        }
        if ((path == echpath) | (path == echallvarpath) | (path == echallvarpath)) {
          if (varname == 'temp2') data <- data - 273.15
          if (varname == 'precip') data <- data * 3600 * 24 * 30
          if (varname == 'slp') data <- data / 100
        }
        if (path == nceppath) {
          if (varname == 'temp2') data <- data - 273.15
          if (varname == 'precip') data <- data * 3600 * 24 * 30
          if (varname == 'slp') data <- data 
        }
        # compute seasonal averages
        #dat <- array(data[,11:(ncol(data)-2)], c(nrow(data), 6, ncol(data)/6 - 2))
        #data <- apply(dat, c(1,3), mean, na.rm=T)
        outdata <- rbind(outdata, data)
        names <- c(names, rep(varname, nrow(data)))
        #rm(data, dat)
      }
    }
    for (varname in setdiff(names(nc$var), c('temp2', 'precip', 'slp'))){
      data <- ncvar_get(nc, varname)[ti]
       #dat <- array(data[11:(length(data) - 2)], c(6, length(data)/6 - 2))
       #data <- apply(dat, 2, mean, na.rm=T)
      outdata <- rbind(outdata, data)
      names <- c(names, varname)
       #rm(data, dat)
    }
    time <- tim[ti]
    tmp[[which(files == f)]] <- outdata
    print(dim(outdata))
    if (f != files[length(files)]) ensmean <- ensmean + outdata
  }
  ensmean <- ensmean/(length(tmp)) #-1) why -1 ???
  if (landonly) {
    data <- list(data=array(unlist(tmp), c(dim(tmp[[1]]), length(files))),
                 ensmean=ensmean, lon=rep(lon[loi], length(lai))[lsm > 0.5],
                 lat=rep(lat[lai], each=length(loi))[lsm > 0.5], 
                 height=alt[lsm > 0.5], lsm.i=which(lsm > 0.5),
                 time=time, names=names)
  } else {
    data <- list(data=array(unlist(tmp), c(dim(tmp[[1]]), length(files))),
                 ensmean=ensmean, lon=rep(lon[loi], length(lai)),
                 lat=rep(lat[lai], each=length(loi)), 
                 height=alt[lsm > 0.5], lsm.i=which(lsm > 0.5),
                 #height=rep(-999,(length(loi)*length(lai))), 
                 #lsm.i=rep(-999,(length(loi)*length(lai))), 
                 time=time, names=names)
  }
  invisible(data) 
}







read_proxy_mxd <- function(syr,eyr){
  files <- list.files(mxdpath, pattern='.csv', full.names=T)
  mxd <- list()
#  ensmean <- 0
  for (f in files){
    # print(f)
    tmp1 <- as.matrix(read.table(f, header=F, sep=" ", dec = ".", na.strings = "NA",
                                 skip=1))
    if (f == files[1]) {
      tmpmxd <- scale(as.numeric(tmp1[3:dim(tmp1)[1],2]))
      tmpts <- ts(tmpmxd,start=as.numeric(tmp1[3,1]),
                  end=as.numeric(tmp1[dim(tmp1)[1],1]))
#      tmpname <- tmp1[1,2]
      tmplon <- tmp1[1,2]
      tmplat <- tmp1[2,2]
    } else {
      tmpmxd <- scale(as.numeric(tmp1[3:dim(tmp1)[1],2]))
      tmpts <- ts.union(tmpts,ts(tmpmxd,start=as.numeric(tmp1[3,1]),
                                 end=as.numeric(tmp1[dim(tmp1)[1],1])))
#      tmpname <- c(tmpname,tmp1[1,2])
      tmplon <- c(tmplon,tmp1[1,2])
      tmplat <- c(tmplat,tmp1[2,2])
    }
  }
  tmptime <- ts(rep(NA,syr,eyr),start=syr,end=eyr)
  tmpts <- ts.union(tmpts,tmptime)[,1:(dim(tmpts)[2])]
  
  regdata <- as.matrix(window(tmpts,1901,1970))
  mxd$data <- as.matrix(window(tmpts,syr,eyr))
#  mxd$name <- tmpname
  mxd$lon <- as.numeric(tmplon)
  mxd$lat <- as.numeric(tmplat)
  mxd$time <- seq(from=tsp(window(tmpts,syr,eyr))[1],to=tsp(window(tmpts,syr,eyr))[2],
                  by=tsp(window(tmpts,syr,eyr))[3])
  # load CRU for same period and locations
  nc=nc_open(paste(crupath,'/cru_allvar_abs_1901-2004.nc',sep=''), write=F)
  #  print.nc(nc)
#  t1=var.get.nc(nc, "temp2") # for CRU temp
  t1=ncvar_get(nc, "temp2") # for CRU temp
#  p1=var.get.nc(nc, "precip") # for CRU precip
#  lonlist=var.get.nc(nc, "lon") 
  lonlist=ncvar_get(nc, "lon") 
  lonlist[lonlist > 180] <- lonlist[lonlist > 180] - 360
#  latlist=var.get.nc(nc, "lat") 
  latlist=ncvar_get(nc, "lat") 
  nc_close(nc)
  t2 <- t1[,,1:840]
#  p2 <- p1[,,1:840] 
  
  load(paste0(echmaskpath,"echam_1911-70.Rdata"))
  t10=echam1901_70$ensmean[echam1901_70$names=='temp2',] # for echam temp
#  p10=echam1901_70$ensmean[echam1901_70$names=='precip',] # for echam precip
  lonlist2=echam1901_70$lon
  lonlist2[lonlist2 > 180] <- lonlist2[lonlist2 > 180] - 360
  latlist2=echam1901_70$lat
  
  for (i in 1:length(mxd$lon)) {
    k=which(abs(lonlist-mxd$lon[i]+0.001)==min(abs(lonlist-mxd$lon[i]+0.001)))
    l=which(abs(latlist-mxd$lat[i])==min(abs(latlist-mxd$lat[i])))
    t3 <- t2[k,l,]
    # make variable for each month
    t4 <- t(array(t3,c(12,length(t3)/12)))
    unab <- t4[,4:9]
    colnames(unab) <-c("t.first","t.second","t.third","t.fourth","t.fifth","t.sixth")
    
    onlytemp<-regression_months[grep("t.",regression_months,fixed=TRUE)]
    unab <- unab[,match(onlytemp,colnames(unab))]
    results <- lm(regdata[,i]~unab)
    corr <- cor(results$coefficients[1]+results$coefficients[2]*t4[,5]+
                   results$coefficients[3]*t4[,6]+results$coefficients[4]*t4[,7]+
                   results$coefficients[5]*t4[,8],
                   regdata[,i])
    # print(corr)
    if (i==1) { 
      mr <- results$coefficients 
      var_residu <- var(results$residuals)
    } else { 
      mr <- rbind(mr,results$coefficients) 
      var_residu <- c(var_residu,var(results$residuals))
    }
    k=which(abs(lonlist2-mxd$lon[i]+0.001)==min(abs(lonlist2-mxd$lon[i]+0.001)))
    l=which(abs(latlist2-mxd$lat[i])==min(abs(latlist2-mxd$lat[i])))
    if (max(match(k,l,nomatch=-99999))==-99999) {
      k=which(abs(lonlist2-mxd$lon[i]+2.001)==min(abs(lonlist2-mxd$lon[i]+2.001)))
      l=which(abs(latlist2-mxd$lat[i])==min(abs(latlist2-mxd$lat[i])))
    }
    if (max(match(k,l,nomatch=-99999))==-99999) {
      k=which(abs(lonlist2-mxd$lon[i]-2.001)==min(abs(lonlist2-mxd$lon[i]-2.001)))
      l=which(abs(latlist2-mxd$lat[i])==min(abs(latlist2-mxd$lat[i])))
    }
    if (max(match(k,l,nomatch=-99999))==-99999) {
      k=which(abs(lonlist2-mxd$lon[i]+0.001)==min(abs(lonlist2-mxd$lon[i]+0.001)))
      l=which(abs(latlist2-mxd$lat[i]+2)==min(abs(latlist2-mxd$lat[i]+2)))
    }
    if (max(match(k,l,nomatch=-99999))==-99999) {
      k=which(abs(lonlist2-mxd$lon[i]+0.001)==min(abs(lonlist2-mxd$lon[i]+0.001)))
      l=which(abs(latlist2-mxd$lat[i]-2)==min(abs(latlist2-mxd$lat[i]-2)))
    }
    t30 <- t10[l[match(k,l)[!is.na(match(k,l))]],]
    # make variable for each month
    t40 <- t(array(t30,c(12,length(t30)/12)))
    unab2 <- cbind(t40[,5:8]) 
    colnames(unab2) <- c('t5','t6','t7','t8') 
    results2 <- lm(regdata[,i]~unab2)
    corr <- cor(results2$coefficients[1]+results2$coefficients[2]*t4[,5]+
                     results2$coefficients[3]*t4[,6]+results2$coefficients[4]*t4[,7]+
                     results2$coefficients[5]*t4[,8], 
                     regdata[,i])
    if (i==1) { 
      mr2 <- results2$coefficients 
    } else { 
      mr2 <- rbind(mr2,results2$coefficients) 
    }    
  }
  
  mrtmp<-matrix(NA, nrow(mr),13)
  colnames(mrtmp) <-c("(Intercept)","unabt.first","unabt.second","unabt.third","unabt.fourth","unabt.fifth","unabt.sixth","unabp.first","unabp.second","unabp.third","unabp.fourth","unabp.fifth","unabp.sixth")
  
  mrtmp[,match(colnames(mr),colnames(mrtmp))]<-mr
  
  
  mxd$mr <- mrtmp
  # print(mxd$mr[,1]-mr2[,1])
#  mxd$mr[,1] <- mr2[,1]
  mxd$var_residu <- var_residu
  invisible(mxd) 
}







read_proxy_schweingr <- function(syr,eyr){
  files <- list.files(schweingrpath, pattern='.ppd', full.names=T)
  mxd <- list()
  #  ensmean <- 0
  for (f in files){
    # print(f)
    tmp1 <- as.matrix(read.table(f, header=F, strip.white=T, dec = ".", 
                                 stringsAsFactors=F, na.strings = "NA", skip=0))
    if (f == files[1]) {
      tmpmxd <- scale(tmp1[4:dim(tmp1)[1],2])
      print(c(as.numeric(tmp1[4,1]),as.numeric(tmp1[dim(tmp1)[1],1])))
      tmpts <- ts(as.numeric(tmpmxd),start=as.numeric(tmp1[4,1]),
                  end=as.numeric(tmp1[dim(tmp1)[1],1]))
      #      tmpname <- tmp1[1,2]
      tmplon <- tmp1[1,2]
      tmplat <- tmp1[2,2]
    } else {
      tmpmxd <- scale(tmp1[4:dim(tmp1)[1],2])
      # print(c(as.numeric(tmp1[4,1]),as.numeric(tmp1[dim(tmp1)[1],1])))
      tmpts <- ts.union(tmpts,ts(as.numeric(tmpmxd),start=as.numeric(tmp1[4,1]),
                                 end=as.numeric(tmp1[dim(tmp1)[1],1])))
      #      tmpname <- c(tmpname,tmp1[1,2])
      tmplon <- c(tmplon,tmp1[1,2])
      tmplat <- c(tmplat,tmp1[2,2])
    }
  }
  tmptime <- ts(rep(NA,syr,eyr),start=syr,end=eyr)
  tmpts <- ts.union(tmpts,tmptime)[,1:(dim(tmpts)[2])]
  
  regdata <- as.matrix(window(tmpts,1901,1960))
  mxd$data <- as.matrix(window(tmpts,syr,eyr))
  #  mxd$name <- tmpname
  mxd$lon <- as.numeric(tmplon)
  mxd$lat <- as.numeric(tmplat)
  mxd$time <- seq(from=tsp(window(tmpts,syr,eyr))[1],to=tsp(window(tmpts,syr,eyr))[2],
                  by=tsp(window(tmpts,syr,eyr))[3])
  # load CRU for same period and locations
#  nc=open.nc(paste(crupath,'/cru_allvar_abs_1901-2004.nc',sep=''), write=F)
  nc=nc_open(paste(crupath,'/cru_allvar_abs_1901-2004.nc',sep=''), write=F)
  #  print.nc(nc)
#  t1=var.get.nc(nc, "temp2") # for CRU temp
  t1=ncvar_get(nc, "temp2") # for CRU temp
  #  p1=var.get.nc(nc, "precip") # for CRU precip
#  lonlist=var.get.nc(nc, "lon") 
  lonlist=ncvar_get(nc, "lon") 
  lonlist[lonlist > 180] <- lonlist[lonlist > 180] - 360
#  latlist=var.get.nc(nc, "lat") 
  latlist=ncvar_get(nc, "lat") 
  nc_close(nc)
  t2 <- t1[,,1:720]
    
  load(paste0(echmaskpath,"echam_1911-70.Rdata"))
  t10=echam1901_70$ensmean[echam1901_70$names=='temp2',1:720] # for echam temp
  #  p10=echam1901_70$ensmean[echam1901_70$names=='precip',] # for echam precip
  lonlist2=echam1901_70$lon
  lonlist2[lonlist2 > 180] <- lonlist2[lonlist2 > 180] - 360
  latlist2=echam1901_70$lat
  
  for (i in 1:length(mxd$lon)) {
    # print(i)
    k=which(abs(lonlist-mxd$lon[i]+0.001)==min(abs(lonlist-mxd$lon[i]+0.001)))
    l=which(abs(latlist-mxd$lat[i])==min(abs(latlist-mxd$lat[i])))
    t3 <- t2[k,l,]
    # make variable for each month
    
    t4 <- t(array(t3,c(12,length(t3)/12)))
    unab <- t4[,4:9]
    colnames(unab) <-c("t.first","t.second","t.third","t.fourth","t.fifth","t.sixth")
    onlytemp<-regression_months[grep("t.",regression_months,fixed=TRUE)] 
    unab <- unab[,match(onlytemp,colnames(unab))]

    results <- lm(regdata[,i]~unab)
    corr <- cor(results$coefficients[1]+results$coefficients[2]*t4[,5]+
                  results$coefficients[3]*t4[,6]+results$coefficients[4]*t4[,7]+
                  results$coefficients[5]*t4[,8],
                regdata[,i])
    # print(corr)
#    tmp <- as.vector(results$coefficients[1]+results$coefficients[2]*t4[,5]+
#              results$coefficients[3]*t4[,6]+results$coefficients[4]*t4[,7]+
#              results$coefficients[5]*t4[,8])
#    plot(ts(tmp,1901,freq=1),ty='l',col='red')
#    lines((regdata[,i]),ty='l')
    if (i==1) { 
      mr <- results$coefficients 
      var_residu <- var(results$residuals)
    } else { 
      mr <- rbind(mr,results$coefficients) 
      var_residu <- c(var_residu,var(results$residuals))
    }
    k=which(abs(lonlist2-mxd$lon[i]+0.001)==min(abs(lonlist2-mxd$lon[i]+0.001)))
    l=which(abs(latlist2-mxd$lat[i])==min(abs(latlist2-mxd$lat[i])))
    if (max(match(k,l,nomatch=-99999))==-99999) {
      k=which(abs(lonlist2-mxd$lon[i]+2.001)==min(abs(lonlist2-mxd$lon[i]+2.001)))
      l=which(abs(latlist2-mxd$lat[i])==min(abs(latlist2-mxd$lat[i])))
    }
    if (max(match(k,l,nomatch=-99999))==-99999) {
      k=which(abs(lonlist2-mxd$lon[i]-2.001)==min(abs(lonlist2-mxd$lon[i]-2.001)))
      l=which(abs(latlist2-mxd$lat[i])==min(abs(latlist2-mxd$lat[i])))
    }
    if (max(match(k,l,nomatch=-99999))==-99999) {
      k=which(abs(lonlist2-mxd$lon[i]+0.001)==min(abs(lonlist2-mxd$lon[i]+0.001)))
      l=which(abs(latlist2-mxd$lat[i]+2)==min(abs(latlist2-mxd$lat[i]+2)))
    }
    if (max(match(k,l,nomatch=-99999))==-99999) {
      k=which(abs(lonlist2-mxd$lon[i]+0.001)==min(abs(lonlist2-mxd$lon[i]+0.001)))
      l=which(abs(latlist2-mxd$lat[i]-2)==min(abs(latlist2-mxd$lat[i]-2)))
    }
    if (max(match(k,l,nomatch=-99999))==-99999) {
      k=which(abs(lonlist2-mxd$lon[i]+2.001)==min(abs(lonlist2-mxd$lon[i]+2.001)))
      l=which(abs(latlist2-mxd$lat[i]+2)==min(abs(latlist2-mxd$lat[i]+2)))
    }
    if (max(match(k,l,nomatch=-99999))==-99999) {
      k=which(abs(lonlist2-mxd$lon[i]-2.001)==min(abs(lonlist2-mxd$lon[i]-2.001)))
      l=which(abs(latlist2-mxd$lat[i]+2)==min(abs(latlist2-mxd$lat[i]+2)))
    }
    if (max(match(k,l,nomatch=-99999))==-99999) {
      k=which(abs(lonlist2-mxd$lon[i]+2.001)==min(abs(lonlist2-mxd$lon[i]+2.001)))
      l=which(abs(latlist2-mxd$lat[i]-2)==min(abs(latlist2-mxd$lat[i]-2)))
    }
    if (max(match(k,l,nomatch=-99999))==-99999) {
      k=which(abs(lonlist2-mxd$lon[i]-2.001)==min(abs(lonlist2-mxd$lon[i]-2.001)))
      l=which(abs(latlist2-mxd$lat[i]-2)==min(abs(latlist2-mxd$lat[i]-2)))
    }
    t30 <- t10[l[match(k,l)[!is.na(match(k,l))]],]
    # make variable for each month
    t40 <- t(array(t30,c(12,length(t30)/12)))
    unab2 <- cbind(t40[,5:8]) 
    colnames(unab2) <- c('t5','t6','t7','t8') 
    results2 <- lm(regdata[,i]~unab2)
    corr <- cor(results2$coefficients[1]+results2$coefficients[2]*t4[,5]+
                  results2$coefficients[3]*t4[,6]+results2$coefficients[4]*t4[,7]+
                  results2$coefficients[5]*t4[,8], 
                regdata[,i])
    if (i==1) { 
      mr2 <- results2$coefficients 
    } else { 
      mr2 <- rbind(mr2,results2$coefficients) 
    }
  }
  
  
  mrtmp<-matrix(NA, nrow(mr),13)
  colnames(mrtmp) <-c("(Intercept)","unabt.first","unabt.second","unabt.third","unabt.fourth","unabt.fifth","unabt.sixth","unabp.first","unabp.second","unabp.third","unabp.fourth","unabp.fifth","unabp.sixth")
  
  mrtmp[,match(colnames(mr),colnames(mrtmp))]<-mr
  
  mxd$mr <- mrtmp # gives 25 possible coefficients
  
  # print(mxd$mr[,1]-mr2[,1])
#  mxd$mr[,1] <- mr2[,1]
  mxd$var_residu <- var_residu
  schweingr <- mxd
  invisible(schweingr) 
}








# read in Petra's TRW tree-ring data
# to use multiple-regression approach for assimilation
# based on monthly data 
# monthly, seasonal or annual plots???
read_proxy2 <- function(syr,eyr){

  load(paste(proxypath,'t35.RData',sep=''))
  load(paste(proxypath,'vslt35.RData',sep=''))
#  data period 1600-1970, cut CRU overlap 1901-1970 
  trw <- t35
  trw$chronologies <- NULL
  t35$chronologies <- scale(t35$chronologies)
  trw$data <- t35$chronologies[302:371,]
  t35$time <- seq(from=1600,to=1600-1+nrow(t35$chronologies))
  ti <- which((t35$time >= syr) & (t35$time <= eyr))
  t35$data <- t35$chronologies[ti,] 
  t35$time <- t35$time[ti]
  t35$lon <- t35$lonlat[1,]
  t35$lat <- t35$lonlat[2,]
  trw$lonlat <- NULL
  t35lon <- t35$lonlat[1,]
  t35lat <- t35$lonlat[2,]
# load CRU for same period and locations
  nc=nc_open(paste(crupath,'/cru_allvar_abs_1901-2004.nc',sep=''), write=F)
#  print.nc(nc)
#  t1=var.get.nc(nc, "temp2") # for CRU temp
#  p1=var.get.nc(nc, "precip") # for CRU precip
#  lonlist=var.get.nc(nc, "lon") 
  t1=ncvar_get(nc, "temp2") # for CRU temp
  p1=ncvar_get(nc, "precip") # for CRU precip
  lonlist=ncvar_get(nc, "lon") 
  lonlist[lonlist > 180] <- lonlist[lonlist > 180] - 360
#  latlist=var.get.nc(nc, "lat") 
  latlist=ncvar_get(nc, "lat") 
  nc_close(nc)
  t2 <- t1[,,1:840]
  p2 <- p1[,,1:840] 
  
  load(paste0(echmaskpath,"echam_1911-70.Rdata"))
  t10=echam1901_70$ensmean[echam1901_70$names=='temp2',] # for echam temp
  p10=echam1901_70$ensmean[echam1901_70$names=='precip',] # for echam precip
  lonlist2=echam1901_70$lon
  lonlist2[lonlist2 > 180] <- lonlist2[lonlist2 > 180] - 360
  latlist2=echam1901_70$lat
  
  for (i in 1:length(t35lon)) {
    k=which(abs(lonlist-t35lon[i]+0.001)==min(abs(lonlist-t35lon[i]+0.001)))
    l=which(abs(latlist-t35lat[i])==min(abs(latlist-t35lat[i])))
    t3 <- t2[k,l,]
    p3 <- p2[k,l,]
    # make variable for each month
    t4 <- t(array(t3,c(12,length(t3)/12)))
    p4 <- t(array(p3,c(12,length(p3)/12)))
    unab <- cbind(t4[,4:9],p4[,4:9]) # from April to September

    colnames(unab) <-c("t.first","t.second","t.third","t.fourth","t.fifth","t.sixth","p.first","p.second","p.third","p.fourth","p.fifth","p.sixth")

    unab <- unab[,match(regression_months,colnames(unab))]
    results <- lm(trw$data[,i]~unab)
    corr <- cor.test((results$coefficients[1]+results$coefficients[2]*t4[,5]+
               results$coefficients[3]*t4[,6]+results$coefficients[4]*t4[,7]+
               results$coefficients[5]*t4[,8]+results$coefficients[6]*p4[,4]+
               results$coefficients[7]*p4[,5]+results$coefficients[8]*p4[,6]),
               trw$data[,i])
    # print(corr[4])
    if (corr[3] < 0.05) {
      if (i==1) { 
        mr <- results$coefficients 
        var_residu <- var(results$residuals)
      } else { 
        mr <- rbind(mr,results$coefficients) 
        var_residu <- c(var_residu,var(results$residuals))
      }
#     print(paste('Proxy:',i))
#     print(summary(results))
#     vslgt <- colMeans(t(array(vslt35$grow_temp[,i],c(12,length(vslt35$grow_temp[,i])/12))))
#     vslgp <- colMeans(t(array(vslt35$grow_moist[,i],c(12,length(vslt35$grow_moist[,i])/12))))
#     plot(vslgt,ty='l',col='red') 
#     par(new=T)
#     plot(vslgp,ty='l',col='blue',axes=F)
#     par(ask=T) 
    

      k=which(abs(lonlist2-t35lon[i]+0.001)==min(abs(lonlist2-t35lon[i]+0.001)))
      l=which(abs(latlist2-t35lat[i])==min(abs(latlist2-t35lat[i])))
      if (max(match(k,l,nomatch=-99999))==-99999) {
        k=which(abs(lonlist2-t35lon[i]+1.001)==min(abs(lonlist2-t35lon[i]+1.001)))
        l=which(abs(latlist2-t35lat[i])==min(abs(latlist2-t35lat[i])))
      } 
      if (max(match(k,l,nomatch=-99999))==-99999) {
        k=which(abs(lonlist2-t35lon[i]-1.001)==min(abs(lonlist2-t35lon[i]-1.001)))
        l=which(abs(latlist2-t35lat[i])==min(abs(latlist2-t35lat[i])))
      }
      if (max(match(k,l,nomatch=-99999))==-99999) {
        k=which(abs(lonlist2-t35lon[i]+0.001)==min(abs(lonlist2-t35lon[i]+0.001)))
        l=which(abs(latlist2-t35lat[i]+1)==min(abs(latlist2-t35lat[i]+1)))
      }
      if (max(match(k,l,nomatch=-99999))==-99999) {
        k=which(abs(lonlist2-t35lon[i]+0.001)==min(abs(lonlist2-t35lon[i]+0.001)))
        l=which(abs(latlist2-t35lat[i]-1)==min(abs(latlist2-t35lat[i]-1)))
      }
      t30 <- t10[l[match(k,l)[!is.na(match(k,l))]],]
      p30 <- p10[l[match(k,l)[!is.na(match(k,l))]],]
      # make variable for each month
      t40 <- t(array(t30,c(12,length(t30)/12)))
      p40 <- t(array(p30,c(12,length(p30)/12)))
      unab2 <- cbind(t40[,5:8],p40[,4:6])
      colnames(unab2) <- c('t5','t6','t7','t8','p4','p5','p6')
      results2 <- lm(trw$data[,i]~unab2)
      if (!exists("mr2")) { 
        mr2 <- results2$coefficients 
      } else { 
        mr2 <- rbind(mr2,results2$coefficients) 
      }    
    } else {
      if (i==1) { 
        mr <- rep(NA, length(results$coefficients)) 
        var_residu <- NA
      } else { 
        mr <- rbind(mr,rep(NA, length(results$coefficients))) 
        var_residu <- c(var_residu,NA)
      }
    }
  }

  mrtmp<-matrix(NA, nrow(mr),13)
  colnames(mrtmp) <-c("(Intercept)","unabt.first","unabt.second","unabt.third","unabt.fourth","unabt.fifth","unabt.sixth","unabp.first","unabp.second","unabp.third","unabp.fourth","unabp.fifth","unabp.sixth")
  mrtmp[,match(colnames(mr),colnames(mrtmp))]<-mr
  
  t35$mr <- mrtmp 
  # print(t35$mr[,1]-mr2[,1])
#  t35$mr[,1] <- mr2[,1]  
  t35$var_residu <- var_residu
  invisible(t35)
}




# PAGES proxies are TEMPERATURE sensitive
# PAGES year is from April to March
# tree: regression calculated for NH from April till September, for SH from October till March
# coral: regression is calculated by using all months (12 reg coeff-s from April till March                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   )
# docu: save docu data as an RData file
# inst: save inst data as an Rdata file
# PATHES SHOULD BE MAKE MORE GENERAL!!!
read_pages = function(fsyr,feyr,archivetype, validate) {
  # load(paste0(paste0(workdir,'/../pages_proxies.RData', sep='')))
  load("/scratch3/veronika/reuse/pages_proxies.RData")

  if(archivetype == "tree") {
    if(exists("mrNH")){rm(mrNH)}
    if(exists("mrSH")){rm(mrSH)}
    mylist.names = c("data","time","lon","lat","archivetype","elevation")
    p_tree = setNames(vector("list", length(mylist.names)), mylist.names)
    # keep the data set by fsyr and feyr
    ti <- which((pages_proxies$year >= fsyr) & (pages_proxies$year <= feyr))
    p_tree$data = scale(pages_proxies$chronologies[ti,which(pages_proxies$archivetype =="tree")])
    p_tree$time = pages_proxies$year[ti]
    p_tree$lon = pages_proxies$lonlat[1, which(pages_proxies$archivetype =="tree")]
    p_tree$lat = pages_proxies$lonlat[2, which(pages_proxies$archivetype =="tree")]
    p_tree$archivetype = pages_proxies$archivetype[which(pages_proxies$archivetype =="tree")]
    p_tree$elevation = pages_proxies$elevation[which(pages_proxies$archivetype =="tree")]
    # Create a new variable for 1901-1970
    # the period on which we will calculate the multiple lin regression -> cut out 1901-1970
    start_yr = which(p_tree$time == "1901")
    end_yr = which(p_tree$time == "1970")
    mylist.names = c("data","time")
    p_tree_1901_1970 = setNames(vector("list", length(mylist.names)), mylist.names)
    # p_tree_1901_1970$data = scale(pages_proxies$chronologies[,which(pages_proxies$archivetype =="tree")])
    p_tree_1901_1970$data =  p_tree$data[start_yr:end_yr,]
    p_tree_1901_1970$time =p_tree$time[start_yr:end_yr]
    if (validate == "CRU") {
      nc=nc_open(paste(crupath,'/cru_allvar_abs_1901-2004.nc',sep=''), write=F)
      t1=ncvar_get(nc, "temp2") # for CRU temp
      t2 <- t1[,,1:852] # from year 1901 till 1971
    } else if (validate == "GISS") {
      nc=nc_open(paste(gisspath,'/air.2x2.250_anom_echamgrid.nc',sep=''), write=F) 
      t1=ncvar_get(nc, "air") # for GISS air temp and SST
      t2 <- t1[,,253:1104]
    }
    lonlist=ncvar_get(nc, "lon") 
    lonlist[lonlist > 180] <- lonlist[lonlist > 180] - 360
    latlist=ncvar_get(nc, "lat") 
    nc_close(nc)
    
    # start calculation the regression for each tree data  

    for (i in 1:length(p_tree$lon)) {
      k=which(abs(lonlist-p_tree$lon[i]+0.001)==min(abs(lonlist-p_tree$lon[i]+0.001)))
      l=which(abs(latlist-p_tree$lat[i])==min(abs(latlist-p_tree$lat[i])))
      t3 <- t2[k,l,]
      if (all(is.na(t3))) { # There is no observation data in the CRU or GISS -> the regression cannot be calculated
        if (p_tree$lat[i] > 0){
          if (!exists("mrNH")) {
            mrNH <- rep(NA,(length(grep("t.",regression_months,fixed=TRUE))+1)) #nr of coeff + intercept
            var_residuNH <- NA
          } else {
            mrNH <- rbind(mrNH,rep(NA,(length(grep("t.",regression_months,fixed=TRUE))+1)))#nr of coeff + intercept
            var_residuNH <- c(var_residuNH,rep(NA,1)) }
        }else if (!exists("mrSH")){
            mrSH <- rep(NA,(length(grep("t.",regression_months,fixed=TRUE))+1)) # nr of coeff + intercept
            var_residuSH <- NA
        }else{
            mrSH <- rbind(mrSH,rep(NA,(length(grep("t.",regression_months,fixed=TRUE))+1)))
            var_residuSH <- c(var_residuSH,rep(NA,1)) 
          }
      }else{ # we can calculate the regression
        # make variable for each month
        t5 = t3[c(4:(length(t3)-3))]
        t4 = t(array(t5,c(6,length(t5)/6))) #makes half year bunches
        tSH = t4[seq(2,nrow(t4),2),] #takes the months from Oct to March
        tNH = t4[seq(1,nrow(t4),2),][-c(1),] #takes the months from April to September
        # t4 = t(array(t5,c(12,length(t5)/12))) # 70 years from Apr-March, cuts the last half year without warning
        
        if (p_tree$lat[i] > 0) { # which half a year we want to use for calculating the regression
          unab <- tNH
        } else {
          unab <- tSH
        }                  
          colnames(unab) <- c('t.first','t.second','t.third','t.fourth','t.fifth','t.sixth')
          onlytemp<-regression_months[grep("t.",regression_months,fixed=TRUE)]
          unab <- unab[,match(onlytemp,colnames(unab))]
        
        # multiple linear regression
        results <- lm(p_tree_1901_1970$data[,i]~unab,  na.action=na.exclude)
        corr = cor.test(fitted.values(results),p_tree_1901_1970$data[,i]) 
        # print(corr[4]) # maybe under a certain corr value we could just set it to NA?
        if (p_tree$lat[i] > 0){ 
        if (!exists("mrNH")) { 
          mrNH <- results$coefficients
          var_residuNH <- var(results$residuals)
        } else { 
          mrNH <- rbind(mrNH,results$coefficients)
          var_residuNH <- c(var_residuNH,var(results$residuals))
        }
        }else{
          if (!exists("mrSH")) { 
            mrSH <- results$coefficients
            var_residuSH <- var(results$residuals)
          } else { 
            mrSH <- rbind(mrSH,results$coefficients)
            var_residuSH <- c(var_residuSH,var(results$residuals))
          }
        }
      }
    }
    
    mrtmp<-matrix(NA, nrow(mrNH)+nrow(mrSH),13)
    colnames(mrtmp) <- c("(Intercept)","unabt.first","unabt.second","unabt.third","unabt.fourth","unabt.fifth","unabt.sixth","unabp.first","unabp.second","unabp.third","unabp.fourth","unabp.fifth","unabp.sixth")
    mrtmp[1:nrow(mrNH),match(colnames(mrNH),colnames(mrtmp))]<-mrNH 
    mrtmp[(nrow(mrNH)+1):nrow(mrtmp),match(colnames(mrSH),colnames(mrtmp))]<-mrSH
    
    p_tree$mr <- mrtmp 
    p_tree$var_residu <- c(var_residuNH,var_residuSH)
    invisible(p_tree)
    
    
  } else if (archivetype == "coral") {
    if(exists("mrNH")){rm(mrNH)}
    if(exists("mrSH")){rm(mrSH)}
    mylist.names = c("data","time","lon","lat","archivetype","elevation")
    p_coral = setNames(vector("list", length(mylist.names)), mylist.names)
    # keep the data set by fsyr and feyr
    ti <- which((pages_proxies$year>= fsyr) & (pages_proxies$year <= feyr))
    p_coral$data = scale(pages_proxies$chronologies[ti,which(pages_proxies$archivetype =="coral")])
    p_coral$time = pages_proxies$year[ti]
    p_coral$lon = pages_proxies$lonlat[1, which(pages_proxies$archivetype =="coral")]
    p_coral$lat = pages_proxies$lonlat[2, which(pages_proxies$archivetype =="coral")]
    p_coral$archivetype = pages_proxies$archivetype[which(pages_proxies$archivetype =="coral")]
    p_coral$elevation = pages_proxies$elevation[which(pages_proxies$archivetype =="coral")]
    # Create a new variable for 1901-1970
    # the period on which we will calculate the multiple lin regression -> cut out 1901-1970
    start_yr = which(p_coral$time == "1901")
    end_yr = which(p_coral$time == "1970")
    mylist.names = c("data","time")
    p_coral_1901_1970 = setNames(vector("list", length(mylist.names)), mylist.names)
    # p_coral_1901_1970$data = scale(pages_proxies$chronologies[,which(pages_proxies$archivetype =="coral")])
    p_coral_1901_1970$data =  p_coral$data[start_yr:end_yr,]
    p_coral_1901_1970$time = p_coral$time[start_yr:end_yr]
    if (validate == "CRU") {
      nc=nc_open(paste(crupath,'/cru_allvar_abs_1901-2004.nc',sep=''), write=F)
      t1=ncvar_get(nc, "temp2") # for CRU temp
      t2 <- t1[,,1:852] # from year 1901 till 1971
    } else if (validate == "GISS") {
      nc=nc_open(paste(gisspath,'/air.2x2.250_anom_echamgrid.nc',sep=''), write=F) 
      t1=ncvar_get(nc, "air") # for GISS air temp and SST
      t2 <- t1[,,253:1104]
    }
    lonlist=ncvar_get(nc, "lon") 
    lonlist[lonlist > 180] <- lonlist[lonlist > 180] - 360
    latlist=ncvar_get(nc, "lat") 
    nc_close(nc)
    
    # start calculation the regression for each tree data  
    for (i in 1:length(p_coral$lon)) {
      if (all(is.na(p_coral_1901_1970$data[,i]))) {
        if (p_coral$lat[i] > 0){
          if (!exists("mrNH")) {
            mrNH <- rep(NA,(length(grep("t.",regression_months,fixed=TRUE))+1)) #nr of coeff + intercept
            var_residuNH <- NA
          } else {
            mrNH <- rbind(mrNH,rep(NA,(length(grep("t.",regression_months,fixed=TRUE))+1)))#nr of coeff + intercept
            var_residuNH <- c(var_residuNH,rep(NA,1)) }
        }else if (!exists("mrSH")){
          mrSH <- rep(NA,(length(grep("t.",regression_months,fixed=TRUE))+1)) # nr of coeff + intercept
          var_residuSH <- NA
        }else{
          mrSH <- rbind(mrSH,rep(NA,(length(grep("t.",regression_months,fixed=TRUE))+1)))
          var_residuSH <- c(var_residuSH,rep(NA,1)) 
        }
      } else {
        k=which(abs(lonlist-p_coral$lon[i]+0.001)==min(abs(lonlist-p_coral$lon[i]+0.001)))
        l=which(abs(latlist-p_coral$lat[i])==min(abs(latlist-p_coral$lat[i])))
        t3 <- t2[k,l,]
        if (all(is.na(t3))) { # There is no observation data in the CRU or GISS -> the regression cannot be calculated
          if (p_coral$lat[i] > 0){
            if (!exists("mrNH")) {
              mrNH <- rep(NA,(length(grep("t.",regression_months,fixed=TRUE))+1)) #nr of coeff + intercept
              var_residuNH <- NA
            } else {
              mrNH <- rbind(mrNH,rep(NA,(length(grep("t.",regression_months,fixed=TRUE))+1)))#nr of coeff + intercept
              var_residuNH <- c(var_residuNH,rep(NA,1)) }
          }else if (!exists("mrSH")){
            mrSH <- rep(NA,(length(grep("t.",regression_months,fixed=TRUE))+1)) # nr of coeff + intercept
            var_residuSH <- NA
          }else{
            mrSH <- rbind(mrSH,rep(NA,(length(grep("t.",regression_months,fixed=TRUE))+1)))
            var_residuSH <- c(var_residuSH,rep(NA,1)) 
          }
        } else { # we can calculate the regression
          # make variable for each month
          t5 = t3[c(4:(length(t3)-3))]
          t4 = t(array(t5,c(6,length(t5)/6))) #makes half year bunches
          tSH = t4[seq(2,nrow(t4),2),] #takes the months from Oct to March
          tNH = t4[seq(1,nrow(t4),2),][-c(1),] #takes the months from April to September
          # t4 = t(array(t5,c(12,length(t5)/12))) # 70 years from Apr-March
          
          if (p_coral$lat[i] > 0) { # which half a year we want to use for calculating the regression
            unab <- tNH
          } else {
            unab <- tSH
          }                  
          colnames(unab) <- c('t.first','t.second','t.third','t.fourth','t.fifth','t.sixth')
          onlytemp<-regression_months[grep("t.",regression_months,fixed=TRUE)]
          unab <- unab[,match(onlytemp,colnames(unab))]
          
          # multiple linear regression
          results <- lm(p_coral_1901_1970$data[,i]~unab,  na.action=na.exclude)
          corr = cor.test(fitted.values(results),p_coral_1901_1970$data[,i]) 
          # print(corr[4]) # correlations are very small
          if (p_coral$lat[i] > 0){ 
            if (!exists("mrNH")) { 
              mrNH <- results$coefficients
              var_residuNH <- var(results$residuals)
            } else { 
              mrNH <- rbind(mrNH,results$coefficients)
              var_residuNH <- c(var_residuNH,var(results$residuals))
            }
          }else{
            if (!exists("mrSH")) { 
              mrSH <- results$coefficients
              var_residuSH <- var(results$residuals)
            } else { 
              mrSH <- rbind(mrSH,results$coefficients)
              var_residuSH <- c(var_residuSH,var(results$residuals))
            }
          }
        }
      }
    }
    
    mrtmp<-matrix(NA, nrow(mrNH)+nrow(mrSH),13)
    colnames(mrtmp) <- c("(Intercept)","unabt.first","unabt.second","unabt.third","unabt.fourth","unabt.fifth","unabt.sixth","unabp.first","unabp.second","unabp.third","unabp.fourth","unabp.fifth","unabp.sixth")
    mrtmp[1:nrow(mrNH),match(colnames(mrNH),colnames(mrtmp))]<-mrNH 
    mrtmp[(nrow(mrNH)+1):nrow(mrtmp),match(colnames(mrSH),colnames(mrtmp))]<-mrSH
    
    p_coral$mr <- mrtmp
    p_coral$var_residu <- c(var_residuNH,var_residuSH)
    invisible(p_coral)
  } else if (archivetype == "documents") {
    mylist.names = c("data","time","lon","lat","archivetype","elevation")
    p_docu = setNames(vector("list", length(mylist.names)), mylist.names)
    # keep the data set by fsyr and feyr
    ti <- which((pages_proxies$year >= fsyr) & (pages_proxies$year <= feyr))
    p_docu$data = pages_proxies$chronologies[ti,which(pages_proxies$archivetype =="documents")]
    p_docu$time = pages_proxies$year[ti]
    p_docu$lon = pages_proxies$lonlat[1, which(pages_proxies$archivetype =="documents")]
    p_docu$lat = pages_proxies$lonlat[2, which(pages_proxies$archivetype =="documents")]
    p_docu$archivetype = pages_proxies$archivetype[which(pages_proxies$archivetype =="documents")]
    p_docu$elevation = pages_proxies$elevation[which(pages_proxies$archivetype =="documents")]
    invisible(p_docu)
  } else if (archivetype == "instrumental") {
    mylist.names = c("data","time","lon","lat","archivetype","elevation")
    p_inst = setNames(vector("list", length(mylist.names)), mylist.names)
    # keep the data set by fsyr and feyr
    ti <- which((pages_proxies$year>= fsyr) & (pages_proxies$year <= feyr))
    p_inst$data = pages_proxies$chronologies[ti,which(pages_proxies$archivetype =="instrumental")]
    p_inst$time = pages_proxies$year[ti]
    p_inst$lon = pages_proxies$lonlat[1, which(pages_proxies$archivetype =="instrumental")]
    p_inst$lat = pages_proxies$lonlat[2, which(pages_proxies$archivetype =="instrumental")]
    p_inst$archivetype = pages_proxies$archivetype[which(pages_proxies$archivetype =="instrumental")]
    p_inst$elevation = pages_proxies$elevation[which(pages_proxies$archivetype =="instrumental")]
    invisible(p_inst)
  }
}


# N-TREND proxies are TEMPERATURE sensitive 
# NH: from May till August as in the paper of Anchukaitis et al, 2017
# Maybe we could use the whole half year
read_ntrend = function(fsyr,feyr, validate) {
  load(paste0(paste0(workdir,'/../ntrend_proxies.RData', sep='')))
  mylist.names = c("data","time","lon","lat","archivetype","elevation","parameter")
  ntrend = setNames(vector("list", length(mylist.names)), mylist.names)
  ti <- which((ntrend_proxies$year >= fsyr) & (ntrend_proxies$year <= feyr))
  ntrend$data =  scale(ntrend_proxies$chronologies[ti,])
  ntrend$time = ntrend_proxies$year[ti]
  ntrend$lon =  ntrend_proxies$lonlat[1,]
  ntrend$lat =  ntrend_proxies$lonlat[2,]
  ntrend$archivetype =  rep("tree", length(ntrend_proxies$parameter))
  ntrend$elevation =  ntrend_proxies$elevation
  ntrend$parameter = ntrend_proxies$parameter
  # Create a new variable for 1901-1970
  # the period on which we will calculate the multiple lin regression -> cut out 1901-1970
  start_yr = which(ntrend$time == "1901")
  end_yr = which(ntrend$time == "1970")
  mylist.names = c("data","time")
  ntrend_1901_1970 = setNames(vector("list", length(mylist.names)), mylist.names)
  ntrend_1901_1970$data = ntrend$data[start_yr:end_yr,]
  ntrend_1901_1970$time = ntrend$time[start_yr:end_yr]
  if (validate == "CRU") {
    nc=nc_open(paste(crupath,'/cru_allvar_abs_1901-2004.nc',sep=''), write=F)
    t1=ncvar_get(nc, "temp2") # for CRU temp
    t2 <- t1[,,1:840] # from year 1901 till 1970
  } else if (validate == "GISS") {
    nc=nc_open(paste(gisspath,'/air.2x2.250_anom_echamgrid.nc',sep=''), write=F) 
    t1=ncvar_get(nc, "air") # for GISS air temp and SST
    t2 <- t1[,,253:1092]
  }
  lonlist=ncvar_get(nc, "lon") 
  lonlist[lonlist > 180] <- lonlist[lonlist > 180] - 360
  latlist=ncvar_get(nc, "lat") 
  nc_close(nc)
  
  # start calculation the regression for each tree data  
  for (i in 1:length(ntrend$lon)) {
    k=which(abs(lonlist-ntrend$lon[i]+0.001)==min(abs(lonlist-ntrend$lon[i]+0.001)))
    l=which(abs(latlist-ntrend$lat[i])==min(abs(latlist-ntrend$lat[i])))
    t3 <- t2[k,l,]
    if (all(is.na(t3))) { # There is no observation data in the CRU or GISS -> the regression cannot be calculated
      if (i==1) {
        mr <- rep(NA,(length(grep("t.",regression_months,fixed=TRUE))+1)) # coeff +intercept
        var_residu <- NA
      } else {
        mr <- rbind(mr,rep(NA,(length(grep("t.",regression_months,fixed=TRUE))+1)))
        var_residu <- c(var_residu,rep(NA,1)) 
      }
    } else { # we can calculate the regression
      # make variable for each month
      t4 <- t(array(t3,c(12,length(t3)/12))) # 70 years from Jan till December
      unab <- t4[,4:9] # 70 years from April to September
      
      colnames(unab) <- c("t.first","t.second","t.third","t.fourth","t.fifth","t.sixth")
      onlytemp<-regression_months[grep("t.",regression_months,fixed=TRUE)]
      unab <- unab[,match(onlytemp,colnames(unab))]
      # multiple linear regression
      results <- lm(ntrend_1901_1970$data[,i]~unab,  na.action=na.exclude)
      corr = cor.test(fitted.values(results),ntrend_1901_1970$data[,i]) 
      # print(corr[4]) # maybe under a certain corr value we could just set it to NA?
      if (i==1) { 
        mr <- results$coefficients
        var_residu <- var(results$residuals)
      } else { 
        mr <- rbind(mr,results$coefficients)
        var_residu <- c(var_residu,var(results$residuals))
      }
    }
  }
  
  mrtmp<-matrix(NA, nrow(mr),13)
  colnames(mrtmp) <-c("(Intercept)","unabt.first","unabt.second","unabt.third","unabt.fourth","unabt.fifth","unabt.sixth","unabp.first","unabp.second","unabp.third","unabp.fourth","unabp.fifth","unabp.sixth")
  mrtmp[,match(colnames(mr),colnames(mrtmp))]<-mr
  
  
  
  ntrend$mr <- mrtmp
  ntrend$var_residu <- var_residu
  invisible(ntrend)
}








compute_Hi_Hredux_proxy <- function(stations, echam, regcoef=NULL, threshold=700){
  H <- array(0, c(nrow(stations$data),(ncol(stations$mr)-1)*2))
  nech <- length(which(echam$names=='temp2'))/6
  # ACHTUNG: works only if all echam data are equally sized fields, i.e. no_stream=T
  for (i in seq(stations$lon)){
    if (!is.na(stations$lon[i])) { 
      dist <- compute_dist(echam$lon, echam$lat, stations$lon[i], stations$lat[i])
      # state vector contains t_apr, t_may, ..., t_sep, p_apr, ..., s_sep, ...
      # multiple regession coeff. must be in right month of t and p
      # H[i, which.min(dist)] <- if (min(dist) < threshold) regcoef[i,1] else 0 # NO coeff for t_apr
      # May T
      if (min(dist) < threshold) {
        #Temperature
        H[i,1] <- if (!is.na(regcoef[i,2])) which.min(dist) else NA 
        H[i,2] <- if (!is.na(regcoef[i,2])) regcoef[i,2]  else NA
        H[i,3] <- if (!is.na(regcoef[i,3])) which.min(dist)+((dim(echam$data)[1]/6)) else NA
        H[i,4] <- if (!is.na(regcoef[i,3])) regcoef[i,3] else NA
        H[i,5] <- if (!is.na(regcoef[i,4])) which.min(dist)+(2*(dim(echam$data)[1]/6)) else NA
        H[i,6] <- if (!is.na(regcoef[i,4])) regcoef[i,4] else NA
        H[i,7] <- if (!is.na(regcoef[i,5])) which.min(dist)+(3*(dim(echam$data)[1]/6)) else NA
        H[i,8] <- if (!is.na(regcoef[i,5])) regcoef[i,5] else NA
        H[i,9] <- if (!is.na(regcoef[i,6])) which.min(dist)+(4*(dim(echam$data)[1]/6)) else NA
        H[i,10] <- if (!is.na(regcoef[i,6])) regcoef[i,6] else NA
        H[i,11] <- if (!is.na(regcoef[i,7])) which.min(dist)+(5*(dim(echam$data)[1]/6))+nech else NA
        H[i,12] <- if (!is.na(regcoef[i,7])) regcoef[i,7] else NA
        #Precipitation
        H[i,13] <- if (!is.na(regcoef[i,8])) which.min(dist)+nech else NA
        H[i,14] <- if (!is.na(regcoef[i,8])) regcoef[i,8] else NA
        H[i,15] <- if (!is.na(regcoef[i,9])) which.min(dist)+((dim(echam$data)[1]/6))+nech else NA
        H[i,16] <- if (!is.na(regcoef[i,9])) regcoef[i,9] else NA
        H[i,17] <- if (!is.na(regcoef[i,10])) which.min(dist)+(2*(dim(echam$data)[1]/6))+nech else NA
        H[i,18] <- if (!is.na(regcoef[i,10])) regcoef[i,10] else NA
        H[i,19] <- if (!is.na(regcoef[i,11])) which.min(dist)+(3*(dim(echam$data)[1]/6))+nech else NA
        H[i,20] <- if (!is.na(regcoef[i,11])) regcoef[i,11] else NA
        H[i,21] <- if (!is.na(regcoef[i,12])) which.min(dist)+(4*(dim(echam$data)[1]/6))+nech else NA
        H[i,22] <- if (!is.na(regcoef[i,12])) regcoef[i,12] else NA
        H[i,23] <- if (!is.na(regcoef[i,13])) which.min(dist)+(5*(dim(echam$data)[1]/6))+nech else NA
        H[i,24] <- if (!is.na(regcoef[i,13])) regcoef[i,13] else NA
       

      }
    }   
  }
  H[which(is.na(H))] <- 0
  H <- H[,abs(colSums(H)) > 0]
  return(H)
}







 compute_H_proxy <- function(stations, echam, regcoef=NULL, threshold=700){
   H <- array(0, c(nrow(stations$data), nrow(echam$data)))
   for (i in seq(stations$lon)){
     if (!is.na(stations$lon[i])) { 
       dist <- compute_dist(echam$lon, echam$lat, stations$lon[i], stations$lat[i])
       # state vector contains t_apr, t_may, ..., t_sep, p_apr, ..., s_sep, ...
       # multiple regession coeff. must be in right month of t and p
       # H[i, which.min(dist)] <- if (min(dist) < threshold) regcoef[i,1] else 0 # NO coeff for t_apr
       H[i, which.min(dist)+(dim(echam$data)[1]/6)] <- if (min(dist) < threshold)
         regcoef[i,2] else 0
       H[i, which.min(dist)+(2*(dim(echam$data)[1]/6))] <- if (min(dist) < threshold)
         regcoef[i,3] else 0
       H[i, which.min(dist)+(3*(dim(echam$data)[1]/6))] <- if (min(dist) < threshold)
         regcoef[i,4] else 0
       H[i, which.min(dist)+(4*(dim(echam$data)[1]/6))] <- if (min(dist) < threshold)
         regcoef[i,5] else 0
       H[i, which.min(dist)+(one_var_dim)] <- if (min(dist) < threshold) 
         regcoef[i,6] else 0     
       H[i, which.min(dist)+(dim(echam$data)[1]/6)+one_var_dim] <- if 
         (min(dist) < threshold) regcoef[i,7] else 0     
       H[i, which.min(dist)+(2*(dim(echam$data)[1]/6))+one_var_dim] <- if 
         (min(dist) < threshold) regcoef[i,8] else 0     
     }
   }   
   H[which(is.na(H))] <- 0
   return(H)
 }






compute_Hi_Hredux_sixmonstatevector <- function(stations, echam, threshold=700){
  H <- array(0, c(nrow(stations$data), 2))
  nstat <- length(stations$names)/6
  nech <- length(which(echam$names=='temp2'))/6
  # ACHTUNG: works only if all echam data are equally sized fields, i.e. no_stream=T
  for (i in seq(stations$lon)){
    if ((!is.na(stations$lon[i])) & (!is.na(stations$lat[i]))) {
      dist <- compute_dist(echam$lon, echam$lat, stations$lon[i], stations$lat[i])
      # state vector contains 6 mon oct to mar or apr to sep
      if (stations$names[i] == 'temp2') {
        if (i < nstat+1) {
          if (min(dist) < threshold) {
            H[i,1] <- which.min(dist)
            H[i,2] <- 1
          } 
        }
        if ((i > nstat) & (i < (2*nstat+1))) {
          if (min(dist) < threshold) {
            H[i,1] <- which.min(dist)+(dim(echam$data)[1]/6)
            H[i,2] <- 1
          } 
        }
        if ((i > (2*nstat)) & (i < (3*nstat+1))) {
          if (min(dist) < threshold) {
            H[i,1] <- which.min(dist)+(2*(dim(echam$data)[1]/6))
            H[i,2] <- 1
          } 
        }
        if ((i > (3*nstat)) & (i < (4*nstat+1))) {
          if (min(dist) < threshold) {
            H[i,1] <- which.min(dist)+(3*(dim(echam$data)[1]/6))
            H[i,2] <- 1
          } 
        }
        if ((i > (4*nstat)) & (i < (5*nstat+1))) {
          if (min(dist) < threshold) {
            H[i,1] <- which.min(dist)+(4*(dim(echam$data)[1]/6))
            H[i,2] <- 1
          } 
        }
        if ((i > (5*nstat)) & (i < (6*nstat+1))) {
          if (min(dist) < threshold) {
            H[i,1] <- which.min(dist)+(5*(dim(echam$data)[1]/6))
            H[i,2] <- 1
          } 
        }
      }
      if (stations$names[i] == 'precip') {
        if (i < nstat+1) {
          if (min(dist) < threshold) {
            H[i,1] <- which.min(dist)+(0*(dim(echam$data)[1]/6)+nech)
            H[i,2] <- 1
          } 
        }
        if ((i > nstat) & (i < (2*nstat+1))) {
          if (min(dist) < threshold) {
            H[i,1] <- which.min(dist)+(1*(dim(echam$data)[1]/6)+nech)
            H[i,2] <- 1
          } 
        }
        if ((i > (2*nstat)) & (i < (3*nstat+1))) {
          if (min(dist) < threshold) {
            H[i,1] <- which.min(dist)+(2*(dim(echam$data)[1]/6)+nech)
            H[i,2] <- 1
          } 
        }
        if ((i > (3*nstat)) & (i < (4*nstat+1))) {
          if (min(dist) < threshold) {
            H[i,1] <- which.min(dist)+(3*(dim(echam$data)[1]/6)+nech)
            H[i,2] <- 1
          } 
        }
        if ((i > (4*nstat)) & (i < (5*nstat+1))) {
          if (min(dist) < threshold) {
            H[i,1]  <- which.min(dist)+(4*(dim(echam$data)[1]/6)+nech)
            H[i,2] <- 1
          } 
        }
        if ((i > (5*nstat)) & (i < (6*nstat+1))) {
          if (min(dist) < threshold) {
            H[i,1] <- which.min(dist)+(5*(dim(echam$data)[1]/6)+nech)
            H[i,2] <- 1
          } 
        }
      }
      if (stations$names[i] == 'slp') {
        if (i < nstat+1) {
          if (min(dist) < threshold) {
            H[i,1] <- which.min(dist)+(0*(dim(echam$data)[1]/6)+(2*nech))
            H[i,2] <- 1
          } 
        }
        if ((i > nstat) & (i < (2*nstat+1))) {
          if (min(dist) < threshold) {
            H[i,1] <- which.min(dist)+(1*(dim(echam$data)[1]/6)+(2*nech))
            H[i,2] <- 1
          } 
        }
        if ((i > (2*nstat)) & (i < (3*nstat+1))) {
          if (min(dist) < threshold) {
            H[i,1] <- which.min(dist)+(2*(dim(echam$data)[1]/6)+(2*nech))
            H[i,2] <- 1
          } 
        }
        if ((i > (3*nstat)) & (i < (4*nstat+1))) {
          if (min(dist) < threshold) {
            H[i,1] <- which.min(dist)+(3*(dim(echam$data)[1]/6)+(2*nech))
            H[i,2] <- 1
          } 
        }
        if ((i > (4*nstat)) & (i < (5*nstat+1))) {
          if (min(dist) < threshold) {
            H[i,1] <- which.min(dist)+(4*(dim(echam$data)[1]/6)+(2*nech))
            H[i,2] <- 1
          } 
        }
        if ((i > (5*nstat)) & (i < (6*nstat+1))) {
          if (min(dist) < threshold) {
            H[i,1] <- which.min(dist)+(5*(dim(echam$data)[1]/6)+(2*nech))
            H[i,2] <- 1
          } 
        }
      }
    }
  }
  return(H)
}











  

compute_H_sixmonstatevector <- function(stations, echam, threshold=700){
  H <- array(0, c(nrow(stations$data), nrow(echam$data)))
  for (i in seq(stations$lon)){
   if ((!is.na(stations$lon[i])) & (!is.na(stations$lat[i]))) {
    dist <- compute_dist(echam$lon, echam$lat, stations$lon[i], stations$lat[i])
    # state vector contains 6 mon oct to mar or apr to sep
    nstat <- length(stations$names)/6
    if (stations$names[i] == 'temp2') {
      if (i < nstat+1) {
        H[i, which.min(dist)] <- if (min(dist) < threshold) 1 else 0
      }
      if ((i > nstat) & (i < (2*nstat+1))) {
        H[i, which.min(dist)+(dim(echam$data)[1]/6)] <- if (min(dist) < threshold) 
          1 else 0
      }
      if ((i > (2*nstat)) & (i < (3*nstat+1))) {
        H[i, which.min(dist)+(2*(dim(echam$data)[1]/6))] <- if (min(dist) < threshold)
          1 else 0
      }
      if ((i > (3*nstat)) & (i < (4*nstat+1))) {
        H[i, which.min(dist)+(3*(dim(echam$data)[1]/6))] <- if (min(dist) < threshold)
          1 else 0
      }
      if ((i > (4*nstat)) & (i < (5*nstat+1))) {
        H[i, which.min(dist)+(4*(dim(echam$data)[1]/6))] <- if (min(dist) < threshold)
          1 else 0
      }
      if ((i > (5*nstat)) & (i < (6*nstat+1))) {
        H[i, which.min(dist)+(5*(dim(echam$data)[1]/6))] <- if (min(dist) < threshold)
          1 else 0
      }
    }
    if (stations$names[i] == 'precip') {
      if (i < nstat+1) {
        H[i, which.min(dist)+(0*(dim(echam$data)[1]/6)+one_var_dim)] <- if 
          (min(dist) < threshold) 1 else 0
      }
      if ((i > nstat) & (i < (2*nstat+1))) {
        H[i, which.min(dist)+(1*(dim(echam$data)[1]/6)+one_var_dim)] <- if 
          (min(dist) < threshold) 1 else 0
      }
      if ((i > (2*nstat)) & (i < (3*nstat+1))) {
        H[i, which.min(dist)+(2*(dim(echam$data)[1]/6)+one_var_dim)] <- if 
          (min(dist) < threshold) 1 else 0
      }
      if ((i > (3*nstat)) & (i < (4*nstat+1))) {
        H[i, which.min(dist)+(3*(dim(echam$data)[1]/6)+one_var_dim)] <- if 
          (min(dist) < threshold) 1 else 0
      }
      if ((i > (4*nstat)) & (i < (5*nstat+1))) {
        H[i, which.min(dist)+(4*(dim(echam$data)[1]/6)+one_var_dim)] <- if 
          (min(dist) < threshold) 1 else 0
      }
      if ((i > (5*nstat)) & (i < (6*nstat+1))) {
        H[i, which.min(dist)+(5*(dim(echam$data)[1]/6)+one_var_dim)] <- if 
          (min(dist) < threshold) 1 else 0
      }
    }
    if (stations$names[i] == 'slp') {
      if (i < nstat+1) {
        H[i, which.min(dist)+(0*(dim(echam$data)[1]/6)+(2*one_var_dim))] <- if 
        (min(dist) < threshold) 1 else 0
      }
      if ((i > nstat) & (i < (2*nstat+1))) {
        H[i, which.min(dist)+(1*(dim(echam$data)[1]/6)+(2*one_var_dim))] <- if 
        (min(dist) < threshold) 1 else 0
      }
      if ((i > (2*nstat)) & (i < (3*nstat+1))) {
        H[i, which.min(dist)+(2*(dim(echam$data)[1]/6)+(2*one_var_dim))] <- if 
        (min(dist) < threshold) 1 else 0
      }
      if ((i > (3*nstat)) & (i < (4*nstat+1))) {
        H[i, which.min(dist)+(3*(dim(echam$data)[1]/6)+(2*one_var_dim))] <- if 
        (min(dist) < threshold) 1 else 0
      }
      if ((i > (4*nstat)) & (i < (5*nstat+1))) {
        H[i, which.min(dist)+(4*(dim(echam$data)[1]/6)+(2*one_var_dim))] <- if 
        (min(dist) < threshold) 1 else 0
      }
      if ((i > (5*nstat)) & (i < (6*nstat+1))) {
        H[i, which.min(dist)+(5*(dim(echam$data)[1]/6)+(2*one_var_dim))] <- if 
        (min(dist) < threshold) 1 else 0
      }
    }
   }
  }
  return(H)
}

  



compute_H_sixmonstatevector_docu <- function(stations, echam, threshold=700){
  H <- array(0, c(nrow(stations$data), nrow(echam$data)))
  for (i in seq(stations$lon)){
#    print(i)
    dist <- compute_dist(echam$lon, echam$lat, stations$lon[i], stations$lat[i])
    # state vector contains 6 mon oct to mar or apr to sep
    if (stations$des[i]=='mon') {
#      print("mon")
      nstat <- length(stations$des[stations$des=='mon'])/6
      if (stations$names[i] == 'temp2') {
        if (i < nstat+1) {
          H[i, which.min(dist)] <- if (min(dist) < threshold) 1 else 0
        }
        if ((i > nstat) & (i < (2*nstat+1))) {
          H[i, which.min(dist)+(dim(echam$data)[1]/6)] <- if (min(dist) < threshold) 
            1 else 0
        }
        if ((i > (2*nstat)) & (i < (3*nstat+1))) {
          H[i, which.min(dist)+(2*(dim(echam$data)[1]/6))] <- if (min(dist) < threshold)
            1 else 0
        }
        if ((i > (3*nstat)) & (i < (4*nstat+1))) {
          H[i, which.min(dist)+(3*(dim(echam$data)[1]/6))] <- if (min(dist) < threshold)
            1 else 0
        }
        if ((i > (4*nstat)) & (i < (5*nstat+1))) {
          H[i, which.min(dist)+(4*(dim(echam$data)[1]/6))] <- if (min(dist) < threshold)
            1 else 0
        }
        if ((i > (5*nstat)) & (i < (6*nstat+1))) {
          H[i, which.min(dist)+(5*(dim(echam$data)[1]/6))] <- if (min(dist) < threshold)
            1 else 0
        }
      }
      if (stations$names[i] == 'precip') {
        if (i < nstat+1) {
          H[i, which.min(dist)+(0*(dim(echam$data)[1]/6)+one_var_dim)] <- if 
          (min(dist) < threshold) 1 else 0
        }
        if ((i > nstat) & (i < (2*nstat+1))) {
          H[i, which.min(dist)+(1*(dim(echam$data)[1]/6)+one_var_dim)] <- if 
          (min(dist) < threshold) 1 else 0
        }
        if ((i > (2*nstat)) & (i < (3*nstat+1))) {
          H[i, which.min(dist)+(2*(dim(echam$data)[1]/6)+one_var_dim)] <- if 
          (min(dist) < threshold) 1 else 0
        }
        if ((i > (3*nstat)) & (i < (4*nstat+1))) {
          H[i, which.min(dist)+(3*(dim(echam$data)[1]/6)+one_var_dim)] <- if 
          (min(dist) < threshold) 1 else 0
        }
        if ((i > (4*nstat)) & (i < (5*nstat+1))) {
          H[i, which.min(dist)+(4*(dim(echam$data)[1]/6)+one_var_dim)] <- if 
          (min(dist) < threshold) 1 else 0
        }
        if ((i > (5*nstat)) & (i < (6*nstat+1))) {
          H[i, which.min(dist)+(5*(dim(echam$data)[1]/6)+one_var_dim)] <- if 
          (min(dist) < threshold) 1 else 0
        }
      }
      if (stations$names[i] == 'slp') {
        if (i < nstat+1) {
          H[i, which.min(dist)+(0*(dim(echam$data)[1]/6)+(2*one_var_dim))] <- if 
          (min(dist) < threshold) 1 else 0
        }
        if ((i > nstat) & (i < (2*nstat+1))) {
          H[i, which.min(dist)+(1*(dim(echam$data)[1]/6)+(2*one_var_dim))] <- if 
          (min(dist) < threshold) 1 else 0
        }
        if ((i > (2*nstat)) & (i < (3*nstat+1))) {
          H[i, which.min(dist)+(2*(dim(echam$data)[1]/6)+(2*one_var_dim))] <- if 
          (min(dist) < threshold) 1 else 0
        }
        if ((i > (3*nstat)) & (i < (4*nstat+1))) {
          H[i, which.min(dist)+(3*(dim(echam$data)[1]/6)+(2*one_var_dim))] <- if 
          (min(dist) < threshold) 1 else 0
        }
        if ((i > (4*nstat)) & (i < (5*nstat+1))) {
          H[i, which.min(dist)+(4*(dim(echam$data)[1]/6)+(2*one_var_dim))] <- if 
          (min(dist) < threshold) 1 else 0
        }  
        if ((i > (5*nstat)) & (i < (6*nstat+1))) {
          H[i, which.min(dist)+(5*(dim(echam$data)[1]/6)+(2*one_var_dim))] <- if 
          (min(dist) < threshold) 1 else 0
        }
      }
    } else if (stations$des[i]=='seas') {
#      print("seas")
      # set DJF for winter = 1/3 and JJA for summer
       nstat <- length(stations$des[stations$des=='seas'])
       if (stations$names[i] == 'temp2') {
       # with seasons set to april to sept and oct to march, JJA/DJF are months 3 to 5
#           H[i, which.min(dist)] <- if (min(dist) < threshold) 1 else 0
#           H[i, which.min(dist)+(dim(echam$data)[1]/6)] <- if (min(dist) < threshold) 1 else 0
           H[i, which.min(dist)+(2*(dim(echam$data)[1]/6))] <- if (min(dist) < threshold) 1/3 else 0
           H[i, which.min(dist)+(3*(dim(echam$data)[1]/6))] <- if (min(dist) < threshold) 1/3 else 0
           H[i, which.min(dist)+(4*(dim(echam$data)[1]/6))] <- if (min(dist) < threshold) 1/3 else 0
#           H[i, which.min(dist)+(5*(dim(echam$data)[1]/6))] <- if (min(dist) < threshold) 1 else 0     
       }
    } else if (stations$des[i]=='JFMA') {
#      print("JFMA")
      # set JFMA for winter = 1/4: NOT possible because april in other season
      # solution: set JFM =1/3
      nstat <- length(stations$des[stations$des=='JFMA'])
      if (stations$names[i] == 'temp2') {
#        H[i, which.min(dist)] <- if (min(dist) < threshold) 1 else 0
#        H[i, which.min(dist)+(dim(echam$data)[1]/6)] <- if (min(dist) < threshold) 1 else 0
#        H[i, which.min(dist)+(2*(dim(echam$data)[1]/6))] <- if (min(dist) < threshold) 1/3 else 0
        H[i, which.min(dist)+(3*(dim(echam$data)[1]/6))] <- if (min(dist) < threshold) 1/3 else 0
        H[i, which.min(dist)+(4*(dim(echam$data)[1]/6))] <- if (min(dist) < threshold) 1/3 else 0
        H[i, which.min(dist)+(5*(dim(echam$data)[1]/6))] <- if (min(dist) < threshold) 1/3 else 0           
      }
    } else if (stations$des[i]=='AMJJA') {
#      print("AMJJA")
      # set AMJJA for summer to 1/5
      nstat <- length(stations$des[stations$des=='AMJJA'])
      if (stations$names[i] == 'temp2') {
        H[i, which.min(dist)] <- if (min(dist) < threshold) 1/5 else 0
        H[i, which.min(dist)+(dim(echam$data)[1]/6)] <- if (min(dist) < threshold) 1/5 else 0
        H[i, which.min(dist)+(2*(dim(echam$data)[1]/6))] <- if (min(dist) < threshold) 1/5 else 0
        H[i, which.min(dist)+(3*(dim(echam$data)[1]/6))] <- if (min(dist) < threshold) 1/5 else 0
        H[i, which.min(dist)+(4*(dim(echam$data)[1]/6))] <- if (min(dist) < threshold) 1/5 else 0
#        H[i, which.min(dist)+(5*(dim(echam$data)[1]/6))] <- if (min(dist) < threshold) 1 else 0           
      }
    } else {
        print("UNKNOWN season in documentary data")
    }
  } 
  return(H)
}












# read in annual proxy data old
read_proxy <- function(timlim=c(syr=1970,eyr=1999),anomlim=c(syr=1970,eyr=1999)){
  
# >500yr (1500-) screened annual temp records from spectra paper
t_prox=c('arge091','az510','buentgen_2005','buentgen_2006','buentgen_2009','buntgen_science_JJA_temp','ca529','ca630','ca631','cana175','carpathian','corona_2009','dolomites','forfjord','gisp2o18','grudd','jamtland','junipars','lapland','mill_pyrenees_MJJAS','norw010','nv513','orokonztr','recjj_yy2','salzer2005','tan_2003_recontemp','thompson_1992_quelccao18','tornetrask','tyrol_mxd','tyrol_trw','vallee_merveille','vinther_2004_scgreenland','ak020','ak021','arge005','arge013','arge070','arge073','ca082','cana036','cana110','chil006','finl021','mexi027','mexi027e','mexi027l','norw009','russ017','schweingruber_mxdabd_grid1','schweingruber_mxdabd_grid10','schweingruber_mxdabd_grid100','schweingruber_mxdabd_grid11','schweingruber_mxdabd_grid110','schweingruber_mxdabd_grid111','schweingruber_mxdabd_grid115','schweingruber_mxdabd_grid12','schweingruber_mxdabd_grid16','schweingruber_mxdabd_grid18','schweingruber_mxdabd_grid19','schweingruber_mxdabd_grid2','schweingruber_mxdabd_grid20','schweingruber_mxdabd_grid21','schweingruber_mxdabd_grid22','schweingruber_mxdabd_grid3','schweingruber_mxdabd_grid30','schweingruber_mxdabd_grid42','schweingruber_mxdabd_grid44','schweingruber_mxdabd_grid45','schweingruber_mxdabd_grid70','schweingruber_mxdabd_grid87','schweingruber_mxdabd_grid88','schweingruber_mxdabd_grid89','wa039','wa048','wa056','wa057','wa064')

## >500yr (1500-) screened annual precip records from spectra paper
#$p_prox=c('ar049','ar050','ar052','ar053','az550','buntgen_science_AMJ_precip','ca051','ca087','ca528','ca529','ca531','ca532','ca533','ca535','ca629','ca632','ca633','cana194','chil002','co067','co511','co555','co556','co580','fisher_1996_cgreenland','fl001','ga002','ga003','ga004','helama_sweden_MJ_precip','id009','id010','la001','mo037','morc001','nc008','nm560','nm572','nv060','nv510','nv514','nv515','nv516','nv518','or012','or015','or060','or061','or062','or063','pola006','qian_2003_yriver','sc004','sd017','thompson_1992_quelccao18_f','treydte_f','tx040','ut024','ut508','va021','w3crn','w42crn','ak020','ar048','arge018','arge070','arge073','az084','az086','az106','az129','az144','az520','az547','az557','ca065','ca073','ca612','ca628','cana036','cana135','cana136','cana137','co066','co076','co509l','co509x','co535','co570','co579','id006','il016','jord001','mexi001','mexi022','mexi023','mexi023e','mexi023l','ms002','nm025','nm026','nm030','nm031','nm035','nm559','nm564','nm565','nv049','nv052','nv053','nv055','nv056','nv058','nv061','nv507','or006','or009','or018','or033','or081','spai011','turk001','tx042','tx042e','tx042l','ut018','ut020','ut022','wimmer_wien_JJA_precip','wy002','wy006','wy019','wy026')
  
for(i in 1:length(t_prox)) {
# print(c(t_prox[i],i))
 file_prox=paste0(dataintdir,"proxy/all_prox/",t_prox[i],".ppd")  
 tmp10 <- read.table(file_prox,header=FALSE)
 if (tmp10[nrow(tmp10),1] > timlim[1]) { 
  lon=tmp10[1,1]
  lat=tmp10[2,1]
  code=tmp10[3,1]
  tmp11=window(ts(tmp10[4:nrow(tmp10),2], tmp10[4,1], deltat=1), start=timlim[1], end=timlim[2], deltat=1)
  #standardize UNsmoothed series over calibration period version 2
  tmp_mean<-mean(window(tmp11,anomlim[1],anomlim[2]),na.rm=F)
  tmp_sd<-sd(window(tmp11,anomlim[1],anomlim[2]),na.rm=F)
  prox<-ts(scale(tmp11, center=tmp_mean,scale=tmp_sd),index(tmp11)[1],deltat=1)
  
  if(i==1){
    prox_all=prox
    lon_all=lon
    lat_all=lat
    name_all=t_prox[i]
    } else {
    prox_all=ts.union(prox_all,prox)
    lon_all=c(lon_all,lon)
    lat_all=c(lat_all,lat)
    name_all=c(name_all,t_prox[i])
  }
 } else { print(paste(t_prox[i],"not in time period",sep=" ")) }
}
prox_all=apply(prox_all,2,rep,each=12)
colnames(prox_all)=c(1:ncol(prox_all))
prox.data=array(as.numeric(prox_all),c(dim(prox_all)[1],dim(prox_all)[2]))

real_prox <- list(data=prox.data, lon=lon_all, lat=lat_all, names=name_all, time=seq(min(timlim) + 1/24, by=1/12, length=nrow(prox.data)))
}








# Yuri's read early instr SLP data
read_slp <- function(){
year_min <- 1749
year_max <- 2012
n_stations <- 200

x.data <- array(dim=c((year_max-year_min+1)*12,n_stations))
x.lon <- array(dim=c(n_stations))
x.lat <- array(dim=c(n_stations))
x.name <- array(data="slp",dim=c(n_stations))
x.height <- array(dim=c(n_stations))
x.time <- seq(year_min + 1/24, by=1/12, length=nrow(x.data))
i_staz <- 0


### READ KUETTEL
coords <- read.table("../assimil_data/data_yuri/slp/Stations_PP_1722_coordinates.csv",skip=1)
data <- read.table("../assimil_data/data_yuri/slp/Stations_PP_1722_monthly_final.csv",skip=1)
for (i_staz in 1:dim(coords)[1]) {
  x.lon[i_staz] <- coords[i_staz,3]
  x.lat[i_staz] <- coords[i_staz,2]
  i_data <- 12
  for (i in 1:dim(data)[1]) {
    x.data[i_data,i_staz] <- data[i,i_staz+2]
    i_data <- i_data+1
  }
# cannot convert surface pressure to SLP because no height information for stations  
#   tmp <- x.data[,i_staz]
#   # convert surface pressure to slp
#   if (mean(x.data[,i_staz],na.rm=T)<1000){
#     slpnew <- 1013.25*(1-((0.0065*x.height[i_staz])/288.15))^5.255
#   }
}

### READ HISTALP
files <- read.table("../assimil_data/data_yuri/slp/list1")
for (filename in t(files)) {
 filename <- paste('../assimil_data/data_yuri/slp/',as.character(filename),sep='')
 i_staz <- i_staz+1
 header <- readLines(filename,n=3)
 x.lon[i_staz] <- as.numeric(substr(header[3],36,40))
 x.lat[i_staz] <- as.numeric(substr(header[3],42,46))
 x.height[i_staz] <- as.integer(substr(header[3],49,52))
 data <- read.table(filename,na.strings="99999",skip=4)
 i_data <- (data[1,1]-year_min)*12+1
 for (i in 1:dim(data)[1]) {
   for (j in 2:13) {
     x.data[i_data,i_staz] <- data[i,j]/10
     i_data <- i_data+1
   }
 }
}

### READ CRU STATIONS
files <- read.table("../assimil_data/data_yuri/slp/list2")
for (filename in t(files)) {
  filename <- paste('../assimil_data/data_yuri/slp/',as.character(filename),sep='')
  i_staz <- i_staz+1
  header <- readLines(filename,n=2)
  x.lon[i_staz] <- as.numeric(substr(header[2],1,6))
  x.lat[i_staz] <- as.numeric(substr(header[2],8,12))
  data <- read.table(filename,na.strings="-999",skip=2)
  i_data <- (data[1,1]-year_min)*12+1
  for (i in 1:dim(data)[1]) {
    for (j in 2:13) {
      x.data[i_data,i_staz] <- data[i,j]/10
      i_data <- i_data+1
    }
  }
}

### READ GHCN
files <- read.table("../assimil_data/data_yuri/slp/list3")
for (filename in t(files)) {
  filename <- paste('../assimil_data/data_yuri/slp/',as.character(filename),sep='')
  i_staz <- i_staz+1
  header <- readLines(filename,n=5)
  x.lon[i_staz] <- as.numeric(substr(header[3],25,31))
  x.lat[i_staz] <- as.numeric(substr(header[3],16,21))
  x.height[i_staz] <- as.integer(substr(header[3],35,38))
  data <- read.table(filename,na.strings="-999.9",skip=5)
  i_data <- (data[1,1]-year_min)*12+1
  for (i in 1:dim(data)[1]) {
    for (j in 2:13) {
      x.data[i_data,i_staz] <- data[i,j]
      i_data <- i_data+1
    }
  }
}

### READ OTHERS
files <- read.table("../assimil_data/slp/list4")
for (filename in t(files)) {
  filename <- paste('../assimil_data/slp/',as.character(filename),sep='')
  i_staz <- i_staz+1
  header <- readLines(filename,n=1)
  x.lon[i_staz] <- as.numeric(substr(header,1,7))
  x.lat[i_staz] <- as.numeric(substr(header,9,14))
  #  x.height[i_staz] <- as.integer(substr(header,16,19))
  data <- read.table(filename,na.strings="-999.9",skip=1)
  i_data <- (data[1,1]-year_min)*12+1
  for (i in 1:dim(data)[1]) {
    for (j in 2:13) {
      x.data[i_data,i_staz] <- data[i,j]
      i_data <- i_data+1
    }
  }
}
# correct non SLP and missing data
notna <- !is.na(colMeans(x.data,na.rm=T))
x.data <- x.data[,notna]
x.lon <- x.lon[notna]
x.lat <- x.lat[notna] 
x.name <- x.name[notna]
x.height <-  x.height[notna]
timesten <- colMeans(x.data,na.rm=T)>1100
x.data[,timesten] = x.data[,timesten]/10
#trueslp <- colMeans(x.data,na.rm=T)>1000
#x.data = x.data[,trueslp]
#x.lon <- x.lon[trueslp]
#x.lat <- x.lat[trueslp] 
#x.name <- x.name[trueslp]
#x.height <-  x.height[trueslp]

slp <- list(data=x.data, lon=x.lon, lat=x.lat, names=x.name, height=x.height, time=x.time)
save(slp,file="../assim_data/slp.Rdata")
write.table(cbind(x.lon,x.lat),file="../assimil_data/data_yuri/slp/list_coord",sep="\t",col.names=FALSE,row.names=FALSE)
#system("./plot_stations.gmt")
}












# read in GHCN version 3 (temperature)
read_ghcn <- function(timlim=c(syr,eyr)){

years <- timlim[1]:timlim[2]

stationfile <- paste(ghcntemppath,'ghcn_mm_t2m_adj.inv',sep='')

ww <- c(11, -1 , 8, -1, 9, -1, 6, -1, 30, -1, 4, 1, -1, 4, 2, 2, 2, 2, 1, 2, 16, 1)
stations <- read.fwf(stationfile, ww, header=F, fill=T, na.string='-999')
colnames(stations) <- c('ID', 'LATITUDE', 'LONGITUDE', 'STNELEV', 'NAME', 'GRELEV', 'POPCLS', 'POPSIZ', 'TOPO', 'STVEG', 'STLOC', 'OCNDIS', 'AIRSTN', 'TOWNDIS', 'GRVEG', 'POPCSS')
stations[stations == -9999 ] <- NA

# problem with orig ghcn v3 because data file has ids that are not in station header file

ww2 <- c(11, 4, 4, rep(c(5,1,1,1),12)) #5, 1, 1, 1, 82, 5, 1, 1, 1)
tmp <- as.matrix(read.fwf(paste(ghcntemppath,'ghcn_mm_t2m_adj.dat',sep=''), ww2, header=F, fill=T, na.string='-9999'))
colnames(tmp) <- c('ID', 'YEAR', 'ELEMENT', 'VALUE1', 'DMFLAG1', 'QCFLAG1', 'DSFLAG1', 'VALUE2', 'DMFLAG2', 'QCFLAG2', 'DSFLAG2', 'VALUE3', 'DMFLAG3', 'QCFLAG3', 'DSFLAG3', 'VALUE4', 'DMFLAG4', 'QCFLAG4', 'DSFLAG4', 'VALUE5', 'DMFLAG5', 'QCFLAG5', 'DSFLAG5', 'VALUE6', 'DMFLAG6', 'QCFLAG6', 'DSFLAG6', 'VALUE7', 'DMFLAG7', 'QCFLAG7', 'DSFLAG7', 'VALUE8', 'DMFLAG8', 'QCFLAG8', 'DSFLAG8', 'VALUE9', 'DMFLAG9', 'QCFLAG9', 'DSFLAG9', 'VALUE10', 'DMFLAG10', 'QCFLAG10', 'DSFLAG10', 'VALUE11', 'DMFLAG11', 'QCFLAG11', 'DSFLAG11', 'VALUE12', 'DMFLAG12', 'QCFLAG12', 'DSFLAG12')

# ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/v3/README

tmp2 <- sort(unique(tmp[,1]))
tmp3 <- sort(unique(stations[,1]))
tmp4 <- match(tmp2,tmp3)
tmp5 <- which(tmp4>0)
stat.i <- tmp2[tmp5]
tmp.arr <- array(NA,c(length(years), 12, length(stat.i)))

for (i in 1:length(seq(stat.i))){
   s.i <- as.numeric(which(tmp[,1] == stat.i[i]))
   y.i <- which(years %in% tmp[s.i,2])
   r.i <- apply(as.matrix(y.i), 1, function(x) min(which(tmp[s.i,2] == years[x])))
   tmp.arr[y.i,,i] <- tmp[s.i[r.i],c(4,8,12,16,20,24,28,32,36,40,44,48)]
 }

stations.new <- stations[which(stations[,'ID'] %in% stat.i),]

# check if at least 80% of data is available and NOT missing (NA)
mask <- apply(!is.na(tmp.arr), 3, sum) > 0.8 * prod(dim(tmp.arr)[1:2])
## check if at least 1% of data is available and NOT missing (NA)
#mask <- apply(!is.na(tmp.arr), 3, sum) > 0.01 * prod(dim(tmp.arr)[1:2])

ghcn.data.tmp <- array(aperm(tmp.arr[,,mask], c(2,1,3)), c(prod(dim(tmp.arr)[1:2]), sum(mask)))

ghcn.data=array(as.numeric(ghcn.data.tmp),c(dim(ghcn.data.tmp)[1],dim(ghcn.data.tmp)[2]))

ghcn <- list(data=ghcn.data/100, lon=stations.new[mask,'LONGITUDE'], lat=stations.new[mask,'LATITUDE'], names=gsub(' *$', '', as.character(stations.new[mask,'NAME'])), height=stations.new[mask,'STNELEV'], time=seq(min(years) + 1/24, by=1/12, length=nrow(ghcn.data)))
}




# read in GHCN version 3 (temperature)
read_ghcn_refyr <- function(syr,eyr,refsyr,refeyr){
  years <- syr:eyr
  ryrs <- refsyr:refeyr
  stationfile <- paste(ghcntemppath,'ghcn_mm_t2m_adj.inv',sep='')
  ww <- c(11, -1 , 8, -1, 9, -1, 6, -1, 30, -1, 4, 1, -1, 4, 2, 2, 2, 2, 1, 2, 16, 1)
  stations <- read.fwf(stationfile, ww, header=F, fill=T, na.string='-999')
  colnames(stations) <- c('ID', 'LATITUDE', 'LONGITUDE', 'STNELEV', 'NAME', 'GRELEV', 'POPCLS', 'POPSIZ', 'TOPO', 'STVEG', 'STLOC', 'OCNDIS', 'AIRSTN', 'TOWNDIS', 'GRVEG', 'POPCSS')
  stations[stations == -9999 ] <- NA
  # problem with orig ghcn v3 because data file has ids that are not in station header file
  ww2 <- c(11, 4, 4, rep(c(5,1,1,1),12)) #5, 1, 1, 1, 82, 5, 1, 1, 1)
  tmp <- as.matrix(read.fwf(paste(ghcntemppath,'ghcn_mm_t2m_adj.dat',sep=''), ww2, header=F, fill=T, na.string='-9999'))
  colnames(tmp) <- c('ID', 'YEAR', 'ELEMENT', 'VALUE1', 'DMFLAG1', 'QCFLAG1', 'DSFLAG1', 'VALUE2', 'DMFLAG2', 'QCFLAG2', 'DSFLAG2', 'VALUE3', 'DMFLAG3', 'QCFLAG3', 'DSFLAG3', 'VALUE4', 'DMFLAG4', 'QCFLAG4', 'DSFLAG4', 'VALUE5', 'DMFLAG5', 'QCFLAG5', 'DSFLAG5', 'VALUE6', 'DMFLAG6', 'QCFLAG6', 'DSFLAG6', 'VALUE7', 'DMFLAG7', 'QCFLAG7', 'DSFLAG7', 'VALUE8', 'DMFLAG8', 'QCFLAG8', 'DSFLAG8', 'VALUE9', 'DMFLAG9', 'QCFLAG9', 'DSFLAG9', 'VALUE10', 'DMFLAG10', 'QCFLAG10', 'DSFLAG10', 'VALUE11', 'DMFLAG11', 'QCFLAG11', 'DSFLAG11', 'VALUE12', 'DMFLAG12', 'QCFLAG12', 'DSFLAG12')
  # ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/v3/README
  tmp2 <- sort(unique(tmp[,1]))
  tmp3 <- sort(unique(stations[,1]))
  tmp4 <- match(tmp2,tmp3)
  tmp5 <- which(tmp4>0)
  stat.i <- tmp2[tmp5]
  tmp.arr <- array(NA,c(length(years), 12, length(stat.i)))
  for (i in 1:length(seq(stat.i))){
    s.i <- as.numeric(which(tmp[,1] == stat.i[i]))
    y.i <- which(years %in% tmp[s.i,2])
    r.i <- apply(as.matrix(y.i), 1, function(x) min(which(tmp[s.i,2] == years[x])))
    ry.i <- which(ryrs %in% tmp[s.i,2])
    if (length(ry.i) > 0) {
      print(paste(i,'of',length(seq(stat.i))))
      tmp.arr[y.i,,i] <- tmp[s.i[r.i],c(4,8,12,16,20,24,28,32,36,40,44,48)]
    } else {
#      print(i)
      if (length(y.i)==1) {
        tmp.arr[y.i,,i] <- array(NA,length(tmp[s.i[r.i],c(4,8,12,16,20,24,28,32,36,40,44,48)]))
      } else {
        tmp.arr[y.i,,i] <- array(NA,dim(tmp[s.i[r.i],c(4,8,12,16,20,24,28,32,36,40,44,48)]))
      }
    }  
  }
  stations.new <- stations[which(stations[,'ID'] %in% stat.i),]
  # check if at least 80% of data is available and NOT missing (NA)
  #mask <- apply(!is.na(tmp.arr), 3, sum) > 0.8 * prod(dim(tmp.arr)[1:2])
  ## check if at least 1% of data is available and NOT missing (NA)
  mask <- apply(!is.na(tmp.arr), 3, sum) > 0.01 * prod(dim(tmp.arr)[1:2])
  ghcn.data.tmp <- array(aperm(tmp.arr[,,mask], c(2,1,3)), c(prod(dim(tmp.arr)[1:2]), sum(mask)))
  ghcn.data=array(as.numeric(ghcn.data.tmp),c(dim(ghcn.data.tmp)[1],dim(ghcn.data.tmp)[2]))
  ghcn <- list(data=ghcn.data/100, lon=stations.new[mask,'LONGITUDE'], lat=stations.new[mask,'LATITUDE'], names=gsub(' *$', '', as.character(stations.new[mask,'NAME'])), height=stations.new[mask,'STNELEV'], time=seq(min(years) + 1/24, by=1/12, length=nrow(ghcn.data)))
}
# 
# # read in GHCN version 3 (temperature)
# read_ghcn_refyr <- function(syr,eyr,refsyr,refeyr){
#   years <- syr:eyr
#   ryrs <- refsyr:refeyr
#   stationfile <- paste(ghcntemppath,'ghcn_mm_t2m_adj.inv',sep='')
#   ww <- c(11, -1 , 8, -1, 9, -1, 6, -1, 30, -1, 4, 1, -1, 4, 2, 2, 2, 2, 1, 2, 16, 1)
#   stations <- read.fwf(stationfile, ww, header=F, fill=T, na.string='-999')
#   colnames(stations) <- c('ID', 'LATITUDE', 'LONGITUDE', 'STNELEV', 'NAME', 'GRELEV', 'POPCLS', 'POPSIZ', 'TOPO', 'STVEG', 'STLOC', 'OCNDIS', 'AIRSTN', 'TOWNDIS', 'GRVEG', 'POPCSS')
#   stations[stations == -9999 ] <- NA
#   # problem with orig ghcn v3 because data file has ids that are not in station header file
#   ww2 <- c(11, 4, 4, rep(c(5,1,1,1),12)) #5, 1, 1, 1, 82, 5, 1, 1, 1)
#   tmp <- as.matrix(read.fwf(paste(ghcntemppath,'ghcn_mm_t2m_adj.dat',sep=''), ww2, header=F, fill=T, na.string='-9999'))
#   colnames(tmp) <- c('ID', 'YEAR', 'ELEMENT', 'VALUE1', 'DMFLAG1', 'QCFLAG1', 'DSFLAG1', 'VALUE2', 'DMFLAG2', 'QCFLAG2', 'DSFLAG2', 'VALUE3', 'DMFLAG3', 'QCFLAG3', 'DSFLAG3', 'VALUE4', 'DMFLAG4', 'QCFLAG4', 'DSFLAG4', 'VALUE5', 'DMFLAG5', 'QCFLAG5', 'DSFLAG5', 'VALUE6', 'DMFLAG6', 'QCFLAG6', 'DSFLAG6', 'VALUE7', 'DMFLAG7', 'QCFLAG7', 'DSFLAG7', 'VALUE8', 'DMFLAG8', 'QCFLAG8', 'DSFLAG8', 'VALUE9', 'DMFLAG9', 'QCFLAG9', 'DSFLAG9', 'VALUE10', 'DMFLAG10', 'QCFLAG10', 'DSFLAG10', 'VALUE11', 'DMFLAG11', 'QCFLAG11', 'DSFLAG11', 'VALUE12', 'DMFLAG12', 'QCFLAG12', 'DSFLAG12')
#   # ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/v3/README
#   tmp2 <- sort(unique(tmp[,1]))
#   tmp3 <- sort(unique(stations[,1]))
#   tmp4 <- match(tmp2,tmp3)
#   tmp5 <- which(tmp4>0)
#   stat.i <- tmp2[tmp5]
#   tmp.arr <- array(NA,c(length(years), 12, length(stat.i)))
#   for (i in 1:length(seq(stat.i))){
#     s.i <- as.numeric(which(tmp[,1] == stat.i[i]))
#     y.i <- which(years %in% tmp[s.i,2])
#     r.i <- apply(as.matrix(y.i), 1, function(x) min(which(tmp[s.i,2] == years[x])))
#     ry.i <- which(ryrs %in% tmp[s.i,2])
#     if (length(ry.i) > 0) {
#       print(paste(i,'of',length(seq(stat.i))))
#       tmp.arr[y.i,,i] <- tmp[s.i[r.i],c(4,8,12,16,20,24,28,32,36,40,44,48)]
#     } else {
#       tmp.arr[y.i,,i] <- array(NA,dim(tmp[s.i[r.i],c(4,8,12,16,20,24,28,32,36,40,44,48)]))
#     }  
#   }
#   stations.new <- stations[which(stations[,'ID'] %in% stat.i),]
#   # check if at least 80% of data is available and NOT missing (NA)
#   #mask <- apply(!is.na(tmp.arr), 3, sum) > 0.8 * prod(dim(tmp.arr)[1:2])
#   ## check if at least 1% of data is available and NOT missing (NA)
#   mask <- apply(!is.na(tmp.arr), 3, sum) > 0.01 * prod(dim(tmp.arr)[1:2])
#   ghcn.data.tmp <- array(aperm(tmp.arr[,,mask], c(2,1,3)), c(prod(dim(tmp.arr)[1:2]), sum(mask)))
#   ghcn.data=array(as.numeric(ghcn.data.tmp),c(dim(ghcn.data.tmp)[1],dim(ghcn.data.tmp)[2]))
#   ghcn <- list(data=ghcn.data/100, lon=stations.new[mask,'LONGITUDE'], lat=stations.new[mask,'LATITUDE'], names=gsub(' *$', '', as.character(stations.new[mask,'NAME'])), height=stations.new[mask,'STNELEV'], time=seq(min(years) + 1/24, by=1/12, length=nrow(ghcn.data)))
# }






# read in GHCN version 2 (precipitation)
read_ghcn_refyr_precip <- function(syr,eyr,prsyr,preyr){
  years <- syr:eyr
  ryrs <- prsyr:preyr
  stationfile <- paste(ghcnprecippath,'v2.prcp.inv',sep='')
  ww <- c(11, -1, 19, -1, 10, -1, 6, -1, 7, -1, 5) 
  stations <- read.fwf(stationfile, ww, header=F, fill=T)
  colnames(stations) <- c('number', 'name', 'country', 'lat', 'lon', 'height')
  stations[stations == -9999 ] <- NA
  ww2 <- c(11,1, 4, rep(5,12))
  # read adjusted precip data 
  #tmp <- as.matrix(read.fwf(paste(ghcnprecippath,'v2.prcp_adj',sep=''), ww2, header=F, fill=T, na.string=c('-8888','-9999')))
  # read unadjusted precip data
  tmp <- as.matrix(read.fwf(paste(ghcnprecippath,'v2.prcp',sep=''), ww2, header=F, fill=T, na.string=c('-8888','-9999')))
  colnames(tmp) <- c('ID', 'TMP', 'YEAR', 'VALUE1', 'VALUE2', 'VALUE3', 'VALUE4', 'VALUE5', 'VALUE6', 'VALUE7', 'VALUE8', 'VALUE9', 'VALUE10', 'VALUE11', 'VALUE12')
  # ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/v2/v2.prcp.readme
  stat.i <- sort(unique(tmp[tmp[,'YEAR'] %in% ryrs,1])) ## speed up computation
  tmp.arr <- array(NA,c(length(years), 12, length(stat.i)),dimnames=c('yr','mon','stat'))
  stations.new <- stations[which(stations[,'number'] %in% stat.i),]
  stat.i <- stations.new[,1]
  # remove gaps in time series by filling with NA
  # search for longest contiguous piece with na.contiguous
  for (i in 1:length(seq(stat.i))){
    s.i <- as.numeric(which(tmp[,1] == stat.i[i]))
    y.i <- which(years %in% tmp[s.i,'YEAR'])
    r.i <- apply(as.matrix(y.i), 1, function(x) min(which(tmp[s.i,'YEAR'] == years[x])))
    ry.i <- which(ryrs %in% tmp[s.i,'YEAR'])
    if (length(ry.i) > 0) {
      tmp.arr[y.i,,i] <- tmp[s.i[r.i],grep('VALUE', colnames(tmp))]
    } 
  }
  # check if at least 80% of data is available and NOT missing (NA)
  #  mask <- apply(!is.na(tmp.arr), 3, sum) > 0.8 * prod(dim(tmp.arr)[1:2])
  ## check if at least 1% of data is available and NOT missing (NA)
  mask <- apply(!is.na(tmp.arr), 3, sum) > 0.01 * prod(dim(tmp.arr)[1:2])
  ghcn.data <- array(aperm(tmp.arr[,,mask], c(2,1,3)), c(prod(dim(tmp.arr)[1:2]), sum(mask)))
  ghcn_precip <- list(data=ghcn.data, lon=stations.new[mask,'lon'], lat=stations.new[mask,'lat'], names=gsub(' *$', '', as.character(stations.new[mask,'name'])), height=stations.new[mask,'height'], time=seq(min(years) + 1/24, by=1/12, length=nrow(ghcn.data)))
  return(ghcn_precip)
}
# read_ghcn_refyr_precip <- function(syr,eyr,prsyr,preyr){
#   years <- syr:eyr
#   ryrs <- prsyr:preyr
#   stationfile <- paste(ghcnprecippath,'v2.prcp.inv',sep='')
#   ww <- c(11, -1, 19, -1, 10, -1, 6, -1, 7, -1, 5)  
#   stations <- read.fwf(stationfile, ww, header=F, fill=T)
#   colnames(stations) <- c('number', 'name', 'country', 'lat', 'lon', 'height')
#   stations[stations == -9999 ] <- NA
#   ww2 <- c(11,1, 4, rep(5,12))
# # read adjusted precip data  
#   #tmp <- as.matrix(read.fwf(paste(ghcnprecippath,'v2.prcp_adj',sep=''), ww2, header=F, fill=T, na.string=c('-8888','-9999')))
# # read unadjusted precip data
#   tmp <- as.matrix(read.fwf(paste(ghcnprecippath,'v2.prcp',sep=''), ww2, header=F, fill=T, na.string=c('-8888','-9999')))
#   colnames(tmp) <- c('ID', 'TMP', 'YEAR', 'VALUE1', 'VALUE2', 'VALUE3', 'VALUE4', 'VALUE5', 'VALUE6', 'VALUE7', 'VALUE8', 'VALUE9', 'VALUE10', 'VALUE11', 'VALUE12')
# # ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/v2/v2.prcp.readme
#   stat.i <- sort(unique(tmp[,1]))
#   tmp.arr <- array(NA,c(length(years), 12, length(stat.i)),dimnames=c('yr','mon','stat'))
#   stations.new <- stations[which(stations[,'number'] %in% stat.i),]
#   stat.i <- stations.new[,1]
# # remove gaps in time series by filling with NA
# # search for longest contiguous piece with na.contiguous
#   for (i in 1:length(seq(stat.i))){
#     s.i <- as.numeric(which(tmp[,1] == stat.i[i]))
#     y.i <- which(years %in% tmp[s.i,3])
#     r.i <- apply(as.matrix(y.i), 1, function(x) min(which(tmp[s.i,3] == years[x])))
# #    tmp.arr[y.i,,i] <- tmp[s.i[r.i],seq(from=3, to=14)]
#     ry.i <- which(ryrs %in% tmp[s.i,3])
#     if (length(ry.i) > 0) {
#       tmp.arr[y.i,,i] <- tmp[s.i[r.i],seq(from=3, to=14)]
#     } else {
#       tmp.arr[y.i,,i] <- array(NA,dim(tmp[s.i[r.i],seq(from=3, to=14)]))
#     }  
#   } 
# # check if at least 80% of data is available and NOT missing (NA)
# #  mask <- apply(!is.na(tmp.arr), 3, sum) > 0.8 * prod(dim(tmp.arr)[1:2])
# ## check if at least 1% of data is available and NOT missing (NA)
#   mask <- apply(!is.na(tmp.arr), 3, sum) > 0.01 * prod(dim(tmp.arr)[1:2])
#   ghcn.data <- array(aperm(tmp.arr[,,mask], c(2,1,3)), c(prod(dim(tmp.arr)[1:2]), sum(mask)))
#   ghcn_precip <- list(data=ghcn.data, lon=stations.new[mask,'lon'], lat=stations.new[mask,'lat'], name=gsub(' *$', '', as.character(stations.new[mask,'name'])), height=stations.new[mask,'height'], time=seq(min(years) + 1/24, by=1/12, length=nrow(ghcn.data)))
# }





# read in GHCN version 2 (precipitation)
read_ghcn_precip <- function(timlim=c(syr,eyr)){
  years <- timlim[1]:timlim[2]
  stationfile <- paste(ghcnprecippath,'v2.prcp.inv',sep='')
  ww <- c(11, -1, 19, -1, 10, -1, 6, -1, 7, -1, 5)  
  stations <- read.fwf(stationfile, ww, header=F, fill=T)
  colnames(stations) <- c('number', 'name', 'country', 'lat', 'lon', 'height')
  stations[stations == -9999 ] <- NA
  ww2 <- c(11,1, 4, rep(5,12))
  tmp <- as.matrix(read.fwf(paste(ghcnprecippath,'v2.prcp_adj',sep=''), ww2, header=F, fill=T, na.string=c('-8888','-9999')))
  colnames(tmp) <- c('ID', 'TMP', 'YEAR', 'VALUE1', 'VALUE2', 'VALUE3', 'VALUE4', 'VALUE5', 'VALUE6', 'VALUE7', 'VALUE8', 'VALUE9', 'VALUE10', 'VALUE11', 'VALUE12')
  # ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/v2/v2.prcp.readme
  stat.i <- sort(unique(tmp[,1]))
  tmp.arr <- array(NA,c(length(years), 12, length(stat.i)),dimnames=c('yr','mon','stat'))
  stations.new <- stations[which(stations[,'number'] %in% stat.i),]
  stat.i <- stations.new[,1]
  
  # remove gaps in time series by filling with NA
  # search for longest contiguous piece with na.contiguous
  for (i in 1:length(seq(stat.i))){
    s.i <- as.numeric(which(tmp[,1] == stat.i[i]))
    y.i <- which(years %in% tmp[s.i,3])
    r.i <- apply(as.matrix(y.i), 1, function(x) min(which(tmp[s.i,3] == years[x])))
    tmp.arr[y.i,,i] <- tmp[s.i[r.i],seq(from=3, to=14)]
  } 
  
  # check if at least 1% of data is available and NOT missing (NA)
  mask <- apply(!is.na(tmp.arr), 3, sum) > 0.01 * prod(dim(tmp.arr)[1:2])
  
  ghcn.data <- array(aperm(tmp.arr[,,mask], c(2,1,3)), c(prod(dim(tmp.arr)[1:2]), sum(mask)))
  ghcn_precip <- list(data=ghcn.data, lon=stations.new[mask,'lon'], lat=stations.new[mask,'lat'], names=gsub(' *$', '', as.character(stations.new[mask,'name'])), height=stations.new[mask,'height'], time=seq(min(years) + 1/24, by=1/12, length=nrow(ghcn.data)))
}







# read in HISTALP
read_histalp <- function(timlim=c(syr,eyr)){

years <- timlim[1]:timlim[2]

stationfile <- paste(histalppath,'histalp_temp.csv',sep='')

tmp=read.table(stationfile,header=T, na.string='999999')
colnames(tmp)=c("id","k_nam","name","l_code","lon","lat","hoehe","datum","jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec","spr","sum","aut","win","apr-sep","oct-mar","ann")

stat.i <- sort(unique(tmp[,1]))
tmp.arr <- array(NA,c(length(years), 12, length(stat.i)))
lon.arr <- array(NA,length(stat.i))
lat.arr <- array(NA,length(stat.i))
elev.arr <- array(NA,length(stat.i))
name.arr <- array(NA,length(stat.i))

for (i in 1:length(seq(stat.i))){
   s.i <- as.numeric(which(tmp[,1] == stat.i[i]))
   y.i <- which(years %in% tmp[s.i,'datum'])
   if (length(y.i) > 0) {
     r.i <- apply(as.matrix(y.i), 1, function(x) min(which(tmp[s.i,'datum'] == years[x])))
     tmp.arr[y.i,,i] <- as.matrix(tmp[s.i[r.i],seq(9,20)])
     lon.arr[i] <- tmp[s.i[1],'lon']
     lat.arr[i] <- tmp[s.i[1],'lat']
     elev.arr[i] <- tmp[s.i[1],'hoehe']
     name.arr[i] <- as.character(tmp[s.i[1],'name'])
   }   
}

## check if at least 80% of data is available and NOT missing (NA)
#mask <- apply(!is.na(tmp.arr), 3, sum) > 0.8 * prod(dim(tmp.arr)[1:2])
# check if at least 1% of data is available and NOT missing (NA)
mask <- apply(!is.na(tmp.arr), 3, sum) > 0.01 * prod(dim(tmp.arr)[1:2])


histalp.data.tmp <- array(aperm(tmp.arr[,,mask], c(2,1,3)), c(prod(dim(tmp.arr)[1:2]), sum(mask)))
histalp.data=array(as.numeric(histalp.data.tmp),c(dim(histalp.data.tmp)[1],dim(histalp.data.tmp)[2]))

histalp <- list(data=histalp.data/10, lon=lon.arr[mask], lat=lat.arr[mask], names=name.arr[mask], height=elev.arr[mask], time=seq(min(years) + 1/24, by=1/12, length=nrow(histalp.data)))
}






# read in OTHER_INSTR
# read_other_instr <- function(timlim=c(syr,eyr)){
# 
#   /data/climdata/instr/a_moberg_cruts2/brussels-1794-2001.dat
#   /data/climdata/instr/early_instr_collection/temp_pre1850/cent_engl_1659-2012_monthly.dat  
#   /data/climdata/instr/auchmann/...
#   
# }









extract_lonlat <- function(x){
    if (!is.na(x)){
        xout <- gsub('deg.*', '', x)
        if (length(grep('-', xout)) == 1){
            xtmp <- strsplit(xout, '-')
            xout <- mean(as.numeric(unlist(xtmp)))
        }
        xout <- as.numeric(xout)
        if (length(grep('min', x)) == 1){
            xadd <- as.numeric(gsub('min.*', '', gsub('.*deg *', '', x)))
            xout <- xout + xadd/60
        }
    } else {    
        xout <- NA
    }
    xout
}        

compute_lonlat <- function(x){
    dir <- gsub('[-. a-z0-9]*', '', x)
    n.i <- dir[,1] %in% c('N', 'S')
    n.i <- cbind(n.i + 1, (!n.i) + 1)
    s.i <- (dir %in% c('E', 'N'))*2 - 1
    xout <- s.i * apply(x, 1:2, extract_lonlat)
    
    # swap x accordingly
    x2 <- cbind(xout[cbind(1:84, n.i[,1])], xout[cbind(1:84, n.i[,2])])
    # remove lons > 180
    lon.i <- x2[,1] > 180 & !is.na(x2[,1])
    x2[lon.i,1] <- x2[lon.i,1] - 360
    xout <- data.frame(lon=x2[,1], lat=x2[,2])    
    return(xout)
}

extra_lonlat <- function(x){
    if (all(!is.na(x))){
        xout <- gsub('deg.*', '', x)
        if (length(grep('min', x)) > 0){
            xadd <- as.numeric(gsub('min.*', '', gsub('.*deg *', '', x)))
            xout <- as.numeric(xout) + xadd/60
        }
        if (length(grep('-', xout)) > 0){
            xtmp <- strsplit(xout, '-')
            if (is.na(xtmp[[1]][2])) { xtmp[[1]][2]=xtmp[[1]][1] }
            if (is.na(xtmp[[2]][2])) { xtmp[[2]][2]=xtmp[[2]][1] }
            xout <- rbind(xtmp[[1]][c(1,1,2,2,1)], xtmp[[2]][c(1,2,2,1,1)])
        }
        if (is.matrix(xout)) {
            xout <- array(as.numeric(xout), dim(xout))
        } else {
            xout <- as.numeric(xout)
        }
        dir <- gsub(' ', '', gsub('[-. a-z0-9]*', '', x))
        n.i <- dir %in% c('N', 'S') + 1
        s.i <- (dir %in% c('E', 'N'))*2 - 1
        xout <- s.i * xout
        if (is.matrix(xout)){
            xout <- xout[n.i,]
        } else {
            xout <- xout[n.i]
        }
    } else {    
        xout <- c(NA,NA)
    }
    xout
}        


# read in proxy locations
read_proxy_loc <- function(){
tmp <- read.table(paste0(dataintdir,'proxies_EnSRF/temperature_proxies.csv'),
                  sep=',', header=TRUE, na.string=c(' -', '-'))
lon <- as.character(tmp[[grep('Longitude', names(tmp))]])
lat <- as.character(tmp[[grep('Latitude', names(tmp))]])
lola <- list()
for (i in seq(along=lon)) lola[[i]] <- c(lon[i], lat[i])
proxy <- lapply(lola, extra_lonlat)
proxies <- compute_lonlat(cbind(lon, lat))
# remove duplicates
proxies <- proxies[!duplicated(proxies) & !is.na(proxies[,1]),]
}



# function to create climatologies for separated by seasons or months
multiyear_seas_average <- function(x, seas){
  # seas e.g. 2 for half year means nov-apr and may-oct
  time    <- x$time
  year    <- floor(time)
  month   <- round((time %% 1 + 1/24)*12)
  # leave only same month/season in 2. dim and create 3. dim for each months/seasons
  sdata   <- aperm(array(x$data,dim=c(dim(x$data)[1],seas,(dim(x$data)[2]/seas))),c(1,3,2))
  sdatamean <- apply(sdata,c(1,3),mean,na.rm=T)
  out <- x
  out$data <- sdatamean
  out$time <- month[1:seas]
  invisible(out)
}









seasonal_average <- function(x, seas){
    # simplified using SMA (moving average)
    # assuming that no data is missing
    # seas e.g. c(4,10) means half year means nov-apr and may-oct
    time    <- x$time
    year    <- floor(time)
    month   <- round((time %% 1 + 1/24)*12)
    nseas   <- 12/length(seas)
    sseas   <- min(which(month == seas[1])[which(month == seas[1]) > nseas])
    ind     <- seq(sseas, max(which(month == seas[length(seas)])), by=nseas)
    index   <- rep(FALSE, length(time))
    index[ind] <- TRUE
    xout    <- lapply(x, seas_ave, n=nseas, index=index)
    xout$index <- which(index)
    invisible(xout)
  }



## seas_ave <- function(x, n=mseas, index=index){
##   dimx    <- dim(x)
##   if (length(dimx) <= 1) {
##     dimx  <- length(x)
##   }
##   if (any(dimx == length(index))){
##     if (length(dimx) <= 1){
##       out     <- SMA(x, n=n)[index]
##     } else {
##       timei   <- max(which(dimx == length(index)))
##       out     <- apply(x, setdiff(seq(along=dimx),timei), SMA, n=n)
##       out     <- array(out[rep(index, prod(dimx[setdiff(seq(dimx),timei)]))],
##                        c(sum(index), dimx[setdiff(seq(dimx),timei)]))
##       outperm <- rep(NA, length(dim(out)))
##       outperm[timei] <- 1
##       outperm[is.na(outperm)] <- 2:length(dim(out))
##       out     <- aperm(out, outperm)
##     }
##   } else {
##     out <- x
##   }
##   invisible(out)
## }

seas_ave <- function(x, n=mseas, index=index){
  dimx    <- dim(x)
  if (length(dimx) <= 1) {
    dimx  <- length(x)
  }
  if (any(dimx == length(index))){
    if (length(dimx) <= 1){
      if (any(is.na(x))){
        index2 <- rep(FALSE, length(index))
        index2[seq(1,length(index)-floor(n/2))] <- index[seq(floor(n/2)+1, length(index))]
        out <- filter(x, rep(1/n, n))[index2]
      } else {
        out     <- SMA(x, n=n)[index]
      }
    } else {
      timei   <- max(which(dimx == length(index)))
      if (any(is.na(x))) {
        index2 <- rep(FALSE, length(index))
        index2[seq(1,length(index)-floor(n/2))] <- index[seq(floor(n/2)+1, length(index))]
        out     <- apply(x, setdiff(seq(along=dimx),timei), filter, rep(1/n,n))
        out     <- array(out[rep(index2, prod(dimx[setdiff(seq(dimx),timei)]))],
                       c(sum(index2), dimx[setdiff(seq(dimx),timei)]))
      } else {
        out     <- apply(x, setdiff(seq(along=dimx),timei), SMA, n=n)
        out     <- array(out[rep(index, prod(dimx[setdiff(seq(dimx),timei)]))],
                       c(sum(index), dimx[setdiff(seq(dimx),timei)]))
      }
      outperm <- rep(NA, length(dim(out)))
      outperm[timei] <- 1
      outperm[is.na(outperm)] <- 2:length(dim(out))
      out     <- aperm(out, outperm)
    }
  } else {
    out <- x
  }
  invisible(out)
}




plot_echam <- function(x, levs, varname='temp2', type='data',  ti=1, 
                       stations=NULL,statpanel=NULL, #centercol=NULL,
                       symmetric=T, siglev=0.05, names=NULL, main='', units='', 
                       cex.pt=2, st.cex=cex.pt*0.25, st.col=NULL, add=FALSE, 
                       cols=NULL, lwd=3, lty=1, seas=c(1,2), xlab='', ylab=NULL, 
                       rownames=NULL, colnames=NULL, latlim=c(-90,90),lonlim=c(-180,180),
                       addcontours=F, contvarname='gph500', conttype='data', contcol='black',
                       contlev=levs,wcol='black',colorbar=T){
  oldpar <- par(no.readonly=TRUE)
  # plot specific lon/lat range (Joerg 2016/08)
  for (i in 1:4) {
    if (i==1) {
      pos <- which(x$lat>latlim[1])
    } else if (i==2){
      pos <- which(x$lat<latlim[2])
    } else if (i==3){
      pos <- which(x$lon>lonlim[1])
    } else if (i==4){
      pos <- which(x$lon<lonlim[2])
    } 
    if (length(dim(x$data)) == 3){
      x$data <- x$data[pos,,,drop=F]
    } else {
      x$data <- x$data[pos,,drop=F]
    }
    x$ensmean <- x$ensmean[pos,,drop=F] 
    x$lon <- x$lon[pos]
    x$lat <- x$lat[pos]
    x$names <- x$names[pos]
  }
  # plot specific lat range  
  # latpos <- match(which(x$lat>latlim[1]),which(x$lat<latlim[2]))
  # if (length(dim(x$data)) == 3){
  #   x$data <- x$data[latpos,,,drop=F]
  # } else {
  #   x$data <- x$data[latpos,,drop=F]
  # }
  # x$ensmean <- x$ensmean[latpos,,drop=F] 
  # x$lon <- x$lon[latpos]
  # x$lat <- x$lat[latpos]
  # x$names <- x$names[latpos]
  dat.i <- which(x$names == varname)
  # Jonas Jan. 2014 for new ECHAM format 
  lons <- x$lon[x$names == varname]
  lats <- x$lat[x$names == varname]
  if (length(dat.i) == 1) ti <- 1:ncol(x$data)
  if (type == 'data'){
    if (length(dim(x$data)) == 3){
      if (dim(x$data)[3]==1) {
        plotdata <- array(x$data[dat.i,ti,],c(length(dat.i),1)) 
      } else {
        plotdata <- x$data[dat.i,ti,] # changed back to original 07/2015 
      }
    } else {
      #        plotdata <- t(x$data[dat.i,ti,drop=F]) 
      plotdata <- x$data[dat.i,ti,drop=F]
    }
  } else if (type == 'ensmean') {
    if (is.matrix(x$ensmean)){
      plotdata <- as.matrix(x$ensmean[dat.i,ti])
    } else {
      plotdata <- as.matrix(x$ensmean[dat.i])
    }
  } else {
    stop("Only 'data' and 'ensmean' are implemented as type")
  }
  if (length(dat.i) == 1){
    if (length(seas) > 1) {
      if (is.null(cols)){
        se.col <- rbfun(14)[c(3,12)]
      } else {
        se.col <- cols[seq(along=seas)]
      }
    } else {
      if (is.null(cols)){
        se.col <- rep(1,max(seas))
      } else {
        se.col <- rep(cols[1], max(seas))
      }
    }
    tis <- as.vector(apply(as.matrix(seas), 1, function(x,y) seq(x, length(y$time),2), y=x))
    if (ncol(plotdata) == 1) {
      if (!add) plot(0, type='n', xlim=range(x$time), ylim=range(plotdata[tis]),
                     xlab=xlab, ylab=if (is.null(ylab)) varname else ylab)
      for (se in 1:2) lines(x$time[seq(se,length(plotdata),2)], 
                        plotdata[seq(se,length(plotdata),2)], col=se.col[se], lwd=lwd, lty=lty)
    } else {
      d.range <- apply(plotdata, 1, range)
      if (!add) plot(0, type='n', xlim=range(x$time), ylim=range(d.range[,tis]),
                     xlab='', ylab=if (is.null(ylab)) varname else ylab)
      for (se in seas){
        for (dr in 1:2){
          lines(x$time[seq(se, ncol(d.range), 2)], d.range[dr,seq(se, ncol(d.range), 2)], 
                col=se.col[se], lwd=lwd/2, lty=lty)
        }
      }
    }
  } else {
    nmod <- ncol(plotdata)
    nrow <- ceiling(nmod/4)
    ncol <- ceiling(nmod/nrow)+1
    nn <- nrow*(ncol-1)
    if (!add){
      dd <- lcm(2)
      if (nmod == 1) dd <- lcm(2.8)
      layout(matrix(c(1:nmod, rep(nmod+2, nn-nmod), rep(nmod+1, nrow)),
                    ncol, nrow, byrow=T), height=c(rep(5, ncol-1),dd))
      oma <- rep(0, 4)
      if (!is.null(rownames)) oma[3] <- 3
      if (!is.null(colnames)) oma[2] <- 3
      par(oma=oma)
    }
    par(mar=rep(0.5,4), cex.axis=1.4, cex.lab=1.4)
    
    if (missing(levs)) {
      if (symmetric){
        levs <- pretty(c(plotdata, -plotdata), 15)
      } else {
        levs <- pretty(plotdata, 15)
      }
    }
    if (is.null(cols)){
      if (varname == 'precip'){
        cols <- gbfun(length(levs) - 1)
      } else { 
        cols <- rbfun(length(levs) - 1)
      }
    } else {
      cols <- rep(cols, length.out=length(levs) - 1)
    }
    for (i in 1:nmod){
      plot(0, type='n', xlim=lonlim, ylim=latlim,
           axes=F, xlab='', ylab='')
      li <- as.numeric(cut(plotdata[,i], breaks=levs))
      points(x$lon, x$lat, pch=15, col=cols[li], cex=cex.pt)
      if (!is.null(colnames)) {
        if (i <= length(colnames)) axis(3, at=mean(range(x$lon)), label=colnames[i], 
                                        tick=FALSE, cex=1.4)
      }
      if (!is.null(rownames)) { 
        if (i%%length(colnames) == 1) {
          axis(2, at=mean(range(x$lat)), 
               label=rownames[ceiling(i/length(colnames))], tick=FALSE, cex=1.4)
        }
      }
      world <- map(interior=F, plot=F)
      world$x[abs(diff(world$x)) > 180] <- NA
      lines(world,col=wcol)
      #map(add=T, region=c('Caspian','Great Lakes'),col=wcol) 
      if (addcontours) {
        add_contour(x, cvarname=contvarname, ctype=conttype, ccol=contcol, 
                    clev=contlev, cont_i=i)
      }
      if ((!is.null(stations)) & (any(statpanel==i))) {
        statlon <- stations$lon[!is.na(stations$data[,i])]
        statlat <- stations$lat[!is.na(stations$data[,i])]
        if (is.null(st.col)){
          st <- stations$data[,i]
          lli <- as.numeric(cut(st, breaks=levs))
          try(points(stations$lon[stations$names=='temp2' & !is.na(stations$data[,i])], 
                     stations$lat[stations$names=='temp2' & !is.na(stations$data[,i])], 
                     pch=1, col=rgb(10,0,0,8,maxColorValue=10), cex=st.cex*2))
          try(points(stations$lon[stations$names=='precip' & !is.na(stations$data[,i])], 
                     stations$lat[stations$names=='precip' & !is.na(stations$data[,i])], 
                     pch=3, col=rgb(10,0,10,8,maxColorValue=10), cex=st.cex))
          try(points(stations$lon[stations$names=='slp' & !is.na(stations$data[,i])], 
                     stations$lat[stations$names=='slp' & !is.na(stations$data[,i])], 
                     pch=4, col=rgb(0,0,10,8,maxColorValue=10), cex=st.cex))
          try(points(stations$lon[stations$names=='prox' & !is.na(stations$data[,i])], 
                     stations$lat[stations$names=='prox' & !is.na(stations$data[,i])], 
                     pch=3, col=rgb(0,10,0,8,maxColorValue=10), cex=st.cex))
        } else {
          points(statlon, statlat, pch=1, col=st.col, cex=st.cex)
        }
      }
      if (!is.null(x$p.value)){
        signif <- x$p.value[,ti,i] > siglev
        points(x$lon[signif], x$lat[signif], pch=4, cex=max(2, cex.pt*0.75))
      }
      if (!is.null(names)) text(min(lonlim), max(latlim), names[i], adj=c(0.5,1), cex=1.4)
      box()
    }
    if (colorbar) {
      if (units == ''){
        par(mar=c(3,0.5,1.5,0.5), cex.axis=1.4, cex.lab=1.4)    
        plot_colourbar(col=cols, lev=levs, cex=1.4, xlab=units)
      } else {
        par(mar=c(3,8,1.5,8), cex.axis=1.4, cex.lab=1.4)  
        tmp <- list(col=cols, lev=levs, units=units)
        class(tmp) <- 'plotmap'  
        plot_colourbar(tmp, cex=1.4)
      }
    }
    mtext(main, side=3, line=2, outer=TRUE, font=1, cex=1.4)
    if (!add) par(oldpar)
  }
}




# plot_echam3 based on plot_echam (NOT version 2!) to overlay wind vectors
plot_echam3 <- function(x, levs, varname='temp2', type='data',  ti=1, 
                       stations=NULL,statpanel=NULL, #centercol=NULL,
                       symmetric=T, siglev=0.05, names=NULL, main='', units='', 
                       cex.pt=2, st.cex=cex.pt*0.25, st.col=NULL, add=FALSE, 
                       cols=NULL, lwd=3, lty=1, seas=c(1,2), xlab='', ylab=NULL, 
                       rownames=NULL, colnames=NULL, latlim=c(-90,90),lonlim=c(-180,180),
                       addcontours=F, contvarname='gph500', conttype='data', contcol='black',
                       contlev=levs, addvectors=F, vecnames=NULL, #vectortype='data',
                       veccol='black', veclen=0.03, vecscale=0.3, vecwd=0.75, every_x_vec=4,
                       wcol='black', colorbar=T){
  oldpar <- par(no.readonly=TRUE)
  # plot specific lat range  
  # latpos <- match(which(x$lat>latlim[1]),which(x$lat<latlim[2]))
  for (i in 1:4) {
    if (i==1) {
      pos <- which(x$lat>latlim[1])
    } else if (i==2){
      pos <- which(x$lat<latlim[2])
    } else if (i==3){
      pos <- which(x$lon>lonlim[1])
    } else if (i==4){
      pos <- which(x$lon<lonlim[2])
    } 
    if (length(dim(x$data)) == 3){
      x$data <- x$data[pos,,,drop=F]
    } else {
      x$data <- x$data[pos,,drop=F]
    }
    x$ensmean <- x$ensmean[pos,,drop=F] 
    x$lon <- x$lon[pos]
    x$lat <- x$lat[pos]
    x$names <- x$names[pos]
  }
#   lonpos <- match(which(x$lon>lonlim[1]),which(x$lon<lonlim[2]))
#   if (length(dim(x$data)) == 3){
#     x$data <- x$data[lonpos,,,drop=F]
#   } else {
#     x$data <- x$data[lonpos,,drop=F]
#   }
#   x$ensmean <- x$ensmean[lonpos,,drop=F] 
#   x$lon <- x$lon[lonpos]
#   x$lat <- x$lat[lonpos]
#   x$names <- x$names[lonpos]
  dat.i <- which(x$names == varname)
  # Jonas Jan. 2014 for new ECHAM format 
  lons <- x$lon[x$names == varname]
  lats <- x$lat[x$names == varname]
  if (length(dat.i) == 1) ti <- 1:ncol(x$data)
  if (type == 'data'){
    if (length(dim(x$data)) == 3){
      if (dim(x$data)[3]==1) {
        plotdata <- array(x$data[dat.i,ti,],c(length(dat.i),1)) 
      } else {
        plotdata <- x$data[dat.i,ti,] # changed back to original 07/2015 
      }
    } else {
      #      plotdata <- t(x$data[dat.i,ti,drop=F]) 
      plotdata <- x$data[dat.i,ti,drop=F]
    }
  } else if (type == 'ensmean') {
    if (is.matrix(x$ensmean)){
      plotdata <- as.matrix(x$ensmean[dat.i,ti])
    } else {
      plotdata <- as.matrix(x$ensmean[dat.i])
    }
  } else {
    stop("Only 'data' and 'ensmean' are implemented as type")
  }
  if (length(dat.i) == 1){
    if (length(seas) > 1) {
      if (is.null(cols)){
        se.col <- rbfun(14)[c(3,12)]
      } else {
        se.col <- cols[seq(along=seas)]
      }
    } else {
      if (is.null(cols)){
        se.col <- rep(1,max(seas))
      } else {
        se.col <- rep(cols[1], max(seas))
      }
    }
    tis <- as.vector(apply(as.matrix(seas), 1, function(x,y) seq(x, length(y$time),2), y=x))
    if (ncol(plotdata) == 1) {
      if (!add) plot(0, type='n', xlim=range(x$time), ylim=range(plotdata[tis]),
                     xlab=xlab, ylab=if (is.null(ylab)) varname else ylab)
      for (se in 1:2) lines(x$time[seq(se,length(plotdata),2)], 
                            plotdata[seq(se,length(plotdata),2)], col=se.col[se], lwd=lwd, lty=lty)
    } else {
      d.range <- apply(plotdata, 1, range)
      if (!add) plot(0, type='n', xlim=range(x$time), ylim=range(d.range[,tis]),
                     xlab='', ylab=if (is.null(ylab)) varname else ylab)
      for (se in seas){
        for (dr in 1:2){
          lines(x$time[seq(se, ncol(d.range), 2)], d.range[dr,seq(se, ncol(d.range), 2)], 
                col=se.col[se], lwd=lwd/2, lty=lty)
        }
      }
    }
  } else {
    nmod <- ncol(plotdata)
    nrow <- ceiling(nmod/4)
    ncol <- ceiling(nmod/nrow)+1
    nn <- nrow*(ncol-1)
    if (!add){
      dd <- lcm(2)
      if (nmod == 1) dd <- lcm(2.8)
      layout(matrix(c(1:nmod, rep(nmod+2, nn-nmod), rep(nmod+1, nrow)),
                    ncol, nrow, byrow=T), height=c(rep(5, ncol-1),dd))
      oma <- rep(0, 4)
      if (!is.null(rownames)) oma[3] <- 3
      if (!is.null(colnames)) oma[2] <- 3
      par(oma=oma)
    }
    par(mar=rep(0.5,4), cex.axis=1.4, cex.lab=1.4)
    
    if (missing(levs)) {
      if (symmetric){
        levs <- pretty(c(plotdata, -plotdata), 15)
      } else {
        levs <- pretty(plotdata, 15)
      }
    }
    if (is.null(cols)){
      if (varname == 'precip'){
        cols <- gbfun(length(levs) - 1)
      } else { 
        cols <- rbfun(length(levs) - 1)
      }
    } else {
      cols <- rep(cols, length.out=length(levs) - 1)
    }
    #    if (!is.null(centercol)) {
    #      if (is.even(length(cols)/2)) {
    #        cols[(length(cols)/2):((length(cols)/2)+1)] <- centercol
    #      } else {
    #        cols[ceiling(length(cols)/2)] <- centercol
    #      }
    #    }  
    for (i in 1:nmod){
      #      plot(0, type='n', xlim=range(x$lon), ylim=range(x$lat),
      #           axes=F, xlab='', ylab='')      
      plot(0, type='n', xlim=lonlim, ylim=latlim,
           axes=F, xlab='', ylab='')
      li <- as.numeric(cut(plotdata[,i], breaks=levs))
      points(x$lon, x$lat, pch=15, col=cols[li], cex=cex.pt)
      if (!is.null(colnames)) {
        if (i <= length(colnames)) axis(3, at=mean(range(x$lon)), label=colnames[i], 
                                        tick=FALSE, cex=1.4)
      }
      if (!is.null(rownames)) { 
        if (i%%length(colnames) == 1) {
          axis(2, at=mean(range(x$lat)), 
               label=rownames[ceiling(i/length(colnames))], tick=FALSE, cex=1.4)
        }
      }
      world <- map(interior=F, plot=F)
      world$x[abs(diff(world$x)) > 180] <- NA
      lines(world,col=wcol)
      #map(add=T, region=c('Caspian','Great Lakes'),col=wcol) 
      if (addcontours) {
        add_contour(x, cvarname=contvarname, ctype=conttype, ccol=contcol, clev=contlev, cont_i=i)
      }
      if (addvectors) {
        if (type == 'ensmean') {
          uw <- x$ensmean[which(x$names==vecnames[1]),i]
          vw <- x$ensmean[which(x$names==vecnames[2]),i]
        } else if (type == 'data') {
          uw <- x$data[which(x$names==vecnames[1]),1,i]
          vw <- x$data[which(x$names==vecnames[2]),1,i]
        }
        uvlon <- x$lon[which(x$names==vecnames[1])]
        uvlat <- x$lat[which(x$names==vecnames[2])]
        uvpos <- !is.na(match(uvlon,unique(uvlon)[seq(1,length(unique(uvlon)),every_x_vec)]))&
                 !is.na(match(uvlat,unique(uvlat)[seq(1,length(unique(uvlat)),every_x_vec)]))
#        plot(0, type='n', xlim=lonlim, ylim=latlim, axes=F, xlab='', ylab='')
#        points(x$lon, x$lat, pch=15, col=rep("white",length(x$lon)), cex=1)
        quiver(x=uvlon[uvpos], y=uvlat[uvpos], u=uw[uvpos], v=vw[uvpos],
             scale=vecscale,length=veclen,col=veccol,lwd=vecwd,lend=1)
#        quiver(x=x$lon, y=x$lat, u=uw, v=vw,
#               scale=vecscale,length=veclen,col=veccol,lwd=vecwd,lend=1)
      }
      if ((!is.null(stations)) & (any(statpanel==i))) {
        statlon <- stations$lon[!is.na(stations$data[,i])]
        statlat <- stations$lat[!is.na(stations$data[,i])]
        if (is.null(st.col)){
          st <- stations$data[,i]
          lli <- as.numeric(cut(st, breaks=levs))
          #          points(statlon, statlat, pch='.', bg=cols[lli], cex=st.cex*0.01)
          #          tst <- stations$data[stations$name=='temp',i]
          #          pst <- stations$data[stations$name=='precip',i]
          #          sst <- stations$data[stations$name=='slp',i]
          try(points(stations$lon[stations$names=='temp2' & !is.na(stations$data[,i])], 
                     stations$lat[stations$names=='temp2' & !is.na(stations$data[,i])], 
                     pch=1, col=rgb(10,0,0,8,maxColorValue=10), cex=st.cex*2))
          try(points(stations$lon[stations$names=='precip' & !is.na(stations$data[,i])], 
                     stations$lat[stations$names=='precip' & !is.na(stations$data[,i])], 
                     pch=3, col=rgb(10,0,10,8,maxColorValue=10), cex=st.cex))
          try(points(stations$lon[stations$names=='slp' & !is.na(stations$data[,i])], 
                     stations$lat[stations$names=='slp' & !is.na(stations$data[,i])], 
                     pch=4, col=rgb(0,0,10,8,maxColorValue=10), cex=st.cex))
          try(points(stations$lon[stations$names=='prox' & !is.na(stations$data[,i])], 
                     stations$lat[stations$names=='prox' & !is.na(stations$data[,i])], 
                     pch=3, col=rgb(0,10,0,8,maxColorValue=10), cex=st.cex))
        } else {
          points(statlon, statlat, pch=1, col=st.col, cex=st.cex)
        }
      }
      if (!is.null(x$p.value)){
        signif <- x$p.value[,ti,i] > siglev
        points(x$lon[signif], x$lat[signif], pch=4, cex=max(2, cex.pt*0.75))
      }
      #      if (!is.null(names)) text(min(x$lon), max(x$lat), names[i], adj=c(0.5,1), cex=1.4)
      if (!is.null(names)) text(min(lonlim), max(latlim), names[i], adj=c(0.5,1), cex=1.4)
      box()
    }
    if (colorbar) {
      if (units == ''){
        par(mar=c(3,0.5,1.5,0.5), cex.axis=1.4, cex.lab=1.4)    
        plot_colourbar(col=cols, lev=levs, cex=1.4, xlab=units)
      } else {
        par(mar=c(3,8,1.5,8), cex.axis=1.4, cex.lab=1.4)  
        tmp <- list(col=cols, lev=levs, units=units)
        class(tmp) <- 'plotmap'  
        plot_colourbar(tmp, cex=1.4)
      }
    }
    mtext(main, side=3, line=2, outer=TRUE, font=1, cex=1.4)
    if (!add) par(oldpar)
  }
}




# plot_echam4 based on plot_echam3 to overlay wind vectors
# + additional zonal mean on side of plot, layout need to be adjusted
plot_echam4 <- function(x, levs, varname='temp2', type='data',  ti=1, 
                        stations=NULL,statpanel=NULL, #centercol=NULL,
                        symmetric=T, siglev=0.05, names=NULL, main='', units='', 
                        cex.pt=2, st.cex=cex.pt*0.25, st.col=NULL, add=FALSE, 
                        cols=NULL, lwd=3, lty=1, seas=c(1,2), xlab='', ylab=NULL, 
                        rownames=NULL, colnames=NULL, latlim=c(-90,90),lonlim=c(-180,180),
                        addcontours=F, contvarname='gph500', conttype='data', contcol='black',
                        contlev=levs, addvectors=F, vecnames=NULL, #vectortype='data',
                        veccol='black', veclen=0.03, vecscale=0.3, vecwd=0.75, every_x_vec=4,
                        wcol='black', zonalmean=F, zmvarname='gph500', colorbar=T){
  oldpar <- par(no.readonly=TRUE)
  # plot specific lat range  
  # latpos <- match(which(x$lat>latlim[1]),which(x$lat<latlim[2]))
  for (i in 1:4) {
    if (i==1) {
      pos <- which(x$lat>latlim[1])
    } else if (i==2){
      pos <- which(x$lat<latlim[2])
    } else if (i==3){
      pos <- which(x$lon>lonlim[1])
    } else if (i==4){
      pos <- which(x$lon<lonlim[2])
    } 
    if (length(dim(x$data)) == 3){
      x$data <- x$data[pos,,,drop=F]
    } else {
      x$data <- x$data[pos,,drop=F]
    }
    x$ensmean <- x$ensmean[pos,,drop=F] 
    x$lon <- x$lon[pos]
    x$lat <- x$lat[pos]
    x$names <- x$names[pos]
  }
  #   lonpos <- match(which(x$lon>lonlim[1]),which(x$lon<lonlim[2]))
  #   if (length(dim(x$data)) == 3){
  #     x$data <- x$data[lonpos,,,drop=F]
  #   } else {
  #     x$data <- x$data[lonpos,,drop=F]
  #   }
  #   x$ensmean <- x$ensmean[lonpos,,drop=F] 
  #   x$lon <- x$lon[lonpos]
  #   x$lat <- x$lat[lonpos]
  #   x$names <- x$names[lonpos]
  dat.i <- which(x$names == varname)
  # Jonas Jan. 2014 for new ECHAM format 
  lons <- x$lon[x$names == varname]
  lats <- x$lat[x$names == varname]
  if (length(dat.i) == 1) ti <- 1:ncol(x$data)
  if (type == 'data'){
    if (length(dim(x$data)) == 3){
      if (dim(x$data)[3]==1) {
        plotdata <- array(x$data[dat.i,ti,],c(length(dat.i),1)) 
      } else {
        plotdata <- x$data[dat.i,ti,] # changed back to original 07/2015 
      }
    } else {
      #      plotdata <- t(x$data[dat.i,ti,drop=F]) 
      plotdata <- x$data[dat.i,ti,drop=F]
    }
  } else if (type == 'ensmean') {
    if (is.matrix(x$ensmean)){
      plotdata <- as.matrix(x$ensmean[dat.i,ti])
    } else {
      plotdata <- as.matrix(x$ensmean[dat.i])
    }
  } else {
    stop("Only 'data' and 'ensmean' are implemented as type")
  }
  if (length(dat.i) == 1){
    if (length(seas) > 1) {
      if (is.null(cols)){
        se.col <- rbfun(14)[c(3,12)]
      } else {
        se.col <- cols[seq(along=seas)]
      }
    } else {
      if (is.null(cols)){
        se.col <- rep(1,max(seas))
      } else {
        se.col <- rep(cols[1], max(seas))
      }
    }
    tis <- as.vector(apply(as.matrix(seas), 1, function(x,y) seq(x, length(y$time),2), y=x))
    if (ncol(plotdata) == 1) {
      if (!add) plot(0, type='n', xlim=range(x$time), ylim=range(plotdata[tis]),
                     xlab=xlab, ylab=if (is.null(ylab)) varname else ylab)
      for (se in 1:2) lines(x$time[seq(se,length(plotdata),2)], 
                     plotdata[seq(se,length(plotdata),2)], col=se.col[se], lwd=lwd, lty=lty)
    } else {
      d.range <- apply(plotdata, 1, range)
      if (!add) plot(0, type='n', xlim=range(x$time), ylim=range(d.range[,tis]),
                     xlab='', ylab=if (is.null(ylab)) varname else ylab)
      for (se in seas){
        for (dr in 1:2){
          lines(x$time[seq(se, ncol(d.range), 2)], d.range[dr,seq(se, ncol(d.range), 2)], 
                col=se.col[se], lwd=lwd/2, lty=lty)
        }
      }
    }
  } else {
    nmod <- ncol(plotdata)
    nrow <- ceiling(nmod/4)
    ncol <- ceiling(nmod/nrow)+1
    nn <- nrow*(ncol-1)
    if (!add){
      dd <- lcm(2)
      if (nmod == 1) dd <- lcm(2.8)
      layout(matrix(c(1:nmod, rep(nmod+2, nn-nmod), rep(nmod+1, nrow)),
                    ncol, nrow, byrow=T), height=c(rep(5, ncol-1),dd))
      oma <- rep(0, 4)
      if (!is.null(rownames)) oma[3] <- 3
      if (!is.null(colnames)) oma[2] <- 3
      par(oma=oma)
    }
    par(mar=rep(0.5,4), cex.axis=1.4, cex.lab=1.4)
    
    if (missing(levs)) {
      if (symmetric){
        levs <- pretty(c(plotdata, -plotdata), 15)
      } else {
        levs <- pretty(plotdata, 15)
      }
    }
    if (is.null(cols)){
      if (varname == 'precip'){
        cols <- gbfun(length(levs) - 1)
      } else { 
        cols <- rbfun(length(levs) - 1)
      }
    } else {
      cols <- rep(cols, length.out=length(levs) - 1)
    }
    #    if (!is.null(centercol)) {
    #      if (is.even(length(cols)/2)) {
    #        cols[(length(cols)/2):((length(cols)/2)+1)] <- centercol
    #      } else {
    #        cols[ceiling(length(cols)/2)] <- centercol
    #      }
    #    }  
    for (i in 1:nmod){
      #      plot(0, type='n', xlim=range(x$lon), ylim=range(x$lat),
      #           axes=F, xlab='', ylab='')      
      plot(0, type='n', xlim=lonlim, ylim=latlim,
           axes=F, xlab='', ylab='')
      li <- as.numeric(cut(plotdata[,i], breaks=levs))
      points(x$lon, x$lat, pch=15, col=cols[li], cex=cex.pt)
      if (!is.null(colnames)) {
        if (i <= length(colnames)) axis(3, at=mean(range(x$lon)), label=colnames[i], 
                                        tick=FALSE, cex=1.4)
      }
      if (!is.null(rownames)) { 
        if (i%%length(colnames) == 1) {
          axis(2, at=mean(range(x$lat)), 
               label=rownames[ceiling(i/length(colnames))], tick=FALSE, cex=1.4)
        }
      }
      world <- map(interior=F, plot=F)
      world$x[abs(diff(world$x)) > 180] <- NA
      lines(world,col=wcol)
      #map(add=T, region=c('Caspian','Great Lakes'),col=wcol) 
      if (addcontours) {
        add_contour(x, cvarname=contvarname, ctype=conttype, ccol=contcol, clev=contlev, cont_i=i)
      }
      if (addvectors) {
        if (type == 'ensmean') {
          uw <- x$ensmean[which(x$names==vecnames[1]),i]
          vw <- x$ensmean[which(x$names==vecnames[2]),i]
        } else if (type == 'data') {
          uw <- x$data[which(x$names==vecnames[1]),1,i]
          vw <- x$data[which(x$names==vecnames[2]),1,i]
        }
        uvlon <- x$lon[which(x$names==vecnames[1])]
        uvlat <- x$lat[which(x$names==vecnames[2])]
        uvpos <- !is.na(match(uvlon,unique(uvlon)[seq(1,length(unique(uvlon)),every_x_vec)]))&
          !is.na(match(uvlat,unique(uvlat)[seq(1,length(unique(uvlat)),every_x_vec)]))
        #        plot(0, type='n', xlim=lonlim, ylim=latlim, axes=F, xlab='', ylab='')
        #        points(x$lon, x$lat, pch=15, col=rep("white",length(x$lon)), cex=1)
        quiver(x=uvlon[uvpos], y=uvlat[uvpos], u=uw[uvpos], v=vw[uvpos],
               scale=vecscale,length=veclen,col=veccol,lwd=vecwd,lend=1)
        #        quiver(x=x$lon, y=x$lat, u=uw, v=vw,
        #               scale=vecscale,length=veclen,col=veccol,lwd=vecwd,lend=1)
      }
      if ((!is.null(stations)) & (any(statpanel==i))) {
        statlon <- stations$lon[!is.na(stations$data[,i])]
        statlat <- stations$lat[!is.na(stations$data[,i])]
        if (is.null(st.col)){
          st <- stations$data[,i]
          lli <- as.numeric(cut(st, breaks=levs))
          #          points(statlon, statlat, pch='.', bg=cols[lli], cex=st.cex*0.01)
          #          tst <- stations$data[stations$name=='temp',i]
          #          pst <- stations$data[stations$name=='precip',i]
          #          sst <- stations$data[stations$name=='slp',i]
          try(points(stations$lon[stations$names=='temp2' & !is.na(stations$data[,i])], 
                     stations$lat[stations$names=='temp2' & !is.na(stations$data[,i])], 
                     pch=1, col=rgb(10,0,0,8,maxColorValue=10), cex=st.cex*2))
          try(points(stations$lon[stations$names=='precip' & !is.na(stations$data[,i])], 
                     stations$lat[stations$names=='precip' & !is.na(stations$data[,i])], 
                     pch=3, col=rgb(10,0,10,8,maxColorValue=10), cex=st.cex))
          try(points(stations$lon[stations$names=='slp' & !is.na(stations$data[,i])], 
                     stations$lat[stations$names=='slp' & !is.na(stations$data[,i])], 
                     pch=4, col=rgb(0,0,10,8,maxColorValue=10), cex=st.cex))
          try(points(stations$lon[stations$names=='prox' & !is.na(stations$data[,i])], 
                     stations$lat[stations$names=='prox' & !is.na(stations$data[,i])], 
                     pch=3, col=rgb(0,10,0,8,maxColorValue=10), cex=st.cex))
        } else {
          points(statlon, statlat, pch=1, col=st.col, cex=st.cex)
        }
      }
      if (!is.null(x$p.value)){
        signif <- x$p.value[,ti,i] > siglev
        points(x$lon[signif], x$lat[signif], pch=4, cex=max(2, cex.pt*0.75))
      }
      # if (!is.null(names)) text(min(x$lon), max(x$lat), names[i], adj=c(0.5,1), cex=1.4)
      if (!is.null(names)) text(min(lonlim), max(latlim), names[i], adj=c(0.5,1), cex=1.4)
      box()
      if (zonalmean) {
        zm=aggregate(x$data[which(x$names==zmvarname),1,i],
                     by=list(x$lat[which(x$names==zmvarname)]),mean)
        li <- max(abs(range(zm[,2])))
#        print(li)
        if (is.na(li)) {
          plot(1,1,ty='n',axes=F,xlab="",ylab="")
        } else {
          plot(zm[,2],zm[,1],ty='l',xlim=c((li*(-1)),li),xlab='',xaxt='n',ylab='',yaxt='n')
          abline(v=0,col='grey')  
        }
      }
    }
    if (colorbar) {
      if (units == ''){
        par(mar=c(3,0.5,1.5,0.5), cex.axis=1.4, cex.lab=1.4)    
        plot_colourbar(col=cols, lev=levs, cex=1.4, xlab=units)
      } else {
        par(mar=c(3,8,1.5,8), cex.axis=1.4, cex.lab=1.4)  
        tmp <- list(col=cols, lev=levs, units=units)
        class(tmp) <- 'plotmap'  
        plot_colourbar(tmp, cex=1.4)
      }
    }
    mtext(main, side=3, line=2, outer=TRUE, font=1, cex=1.4)
    if (!add) par(oldpar)
  }
}












select_stations <- function(x, x.i){
  lapply(x, function(x) if (is.matrix(x)) x[x.i,] else x[x.i])
}



# matrix with dimension: number of proxies times vector of model variables at locations
# for each proxy a matrix which =1 at one location
compute_H <- function(stations, echam, threshold=700){
  H <- array(0, c(nrow(stations$data), nrow(echam$data)))
  for (i in seq(stations$lon)){
    dist <- compute_dist(echam$lon, echam$lat, stations$lon[i], stations$lat[i])
    H[i, which.min(dist)] <- if (min(dist) < threshold) 1 else 0
  }
  H
}

compute_H_temp_precip <- function(stations, echam, threshold=700){
  H <- array(0, c(nrow(stations$data), nrow(echam$data)))
  for (i in seq(stations$lon)){
    dist <- compute_dist(echam$lon, echam$lat, stations$lon[i], stations$lat[i])
    if (stations$names[i] == 'temp2') {
      H[i, which.min(dist)] <- if (min(dist) < threshold) 1 else 0
    }
    if (stations$names[i] == 'precip') {
      H[i, which.min(dist)+one_var_dim] <- if (min(dist) < threshold) 1 else 0
    }
  }
  H
}

compute_H_temp_precip_slp_redux <- function(stations, echam, threshold=700){
  H <- array(0, c((nrow(stations$data)/6), nrow(echam$data)))
  nprox <- length(stations$lon)/6
  for (i in seq(nprox)){
    dist <- compute_dist(echam$lon, echam$lat, stations$lon[i], stations$lat[i])
    if (stations$names[i] == 'temp2') {
      H[i, which.min(dist)] <- if (min(dist) < threshold) 1/6 else 0
      H[i, (which.min(dist)+((dim(echam$data)[1]/6)))] <- if 
        (min(dist) < threshold) 1/6 else 0
      H[i, (which.min(dist)+(2*(dim(echam$data)[1]/6)))] <- if 
        (min(dist) < threshold) 1/6 else 0
      H[i, (which.min(dist)+(3*(dim(echam$data)[1]/6)))] <- if 
        (min(dist) < threshold) 1/6 else 0
      H[i, (which.min(dist)+(4*(dim(echam$data)[1]/6)))] <- if 
        (min(dist) < threshold) 1/6 else 0
      H[i, (which.min(dist)+(5*(dim(echam$data)[1]/6)))] <- if 
        (min(dist) < threshold) 1/6 else 0
    }
    if (stations$names[i] == 'precip') {
      H[i, (which.min(dist)+
        (length(which(echam$names=='temp2'))/6))] <- if 
        (min(dist) < threshold) 1/6 else 0
      H[i, (which.min(dist)+(1*(dim(echam$data)[1]/6)+
        (length(which(echam$names=='temp2'))/6)))] <- if 
        (min(dist) < threshold) 1/6 else 0
      H[i, (which.min(dist)+(2*(dim(echam$data)[1]/6)+
        (length(which(echam$names=='temp2'))/6)))] <- if 
        (min(dist) < threshold) 1/6 else 0
      H[i, (which.min(dist)+(3*(dim(echam$data)[1]/6)+
        (length(which(echam$names=='temp2'))/6)))] <- if 
        (min(dist) < threshold) 1/6 else 0
      H[i, (which.min(dist)+(4*(dim(echam$data)[1]/6)+
        (length(which(echam$names=='temp2'))/6)))] <- if 
        (min(dist) < threshold) 1/6 else 0
      H[i, (which.min(dist)+(5*(dim(echam$data)[1]/6)+
        (length(which(echam$names=='temp2'))/6)))] <- if 
        (min(dist) < threshold) 1/6 else 0
    }
    if (stations$names[i] == 'slp') {
      H[i, (which.min(dist)+
        (2*length(which(echam$names=='temp2'))/6))] <- if 
        (min(dist) < threshold) 1/6 else 0
      H[i, (which.min(dist)+(1*(dim(echam$data)[1]/6)+
        (2*length(which(echam$names=='temp2'))/6)))] <- if 
        (min(dist) < threshold) 1/6 else 0
      H[i, (which.min(dist)+(2*(dim(echam$data)[1]/6)+
        (2*length(which(echam$names=='temp2'))/6)))] <- if 
        (min(dist) < threshold) 1/6 else 0
      H[i, (which.min(dist)+(3*(dim(echam$data)[1]/6)+
        (2*length(which(echam$names=='temp2'))/6)))] <- if 
        (min(dist) < threshold) 1/6 else 0
      H[i, (which.min(dist)+(4*(dim(echam$data)[1]/6)+
        (2*length(which(echam$names=='temp2'))/6)))] <- if 
        (min(dist) < threshold) 1/6 else 0
      H[i, (which.min(dist)+(5*(dim(echam$data)[1]/6)+
        (length(which(echam$names=='temp2'))/6)))] <- if 
        (min(dist) < threshold) 1/6 else 0
    }
  }
  H
}

compute_H_temp_precip_slp <- function(stations, echam, threshold=700){
  H <- array(0, c(nrow(stations$data), nrow(echam$data)))
  for (i in seq(stations$lon)){
    dist <- compute_dist(echam$lon, echam$lat, stations$lon[i], stations$lat[i])
    if (stations$names[i] == 'temp2') {
      H[i, which.min(dist)] <- if (min(dist) < threshold) 1 else 0
    }
    if (stations$names[i] == 'precip') {
      H[i, which.min(dist)+one_var_dim] <- if (min(dist) < threshold) 1 else 0
    }
    if (stations$names[i] == 'slp') {
      H[i, which.min(dist)+(2*one_var_dim)] <- if (min(dist) < threshold) 1 else 0
    }
  }
  H
}

compute_H_bias_correct <- function(stations, echam, echamatprox, seas=1, threshold=700){       
  stat.seas=array(stations$data,c(dim(stations$data)[1],seas,(dim(stations$data)[2])/seas))
  # stat.seas[stat.seas < -100] = NA
  # stat.seas[stat.seas > 100] = NA 
  echam.seas=array(echamatprox,c(dim(echamatprox)[1],seas,(dim(echamatprox)[2])/seas))
  diff.mn=(echam.seas - stat.seas)
#  diff.mean=apply(diff.mn,1:2,mean,na.rm=T)
  diff.mean1=apply(diff.mn[,1,],1,mean,na.rm=T)
  diff.mean2=apply(diff.mn[,2,],1,mean,na.rm=T)
  diff.mean=apply(cbind(diff.mean1,diff.mean2),1,mean,na.rm=T)
  H <- array(0, c(nrow(stations$data), nrow(echam$data)))
  for (i in seq(stations$lon)){
    dist <- compute_dist(echam$lon, echam$lat, stations$lon[i], stations$lat[i])
    H[i, which.min(dist)] <- if (min(dist) < threshold) diff.mean[i] else 0
  }  
  H
# H has to have dim 380x30 (# of stations x ensemble size)
}

compute_dist_2d <- function(lon, lat, lon0, lat0){
  lo <- lon/180*pi
  la <- lat/180*pi
  lo0 <- lon0/180*pi
  la0 <- lat0/180*pi
  tmp <- outer(sin(la), sin(la0), '*') + outer(cos(la), cos(la0), '*') * outer(lo, lo0, function(x,y) cos(x - y))
  ##tmp <- sin(la)*sin(la0) + cos(la) * cos(la0) * cos(lo - lo0)
  # deal with rounding errors
  tmp[abs(tmp) >= 1] <- round(tmp[abs(tmp) >= 1])
  out <- 6375*acos(tmp)
  out
}

compute_dist <- function(lon, lat, lon0, lat0){
  lo <- lon/180*pi
  la <- lat/180*pi
  lo0 <- lon0/180*pi
  la0 <- lat0/180*pi
  tmp <- sin(la)*sin(la0) + cos(la) * cos(la0) * cos(lo - lo0)
  # deal with rounding errors
  tmp[abs(tmp) >= 1] <- round(tmp[abs(tmp) >= 1])
  out <- 6375*acos(tmp)
  out
}

compute_giorgi_H_sixmon <- function(giorgi, echam) { #}, numvar){
  nll <- length(echam$lon)/(length(unique(echam$names))*6)
  giorgi.arr <- t(sapply(giorgi, function(x, y){
  lons <- y$lon[1:nll]
  lats <- y$lat[1:nll]
  ind <- c(lons > x$edges[1] & lons < x$edges[2] & lats > x$edges[3] & lats < x$edges[4])
  weights <- ind * cos(lats/180*pi)
  weights <- weights/sum(weights)
  weights}, y=echam))
  # stitch together to form the full giorgi_H
  nreg <- nrow(giorgi.arr)
  #  ndim <- one_var_dim
  ndim <- length(which(echam$names=='temp2'))
  nvar <- floor(nrow(echam$data)/ndim)
#  nvar <- numvar
#  H.giorgi <- array(0, c(nreg*nvar, nrow(echam$data)))
  H.giorgi <- array(0, c(nreg, nrow(echam$data)))
  for (i in 1:nvar) {
#    H.giorgi[(i-1)*nreg + 1:nreg,(i-1)*ndim + 1:ndim] <- giorgi.arr
    H.giorgi[,(i-1)*ndim + 1:ndim] <- giorgi.arr
  }
  return(H.giorgi)
}


compute_giorgi_H_numvar <- function(giorgi, echam, numvar){
  giorgi.arr <- t(sapply(giorgi, function(x, y){
    lons <- y$lon
    lats <- y$lat
#    print(cbind(lons,lats))
    ind <- c(lons > x$edges[1] & lons < x$edges[2] & lats > x$edges[3] & lats < x$edges[4])
    weights <- ind * cos(lats/180*pi)
    weights <- weights/sum(weights)
    weights}, y=echam))
  # stitch together to form the full giorgi_H
  nreg <- nrow(giorgi.arr)
  ndim <- one_var_dim
  #  nvar <- floor(nrow(echam$data)/ndim)
  nvar <- numvar
  H.giorgi <- array(0, c(nreg*nvar, nrow(echam$data)))
  for (i in 1:nvar) {
    H.giorgi[(i-1)*nreg + 1:nreg,(i-1)*ndim + 1:ndim] <- giorgi.arr
  }
  return(H.giorgi)
}


compute_giorgi_H <- function(giorgi, echam){
  giorgi.arr <- t(sapply(giorgi, function(x, y){
    lons <- y$lon
    lats <- y$lat
    ind <- c(lons > x$edges[1] & lons < x$edges[2] & lats > x$edges[3] & lats < x$edges[4])
    weights <- ind * cos(lats/180*pi)
    weights <- weights/sum(weights)
    weights}, y=echam))
  # stitch together to form the full giorgi_H
  nreg <- nrow(giorgi.arr)
  ndim <- length(echam$lon)
  nvar <- floor(nrow(echam$data)/ndim)
  H.giorgi <- array(0, c(nreg*nvar, nrow(echam$data)))
  for (i in 1:nvar) {
    H.giorgi[(i-1)*nreg + 1:nreg,(i-1)*ndim + 1:ndim] <- giorgi.arr
  }
  return(H.giorgi)
}

compute_giorgi_H_v2 <- function(giorgi, echam){
  nll <- length(echam$lon)/(length(unique(echam$names)))
  giorgi.arr <- t(sapply(giorgi, function(x, y){
    lons <- y$lon[1:nll]
    lats <- y$lat[1:nll]
    ind <- c(lons > x$edges[1] & lons < x$edges[2] & lats > x$edges[3] & lats < x$edges[4])
    weights <- ind * cos(lats/180*pi)
    weights <- weights/sum(weights)
    weights}, y=echam))
  # stitch together to form the full giorgi_H
  nreg <- nrow(giorgi.arr)
  ndim <- length(echam$lon)
  nvar <- floor(nrow(echam$data)/ndim)
  H.giorgi <- array(0, c(nreg*nvar, nrow(echam$data)))
  for (i in 1:nvar) {
    H.giorgi[(i-1)*nreg + 1:nreg,(i-1)*ndim + 1:ndim] <- giorgi.arr
  }
  return(H.giorgi)
}


height_correction <- function(stations, echam){
  H <- compute_H(stations, echam)
  dT <- 0.0065 * (H %*% echam$height - stations$height)
  dT[is.na(dT)] <- 0
  stations$data <- stations$data - as.vector(dT)
  invisible(stations)
}

reorganize <- function(x, ind, seas=6){
  dims <- dim(x)
  if (length(dims) == 0) dims <- length(x)
  if (any(dims == length(ind))){
    dim.i <- max(which(dims == length(ind)))
    if (length(dims) > 1){
      xx <- array(aperm(x, c(dim.i, setdiff(seq(along=dims), dim.i))), c(length(ind), length(x)/length(ind)))
    } else {
      xx <- as.matrix(x)
    }
    xout <- array(xx[ind,], c(seas, 12/seas,length(ind)/12,length(xx)/length(ind)))
    xout <- array(aperm(xout, c(2,1,3,4)), c(12/seas, length(ind)/12*seas, dims[setdiff(seq(along=dims), dim.i)]))
    if (length(dims) > 1){
      perm.i <- c(dim.i + 1:(dim.i-1), 1:dim.i, (2*dim.i):length(dim(xout)))
      perm.i <- perm.i[1:length(dim(xout))]
      xout <- aperm(xout, perm.i)
    } else {
      xout <- as.matrix(xout)
    }
    return(xout)
  } else {
    return(x)
  }
}

corr_function <- function(dist, L){
  exp(-1/2 * (dist/L)**2)
}

rmse_fun <- function(x, y, seas=2, center=F){
  if (!is.null(x$data)) {
    anadat <- array(x$data, c(nrow(x$data), seas, ncol(x$data)/seas, length(x$data)/prod(dim(x$data)[1:2])))
    obsdat <- array(y$data, c(nrow(y$data), seas, ncol(y$data)/seas, length(y$data)/prod(dim(y$data)[1:2])))
    if (center){
      anadat <- aperm(apply(anadat, c(1,2,4), function(x) x - mean(x, na.rm=T)), c(2,3,1,4))
      obsdat <- aperm(apply(obsdat, c(1,2,4), function(x) x - mean(x, na.rm=T)), c(2,3,1,4))
    }
    x$data <- sqrt(apply((anadat - as.vector(obsdat))**2, c(1,2,4), mean, na.rm=T))[,,1:dim(anadat)[4]]
  }
  if (!is.null(x$ensmean)) {
    anadat <- array(x$ensmean, c(nrow(x$ensmean), seas, ncol(x$ensmean)/seas))
    obsdat <- array(y$data, dim(anadat))
    if (center){
      anadat <- aperm(apply(anadat, 1:2, function(x) x-mean(x,na.rm=T)), c(2,3,1))
      obsdat <- aperm(apply(obsdat, 1:2, function(x) x-mean(x,na.rm=T)), c(2,3,1))
    }
    x$ensmean <- sqrt(apply((anadat - as.vector(obsdat))**2, 1:2, mean,na.rm=T))
  }
  return(x)
}


# average bias over entire time period of data
bias_fun <- function(x, y, seas=2){
  if (!is.null(x$data)) {
    anadat <- array(x$data, c(nrow(x$data), seas, ncol(x$data)/seas, length(x$data)/prod(dim(x$data)[1:2])))
    obsdat <- array(y$data, c(nrow(y$data), seas, ncol(y$data)/seas, length(y$data)/prod(dim(y$data)[1:2])))
    x$data <- apply(anadat - as.vector(obsdat), c(1,2,4), mean, na.rm=T)[,,1:dim(anadat)[4], drop=F]
  }
  if (!is.null(x$ensmean)) {
    anadat <- array(x$ensmean, c(nrow(x$ensmean), seas, ncol(x$ensmean)/seas))
    if (is.null(y$ensmean)) {
      obsdat <- array(y$data, c(nrow(y$data), seas, ncol(y$data)/seas))
    } else {
      obsdat <- array(y$ensmean, c(nrow(y$ensmean), seas, ncol(y$ensmean)/seas))
    }
    x$ensmean <- apply(anadat - as.vector(obsdat), 1:2, mean,na.rm=T)
  }
  return(x)
}


# difference at each time step, NO time average. For observation error calculation
diff_fun <- function(x, y, seas=2){
  if (!is.null(x$data)) {
    anadat <- array(x$data, c(nrow(x$data), seas, ncol(x$data)/seas, length(x$data)/prod(dim(x$data)[1:2])))
    obsdat <- array(y$data, c(nrow(y$data), seas, ncol(y$data)/seas, length(y$data)/prod(dim(y$data)[1:2])))
    #x$data <- apply(anadat - as.vector(obsdat), c(1,2,4), mean, na.rm=T)[,,1:dim(anadat)[4], drop=F]
    x$data <- anadat - as.vector(obsdat)
  }
  if (!is.null(x$ensmean)) {
    anadat <- array(x$ensmean, c(nrow(x$ensmean), seas, ncol(x$ensmean)/seas))
    if (is.null(y$ensmean)) {
      obsdat <- array(y$data, c(nrow(y$data), seas, ncol(y$data)/seas))
    } else {
      obsdat <- array(y$ensmean, c(nrow(y$ensmean), seas, ncol(y$ensmean)/seas))
    }
    #x$ensmean <- apply(anadat - as.vector(obsdat), 1:2, mean,na.rm=T)
    x$ensmean <- anadat - as.vector(obsdat)
  }
  return(x)
}


rank_fun <- function(x, y, seas=2){
  anadat <- array(x$data, c(nrow(x$data), seas, ncol(x$data)/seas, length(x$data)/prod(dim(x$data)[1:2])))
  obsdat <- array(y$data, c(nrow(y$data), seas, ncol(y$data)/seas, length(y$data)/prod(dim(y$data)[1:2])))
  x$data <- NULL
  x$ensmean <- apply(anadat > as.vector(obsdat), c(1,2,4), sum)[,,1:dim(anadat)[4], drop=F]
  return(x)
}

RE_fun <- function(x, y){
  if (!is.null(x$data)) {
    if (!is.null(y$data)){
      x$data <- 1 - (x$data/as.vector(y$data))**2
    } else {
      x$data <- 1 - (x$data/as.vector(y$ensmean))**2
    }
  }
  if (!is.null(x$ensmean)){
    if (!is.null(y$ensmean)){
      x$ensmean <- 1- (x$ensmean/as.vector(y$ensmean))**2
    } else {
      x$ensmean <- 1- (x$ensmean/as.vector(y$data))**2
    }
  }
  return(x)
}

corr_fun <- function(x, y, seas=2){
  if (!is.null(x$data)) {
    xdat <- array(x$data, c(nrow(x$data), seas, ncol(x$data)/seas, length(x$data)/prod(dim(x$data)[1:2])))
    ydat <- array(y$data, c(nrow(x$data), seas, ncol(x$data)/seas, 1))[,,,rep(1, dim(xdat)[4]), drop=F]
    perm.i <- c(3,setdiff(seq(along=dim(xdat)), 3))
    xxdat <- array(aperm(xdat, perm.i), c(dim(xdat)[3], length(xdat)/dim(xdat)[3]))
    yydat <- array(aperm(ydat, perm.i), c(dim(ydat)[3], length(ydat)/dim(ydat)[3]))
    outcor <- apply(as.matrix(1:ncol(xxdat)), 1, function(x) cor(xxdat[,x], yydat[,x],use="pairwise.complete.obs"))
    outdim <- dim(xdat)[perm.i[2:length(perm.i)]]
    x$data <- array(outcor, outdim[outdim > 1])
  }
  if (!is.null(x$ensmean)) {
    xdat <- array(x$ensmean, c(nrow(x$ensmean), seas, ncol(x$ensmean)/seas))
    ydat <- array(y$data, c(nrow(x$ensmean), seas, ncol(x$ensmean)/seas))
    perm.i <- c(3,setdiff(seq(along=dim(xdat)), 3))
    xxdat <- array(aperm(xdat, perm.i), c(dim(xdat)[3], length(xdat)/dim(xdat)[3]))
    yydat <- array(aperm(ydat, perm.i), c(dim(ydat)[3], length(ydat)/dim(ydat)[3]))
    outcor <- apply(as.matrix(1:ncol(xxdat)), 1, function(x) cor(xxdat[,x], yydat[,x],use="pairwise.complete.obs"))
    outdim <- dim(xdat)[perm.i[2:length(perm.i)]]
    x$ensmean <- array(outcor, outdim[outdim > 1])
  }
  return(x)
}

compute_mean <- function(x, H.mn) {
  xdat <- array(NA, c(nrow(H.mn), dim(x$data)[2:3]))
  for (i in 1:dim(x$data)[3]) xdat[,,i] <- H.mn %*% x$data[,,i]
  ensdat <- H.mn %*% x$ensmean
  x$data <- xdat
  x$ensmean <- ensdat
  return(x)}

compute_avg_RE_pseudoproxy <- function(H, ech,ana,val){
    valmn <- val
    valmn$data <- H %*% val$data
    ech.mn <- compute_mean(ech, H)
    ana.mn <- compute_mean(ana, H)
    ech.rmse <- rmse_fun(ech.mn, valmn, seas=2)
    ana.rmse <- rmse_fun(ana.mn, valmn, seas=2)
    RE <- RE_fun(ana.rmse, ech.rmse)
    if (no_forc_big_ens) {
      n<-n_no_forc
    } else { 
      n=nmem
    }
    RE.mat <- array(0, c(nrow(H), 2, n))
    RE.mat[,,1:(n-1)] <- RE$data
    RE.mat[,,n] <- RE$ensmean
    return(RE.mat)
}

compute_avg_RE <- function(H, ech,ana,val){
    valmn <- val
##    tpsdim <- length(c(which(val$name=='temp2'),which(val$name=='precip'),which(val$name=='slp')))
##    Hvalsh <- H[,1:tpsdim]
#    valmn$data[is.na(val$data)]=9e36
    H <- H[,1:(length(which(ech$names=="temp2")))] #*3)]
##   valmn$data <- H %*% val$data
# land sea mask to take care of missing validata in the ocean
    vpos <- match(which(!is.na(valmn$data[,1])), (which(val$names=="temp2")))
    vpos <- vpos[!is.na(vpos)]
    H <- H[,vpos]
    valmn$data <- valmn$data[vpos,]
    valmn$data <- H %*% valmn$data #[1:(length(which(ech$names=="temp2"))*3),]
    ech$data <- ech$data[vpos,,] #[1:(length(which(ech$names=="temp2"))*3),,]
    ech$ensmean <- ech$ensmean[vpos,] #[1:(length(which(ech$names=="temp2"))*3),]
    ana$data <- ana$data[vpos,,] #[1:(length(which(ana$names=="temp2"))*3),,]
    ana$ensmean <- ana$ensmean[vpos,] #[1:(length(which(ana$names=="temp2"))*3),]
    ech.mn <- compute_mean(ech, H)
    ana.mn <- compute_mean(ana, H)
    ech.rmse <- rmse_fun(ech.mn, valmn, seas=2)
    ana.rmse <- rmse_fun(ana.mn, valmn, seas=2)
    RE <- RE_fun(ana.rmse, ech.rmse)
    if (no_forc_big_ens) {
      n<-n_no_forc
    } else { 
      n=nmem
    }
    RE.mat <- array(0, c(nrow(H), 2, (n+1)))
    RE.mat[,,1:n] <- RE$data
    RE.mat[,,(n+1)] <- RE$ensmean
    return(RE.mat)
}


compute_avg_corr_pseudoprox <- function(H, ech,ana,val){
    valmn <- val
    valmn$data <- H %*% val$data
    ech.mn <- compute_mean(ech, H)
    ana.mn <- compute_mean(ana, H)
    ech.corr <- corr_fun(ech.mn, valmn, seas=2)
    ana.corr <- corr_fun(ana.mn, valmn, seas=2)
    corr.mat <- array(0, c(nrow(H)*2, 2, 30))
    corr.mat[seq(1,nrow(H)*2, 2),,1:29] <- ech.corr$data
    corr.mat[seq(1,nrow(H)*2, 2),,30] <- ech.corr$ensmean
    corr.mat[seq(2,nrow(H)*2, 2),,1:29] <- ana.corr$data
    corr.mat[seq(2,nrow(H)*2, 2),,30] <- ana.corr$ensmean
    return(corr.mat)
}

compute_avg_corr <- function(H, ech,ana,val){
    valmn <- val
     H <- H[,1:length(which(ech$names=="temp2"))*3]
     val$data <- val$data[1:length(which(ech$names=="temp2"))*3,]
    val$data[is.na(val$data)] = 0
    valmn$data <- H %*% val$data
     ech$data <- ech$data[1:length(which(ech$names=="temp2"))*3,,]
     ech$ensmean <- ech$ensmean[1:length(which(ech$names=="temp2"))*3,]
     ana$data <- ana$data[1:length(which(ana$names=="temp2"))*3,,]
     ana$ensmean <- ana$ensmean[1:length(which(ana$names=="temp2"))*3,]
    ech.mn <- compute_mean(ech, H)
    ana.mn <- compute_mean(ana, H)
    ech.corr <- corr_fun(ech.mn, valmn, seas=2)
    ana.corr <- corr_fun(ana.mn, valmn, seas=2)
    corr.mat <- array(0, c(nrow(H)*2, 2, 31))
    corr.mat[seq(1,nrow(H)*2, 2),,1:30] <- ech.corr$data
    corr.mat[seq(1,nrow(H)*2, 2),,31] <- ech.corr$ensmean
    corr.mat[seq(2,nrow(H)*2, 2),,1:30] <- ana.corr$data
    corr.mat[seq(2,nrow(H)*2, 2),,31] <- ana.corr$ensmean
    return(corr.mat)
}


# dyn.load('EnSRF.so')
# 
# EnSRF <- function(echam, calibrate, Hcal=NULL, weights=NULL, R=NULL) {
#   
#   x <- echam$data
#   y <- calibrate$data
#   if (is.null(Hcal)) Hcal <- compute_H(calibrate, echam)
#   ndim <- nrow(x)
#   nens <- dim(x)[3]
#   ntim <- ncol(x)
#   nobs <- nrow(y)
#   # JF May 2013
# #   nh <- length(Hcal[1,Hcal[1,]!=0,drop=FALSE])
# #   if (nh < 2) {
# #     print("WARNING: < 2 values in short H operator")
# #     nh=2
# #   }
# #  print(nh)
#   #
#   if (is.null(weights)) weights <- 1
#   if (length(weights) != ndim**2){
#     print('WARNING: weights matrix has wrong dimensions, set to 1')
#     weights <- array(1, c(ndim, ndim))
#   }
#   if (is.null(R)) R <- rep(0, nobs)
#   R <- array(R, c(nobs, ntim))
#   
#   tmp <- .Fortran('EnSRF', 
#                   x=as.double(x),
#                   y=as.double(y),
#                   Hcal=as.double(Hcal),
#                   R=as.double(R),
#                   weights=as.double(weights),
#                   ndim=as.integer(ndim),
#                   ntim=as.integer(ntim),
#                   nens=as.integer(nens),
#                   nobs=as.integer(nobs),
#                   xout=as.double(x),
#                   xens=as.double(x[,,1]))
# 
#   out <- echam
#   out$data <- array(tmp$xout, dim(echam$data))
#   out$ensmean <- array(tmp$xens, dim(echam$data)[1:2])
#   invisible(out)
#   
# }



# dyn.load('EnSRF_new.so')
# 
# EnSRF_new <- function(echam, calibrate, Hcal=NULL, weights=NULL, R=NULL) {
#   
#   x <- echam$data
#   y <- calibrate$data
#   if (is.null(Hcal)) Hcal <- compute_H(calibrate, echam)
#   ndim <- nrow(x)
#   nens <- dim(x)[3]
#   ntim <- ncol(x)
#   nobs <- nrow(y)
#   # JF May 2013
#     nh <- length(Hcal[1,Hcal[1,]!=0,drop=FALSE])
#     if ((nh < 2) & (addbias==T)) {
#       print("WARNING: < 2 values in short H operator")
#       nh=2
#     }
# #   print(nh)
#   
#   if (is.null(weights)) {
#     print('ATTENTION: distance weigths are zero')
#     weights <- 1
#   }
#   if (length(weights) != ndim**2){
#     print('ATTENTION: distance weigths have wrong dimension')
#     weights <- array(1, c(ndim, ndim))
#   }
#   if (is.null(R)) {
#     print('ATTENTION: R zero')
#     R <- rep(0, nobs)
#   }
#   R <- array(R, c(nobs, ntim))
#   
#   tmp <- .Fortran('EnSRF_new', 
#                   x=as.double(x),
#                   y=as.double(y),
#                   Hcal=as.double(Hcal),
#                   R=as.double(R),
#                   weights=as.double(weights),
#                   ndim=as.integer(ndim),
#                   ntim=as.integer(ntim),
#                   nens=as.integer(nens),
#                   nobs=as.integer(nobs),
#                   nh=as.integer(nh),
#                   xout=as.double(x),
#                   xens=as.double(x[,,1]))
#   
#   out <- echam
#   out$data <- array(tmp$xout, dim(echam$data))
#   out$ensmean <- array(tmp$xens, dim(echam$data)[1:2])
#   invisible(out)
#   
# }



EnSRF2 <- function(echam, calibrate, R=NULL, Hcal=NULL, weights=1){

   if (is.list(echam)){
     x <- echam$data
     if (is.null(echam$ensmean)) {
       xmn <- apply(x, 1:2, mean, na.rm=T)
     } else {
       xmn <- echam$ensmean
     }
   } else {
     x <- echam
     xmn <- apply(echam, 1:2, mean, na.rm=T)
   }
   if (is.list(calibrate)){
     y <- calibrate$data
   } else {
     if (is.null(Hcal)) stop('Cannot figure out forward operator')
   }
   if (is.null(Hcal)){
     if (is.list(echam) & is.list(calibrate)) {
       Hcal <- compute_H(calibrate, echam)  
     } else {
       stop('Cannot figure out forward operator')
     }
   }

   nens <- dim(x)[3]
   ntim <- dim(x)[2]
   ndim <- dim(x)[1]
   nobs <- nrow(Hcal)

   xout <- array(NA, dim(x))
   xoutmn <- array(NA, dim(xmn))
   if (is.null(R)) R <- rep(0, nobs)
   R <- array(R, c(nobs, ntim))
   print('Processing time step')
   for (i in 1:ntim) {
     print(i)
     xb <- xmn[,i]
     xb1 <- x[,i,] - xb
     # Covariance maxtrix P for 1 time step, full state vector and 30 members 
     P <- xb1 %*% t(xb1) / (nens - 1) * weights
     for (j in which(!is.na(y[,i]))) {
       H <- Hcal[j,, drop=F]
       HPHR <-  as.vector(H %*% P %*% t(H) + R[j,i])
       K <- P %*% t(H) / HPHR
       Ktilde <- K / (1 + sqrt(R[j,i]/HPHR))
       xa <- xb + K %*% (y[j,i] - H %*% xb)
       # bias correction in H with state vector extension of extra 1 value at end
       # now H at this position in vector has to include bias
       # e.g. +3 if model is 3 degree colder than observations!
       xa1 <- xb1 - Ktilde %*% H %*% xb1
       # redefine xb and xb1 for next iteration
       xb <- xa
       xb1 <- xa1
     }
     xout[,i,] <- xb1 + as.vector(xb)
     xoutmn[,i] <- xb   
   }
   if (is.list(echam)){
     out <- echam
     out$data <- xout
     out$ensmean <- xoutmn
   } else {
     out <- list(data=xout, ensmean=xoutmn)
   }

   invisible(out)
 }




EnSRF3 <- function(echam, calibrate, R=NULL, Hcal=NULL, Hadd=NULL, weights=1){

   if (is.list(echam)){
     x <- echam$data
     if (is.null(echam$ensmean)) {
       xmn <- apply(x, 1:2, mean, na.rm=T)
     } else {
       xmn <- echam$ensmean
     }
   } else {
     x <- echam
     xmn <- apply(echam, 1:2, mean, na.rm=T)
   }
   if (is.list(calibrate)){
     y <- calibrate$data
   } else {
     if (is.null(Hcal)) stop('Cannot figure out forward operator')
   }
   if (is.null(Hcal)){
     if (is.list(echam) & is.list(calibrate)) {
       Hcal <- compute_H(calibrate, echam)  
     } else {
       stop('Cannot figure out forward operator')
     }
   }
   if (is.null(Hadd)){
     Hadd <- array(0, dim(Hcal)[1])
   }

   nens <- dim(x)[3]
   ntim <- dim(x)[2]
   ndim <- dim(x)[1]
   nobs <- nrow(Hcal)

   xout <- array(NA, dim(x))
   xoutmn <- array(NA, dim(xmn))
   if (is.null(R)) R <- rep(0, nobs)
   R <- array(R, c(nobs, ntim))
   print('Processing time step')
   for (i in 1:ntim) {
     print(i)
     xb <- xmn[,i]
     xb1 <- x[,i,] - xb
     #P <- xb1H %*% xb11) / (nens - 1) * weights
     for (j in which(!is.na(y[,i]))) {
       H <- Hcal[j,, drop=F]
       Hbx <- H %*% xb1 + Hadd[j,, drop=F]
       #Hbx <- H %*% xb1 + dim(array(rep(Hadd[j, drop=F],30),c(30,1)))
       ##HPHR <-  as.vector(H %*% P %*% t(H) + R[j,i])
       HPHR <- as.vector(Hbx %*% t(Hbx) + R[j,i])
       ##K <- P %*% t(H) / HPHR
       K <- xb1 %*% t(Hbx) / (nens - 1) / HPHR 
       Ktilde <- K / (1 + sqrt(R[j,i]/HPHR))
       ##xa <- xb + K %*% (y[j,i] - H %*% xb)
       ##xa1 <- xb1 - Ktilde %*% H %*% xb1
       #xa <- xb + K %*% (y[j,i] - H %*% xb - Hadd[j,, drop=F])
       xa <- xb + K %*% (y[j,i] - Hbx)
       xa1 <- xb1 - Ktilde %*% Hbx
       # redefine xb and xb1 for next iteration
       xb <- xa
       xb1 <- xa1
     }
     xout[,i,] <- xb1 + as.vector(xb)
     xoutmn[,i] <- xb   
   }
   if (is.list(echam)){
     out <- echam
     out$data <- xout
     out$ensmean <- xoutmn
   } else {
     out <- list(data=xout, ensmean=xoutmn)
   }

   invisible(out)
 }



# smaller covariance matrix to reduce memory
# Fortran Version compilation: R CMD SHLIB EnSRF.f95
EnSRF4 <- function(echam, calibrate, R=NULL, Hcal=NULL, weights=1){
  
  if (is.list(echam)){
    x <- echam$data
    if (is.null(echam$ensmean)) {
      xmn <- apply(x, 1:2, mean, na.rm=T)
    } else {
      xmn <- echam$ensmean
    }
  } else {
    x <- echam
    xmn <- apply(echam, 1:2, mean, na.rm=T)
  }
  if (is.list(calibrate)){
    y <- calibrate$data
  } else {
    if (is.null(Hcal)) stop('Cannot figure out forward operator')
  }
  if (is.null(Hcal)){
    if (is.list(echam) & is.list(calibrate)) {
      Hcal <- compute_H(calibrate, echam)  
    } else {
      stop('Cannot figure out forward operator')
    }
  }
  
  nens <- dim(x)[3]
  ntim <- dim(x)[2]
  ndim <- dim(x)[1]
  nobs <- nrow(Hcal)
  
  xout <- array(NA, dim(x))
  xoutmn <- array(NA, dim(xmn))
  if (is.null(R)) R <- rep(0, nobs)
  R <- array(R, c(nobs, ntim))
  print('Processing time step')
  for (i in 1:ntim) {
    print(i)
    xb <- xmn[,i]
    xb1 <- x[,i,] - xb
    ## Covariance maxtrix P for 1 time step, full state vector and 30 members 
    #P <- xb1 %*% t(xb1) / (nens - 1) * weights
    # NEW J. Franke 05/2013: xb NO anomalies
#   xb2 <- x[,i,]
    for (j in which(!is.na(y[,i]))) {
      print(j)
      H <- Hcal[j,, drop=F]
      # NEW J. Franke 05/2013: P in loop
      # P and P2 should be identical but there are differences after 3. decimal
      # position. NO clue what COV does (NOT nens-1). 
      # We use correct matrix calculation!
      #P2 <- t(cov(t(xb2[H!=0,]),t(xb2)))*weights[,H!=0]
      P <- xb1 %*% t(xb1[H!=0,,drop=FALSE]) / (nens - 1) * weights[,H!=0,drop=FALSE]
      #xb4 <- xb3[H!=0,]
      #P2 <- xb4 %*% t(xb3) / (nens - 1) * weights
      H1 <- H[1,H!=0,drop=FALSE] 
      HPHR <-  as.vector(H %*% P %*% t(H1) + R[j,i])
      K <- P %*% t(H1) / HPHR
      Ktilde <- K / (1 + sqrt(R[j,i]/HPHR))
      xa <- xb + K %*% (y[j,i] - H %*% xb)
      xa1 <- xb1 - Ktilde %*% H %*% xb1
      # working version with small cov matrix
      #  P2 <- t(cov(t(xb2[H!=0,]),t(xb2)))
      #  H2 <- H[1,H!=0,drop=FALSE] 
      #  HPHR2 <-  as.vector(H %*% P2 %*% t(H2) + R[j,i])
      #  K2 <- P2 %*% t(H2) / HPHR2
      #  Ktilde2 <- K / (1 + sqrt(R[j,i]/HPHR))
      #  xa2 <- xb + K %*% (y[j,i] - H %*% xb)
      #  xa12 <- xb1 - Ktilde %*% H %*% xb1
      # original version from Jonas (P calculated before observation loop)
      #HPHR <-  as.vector(H %*% P %*% t(H) + R[j,i])
      #K <- P %*% t(H) / HPHR
      #Ktilde <- K / (1 + sqrt(R[j,i]/HPHR))
      #xa <- xb + K %*% (y[j,i] - H %*% xb)
      ## bias correction in H with state vector extension of extra 1 value at end
      ## now H at this position in vector has to include bias
      ## e.g. +3 if model is 3 degree colder than observations!
      #xa1 <- xb1 - Ktilde %*% H %*% xb1
      ## redefine xb and xb1 for next iteration
      xb <- xa
      xb1 <- xa1
    }
    xout[,i,] <- xb1 + as.vector(xb)
    xoutmn[,i] <- xb   
  }
  if (is.list(echam)){
    out <- echam
    out$data <- xout
    out$ensmean <- xoutmn
  } else {
    out <- list(data=xout, ensmean=xoutmn)
  }
  invisible(out)
}



# ------------------------------------------------------------------------------
# compute localization weights for covariance matrix P
# ------------------------------------------------------------------------------
compute_weights <- function(l_dist_slp=l_dist_slp, l_dist_precip=l_dist_precip, 
                            l_dist=l_dist){
#
if (addbias) {
  load(paste("../data/echam/echam_addbias_",syr,"-",eyr,".Rdata",sep=""))
} else {
  load(paste("../data/echam/echam_",syr,"-",eyr,".Rdata",sep=""))
}  
dist.outer <- array(0, c(length(echam$lon),length(echam$lon)))
dist.outer.stream <- array(0, c(length(echam$lonstream),length(echam$lon)))
dist.outer.stream2 <- array(0, c(length(echam$lonstream),length(echam$lonstream)))
dist.outer.stream_gph300 <- array(0, c(length(echam$lonstream),length(echam$longph300)))
dist.outer.stream_gph100 <- array(0, c(length(echam$lonstream),length(echam$longph100)))
dist.outer.stream_gph100_2 <- array(0, c(length(echam$lonstream),length(echam$longph100_2)))
dist.outer.stream_u200 <- array(0, c(length(echam$lonstream),length(echam$lonu200)))
dist.outer.stream_omega <- array(0, c(length(echam$lonstream),length(echam$lonomega)))
dist.outer.stream_omega_2 <- array(0, c(length(echam$lonstream),length(echam$lonomega_2)))
dist.outer.gph300 <- array(0, c(length(echam$longph300),length(echam$lon)))
dist.outer.gph3002 <- array(0, c(length(echam$longph300),length(echam$longph300)))
dist.outer.gph300_gph100 <- array(0, c(length(echam$longph300),length(echam$longph100)))
dist.outer.gph300_gph100_2 <- array(0, c(length(echam$longph300),length(echam$longph100_2)))
dist.outer.gph300_u200 <- array(0, c(length(echam$longph300),length(echam$lonu200)))
dist.outer.gph300_omega <- array(0, c(length(echam$longph300),length(echam$lonomega)))
dist.outer.gph300_omega_2 <- array(0, c(length(echam$longph300),length(echam$lonomega_2)))
dist.outer.gph100 <- array(0, c(length(echam$longph100),length(echam$lon)))
dist.outer.gph1002 <- array(0, c(length(echam$longph100),length(echam$longph100)))
dist.outer.gph100_gph100_2 <- array(0, c(length(echam$longph100),length(echam$longph100_2)))
dist.outer.gph100_u200 <- array(0, c(length(echam$longph100),length(echam$lonu200)))
dist.outer.gph100_omega <- array(0, c(length(echam$longph100),length(echam$lonomega)))
dist.outer.gph100_omega_2 <- array(0, c(length(echam$longph100),length(echam$lonomega_2)))
dist.outer.gph100_2 <- array(0, c(length(echam$longph100_2),length(echam$lon)))
dist.outer.gph100_22 <- array(0, c(length(echam$longph100_2),length(echam$longph100_2)))
dist.outer.gph100_2_u200 <- array(0, c(length(echam$longph100_2),length(echam$lonu200)))
dist.outer.gph100_2_omega <- array(0, c(length(echam$longph100_2),length(echam$lonomega)))
dist.outer.gph100_2_omega_2 <- array(0, c(length(echam$longph100_2),length(echam$lonomega_2)))
dist.outer.u200 <- array(0, c(length(echam$lonu200),length(echam$lon)))
dist.outer.u2002 <- array(0, c(length(echam$lonu200),length(echam$lonu200)))
dist.outer.u200_omega <- array(0, c(length(echam$lonu200),length(echam$lonomega)))
dist.outer.u200_omega_2 <- array(0, c(length(echam$lonu200),length(echam$lonomega_2)))
dist.outer.omega <- array(0, c(length(echam$lonomega),length(echam$lon)))
dist.outer.omega2 <- array(0, c(length(echam$lonomega),length(echam$lonomega)))
dist.outer.omega_omega_2 <- array(0, c(length(echam$lonomega),length(echam$lonomega_2)))
dist.outer.omega_2 <- array(0, c(length(echam$lonomega_2),length(echam$lon)))
dist.outer.omega_22 <- array(0, c(length(echam$lonomega_2),length(echam$lonomega_2)))

for (i in 1:length(echam$lon)) {
  dist.outer[i,1:length(echam$lon)] <- compute_dist(echam$lon, echam$lat, echam$lon[i], echam$lat[i])
}
for (i in 1:length(echam$lonstream)) {
  dist.outer.stream[i,1:length(echam$lon)] <- compute_dist(echam$lon, echam$lat, echam$lonstream[i], echam$latstream[i])
  dist.outer.stream2[i,1:length(echam$lonstream)] <- compute_dist(echam$lonstream, echam$latstream, echam$lonstream[i], echam$latstream[i])
  dist.outer.stream_gph100[i,1:length(echam$longph100)] <- compute_dist(echam$longph100, echam$latgph100, echam$lonstream[i], echam$latstream[i])
  dist.outer.stream_gph100_2[i,1:length(echam$longph100_2)] <- compute_dist(echam$longph100_2, echam$latgph100_2, echam$lonstream[i], echam$latstream[i])
  dist.outer.stream_gph300[i,1:length(echam$longph300)] <- compute_dist(echam$longph300, echam$latgph300, echam$lonstream[i], echam$latstream[i])
  dist.outer.stream_u200[i,1:length(echam$lonu200)] <- compute_dist(echam$lonu200, echam$latu200, echam$lonstream[i], echam$latstream[i])
  dist.outer.stream_omega[i,1:length(echam$lonomega)] <- compute_dist(echam$lonomega, echam$latomega, echam$lonstream[i], echam$latstream[i])
  dist.outer.stream_omega_2[i,1:length(echam$lonomega_2)] <- compute_dist(echam$lonomega_2, echam$latomega_2, echam$lonstream[i], echam$latstream[i])
}
for (i in 1:length(echam$longph300)) {
  dist.outer.gph300[i,1:length(echam$lon)] <- compute_dist(echam$lon, echam$lat, echam$longph300[i], echam$latgph300[i])
  dist.outer.gph3002[i,1:length(echam$longph300)] <- compute_dist(echam$longph300, echam$latgph300, echam$longph300[i], echam$latgph300[i])
  dist.outer.gph300_gph100[i,1:length(echam$longph100)] <- compute_dist(echam$longph100, echam$latgph100, echam$longph300[i], echam$latgph300[i])
  dist.outer.gph300_gph100_2[i,1:length(echam$longph100_2)] <- compute_dist(echam$longph100_2, echam$latgph100_2, echam$longph300[i], echam$latgph300[i])
  dist.outer.gph300_u200[i,1:length(echam$lonu200)] <- compute_dist(echam$lonu200, echam$latu200, echam$longph300[i], echam$latgph300[i])
  dist.outer.gph300_omega[i,1:length(echam$lonomega)] <- compute_dist(echam$lonomega, echam$latomega, echam$longph300[i], echam$latgph300[i])
  dist.outer.gph300_omega_2[i,1:length(echam$lonomega_2)] <- compute_dist(echam$lonomega_2, echam$latomega_2, echam$longph300[i], echam$latgph300[i])
}  
for (i in 1:length(echam$longph100)) {
  dist.outer.gph100[i,1:length(echam$lon)] <- compute_dist(echam$lon, echam$lat, echam$longph100[i], echam$latgph100[i])
  dist.outer.gph1002[i,1:length(echam$longph100)] <- compute_dist(echam$longph100, echam$latgph100, echam$longph100[i], echam$latgph100[i])
  dist.outer.gph100_gph100_2[i,1:length(echam$longph100_2)] <- compute_dist(echam$longph100_2, echam$latgph100_2, echam$longph100[i], echam$latgph100[i])
  dist.outer.gph100_u200[i,1:length(echam$lonu200)] <- compute_dist(echam$lonu200, echam$latu200, echam$longph100[i], echam$latgph100[i])
  dist.outer.gph100_omega[i,1:length(echam$lonomega)] <- compute_dist(echam$lonomega, echam$latomega, echam$longph100[i], echam$latgph100[i])
  dist.outer.gph100_omega_2[i,1:length(echam$lonomega_2)] <- compute_dist(echam$lonomega_2, echam$latomega_2, echam$longph100[i], echam$latgph100[i])
}  
for (i in 1:length(echam$longph100_2)) {
  dist.outer.gph100_2[i,1:length(echam$lon)] <- compute_dist(echam$lon, echam$lat, echam$longph100_2[i], echam$latgph100_2[i])
  dist.outer.gph100_22[i,1:length(echam$longph100_2)] <- compute_dist(echam$longph100_2, echam$latgph100_2, echam$longph100_2[i], echam$latgph100_2[i])
  dist.outer.gph100_2_u200[i,1:length(echam$lonu200)] <- compute_dist(echam$lonu200, echam$latu200, echam$longph100_2[i], echam$latgph100_2[i])
  dist.outer.gph100_2_omega[i,1:length(echam$lonomega)] <- compute_dist(echam$lonomega, echam$latomega, echam$longph100_2[i], echam$latgph100_2[i])
  dist.outer.gph100_2_omega_2[i,1:length(echam$lonomega_2)] <- compute_dist(echam$lonomega_2, echam$latomega_2, echam$longph100_2[i], echam$latgph100_2[i])
}  
for (i in 1:length(echam$lonu200)) {
  dist.outer.u200[i,1:length(echam$lon)] <- compute_dist(echam$lon, echam$lat, echam$lonu200[i], echam$latu200[i])
  dist.outer.u2002[i,1:length(echam$lonu200)] <- compute_dist(echam$lonu200, echam$latu200, echam$lonu200[i], echam$latu200[i])
  dist.outer.u200_omega[i,1:length(echam$lonomega)] <- compute_dist(echam$lonomega, echam$latomega, echam$lonu200[i], echam$latu200[i])
  dist.outer.u200_omega_2[i,1:length(echam$lonomega_2)] <- compute_dist(echam$lonomega_2, echam$latomega_2, echam$lonu200[i], echam$latu200[i])
}
for (i in 1:length(echam$lonomega)) {
  dist.outer.omega[i,1:length(echam$lon)] <- compute_dist(echam$lon, echam$lat, echam$lonomega[i], echam$latomega[i])
  dist.outer.omega2[i,1:length(echam$lonomega)] <- compute_dist(echam$lonomega, echam$latomega, echam$lonomega[i], echam$latomega[i])
  dist.outer.omega_omega_2[i,1:length(echam$lonomega_2)] <- compute_dist(echam$lonomega_2, echam$latomega_2, echam$lonomega[i], echam$latomega[i])
}  
for (i in 1:length(echam$lonomega_2)) {
  dist.outer.omega_2[i,1:length(echam$lon)] <- compute_dist(echam$lon, echam$lat, echam$lonomega_2[i], echam$latomega_2[i])
  dist.outer.omega_22[i,1:length(echam$lonomega_2)] <- compute_dist(echam$lonomega_2, echam$latomega_2, echam$lonomega_2[i], echam$latomega_2[i])
}
# NOT sure if one matrix with all is needed later?
# <- rbind(dist.outer, dist.outer.stream, dist.outer.gph100, dist.outer.gph300, dist.outer.u200, dist.outer.omega, dist.outer.omega_2)


# localization only for first variable (d.weights)
d.weights.zero <- array(0, c(length(echam$lon),length(echam$lon)))
d.weights <- corr_function(dist.outer, L=l_dist)
d.weights.precip <- corr_function(dist.outer, L=l_dist_precip)
d.weights.slp <- corr_function(dist.outer, L=l_dist_slp)

if (addbias) {
  nvarind <- 6+1 # manually entered for stream, gph100, gph300, u200, omega, omega_2 and bias
} else {
  nvarind <- 6 # manually entered for stream, gph100, gph300, u200, omega, omega_2
}  
nvar <- length(unique(echam$names))-nvarind
d.weights.stream <- corr_function(dist.outer.stream, L=l_dist)
d.weights.stream2 <- corr_function(dist.outer.stream2, L=l_dist)
d.weights.stream_gph300 <- corr_function(dist.outer.stream_gph300, L=l_dist)
d.weights.stream_gph100 <- corr_function(dist.outer.stream_gph100, L=l_dist)
d.weights.stream_gph100_2 <- corr_function(dist.outer.stream_gph100_2, L=l_dist)
d.weights.stream_u200 <- corr_function(dist.outer.stream_u200, L=l_dist)
d.weights.stream_omega <- corr_function(dist.outer.stream_omega, L=l_dist)
d.weights.stream_omega_2 <- corr_function(dist.outer.stream_omega_2, L=l_dist)
d.weights.gph300 <- corr_function(dist.outer.gph300, L=l_dist)
d.weights.gph3002 <- corr_function(dist.outer.gph3002, L=l_dist)
d.weights.gph300_gph100 <- corr_function(dist.outer.gph300_gph100, L=l_dist)
d.weights.gph300_gph100_2 <- corr_function(dist.outer.gph300_gph100_2, L=l_dist)
d.weights.gph300_u200 <- corr_function(dist.outer.gph300_u200, L=l_dist)
d.weights.gph300_omega <- corr_function(dist.outer.gph300_omega, L=l_dist)
d.weights.gph300_omega_2 <- corr_function(dist.outer.gph300_omega_2, L=l_dist)
d.weights.gph100 <- corr_function(dist.outer.gph100, L=l_dist)
d.weights.gph1002 <- corr_function(dist.outer.gph1002, L=l_dist)
d.weights.gph100_gph100_2 <- corr_function(dist.outer.gph100_gph100_2, L=l_dist)
d.weights.gph100_u200 <- corr_function(dist.outer.gph100_u200, L=l_dist)
d.weights.gph100_omega <- corr_function(dist.outer.gph100_omega, L=l_dist)
d.weights.gph100_omega_2 <- corr_function(dist.outer.gph100_omega_2, L=l_dist)
d.weights.gph100_2 <- corr_function(dist.outer.gph100_2, L=l_dist)
d.weights.gph100_22 <- corr_function(dist.outer.gph100_22, L=l_dist)
d.weights.gph100_2_u200 <- corr_function(dist.outer.gph100_2_u200, L=l_dist)
d.weights.gph100_2_omega <- corr_function(dist.outer.gph100_2_omega, L=l_dist)
d.weights.gph100_2_omega_2 <- corr_function(dist.outer.gph100_2_omega_2, L=l_dist)
d.weights.u200 <- corr_function(dist.outer.u200, L=l_dist)
d.weights.u2002 <- corr_function(dist.outer.u2002, L=l_dist)
d.weights.u200_omega <- corr_function(dist.outer.u200_omega, L=l_dist)
d.weights.u200_omega_2 <- corr_function(dist.outer.u200_omega_2, L=l_dist)
d.weights.omega <- corr_function(dist.outer.omega, L=l_dist)
d.weights.omega2 <- corr_function(dist.outer.omega2, L=l_dist)
d.weights.omega_omega_2 <- corr_function(dist.outer.omega_omega_2, L=l_dist)
d.weights.omega_2 <- corr_function(dist.outer.omega_2, L=l_dist)
d.weights.omega_22 <- corr_function(dist.outer.omega_22, L=l_dist)
# localization for all variables except indices
d.weights1 <- d.weights
# 1. column
d.weights1 <- rbind(d.weights,d.weights.zero,d.weights.zero,d.weights.stream,d.weights.gph300,d.weights.gph100,d.weights.gph100_2,d.weights.u200,d.weights.omega,d.weights.omega_2)
# 2. column
d.weights2 <- rbind(d.weights.zero,d.weights.precip,d.weights.zero,d.weights.stream,d.weights.gph300,d.weights.gph100,d.weights.gph100_2,d.weights.u200,d.weights.omega,d.weights.omega_2)
# 3. column
d.weights3 <- rbind(d.weights.zero,d.weights.zero,d.weights.slp,d.weights.stream,d.weights.gph300,d.weights.gph100,d.weights.gph100_2,d.weights.u200,d.weights.omega,d.weights.omega_2)
## 1. column
#d.weights1 <- rbind(d.weights,d.weights,d.weights,d.weights.stream,d.weights.gph300,d.weights.gph100,d.weights.gph100_2,d.weights.u200,d.weights.omega,d.weights.omega_2)
## 2. column
#d.weights2 <- rbind(d.weights,d.weights.precip,d.weights,d.weights.stream,d.weights.gph300,d.weights.gph100,d.weights.gph100_2,d.weights.u200,d.weights.omega,d.weights.omega_2)
## 3. column
#d.weights3 <- rbind(d.weights,d.weights,d.weights,d.weights.stream,d.weights.gph300,d.weights.gph100,d.weights.gph100_2,d.weights.u200,d.weights.omega,d.weights.omega_2)
# 4. column
d.weights4 <- rbind(t(d.weights.stream),t(d.weights.stream),t(d.weights.stream),t(d.weights.stream2),t(d.weights.stream_gph300),t(d.weights.stream_gph100),t(d.weights.stream_gph100_2),t(d.weights.stream_u200),t(d.weights.stream_omega),t(d.weights.stream_omega_2))
# 5. column
d.weights5 <- rbind(t(d.weights.gph300),t(d.weights.gph300),t(d.weights.gph300),d.weights.stream_gph300,t(d.weights.gph3002),t(d.weights.gph300_gph100),t(d.weights.gph300_gph100_2),t(d.weights.gph300_u200),t(d.weights.gph300_omega),t(d.weights.gph300_omega_2))
# 6. column
d.weights6 <- rbind(t(d.weights.gph100),t(d.weights.gph100),t(d.weights.gph100),d.weights.stream_gph100,d.weights.gph300_gph100,t(d.weights.gph1002),t(d.weights.gph100_gph100_2),t(d.weights.gph100_u200),t(d.weights.gph100_omega),t(d.weights.gph100_omega_2))
# 7. column
d.weights7 <- rbind(t(d.weights.gph100_2),t(d.weights.gph100_2),t(d.weights.gph100_2),d.weights.stream_gph100_2,d.weights.gph300_gph100_2,d.weights.gph100_gph100_2,t(d.weights.gph100_22),t(d.weights.gph100_2_u200),t(d.weights.gph100_2_omega),t(d.weights.gph100_2_omega_2))
# 8. column
d.weights8 <- rbind(t(d.weights.u200),t(d.weights.u200),t(d.weights.u200),d.weights.stream_u200,d.weights.gph300_u200,d.weights.gph100_u200,d.weights.gph100_2_u200,t(d.weights.u2002),t(d.weights.u200_omega),t(d.weights.u200_omega_2))
# 9. column
d.weights9 <- rbind(t(d.weights.omega),t(d.weights.omega),t(d.weights.omega),d.weights.stream_omega,d.weights.gph300_omega,d.weights.gph100_omega,d.weights.gph100_2_omega,d.weights.u200_omega,t(d.weights.omega2),t(d.weights.omega_omega_2))
# 10. column
d.weights10 <- rbind(t(d.weights.omega_2),t(d.weights.omega_2),t(d.weights.omega_2),d.weights.stream_omega_2,d.weights.gph300_omega_2,d.weights.gph100_omega_2,d.weights.gph100_2_omega_2,d.weights.u200_omega_2,d.weights.omega_omega_2,t(d.weights.omega_22))
if (addbias) {
  d.weights11 <- rep(1,dim(echam$data)[1]-1)
  d.weights12 <- cbind(d.weights1,d.weights2,d.weights3,d.weights4,d.weights5,d.weights6,d.weights7,d.weights8,d.weights9,d.weights10,d.weights11)
  d.weights_all <- rbind(d.weights12,t(rep(1,dim(echam$data)[1])))
} else {
  d.weights11 <- NULL
  d.weights12 <- NULL
  d.weights_all <- cbind(d.weights1,d.weights2,d.weights3,d.weights4,d.weights5,d.weights6,d.weights7,d.weights8,d.weights9,d.weights10) 
}  

save(d.weights_all, file='../data/weights/distance_weights.Rdata')
rm(d.weights1,d.weights2,d.weights3,d.weights4,d.weights5,d.weights6,d.weights7,d.weights8,d.weights9,d.weights10,d.weights11,d.weights12) 

}













# prox2grid <- function(prox, gridfile='~/unibe/data/cru/cru_echamgrid_temp2.nc'){
# 
# nc=open.nc(gridfile, write=F)
# print.nc(nc)
# tmp1=var.get.nc(nc, "temp2") # for CRU temp
# lon=dim(tmp1)[1]
# lat=dim(tmp1)[2]
# ts=dim(tmp1)[3]
# anocru=c(1:(lon*lat*ts))
# anoprox=c(1:(lon*lat*ts))
# for(i in 1:length(anocru)){anocru[i]=NA}
# for(i in 1:length(anoprox)){anoprox[i]=NA}
# dim(anocru)=c(lon,lat,ts)
# dim(anoprox)=c(lon,lat,ts)
# lonlist=var.get.nc(nc, "lon") 
# latlist=var.get.nc(nc, "lat") 
# #write("correlation log file", file = "corr.log", append = FALSE, sep = " ")
# 
# #read proxy data into same grid as cru
#   lon=prox$lon
#   lat=prox$lat
#   k=which(abs(lonlist-lon+0.001)==min(abs(lonlist-lon+0.001))) # +0.001 to avoid find 2 locations with same distance
#   l=which(abs(latlist-lat)==min(abs(latlist-lat)))
#   # cru and prox dates
#   tmp34=anocru[k,l,1:length(tmp21)]
#   tmp35=cbind(tmp21[tmp13],tmp34[tmp13])
#   tmp40=try(cor.test(tmp35[,2],prox[,2],use="pairwise.complete.obs"))
#   df=try(round(as.numeric(tmp40[2]),digits=2))
#   p=try(round(as.numeric(tmp40[3]),digits=2))
#   r=try(round(as.numeric(tmp40[4]),digits=2))
#       if (is.na(p)) {
#       write(as.character(filenames[i]), file = "no_corr.log", append = TRUE, sep = " ")
#       } else if (p < 0.05) {
#       write(paste(as.character(filenames[i]),"p:",p,"cor:",r,"df:",df), file = "corr.log", append = TRUE, sep = " ")
#     } else if (p > 0.05) { 
#       write(paste(as.character(filenames[i]),"p:",p,"cor:",r,"df:",df), file = "no_corr.log", append = TRUE, sep = " ")
#     }
#  }
#  }
# }
# 





# ------------------------------------------------------------------------------
# define giorgi regions
# ------------------------------------------------------------------------------
giorgi.names <- c(
  'Southern Hemisphere',
  'Extratropical Southern Hemisphere',
  'Northern Hemisphere',
  'North America',
  'South America',
  'Asia',
  'Arctic',
  'Antarctica',
  'Africa',
  'Australia',
  'Amazon Basin',
  'Southern South America',
  'Central America',
  'Western North America',
  'Central North America',
  'Eastern North America',
  'Alaska',
  'Greenland',
  'Mediterranean Basin',
  'Northern Europe',
  'Western Africa',
  'Eastern Africa',
  'Southern Africa',
  'Sahara',
  'Southeast Asia',
  'East Asia',
  'South Asia',
  'Central Asia',
  'Tibet',
  'North Asia',
  'Europe',
  'Extratropical Northern Hemisphere',
  'Earth',
  'Hadley cell strength',
  'Hadley cell poleward extend',
  'ITCZ location',
  'Midlatitude circulation',
  'Subtropical jet',
  'Stratospheric polar vortex region 1',
  'Stratospheric polar vortex region 2',
  'Pacific walker circulation region 1',
  'Pacific walker circulation region 2',
  'Dynamic Indian Monsoon index region 1',
  'Dynamic Indian Monsoon index region 2',
  'North Atlantic Oscillation location 1',
  'North Atlantic Oscillation location 2',
  'Pacific North American Pattern location 1',
  'Pacific North American Pattern location 2',
  'Pacific North American Pattern location 3',
  'Pacific North American Pattern location 4'
  )

giorgi.short <- c('SH','ESH','NH','NAM','SAM','ASI','ARC','ANT','AFR','AUS', 'AMZ', 'SSA', 'CAM', 'WNA', 'CNA', 'ENA', 'ALA', 'GRL', 'MED', 'NEU', 'WAF', 'EAF', 'SAF', 'SAH', 'SEA', 'EAS', 'SAS', 'CAS', 'TIB', 'NAS', 'EU', 'ENH','GLO','HC','HCL', 'ITCZ','MC','SJ','PV1','PV2','PWC1','PWC2','DIMI1','DIMI2','NAO1','NAO2','PNA1','PNA2','PNA3','PNA4')

giorgi.edges <- cbind(
  c( -180, 180, -90, 0),  
  c( -180, 180, -90, -20),  
  c( -180, 180, 0, 90),  
  c( -130, -60, 25, 60),
  c( -82, -34, -56, 12),
  c( 40, 180, 20, 70),
  c( -180, 180, 67, 90),
  c( -180, 180, -90, -67),
  c( -20, 65, -35, 30),
  c( 110, 155, -45, -11),
  c( -82, -34, -20, 12),
  c( -76, -40, -56, -20),
  c( -116, -83, 10, 30),
  c( -130, -103 , 30, 60),
  c( -103, -85, 30 , 50),
  c( -85, -60, 25, 50),
  c( -170, -103, 60, 72),
  c( -103, -10, 50, 85),
  c( -10, 40, 30, 48),
  c( -10, 40, 48, 75),
  c( -20, 22, -12, 18),
  c( 22, 52, -12, 18),
  c( -10, 52, -35, -12),
  c( -20, 65, 18, 30),
  c( 95 , 155, -11, 20),
  c( 100, 145, 20, 50),
  c( 65, 100, 5, 30),
  c( 40, 75, 30, 50),
  c( 75, 100, 30, 50),
  c( 40, 180, 50, 70),
  c( -10, 40, 30, 75),
  c( -180, 180, 20, 90),
  c( -180, 180, -90, 90),
  c(-180,180,0,30),       # 'Hadley cell strength'
  c(-180,180,0,30),       # 'Hadley cell poleward extend'
#  The poleward extent (HCL) of the Hadley circulation
#  is defined as the latitude at which the 850 hPa zonal mean meridional wind becomes poleward when 
#  moving northward from the latitude of maximum ??.
  c(-180,180,0,30),      # 'ITCZ location',
  c(-180,180,30,60),     # 'Midlatitude circulation',
  c(-180,180,0,50),      # 'Subtropical jet',
#  c(-180,180,75,90),     #'Stratospheric polar vortex',
  c(-180,180,75,90),    #'Stratospheric polar vortex 1',  
  c(-180,180,40,55),    #'Stratospheric polar vortex 2',  
#  c(-180,-100,-10,10),   # 'Pacific walker circulation',
  c(-180,-100,-10,10),     # 'Pacific walker circulation 1',
  c(100,150,-10,10),     # 'Pacific walker circulation 2',
#  c(40,80,-15,-5),       # 'Dynamic Indian Monsoon index',
  c(40,80,5,15),         # 'Dynamic Indian Monsoon index',
  c(70,90,20,30),        # 'Dynamic Indian Monsoon index',
  c(-27,-24,36,39),      # 'North Atlantic Oscillation 1',
  c(-24,-20,62,66),      # 'North Atlantic Oscillation 2',
#  difference in the standardised monthly SLP anomalies at Ponta Delgada (Azores, 37.7N, 25.6W) 
#  and Reykjavik (Iceland,64N, 22W).
  c(-162,-158,18,22), #'Pacific North American Pattern 1'
  c(-167,-163,43,47), #'Pacific North American Pattern 2'
  c(-117,-113,53,57), #'Pacific North American Pattern 3'
  c(-87,-83,28,32)    #'Pacific North American Pattern 4'
#  PNA = 0.25*(Z[20N,160W] ??? Z[45N,165W] + Z[55N,115W] ??? Z[30N,85W]) 
#  where Z is standardised 500 hPa GPH
  )

giorgi <- list()
for (i in 1:ncol(giorgi.edges)) giorgi[[i]] <- list(names=giorgi.names[i], edges=giorgi.edges[,i])
names(giorgi) <- giorgi.short

# H.giorgi <- compute_giorgi_H(giorgi, echam)
# rm(giorgi.edges)







getgridboxnum <- function(stat, echam) {
  m <- rep(NA,length(stat$lon))
  echamatprox.arr <- array(NA,c(length(stat$lon), length(stat$time)))
  for(i in 1:(length(stat$lon))){
    if ((i %% 100 == 0) | (i==(length(stat$lon)))) {
      print(paste("station number:",i))
    }
    plon <- stat$lon[i]
    plat <- stat$lat[i]
    clon <- echam$lon[!is.na(echam$lon)]
    clat <- echam$lat[!is.na(echam$lat)]
    k=which(abs(clon-plon+0.001)==min(abs(clon-plon+0.001))) 
    l=which(abs(clat-plat)==min(abs(clat-plat)))
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
    mtmp=k[which(match(k,l)>0)]
#    print(mtmp)
    if (length(mtmp) > 1) {
      print("ACHTUNG: more than 1 grid box identified. Only first one will be saved. 
            Probably due to 6 months state vector")
      mtmp <- mtmp[1]
    }
    if (length(mtmp) > 0) {
#      print("mtmp > 0")
      m[i] <- mtmp
      try(if (stat$names[i]=="precip") {
        m[i] <- m[i]+one_var_dim
      } else if (stat$names[i]=="slp") {
        m[i] <- m[i]+(2*one_var_dim)
      },silent = T)
    } else { m[i] <- NA }
  }
  invisible(m)
}  




################################################################################
# netcdfwrite function
################################################################################
netcdfwrite <- function(lon,lat,data,filename="data.nc",time=1,mv=-999) {
  
  # Writes two-dimensional map or three-dimensional array of data, longitude vector,
  # latitude vector and time vector into a netcdf file.
  #
  # Description:
  #
  #      Returns a netcdf file that contains a two-dimensional map or a
  #      three-dimensional array of data, a longitude vector, a latitude vector
  #      and a time vector
  #
  # Usage:
  #
  #      netcdfwrite(lon,lat,data,filename,time,mv)
  #
  # Arguments:
  #
  #        data:     two-dimensional map with first dimension longitude and
  #                  second dimension latitude or three-dimensional array
  #                  with first dimension longitude, second dimension latitude and
  #                  third dimension time, containing the data to be written
  #                  to the netcdf file
  #     lon:       a vector containing the longitudes of data
  #        lat:      a vector containing the latitudes of data
  #     time:       a vector containing the times of data in hours
  #                  since 1900-1-1 00:00:0.0 . Default is 1, which is
  #                  appropriate for writing a two-dimensional map
  #     filename: string with the name of the netcdf file to be created.
  #                  Default is "data.nc"
  #        mv:       Value attributed to missing data. Default is -999
  #
  # Output:
  #
  #     netcdf-file containing the data
  #
  #
  # Example:
  #
  #     x <- seq(-20, 20, 5)
  #     y <- seq(30, 60, 5)
  #     dims <- c(length(x), length(y), 100)
  #     data <- array(rnorm(prod(dims)), dims)
  #     time <- 1:100
  #     # write a two-dimensional map of data to the netcdf file
  #     netcdfwrite(x,y,data[,,3],filename="data.nc")
  #     # write a three-dimensional array of data to the netcdf file
  #     netcdfwrite(x,y,data,filename="data.nc",time)
  
  
  
  if(length(time)==1){data<-array(data,c(dim(data)[1],dim(data)[2],1))}
  
  #create netcdf file
  #ncfile <- create.nc(filename,clobber)
  ncfile <- create.nc(filename)
  
  
  #give dimensions to variables
  dim.def.nc(ncfile,"lon",length(lon))
  dim.def.nc(ncfile,"lat",length(lat))
  dim.def.nc(ncfile,"time",length(time))
  
  
  # create variables
  var.def.nc(ncfile,"lon","NC_FLOAT","lon")
  var.def.nc(ncfile,"lat","NC_FLOAT","lat")
  var.def.nc(ncfile,"time","NC_DOUBLE","time")
  var.def.nc(ncfile,"data","NC_FLOAT",c("lon","lat","time"))
  #var.def.nc(ncfile,"data","NC_DOUBLE",c("lon","lat"))
  
  # add description information to variables
  att.put.nc(ncfile,"data","missing_value","NC_FLOAT",mv)
  
  att.put.nc(ncfile,"lon","long_name","NC_FLOAT","Longitude")
  att.put.nc(ncfile,"lon","units","NC_FLOAT","degrees_east")
  att.put.nc(ncfile,"lon","axis","NC_FLOAT","X")
  
  att.put.nc(ncfile,"lat","long_name","NC_FLOAT","Latitude")
  att.put.nc(ncfile,"lat","units","NC_FLOAT","degrees_north")
  att.put.nc(ncfile,"lat","axis","NC_FLOAT","Y")
  
  att.put.nc(ncfile,"time","long_name","NC_DOUBLE","Time")
  att.put.nc(ncfile,"time","units","NC_DOUBLE","hours since 1900-1-1 00:00:0.0")
  att.put.nc(ncfile,"time","delta_t","NC_DOUBLE","0000-01-00 00:00:00")
  
  # add data to variables
  var.put.nc(ncfile,"lat",lat)
  var.put.nc(ncfile,"lon",lon)
  var.put.nc(ncfile,"time",time)
  var.put.nc(ncfile,"data",data)
  
  # syncronize and close the netcdf file
  sync.nc(ncfile)
  close.nc(ncfile) 
}





add_contour <- function(x, cvarname = "gph500", ctype="ensmean", ti=1, ccol='black', 
                        clev=c(-100,100), cont_i=1, ... ) {
  dat.i <- which(x$names == cvarname)
  if (length(dat.i) == 1) ti <- 1:ncol(x$data)
  if (ctype == 'data'){
    if (length(dim(x$data)) == 3){
      plotdata <- x$data[dat.i,ti,cont_i]
    } else {
      plotdata <- t(x$data[dat.i,ti, drop=F])
    }
  } else if (ctype == 'ensmean') {
    if (is.matrix(x$ensmean)){
      plotdata <- as.matrix(x$ensmean[dat.i,ti])
    } else {
      plotdata <- as.matrix(x$ensmean[dat.i])
    }
  } else {
    stop("Only 'data' and 'ensmean' are implemented as type")
  }
  
  ## sort data into array
  lon <- x$lon[dat.i]
  lat <- x$lat[dat.i]
  ulon <- sort(unique(x$lon[dat.i]))
  ulat <- sort(unique(x$lat[dat.i]))
  
  dat.arr <- array(NA, c(length(ulon), length(ulat)))
  dat.arr[cbind(match(lon, ulon), match(lat, ulat))] <- plotdata
  
#  contour(ulon, ulat, dat.arr, col=ccol, add=T)
#  contour(ulon, ulat, dat.arr, levels = head(clev[clev>0],-1), lty=1, col=ccol, add=T)
#  contour(ulon, ulat, dat.arr, levels = tail(clev[clev<0],-1), lty=2, col=ccol, add=T)
  if (length(which(contlevs>0))>0) {
    contour(ulon, ulat, dat.arr, levels = clev[clev>0], lty=1, col=ccol, add=T)}
  if (length(which(contlevs<0))>0) {
    contour(ulon, ulat, dat.arr, levels = clev[clev<0], lty=2, col=ccol, add=T)}
  
  
#   lon <- x$lon	
#   lat <- x$lat
#   ulondtmp <- sort(unique(diff(lon)))
#   ulond <- abs(median(ulondtmp[ulondtmp<10 & ulondtmp > -10]))
#   ulatdtmp <- sort(unique(diff(lat)))
#   ulatd <- abs(median(ulatdtmp[ulatdtmp<10 & ulatdtmp > -10]))
#   ulon <- seq(min(lon), max(lon), ulond)
#   ulat <- seq(min(lat), max(lat), ulatd)
  
#  dat.arr <- array(NA, c(length(ulon), length(ulat)))
#  dat.arr[cbind(match(lon, ulon), match(lat, ulat))] <- plotdata
  
#  contour(ulon, ulat, dat.arr, add=T, ...)
}





plot_pages <- function(wl=31, reg, xa='n',box='n',seas='yrmean',anom=T,anomper=(1901:1980),lw=1,
                       sca=F,yl="Temperature [??C]",xl='Year',plotbem=T,plotmem=T,plotech=T,
                       pag=T,use18=T,title='',plotts=NULL){ 
  if (title=='') { title <- paste(eind.allts$names[reg],"anomalies") }
  # set memb 18 to NA because of SH bias
  if (!use18) {
    eind.allts$data[,,18]=eind.allts$data[,,18]*NA
    aind.allts$data[,,18]=aind.allts$data[,,18]*NA
  }
  # wl = running mean window length
  # reg = index of giorgi region
  # seas = yrmean or sum or win
  # lw = line width
  period <- eind.allts$time #seq(syr,eyr)
  if (pag) {
    pages <- read.table('../comparison_data/pages_continental_buch_stefan.txt',header=T,sep='\t')
    periodv <- pages[,1] #vind.allts$time #seq(1751,1980)
    pos <- match(period,periodv)
    pages <- pages[pos,-1]
    periodv <- periodv[pos]
    cnpages <- colnames(pages)
  }
  if (seas=='yrmean') {
    years <- rep(eind.allts$time,each=2)  
  } else {
    wincol <- is.odd(seq(1,ncol(aind.allts$ensmean)))
    somcol <- is.even(seq(1,ncol(aind.allts$ensmean)))  
    wincolv <- is.odd(seq(1,ncol(vind.allts$data)))
    somcolv <- is.even(seq(1,ncol(vind.allts$data)))    
  }
  if (seas=='yrmean') {
    print("annual mean")
    eindmean <- aggregate(eind.allts$ensmean[reg,],list(years),mean,na.rm=T)[,2]
    aindmean <- aggregate(aind.allts$ensmean[reg,],list(years),mean,na.rm=T)[,2]
    eindbem <- aggregate(eind.allts$bem[reg,],list(years),mean,na.rm=T)[,2]
    aindbem <- aggregate(aind.allts$bem[reg,],list(years),mean,na.rm=T)[,2]
    eindmin <- apply(aggregate(eind.allts$data[reg,,],list(years),mean,na.rm=T)[2:31],1,min,na.rm=T)
    eindmax <- apply(aggregate(eind.allts$data[reg,,],list(years),mean,na.rm=T)[2:31],1,max,na.rm=T)
    aindmin <- apply(aggregate(aind.allts$data[reg,,],list(years),mean,na.rm=T)[2:31],1,min,na.rm=T)
    aindmax <- apply(aggregate(aind.allts$data[reg,,],list(years),mean,na.rm=T)[2:31],1,max,na.rm=T)
    #  vinddat <- runmean(vind.allts$data[1,somcolv],wl)
    #  plot(period,eindbem,ty='l',col=rgb(0,0,0,10,maxColorValue=10), ylab="Temperature [??C]",
    #     xlab="Year", ylim=c(min(eindmin),max(eindmax)))
    #     tmp <- aggregate(eind.allts$data[reg,,],list(years),mean)
    #     plot(tmp[,1],tmp[,2],ty='l',lwd='2',col=rgb(5,5,5,3,maxColorValue=10),ylab="",xlab="",ylim=c(1,8))
    #     for (i in 1:(nmem-1)) {
    #       lines(tmp[,1],tmp[,(i+2)],ty='l',lwd='2',col=rgb(5,5,5,3,maxColorValue=10),ylab="",xlab="")
    #     }  
    #     lines(period,allind[,1],ty='l',lwd=lw,lty=1,col=rgb(0,0,0,10,maxColorValue=10), ylab="",xlab="")
    #   for (i in 1:nmem) {
    #     lines(period,eind.allts$data[1,somcol,i],ty='l',lwd='2',col=rgb(5,5,5,1,maxColorValue=10),
    #           ylab="",xlab="")
    #   }  
    #   for (i in 1:nmem) {
    #     lines(period,aind.allts$data[1,somcol,i],ty='l',lwd='1',lty=2,
    #     col=rgb(10,2,0,1,maxColorValue=10), ylab="",xlab="")
    #   }
  } else if (seas=='sum'){  
    print("summer season")
    eindmean <- eind.allts$ensmean[reg,somcol]  
    aindmean <- aind.allts$ensmean[reg,somcol]  
    eindbem <- eind.allts$bem[reg,somcol]
    aindbem <- aind.allts$bem[reg,somcol]
    eindmin <- apply(eind.allts$data[reg,somcol,],1,min,na.rm=T)
    eindmax <- apply(eind.allts$data[reg,somcol,],1,max,na.rm=T)
    aindmin <- apply(aind.allts$data[reg,somcol,],1,min,na.rm=T)
    aindmax <- apply(aind.allts$data[reg,somcol,],1,max,na.rm=T)
  } else if (seas=='win'){
    print("winter season")
    eindmean <- eind.allts$ensmean[reg,wincol]  
    aindmean <- aind.allts$ensmean[reg,wincol]  
    eindbem <- eind.allts$bem[reg,wincol]
    aindbem <- aind.allts$bem[reg,wincol]
    eindmin <- apply(eind.allts$data[reg,wincol,],1,min,na.rm=T)
    eindmax <- apply(eind.allts$data[reg,wincol,],1,max,na.rm=T)
    aindmin <- apply(aind.allts$data[reg,wincol,],1,min,na.rm=T)
    aindmax <- apply(aind.allts$data[reg,wincol,],1,max,na.rm=T)
  }
  print("join indices to allind")
  allind <- cbind(eindmean,aindmean,eindbem,aindbem,eindmin,eindmax,aindmin,aindmax)
  #   allind <- allind2
  #   # check errors
  #   std <- apply(allind2,2,sd,na.rm=T)
  #   mea <- apply(allind2,2,mean,na.rm=T)
  #   for (c in 1:ncol(allind2)) {
  #     for (r in 1:nrow(allind2)) {
  # #      if (!is.na(allind[r,c])) {if (allind[r,c] > mea[c]+2*std[c]) {allind[r,c]=NA}}
  # #      if (!is.na(allind[r,c])) {if (allind[r,c] < mea[c]-2*std[c]) {allind[r,c]=NA}}
  #       if (!is.na(allind2[r,c])) {if (allind2[r,c] > mea[c]+2*std[c]) {allind2[r,]=NA}}
  #       if (!is.na(allind2[r,c])) {if (allind2[r,c] < mea[c]-2*std[c]) {allind2[r,]=NA}}
  #     }
  #   } 
  #   allind2 <- allind2[which(!is.na(allind2[,1])),]
  if (anom){
    anopos <- period %in% anomper
    centerfac <- apply(allind[anopos,],2,mean)
    centerfac[5] <- centerfac[6] <- centerfac[3] <- centerfac[1]
    centerfac[7] <- centerfac[8] <- centerfac[4] <- centerfac[2]
    if (sca) {
      scalefac <- apply(allind[anopos,],2,sd)
      scalefac[5] <- scalefac[6] <- scalefac[3] <- scalefac[1]
      scalefac[7] <- scalefac[8] <- scalefac[4] <- scalefac[2]
    } else { scalefac <- rep(1,8) }  
    allindann <- allind
    if (pag) {
      pagesann <- pages 
    }
    print("scale allind")
    if (wl > 1) {
      allind <- runmean(scale(allind,center=centerfac,scale=scalefac),wl)
    } else {
      allind <- scale(allind,center=centerfac,scale=scalefac)
    }
    if (pag) {
      centerfacp <- apply(pages[anopos,],2,mean,na.rm=T)
      # min and max col center factors to pages mean 
      centerfacp[3:4] <- centerfacp[2] 
      centerfacp[8:9] <- centerfacp[7]
      centerfacp[13:14] <- centerfacp[12]
      centerfacp[16:17] <- centerfacp[15]
      centerfacp[20:21] <- centerfacp[19]
      centerfacp[24:25] <- centerfacp[23]
      centerfacp[28:29] <- centerfacp[27]
      if (sca) {
        scalefacp <- apply(pages[anopos,],2,sd,na.rm=T)
      } else { scalefacp <- rep(1,ncol(pages))}    
      print(wl)
      if (wl > 1) {
        pages <- runmean(scale(pages,center=centerfacp,scale=scalefacp),wl)
      } else {
        pages <- scale(pages,center=centerfacp,scale=scalefacp)  
      }
      colnames(pages) <- cnpages 
    }
  }
  print("create axis")
  print(period)
  print(allind[,2])
  if (plotmem) {
    plot(period,allind[,2],ty='l',col=rgb(0,0,0,0,maxColorValue=10), ylab=yl, xlab=xl,
         ylim=c(min(allind[,5],na.rm=T)-0.5,max(allind[,6],na.rm=T)+0.5),
         main=title,xaxt=xa,bty=box)
  } else {
    if (is.null(plotts)){
      plot(period,allind[,2],ty='l',col=rgb(0,0,0,0,maxColorValue=10), ylab=yl, xlab=xl,
           ylim=c(min(allind[,2],na.rm=T)-0.5,max(allind[,2],na.rm=T)+0.5),
           main=title,xaxt=xa,bty=box)
    } else {
      plot(period[(plotts[1]:plotts[2])],allind[(plotts[1]:plotts[2]),2],ty='l',
           col=rgb(0,0,0,0,maxColorValue=10), ylab=yl, xlab=xl,
           ylim=c(min(allind[,2],na.rm=T)-0.5,max(allind[,2],na.rm=T)+0.5),
           main=title,xaxt=xa,bty=box)
    }  
  }
  #  plot(period,allind[,1],ty='l',col=rgb(0,0,0,10,maxColorValue=10), ylab="Temperature [??C]",
  #       xlab="Year", ylim=c(min(allind[,5]),max(allind[,6])),main=title,
  #       'anomalies'),xaxt=xa,bty=box)
  if (plotmem) {
    print("plotmem")
    if (plotech) {
      polygon(c(period,rev(period)),c(allind[,5],rev(allind[,6])),density=NA,col=rgb(1,1,1,3,maxColorValue=10))
    }
    polygon(c(period,rev(period)),c(allind[,7],rev(allind[,8])),density=NA,col=rgb(10,0,0,3,maxColorValue=10))
  }
  #  lines(period,allind[,5],ty='l',lwd='2',lty=3,col=rgb(3,3,3,10,maxColorValue=10),xlab='',ylab='')  
  #  lines(period,allind[,6],ty='l',lwd='2',lty=3,col=rgb(3,3,3,10,maxColorValue=10),xlab='',ylab='')  
  #  lines(period,allind[,7],ty='l',lwd='2',lty=3,col=rgb(10,5,0,10,maxColorValue=10),xlab='',ylab='')  
  #  lines(period,allind[,8],ty='l',lwd='2',lty=3,col=rgb(10,5,0,10,maxColorValue=10),xlab='',ylab='')  
  if (plotbem) {
    print("plotbem")
    #    lines(period,allind[,3],ty='l',lwd=lw,lty=1,col=rgb(3,3,3,10,maxColorValue=10), ylab="",xlab="")
    lines(period,allind[,4],ty='l',lwd=lw,lty=1,col=rgb(4,0,4,10,maxColorValue=10), ylab="",xlab="")
  }
  print("plot ens. mean")
  if (plotech) {
    lines(period,allind[,1],ty='l',lwd=lw,lty=1,col=rgb(0,0,0,10,maxColorValue=10), ylab="",xlab="")
  }
  lines(period,allind[,2],ty='l',lwd=lw,lty=1,col=rgb(10,0,0,10,maxColorValue=10), ylab="",xlab="")
  if (pag) {
    return(list(period=period,periodv=periodv,allindann=allindann,pagesann=pagesann,
                allind=allind,pages=pages))
  } else {
    return(list(period=period,allindann=allindann, allind=allind))
  }
}






is.even <- function(x) x %% 2 == 0
is.odd <- function(x) x %% 2 != 0 









# 
# # NEW read in year by year more efficiently, therefore put ensemble together later
# read_pdsi <- function(filehead, path=echallvarpath, xlim=c(-180,180), ylim=c(-90,90), 
#                         timlim=c(1600, 2005), small=F, landonly=F, anom=F, clim=F, std=F){
#   # read in the land-sea mask of echam
#   # mask out sea grid boxes (no variability)
#   if (landonly){
#     nc <- nc_open(paste(echmaskpath, 'landseamask.nc', sep='/'))
#     nc2 <-nc_open(paste(echmaskpath, 'orography.nc', sep="/"))
#     lon <- nc$dim$lon$vals
#     lat <- nc$dim$lat$vals
#     lon[lon > 180] <- lon[lon > 180] - 360
#     loi <- which(lon >= xlim[1] & lon <= xlim[2])
#     lai <- which(lat >= ylim[1] & lat <= ylim[2])
#     if (small==T) {
#       # mulc for reading each 3rd grid cell to avoid memory problems
#       mulc <- floor(length(loi)/96)
#       loi <- loi[seq(ceiling(mulc/2),length(loi),mulc)]
#       lai <- lai[seq(ceiling(mulc/2), length(lai),mulc)]
#     }
#     lsm <- ncvar_get(nc)[loi, lai]
#     alt <- ncvar_get(nc2)[loi, lai]
#     nc_close(nc)
#     nc_close(nc2)
#     lon2 <- lon
#     lat2 <- lat
#   }
#   
#   files <- list.files(path, pattern=paste('^', filehead, sep=''), full.names=T)
#   tmp <- list()
#   ensmean <- 0
#   num <- 0
#   for (f in files){
#     num <- num + 1
#     print(f)
#     nc <- nc_open(f)
#     lon <- nc$dim$lon$vals
#     lat <- nc$dim$lat$vals
#     lon[lon > 180] <- lon[lon > 180] - 360
#     loi <- which(lon >= xlim[1] & lon <= xlim[2])
#     lai <- which(lat >= ylim[1] & lat <= ylim[2])
#     if (small==T) {
#       mulc <- floor(length(loi)/96)
#       loi <- loi[seq(ceiling(mulc/2),length(loi),mulc)]
#       lai <- lai[seq(ceiling(mulc/2), length(lai),mulc)]
#     }
#     addc <-0
#     year <- as.numeric(format(ncdf_times(nc) + addc, "%Y"))
#     month <- as.numeric(format(ncdf_times(nc) + addc, "%m"))
#     tim <- year + (month-0.5)/12
#     for (t in timlim[1]:(timlim[2]-1)) {
#       ti <- which(tim >= t & tim < (t[1]+1)) #+2))
#       #    }
#       #    ti <- which(tim >= timlim[1] & tim < (timlim[2]+1))
#       outdata <- NULL
#       names <- NULL
#       # T2m (land only), SLP, Precip, Stream (should be zonal mean),GPH500, GPH100, u200, omega500, u850, v850, T850
#       varname <- 'PDSI'
#       if (varname %in% names(nc$var)){
#         data <- ncvar_get(nc, varname, start=c(1,1,min(ti)),
#                              count=c(-1,-1,length(ti)))[loi, lai] #,]
#         if (landonly){
#           data <- array(data, c(length(lsm),1))[lsm>0.5]  #, dim(data)[3]))[lsm > 0.5] #, ]
#         } else {
#           data <- array(data, c(dim(data)[1]*dim(data)[2], dim(data)[3]))
#         }
#         print(dim(data))
#         outdata <- rbind(outdata, data)
#         names <- c(names, rep(varname, nrow(data)))
#       }
#     }  
#       if (small) {
#         save(outdata, file=paste("../data/pdsi/pdsi_",num,"_",t,"-",(t+1),
#                                    "_2ndgrid.Rdata",sep=""))
#       } else {
#         save(outdata, file=paste("../data/pdsi/pdsi_",num,"_",t,"-",(t+1),".Rdata",sep=""))
#       }
#     } # end time loop
#   } # end files loop
# 
#     if (landonly) {
#       echam <- list(data=array(unlist(tmp), c(dim(tmp[[1]]), length(files))),
#                     ensmean=ensmean, lon=rep(lon[loi], length(lai))[lsm > 0.5],
#                     lat=rep(lat[lai], each=length(loi))[lsm > 0.5],
#                     height=alt[lsm > 0.5], lsm.i=which(lsm > 0.5),
#                     time=time, names=names)
#     } else {
#       echam <- list(data=array(unlist(tmp), c(dim(tmp[[1]]), length(files))),
#                     ensmean=ensmean, lon=rep(lon[loi], length(lai)),
#                     lat=rep(lat[lai], each=length(loi)),
#                     height=rep(-999,(length(loi)*length(lai))),
#                     lsm.i=rep(-999,(length(loi)*length(lai))), time=time, names=names)
#     }
#     if (small) {
#       save(echam, file=paste("../data/pdsi/pdsi_",t,"-",(t+1),
#                                "_2ndgrid.Rdata",sep=""))
#     } else {
#       save(echam, file=paste("../data/pdsi/pdsi_",t,"-",(t+1),".Rdata",sep=""))
#     }
#   }
# }







read_pdsi <- function(filehead, path=echpath, xlim=c(-180,180), ylim=c(-90,90), timlim=c(1600, 1630), small=F, landonly=F){
  
  # read in the land-sea mask of echam
  # mask out sea grid boxes (no variability)
  if (landonly){
    nc <- nc_open(paste(echmaskpath, 'landseamask.nc', sep='/'))
    nc2 <-nc_open(paste(echmaskpath, 'orography.nc', sep="/"))
    lon <- nc$dim$lon$vals
    lat <- nc$dim$lat$vals
    lon[lon > 180] <- lon[lon > 180] - 360
    loi <- which(lon >= xlim[1] & lon <= xlim[2])
    lai <- which(lat >= ylim[1] & lat <= ylim[2])
    if (small==T) {
      # mulc for reading each 3rd grid cell to avoid memory problems
      mulc <- floor(length(loi)/96)
      loi <- loi[seq(ceiling(mulc/2),length(loi),mulc)]
      lai <- lai[seq(ceiling(mulc/2), length(lai),mulc)]
    }  
    lsm <- ncvar_get(nc)[loi, lai]
    alt <- ncvar_get(nc2)[loi,lai]
    nc_close(nc)
    nc_close(nc2)
    lon2 <- lon
    lat2 <- lat
  } 
  files <- list.files(path, pattern=paste('^', filehead, sep=''), full.names=T)
  tmp <- list()
  ensmean <- 0
  for (f in files){
    print(f)
    nc <- nc_open(f)
    lon <- nc$dim$lon$vals
    lat <- nc$dim$lat$vals
    lon[lon > 180] <- lon[lon > 180] - 360
    loi <- which(lon >= xlim[1] & lon <= xlim[2])
    lai <- which(lat >= ylim[1] & lat <= ylim[2])
    if (small==T) {
      mulc <- floor(length(loi)/96)
      loi <- loi[seq(ceiling(mulc/2),length(loi),mulc)]
      lai <- lai[seq(ceiling(mulc/2), length(lai),mulc)]
    }
    # change processing of the time coordinate
    # tim <- seq(1600 + 1/24, by=1/12, length=nc$dim$time$len)
    # ti  <- which(tim >= 1600 & tim < 1631)
    # glitches in time representation
    #    if (timlim[1] < 1900){
    #      addc <- -15
    #    } else {
    #      addc <- 3
    #    }
    addc <-0
    year <- as.numeric(format(ncdf_times(nc) + addc, "%Y"))
    month <- as.numeric(format(ncdf_times(nc) + addc, "%m"))
    tim <- year + (month-0.5)/12
    ti <- which(tim >= timlim[1] & tim < (timlim[2]+1))
    outdata <- NULL
    names <- NULL
    for (varname in c('PDSI')) { #temp2', 'precip', 'slp')){
      if (varname %in% names(nc$var)){
        if (length(nc$var[[varname]]$dim) == 3){
          data <- ncvar_get(nc, varname, start=c(1,1,min(ti)), 
                               count=c(-1,-1,length(ti)))[loi, lai,]
        } else if (length(nc$var[[varname]]$dim) == 4){
          data <- ncvar_get(nc, varname, start=c(1,1,1,min(ti)), 
                               count=c(-1,-1,1,length(ti)))[loi, lai,]
        }
        if (landonly){
          data <- array(data, c(length(lsm), dim(data)[3]))[lsm > 0.5, ]
        } else {
          data <- array(data, c(dim(data)[1]*dim(data)[2], dim(data)[3]))
        }
        if ((path == echpath) | (path == echallvarpath) | (path == echallvarpath)) {
          if (varname == 'temp2') data <- data - 273.15
          if (varname == 'precip') data <- data * 3600 * 24 * 30
          if (varname == 'slp') data <- data / 100
        }
        if (path == nceppath) {
          if (varname == 'temp2') data <- data - 273.15
          if (varname == 'precip') data <- data * 3600 * 24 * 30
          if (varname == 'slp') data <- data 
        }
        # compute seasonal averages
        #dat <- array(data[,11:(ncol(data)-2)], c(nrow(data), 6, ncol(data)/6 - 2))
        #data <- apply(dat, c(1,3), mean, na.rm=T)
        outdata <- rbind(outdata, data)
        names <- c(names, rep(varname, nrow(data)))
        #rm(data, dat)
      }
    }
#     for (varname in setdiff(names(nc$var), c('temp2', 'precip', 'slp'))){
#       data <- ncvar_get(nc, varname)[ti]
#       #dat <- array(data[11:(length(data) - 2)], c(6, length(data)/6 - 2))
#       #data <- apply(dat, 2, mean, na.rm=T)
#       outdata <- rbind(outdata, data)
#       names <- c(names, varname)
#       #rm(data, dat)
#     }
    # add in the global and european mean of all 3-d datasets
    #    lons <- rep(lon[loi], length(lai))[lsm > 0.5]
    #    lats <- rep(lat[lai], each=length(loi))[lsm > 0.5]
    # compute additional operators for extraction of global and European means
    #    nh.mn <- (lats > 20) * cos(lats/180*pi)
    #    nh.mn <- nh.mn/sum(nh.mn)
    #    neu.mn <- (lats > 48 & lats < 75 & lons > -10 & lons < 40)*cos(lats/180*pi)
    #    neu.mn <- neu.mn/sum(neu.mn)
    #    n.mn <- nrow(outdata) %/% length(lons) 
    #    H.mn <- array(0, c(2*n.mn, nrow(outdata)))
    #    for (i in 1:n.mn){
    #      H.mn[i,i*length(lons) - length(lons):1 + 1] <- nh.mn
    #      H.mn[n.mn+i,i*length(lons) - length(lons):1 + 1] <- neu.mn
    #    }
    #    outdata <- rbind(outdata, H.mn %*% outdata)
    # check whether all of temp2, precip, and slp are available
    #    add.names <- as.vector(t(outer(c('NH', 'NEU'), unique(names[seq(1,n.mn*length(lons))]), paste, sep='.')))
    #    names <- c(names, add.names)
    #time <- apply(array((tim[ti])[11:(length(ti)-2)], c(6, length(ti)/6 -2)), 2, mean, na.rm=T)
    time <- tim[ti]
    tmp[[which(files == f)]] <- outdata
    print(dim(outdata))
    if (f != files[length(files)]) ensmean <- ensmean + outdata
  }
  ensmean <- ensmean/(length(tmp)) #-1) why -1 ???
  if (landonly) {
    data <- list(data=array(unlist(tmp), c(dim(tmp[[1]]), length(files))),
                 ensmean=ensmean, lon=rep(lon[loi], length(lai))[lsm > 0.5],
                 lat=rep(lat[lai], each=length(loi))[lsm > 0.5], 
                 height=alt[lsm > 0.5], lsm.i=which(lsm > 0.5),
                 time=time, names=names)
  } else {
    data <- list(data=array(unlist(tmp), c(dim(tmp[[1]]), length(files))),
                 ensmean=ensmean, lon=rep(lon[loi], length(lai)),
                 lat=rep(lat[lai], each=length(loi)), 
                 height=rep(-999,(length(loi)*length(lai))), 
                 lsm.i=rep(-999,(length(loi)*length(lai))), time=time, names=names)
  }
  #  data <- list(data=array(unlist(tmp), c(dim(tmp[[1]]), length(files))),
  #               ensmean=ensmean, lon=rep(lon[loi], length(lai))[lsm > 0.5],
  #               lat=rep(lat[lai], each=length(loi))[lsm > 0.5],
  #               height=alt[lsm > 0.5], lsm.i=which(lsm > 0.5),
  #               time=time, names=names)
  invisible(data) 
}






read_20cr <- function(filehead, path=twentycrpath, xlim=c(-180,180), ylim=c(-90,90), timlim=c(1980, 2000), small=F, landonly=F){
  
  # read in the land-sea mask of echam
  # mask out sea grid boxes (no variability)
  if (landonly){
    nc <- nc_open(paste(echmaskpath, 'landseamask.nc', sep='/'))
    nc2 <-nc_open(paste(echmaskpath, 'orography.nc', sep="/"))
    lon <- nc$dim$lon$vals
    lat <- nc$dim$lat$vals
    lon[lon > 180] <- lon[lon > 180] - 360
    loi <- which(lon >= xlim[1] & lon <= xlim[2])
    lai <- which(lat >= ylim[1] & lat <= ylim[2])
    if (small==T) {
      # mulc for reading each 3rd grid cell to avoid memory problems
      mulc <- floor(length(loi)/96)
      loi <- loi[seq(ceiling(mulc/2),length(loi),mulc)]
      lai <- lai[seq(ceiling(mulc/2), length(lai),mulc)]
    }  
    lsm <- ncvar_get(nc)[loi, lai]
    alt <- ncvar_get(nc2)[loi,lai]
    nc_close(nc)
    nc_close(nc2)
    lon2 <- lon
    lat2 <- lat
  } 
  files <- list.files(path, pattern=paste('^', filehead, sep=''), full.names=T)
  tmp <- list()
  ensmean <- 0
  for (f in files){
    print(f)
    nc <- nc_open(f)
    lon <- nc$dim$lon$vals
    lat <- nc$dim$lat$vals
    lon[lon > 180] <- lon[lon > 180] - 360
    loi <- which(lon >= xlim[1] & lon <= xlim[2])
    lai <- which(lat >= ylim[1] & lat <= ylim[2])
    if (small==T) {
      mulc <- floor(length(loi)/96)
      loi <- loi[seq(ceiling(mulc/2),length(loi),mulc)]
      lai <- lai[seq(ceiling(mulc/2), length(lai),mulc)]
    }
    # change processing of the time coordinate
    # tim <- seq(1600 + 1/24, by=1/12, length=nc$dim$time$len)
    # ti  <- which(tim >= 1600 & tim < 1631)
    # glitches in time representation
    #    if (timlim[1] < 1900){
    #      addc <- -15
    #    } else {
    #      addc <- 3
    #    }
    addc <-0
    year <- as.numeric(format(ncdf_times(nc) + addc, "%Y"))
    month <- as.numeric(format(ncdf_times(nc) + addc, "%m"))
    tim <- year + (month-0.5)/12
    ti <- which(tim >= timlim[1] & tim < (timlim[2]+1))
    outdata <- NULL
    names <- NULL
    for (varname in c('temp2', 'precip', 'slp')){
      if (varname %in% names(nc$var)){
        if (length(nc$var[[varname]]$dim) == 3){
          if (length(ti)==1) {
            data <- ncvar_get(nc, varname, start=c(1,1,min(ti)), 
                                 count=c(-1,-1,length(ti)))[loi, lai]
          } else {
            data <- ncvar_get(nc, varname, start=c(1,1,min(ti)), 
                               count=c(-1,-1,length(ti)))[loi, lai,]
          }
        } else if (length(nc$var[[varname]]$dim) == 4){
          if (length(ti)==1) {
            data <- ncvar_get(nc, varname, start=c(1,1,1,min(ti)), 
                                 count=c(-1,-1,1,length(ti)))[loi, lai]
          } else {  
            data <- ncvar_get(nc, varname, start=c(1,1,1,min(ti)), 
                               count=c(-1,-1,1,length(ti)))[loi, lai,]
          }
        }
        if (landonly){
          if (length(ti)==1) {
            data <- array(data, c(length(lsm)))[lsm > 0.5, ]
          } else {  
            data <- array(data, c(length(lsm), dim(data)[3]))[lsm > 0.5, ]
          }
        } else {
          if (length(ti)==1) {
            data <- array(data, c(dim(data)[1]*dim(data)[2]))
          } else {
            data <- array(data, c(dim(data)[1]*dim(data)[2], dim(data)[3]))
          }
        }
        if (path == twentycrpath) {
          if (varname == 'precip') data <- data * 3600 * 24 * 30
          if (varname == 'slp') data <- data / 100
        }
        # compute seasonal averages
        #dat <- array(data[,11:(ncol(data)-2)], c(nrow(data), 6, ncol(data)/6 - 2))
        #data <- apply(dat, c(1,3), mean, na.rm=T)
        if (length(ti)==1) {
          outdata <- c(outdata, data)  
        } else {
          outdata <- rbind(outdata, data)
        }
        names <- c(names, rep(varname, nrow(data)))
        #rm(data, dat)
      }
    }
    time <- tim[ti]
    tmp[[which(files == f)]] <- outdata
    if (length(ti)==1) {
      print(length(outdata))
    } else {  
      print(dim(outdata))
    }
    if (f != files[length(files)]) ensmean <- ensmean + outdata
  }
  ensmean <- ensmean/(length(tmp)) #-1) why -1 ???
  if (landonly) {
    if (length(ti)==1) {
      data <- list(data=unlist(tmp),
                   ensmean=ensmean, lon=rep(lon[loi], length(lai))[lsm > 0.5],
                   lat=rep(lat[lai], each=length(loi))[lsm > 0.5], 
                   height=alt[lsm > 0.5], lsm.i=which(lsm > 0.5),
                   time=time, names=names)
    } else {  
      data <- list(data=array(unlist(tmp), c(dim(tmp[[1]]), length(files))),
                 ensmean=ensmean, lon=rep(lon[loi], length(lai))[lsm > 0.5],
                 lat=rep(lat[lai], each=length(loi))[lsm > 0.5], 
                 height=alt[lsm > 0.5], lsm.i=which(lsm > 0.5),
                 time=time, names=names)
    }
  } else {
    if (length(ti)==1) {
      data <- list(data=unlist(tmp),
                   ensmean=ensmean, lon=rep(lon[loi], length(lai)),
                   lat=rep(lat[lai], each=length(loi)), 
                   height=rep(-999,(length(loi)*length(lai))), 
                   lsm.i=rep(-999,(length(loi)*length(lai))), time=time, names=names) 
    } else {  
      data <- list(data=array(unlist(tmp), c(dim(tmp[[1]]), length(files))),
                 ensmean=ensmean, lon=rep(lon[loi], length(lai)),
                 lat=rep(lat[lai], each=length(loi)), 
                 height=rep(-999,(length(loi)*length(lai))), 
                 lsm.i=rep(-999,(length(loi)*length(lai))), time=time, names=names)
    }
  }
  invisible(data) 
}


# plot_echam2 should also allow for single plots
plot_echam2 <- function(x, levs, varname='temp2', type='data',  ti=1, stations=NULL,
                       symmetric=T, siglev=0.05, names=NULL, main='', units='', 
                       cex.pt=2, st.cex=cex.pt*0.25, st.col=NULL, add=FALSE, 
                       cols=NULL, lwd=3, lty=1, seas=c(1,2), x_lab='', y_lab=NULL, 
                       rownames=NULL, colnames=NULL){
  oldpar <- par(no.readonly=TRUE)
  dat.i <- which(x$names == varname)
  # Jonas Jan. 2014 for new ECHAM format 
  lons <- x$lon[x$names == varname]
  lats <- x$lat[x$names == varname]
  if (length(dat.i) == 1) ti <- 1:ncol(x$data)
  if (type == 'data'){
    if (length(dim(x$data)) == 3){
      plotdata <- x$data[dat.i,ti,drop=F]
    } else if (length(dim(x$data)) == 2){
      plotdata <- x$data[dat.i,ti,drop=F]  
    } else {
      plotdata <- t(x$data[dat.i,ti,drop=F])
    }
  } else if (type == 'ensmean') {
    if (is.matrix(x$ensmean)){
      plotdata <- as.matrix(x$ensmean[dat.i,ti])
    } else {
      plotdata <- as.matrix(x$ensmean[dat.i])
    }
  } else {
    stop("Only 'data' and 'ensmean' are implemented as type")
  }
  if (length(dat.i) == 1){
    if (length(seas) > 1) {
      if (is.null(cols)){
        se.col <- rbfun(14)[c(3,12)]
      } else {
        se.col <- cols[seq(along=seas)]
      }
    } else {
      if (is.null(cols)){
        se.col <- rep(1,max(seas))
      } else {
        se.col <- rep(cols[1], max(seas))
      }
    }
    tis <- as.vector(apply(as.matrix(seas), 1, function(x,y) seq(x, length(y$time),2), y=x))
    if (ncol(plotdata) == 1) {
      if (!add) plot(0, type='n', xlim=range(x$time), ylim=range(plotdata[tis]),
                     xlab=xlab, ylab=if (is.null(ylab)) varname else ylab)
      for (se in 1:2) lines(x$time[seq(se,length(plotdata),2)], plotdata[seq(se,length(plotdata),2)], col=se.col[se], lwd=lwd, lty=lty)
    } else {
      d.range <- apply(plotdata, 1, range)
      if (!add) plot(0, type='n', xlim=range(x$time), ylim=range(d.range[,tis]),
                     xlab='', ylab=if (is.null(ylab)) varname else ylab)
      for (se in seas){
        for (dr in 1:2){
          lines(x$time[seq(se, ncol(d.range), 2)], d.range[dr,seq(se, ncol(d.range), 2)], col=se.col[se], lwd=lwd/2, lty=lty)
        }
      }
    }
  } else {
    nmod <- ncol(plotdata)
    nrow <- ceiling(nmod/4)
    ncol <- ceiling(nmod/nrow)+1
    nn <- nrow*(ncol-1)
    if (!add){
      dd <- lcm(2)
      if (nmod == 1) dd <- lcm(2.8)
      layout(matrix(c(1:nmod, rep(nmod+2, nn-nmod), rep(nmod+1, nrow)),
                    ncol, nrow, byrow=T), height=c(rep(5, ncol-1),dd))
      oma <- rep(0, 4)
      if (!is.null(rownames)) oma[3] <- 3
      if (!is.null(colnames)) oma[2] <- 3
      par(oma=oma)
    }
    par(mar=rep(0.5,4), cex.axis=1.4, cex.lab=1.4)
    
    if (missing(levs)) {
      if (symmetric){
        levs <- pretty(c(plotdata, -plotdata), 15)
      } else {
        levs <- pretty(plotdata, 15)
      }
    }
    if (is.null(cols)){
      if (varname == 'precip'){
        cols <- gbfun(length(levs) - 1)
      } else { 
        cols <- rbfun(length(levs) - 1)
      }
    } else {
      cols <- rep(cols, length.out=length(levs) - 1)
    }
    for (i in 1:nmod){
      plot(0, type='n', xlim=range(x$lon), ylim=range(x$lat),
           axes=F, xlab='', ylab='')
      li <- as.numeric(cut(plotdata[,i], breaks=levs))
      points(x$lon, x$lat, pch=15, col=cols[li], cex=cex.pt)
      if (!is.null(colnames)) {
        if (i <= length(colnames)) axis(3, at=mean(range(x$lon)), label=colnames[i], tick=FALSE, cex=1.4)
      }
      if (!is.null(rownames)) { 
        if (i%%length(colnames) == 1) axis(2, at=mean(range(x$lat)), label=rownames[ceiling(i/length(colnames))], tick=FALSE, cex=1.4)
      }
      world <- map(interior=F, plot=F)
      world$x[abs(diff(world$x)) > 180] <- NA
      lines(world)
#      map(add=T, region=c('Caspian','Great Lakes')) 
      if (!is.null(stations)){
        if (is.null(st.col)){
          st <- stations$data[,i]
          lli <- as.numeric(cut(st, breaks=levs))
          points(stations$lon, stations$lat, pch=21, bg=cols[lli], cex=st.cex*0.5)
          #          tst <- stations$data[stations$name=='temp',i]
          #          pst <- stations$data[stations$name=='precip',i]
          #          sst <- stations$data[stations$name=='slp',i]
          try(points(stations$lon[stations$names=='temp2'], 
                     stations$lat[stations$names=='temp2'], pch=2, col="red", cex=st.cex*2))
          try(points(stations$lon[stations$names=='precip'], 
                     stations$lat[stations$names=='precip'], pch=1, col="blue", cex=st.cex))
          try(points(stations$lon[stations$names=='slp'], 
                     stations$lat[stations$names=='slp'], pch=5, col="magenta", cex=st.cex*3))
        } else {
          points(stations$lon, stations$lat, pch=1, col=st.col, cex=st.cex)
        }
      }
      if (!is.null(x$p.value)){
        signif <- x$p.value[,ti,i] > siglev
        points(x$lon[signif], x$lat[signif], pch=4, cex=max(2, cex.pt*0.75))
      }
      if (!is.null(names)) text(min(x$lon), max(x$lat), names[i], adj=c(0.5,1), cex=1.4)
      box()
    }
    if (units == ''){
      par(mar=c(3,0.5,1.5,0.5), cex.axis=1.4, cex.lab=1.4)    
      plot_colourbar(col=cols, lev=levs, cex=1.4, xlab=units)
    } else {
      par(mar=c(3,8,1.5,8), cex.axis=1.4, cex.lab=1.4)  
      tmp <- list(col=cols, lev=levs, units=units)
      class(tmp) <- 'plotmap'  
      plot_colourbar(tmp, cex=1.4)
    }
    mtext(main, side=3, line=2, outer=TRUE, font=1, cex=1.4)
    if (!add) par(oldpar)
  }
}




FairSprErr2 <- function (ens, obs, obs.sd=NULL,tmean=T) {
  stopifnot(is.matrix(ens), is.vector(obs), nrow(ens) == length(obs))
  xmask <- apply(!is.na(ens), 1, any) & !is.na(obs)
  nens <- ncol(ens)
  spread <- mean(apply(ens[xmask, , drop = F], 1, sd, na.rm = T)^2,na.rm = T)
  if (is.null(obs.sd)) {
    error <- mean((obs - rowMeans(ens))^2, na.rm = T)
  } else if (!is.null(obs.sd) & tmean) {    
    error <- mean((obs - rowMeans(ens))^2, na.rm=T) - mean(obs.sd^2, na.rm=T)
  } else if (!is.null(obs.sd) & !tmean) {    
    error <- mean((obs - rowMeans(ens))^2 - (obs.sd**2), na.rm = T)
  } else {
    stop("ERROR in spread/error ratio calculation!")
  }
  return(sqrt((nens + 1)/nens * spread/error))
}




# calculated echam covariance over ensemble and all time step instead of just current year
echam_covar <- function(syr=1603,eyr=2004){
  if (syr<1603) {syr=1603; print("syr set to 1603")}
  for (cyr in (syr:eyr)) {
    print(cyr)
    if (every2grid){
      load(paste0(dataintdir,"analysis/analysis_",cyr,"_2ndgrid.Rdata"))
    } else {
      load(paste0(dataintdir,"analysis/analysis_",cyr,".Rdata"))
    }
    if (cyr == syr) {
      echanomallts <- echam.anom
      echanomallts$data <- echam.anom$data[,,]
    } else {
      echanomallts$data <- abind(echanomallts$data,echam.anom$data[,,],along=3)
    }
  }
  save(echanomallts,file="../data/echam/echallts_for_covar.Rdata")
  # calc covar in analysis code just for H.i, otherwise too large
}


read_PAGES <- function(type){
  get("fsyr")
  get("feyr")
  
  if (any(!is.na(match(type, "tree"))) & all(is.na(match(type, "coral"))) ) {
    print("generate_PAGES_tree")
    p_tree = read_pages(fsyr,feyr, archivetype ="tree", validate=pages_lm_fit)
    pagesprox = p_tree
    save(pagesprox, file=paste0(workdir,"../pages_tree_",fsyr,"-",feyr,"_",pages_lm_fit,".Rdata"))
  } 
  if (any(!is.na(match(type, "coral"))) & all(is.na(match(type, "tree"))) ) {
    print("generate_PAGES_coral")
    p_coral = read_pages(fsyr,feyr, archivetype ="coral",validate=pages_lm_fit)
    pagesprox = p_coral
    save(pagesprox, file=paste0(workdir,"../pages_coral_",fsyr,"-",feyr,"_",pages_lm_fit,".Rdata"))
  } 
  if (any(!is.na(match(type, "documents")))) {
    print("generate_PAGES_docu")
    p_docu = read_pages(fsyr,feyr, archivetype ="documents",validate=pages_lm_fit) # validate doesnt matter
    save(p_docu, file=paste0(workdir,"../pages_docu_",fsyr,"-",feyr,".Rdata"))
  } 
  if (any(!is.na(match(type, "instrumental")))) {
    print("generate_PAGES_inst")
    p_inst = read_pages(fsyr,feyr, archivetype ="instrumental", validate=pages_lm_fit) # validate doesnt matter
    save(p_inst, file=paste0(workdir,"../pages_inst_",fsyr,"-",feyr,".Rdata"))
  } 
  if (any(!is.na(match(type, "tree"))) & any(!is.na(match(type, "coral"))) )  {
    print("generate_PAGES_tree_&_coral")
    p_tree = read_pages(fsyr,feyr, archivetype ="tree", validate=pages_lm_fit)
    p_coral = read_pages(fsyr,feyr, archivetype ="coral", validate=pages_lm_fit)
    pagesprox <- list()
    pagesprox$data <- cbind(p_tree$data, p_coral$data)
    pagesprox$lon <- c(p_tree$lon, p_coral$lon)
    pagesprox$lat <- c(p_tree$lat, p_coral$lat)
    pagesprox$time <-p_tree$time
    pagesprox$mr <- rbind(p_tree$mr, p_coral$mr) 
    pagesprox$var_residu <- c(p_tree$var_residu, p_coral$var_residu)
    save(pagesprox, file=paste0(workdir,"../pages_tree_&_coral_",fsyr,"-",feyr,"_",pages_lm_fit,".Rdata"))
  }
  return(pagesprox)
}


## function to calculate 71 year climatology of list x

calculate_climatology <- function(x,cyr,subtracted,added,source){
  
  ti=which(floor(x$time)>=(cyr-(subtracted)) & floor(x$time)<=(cyr+added))
  
  if (source=="proxy"){
    tiv=which(floor(x$time)==cyr)
    sts=ti[1]
    ets=ti[length(ti)]
    x.clim = x
    x.clim$data=t(x$data[sts:ets,])
    x.clim$time=x$time[sts:ets]
  }
  # transform it to Oct-Septyears
  if (source=="inst"){
    x.clim = x
    x.clim$data=t(x$data)
    if (added != 0 ) {
      if (cyr < 1636) {
        y_start = agrep(paste0(1600,'.792'),as.character(x.clim$time))
      } else {
        y_start =
          agrep(paste0(cyr-subtracted,'.792'),as.character(x.clim$time))
      }
      if (cyr > 1970) {
        y_end = grep(paste0(2005,'.708'),as.character(x.clim$time)) #
        # 2012 should be changed to 2005, everywhere
      } else {
        y_end = grep(paste0(cyr+added,'.708'),as.character(x.clim$time))
      }
    } else {
      y_start =
        agrep(paste0(cyr-subtracted,'.792'),as.character(x.clim$time))
      y_end = grep(paste0(cyr+added,'.708'),as.character(x.clim$time))
    }
    x.clim$data = x.clim$data[,y_start:y_end]
    x.clim$time = x.clim$time[y_start:y_end]
    x.clim$data = t(x.clim$data)
  }
  return(x.clim)
}


## function to screen realprox data (+/- 5 std.)

screenstd <- function (x,cyr,source) { #sources: proxy/inst
  
  
  
  if (source=="proxy"){
    
    x.clim<-calculate_climatology(x,cyr,35,35,source)
    
    tiv=which(floor(x$time)==cyr)  
    x$data<-t(x$data)
    for (i in 1:length(x$lon)) {
      
      rpmean<- mean(x.clim$data[i,],na.rm=T)
      rpsd   <- sd(x.clim$data[i,],na.rm=T)
      #rpmean <- mean(realprox$data[,i],na.rm=T)
      #rpsd   <- sd(realprox$data[,i],na.rm=T)
      if ((!is.na(rpmean)) & (!is.na(rpsd))) {
        if ((!is.na(x$data[i,tiv])) & ((x$data[i,tiv] < rpmean-5*rpsd)
                                       | (x$data[i,tiv] > rpmean+5*rpsd))) {
          x$data[i,tiv] <- NA
          print(paste('proxy data', i, 'out of range'))
          write(paste('proxy data', i, 'out of range'),file=paste0('../log/',logfn),append=T)
          
        }
      }
    }
  }
  
  if (source=="inst"){
    
    x.clim<-calculate_climatology(x,cyr,36,35,source)
    jecham.sd<-get("echam.sd")
    echam_clim_mon_ensmean <- get("echam_clim_mon_ensmean")
    
    gpos <- getgridboxnum(x,echam.sd)
    
    # year from Oct to Sept for the climatology
    
    vtmp <- array(x.clim$data,c(12, nrow(x.clim$data)/12, dim(x.clim$data)[2]))
    stsv = round(dim(vtmp)[2]/2)
    tiv=which(floor(x$time)==cyr)  
    sts=tiv[1]
    
    for (i in 1:length(x$lon)) { 
      if (!is.na(gpos[i])) {
        m <- gpos[i]
        if (!is.na(m)) {
          
          
          for (j in 1:12) {
            # if bias corrected proxy/inst is outside echam ens range +- 5SD,
            # data point will not be assimilated at this time step
            if (!all(is.na(vtmp[j,,i])) ) { 
              biasm <- echam_clim_mon_ensmean[m,j] - mean(vtmp[j,,i],na.rm=T) 
              if (!is.na(biasm) & !is.na(vtmp[j,stsv,i]) ) {  # Veronika added the second term of the if
                #  if (((vtmp[j,((stsv-1)/12+1),i]+ biasm) < echam$ensmean[m,(j+12)]-5*echam.sd$data[m,j]) |
                #     ((vtmp[j,((stsv-1)/12+1),i]+ biasm) > echam$ensmean[m,(j+12)]+5*echam.sd$data[m,j])) {
                if (((vtmp[j,stsv,i]+biasm) < echam_clim_mon_ensmean[m,j]-5*echam.sd$data[m,j]) |
                    ((vtmp[j,stsv,i]+ biasm) > echam_clim_mon_ensmean[m,j]+5*echam.sd$data[m,j])) {
                  
                  x$data[sts-1+j,i]<-NA
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
                  
                }
              }
            }
          }
        }
      }
    }
  }
  return(x)
}



screendistance <- function (echam.sd,varname) {
  var <- get(varname)
  gpos <- getgridboxnum(var,echam.sd)
  for (i in 1:length(var$lon)) { # order of this and next if statement changed 2017-07-25
    if (!is.na(gpos[i])) {
      m <- gpos[i]
      d <- compute_dist(var$lon[i],var$lat[i],echam.sd$lon[m],echam.sd$lat[m])
      if ((!is.na(d)) & (d > 600)) { # check distance of assim data to next model grid box
        m=NA
        print(paste('inst data', varname, i, '>600km from echam grid box; set to NA'))
        write(paste('inst data', varname, i, '>600km from echam grid box; set to NA'),
              file=paste0('../log/',logfn),append=T)
        var$data[,i]<-NA
        
        
      }
    }
  }
  return(var)
}

# change array to have 6 months in state vector for winter and summer
# first winter starts in oct of syr
# 6 mon stat vectors for oct-mar and apr and sep
# works for echam format
convert_to_sixmonstatevector <- function(x,cyr){
  tmp21 <- array(x$data,c(dim(x$data)[1]*dim(x$data)[2],
                          dim(x$data)[3]))
  tmp22 <- tmp21[((9*dim(x$data)[1]+1):(dim(tmp21)[1]-(3*dim(x$data)[1]))),]
  x$data <- array(tmp22,c(dim(tmp22)[1]/(((dim(x$data)[2]/12)-1)*2),
                          (((dim(x$data)[2]/12)-1)*2),dim(tmp22)[2]))
  tmp31 <- array(x$ensmean,c(dim(x$ensmean)[1]*dim(x$ensmean)[2]))
  tmp32 <- tmp31[((9*dim(x$ensmean)[1]+1):(dim(tmp31)[1]-
                                             (3*dim(x$ensmean)[1])))]
  x$ensmean <- array(tmp32,c(dim(tmp32)[1]/(((dim(x$ensmean)[2]/12)-1)*2),
                             (((dim(x$ensmean)[2]/12)-1)*2)))
  x$time <- c(cyr,cyr+0.5)
  return(x)
}

# takes realprox-data and converts it into a [x,1:2] matrix for the two seasons
convert_to_2_seasons <- function(x,source){
  fsyr<-get("fsyr")
  feyr<-get("feyr")
  cyr<-get("cyr")
  if (source=="proxy"){
    x.allts <- x
    tmp1=t(x$data)
    points.on.SH <- which(x$lat<0)
    tmp2=array(NA,c(dim(tmp1)[1],2,dim(tmp1)[2]))
    tmp2[,2,]=tmp1
    if (!isempty(points.on.SH)){
    tmp2[points.on.SH,1,]<-tmp2[points.on.SH,2,]
    tmp2[points.on.SH,2,]<-NA
    }
    x.allts$data=array(tmp2,c(dim(tmp2)[1],dim(tmp2)[2]*dim(tmp2)[3]))
    x.allts$time=seq(fsyr,feyr+0.5,0.5) 
    ti=which(floor(x.allts$time)==cyr)
    sts=ti[1]
    ets=ti[length(ti)]      
    x$data=x.allts$data[,sts:ets]
    x$time=x.allts$time[sts:ets]
    x$names=rep("prox",dim(x.allts$data)[1])
  }
  if (source=="doc"){
    tmp1 <- array(t(x$data),c(dim(x$data)[1] *  dim(x$data)[2]))
    x$data <- array(tmp1,c(dim(tmp1)[1]/(dim(x$data)[1]/6),2))
    rm(tmp1)
    x$time <- c(cyr,cyr+0.5)
    x$names <- rep(x$names,6)
    x$lon <- rep(x$lon,6)
    x$lat <- rep(x$lat,6)
    x$des=rep('mon',nrow(x$data))
  }
  return(list(x=x,x.allts=x.allts))
}

convert_to_monthly <- function(dataset) {
  # converts dataset from sixmonstatevector back to twelve months (oktober-september)
  syr<-get("syr")
  eyr<-get("eyr")
  s<-12 # set back to 12 months for plotting
  tmptime <-  get("tmptime")
  season<-get("season")
  #reshape dataset
  
  if (length(dim(dataset$data))==3){ 
    dataset$data <- array(dataset$data, c((dim(dataset$data)[1]/6), dim(dataset$data)[2]*6, dim(dataset$data)[3]))
  }
  else{ #validate and calibrate do not have multiple ensemble members (one dim less)
    dataset$data <- array(dataset$data, c((dim(dataset$data)[1]/6), dim(dataset$data)[2]*6))  
  }
  if ("ensmean" %in% names(dataset)){  #some datasets do not have ensmean
  dataset$ensmean <- array(dataset$ensmean, c((dim(dataset$ensmean)[1]/6), dim(dataset$ensmean)[2]*6))
  }
  dataset$time <- tmptime[(season[2]+1):(length(tmptime)-4)]
  dataset$names <- dataset$names[1:dim(dataset$data)[1]]
  dataset$lon <- dataset$lon[1:(dim(dataset$data)[1])] #/length(unique(dataset$names)))]
  dataset$lat <- dataset$lat[1:(dim(dataset$data)[1])] #/length(unique(dataset$names)))]
  return(dataset)
}



