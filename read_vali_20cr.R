rm(list=ls())

user <- system("echo $USER",intern=T)
print(paste('User:',user))
if (user=="veronika") {
  # workdir('/scratch/veronika/rerun/r_code')
  workdir ='/scratch3/veronika/reuse/reuse_git/' # where are the scripts from github
} else if (user=="lucaf") {
  workdir='/scratch3/lucaf/reuse/reuse_git/'
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


#air.2m.mon: NOAA-CIRES 20th Century Reanalysis version 2c Monthly Averages monthly ensemble mean of 2m from 1851 to 2014
#https://www.esrl.noaa.gov/psd/cgi-bin/db_search/DBSearch.pl?Dataset=NOAA-CIRES+20th+Century+Reanalysis+Version+2c+Monthly+Averages&Variable=Air+Temperature&group=1&submit=Search

#prate.mon : NOAA-CIRES 20th Century Reanalysis version 2 Monthly Averages: Ensemble mean from 1871-2012

#GeopotH : https://www.esrl.noaa.gov/psd/cgi-bin/db_search/DBSearch.pl?Dataset=NOAA-CIRES+20th+Century+Reanalysis+Version+2c+Monthly+Averages&Variable=Geopotential+Height&group=1&submit=Search
#          various pressure levels. NOAA-CIRES 20th Century Reanalysis version 2 Monthly Averages: Ensemble mean from 1851-2012

#Omega : https://www.esrl.noaa.gov/psd/cgi-bin/db_search/DBSearch.pl?Dataset=NOAA-CIRES+20th+Century+Reanalysis+Version+2c+Monthly+Averages&Variable=Omega+%28dp%2Fdt%29&group=1&submit=Search 
#        pressure levels. NOAA-CIRES 20th Century Reanalysis version 2 Monthly Averages: Ensemble mean from 1851-2012


# SOMEHOW COULDN'T FIND ANY SLP DATA ON NOOA

# var name anzeigen
# ncdump -h filename.nc
# oder
# cdo vardes filename
# cdo sinfo filename

# GEOPOTH
## Time_bnds variable löschen und 500 level auswählen
system(paste0("cdo -sellevel,500 -delname,time_bnds ",twentycrpath,"hgt.mon.mean.nc ",twentycrpath,"new.hgt.nc"))
## Interpolate onto echam grid
system(paste0("cdo -r remapcon,t63grid ",twentycrpath,"new.hgt.nc ",twentycrpath,"newnew.hgt.nc"))
## change variablename
system(paste0("cdo chname,hgt,gph500 ",twentycrpath,"newnew.hgt.nc ",twentycrpath,"newnewnew.hgt.nc"))

# TEMPERATURE
system(paste0("cdo delname,time_bnds ",twentycrpath,"air.2m.mon.mean.nc ",twentycrpath,"new.air.nc"))
system(paste0("cdo -r remapcon,t63grid ",twentycrpath,"new.air.nc ",twentycrpath,"newnew.air.nc"))
system(paste0("cdo chname,air,temp2 ",twentycrpath,"newnew.air.nc ",twentycrpath,"newnewnew.air.nc"))

# PRECIPITATION -> needs to be changed from rate to normal
system(paste0("cdo delname,time_bnds ",twentycrpath,"prate.mon.mean.nc ",twentycrpath,"new.prate.nc"))
system(paste0("cdo -r remapcon,t63grid ",twentycrpath,"new.prate.nc ",twentycrpath,"newnew.prate.nc"))
system(paste0("cdo chname,prate,precip ",twentycrpath,"newnew.prate.nc ",twentycrpath,"newnewnew.prate.nc"))

# SLP
system(paste0("cdo delname,time_bnds ",twentycrpath,"prmsl.mon.mean.nc ",twentycrpath,"new.prmsl.nc"))
system(paste0("cdo -r remapcon,t63grid ",twentycrpath,"new.prmsl.nc ",twentycrpath,"newnew.prmsl.nc"))
system(paste0("cdo chname,prmsl,slp ",twentycrpath,"newnew.prmsl.nc ",twentycrpath,"newnewnew.prmsl.nc"))

# OMEGA
system(paste0("cdo -sellevel,500 -delname,time_bnds ",twentycrpath,"omega.mon.mean.nc ",twentycrpath,"new.omega.nc"))
system(paste0("cdo -r remapcon,t63grid ",twentycrpath,"new.omega.nc ",twentycrpath,"newnew.omega.nc"))
system(paste0("cdo chname,omega,omega500 ",twentycrpath,"newnew.omega.nc ",twentycrpath,"newnewnew.omega.nc"))
#already has the right variable name

## select time range
system(paste0("cdo seldate,1901-01-01,2004-12-31 ",twentycrpath,"newnewnew.hgt.nc ",twentycrpath,"geopoth.nc"))
system(paste0("cdo seldate,1901-01-01,2004-12-31 ",twentycrpath,"newnewnew.air.nc ",twentycrpath,"temp2.nc"))
system(paste0("cdo seldate,1901-01-01,2004-12-31 ",twentycrpath,"newnewnew.prate.nc ",twentycrpath,"precip.nc"))
system(paste0("cdo seldate,1901-01-01,2004-12-31 ",twentycrpath,"newnewnew.prmsl.nc ",twentycrpath,"slp.nc"))
system(paste0("cdo seldate,1901-01-01,2004-12-31 ",twentycrpath,"newnewnew.omega.nc ",twentycrpath,"omega.nc"))

### MERGE ALL VARS 1901-2004
system(paste0("cdo merge ",twentycrpath,"temp2.nc ",twentycrpath,"precip.nc ",twentycrpath,"slp.nc ",twentycrpath,"geopoth.nc ",twentycrpath,"omega.nc ",twentycrpath,"twentycr_allvar_1901-2004.nc"))

## check data
system(paste0("cdo vardes ",twentycrpath,"twentycr_allvar_1901-2004.nc"))

### LOAD DATA 
twentycr.all <- read_20cr(filehead="twentycr_allvar_1901-2004",path=twentycrpath,xlim=c(-180,180), ylim=c(-90,90), timlim=c(1871, 2004), small=T, landonly=F,calc_ensmean=FALSE)

## check data
length(which(twentycr.all$names=="temp2"))
length(which(twentycr.all$names=="precip"))
length(which(twentycr.all$names=="slp"))
length(which(twentycr.all$names=="omega500"))
length(which(twentycr.all$names=="gph5"))

save(twentycr.all,file=paste0(twentycrpath,"twentycr_allvar_1901-2004_2ndgrid.Rdata"))





