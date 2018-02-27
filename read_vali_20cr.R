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

# GEOPOTH 500
## Time_bnds variable löschen und 500 level auswählen
system(paste0("cdo -sellevel,500 -delname,time_bnds ",twentycrpath,"hgt.mon.mean.nc ",twentycrpath,"new.hgt.nc"))
## Interpolate onto echam grid
system(paste0("cdo -r remapcon,t63grid ",twentycrpath,"new.hgt.nc ",twentycrpath,"newnew.hgt.nc"))
## change variablename
system(paste0("cdo chname,hgt,gph500 ",twentycrpath,"newnew.hgt.nc ",twentycrpath,"newnewnew.hgt.nc"))

# GEOPOTH 100
system(paste0("cdo -sellevel,100 -delname,time_bnds ",twentycrpath,"hgt.mon.mean.nc ",twentycrpath,"new.hgt100.nc"))
system(paste0("cdo -r remapcon,t63grid ",twentycrpath,"new.hgt100.nc ",twentycrpath,"newnew.hgt100.nc"))
system(paste0("cdo chname,hgt,gph100 ",twentycrpath,"newnew.hgt100.nc ",twentycrpath,"newnewnew.hgt100.nc"))

# TEMPERATURE 2m
system(paste0("cdo delname,time_bnds ",twentycrpath,"air.2m.mon.mean.nc ",twentycrpath,"new.air.nc"))
system(paste0("cdo -r remapcon,t63grid ",twentycrpath,"new.air.nc ",twentycrpath,"newnew.air.nc"))
system(paste0("cdo chname,air,temp2 ",twentycrpath,"newnew.air.nc ",twentycrpath,"newnewnew.air.nc"))

# TEMPERATURE 500
system(paste0("cdo -sellevel,500 -delname,time_bnds ",twentycrpath,"air.mon.mean.nc ",twentycrpath,"new.air500.nc"))
system(paste0("cdo -r remapcon,t63grid ",twentycrpath,"new.air500.nc ",twentycrpath,"newnew.air500.nc"))
system(paste0("cdo chname,air,t500 ",twentycrpath,"newnew.air500.nc ",twentycrpath,"newnewnew.air500.nc"))

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

# u850
system(paste0("cdo -sellevel,850 -delname,time_bnds ",twentycrpath,"uwnd.mon.mean.nc ",twentycrpath,"new.u850.nc"))
system(paste0("cdo -r remapcon,t63grid ",twentycrpath,"new.u850.nc ",twentycrpath,"newnew.u850.nc"))
system(paste0("cdo chname,uwnd,u850 ",twentycrpath,"newnew.u850.nc ",twentycrpath,"newnewnew.u850.nc"))

# u200
system(paste0("cdo -sellevel,200 -delname,time_bnds ",twentycrpath,"uwnd.mon.mean.nc ",twentycrpath,"new.u200.nc"))
system(paste0("cdo -r remapcon,t63grid ",twentycrpath,"new.u200.nc ",twentycrpath,"newnew.u200.nc"))
system(paste0("cdo chname,uwnd,u200 ",twentycrpath,"newnew.u200.nc ",twentycrpath,"newnewnew.u200.nc"))

# v850
system(paste0("cdo -sellevel,850 -delname,time_bnds ",twentycrpath,"vwnd.mon.mean.nc ",twentycrpath,"new.v850.nc"))
system(paste0("cdo -r remapcon,t63grid ",twentycrpath,"new.v850.nc ",twentycrpath,"newnew.v850.nc"))
system(paste0("cdo chname,vwnd,v850 ",twentycrpath,"newnew.v850.nc ",twentycrpath,"newnewnew.v850.nc"))

# v200
system(paste0("cdo -sellevel,200 -delname,time_bnds ",twentycrpath,"vwnd.mon.mean.nc ",twentycrpath,"new.v200.nc"))
system(paste0("cdo -r remapcon,t63grid ",twentycrpath,"new.v200.nc ",twentycrpath,"newnew.v200.nc"))
system(paste0("cdo chname,vwnd,v200 ",twentycrpath,"newnew.v200.nc ",twentycrpath,"newnewnew.v200.nc"))


# TO DO:
# merge all and save


## select time range
system(paste0("cdo seldate,1901-01-01,2004-12-31 ",twentycrpath,"newnewnew.hgt.nc ",twentycrpath,"gph500.nc"))
system(paste0("cdo seldate,1901-01-01,2004-12-31 ",twentycrpath,"newnewnew.hgt100.nc ",twentycrpath,"gph100.nc"))
system(paste0("cdo seldate,1901-01-01,2004-12-31 ",twentycrpath,"newnewnew.air.nc ",twentycrpath,"temp2.nc"))
system(paste0("cdo seldate,1901-01-01,2004-12-31 ",twentycrpath,"newnewnew.air500.nc ",twentycrpath,"t500.nc"))
system(paste0("cdo seldate,1901-01-01,2004-12-31 ",twentycrpath,"newnewnew.prate.nc ",twentycrpath,"precip.nc"))
system(paste0("cdo seldate,1901-01-01,2004-12-31 ",twentycrpath,"newnewnew.prmsl.nc ",twentycrpath,"slp.nc"))
system(paste0("cdo seldate,1901-01-01,2004-12-31 ",twentycrpath,"newnewnew.omega.nc ",twentycrpath,"omega500.nc"))
system(paste0("cdo seldate,1901-01-01,2004-12-31 ",twentycrpath,"newnewnew.u850.nc ",twentycrpath,"u850.nc"))
system(paste0("cdo seldate,1901-01-01,2004-12-31 ",twentycrpath,"newnewnew.u200.nc ",twentycrpath,"u200.nc"))
system(paste0("cdo seldate,1901-01-01,2004-12-31 ",twentycrpath,"newnewnew.v850.nc ",twentycrpath,"v850.nc"))
system(paste0("cdo seldate,1901-01-01,2004-12-31 ",twentycrpath,"newnewnew.v200.nc ",twentycrpath,"v200.nc"))




### MERGE ALL VARS 1901-2004
system(paste0("cdo merge ",twentycrpath,"temp2.nc ",twentycrpath,"precip.nc ",twentycrpath,"slp.nc ",twentycrpath,"gph500.nc ",twentycrpath,"gph100.nc ",twentycrpath,"u850.nc ",twentycrpath,"u200.nc ",twentycrpath,"v850.nc ",twentycrpath,"v200.nc ",twentycrpath,"omega500.nc ",twentycrpath,"t500.nc ",twentycrpath,"twentycr_allvar_1901-2004.nc"))

## check data
system(paste0("cdo vardes ",twentycrpath,"twentycr_allvar_1901-2004.nc"))

### LOAD DATA 
twentycr.all <- read_20cr(filehead="twentycr_allvar_1901-2004",path=twentycrpath,xlim=c(-180,180), ylim=c(-90,90), timlim=c(1871, 2004), small=T, landonly=F,calc_ensmean=FALSE)

## check data
length(which(twentycr.all$names=="temp2"))
length(which(twentycr.all$names=="precip"))
length(which(twentycr.all$names=="slp"))
length(which(twentycr.all$names=="gph500"))
length(which(twentycr.all$names=="gph100"))
length(which(twentycr.all$names=="u850"))
length(which(twentycr.all$names=="u200"))
length(which(twentycr.all$names=="v850"))
length(which(twentycr.all$names=="v200"))
length(which(twentycr.all$names=="omega500"))
length(which(twentycr.all$names=="t500"))

save(twentycr.all,file=paste0(twentycrpath,"twentycr_allvar_1901-2004_2ndgrid.Rdata"))





