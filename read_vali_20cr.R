source("EnSRF_functions.R")





#air.2m.mon: NOAA-CIRES 20th Century Reanalysis version 2c Monthly Averages monthly ensemble mean of 2m from 1851 to 2014
#https://www.esrl.noaa.gov/psd/cgi-bin/db_search/DBSearch.pl?Dataset=NOAA-CIRES+20th+Century+Reanalysis+Version+2c+Monthly+Averages&Variable=Air+Temperature&group=1&submit=Search
#prate.mon : NOAA-CIRES 20th Century Reanalysis version 2 Monthly Averages: Ensemble mean from 1871-2012
twentycr.all <- read_20cr(filehead="",path=paste0(workdir,"../data/20cr/"),xlim=c(-180,180), ylim=c(-90,90), timlim=c(1871, 2004), small=F, landonly=F,calc_ensmean=FALSE)

save(twentycr.all,file=)
