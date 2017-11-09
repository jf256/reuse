expname="EKF400_v1.3_test_merged_code" # "EKF400_v1.3_full_res" #
# TODO
#  "mon_from_seas"               # can we get monthly res from seasonal proxies, 
                                 # maybe idealized pseudoproxy experiment
#  "EKF400_v1.4_new_error"       # run with new error covar estimates from Luca 
#  "PAGES_2_0_full"              # complete PAGES proxy database version 2.0
#  "PAGES_2_0_screened"          # screened version of PAGES DB 2.0, just good asian rec.
#  "all_data_2017"               # include ISTI, Petra's TRW and more docu and ... data
#
# DONE
#  "EKF400_v1.3_full_res"        # same but tps only and EVERY grid cell, 1901-2004
#  "EKF400_v1.3_no_forc_30_ens" # use all years as one big 30 member ensemble regardless of forcing, 1901-2004
#  "EKF400_v1.3_no_forc_100_ens" # use all years as one big 100 member ensemble regardless of forcing, 1901-2004
#  "EKF400_v1.3_no_forc_500_ens" # use all years as one big 500 member ensemble regardless of forcing, 1901-2004
#  "clim_covar_25_NO_loc_250mem" # echam covar 75% current year 25% climatology
                                 # but clim. calc. from n=200 instead of 500; code version 1.3
#  "EKF400_v1.3_corr_echam_clim" # echam climatology had wrong dates, 35 years shift
                                 # validation data had bug, full period
#  "clim_covar_75_NO_loc"        # echam covar 25% current year 75% climatology
#  "clim_covar_75_NO_with"       # code version 1.2
#  "clim_covar_25_NO_loc"        # echam covar 75% current year 25% climatology
#  "clim_covar_25_with_loc"      # code version 1.2
#  "clim_covar_50_NO_loc_250mem" # echam covar 50% current year 50% climatology
                                 # but clim. calc. from n=200 instead of 500
                                 # REPEAT: obviously something wrong (see spread recuction)
#  "clim_covar_50_NO_loc"        # echam covar 50% current year 50% climatology
#  "clim_covar_50_with_loc"      # run with version 1.1???
#  "clim_covar_100_with_loc"     # echam covar 0% current year 100% climatology
#  "clim_covar_100_NO_loc"       # run with version 1.1???

#  "EKF400_v1.2_corr_inst_screen"  # outlier screening of instr. data corrected,
                                   # now 71-yr window instead of current year only
                                   # proper vali data, too. Hopefully! NO!!!
#  "EKF400_v1.1_correct_71yr_anom" # corrected 71-yr anomaly in assim data 08-2017
#expname=paste(expname,format(Sys.time(),"%Y%m%d%H%M"),sep="_")




# Experiments/Runs
# 2017-07-27 climcal4/scratch4 1600-2005 with corrected 71-year anom. input data
#   DETAILS THAT WERE SOLVED:
#   - ti=which(floor(realprox.allts$time)>=(cyr-35) &
#       floor(realprox.allts$time)<=(cyr+35))
#     has probably no effect because proxy data are already cut in time before if code is 
#     run for short periods (syr to eyr). I.e. all my results are with assimilated ~20-year 
#     anomalies instead of the desired 70-year anomalies
#   - check syr and eyr in "mean(vtmp[j,,i])"
#   - calibrate.allts AND proxies.allts COMMENTED AS NOT NEEDED AND ALL DATA COVER DIFF. PERIODS
#   - pspline package commented on climcal4 because not available for installed r version
#   - check difference between old and new version!!!

# 2017-07-31 climcal4/scratch4 1941-1970 check error covar from climatology
#   instead of current year only to have less spurious correlation that were due to 
#   small 30 member ensemble. 
# - 50% and 100% error covar from climatology
# - WITH and WITHOUT localization
# - Now very 7th field from 1600-1900 period and 30 member ens.
#   which roughly equals a 500 member ensemble
# - run tps_only because super slow otherwise (1 year in ~2 days)

# Planned Experiments
# 2. - Paleodata only, NO instr. (after adding N-TREND, PAGES CORALS and more hist. DOCUMENTS)
# 3. - all year as one big ensemble (Hakim style)
# 4. - 60 ens. mem. vs 30 mem.
# 5. - estimate obs. err. (see below)
# 6. - check if assim. data out of ens. spread!
# 7. - check proposal!

# - change to 12-months state vector to assimilate corals???
# - or use n-trend instead of pages, 
# - or somehow only corals from pages as ice cores do not help
# FAVORIT: N-TREND, PAGES CORALS and more hist. DOCUMENTS

# estimate observation errors using "Simultaneous estimation of covariance inflation and 
# observation errors within an ensemble Kalman filter"
# Hong Li, Eugenia Kalnay and Takemasa Miyoshi

# change every2grid  to lowres, i.e. interpolation
# library(akima)
# interpp(x=x-coor_verctor,y=y-coor_vector,z=z-values,xo=x_out-coor_vector,yo=y_out-coor_vector)

#FIX:
# strange no of input data 1641-1643  
# [37,] 1639     0   30   164   35  129
# [38,] 1640     0   30   164   35  129
# [39,] 1641  1668   30    50   19   31
# [40,] 1642  1668   30    50   19   31
# [41,] 1643  1668   30    50   19   31
# [42,] 1644     0   30   164   35  129
# [43,] 1645     0   30   164   35  129

# why do mxd input data drop around 1650?

# JONAS IDEEN:
# instr. mitteln, aber fehler nicht reduzieren
# proxies und instr. assimilieren damit wenn einzelne monate NAs haben 
# immner noch proxy infos eingehen

# Olivia Nov. 2015: 
# Z500 and temp. anom. at same place = radiation is cause
# Z500 and temp. anom. shifted = advection is cause
# precip anom. in % instead of mm, otherwise just tropics visible


# current code in R only without Fortran!!!
# ATTENTION: to use the faster fortran EnSRF routine execute on the command line:
# R CMD SHLIB EnSRF.f95 uncomment/comment related code in EnSRF_functions.R
# R --arch=x86_64 CMD SHLIB EnSRF.f95 # to compile EnSRF.f95 on 64bit MAC


#####################################################################################
# general switches 
#####################################################################################
fsyr=1602 # full period start year for data generation that allows 71-yr anom. calc.
feyr=2004
syr_cru=1901
eyr_cru=2004
syr_recon=1750
eyr_recon=1900
syr_ncep=1948
eyr_ncep=2009
#syr_ind=1901
#eyr_ind=2004

#!!!ATTENTION: sixmonstatevector year starts in October of previous year (cyr-1)
sixmonstatevector=T    # 6 months of data in state vector for real proxy multiple 
# regression approach. ATTENTION: monthly has to be TRUE
if (sixmonstatevector) {
  s <- 2 # only 2 seasons to calculate but still monthly results in long state vector
} else {
  s <- 12
}
season=c(3,9)  # 3,9 = apr-sep and oct-mar, num=end month of season; # season=c(2,5,8,11)
nmem=30        # number of ensemble members

#####################################################################################
# EnSRF_data switches
#####################################################################################
# load or generate data from scratch
generate_ECHAM=F       # if TRUE -> orig. echam data is read
generate_ECHAM_1901_70=F # ECHAM data for bias calc with real_prox data 
generate_ECHAM_103=F   # ECHAM ens. mem. 103 with corrected land use forcing
generate_ECHAM_covar=F # generate_ECHAM all time step array long-term covariance")
generate_ECHAM_anom=F  # read echam anom, clim and sd from cdo
# ATTENTION if echam=T the proxy/instr. data have to be generated, too
generate_NCEP=F        # generate NCEP/NCAR reanalysis for independent verification
# next line not included yet: 
generate_20CR=F        # generate 20CR reanalysis for independent verification
#if ((syr > 1900) & (eyr < 2005)) {
# generate_CRUALLVAR=T  # if TRUE -> read CRUTEM3 temp, TS3 precip and 
#} else {               # HADSLP2 gridded instrumentals
generate_CRUALLVAR=F   # if FALSE -> cru_allvar.Rdata 
generate_HadCRU4=F     # HadCRU ens. SD for instr. uncertainty and error-spread ratio
#}
#if ((syr < 1901) & (eyr > 1749)) {
#  generate_LUTPAULKUT=T # if TRUE -> read Luterbacher, Pauling, 
#} else {               # Kuettel's temp. precip. and SLP 
generate_LUTPAULKUT=F # gridded seasonal recons (1750-1999)
#} 
#generate_ind_recon=F   # read Stefan's indices 1900-2000 from .txt to .RData
# use scripts in data_yuri to generate .Rdata files 
generate_t_yuri=F      # if TRUE -> yuri's temp. data collection including HISTALP is read
generate_slp_yuri=F    # if TRUE -> yuri's slp data collection is read
generate_GHCN=F        # if TRUE -> orig. GHCN stat. data is read; 
generate_GHCN_precip=F # if FALSE -> ghcn.Rdata 
generate_DOCUM=F        # if TRUE -> yuri's docu. data collection is read 
generate_PROXIES=F

yuri_temp=T          # yuri's data compilation, SLP always loaded
yuri_slp=T
ghcn_temp=T
ghcn_prec=F
trw_only=F           # Petra's TRW only
mxd_only=F           # Use only MXD tree ring proxies, NOT Petra's TRW
schweingr_only=F     # Use Schweingruber MXD grid only

# all available data selected above are automatically switched on when available in EnSRF_data


loc=T      # T = WITH localization, F without
covarclim=0 # set 50 or 100 [%] how much echam climatology covariance should be used
# default=0, i.e. current year covar from ECHAM ensemble
n_covar=500  # set sample size for covar calc, e.g. 250 or 500

calc_decorr_dist=F     # calculate decorrelation distance for each variable from ECHAM to set L
if (loc) {
  l_dist_temp2=1000*1.5  # factor *1.5 after stefans recommendation
  l_dist_slp=1800*1.5
  l_dist_precip=300*1.5
  l_dist_gph500=1800*1.5
  l_dist_gph100=2500*1.5
  l_dist_u850=1200*1.5
  l_dist_u200=1200*1.5
  l_dist_v850=1000*1.5
  l_dist_v200=1000*1.5
  l_dist_omega500=300*1.5
  l_dist_t850=1000*1.5
  l_dist_ind=999999 # precalculated indices should be removed
} else {
  l_dist_temp2=999999
  l_dist_slp=999999
  l_dist_precip=999999
  l_dist_gph500=999999
  l_dist_gph100=999999
  l_dist_u850=999999
  l_dist_u200=999999
  l_dist_v850=999999
  l_dist_v200=999999
  l_dist_omega500=999999
  l_dist_t850=999999
  l_dist_ind=999999 
}

# ATTENTION: landcorrected only works with anomaly_assim==T and every2grid==T!!!
  landcorr = F      # use simulation WITHOUT land use bug if TRUE
# how to treat multiple input series in same grid box
first_prox_per_grid=F  # first proxy per echam grid box ATTENTION: only this 
# or second next option (avg_prox_per_grid) can be TRUE
  firstproxres=10      # grid resolution for instr. stations (5 = echamgrid/5)
avg_prox_per_grid=T    # average more than one proxy per echam grid box 
                       # and calc proxy vs echam correlation
instmaskprox=F         # remove proxy data from grid boxes that have instr. data
reduced_proxies=F      # use every ??th (see code below) proxy record
every2grid=T           # only use every third grid cell of ECHAM, CRU validation, ...
land_only=F            # calc on land only
fasttest=F             # use even less data
tps_only=F             # only use temp, precip and slp in state vector, remove other vars
no_stream=T            # all echam vars but stream function as there is problem with 
#                       # 5/9 levels, which are in lat dimension before and after 1880
loo=F                  # leave-one-out validation 
if (loo) {tps_only=T;no_stream=F}  # reduce state vector for faster validation
no_forc_big_ens=F      # use all years as one big ensemble regardless of forcing like LMR
                       # ONLY works with next option load_71yr_anom=T
  n_no_forc=500         # ensemble size for no_forc LMR like experiment
#load_71yr_anom=T       # load 71yr echam anomalies calculated with cdo
#anom_reload=F          # reload anom calculated in R (next option)
#anom_save=F            # save anom calculated in R to reload next time
#if (load_71yr_anom==T) {
#  anom_reload=F
#  anom_save=F}
check_assimdata=T      # screen assimilation data before using it

if (no_stream & tps_only) {
  tps_only = F
  print('ACHTUNG: tps_only was set to FALSE')
}


# other options
scaleprox=T            # scale standardized docu and prox data the echam variance at location
anomaly_assim=T        # work with anomalies to avoid reg. const in state vector
nseas <- 12            # year with 12 months
check_dist=F           # test for ideal cut-off distance of spatial correlations
#H_non_lin=F           # new H operator that also allows non-linear functions
ana.enssize=F
NCEP_SOCOL=F

# choose validation data set
# ONLY one can be TRUE
# # next line not included yet: 
# if (eyr < 1750) {
#   vali=F                 # switch off prepplot if no vali data selected
# } else {
#   vali=T
# }
twcr_vali=F            # 20CR reanalysis data for validation
ncep_vali=F            # NCEP/NCAR reanalysis data for validation
# if ((syr > 1900) & (eyr < 2006)) {
#   cru_vali=T             # monthly CRU TS3 temp, precip and HADSLP2 gridded instrumentals (1901-2004)
# #  ind_recon=T            # Stefan's reconstructed indices until 1948 and NCAR reanalysis later added to CRU and NCEP
# } else {
#   cru_vali=F 
# #  ind_recon=F
# }
# #ind_recon=F
# if ((syr < 1901) & (eyr > 1749)) {
#   recon_vali=T           # seasonal luterbacher, pauling, kuettel recons (1750-1999)
# } else {
#   recon_vali=F
# }

#####################################################################################
# prepare plot switches
#####################################################################################
monthly_out=F    # if sixmonstatevector=T output is backtransformed to seasonal 
                 # average or monthly data if monthly_out=T 
calc_prepplot=F  # save half year averages calc from monthly data into /prepplot folder
  write_coor=F     # write ascii files with assimilated stations and data per ts
# maybe change files names for new EKF400 version "1.0" to "1.1"
# write_netcdf requires to run calc_prepplot before 
# best set load_prepplot=F
write_netcdf=F   # write entire EKF400 to NetCDF files
version="v1.3"   # set version number for netcdf file name
# v1.1 at DKRZ is experiment 1.2 here!
if (!monthly_out & write_netcdf) {
  write_netcdf=F
  print('ACHTUNG: write_netcdf set to FALSE because monthly_out=F')
}
# run next option "load_prepplot" for entire validation period, usually 
# 1902-2003, because it creates time series
load_prepplot=F  # ATTENTION check if folder prepplot on scratch contains monthly or seasonal data!
                 # saves image and only needs to be run once, afterward set "load_image=T" 
statyr=1904      # 1941 1850/69 year, when station network is kept constant
load_image=F     # directly load image for syr-eyr period: 1902-2001 or 1651-1750 image
calc_vali_stat=F # calculate validation statistics after preparation (set "load_image=T")
vali_plots=F     # source EnSRF_plots.R script 
ind_ECHAM=F      # delete/comment code in prepplot script and then delete switches here
ind_recon=F      # delete/comment code in prepplot script and then delete switches here


#####################################################################################
# plot switches
#####################################################################################
monthly=F
pseudoproxy=F
plot_dweights=F
write_nc=F
recalc <- F
reload <- F
plstat <- NULL #calibrate # NULL or calibrate
countseries <- F
#PAGES <- F          # write output for PAGES paper



# -----------------------------------------------------------------------------------------
