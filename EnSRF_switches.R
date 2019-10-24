######################################################################################
################################ EnSRF_switches.R ####################################
######################################################################################

expname="allinput_fullres_50climcovar_alloldvar" # "test_v2_oct17" #

version="v1.4"    # set code version number for experiment and netcdf file names
# v1.1 at DKRZ is experiment 1.2 here!
# v1.4 after merging JF proxy fixes and VV precip assimilation
#####################################################################################
# general switches 
#####################################################################################

fsyr=1602 # full period start year for data generation that allows 71-yr anom. calc.
feyr=2004
syr_cru=1901
eyr_cru=2004
syr_recon=1750
eyr_recon=1900
syr_twentycr=1901 #1836 # 20CR data start from 1836, currently only TPS, validation stats need to be adapted
  # currently statistics only work for same periods of CRU and 20CR (20CR data actually start from 1850)
eyr_twentycr=2004
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
# generate MODEL data in .Rdata format
generate_ECHAM=F                # if TRUE -> orig. echam data is read
generate_ECHAM_103=F            # ECHAM ens. mem. 103 with corrected land use forcing
generate_ECHAM_covar=F          # generate_ECHAM all time step array long-term covariance")
generate_ECHAM_anom=F           # read echam anom, clim and sd from cdo
generate_CCSM_last_mill_ens=F   # CCSM last millinnium ensemble for DAPS no forc./stat. offline pseudoproxy exp.

# ATTENTION if echam=T the proxy/instr. data have to be generated, too

# old switches, not working anymore
#generate_ECHAM_1901_70=F        # ECHAM data for bias calc with real_prox data 
#generate_NCEP=F                 # generate NCEP/NCAR reanalysis for independent verification

# generate VALIDATION data in .Rdata format
# next line not included yet: 
generate_20CR=F                # generate 20CR v3 reanalysis for independent verification
generate_CRUALLVAR=F             # if FALSE -> cru_allvar.Rdata 
generate_HadCRU4SD=F             # HadCRU ens. SD for instr. uncertainty and error-spread ratio
#generate_LUTPAULKUT=F           # gridded seasonal recons (1750-1999) NOT tested forever
generate_ind_recon=F             # read Stefan's indices 1900-2000 from .txt to .RData

# generate INPUT/ASSIMILATION data in .Rdata format
# use scripts in data_yuri to generate .Rdata files 
generate_instr_yuri=F            # if TRUE -> yuri's temp. data collection .Rdata file is created
statyr=1880                      # year, when GHCN/ISTI station network is kept constant
 generate_GHCN=F                 # if TRUE -> orig. GHCN stat. data .Rdata file is created 
 generate_GHCN_precip=F          
 generate_ISTI=F
generate_DOCUM=F                 # if TRUE -> angies's + yuri's docu. data collection .Rdata file is created 
# do we still need next 3 lines? why after generate_proxies in EnSRF_generate?
#generate_PAGES=F                # using the screened PAGES proxy dataset
                    # CODE from Roni missing to convert csv to RData 
#generate_NTREND=F   # CODE from Roni missing to convert csv to RData 
# REMEMBER to use 12 months temp. forward model and yearly output for next which
generate_PSEUDO=F               # DAPS PSEUDO PROXY EXPERIMENT: Set pseudo_prox to T further down (approx Line 145)
                                # if F but pseudo_prox below is T, then .RData file is loaded


# At the moment: any combination of read.these should be possible
# To do: noch ntrend machen mit regressionmonth etc - gemacht aber ntrend proxies nicht mehr im verz.
#        Testen: verschiedene regression_months mit und ohne t/p usw.
#        bei read_pages ist das directory von pagesproxies noch zu verallgemeinern
#        Abklären: ist in den jeweiligen  read_proxy_mxd zb t4 immer von t1-t12, weil sonst machen 
#                  die colnames nicht Sinn!!!!! Zb Zeile:
#                  Weil in Pages sind es nicht dieselben

generate_PROXIES=F
# for multiple parallel runs generate_PROXIES only has to be true for the first run
# once they are create they can be loaded with this switch=F
# NOTE: The proxy data set switches need to remain to automatically set real_prox=T in EnSRF data script
#if (generate_PROXIES==T) {
  # You can choose any combination of months and variable (T&P) for regression_months.
  # Then you can choose for each source whether it to be included or not. 
  # The resulting realprox$mr is a matrix of dimension [1:x,1:25], where x depends on your chosen sources.
  # The 25 originates from 12 months for T and 12 for P plus the intercept. 
  # The respective columns that were not chosen remain NA. 
  # For MXD and SCHWEINGR it only takes the temperature and leaves the precip. months NA.
  # PAGES_tree data also consists of location on the SH: if for ex. t4 (is chosen), it takes t10 (t4+6) 
  # for any locations with lat<0. 
  
  #regression_months = c('t.first', 't.second','t.third','t.fourth','t.fifth','t.sixth') 
  regression_months = c('t.first', 't.second','t.third','t.fourth','t.fifth','t.sixth', 
                        'p.first', 'p.second','p.third','p.fourth','p.fifth','p.sixth')
  
  #### PROXIES ####
  NTREND=T       # NTREND best tree data version 2018 (identical with 2015 paper)
  PAGES=T        # PAGESdata base version 2 from 01/2018
                 # works only with t and p regression months
  TRW_PETRA=T    # all TRW series from ITRDB and recalibrated by Petra
  # Following switches have not been adapted to code version 08/2019, 2 col var_residu etc.
    TRW=F          # 35 best TRW records from Petra's collection
    MXD=F          # additional MXD, not gridded Schweizgruber data
    SCHWEINGR=F    # Schweingruber/Briffa MXD grid
  #################



  # ATTENTION: precip. forward model only works with CRU and NOT GISS
  lm_fit_data = "BEST"      # can be BEST, CRU or GISS to calculate the reg coeff-s 
                            # (BEST has CRU SLP and prec, GISS is land+ocean temp only)
  type = c("coral","tree")  # only works with tree and coral (and both indiviually as well)
  #}  
  # END if generate_PROXIES


  # Nevin: May 2018 
  ####### SCREENING FOR PROXIES ##########
  # (For now only works for temperature)
  # Only either AIC or PVALUE can be TRUE, if both are set to TRUE only the AIC part will be run
  # if neither of AIC and PVALUE are TRUE then the full model is used without screening

  AIC=T                           # calculates linear regression models for different continuous 
  # subperiods and takes the best one according to the AIC value.
  # furthermore, it only keeps the signifcant models with pvalue>alpha (alpha set below)
  # Example T1-T6 AIC = 9, T2-T4 AIC=-1 (all combinations are respected)
  # => smallest AIC => best model=> if not significant => tree excluded

  PVALUE=T                        # calculates only the full regression model and only keeps the significant ones (pval>alpha)
  alpha=0.05                      # Significance level default: 0.05 or 0.01

  avg_realprox_per_grid=F         # if more than one tree is situated in one Echam-Gridcell an average 
  # of all treeringwidth is calculated before makeing the regression model. Because of the independed 
  # loading of the different datasets (ntrend, pages, petra), the average is only calculated taken 
  # from trees of the same dataset ->if all 3 datasets are used it can occur, that still 3 averaged 
  # trees are in one gridbox.


#### Instrumental Data ####
old_statvec = T
new_statvec = F          # has +: wetdays, block, cycfreq; -: v200, t500

yuri_temp = T            # yuri's data compilation, SLP always loaded
yuri_slp = T
yuri_prec = T
 inst_slp_err = sqrt(10) # instrumental slp error (10 is the variance of slp error)
ghcn_temp = T
 inst_t_err = sqrt(0.9)  # instrumental temp error (0.9 is the variance of temp error)
 isti_instead_ghcn = T    # switch from ghcn to isti (ghcn_temp must still be set to TRUE)
ghcn_prec = T
  ghcn_p_err = 0.3       # error in percent (based on US stations estimation should be 30%)
  ghcn_p_min = 10        # minimum error 10 mm
  precip_ratio = F       # if T assimilating ratio, if F assimilating the difference
  gauss_ana = F          # use Gaussian anamorphosis for precipitation ratio
  check_norm = F         # check whether the GA transformed values normally distributed and use only those that are
  ghcn_wday = F          # assimilating wetdays calculated from daily precip ghcn data
  ghnc_w_err = 2         # error number of days (based on US stations estimation should be 2 days)


#### Documentary Data ####
assim_docu = T           # use docu series from version 1 and series from angie 2019, all monthly resolution only! 
  docu_err = sqrt(0.25)  # equals 0.5 std. dev.


#### Pseudo Proxy Data ####
pseudo_prox=F                    # use DAPS Pseudo-Proxies and annual resolution Jan-Dez
  last_mill_prior=F              # use NCAR last Millennium ensemble as prior in stat. offline 
                                 # DAPS pseudoproxy experiment instead of ECHAM CCC400
                                 # only works with landonly=F
if (pseudo_prox) {
  if (generate_PROXIES==T) {stop("pseudo_prox and generate_PROXIES cannot be used together")}
  sixmonstatevector=F    # using annual mean
  s <- 1
  season=12              # year from jan-dec
  print("ATTENTION: sixmonstatevector has been set to FALSE for annual resolution pseudoproxy experiment!")
}
# all available data selected above are automatically switched on when available in EnSRF_data

# if (generate_PROXIES){
#   if ((generate_PAGES & PAGES) | (generate_NTREND & NTREND) | (trw_only) | (mxd_only) | (schweingr_only)){
#     stop("WARNING! These switches should not be set to TRUE simultaneously: 
#          generate_PAGES   & PAGES
#          generate_NTREND  & NTREND
#          generate_PROXIES & trw_only
#          generate_PROXIES & mxd_only
#          generate_PROXIES & schweingr_only
#          ")
#     }
#   }

#################
# To use a bigger ensemble for the background
no_forc_big_ens= F      # use all years as one big ensemble regardless of forcing like LMR
                        # ONLY works with next option load_71yr_anom=T
covarclim=50            # set 50 or 100 [%] how much echam climatology covariance should be used
                            # default=0, i.e. current year covar from ECHAM ensemble
cov_inflate = F         # inflate the PB matrix
  inflate_fac = 1.02      # the factor of covariance inflation
# Only used if no_forc_big_ens=T or covarclim>0
state = "changing"      # can be "static" or "changing" (static = the same big ens used for all year, changing = it is recalculated for every year)
n_covar=100             # set sample size for covar calc or for no_forc LMR like experiment, e.g. 250 or 500
PHclim_loc = T          # whether we want to localize the PHclim, only works if covarclim > 0
PHclim_lvec_factor = 2  # if PHclim_loc=T, we can use eg. 2times the distances as in the 30 ensemble member, at the moment only works for shape_wgt= "circle"
mixed_loc = F           # first combining Pb and Pclim then localizing
update_PHclim = T       # whether PHclim should be updated assimilating observation-by-observation
save_ananomallts = F    # in the covarclim exps if we update the climatology part -> whether to save the "climatological" analysis or not

if (covarclim==0) {
  update_PHclim = F
}

if((covarclim>0)&no_forc_big_ens){
  stop("Warning: covarclim>0 and no_forc_big_ens are both TRUE: These two experiments cannot be combined")
}


# Calculate decorr length -> was done already
calc_decorr_dist=F              # calculate decorrelation distance for each variable from ECHAM to set L
region = "global"               # region: where the decorrelation length should be calculated
# default = "global"
# can select: "golbal", "ENH", "ESH", "tropics", "lat_band", "lon_band"
cor_length_period = "annual"    # period: over which the decorrelation length should be calculated
# default = "annual
# can select: "annual", "summer", "winter"
# for corr_over_region function both region and cor_length_period is needed
# for compute_dist_2d function region is needed

# Localizing the 30 ensemble members: distance and shape
loc=T                           # T = WITH localization, F without
if (loc) {
  l_dist_temp2=1000*1.5         # factor *1.5 after stefans recommendation
  l_dist_slp=1800*1.5
  l_dist_precip=300*1.5
  #l_dist_precip=900
  l_dist_gph500=1800*1.5
  l_dist_gph100=2500*1.5
  l_dist_u850=1200*1.5
  l_dist_u200=1200*1.5
  l_dist_v850=1000*1.5
  l_dist_v200=1000*1.5
  l_dist_omega500=300*1.5
  # l_dist_t850=1000*1.5 #Roni: in the echam there is no t850 but t500. They refer to the same level but I forgot which one is the correct one
  l_dist_t500=1000*1.5
  l_dist_ind=999999 # precalculated indices should be removed
  l_dist_wdays = 300*1.5
  #l_dist_wdays = 900
  l_dist_blocks = 1800*1.5
  l_dist_cycfreq = 1800*1.5
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
  l_dist_t500=999999
  # l_dist_t850=999999 # see above
  l_dist_ind=999999 
}
shape_wgt = "circle"          # can be "circle" or "ellipse" depends on how we want to do the localization
# default is "circle"


# ATTENTION: landcorrected only works with anomaly_assim==T and every2grid==T!!!
landcorr = F                  # use simulation WITHOUT land use bug if TRUE

# how to treat multiple input series in same grid box
best2worst=T           # order proxies from best first to worst last based on residuals
                       # if FALSE then worst first and best last
best_tree_per_grid=T   # use only one real tree proxy record with lowest residuals per grid cell 
  besttreeres=0.1      # generates lon/lat grid with set resolution 0.1 mean only best proxy in ~10km2 grid
                       # NOT for corals that have two values per year because 2nd val would be deleted
first_inst_per_grid=F  # first instrumental observations per echam grid box 
                       # ATTENTION: only this or second next option (avg_prox_per_grid) can be TRUE
                       # "DON'T USE FIRST_INST_PER_GRID IF YOU WANT TO INCLUDE REAL PROXY DATA AND 
                       # INSTRUMENTALS AT THE SAME TIME")
  firstinstres=1       # grid resolution for instr. stations (5 = echamgrid/5)
avg_obs_per_grid=T     # average more than one observation per echam grid box 
                       # and calc proxy vs echam correlation
ins_tim_loc = T        # whether the instrumental obs-s should be localized in time or not
instmaskprox=F         # remove proxy data from grid boxes that have instr. data
every2grid=F           # only use every third grid cell of ECHAM, CRU validation, ...
land_only=F            # calc on land only
tps_only=F             # only use temp, precip and slp in state vector, remove other vars
tpsw_only=F            # only use temp, precip, slp and wetdays in state vector, remove other vars
no_stream=T            # all echam vars but stream function as there is problem with 
#                       # 5/9 levels, which are in lat dimension before and after 1880
loo=F                  # leave-one-out validation 

if (loo) {tps_only=T;no_stream=F}  # reduce state vector for faster validation
#load_71yr_anom=T               # load 71yr echam anomalies calculated with cdo
#anom_reload=F                  # reload anom calculated in R (next option)
#anom_save=F                    # save anom calculated in R to reload next time
#if (load_71yr_anom==T) {
#  anom_reload=F
#  anom_save=F}
check_assimdata=T               # screen assimilation data before using it

if (no_stream & tps_only) {
  stop("Either no_stream or tps_only has to be TRUE but not both!")
}else if(!tps_only &!no_stream){
  stop("Either no_stream or tps_only has to be TRUE but not both!")
}


# other options
scaleprox=T                     # scale standardized docu to echam variance at location
                                # not necessary for proxy data because of regression model
anomaly_assim=T                 # work with anomalies to avoid reg. const in state vector
if (anomaly_assim==F) {
  stop("Only working with anomaly_assim=T at the moment!")
}
# nseas <- 12                   # year with 12 months
check_dist=F                    # test for ideal cut-off distance of spatial correlations
#H_non_lin=F                    # new H operator that also allows non-linear functions
ana.enssize=F
#NCEP_SOCOL=F

# choose validation data sets saved in the analysis step (EnSRF_data):
# (all three can be selected simultaneously)
vali_cru=T
vali_twentycr=T # now 20CRv3
vali_recon=F
#####################################################################################







#####################################################################################
# prepare plot switches
#####################################################################################

tps_only_postproc=F             # due validation only with tps
validation_set=c("cru_vali") #,"twentycr_vali")    #can be set to cru_vali, or twentycr_vali or 
# both together c("cru_vali","twentycr_vali")
# choses which validation set should be used in the postprocessing and plots
monthly_out=T                   # if sixmonstatevector=T output is backtransformed to seasonal 
# yearly out is old switch from nevin. annual output automatically if pseudo_prox=T
  yearly_out=F                    # if both false, plots seasonal averages are calculated
                                # average or monthly data if monthly_out=T 
temporal_postproc=T             # save half year averages calc from monthly data into /prepplot folder
mergetime_indices=F             # if TRUE: indices are combined to allts variables for whole period 
                                # (e.g. also for 1604-2004) and saved into image folder for TS-plots
# run next option "load_prepplot" for entire validation period, usually 
# 1902-2003, because it creates time series
mergetime_fields=F              # if calc_prepplot has been run, load_prepplot can be used
# saves image and only needs to be run once, afterward set "load_images=T" 
load_images=F                   # directly load image for syr-eyr period: 1902-2001 or 1651-1750 image
                                # of merged indices and fields. NOT need if just calculated before
calc_vali_stat=F                # calculate validation statistics after preparation (set "load_image=T")
CRPS=F                          # calculate Continuous Ranked Probability Score
ind_anom=T                      # calculate indices from anomaly data

vali_plots=F                    # source EnSRF_plots.R script 

#ind_ECHAM=T                     # delete/comment code in prepplot script and then delete switches here
#ind_recon=F                     # delete/comment code in prepplot script and then delete switches here


write_coor=F                    # write ascii files with assimilated stations and data per ts
# maybe change files names for new EKF400 version "1.0" to "1.1"
# write_netcdf requires to run calc_postproc before 
# best set mergetime_*=F and load_image=F
write_netcdf=T                  # write entire EKF400 to NetCDF files
if (!monthly_out & write_netcdf) {
  write_netcdf=F
  print('ACHTUNG: write_netcdf set to FALSE because monthly_out=F')
}

#####################################################################################



#####################################################################################
# plot switches
#####################################################################################

monthly=F
#station_yr=1905                # specific year to plot the station locations (syr>=station_yr<=eyr)
#pseudoproxy=F                 # old from Nevin pseudoproxy try
#plot_dweights=F
#write_nc=F
#recalc <- F
#reload <- F
plstat <- NULL                  #calibrate # NULL or calibrate
countseries <- T
#PAGES <- F                     # write output for PAGES paper

# -----------------------------------------------------------------------------------------
