###############################################################################
###############################################################################
## MBG Combine Script
##
## Author: USERNAME
## Indicator: Combines dpt3_cov, dpt2_cond, and dpt1_cond to create dpt1_cov, dpt2_cov, dropout metrics, etc.
## Date: DATE
##
## Source:
##   source("FILEPATH")
###############################################################################
###############################################################################

###############################################################################
## SETUP
###############################################################################

## clear environment
rm(list=ls())

## Set repo location and indicator group
user               <- 'USERNAME'
core_repo          <- sprintf('FILEPATH')
indic_repo         <- sprintf('FILEPATH')

## sort some directory stuff
setwd(core_repo)
commondir      <- sprintf('FILEPATH')
package_list <- c(t(read.csv(sprintf('FILEPATH'),header=FALSE)))

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0('FILEPATH'))
mbg_setup(package_list = package_list, repos = core_repo)

# Custom load indicator-specific functions
source(paste0('FILEPATH'))

## Script-specific stuff begins here ##########################################

# Load objects from qsub 
#   - variable: region, holdout
#   - static: save_2_cov, save_doses, vaccine, doses, run_date, indicator_group

load_from_parallelize()

## Read config file (from sharedir) and save all parameters in memory
config <- load_config(repo            = indic_repo,
                      indicator_group = indicator_group,
                      indicator       = paste0(vaccine, doses, "_cov"),
                      post_est_only   = TRUE,           
                      run_date        = run_date)

## Ensure you have defined all necessary settings in your config
check_config()

## Create a few objects from options above
if (class(Regions) == "character" & length(Regions) == 1) Regions <- eval(parse(text=Regions))
if (class(summstats) == "character" & length(summstats) == 1) summstats <- eval(parse(text=summstats))
gaul_list <- get_adm0_codes(Regions, shapefile_version=modeling_shapefile_version, subnational_raking=subnational_raking)
if (class(year_list) == "character") year_list <- eval(parse(text=year_list))
test <- as.logical(test)

# Check if correct # doses for now
if (doses != 3) stop("Only set up for 3-dose vaccines for now.")

#############################################################
## Load cell pred objects ###################################
#############################################################
ig_dir          <- paste0("FILEPATH") 

if (doses == 3) {

  message("\nLoading cell pred objects...")

  message(paste0("  ", vaccine, "1_cond..."))
  load(paste0(ig_dir, vaccine, "1_cond/output/", run_date, "/", vaccine, "1_cond_cell_draws_eb_bin0_", region, "_", holdout, ".RData"))
  cpred_1 <- cell_pred
  rm(cell_pred)

  message(paste0("  ", vaccine, "2_cond..."))
  load(paste0(ig_dir, vaccine, "2_cond/output/", run_date, "/", vaccine, "2_cond_cell_draws_eb_bin0_", region, "_", holdout, ".RData"))
  cpred_2 <- cell_pred
  rm(cell_pred)

  message(paste0("  ", vaccine, "3_cov..."))
  load(paste0(ig_dir, vaccine, "3_cov/output/", run_date, "/", vaccine, "3_cov_cell_draws_eb_bin0_", region, "_", holdout, ".RData"))
  cpred_3 <- cell_pred
  rm(cell_pred)

}

#############################################################
## Apply continuation ratio method to calculate doses #######
#############################################################

message("\nApplying continuation ratio calculations & saving cell_preds... ")

if (doses == 3 & save_doses) {

  # Math double-checked here: http://bit.ly/2lvzdzG

  # Check to make sure dimensions equal
  if (!all.equal(dim(cpred_1), dim(cpred_2), dim(cpred_3))) stop("Dimensions of cell preds are not equal.")

  # vax3_dose = directly modeled, so nothing to do
  
  ### vax2_dose #############################################

  message(paste0("  ", vaccine, "2_dose..."))

  # The math: vax2_dose = vax2_cond - vax2_cond*vax3_cond
  cell_pred <- cpred_2 - (cpred_2*cpred_3) 


  # Save output in standard directory
  vax2_dose_out_dir <- paste0(ig_dir, vaccine, "2_dose/output/", run_date, "/")
  dir.create(vax2_dose_out_dir, recursive = T, showWarnings=F)
  save(cell_pred, file = paste0(vax2_dose_out_dir, vaccine, "2_dose_cell_draws_eb_bin0_", region, "_", holdout, ".RData"))
  
  rm(cell_pred)

  ### vax1_dose #############################################

  message(paste0("  ", vaccine, "1_dose..."))

  # The math: vax1_dose = vax1_cond(vax2_cond-1)(vax3_cond-1)
  cell_pred <- (cpred_1)*(cpred_2 - 1)*(cpred_3 - 1) 

  # Save output in standard directory
  vax1_dose_out_dir <- paste0(ig_dir, vaccine, "1_dose/output/", run_date, "/")
  dir.create(vax1_dose_out_dir, recursive = T, showWarnings=F)
  save(cell_pred, file = paste0(vax1_dose_out_dir, vaccine, "1_dose_cell_draws_eb_bin0_", region, "_", holdout, ".RData"))
  
  rm(cell_pred)

  ### vax0_dose #############################################

  message(paste0("  ", vaccine, "0_dose..."))

  # The math: vax0_dose = vax1_cond(vax2_cond-1)(vax3_cond-1)
  cell_pred <- -1*(cpred_1-1)*(cpred_2 - 1)*(cpred_3 -1)

  # Save output in standard directory
  vax0_dose_out_dir <- paste0(ig_dir, vaccine, "0_dose/output/", run_date, "/")
  dir.create(vax0_dose_out_dir, recursive = T, showWarnings=F)
  save(cell_pred, file = paste0(vax0_dose_out_dir, vaccine, "0_dose_cell_draws_eb_bin0_", region, "_", holdout, ".RData"))
  
  rm(cell_pred)

}

#############################################################
## Now calculate the coverage metrics (P>= doses) ###########
#############################################################

if (doses == 3) {

  ### vax1_cov ##############################################
  
  message(paste0("  ", vaccine, "1_cov..."))

  # The math: vax1_cov = vax1_dose + vax2_dose + vax3_dose
  # vax3_dose = vax3_cond by definition

  cell_pred <- ((cpred_1)*(cpred_2 - 1)*(cpred_3 - 1)) + 
               (cpred_2 - (cpred_2*cpred_3)) +
               (cpred_3)

  vax1_cov_out_dir <- paste0(ig_dir, vaccine, "1_cov/output/", run_date, "/")
  dir.create(vax1_cov_out_dir, recursive = T, showWarnings=F)
  save(cell_pred, file = paste0(vax1_cov_out_dir, vaccine, "1_cov_cell_draws_eb_bin0_", region, "_", holdout, ".RData"))
  # if (holdout == 0) save_mean_raster(cell_pred, paste0(vaccine, "1_cov"), region)
  rm(cell_pred)

  ### vax2_cov ##############################################
  
  if (save_2_cov) {
    message(paste0("  ", vaccine, "2_cov..."))

    # The math: vax2_cov = vax2_dose + vax3_dose
    # vax3_dose = vax3_cond by definition

    cell_pred <- (cpred_2 - (cpred_2*cpred_3)) +
                 (cpred_3)

    vax2_cov_out_dir <- paste0(ig_dir, vaccine, "2_cov/output/", run_date, "/")
    dir.create(vax2_cov_out_dir, recursive = T, showWarnings=F)
    save(cell_pred, file = paste0(vax2_cov_out_dir, vaccine, "2_cov_cell_draws_eb_bin0_", region, "_", holdout, ".RData"))
    # if (holdout == 0) save_mean_raster(cell_pred, paste0(vaccine, "2_cov"), region)
    rm(cell_pred)
  }

}

#############################################################
## Now calculate dropout ####################################
#############################################################

if (doses == 3) {

  ### vax1_3_abs_dropout ####################################
  message(paste0("  ", vaccine, "1_3_abs_dropout..."))

  # The math: vax1_3_abs_dropout = vax1_dose + vax2_dose

  cell_pred <- ((cpred_1)*(cpred_2 - 1)*(cpred_3 - 1)) + 
               (cpred_2 - (cpred_2*cpred_3))

  vax1_3_abs_dropout_out_dir <- paste0(ig_dir, vaccine, "1_3_abs_dropout/output/", run_date, "/")
  dir.create(vax1_3_abs_dropout_out_dir, recursive = T, showWarnings = F)
  save(cell_pred, file = paste0(vax1_3_abs_dropout_out_dir, vaccine, "1_3_abs_dropout_cell_draws_eb_bin0_", region, "_", holdout, ".RData"))
  # if (holdout == 0) save_mean_raster(cell_pred, paste0(vaccine, "1_3_abs_dropout"), region)
  rm(cell_pred)

  ### vax1_3_rel_dropout #####################################
  message(paste0("  ", vaccine, "1_3_rel_dropout..."))

  # The math: vax1_3_rel_dropout = (vax1_dose + vax2_dose)/(1-vax0_dose)

  cell_pred <- (((cpred_1)*(cpred_2 - 1)*(cpred_3 - 1)) + 
               (cpred_2 - (cpred_2*cpred_3))) /
               (1-(-1*(cpred_1-1)*(cpred_2 - 1)*(cpred_3 -1)))

  vax1_3_rel_dropout_out_dir <- paste0(ig_dir, vaccine, "1_3_rel_dropout/output/", run_date, "/")
  dir.create(vax1_3_rel_dropout_out_dir, recursive = T, showWarnings = F)
  save(cell_pred, file = paste0(vax1_3_rel_dropout_out_dir, vaccine, "1_3_rel_dropout_cell_draws_eb_bin0_", region, "_", holdout, ".RData"))
  # if (holdout == 0) save_mean_raster(cell_pred, paste0(vaccine, "1_3_rel_dropout"), region)
  rm(cell_pred)

}

message("\nCombine script finished.")

## END OF FILE
###############################################################################