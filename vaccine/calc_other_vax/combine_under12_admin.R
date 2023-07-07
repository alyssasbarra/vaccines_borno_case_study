

## clear environment
rm(list=ls())

## Set repo location and indicator group
user               <- 'USERNAME'
core_repo          <- sprintf('FILEPATH')
indic_repo         <- sprintf('FILEPATH')

## sort some directory stuff
commondir      <- sprintf('FILEPATH')
package_list <- c(t(read.csv(sprintf('FILEPATH'),header=FALSE)))

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0('FILEPATH'))
mbg_setup(package_list = package_list, repos = core_repo)

# Custom load indicator-specific functions
source(paste0('FILEPATH'))

## Script-specific stuff begins here ##########################################
indicator_group <- 'vaccine'
indicator <- 'dpt3_cov_under12'

# Load objects from qsub
load_from_parallelize() # run_date, indicator, vaccine, indicator_group, reg
run_date_ratio <- run_date
raking_shapefile_version <- shapefile_version
reg <- region 

raking_shapefile_version <- modeling_shapefile_version <- shapefile_version

# Define directories
dpt3_dir <- paste0("FILEPATH")
ratio_dir <- paste0("FILEPATH")
under12_dir <- paste0("FILEPATH")
dir.create(under12_dir)
sharedir <- under12_dir
## LOAD CELL PREDS ####################




load(paste0('FILEPATH'))

ratio_ad0 <- copy(admin_0)
ratio_ad1 <- copy(admin_1)
ratio_ad2 <- copy(admin_2)


# Now, load in the cell_preds for the ratios
load(paste0(dpt3_dir, 'dpt3_cov_raked_admin_draws_eb_bin0_', reg, '_0.RData'))

dpt3_ad0 <- copy(admin_0)
dpt3_ad1 <- copy(admin_1)
dpt3_ad2 <- copy(admin_2)


ndraws <- ncol(dpt3_ad0) -3
overs <- paste0("V", 1:ndraws)

### getting in same order?
dpt3_ad0 <- dpt3_ad0[order(year),] 
dpt3_ad0 <- dpt3_ad0[order(ADM0_CODE),] 
ad0_keeps <- dpt3_ad0[,!..overs]
dpt3_ad0 <- dpt3_ad0[,..overs]

dpt3_ad1 <- dpt3_ad1[order(year),] 
dpt3_ad1 <- dpt3_ad1[order(ADM1_CODE),] 
ad1_keeps <- dpt3_ad1[,!..overs]
dpt3_ad1 <- dpt3_ad1[,..overs]

dpt3_ad2 <- dpt3_ad2[order(year),] 
dpt3_ad2 <- dpt3_ad2[order(ADM2_CODE),] 
ad2_keeps <- dpt3_ad2[,!..overs]
dpt3_ad2 <- dpt3_ad2[,..overs]

ratio_ad0 <- ratio_ad0[order(year),] 
ratio_ad0 <- ratio_ad0[order(ADM0_CODE),] 
ratio_ad0 <- ratio_ad0[,..overs]

ratio_ad1 <- ratio_ad1[order(year),] 
ratio_ad1 <- ratio_ad1[order(ADM1_CODE),] 
ratio_ad1 <- ratio_ad1[,..overs]

ratio_ad2 <- ratio_ad2[order(year),] 
ratio_ad2 <- ratio_ad2[order(ADM2_CODE),] 
ratio_ad2 <- ratio_ad2[,..overs]



## COMBINE CELL PREDS ######################################################################################

cell_pred_ad0 <- ratio_ad0 * dpt3_ad0
cell_pred_ad1 <- ratio_ad1 * dpt3_ad1
cell_pred_ad2 <- ratio_ad2 * dpt3_ad2

admin_0 <- cbind(ad0_keeps[,c(1,2)], cell_pred_ad0, ad0_keeps[,c(3)])
admin_1 <- cbind(ad1_keeps[,c(1,2)], cell_pred_ad1, ad1_keeps[,c(3)])
admin_2 <- cbind(ad2_keeps[,c(1,2)], cell_pred_ad2, ad2_keeps[,c(3)])



## SAVE CELL PRED IN NEW RUN_DATE FOLDER ######################################################################################

## save raked cell preds
save(admin_0, admin_1, admin_2, sp_hierarchy_list, file = paste0(
  under12_dir, "dpt3_cov_under12_raked_admin_draws_eb_bin0_", reg, "_0.RData"
))




