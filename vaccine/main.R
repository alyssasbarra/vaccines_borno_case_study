###############################################################################
###############################################################################
## MBG main Launch Script
##
## Author: USERNAME
## Vaccine: DPT
## Date: DATE
##
## This is the main launch script.  Will run individual models for each needed
## vaccine or conditional vaccine data set for ordinal regression using a
## continuation-ratio model, then combine these arithmatically in order to
## the desired vaccine coverage metrics. 
##
## Source:
##   source("FILEPATH")
## Prod cluster as qsub:
##   qsub -e FILEPATH -o FILEPATH -cwd -pe multi_slot 16 -P proj_geospatial -N main_dpt FILEPATH FILEPATH
## Geo cluster as qsub:
##   qsub -e FILEPATH -o FILEPATH -N main_dpt -P proj_geospatial -l gn=TRUE FILEPATH FILEPATH
###############################################################################
###############################################################################

###############################################################################
## SETUP
###############################################################################

## clear environment
rm(list=ls())

## Set repo location and indicator group
user               <- 'USERNAME'
core_repo          <- sprintf('FILEPATH',user)
indic_repo         <- sprintf('FILEPATH',user)
remote             <- 'origin'
branch             <- 'develop'
pullgit            <- FALSE

## sort some directory stuff
commondir      <- sprintf('FILEPATH')
package_list <- c(t(read.csv(sprintf('FILEPATH',commondir),header=FALSE)))

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0('FILEPATH'))
mbg_setup(package_list = package_list, repos = core_repo)

# Custom load indicator-specific functions
source(paste0('FILEPATH'))

'%!in%' <- function(x,y)!('%in%'(x,y))

## Script-specific stuff begins here ##########################################
indicator_group <- "vaccine"
coverage <- FALSE # if false, modeling ratio! i.e. the timeliness ratio
vacc <- 'timeliness'

if(coverage == TRUE){
  vaccine <- vacc
  set_up_indicators(stem = vaccine, 
                  doses = 3, # if modeling mcv, this should be 1
                  single_dose = F, # if modeling mcv, this should be TRUE
                  save_doses = F,
                  save_2_cov = F)
}


if(vacc == "timeliness"){  
  vaccine <- "dpt3_timeliness_ratio"
  model_indicators <- vaccine
  dose_indicators <- vaccine
  rake_indicators <- vaccine
  postest_indicators <- vaccine
  all_indicators <- vaccine
  doses <- ''
}

slots              <- 4

## Create run date in correct format
#run_date <- make_time_stamp(TRUE)
run_date <- 'RUN_DATE'

# Define a log directory and clean out any files that haven't been touched in the last week
log_dir <- paste0("FILEPATH")
dir.create(log_dir, recursive = T)
system(paste0("find FILEPATH", user, "/logs -type f -atime +7 -delete"), intern=T)

if(vaccine %in% c("dpt","polio", "bcg", "mcv")){
  config <- load_config(repo            = indic_repo,
                    indicator_group = "",
                    indicator       = NULL,
                    config_name     = paste0("config_", vaccine), 
                    covs_name       = paste0(vaccine, "_covs"),
                    run_test        = FALSE) 
}
if(vaccine == "dpt3_timeliness_ratio"){
  config <- load_config(repo            = indic_repo,
                      indicator_group = "",
                      indicator       = NULL,
                      config_name     = paste0("config_dpt"), 
                      covs_name       = paste0("dpt_covs"),
                      run_test        = FALSE)
  # config[V1=='gbd_date', V2:='DATE']
  # gbd_date <- 'DATE'
  dpt3_run_date <- 'RUN_DATE'
}
#if(vaccine == "mcv"){
#  config <- load_config(repo            = indic_repo,
#                      indicator_group = "",
##                      indicator       = NULL,
#                      config_name     = paste0("config_", vaccine), 
#                      covs_name       = paste0(vaccine, doses, "_covs"),
#                      run_test        = FALSE)
#}

## Ensure you have defined all necessary settings in your config
check_config()

# distribute config to all indicator directories
distribute_config(cfg = config, indicators = all_indicators)

# parse formatting
if (class(Regions) == "character" & length(Regions) == 1) Regions <- eval(parse(text=Regions))
if (class(summstats) == "character" & length(summstats) == 1) summstats <- eval(parse(text=summstats))
if (class(year_list) == "character") year_list <- eval(parse(text=year_list))

###############################################################################
## Make Holdouts
###############################################################################
if(as.logical(makeholdouts)){

  if(as.logical(load_other_holdouts)){
    message(paste0("Loading holdouts from ", holdout_rundate, "..."))

    stratum_filename <- paste0("FILEPATH") 
    
    stratum_ho <- readRDS(stratum_filename)

    saveRDS(stratum_ho,
            file = paste0("FILEPATH"))

  } else {
    # load the full input data for last dose of vaccine in series
    df <- load_input_data(indicator   = paste0(vaccine, doses, "_cov"),
                          simple      = NULL,
                          removeyemen = TRUE,
                          years       = yearload,
                          withtag     = as.logical(withtag),
                          datatag     = datatag,
                          use_share   = as.logical(use_share),
                          yl          = year_list)

    # add in location information
    df <- merge_with_ihme_loc(df, shapefile_version = modeling_shapefile_version)

    indicator <- paste0(vaccine, doses, "_cov") # referenced internally in make_folds

  # Load simple polygon template by region
    gaul_list           <- get_adm0_codes(Regions, shapefile_version = modeling_shapefile_version)
    simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 1, tolerance = 0.4, shapefile_version = modeling_shapefile_version)
    subset_shape        <- simple_polygon_list[[1]]
    # Load simple raster by region
    raster_list        <- build_simple_raster_pop(subset_shape)
    simple_raster      <- raster_list[['simple_raster']]

    # load admin 2 shapefile and crop by region
    shapefile_admin <- shapefile(get_admin_shapefile(admin_level = 1, version = modeling_shapefile_version))
    shapefile_admin <- shapefile_admin[which(shapefile_admin$ADM0_CODE %in% gaul_list),]
    # load mask raster and crop by region
    raster_mask <- raster('FILEPATH')
    raster_mask <- crop(raster_mask, extent(simple_raster))
    # make a list of dfs for each region, with 5 qt folds identified in each    
    lat_col = 'latitude'
    long_col = 'longitude';
    ss_col = 'N'
    yr_col = 'year'
    admin_raster <- simple_raster


    stratum_ho <- make_folds(data = df, n_folds = 5, spat_strat = 'poly',
                              temp_strat = NULL, strat_cols = 'region',
                              #admin_raster=simple_raster,
                              shape_ident="ADM1_CODE",
                              ts         = as.numeric(ho_ts),
                              mb         = as.numeric(ho_mb),
                              admin_shps=shapefile_admin,
                              admin_raster=simple_raster,
                              mask_shape=subset_shape,
                              mask_raster=raster_mask,
                              lat_col = 'latitude', long_col = 'longitude',
                              ss_col = 'N', yr_col = 'year', seed = 98112)
    # Recreate & distribute holdouts
  } 

  indicators <- indicator
  invisible(lapply((unique(c(indicators, postest_indicators))), function(ind) {

    sharedir <- sprintf('FILEPATH')

    ## Create holdouts if not already present
    if(as.logical(makeholdouts) & !(file.exists(paste0(sharedir,'/output/',run_date,"/stratum.rds")))){

      message(paste0("Recreating holdouts for ", ind, "..."))
      # load the full input data
      df <- load_input_data(indicator   = ind,
                            simple      = NULL,
                            removeyemen = TRUE,
                            years       = yearload,
                            withtag     = as.logical(withtag),
                            datatag     = datatag,
                            use_share   = as.logical(use_share),
                            yl          = year_list)

      # add in location information
      df <- merge_with_ihme_loc(df, shapefile_version = modeling_shapefile_version)

      # make a list of dfs for each region, with 5 qt folds identified in each

      # New way: use the loaded data and the dpt3_cov holdouts to assign each
      #          row in the data to the same data across modeled indicators
      stratum_ho <- recreate_holdouts(data = df,
                                      row_id_col = "row_id",
                                      load_from_indic = paste0(vaccine, "1_cov"),
                                      rd = run_date,
                                      ig = indicator_group)

      # Save stratum_ho object for reference
      save_stratum_ho(indic       = ind,
                      ig          = indicator_group,
                      rd          = run_date,
                      stratum_obj = stratum_ho)

    } # if(as.logical(makeholdouts) & !(file.exists(...
  })) # invisible(lapply...
} #if(as.logical(makeholdouts))

###############################################################################
## Launch parallel models
###############################################################################

# launch_lv = list(indicator = model_indicators,
#                  use_gn = TRUE)
# 
# launch_para_output <- parallelize(script = "launch",
#                                   log_location = paste0(log_dir, "launch/"),
#                                   expand_vars = launch_lv,
#                                   save_objs = c("run_date", "indicator_group", "vaccine", "core_repo", "shapefile_version"),
#                                   prefix = "launch",
#                                   slots = 1,
#                                   memory = 10,
#                                   script_dir = indic_repo,
#                                   geo_nodes = FALSE,
#                                   singularity = 'default',
#                                   queue             = 'long.q',
#                                   proj              = 'proj_geospatial')
# stop('stop here for now')
# 
# monitor_jobs(launch_para_output, notification = "pushover")


###############################################################################
## Rake indicators
###############################################################################
 # pushover_notify("main script: starting raking", 
  #                title = "Raking")

multiple_raking_indicators <-rake_indicators  #Dealing with multiple indicators w DPT
#   
# for(r in multiple_raking_indicators){
#     rake_indicators <- r
# 
# 
#   rake_lv <- list(region = Regions,
#                   holdout = 0)
#   
# 
#   rake_para_output <- parallelize(script = "rake_fractional",
#                                   log_location = paste0("FILEPATH"),
#                                   expand_vars = rake_lv,
#                                   save_objs = c("indicator_group", "run_date", "vaccine", 
#                                                 "rake_indicators", "doses", "shapefile_version"),
#                                   prefix = "rake_fractional",
#                                   slots = 8, # geo_nodes
#                                   threads = 4,
#                                   use_c2_nodes = FALSE,
#                                   script_dir = indic_repo, 
#                                   memory = 300,
#                                   geo_nodes = F,
#                                   singularity = 'default')
#   
#  #if(length(multiple_raking_indicators)==1) monitor_jobs(rake_para_output, notification = "pushover")
# }
# 
# stop('stop here for now')

####################################################################################
#Unraked aggregation
####################################################################################

### Run in parallel #######################

# for(r in multiple_raking_indicators){
#   indicator <- r
# expand_lv <- list(region = Regions)
# # Define a log directory and clean out any files that haven't been touched in the last week
# #log_dir <- paste0("FILEPATH")
# dir.create(log_dir, recursive = T, showWarnings = F)
# expand_para_output <- parallelize(script = "unraked_aggregate_cell_draws",
#                                   log_location = paste0("FILEPATH"),
#                                   expand_vars = expand_lv,
#                                   save_objs = c("Regions", "indicator_group",
#                                                 "coverage", "vacc", "doses", "indicator",
#                                                 "run_date","summstats","shapefile_version","year_list",
#                                                 "pop_measure","interval_mo","pop_release","countries_not_to_subnat_rake"),
#                                   prefix = "unraked_agg",
#                                   slots = 5,
#                                   threads = 4,
#                                   script_dir = indic_repo,
#                                   memory = 200,
#                                   geo_nodes = FALSE,
#                                   use_c2_nodes = FALSE,
#                                   singularity = "default")
# 
# 
# }
# #stop('stop here for now')
# monitor_jobs(expand_para_output, notification = "pushover")
####################################################################################
#processing of conditional DTP indicators (DTP1 & DTP3<12)
####################################################################################

### Run in parallel #######################

# if(vaccine == 'dpt'){
# 
#   # Define a log directory and clean out any files that haven't been touched in the last week
#   log_dir <- paste0("FILEPATH")
#   dir.create(log_dir, recursive = T, showWarnings = F)
#   expand_para_output <- parallelize(script = "post_hoc_combine_dpt1",
#                                     log_location = log_dir,
#                                     expand_vars = run_date,         
#                                     save_objs = c("Regions", "indicator_group",
#                                                   "coverage", "vacc", "doses",  
#                                                   "run_date", "shapefile_version", "summstats"),
#                                     prefix = "combine_dpt1",
#                                     slots = 5,
#                                     threads = 8,
#                                     script_dir = indic_repo,
#                                     memory = 100,
#                                     geo_nodes = FALSE,
#                                     use_c2_nodes = FALSE,
#                                     singularity = "default")
#   monitor_jobs(expand_para_output, notification = "pushover")
# }

# if(vaccine == "dpt3_timeliness_ratio"){
#   expand_lv <- list(region = Regions)
#   # Define a log directory and clean out any files that haven't been touched in the last week
#   log_dir <- paste0("FILEPATH")
#   dir.create(log_dir, recursive = T, showWarnings = F)
#   
#   
#   expand_para_output <- parallelize(script = "combine_under12_cells_rasters",
#                                     log_location = paste0("FILEPATH"),
#                                     expand_vars = expand_lv,  
#                                     save_objs = c("Regions", "indicator_group",
#                                                   "coverage", "vacc", "doses",  
#                                                   "run_date", "dpt3_run_date", "shapefile_version", "summstats"),
#                                     prefix = "combine_dpt_under12",
#                                     slots = 5,
#                                     threads = 8,
#                                     script_dir = paste0(indic_repo, "/calc_other_vax/"),
#                                     memory = 400,
#                                     geo_nodes = FALSE,
#                                     use_c2_nodes = FALSE,
#                                     singularity = "default")
#   
#   expand_para_output <- parallelize(script = "combine_under12_admin",
#                                     log_location = paste0("FILEPATH"),
#                                     expand_vars = expand_lv,  
#                                     save_objs = c("Regions", "indicator_group",
#                                                   "coverage", "vacc", "doses",  
#                                                   "run_date", "dpt3_run_date", "shapefile_version", "summstats"),
#                                     prefix = "combine_dpt_under12",
#                                     slots = 5,
#                                     threads = 8,
#                                     script_dir = paste0(indic_repo, "/calc_other_vax/"),
#                                     memory = 50,
#                                     geo_nodes = FALSE,
#                                     use_c2_nodes = FALSE,
#                                     singularity = "default")
#   
#   monitor_jobs(expand_para_output, notification = "pushover")
# }


###############################################################################
## Post predict (summarize, combine, save files)
###############################################################################

## Convenience
if(vaccine=='dpt') all_indicators <- c(all_indicators, 'dpt1_cov')
if(vaccine=='dpt3_timeliness_ratio') all_indicators <- c(all_indicators, 'dpt3_cov_under12')


for(indicator in all_indicators[2]){
  
  sharedir <- sprintf('FILEPATH')
  
  #Regions <- c("vax_eaas","vax_name", "vax_trsa", "vax_ansa", "vax_caeu", "vax_crbn", "vax_cssa", "vax_ctam", "vax_essa", "vax_seas", "vax_soas", "vax_sssa", "vax_wssa")
  
  
  #Regions <- c("vax_ansa_pcv", "vax_caeu_pcv", "vax_crct_pcv", "vax_cssa_pcv", "vax_eaas_pcv", "vax_essa_pcv", "vax_seas_pcv", "vax_some_pcv", "vax_sssa_pcv", "vax_trsa_pcv", "vax_wssa_pcv")
  
  
  strata <- Regions
  #######################################################################################
  ## Merge!! 
  #######################################################################################
  
  if(indicator %!in% c('dpt1_cov', 'dpt3_cov_under12')) rr <- c("raked", "unraked") else rr <- 'raked'
  if(indicator %!in% c('dpt1_cov', 'dpt3_cov_under12')) rr_tf <- c(T,F) else rr_tf <- T
  if(indicator %!in% c('dpt1_cov', 'dpt3_cov_under12')) rf_tab <- T else rf_tab <- F
  run_summary <- T

if (("p_below" %in% summstats) & !(indicator %in% c(paste0(vaccine, "1_cov"), paste0(vaccine, "3_cov")))) {
  # Only run this for 1st and 3rd doses
  summstats <- summstats[summstats != "p_below"]
}

post_load_combine_save(regions    = Regions,
                       summstats  = summstats,
                       raked      = rr,
                       rf_table   = rf_tab,
                       run_summ   = run_summary,
                       indic      = indicator,
                       ig         = indicator_group,
                       sdir       = sharedir)

# Clean up / delete unnecessary files
clean_after_postest(indicator             = indicator,
                    indicator_group       = indicator_group,
                    run_date              = run_date,
                    strata                = strata,
                    delete_region_rasters = F)


#######################################################################################
## Combine aggreagations
#######################################################################################
holdouts <- 0
combine_aggregation(rd = run_date, indic = indicator, ig = indicator_group,
                    ages = 0, 
                    regions = strata,
                    holdouts = holdouts,
                    raked = rr_tf,
                    delete_region_files = F)


#######################################################################################
## Summarize
#######################################################################################


summarize_admins(summstats = c("mean", "lower", "upper", "cirange", "cfb"), 
                 ad_levels = c(0,1,2), 
                 raked = rr_tf)

message("Summarize p_above")
summarize_admins(indicator, indicator_group,
                 summstats = c("p_above"),
                 raked = rr_tf,
                 ad_levels = c(0,1,2),
                 file_addin = "p_0.8_or_better",
                 value = 0.8,
                 equal_to = T)
}


stop('stop here until mapping')
########################################################################################
##Save nicely for mapping
########################################################################################
postest_indicators <- c('dpt3_cov_under12') #indicator
#summstats <- summstats[-5]

tifs_to_copy <- data.table(expand.grid(summstats, c(T), postest_indicators, stringsAsFactors = F))
names(tifs_to_copy) <- c("measure", "raked", "indicator") 

for (i in 1:nrow(tifs_to_copy)) {
  copy_tif_to_map_input_dir(ind = tifs_to_copy[i, "indicator"],
                            ig = indicator_group,
                            measure = tifs_to_copy[i, "measure"],
                            rd = run_date,
                            raked = tifs_to_copy[i, "raked"],
                            yl = year_list)
}

# Save ads
admin_lvs <- data.table(expand.grid(summstats, c(T), c(0,1,2), postest_indicators, stringsAsFactors = F))
#admin_lvs <- rbind(admin_lvs,
#                   data.table(expand.grid("psup80", 
#                                          c(T,F), 
#                                          c(0,1,2),
#                                          postest_indicators[postest_indicators != paste0(vaccine, "1_3_rel_dropout")], 
#                                          stringsAsFactors = F)))

names(admin_lvs) <- c("measure", "raked", "ad_level", "indicator") 

for (i in 1:nrow(admin_lvs)) {
  copy_admins_to_map_input_dir(ind = admin_lvs[i, "indicator"],
                               ig = indicator_group, 
                               measure = as.character(admin_lvs[i, "measure"]), 
                               rd = run_date, 
                               raked = admin_lvs[i, "raked"], 
                               yl = year_list, 
                               ad_level = admin_lvs[i, "ad_level"]) 
}







summstats <- c("mean", "cfb")
summstats <- 'mean'


#### gather list of admin codes to keep in csvs



ad0 <- fread(paste0('FILEPATH'))
ad1 <- fread(paste0('FILEPATH'))
ad2 <- fread(paste0('FILEPATH'))



ad0 <- unique(ad0$ADM0_CODE)
ad1 <- unique(ad1$ADM1_CODE)
ad2 <- unique(ad2$ADM2_CODE)



scope <- raster(paste0('FILEPATH'))


for(indicator in postest_indicators){
  message(indicator)
dir.create(paste0('FILEPATH'))
dir.create(paste0('FILEPATH'))



for(stat in summstats){
  message(stat)
  
  #### first, re-save csv's
  ad0_csv <- fread(paste0('FILEPATH'))
  ad1_csv <- fread(paste0('FILEPATH'))
  ad2_csv <- fread(paste0('FILEPATH'))
  
  ad0_csv <- subset(ad0_csv, ADM0_CODE %in% ad0)
  ad1_csv <- subset(ad1_csv, ADM1_CODE %in% ad1)
  ad2_csv <- subset(ad2_csv, ADM2_CODE %in% ad2)
  
  write.csv(ad0_csv, paste0('FILEPATH'))
  write.csv(ad1_csv, paste0('FILEPATH'))
  write.csv(ad2_csv, paste0('FILEPATH'))
  
  
  
  
  #### next, re-save tifs
  tif_yr <- brick(paste0('FILEPATH'))
  
  for(yr in 1:20){
    message((yr + 1999))
    year <- yr+1999
    
    single_yr <- tif_yr[[yr]]     
    single_yr <- crop(single_yr, scope)
    single_yr <- mask(single_yr, scope)
    
    writeRaster(single_yr, paste0('FILEPATH'))
  }
  
  
}  
}

########################################################################################
message('Done with main script!')
## END OF FILE
###############################################################################