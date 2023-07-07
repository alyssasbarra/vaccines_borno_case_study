

## clear environment
rm(list=ls())

## Set repo location and indicator group
user               <- 'USERNAME'
core_repo          <- sprintf('FILEPATH')
indic_repo         <- sprintf('FILEPATH')
remote             <- 'origin'
branch             <- 'develop'
pullgit            <- FALSE

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

load_from_parallelize()


if(coverage == TRUE){
  vaccine <- vacc
  set_up_indicators(stem = vaccine, 
                    doses = doses, 
                    single_dose = T,
                    save_doses = F,
                    save_2_cov = F)
  # set_up_indicators(stem = "mcv", doses = 1, single_dose = T) # just mcv1_cov
  # set_up_indicators(stem = "dpt", doses = 3, single_dose = T) # just dpt3_cov
}



## Read config file (from sharedir) and save all parameters in memory
#config <- load_config(repo            = indic_repo,
#                      indicator_group = indicator_group,
#                      indicator       = indicator,
#                      post_est_only   = TRUE,   
#                      run_tests       =FALSE,        
#                      run_date        = run_date)


# parse formatting
#if (class(Regions) == "character" & length(Regions) == 1) Regions <- eval(parse(text=Regions))
if (class(summstats) == "character" & length(summstats) == 1) summstats <- eval(parse(text=summstats))
if (class(year_list) == "character") year_list <- eval(parse(text=year_list))

raking_shapefile_version <- modeling_shapefile_version <-shapefile_version


# Define a log directory and clean out any files that haven't been touched in the last week
log_dir <- paste0("/FILEPATH")
dir.create(log_dir, recursive = T, showWarnings = F)
system(paste0("find FILEPATH -type f -atime +7 -delete"), intern=T)



load(paste0('/FILEPATH'))

main_dir <- paste0('FILEPATH')



  # setting a reference for the number of draws
  ndraws <- ncol(cell_pred)
  overs <- paste0("V", 1:ndraws)
  
  # setup the output folder
  output_dir <- main_dir
  
  # setup values to test against
  cell_pred_dims <- dim(cell_pred)

  reg <- region
  
  
  '%!in%' <- function(x,y)!('%in%'(x,y))
  
  ##### Prep input data into raking:
  if(indicator %!in% c("hib3_cov", "hepb3_cov", "pcv3_cov", "rotac_cov", "mcv2_cov")){
    setwd(paste0("FILEPATH"))
    load(paste0(run_date, "_bin0_", reg, "_0.RData"))
  }
  
  simple_polygon <- simple_raster
  new_simple_polygon <- simple_raster
  
  pixel_id <- pixel_id <- seegSDM:::notMissingIdx(simple_raster)
  
  #####################################################################
  # collect and load the population data from the WorldPop rasters
  covdt <- load_populations_cov(reg, pop_measure, measure = 'count', simple_polygon, simple_raster, year_list, interval_mo, pixel_id = pixel_id)
  
  #####################################################################
  #load the cell id to admin units link
  link_table <- get_link_table(simple_raster, shapefile_version = modeling_shapefile_version)
  
  
  
  #####################################################################
  # Prepping the cell_pred and link table to be linked by making sure they have the appropriate identifiers.  Also performs a 
  # zippering at the region boundary where cells that have a portion of their area outside of the modeling region are reintegrated
  # as a whole cell and considered to be only in one region.  This works becasue a given cell is only modeled in one region.
  link <- prep_link_table(link_table = link_table,
                          simple_raster = simple_raster,
                          pixel_id = pixel_id)
  
  cell_ids <- link_table[[2]]
  
  # getting the connector for sub-national or national raking, This connector gets the IHME location_code for our
  # gbd targets and connects that to the ADM0_CODE or ADM1_CODE as nessecary 
  
  connector <- get_gbd_locs(rake_subnational = TRUE,
                            reg = reg,
                            shapefile_version = modeling_shapefile_version)
  
  # getting the connector for sub-national raking - used to implement countries_not_to_subnat_rake
  nat_connector <- get_gbd_locs(rake_subnational = F,
                                reg = reg,
                                shapefile_version = modeling_shapefile_version)
  
  # merge the connectors on to the link table
  link <- sub_nat_link_merge(rake_subnational=T,
                             link,
                             connector,
                             nat_connector,
                             countries_not_to_subnat_rake)
  
  # set cell pred as a data table, and rename things
  cell_pred <- prep_cell_pred(cell_pred = cell_pred,
                              cell_ids  = cell_ids,
                              pixel_id  = pixel_id,
                              covdt     = covdt)
  
  # merge cell_pred on the link
  cell_pred = merge(link, cell_pred, by.x = 'ID', by.y = 'cell_id',allow.cartesian = TRUE)
  
  # space
  #link <- NULL
  
  ## Raking Population ###################################################################
  # This is done to ensure that the total pop in each raking geography is the same as GBD
  message("raking population")
  
  #convert to fractional population 
  cell_pred = cell_pred[,pop := pop * area_fraction] 
  
  #NA out population where the pixel value is NA (to prevent weirdness with denominators)
  cell_pred = cell_pred[is.na(V1), pop := NA]
  
  
  ##### Using fractional raking #####
  scalars <- fread(paste0(main_dir, indicator, "_", reg, "_pop_rf.csv"))
  
  cell_pred <- merge(cell_pred, scalars, by = c("location_id", "year"))
  
  # rake fractional populations
  cell_pred$pop_raked <- 0
  cell_pred = cell_pred[,pop := pop * pop_scalar] #!# This was formerly creating pop_raked instead of rewriting pop, but this conflicts with the summarize admins function (which calls 'pop' instead)
  
  ##Unraked Counts#########################################
  # creating unraked counts aggregations 
  message("creating a unraked counts aggregations")
  cell_pred$V1 <- cell_pred$V1.x
  #convert cell_pred to counts
  cell_pred <- cell_pred[, (overs) := lapply(.SD, function(x) x * pop), .SDcols=overs]
  
  ## Create unraked counts agregation
  admin_2 <- cell_pred[, lapply(.SD, sum, na.rm=T), .SDcols=c(overs, 'pop'), by = .(year, ADM2_CODE)]
  admin_1 <- cell_pred[, lapply(.SD, sum, na.rm=T), .SDcols=c(overs, 'pop'), by = .(year, ADM1_CODE)]
  admin_0 <- cell_pred[, lapply(.SD, sum, na.rm=T), .SDcols=c(overs, 'pop'), by = .(year, ADM0_CODE)]
  
  sp_hierarchy_list <- 
    link[ADM0_CODE %in% unique(admin_0[, ADM0_CODE]), 
         .(ADM0_CODE, ADM1_CODE, ADM2_CODE, ADM0_NAME, ADM1_NAME, ADM2_NAME)] %>% 
    unique %>% 
    .[, region := reg]
  
  ## save unraked counts aggregations
    paste0(output_dir, "/", indicator, "_unraked_", "_c_admin_draws_eb_bin", 
           "0_", reg, "_0.RData") %>%  
      save(admin_0, admin_1, admin_2, sp_hierarchy_list, file = .)

  
  ##Raked Rates Aggregation#########################################
  # creating unraked rates aggregations
  message("creating a unraked rates aggregations")
  
  # convert back to rates
  admin_0 <- admin_0[, (overs) := lapply(.SD, function(x) x / pop), .SDcols=overs]
  admin_1 <- admin_1[, (overs) := lapply(.SD, function(x) x / pop), .SDcols=overs]
  admin_2 <- admin_2[, (overs) := lapply(.SD, function(x) x / pop), .SDcols=overs]
  
  ## save unraked rates aggregations
    paste0(output_dir, "/", indicator, "_unraked_admin_draws_eb_bin0", 
           "_", reg, "_0.RData") %>%
      save(admin_0, admin_1, admin_2, sp_hierarchy_list, file = .)

  ######## END