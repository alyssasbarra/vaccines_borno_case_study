fractional_rake_rates_pre_2000 <- function(cell_pred = cell_pred,
                                  simple_raster = simple_raster,
                                  simple_polygon = simple_polygon,
                                  pixel_id = pixel_id,
                                  shapefile_version = shapefile_version,
                                  reg = reg,
                                  pop_measure = pop_measure,
                                  year_list = year_list,
                                  use_intermediate_years = TRUE,
                                  interval_mo = interval_mo,
                                  rake_subnational = rake_subnational,
                                  age_group = age_group,
                                  sex_id = sex_id,
                                  sharedir = sharedir,
                                  run_date = run_date,
                                  indicator = indicator,
                                  gbd = gbd,
                                  rake_method = "linear",
                                  gbd_pops = gbd_pops,
                                  countries_not_to_rake = NULL,
                                  countries_not_to_subnat_rake = NULL,
                                  custom_output_folder = NULL){
  # setting a reference for the number of draws
  ndraws = ncol(cell_pred)
  
  #####################################################################
  # collect and load the population data from the WorldPop rasters
  covdt <- load_populations_cov_vax_custom(reg, vax_age=NULL, pop_measure=pop_measure, measure = 'count', simple_polygon, simple_raster, year_list, interval_mo, pixel_id = pixel_id, age_is_gbd_age = FALSE)
  
  if(!use_intermediate_years) {
    print("Subsetting to supplied years only")
    covdt <- covdt[year %in% year_list]
  }
  
  
  #####################################################################
  #load the cell id to admin units link
  link_table <- get_link_table(simple_raster, shapefile_version = shapefile_version)
  
  
  
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
  
  connector <- get_gbd_locs(rake_subnational = rake_subnational,
                            reg = reg,
                            shapefile_version = shapefile_version)
  
  # getting the connector for sub-national raking - used to implement countries_not_to_subnat_rake
  nat_connector <- get_gbd_locs(rake_subnational = F,
                                reg = reg,
                                shapefile_version = shapefile_version)
  
  # merge the connectors on to the link table
  link <- sub_nat_link_merge(rake_subnational,
                             link,
                             connector,
                             nat_connector,
                             countries_not_to_subnat_rake)
  
  
  #set cell pred as a data table, and rename things
  if (!inherits(x = cell_pred, 'data.table')) {
    cell_pred = as.data.table(cell_pred)
    names(cell_pred) = paste0('V',1:ncol(cell_pred))
  }
  
  ####################################
  #Link cell pred 2000 and 1980-1999 population by hand
  cell_pred[, cell_pred_id := .I] #cell_pred object ID
  cell_pred[,cell_id := rep(cell_ids, times = nrow(cell_pred) / length(cell_ids))]  #cell id references the africa map
  cell_pred[,pixel_id := rep(pixel_id, times = nrow(cell_pred) / length(pixel_id))] #pixel id references the regional map
  cell_pred[,year_id := rep(1:20, each = length(cell_ids))] #pixel id references the regional map
  
  cell_pred <- cell_pred[year_id==1,]
  cell_pred[,year_id:=NULL]

  #Copy 2000 cell pred data over for each pre-2000 year
  cell_pred <- do.call("rbind", replicate(length(year_list), cell_pred, simplify = FALSE))

  #add population, year and potentially the stackers
  cell_pred = cbind(cell_pred, covdt[,c('year', 'pop'),with = F])
  
  
  ############################################################################################################
  # merge cell_pred on the link
  cell_pred = merge(link, cell_pred, by.x = 'ID', by.y = 'cell_id',allow.cartesian = TRUE)
  
  # space
  link <- NULL
  
  ## Raking Population ###################################################################
  # This is done to ensure that the total pop in each raking geography is the same as GBD
  message("raking population")
  
  #convert to fractional population 
  cell_pred = cell_pred[,pop := pop * area_fraction] 
  
  #NA out population where the pixel value is NA (to prevent weirdness with denominators)
  cell_pred = cell_pred[is.na(V1), pop := NA]
  
  scalars <- calculate_pop_scalars(cell_pred = cell_pred,
                                   age_group = age_group,
                                   connector = connector,
                                   sex       = sex_id,
                                   sharedir  = sharedir,
                                   run_date  = run_date,
                                   indicator = indicator,
                                   stratum   = paste0('bin', age_bin, '_', reg),
                                   gbd_pops  = gbd_pops,
                                   custom_output_folder = custom_output_folder)
  # add back to the cell_pred as a population rf
  cell_pred <- merge(cell_pred, scalars, by = c("location_id", "year"))
  
  # rake fractional populations
  cell_pred$pop_raked <- 0
  cell_pred = cell_pred[,pop_raked := pop * pop_scalar]
  cell_pred$pop <- NULL
  
  ## Raking actual data ###################################################################
  message("raking rates")
  # Calculate Fractional Raking Factors
  fractional_rf <- calculate_fractional_rfs(ndraws    = ndraws,
                                            cell_pred = cell_pred,
                                            gbd       = gbd,
                                            sharedir  = sharedir,
                                            run_date  = run_date,
                                            indicator = indicator,
                                            shapefile_version = shapefile_version,
                                            stratum   = paste0('bin', age_group, '_', reg),
                                            countries_not_to_rake = countries_not_to_rake,
                                            custom_output_folder = custom_output_folder,
                                            countries_not_to_subnat_rake = countries_not_to_subnat_rake,
                                            rake_method = rake_method)
}
