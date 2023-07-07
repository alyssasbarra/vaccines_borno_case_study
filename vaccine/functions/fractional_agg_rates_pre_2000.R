fractional_agg_rates_pre_2000 <- function(cell_pred = cell_pred,
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
                                 sharedir = sharedir,
                                 run_date = run_date,
                                 indicator = indicator,
                                 main_dir = main_dir,
                                 rake_method = "linear",
                                 age = 0,
                                 holdout = 0,
                                 return_objects = FALSE,
                                 countries_not_to_subnat_rake = countries_not_to_subnat_rake,
                                 custom_output_folder = NULL,
                                 stratum=reg) {
  # setting a reference for the number of draws
  ndraws <- ncol(cell_pred)
  region <- reg
  #####################################################################
  # load the cell id to admin units link
  
  link_table <- get_link_table(simple_raster, shapefile_version = shapefile_version)
  
  #####################################################################
  # collect and load the population data
  covdt <- load_populations_cov_vax_custom(reg, vax_age=NULL, pop_measure=pop_measure, measure = 'count', simple_polygon, simple_raster, year_list, interval_mo, pixel_id = pixel_id, age_is_gbd_age = FALSE)
  
  if(!use_intermediate_years) {
    print("Subsetting to supplied years only")
    covdt <- covdt[year %in% year_list]
  }
  
  #####################################################################
  # Prepping the cell_pred and link table to be linked and then merging them
  link <- prep_link_table(
    link_table = link_table,
    simple_raster = simple_raster,
    pixel_id = pixel_id
  )
  
  cell_ids <- link_table[[2]]
  
  # getting the connector for sub-national raking
  connector <- get_gbd_locs(
    rake_subnational = rake_subnational,
    reg = reg,
    shapefile_version = shapefile_version
  )
  
  # getting the connector for sub-national raking
  nat_connector <- get_gbd_locs(
    rake_subnational = F,
    reg = reg,
    shapefile_version = shapefile_version
  )
  
  # merge the connector on to the link table
  link <- sub_nat_link_merge(
    rake_subnational,
    link,
    connector,
    nat_connector,
    countries_not_to_subnat_rake
  )
  
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
  cell_pred[, cell_pred_id := .I] #cell_pred object ID #to help with deduplication later
  
  #add population, year and potentially the stackers
  cell_pred = cbind(cell_pred, covdt[,c('year', 'pop'),with = F])
  
  
  ############################################################################################################
  
  # merge on the link
  cell_pred <- merge(link, cell_pred, all.y=T, by.x = "ID", by.y = "cell_id", allow.cartesian = TRUE)
  
  
  ############################################################
  # adding the raking factors and scaling the populations
  
  message("adding raking factors")
  # convert to fractional population
  cell_pred <- cell_pred[, pop := pop * area_fraction]
  
  scalars <- read.csv(file = ifelse(is.null(custom_output_folder),
                                    paste0(sharedir, "/output/", run_date, "/", indicator, "_", stratum, "_pop_rf.csv"),
                                    paste0(custom_output_folder, "/", indicator, "_", stratum, "_pop_rf.csv")
  ))
  
  # add back to the cell_pred as a population rf
  cell_pred <- merge(cell_pred, scalars, all.x=T, by = c("location_id", "year"))
  
  ## load the raking factors
  fractional_rf <- read.csv(file = ifelse(is.null(custom_output_folder),
                                          paste0(sharedir, "/output/", run_date, "/", indicator, "_", stratum, "_rf.csv"),
                                          paste0(custom_output_folder, "/", indicator, "_", stratum, "_rf.csv")
  ))
  
  ## merge them onto the cell_pred
  cell_pred <- merge(cell_pred, fractional_rf, all.x=T, by = c("location_id", "year"))
  cell_pred$rf <- as.numeric(as.character(cell_pred$rf))
  
  ############################################################
  # creating a raked rates cell_pred (this happens first becasue once we go to counts in the cell_pred we can't do back without loading a fresh copy)
  message("creating a raked rates cell_pred")
  # rake rates
  overs <- paste0("V", 1:ndraws)
  
  if(rake_method == "linear"){
    cell_pred <- cell_pred[, (overs) := lapply(overs, function(x) get(x) * rf)]
  }else{
    cell_pred <- cell_pred[, (overs) := lapply(overs, function(x) invlogit(logit(get(x)) + rf))]
  }
  
  # multiply the cell_pred by the area fraction for the dedupe function (so that each cell will add to 1 and the constituent rates are weighted by area)
  cell_pred <- cell_pred[, (overs) := lapply(overs, function(x) get(x) * area_fraction) ]
  
  raked_cell_pred <- dedupe_linked_cell_pred(cell_pred, overs)
  
  # Save raked rates cell_pred object
  save(raked_cell_pred, file = ifelse(is.null(custom_output_folder),
                                      paste0(
                                        sharedir, "/output/", run_date, "/",
                                        indicator, "_raked_cell_draws_eb_bin", age, "_", reg, "_", holdout, ".RData"
                                      ),
                                      paste0(
                                        custom_output_folder, "/",
                                        indicator, "_raked_cell_draws_eb_bin", age, "_", reg, "_", holdout, ".RData"
                                      )
  ))
  
  
  # SPACE!!!!!
  if (!return_objects) {
    raked_cell_pred <- NULL
    gc()
  }
  
  # un do the area fraction (so that you can use this cell pred again)
  overs <- paste0("V", 1:ndraws)
  cell_pred <- cell_pred[, (overs) := lapply(overs, function(x) get(x) / area_fraction) ]
  
  ############################################################
  # creating a raked counts cell_pred
  message("creating a raked counts cell_pred")
  # rake fractional populations
  cell_pred$pop_raked <- 0
  cell_pred <- cell_pred[, pop_raked := pop * pop_scalar]
  cell_pred$pop <- NULL
  
  message("converting from prevalence to counts")
  # set the variables to aggregate
  overs <- paste0("V", 1:ndraws)
  cell_pred <- cell_pred[, (overs) := lapply(overs, function(x) get(x) * pop_raked)]
  
  # NA out population where the pixel value is NA (to prevent weirdness with denominators)
  cell_pred <- cell_pred[is.na(V1), pop_raked := NA]
  raked_cell_pred_c <- dedupe_linked_cell_pred(cell_pred, overs)
  save(raked_cell_pred_c, file = ifelse(is.null(custom_output_folder),
                                        paste0(
                                          sharedir, "/output/", run_date, "/",
                                          indicator, "_raked_c_cell_draws_eb_bin", age, "_", reg, "_", holdout, ".RData"
                                        ), paste0(
                                          custom_output_folder, "/",
                                          indicator, "_raked_c_cell_draws_eb_bin", age, "_", reg, "_", holdout, ".RData"
                                        )
  ))
  
  # SPACE
  raked_cell_pred_c <- NULL
  
  ############################################################
  # creating a raked counts aggregations
  message("creating a raked counts aggregations")
  
  # do calculations!
  admin_2 <- cell_pred[, lapply(c(overs, "pop_raked"), function(x) sum(get(x), na.rm = T)), by = c("year", "ADM2_CODE", "ADM0_CODE")]
  admin_1 <- cell_pred[, lapply(c(overs, "pop_raked"), function(x) sum(get(x), na.rm = T)), by = c("year", "ADM1_CODE", "ADM0_CODE")]
  admin_0 <- cell_pred[, lapply(c(overs, "pop_raked"), function(x) sum(get(x), na.rm = T)), by = c("year", "ADM0_CODE")]
  
  
  setnames(admin_2, grep("V[0-9]", names(admin_2), value = T), c(overs, "pop_raked"))
  setnames(admin_1, grep("V[0-9]", names(admin_1), value = T), c(overs, "pop_raked"))
  setnames(admin_0, grep("V[0-9]", names(admin_0), value = T), c(overs, "pop_raked"))
  
  # create the spatial hierarchy
  sp_hierarchy_list <- unique(link[ADM0_CODE %in% unique(admin_0[, ADM0_CODE]), .(ADM0_CODE, ADM1_CODE, ADM2_CODE, ADM0_NAME, ADM1_NAME, ADM2_NAME, region)])
  sp_hierarchy_list[, region := reg]
  
  # cleaning raked admin draws in count space
  admin_0$pop <- admin_0$pop_raked
  admin_0$pop_raked <- NULL
  admin_1$ADM0_CODE <- NULL
  admin_1$pop <- admin_1$pop_raked
  admin_1$pop_raked <- NULL
  admin_2$ADM0_CODE <- NULL
  admin_2$pop <- admin_2$pop_raked
  admin_2$pop_raked <- NULL
  
  ## save raked counts aggregations
  save(admin_0, admin_1, admin_2, sp_hierarchy_list,
       file = ifelse(is.null(custom_output_folder),
                     paste0(main_dir, "/", indicator, "_", "raked", "_c", "_admin_draws_eb_bin", age, "_", reg, "_", holdout, ".RData"),
                     paste0(custom_output_folder, "/", indicator, "_", "raked", "_c", "_admin_draws_eb_bin", age, "_", reg, "_", holdout, ".RData")
       )
  )
  
  ############################################################
  # creating raked rates aggregations (you can work back at the admin level because there are no admin's with pop = 0)
  message("creating a raked rates aggregations")
  
  # convert back to rates/prevalence
  overs <- paste0("V", 1:ndraws)
  admin_0 <- admin_0[, (overs) := lapply(overs, function(x) get(x) / pop) ]
  admin_1 <- admin_1[, (overs) := lapply(overs, function(x) get(x) / pop) ]
  admin_2 <- admin_2[, (overs) := lapply(overs, function(x) get(x) / pop) ]
  
  ## save raked rates aggregations
  save(admin_0, admin_1, admin_2, sp_hierarchy_list,
       file = ifelse(is.null(custom_output_folder),
                     paste0(main_dir, "/", indicator, "_", "raked", "_admin_draws_eb_bin", age, "_", reg, "_", holdout, ".RData"),
                     paste0(custom_output_folder, "/", indicator, "_", "raked", "_admin_draws_eb_bin", age, "_", reg, "_", holdout, ".RData")
       )
  )
  
  
  ############################################################
  # creating unraked counts aggregations (can be done two ways, un raking the counts cell_pred or reloading the unraked cell_pred and multiplying by the pop.  I chose this to avoid reloading cell_preds)
  message("creating a unraked counts aggregations")
  
  # unrake all draws
  overs <- paste0("V", 1:ndraws)
  cell_pred <- cell_pred[, (overs) := lapply(overs, function(x) get(x) / rf) ]
  
  ## Create unraked counts agregation
  admin_2 <- cell_pred[, lapply(c(overs, "pop_raked"), function(x) sum(get(x), na.rm = T)), by = c("year", "ADM2_CODE", "ADM0_CODE")]
  admin_1 <- cell_pred[, lapply(c(overs, "pop_raked"), function(x) sum(get(x), na.rm = T)), by = c("year", "ADM1_CODE", "ADM0_CODE")]
  admin_0 <- cell_pred[, lapply(c(overs, "pop_raked"), function(x) sum(get(x), na.rm = T)), by = c("year", "ADM0_CODE")]
  
  setnames(admin_2, grep("V[0-9]", names(admin_2), value = T), c(overs, "pop_raked"))
  setnames(admin_1, grep("V[0-9]", names(admin_1), value = T), c(overs, "pop_raked"))
  setnames(admin_0, grep("V[0-9]", names(admin_0), value = T), c(overs, "pop_raked"))
  
  admin_0$pop <- admin_0$pop_raked
  admin_0$pop_raked <- NULL
  admin_1$pop <- admin_1$pop_raked
  admin_1$pop_raked <- NULL
  admin_1$ADM0_CODE <- NULL
  admin_2$pop <- admin_2$pop_raked
  admin_2$pop_raked <- NULL
  admin_2$ADM0_CODE <- NULL
  gc()
  
  ## save unraked counts aggregations
  save(admin_0, admin_1, admin_2, sp_hierarchy_list,
       file = ifelse(is.null(custom_output_folder),
                     paste0(main_dir, "/", indicator, "_", "unraked", "_c", "_admin_draws_eb_bin", age, "_", reg, "_", holdout, ".RData"),
                     paste0(custom_output_folder, "/", indicator, "_", "unraked", "_c", "_admin_draws_eb_bin", age, "_", reg, "_", holdout, ".RData")
       )
  )
  
  
  ############################################################
  # creating unraked rates aggregations
  message("creating a unraked rates aggregations")
  
  # convert back to rates
  overs <- paste0("V", 1:ndraws)
  admin_0 <- admin_0[, (overs) := lapply(overs, function(x) get(x) / pop) ]
  admin_1 <- admin_1[, (overs) := lapply(overs, function(x) get(x) / pop) ]
  admin_2 <- admin_2[, (overs) := lapply(overs, function(x) get(x) / pop) ]
  
  ## save unraked rates aggregations
  save(admin_0, admin_1, admin_2, sp_hierarchy_list,
       file = ifelse(is.null(custom_output_folder),
                     file.path(main_dir,
                               sprintf('%s_unraked_admin_draws_eb_bin%i_%s_%i.RData',
                                       indicator, age, reg, holdout)),
                     file.path(custom_output_folder,
                               sprintf('%s_unraked_admin_draws_eb_bin%i_%s_%i.RData',
                                       indicator, age, reg, holdout))
       )
  )
  
  ## Return cell pred and raking factors if desired
  if (return_objects) {
    output_list <- list()
    output_list[["rf"]] <- data.table(fractional_rf)
    output_list[["raked_cell_pred"]] <- raked_cell_pred
    return(output_list)
  }
}
