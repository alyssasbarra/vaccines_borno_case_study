fr_aggregate_pop_custom_shapefile <- function(shapefile_path,
                                              field,
                                              covs,
                                              measures,
                                              year_list,
                                              interval_mo,
                                              cores,
                                              link_table = NULL,
                                              id_raster = NULL) {
  #   #build link table if not passed in to function call
  if(is.null(link_table) & is.null(id_raster)) {
    message("Building link table for provided shapefile. Region or Global shapefiles may take several hours to run.")
    link_output <- build_link_table(shapefile_path=NULL, field, cores)
    link_table <- link_output$link_table
    id_raster <- link_output$id_raster
  }
  #   #loop through covariates, because load_and_crop can't do multiple measures of the same covariate
  agg_list = list()
  for(i in 1:length(covs)) {
    pop_raster_annual <- load_and_crop_covariates_annual(covs           = covs[i],                
                                                         measures       = measures[i],      
                                                         simple_polygon = id_raster,
                                                         start_year     = min(year_list),
                                                         end_year       = max(year_list),
                                                         interval_mo    = as.numeric(interval_mo),
                                                         agebin=1)[[1]]
    #get pixel_id, year, and covariate values into vectors of the same length
    id_raster_extract <- raster::extract(id_raster, extent(id_raster))
    year_extract <- rep(year_list, each = length(id_raster_extract))
    id_raster_extract <- rep(id_raster_extract, times = length(year_list))
    pop_raster_extract <- as.vector(raster::extract(pop_raster_annual, extent(pop_raster_annual)))
    if(length(id_raster_extract) != length(pop_raster_extract)) {
      message("cov raster and id raster have different length, this should never happen")
    }
    #make a table of pixel_id, year, and covariate values
    pop_table <- data.table(cbind(id_raster_extract, year_extract, pop_raster_extract))
    #merge onto link table by pixel_id (this duplicates rows for each year)
    link_pop <- merge(link_table, pop_table, by.x = "pixel_id", by.y = "id_raster_extract", all.x = T)
    #check the duplication didn't go wrong
    if(nrow(link_table) != (nrow(link_pop) / length(year_list))) {
      message("problem merging covariate and link table, extra rows added")
    }
    #calculate fractional covariate value
    link_pop[, frac_pop := pop_raster_extract * area_fraction]
    #aggregate to field and year
    agg_table <- link_pop[, .(pop = sum(frac_pop, na.rm = T)), by = c(field, "year_extract")]
    setnames(agg_table, "pop", paste0(covs[i], "_", measures[i]))
    agg_list[[paste0(covs[i], "_", measures[i])]] <- agg_table
  }
  #combine covariates into a single table
  if(length(agg_list) > 1){
    agg_table <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = c(field, "year_extract")), agg_list)
  } else {
    agg_table <- agg_list[[1]]
  }
  setnames(agg_table, c("year_extract"), "year")
  return(list("agg_table" = agg_table, "link_table" = link_table, "id_raster" = id_raster))
}