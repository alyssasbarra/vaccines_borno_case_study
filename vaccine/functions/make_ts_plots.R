make_ts_plots <- function(indicator, indicator_group, run_date, regions, output_dir, shapefile_version = 'current', counts = F, plot_levels, tol=0.1) {
  #indicator
  #indicator_group
  #run_date
  #regions <- strata
  #output_dir <- paste0(outputdir, '/diagnostic_plots/')
  
  version = shapefile_version
  
  #Set up for counts file pulls and naming
  if(counts == T) cts <- 'c_' else cts <- ''
  
  # set regions names
  region_names <- regions
  
  # load collapsed input data
  in_dir <- paste("FILEPATH")
  input_data <- fread(paste0('FILEPATH'))
  input_data[,source := 'Data']
  input_data[,nid:=svy_id]
  input_data[,svy_id:=NULL]
  
  input_data[, nid := as.double(nid)]
  input_data[nid %in% input_data[, uniqueN(year), 'nid'][V1 > 2, nid], nid := 10000*nid + year]
  input_data[,N_copy:=N]
  
  # get aggregated estimates
  admin_data <- input_aggregate_admin(indicator, indicator_group, input_data = input_data,
                                      regions = regions, shapefile_version = shapefile_version, sample_column="N_copy")
  
  # reset the polygon variable based on original input data (the function tries to guess based on
  # resampling weights, but gets this wrong if these are used to adjust N rather than as analytic
  # weights)
  input_data_poly <- input_data[, list(polygon = as.numeric(mean(point) < 1)), 'nid,country']
  input_data_poly[, polygon := as.factor(polygon)]
  input_data_poly <- merge(input_data_poly, get_location_code_mapping(shapefile_version = shapefile_version), by.x = 'country', by.y = 'ihme_lc_id')
  input_data_poly <- input_data_poly[, list(svy_id = nid, ADM0_CODE = ADM_CODE, polygon)]
  
  admin_data <- lapply(admin_data, function(d) {
    d[, polygon := NULL]
    d <- merge(d, input_data_poly)
    return(d)
  })
  
  # get aggregated preds
  admin_preds <- lapply(0:2, function(d) {
    raked<-fread(paste0(in_dir, indicator, "_admin_", d, "_", cts, "raked_summary.csv"))
    raked[,run:='raked']
    unraked<-fread(paste0(in_dir, indicator, "_admin_", d, "_", cts, "unraked_summary.csv"))
    unraked[,run:='unraked']
    df<-rbind(raked,unraked)
    return(df)
  })
  names(admin_preds) <- paste0('ad', 0:2)
  # run subnational_ts_plots() function
  subnational_ts_plots(ad0_df = admin_preds$ad0, ad1_df = admin_preds$ad1, ad2_df = admin_preds$ad2,
                       ad0_data = admin_data$ad0, ad1_data = admin_data$ad1, ad2_data = admin_data$ad2,
                       ind_title = 'Coverage', highisbad = F, val_range = c(0, 1),
                       plot_levels = plot_levels, multiple_runs = T,
                       ad0_map_regions = regions, ad0_map_region_titles = region_names,
                       out_dir = output_dir, plot_data = T, verbose = T, shapefile_version = shapefile_version, tol=tol, 
                       out_filename_format = paste0("subnational_ts_plots_", regions, "_%s.pdf"))
  
  return(paste('Plots saved to:', output_dir))
}
