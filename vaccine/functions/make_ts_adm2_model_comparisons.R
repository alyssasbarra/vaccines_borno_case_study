make_ts_adm2_model_comparisons <- function(indicator_group = indicator_group,
                                           indicator       = indicator,
                                           run_date        = run_date,
                                           shapefile_version = 'current',
                                           regions=Regions,
                                           tol=0.01,
                                           version= 'current',
                                           z_list,
                                           plot_agebin_faceted = TRUE,
                                           plot_year_faceted = TRUE,
                                           models,
                                           model_dates) {
  
  
  #Functions############
  # Title generate makes the title Grob that goes in the right hand corner of the pdf
  title_generate <- function(ind_title,
                             year_list,
                             admin = 0,
                             ad0_reg_title = "",
                             ctry = "",
                             ad1_df_ad1 = ""){
    arrangeGrob(textGrob("", gp = gpar(fontsize = 30)),
                textGrob(str_wrap(ind_title, 18),
                         gp = gpar(fontsize = title_grob_size, fontface = "bold")),
                textGrob("", gp = gpar(fontsize = 10)),
                textGrob(if (admin == 0) ifelse(ad0_reg_title == "All countries", NULL, ad0_reg_title) else ctry,
                         gp = gpar(fontsize = 20)),
                textGrob("", gp = gpar(fontsize = 10)),
                textGrob(if (admin == 0) "National estimates" else if (admin == 1) "By First-level Administrative Unit" else paste0("Admin 1: ", unique(ad1_df_ad1$ADM1_NAME)),
                         gp = gpar(fontsize = 15)),
                textGrob("", gp = gpar(fontsize = 10)),
                textGrob(paste0(min(year_list), " - ", max(year_list)),
                         gp = gpar(fontsize = 18)),
                ncol = 1,
                heights = c(30, 25, 10, 25, 10, 15, 10, 20))
  }
  
  # Plot overlay defines where each object (ggplot and title) are placed on the pdf.
  # This depends on if data is being plotted, and if multiple runs are used
  plot_overlay <- function(plot_data, multiple_runs){
    if (plot_data == T) {
      if (multiple_runs == T) {
        lay <- rbind(
          c(1, 1, 1, 1, 1, 2, 2, 3, 3),
          c(1, 1, 1, 1, 1, 2, 2, 3, 3),
          c(1, 1, 1, 1, 1, 4, 4, 4, NA),
          c(1, 1, 1, 1, 1, 4, 4, 4, NA),
          c(1, 1, 1, 1, 1, 4, 4, 4, 5)
        )
      } else {
        lay <- rbind(
          c(1, 1, 1, 1, 1, 2, 2, 3, 3),
          c(1, 1, 1, 1, 1, 2, 2, 3, 3),
          c(1, 1, 1, 1, 1, 4, 4, 4, 4),
          c(1, 1, 1, 1, 1, 4, 4, 4, 4),
          c(1, 1, 1, 1, 1, 4, 4, 4, 4)
        )
      }
    } else {
      lay <- rbind(
        c(1, 1, 1, 1, 2, 2, 3, 3),
        c(1, 1, 1, 1, 2, 2, 3, 3),
        c(1, 1, 1, 1, 4, 4, 4, 4),
        c(1, 1, 1, 1, 4, 4, 4, 4),
        c(1, 1, 1, 1, 4, 4, 4, 4),
        c(1, 1, 1, 1, NA, NA, NA, 5)
      )
    }
    return(lay)
  }
  
  ####################
  
  ## Set repo
  #core_repo  <- paste0("FILEPATH")
  #if (!dir.exists(core_repo)) core_repo <- "FILEPATH"
  #indic_repo <- paste0("FILEPATH")
  
  ## Load libraries and  MBG project functions.
  #source(paste0('FILEPATH))
  #package_list <- readLines(paste0("FILEPATH"))
  #mbg_setup(package_list = package_list, repos = core_repo)
  # Set run_date, indicator, indicator_group, out_dir per your preferences
  
  
  #config <- set_up_config(repo            = indic_repo,
  #                         core_repo       = core_repo,
  #                         indicator       = "",
  #                         indicator_group = "",
  #                         config_file     = paste0(indic_repo, 'mbg/', indicator, '/2_modeling/hiv_config.csv'),
  #                         covs_file       = paste0(indic_repo, 'mbg/', indicator, '/2_modeling/cov_list_std_and_hiv.csv'),
  #                         run_tests       = TRUE)
  
  #if (class(Regions) == "character" & length(Regions) == 1) Regions <- eval(parse(text=Regions))
  if (class(year_list) == "character") year_list <- eval(parse(text=year_list))
  # if (length(summstats) == 1 & grepl("c(", summstats, fixed =T)) summstats <- eval(parse(text=summstats))
  if (class(z_list) == "character") z_list <- eval(parse(text=z_list))
  
  
  ad0_df_mms<-NULL
  ad1_df_mms<-NULL
  ad2_df_mms<-NULL  
  
  
  for(i in rev(1:length(models))) {
    share_dir <- paste0("FILEPATH")
    in_dir <- paste0("FILEPATH")
    out_dir<-paste0('FILEPATH')
    dir.create(out_dir)
    out_filename_format = "subnational_model_comparison_w_uncertainty_%s"

    
    ad0_df<-NULL
    ad1_df<-NULL
    ad2_df<-NULL
    for(rake_lev in c('unraked','raked')){
      for (age in z_list){
        in_file_ad0 <- paste0(in_dir, indicator, "_admin_0_", rake_lev, "_bin", age,  "_summary.csv")
        in_file_ad1 <- paste0(in_dir, indicator, "_admin_1_", rake_lev, "_bin", age,  "_summary.csv")
        in_file_ad2 <- paste0(in_dir, indicator, "_admin_2_", rake_lev, "_bin", age,  "_summary.csv")
        
        # Prepare inputs #######################################################
        
        ad0_df_bin <- fread(in_file_ad0)
        ad1_df_bin <- fread(in_file_ad1)
        ad2_df_bin <- fread(in_file_ad2)
        
        ad0_df_bin$agebin<-age
        ad1_df_bin$agebin<-age
        ad2_df_bin$agebin<-age
        
        
        ad0_df_bin$run<-rake_lev
        ad1_df_bin$run<-rake_lev
        ad2_df_bin$run<-rake_lev
        
        
        ad0_df<-rbind(ad0_df, ad0_df_bin)
        ad1_df<-rbind(ad1_df, ad1_df_bin)
        ad2_df<-rbind(ad2_df, ad2_df_bin)
        
        rm(ad0_df_bin)
        rm(ad1_df_bin)
        rm(ad2_df_bin)
      }
    }
    # Drop Ma'tan al-Sarra if present
    ad0_df <- subset(ad0_df, ADM0_CODE != 40762)
    ad1_df <- subset(ad1_df, ADM0_CODE != 40762)
    ad2_df <- subset(ad2_df, ADM0_CODE != 40762)
    
    ad0_df$model<-models[i]
    ad1_df$model<-models[i]
    ad2_df$model<-models[i]
    
    ad0_df_mms<-rbind(ad0_df_mms, ad0_df)
    ad1_df_mms<-rbind(ad1_df_mms, ad1_df)
    ad2_df_mms<-rbind(ad2_df_mms, ad2_df)
  }
  
  #Make sure writing to designated run date
  share_dir <- paste0("FILEPATH")
  in_dir <- paste0("FILEPATH")
  out_dir<-paste0('FILEPATH')
  
  # Load GAUL shapefiles #################################################
  
  ad0_shape <- readOGR(get_admin_shapefile(admin_level=0))
  ad1_shape <- readOGR(get_admin_shapefile(admin_level=1))
  ad2_shape <- readOGR(get_admin_shapefile(admin_level=2))
  # Input data
  # set regions names
  region_names <- c(cssa = 'Central SSA', essa_sdn = 'Eastern SSA', sssa = 'Southern SSA', wssa = 'Western SSA')[regions]
  region_names[is.na(region_names)] <- regions[is.na(region_names)]
  message('made it to 156')
  # load collapsed input data
  in_dir <- paste("FILEPATH", sep='/')
  input_data <- fread(paste0('FILEPATH'))
  message('made it to 162') 

  input_data[, nid := as.double(svy_id)]
  input_data[nid %in% input_data[, uniqueN(year), 'nid'][V1 > 2, nid], nid := 10000*nid + year]
  input_data[,V1:=NULL]
  input_data[,svy_id:=NULL]
  
  ad0_data<-NULL
  ad1_data<-NULL
  ad2_data<-NULL
  for (age in z_list) {
    id<-input_data[age_bin==age & !is.na(latitude)]
    # get aggregated estimates
    admin_data <- input_aggregate_admin(indicator, indicator_group, input_data = id, 
                                        regions = regions, shapefile_version = shapefile_version,
                                        sample_column = 'weighted_n')
    
    ad0_data_bin<-admin_data$ad0
    ad1_data_bin<-admin_data$ad1
    ad2_data_bin<-admin_data$ad2
    
    ad0_data_bin$agebin<-age
    ad1_data_bin$agebin<-age
    ad2_data_bin$agebin<-age
    
    
    ad0_data<-rbind(ad0_data, ad0_data_bin)
    ad1_data<-rbind(ad1_data, ad1_data_bin)
    ad2_data<-rbind(ad2_data, ad2_data_bin)
    
    rm(ad0_data_bin)
    rm(ad1_data_bin)
    rm(ad2_data_bin)
  }
  #Plot
  message('made it to 203')
  library(ggrepel)
  library(gridExtra)
  library(ggplot2)
  library(stringr)
  
  ad0_map_regions = regions
  
  subset_codes <- get_adm0_codes(ad0_map_regions, shapefile_version = shapefile_version)
  
  # Wrapper for gSimplify() to allow it to work with SPDFs
  simplify_spdf <- function(the_spdf, ...) {
    simple_spdf <- gSimplify(the_spdf, topologyPreserve = T, ...)
    return(SpatialPolygonsDataFrame(simple_spdf, the_spdf@data))
  }
  
  # Use simplified GAUL shapefiles as default if none provided, if provided simplify the shapefile for speed
  
  ad0_shape_simple <- readOGR(get_admin_shapefile(admin_level = 0, raking = F, version = shapefile_version)) %>% subset(ADM0_CODE %in% subset_codes)
  ad1_shape_simple <- readOGR(get_admin_shapefile(admin_level = 1, raking = F, version = shapefile_version)) %>% subset(ADM0_CODE %in% subset_codes)
  ad2_shape_simple <- readOGR(get_admin_shapefile(admin_level = 2, raking = F, version = shapefile_version)) %>% subset(ADM0_CODE %in% subset_codes)
  
  ad0_shape_simple <- simplify_spdf(ad0_shape_simple, tol)
  ad1_shape_simple <- simplify_spdf(ad1_shape_simple, tol)
  ad2_shape_simple <- simplify_spdf(ad2_shape_simple, tol)
  
  
  # Merge on ihme_lc_ids
  gaul_to_loc_id <- get_location_code_mapping(shapefile_version = shapefile_version)
  gaul_to_loc_id <- subset(gaul_to_loc_id, select = c("GAUL_CODE", "ihme_lc_id"))
  setnames(gaul_to_loc_id,"GAUL_CODE", "ADM0_CODE")
  
  #renaming the ADMX_NAME columns using the names in ad2_shape, because this functions matches
  #everything on names and things get dropped if the names don't line up perfectly
  #This is an alternative to rewriting the entire function to match on coes
  admin_shp_data_adm0 <- as.data.table(ad2_shape_simple)
  admin_shp_data_adm0 <- unique(admin_shp_data_adm0[,c("ADM0_CODE", "ADM0_NAME")])
  admin_shp_data_adm0$ADM0_CODE <- as.integer(as.character(admin_shp_data_adm0$ADM0_CODE))
  admin_shp_data_adm0$ADM0_NAME <- as.character(admin_shp_data_adm0$ADM0_NAME)
  gaul_to_loc_id_adm0 <- merge(gaul_to_loc_id, admin_shp_data_adm0, by = "ADM0_CODE")
  
  admin_shp_data_adm1 <- as.data.table(ad2_shape_simple)
  admin_shp_data_adm1 <- unique(admin_shp_data_adm1[,c("ADM0_CODE", "ADM0_NAME", "ADM1_CODE", "ADM1_NAME")])
  admin_shp_data_adm1$ADM1_CODE <- as.integer(as.character(admin_shp_data_adm1$ADM1_CODE))
  admin_shp_data_adm1$ADM0_CODE <- as.integer(as.character(admin_shp_data_adm1$ADM0_CODE))
  admin_shp_data_adm1$ADM1_NAME <- as.character(admin_shp_data_adm1$ADM1_NAME)
  admin_shp_data_adm1$ADM0_NAME <- as.character(admin_shp_data_adm1$ADM0_NAME)
  gaul_to_loc_id_adm1 <- merge(gaul_to_loc_id, admin_shp_data_adm1, by = "ADM0_CODE")
  
  admin_shp_data_adm2 <- as.data.table(ad2_shape_simple)
  admin_shp_data_adm2 <- unique(admin_shp_data_adm2[,c("ADM0_CODE", "ADM0_NAME", "ADM1_CODE", "ADM1_NAME", "ADM2_CODE", "ADM2_NAME")])
  admin_shp_data_adm2$ADM2_CODE <- as.integer(as.character(admin_shp_data_adm2$ADM2_CODE))
  admin_shp_data_adm2$ADM1_CODE <- as.integer(as.character(admin_shp_data_adm2$ADM1_CODE))
  admin_shp_data_adm2$ADM0_CODE <- as.integer(as.character(admin_shp_data_adm2$ADM0_CODE))
  admin_shp_data_adm2$ADM2_NAME <- as.character(admin_shp_data_adm2$ADM2_NAME)
  admin_shp_data_adm2$ADM1_NAME <- as.character(admin_shp_data_adm2$ADM1_NAME)
  admin_shp_data_adm2$ADM0_NAME <- as.character(admin_shp_data_adm2$ADM0_NAME)
  gaul_to_loc_id_adm2 <- merge(gaul_to_loc_id, admin_shp_data_adm2, by = "ADM0_CODE")
  
  ad0_df_mms[, c("ADM0_NAME")] <- NULL
  ad1_df_mms[, c("ADM0_NAME", "ADM1_NAME", "ADM0_CODE")] <- NULL
  ad2_df_mms[, c("ADM0_NAME", "ADM1_NAME", "ADM2_NAME", "ADM0_CODE", "ADM1_CODE")] <- NULL
  
  ad0_df <- merge(copy(as.data.table(ad0_df_mms)), gaul_to_loc_id_adm0, all.x=T, all.y=F, by = 'ADM0_CODE')
  ad1_df <- merge(copy(as.data.table(ad1_df_mms)), gaul_to_loc_id_adm1, all.x=T, all.y=F, by = 'ADM1_CODE')
  ad2_df <- merge(copy(as.data.table(ad2_df_mms)), gaul_to_loc_id_adm2, all.x=T, all.y=F, by = 'ADM2_CODE')
  
  
  ad0_data <- merge(copy(as.data.table(ad0_data)), gaul_to_loc_id, all.x=T, all.y=F, by = 'ADM0_CODE')
  ad1_data <- merge(copy(as.data.table(ad1_data)), gaul_to_loc_id, all.x=T, all.y=F, by = 'ADM0_CODE')
  ad2_data <- merge(copy(as.data.table(ad2_data)), gaul_to_loc_id, all.x=T, all.y=F, by = 'ADM0_CODE')
  
  plot_levels<-c('ad0', 'ad1', 'ad2')
  message("Simplifying shapes")
  
  # Subset and simplify shapefiles for memory & speed
  ad0_codes <- unique(ad0_df$ADM0_CODE)
  # If any region of Africa is being used, use full map of Africa so it does not look disjointed
  if (any(ad0_map_regions %in% c("cssa", "wssa", "essa", "sssa", "name"))) ad0_codes <- union(ad0_codes, get_adm0_codes("Africa", shapefile_version = shapefile_version))
  subset_codes <- ad0_codes
  
  # Make sure that the admin codes are converted to numeric if they are factors
  if (is.factor(ad0_shape_simple@data$ADM0_CODE)) ad0_shape_simple@data$ADM0_CODE <- as.numeric(levels(ad0_shape_simple@data$ADM0_CODE))[ad0_shape_simple@data$ADM0_CODE]
  if (is.factor(ad1_shape_simple@data$ADM1_CODE)) ad1_shape_simple@data$ADM1_CODE <- as.numeric(levels(ad1_shape_simple@data$ADM1_CODE))[ad1_shape_simple@data$ADM1_CODE]
  if (is.factor(ad2_shape_simple@data$ADM2_CODE)) ad2_shape_simple@data$ADM2_CODE <- as.numeric(levels(ad2_shape_simple@data$ADM2_CODE))[ad2_shape_simple@data$ADM2_CODE]
  
  # Add simple world shapefile if region is specified as Africa
  if ("ad0" %in% plot_levels & "africa" %in% tolower(ad0_map_regions)) world_shape_simple <- readRDS('FILEPATH')
  
  
  
  ad0_code_list <- lapply(ad0_map_regions, get_adm0_codes, shapefile_version = shapefile_version)
  ad0_map_region_titles = region_names
  names(ad0_code_list) <- ad0_map_region_titles
  
  ad0_code_list<-lapply(1:length(regions), function(i){
    ad0_code_list[[i]] <- as.data.table(ad0_code_list[[i]])
  })
  
  ad0_codes<-rbindlist(ad0_code_list)
  
  
  
  
  #ad0_reg_codes <- ad0_code_list[[i]]
  #ad0_reg_title <- names(ad0_code_list)[i]
  
  ad0_df_country <- ad0_df
  ad0_data_country <- ad0_data
  
  val_range = c(0,1)
  title_plot_size <- 20
  ind_title = "MCV1 Coverage"
  
  # First, time series plot of all admin0s in the region
  n_ad0s <- length(unique(ad0_df_country$ADM0_CODE))
  
  # Try to wrap names of countries as needed
  if (n_ad0s >  0 & n_ad0s <= 16) wrap_width <- 18
  if (n_ad0s > 16 & n_ad0s <= 25) wrap_width <- 12
  if (n_ad0s > 25)                wrap_width <- 9
  
  if (n_ad0s < 25) {
    ad0_df_country[, plot_name := stringr::str_wrap(ADM0_NAME, width = wrap_width)]
  } else {
    ad0_df_country[, plot_name := stringr::str_wrap(ADM0_NAME, width = wrap_width)]
  }
  
  # When plotting data, aligning plot names depending on number of admin0 in region
  if (n_ad0s < 25) {
    ad0_data_country[, plot_name := stringr::str_wrap(ADM0_NAME, width = wrap_width)]
  } else {
    ad0_data_country[, plot_name := stringr::str_wrap(ADM0_NAME, width = wrap_width)]
  }
  #####
  
  
  
  lay<-rbind(
    c(1, 1, 1, 1, 1, 2, 2),
    c(1, 1, 1, 1, 1, 2, 2),
    c(1, 1, 1, 1, 1, 2, 2),
    c(1, 1, 1, 1, 1, 2, 2),
    c(1, 1, 1, 1, 1, 2, 2))
  
  
  
  
  
  #------------------------------------------------------------------------------------------------
  if (plot_agebin_faceted==TRUE){
    pdf_filename <- paste0(out_dir, sprintf(out_filename_format, regions), "_by_admin2_x_as_year.pdf")
    
    
    message("\n############################################################")
    message("############################################################")
    message(paste0("Plotting Admin-2-level estimates by Admin 1\n"))
    message(paste0("  Writing output to ", pdf_filename, "\n"))
    
    
    
    pdf(file = pdf_filename,
        height = 10,
        width = 18)
    
    n_plots<-0
    
    
    #ad0_shape_simple$ord<- 1:length(ad0_shape_simple)
    #ad0_shape_simple[ad0_shape_simple$ADM0_NAME=='Mauritania',]$ord #26
    
    
    #----------------------------------------------------------------------------------------------
    for (i in ad0_shape_simple$ADM0_CODE) {
      print(ad1_shape_simple[ad1_shape_simple$ADM0_CODE==i,]$ADM0_NAME[1])
      for (j in ad1_shape_simple[ad1_shape_simple$ADM0_CODE==i,]$ADM1_CODE) {
        print(j)
        for (k in ad2_shape_simple[ad2_shape_simple$ADM1_CODE==j,]$ADM2_CODE) {
        model_data=ad2_df[ADM0_CODE==i & ADM1_CODE==j & ADM2_CODE==k]
        input_data=ad2_data[ADM0_CODE==i & ADM1_CODE==j & ADM2_CODE==k]
        
        if (nrow(input_data) == 0 & nrow(model_data) == 0){
          next
        }
        
        if (any(is.na(input_data$mean)) | any(is.na(model_data$mean))){
          message('skipping ^this one due to NAs in estimates')
          next
        }
        
        admin = 0
        
        
        title_plot_size=title_plot_size
        ind_title='MCV1 coverage'
        
        carto_discrete <- rep(c("#7F3C8D","#11A579","#F2B701","#E73F74",
                                "#3969AC","#80BA5A","#E68310","#008695",
                                "#CF1C90","#f97b72","#4b4b8f","#A5AA99"), 3)
        model_data$agebin<-as.factor(model_data$agebin)
        input_data$agebin<-as.factor(input_data$agebin)
        levels(model_data$agebin)<-c(#'6-8m',
          '9-11m',
          '1y',
          '2y',
          '3y',
          '4y')
        levels(input_data$agebin)<-c(#'6-8m',
          '9-11m',
          '1y',
          '2y',
          '3y',
          '4y')
        
        
        ####################################
        gg_ad0_ts <-
          ggplot() +
          geom_ribbon(data = model_data[run=='raked',], aes(x = year, ymin = lower, ymax = upper, fill=as.factor(model)), alpha = 0.3) +
          geom_line(data = model_data, aes(x = year, y = mean, color=as.factor(model), lty=run)) +
          theme_bw(base_size = 16)
        if (nrow(input_data) > 0){
          gg_ad0_ts <- gg_ad0_ts + 
            geom_point(data = input_data,
                       aes(x = year, y = outcome, size = N, shape = point), alpha = 0.7) 
        }
        gg_ad0_ts <- gg_ad0_ts +
          facet_wrap( ~ agebin, ncol=2) +
          coord_cartesian(ylim = c(0,1)) +
          scale_shape_manual("Data type", values = c(17,16), label = c("polygon", "point"), drop = F) +
          scale_color_manual("", values = carto_discrete) +
          scale_fill_manual("", values = carto_discrete) +
          scale_size(range = c(1,5)) +
          theme(
            strip.background = element_blank(),
            plot.caption = element_text(hjust = 0.5),
            plot.title = element_text(size = 16, face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1),
            # legend.title = element_text(size = 10),
            # legend.text = element_text(size = 7),
            legend.justification = "top"
          ) +
          guides( fill = guide_legend(order = 1, override.aes = list(size = 5, shape = 21, lty='solid'))) +
          labs(y = ind_title, size='Sample size',
               title = paste0(model_data[1]$ADM1_NAME, ": ", model_data[1]$ADM2_NAME),
               x='Year')
        
        
        ###############################################################
        ad0_shape_country <- subset(ad0_shape_simple, ADM0_CODE==i)
        ad1_shape_country <- subset(ad1_shape_simple, ADM0_CODE==i)
        ad2_shape_country <- subset(ad2_shape_simple, ADM0_CODE==i)
        #   ad1_shape_country <- st_as_sf(subset(ad1_shape_simple, ADM0_CODE==i))
        ad1_shape_adm1    <- subset(ad1_shape_country, ADM1_CODE==j)
        ad2_shape_adm2    <- subset(ad2_shape_country, ADM2_CODE==k)
        
        gg_0 <- 
          ggplot() +
          geom_polygon_quiet(data = ad1_shape_adm1,
                             aes(x = long, y = lat, group = group),
                             fill = "light gray") +
          geom_polygon_quiet(data = ad2_shape_adm2,
                             aes(x = long, y = lat, group = group),
                             fill = "red") +
          geom_path_quiet(data = ad1_shape_country,
                          aes(x=long, y=lat, group=group),
                          size = 0.2) +
          geom_path_quiet(data = ad0_shape_country,
                          aes(x=long, y = lat, group=group),
                          size = 0.5) +
          coord_equal()+
          theme_void()
        
        
        
        
        
        main_plot <- arrangeGrob(gg_ad0_ts, gg_0, layout_matrix = lay)
        grid.draw(main_plot)
        n_plots<- n_plots + 1
        if (n_plots <713 ) plot.new()
      }
      }
    }
    dev.off()
    
  }
  #--------------------------------------------------------------------------
  ##Same thing but switch age & year
  if (plot_year_faceted==TRUE) {
    pdf_filename <- paste0(out_dir, sprintf(out_filename_format, regions), "_by_admin2_x_as_age.pdf")
    
    
    message("\n############################################################")
    message("############################################################")
    message(paste0("Plotting Admin-2-level estimates by Admin 1\n"))
    message(paste0("  Writing output to ", pdf_filename, "\n"))
    
    
    pdf(file = pdf_filename,
        height = 10,
        width = 18)
    
    n_plots<-0
    
    mround <- function(x,base){ 
      base*round(x/base) 
    } 
    
    
    
    for (i in ad0_shape_simple$ADM0_CODE) {
      print(ad1_shape_simple[ad1_shape_simple$ADM0_CODE==i,]$ADM0_NAME[1])
      for (j in ad1_shape_simple[ad1_shape_simple$ADM0_CODE==i,]$ADM1_CODE) {
        print(j)
        for (k in ad2_shape_simple[ad2_shape_simple$ADM1_CODE==j,]$ADM2_CODE) {
          model_data=ad2_df[ADM0_CODE==i & ADM1_CODE==j & ADM2_CODE==k]
          input_data=ad2_data[ADM0_CODE==i & ADM1_CODE==j & ADM2_CODE==k]
        
        if (nrow(input_data) == 0 & nrow(model_data) == 0){
          next
        }
        
        admin = 0
        
        title_plot_size=title_plot_size
        ind_title='MCV1 coverage'
        
        carto_discrete <- rep(c("#7F3C8D","#11A579","#F2B701","#E73F74",
                                "#3969AC","#80BA5A","#E68310","#008695",
                                "#CF1C90","#f97b72","#4b4b8f","#A5AA99"), 3)
        
        
        gg_ad0_ts<-
          ggplot() +
          geom_ribbon(data = model_data[run=='raked',], aes(x = as.numeric(agebin), ymin = lower, ymax = upper, fill=as.factor(model)), alpha = 0.3) +
          geom_line(data = model_data, aes(x = as.numeric(agebin), y = mean, color=as.factor(model), lty=as.factor(run))) +
          theme_bw(base_size = 16)+
          #   if (nrow(input_data[year ==2000 | year ==2005 | year ==2010 | year ==2018]) > 0){
          #     gg_ad0_ts <- gg_ad0_ts + 
          geom_point(data = input_data,
                     aes(x = as.numeric(agebin), y = outcome, size = N, shape = as.factor(point))) +
          
          scale_x_continuous(labels = c(#'6-8m', 
            '9-11m', '1y', '2y', '3y', '4y'), breaks=1:5) +
          #  }
          #  gg_ad0_ts <- gg_ad0_ts +
          facet_wrap( ~ year) +
          coord_cartesian(ylim = c(0,1)) +
          scale_shape_manual("Data type", values = c(17,16), label = c("polygon", "point"), drop = F) +
          scale_size(range = c(1,5)) +
          scale_color_manual("", values = carto_discrete) +
          scale_fill_manual("", values = carto_discrete) +
          theme(
            strip.background = element_blank(),
            plot.caption = element_text(hjust = 0.5),
            plot.title = element_text(size = 16, face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.justification = "top"
          ) +
          guides( fill = guide_legend(order = 1, override.aes = list(size = 5, shape = 21, lty='solid')),
                  alpha = FALSE) +
          labs(y = ind_title,
               title = paste0(model_data[1]$ADM0_NAME, ": ", model_data[1]$ADM1_NAME))
        
        
        ###############################################################
        ad0_shape_country <- subset(ad0_shape_simple, ADM0_CODE==i)
        ad1_shape_country <- subset(ad1_shape_simple, ADM0_CODE==i)
        ad2_shape_country <- subset(ad2_shape_simple, ADM0_CODE==i)
        #   ad1_shape_country <- st_as_sf(subset(ad1_shape_simple, ADM0_CODE==i))
        ad1_shape_adm1    <- subset(ad1_shape_country, ADM1_CODE==j)
        ad2_shape_adm2    <- subset(ad2_shape_country, ADM2_CODE==k)
        
        gg_0 <- 
          ggplot() +
          geom_polygon_quiet(data = ad1_shape_adm1,
                             aes(x = long, y = lat, group = group),
                             fill = "light gray") +
          geom_polygon_quiet(data = ad2_shape_adm2,
                             aes(x = long, y = lat, group = group),
                             fill = "red") +
          geom_path_quiet(data = ad1_shape_country,
                          aes(x=long, y=lat, group=group),
                          size = 0.2) +
          geom_path_quiet(data = ad0_shape_country,
                          aes(x=long, y = lat, group=group),
                          size = 0.5) +
          coord_equal()+
          theme_void()
        
        main_plot <- arrangeGrob(gg_ad0_ts, gg_0, layout_matrix = lay)
        grid.draw(main_plot)
        n_plots<- n_plots + 1
        if (n_plots <713 ) plot.new()
      }
    }
  }
    dev.off()
  }

}