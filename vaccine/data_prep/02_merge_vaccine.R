# HEADER ------------------------------------------------------------------
# Author: USERNAME
# Date: DATE
# Project: Vaccines: Data Prep
# Purpose: Merge vaccine data together, clean, apply case definitions
# Details: 
# source("FILEPATH")
#************************************************************************** 

# Note: run from 00_main.R, which sets up all directories & files

### I: LOAD DATA ##################################################################

# Load in the Hib introduction data

intro_dt <- readRDS(intro_file) %>% 
  as.data.table %>%
  subset(., select = c("ihme_loc_id", "cv_intro", "me_name")) %>%
  unique

# Load all vaccine data
message("Loading vaccine data... \n")
vaccine_data <- readRDS(paste0(in_dir, in_run_date, '.Rda')) %>% 
  as.data.table

#### Load in and apply outliers ----------------------------------------------------

# Read in shared outlier table
df_outlier <- fread("FILEPATH")

# 1. Batch outliers
nids_to_outlier <- df_outlier[batch_outlier == 1 & lbd == 1, nid]

# 2. Vaccine-specific outliers
if (vaccine_to_run == "mcv") {
  gbd_names <- "vacc_mcv1"
} else if (vaccine_to_run == "dpt") {
  gbd_names <- c("vacc_dpt3", "vacc_dpt2", "vacc_dpt1")  
} else {
  message(paste0("Need to add outlier code for ", vaccine_to_run))
  warning("No outliers applied")
}

nids_to_outlier <- c(nids_to_outlier, 
                     df_outlier[(me_name %in% gbd_names & lbd == 1), nid])
nids_to_outlier <- unique(nids_to_outlier)
message(paste0("Outliering ", length(nids_to_outlier), " nids and saving table to log directory..."))

# Save an outlier csv to the log directory
df_outliered <- subset(df_outlier, nid %in% nids_to_outlier & (batch_outlier == 1 | me_name %in% gbd_names))
fwrite(df_outliered, paste0(log_dir, vaccine_to_run, "_outliered_nids.csv"))

# Outlier
nrow_before <- nrow(vaccine_data)
vaccine_data <- subset(vaccine_data, !(nid %in% nids_to_outlier))
nrow_after <- nrow(vaccine_data)

message(paste0("Dropped ", nrow_before - nrow_after, " out of ", nrow_before, 
               " (", round(((nrow_before - nrow_after)/nrow_before)*100, 2), 
               "%) rows due to outliering"))

### II: DATA CLEANING ##############################################################
# General data cleaning (applies to all vaccines)

message("cleaning data... \n")
# combine age categories
if ("child_age_month" %in% names(vaccine_data)) {
  vaccine_data[!is.na(child_age_month) & is.na(age_month), age_month := child_age_month]
  vaccine_data[, child_age_month := NULL]
}

if ("child_age_year" %in% names(vaccine_data)) {
  vaccine_data[!is.na(child_age_year) & is.na(age_year), age_year := child_age_year]
  vaccine_data[, child_age_year := NULL]
}

if ("child_sex_id" %in% names(vaccine_data)) {
  vaccine_data[, child_sex_id := NULL]
  vaccine_data[!is.na(child_sex_id) & is.na(sex_id), sex_id := child_sex_id]
}

# drop some variables
drop_vars <- c("smaller_site_unit", "line_id", "hh_id", 
               "admin_2_id", "admin_2_mapped", "admin_2",
               "admin_3_id", "admin_3_mapped", "admin_3",
               "admin_4_id", "admin_4_mapped", "admin_4",
               "admin_1_id", "admin_1_mapped", "file_path", 
               "nid.x", "iso3")

vaccine_data <- vaccine_data[,!(names(vaccine_data) %in% drop_vars), with = F]

# variable for year midpoint
# now switched to year_start
vaccine_data[, year_mid := (year_start + year_end) / 2]

#### calculate pweight from hhweight if pweight is missing
vaccine_data[is.na(pweight) & is.na(latitude) & !(is.na(hhweight)), pweight := hhweight]

# drop if not successfully geopositioned
vaccine_data[shapefile == "", shapefile := NA]
nrow_all <- nrow(vaccine_data)
geopos_pre <- vaccine_data[,.N, by=.(nid)]
names(geopos_pre) <- c("nid","nrow_pre_geoposition")

vaccine_data <- vaccine_data[!(is.na(lat) & is.na(shapefile))]
nrow_post_drop <- nrow(vaccine_data)
geopos_post <- vaccine_data[,.N, by=.(nid)]
names(geopos_post) <- c("nid","nrow_post_geoposition")

geopos_drop_log <- merge(geopos_pre, geopos_post, by="nid", all.x=T)
geopos_drop_log[is.na(nrow_post_geoposition), nrow_post_geoposition:=0]
message(paste0("Dropping ", nrow_all-nrow_post_drop,
               " rows not geopositioned (",
               round(((nrow_all-nrow_post_drop)/nrow_all)*100,2),
               "%)"))

geopos_drop_log[, percent_lost := ((nrow_pre_geoposition - nrow_post_geoposition)/nrow_pre_geoposition)*100]

message(paste0("Saving log file for rows dropped due to lack of geopositioning at: ", log_dir, "geopositioning_drop.csv"))
write.csv(geopos_drop_log, paste0(log_dir, "geopositioning_drop.csv"))

if(!"psu" %in% names(vaccine_data)){
  vaccine_data$psu <- vaccine_data$psu.y
}

### III. VACCINE-SPECIFIC PROCESSING ###########################################

# Loop over this code for each vaccine
# Generates a series of flat files (data frames)

for (vax in vaccines) {

  ### LOAD VACCINE-SPECIFIC PARAMETERS -------------------------------------------
  vax_df <- copy(vaccine_data)

  vax_prefix <- vaccine_list[[vax]][["prefix"]]
  vax_title  <- vaccine_list[[vax]][["title"]]
  vax_doses  <- vaccine_list[[vax]][["all_doses"]]
  age_range  <- vaccine_list[[vax]][["age_range"]]
  age_min <- as.numeric(min(age_range))
  age_max <- as.numeric(max(age_range))

  message("\n################################################")
  message(paste0("Working on ", vax_title, "...\n"))

  ### ELIGIBILITY AND AGE COHORTS ------------------------------------------------

  # Process age eligibility 
  message("Processing age eligibility ...")

  # Age information present?
  vax_df[!is.na(age_month), age_info := 1]
  vax_df[!is.na(age_year), age_info := 1]
  vax_df[!is.na(age_categorical), age_info := 1]

  total_ages <- vax_df[,.N,by="nid"]
  names(total_ages) <- c("nid","total_rows")
  no_age <- vax_df[!(is.na(age_info)),.N,by="nid"]
  names(no_age) <- c("nid","no_age_rows")

  no_age <- merge(no_age, total_ages, by="nid", all.y=T)
  no_age[is.na(no_age_rows),no_age_rows:=0]
  no_age[,percent_missing_ages := (total_rows - no_age_rows)/total_rows*100]
  # Drop all without ages
  message(paste0("Dropping ", nrow(vax_df[is.na(age_info)]), 
                 " out of ", nrow(vax_df), " rows (",
                 round(nrow(vax_df[is.na(age_info)])/nrow(vax_df) * 100, 1),
                 "%) without age information."))

  vax_df <- vax_df[!(is.na(age_info))]

  message(paste0("Saving log file for rows dropped due to lack of age information at: ", log_dir, "age_data_drop.csv"))
  write.csv(no_age, paste0(log_dir, "age_data_drop.csv"))

  # Check if there are any rows with age_year but not age_month
  if (nrow(vax_df[!is.na(age_month) & is.na(age_year)]) > 0) {
    warning(paste0("Warning: some rows have age in years, but not in months",
                   "\nNot currently equipped to handle this - need to modify code."))
  }

  # Note: This section uses age_month preferentially, then age_year if age_month not available
  # Also treating age_year == 1 as 12-23 months inclusive, etc.

  # Placeholder: Check if there are any categorical ages (not currently)
  if ("age_categorical" %in% names(vax_df)) {
  if (nrow(vax_df[!is.na(age_categorical) & is.na(age_month) & is.na(age_year)]) > 0) {
    warning("Dropping data without categorical ages - need to implement")
    message("#########################################################")
    message("Warning! Need to implement categorical ages!")
    message("You now have data that has this type of value - ")
    message("program this now!")
    message("#########################################################")
   }
  }

  # Drop those outside of age eligibility ranges
  n_before <- nrow(vax_df)
  age_range_pre <- vax_df[,.N, by=.(nid)]
  names(age_range_pre)<-c("nid","pre_age_range_drop")
  if (age_min == 12 & age_max == 59) {
    vax_df <- vax_df[(!is.na(age_month) & age_month >= 12) | (is.na(age_month) & age_year >= 1)]
    vax_df <- vax_df[(!is.na(age_month) & age_month <= 59) | (is.na(age_month) & age_year <= 4)]

    # Also drop a small number of rows with age_month < 60 but age_year > 4 
    vax_df <- vax_df[!(age_year > 4)]

  } else {
    stop("Not currently set up to handle age ranges aside from 12-59 mo.\nNeed to implement in code.")
  }

  message(paste0("Dropping ", n_before - nrow(vax_df),
                 " out of ", n_before, " rows (", 
                 round((n_before - nrow(vax_df))/n_before * 100, 2), 
                 "%) under ", age_min, " or over ", age_max, " months."))
  age_range_post <- vax_df[,.N, by=.(nid)]
  names(age_range_post)<-c("nid","post_age_range_drop")
  age_range_post <- merge(age_range_post, age_range_pre, by="nid", all.y=T)
  age_range_post[is.na(post_age_range_drop),post_age_range_drop := 0]
  age_range_post[, percent_dropped_age_range := (pre_age_range_drop - post_age_range_drop)/pre_age_range_drop*100]
  
  message(paste0("Saving log file for rows dropped from being out of the age range: ", log_dir, "age_range_drop.csv"))
  write.csv(age_range_post, paste0(log_dir, "age_range_drop.csv"))

  ### Generate birth cohorts ------------------------------------------------

  message("\nGenerating birth cohorts ...")

  if (age_min == 12 & age_max == 59) {
    vax_df[!is.na(age_month), birth_cohort := floor(age_month / 12)]
    vax_df[is.na(age_month), birth_cohort := as.numeric(age_year)]

    # Note: add one to year because estimating 12-23 month as target age bin
    vax_df[, year_start_cohort := year_start - birth_cohort + 1] 

  } else {
    stop("Not currently set up to handle age ranges aside from 12-59 mo.\nNeed to implement in code.")
  }

  ### Exclude children in a cohort before 2000 -------------------------------------------------
  if (drop_pre_2000 == T){
    #drop all children with a year_start_cohort before 2000
    message(paste0("\nDropping ", vax_df[year_start_cohort<2000,.N], " rows out of ",vax_df[,.N]," that start before 2000"))
    vax_df <- vax_df[year_start_cohort >= 2000]
  }

  ### Assign weights based on age cohorts ---------------------------------------

  ### Generate a categorical dummy variable for doses of vaccine given

  dose_col <- paste0(vax_prefix, "_dose")
  
  # Convert to numeric prior to vacc_ever check
  vax_df[[dose_col]] <- as.numeric(vax_df[[dose_col]])
#  vax_df[[paste0(vax_prefix, "_ever")]] <- as.numeric(vax_df[[paste0(vax_prefix, "_ever")]])

  # Ensure numeric
  vax_df[[dose_col]] <- as.numeric(vax_df[[dose_col]])

  for (i in 0:(max(vax_doses)-1)) {
    vax_df[get(dose_col) == i & !is.na(get(dose_col)), paste0(vax_prefix, "_dose_", i) := 1]
    vax_df[get(dose_col) != i & !is.na(get(dose_col)), paste0(vax_prefix, "_dose_", i) := 0]
  }

  vax_df[get(dose_col) >= max(vax_doses) & !is.na(get(dose_col)), paste0(vax_prefix, "_dose_", max(vax_doses)) := 1]
  vax_df[get(dose_col) <  max(vax_doses) & !is.na(get(dose_col)), paste0(vax_prefix, "_dose_", max(vax_doses)) := 0]

  ### column for missingness
  #vax_df[,missing:=0]
  #ax_df[is.na(get(dose_col)), missing := 1]
  # List for missing shapefiles
  missing_shapefiles <- list()

  ### Collapse 
    
  message(paste0("Collapsing ", vax_title))
  
  # Only keep those with vaccine information
  vax_df <- vax_df[!is.na(get(dose_col))]
  
  # Keep only useful Variables    
  keep_vars <- c("survey_name", "ihme_loc_id", "year_start", "nid", "geospatial_id", "point",
                 "lat", "long", "location_code", "shapefile", "year_start_cohort", "shapefile1","loc_name1","loc_code1",
                 "admin_level1", "pweight", "strata","psu",paste0(vax_prefix, "_dose_", vax_doses))
  vax_df <- subset(vax_df, select=keep_vars)
  
  # Rename columns
  rename_table <- data.table(rbind(c("survey_name", "source"), 
                                   c("ihme_loc_id", "country"), 
                                   c("year_start", "svy_year"), 
                                   c("nid", "svy_id"), 
                                   c("geospatial_id", "geo_id"), 
                                   c("lat", "latitude"),
                                   c("long", "longitude"), 
                                   c("year_start_cohort", "year"),
                                   c("shapefile1","admin1_shapefile"),
                                   c("loc_name1","admin1_name"),
                                   c("loc_code1","loc_code_admin1"),
                                   c("admin_level1","admin_level1"),
                                   c("point","point"),
                                   c("pweight","pweight"), 
                                   c("strata", "strata"),
                                   c("psu","psu")))
  
  names(rename_table) <- c("old", "new")
  setnames(vax_df, rename_table$old, rename_table$new)
  
  # recode shapefile
  vax_df[shapefile == "", shapefile := NA]
  vax_df[admin1_shapefile == "", admin1_shapefile := NA]
  vax_df[admin1_name == "", admin1_name := NA]
  vax_df[loc_code_admin1 == "", loc_code_admin1 := NA]
  vax_df[admin_level1 == "", admin_level1 := NA]
  vax_df[,admin1_shapefile := as.character(admin1_shapefile)]

  # Placeholder value to sum over
  vax_df[,N := 1]
  
  # Ensure column formats correct
  if (class(vax_df$latitude) != "numeric") vax_df[, latitude := as.numeric(as.character(latitude))]
  if (class(vax_df$longitude) != "numeric") vax_df[, longitude := as.numeric(as.character(longitude))]

  ### COLLAPSE POLYS AND POINTS ###########################################

  # Get shapefiles into character format (in case factor)
  vax_df$shapefile <- as.character(vax_df$shapefile)

  # Recode shapefile for points

  # Figure out what shapefiles are missing and exclude a priori
  # (Otherwise will break resample_polygons)

  shapefile_list <- unique(vax_df$shapefile[!is.na(vax_df$shapefile)])
  shapefile_dir <- paste0("FILEPATH")
  shapefiles_available <- list.files(shapefile_dir, ".shp$") %>% 
                            gsub(".shp", "", .)

  no_shapefile <- shapefile_list[!(shapefile_list %in% shapefiles_available)]

  if (length(no_shapefile > 0)) {
    message("The following shapefiles are missing from the shapefile directory. Associated data will be dropped:")
    message(paste(paste0("  ", no_shapefile, ".shp"), collapse = " \n"))
    vax_df <- vax_df[!(shapefile %in% no_shapefile)]
    missing_shapefiles[[vax]] <- no_shapefile
  }

  # First, check if there are any polygons

  if (nrow(vax_df[point == 0]) == 0) {

    #If so, skip ahead
    message("No polygon data found - moving to next cycle")
    df_pointpoly <- vax_df
    keep_vars <- c("svy_id", "source", "country", "point", "svy_year", 
                   "location_code", "shapefile", "year", 
                   "psu", "latitude", "longitude", paste0(vax_prefix, "_dose_", vax_doses), "N")
    df_pointpoly <- subset(df_pointpoly, select = names(df_pointpoly %in% keep_vars))

  } else {

    # Collapse the polygon data
    vax_df[ ,cluster_id := NULL]
    vax_df[, strata := as.numeric(strata)]

    DHS_subset <- copy(vax_df)
    DHS_subset <- DHS_subset[source=="MACRO_DHS" & !is.na(admin1_shapefile),]

    df_point <- vax_df[point == 1]
    df_poly <- vax_df[point == 0]

    pre_pweight_drop <- df_poly[,.N, .(svy_id)]
    names(pre_pweight_drop) <- c("nid","pre_pweight_drop")
    # Drop all rows without pweight
    warning(paste0('Dropping ', nrow(df_poly[is.na(pweight) & is.na(longitude),]), 
                   ' rows of ', nrow(df_poly), ' (', 
                   round(nrow(df_poly[is.na(pweight) & is.na(longitude),])/nrow(df_poly)*100,2),
                   '%) due to missing pweights'))
    
    df_poly <- df_poly[!is.na(longitude) | is.na(longitude) & !is.na(pweight),]
    post_pweight_drop <- df_poly[,.N, .(svy_id)]
    names(post_pweight_drop) <- c("nid","post_pweight_drop")
    post_pweight_drop <- merge(pre_pweight_drop, post_pweight_drop,by="nid", all.x=T)
    post_pweight_drop[is.na(post_pweight_drop), post_pweight_drop:=0]

    post_pweight_drop[,percent_pweight_missing := (pre_pweight_drop - post_pweight_drop)/pre_pweight_drop*100]

    message(paste0("Saving log file for rows dropped from having no pweight: ", log_dir, "pweight_drop.csv"))
    write.csv(post_pweight_drop, paste0(log_dir, "pweight_drop.csv"))

    # Collapse polygons -------------------------------------------------------------

    if (collapse_method == "survey") {
      # Collapse using the original method - survey-package based

      # create "lonely" column to indicate whether a polygon_year only has 1 individual
      by_vars <- c("svy_id", "source", "country", "point", 
                   "svy_year", "location_code", "shapefile", "year")
      df_poly[is.na(latitude) & !is.na(shapefile), lonely := lapply(.SD, length), .SDcols="N", by=by_vars]
      
      warning(paste0('Dropping ', nrow(df_poly[lonely==1,]), 
                     " out of ", uniqueN(df_poly[!is.na(lonely),], by=c("location_code", "shapefile", "year")),
                     " polygon-years that contain just 1 individual"))
      
      # Drop rows in which lonely is 1.  These are polygon_years that just have 1 individual 
      df_poly<- df_poly[lonely!=1,]

      ## Collapse each NID separately to most granular shapefiles
      collapse_each_nid <- function(this_nid) {
          
          message(paste0('Collapsing NID: ', this_nid))
          test_poly_data <- df_poly[svy_id==this_nid,]
          test_poly_data$strata <- 0
          
          # Check for missings
          if(length(test_poly_data$pweight[is.na(test_poly_data$pweight)])>0) {
              message(paste0(length(test_poly_data$pweight[is.na(test_poly_data$pweight)]), ' / ', length(test_poly_data$pweight), ' are missing pweight'))
              return(NULL)
          } else {
              collapse_polys <- function(x) {
                  setup_design(df = test_poly_data, var = x)
                  by_vars <- c('year', 'country', 'location_code', 'shapefile', 'source','svy_id', 'svy_year')
                  poly <- collapse_by(df = test_poly_data,
                                      var = x,
                                      by_vars = by_vars)
                  collapsed <- poly[, c(by_vars, 'mean', 'ss')]
                  names(collapsed)[names(collapsed)=='mean'] <- x
                  names(collapsed)[names(collapsed)=='ss'] <- 'N'
                  return(collapsed)
                  }
              doses <- c(paste0(vax_prefix, "_dose_", vax_doses))
              polys <- lapply(doses, collapse_wrapper)
              merged_polys <- Reduce(function(...) merge(..., all=T), polys)
              return(merged_polys)
          }
      }

       collapse_wrapper <- function(nid) {
        output <- tryCatch({
          collapse_each_nid(nid)
        }, error = function(e) NULL)
        return(output)
      }

      poly_nids <- unique(df_poly$svy_id)
      poly_list <- mclapply(poly_nids, collapse_wrapper, mc.cores=max(cores, 10))
      names(poly_list) <- poly_nids

      # Identify broken NIDs
      find_error <- function(p) {
        if ("error" %in% class(p)) {
          return(T)
        } else {
          return(F)
        }
      }
      
      broken_nids <- poly_nids[which(sapply(poly_list, find_error))]
      message(paste0("Unable to successfully collapse polys for the following NIDs:\n  ",
                     paste(broken_nids, collapse = "\n  "),
                     "\nThese NIDs will be dropped..."))
      
      # Drop broken NIDs
      poly_list<-poly_list[which(!(names(poly_list) %in% broken_nids))]

      df_poly <- rbindlist(poly_nids)

      for (var in paste0(vax_prefix, "_dose_", vax_doses)) {
        df_poly[, eval(var) := get(var)*N]
      }

      # Mark as poly 
      df_poly[, point := 0]

    } else if (collapse_method == "kish") {

      # Collapse using the Kish approximation to generate an effective sample size (N)
      # Keeps the original sample size in N_obs

      sum_cols <- c(paste0(vax_prefix, "_dose_", vax_doses))

      df_poly <-   df_poly[, c(list(N_obs = sum(N),
                                    N = sum(pweight)^2/sum(pweight^2)),
                               lapply(.SD, FUN = weighted.mean, w = pweight)), 
                             by = list(svy_id, source, country, point, svy_year, 
                                    location_code, shapefile, year),
                             .SDcols = sum_cols]

      # Convert these back to counts
      df_poly[, (sum_cols) := lapply(.SD, '*', N), .SDcols = sum_cols]                            

    }

    # Collapse the point data to clusters
    sum_cols <- c("N", paste0(vax_prefix, "_dose_", vax_doses))
    df_point <- df_point[, lapply(.SD, sum), 
                           by = .(svy_id, source, country, point, svy_year, 
                                  location_code, shapefile, year, pweight, 
                                  strata, psu, latitude, longitude),
                           .SDcols = sum_cols]

    # Harmonize names & combine
    drop_vars <- c("strata", "pweight")
    df_point <- subset(df_point, select = !(names(df_point) %in% drop_vars))

    df_pointpoly <- rbind(df_point, df_poly, fill = T)

  }

  # Add a weights column - this will be resampled in resample_polygons
  df_pointpoly$weight <- 1

  ### Final renaming stuff
  setnames(df_pointpoly, "svy_year", "original_year")

  df_pointpoly[,country:=tstrsplit(country,"_")[1]]

  ### Saving DHS files for plotting

  # Create a DHS subset if applicable
  if (collapse_dhs == T) {

    #####subset to just DHS for validation
    ### Subsetting to 12-23 months
    collapse_each_nid_admin1 <- function(this_nid) {
          
        message(paste0('Collapsing NID: ', this_nid))
        test_poly_data <- DHS_subset[svy_id==this_nid,]
        test_poly_data$strata <- 0
        
        # Check for missings
        if(length(test_poly_data$pweight[is.na(test_poly_data$pweight)])>0) {
          message(paste0(length(test_poly_data$pweight[is.na(test_poly_data$pweight)]), ' / ', length(test_poly_data$pweight), ' are missing pweight'))
          return(NULL)
        } else {
          collapse_polys <- function(x) {
            setup_design(df = test_poly_data, var = x)
            by_vars <- c('year', 'country', 'loc_code_admin1', 'admin1_shapefile', 'source','svy_id', 'svy_year')
            poly <- collapse_by(df = test_poly_data,
                                var = x,
                                by_vars = by_vars)
            collapsed <- poly[, c(by_vars, 'mean', 'ss')]
            names(collapsed)[names(collapsed)=='mean'] <- x
            names(collapsed)[names(collapsed)=='ss'] <- 'N'
            return(collapsed)
          }

          polys <- c(paste0(vax_prefix, "_dose_", vax_doses))
          polys <- lapply(polys, collapse_polys)
          merged_polys <- Reduce(function(...) merge(..., all=T), polys)
          return(merged_polys)
        }
      }
      #### Problem surveys with merging, drop for now
      bad_svys <- c(20120, 218593)
      DHS_subset <- DHS_subset[!(svy_id %in% bad_svys),]
      poly_nids_dhs <- unique(DHS_subset$svy_id)
      poly_nids_dhs <- mclapply(poly_nids_dhs, collapse_each_nid_admin1,mc.cores=10)

      DHS_sub<- rbindlist(poly_nids_dhs)

      fast_shapefiles <- TRUE

      ## Create DHS shapefile spdf
      df_shape_loc <- unique(DHS_sub[, c("admin1_shapefile", "loc_code_admin1")])
      poly_list <- pull_polys_in_parallel(shape_loc_list = df_shape_loc,
                                          shapefile_col = "admin1_shapefile",
                                          location_code_col = "loc_code_admin1",
                                          cores = cores,
                                          fast_shapefiles = fast_shapefiles)
      DHS_polys <- poly_list$poly_shapes_all

      for (var in paste0(vax_prefix, "_dose_", vax_doses)) {
        DHS_sub[, eval(var) := get(var)*N]
      }

      ## Fixes for name and vars for DHS_MBG compare
      names <- names(DHS_sub)
      names[7]<- "nid"
      names[5] <-"shapefile"
      names[4] <- "location_code"
      names[3] <- "iso3"
      names(DHS_sub) <- names
      DHS_sub[,outcome:=dpt_dose_3/N]


      message(paste0("Saving cleaned data for ", vax_title, " here: "))
      message(paste0("  ", vaccine_cleaned_file_prefix, vax_prefix, "_DHS_subset.csv"))
      write.csv(DHS_sub, paste0(vaccine_cleaned_file_prefix, vax_prefix, "_DHS_subset.csv"))

      ### Shapefile SPDF
      message(paste0("Saving DHS SPDF for ", vax_title, " here: "))
      message(paste0("  ", vaccine_cleaned_file_prefix, vax_prefix, "_DHS_polys.Rds"))
      saveRDS(DHS_polys, paste0(vaccine_cleaned_file_prefix, vax_prefix, "_DHS_polys.Rds"))

    }

  ### Save an RDS file  of the cleaned data
  message(paste0("Saving cleaned data for ", vax_title, " here: "))
  message(paste0("  ", vaccine_cleaned_file_prefix, vax_prefix, ".rds"))
  saveRDS(df_pointpoly, paste0(vaccine_cleaned_file_prefix, vax_prefix, ".rds"))
    
} # close (for vax in vaccines)
