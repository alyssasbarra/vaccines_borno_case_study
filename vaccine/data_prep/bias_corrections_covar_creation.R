

rm(list=ls())

# Set up environment and directories
## Set repo location and indicator group
user            <- 'USERNAME'
repo            <- sprintf('FILEPATH')

## drive locations
root           <- ifelse(Sys.info()[1]=='Windows', 'FILEPATH', 'FILEPATH')
commondir      <- sprintf('FILEPATH')
package_list <- c(t(read.csv(sprintf('FILEPATH',commondir),header=FALSE)))
setwd(repo)
core_repo <- repo

for (p in package_list) {
  try(library(p, character.only = T))
}

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0('FILEPATH'))

mbg_setup(package_list = package_list, repos = repo)


vaccines = c("vacc_bcg", "vacc_polio3", "vacc_dpt3","vacc_mcv1")
year_ids=c(2000:2020)
gbd_date = "DATE"
gbd_date_mbg_format = "DATE"
gbd_yr = 2020

# save_biascorr_admin_est <- function(vaccines, 
#                                     year_ids, 
#                                     gbd_date = "DATE",
#                                     gbd_date_mbg_format = "DATE",
#                                     gbd_yr = 2020) {
  
  # -------------------------------------------------------------------------------------------------------
  # Set up directories and files --------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------
  message("\nSetting up directories and files")
  
  std_cov_dir <- "FILEPATH"
  bias_corr_file <- paste0("FILEPATH")
  
  # Load in bias correction file
  df_biascorr <- readRDS(bias_corr_file)
  
  # Subset to bias corrected admin data and vaccines of interest only
  df_biascorr <- subset(df_biascorr, me_name %in% vaccines & cv_admin == 1)
  
  # Drop any outliered values
  if (nrow(subset(df_biascorr, !is.na(cv_outlier))) > 0) {
    message(paste0("Dropping ", nrow(subset(df_biascorr, !is.na(cv_outlier))), 
                   " outliered rows out of ", nrow(df_biascorr), " total rows."))
    df_biascorr <- subset(df_biascorr, is.na(cv_outlier))
  }
  
  # -------------------------------------------------------------------------------------------------------
  # Check to see which rows are misisng --------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------
  
  # Create a list of all country-year combinations
  ihme_loc_ids <- unique(df_biascorr$ihme_loc_id)
  df_year_loc_vax <- expand.grid(ihme_loc_ids, year_ids, vaccines, stringsAsFactors = F) %>% as.data.table
  names(df_year_loc_vax) <- c("ihme_loc_id", "year_id", "me_name")
  
  # Merge on to df_biascorr to create empty rows where needed
  df_biascorr <- merge(df_biascorr, df_year_loc_vax, all.x = T, all.y = T, by = c("year_id", "ihme_loc_id", "me_name"))
  
  
  # -------------------------------------------------------------------------------------------------------
  # Fill in missing data rows via linear interpolation or forwards / backwards fill where needed ----------
  # -------------------------------------------------------------------------------------------------------
  
  # Determine whether there are any missing data rows
  n_missing <- nrow(subset(df_biascorr, year_id %in% year_ids & is.na(data)))
  
  if (n_missing > 0) {
    message(paste0("\nMissing values found for ", n_missing, " rows - attempting to fill via linear interpolation"))
    
    # A function to fill in missing data via linear interporation or back-fill / forwards-fill
    
    interpolate_missing <- function(loc_id, yr, vax) {
      
      # Find last non-NA value
      last_row <- df_biascorr[ihme_loc_id == loc_id & me_name == vax & year_id < yr & !is.na(data)]
      last_row <- suppressWarnings(last_row[year_id == max(year_id)])
      last_year <- last_row$year_id
      last_data <- last_row$data
      
      if(length(last_data)>1) {
        last_data = last_data[1]
      #  message(paste('>1 last data option for loc_id and year:', loc_id, yr, sep=' '))
      }
      
      # Find next non-NA value
      next_row <- df_biascorr[ihme_loc_id == loc_id & me_name == vax & year_id > yr & !is.na(data)]
      next_row <- suppressWarnings(next_row[year_id == min(year_id)])
      next_year <- next_row$year_id
      next_data <- next_row$data
      
      if(length(next_data>1)) {
        next_data = next_data[1]
       # message(paste('>1 next data option for loc_id and year:', loc_id, yr, sep=' '))
      }
      
      # If not able to locate these, fill in with the last or next non-missing value
      # This is not ideal but follows what is routinely done with covariates 
      
      if (length(next_data) == 0) {
        data_val <- last_data
        notes <- data.table(ihme_loc_id = loc_id,
                            year = yr,
                            me_name = vax,
                            created_value = data_val,
                            method = paste0("Filled with last non-NA value from ", last_row$ihme_loc_id, " - ", last_year))
      } else if (length(last_data) == 0) {
        data_val <- next_data
        notes <- data.table(ihme_loc_id = loc_id,
                            year = yr,
                            me_name = vax,
                            created_value = data_val,
                            method = paste0("Filled with next non-NA value from ", next_row$ihme_loc_id, " - ", next_year))
      } else {
        # Calculate linear interpolation value
        annual_change <- (next_data - last_data) / (next_year - last_year)
        data_val <- last_data + annual_change * (yr - last_year)
        notes <- data.table(ihme_loc_id = loc_id,
                            year = yr,
                            me_name = vax,
                            created_value = data_val,
                            method = paste0("Linear interpolation between ", last_row$ihme_loc_id, " - ", last_year, 
                                            " and ", next_row$ihme_loc_id, " - ", next_year))
        
      }
      
      return(list(data_val = data_val,
                  notes = notes))
    }
    
    # Generate an ID variable and an empty interpolated data column
    df_biascorr[, id := .I]
    df_biascorr[, data_interp := NULL]
    df_biascorr[, data_interp := as.numeric()]
    
    # Initialize a notes list to capture all of the notes about how things were interpolated
    notes_list <- list()
    
    # Do the interpolation row-by-row
    missing_ids <- df_biascorr[is.na(data) & year_id %in% year_ids, id]
    for (i in 1:length(missing_ids)) {
      print(i)
      interpolation <- interpolate_missing(loc_id = df_biascorr[id == missing_ids[i], ihme_loc_id],
                                           yr = df_biascorr[id == missing_ids[i], year_id],
                                           vax = df_biascorr[id == missing_ids[i], me_name])
      if(length(unique(interpolation$data_val))>0){
        df_biascorr[id == missing_ids[i],]$data_interp <-  unique(interpolation$data_val)
      }else{
        df_biascorr[id == missing_ids[i],]$data_interp <- NA
      }
      notes_list[[i]] <- interpolation$notes
    }
    
    # Fill in NAs with interpolated values and create a table of notes to save later
    df_biascorr[is.na(data) & !is.na(data_interp), data := data_interp]
    interpolation_notes <- rbindlist(notes_list) %>% setorderv(., c("me_name", "ihme_loc_id", "year"))
    rm(notes_list)
    
    # check to ensure that this worked 
    if (nrow(df_biascorr[is.na(data)]) == 0) {
      message(paste0("Succesfully filled ", n_missing, " rows.  See interpolation.csv in covariate directory for details")) 
    } else {
      warning(paste0("Unable to fill ", nrow(df_biascorr[is.na(data)]), " rows via interpolation or back/forward filling."))
      message(df_biascorr[is.na(data)]$ihme_loc_id)
    }
    
  } # end `if (n_missing > 0)
  
  #save an intermediate file if something goes wrong in raster creation
  saveRDS(df_biascorr,paste0('FILEPATH'))
  # -------------------------------------------------------------------------------------------------------
  # Assign data spatially to make a covariate in the proper format-----------------------------------------
  # -------------------------------------------------------------------------------------------------------
  
  message('\nLoading global raster template for GBD covs')
  template <- raster(paste0("FILEPATH"))
  
  # Load the raking shapefile
  shapefile_path <- "FILEPATH"
  
  location_metadata <- get_location_code_mapping( shapefile_version <- 'DATE')
  location_metadata <- location_metadata[,c("GAUL_CODE", "loc_id", "ihme_lc_id")]
  setnames(location_metadata, c("GAUL_CODE", "ihme_lc_id"), c("ADM0_CODE", "ihme_loc_id"))
  
  ad0_shape <- sf::st_read(shapefile_path, stringsAsFactors = F)
  
  
  regs<- get_location_code_mapping(shapefile_version = 'DATE')
  
  # Merge on loc_ids for making global raster
  ad0_shape <- merge(ad0_shape, regs, by.x = "ADM0_CODE", by.y="ADM_CODE")
  
  # Merge loc_id with bias corr data
  df_biascorr <- merge(df_biascorr, regs, by.x = c("ihme_loc_id"), by.y = c("ihme_lc_id"))
  
  # -------------------------------------------------------------------------------------------------------
  # Define a function to quickly make year rasters from data ----------------------------------------------
  # -------------------------------------------------------------------------------------------------------
  
  make_year_raster <- function(ad0_shape, template, df_biascorr, year, vax) {
    
    # Subset to year of interest an dmerge with admin0 shape
    df_yr <- subset(df_biascorr, year_id == year & me_name == vax, select = c("loc_id", "data", "year_id"))
    ad0_shape_yr <- merge(ad0_shape, df_yr)
    
    # Rasterize using loc_id
    ad0_ras <- fasterize::fasterize(ad0_shape_yr, template, "data")
    
    # Ensure that extents match
    ad0_ras <- extend(ad0_ras, y = template, value = NA)
    
    # Check to ensure that rasters match
    if (!compareRaster(ad0_ras, template)) {
      stop("Your admin raster resolution, extent, or projection does not match your template - please check this!")
    }    
    return(ad0_ras)
  }
  
  # -------------------------------------------------------------------------------------------------------
  # Loop over vaccines & write to standard covariate directory --------------------------------------------
  # -------------------------------------------------------------------------------------------------------
  
  for (vax in vaccines) {
    
    message(paste0("\nWriting outputs for ", vax))
    
    message("  Generating raster brick")
    # Generate a raster brick for the vaccine of interest
    ras_brick <- lapply(year_ids, function(yyy) make_year_raster(ad0_shape = ad0_shape,
                                                                 template = template,
                                                                 df_biascorr = df_biascorr,
                                                                 year = yyy,
                                                                 vax = vax)) %>% brick
    
    output_dir <- paste0(std_cov_dir, "bias_corr_", vax, '/mean/',gbd_date_mbg_format)
    dir.create(output_dir, recursive = T)
    output_dir <- paste0(output_dir,'/1y/')
    dir.create(output_dir, recursive = T, showWarnings = F)
    
    message("  Writing .tif for year:")
    # Write raster layer for each year, raked and unraked
    for(i in 1:length(year_ids)) {
      subset_r_brick <- ras_brick[[i]]
      year <- year_ids[i]
      message(paste0("  -- ", year))
      writeRaster(subset_r_brick,
                  filename = paste0(output_dir, "bias_corr_", vax, "_mean_1y_", year, "_00_00.tif"),
                  format = "GTiff",
                  overwrite = TRUE)
    }
    
    # Write an .rds file containing the bias corrected data table
    message("  Writing .rds file of bias-corrected admin data used for covariate")
    df_biascorr_output <- subset(df_biascorr, 
                                 me_name == vax,
                                 select = c("ihme_loc_id", "loc_id", "year_id", "me_name", "data", "data_interp"))
    df_biascorr_output[!is.na(data_interp), interpolated := T]
    df_biascorr_output[is.na(data_interp), interpolated := F]
    df_biascorr_output[, data_interp := NULL]
    
    saveRDS(df_biascorr_output, file = paste0(output_dir, vax, "_bias_corr_values.rds"))
    
    # Write various additional notes to output directory  ---------------------------------------------------
    message("  Writing notes:")
    
    # Interpolation notes
    message("  -- interpolation_notes.csv")
    fwrite(interpolation_notes[me_name == vax,], file = paste0(output_dir, "interpolation_notes.csv"))
    
    # Save a readme file 
    message("  -- readme.txt")
    readme_file <- paste0(output_dir, "readme.txt")
    fileConn <- file(readme_file)
    writeLines(c("## MBG COVARIATE README #############################",
                 "",
                 paste0("Covariate generated by ", Sys.info()['user'], " on ", Sys.time(), "."),
                 "",
                 paste0("Covariate name: bias_corr_", vax),
                 "",
                 paste0("This represents bias-corrected administrative data for ", vax, " from GBD"),
                 "",
                 "####################################################"),
               fileConn)
    close(fileConn)
    
    # Find missing country_years
    all_country_years <- expand.grid(loc_id = unique(ad0_shape$loc_id),
                                     year = year_ids,
                                     stringsAsFactors = F) %>% as.data.table
    all_country_years <- merge(all_country_years, location_metadata) %>% 
      subset(., select = c("loc_id", "ihme_loc_id", "year"))
    
    all_country_years <- merge(all_country_years, subset(df_biascorr, select = c("loc_id", "year_id", "data")),
                               all.x = T, all.y = F, by.x = c("loc_id", "year"), by.y = c("loc_id", "year_id"))
    
    missing_country_years <- all_country_years[is.na(data),]
    
    if (nrow(missing_country_years) < 0) {
      warning(paste0("You have ", nrow(missing_country_years), " country-years with no bias corrected administrative estimate available.\n",
                     "Please see the missing_country_years.csv file in the output directory for details"))
      message("  -- missing_country_years.csv")
      fwrite(missing_country_years, file = paste0(output_dir, "missing_country_years.csv"))
    }
    message(paste0("  Completed ", vax, " succesfully.\n",
                   "  Output directory: ", output_dir))
  } # Close vaccines loop  
# }



# save_biascorr_admin_est(vaccines = c("vacc_bcg", "vacc_polio3", "vacc_dpt3","vacc_mcv1"), 
#                         year_ids=c(2000:2020), 
#                         gbd_date = "DATE",
#                         gbd_date_mbg_format = "DATE",
#                         gbd_yr = 2020)





