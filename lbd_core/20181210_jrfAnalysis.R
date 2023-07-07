# HEADER ------------------------------------------------------------------
# Author: USERNAME
# Date: DATE
# Project: Custom project - JRF admin analysis
# Purpose: Combine all geo-matched JRF admins into a single big spdf and save
#************************************************************************** 

###########################################################################
# SETUP 
###########################################################################

## clear environment
rm(list=ls())

## Set repo location and indicator group
user               <- Sys.info()['user']
core_repo          <- sprintf('FILEPATH')
indic_repo         <- sprintf('FILEPATH')

## sort some directory stuff
setwd(core_repo)
commondir      <- sprintf('FILEPATH')
package_list <- c(t(read.csv(sprintf('FILEPATH',commondir),header=FALSE)))

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0('FILEPATH'))
mbg_setup(package_list = package_list, repos = core_repo)

# Custom load indicator-specific functions
source(paste0('FILEPATH'))

###########################################################################
# OPTIONS 
###########################################################################

year <- 2017
jrf_file <- "FILEPATH"
date_stamp <- gsub("-|:| ","_",Sys.time())
adm0_list <- get_gaul_codes("africa") #adm0 codes to match
vax <- "DTP3" # in jrf parlance

output_dir <- paste0("FILEPATH")
dir.create(output_dir, recursive = T)

###########################################################################
# READ JRF FILE AND FORMAT
###########################################################################

# Archive JRF lookup file
file.copy(jrf_file, paste0(output_dir, basename(jrf_file)))

jrf_match <- fread(jrf_file)
setnames(jrf_match, 
         c("NAME.who_adm1", "NAME.who", "NAME.shpfile"),
         c("who_name_adm1", "who_name", "shpfile_name"))

# Recode NAs
jrf_match[who_coverage == -2222, who_coverage := NA]

# Subset to only countries in adm0_list
# Note odd syntax currently needed because get_adm0_codes() is not properly vectorized
# ideally this would be as simple as jrf_match[, ADM0_CODE := get_adm0_codes(iso3)]

adm0_table <- data.table(iso3 = unique(jrf_match$iso3))
adm0_table[, ad0 := sapply(iso3, function(x) suppressMessages(get_adm0_codes(x)))]
adm0_table$ad0[sapply(adm0_table$ad0, function(x) (length(x)==0))] <- NA
adm0_table$ADM0_CODE <- unlist(adm0_table$ad0)
adm0_table[, ad0 := NULL]
jrf_match <- merge(jrf_match, adm0_table, by="iso3", all.x=T, all.y=F)

###########################################################################
# PULL POLYGONS
###########################################################################

# Generate a list of shapefiles
shapefile_list <- unique(subset(jrf_match, !is.na(shapefile), select = c(shapefile)))

poly_list <- lapply(shapefile_list$shapefile, function(shpfile, 
                                                       jrf_df = jrf_match,
                                                       yr = year,
                                                       vx = vax) {
  
  message("\n########################################")
  message(shpfile)
  
  # Figure out where the shapefile lives
  df_shapefile_lookup <- fread("FILEPATH")
  
  # Check to make sure that there's an entry in the lookup table
  if (!(shpfile %in% df_shapefile_lookup$shp)) {
    stop(paste0("Need to add a row to lookup table for this shapefile: ", shpfile))
  }
  
  # Find the file path 
  shpfile_path <- as.character(df_shapefile_lookup[shp == shpfile, shapefile_directory])
  
  # Load the shapefile using fast loading from sf library
  the_shp <- sf::st_read(paste0(shpfile_path, shpfile, ".shp"))
  
  # Merge JRF data and harmonize names ------------------------------------
  jrf_subset <- subset(jrf_df, shapefile == shpfile & vaccine == vax & !is.na(location_code))
  
  # Grab the field name
  field_nm <- unique(jrf_subset$field_name)
  field_nm <- as.character(field_nm[!is.na(field_nm)])
  
  # Rename the field name and set up for merge
  setnames(the_shp, field_nm, "the_field_name")
  jrf_to_merge <- subset(jrf_subset, 
                         !is.na(location_code) & !is.na(who_coverage),
                         select = c("who_coverage", "location_code"))
  
  setnames(jrf_to_merge, c("who_coverage", "location_code"), c("jrf_coverage", "the_field_name"))
  jrf_to_merge[, jrf_coverage := jrf_coverage / 100]
  
  # merge on JRF coverage and set names
  the_shp <- merge(the_shp, jrf_to_merge, by="the_field_name")
  the_shp$original_field_name <- field_nm
  the_shp$original_shpfile <- shpfile
  the_shp$vaccine <- vx
  setnames(the_shp, "the_field_name", "location_code")
  
  the_shp <- subset(the_shp, select = c("location_code", "original_field_name", "original_shpfile", "jrf_coverage", "vaccine"))
  
})

names(poly_list) <- shapefile_list$shapefile

# Adjust CRS of each sf object to conform with gadm's
main_crs <- sf::st_crs(poly_list[["gadm_36_ad2"]])

poly_list <- lapply(poly_list, function(x, mcrs = main_crs) {
  the_crs <- sf::st_crs(x)
  if (is.na(the_crs)) {
    # Unprojected, so go ahead and add gadm crs (unprojected too)
    sf::st_crs(x) <- mcrs
  } else if (the_crs != mcrs & !is.na(the_crs)) {
    x <- sf::st_transform(x, mcrs)
  }
  return(x)
})

# Merge everything together and write to disk
main_jrf_shapefile <- Reduce(rbind, poly_list)

# Delete if already present
suppressWarnings(file.remove(paste0(output_dir, "who_jrf_shapefile_", vax, ".shp")))

sf::st_write(obj = main_jrf_shapefile, 
             dsn = paste0(output_dir, "who_jrf_shapefile_", vax, ".shp"))
