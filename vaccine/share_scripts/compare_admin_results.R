# Script to perform custom aggregations, then compare to JRF data

# -------------------------------------------------------------------------------------------------
# Setup -------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

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
package_list <- c(t(read.csv(sprintf('FILEPATH',commondir),header=FALSE)))

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0('FILEPATH'))
mbg_setup(package_list = package_list, repos = core_repo)

# Custom load indicator-specific functions
source(paste0('FILEPATH'))

#####################################################################################################

# -------------------------------------------------------------------------------------------------
# Create table to match ---------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

year <- 2017
indicator <- "dpt3_cov"
indicator_group <- "vaccine"
run_date <- "RUN_DATE"
raked <- TRUE
pop_measure <- "a0004t"
overwrite <- TRUE
age <- 0
holdout <- 0
crop_shapefile <- T
use_lookup <- T
jrf_file <- "FILEPATH"
date_stamp <- gsub("-|:| ","_",Sys.time())

plot_output_dir <- paste0("FILEPATH")
dir.create(plot_output_dir, recursive = T)

# Archive JRF lookup file
file.copy(jrf_file, paste0(plot_output_dir, basename(jrf_file)))

# Set up a readme
fileConn <- file(paste0(plot_output_dir, "readme.txt"))
writeLines(c(paste0("Indicator: ", indicator),
             paste0("Indicator Group: ", indicator_group),
             paste0("Run date for indicator: ", run_date),
             paste0("JRF file path: ", jrf_file), 
             paste0("User: ", Sys.info()["user"]),
             paste0("This comparison launched ", Sys.time())),
          fileConn)
close(fileConn)

jrf_match <- fread(jrf_file)
setnames(jrf_match, 
         c("NAME.who_adm1", "NAME.who", "NAME.shpfile"),
         c("who_name_adm1", "who_name", "shpfile_name"))

if (indicator == "dpt3_cov") {
  jrf_match <- subset(jrf_match, vaccine == "DTP3")
} else {
  stop("Need to define indicator / vaccine name relationship.")
}

countries_to_match <- unique(subset(jrf_match, ADM0_NAME != "", select = c("ADM0_NAME", "iso3", "field_name", "shapefile")))
countries_to_match[, gaul_code := sapply(iso3, function(x) suppressMessages(get_gaul_codes(x)))]
countries_to_match <- countries_to_match[!is.na(gaul_code) & gaul_code != 0,]
countries_to_match$gaul_code <- unlist(countries_to_match$gaul_code)

region_list <- get_output_regions(in_dir = paste0("FILEPATH"))

identify_region <- function(gaul_list, region_list) {
	
	regions <- c('vax_soas','vax_seas','vax_eaas','vax_caeu','vax_crbn','vax_ctam','vax_ansa','vax_trsa','vax_name','vax_cssa','vax_essa','vax_sssa','vax_wssa')

	reg_table <- lapply(regions, function(rr) {
		return(data.table(region = rr, gaul_code = suppressMessages(get_gaul_codes(rr))))
	}) %>% rbindlist

	gl <- data.table(gaul_code = gaul_list)
	gl <- merge(gl, reg_table, by = "gaul_code", all.x = T, all.y = F, sort = F)

	return(gl$region)	
}

countries_to_match[, region := identify_region(gaul_code, region_list)]
countries_to_match[, filename_addin := paste0("_", iso3)]

# Look only at countries with assigned regions
countries_to_match <- subset(countries_to_match, !is.na(region))

# And subset to those only with shapefiles
countries_to_match <- subset(countries_to_match, !is.na(shapefile))

# Set up subsetting codes
countries_to_match[shapefile %in% c("gadm_36_ad1", "gadm_36_ad2"), 
                   subset_field_codes := sapply(iso3, function(x) suppressMessages(get_adm0_codes(x, adm0_type = "gadm")))]

countries_to_match[shapefile %in% c("gadm_36_ad1", "gadm_36_ad2"), 
                   subset_field := "ADM0_CODE"]

countries_to_match[shapefile %in% c("lf_g2015_2014_1", "lf_g2015_2014_2"), 
                   subset_field_codes := sapply(iso3, function(x) suppressMessages(get_adm0_codes(x, adm0_type = "gaul")))]

countries_to_match[shapefile %in% c("lf_g2015_2014_1", "lf_g2015_2014_2"), 
                   subset_field := "ADM0_CODE"]

# Archive this file
fwrite(countries_to_match, file = paste0(plot_output_dir, "countries_to_match.csv"))                   

# -------------------------------------------------------------------------------------------------
# Run aggregation in parallel ---------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

agg_lv <- subset(countries_to_match, select = c("shapefile", "filename_addin", "region", "subset_field", "subset_field_codes", "field_name"))

agg_para_output <- parallelize(script = "custom_aggregate_results",
                               log_location = "sgeoutput",
                               lv_table = agg_lv,
                               save_objs = c("indicator", "indicator_group", "run_date", "raked", "pop_measure",
                                      			 "overwrite", "age", "holdout", "crop_shapefile", "use_lookup",
                                				     "core_repo"),
                               prefix = "agg",
                               slots = 8,
                               memory = 64,
                               script_dir = paste0(indic_repo, "share_scripts/"),
                               geo_nodes = F,
                               use_c2_nodes = T,
                               singularity = 'default')

monitor_jobs(agg_para_output, notification = "pushover", max_tries = 1)

# -------------------------------------------------------------------------------------------------
# Create plots for each country -------------------------------------------------------------------
plot_lv <- subset(countries_to_match, select = c("shapefile", "filename_addin", "subset_field", "subset_field_codes", "field_name"))

plot_para_output <- parallelize(script = "plot_admin_results",
                                log_location = "sgeoutput",
                                lv_table = plot_lv,
                                save_objs = c("indicator", "indicator_group", "run_date", "raked", "core_repo"),
                                prefix = "plot",
                                slots = 1,
                                memory = 1,
                                script_dir = paste0(indic_repo, "share_scripts/"),
                                geo_nodes = T,
                                use_c2_nodes = F,
                                singularity = 'default')

monitor_jobs(plot_para_output, notification = "pushover", max_tries = 1)

# -------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------------
# Conduct the overlap analysis --------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

all_matches <- lapply(1:nrow(countries_to_match), function(i) {

  country <- countries_to_match[i, Country]
  iso3 <- countries_to_match[i, iso3]
  field_name <- countries_to_match[i, field_name]
  shapefile <- countries_to_match[i, shapefile]
  filename_addin <- countries_to_match[i, filename_addin]
  the_iso3 <- countries_to_match[i, iso3]
  the_year <- year
  the_field_name <- countries_to_match[i, field_name]

  # Read in the file
  df_agg <- fread(paste0("FILEPATH"))

  df_merge <- merge(subset(df_agg, year == the_year),
                    subset(jrf_match, iso3 == the_iso3, select = c("NAME.who", 
                                                                   "who_coverage", 
                                                                   "location_code", 
                                                                   "shapefile")),
                    by.x = the_field_name, by.y = "location_code")

  setnames(df_merge, "NAME.who", "who_name")
  df_merge[, who_coverage := who_coverage / 100]

})

all_matches <- rbindlist(all_matches, fill = T)





