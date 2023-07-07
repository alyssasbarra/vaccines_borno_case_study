RUN_DATE
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
setwd(core_repo)
commondir      <- sprintf('FILEPATH')
package_list <- c(t(read.csv(sprintf('FILEPATH'),header=FALSE)))

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0('FILEPATH'))
mbg_setup(package_list = package_list, repos = core_repo)

# Custom load indicator-specific functions
source(paste0('FILEPATH'))

## Options ########################################

indicators <- c("dpt3_cov", "dpt1_cov", "dpt1_3_abs_dropout")
measures <- c("mean", "upper", "lower", "cfb")
indicator_group <- "vaccine"
run_date <- "RUN_DATE"
raked <- T
year_list <- c(2000:2016)

copy_table <- expand.grid(indicator = indicators,
                          measure = measures,
                          stringsAsFactors = F) %>%
							as.data.table

## Copy tifs #####################################

for (i in 1:nrow(copy_table)) {
	copy_tif_to_geoviz_dir(ind = copy_table[i, indicator],
	                       ig = indicator_group,
	                       measure = copy_table[i, measure],
	                       raked = TRUE,
	                       yl = year_list,
	                       rd = run_date)	
}

# Manually copy over the "diff" rasters
for (ind in indicators) {
	for (measure in measures[measures != "cfb"]) {

		out_dir <- paste0("FILEPATH")
		dir.create(out_dir, recursive = T)

	  # Define input dirs / files
	  in_dir <- paste0("FILEPATH")
	  in_file <- paste0("FILEPATH")

	  # Load raster and save with standard filename
	  ras <- raster(in_file)

		# Adjust so that values are scaled to 0-100 for percents
		ras <- ras*100

		# Replace NA values with -999999 for viz tool
		# ras[is.na(ras)] <- -999999
		NAvalue(ras) <- -999999

		ind_name <- ifelse(ind=="dpt1_3_abs_dropout", "dropout", ind)
    out_filename <- paste0(ind_name, "_", measure, "_diff")
    out_file <- paste0(out_dir, out_filename)

    writeRaster(ras,
                file = out_file,
                overwrite = TRUE,
                format='GTiff',
			          datatype='FLT4S', 
			          NAflag = -999999)
	}
}


## Copy admins #####################################

for (i in 1:nrow(copy_table)) {
	copy_admins_to_geoviz_dir(ind = copy_table[i, indicator],
			                      ig = indicator_group,
			                      measure = copy_table[i, measure],
			                      raked = TRUE,
			                      yl = year_list,
			                      rd = run_date)	
}

# Manually copy "diff" admin csvs

for (ind in indicators) {
	out_dir <- paste0("/FILEPATH")
	dir.create(out_dir, recursive = T)
	
	all_admins <- lapply(c(0,1,2), function(ad_level) {
									in_file <- paste0("FILEPATH")
									in_df <- fread(in_file)
									setnames(in_df, 
									         paste0("ADM", ad_level, "_CODE"),
									         "gaul_code")
									return(subset(in_df, select = c("gaul_code", measures[measures != "cfb"])))
								})

	all_admins <- rbindlist(all_admins)

	for (meas in measures[measures != "cfb"]) {
		meas_df <- subset(all_admins, select = c("gaul_code", meas))
		setnames(meas_df, meas, "value")
		meas_df <- meas_df[!is.na(value),]
		meas_df[, value := value*100] # Adjust so that values are scaled 0-100 for percent

		ind_name <- ifelse(ind=="dpt1_3_abs_dropout", "dropout", ind)
		
		fwrite(meas_df, file = paste0(out_dir, ind_name, "_", meas, "_diff.csv"))
	}
	
}

