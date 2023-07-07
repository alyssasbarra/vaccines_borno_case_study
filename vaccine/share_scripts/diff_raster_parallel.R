
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
package_list <- c(t(read.csv(sprintf('FILEPATH'),header=FALSE)))

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0('FILEPATH'))
mbg_setup(package_list = package_list, repos = core_repo)

# Custom load indicator-specific functions
source(paste0('FILEPATH'))

### OPTIONS DEFINED HERE ################################################################

## Options

# Pass: indicator, indicator_group, run_date, region, year_list, raked, summstats
load_from_parallelize()

message(indicator, " | ", region, " | ", run_date)

## Set up directories, etc.
rake_addin <- ifelse(raked, "_raked", "")

sharedir <- paste0("FILEPATH")

### GENERATE A DIFFERENCE CELL PRED #######################################################
## Load cell pred
message("Loading cell pred object")
cell_pred_file <- paste0(sharedir, indicator, rake_addin, "_cell_draws_eb_bin0_", region, "_0.RData")

load(cell_pred_file, verbose = TRUE)

if (raked) {
	cell_pred <- raked_cell_pred
	rm(raked_cell_pred)
}

# Check to ensure lengths conformable
if (!(dim(cell_pred)[1] %% length(year_list) == 0)) {
	stop("The length of the cell pred is not divisible by the number of years")
}

# Create a year index for cell pred rows
n_rows_per_year <- dim(cell_pred)[1] / length(year_list)
year_idx <- rep(year_list, each = n_rows_per_year)

# Split out the first and last matrices
first_year <- cell_pred[which(year_idx == min(year_list)),]
last_year  <- cell_pred[which(year_idx == max(year_list)),]

if (!identical(dim(first_year), dim(last_year))) {
	stop("First and last year cell_pred chunks are not identical")
}

# Free up some memory
rm(cell_pred); gc()

# Calculate a diff cell pred
diff_cell_pred <- last_year - first_year

# Free up more memory
rm(first_year); rm(last_year); gc()

### SUMMARIZE THE DIFFERENCE CELL PRED ####################################################

# Load simple raster
message("Loading simple raster...")
simple_polygon_list <- load_simple_polygon(gaul_list = get_gaul_codes(region), buffer = 0.4, subset_only = FALSE)
subset_shape   <- simple_polygon_list[[1]]
raster_list    <- build_simple_raster_pop(subset_shape)
simple_raster  <- raster_list[['simple_raster']]

for (ss in summstats) {
	message(paste0("Summarizing ", ss, "..."))
	ras <- make_cell_pred_summary(draw_level_cell_pred = diff_cell_pred,
	                              mask = simple_raster,
	                              return_as_raster = TRUE,
	                              summary_stat = ss)
	output_file <- paste0(sharedir, indicator, "_", region, rake_addin, 
	                      "_diff_", min(year_list),"-", max(year_list), "_", ss, "_raster.tif")

    writeRaster(ras,
	            file = output_file,
	            format='GTiff',
	            overwrite = TRUE)

}


message("Done with script")