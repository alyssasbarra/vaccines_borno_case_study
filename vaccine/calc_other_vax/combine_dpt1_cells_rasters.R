

## clear environment
rm(list=ls())

## Set repo location and indicator group
user               <- 'USERNAME'
core_repo          <- sprintf('FILEPATH')
indic_repo         <- sprintf('FILEPATH')

## sort some directory stuff
commondir      <- sprintf('FILEPATH')
package_list <- c(t(read.csv(sprintf('FILEPATH'),header=FALSE)))

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0('FILEPATH'))
mbg_setup(package_list = package_list, repos = core_repo)

# Custom load indicator-specific functions
source(paste0('FILEPATH'))

## Script-specific stuff begins here ##########################################
indicator_group <- 'vaccine'

# Load objects from qsub
load_from_parallelize() # run_date, indicator, vaccine, indicator_group, reg

reg <- region

# Define directories
dpt3_dir <- paste0("FILEPATH")
dpt12_dir <- paste0("FILEPATH")
dpt1_dir <- paste0("FILEPATH")

raking_shapefile_version <- modeling_shapefile_version <- shapefile_version
## LOAD CELL PREDS ######################################################################################

# First, load in the cell_preds for the base vaccine
message("Loading Data...")
load(paste0(dpt3_dir, "dpt3_cov_raked_cell_draws_eb_bin0_", reg,"_0.RData"))

dpt3_cell_pred <- copy(raked_cell_pred)


# Now, load in the cell_preds for the ratios
message("Loading Data...")
load(paste0(dpt12_dir, "dpt12_cond_raked_cell_draws_eb_bin0_", reg, "_0.RData"))

dpt12_cell_pred <- copy(raked_cell_pred)


## COMBINE CELL PREDS ######################################################################################

cell_pred <- ( dpt12_cell_pred * (1 - dpt3_cell_pred) ) + dpt3_cell_pred



## SAVE CELL PRED IN NEW RUN_DATE FOLDER ######################################################################################

## save raked cell preds
save(cell_pred, file = paste0(
  dpt1_dir, "dpt1_cov_raked_cell_draws_eb_bin0_", reg, "_0.RData"
))

#load(paste0(dpt1_dir, "dpt1_cov_raked_cell_draws_eb_bin0_", reg, "_0.RData"))


## SAVE THE RESULTS #####################################################################
message("Saving results...")
# make and save summaries

save_cell_pred_summary <- function(summstat, raked, ...) {
  message(paste0("Making summmary raster for: ", summstat, " (", raked, ")"))
  
  if (raked == "unraked") {
    cpred <- "cell_pred"
    mask_raster <- "simple_raster"
  }
  if (raked == "raked") {
    cpred <- "cell_pred"
    mask_raster <- "raked_simple_raster"
  }
  if (raked == "raked_c") {
    cpred <- "raked_cell_pred_c"
    load(paste0(sharedir, "/output/", run_date, "/", indicator, "_raked_c_cell_draws_eb_bin0_", reg, "_0.RData" ))
    mask_raster <- "raked_simple_raster"
  }
  ras <- make_cell_pred_summary(
    draw_level_cell_pred = get(cpred),
    mask = get(mask_raster),
    return_as_raster = TRUE,
    summary_stat = summstat,
    ...
  )
  save_post_est(ras,'raster',paste0(reg, ifelse(raked == "raked", "_raked", ifelse(raked == 'raked_c', '_raked_c', '')), '_', summstat, '_raster'))
}

# Do this as lapply to not fill up memory in global env with big obs
rake_list <- 'raked'


summ_list <- expand.grid(summstats[summstats != "p_below"], rake_list)

# lapply(1:nrow(summ_list), function(i) {
#   summstat <- as.character(summ_list[i, 1])
#   raked <- as.character(summ_list[i, 2])
#   save_cell_pred_summary(summstat, raked)
# })

indicator <- 'dpt1_cov'

  ## Get the simple and new_simple rasters prepped up for us
  setwd("FILEPATH")
  load(paste0(run_date, "_bin0_", reg, "_0.RData"))
  # setwd("FILEPATH")
  # load(paste0("dpt3_cov_cell_draws_eb_bin0_", reg,"_0.RData"))
  
  ## Take out the objects from the list that actually matters to us:
  # simple_raster <- raster_outputs[["simple_raster"]]
  # new_simple_raster <- raster_outputs[["new_simple_raster"]]
  
  simple_polygon <- simple_raster
  new_simple_polygon <- simple_raster
  
  pixel_id <- pixel_id <- seegSDM:::notMissingIdx(simple_raster)
  
  raked_simple_raster <- simple_raster

## Can't pass additional params in the above framework, so will code by hand here
for (r in rake_list) {
  if ("mean" %in% summstats) {
    save_cell_pred_summary(
      summstat = "mean",
      raked = r
    )
  }
  if ("cirange" %in% summstats) {
    save_cell_pred_summary(
      summstat = "cirange",
      raked = r
  )
  }
  if ("lower" %in% summstats) {
    save_cell_pred_summary(
      summstat = "lower",
      raked = r
    )
  }
  if ("upper" %in% summstats) {
    save_cell_pred_summary(
      summstat = "upper",
      raked = r
    )
  }
  if ("cfb" %in% summstats) {
    save_cell_pred_summary(
      summstat = "cfb",
      raked = r
    )

  }
}


# All done
message(paste0("Done with ", reg, "!"))








