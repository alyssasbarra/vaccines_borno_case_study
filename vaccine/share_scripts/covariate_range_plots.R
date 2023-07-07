
# clear memory
rm(list=ls())

### I: SETUP #######################################################################

repo <- 'FILEPATH'

## Set repo location and indicator group
user               <- 'USERNAME'
core_repo          <- sprintf('FILEPATH')
indic_repo         <- sprintf('FILEPATH')

root <- ifelse(Sys.info()[1]=="Windows", "FILEPATH", "FILEPATH")
j_root <- "FILEPATH"
## Load libraries and miscellaneous MBG project functions.
setwd(core_repo)
root <- ifelse(Sys.info()[1]=="Windows", "FILEPATH", "FILEPATH")
  package_lib    <- ifelse(grepl('geos', Sys.info()[4]),
                          paste0(j_root,'FILEPATH',
                          paste0(j_root,'FILEPATH'))
.libPaths(package_lib)                                  # Ensures packages look for dependencies here when called with library(). 
                                                        #    Necessary for seeg libraries.

source('FILEPATH')                   # Functions to run MBG model.
source('FILEPATH')                  # Functions to setup MBG run
source('FILEPATH')             # Functions to prep and transform 5*5 covariates
source('FILEPATH')                  # Other computational MBG-related functions.
source('FILEPATH')
source('FILEPATH')
source('FILEPATH')
source('FILEPATH')
source('FILEPATH')
source('FILEPATH')
source(paste0('FILEPATH'))
source(paste0("FILEPATH"))
source('FILEPATH')     # Using USERNAME's edit for now that can take temporally varying covariates

package_list <- c('survey', 'magrittr', 'foreign', 'rgeos', 'data.table',
                  'raster','rgdal','INLA','seegSDM','seegMBG','plyr','dplyr', 
                  'foreach', 'snow', 'parallel', 'ggplot2')
for(package in package_list) {
  library(package, lib.loc = package_lib, character.only=TRUE)
}

# Custom load indicator-specific functions
source(paste0('FILEPATH'))

# Setup run options ##########################################################

rd <- "RUN_DATE"
ind <- "dpt3_cov"
reg <-"cssa"
ig <- "vaccine"
stacker_names <- c("gam", "gbm", "lasso")
year_list <- c(2000:2016)
out_dir <- paste0("FILEPATH")
dir.create(out_dir, showWarnings = F)

for (reg in c("cssa", "name", "sssa", "essa", "wssa")) {

  message(paste0("Working on region ", reg, "..."))
  
  # Load stackers
  load(paste0("FILEPATH"))  
  
  stacker_list <- cov_list[which(names(cov_list) %in% stacker_names)]
  
  # Load data
  df <- fread(paste0("FILEPATH"))
  
  extract_values <- function(stacker, yr, stack_list = stacker_list, yl = year_list, the_df = df) {
    df_year <- subset(the_df, year == yr)
    if (nrow(df_year) == 0) return(NULL)
    ras_year <- stack_list[[stacker]][[which(year_list == yr)]]
     #vals <- raster::extract(ras_year, cbind(df_year$longitude, df_year$latitude))
    # range_vals <- range(vals, na.rm = T)
    # min_vals <- min(range_vals)
    # max_vals <- max(range_vals)
    
    # Figure out where we have data
    ras_idx <- matrix(1:(nrow(ras_year)*ncol(ras_year)), nrow=nrow(ras_year), ncol=ncol(ras_year))
    ras_idx <- raster(x = ras_idx)
    crs(ras_idx) <- crs(ras_year)
    extent(ras_idx) <- extent(ras_year)  
    
    vals_idx <- na.omit(raster::extract(ras_idx, cbind(df_year$longitude, df_year$latitude)))
    cells_with_data <- as.vector(as.matrix(ras_year))[vals_idx]
    if (length(na.omit(cells_with_data)) == 0) return(NULL)
    cells_without_data <- as.vector(as.matrix(ras_year))[-na.omit(vals_idx)]
    
    range_with_data <- range(cells_with_data, na.rm = T)
    
    df_with_data <- data.table(val = na.omit(cells_with_data), 
                               data = "with_data",
                               stacker = stacker,
                               year = yr)
    
    df_without_data <- data.table(val = na.omit(cells_without_data),
                                  data = "without_data",
                                  stacker = stacker,
                                  year = yr)
    
    df_return <- rbind(df_with_data, df_without_data, use.names=T)
   
  }
  
  out_list <- list()
  
  for (ss in stacker_names) {
    for (yy in year_list) {
      #message(ss, "_", yy)
      out_list[[paste0(ss, "_", yy)]] <- extract_values(stacker = ss, yr = yy)
    }
  }
  
  out_df <- rbindlist(out_list)
  out_df[stacker == "gam", stacker := "GAM"]
  out_df[stacker == "gbm", stacker := "GBM"]
  out_df[stacker == "lasso", stacker := "Lasso"]
  out_df[data == "with_data", data := "Observed coverage data"]
  out_df[data == "without_data", data := "No observed coverage data"]
  
  gg_all <- ggplot(data = out_df, aes(x = val, color = data)) +
    geom_density() +
    theme_bw() +
    facet_wrap(~stacker) +
    labs(x = "Value", y = "Density", color = "Pixel type") +
    theme(legend.position = "bottom")
  
  png(file = paste0(out_dir, "cov_range_plot_all_", reg, ".png"),
      height = 3, width = 8,
      units = "in", 
      res = 300)
  
  print(gg_all)
  dev.off()

}
