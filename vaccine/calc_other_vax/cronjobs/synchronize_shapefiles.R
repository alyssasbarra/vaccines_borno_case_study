# HEADER ------------------------------------------------------------------
# Author: USERNAME
# Date: DATE
# Project: Geospatial - General
# Purpose: Synchronize shapefile directories & report results
# Details: To be run as a cronjob
#************************************************************************** 

## clear environment
rm(list=ls())

## Set repo location and indicator group
user               <- 'USERNAME'
core_repo          <- sprintf('FILEPATH')
indic_repo         <- sprintf('FILEPATH')

setwd(core_repo)
commondir      <- sprintf('FILEPATH')
package_list <- c(t(read.csv(sprintf('FILEPATH',commondir),header=FALSE)))

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0('FILEPATH'))
mbg_setup(package_list = package_list, repos = core_repo)

# Custom load indicator-specific functions
source(paste0('FILEPATH'))

# Set up directories
log_dir        <- "FILEPATH"
time_stamp     <- gsub("[-|:| ]", "_", Sys.time())
out_file       <- paste0(log_dir, "cron_log_", time_stamp, ".csv")

# Run
output <- synchronize_shapefile_directories(verbose = T, cores = 20)

# Write response
fwrite(output, file = out_file)