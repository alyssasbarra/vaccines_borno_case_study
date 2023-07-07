# HEADER ------------------------------------------------------------------
# Author: USERNAME
# Date: DATE
# Project: Vaccines: Data Prep
# Purpose: Prepare vaccine data
# Details: 
# qlogin -pe multi_slot 32 -l mem_free=100g -now no -P proj_geospatial
# source("FILEPATH")
# source("FILEPATH")
#************************************************************************** 

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
package_list <- c(t(read.csv(sprintf('FILEPATH',commondir),header=FALSE)))

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0('FILEPATH'))
mbg_setup(package_list = package_list, repos = core_repo)

# Custom load indicator-specific functions
source(paste0('FILEPATH'))

## drive locations
root           <- ifelse(Sys.info()[1]=='Windows', 'FILEPATH', 'FILEPATH')

# Custom load indicator-specific functions
source(paste0('FILEPATH'))
source(paste0("FILEPATH"))
source(paste0("FILEPATH"))

# Specific function loads (avoid namespace conflicts)
str_match <- stringr::str_match

### II: DEFINE VARIABLES, ETC. #####################################################

# Which vaccine to run? ------------------------------------------------------------
  vaccine_to_run     <- "dpt" # "dpt" or "mcv" for now

# Variables ------------------------------------------------------------------------
  run_date <- NULL
  #run_date     <- NULL # Specify here as "YYYY_MM_DD" if running on specific date
                       # Otherwise will use today if re-reading data
                       # Or the most recent existing if not
  cores 	 <- 10                   

# Vaccines & vaccine titles --------------------------------------------------------

  # Set up the vaccine depending on the "vaccine_to_run" parameter

  if (vaccine_to_run == "mcv") {

      vaccine_list <- add_vaccine(prefix = "mcv", 
                              title = "MCV",
                              doses = 1,
                              age_min = 12,
                              age_max = 59)

        vaccines <- names(vaccine_list)

        graph_vaccines <- c("mcv1_cov")
        final_model_vaccines <- c("mcv1_cov")
        graph_vaccine_titles <- c("MCV1")

  } else if (vaccine_to_run == "dpt") {

      vaccine_list <- add_vaccine(prefix = "dpt", 
                              title = "DPT",
                              doses = 3,
                              age_min = 12,
                              age_max = 59)

      vaccines <- names(vaccine_list)

      graph_vaccines <- c("dpt3_cov")
      final_model_vaccines <- c("dpt3_cov", "dpt1_cov", "dpt1_3_abs_dropout", "dpt1_3_rel_dropout")
      graph_vaccine_titles <- c("DPT3")
  }

# Toggles --------------------------------------------------------------------------

  skip_IDN        <- F        # Skip Indonesia (big files)
  drop_pre_2000   <- T        # Drop surveys and individuals before 2000

  read_data       <- F        # Read in data from FILEPATH and combine
  merge_data			<- T        # Merge all the files together to vaccine_clean_list_file (below)
  make_graphs			<- T        # Make data coverage graphs
  do_pointpoly 		<- T        # Run point/poly processing
  make_csvs       <- T        # Generate CSVs
  collapse_dhs    <- F        # Collapse a DHS dataset for comparisons
  collapse_method <- "kish"   # `survey` or `kish`
  crop_to_region  <- "all"    # Region to crop data to ("africa", etc); all = no crop

# Directories & files --------------------------------------------------------------
  
  load_data_date  <- NULL # Specify if you want to load data from a specific date
                                  # If NULL will load same as run_date

  #file for Hib introduction dates
  intro_file <- 'FILEPATH'
  
  script_dir <- paste0("FILEPATH")
  in_dir_stem <- paste0('FILEPATH')
  out_dir_stem <- paste0("FILEPATH")

# Prep a run date-------------------------------------------------------------------
  run_date <- Sys.Date() #date-stamp 
  run_date <- gsub("-", "_", run_date)

 if (read_data == F & is.null(load_data_date)) {
  # make list of dates already  generated
    dates <- list.files(in_dir_stem)
    dates <- sort(dates[grep(".*_.*_.*", dates)], decreasing = T)
  # pull most recent date folder & use that
    in_run_date <- dates[1]
 } else if (read_data == T & is.null(load_data_date)) {
  # just use the most recent run date
    in_run_date <- run_date
 } else if (read_data == F & !is.null(load_data_date)) {
  # use the specified date
    in_run_date <- load_data_date
 }

 # Set up input dir
 in_dir <- paste0(in_dir_stem, in_run_date, "/")
 dir.create(in_dir, showWarnings = F)

 out_dir <- paste0(out_dir_stem, run_date, "/")
 dir.create(out_dir, showWarnings = F)  

# File for temporary, cleaned data set of processed vaccines - in same data as input data
# Will have the following added to it:  "..._[vaccine_prefix].rds"
 vaccine_cleaned_file_prefix <- paste0(in_dir, vaccine_to_run, "_combined_vaccine_clean_")

# File for temporary, *resampled* data set of processed vaccines - in same data as input data
# Will have the following added to it:  "..._[vaccine_prefix].rds"
 vaccine_resampled_file_prefix <- paste0(in_dir, vaccine_to_run, "_combined_vaccine_resampled_")

 # Log directory - now in output directory
 log_dir <- paste0(out_dir, "logs_", vaccine_to_run, "/")
 dir.create(log_dir, showWarnings = F, recursive = T)

### III: RUN BITS PER TOGGLES ######################################################

if (read_data == T) {
  # This will read in all data and geoposition it
  # Output: combined CSV files, one per survey series, in "FILEPATH"
  print_header()
  message("\nReading data...\n")
  source(paste0(script_dir, "01_highspeed_cluster_post_extraction_v2.R"))
}

if (merge_data == T) {
 # This will combine the survey-specific CSV files, then merge and apply case definitions
 # Output: an RData file (one line per individual) in "FILEPATH"
 print_header()
 message("\nMerging vaccine data...\n")
 source(paste0(script_dir, "02_merge_vaccine.R"))
}

if (make_graphs == T) {
 # Creates graphs from the file made in merge_data
 # Outputs: vaccine-specific graphs in "FILEPATH"
 print_header()
 message("\nGraphing coverage...\n")
 source(paste0(script_dir, "03_graph_coverage.R"))
}

if (do_pointpoly == T) {
 # Creates graphs from the file made in merge_data
 # Outputs: pointpoly-ed individual csvs in out_dir
 print_header()
 message("\nRunning pointpoly...\n")
 source(paste0(script_dir, "04_pointpoly.R"))
}

if (make_csvs == T) {
 # Creates graphs from the file made in merge_data
 # Outputs: pointpoly-ed individual csvs in out_dir
 print_header()
 message("\nMaking CSVs...\n")
 source(paste0(script_dir, "05_generate_vaccine_csvs.R"))
}
