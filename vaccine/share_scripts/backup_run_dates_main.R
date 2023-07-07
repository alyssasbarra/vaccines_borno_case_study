###############################################################################
## SETUP
###############################################################################

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

# Load `tools` for `md5sum()` function
library(tools)
library(lubridate)


###############################################################################
## OPTIONS
###############################################################################

indicators <- c("dpt1_cov", "dpt2_cov", "dpt3_cov", "dpt0_dose", "dpt1_dose", "dpt2_dose", "dpt1_cond", "dpt2_cond", "dpt1_3_abs_dropout", "dpt1_3_rel_dropout")
run_dates <- c("RUN_DATE","RUN_DATE", "RUN_DATE", "RUN_DATE", "RUN_DATE")
indicator_group <- "vaccine"

### As a function #######################

backup_lv <- data.table(run_dates = run_dates) 

backup_para_output <- parallelize(script = "backup_run_dates_parallel",
	                                log_location = "sgeoutput",
	                                expand_vars = backup_lv,
	                                save_objs = c("indicators", "indicator_group"),
	                                prefix = "backup",
	                                slots = 2,
	                                script_dir = paste0(indic_repo,'share_scripts/'), 
	                                memory = 4,
	                                geo_nodes = FALSE,
	                                use_c2_nodes = FALSE,
	                                singularity = "default")

monitor_jobs(backup_para_output)