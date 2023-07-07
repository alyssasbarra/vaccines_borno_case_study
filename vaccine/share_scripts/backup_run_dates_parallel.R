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

# Load options from parallelize()
#  - indicators
#  - run_dates
#  - indicator_group

load_from_parallelize()

message(paste0("Run Dates: ", paste(run_dates, collapse = ", ")))
message(paste0("Indicators: ", paste(indicators, collapse = ", ")))
message(paste0("Indicator Group: ", indicator_group))        

###############################################################################
## OPTIONS
###############################################################################

### As a function #######################

## Create a table of directories
all_dirs <- expand.grid(indicators, run_dates, stringsAsFactors = F) %>% 
							as.data.table %>%
							setnames(c("Var1", "Var2"), c("indicator", "run_date"))

archive_run(indicator = all_dirs$indicator, 
            indicator_group = indicator_group, 
            run_date = all_dirs$run_date, 
            remove_source_files = F, 
            direction = "share_to_j")
