###############################################################################
## SETUP
###############################################################################

## clear environment
rm(list=ls())

## Set repo location and indicator group
user               <- 'USERNAME'
core_repo          <- sprintf('FILEPATH',user)
indic_repo         <- sprintf('FILEPATH',user)
remote             <- 'origin'
branch             <- 'develop'
pullgit            <- FALSE

## sort some directory stuff and pull newest code into share
setwd(core_repo)
if(pullgit) system(sprintf('cd %s\ngit pull %s %s', indic_repo, remote, branch))

## drive locations
root           <- ifelse(Sys.info()[1]=='Windows', 'FILEPATH', 'FILEPATH')
package_lib    <- ifelse(grepl('geos', Sys.info()[4]),
                          sprintf('FILEPATH',root),
                          sprintf('FILEPATH',root))
commondir      <- sprintf('FILEPATH')

## Load libraries and  MBG project functions.
.libPaths(package_lib)
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))

if(Sys.info()[1] == 'Windows'){
  stop('STOP! you will overwrite these packages if you run from windows\n
        STOP! also, lots of this functions wont work so get on the cluster!')
} else {
  for(package in package_list)
    require(package, lib.loc = package_lib, character.only=TRUE)
  for (funk in list.files(core_repo, recursive=TRUE, pattern='functions')) {
    message(paste("loading from core repo:", funk))
    source(paste0(core_repo, funk))
  }
}

# Custom load indicator-specific functions
source(paste0(indic_repo,'functions/misc_vaccine_functions.R'))

## ##################################################################################################
## ##################################################################################################

## setup some required objects and directories
# Varying: indicator, run_date
# Static:indicator_group, samples

load_from_qsub() 
ind <- indicator
rd <- run_date

mod.dir <- sprintf('FILEPATH')
si.fig.dir <- paste0('FILEPATH')
dir.create(si.fig.dir, showWarnings = F)

## for admin0
draws.df <- fread(sprintf("FILEPATH"))

samples <- length(grep('draw', colnames(draws.df)))

## ad0
message("Admin 0")
ad0.pvtable <- get_pv_table(d = draws.df,
                            indicator_group = indicator_group,
                            rd = run_date,
                            plot=F,
                            indicator=indicator,
                            aggregate_on='country',
                            result_agg_over="oos",
                            draws = as.numeric(samples),
                            plot_ci = TRUE,
                            out.dir = si.fig.dir,
                            save_csv = F)
cols <- names(ad0.pvtable)[-(1:2)]
ad0.pvtable[,(cols) := round(.SD,3), .SDcols=cols]
ad0.pvtable[, c(4, 5)] <- NULL
ad0.pvtable <- ad0.pvtable[, c(1, 2, 5, 3, 4, 6, 7), with = FALSE]
write.csv(ad0.pvtable,
          file = sprintf("%s/ad0_metrics.csv",
                         si.fig.dir),
          row.names = FALSE)

## ad1
message("Admin 1")
ad1.pvtable <- get_pv_table(d = draws.df,
                            indicator_group = indicator_group,
                            rd = run_date,
                            plot=F,
                            indicator=indicator,
                            aggregate_on='ad1',
                            result_agg_over="oos",
                            draws = as.numeric(samples),
                            plot_ci = TRUE,
                            out.dir = si.fig.dir,
                            save_csv = F)
cols <- names(ad1.pvtable)[-(1:2)]
ad1.pvtable[,(cols) := round(.SD,3), .SDcols=cols]
ad1.pvtable[, c(4, 5)] <- NULL
ad1.pvtable <- ad1.pvtable[, c(1, 2, 5, 3, 4, 6, 7), with = FALSE]
write.csv(ad1.pvtable,
          file = sprintf("%s/ad1_metrics.csv",
                         si.fig.dir),
          row.names = FALSE)

## ad2
message("Admin 2")
ad2.pvtable <- get_pv_table(d = draws.df,
                            indicator_group = indicator_group,
                            rd = run_date,
                            plot=F,
                            indicator=indicator,
                            aggregate_on='ad2',
                            result_agg_over="oos",
                            draws = as.numeric(samples),
                            plot_ci = TRUE,
                            out.dir = si.fig.dir,
                            save_csv = F)

cols <- names(ad2.pvtable)[-(1:2)]
ad2.pvtable[,(cols) := round(.SD,3), .SDcols=cols]
ad2.pvtable[, c(4, 5)] <- NULL
ad2.pvtable <- ad2.pvtable[, c(1, 2, 5, 3, 4, 6, 7), with = FALSE]
write.csv(ad2.pvtable,
          file = sprintf("%s/ad2_metrics.csv",
                         si.fig.dir),
          row.names = FALSE)