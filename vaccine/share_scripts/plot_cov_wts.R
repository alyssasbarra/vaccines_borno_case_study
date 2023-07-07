## clear environment
rm(list=ls())

## Set repo location and indicator group
user               <- Sys.info()['user']
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

## Options ---------------------------------------

run_date <- "RUN_DATE"
indicator <- "dpt3_cov"
indicator_group <- "vaccine"
regions <- c("cssa", "essa", "name", "sssa", "wssa")c

output_plot_dir <- paste0('FILEPATH')
dir.create(output_plot_dir)

## Make covariate importance heat map
if(length(list.files(paste0('FILEPATH'), pattern = 'child_model_list')) != 0) {

  for (rr in regions) {

    message(reg)
    
    cov_importance <- get.cov.wts(rd = run_date,
                                  ind = indicator, ## indicator
                                  ind_gp = indicator_group, ## indicator_group
                                  reg = rr,
                                  age = 0,
                                  holdout = 0)

    cov_gg <- plot.cov.wts(rd = run_date,
                           ind = indicator, ## indicator
                           ind_gp = indicator_group, ## indicator_group
                           reg = rr,
                           age = 0,
                           holdout = 0,
                           plot.inla.col = F)

    png(file = paste0(output_plot_dir, indicator, "_rel_cov_imp_", reg, ".png"),
        height = 4,
        width = 8,
        units = "in", 
        res = 300)

    print(cov_gg)
    dev.off()
  }
  
}
