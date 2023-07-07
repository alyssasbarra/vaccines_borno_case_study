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

## sort some directory stuff and pull newest code into share
setwd(core_repo)
if(pullgit) system(sprintf('cd %s\ngit pull %s %s', indic_repo, remote, branch))

## drive locations
root           <- ifelse(Sys.info()[1]=='Windows', 'FILEPATH', 'FILEPATH')
package_lib    <- ifelse(grepl('geos', Sys.info()[4]),
                          sprintf('FILEPATH'),
                          sprintf('FILEPATH'))
commondir      <- sprintf('FILEPATH')

## Load libraries and  MBG project functions.
.libPaths(package_lib)
package_list <- c(t(read.csv(sprintf('FILEPATH',commondir),header=FALSE)))

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
source(paste0('FILEPATH'))

###############################################################################
indicator_group <- 'vaccine'
indicator <- "dpt3_cov"
run_date <- 'RUN_DATE'
input_date <- 'DATE'  

indicators <- c("dpt3_cov", "dpt1_cov", "dpt1_3_abs_dropout")
regions <- c("wssa", "cssa", "name", "sssa", "essa")
raked_vals <- c("_raked", "")

dhs_lv <- expand.grid(indicators,
                      regions,
                      raked_vals,
                      stringsAsFactors = F)
dhs_lv <- as.data.table(dhs_lv)
names(dhs_lv) <- c("indicator", "region", "use_raked")
 
dhs_qsub_output <- make_qsub2(code = "dhs_parallel_script",
                              script_dir = "share_scripts",
                              script_repo = indic_repo, 
                              log_location = "sgeoutput",
                              lv_table = dhs_lv,
                              save_objs = c('indicator_group', 'run_date', 'input_date'),
                              prefix = "dhs",
                              cores = 4,
                              memory = 20,
                              geo_nodes = TRUE)

waitformodelstofinish2(dhs_qsub_output)

for (ind in indicators) {
  # Set up the directory with DHS compare data
  in_dir <- paste0("FILEPATH")
  
  # Do this once each for raked & unraked 
  for (rake in raked_vals) {

    message(ind, " | ", rake)
  
    # Grab all data for raked/unraked
    all_df <- lapply(regions, function(reg) {
      reg_df <- fread(paste0(in_dir, reg, "_", rake, ".csv"))
      reg_df[, region := reg]
      return(reg_df)                      
    })
    all_df <- rbindlist(all_df)

  if (ind == "dpt3_cov") title_ind <- "DPT3 Coverage"
  if (ind == "dpt1_cov") title_ind <- "DPT1 Coverage"
  if (ind == "dpt1_3_abs_dropout") title_ind <- "DPT 1-3 Absolute Dropout"

  plot_title <- paste0(title_ind, " at the Admin 1 Level")

  # carto_colors <- c("#5F4690","#1D6996","#38A6A5","#0F8554",
  #                   "#73AF48","#EDAD08","#E17C05","#CC503E",
  #                   "#94346E","#6F4070","#994E95","#666666")
  
  all_df[region == "wssa", region := "Western Sub-Saharan Africa"]
  all_df[region == "cssa", region := "Central Sub-Saharan Africa"]
  all_df[region == "sssa", region := "Southern Sub-Saharan Africa"]
  all_df[region == "essa", region := "Eastern Sub-Saharan Africa"]
  all_df[region == "name", region := "Northern Africa"]

  gg <- ggplot(all_df, aes(x=outcome,y=geo_mean))+
    geom_abline(intercept=0,slope=1,colour='red')+
    geom_point(aes(size = N), colour='black', alpha = 0.2, pch = 16)+
    theme_bw()+
    xlab('Data Estimate (DHS)') +
    ylab('Mean Prediction (MBG)')+
    theme(strip.background = element_rect(fill="white"))+
    geom_errorbar(aes(ymin=geo_lower, ymax=geo_upper), colour="black", width=0, size=.5, alpha = 0.2) +
    geom_abline(intercept=0,slope=1,colour='red') +
    scale_size_area() +
    coord_equal() +
    scale_x_continuous(expand = c(0, 0), limits = c(0,1)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
    labs(title = plot_title,
         subtitle = paste0("DHS vs ", 
                           ifelse(rake == "", 
                                  "uncalibrated", 
                                  "calibrated"), 
                           " MBG estimates"),
         size = "N",
         color = "Modeling region")

  png(file = paste0(in_dir, "compare_dhs_mbg_", ind, rake, ".png"),
      width = 8, 
      height = 7,
      units = "in",
      res = 300)
  print(gg)
  dev.off()

  # By region
  gg <- gg + facet_wrap(~region)
  png(file = paste0(in_dir, "compare_dhs_mbg_", ind, rake, "_by_region.png"),
      width = 8, 
      height = 7,
      units = "in",
      res = 300)
  print(gg)
  dev.off()
  }
}
