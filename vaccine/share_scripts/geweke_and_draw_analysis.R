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

## sort some directory stuff
commondir      <- sprintf('FILEPATH')
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0('FILEPATH'))
mbg_setup(package_list = package_list, repos = core_repo)

# Custom load indicator-specific functions
source(paste0('FILEPATH'))

# Options
age <- 0
holdout <- 0
run_date <- "RUN_DATE"
indicators <- c("dpt3_cov", "dpt1_cov", "dpt1_3_abs_dropout")
regions <- c("cssa", "essa", "wssa", "name", "sssa")

out_dir <- "FILEPATH"
dir.create(out_dir)

# Launch jobs in parallel
geweke_lv <- data.table(expand.grid(ind = indicators, reg = regions, stringsAsFactors = F))

geweke_para_output <- parallelize(script = "geweke_and_draw_analysis_parallel",
                                  log_location = "sgeoutput",
                                  lv_table = geweke_lv,
                                  save_objs = c("run_date", "out_dir"),
                                  prefix = "geweke",
                                  rd = NULL,
                                  slots = 12,
                                  memory = 64,
                                  script_dir = paste0(indic_repo, "share_scripts/"),
                                  geo_nodes = F,
                                  use_c2_nodes = T,
                                  singularity = 'default')

monitor_jobs(geweke_para_output, notification = "pushover", max_tries = 1)

# Print combined geweke statistics after loading from disk

geweke_list <- lapply(1:nrow(geweke_lv), function(i) {
  reg <- geweke_lv[i, reg]
  ind <- geweke_lv[i, ind]
  df_geweke <- readRDS(paste0(out_dir, ind, "_", reg, "_", "geweke_vals.rds"))
})

geweke_plot <- rbindlist(geweke_list)
geweke_plot[indicator == "dpt3_cov", ind_title := "DPT3 Coverage"]
geweke_plot[indicator == "dpt1_cov", ind_title := "DPT1 Coverage"]
geweke_plot[indicator == "dpt1_3_abs_dropout", ind_title := "DPT1-3 Absolute Dropout"]

geweke_plot[region == "cssa", reg_title := "Central Sub-Saharan Africa"]
geweke_plot[region == "essa", reg_title := "Eastern Sub-Saharan Africa"]
geweke_plot[region == "name", reg_title := "North Africa"]
geweke_plot[region == "sssa", reg_title := "Southern Sub-Saharan Africa"]
geweke_plot[region == "wssa", reg_title := "Western Sub-Saharan Africa"]

color_scale <- c("Standard Normal Distribution" = "black",
                 "Central Sub-Saharan Africa" = "#7F3C8D",
                 "Eastern Sub-Saharan Africa" = "#11A579",
                 "North Africa" = "#3969AC",
                 "Southern Sub-Saharan Africa" = "#F2B701",
                 "Western Sub-Saharan Africa" = "#E73F74")

gg <- ggplot(data = geweke_plot, aes(x = geweke)) +
    geom_line(aes(color = reg_title), stat = "density", alpha = 0.5) +
    stat_function(aes(color = "Standard Normal Distribution"), 
                  fun = dnorm, n = 10000, args = list(mean=0, sd=1), 
                  alpha = 0.5) +
    theme_classic() +
    labs(x = "Geweke statistic",
         y = "Density",
         title = "Distributions of Geweke statistics by indicator and region",
         subtitle = "Distribution of Geweke statistics (for 10,000 randomly sampled pixels)") +
    scale_color_manual(name = NULL, 
                       values = color_scale, 
                       breaks = names(color_scale)) +
    facet_wrap(~ind_title)
  
  png(file=paste0(out_dir, "geweke_all.png"), 
      height = 4, 
      width = 8,
      res = 200, 
      units = "in")
  
  print(gg)
  dev.off()
