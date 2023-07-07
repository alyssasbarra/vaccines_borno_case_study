# clear memory
rm(list=ls())

### I: SETUP #######################################################################

repo <- 'FILEPATH'

## Set repo location and indicator group
user               <- 'USERNAME'
core_repo          <- sprintf('FILEPATH',user)
indic_repo         <- sprintf('FILEPATH',user)

root <- ifelse(Sys.info()[1]=="Windows", "FILEPATH", "FILEPATH")
j_root <- "FILEPATH"
## Load libraries and miscellaneous MBG project functions.
setwd(core_repo)
root <- ifelse(Sys.info()[1]=="Windows", "FILEPATH", "FILEPATH")
  package_lib    <- ifelse(grepl('geos', Sys.info()[4]),
                          paste0('FILEPATH'),
                          paste0('FILEPATH'))
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
source(FILEPATH("FILEPATH"))
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

# Plot error
indicator_group <- "vaccine"
run_date <- "RUN_DATE"
regions <- c("wssa", "cssa", "essa", "sssa", "name")
ig <- indicator_group

# Background shapes
background_shape <- readRDS("FILEPATH")

background_shape@data$id <- rownames(background_shape@data)
background_shape_points <- fortify(background_shape, region = "id")
background_shape_df <- join(background_shape_points, background_shape@data, by = "id")
background_shape_df <- as.data.table(background_shape_df)

theme_empty <- theme_classic() +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        strip.background = element_blank())

color_scheme <- c("#990021", "#CC51AB", "#CE89E5", "#6151CC", "#004999")

for (indicator in c("dpt3_cov", "dpt1_cov", "dpt1_3_abs_dropout")) {

  if (indicator == "dpt3_cov") indicator_title <- "DPT3 coverage"
  if (indicator == "dpt1_cov") indicator_title <- "DPT1 coverage"
  if (indicator == "dpt1_3_abs_dropout") indicator_title <- "DPT 1-3 absolute dropout"
  
  out_dir <- paste0("FILEPATH")
  dir.create(out_dir, showWarnings = F)

  message(paste0("Working on ", indicator, "..."))
  
  sharedir <- paste0("FILEPATH")
  
  message("   loading data...")
  input_df <- fread(paste0(sharedir, "output_draws_data.csv"))
  
  draw_cols <- names(input_df)[grepl("draw", names(input_df))]
  
  input_df[, pred_mean := rowMeans(.SD, na.rm = T), .SDcols = draw_cols]
  input_df[, me := pred_mean - (get(indicator) / N)]
  
  # Drop draws for memory
  input_df <- subset(input_df, select = !(names(input_df) %in% draw_cols))
  
  # Set years
  input_df[year > 1999 & year <= 2003, plot_year := "2000-2003"]
  input_df[year > 2003 & year <= 2007, plot_year := "2004-2007"]
  input_df[year > 2007 & year <= 2011, plot_year := "2008-2011"]
  input_df[year > 2011 & year <= 2016, plot_year := "2012-2016"]
  
  # Basic plot
  gg <- ggplot() +
    geom_polygon(data = background_shape_df,
                 aes(x = long, y = lat, group = group),
                 fill="white") +
    geom_path(data = background_shape_df,
              aes(x = long, y = lat, group = group),
              color="black", size = 0.2) +
    geom_point(data = input_df[weight != 1,], 
               aes(x = longitude, y = latitude, size = N, alpha = weight, color = me),
               pch = 16) +
    geom_point(data = input_df[weight == 1,], 
               aes(x = longitude, y = latitude, size = N, alpha = weight, color = me),
               pch = 16) +
    facet_wrap(~plot_year, ncol = 2) +
    scale_size_area(max_size = 1) +
    #scale_color_distiller(palette = "RdBu", limits = c(-1,1)) +
    scale_color_gradientn(colors = color_scheme, limits = c(-1,1)) +
    coord_equal() +
    labs(title = paste0("Residual error: ", indicator_title), 
         color = "Residual error", alpha = "Weight", size = "N") +
    theme_empty
  
  message("    plotting...")
  png(file = paste0(out_dir, "error_plot_", indicator, ".png"),
      height = 7,
      width = 8,
      units = "in",
      res = 300)
  print(gg)
  dev.off()
  
}