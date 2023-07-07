# USERNAME
# DATE <- DATE <- DATE
# Make Admin 1 and 2 level summary draw-level objects
# Based on USERNAME code: FILEPATH
# Writing own because not wanting to wait for PR

# indicator <- "dpt3_cov"
# indicator_group <- "vaccine"
# run_date <- "RUN_DATE"
# raked <- T
# pop_measure <- "a0004t"
# overwrite <- T
# age <- 0
# holdout <- 0
# region <- "vax_essa"
# rasterize_field <- "GAUL_CODE"
# subset_field <- "ISO"
# subset_field_codes <- c("KEN")
# filename_addin <- "_KEN"
# shapefile_path <- "FILEPATH"
#     Can also just pass a shapefile (without .shp) at the end
# sp_hierarchy_columns <- c("NAME_0", "NAME_1", "GAUL_CODE")
#     Has some defaults in there if needed
# crop_shapefile <- TRUE

# Setup
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## clear environment
rm(list=ls())

# Use parallelize infrastructure
source("FILEPATH")

# What to load through parallelize: 
#  - varying: c("shapefile", "filename_addin", "region", "subset_field", "subset_field_codes")
#  - static: c("indicator", "indicator_group", "run_date", "raked", "pop_measure",
#              "overwrite", "age", "holdout", "crop_shapefile", "use_lookup", core_repo")

load_from_parallelize()  

## sort some directory stuff
commondir      <- sprintf('FILEPATH')
package_list <- c(t(read.csv(sprintf('FILEPATH',commondir),header=FALSE)))

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0('FILEPATH'))
mbg_setup(package_list = package_list, repos = core_repo)

# Custom function
# library(HatchedPolygons, lib.loc = "FILEPATH")

## -------------------------------------------------------------------------------------------------
## Start plotting script here -------------------------------------------------------------------
## -------------------------------------------------------------------------------------------------

# Define modeling directory
sharedir <- paste0("FILEPATH")

# Define output directory
output_dir <- paste0("FILEPATH")
dir.create(output_dir, showWarnings = F)

# Load the JRF data frame
df_jrf_all <- fread("FILEPATH")
setnames(df_jrf_all, 
         c("NAME.who_adm1", "NAME.who", "NAME.shpfile"),
         c("who_name_adm1", "who_name", "shpfile_name"))

if (indicator == "dpt3_cov") {
  df_jrf_all <- subset(df_jrf_all, vaccine == "DTP3")
} else {
  stop("Need to define indicator / vaccine name relationship.")
}

# Determine if ad1 or ad2
if (length(unique(df_jrf_all) > 1)) {
  admin_level <- 2
} else {
  admin_level <- 1
}

# Define some objects
if (!is.null(subset_field)) {
  if (is.na(subset_field) | subset_field == "NA") subset_field <- NULL
  if (is.na(subset_field_codes) | subset_field_codes == "NA") subset_field_codes <- NULL
}

# Pull the country name
country_name <- get_location_code_mapping(shapefile_version = "current") %>%
                .[ihme_lc_id == stringr::str_match(filename_addin, "_(.*)")[,2], loc_name]

# Load the aggregates
df_agg_all <- fread(paste0(sharedir, "custom_aggregation/", shapefile, "/", indicator, "_admin_raked_summary", filename_addin, ".csv"))

# Subset JRF to the country of interest
df_jrf_all <- subset(df_jrf_all, ADM0_NAME == country_name)

# Replace -2222 with NA
df_jrf_all[who_coverage == -2222, who_coverage := NA]

# Find JRF duplicates and set to NA for now
dups <- df_jrf_all[duplicated(df_jrf_all$location_code), location_code] %>% .[!is.na(.)] %>% unique
df_jrf_all[location_code %in% dups, shpfile_name := NA]
df_jrf_all[location_code %in% dups, location_code := NA]
df_jrf_all[, location_code := as.numeric(location_code)]

# Use lookup table to get shapefile attributes
df_shapefile_lookup <- fread("FILEPATH")
if (!(shapefile %in% df_shapefile_lookup$shp)) stop("Need to add a row to lookup table for this shapefile!")
assign("shapefile_path", paste0(df_shapefile_lookup[shp == shapefile]$shapefile_directory, 
                              df_shapefile_lookup[shp == shapefile]$shp, ".shp"))
assign("sp_hierarchy_columns", eval(parse(text = df_shapefile_lookup[shp == shapefile]$sp_hierarchy_columns)))

# Merge
df_all <- merge(subset(df_agg_all, year == 2017, select = c(sp_hierarchy_columns, "mean", "upper", "lower", "pop")),
                subset(df_jrf_all, select = c("iso3", "vaccine", "who_name_adm1", "who_name", "who_coverage", "shpfile_name", "merge_status", "location_code")),
                by.x = field_name, by.y = "location_code", all.x = T, all.y = T)

n_missing_jrf <- nrow(df_all[is.na(who_coverage)])
n_missing_agg <- nrow(df_all[is.na(mean)])
       
df_all[, who_coverage := who_coverage / 100]

df_plot <- subset(df_all, !is.na(mean) & !is.na(who_coverage))

setorderv(df_plot, "mean")
df_plot[, mbg_rank := .I]

# Reshape to long
df_mbg <- subset(df_plot, select = c(sp_hierarchy_columns, "mean", "upper", "lower", "mbg_rank"))
df_jrf <- subset(df_plot, select = c(sp_hierarchy_columns, "mbg_rank", "who_coverage"))
setnames(df_jrf, "who_coverage", "mean")
df_jrf[, source := "JRF"]
df_mbg[, source := "MBG"]
df_plot_long <- rbind(df_jrf, df_mbg, fill = T)

plot_ordered_admins <- function(plot_df, country_title, error_bars = "yaxis", admin_level, truncate_above = 200, n_missing_jrf) {
  # error_bars: none, yaxis, xaxis, or both
  plot_df <- copy(plot_df)
  
  if (admin_level == 1) {
    xaxis_title <- "First-level administrative unit (ordered by mean IHME estimate)"
  }
  if (admin_level == 2) {
    xaxis_title <- "Second-level administrative unit (ordered by mean IHME estimate)"
  }
  
  plot_df[mean <= truncate_above/100, actual_value := T]
  plot_df[mean > truncate_above/100, actual_value := F]
  plot_df[mean > truncate_above/100, mean := truncate_above/100]
  
  # Rename for clarity
  plot_df[source == "JRF", source := "JRF (country-reported)"]
  plot_df[source == "MBG", source := "IHME (MBG estimate)"]
  
  gg_p <- ggplot(data = plot_df[actual_value == T],
                 aes(x = mbg_rank, 
                     y = mean,
                     shape = source,
                     color = source)) 
  
  if (error_bars == "yaxis" | error_bars == "both") {
    gg_p <- gg_p + geom_errorbar(data = plot_df[actual_value == T & source == "IHME (MBG estimate)"], 
                                 aes(x = mbg_rank, ymin = lower, ymax = upper), width = 0, alpha = 0.2)
  }
  
  gg_p <- gg_p + 
    geom_point(alpha = 0.7, size = 2) +
    geom_point(data = plot_df[actual_value == F], aes(x = mbg_rank, y = mean), 
               size = 6, color = "black", alpha = 0.7, shape = 94) +
    geom_abline(slope = 0, intercept = 1, color = "black") +
    geom_abline(slope = 0, intercept = 0, color = "black") +
    #geom_abline(slope = 0, intercept = 0.8, color = "black", linetype = "dotted") +
    #geom_smooth(method = "loess", color = "black") +
    theme_classic() +
    labs(y = paste0("Mean DPT3 Coverage: 2017"),
         x = xaxis_title, 
         title = country_title,
         shape = "Source",
         color = "Source",
         caption = paste0(ifelse(!all(plot_df$actual_value == T), paste0("Reported administrative coverage values > ", 
                                                                         truncate_above, "% represented by the ^ symbol at top of graph"), ""),
                          ifelse(n_missing_jrf > 0, paste0("\nExcludes ", n_missing_jrf, 
                                                           " administrative unit", ifelse(n_missing_jrf == 1, "", "s"), 
                                                           " reported through JRF for which no geographical boundaries were available"), ""))) +
    scale_y_continuous(expand = c(0, 0), limits = c(0,(truncate_above/100)+0.1), labels = scales::percent) +
    scale_color_manual(values = c("#630758", "black")) +
    scale_shape_manual(values = c(19, 4)) +
    theme(legend.position="right") +
    theme(axis.text.x=element_blank(),
          axis.ticks.x = element_blank(),
          plot.caption = element_text(hjust = 0.5)) +
    guides(color = guide_legend(override.aes = list(linetype = 0)))
  
  return(gg_p)
  
}

png(file = paste0(output_dir, "ordered_admin_comparison_2017", filename_addin, ".png"),
    height = 6,
    width = 12,
    units = "in",
    res = 300)
plot_ordered_admins(plot_df = df_plot_long, country_title = country_name, admin_level = admin_level, n_missing_jrf = n_missing_jrf)
dev.off()

max_val <- max(max(df_plot$mean), max(df_plot$who_coverage), 1)

### Make a scatter-plot
gg_scatter <- ggplot(data = df_plot, aes(x = who_coverage, y = mean)) +
                geom_abline(slope = 1, intercept = 0) +
                geom_vline(xintercept = 1, color = "red", linetype = "dashed") +
                geom_point(aes(size = pop), alpha = 0.5) +
                geom_errorbar(aes(ymin = lower, ymax = upper), alpha = 0.2) +
                theme_classic() +
                scale_size_area() +
                coord_equal() +
                labs(x = "JRF (country-reported)",
                     y = "IHME (MBG estimate)",
                     title = country_name,
                     size = "Number of children\nunder 5 years\n(2017)",
                     caption = ifelse(n_missing_jrf > 0, paste0("\nExcludes ", n_missing_jrf, 
                                                           " administrative units reported through JRF for which no geographical boundaries were available"), "")) +
                scale_x_continuous(labels = scales::percent, limits = c(0, max_val)) +
                scale_y_continuous(labels = scales::percent, limits = c(0, max_val)) +
                theme(legend.position = "right",
                      plot.caption = element_text(hjust = 0.5))

png(file = paste0(output_dir, "scattered_admin_comparison_2017", filename_addin, ".png"),
    height = 6,
    width = 12,
    units = "in",
    res = 300)
print(gg_scatter)
dev.off()

### Make some maps
ad2_spdf <- rgdal::readOGR(shapefile_path, stringsAsFactors = FALSE, GDAL1_integer64_policy=TRUE)

if (!is.null(subset_field) & !is.null(subset_field_codes)) {
    message("Subsetting to `subset_field_codes`")
    ad2_spdf <- subset(ad2_spdf, get(subset_field) %in% subset_field_codes)
}

plot_est <- function(plot_df, src, the_spdf, loc_field, title, subtitle, legend_title, caption) {
  
  plot_df <- copy(plot_df)
  plot_df <- subset(plot_df, source == src)
  
  # Format spdf
  the_spdf <- merge(the_spdf, subset(plot_df, select = c(loc_field, "mean")), by = loc_field)
  the_spdf@data$id <- rownames(the_spdf@data)
  the_spdf.points <- fortify(the_spdf, region = "id")
  the_spdf.df <- join(the_spdf.points, the_spdf@data, by = "id") %>% as.data.table
  
  # if (src == "JRF") {
  #   hatch_me <- (subset(the_spdf, mean > 1 & !is.na(mean)))    
  #   hatched <- lapply(1:nrow(hatch_me), function(i) {try(suppressWarnings(suppressMessages(hatched.SpatialPolygons(the_spdf[i,], angle = 135, density = 8))))})
    
  #   # Remove small geometries
  #   hatched_classes <- sapply(hatched, "class")
  #   hatched <- hatched[which(hatched_classes == "SpatialLinesDataFrame")]
    
  #   hatched_sldf <- do.call(rbind, hatched)
    
  #   hatched_df <- fortify(hatched_sldf) %>% as.data.table
    
  # }

  # Custom vaccines color scale -----------------------------------------------------
  vals <- c( 1,       0.8,        0.6,      0.496094,   0.4,       0.2,       0.164063,  0.000000)
  cols <- c("#3D649D", "#91BEDC",  "#DEF3F8", "#FFE291", "#FEE191", "#FC8D58", "#FC8C58", "#B03027")

  geom_polygon_quiet <- function(...) {suppressMessages(ggplot2::geom_polygon(...))}
  geom_path_quiet     <- function(...) {suppressMessages(ggplot2::geom_path(...))}

  gg <- ggplot() +
    geom_polygon_quiet(data = subset(the_spdf.df, mean > 1),
                       aes(x=long, y=lat, group=group),
                       fill = "#8436A8", color = "black", size = 0.3)

  # if (src == "JRF") {
  #   gg <- gg + geom_line(data = hatched_df, 
  #                        aes(x=long, y=lat, group=group), 
  #                        color = "black", size = 0.2) 
  # }

    gg <- gg + geom_polygon_quiet(data = subset(the_spdf.df, mean <= 1 | is.na(mean)),
                                  aes(x=long, y=lat, group=group, fill=mean),
                                  size = 0.3, color = "black") +
    scale_fill_gradientn(colors = cols, values = vals, 
                         breaks = c(0,0.2,0.4,0.6,0.8,1),
                         limits = c(0,1),
                         labels = scales::percent) +
    theme_empty() +
    coord_equal() +
    labs(title = title,
         subtitle = subtitle,
         fill = legend_title,
         caption = caption) +
    theme(plot.caption = element_text(hjust = 0.5),
          plot.title = element_text(hjust = 0.5))

  return(gg)
}

gg_jrf <- plot_est(plot_df = df_plot_long, 
                   src = "JRF", 
                   the_spdf = ad2_spdf, 
                   loc_field = field_name, 
                   title = "JRF-Reported Subnational Administrative Data", 
                   subtitle = NULL, 
                   legend_title = "DPT3 Coverage", 
                   caption = paste0("Administrative units with no matches found in JRF data are displayed in gray.\n",
                                    "Administrative units with JRF reported coverage > 100% are displayed in purple"))


gg_mbg <- plot_est(plot_df = df_plot_long, 
                   src = "MBG", 
                   the_spdf = ad2_spdf, 
                   loc_field = field_name, 
                   title = "IHME Estimates (MBG Models)", 
                   subtitle = NULL, 
                   legend_title = "DPT3 Coverage", 
                   caption = paste0("Administrative units with no matches found in JRF data are displayed in gray.\n",
                                    " "))

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

the_legend <- g_legend(gg_jrf)

title_grob <- textGrob(paste0(country_name, ": Subnational DPT3 Coverage, 2017"), gp = gpar(fontsize = 24, fontface = "bold"))

gg_jrf <- gg_jrf + theme(legend.position = "none")
gg_mbg <- gg_mbg + theme(legend.position = "none")

lay <- rbind(c(1,1,1,1,1,1,1),
             c(2,2,2,3,3,3,4),
             c(2,2,2,3,3,3,4),
             c(2,2,2,3,3,3,4))

plot_all <- arrangeGrob(title_grob, gg_jrf, gg_mbg, the_legend, 
                        layout_matrix = lay,
                        heights = c(0.2,1,1,1))

png(file = paste0(output_dir, "map_admin_comparison_2017", filename_addin, ".png"),
    width = 14,
    height = 8,
    units = "in", 
    res = 300)
grid.draw(plot_all)
dev.off()

