####################################################################################################
## Description:   Make plots of national-level temporal trends and data.
####################################################################################################

require(data.table)
require(ggplot2)
require(grid)
require(gridExtra)
require(RColorBrewer)

admin0_data_and_estimates_plots <- function(run_date, indicator) {
  
  ## Load data and estimates -------------------------------------------------------------------------
  
  pred_dir <- paste0("FILEPATH")
  
  inla <- rbindlist(lapply(c("raked", "unraked"), function(x) {
    pred <- fread(paste0(pred_dir,indicator,"_admin_0_", x, "_summary.csv"))
    pred[, list(type = x, ADM0_CODE, ADM0_NAME, year, mean, lower, upper)]
  }))
  
  if (file.exists(paste0(pred_dir,indicator,"_admin_0_stackers.csv"))) {
    stackers <- fread(paste0(pred_dir,indicator,"_admin_0_stackers.csv"))
    stackers <- melt(stackers, id.vars = c("ADM0_CODE", "year"), value.name = "mean", variable.name = "type")
    levels(stackers$type) <- paste("stackers:", levels(stackers$type))
  } else {
    stackers <- inla[type == "stacker",] # empty data frame, but with the correct columns
    stackers[, type := factor(type)]
  }
  
  pred <- rbind(inla, stackers, fill=T)
  pred[, ADM0_NAME := ADM0_NAME[1], ADM0_CODE]
  pred[, type := factor(type, levels = c("raked", "unraked", rev(levels(stackers$type))))]
  rm(pred_dir, inla, stackers)
  
  # input data aggregated to admin0
  input_data <- fread(paste0("FILEPATH"))
  input_data <- input_data[, list(svy_id, country, source, year, prev = get(indicator) / N, weight)]
  input_data[, svy_id := as.numeric(svy_id)]
  input_data <- unique(input_data[, list(year,
                        prev = weighted.mean(prev, weight)),
                 by = 'svy_id,country,source'])
  
  ## Make plots --------------------------------------------------------------------------------------
  
  # get GAUL to iso mapping
  loc_codes <- get_location_code_mapping(shapefile_version=shapefile_version)
  loc_codes <- loc_codes[, list(location_id = loc_id, ADM_CODE)]
  
  source("FILEPATH")
  loc <- get_location_metadata(location_set_id = 1, gbd_round_id = 5)
  loc <- merge(loc[location_type ==  "admin0",], loc_codes)
  
  # loop over countries and make plots
  dir.create(paste0('FILEPATH'), showWarnings = F)
  pdf(paste0("FILEPATH"), width = 14, height = 8)
  for (cc in na.omit(pred[order(ADM0_NAME), unique(ADM0_CODE)])) {
    iso <- loc[ADM_CODE ==  cc, ihme_loc_id]
    name <- loc[ADM_CODE ==  cc, location_name]
    if(length(iso)==0)next()
    
    # main data and estimates plot
    plot_colors <- brewer.pal(nlevels(pred$type), "Set2")
    p1 <- ggplot() +
      geom_line(data = pred[ADM0_CODE ==  cc,], aes(x = year, y = mean, color = type, size = type)) +
      scale_color_manual(values = plot_colors) +
      scale_size_manual(values = c(1.5, 1.5, rep(0.5, length(plot_colors) - 2))) +
      scale_x_continuous(limits = range(pred$year) + c(-0.5, 0.5), expand = c(0, 0)) +
      scale_y_continuous(labels = function(x) format(x, nsmall = 3)) +
      labs(y = paste0(indicator, 'coverage'), title = name) +
      theme_bw() + theme(legend.position = "bottom", legend.direction = "horizontal",
                         legend.title = element_blank(), legend.margin = margin(0, 0, 0, 0, "cm"),
                         axis.title.x = element_blank(), plot.margin = unit(rep(0.5, 4), "cm"))
    
    p2 <- ggplot() +
      geom_ribbon(data = pred[ADM0_CODE ==  cc & type == "raked",], aes(x = year, ymin = lower, ymax = upper, fill = type), color = NA, alpha = 0.2, show.legend = F) +
      geom_line(data = pred[ADM0_CODE ==  cc & type %in% c("raked", "unraked"),], aes(x = year, y = mean, color = type), size = 1.5, show.legend = F) +
      geom_point(data = input_data[country ==  iso,], aes(x = year, y = prev), color = "gray40", size = 2, position = position_jitter(width = 0.1)) +
      scale_color_manual(values = plot_colors[1:2]) +
      scale_fill_manual(values = plot_colors[1]) +
      scale_x_continuous(limits = range(pred$year) + c(-0.5, 0.5), expand = c(0, 0)) +
      scale_y_continuous(labels = function(x) format(x, nsmall = 3)) +
      guides(shape = guide_legend(override.aes = list(cex = 2))) +
      labs(x = "Year", y = paste0(indicator, 'coverage')) +
      theme_bw() + theme(legend.position = "top", legend.direction = "horizontal",
                         legend.title = element_blank(), legend.margin = margin(0, 0, 0, 0, "cm"),
                         plot.margin = unit(rep(0.5, 4), "cm"))
    
    # raking factors plot
    fdata <- pred[ADM0_CODE == cc, list(raking_factor = mean[type == "raked"] / mean[type == "unraked"]), year]
    p3 <- ggplot(fdata, aes(x = year, y = raking_factor)) +
      geom_hline(yintercept = 1) +
      geom_hline(yintercept = c(0.75, 1.25), linetype = 2) +
      geom_line(color = "red") +
      scale_x_continuous(limits = range(pred$year) + c(-0.5, 0.5), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, 2), expand = c(0, 0)) +
      labs(x = "Year", y = "Raking factor") +
      theme_bw() + theme(plot.margin = unit(c(0.5, 1, 0.5, 0.5), "cm"))
    
    # data table
    find_year_ranges <- function(run){
      rundiff <- c(1, diff(run))
      difflist <- split(run, cumsum(rundiff != 1))
      ranges <- unlist(lapply(difflist, function(x){
        if (length(x) == 1) as.character(x) else paste0(x[1], "-", substr(x[length(x)], 3, 4))
      }), use.names = FALSE)
      paste(ranges, collapse = ", ")
    }
    
    tab <- unique(input_data[country == iso, list(svy_id, source, year)])[order(svy_id, year)]
    tab <- tab[, list(years = paste(strwrap(find_year_ranges(year), width = 15), collapse = "\n")), by = 'svy_id,source']
    if (nrow(tab) ==  0) tab <- data.table(svy_id = "", source = "", years = "")
    tab <- tableGrob(tab, rows = NULL,
                     theme = ttheme_default(core = list(fg_params = list(hjust = 0, x = 0.05, fontsize = 6)),
                                            colhead = list(fg_params = list(hjust = 0, x = 0.05, fontsize = 6, fontface = "bold"))))
    
    # plot chart and table into one object
    grid.newpage()
    grid.draw(arrangeGrob(p1, p2, tab, p3, layout_matrix = matrix(c(rep(c(1, 1, 2, 2, 2), 3), 3, 3, 4, 4, 4), nrow = 5)))
  }
  dev.off()
  
  return("Plots saved!")
}
