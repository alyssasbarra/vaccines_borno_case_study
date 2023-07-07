####################################################################################################
## Description:   Make plots of national-level temporal trends and data.
##
## Inputs:        
##                Collapsed MBG input data 
##                Unraked admin0 predictions 
##                Raked admin0 predictions 
##
## Output:        PDF of plots 
##############################

require(data.table)
require(ggplot2)
require(grid)
require(gridExtra)
require(RColorBrewer)

admin0_data_and_estimates_model_comparisons <- function(run_date, ages, 
                                            plot_year_as_x = FALSE, 
                                            plot_age_as_x = FALSE,
                                            cohort_plots=FALSE,
                                            models,
                                            model_dates,
                                            countries = NULL) {
  
  ## Load data and estimates -------------------------------------------------------------------------
  fit_all<-NULL
  stackers_all <-NULL
  for (i in 1:length(models)){
  # national estimates
  pred_dir <- paste0("FILEPATH")
  
  rake_levels<- c('raked', 'unraked')
  
  fit <- rbindlist(lapply(rake_levels, function(x) {
    fit_age <- rbindlist(lapply(ages, function(age) {
      pred <- fread(paste0(pred_dir, "mcv1_disagg_admin_0_", x, "_bin", age, "_summary.csv"))
      pred[, list(ADM0_CODE, ADM0_NAME, year, mean, lower, upper)]
      pred[, agebin := age]
    }))
    fit_age[, type := x]
    #  return(fit_age)
  }))    
  fit[,model := factor(models[i])]
  fit_all<-rbind(fit_all, fit)
  
  
  stackers <- rbindlist(lapply(ages, function(age) {
    if (file.exists(paste0(pred_dir, "mcv1_disagg_admin_0_bin", age, "_stackers.csv"))) {
      stackers <- fread(paste0(pred_dir, "mcv1_disagg_admin_0_bin", age, "_stackers.csv"))
      stackers <- melt(stackers, id.vars = c("ADM0_CODE", "year"), value.name = "mean", variable.name = "type")
      levels(stackers$type) <- paste("stackers:", levels(stackers$type))
      stackers[, agebin := age]
    } else {
      stackers <- fit[type == "stacker",] # empty data frame, but with the correct columns
      stackers[, type := factor(type)]
      stackers[, agebin := age]
    }
  }))
  stackers[, model:=factor(models[i])]
  stackers_all<-rbind(stackers_all, stackers)
  }
  
  pred <- rbind(fit_all, stackers, fill = T)
  pred[, ADM0_NAME := ADM0_NAME[1], ADM0_CODE]
  pred[, type := factor(type, levels = c("raked", "unraked", rev(levels(stackers$type))))]
  rm(pred_dir, fit, stackers)
  
  #Limit to countries with both models
  if(!is.null(countries)){
    pred <- pred[ADM0_CODE %in% countries]
  }
  
  # input data
  try(input_data <- readRDS(paste0("FILEPATH")))
  try(input_data <- read.csv(paste0("FILEPATH")))
  input_data<-as.data.table(input_data)
  
  input_data <- input_data[, list(nid=svy_id, country, source, year, agebin=age_bin,  mcv1_disagg = mcv1_disagg / N,  agg_weight, N)]
  # input_data[, type := factor(type, levels = c( "ANC","Survey microdata", "Survey report", "Literature review"))]
  input_data[, nid := as.numeric(nid)]
  #  input_data[type == "ANC", nid := floor(nid / 10000)]
  
  input_data <- input_data[, list(year=floor(median(year, na.rm = T)),
                                  mcv1_disagg = weighted.mean(mcv1_disagg, agg_weight*N),
                                  N=sum(N*agg_weight)),
                           by = 'nid,source,country,agebin']
  
  
  
  ## Make plots --------------------------------------------------------------------------------------
  
  # get GAUL to iso mapping
  loc_codes <- get_location_code_mapping(shapefile_version=shapefile_version)
  loc_codes <- loc_codes[, list(location_id = loc_id, ADM_CODE)]
  
  source("FILEPATH")
  loc <- get_location_metadata(location_set_id = 1, gbd_round_id = 5)
  loc <- merge(loc[location_type ==  "admin0",], loc_codes)
  
  
  
  #  micro[, sex_id := as.factor(sex_id)]
  #  levels(micro$sex_id) <- c('men', 'women')
  
  
  
  
  # if(length(z_list)==6){
  input_data[,agebin:=as.factor(agebin)]
  levels(input_data$agebin) <- c(#'6-8m',
    '9-11m','1y','2y','3y','4y')
  pred[,agebin:=as.factor(agebin)]
  levels(pred$agebin) <- c(#'6-8m',
    '9-11m','1y','2y','3y','4y')
  #}
  
  if(ages == 0) input_data$agebin <- 0
  pred <- pred[ADM0_CODE %!in% c(94,69,241),]
  
  # main data and estimates plot
  plot_colors <- brewer.pal(nlevels(pred$type), "Set2")
  
  #Plot year as X ------------------------------------------------------------------------------------
  if (plot_year_as_x == TRUE) {
    message('Doing plot year as x plots now')
    # loop over countries and make plots
    dir.create(paste0('FILEPATH'), showWarnings = F)
    pdf(paste0("FILEPATH"), width = 14, height = 8)
    for (cc in pred[order(ADM0_NAME), unique(ADM0_CODE)]) {
      message(paste0('Plotting country ', cc))
      iso <- loc[ADM_CODE ==  cc, ihme_loc_id]
      name <- loc[ADM_CODE ==  cc, location_name]
      
      
      p2 <- 
        ggplot() +
        geom_line(data = pred[ADM0_CODE ==  cc & type %in% c("raked", "unraked"),],
                  aes(x = year, y = mean, color = type,lty=model), size = 1.5) +
        geom_point(data = input_data[country ==  iso,], 
                   aes(x = year, y = mcv1_disagg, size=N)) +
        
        facet_wrap(.~agebin)+
        #geom_point(data = report[country ==  iso,], aes(x = int_year, y = hiv_prev_disagg, shape = type), color = "black", size = 5) +
        scale_color_manual(values = plot_colors[1:2]) +
        scale_fill_manual(values = plot_colors[1:2]) +
        #   scale_shape_manual(values = c(16, 18, 15, 8), drop = F) +
        scale_x_continuous(limits = range(pred$year) + c(-0.5, 0.5), expand = c(0, 0)) +
        #   guides(shape = guide_legend(override.aes = list(cex = 2))) +
        labs(x = "Year", y = "MCV1 coverage", size='Sample size',title = paste0(name), color='', fill='') +
        theme_bw() + theme(legend.position = "bottom", legend.direction = "horizontal",
                           legend.margin = margin(0, 0, 0, 0, "cm"),
                           plot.margin = unit(rep(0.5, 4), "cm"))
      
      print(p2)
    }
    dev.off()
    
    return("Year as x plots saved!")
  }
  
  #Plot age as X ------------------------------------------------------------------------------------  
  if (plot_age_as_x == TRUE) {
    
    message('Doing plot age as x plots now')
    # loop over countries and make plots
    dir.create(paste0('FILEPATH'), showWarnings = F)
    pdf(paste0("FILEPATH"), width = 14, height = 8)
    for (cc in pred[order(ADM0_NAME), unique(ADM0_CODE)]) {
      message(paste0('Plotting country ', cc)) 
      iso <- loc[ADM_CODE ==  cc, ihme_loc_id]
      name <- loc[ADM_CODE ==  cc, location_name]
      
      p2 <- 
        ggplot() +
        geom_line(data = pred[ADM0_CODE ==  cc & type %in% c("raked", "unraked"),], aes(x = as.numeric(agebin), y = mean, color = type, lty=model), size = 1.5) +
        geom_point(data = input_data[country ==  iso,],
                   aes(x = as.numeric(agebin), y = mcv1_disagg, size=N), 
                   color = "gray40",  position = position_jitter(width = 0.1))+
        facet_wrap(.~year) +
        scale_color_manual(values = plot_colors[1:2]) +
        scale_fill_manual(values = plot_colors[1:2]) +
        # scale_shape_manual(values = c(16, 18, 15, 8), drop = F) +
        scale_x_continuous(limits = range(as.numeric(pred$agebin)+ c(-.2,.2)) , expand = c(0, 0), 
                           labels = unique(pred$agebin))+
        guides(alpha = FALSE) +
        labs(x = "Age bin", y = "MCV1 coverage", size='Sample size',title = paste0(name), color='', fill='') +
        theme_bw() + theme(legend.position = "bottom", legend.direction = "horizontal",
                           legend.margin = margin(0, 0, 0, 0, "cm"),
                           plot.margin = unit(rep(0.5, 4), "cm"))
      
      print(p2)
    }
    dev.off()
    
    return("Age as x plots saved!")
    
  }  
  
  #Cohort plots---------------------------------------------------------------------------------
  
  
  
  
  
  
  if(cohort_plots == TRUE){
    pred[,birth_year:=as.numeric(agebin)]
    pred[birth_year == 1, birth_year:=year]
    pred[birth_year == 2, birth_year:=year-1]
    pred[birth_year == 3,     birth_year:=year-2]
    pred[birth_year == 4,     birth_year:=year-3]
    pred[birth_year == 5,     birth_year:=year-4]
    
    input_data[,birth_year:=as.numeric(agebin)]
    input_data[birth_year == 1, birth_year:=year]
    input_data[birth_year == 2, birth_year:=year-1]
    input_data[birth_year == 3,     birth_year:=year-2]
    input_data[birth_year == 4,     birth_year:=year-3]
    input_data[birth_year == 5,     birth_year:=year-4]
    
    dir.create(paste0('FILEPATH'), showWarnings = F)
    pdf(paste0("FILEPATH"), width = 14, height = 8)
    for (cc in pred[order(ADM0_NAME), unique(ADM0_CODE)]) {
      message(paste0('Plotting country ', cc))
      iso <- loc[ADM_CODE ==  cc, ihme_loc_id]
      name <- loc[ADM_CODE ==  cc, location_name]
      
      p2 <- 
        ggplot() +
        geom_line(data = pred[ADM0_CODE ==  cc & type %in% c("raked", "unraked"),], aes(x = as.numeric(agebin), y = mean, color = type,lty=model), size = 1.5) +
        geom_point(data = input_data[country ==  iso,],
                   aes(x = as.numeric(agebin), y = mcv1_disagg, size=N), 
                   color = "gray40",  position = position_jitter(width = 0.1))+
        facet_wrap(.~birth_year) +
        scale_color_manual(values = plot_colors[1:2]) +
        scale_fill_manual(values = plot_colors[1:2]) +
        # scale_shape_manual(values = c(16, 18, 15, 8), drop = F) +
        scale_x_continuous(limits = range(as.numeric(pred$agebin)+ c(-.2,.2)) , expand = c(0, 0), 
                           labels = unique(pred$agebin))+
        guides(alpha = FALSE) +
        labs(x = "Age bin", y = "MCV1 coverage", size='Sample size',title = paste0(name, ' - cohorts by birth year'), color='', fill='') +
        theme_bw() + theme(legend.position = "bottom", legend.direction = "horizontal",
                           legend.margin = margin(0, 0, 0, 0, "cm"),
                           plot.margin = unit(rep(0.5, 4), "cm"))
      
      print(p2)
    }
    dev.off()
    
    return("Cohort plots saved!")
  }
  
  
}
