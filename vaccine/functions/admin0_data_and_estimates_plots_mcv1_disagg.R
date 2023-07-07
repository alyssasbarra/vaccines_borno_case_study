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
# ages = z_list
# plot_age_as_x = FALSE
# plot_year_as_x = TRUE
# cohort_plots=FALSE



require(data.table)
require(ggplot2)
require(grid)
require(gridExtra)
require(RColorBrewer)

admin0_data_and_estimates_plots <- function(run_date, ages,
                                            plot_year_as_x = FALSE,
                                            plot_age_as_x = FALSE,
                                            cohort_plots=FALSE,
                                            Regions) {
  
  ## Load data and estimates -------------------------------------------------------------------------
  
  # national estimates
  pred_dir <- paste0("FILEPATH")
  
  rake_levels<- c('raked', 'unraked')
  
  fit <- rbindlist(lapply(rake_levels, function(x) {
    fit_age <- rbindlist(lapply(ages, function(age) {
        pred <- fread(paste0(pred_dir, "mcv1_disagg_admin_0_", x, "_bin", age, "_summary.csv"))
        pred[, list(ADM0_CODE, ADM0_NAME, year, mean, lower, upper)]
        pred[, agebin := age]
    }))
    fit_age[, type := paste0('0',x)]
    #  return(fit_age)
  }))
  
stackers_all<-data.table()  
for(age in ages){
    if (file.exists(paste0(pred_dir, "mcv1_disagg_admin_0_stackers1.csv"))) {
      stackers <- fread(paste0(pred_dir, "mcv1_disagg_admin_0_stackers1.csv"))
      stackers[,V1:=NULL]
      stackers[,region:=NULL]
      if('agebin' %in% colnames(stackers)){
        stackers <- melt(stackers, id.vars = c("ADM0_CODE", "year","agebin"), value.name = "mean", variable.name = "type")
        levels(stackers$type) <- paste("0stackers:", levels(stackers$type))
        stackers_all <-rbind(stackers_all, stackers)
        break
      }else{
      stackers <- melt(stackers, id.vars = c("ADM0_CODE", "year"), value.name = "mean", variable.name = "type")
      levels(stackers$type) <- paste("1stackers:", levels(stackers$type))
      stackers[, agebin := age]
      }
    } else {
      stackers <- fit[type == "stacker",] # empty data frame, but with the correct columns
      stackers[, type := factor(type)]
      stackers[, agebin := age]
    }
  stackers_all <-rbind(stackers_all, stackers)
}
stackers<-copy(stackers_all)
remove(stackers_all)
stackers<-stackers[ADM0_CODE %in% unique(fit$ADM0_CODE)]
  
  pred <- rbind(fit, stackers, fill = T)
  pred[, ADM0_NAME := ADM0_NAME[1], ADM0_CODE]
  pred[, type := factor(type, levels = c("0raked", "0unraked", rev(levels(stackers$type))))]
  pred[type %in% c('0raked', '0unraked'),stacker:=FALSE]
  pred[type %!in% c('0raked', '0unraked'),stacker:=TRUE]
  rm(pred_dir, fit, stackers)
  
  # input data
  try(input_data <- readRDS(paste0("FILEPATH")))
  try(input_data <- read.csv(paste0("FILEPATH")))
  input_data<-as.data.table(input_data)

    input_data <- input_data[, list(nid=svy_id, country, source, year, agebin=age_bin,  mcv1_disagg = mcv1_disagg / N,  agg_weight, N,region)]
 # input_data[, type := factor(type, levels = c( "ANC","Survey microdata", "Survey report", "Literature review"))]
  input_data[, nid := as.numeric(nid)]
#  input_data[type == "ANC", nid := floor(nid / 10000)]
  
  input_data <- input_data[, list(year=floor(median(year, na.rm = T)),
                                  mcv1_disagg = weighted.mean(mcv1_disagg, agg_weight*N),
                                  N=sum(N*agg_weight)),
                           by = 'nid,source,country,agebin,region']
  ########Grab country shapes
  #Adm0 outline
  shape_3<-read_sf('FILEPATH')
  shape_12<-read_sf('FILEPATH')
  shape<-rbind(shape_3, shape_12)
  shape<-shape[shape$ADM0_CODE %in% unique(pred$ADM0_CODE),]

  
  ## Make plots --------------------------------------------------------------------------------------
  
  # get GAUL to iso mapping
  loc_codes <- get_location_code_mapping(shapefile_version=shapefile_version)
  loc_codes <- loc_codes[, list(location_id = loc_id, ADM_CODE)]
  
  source("FILEPATH")
  loc <- get_location_metadata(location_set_id = 1, gbd_round_id = 5)
  loc <- merge(loc[location_type ==  "admin0",], loc_codes)
  
  
  
  #  micro[, sex_id := as.factor(sex_id)]
  #  levels(micro$sex_id) <- c('men', 'women')
  
  #For each subregion, make lists of all the countries that fell into the data for modeling that subregion
 # subregions<-fread(paste0('FILEPATH'))
  region_assignments <- data.table()
  for(s in Regions ){
    gadm_codes <- get_adm0_codes(s,shapefile_version = shapefile_version,adm0_type='gadm')
  #  if(is.null(gadm_codes)) stop('empty gadm codes')
    region_assignments <- rbind(region_assignments, data.table(ADM0_CODE=gadm_codes, region=s))
  }

  
 # if(length(z_list)==6){
  input_data[,agebin:=as.factor(agebin)]
  levels(input_data$agebin) <- c(#'6-8m',
                                 '9-11m','1y','2y','3y','4y')
  pred[,agebin:=as.factor(agebin)]
  levels(pred$agebin) <- c(#'6-8m',
                           '9-11m','1y','2y','3y','4y')
  #}

  if(ages == 0) input_data$agebin <- 0
  pred <- pred[ADM0_CODE %!in% c(94,69,241),] #Drop some countries breaking the function
  
  # main data and estimates plot
  plot_colors <- c(brewer.pal(8, "Set2")[1:2],brewer.pal(8, "Set1"),brewer.pal(8, "Set2")[3:8],brewer.pal(8, "Set3"))
  
  pred<-pred[!is.na(ADM0_CODE)]
  
  #Plot year as X ------------------------------------------------------------------------------------
  if (plot_year_as_x == TRUE) {
    message('Doing plot year as x plots now')
    # loop over countries and make plots
    dir.create(paste0('FILEPATH'), showWarnings = F)
    pdf(paste0("FILEPATH"), width = 14, height = 8)
    for (cc in pred[order(ADM0_NAME), unique(ADM0_CODE)]) {
      if(cc==105|cc==214) next
      message(paste0('Plotting country ', cc))
        iso <- loc[ADM_CODE ==  cc, ihme_loc_id]
        name <- loc[ADM_CODE ==  cc, location_name]
        
        
        p2 <- 
          ggplot() +
          geom_ribbon(data = pred[ADM0_CODE ==  cc & type %in% c("0raked", "0unraked"),], 
                      aes(x = year, ymin = lower, ymax = upper, fill = type), 
                      color = NA, alpha = 0.2) +
          geom_line(data = pred[ADM0_CODE ==  cc,],
                    aes(x = year, y = mean, color = type, lty=stacker), size = 1) +
          geom_point(data = input_data[region == region_assignments[ADM0_CODE==cc]$region,], 
                     aes(x = year, y = mcv1_disagg, size=N, color=country), alpha=0.8) +
          geom_point(data = input_data[country ==  iso & region == region_assignments[ADM0_CODE==cc]$region,], 
                     aes(x = year, y = mcv1_disagg, size=N)) +

          facet_wrap(.~agebin)+
          scale_color_manual(values = plot_colors) +
          scale_fill_manual(values = plot_colors) +
          scale_x_continuous(limits = range(pred$year) + c(-0.5, 0.5), expand = c(0, 0)) +
       #   guides(shape = guide_legend(override.aes = list(cex = 2))) +
          labs(x = "Year", y = "MCV1 coverage", size='Sample size',title = paste0(name), color='', fill='') +
          theme_bw() + theme(legend.position = "bottom", legend.direction = "horizontal",
                             legend.margin = margin(0, 0, 0, 0, "cm"),
                             plot.margin = unit(rep(0.5, 4), "cm"))
        
        #print(p2)
        p3 <-ggplot()+geom_sf(data=shape[shape$iso3 %in% unique(input_data[region == region_assignments[ADM0_CODE==cc]$region,]$country),],
                              aes(fill=iso3))+
          geom_sf(data=shape[shape$ADM0_CODE %in% unique(pred[region == region_assignments[ADM0_CODE==cc]$region,]$ADM0_CODE),],
                  fill=NA,size=1.5)+
          geom_sf(data=shape[shape$ADM0_CODE ==cc,],
                  fill='black',size=1.5)+
          theme_void()+
          labs(color='', fill='')+
          theme(legend.position = "none")+
          scale_fill_manual(values=plot_colors[(length(unique(pred$type))+1):
                                                 (length(unique(pred$type))+length(unique(input_data[region == region_assignments[ADM0_CODE==cc]$region,]$country)))])
    
        
        gp <- ggplotGrob(p2)
        
        library(gtable)
        
        # visual check of gp's layout (in this case, it has 21 rows, 15 columns)
     #   gtable_show_layout(gp)
        
        empty.area <- gtable_filter(gp, "panel", trim = F)
        empty.area <- empty.area$layout[sapply(empty.area$grob,
                                               function(x){class(x)[[1]]=="zeroGrob"}),]
        
        empty.area$t <- empty.area$t - 1 #extend up by 1 cell to cover facet header
        empty.area$b <- empty.area$b + 1 #extend down by 1 cell to cover x-axis
        
        gp0 <- gtable_add_grob(x = gp,
                               grobs = ggplotGrob(p3),
                               t = min(empty.area$t), #16 in this case
                               l = min(empty.area$l), #8
                               b = max(empty.area$b), #18
                               r = max(empty.area$r))
        grid::grid.draw(gp0)
        plot.new()
        
        
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
      if(cc==105|cc==214) next
      message(paste0('Plotting country ', cc)) 
      iso <- loc[ADM_CODE ==  cc, ihme_loc_id]
      name <- loc[ADM_CODE ==  cc, location_name]
      
      p2 <- 
        ggplot() +
        geom_ribbon(data = pred[ADM0_CODE ==  cc & type %in% c("0raked", "0unraked"),], aes(x = as.numeric(agebin), ymin = lower, ymax = upper, fill = type), color = NA, alpha = 0.2) +
        geom_line(data = pred[ADM0_CODE ==  cc,], aes(x = as.numeric(agebin), y = mean, color = type, lty=stacker), size = 1) +
        geom_point(data = input_data[region == region_assignments[ADM0_CODE==cc]$region,], 
                   aes(x = as.numeric(agebin), y = mcv1_disagg, size=N, color=country), alpha=0.8) +
        geom_point(data = input_data[country ==  iso & region == region_assignments[ADM0_CODE==cc]$region,], 
                   aes(x = as.numeric(agebin), y = mcv1_disagg, size=N)) +
        facet_wrap(.~year) +
        scale_color_manual(values = plot_colors) +
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
      if(cc==105|cc==214) next
      message(paste0('Plotting country ', cc))
      iso <- loc[ADM_CODE ==  cc, ihme_loc_id]
      name <- loc[ADM_CODE ==  cc, location_name]
      
      p2 <- 
        ggplot() +
        geom_ribbon(data = pred[ADM0_CODE ==  cc & type %in% c("0raked", "0unraked"),], aes(x = as.numeric(agebin), ymin = lower, ymax = upper, fill = type), color = NA, alpha = 0.2) +
        geom_line(data = pred[ADM0_CODE ==  cc,], aes(x = as.numeric(agebin), y = mean, color = type,lty=stacker), size = 1) +
        geom_point(data = input_data[region == region_assignments[ADM0_CODE==cc]$region,], 
                   aes(x = as.numeric(agebin), y = mcv1_disagg, size=N, color=country), alpha=0.8) +
        geom_point(data = input_data[country ==  iso & region == region_assignments[ADM0_CODE==cc]$region,], 
                   aes(x = as.numeric(agebin), y = mcv1_disagg, size=N)) +
        
        facet_wrap(.~birth_year) +
        scale_color_manual(values = plot_colors) +
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
