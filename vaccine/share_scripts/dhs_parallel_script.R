# Code to generate plots of admin1 dhs estimates vs admin1 aggregates from our models
# Written by: USERNAME
# Adapted by USERNAME
# DATE
# 
# Must pass in 2 files
# 1) SPDF of shapefiles corresponding to DHS Data
# 2)  Collapsed .CSV of just MACRO DHS data

## Pull draws by country and save template rasters by country in memory
###############################################################################
## SETUP
###############################################################################

## clear environment
rm(list=ls())

## Set repo location and indicator group
user               <- Sys.info()['user']
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
                          sprintf('FILEPATH',root),
                          sprintf('FILEPATH',root))
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
source(paste0(indic_repo,'FILEPATH'))

package_list <- c('SDMTools', 'survey', 'pbapply', 'readstata13', 'foreign', 'rgeos', 'data.table','raster','rgdal','INLA','seegSDM','seegMBG','plyr','dplyr', 'foreach', 'doParallel','ggplot2')
for(package in package_list) {
  library(package, lib.loc = package_lib, character.only=TRUE)
}

## Manually set indicator and run_date for manuscript prep
load_from_qsub() #indicator, reg

reg <- region
indicator_name <- indicator
#indicator_name <- 'dpt3_cov'
#indicator <- indicator_name
#indicator_group <- 'vaccine'
pop_measure<-'a0004t'
#run_date <- 'RUN_DATE'
#input_date <- 'DATE'                #date to read from FILEPATH input data
#reg <- 'wssa'
#use_raked='_raked'
#use_raked='_raked'
#cohort <- "12_23_months" to subset to just 1-year old
cohort <- "12_23_months"
year_list <- c(2000:2016)

message(paste0("input_date: ", input_date))
message(paste0("indicator: ", indicator))
message(paste0("reg: ", reg))

pull_country_draws <- function(reg, periods, raked = "", in_dir, pop_measure, start_year, end_year, admin2_shapes, all_region_pops, 
                                shapefile_version=modeling_shapefile_version,
                                subnational_raking=subnational_raking) {
  
  ## Load simple_raster and pop_raster into memory for this GAUL_CODE
  gaul_list <- get_adm0_codes(reg, shapefile_version=modeling_shapefile_version, subnational_raking=subnational_raking)
  simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list,
                                             buffer = 0.4,
                                             subset_only = TRUE, 
                                             shapefile_version=modeling_shapefile_version)
  subset_shape   <- simple_polygon_list[[1]]
  simple_polygon <- simple_polygon_list[[2]]
  raster_list    <- build_simple_raster_pop(subset_shape)
  simple_raster  <- raster_list[['simple_raster']]
  pop_raster     <- raster_list[['pop_raster']]
  
  ## Load draws into memory for this GAUL_CODE.
  message(paste0('Loading draws for ', reg,'...'))
  load(paste0(in_dir, '/', indicator, raked, '_cell_draws_eb_bin0_',reg,'_0.RData'))
  if (raked == ''){
    draws <- as.data.table(cell_pred)
    rm(cell_pred)
  }
  else{
    draws <- as.data.table(raked_cell_pred)
    rm(raked_cell_pred)
  }
  
  ## Get dt or draws for this region indexed by GAUL_CODE.
  cell_idx <- seegSDM:::notMissingIdx(simple_raster)
  coords <- xyFromCell(simple_raster, seegSDM:::notMissingIdx(simple_raster))
  template <- raster::extract(simple_raster,coords)
  draws_by_gaul <- draws[, GAUL_CODE := rep(template,periods)]
  
  ## Pull out draws/templates for each GAUL_CODE in this region.
  message('Pulling draws by country for ', reg, '...')
  make_country_list <- function(gaul, reg_simple_raster) {
    message(gaul)
    # Save draws
    country_draws <- draws_by_gaul[GAUL_CODE==gaul,]
    country_draws <- country_draws[, GAUL_CODE := NULL]
    # Save a template raster and population raster
    country_simple_raster <- reg_simple_raster
    country_simple_raster[country_simple_raster!=gaul] <- NA
    # gaul_list <- gaul
    # country_shape  <- subset_shape[subset_shape@data$GAUL_CODE %in% gaul_list, ]
    # raster_list    <- build_simple_raster_pop(country_shape)
    # simple_raster  <- raster_list[['simple_raster']]
    # pop_raster     <- raster_list[['pop_raster']]
    # Get population brick for all periods
    pop_raster_annual  <- all_region_pops
    pop_raster_annual  <- crop(pop_raster_annual, extent(country_simple_raster))
    pop_raster_annual  <- setExtent(pop_raster_annual, country_simple_raster)
    pop_raster_annual  <- mask(pop_raster_annual, country_simple_raster)
    # Get raster of admin2 codes
    admin_level <- 2
    shapes <- admin2_shapes
    cropped_shapes <- crop(shapes, extent(country_simple_raster), snap="out")
    ## Fix rasterize
    initial_raster <- rasterize(cropped_shapes, country_simple_raster, field = paste0('ADM', admin_level,'_CODE'))
    if(length(cropped_shapes[!cropped_shapes[[paste0('ADM', admin_level,'_CODE')]]%in%unique(initial_raster),])!=0) {
      rasterized_shape <- merge(rasterize(cropped_shapes[!cropped_shapes[[paste0('ADM', admin_level,'_CODE')]]%in%unique(initial_raster),], country_simple_raster, field = paste0('ADM', admin_level,'_CODE')), initial_raster)
    }
    if(length(cropped_shapes[!cropped_shapes[[paste0('ADM', admin_level,'_CODE')]]%in%unique(initial_raster),])==0) {
      rasterized_shape <- initial_raster
    }
    masked_shapes <- mask(x=rasterized_shape, mask=country_simple_raster)
    # Add items to list
    return_list <- list(country_draws,
                        pop_raster_annual,
                        country_simple_raster,
                        masked_shapes)
    names(return_list) <- c(paste0('draws_', gaul), paste0('pops_', gaul), paste0('simple_', gaul), paste0('admin2_', gaul))
    return(return_list)
  }
  reg_gaul_list <- get_adm0_codes(reg, shapefile_version=modeling_shapefile_version, subnational_raking=subnational_raking)
  reg_gaul_list <- reg_gaul_list[!(reg_gaul_list %in% 40762)]
  return_list <- lapply(reg_gaul_list, make_country_list, 
                        reg_simple_raster = simple_raster)
  names(return_list) <- paste0('list_', reg_gaul_list)
  
  ## Return list of draws, pops, templates
  return(return_list)
  
}

#indicator_name <- as.character(commandArgs()[4])
#run_date <- as.character(commandArgs()[5])
#regions <- as.character(commandArgs()[6])
#use_raked <- as.character(commandArgs()[7])
#pop_measure <- as.character(commandArgs()[8])
#indicator_group  <- as.character(commandArgs()[10])

in_dir  <- paste0('FILEPATH')
default_rf_path <- paste0('FILEPATH')
all_rfs <- fread(default_rf_path)
## Define path to .tif of results and raked results, load in raster bricks.
default_raked_results_path <- paste0(in_dir, '/', indicator, '_', reg,'_raked_mean_raster.tif')
results_raked <- brick(default_raked_results_path)
default_results_path <- paste0(in_dir, '/', indicator, '_', reg, '_mean_raster.tif')
results <- brick(default_results_path)

## Load regional pops
simple_polygon_list <- load_simple_polygon(gaul_list = get_adm0_codes(reg, shapefile_version=modeling_shapefile_version, subnational_raking=subnational_raking),
                                           buffer = 0.4,
                                           subset_only = FALSE,
                                           shapefile_version=modeling_shapefile_version)
subset_shape   <- simple_polygon_list[[1]]
simple_polygon <- simple_polygon_list[[2]]
pop_raster_annual <- suppressMessages(suppressWarnings(load_and_crop_covariates_annual(covs = 'worldpop_raked',
                                                                                       measures = pop_measure,
                                                                                       simple_polygon = simple_polygon,
                                                                                       start_year  = min(year_list),
                                                                                       end_year    = max(year_list),
                                                                                       interval_mo = 12,
                                                                                       agebin      = 1)))

## Make main list
total_periods <- length(names(results))
#regions <- get_output_regions(in_dir)

admin_level <- 0
admin0 <- shapefile(paste0("FILEPATH"))

## Load admin2 raster
admin_level <- 2
shapes <- shapefile(paste0("FILEPATH"))
## Fix several african disputed territories (assign them to a country so we can include them in our estimates and they dont drop out)
shapes[[paste0('ADM', admin_level,'_CODE')]][shapes[[paste0('ADM', admin_level,'_CODE')]]==61013]=133
shapes[[paste0('ADM', admin_level,'_CODE')]][shapes[[paste0('ADM', admin_level,'_CODE')]]==40760]=40765
shapes[[paste0('ADM', admin_level,'_CODE')]][shapes[[paste0('ADM', admin_level,'_CODE')]]==40762]=145
shapes[[paste0('ADM', admin_level,'_CODE')]][shapes[[paste0('ADM', admin_level,'_CODE')]]==102  ]=6

main_list <- lapply(reg, pull_country_draws,
                      periods = total_periods,
                      in_dir = in_dir,
                      pop_measure = pop_measure,
                      start_year =  min(year_list),
                      end_year =  max(year_list),
                      admin2_shapes = shapes,
                      all_region_pops = pop_raster_annual[[1]],
                      raked = use_raked)
main_list <- do.call(c, unlist(main_list, recursive=FALSE))

## Load premade .csv of DHS admin1 aggregates for edu_mean and an SPDF with all the relevant location_code/shapefile entries.
# load('FILEPATH')
# dhs_admin1_data <- fread(paste0('FILEPATH'))

#### Load
poly_shapes_all <- readRDS(paste0("FILEPATH"))
DHS_data <- fread(paste0("FILEPATH"))

## subset to 12_23_month cohort

# ##rename
names <- names(DHS_data)
names[7]<- "nid"
names[5] <-"shapefile"
names[4] <- "location_code"
names[6] <- "survey_series"
names[3] <- "iso3"
names[2] <- "cohort_year"
names[8] <- "year"
names(DHS_data) <- names

# Subset to 12-23 months
DHS_data <- DHS_data[cohort_year == year,]

#DHS_data[,outcome:=dpt_dose_3/N]

## Choose which NIDs to compare (all for this region)
loc_names <- fread("FILEPATH")
region_iso3s <- loc_names[GAUL_CODE %in% get_gaul_codes(reg), ihme_lc_id]
region_nid_list <- unique(DHS_data[iso3 %in% region_iso3s, nid])

pull_nid_draws <- function(this_nid, iso3_results = iso3_results, compare_spdf = full_spdf, simple = simple, all_draws, all_pops) {
  #### Need to re-write to get more than just 12-23 month cohort
  message(this_nid)
  nid_results <- iso3_results[year >= min(year_list) & year <= max(year_list) & nid == this_nid, ]
  #all_data <- merge(compare_spdf, nid_results, by=c('location_code','shapefile'))
  
  ### Testing how to create 1 data frame for all years
  years <- unique(nid_results$year)
  
  if(length(years) > 1){
    all_data <- lapply(years, function(x) {merge(compare_spdf, nid_results[year==x,], by=c('location_code','shapefile'))})
    all_data <- do.call(rbind,all_data)
  }
  else{
    all_data <- merge(compare_spdf, nid_results, by=c('location_code','shapefile'))
  }
  
  
  all_data <- all_data[!is.na(all_data@data$outcome),]
  all_data$geo_mean <- 0
  all_data$geo_upper <- 0
  all_data$geo_lower <- 0
  for(shape_x in unique(all_data$location_code)) {
    
    this_poly <- all_data[all_data$location_code == shape_x, ]
    masked_simple <- mask(simple, this_poly)
    poly_idx <- cellIdx(masked_simple)
    simple_idx <- cellIdx(simple)
    simple_idx[simple_idx %in% poly_idx] <- -999
    simple_idx[!(simple_idx %in% -999)] <- 0
    
    for (y in unique(all_data$year)){
      nid_results_year <- nid_results[year==y,]
      period <- unique(nid_results_year[, year]) - 2000 + 1
      poly_draws <- all_draws[period_index == paste0('period_', period), ]
      poly_draws <- poly_draws[, pops := raster::extract(all_pops[[period]], cellIdx(simple))]
      poly_draws <- poly_draws[, poly_id := simple_idx]
      poly_draws <- poly_draws[poly_id == -999, ]
      poly_means <- poly_draws[, lapply(.SD, weighted.mean, w=pops, na.rm=TRUE), .SDcols=grep("^V", names(poly_draws)) ]
      poly_means <- poly_means[, lower := apply(.SD, 1, quantile, c(.025), na.rm=T), .SDcols=grep("^V", names(poly_means))]
      poly_means <- poly_means[, mean := apply(.SD, 1, mean), .SDcols=grep("^V", names(poly_means))]
      poly_means <- poly_means[, upper := apply(.SD, 1, quantile, c(.975), na.rm=T), .SDcols=grep("^V", names(poly_means))]
      all_data$geo_mean[all_data$location_code == shape_x & all_data$year == y] <- poly_means[, mean]
      all_data$geo_upper[all_data$location_code == shape_x & all_data$year == y] <- poly_means[, upper]
      all_data$geo_lower[all_data$location_code == shape_x & all_data$year == y] <- poly_means[, lower]
    }
  }
  this_data <- as.data.table(all_data)
  return(this_data)
  
}

pull_iso3_results <- function(this_iso3, main_list, compare_data, nperiod, full_spdf) {
  
  message(this_iso3)
  iso3_results <- compare_data[year >= min(year_list) & year <= max(year_list) & iso3 == this_iso3, ]
  gaul <- gaul_convert(this_iso3, shapefile_version=modeling_shapefile_version)
  country_pops <- main_list[[paste0('list_', gaul, '.pops_', gaul)]]
  country_pops <- crop(country_pops, extent(main_list[[paste0('list_', gaul, '.simple_', gaul)]]))
  country_pops <- setExtent(country_pops, main_list[[paste0('list_', gaul, '.simple_', gaul)]])
  country_pops <- mask(country_pops, main_list[[paste0('list_', gaul, '.simple_', gaul)]])
  message('get country draws...')
  cell_pred.dt <- as.data.table(main_list[[paste0('list_', gaul, '.draws_', gaul)]])
  cell_pred.dt[cell_pred.dt < 0] <- 0
  names <- names(cell_pred.dt)
  cols <- names(cell_pred.dt)
  period_index <- c()
  for(i in 1:nperiod) {
    period_index <- c(period_index, rep(paste0("period_", i), length(cell_pred.dt$V1)/nperiod))
  }
  pixel_id <- c(rep(1:(length(cell_pred.dt$V1)/nperiod),
                    nperiod))
  cell_pred.dt <- cbind(cell_pred.dt, period_index, pixel_id)
  simple <- main_list[[paste0('list_', gaul, '.simple_', gaul)]]
  all_nid_compare <- rbindlist(lapply(unique(iso3_results[nid!="157025" & nid!="218593", nid]), pull_nid_draws,
                                      iso3_results = iso3_results,
                                      compare_spdf = full_spdf,
                                      simple = simple,
                                      all_draws = cell_pred.dt,
                                      all_pops = country_pops))
  return(all_nid_compare)
  
}

all_data <- rbindlist(lapply(unique(DHS_data[year >= min(year_list) & nid %in% region_nid_list, iso3]), pull_iso3_results,
                             main_list = main_list,
                             compare_data = DHS_data,
                             nperiod = length(year_list),
                             full_spdf = poly_shapes_all))



## Calculate super simple bias estimate
## Fit simple lm
all_data <- all_data[!is.na(outcome) & !is.na(geo_mean), ]
bias_model <- lm(outcome ~ geo_mean, data = all_data)
## Intercept
intercept <- summary(bias_model)$coefficients[1, 1]
intercept_p <- summary(bias_model)$coefficients[1, 4]
if(intercept_p <= 0.0005) intercept_p <- 0.0005
## Slope
slope <- summary(bias_model)$coefficients[2, 1]
slope_p <- summary(bias_model)$coefficients[2, 4]
if(slope_p <= 0.0005) slope_p <- 0.0005
summary(bias_model)
## Predict fitted values for line
all_data$lm_pred <- predict(bias_model)

compare_source_title <- "test"

## Make title for plot depending on if comparing a single source (grab year/iso3) or many sources
if(length(unique(all_data[, nid]))==1) add_title <- paste0('\n', all_data[, iso3][1], ' ', min(unique(all_data[, year])),'-', max(unique(all_data[, year])), ' ', compare_source_title, '\n')
if(length(unique(all_data[, nid]))!=1) add_title <- '\n'
if(use_raked == '_raked') raked_title <- 'raked'
if(use_raked == '') raked_title <- 'unraked'
main_gg_title <- paste0(indicator,'_','DHS Admin1 vs. aggregated ', raked_title, ' MBG estimates', add_title, 'Intercept: ', round(intercept, 2), ' (p <= ', round(intercept_p, 4),')\nSlope: ', round(slope, 2), ' (p <= ', round(slope_p, 4),')')

## Plot gg for unraked
comparison_gg=ggplot(all_data, aes(x=outcome,y=geo_mean))+
  geom_abline(intercept=0,slope=1,colour='red')+
  geom_point(aes(size = N), colour='grey')+
  theme_bw()+
  xlab('Data Estimate') +
  ylab('Mean Prediction')+
  theme(strip.background = element_rect(fill="white"))+
  ggtitle(main_gg_title) +
  geom_errorbar(aes(ymin=geo_lower, ymax=geo_upper), colour="black", width=0, size=.5, alpha = 0.2) +
  geom_abline(intercept=0,slope=1,colour='red') +
  scale_size_area() +
  xlim(c(0,1)) +
  ylim(c(0,1))+
  coord_equal()

pdf(paste0('FILEPATH'))
comparison_gg
dev.off()

dir.create(paste0('FILEPATH'))
write.csv(all_data, paste0('FILEPATH'))



