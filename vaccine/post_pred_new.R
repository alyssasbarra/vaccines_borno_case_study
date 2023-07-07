
## Convenience

indicator <- rake_indicators

sharedir <- sprintf('FILEPATH')

#Regions <- c("vax_eaas","vax_name", "vax_trsa", "vax_ansa", "vax_caeu", "vax_crbn", "vax_cssa", "vax_ctam", "vax_essa", "vax_seas", "vax_soas", "vax_sssa", "vax_wssa")


#Regions <- c("vax_ansa_pcv", "vax_caeu_pcv", "vax_crct_pcv", "vax_cssa_pcv", "vax_eaas_pcv", "vax_essa_pcv", "vax_seas_pcv", "vax_some_pcv", "vax_sssa_pcv", "vax_trsa_pcv", "vax_wssa_pcv")


strata <- Regions
#######################################################################################
## Merge!! 
#######################################################################################

  rr <- c("raked", "unraked")
  rf_tab <- T
  run_summary <- T

if (("p_below" %in% summstats) & !(indicator %in% c(paste0(vaccine, "1_cov"), paste0(vaccine, "3_cov")))) {
  # Only run this for 1st and 3rd doses
  summstats <- summstats[summstats != "p_below"]
}

post_load_combine_save(regions    = Regions,
                       summstats  = summstats,
                       raked      = rr,
                       rf_table   = rf_tab,
                       run_summ   = run_summary,
                       indic      = indicator,
                       ig         = indicator_group,
                       sdir       = sharedir)

# Clean up / delete unnecessary files
clean_after_postest(indicator             = indicator,
                    indicator_group       = indicator_group,
                    run_date              = run_date,
                    strata                = strata,
                    delete_region_rasters = F)


#######################################################################################
## Combine aggreagations
#######################################################################################
holdouts <- 0
combine_aggregation(rd = run_date, indic = indicator, ig = indicator_group,
                    ages = 0, 
                    regions = strata,
                    holdouts = holdouts,
                    raked = c(T, F),
                    delete_region_files = F)


#######################################################################################
## Summarize
#######################################################################################


summarize_admins(summstats = c("mean", "lower", "upper", "cirange", "cfb"), 
                 ad_levels = c(0,1,2), 
                 raked = c(T,F))

message("Summarize p_above")
summarize_admins(indicator, indicator_group,
                 summstats = c("p_above"),
                 raked = c(T,F),
                 ad_levels = c(0,1,2),
                 file_addin = "p_0.8_or_better",
                 value = 0.8,
                 equal_to = T)



#######################################################################################
## Save everything nicely for mapping
#######################################################################################

postest_indicators <- indicator
#summstats <- summstats[-5]

  tifs_to_copy <- data.table(expand.grid(summstats, c(T), postest_indicators, stringsAsFactors = F))
  names(tifs_to_copy) <- c("measure", "raked", "indicator") 

  for (i in 1:nrow(tifs_to_copy)) {
    copy_tif_to_map_input_dir(ind = tifs_to_copy[i, "indicator"],
                              ig = indicator_group,
                              measure = tifs_to_copy[i, "measure"],
                              rd = run_date,
                              raked = tifs_to_copy[i, "raked"],
                              yl = year_list)
  }

  # Save ads
  admin_lvs <- data.table(expand.grid(summstats, c(T), c(0,1,2), postest_indicators, stringsAsFactors = F))
  admin_lvs <- rbind(admin_lvs,
                     data.table(expand.grid("psup80", 
                                            c(T,F), 
                                            c(0,1,2),
                                            postest_indicators[postest_indicators != paste0(vaccine, "1_3_rel_dropout")], 
                                stringsAsFactors = F)))

  names(admin_lvs) <- c("measure", "raked", "ad_level", "indicator") 

  for (i in 1:nrow(admin_lvs)) {
    copy_admins_to_map_input_dir(ind = admin_lvs[i, "indicator"],
                                 ig = indicator_group, 
                                 measure = as.character(admin_lvs[i, "measure"]), 
                                 rd = run_date, 
                                 raked = admin_lvs[i, "raked"], 
                                 yl = year_list, 
                                 ad_level = admin_lvs[i, "ad_level"]) 
  }







summstats <- c("mean", "cfb")
summstats <- 'mean'


#### gather list of admin codes to keep in csvs


ad0 <- fread('FILEPATH')
ad1 <- fread('FILEPATH')
ad2 <- fread('FILEPATH')



ad0 <- unique(ad0$ADM0_CODE)
ad1 <- unique(ad1$ADM1_CODE)
ad2 <- unique(ad2$ADM2_CODE)



africa <- raster('FILEPATH')


dir.create(paste0('FILEPATH'))
dir.create(paste0('FILEPATH'))



for(stat in summstats){
  message(stat)
  
  #### first, re-save csv's
  ad0_csv <- fread(paste0('FILEPATH'))
  ad1_csv <- fread(paste0('FILEPATH'))
  ad2_csv <- fread(paste0('FILEPATH'))
  
  ad0_csv <- subset(ad0_csv, ADM0_CODE %in% ad0)
  ad1_csv <- subset(ad1_csv, ADM1_CODE %in% ad1)
  ad2_csv <- subset(ad2_csv, ADM2_CODE %in% ad2)
  
  write.csv(ad0_csv, paste0('FILEPATH'))
  write.csv(ad1_csv, paste0('FILEPATH'))
  write.csv(ad2_csv, paste0('FILEPATH'))
  
  
  
  
  #### next, re-save tifs
  tif_yr <- brick(paste0('FILEPATH'))
  
  for(yr in 1:19){
    message((yr + 1999))
    year <- yr+1999
    
    single_yr <- tif_yr[[yr]]     
    single_yr <- crop(single_yr, africa)
    single_yr <- mask(single_yr, africa)
    
    writeRaster(single_yr, paste0('FILEPATH'))
  }
                    
  
  
}





#######################################################################################
## In sample metric validation
#######################################################################################

# Combine csv files only if none present
csvs <- list.files(paste0("FILEPATH"), 
                   full.names = T)

if (sum(grepl("input_data.csv", csvs)) == 0) {
  csv_main <- lapply(csvs, read.csv, stringsAsFactors = F) %>% 
                rbindlist %>%
                subset(., select = !(names(.) %in% c("X.2", "X.1"))) %>%
                as.data.table
  write.csv(csv_main, file=paste0('FILEPATH'))
}

# Get in and out of sample draws
run_in_oos <- get_is_oos_draws(ind_gp = indicator_group,
                               ind = indicator,
                               rd = run_date,
                               ind_fm = 'binomial',
                               model_domain = Regions,
                               age = 0,
                               nperiod = length(year_list),
                               yrs = year_list,
                               get.oos = as.logical(makeholdouts),
                               write.to.file = TRUE,
                               year_col = "year",
                               shapefile_version = modeling_shapefile_version) # uses unraked

## set out_dir
out_dir <- paste0("FILEPATH")
dir.create(out_dir, recursive = T, showWarnings = F)

## set up titles
if (indicator == "dpt3_cov") plot_title <- "DPT3 Coverage"
if (indicator == "dpt1_cov") plot_title <- "DPT1 Coverage"
if (indicator == "dpt1_3_abs_dropout") plot_title <- "DPT1-3 Absolute Dropout"
if (indicator == "dpt1_3_rel_dropout") plot_title <- "DPT1-3 Relative Dropout"

if (indicator == "mcv1_cov") plot_title <- "MCV1 Coverage"

draws.df <- fread(sprintf("FILEPATH",
                          indicator_group, indicator, run_date))

# By region
pvtable.reg <- get_pv_table(d = draws.df,
                            indicator_group = indicator_group,
                            rd = run_date,
                            indicator=indicator,
                            result_agg_over = c("year", "oos", "region"),
                            coverage_probs = seq(from=5, to=95, by=10),
                            aggregate_on=c("country", "ad1", "ad2"),
                            draws = sum(grepl("draw[0-9]+",names(draws.df))),
                            out.dir = out_dir,
                            plot = TRUE,
                            plot_by = "region",
                            plot_by_title = "Region",
                            plot_ci = TRUE,
                            point_alpha = 0.5,
                            point_color = "black",
                            ci_color = "gray",
                            plot_title = plot_title)

# All regions
pvtable.all <- get_pv_table(d = draws.df,
                            indicator_group = indicator_group,
                            rd = run_date,
                            indicator=indicator,
                            result_agg_over = c("year", "oos"),
                            coverage_probs = seq(from=5, to=95, by=10),
                            aggregate_on=c("country", "ad1", "ad2"),
                            draws = sum(grepl("draw[0-9]+",names(draws.df))),
                            out.dir = out_dir,
                            plot = TRUE,
                            plot_ci = TRUE,
                            point_alpha = 0.1,
                            point_color = "black",
                            ci_color = "gray",
                            plot_title = plot_title)



Regions <- c("vax_eaas", "vax_name", "vax_trsa", "vax_ansa", "vax_caeu", "vax_crbn", "vax_cssa", "vax_ctam", "vax_essa", "vax_seas", "vax_soas", "vax_sssa", "vax_wssa")
