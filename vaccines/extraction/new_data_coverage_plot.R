
#Function used in outliering process
parse_nid <- function(arow) {
  
  if (is.na(arow$year_id)) {
    
    # NAs: just return as NA    
    return(arow)
    
  } else if (grepl("-", arow$year_id)) {
    
    # dashes: range inclusive
    year_start <- as.numeric(stringr::str_match(arow$year_id, "(.*)-.*")[,2])
    year_end <- as.numeric(stringr::str_match(arow$year_id, ".*-(.*)")[,2])
    return(data.table(nid = arow$nid,
                      year_id = year_start:year_end))
    
  } else if (grepl("/", arow$year_id)) {
    
    # slashes: just the specified elements
    years <- as.numeric(unlist(strsplit(arow$year_id, "/")))
    return(data.table(nid = arow$nid,
                      year_id =years))
    
  } else {
    
    return(arow)
    
  }
}
###


list = as.numeric(gsub(".csv", "", list.files(paste0("FILEPATH"))))
coverage_data <- data.table()
for(nid in list){
  fn <- paste0("FILEPATH")
  
  if (file.exists(fn)) {
    
    df_pointpoly <- fread(fn)
    if (nrow(df_pointpoly) > 0) {
      
      # Specific function loads (avoid namespace conflicts)
      str_match <- stringr::str_match
      
      # Replace some functions to be quieter
      geom_polygon_quiet <- function(...) {suppressMessages(ggplot2::geom_polygon(...))}
      geom_path_quiet    <- function(...) {suppressMessages(ggplot2::geom_path(...))}
      
      log <- fread(file.path(extraction_root, "log/details", paste0(nid, ".csv")))
      
      #########################################
      ####Run prep from 03_plot_coverage#######
      #########################################
      #df_pointpoly <- df_pointpoly[me_name != "vacc_mcv2",] #don;t want mcv2
      
      # Prep data for graphing
      df_pointpoly[, latitude := as.numeric(latitude)]
      df_pointpoly[, longitude := as.numeric(longitude)]
      
      df_pointpoly$survey_name <- gsub("/", "_", df_pointpoly$survey_name)
      
      # For now, don't use age-specific "pseudoclusters"
      # First, group everyone into one year
      max_y <- max(df_pointpoly$year_id) + 1
      df_pointpoly[, svy_year := max_y]
      df_pointpoly[ svy_year >= 2000, ]
      coverage_data <- rbind(coverage_data, df_pointpoly, fill=T)
      
    }
  }  
}

#names(coverage_data)[names(coverage_data == "svy_id")] <- "nid"
#names(coverage_data)[names(coverage_data == "survey_name")] <- "source"
#names(coverage_data)[names(coverage_data == "year_id")] <- "year"

covs = c("dpt3_cov")

for(c in covs){
  
  df <- copy(coverage_data)
  vax_prefix <- str_match(vax, "^([a-z]*)[0-9]?_cov")[2]
  if(vax == "rotac_cov"){ vax_prefix <- "rota" }
  vax_dose <- str_match(vax, "^[a-z]*([0-9]?)_cov")[2] %>% as.numeric
  
  if(vax_prefix %in% c("bcg", "yfv")){ vax_dose <- ""  }
  if(vax_prefix == "rota"){ vax_dose = "c" }
  
  # Set Output directory
  #out_dir <- paste0(path, "/vacc_", gsub("_cov", "", vax), "/")
  #if (!dir.exists(out_dir)) dir.create(out_dir, showWarnings = F, recursive = T)
  
  ### Create vaccine variable
  # Figure out which doses are present in the data & grab max
  dose_vars <- unique(df$me_name)[grepl(vax_prefix, unique(df$me_name))]
  if(vax == "mcv1_cov" & length(dose_vars) > 0){ dose_vars <- "vacc_mcv1" } ##Ignore mcv2 for now
  if (length(dose_vars) > 0) {
    #max_doses <- str_match(dose_vars, paste0("vacc_", vax_prefix, "([0-9])"))[,2] %>%
    #  as.numeric %>% max
    
    #subset vaccine of interest
    df <- df[me_name %in% dose_vars]
    #reshape wide
    idvars <- c(names(df)[!names(df) %in% c("me_name", "value")])
    #test <- reshape(df, timevar = "me_name", idvar = idvars, direction="wide" )
    df <- data.table(data.table::dcast(df, svy_id + country + point + latitude + longitude + location_code +
                                         shapefile + survey_name + year_id + svy_year +  psu + weight ~ me_name, sum,value.var = c("value", "N")))
    
    val_vars <- paste0("value_", dose_vars)
    N_vars <- paste0("N_", dose_vars)
    names(df)[names(df) %in% val_vars] <- dose_vars
    names(df)[names(df) %in% N_vars] <- c(rep("drop", length(N_vars)-1), "N")
    if(is.na(vax_dose)) { cov_var <- paste0("vacc_",vax_prefix) } else { cov_var <-  paste0("vacc_",vax_prefix, vax_dose) }
    
    if (cov_var %in% names(df)) {
      
      df[, tmp := get(cov_var)] # select the dose of interest
      
      # Remove other vaccine dose variables
      df <- subset(df, select = !(names(df) %in% dose_vars))
      df <- subset(df, select = !(names(df) == "drop"))
      # Rename
      setnames(df, "tmp", vax)
      
      # Next, sum up over all variables but outcome & N
      df <- df[, list(N = sum(N), outcome = sum(get(vax))),
               by = setdiff(names(df), c("N", vax))]
      
      # Restore names
      setnames(df, "outcome", vax)
      
      # Finally, convert coverage to percents
      df[, eval(vax) := get(vax) / N]
      
      # Truncate long names
      df[survey_name == "GLOBAL_FUND_HOUSEHOLD_SURVEY", survey_name := "GLOBAL_FUND"]
      df[survey_name == "ARAB_LEAGUE_PAPFAM", survey_name := "ARAB_LG_PAPFAM"]
      df[survey_name == "MACRO_DHS_SP", survey_name := "MACRO_DHS"]
      df[survey_name == "WB_LSMS_ISA", survey_name := "WB_LSMS"]
      
      if("original_year" %in% names(df)) {
        setnames(df, "original_year", "year")
      }
      setnames(df, "svy_id", "nid")
      
      
      # Drop troublesome shapefiles 
      drop_shapefiles <- c("CRI_Pavas", 
                           "CRI_Tambor",
                           "CRI_Ulloa") 
      df <- subset(df, !(shapefile %in% drop_shapefiles))
      
      #########################################
      #########################################
      
      df[ ,source_label := "Other"]
      df[survey_name == "MACRO_DHS", source_label := "Demographic and Health Survey"]
      df[survey_name == "UNICEF_MICS", source_label := "UNICEF Multiple Indicator Cluster Survey"]
      df[survey_name == "WB_CWIQ", source_label := "World Bank Core Welfare Indicator Questionnaire Survey"]
      #df[is.na(source_label), source_label := survey_name]
      df[, source := source_label]
      
   
      ###Report and Literature Extractions
      
      df_lbd <- fread(paste0("FILEPATH", 
                               vax_prefix, "_custom_presample.csv"))
      setnames(df_lbd, "svy_id", "nid")
      df_lbd[ ,svy_year:=floor((year_start+year_end)/2)]
      df_lbd[ ,(vax) := get(paste0(vax_prefix, "_dose_", max_doses)) / N ]
      df_lbd[ ,point := 0]
      df_lbd[ ,source := "Custom"]
      df_lbd[ ,source_label := "Custom"]
      df <- data.table(rbind(df, df_lbd, fill=T))
      
      
      # A few more fixes for the coverage plotting code
      df[, cluster_id := .I]
      df[location_code == "" | location_code == "#N/A", location_code := NA]
      df[, location_code := as.numeric(location_code)]
      
      
      ###Apply Outliering
  
      outlier_names <- dose_vars
      
      outlier_fn <- "FILEPATH"
      outlier <- data.table(read.csv(outlier_fn))
      outlier <- outlier[lbd == 1 & (me_name == "" | me_name %in% outlier_names),] #keep only lbd outliers specific to mcv or all vaccines
      outlier$year_id <- as.character(outlier$year_id)
      outlier_yrs <- outlier[!is.na(as.numeric(year_id)),] #identify year ranges we need to split
      outlier_nid <- outlier[is.na(year_id) | year_id == "" | batch_outlier==1,]  #outlier the entire nid if batch = 1 or year_id is blank
      outlier_y2  <- outlier_yrs[year_id != "",] #ignore missing years
      outlier_yf <- lapply(1:nrow(outlier_y2), function(i)parse_nid(outlier_y2[i,])) %>% rbindlist #parse years 
      outlier_yrs <- rbind(outlier_yrs, outlier_yf, fill=T) #bind split years back into year level data
      outlier_ids <- paste(outlier_yrs$nid, outlier_yrs$year_id, sep = "_")
      df[ , outlier := 0]
      df[paste(nid, year_id, sep="_") %in% outlier_ids , outlier := 1]
      df[(as.character(nid)) %in% unique(as.character(outlier_nid$nid)) , outlier := 1]

    }
  }
}


names(df)[names(df) == "year_id"] <- "year"

reg = "sssa"

coverage_maps <- graph_data_coverage_values(df = df,
                                            var = vax,
                                            title = vax_title,
                                            year_min = 1997,
                                            year_max = 2016,
                                            year_var = 'svy_year',
                                            region = reg,
                                            cores = cores,
                                            indicator = vax,
                                            since_date = 'DATE',
                                            high_is_bad = FALSE,
                                            return_maps = TRUE,
                                            legend_title = "Vaccine \nCoverage",
                                            endemic_gauls = NULL,
                                            map_point_size = 0.8,
                                            fast_shapefiles = T,
                                            simplify_polys = T,
                                            tolerance = 0.01,
                                            save_on_share = F,
                                            base_font_size = 18,
                                            core_repo = core_repo,
                                            #out_dir = paste0("FILEPATH"),
                                            out_dir = "FILEPATH",
                                            prep_shiny = F,
                                            color_scheme = "classic")

