##################################
# POST UBCOV EXTRACTION DATA CLEANING FOR GEOSPATIAL DATA EXTRACTIONS & GEOGRAPHY MATCHING
# PIONEERED BY USERNAME
# UPDATED & OVERHAULED BY USERNAME
# EMAIL USERNAME
# EMAIL USERNAME

# INSTRUCTIONS: 
# UBCOV OUTPUTS MUST BE SAVED IN LIMITED USE DIRECTORY
###################################

if (as.character(Sys.info()[1]) == "Windows"){
  j <- "FILEPATH"
} else{
  j <- "FILEPATH"
}

folder_in <- paste0("FILEPATH") #where your extractions are stored
folder_out <- paste0("FILEPATH") #where you want to save the big csv of all your extractions together
topic <- "vaccination"
cores <- 30 #qlogin -pe multi_slot 30 -P proj_geospatial -now no
#source('FILEPATH')


##Get all geog codebooks and package them together
files <- list.files(paste0("FILEPATH"), pattern=".csv$", ignore.case = T, full.names = T)
files <- grep("Copy|together|linkage|IPUMS", files, value = T, invert = T) #IPUMS is handled separately

read_add_name_col <- function(file){
  rn <- gsub(".csv", "", file, ignore.case=T)
  spl <- strsplit(rn, "/") %>% unlist()
  svy <- spl[length(spl)]
  df <- fread(file)
  df$survey_series <- svy
  return(df)
}

geo <- lapply(files, read_add_name_col) %>% rbindlist(fill=T, use.names=T)

extractions <- list.files(folder_in, full.names=T)
extractions <- grep("IPUMS_CENSUS", extractions, invert=T, value = T) #IPUMS is handled separately

say_file_and_read_dta <- function(file){
  print(file)
  dta <- fread(file, encoding="Latin-1", data.table=FALSE)
  bye <- grep("map", names(dta), value=T)
  bye <- c(bye, "nid_n", "year_n", "end_year_n")
  dta <- dta[, !(names(dta) %in% bye)]
  return(dta)
}
topics <- mclapply(extractions, say_file_and_read_dta, mc.cores=cores) %>% rbindlist(fill=T, use.names=T)


#make everything the same data type to merge
geo[,nid <- as.character(nid)]
#geo$psu <- as.character(geo$psu)
topics[,nid := as.character(nid)]
topics$geospatial_id[is.na(topics$geospatial_id)] <- topics$psu[is.na(topics$geospatial_id)]
topics[,geospatial_id := as.character(geospatial_id)]

### Remove any characters after ISO3 code
#geo$iso3 <- as.character(geo$iso3)
#topics$ihme_loc_id <- as.character(topics$ihme_loc_id)
#geo$iso3 <- substr(geo$iso3, 0, 3)
#topics$ihme_loc_id <- substr(topics$ihme_loc_id, 0, 3)

all <- merge(geo, topics, by.x=c("nid", "geospatial_id", "iso3"), by.y=c("nid", "geospatial_id","ihme_loc_id"), all.x=F, all.y=T, allow.cartesian=TRUE)
all[,ihme_loc_id := iso3]

###SPECIFIC RECODES

match_on_admin1_string <- function(nid_vec){
  nid_vec <- as.character(nid_vec)
  for (nid_i in nid_vec){
    nid_geo <- subset(geo, nid == nid_i)
    nid_dta <- subset(topics, nid == nid_i)
    nid_geo[, "svy_area1"] <- as.character(nid_geo[, "svy_area1"])
    nid_dta[, "admin_1"] <- str_trim(nid_dta[, "admin_1"])
    nid_all <- merge(nid_geo, nid_dta, by.x = c("nid", "svy_area1"), by.y=c("nid", "admin_1"), all.x=F, all.y = T)
    all <<- all[all$nid != nid_i, ] #clear from all
    all <<- rbind.fill(all, nid_all) #rbind to all
  }
}

drop_x <- grep("\\.x$", names(all), value=T)


keep <- grep(paste0(drop_x, collapse="|"), names(all), value=T, invert=T)

all <- subset(all, select=keep)

fix <- all
fix[, num_clusters:= uniqueN(geospatial_id), by = c("nid", "ihme_loc_id", "year_start", "survey_name", "file_path.y")]

fix <- data.table(subset(fix, is.na(shapefile) & (is.na(lat) | is.na(long)) ))
fix_collapse <- fix[, uniqueN(geospatial_id), by = c("nid", "ihme_loc_id", "year_start", "survey_name", "file_path.y","num_clusters")]
names(fix_collapse)[names(fix_collapse) == "V1"] = "num_unmatched_clusters"

fix_total <- fix[, unique(geospatial_id), by=.(nid, ihme_loc_id, year_start, survey_name, file_path.y)]
names(fix_total)[names(fix_total) == "V1"] = "geospatial_id"


vaccine_data <- data.table(all)
rm(all)
module_date <- Sys.Date()
module_date <- gsub("-", "_", module_date)
  folder_out <- paste0(folder_out, "/", module_date)
  #SAVE .Rda OF EXTRACTIONS MATCHED TO GEOGRAPHIES
  saveRDS(vaccine_data, file=paste0(folder_out, "/", module_date, ".Rda"))
  ##WRITE .csv OF SURVEYS TO MATCH TO GEOGRAPHIES
  write.csv(fix_collapse, file=paste0(folder_out, "/", module_date, "_", topic, "_", "geography_matching_to_do", ".csv"), row.names = F)
  ##WRITE .csv of surveys + cluster numbers to match(gives an idea of how much is remaining)
  write.csv(fix_total, file=paste0(folder_out, "/", module_date, "_", topic, "_", "cluster_matching_to_do", ".csv"), row.names = F)
  ##WRITE .csv OF EXTRACTIONS MATCHED TO GEOGRAPHIES
  ##write.csv(vaccine_data, file=paste0(folder_out, "/", module_date, "_", topic, ".csv"), row.names = F)
