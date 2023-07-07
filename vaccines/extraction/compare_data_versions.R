#----HEADER-------------------------------------------------------------------------------------------------------------
# Author:  USERNAME
# Date:    DATE
# Purpose: Compare prepared GBD input data (vaccination.rds files) to ensure no data drops, etc. after a new data lock
# Run: Interactive: file.path("FILEPATH")
#***********************************************************************************************************************


#----SET UP-------------------------------------------------------------------------------------------------------------
# clean work space
rm(list=ls())

# set os flexibility
os <- .Platform$OS.type
if (os == "windows") {
  j_root <- "FILEPATH"
} else {
  j_root <- "FILEPATH"
  username <- "USERNAME"
}

# load packages
library(data.table)  
library(magrittr)

# source paths and functions
source(paste0("FILEPATH"))  
source("FILEPATH")  
source("FILEPATH")
source("FILEPATH")
source("FILEPATH")
'%!in%'  <- function(x,y)!('%in%'(x,y))

### key comparison objects
# versions to compare
date_new_ <- "DATE"
cycle_new_ <- "gbd2020"
date_old_ <- "DATE"  # make sure these are the dates pointing to the to_model_dir, not "-covidfree" post-processing
cycle_old_  <- "gbd2020"
# columns to compare
cols_to_compare <- c("nid", "survey_name", "ihme_loc_id", "year_id", "age_year", "me_name", "data",
                     "cv_lit", "cv_admin", "cv_survey", "cv_outlier", "cv_admin_orig")
exact_mes <- c("vacc_bcg", "vacc_dpt1", "vacc_dpt3", "vacc_hepb3", "vacc_hib3", "vacc_mcv1", "vacc_mcv2",
               "vacc_pcv3", "vacc_polio3", "vacc_rcv1", "vacc_rotac")
mes_to_compare <- c(exact_mes, "vacc_dpt12_cond", "vacc_dpt3_timeliness_ratio")

### output directory 
to_model_dir_data <- file.path(data_root, "exp/to_model", cycle_new_, date_new_, "input_data_comparison")
ifelse(!dir.exists(to_model_dir_data), dir.create(to_model_dir_data), FALSE)

### write an output file of the two data date-versions being compared
cat(paste0(date_new_, " vs ", date_old_, " 'vaccination.rds' input data comparison"), file=file.path("FILEPATH"))
#***********************************************************************************************************************


#----FUNCTIONS TO CREATE DATA COMPARISON DF OBJECTS---------------------------------------------------------------------
# function to remove different cohorts from "data_old_only" object (e.g. if Decider decisions changed)
remove_non_dropped <- function(me, dt, dt_new) {
  
  dt_me <- dt[me_name==me]
  dt_new_me <- dt_new[me_name==me]
  dt_me[nid %in% unique(dt_new_me$nid), remove := 1]
  return(dt_me)
  
}

# function to create data comparison dfs
create_data_comparisons <- function(date_new=date_new_, cycle_new=cycle_new_, 
                                    date_old=date_old_, cycle_old=cycle_old_, 
                                    cols=cols_to_compare, 
                                    mes=mes_to_compare) {
  
  # read in and subset data
  message("Reading in and subsetting the new and old data versions")
  data_new <- readRDS(file.path(data_root, "exp/to_model", cycle_new, date_new, "vaccination.rds")) %>% 
    .[, ..cols] %>% .[me_name %in% mes]
  data_old <- readRDS(file.path(data_root, "exp/to_model", cycle_old, date_old, "vaccination.rds")) %>% 
    .[, ..cols] %>% .[me_name %in% mes]
  
  # add identifying column to both dfs
  data_new[, new := 1]
  data_old[, old := 1]
  setnames(data_new, c("data", "cv_outlier"), c("data_new", "cv_outlier_new"))
  setnames(data_old, c("data", "cv_outlier"), c("data_old", "cv_outlier_old"))
  
  # merge 
  message("Merging data versions together")
  data_cp <- merge(data_new, data_old, by=names(data_new)[!names(data_new) %in% c("data_new", "cv_outlier_new", "new")], all.x=T, all.y=T)
  
  # simple clean
  message("Making the comparisons")
  data_cp <- data_cp[!(is.na(data_new) & is.na(cv_outlier_new) & is.na(data_old) & is.na(cv_outlier_old))]
  
  ### start making comparison objects to return
  data_cp_all <- copy(data_cp)
  
  # subset out the data that's the same
  data_cp_same <- data_cp[new==1 & old==1]  
  
  # subset out the data that's only in one or the other, by location-nid-antigen
  data_new_only <- data_cp[new==1 & is.na(old)]  # includes new/alternative cohorts
  data_old_only <- data_cp[is.na(new) & old==1]
  
  # filter out alternative cohorts by me_name for data_old_only (e.g., don't highlight Decider diffs here)
  data_old_only <- lapply(mes[mes !="vacc_dpt3_timeliness_ratio"], function(x) {
    remove_non_dropped(me=x, dt=data_old_only, dt_new=data_new)
  }) %>% rbindlist()
  
  data_old_only <- data_old_only[is.na(remove)]
  
  # return a list object of the various resultant comparison dfs
  data <- list(data_cp_all,  
               data_cp_same, 
               data_new_only, 
               data_old_only)
  
  message("Done")
  return(data)
  
  
}
#***********************************************************************************************************************


#----CALL FUNCTIONS-----------------------------------------------------------------------------------------------------
#' Returned list object contains the following dfs:
#' [[1]] : the full merged comparison table -- inclusive of all new + old data 
#' [[2]] : all of the data that's the same between new and old
#' [[3]] : the data in the new file only (inclusive of new, alternative cohorts for me-NIDs that might have been used previously)
#' [[4]] : the data in the old file only, supposedly "dropped", needing more review

data_list <- create_data_comparisons()
#***********************************************************************************************************************


#----SAVE DATA OBJECTS--------------------------------------------------------------------------------------------------
saveRDS(data_list[[1]], file.path(to_model_dir_data, "1_all_new_and_old_data.rds"))
saveRDS(data_list[[2]], file.path(to_model_dir_data, "2_same_new_and_old_data.rds"))
saveRDS(data_list[[3]], file.path(to_model_dir_data, "3_new_data_or_cohorts.rds"))
saveRDS(data_list[[4]], file.path(to_model_dir_data, "4_old_or_dropped_data.rds"))
#***********************************************************************************************************************

# Done - time to custom review "4_old_or_dropped_data.csv"! (Use outliering_prep.R to help subset data vaccination.rds objects)

#----PLOT---------------------------------------------------------------------------------------------------------------
library(ggplot)
library(ggrepel)

# Get all data
data <- data_list[[1]]

# Drop subnationals and conditional antigens
data <- data[!grepl("_", ihme_loc_id) & !grepl("cond|ratio", me_name)]

# Hotfixes
data <- data[nid != 12389, ]# Remove Syria nid 12389, which is mistakenly being compared against nid 10023 (newly outliered)

# Subset to valid data present in both versions
data <- data[new == 1 & old == 1 & is.na(cv_outlier_old) & is.na(cv_outlier_new), ]

# Get difference
data[, diff := round(data_new - data_old, digits = 2)]

# Make useful labels: for each nid, only show label once
threshold <- .04
data[, row_id := 1:nrow(data)]

row_ids_to_label <- data[abs(diff) > threshold, .SD[abs(diff) == max(abs(diff)), .(row_ids_to_label = row_id)], by = c("nid", "ihme_loc_id")]$row_ids_to_label
row_ids_to_label <- data[row_id %in% row_ids_to_label, .SD[1, .(row_ids_to_label = row_id)], by = c("nid", "ihme_loc_id")]$row_ids_to_label

# row_ids_to_label           <- c(admin_row_ids_to_label, non_admin_row_ids_to_label)
data[row_id %in% row_ids_to_label, lab := paste0(nid, "_", ihme_loc_id)]

# For each nid/iso3 combination > threshold, color code uniquely
data[abs(diff) > threshold, color_identifier := paste0(nid, "_", ihme_loc_id)]

# Plot all antigens
gg_all <- ggplot(data) +
  geom_abline() +
  geom_point(mapping = aes(x = data_old, y = data_new, color = color_identifier), 
             alpha = .5) +
  geom_text_repel(data = data[!is.na(lab), ], 
                  mapping = aes(x = data_old, y = data_new,  label = lab),
                  box.padding = 0.1, max.overlaps = Inf) +
  xlab("DATE") +
  ylab("DATE") +
  ggtitle("Data Comparison: All Antigens") +
  theme_bw()

# Facet antigen-specifically
gg_antigen <- ggplot(data) +
  geom_abline() +
  geom_point(mapping = aes(x = data_old, y = data_new, color = color_identifier), 
             alpha = .5) +
  geom_text_repel(data = data[!is.na(lab), ], 
                  mapping = aes(x = data_old, y = data_new,  label = lab),
                  box.padding = 0.1, max.overlaps = Inf) +
  facet_wrap(~ me_name) +
  xlab("DATE") +
  ylab("DATE") +
  ggtitle("Data Comparison: All Antigens") +
  theme_bw()




