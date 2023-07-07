# LOAD PACKAGES ---------------------------------------------------------------
library(data.table)
library(reshape2)
library(ggplot2)


# HELPER FUNCTIONS ------------------------------------------------------------

# Given a filepath of tabulated data, read and rbind all data together
load_data <- function(path) {
  files <- list.files(path)
  files <- files[grepl(pattern = ".csv", files)]
  
  loaded_data <- data.table()
  
  for (file in files) {
    print(file)
    data <- fread(paste0(path, file))
    loaded_data <- rbind(loaded_data, 
                         data,
                         fill = TRUE)
  }
  loaded_data
}


collapse_data <- function(data) {
  collapsed_data <- unique(data[, .(survey_name, country, value = sum(value), N = sum(N)), by = .(me_name, year_id, svy_id)])
  collapsed_data[, value := value / N]
  collapsed_data$N <- NULL
  collapsed_data
}


# SET VARIABLES ---------------------------------------------------------------

new_data_path <- "FILEPATH"
old_data_path <- "FILEPATH"

outdir   <- "FILEPATH"

new_data <- load_data(new_data_path)
old_data <- load_data(old_data_path)

new_data <- collapse_data(new_data)
old_data <- collapse_data(old_data)



# PREPARE FOR PLOTTING --------------------------------------------------------


# merge gbd2019 and tabs together
data <- merge(old_data, new_data, by=c("svy_id", "survey_name", "country", "me_name", "year_id"), all.x = T, all.y = T)
setnames(data, "value.x", "new")
setnames(data, "value.y", "old")
# both <- both[complete.cases(both)]
# both <- both[, tabulation := tabulation/100]

# add a 'difference' column to 'both' dataset
data <- data[, difference := new - old]
data <- data[, difference := abs(difference)]



# PLOT ------------------------------------------------------------------------

antigens_to_model <- c("MCV" = "vacc_mcv1", "DPT" = "vacc_dpt3", "Polio3" = "vacc_polio3", "BCG" = "vacc_bcg")

for (i in 1:length(antigens_to_model)){
  
  plot <- ggplot(data = data[me_name == antigens_to_model[i]], aes(x=old, y=new)) +
    geom_point() +
    theme_bw() +
    # geom_smooth(method=lm, se=FALSE) +
    labs(title = paste0(names(antigens_to_model)[i], " - LBD Tabulated: New vs Prior"),
         x = "Prior LBD Tabulated Estimates",
         y = paste0("New LBD Tabulated Estimates")) +
    geom_abline(slope = 1, intercept = 0) +
    coord_equal() +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    geom_text(aes(label=ifelse(abs(difference) > 0.01, as.character(paste(nid, "|", me_name, "|", year_id)), "")),
              hjust=0, vjust=0, size = 3)
  
  ggsave(plot, file=paste0(outdir, names(antigens_to_model)[i], "_LBD_tabulate_comparison.png"), width = 12.6, height = 6)
  
}














