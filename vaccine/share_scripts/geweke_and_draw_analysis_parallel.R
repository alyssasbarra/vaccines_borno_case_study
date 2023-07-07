###############################################################################
## SETUP
###############################################################################

## clear environment
rm(list=ls())

## Set repo location and indicator group
user               <- 'USERNAME'
core_repo          <- sprintf('FILEPATH')
indic_repo         <- sprintf('FILEPATH')
remote             <- 'origin'
branch             <- 'develop'
pullgit            <- FALSE

## sort some directory stuff
commondir      <- sprintf('FILEPATH')
package_list <- c(t(read.csv(sprintf('FILEPATH'),header=FALSE)))

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0('FILEPATH'))
mbg_setup(package_list = package_list, repos = core_repo)

library(coda)

# Custom load indicator-specific functions
source(paste0('FILEPATH'))

load_from_parallelize() 
# varying: ind, reg
# static: rd, out_dir

# Options
age <- 0
holdout <- 0

message(paste0("Run date: ", run_date))
message(paste0("Output directory: ", out_dir))

# Names
if (ind == "dpt3_cov") ind_title <- "DPT3 Coverage"
if (ind == "dpt1_3_abs_dropout") ind_title <- "DPT1-3 Absolute Dropout"
if (ind == "dpt1_cov") ind_title <- "DPT1 Coverage"

if (reg == "cssa") reg_title <- "Central Sub-Saharan Africa"
if (reg == "essa") reg_title <- "Eastern Sub-Saharan Africa"
if (reg == "name") reg_title <- "North Africa"
if (reg == "sssa") reg_title <- "Southern Sub-Saharan Africa"
if (reg == "wssa") reg_title <- "Western Sub-Saharan Africa"

message(paste0("\nIndicator: ", ind, " | Region: ", reg))

# Load a cell pred
message("  -- Loading cell pred object...")
main_dir <- paste0('FILEPATH')
load(paste0(main_dir, ind, '_cell_draws_eb_bin0_', reg, '_0.RData'))

draw_cols <- names(cell_pred)[grep("V[0-9]*", names(cell_pred))]

# Define a function to calculte the geweke statistic
geweke_vec <- function(x) {
  if (anyNA(x)) return(NA)
  output <- geweke.diag(mcmc(x))[[1]]
  return(as.numeric(output))
}

# Subset cell pred to 10,000 observations
cell_pred <- cell_pred[sample(nrow(cell_pred), 10000),]

# Calculate geweke values
message("  -- Calculating Geweke statistics...")
geweke_vals <- apply(cell_pred, 1, geweke_vec)  
geweke_vals <- data.table(geweke = geweke_vals,
                          indicator = ind,
                          region = reg)

message("  -- Plotting...")

# Plot distribution of geweke values
gg <- ggplot(data = geweke_vals, aes(x = geweke)) +
  stat_function(aes(color = "red"), fun = dnorm, n = nrow(cell_pred), args = list(mean=0, sd=1), alpha = 0.8) +
  stat_density(aes(color = "blue"), geom = "line", alpha = 0.8) +
  theme_classic() +
  labs(x = "Geweke statistic",
       y = "Density",
       title = paste0(ind_title, " | Region: ", reg_title),
       subtitle = "Distribution of Geweke statistics (for 10,000 randomly sampled pixels)") +
  scale_color_manual(name = NULL, values = c('red' = 'red', 'blue' = 'blue'), labels = c("Geweke statistic distribution", "Standard normal distribution"))

png(file=paste0(out_dir, "geweke_", ind, "_", reg, ".png"), 
    height = 4, 
    width = 8,
    res = 200, 
    units = "in")

print(gg)
dev.off()

# Plot 50 random rows (MCMC-style)
cell_pred2 <- data.table(cell_pred[sample(nrow(cell_pred), 50),])
cell_pred2[, idx := .I]
cell_pred2 = melt(cell_pred2, id.vars = "idx")
cell_pred2[, sample := as.numeric(stringr::str_sub(variable, 2))]

gg <- ggplot(data = cell_pred2,
             aes(x=sample, y=value)) +
  geom_line() +
  facet_wrap(~idx) +
  theme_bw() +
  labs(title = paste0(ind_title, " | Region: ", reg_title), 
       subtitle = "50 randomly-selected pixels",
       x = "Sample", y = "Value")


png(file=paste0(out_dir, "draws_", ind, "_", reg, ".png"), 
    height = 8, 
    width = 12,
    res = 200, 
    units = "in")

print(gg)
dev.off()


saveRDS(geweke_vals, file = paste0(out_dir, ind, "_", reg, "_", "geweke_vals.rds"))