
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
package_list <- c(t(read.csv(sprintf('FILEPATH',commondir),header=FALSE)))

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

# Custom load indicator-specific functions
source(paste0('FILEPATH'))

# Options
indicator_group <- "vaccine"
indicator <- "dpt3_cov"
out_dir <- "FILEPATH"
dir.create(out_dir)

ras_old <- brick(paste0("FILEPATH"))
ras_new <- brick(paste0("FILEPATH"))

# Mask both of these to Africa only
template <- brick("FILEPATH")[[17]]

ras_old <- raster::crop(ras_old, template)
ras_new <- raster::crop(ras_new, template)

ras_old <- raster::mask(ras_old, template)
ras_new <- raster::mask(ras_new, template)

ras_diff <- ras_old - ras_new

pdf(file = "FILEPATH")

plot(ras_diff[[17]])
dev.off()

extract_vals <- function(ras, i) {
  val <- getValues(ras)
  dt_val <- data.table(vals = val[which(!is.na(val))],
                      layer = i)
  dt_val[, idx := .I]
  return(dt_val)
}
vals_all <- lapply(1:nlayers(ras_diff), function(i) extract_vals(ras_diff, i))
vals_all <- rbindlist(vals_all)

vals_all <- merge(vals_all, data.table(layer = 1:17, year = 2000:2016))

sample_vals_all <- sample(unique(vals_all$idx), 10000)
vals_all_plot <- vals_all[idx %in% sample_vals_all,]

gg <- ggplot(data = vals_all_plot, aes(x = vals)) +
        geom_density() +
        theme_bw() +
        xlim(-0.5, 0.5) +
        facet_wrap(~year) +
        labs(title = "Change in estimated DPT3 coverage with changes in country random effect and nugget effect priors",
             subtitle = "10,000 randomly selected pixels for each year",
             y = "Density",
             x = "Difference in estimated DPT3 coverage (old priors - new priors)")

png(file = paste0(out_dir, "change_in_dpt3_coverage_with_priors.png"), 
    height = 6,
    width = 10,
    res = 300,
    units = "in")
print(gg)
dev.off()    

summary(vals_all$vals)