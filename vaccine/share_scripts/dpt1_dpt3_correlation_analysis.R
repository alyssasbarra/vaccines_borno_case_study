# Correlation analysis for DPT1 and DPT3 coverage

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
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0('FILEPATH'))
mbg_setup(package_list = package_list, repos = core_repo)

# Custom load indicator-specific functions
source(paste0('FILEPATH'))

dpt3_brick <- brick(paste0("FILEPATH"))
dpt1_brick <- brick(paste0("FILEPATH"))

compare_sample <- function(dpt3_brick, dpt1_brick, lyr) {
  message(paste0("  ", lyr))
  if (length(dpt3_brick[[lyr]]) != length(dpt1_brick[[lyr]])) stop("Unequal lengths")
  
  # Sample non-NA pixels only
  smp_idx <- sample(cellIdx(dpt3_brick[[lyr]]), 100000)
  dpt1_smp <- dpt1_brick[[lyr]][smp_idx]
  dpt3_smp <- dpt3_brick[[lyr]][smp_idx]
  
  return(data.table(dpt1 = dpt1_smp,
                    dpt3 = dpt3_smp,
                    layer = lyr))
}

# df_sample <- lapply(1:17, function(i) compare_sample(dpt3_brick = dpt3_brick, dpt1_brick = dpt1_brick, lyr = i))
# df_sample <- rbindlist(df_sample)
df_sample <- compare_sample(dpt3_brick, dpt1_brick, lyr = 17)

year_map <- data.table(layer = 1:17,
                       year = c(2000:2016))

df_sample <- merge(df_sample, year_map)

gg <- ggplot(df_sample[layer == 17], aes(x = dpt1, y = dpt3)) +
  geom_point(alpha= 0.05) +
  theme_bw() +
  labs(x = "DPT1 Coverage", 
       y = "DPT3 coverage", 
       title = "DPT1 coverage vs DPT3 coverage", 
       subtitle = "100,000 randomly selected pixels, 2016")

out_dir <- "FILEPATH"
dir.create(out_dir)

png(file = paste0(out_dir, "dpt1_dpt3_correlation.png"),
    height = 6,
    width = 10,
    units = "in", 
    res = 200)

print(gg)
dev.off()
