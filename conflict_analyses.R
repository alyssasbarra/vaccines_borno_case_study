
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
setwd(core_repo)
commondir      <- sprintf('FILEPATH')
package_list <- c(t(read.csv(sprintf('FILEPATH'),header=FALSE)))

library(data.table)
library(raster)
library(sf)

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0('FILEPATH'))
mbg_setup(package_list = package_list, repos = core_repo)

# Custom load indicator-specific functions
source(paste0('FILEPATH'))
source('FILEPATH')
'%!in%' <- function(x,y)!('%in%'(x,y))

source(paste0('FILEPATH'))
source('FILEPATH')

##### set up
run_date <- 'RUN_DATE'
adm_nga <- 164

makefigs <- T

############################################################################################################
################## Load shapefiles to use later

ad0_shp <- prep_ad_shp_for_mapping(ad_level = 0,
                                   shapefile_version = "DATE", 
                                   adm0_list = get_adm0_codes("NGA"))

ad1_shp <- prep_ad_shp_for_mapping(ad_level = 1,
                                   shapefile_version = "DATE", 
                                   adm0_list = get_adm0_codes("NGA"))

borno_shp <- subset(ad1_shp, ADM1_NAME == "Borno")

ad2_shp <- prep_ad_shp_for_mapping(ad_level = 2,
                                   shapefile_version = "DATE", 
                                   adm0_list = get_adm0_codes("NGA"))

borno_ad2_shp <- subset(ad2_shp, ADM1_NAME == "Borno")

############################################################################################################
################## pull GBD for introduction

gbd <-readRDS("FILEPATH")
gbd <- subset(gbd, location_id == 214 & year_id %in% c(2000, 2019))

############################################################################################################
################## conflict data

#### load in conflict data, mean raster, and shapefile
conflict <- fread('FILEPATH')

conflict2016 <- subset(conflict, year == 2016)
conflict2017 <- subset(conflict, year == 2017)
conflict2018 <- subset(conflict, year == 2018)

borno16 <- subset(conflict2016, admin1 == "Borno")
borno17 <- subset(conflict2017, admin1 == "Borno")
borno18 <- subset(conflict2018, admin1 == "Borno")

# conflict-related incidents reported in borno state 
dim(borno16)[1] + dim(borno17)[1] + dim(borno18)[1] 

# conflict-related incidents reported in borno state 
sum(borno16$fatalities) + sum(borno17$fatalities) + sum(borno18$fatalities) 

# prop of conflict in NGA that are in borno state alone
(dim(borno16)[1] + dim(borno17)[1] + dim(borno18)[1]) / (dim(conflict2016)[1] + dim(conflict2017)[1] + dim(conflict2018)[1] )

#### figure 1

if(makefigs){
  
  all <- conflict2018 %>% sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
  
  only_conflict_nga <- ggplot() + 
    geom_sf(aes(shape = event_type, size=fatalities, color=event_type), alpha=0.75, data = all)+ 
    geom_sf(data = ad1_shp, color = "black", lwd = 0.3, fill = NA) +  theme_empty() + 
    scale_colour_manual(values = c("#fb8072", "#fdb462", "#ffffb3", "#8dd3c7", "#80b1d3", "#bebada")) + 
    coord_sf(datum = NA) + labs(color = "Event type", size=NULL, shape=NULL) + scale_shape(guide = 'none') + 
    scale_size(guide = 'none')   + theme(legend.position = "none")
  
  
  png(file = paste0("FILEPATH"),
      width = 10,
      height = 10,
      units = "in",
      res = 400,
      type = "cairo")
  print(only_conflict_nga)
  dev.off()
  
  borno_conflict2018 <- subset(conflict2018, admin1 == "Borno")
  borno_conflict2018 <- borno_conflict2018 %>% sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
  
  only_conflict_b <- ggplot() + 
    geom_sf(data = borno_ad2_shp,color = "lightgrey",lwd = 0.3,fill = NA) + 
    geom_sf(aes(shape = event_type, size=fatalities, color=event_type), alpha=0.75, data = borno_conflict2018)+ 
    geom_sf(data = borno_shp,color = "black",lwd = 0.3,fill = NA) +  theme_empty() + 
    scale_colour_manual(values = c("#fb8072", "#fdb462", "#ffffb3", "#8dd3c7", "#80b1d3", "#bebada")) + 
    coord_sf(datum = NA) + labs(color = "Event type", size=NULL, shape=NULL) + scale_shape(guide = 'none') + 
    scale_size(guide = 'none')  + theme(legend.position = "none")
  
  
  png(file = paste0("FILEPATH"),
      width = 5,
      height = 10,
      units = "in",
      res = 400,
      type = "cairo")
  print(only_conflict_b)
  dev.off()
}

############################################################################################################
################## survey only
############################################################################################################

library(sf)
dhs <- fread('FILEPATH')
head(dhs)
dhs <- subset(dhs, me_name == "vacc_dpt3")
dhs <- subset(dhs, !is.na(latitude))
dhs <- dhs %>% sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

test_dhs <- st_intersection(dhs, borno_ad2_shp)
length(unique(test_dhs$NAME_2)) # number of districts in borno where DHS sampled

mics <- fread('FILEPATH')
head(mics)
#mics <- subset(mics, year_id == 2015)
mics <- subset(mics, me_name == "vacc_dpt3")
mics <- subset(mics, !is.na(latitude))
mics <- mics %>% sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

test_mics <- st_intersection(mics, borno_ad2_shp)
length(unique(test_mics$NAME_2)) # number of districts in borno where DHS sampled

############################################################################################################
################## conflict + survey
############################################################################################################

### first MCIS
df_mics <- copy(test_mics)
data_codes_mics <- unique(df_mics$ADM2_CODE)
data_names_mics <- unique(df_mics$ADM2_NAME)

conflict_space_2017 <- conflict2017 %>% sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
test_conflict <- st_intersection(conflict_space_2017, borno_ad2_shp)
conflict_codes <- unique(test_conflict$ADM2_CODE)

mics_diff <- setdiff(data_codes_mics, conflict_codes)
length(data_codes_mics) - length(mics_diff) # of districts where there was conflict when sampling for MICS


### now DHS
df_dhs <- copy(test_dhs)
data_codes_dhs <- unique(df_dhs$ADM2_CODE)
data_names_dhs <- unique(df_dhs$ADM2_NAME)

conflict_space_2018 <- conflict2018 %>% sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
test_conflict <- st_intersection(conflict_space_2018, borno_ad2_shp)
conflict_codes <- unique(test_conflict$ADM2_CODE)

dhs_diff <- setdiff(data_codes_dhs, conflict_codes)
length(data_codes_dhs) - length(dhs_diff) # of districts where there was conflict when sampling for MICS

if(makefigs){
  
  only_conflict_borno_same_color_2017 <- ggplot() + 
    geom_sf(aes(), color="#9999CC", data = borno_conflict2017)+ 
    geom_sf(data = borno_shp,color = "black",lwd = 0.3,fill = NA) +  
    theme_empty() + coord_sf(datum = NA) + theme(legend.background = element_rect())
  
  borno_mics_points <- only_conflict_borno_same_color_2017 + 
    geom_sf(data = borno_ad2_shp,color = "grey",lwd = 0.75,fill = NA) +
    geom_sf(data = borno_shp,color = "black",lwd = 1,fill = NA) +  
    geom_sf(aes(size=N), color="#66CC99", shape=17, data = test_mics) + theme_empty() + theme(legend.position = "none") + 
    coord_sf(datum = NA) 
  
  png(file = paste0("FILEPATH"),
      width = 5,
      height = 10,
      units = "in",
      res = 400,
      type = "cairo")
  print(borno_mics_points)
  dev.off()
  
  only_conflict_borno_same_color_2018 <- ggplot() + 
    geom_sf(aes(), color="#9999CC", data = borno_conflict2018)+ 
    geom_sf(data = borno_shp,color = "black",lwd = 0.3,fill = NA) +  
    theme_empty() + coord_sf(datum = NA) + theme(legend.background = element_rect())
  
  borno_dhs_points <- only_conflict_borno_same_color_2018 + 
    geom_sf(data = borno_ad2_shp,color = "grey",lwd = 0.75,fill = NA) +
    geom_sf(data = borno_shp,color = "black",lwd = 1,fill = NA) +  
    geom_sf(aes(size=N), color="#66CC99", shape=17, data = test_dhs) + theme_empty() + theme(legend.position = "none") + 
    coord_sf(datum = NA) 
  
  png(file = paste0("FILEPATH"),
      width = 5,
      height = 10,
      units = "in",
      res = 400,
      type = "cairo")
  print(borno_dhs_points)
  dev.off()
}

############################################################################################################
################## all post -2010 surveys

## for the OTHER dhs
dhs_77390 <- fread('FILEPATH')
head(dhs_77390)
dhs_77390 <- subset(dhs_77390, me_name == "vacc_dpt3")
dhs_77390 <- subset(dhs_77390, !is.na(latitude))
dhs_77390 <- dhs_77390 %>% sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

keep_nga_dhs_77390 <- st_intersection(dhs_77390, ad2_shp)
keep_nga_dhs_77390 <- unique(keep_nga_dhs_77390$NAME_2)

test_dhs_77390 <- st_intersection(dhs_77390, borno_ad2_shp)
keep_dhs1_dhs_77390 <- unique(test_dhs_77390$NAME_2)

## for the dhs
keep_nga_dhs <- st_intersection(dhs, ad2_shp)
keep_nga2_dhs <- unique(keep_nga_dhs$NAME_2)
keep_dhs2_dhs <- unique(test_dhs$NAME_2)

## for the mics
keep_nga_mics <- st_intersection(mics, ad2_shp)
keep_nga2_mics <- unique(keep_nga_mics$NAME_2)
keep_dhs2_mics <- unique(test_mics$NAME_2)

## # of lgas that havent been sampled since 2010 in borno
length(unique(borno_ad2_shp$NAME_2)) - length(unique(c(keep_dhs1_dhs_77390, keep_dhs2_dhs, keep_dhs2_mics)))
length(unique(borno_ad2_shp$NAME_2))
(length(unique(borno_ad2_shp$NAME_2)) - length(unique(c(keep_dhs1_dhs_77390, keep_dhs2_dhs, keep_dhs2_mics)))) / length(unique(borno_ad2_shp$NAME_2)) 


## # of lgas that havent been sampled since 2010 in nigeria (there are 774 lgas in nigeria)
774 - length(unique(c(keep_nga_dhs_77390, keep_nga2_dhs, keep_nga2_mics)))
(774 - length(unique(c(keep_nga_dhs_77390, keep_nga2_dhs, keep_nga2_mics))) ) / 774

### test to make sure there aren't other nids!
input_data_file <- fread(paste0('FILEPATH'))
nga_nids <- subset(input_data_file, country == "NGA")
nid_year_combos <- paste0(nga_nids$svy_id, "_", nga_nids$year_end)
unique(nid_year_combos) 

############################################################################################################
################## underlying population represented

file_path <- paste0('FILEPATH')

load(paste0('FILEPATH'))
rm(admin_1)
rm(admin_0)
head(admin_2)

admin_2_codes <- fread(paste0('FILEPATH'))
admin_2_codes <- data.table(admin_2_codes$ADM0_NAME, admin_2_codes$ADM0_CODE, admin_2_codes$ADM1_NAME, admin_2_codes$ADM1_CODE, admin_2_codes$ADM2_NAME, admin_2_codes$ADM2_CODE)
colnames(admin_2_codes) <- c("ADM0_NAME", "ADM0_CODE", "ADM1_NAME", "ADM1_CODE", "ADM2_NAME", "ADM2_CODE")
admin_2_codes <- unique(admin_2_codes)

df <- merge(admin_2, admin_2_codes, by="ADM2_CODE")
df <- subset(df, ADM0_CODE == adm_nga)
df <- subset(df, !is.na(ADM2_CODE))

df <- subset(df, ADM1_NAME == "Borno")
df16 <- subset(df, year == 2016)
tot16 <- sum(df16$pop)
df17 <- subset(df, year == 2017)
tot17 <- sum(df17$pop)
df18 <- subset(df, year == 2018)
tot18 <- sum(df17$pop)

mics_pop <- subset(df17, ADM2_CODE %in% data_codes_mics)
sum(mics_pop$pop) / tot17

dhs_pop <- subset(df18, ADM2_CODE %in% data_codes_dhs)
sum(dhs_pop$pop) / tot18


############# underlying prop from grid sampled
grid <- fread('FILEPATH')
grid <- subset(grid, state == "Borno")
grid$local <- ifelse(grid$local=="Askira/Uba", "Askira/U", grid$local)
grid$local <- ifelse(grid$local=="Maiduguri", "Maidugur", grid$local)
grid <- data.table(cbind(grid$local, grid$mean))
colnames(grid) <- c('ADM2_NAME', 'grid')

grid$grid <- as.numeric(as.character(grid$grid))
nga_tot_grid <- sum(grid$grid)


mics_pop_grid <- subset(grid, ADM2_NAME %in% data_names_mics)
sum(mics_pop_grid$grid) / nga_tot_grid


dhs_pop_grid <- subset(grid, ADM2_NAME %in% data_names_dhs)
sum(dhs_pop_grid$grid) / nga_tot_grid 

############################################################################################################
################## Coverage changes / calculations

file_path <- paste0('FILEPATH')

load(paste0('FILEPATH'))
rm(admin_0)
rm(admin_1)

################## district-level increases with uncertainty (2018 - 2000)
admin_2_codes <- fread(paste0('FILEPATH'))
admin_2_codes <- data.table(admin_2_codes$ADM0_NAME, admin_2_codes$ADM0_CODE, admin_2_codes$ADM1_NAME, admin_2_codes$ADM1_CODE, admin_2_codes$ADM2_NAME, admin_2_codes$ADM2_CODE)
colnames(admin_2_codes) <- c("ADM0_NAME", "ADM0_CODE", "ADM1_NAME", "ADM1_CODE", "ADM2_NAME", "ADM2_CODE")
admin_2_codes <- unique(admin_2_codes)

df <- merge(admin_2, admin_2_codes, by="ADM2_CODE")
df <- subset(df, ADM0_CODE == adm_nga)
df <- subset(df, !is.na(ADM2_CODE))


draws0 <- df[which(df$year==2000),] 
draws18 <- df[which(df$year==2019),] 

draws0 <- draws0[,-c("pop", "ADM2_CODE", "year", "ADM0_NAME", "ADM0_CODE", "ADM1_NAME", "ADM1_CODE", "ADM2_NAME")]
draws18 <- draws18[,-c("pop", "ADM2_CODE", "year", "ADM0_NAME", "ADM0_CODE", "ADM1_NAME", "ADM1_CODE", "ADM2_NAME")]


all <- draws18 - draws0 
all.t <- ifelse(all > 0, 1, 0)

vect.t <- colSums (all.t, na.rm = TRUE, dims = 1)
vect.t <- vect.t / length(which(!is.na(all.t[,1])))

mean(vect.t)
quantile(vect.t, probs=c(0.025, 0.975))

##### borno only
admin_1_means <- fread(paste0('FILEPATH'))
admin_1_means <- subset(admin_1_means,  ADM1_NAME == "Borno" & year %in% c(2000, 2016, 2017, 2018, 2019))
admin_1_means

file_path_dpt3 <- paste0('FILEPATH')
admin_1_means <- fread(paste0('FILEPATH'))
admin_1_means <- subset(admin_1_means,  ADM1_NAME == "Borno" & year %in% c(2000, 2016, 2017, 2018, 2019))
admin_1_means

df <- subset(df, ADM1_NAME == "Borno")
df <- subset(df, !is.na(ADM2_CODE))

draws0 <- df[which(df$year==2000),] 
draws18 <- df[which(df$year==2019),] 

draws0 <- draws0[,-c("pop", "ADM2_CODE", "year", "ADM0_NAME", "ADM0_CODE", "ADM1_NAME", "ADM1_CODE", "ADM2_NAME")]
draws18 <- draws18[,-c("pop", "ADM2_CODE", "year", "ADM0_NAME", "ADM0_CODE", "ADM1_NAME", "ADM1_CODE", "ADM2_NAME")]

all <- draws18 - draws0 
all.t <- ifelse(all > 0, 1, 0)

vect.t <- colSums (all.t, na.rm = TRUE, dims = 1)
vect.t <- vect.t / length(which(!is.na(all.t[,1])))

mean(vect.t)
quantile(vect.t, probs=c(0.025, 0.975))

############################################################################################################
################## GRID re-aggregated

grid <- fread('FILEPATH')
grid <- subset(grid, state == "Borno")
grid$local <- ifelse(grid$local=="Askira/Uba", "Askira/U", grid$local)
grid$local <- ifelse(grid$local=="Maiduguri", "Maidugur", grid$local)
grid <- data.table(cbind(grid$local, grid$mean))
colnames(grid) <- c('ADM2_NAME', 'grid')


file_path <- paste0('FILEPATH')
admin_2_codes <- fread(paste0('FILEPATH'))

ad2 <- subset(admin_2_codes, ADM1_NAME == "Borno")
ad2 <- subset(ad2, year == 2019)

test <- merge(ad2, grid, by="ADM2_NAME")
test$mean <- as.numeric(as.character(test$mean))
test$grid <- as.numeric(as.character(test$grid))

weighted.mean(test$mean, test$grid) 

file_path <- paste0('FILEPATH')
admin_2_codes <- fread(paste0('FILEPATH'))

ad2 <- subset(admin_2_codes, ADM1_NAME == "Borno")
ad2 <- subset(ad2, year == 2019)

test <- merge(ad2, grid, by="ADM2_NAME")
test$mean <- as.numeric(as.character(test$mean))
test$grid <- as.numeric(as.character(test$grid))

weighted.mean(test$mean, test$grid) 

###### extended values
ad1_for_table <- fread(paste0('FILEPATH'))
ad1_for_table_dpt1 <- subset(ad1_for_table, ADM1_NAME == "Borno" & year %in% c(2016, 2017))
ad1_for_table_dpt1

ad1_for_table <- fread(paste0('FILEPATH'))
ad1_for_table_dpt3 <- subset(ad1_for_table, ADM1_NAME == "Borno" & year %in% c(2016, 2017))
ad1_for_table_dpt3

nga <- readRDS("FILEPATH")
nga

#########################################################################
##### Nigeria deep dive / Borno state case study, pt. 2
#########################################################################

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

# Custom load indicator-specific functions
source(paste0('FILEPATH'))
library(fasterize)

#### load in conflict data, mean raster, and shapefile
year <- 2020
i <- year - 1999

conflict <- fread('FILEPATH')
conflict2018 <- subset(conflict, year == 2018)
conflict2017 <- subset(conflict, year == 2017)

mean <- brick('FILEPATH')
mean <- mean[[i]]

ashp <- sf::st_read('FILEPATH') 
ashp <- subset(ashp, ADM0_NAME == "Nigeria")

masked <- mask(mean, ashp)
cropped <- crop(masked, extent(ashp))

plot(cropped)

indicator <- 'dpt3_cov'
indicator_group <- 'vaccine'
run_date <- 'RUN_DATE'

source('FILEPATH')
'%!in%' <- function(x,y)!('%in%'(x,y))

###############################################################################################################
ind = indicator
ig = indicator_group
rd = run_date
shapefile_version = "DATE"
yl = year_list
adm0_list = NULL
stat = "mean"
mask_output = T
lakes_and_rivers = T
ad1_borders = T
ad2_borders = F
year_list <- c(2000:2018)

adm0_list = get_adm0_codes("NGA", shapefile_version='DATE', subnational_raking=T)
plot_title = "IHME geospatial estimates"
legend_title = paste0("DTP3", "\nCoverage")
stat = "mean"
raked=T
mask_output = T
lakes_and_rivers = T

# Directories and files -----------------------------------------------------------------------
sharedir <- paste0("FILEPATH")
ras_file <- paste0(sharedir, ind, "_", stat, ifelse(raked, "_raked", ""), "_raster.tif")
if(indicator=='dpt1_cov'){
  ras_file <- ('FILEPATH')
}
# Checks to ensure everything exists
if (!dir.exists(sharedir)) stop("No directory found for this indicator and run date")
if (!file.exists(ras_file)) stop("Raster file for this statistic and indicator not found")

# Load shapefiles -----------------------------------------------------------------------------
ad0_shp <- prep_ad_shp_for_mapping(ad_level = 0,
                                   shapefile_version = "DATE", 
                                   adm0_list = adm0_list)

if (ad1_borders) {
  ad1_shp <- prep_ad_shp_for_mapping(ad_level = 1,
                                     shapefile_version = "DATE", 
                                     adm0_list = adm0_list)
}

borno_shp <- subset(ad1_shp, ADM1_NAME == "Borno")

ad2_shp <- prep_ad_shp_for_mapping(ad_level = 2,
                                   shapefile_version = "DATE", 
                                   adm0_list = adm0_list)
borno_ad2_shp <- subset(ad2_shp, ADM1_NAME == "Borno")

only_conflict_nga <- ggplot() + 
  geom_sf(aes(shape = event_type, size=fatalities, color=event_type), alpha=0.75, data = all)+ geom_sf(data = ad1_shp,
                                                                                                       color = "black",
                                                                                                       lwd = 0.3,
                                                                                                       fill = NA) +  theme_empty() + 
  scale_colour_manual(values = c("#fb8072", "#fdb462", "#ffffb3", "#8dd3c7", "#80b1d3", "#bebada")) + 
  theme(legend.position = "none") + 
  coord_sf(datum = NA) + labs(color = "Event type", size=NULL, shape=NULL) + scale_shape(guide = 'none') + scale_size(guide = 'none')   

only_conflict_nga

png(file = paste0("FILEPATH"),
    width = 10,
    height = 10,
    units = "in",
    res = 400,
    type = "cairo")
print(only_conflict_nga)
dev.off()

borno_conflict2018 <- subset(conflict2018, admin1 == "Borno")
borno_conflict2018 <- borno_conflict2018 %>% sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

borno_conflict2017 <- subset(conflict2017, admin1 == "Borno")
borno_conflict2017 <- borno_conflict2017 %>% sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

conflict2017 <- conflict2017 %>% sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
conflict2018 <- conflict2018 %>% sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

only_conflict_b <- ggplot() + geom_sf(data = borno_ad2_shp,
                                      color = "grey",
                                      lwd = 0.75,
                                      fill = NA) +
  geom_sf(data = borno_shp,
          color = "black",
          lwd = 1,
          fill = NA)+  theme_empty() + 
  geom_sf(aes(shape = event_type, size=fatalities, color=event_type), alpha=0.75, data = borno_conflict2018) +
  theme(legend.position = "none") + 
  scale_colour_manual(values = c("#fb8072", "#fdb462", "#ffffb3", "#8dd3c7", "#80b1d3", "#bebada")) + 
  coord_sf(datum = NA) + labs(color = "Event type", size=NULL, shape=NULL) + scale_shape(guide = 'none') + scale_size(guide = 'none')   

png(file = paste0("FILEPATH"),
    width = 8,
    height = 8,
    units = "in",
    res = 400,
    type = "cairo")
print(only_conflict_b)
dev.off()

only_conflict_borno_same_color_2018 <- ggplot() + 
  geom_sf(aes(), color="#9999CC", data = borno_conflict2018)+ geom_sf(data = borno_shp,
                                                                      color = "black",
                                                                      lwd = 0.3,
                                                                      fill = NA) +  
  theme_empty() + coord_sf(datum = NA) + theme(legend.background = element_rect())

only_conflict_borno_same_color_2017 <- ggplot() + 
  geom_sf(aes(), color="#9999CC", data = borno_conflict2017)+ geom_sf(data = borno_shp,
                                                                      color = "black",
                                                                      lwd = 0.3,
                                                                      fill = NA) +  
  theme_empty() + coord_sf(datum = NA) + theme(legend.background = element_rect())

only_conflict_nga_same_color_2018 <- ggplot() + 
  geom_sf(aes(), color="#9999CC", data = conflict2018)+ geom_sf(data = ad0_shp,
                                                                color = "black",
                                                                lwd = 0.3,
                                                                fill = NA) +  
  theme_empty() + coord_sf(datum = NA) + theme(legend.background = element_rect())

only_conflict_nga_same_color_2017 <- ggplot() + 
  geom_sf(aes(), color="#9999CC", data = conflict2017)+ geom_sf(data = ad0_shp,
                                                                color = "black",
                                                                lwd = 0.3,
                                                                fill = NA) +  
  theme_empty() + coord_sf(datum = NA) + theme(legend.background = element_rect())

png(file = paste0("FILEPATH"),
    width = 9,
    height = 9,
    units = "in",
    res = 400,
    type = "cairo")
print(only_conflict_borno_same_color_2017)
dev.off()

png(file = paste0("FILEPATH"),
    width = 9,
    height = 9,
    units = "in",
    res = 400,
    type = "cairo")
print(only_conflict_borno_same_color_2018)
dev.off()

png(file = paste0("FILEPATH"),
    width = 9,
    height = 9,
    units = "in",
    res = 400,
    type = "cairo")
print(only_conflict_nga_same_color_2017)
dev.off()

png(file = paste0("FILEPATH"),
    width = 9,
    height = 9,
    units = "in",
    res = 400,
    type = "cairo")
print(only_conflict_nga_same_color_2018)
dev.off()

####################################################################################

library(sf)
dhs <- fread('FILEPATH')
head(dhs)
dhs <- subset(dhs, me_name == "vacc_dpt3")
dhs <- subset(dhs, !is.na(latitude))
dhs <- dhs %>% sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

test <- st_intersection(dhs, borno_ad2_shp)

borno_dhs_points <- only_conflict_borno_same_color_2018 + 
  geom_sf(data = borno_ad2_shp,
          color = "grey",
          lwd = 0.75,
          fill = NA) +
  geom_sf(data = borno_shp,
          color = "black",
          lwd = 1,
          fill = NA) +  
  geom_sf(aes(size=N), color="#66CC99", shape=17, data = test) + theme_empty() + theme(legend.position = "none") + 
  coord_sf(datum = NA) 

png(file = paste0("FILEPATH"),
    width = 8,
    height = 8,
    units = "in",
    res = 400,
    type = "cairo")
print(borno_dhs_points)
dev.off()

############################# mics

mics <- fread('FILEPATH')
head(mics)
mics <- subset(mics, me_name == "vacc_dpt3")
mics <- subset(mics, !is.na(latitude))
mics <- mics %>% sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

test <- st_intersection(mics, borno_ad2_shp)
test_yr_dhs <- subset(test, year_id == "2015")

borno_mics_points <-  only_conflict_borno_same_color_2017 + 
  geom_sf(data = borno_ad2_shp,
          color = "grey",
          lwd = 0.75,
          fill = NA) +
  geom_sf(data = borno_shp,
          color = "black",
          lwd = 1,
          fill = NA) +  
  geom_sf(aes(size=N), color="#66CC99", shape=17, data = test) + theme_empty() + theme(legend.position = "none") + 
  coord_sf(datum = NA) 

png(file = paste0("FILEPATH"),
    width = 9,
    height = 9,
    units = "in",
    res = 400,
    type = "cairo")
print(borno_mics_points)
dev.off()

##### LGA maps
head(ad2_shp)

est <- fread('FILEPATH')
est <- subset(est, ADM0_NAME == "Nigeria")
est <- subset(est, year == 2017)
ad2_dpt3_2017 <- merge(ad2_shp, est, by=c("ADM2_CODE","ADM2_NAME"), all.x=T)
ad2_dpt3_2017$mean <- ifelse(ad2_dpt3_2017$ADM2_NAME== "Lake Chad", NA, ad2_dpt3_2017$mean)

est <- fread('FILEPATH')
est <- subset(est, ADM0_NAME == "Nigeria")
est <- subset(est, year == 2018)
ad2_dpt3_2018 <- merge(ad2_shp, est, by=c("ADM2_CODE","ADM2_NAME"), all.x=T)
ad2_dpt3_2018$mean <- ifelse(ad2_dpt3_2018$ADM2_NAME== "Lake Chad", NA, ad2_dpt3_2018$mean)

est <- fread('FILEPATH')
est <- subset(est, ADM0_NAME == "Nigeria")
est <- subset(est, year == 2017)
ad2_dpt1_2017 <- merge(ad2_shp, est, by=c("ADM2_CODE","ADM2_NAME"), all.x=T)
ad2_dpt1_2017$mean <- ifelse(ad2_dpt1_2017$ADM2_NAME== "Lake Chad", NA, ad2_dpt1_2017$mean)

est <- fread('FILEPATH')
est <- subset(est, ADM0_NAME == "Nigeria")
est <- subset(est, year == 2018)
ad2_dpt1_2018 <- merge(ad2_shp, est, by=c("ADM2_CODE","ADM2_NAME"), all.x=T)
ad2_dpt1_2018$mean <- ifelse(ad2_dpt1_2018$ADM2_NAME== "Lake Chad", NA, ad2_dpt1_2018$mean)

ad2_b_dpt3_2018 <- subset(ad2_dpt3_2018, ADM1_NAME.x == "Borno")
ad2_b_dpt3_2017 <- subset(ad2_dpt3_2017, ADM1_NAME.x == "Borno")
ad2_b_dpt1_2018 <- subset(ad2_dpt1_2018, ADM1_NAME.x == "Borno")
ad2_b_dpt1_2017 <- subset(ad2_dpt1_2017, ADM1_NAME.x == "Borno")

dpt3_lga_18 <- ggplot() + 
  geom_sf(data = ad2_shp,
          color = "black",
          lwd = 0.3,
          fill = NA) + 
  geom_sf(data = ad2_dpt3_2018,
          color = "black",
          lwd = 0.1,aes(fill=mean)) + scale_fill_vaccine() + 
  theme_empty() + coord_sf(datum = NA) + theme(legend.position = "none")

dpt3_lga_17 <- ggplot() + 
  geom_sf(data = ad2_shp,
          color = "black",
          lwd = 0.3,
          fill = NA) + 
  geom_sf(data = ad2_dpt3_2017,
          color = "black",
          lwd = 0.1,aes(fill=mean)) + scale_fill_vaccine() + 
  theme_empty() + coord_sf(datum = NA) + theme(legend.position = "none")

dpt1_lga_18 <- ggplot() + 
  geom_sf(data = ad2_shp,
          color = "black",
          lwd = 0.3,
          fill = NA) + 
  geom_sf(data = ad2_dpt1_2018,
          color = "black",
          lwd = 0.1,aes(fill=mean)) + scale_fill_vaccine() + 
  theme_empty() + coord_sf(datum = NA) + theme(legend.position = "none")

dpt1_lga_17 <- ggplot() + 
  geom_sf(data = ad2_shp,
          color = "black",
          lwd = 0.3,
          fill = NA) + 
  geom_sf(data = ad2_dpt1_2017,
          color = "black",
          lwd = 0.1,aes(fill=mean)) + scale_fill_vaccine() + 
  theme_empty() + coord_sf(datum = NA) + theme(legend.position = 'none')

png(file = paste0("FILEPATH"),
    width = 10,
    height = 10,
    units = "in",
    res = 400,
    type = "cairo")
print(dpt1_lga_17)
dev.off()

png(file = paste0("FILEPATH"),
    width = 10,
    height = 10,
    units = "in",
    res = 400,
    type = "cairo")
print(dpt1_lga_18)
dev.off()

png(file = paste0("FILEPATH"),
    width = 10,
    height = 10,
    units = "in",
    res = 400,
    type = "cairo")
print(dpt3_lga_17)
dev.off()

png(file = paste0("FILEPATH"),
    width = 10,
    height = 10,
    units = "in",
    res = 400,
    type = "cairo")
print(dpt3_lga_18)
dev.off()

dpt3_lga_borno_18 <- ggplot() + 
  geom_sf(data = borno_shp,
          color = "black",
          lwd = 0.3,
          fill = NA) + 
  geom_sf(data = ad2_b_dpt3_2018,
          color = "black",
          lwd = 0.1,aes(fill=mean)) + scale_fill_vaccine() + 
  theme_empty() + coord_sf(datum = NA) + theme(legend.position = 'none')

dpt1_lga_borno_18 <- ggplot() + 
  geom_sf(data = borno_shp,
          color = "black",
          lwd = 0.3,
          fill = NA) + 
  geom_sf(data = ad2_b_dpt1_2018,
          color = "black",
          lwd = 0.1,aes(fill=mean)) + scale_fill_vaccine() + 
  theme_empty() + coord_sf(datum = NA) + theme(legend.position = 'none')

library(sf)
displaced <- fread('FILEPATH')
head(mics)
displaced <- subset(displaced, !is.na(latitude))
displaced <- displaced %>% sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

test <- st_intersection(displaced, borno_ad2_shp)
test_new_mics <- test

dhs <- fread('FILEPATH')
head(dhs)
dhs <- subset(dhs, me_name == "vacc_dpt1")
dhs <- subset(dhs, !is.na(latitude))
dhs <- dhs %>% sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

test <- st_intersection(dhs, borno_ad2_shp)
test_dhs <- test

borno_mics_points_new <- dpt3_lga_borno_18 + 
  geom_sf(aes(), color="black", shape=1, data = test_new_mics, size=5) +   geom_sf(data = borno_ad2_shp,
                                                                                    color = "grey",
                                                                                    lwd = 0.1,
                                                                                    fill = NA) +  theme_empty() + theme(legend.background = element_rect()) + 
  coord_sf(datum = NA)

borno_mics_points_new_dhs_dpt3 <- borno_mics_points_new +
  geom_sf(aes(), color="black", shape=3, data = test_dhs, size=5) +   geom_sf(data = borno_ad2_shp,
                                                                              color = "grey",
                                                                              lwd = 0.1,
                                                                              fill = NA) +  theme_empty() + theme(legend.position = 'none') + 
  coord_sf(datum = NA)

borno_mics_points_new <- dpt1_lga_borno_18 + 
  geom_sf(aes(), color="black", shape=1, data = test_new_mics, size=5) +   geom_sf(data = borno_ad2_shp,
                                                                                    color = "grey",
                                                                                    lwd = 0.1,
                                                                                    fill = NA) +  theme_empty() + theme(legend.background = element_rect()) + 
  coord_sf(datum = NA)

borno_mics_points_new_dhs_dpt1 <- borno_mics_points_new +
  geom_sf(aes(), color="black", shape=3, data = test_dhs, size=5) +   geom_sf(data = borno_ad2_shp,
                                                                              color = "grey",
                                                                              lwd = 0.1,
                                                                              fill = NA) +  theme_empty() + theme(legend.position = 'none') + 
  coord_sf(datum = NA)

dpt3_lga_borno_17 <- ggplot() + 
  geom_sf(data = borno_shp,
          color = "black",
          lwd = 0.3,
          fill = NA) + 
  geom_sf(data = ad2_b_dpt3_2017,
          color = "black",
          lwd = 0.1,aes(fill=mean)) + scale_fill_vaccine() + 
  theme_empty() + coord_sf(datum = NA) + theme(legend.position = "none")

png(file = paste0("FILEPATH"),
    width = 8,
    height = 8,
    units = "in",
    res = 400,
    type = "cairo")
print(dpt3_lga_borno_18)
dev.off()

png(file = paste0("FILEPATH"),
    width = 8,
    height = 8,
    units = "in",
    res = 400,
    type = "cairo")
print(dpt3_lga_borno_17)
dev.off()

dpt1_lga_borno_18 <- ggplot() + 
  geom_sf(data = borno_shp,
          color = "black",
          lwd = 0.3,
          fill = NA) + 
  geom_sf(data = ad2_b_dpt1_2018,
          color = "black",
          lwd = 0.1,aes(fill=mean)) + scale_fill_vaccine() + 
  theme_empty() + coord_sf(datum = NA) + theme(legend.position = 'none')

dpt1_lga_borno_17 <- ggplot() + 
  geom_sf(data = borno_shp,
          color = "black",
          lwd = 0.3,
          fill = NA) + 
  geom_sf(data = ad2_b_dpt1_2017,
          color = "black",
          lwd = 0.1,aes(fill=mean)) + scale_fill_vaccine() + 
  theme_empty() + coord_sf(datum = NA) + theme(legend.position = "none")

png(file = paste0("FILEPATH"),
    width = 8,
    height = 8,
    units = "in",
    res = 400,
    type = "cairo")
print(dpt1_lga_borno_18)
dev.off()

png(file = paste0("FILEPATH"),
    width = 8,
    height = 8,
    units = "in",
    res = 400,
    type = "cairo")
print(dpt1_lga_borno_17)
dev.off()

png(file = paste0("FILEPATH"),
    width = 8,
    height = 8,
    units = "in",
    res = 400,
    type = "cairo")
print(borno_mics_points_new_dhs_dpt1)
dev.off()

png(file = paste0("FILEPATH"),
    width = 8,
    height = 8,
    units = "in",
    res = 400,
    type = "cairo")
print(borno_mics_points_new_dhs_dpt3)
dev.off()
