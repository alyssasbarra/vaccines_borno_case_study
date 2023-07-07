# USERNAME
# DATE <- DATE <- DATE
# Make Admin 1 and 2 level summary draw-level objects
# Based on USERNAME code: FILEPATH
# Writing own because not wanting to wait for PR

# indicator <- "dpt3_cov"
# indicator_group <- "vaccine"
# run_date <- "RUN_DATE"
# raked <- T
# pop_measure <- "a0004t"
# overwrite <- T
# age <- 0
# holdout <- 0
# region <- "vax_essa"
# rasterize_field <- "GAUL_CODE"
# subset_field <- "ISO"
# subset_field_codes <- c("KEN")
# filename_addin <- "_KEN"
# shapefile_path <- "FILEPATH"
#     Can also just pass a shapefile (without .shp) at the end
# sp_hierarchy_columns <- c("NAME_0", "NAME_1", "GAUL_CODE")
#     Has some defaults in there if needed
# crop_shapefile <- TRUE

# Setup
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## clear environment
rm(list=ls())

# Use parallelize infrastructure (only uses 5 commandArgs)
source("FILEPATH")

# What to load through parallelize: 
#  - varying: c("shapefile", "filename_addin", "region", "subset_field", "subset_field_codes")
#  - static: c("indicator", "indicator_group", "run_date", "raked", "pop_measure",
#              "overwrite", "age", "holdout", "crop_shapefile", "use_lookup", core_repo")

load_from_parallelize()  

## sort some directory stuff
commondir      <- sprintf('FILEPATH')
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0('FILEPATH'))
mbg_setup(package_list = package_list, repos = core_repo)

## -------------------------------------------------------------------------------------------------
## Start aggregation script here -------------------------------------------------------------------
## -------------------------------------------------------------------------------------------------

reg <- region #convenience

# Use lookup table to get shapefile attributes
df_shapefile_lookup <- fread("FILEPATH")
if (!(shapefile %in% df_shapefile_lookup$shp)) stop("Need to add a row to lookup table for this shapefile!")
assign("shapefile_path", paste0(df_shapefile_lookup[shp == shapefile]$shapefile_directory, 
                              df_shapefile_lookup[shp == shapefile]$shp, ".shp"))
assign("sp_hierarchy_columns", eval(parse(text = df_shapefile_lookup[shp == shapefile]$sp_hierarchy_columns)))

# A few checks to start; defining some objects
if (!(file.exists(shapefile_path))) stop("Shapefile not found!")

  filename_addin <- ""
} else if (is.null(filename_addin)) {
  filename_addin <- ""
}

if (!is.null(subset_field)) {
  if (is.na(subset_field) | subset_field == "NA") subset_field <- NULL
  if (is.na(subset_field_codes) | subset_field_codes == "NA") subset_field_codes <- NULL
}
rasterize_field <- field_name

# set share_dir
share_dir <- paste0("FILEPATH")

# set output directory
if (!exists("output_dir")) output_dir <- paste0(share_dir, "custom_aggregation2/", stringr::str_match(basename(shapefile_path), "(.*).shp")[,2], "/")
dir.create(output_dir, recursive = T, showWarnings = F)

# get year_list from the config file
this_config <- fread(paste0(share_dir, 'config.csv'))
year_list <- this_config[V1 == 'year_list', V2]
if (is.character(year_list)) year_list <- eval(parse(text = year_list))

# Determinining whether this is a re-run, stopping if overwrite is off and the file already exists.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(overwrite==T){
  if(file.exists(paste0(share_dir, indicator, "_admin_draws", ifelse(raked,"_raked.Rdata", ".RData")))){
    message("This file already exists, this must be a re-run.")
  }
}else{
  if(file.exists(paste0(share_dir, indicator, "_admin_draws", ifelse(raked,"_raked.Rdata", ".RData")))){
    stop("Looks like you have already prepped that admin draw file; skipping this!")
  }
}

# These are where draw-level, admin-level results will be stored
message("Starting empty lists for results")
regions<-list()
admin_list<-list()
sp_hierarchy_list<-list()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Looping through Regions, aggregating data to Admin Units
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
start <-proc.time()

message(paste0("Pulling numbers for indicator: ",indicator," region: ",reg," rundate: ",run_date," using the population measure ",pop_measure," for weighting."))
pathaddin<-paste0('_bin',age,'_',reg,'_',holdout)

# Loading Data
#~~~~~~~~~~~~~~
message("Loading in the draw-level data for this region; this may take a while (but probably not longer than 10 minutes).")
# Load in cell-level results based on that same model run
if(raked){
  # Different file conventions in use; try either rds or rdata until standardized
  filename_rds <- paste0(share_dir,reg,"_raked_cell_pred.RDs")
  filename_rdata <- paste0(share_dir, indicator, "_", ifelse(raked, "raked", "unraked"), "_", "cell_draws_eb_bin0_", reg, "_", holdout, ".RData")
  if (file.exists(filename_rds)) cell_pred<-readRDS(filename_rds)
  if (file.exists(filename_rdata)) {
    load(filename_rdata, verbose = T)
    cell_pred <- raked_cell_pred
    rm(raked_cell_pred)
  }
} else {
  load(paste0(share_dir, indicator, "_cell_draws_eb", pathaddin, ".RData"))
}

# Check if load was successful; stop if not
if (!exists("cell_pred")) stop("Unable to load cell_pred object! Check to make sure that the relevant object exists.")

# Getting the simple polygon and simple raster objects for this region alone
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
message("Getting the spatial objects associated with this region.")
simple_polygon_list <- load_simple_polygon(gaul_list = get_adm0_codes(reg, shapefile_version= 'DATE'), buffer = 0.4, subset_only = FALSE, shapefile_version='DATE')
subset_shape   <- simple_polygon_list[[1]]
simple_polygon <- simple_polygon_list[[2]]
message("Building simple raster from subset_shape")
raster_list    <- build_simple_raster_pop(subset_shape)
simple_raster  <- raster_list[['simple_raster']]
pop_raster     <- raster_list[['pop_raster']]
rm(raster_list,simple_polygon_list,pop_raster);gc()
message("All done loading spatial template files (subset_shape,simple_polygon,simple_raster,pop_raster, etc)")

# Determining a list of the valid pixel indices based on the simple raster template
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pixel_id <- seegSDM:::notMissingIdx(simple_raster)
pixel_spatial<-data.table(pixel_id=pixel_id)
message("Pixel ID that links cell_pred observations to raster locations created from simple raster, checking for sameness in dimensions compared to cell_pred.")
if(!(length(pixel_id)==nrow(cell_pred)/length(year_list))){
  stop("Excuse me, but the number of valid pixels in your simple raster is not the same number of valid pixels in your cell_pred object. Look into this!")
}else{
  message(paste0("Check passed: There are ",length(pixel_id),
                 " pixels in your simple_raster object, the same number of pixels in the cell_pred object for each year."))
}

# Pulling Population
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Pull annual population brick using new covariates function
message(paste0("Pulling ",pop_measure," population raster."))
pop <- load_and_crop_covariates_annual(covs = 'worldpop',
                                       measures = pop_measure, # Defined above
                                       simple_polygon = simple_polygon,
                                       start_year  = min(year_list),
                                       end_year    = max(year_list),
                                       interval_mo = 12,
                                       agebin=1)
message("Ensuring that the spatial extents of the population matches the simple raster")
pop<-crop_set_mask(pop[[1]],simple_raster)
message("Done generating population rasterbrick; extracting relevant pixels based on simple_raster.")
pop <- data.table(extract(pop, pixel_id)) # getting the pixels from the population raster that correspond to the valid pixels in the pop-raster
message("Reshaping population from wide to long, and generating an id variable that corresponds to the simple raster.")
pop[,pixel_id:=pixel_id]
pop<-melt(pop,id.vars="pixel_id") # Melting the dataframe such that it is long () and should match the number of rows in cell_pred
pop[, year := (min(year_list) - 1) + as.numeric(gsub("worldpop.", "", variable))] # Converting "worldpop.1" variables to actual years.
pop<-pop[,list(pixel_id,year,pop=value)] # Subsetting columns
pop[is.na(pop),pop:=0] # Setting values where pop is NA to 0
message("Checking whether the population has the same dimensions as cell_pred")
if(!(nrow(pop)==nrow(cell_pred))){
  stop("Excuse me, but the number of rows in your population table is *not* the same as the number of rows in your cell predictions. Look into this!")
}else{
  message("Check passed: The population file has the same number of rows as the cell predictions after transformation to long format.")
}


# Loading and Assigning Admin Units to Pixels
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
message("Getting the spatial location (admin unit) of each of the pixel locations.")

region_adm0_list<-get_adm0_codes(reg, shapefile_version='DATE') # Getting the adm0 GAUL codes, we can use this to make sure we don't accidentally include countries from buffers that shouldn't be in this region

# Defining a function that will get the raster versions of each Admin level:
GetAdmin<-function(simple_raster, region_adm0_list, rasterize_field, subset_field = NULL, subset_field_codes = NULL){
  
  # load admin shape file
  admin_shp <- sf::st_read(shapefile_path)
  
  if (!is.null(subset_field) & !is.null(subset_field_codes)) {
    message("Subsetting to `subset_field_codes`")
    admin_shp <- subset(admin_shp, get(subset_field) %in% subset_field_codes)
  }
  
  message("Rasterizing...")
  if (!("Shape_Area" %in% names(admin_shp))) admin_shp$Shape_Area <- sf::st_area(admin_shp)
  admin_rast<-fasterize::fasterize(sf = admin_shp[order(admin_shp$Shape_Area),],
                                   raster = simple_raster,
                                   field = rasterize_field, 
                                   fun="first")
  
  message("Converted to raster based on simple_raster template. Cropping and masking:")
  admin_rast  <- crop(admin_rast, extent(simple_raster))
  admin_rast  <- extend(admin_rast, simple_raster)
  admin_rast  <- mask(admin_rast, simple_raster)
  
  message("Subsetting polygon and point objects to only contain the relevant ADM0 codes; calculating centroids.")
  admin_centroids<-SpatialPointsDataFrame(gCentroid(as(admin_shp, 'Spatial'), byid=TRUE), as(admin_shp, 'Spatial')@data, match.ID=FALSE)
  
  message("Compiling and returning results.")
  admin<-list()
  admin[["spdf"]]<-as(admin_shp, 'Spatial')
  admin[["centroids"]]<-admin_centroids
  admin[["rast"]]<-admin_rast
  admin[["attributes"]]<-copy(data.table(as(admin_shp, 'Spatial')@data))
  
  return(admin)
}

message("Rasterizing shapefiles; this may take a while.")
admin_info<-GetAdmin(simple_raster,region_adm0_list, rasterize_field, subset_field = subset_field, subset_field_codes = subset_field_codes)
pixel_spatial[[rasterize_field]]<-extract(admin_info[["rast"]],pixel_spatial$pixel_id) # Generate a field based on the admin boundary unit that has the ID code in the pixel_spatial data.table
if(sum(is.na(pixel_spatial[[rasterize_field]]))>0){ # Check to see if any of the pixels don't have a location assigned
  message("   Caution: there are some pixels that are NA, and have not been assigned a location.\n   This may be normal, however, if you're only aggregating a subset of your cell_pred")
}


# Merging together cell_pred and spatial information, population information.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pop<-merge(pop,pixel_spatial,by="pixel_id",all.x=T) # Merging on the spatial information to the population information.
pop<-pop[order(year,pixel_id)] # Re-ordering the pop object by year such that pixels ascent,a nd years ascend (same as cell_pred)
cell_pred<-data.table(cell_pred) # Converting cell_pred to a data.table
draw_colnames<-names(cell_pred)
cell_pred<-cbind(cell_pred,pop) # Adding on the population and spatial information to cell_pred.

# Loading and Assigning Pixels to Admin Units (that are missing)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# We now need to figure out what admin units don't show up when we pull the locations of individual raster pixels of our results.
# To do this, we need to know what Admin 1 and 2 units exist in the shapefiles we want to use, but don't show up in the results.

message("Generating a raster of the pixel_id values")
pixel_id_raster<-copy(simple_raster)
for (pix in pixel_id){ # For each of the valid pixels, replace the value with the pixel_id, not whatever country ID is used in the simple_polygon/raster
  pixel_id_raster@data@values[pix]<-pix
}

missing_admins<-list() # This will contain the cell_preds for otherwise missing admin units, by level.
fieldname<-rasterize_field

# First, we discover what GAUL codes are missing from the pixel_spatial table compared to the shapefile:
shpfile_attributes<-admin_info[["attributes"]] # Get the attribute table from the admin level's shapefile

attribute_codes<-shpfile_attributes[[fieldname]] # Get list of codes for that level from the admin level's attribute table
pixel_codes<-pixel_spatial[[fieldname]] # Get list of codes based on what's in the pixel_spatial object
missing_codes<-attribute_codes[!(attribute_codes %in% pixel_spatial[[fieldname]])] # Get list of missing codes
missing_table<-shpfile_attributes[shpfile_attributes[[fieldname]]%in%missing_codes,fieldname,with=F]

if(length(missing_codes)==0){
  message("No missing codes")
}else{
  message("Missing codes found")
  # Strategy 1: Assign a pixel location based on centroid location
  # Develop a raster of the pixel-IDs:
  message("  Discovering centroid location")
  points<-admin_info[["centroids"]] # SpatialPointsDataFrame of that level
  missing_points<-points[(points@data[[fieldname]] %in% missing_codes),] # getting only the missing points
  missing_centroid_locs<-data.table(raster::extract(pixel_id_raster,missing_points)) # Extracting the missing location values as centroid points...
  names(missing_centroid_locs)<-"point_method"
  missing_admins_centroids<-data.table(missing_codes,missing_centroid_locs)
  
  # Strategy 2: Assign a pixel location based on polygon extraction
  # Develop a raster of the pixel-IDs:
  message("  Discovering first raster pixel touching polygon")
  polys<-admin_info[["spdf"]] # SpatialPointsDataFrame of that level
  missing_polys<-polys[(polys@data[[fieldname]] %in% missing_codes),] # getting only the missing points
  missing_poly_locs<-data.table(raster::extract(x=pixel_id_raster,y=missing_polys,small=T,fun=function(x,...)first(x))) # Extracting the missing location values as polygons, pulling the first raster cell that it touches...
  names(missing_poly_locs)<-"poly_method"
  missing_admins_polys<-data.table(missing_codes,missing_poly_locs) # Extracting the
  
  # Merging strategies together: centroids and polygons; adding to list.
  missing_locs<-merge(missing_admins_polys,missing_admins_centroids,by="missing_codes")
  setnames(missing_locs,"missing_codes",fieldname)
  missing_locs<-merge(missing_locs,missing_table,by=fieldname) # Add in admin 0, 1 levels if relevant
  missing_locs[,pixel_id:=point_method] # If centroid produced something, go with centroid
  missing_locs[is.na(point_method),pixel_id:=poly_method] # Otherwise, go with the polygon method.
  
  # For those still missing, assign them to nearest non-NA pixel using centroid
  if(NA %in% missing_locs$pixel_id){
    message('  After centroids and polygon methods, there are still some missing admins.')
    message('  Now sampling to find nearest non-NA pixel and using that')
    
    for(rr in which(is.na(missing_locs$pixel_id))){
      message(sprintf('  -- finding nearest pixel to gaul_code: %i',
                      as.numeric(missing_locs[rr, fieldname, with = F])))
      
      ## get centroid of chape
      mp <- polys[polys[[fieldname]] == missing_locs[[fieldname]][rr],]
      cent <- getSpPPolygonsLabptSlots(mp)
      
      ## loop through withh an increasing radius and see if nearby points are non-NA
      found <- 0
      radius <- .005
      while(found != 1 & radius < 1.5){ ## stop for max radius or match found
        ## sample 1000 nearby locs
        near <- matrix(runif(2000, -radius, radius), ncol = 2)
        near[, 1] <- near[, 1] + cent[1, 1]
        near[, 2] <- near[, 2] + cent[1, 2]
        
        ## extract raster pixels
        near <- data.table(raster::extract(pixel_id_raster, near), near)
        colnames(near) <- c('pixel_id', 'x', 'y')
        
        if(mean(is.na(near[,pixel_id])) < 1){ ## then we've found a non-NA neighbor
          found <- 1 ## end while loop
          message(sprintf('  ---- found neighbor using radius: %f', radius))
          
          ## find closest neighbor
          near <- na.omit(near) ## those with non-NA pixels
          dist <- sqrt((near[, x] - cent[1, 1]) ^ 2 +
                         (near[, y] - cent[1, 2]) ^ 2)
          min.ind <- which.min(dist) ## in case of tie, this returns 1st
          
          ## take the pixel id of the nearest sampled neighbor
          missing_locs[rr, pixel_id := near[min.ind, pixel_id] ]
        }
        
        ## increase radius in case we didn't catch anything
        radius <- radius + .005
        
      } ## end while loop
    } ## for each admin with NA pixel_id
  } ## if any NA pixel ids after centroid and poly methods
  
  
  # Check to see if NAs still exist in missing locations, make a warning about missing areas.
  if(NA %in% missing_locs$pixel_id){
    message( "The following admin units appear to be getting lost forever:")
    print(missing_locs[is.na(pixel_id),c("pixel_id",fieldname),with=F])
  }
  
  # Merging on locations with pixel IDs
  missing_locs<-missing_locs[!is.na(pixel_id),c("pixel_id",fieldname),with=F]
  missing_locs<-merge(missing_locs,cell_pred[,c("pixel_id","pop","year",draw_colnames),with=F],by="pixel_id",allow.cartesian=T)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Collapsing down draw information using weighted means, based on population.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

message("Generating data.tables with draws for each admin level")

# Note that in the standard aggregation script there are checks to ensure that
# a region-level draws object is created to ensure that there aren't any overlaps 
# between regions -- here, that's not done since there aren't regions to line up 
# -- just aggregating one-off shapefiles. generating entire shapefile estimates
# here but not currently using them at all.

entire_shapefile<-cell_pred[!is.na(get(fieldname)),lapply(.SD,weighted.mean,w=pop, na.rm=T),by=list(year), .SDcols=draw_colnames]
entire_shapefile[,region_name:="entire shapefile"]

# Adding in the missing admin information
if (exists("missing_locs")) {
  adm <- rbind(cell_pred, missing_locs, fill = T)
} else {
  adm <- copy(cell_pred)
}

setnames(adm, rasterize_field, "the_field_name")
adm <- adm[!is.na(the_field_name), lapply(.SD,weighted.mean,w=pop, na.rm=T),by=list(the_field_name,year), .SDcols=draw_colnames]
setnames(adm, "the_field_name", rasterize_field)

# Generate total populations for each level, make sure to account for missing admins
pop_entire_shapefile<-pop[!is.na(get(rasterize_field)),list(pop=sum(pop)),by=list(year)]
pop_entire_shapefile[,region_name:="entire_shapefile"]

pop_adm <- copy(pop)
setnames(pop_adm, rasterize_field, "the_field_name")
pop_adm <- pop_adm[!is.na(the_field_name), list(pop=sum(pop)),by=list(the_field_name,year)]

if (exists("missing_locs")) {
  missing_adm_pops <- subset(missing_locs, !(get(rasterize_field) %in% pop_adm$the_field_name), select = c(fieldname, "year", "pop"))
  setnames(missing_adm_pops, rasterize_field, "the_field_name")
  pop_adm <- rbind(pop_adm, missing_adm_pops, use.names = T)                      
}

setnames(pop_adm, "the_field_name", rasterize_field)

# Merge results together
admin_results <- merge(adm, pop_adm, by = c("year", rasterize_field))

# Defining the hierarchy of what lives within what:
sp_hierarchy<-unique(subset(admin_info[["attributes"]], select = sp_hierarchy_columns)) 
sp_hierarchy[,region:=reg]

# Merge together
all_results <- merge(sp_hierarchy, admin_results, by=fieldname, all.x=T, all.y=T)

elapsed <- proc.time() - start
print(elapsed)

# Collapsing and Saving Results
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(all_results, file = paste0(output_dir, indicator, "_", ifelse(raked, "raked", "unraked"), "_admin_draws_eb_bin", age, "_", reg, "_", holdout, filename_addin,".rds"))

# Summarize
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (!exists("summ_stats")) summ_stats <- c("mean", "upper", "lower")

all_results[, a_row_id := .I]
draw_cols <-names(all_results)[grep("V[0-9]+", names(all_results))]
draws <- subset(all_results, select = c("a_row_id", draw_cols))

non_draw_cols <- names(all_results)[-grep("V[0-9]+", names(all_results))]
non_draws <- subset(all_results, select = non_draw_cols)

df_summstats <- lapply(summ_stats, function(ss) {
  df_ss <- draws[, list(summ = apply(.SD, 1, ss),
                        a_row_id = a_row_id), 
                 .SDcols = draw_cols]
  setnames(df_ss, "summ", ss)
  return(df_ss)
})

df_summstats <- Reduce(merge, df_summstats)

df_summstats <- merge(non_draws, df_summstats, by = "a_row_id")
df_summstats[, a_row_id := NULL]

fwrite(df_summstats, file = paste0(output_dir, indicator, "_admin_", ifelse(raked, "raked", "unraked"), "_summary", filename_addin, ".csv")) 

# Plot
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

message("Plotting results...")

# Read in mean pixel estimates
mean_ras <- brick(paste0(share_dir, indicator, "_mean", ifelse(raked, "_raked", ""), "_raster.tif"))
mean_ras <- crop(mean_ras, admin_info[["spdf"]])
mean_ras <- mask(mean_ras, admin_info[["spdf"]])

plot_dir <- paste0(output_dir, "plots/")
dir.create(plot_dir)

for (yy in year_list) {
  
  message(paste0("  ", yy, "..."))
  year_ras <- mean_ras[[which(year_list == yy)]]
  ras_df <- as.data.frame(as(year_ras, "SpatialPixelsDataFrame")) %>% as.data.table
  names(ras_df) <- c("mean", "x", "y")

  # Get the year data and convert to spatial format - aggregates
  agg_pred <- subset(df_summstats, year == yy, select = c(fieldname, "mean"))
  the_spdf <- admin_info[["spdf"]]
  the_spdf <- merge(the_spdf, agg_pred, by = fieldname)
  the_spdf@data$id <- rownames(the_spdf@data)
  the_spdf.points <- fortify(the_spdf, region = "id")
  the_spdf.df <- join(the_spdf.points, the_spdf@data, by = "id") %>% as.data.table
  
  # Custom vaccines color scale -----------------------------------------------------
  vals <- c(1.000000,  0.9,       0.8,        0.6,      0.496094,   0.4,       0.2,       0.164063,  0.000000)
  cols <- c("#8436A8", "#3D649D", "#91BEDC",  "#DEF3F8", "#FFE291", "#FEE191", "#FC8D58", "#FC8C58", "#B03027")
  
  geom_polygon_quiet <- function(...) {suppressMessages(ggplot2::geom_polygon(...))}
  geom_path_quiet     <- function(...) {suppressMessages(ggplot2::geom_path(...))}
  
  gg_agg <- ggplot() +
    geom_polygon_quiet(data = the_spdf.df,
                       aes(x=long, y=lat, group=group, fill=mean),
                       size = 0.2) +
    geom_path_quiet(data = the_spdf.df,
                    aes(x=long, y=lat, group=group),
                    size = 0.2,
                    color = "black") +
    scale_fill_gradientn(colors = cols, values = vals, 
                         breaks = c(0,0.2,0.4,0.6,0.8,0.9,1),
                         limits = c(0,1),
                         labels = scales::percent) +
    theme_empty() +
    coord_equal() +
    labs(title = "Aggregated estimates",
         subtitle = ifelse(raked, "(Raked)", "(Unraked)"))
  
  gg_pix <- ggplot() +
    geom_tile(data = ras_df,
              aes(x=x, y=y, fill=mean)) +
    geom_path_quiet(data = the_spdf.df,
                    aes(x=long, y=lat, group=group),
                    size = 0.2,
                    color = "black") +
    scale_fill_gradientn(colors = cols, values = vals, 
                         breaks = c(0,0.2,0.4,0.6,0.8,0.9,1),
                         limits = c(0,1),
                         labels = scales::percent) +
    theme_empty() +
    coord_equal() +
    labs(title = "5x5 km estimate",
         subtitle = ifelse(raked, "(Raked)", "(Unraked)"))
  
  g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}
  
  the_legend <- g_legend(gg_pix)
  
  title_grob <- textGrob(paste0(indicator, ": ", yy), gp = gpar(fontsize = 24, fontface = "bold"))
  
  gg_pix <- gg_pix + theme(legend.position = "none")
  gg_agg <- gg_agg + theme(legend.position = "none")
  
  lay <- rbind(c(1,1,1,1,1,1,1),
               c(2,2,2,3,3,3,4),
               c(2,2,2,3,3,3,4),
               c(2,2,2,3,3,3,4))
  
  plot_all <- arrangeGrob(title_grob, gg_pix, gg_agg, the_legend, 
                          layout_matrix = lay,
                          heights = c(0.2,1,1,1))
  
  png(file = paste0(plot_dir, indicator, ifelse(raked, "_raked","_unraked"), "_aggregated", filename_addin, "_", yy, ".png"),
      width = 14,
      height = 8,
      units = "in", 
      res = 300)
  grid.draw(plot_all)
  dev.off()
  
}


# Finish up
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

message("Done!")
