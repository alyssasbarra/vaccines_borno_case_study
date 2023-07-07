## expand_bbox() ----------------------------------------------
#'
#' @title expand_bbox
#'
#' @description Expand a given country's bounding box to encompass a larger area than
#' the default extent of the country.
#'
#' @param bound_box A BoundingBox object from a country shapefile
#' @param deg_x_neg Decimal degrees in the -x direction to expand the bounding box (default .5).
#' @param deg_x_pos Decimal degrees in the +x direction to expand the bounding box (default .5).
#' @param deg_y_neg Decimal degrees in the -y direction to expand the bounding box (default .5).
#' @param deg_y_pos Decimal degrees in the +y direction to expand the bounding box (default .5).
#' @return BoundingBox object with borders expanded
expand_bbox <- function(bound_box, deg_x_neg=.5, deg_x_pos=.5, deg_y_neg=.5, deg_y_pos=.5) {
  new_bound_box <- copy(bound_box)
  new_bound_box@xmin <- new_bound_box@xmin - deg_x_neg
  new_bound_box@xmax <- new_bound_box@xmax + deg_x_pos
  new_bound_box@ymin <- new_bound_box@ymin - deg_y_neg
  new_bound_box@ymax <- new_bound_box@ymax + deg_y_pos
  return(new_bound_box)
}

## country_data_map() ----------------------------------------------
#'
#' @title Country Data Map
#'
#' @description Produces a subnational plot of a given indicator
#'
#' @param df Input dataframe of values to plot. Must contain 'uid' and a column mapping to 'year', 
#' (usually year_id). Typically this
#' will be provided by the prep_country_vr_data() function output.
#' Here is the structure of a dataframe/data table that will work with the function. Note that
#' instead of  total deaths, you can supply another column to provide your 
#' values to map, as long as they line up with a uid and year id 
#' for a given country.
#'
#' year_id  uid total_deaths     
#' 2000      0        13042
#' 2000      1          168  
#' 2000      2           20   
#' 2000      3           54
#' 2000      4          145
#' ......
#' 2017   1110            9
#' 2017   1111          207
#' 2017   1112          200
#' 2017   1113          251
#' 2017   1114           42
#'
#' @param adm_level Most granular admin level to plot at. Only supported for admin level 2
#' at present
#' @param id_col Column of values in the input dataframe to plot.
#' @param merge_col Column in input dataframe containing the unique codes to match
#' with that country's VR specific shapefile (default 'uid')
#' @param iso3 The corresponding iso3 code matching the country of interest,
#' necessary for cropping plot to the correct area
#' and matching the appropriate cities, if cities are being plotted
#' @param years Vector of years to produce plots for. If provided as 'ALL', no subsetting
#' by year will be done on the input dataframe. 'ALL' should only be used if an already aggregated
#' dataframe is being used as input or undesirable behavior will result.
#' @param shapefile Path to the VR shapefile to use for plotting. 
#' @param cities_overlay Whether to include labeled city points in the country map.
#' At present, city overlays are always included and removing them is not yet implemented.
#' @param title Title string or vector of title strings to override the default map title
#' generation of 'Year: {year}'. If single string provided, same title will be applied
#' across all maps. If a vector of titles is provided, the number of tiles must match the number
#' of years.
#' @param legend_title String to use for legend title in output plot
#' @param expand_bounds If TRUE, will expand the mapping area by two decimal degrees in
#' each direction. Alternatively, custom map boundaries may be set by
#' setting this to TRUE and passing in a vector of custom boundaries.
#' @param min Minimum value to set on legend. If not specified will be automatically
#' determined by ggplot.
#' @param min Maximum value to set on legend. If not specified will be automatically
#' determined by ggplot.
#' @param transform_arg Transformation that should be applied to the legend scaling.
#' (default 'identity'). If set to 'identity', no transformation applied. If 'log' provided
#' will apply a log transform to the legend scale.
#' If 'custom_cv' provided, will provide a custom legend for plotting coefficient of
#' variation scaled from '0' to '3+'
#' If 'custom_cv_0_1' provided, will provide a custom legend for plotting coefficient of
#' variation scaled from '0' to '1+'
#' @param cities_shapefile Path to a specific cities shapefile to use for plotting. 
#' By default will point to 'FILEPATH'. If custom
#' shapefile provided, field names will need to match the field names in the World Cities shapefile.
#' @param num_cities How many of the top cities (by population) to plot on map (default 5)
#' @param text_scale_factor Scaling adjustment to make map text larger or smaller. Useful if
#' placing more/less than 3 maps in a panel. Default is 1 but if text becomes too small try setting
#' to 1.5 or 2. Conversely, if text is too large, try setting to .85 or .6.
#' @param legend_position Position vector (c(x, y)) defining where to place the legend on map
#' By default places legend in bottom right (c(1,0)). Depending on the geometry of a country it
#' may be desirable to move the legend to a different position when plotting.
#' @param color_scheme Vector of color codes to use for plotting. 
#' Also accepts a RColorBrewer palette. By default will use a reverse Spectral 9-color palette.
#' @param custom_bounds Override to use custom defined boundaries instead of the 
#' default shapefile extent if passed a vector of (xmin, xmax, ymin, ymax) values.
#' Unexpected behavior may occur if overriding shapefile boundaries. Must be used with
#' parameter expand_bbox set to TRUE.
#'
#' @example output_maps <- country_data_map(df = vr_results,
#' adm_level = 2,
#' id_col = 'sum_cause_deaths',
#' iso3 = 'BRA',
#' shapefile = paste0("FILEPATH"),
#' years = c(2000,2005,2010,2015),
#' min = 0,
#' max = 50,
#' legend_position = c(1,0),
#' legend_title ='Total\nDeaths',
#' color_scheme = RColorBrewer::brewer.pal("BuPu", n=7),
#' text_scale_factor = .9
#'
#' @return List of ggplot objects, one for each year of data provided. If year is 'ALL', wil
#' return a single map.
country_data_map <- function(df,
                             adm_level = 2,
                             id_col,
                             merge_col = 'uid',
                             iso3,
                             years,
                             shapefile,
                             cities_overlay = TRUE,
                             title = NULL,
                             legend_title = NULL,
                             expand_bounds = TRUE,
                             min = NULL, 
                             max = NULL,
                             transform_arg = 'identity',
                             cities_shapefile = 'FILEPATH',
                             num_cities = 5,
                             text_scale_factor = 1,
                             legend_position = c(1,0),
                             color_scheme = NULL,
                             custom_bounds = NULL) {
  
  package_list <- c('data.table', 'ggplot2', 'rgdal', 'scales', 'grid',
                    'gridExtra', 'units', 'RColorBrewer', 'rgdal', 'sp',
                    'raster', 'sf', 'ggrepel', 'viridis')
  load_R_packages(package_list)
  if(!merge_col %in% names(df)) {
    stop("Missing declared merge column.")
  }

  country_shp <- sf::st_read(shapefile)
  ad0_sf <- suppressWarnings(suppressMessages(
    sf::st_read(get_admin_shapefile(admin_level=0))
  ))
  ad1_sf <- suppressWarnings(suppressMessages(
  sf::st_read(get_admin_shapefile(admin_level=1))
  ))
  lookup_table <- load_adm0_lookup_table()
  lookup_table$iso3 <- lapply(lookup_table$iso3, FUN = function(x){toupper(x)})
  ad1_sf <- merge(ad1_sf, lookup_table, by.x= 'ADM0_CODE', by.y='gadm_geoid', all.x=TRUE)
  ad0_sf <- merge(ad0_sf, lookup_table, by.x= 'ADM0_CODE', by.y='gadm_geoid', all.x=TRUE)
  cities_shp <- sf::st_read(cities_shapefile)
 
  if (is.null(color_scheme)) {
  color_scheme <- rev(RColorBrewer::brewer.pal('Spectral', n=9))
  }
  if(!is.null(custom_bounds)) {
    if (length(custom_bounds != 4)) {
      stop("Custom bounds must be a vector of 4 boundary positions.
           See the function documentation for more details.")
    }
    expanded_bbox <- raster::extent(custom_bounds[1],
                                    custom_bounds[2],
                                    custom_bounds[3],
                                    custom_bounds[4])
  } else {
    bounds <- extent(country_shp)
    bbox <-c(bounds@xmin, bounds@xmax, bounds@ymin, bounds@ymax)
    expanded_bbox <- expand_bbox(raster::extent(country_shp))
  }
  #Crop shapefile to match extent of VR country of interest
  iso3 <- toupper(iso3)
  cropped_ad0 <- sf::st_crop(x = ad0_sf, y=expanded_bbox)
  cropped_ad1 <- sf::st_crop(x = ad1_sf, y=expanded_bbox)
  ad0_subset <- cropped_ad0[which(cropped_ad0$iso3 == iso3),]
  ad0_inverse_subset <- cropped_ad0[which(cropped_ad0$iso3 != iso3),]
  ad1_subset <- cropped_ad1[which(cropped_ad1$iso3 == iso3),]
  merged_shp <- merge(country_shp, df, by.x='GAUL_CODE', by.y=merge_col, all.x=TRUE)
  SCALE_SIZE <- text_scale_factor
  SCALE_TEXT_SIZE <- ((SCALE_SIZE * 7.34))
  
  #Build custom map theme
  theme_custom_map <- theme_set(theme_bw())
  theme_custom_map <- theme_update(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.border = element_blank(),   
    legend.key.height = unit((SCALE_SIZE), "lines"),
    legend.key.width = unit((SCALE_SIZE*.85), "lines"),
    legend.background = element_rect(fill=alpha('#c9c9c9', .5)),
    legend.title = element_text(size = SCALE_TEXT_SIZE),
    legend.text = element_text(size = (SCALE_TEXT_SIZE)),
    legend.justification = legend_position,
    legend.position = legend_position
  )
  
  #Get coordinates of top n cities  for given country
  top_cities <- cities_shp[cities_shp$iso3 == iso3,]
  top_cities <- top_cities[order(-top_cities$POP),][1:num_cities,]
  print(head(top_cities))
  coords <- sf::st_coordinates(top_cities)
  x_coords <- coords[, 'X']
  y_coords <- coords[, 'Y']
  city_text_size <- SCALE_SIZE * 1.67
  merged_shp[[id_col]] <- as.numeric(merged_shp[[id_col]])
  
  plots_list <- list()

  #Iterate through each year provided and create a map to add to the plot list
  for(i in seq(1, length(years))) {
    if(is.null(title)) {
      map_title <- sprintf("Year: %s", years[[i]])


    } else if(length(as.vector(title)) == 1) {
      map_title <- title
    } else {
      if (length(title) != length(years)) {
        stop("Number of titles does not match number of years--please\n
                 ensure you have an equal number of titles and years,\n
                 or a single title across all years")
      } else {
        map_title <-title[[i]]
      }

    }
    
    if(years[[1]] == 'ALL') {
      full_plot <- ggplot() +
        geom_sf(data=merged_shp, 
        aes_string(fill=id_col), 
        lwd=0)
    } else {
      full_plot <- ggplot() +
        geom_sf(data=merged_shp[merged_shp$year_id == years[[i]], ], 
        aes_string(fill=id_col), 
        lwd=0)
      
    }
    
    #Set custom scales depending on arguments provided to function
    if(transform_arg == 'log') {
      breaks <- c(max, max/10, max/100, min)
      labels <- c(max, max/10, max/100, min)
      labels <- as.character(labels)
      full_plot <- full_plot + scale_fill_gradientn(colors=color_scheme, 
      limits= c(min, max), 
      breaks=breaks, 
      labels=labels, 
      trans=transform_arg, 
      oob=squish)
    } else if (transform_arg == 'custom_cv') {
      breaks <- c(.25, .5, 1, 1.5, 2,3)
      labels <- c('', '', '1', '', '2', '3+')
      labels <- as.character(labels)
      full_plot <- full_plot + scale_fill_gradientn(colors=color_scheme, 
      limits= c(min, max), 
      breaks=breaks, 
      labels=labels, 
      oob=squish)
    } else if (transform_arg == 'custom_cv_0_1') {
      breaks <- c(.25, .5, .75, 1)
      labels <- c('.25', '.5', '.75', '1+')
      labels <- as.character(labels)
      full_plot <- full_plot + scale_fill_gradientn(colors=color_scheme, 
      limits= c(min, max), 
      breaks=breaks, 
      labels=labels, 
      oob=squish)
    } else {
      full_plot <- full_plot + scale_fill_gradientn(colors=color_scheme, 
      limits= c(min, max), 
      trans=transform_arg, 
      oob=squish)
    }
    
    full_plot <- full_plot +
      geom_sf(
        data = ad1_subset,
        fill = NA,
        color = '#575757',
        size= .2
      ) +
      geom_sf(
        data = ad0_inverse_subset,
        fill = "#d9d9d9",
        color = '#222222',
        size= .2
      ) +
      geom_sf(
        data = ad0_subset,
        fill = NA,
        color = '#000000',
        size= .3
      ) +
      geom_sf(data=top_cities, 
              size=1.0, 
              color='white', 
              pch=21, 
              fill='black', 
              alpha=.5) +
      geom_label_repel(aes(label=top_cities$CITY_NAME, 
                           x=x_coords, 
                           y=y_coords), 
                       seed = 123, 
                       color= NA, 
                       nudge_y = .4, 
                       hjust=0, 
                       nudge_x = -.1, 
                       size=city_text_size,
                       segment.size = .2,
                       fill='white',
                       label.padding = .1,
                       alpha=.5) +
      geom_label_repel(aes(label=top_cities$CITY_NAME, 
                           x=x_coords, 
                           y=y_coords), 
                       seed = 123, 
                       color='black', 
                       nudge_y = .4, 
                       nudge_x = -.1, 
                       hjust=0,
                       size=city_text_size,
                       segment.size = .2,
                       label.padding = .1,
                       label.size = NA,
                       fill=NA) +
      
      labs(fill=legend_title, title=map_title)
    
    if(expand_bounds) {
      full_plot <- full_plot + 
        coord_sf(xlim=expanded_bbox[1:2],
                 ylim=expanded_bbox[3:4], 
                 expand = FALSE)
    } else {
      full_plot <- full_plot + 
        coord_sf(xlim=bbox[1:2],
                 ylim=bbox[3:4],
                 expand=FALSE)
    }
    #Return list of ggplot objects
    plots_list[[i]] <- full_plot
    full_plot <- NULL
    
  }
  return(plots_list)
}