require(INLA)
require(data.table)
require(maptools)
require(sp)
require(rgeos)
require(raster)
require(ggplot2)
require(RColorBrewer)
require(gridExtra)
require(scales)

plot_predicted_params <- function(indicator, 
                               indicator_group, 
                               run_date, 
                               holdout=0, 
                               shapefile_version='current', 
                               r 
) {
  
  
  # set directories
  modeldir <- paste('FILEPATH', sep='/')
  outputdir <- paste('FILEPATH', sep='/')
  
  # load GAUL codes and region mapping
  #config <- fread(paste0(outputdir, 'config.csv'))
  
#  source("FILEPATH")
  # load ad0 shapefile
  library(maptools)
  ad0 <- readShapePoly(get_admin_shapefile(0, version = shapefile_version))
  ad0 <- ad0[ad0@data$ADM0_CODE %in% get_adm0_codes(r),]

  try(load(paste0(modeldir, run_date, '_bin0_', r, '_0.RData')))
  input_data <- readRDS(paste0(outputdir, 'tmb_data_input_list_', r, '_holdout_0_agebin_0.RDS'))
  model_fit<-readRDS(paste0(outputdir, 'tmb_model_fit_pre_preds_', r, '_holdout_0_agebin_0.RDS'))
  
  
  if (class(mesh_t_knots) == "character") mesh_t_knots <- eval(parse(text=mesh_t_knots))
  # get GP preds ----------------------------------------------------------------------------------------------------------
  outputdir_chunks <- file.path(outputdir, 'prediction_chunks')
  dir.create(outputdir_chunks, showWarnings = FALSE)
  chunks <- rep(50, 1)
  
  pm_names    <- c(if(as.logical(use_gp))          'pred_gp_int_1' else NA,
                   if(as.logical(use_tz_gp))          'pred_gp_tz' else NA,
                   if(as.logical(use_space_only_gp))  'pred_gp_s'   else NA,
                   if(as.logical(use_time_only_gmrf)) 'pred_gp_t'   else NA,
                   if(as.logical(use_age_only_gmrf))  'pred_gp_z'   else NA,
                   if(as.logical(use_sz_gp))         'pred_gp_sz'  else NA,
                   if(as.logical(use_cre_z_gp))          'pred_gp_cre_z'  else NA,
                   if(as.logical(use_cre_t_gp))          'pred_gp_cre_t'  else NA,
                   if(as.logical(use_cre_tz_gp))          'pred_gp_cre_tz'  else NA)
  
  pm_names<- pm_names[!is.na(pm_names)]
  
  chunk<-1   
  for (chunk in 1:length(chunks)) {
    pm <- predict_gps_tmb(samples      = chunks[chunk],
                         seed             = NULL,
                         tmb_input_stack  = input_data,
                         model_fit_object = model_fit,
                         sr               = simple_raster,
                         yl               = year_list,
                         zl               = z_list,
                         int_gp_1_effs    = interacting_gp_1_effects,
                         use_tz_gp = as.logical(use_tz_gp),
                         transform        = 'inverse-logit',
                         int_gp_1_effect = as.logical(use_gp),
                         use_sz_gp = as.logical(use_sz_gp),
                         use_cre_z_gp = as.logical(use_cre_z_gp),
                         use_cre_t_gp = as.logical(use_cre_t_gp),
                         use_cre_tz_gp = as.logical(use_cre_tz_gp),
                         mesh_t = mesh_t)
    for (i in 1:length(pm_names)){
      if(pm_names[i] == 'pred_gp_int_1') {
        if(z_list != 0 & 'age' %in% interacting_gp_1_effects)  zs <- z_list else zs <-1
      } else {
        zs <-1
      }
        for(z in zs) {
          if (length(zs)>1){
              pm_izx <- pm[[i]][[z]] 
              } else {
            pm_izx<- pm[[i]]
          }
          
          
          pathaddin <- paste0('_bin',z, '_chunk',chunk, '_', r,'_',holdout, '_', pm_names[i])
          message(paste0('Saving gp draws for chunk ', chunk, ' for bin ', z,  ' for ', pm_names[i]))
          save(
            pm_izx,
            file     = paste0(outputdir_chunks, '/', indicator,'_gp_draws_eb',pathaddin,'.RData'),
            compress = TRUE
          )
          
          rm(pm_izx)
          gc()
        }
      
    }
    rm(pm)
    gc()
  }
  
  gps <- list()
  
  if(z_list != 0 & 'age' %in% interacting_gp_1_effects)  zs <- z_list else zs <-1
  for (i in 1:length(pm_names)){
    mean_ras <- list()
      #mean_ras <- list()
      for (z in zs) {
        if(pm_names[i] != 'pred_gp_int_1' & z > 1) next
        cell_pred<-NULL
        message('Reload and combine temporary gp chunk files')
        for (chunk in 1:length(chunks)) {
          pathaddin <- paste0('_bin',z, '_chunk',chunk, '_', r,'_',holdout, '_', pm_names[i])
          load(paste0(outputdir_chunks, '/', indicator,'_gp_draws_eb',pathaddin,'.RData'))
          
          # if(pm_names[i]=='pred_gp_s') pm_izx <- pm_izx[1:(nrow(pm_izx)/length(year_list))]
          if(pm_names[i] == 'pred_gp_int_1' & length(zs)==1) {
          cell_pred <- cbind(cell_pred, pm_izx[[1]])
          } else{
            cell_pred <- cbind(cell_pred, pm_izx)
          }
          rm(pm_izx)
          gc()
        }
        
        pathaddin <- paste0('_bin',z, '_', r,'_',holdout, '_', pm_names[i]) # new pathaddin
        
        # make a mean raster
        library(matrixStats)
        if(pm_names[i] == 'pred_gp_int_1') {
            if(length(zs) != 1) {
            mean_ras[[z]]  <- insertRaster(simple_raster,matrix(rowMeans(cell_pred),
                                                                ncol = max(period_map$period)))
          } else {
            mean_ras  <- insertRaster(simple_raster,matrix(rowMeans(cell_pred),
                                                           ncol = max(period_map$period)))
          }
          #sd_ras    <- insertRaster(simple_raster,matrix(  rowSds(cell_pred),ncol = max(period_map$period)))
        } else if(pm_names[i] == 'pred_gp_s'){
          mean_ras  <- insertRaster(simple_raster,matrix(rowMeans(cell_pred),ncol = 1))
          mean_ras  <- stack(mean_ras)
          #sd_ras    <- insertRaster(simple_raster,matrix(  rowSds(cell_pred),ncol = 1))
        } else if(pm_names[i] == 'pred_gp_sz'){
          mean_ras  <- insertRaster(simple_raster,matrix(rowMeans(cell_pred),ncol = length(z_list)))
          mean_ras  <- stack(mean_ras)
          #sd_ras    <- insertRaster(simple_raster,matrix(  rowSds(cell_pred),ncol = 1))
        }  else {
          mean_ras <-as.data.table(cbind(rowMeans(cell_pred), rowSds(cell_pred)))
          names(mean_ras) <- c('mean', 'sd')
          if(pm_names[i] == 'pred_gp_tz') {
              mean_ras[,year := rep(year_list, length(z_list))]
              mean_ras[,agebin := rep(z_list, each=length(year_list))]
          } else if(pm_names[i] == 'pred_gp_cre_z') {
            mean_ras[,country := rep(input_data$cntry_re_map$country, each=length(z_list))]
            mean_ras[,year:= 0]
            mean_ras[,agebin:= rep(z_list, length(input_data$cntry_re_map$country))]
          } else if(pm_names[i] == 'pred_gp_cre_t') {
            mean_ras[,country := rep(input_data$cntry_re_map$country, each=length(year_list))]
            mean_ras[,agebin:= 0]
            mean_ras[,year:= rep(year_list, length(input_data$cntry_re_map$country))]
          } else if(pm_names[i] == 'pred_gp_cre_tz') {
            mean_ras[,country := rep(input_data$cntry_re_map$country, each=length(year_list)*length(z_list))]
            mean_ras[,agebin:= rep(z_list, length(input_data$cntry_re_map$country), each=length(year_list))]
            mean_ras[,year:= rep(year_list, length(input_data$cntry_re_map$country)*(length(z_list)))]
          }
          
        }
        
        # save z specific objects
        if(pm_names[i] %in% c('pred_gp_int_1', 'pred_gp_s', 'pred_gp_sz') ){
          writeRaster(
            mean_ras[[z]],
            file      = paste0(outputdir, '/', indicator,'_prediction_gp',pathaddin),
            overwrite = TRUE
          )
        } else {
          save(
            mean_ras,
            file = (paste0(outputdir, '/', indicator,'_prediction_gp',pathaddin, ".RData")),
            compress = TRUE
          )
        }
        
        #Clean up
        rm(cell_pred)
        #rm(sd_ras)
        gc()
        
      }
    gps[[i]] <- mean_ras
    
    rm(mean_ras)
  }
  names(gps) <- pm_names
  
  #Remove temporary chunk files
  #message('Remove temporary chunk files')
  #for (z in zs) {
  #  for (chunk in 1:length(chunks)) {
  #    pathaddin <- paste0('_bin',z,'_chunk',chunk, '_', r,'_',holdout)
  #    file.remove(paste0(outputdir_chunks, '/', indicator,'_gp_draws_eb',pathaddin,'.RData'))
  #  }
  #}
  #file.remove(outputdir_chunks)
  
  
  
  
  
  #---------------------------------------------------------------------------------------------------------------------------
  
  gps <- lapply(names(gps), function(x){
    if(!(x %in% c( 'pred_gp_int_1', 'pred_gp_sz'))){
      if(!is.data.table(gps[[x]])) {
        gps[[x]] <- data.table(rasterToPoints(gps[[x]])) 
        names(gps[[x]])[3] <- 'mean'
      }
      gps[[x]][, diff := mean - mean(mean)]
    } else return(gps[[x]])
  })
  
  names(gps) <- pm_names
  
  
  #---------------------------------------------------------------------------------------------------------------------------
  # extract country outlines
  poly <- ad0
  poly <- suppressMessages(fortify(poly))
  
  colors <- c('#ffffe0','#ffe4ac','#ffc879','#ffa84c','#ff8725','#ff5c03','#f12861','#cc117d','#a60383','#800080')
  
  if(use_tz_gp==TRUE){
    data=gps[['pred_gp_tz']]
    
    data$agebin <- as.factor(data$agebin)
    
    
    p_int2 <- ggplot() +
      geom_line(data=data, aes(x=year, y=mean, color=agebin)) +
      geom_ribbon(data=data, aes(x=year, ymin=mean-sd, ymax=mean+sd, fill=agebin), alpha=0.1)+
      theme_bw() +
      labs(x='Year', y= 'Interacting GMRF 2 Pred', color='Age', fill='Age')
    
  }
  
  
  if(as.logical(use_cre_z_gp)==TRUE){
    data=gps[['pred_gp_cre_z']]
    
    
    p11 <- ggplot() +
      geom_line(data=data, aes(x=agebin, y=mean, color=country)) +
      geom_ribbon(data=data, aes(x=agebin, ymin=mean-sd, ymax=mean+sd, fill=country), alpha=0.1)+
      theme_bw() +
      facet_wrap(.~country)+
      labs(x='Agebin', y= 'Country-specific Z pred', color='Country', fill='Country')
    
  }
  
  
  if(as.logical(use_cre_t_gp)==TRUE){
    data=gps[['pred_gp_cre_t']]
    
    
    p11_1 <- ggplot() +
      geom_line(data=data, aes(x=year, y=mean, color=country)) +
      geom_ribbon(data=data, aes(x=year, ymin=mean-sd, ymax=mean+sd, fill=country), alpha=0.1)+
      facet_wrap(.~country)+
      theme_bw() +
      labs(x='Year', y= 'Country-specific T pred', color='Country', fill='Country')
    
  }
  
  
  if(as.logical(use_cre_tz_gp)==TRUE){
    data=gps[['pred_gp_cre_tz']]
    
    
    p11_2 <- ggplot() +
      geom_line(data=data, aes(x=year, y=mean, color=as.factor(agebin))) +
      geom_ribbon(data=data, aes(x=year, ymin=mean-sd, ymax=mean+sd, fill=as.factor(agebin)), alpha=0.1)+
      facet_wrap(.~country)+
      theme_bw() +
      labs(x='Year', y= 'Country-specific TZ pred', color='Agebin', fill='Agebin')
    
  }
  
  
  if(use_space_only_gp){
    data=gps[['pred_gp_s']]
    
    p1 <- ggplot() +
      geom_raster(data=data, aes(x=x, y=y, fill=mean))+
      #geom_point(data = fdata, aes(x = x, y = y), size = 0.75, alpha=0.3) +
      #geom_point(data = data, aes(x = longitude, y = latitude, color=100*pmin(0.25, hiv_prev_disagg/N)), size = 0.2) +
      #geom_path(data = poly, aes(x = long, y = lat), size = 0.05) +
      #scale_color_gradientn(colors = colors, values = color_values) +
      scale_fill_gradientn(colors = colors, name = 'Space GP Pred')+
      guides(fill = guide_colorbar(barheight = 7, barwidth = 0.5, nbin = 100)) +
   #   coord_equal(xlim = range(poly$long), ylim = range(poly$lat)) +
      labs( x = NULL, y = NULL, title = NULL) +
      theme_classic(base_size = 12) +
      theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
            legend.title = element_text(angle = 90, hjust = 0.5),
            plot.title = element_text(hjust = 0.5), #plot.margin = unit(c(-.13, -.13, -.13, -.13), "in"), 
            panel.background = element_blank()) +
      theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                            colour = "white"), 
            panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                            colour = "white"))+
      guides(fill = guide_colorbar(barwidth = 0.75, barheight = 5, nbin = 1000, title.position = 'left'))
  }
  
  
  if(use_time_only_gmrf){
    data=gps[['pred_gp_t']]
    data[,year:= 1:nrow(data) + 1999]
    
    p2 <- ggplot() +
      geom_line(data=data, aes(x=year, y=mean)) +
      geom_ribbon(data=data, aes(x=year, ymin=mean-sd, ymax=mean+sd), alpha=0.1)+
      theme_bw() +
      labs(x='Year', y= 'Time GMRF Pred')
    
  }
  
  if(use_age_only_gmrf){
    data=gps[['pred_gp_z']]
    data[,agebin:= 1:nrow(data)]
    
    #    data$agebin<-as.factor(data$agebin)
    #   levels(data$agebin)<-c('15-19',
    #                             '20-24',
    #                             '25-29',
    ##                             '30-34',
    #                             '35-39',
    #                             '40-44',
    #                             '45-49')
    
    p3 <- ggplot() +
      geom_line(data=data, aes(x=agebin, y=mean)) +
      geom_ribbon(data=data, aes(x=agebin, ymin=mean-sd, ymax=mean+sd), alpha=0.1)+
      theme_bw() +
      labs(x='Agebins', y= 'Age GMRF Pred')
  }
  
  
  
  if(as.logical(use_gp)){
    
    rasters<-gps[['pred_gp_int_1']]
    if(length(z_list) > 1 & 'age' %in% interacting_gp_1_effects) {
      rasters <- rbindlist(lapply(z_list, function(z) {
        rasters <- rbindlist(lapply(c(2000:2019), function(year){
          raster_year <- data.table(rasterToPoints(rasters[[z]][[year - 1999]]))
          raster_year[,year := year]
          raster_year[,agebin:=z]
          setnames(raster_year, c("long", 'lat', 'mean', 'year', 'agebin'))
          setkey(raster_year, agebin, long, lat)
        }))
      }))
      
      rasters$agebin<-as.factor(rasters$agebin)
      levels(rasters$agebin)<-c('Ages 9-11mo',
                                'Ages 1y',
                                'Ages 2y',
                                'Ages 3y',
                                'Ages 4y')
      
      #  p5 <- lapply(c(levels(rasters$agebin)[5]), function(z){
      p5 <- ggplot() +
        geom_raster(data=rasters[agebin==levels(rasters$agebin)[2]], aes(x=long, y=lat, fill=mean))+
        #geom_point(data = fdata, aes(x = x, y = y), size = 0.75, alpha=0.3) +
        #geom_point(data = data, aes(x = longitude, y = latitude, color=100*pmin(0.25, hiv_prev_disagg/N)), size = 0.2) +
        #     geom_path(data = poly, aes(x = long, y = lat, group = group), size = 0.05) +
        facet_wrap(.~year)+
        geom_sf(data=st_as_sf(ad0), fill=NA, size=0.1)+
        #scale_color_gradientn(colors = colors, values = color_values) +
        scale_fill_gradientn(colors = colors, name = 'GP Pred')+
        guides(fill = guide_colorbar(barheight = 7, barwidth = 0.5, nbin = 100)) +
        #      coord_equal(xlim = range(poly$long), ylim = range(poly$lat)) +
        labs(x = '', y = '', title = paste0(levels(rasters$agebin)[2])) +
        theme_bw() +
        theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())
      #   return(p5)
      # })
      
      
      p6<- ggplot(data=rasters, aes(x=diff, color=agebin)) + 
        geom_density() + 
        theme_bw()+
        facet_wrap(.~as.factor(year))+
        labs(x = 'De-meaned GP Pred')
      
    } else {
      rasters <- rbindlist(lapply(c(2000, 2005, 2010, 2019), function(year){
        raster_year <- data.table(rasterToPoints(rasters[[year - 1999]]))
        raster_year[,year := year]
        setnames(raster_year, c("long", 'lat', 'mean', 'year'))
        setkey(raster_year, long, lat)
      }))
      
      years<-c(2000, 2005, 2010, 2019)
      
      p5 <- list()
      for(y in 1:4) {
        
        p5[[y]] <- ggplot() +
          geom_raster(data=rasters[year==years[y]], aes(x=long, y=lat, fill=mean))+
          #geom_point(data = fdata, aes(x = x, y = y), size = 0.75, alpha=0.3) +
          #geom_point(data = data, aes(x = longitude, y = latitude, color=100*pmin(0.25, hiv_prev_disagg/N)), size = 0.2) +
          #     geom_path(data = poly, aes(x = long, y = lat, group = group), size = 0.05) +
          #scale_color_gradientn(colors = colors, values = color_values) +
          scale_fill_gradientn(colors = colors, name = 'GP Pred', limits=c(min(rasters$mean),max(rasters$mean)))+
          guides(fill = guide_colorbar(barheight = 7, barwidth = 0.5, nbin = 100)) +
          #      coord_equal(xlim = range(poly$long), ylim = range(poly$lat)) +
          labs(x = '', y = '', title = paste0('Space*Year: ', years[y])) +
          theme_bw() +
          theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())
      }
    }
  }
  
  if(as.logical(use_sz_gp)){
    
    rasters<-gps[['pred_gp_sz']]
    rasters <- rbindlist(lapply(1:length(z_list), function(z){
      raster_age <- data.table(rasterToPoints(rasters[[z]]))
      raster_age[,agebin := z]
      setnames(raster_age, c("long", 'lat', 'mean', 'agebin'))
      setkey(raster_age, long, lat)
    }))
    
    
    rasters$agebin<-as.factor(rasters$agebin)
    #  levels(rasters$agebin)<-c('Ages 15-19',
    ##                            'Ages 20-24',
    #                            'Ages 25-29',
    #                            'Ages 30-34',
    #                            'Ages 35-39',
    ##                            'Ages 40-44',
    #                            'Ages 45-49',
    #                            'Ages 50-54',
    #                            'Ages 55-59')
    p7 <- ggplot() +
      geom_raster(data=rasters, aes(x=long, y=lat, fill=mean))+
      facet_wrap(.~agebin)+
      scale_fill_gradientn(colors = colors, name = 'GP Pred', limits=c(min(rasters$mean),max(rasters$mean)))+
      guides(fill = guide_colorbar(barheight = 7, barwidth = 0.5, nbin = 100)) +
      #   coord_equal(xlim = range(poly$long), ylim = range(poly$lat)) +
      labs(x = '', y = '', title = 'Space*Age') +
      theme_bw() +
      theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())
  }
  
  
  
  
  if(!exists('p1')) p1 <- ggplot() + theme_void()
  if(!exists('p2')) p2 <- ggplot() + theme_void()
  if(!exists('p3')) p3 <- ggplot() + theme_void()
  # if(!exists('p4')) p4 <- ggplot() + theme_void()
  if(!exists('p5')) {
    p5 <-list()
    p5[[1]] <- ggplot() + theme_void()
    p5[[2]] <- ggplot() + theme_void()
    p5[[3]] <- ggplot() + theme_void()
    p5[[4]] <- ggplot() + theme_void()
  }
  if(!exists('p6')) p6 <- ggplot() + theme_void()
  if(!exists('p7')) p7 <- ggplot() + theme_void()
  # if(!exists('p8')) p8 <- ggplot() + theme_void()
  # if(!exists('p9')) p9 <- ggplot() + theme_void()
  # if(!exists('p10')) p10 <- ggplot() + theme_void()
  if(!exists('p_int2')) p_int2 <- ggplot() + theme_void()
  
  
  if(!exists('p11')) p11 <- ggplot() + theme_void()
  if(!exists('p11_1')) p11_1 <- ggplot() + theme_void()
  if(!exists('p11_2')) p11_2 <- ggplot() + theme_void()
  # if(!exists('p12')) p12 <- ggplot() + theme_void()
  file_name <- paste0(outputdir, 'diagnostic_plots/gp_preds_',r,'_',holdout, '.pdf') # new pathaddin)
  
  
  pdf(file_name, width = 10, height = 5)
  
  grid.arrange(p1, p2, p3, ncol=2)
  
  # grid.arrange(p5[[1]], p5[[2]], p5[[3]], p5[[4]], ncol=2)
  try(grid.arrange(p5, ncol=1))
  
  plot(p7)
  
  plot(p_int2)
  
  # grid.arrange(p9, p10, ncol=2)
  plot(p11)
  
  plot(p11_1)
  
  plot(p11_2)
  
  dev.off()
  
  
  
  return(paste0('Plots saved for region: ', r))
  
}


