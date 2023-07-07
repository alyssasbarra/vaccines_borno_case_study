plot_hyperparameters_tmb<-function(indicator, indicator_group, run_date, age=0, holdout, save_file = NULL, regions=Regions) {
  
  # get regions
  outputdir <- paste0('FILEPATH')
  #config <- fread(paste0(outputdir, 'config.csv'))
  #regions <- eval(parse(text = config[V1 == "Regions", V2]))
  
  # extract prior and posterior distributions from INLA model objects
  message("Load models & extract priors and posteriors for hyper-parameters")
  # dist <- rbindlist(lapply(regions, function(r) {
  
  d <- data.table()
  dir <- sprintf('FILEPATH')
  files <- list.files(dir,'fitted_parameter')
  
  for (f in files) {
    tmp <- fread(sprintf('%s/%s',dir,f))
    r  <- gsub('.csv','',gsub('fitted_parameter_summary_table_','',f))
    tmp$region <- r
    d   <- rbind(d,tmp)
  }
  
  # make region names nice
  d[ grepl('cssa'     , region), region := 'Central SSA']
  d[ grepl('essa_sdn-COM'     , region), region := 'Eastern SSA']
  d[ grepl('sssa'     , region), region := 'Southern SSA']
  d[ grepl('wssa-CPV-STP-MRT'     , region), region := 'Western SSA']
  
  ## hyper parameters
  dd <- d[param_name %in% c(#'logtau','logkappa',
    'age_rho','year_rho', 'sex_rho', 
    'age_rho_me','year_rho_me', 'sex_rho_me',
    'age_rho_tz','year_rho_tz', 'sex_rho_tz',
    'age_rho_zx','year_rho_tx', 'sex_rho_tx', 'sex_rho_zx',
    'zrho_sz', 'xrho_sx',
    'age_variance','year_variance', 'sex_variance', 'tz_variance', 'tx_variance', 'zx_variance',
    'country_RE_variance','nugget_variance',
    'gam', 'lasso', 'xgboost','gbm', 'int', 'gbd_est',
    'range', 'nominal_variance',
    'range_me', 'nominal_variance_me',
    'range_sz', 'nominal_variance_sz',
    'range_sx', 'nominal_variance_sx',
    'ANC_01','site_RE_variance', 'site_fe',
    'country-specific-Z_rho','country-specific-T_rho', 'country-specific-X_rho', 'TZ_country-specific-T_rho','TZ_country-specific-Z_rho',
    'country-specific-TZ_variance','country-specific-Z_variance','country-specific-T_variance', 'country-specific-X_variance', 'observation_error_variance',
    'FE_z_level__2','FE_z_level__3','FE_z_level__4','FE_z_level__5','FE_z_level__6')]
  dd$param_name[dd$param_name=='nugget_variance']      <- 'Nugget Var'
  dd$param_name[dd$param_name=='country_RE_variance']  <- 'CRE Var'
  dd$param_name[dd$param_name=='age_rho']        <- 'Z-Rho STZ'
  dd$param_name[dd$param_name=='year_rho']       <- 'T-Rho ST+'
  dd$param_name[dd$param_name=='sex_rho']       <- 'X-Rho STZX'
  dd$param_name[dd$param_name=='age_rho_tz']        <- 'Z-Rho TZ'
  dd$param_name[dd$param_name=='year_rho_tz']       <- 'T-Rho TZ'
  dd$param_name[dd$param_name=='sex_rho_tz']       <- 'X-Rho TZX'
  dd$param_name[dd$param_name=='age_rho_me']        <- 'Z-Rho ME'
  dd$param_name[dd$param_name=='year_rho_me']       <- 'T-Rho ME'
  dd$param_name[dd$param_name=='sex_rho_me']       <- 'X-Rho ME'
  dd$param_name[dd$param_name=='age_rho_zx']        <- 'Z-Rho ZX'
  dd$param_name[dd$param_name=='year_rho_tx']       <- 'T-Rho TX'
  dd$param_name[dd$param_name=='sex_rho_tx']       <- 'X-Rho TX'
  dd$param_name[dd$param_name=='sex_rho_zx']       <- 'X-Rho ZX'
  dd$param_name[dd$param_name=='zrho_sz']       <- 'Z-Rho SZ'
  dd$param_name[dd$param_name=='xrho_sx']       <- 'Z-Rho SX'
  dd$param_name[dd$param_name=='logtau']         <- 'Theta 1'
  dd$param_name[dd$param_name=='logkappa']       <- 'Theta 2'
  dd$param_name[dd$param_name=='gam']            <- 'FE-GAM'
  dd$param_name[dd$param_name=='lasso']          <- 'FE-Lasso'
  dd$param_name[dd$param_name=='xgboost']        <- 'FE-BRT'
  dd$param_name[dd$param_name=='gbm']        <- 'FE-GBM'
  dd$param_name[dd$param_name=='age_variance']   <- 'Z Var'
  dd$param_name[dd$param_name=='year_variance']  <- 'T Var'
  dd$param_name[dd$param_name=='sex_variance']  <- 'X Var'
  dd$param_name[dd$param_name=='tz_variance']  <- 'TZ Var'
  dd$param_name[dd$param_name=='tx_variance']  <- 'TX Var'
  dd$param_name[dd$param_name=='zx_variance']  <- 'ZX Var'
  dd$param_name[dd$param_name=='nominal_variance'] <- 'ST+ Var'
  dd$param_name[dd$param_name=='range']          <- 'Range ST'
  dd$param_name[dd$param_name=='nominal_variance_sz'] <- 'SZ Var'
  dd$param_name[dd$param_name=='range_sz']          <- 'Range SZ'
  dd$param_name[dd$param_name=='nominal_variance_sx'] <- 'SX Var'
  dd$param_name[dd$param_name=='range_sx']          <- 'Range SX'
  dd$param_name[dd$param_name=='nominal_variance_me'] <- 'S Var'
  dd$param_name[dd$param_name=='range_me']          <- 'Range S'
  dd$param_name[dd$param_name=='int']          <- 'Intercept FE'
  dd$param_name[dd$param_name=='site_fe']          <- 'ANC FE'
  dd$param_name[dd$param_name=='gbd_est']          <- 'GBD FE'
  dd$param_name[dd$param_name=='site_RE_variance']          <- 'Site RE Var'
  dd$param_name[dd$param_name=='country-specific-Z_rho']          <- 'Rho CRE Z'
  dd$param_name[dd$param_name=='country-specific-T_rho']          <- 'Rho CRE T'
  
  dd$param_name[dd$param_name=='TZ_country-specific-Z_rho']          <- 'Z-Rho CRE TZ'
  dd$param_name[dd$param_name=='TZ_country-specific-T_rho']          <- 'T-Rho CRE TZ'
  dd$param_name[dd$param_name=='country-specific-X_rho']          <- 'Rho CRE X'
  dd$param_name[dd$param_name=='country-specific-Z_variance']          <- 'CRE Z Var'
  dd$param_name[dd$param_name=='country-specific-T_variance']          <- 'CRE T Var'
  dd$param_name[dd$param_name=='country-specific-TZ_variance']          <- 'CRE TZ Var'
  dd$param_name[dd$param_name=='country-specific-X_variance']          <- 'CRE X Var'
  dd$param_name[dd$param_name=='observation_error_variance']          <- 'Observ RE Var'
  
  #Re-order levels for nicer figure layout
  dd$param_name         <- as.factor(dd$param_name)
  dd$param_type         <- as.factor(dd$param_type)
  
  dd[,type:='other']
  dd[grep('Range', param_name), type := 'Range']
  dd[grep('Var', param_name), type := 'Variance']
  dd[grep('Rho', param_name), type := 'Rho']
  dd[grep('FE', param_name), type := 'Fixed effects']
  
  dd$param_name <-gsub('Range ', '', dd$param_name)
  dd$param_name <-gsub('Var ', '', dd$param_name)
  dd$param_name <-gsub('Var', '', dd$param_name)
  dd$param_name <-gsub('Rho ', '', dd$param_name)
  dd$param_name <-gsub(' ME', '', dd$param_name)
  dd$param_name <-gsub(' FE', '', dd$param_name)
  dd$param_name <-gsub('FE-', '', dd$param_name)
  dd$param_name <-gsub('FE_z_level__', 'Age group', dd$param_name)
  
  # if(use_gp == FALSE){
  #dd$param_name <- factor(dd$param_name,levels(dd$param_name)[c(3,2,6,4,7,9,10,8,5,1
  # 2,5,3,6,8,9,7,1,4
  #9, 2, 4, 3, 5, 6, 8, 1, 7
  #                                                               )])
  # } 
  # make plots
  message("Plotting hyper-parameters")
  
  if (is.null(save_file)) save_file <- paste0(outputdir, "/diagnostic_plots/tmb_hyperparameters.pdf")
  pdf(save_file, width = 15, height = 8)
  
  gg<-ggplot(dd,aes(x=param_name,y=median,color=region)) + 
    ylab('Posterior fitted hyper-parameter with 95% UI') + 
    xlab('') + 
    labs(color='')+
    facet_wrap(~type,scales="free") +
    geom_pointrange(aes(ymin=lower,ymax=upper),position=position_dodge(width=.5)) + 
    theme_minimal() +
    theme(
      axis.text=element_text(size=12),
      axis.title=element_text(size=12,face="bold"),
      legend.text=element_text(size=12),
      strip.text.x = element_text(size = 12))+
    theme(axis.text.x = element_text(angle = 45))
  plot(gg)
  dev.off()
  
  return(dist)
}
