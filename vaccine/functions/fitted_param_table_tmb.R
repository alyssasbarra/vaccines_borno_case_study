#' @title Get fitted model parameters from TMB object
#' @description Make a nice table of fitted model parameters from a geostat TMB object
#' @author USERNAME
#'
#' @param model_fit fitted TMB object that comes out of fit_mbg_tmb()
#' @param exp_fixed_effects Boolean, should the fixed effects be exponentiated (if model was fit with logit link). Defaults to TRUE
#' @param transform_hyperparams Boolean, should hyperparmeters be transformed from fitted to more natural space. Defaults to TRUE
#' @param draws Integer, number of draws to use. Defaults to 1000
#' @param calculate_range Boolean, should we calculate and report range using kappa in draws space. Defaults to TRUE
#' @param calculate_nom_var  Boolean, should we calculate and report nominal variance using kappa and tau in draws space. Defaults to TRUE
#'
#' @return formatted data.table object
#' @export
fitted_param_table_tmb <- function(model_fit,
                                   input_data,
                                   exp_fixed_effects     = TRUE,
                                   transform_hyperparams = TRUE,
                                   draws                 = 1000,
                                   calculate_range       = TRUE,
                                   calculate_nom_var     = TRUE,
                                   cov_constraints = covariate_constraint_vectorize(config),
                                   apply_constraints     = FALSE) {
  
  # get draws of parameter values
  mu <- model_fit$sdrep$par.fixed
  pm <- model_fit$sdrep$jointPrecision[1:length(mu),1:length(mu)]
  set.seed(123)
  draws <- rmvnorm_prec(mu,pm,draws)
  
  # deal with given names of parameters
  fn <- names(model_fit$sdrep$par.fixed)
  fn[fn == 'alpha_j'] <- model_fit$fenames
  
  # Apply constraint transformations
  if(apply_constraints == T) {
    tmb_const <- tmb_cov_constraint(model_fit$fenames, cov_constraints)
    draws[names(model_fit$sdrep$par.fixed) == "alpha_j",] <-
      apply_constraints(tmb_const, draws[names(model_fit$sdrep$par.fixed) == "alpha_j",])
  }
  
  # Transform fixed effects
  if(exp_fixed_effects == TRUE){
    draws[which(fn %in% model_fit$fenames),] <- exp(draws[which(fn %in% model_fit$fenames),])
  }
  
  # Get the range parameter
  if(calculate_range == TRUE){
    # see equation 6.16 (pg 196) of Blangiardo and Cameletti Book 
    ranger <- sqrt(8) / exp(draws[which(fn == 'logkappa'),])
    if(any(ranger!=0)) {
      draws  <- rbind(draws, ranger)
      fn     <- c(fn, 'range')
    }
    
    ranger_me <- sqrt(8) / exp(draws[which(fn == 'logkappa_me'),])
    if(any(ranger_me!=0)) {
      draws  <- rbind(draws, ranger_me)
      fn     <- c(fn, 'range_me')
    }
    
    ranger_sz <- sqrt(8) / exp(draws[which(fn == 'logkappa_sz'),])
    if(any(ranger_sz!=0)) {
      draws  <- rbind(draws, ranger_sz)
      fn     <- c(fn, 'range_sz')
    }
    
    ranger_sx <- sqrt(8) / exp(draws[which(fn == 'logkappa_sx'),])
    if(any(ranger_sx!=0)) {
      draws  <- rbind(draws, ranger_sx)
      fn     <- c(fn, 'range_sx')
    }
  }
  
  # Get the nominal variance parameter
  if(calculate_nom_var == TRUE){
    # see equation 6.17 (pg 196) of Blangiardo and Cameletti Book 
    nominal_variance <- 1 / (4 * pi * (exp(draws[which(fn == 'logkappa'),]))^2 * (exp(draws[which(fn == 'logtau'),]))^2)
    if(any(nominal_variance!=0)) {
      draws <- rbind(draws, nominal_variance)
      fn    <- c(fn, 'nominal_variance')
    }
    
    nominal_variance_me <- 1 / (4 * pi * (exp(draws[which(fn == 'logkappa_me'),]))^2 * (exp(draws[which(fn == 'logtau_me'),]))^2)
    if(any(nominal_variance_me!=0)) {
      draws <- rbind(draws, nominal_variance_me)
      fn    <- c(fn, 'nominal_variance_me')
    }
    
    nominal_variance_sz <- 1 / (4 * pi * (exp(draws[which(fn == 'logkappa_sz'),]))^2 * (exp(draws[which(fn == 'logtau_sz'),]))^2)
    if(any(nominal_variance_sz!=0)) {
      draws <- rbind(draws, nominal_variance_sz)
      fn    <- c(fn, 'nominal_variance_sz')
    }
    
    nominal_variance_sx <- 1 / (4 * pi * (exp(draws[which(fn == 'logkappa_sx'),]))^2 * (exp(draws[which(fn == 'logtau_sx'),]))^2)
    if(any(nominal_variance_sx!=0)) {
      draws <- rbind(draws, nominal_variance_sx)
      fn    <- c(fn, 'nominal_variance_sx')
    }
  }
  
  # transform hyperparmeters
  if(transform_hyperparams == TRUE){
    draws[which(fn == 'logtau'),]             <- draws[which(fn == 'logtau'),]
    draws[which(fn == 'logkappa'),]           <- draws[which(fn == 'logkappa'),]
    draws[which(fn == 'logtau_me'),]             <- draws[which(fn == 'logtau_me'),]
    draws[which(fn == 'logkappa_me'),]           <- draws[which(fn == 'logkappa_me'),]
    draws[which(fn == 'logtau_sz'),]             <- draws[which(fn == 'logtau_sz'),]
    draws[which(fn == 'logkappa_sz'),]           <- draws[which(fn == 'logkappa_sz'),]
    draws[which(fn == 'logtau_sx'),]             <- draws[which(fn == 'logtau_sx'),]
    draws[which(fn == 'logkappa_sx'),]           <- draws[which(fn == 'logkappa_sx'),]
    draws[which(fn == 'log_nugget_prec'),]   <- 1/sqrt(exp(draws[which(fn == 'log_nugget_prec'),]))
    draws[which(fn == 'log_cre_prec'),]      <- 1/sqrt(exp(draws[which(fn == 'log_cre_prec'),]))
    draws[which(fn == 'log_nidre_sigma'),]    <- exp(draws[which(fn == 'log_nidre_sigma'),])
    draws[which(fn == 'log_pixelre_sigma'),]   <- exp(draws[which(fn == 'log_pixelre_sigma'),])^2
    draws[which(fn == 'log_cre_sigma'),]      <- exp(draws[which(fn == 'log_cre_sigma'),])^2
    draws[which(fn == 'log_site_iid_sigma'),]      <- exp(draws[which(fn == 'log_site_iid_sigma'),])^2
    draws[which(fn == 'log_error_iid_sigma'),]      <- exp(draws[which(fn == 'log_error_iid_sigma'),])^2
    draws[which(fn == 'log_tz_sigma'),]      <- exp(draws[which(fn == 'log_tz_sigma'),])^2
    
    
    draws[which(fn == 'log_cre_z_sigma'),]      <- exp(draws[which(fn == 'log_cre_z_sigma'),])^2
    draws[which(fn == 'log_cre_t_sigma'),]      <- exp(draws[which(fn == 'log_cre_t_sigma'),])^2
    draws[which(fn == 'log_cre_x_sigma'),]      <- exp(draws[which(fn == 'log_cre_x_sigma'),])^2
    draws[which(fn == 'log_cre_tz_sigma'),]      <- exp(draws[which(fn == 'log_cre_tz_sigma'),])^2
    
    
    fn[fn == 'logtau']           <- 'logtau'
    fn[fn == 'logkappa']         <- 'logkappa'
    fn[fn == 'logtau_me']           <- 'logtau_me'
    fn[fn == 'logkappa_me']         <- 'logkappa_me'
    fn[fn == 'log_nugget_prec'] <- 'nugget_SD'
    fn[fn == 'log_cre_prec']    <- 'country_RE_SD'
    fn[fn == 'log_nidre_sigma']  <- 'NID_RE_SD'
    fn[fn == 'log_pixelre_sigma'] <- 'nugget_variance'
    fn[fn == 'log_cre_sigma']    <- 'country_RE_variance'
    fn[fn == 'log_site_iid_sigma']    <- 'site_RE_variance'
    fn[fn == 'log_error_iid_sigma']    <- 'observation_error_variance'
    
    
    fn[fn == 'log_cre_z_sigma']    <- 'country-specific-Z_variance'
    fn[fn == 'log_cre_t_sigma']    <- 'country-specific-T_variance'
    
    fn[fn == 'log_cre_tz_sigma']    <- 'country-specific-TZ_variance'
    
    draws[which(fn == 'zrho'),] <- (exp( draws[which(fn == 'zrho'),] ) - 1) / (exp( draws[which(fn == 'zrho'),] ) + 1)
    draws[which(fn == 'trho'),] <- (exp( draws[which(fn == 'trho'),] ) - 1) / (exp( draws[which(fn == 'trho'),] ) + 1)
    
    
    draws[which(fn == 'zrho_me'),] <- (exp( draws[which(fn == 'zrho_me'),] ) - 1) / (exp( draws[which(fn == 'zrho_me'),] ) + 1)
    draws[which(fn == 'trho_me'),] <- (exp( draws[which(fn == 'trho_me'),] ) - 1) / (exp( draws[which(fn == 'trho_me'),] ) + 1)
    
    
    draws[which(fn == 'zrho_tz'),] <- (exp( draws[which(fn == 'zrho_tz'),] ) - 1) / (exp( draws[which(fn == 'zrho_tz'),] ) + 1)
    draws[which(fn == 'trho_tz'),] <- (exp( draws[which(fn == 'trho_tz'),] ) - 1) / (exp( draws[which(fn == 'trho_tz'),] ) + 1)
    
    draws[which(fn == 'zrho_sz'),] <- (exp( draws[which(fn == 'zrho_sz'),] ) - 1) / (exp( draws[which(fn == 'zrho_sz'),] ) + 1)
    draws[which(fn == 'xrho_sx'),] <- (exp( draws[which(fn == 'xrho_sx'),] ) - 1) / (exp( draws[which(fn == 'xrho_sx'),] ) + 1)
    
    
    
    draws[which(fn == 'zrho_cre_z'),] <- (exp( draws[which(fn == 'zrho_cre_z'),] ) - 1) / (exp( draws[which(fn == 'zrho_cre_z'),] ) + 1)
    
    draws[which(fn == 'trho_cre_t'),] <- (exp( draws[which(fn == 'trho_cre_t'),] ) - 1) / (exp( draws[which(fn == 'trho_cre_t'),] ) + 1)
    
    draws[which(fn == 'trho_cre_tz'),] <- (exp( draws[which(fn == 'trho_cre_tz'),] ) - 1) / (exp( draws[which(fn == 'trho_cre_tz'),] ) + 1)
    draws[which(fn == 'zrho_cre_tz'),] <- (exp( draws[which(fn == 'zrho_cre_tz'),] ) - 1) / (exp( draws[which(fn == 'zrho_cre_tz'),] ) + 1)
    
    fn[fn == 'zrho_cre_z'] <- 'country-specific-Z_rho'
    
    fn[fn == 'trho_cre_t'] <- 'country-specific-T_rho'
    
    fn[fn == 'trho_cre_tz'] <- 'TZ_country-specific-T_rho'
    fn[fn == 'zrho_cre_tz'] <- 'TZ_country-specific-Z_rho'
    
    fn[fn == 'zrho'] <- 'age_rho'
    fn[fn == 'trho'] <- 'year_rho'
    
    fn[fn == 'zrho_me'] <- 'age_rho_me'
    fn[fn == 'trho_me'] <- 'year_rho_me'
    fn[fn == 'xrho_me'] <- 'sex_rho_me'
    
    fn[fn == 'zrho_tz'] <- 'age_rho_tz'
    fn[fn == 'trho_tz'] <- 'year_rho_tz'
    fn[fn == 'xrho_int2'] <- 'sex_rho_int2'
    
    
    
    draws[which(fn == 'log_t_sigma'),] <- exp(draws[which(fn == 'log_t_sigma'),])^2
    draws[which(fn == 'log_z_sigma'),] <- exp(draws[which(fn == 'log_z_sigma'),])^2
    draws[which(fn == 'log_x_sigma'),] <- exp(draws[which(fn == 'log_x_sigma'),])^2
    fn[fn == 'log_z_sigma'] <- 'age_variance'
    fn[fn == 'log_t_sigma'] <- 'year_variance'
    
    draws[which(fn == 'log_tz_sigma'),] <- exp(draws[which(fn == 'log_tz_sigma'),])^2
    fn[fn == 'log_tz_sigma'] <- 'tz_variance'
    
    
    
    fn[fn == 'site_fe'] <- 'site_fe'
    
  }
  
  # summarize draws and clean up data table
  su <- data.table(t(apply(draws,1,quantile,c(0.025,0.500,0.975))))
  su[, fn := fn]
  colnames(su) <- c('lower','median','upper','param_name')
  su <- su[,c('param_name','median','lower','upper'),with=FALSE]


  # return the final table
  return(su)
}
