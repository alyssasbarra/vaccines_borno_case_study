#' @title Predict MBG from a TMB Model
#'
#' @description Project out to a full sample space defined by the sr argument
#' 
#' @author USERNAME
#' 
#' @param samples Number of draws to take
#' @param seed Seed to set for RNG
#' @param model_fit_object object output from function fit_mbg_tmb()
#' @param tmb_input_data input_data output of build_mbg_data_stack_tmb
#' @param fes a string of model fixed effects of form 'var1 + var2 + var3' or as string vector
#'  corresponding to column names in d
#' @param sr simple raster object
#' @param yl a vector of years for analysis (i.e. c(2001,2002,2003))
#' @param zl a vector of z for analysis (i.e. c(1,2,3,4,5)). Config item z_list. No z-dimension it should just be zero
#' @param covs_list named list of covariate bricks, each should be length(yl) long
#'
#' @return a cell_preds object
#'
#' @export
predict_gps_tmb <- function(samples,
                           seed             = NULL,
                           tmb_input_stack  = input_data,
                           model_fit_object = model_fit,
                           sr               = simple_raster,
                           yl               = year_list,
                           zl               = z_list,
                           int_gp_1_effs    = c('space', 'time'),
                           use_tz_gp = as.logical(use_tz_gp),
                           int_gp_2_effs    = '',
                           transform        = 'inverse-logit',
                           int_gp_1_effect = FALSE,
                           use_sz_gp = as.logical(use_sz_gp),
                           use_cre_z_gp = as.logical(use_cre_z_gp),
                           use_cre_t_gp = as.logical(use_cre_t_gp),
                           use_cre_tz_gp = as.logical(use_cre_tz_gp),
                           mesh_t = mesh_t) {
  
  # libs
  require(raster)
  require(sp)
  
  #
  sdrep     <- model_fit_object$sdrep
  
  # pull a few useful things from the input data stack
  cs_transform         <- tmb_input_stack$cs_df
  mesh_int                 <- tmb_input_stack$mesh_int
  mesh_s                 <- tmb_input_stack$mesh_s
  cntry_re_map         <- tmb_input_stack$cntry_re_map
  if(clamp_covs == TRUE) {
    clamper            <- tmb_input_stack$clamper
  } else {
    clamper            <- NULL
  }
  
  
  ## dummy raster
  cell_idx <- seegSDM:::notMissingIdx(sr)
  
  # set seed if it is requested
  if(!is.null(seed)) set.seed(seed)
  
  # vector of means
  mu    <- c(sdrep$par.fixed,sdrep$par.random)
  
  # simulate draws
  #if(any(use_gp,use_space_only_gp, use_time_only_gmrf, use_age_only_gmrf)){
  draws <- rmvnorm_prec(mu = mu , prec = sdrep$jointPrecision, n.sims = samples)
  
  ## separate out the draws
  parnames      <- c(names(sdrep$par.fixed), names(sdrep$par.random))
  
  if (use_gp == T){
    epsilon_int_1_draws <- draws[parnames=='Epsilon_int_1',]
  }
  if (use_tz_gp==TRUE){
    epsilon_tz_draws <- draws[parnames=='Epsilon_tz',]
  }  
  if (use_space_only_gp == T){
    epsilon_s_draws <- draws[parnames=='Epsilon_s',]
  }
  if (use_time_only_gmrf == T){
    epsilon_t_draws <- draws[parnames=='Epsilon_t',]
  }
  if (use_age_only_gmrf == T){
    epsilon_z_draws <- draws[parnames=='Epsilon_z',]
  }
  if (use_sz_gp == T){
    epsilon_sz_draws <- draws[parnames=='Epsilon_sz',]
  }
  
  

  
  if (use_cre_z_gp == T){
    epsilon_cre_z_draws <- draws[parnames=='Epsilon_cre_z',]
  }
  
  
  if (use_cre_t_gp == T){
    epsilon_cre_t_draws <- draws[parnames=='Epsilon_cre_t',]
  }
  
  
  if (use_cre_tz_gp == T){
    epsilon_cre_tz_draws <- draws[parnames=='Epsilon_cre_tz',]
  }

  
  
  # get coordinates of full projection space
  # Extract admin0 code
  f_orig <- data.table(cbind(xyFromCell(sr, seegSDM:::notMissingIdx(sr)), adm_code=as.vector(sr[seegSDM:::notMissingIdx(sr)])))
  f_orig$t <- f_orig$z <- 1 # set initial time and Z
  f_orig[,tmpord:=1:.N]
  
  # use the country code dt from input_data to map admin0 code to RE values 
  f_orig <- merge(f_orig,cntry_re_map[,c('adm_code','re_id'),with=FALSE],by='adm_code',all.x=TRUE)
  f_orig <- f_orig[order(tmpord)] # make 100% sure everything is correctly ordered after the merge. 
  f_orig[, re_id := re_id+1 ]  # to deal with indexing which started at 0 in the cpp                        
  f_orig$re_id[is.na(f_orig$re_id)] <- 0 # to deal with countries not in the data
  
  # add time periods and z periods as needed
  message('Creating full sample space index')
  grp <- setDT(expand.grid(1:length(yl), 1:length(zl)))
  setnames(grp,c('Var1','Var2'),c('t','z'))
  grp[,group := 1:.N]
  fullsamplespace <- data.table()
  for(g in 1:max(grp$group)){
    tmp <- f_orig
    tmp[,z  := grp$z[grp$group==g]]
    tmp[,t  := grp$t[grp$group==g]]
    tmp[,gp := g]
    fullsamplespace <- rbind(fullsamplespace,tmp)
  }
  fullsamplespace[,idx := 1:.N]
  
  
  #Create time/age mesh knots for the projection matrix grouping
  mesh_tz_knots <- grp[t %in% 1:20]$group
 # mesh_tz_knots <- 1:20
  
  ## add time periods and z periods as needed -- interaction #1 sample space
  message('Creating interacting gp 1 sample space index')
  if('year' %in% int_gp_1_effs) t <- length(yl) else if(int_gp_1_effect == TRUE) stop('interacting gp without year is not yet supported') else t <- 1
  if('age' %in% int_gp_1_effs)  z <- length(zl) else z <- 1 
  
  if(use_gp == TRUE & !('space' %in% int_gp_1_effs)) stop('interacting gp without space is not yet supported')
  
  grp <- setDT(expand.grid(1:t, 1:z))
  setnames(grp,c('Var1','Var2'),c('t','z'))
  grp[,group := 1:.N]
  intsamplespace <- data.table()
  for(g in 1:max(grp$group)){
    tmp <- f_orig
    tmp[,z  := grp$z[grp$group==g]]
    tmp[,t  := grp$t[grp$group==g]]
    tmp[,gp := g]
    intsamplespace <- rbind(intsamplespace,tmp)
  }
  intsamplespace[,idx := 1:.N]
  
  #Sample space for SZ gp
  message('Creating sz gp sample space index')
  z=length(zl)
  grp <- setDT(expand.grid(1:z))
  setnames(grp,'z')
  grp[,group := 1:.N]
  szsamplespace <- data.table()
  for(g in 1:max(grp$group)){
    tmp <- f_orig
    tmp[,z  := grp$z[grp$group==g]]
    tmp[,gp := g]
    szsamplespace <- rbind(szsamplespace,tmp)
  }
  szsamplespace[,idx := 1:.N]
  
  message('line 166')
  # get surface locs to project on to
  pcoords        <- cbind(x=fullsamplespace$x, y=fullsamplespace$y) ## used for cov raster extract
  intcoords        <- cbind(x=intsamplespace$x, y=intsamplespace$y) ## used for cov raster extract
  pcoords_s      <- cbind(x=fullsamplespace[gp == 1,]$x, y=fullsamplespace[gp == 1,]$y) ## used for cov raster extract
  szcoords        <- cbind(x=szsamplespace$x, y=szsamplespace$y) ## used for cov raster extract
  
  if(mesh_int$manifold == "S2"){
    gp_coords <- lonlat3D(pcoords[, 1], pcoords[, 2])
    gp_coords_s <- lonlat3D(pcoords_s[, 1], pcoords_s[, 2])
    gp_coords_int <- lonlat3D(intcoords[, 1], intcoords[, 2])
    gp_coords_sz <- lonlat3D(szcoords[, 1], szcoords[, 2])
  } else {
    gp_coords <- pcoords
    gp_coords_s <- pcoords_s
    gp_coords_int <- intcoords
    gp_coords_sz <- szcoords
  }
  
  ## define grouping across periods
  groups_periods <- fullsamplespace$gp
  int_periods    <- intsamplespace$gp
  sz_periods    <- szsamplespace$gp
  message('line 189')
  ## use inla helper functions to project the spatial effect.
  if(use_gp){
    if(length(zl) > 1 & 'age' %in% interacting_gp_1_effects){
      message('using mesh tz in projection matrix: ')
      print(mesh_tz_knots)
      mesh_tz <- build_time_mesh(periods = mesh_tz_knots)
      print(mesh_tz)
      A.pred_int_1 <- inla.spde.make.A(
        mesh  = tmb_input_stack$mesh_int,
        loc   = gp_coords_int,
        group = int_periods,   
        group.mesh = mesh_tz)
    }else{
      message('using mesh t in projection matrix: ')
      print(mesh_t)
      A.pred_int_1 <- inla.spde.make.A(
        mesh  = tmb_input_stack$mesh_int,
        loc   = gp_coords_int,
        group = int_periods,                                 
        group.mesh = mesh_t)
    }
  }
  
  #if(use_space_only_gp){
  # make a projection matrix from data to s mesh
  #  A.pred_s <- inla.spde.make.A(
  #    mesh  = mesh_int,
  #    loc   = gp_coords)
  #}
  
  ### values of GP ST surface at each cell (long by nperiods)
  # if we have multiple zs then do this by z since its possible to throw a SuiteSparse 'Problem too large' error here. 
  if(use_gp){
    if(length(zl) > 1 & 'age' %in% int_gp_1_effs){
      cell_int_1 <- list()
        for(zz in zl){
          cell_int_1[[zz]] <- as.matrix(A.pred_int_1[(which(intsamplespace$z==zz)),] %*% epsilon_int_1_draws)
        }
      cell_int_1 <- do.call('rbind',cell_int_1)
    } else{
      cell_int_1 <- as.matrix(A.pred_int_1 %*% epsilon_int_1_draws)
    }
  }
  
  if(use_sz_gp){  
    A.pred_sz <- inla.spde.make.A(
      mesh  = tmb_input_stack$mesh_int,
      loc   = gp_coords_sz,
      group = sz_periods)
    cell_sz <- as.matrix(A.pred_sz %*% epsilon_sz_draws)
  }
  

  
  ### values of GP S surface at each cell (long by nperiods)
  
  
  if(use_space_only_gp)  {
    A.pred_s <- inla.spde.make.A(
      mesh  = mesh_s,
      loc   = gp_coords_s)
    
    cell_s <- as.matrix(A.pred_s %*% epsilon_s_draws)
  }
  
  # if(use_gp)             pred_gp_int_1 <- plogis(as.matrix(cell_int_1))
  #  if(use_space_only_gp)  pred_gp_s <- plogis(as.matrix(cell_s))
  #  if(use_time_only_gmrf) pred_gp_t <- plogis(epsilon_t_draws)
  #  if(use_age_only_gmrf)  pred_gp_z <- plogis(epsilon_z_draws)
  
  preds <- list(if(use_gp)             list(as.matrix(cell_int_1)) else NULL,
                if(use_tz_gp) epsilon_tz_draws else NULL,
                if(use_space_only_gp)  list(as.matrix(cell_s)) else NULL,
                if(use_time_only_gmrf) epsilon_t_draws else NULL,
                if(use_age_only_gmrf)  epsilon_z_draws else NULL,
                if(use_sz_gp)  list(as.matrix(cell_sz)) else NULL,
                if(use_cre_z_gp)  epsilon_cre_z_draws else NULL,
                if(use_cre_t_gp)  epsilon_cre_t_draws else NULL,
                if(use_cre_tz_gp)  epsilon_cre_tz_draws else NULL)
  
  for(i in 1:length(preds)){
    try( for(i in 1:length(preds)) {
      if (is.null(preds[[i]])) preds[[i]] <- NULL
    }, silent = T)
  }
  
  names(preds) <-c(if(use_gp)             'pred_gp_int_1',
                   if(use_tz_gp) 'pred_tz_gp',
                   if(use_space_only_gp)  'pred_gp_s',
                   if(use_time_only_gmrf) 'pred_gp_t',
                   if(use_age_only_gmrf)  'pred_gp_z',
                   if(use_sz_gp)  'pred_gp_sz',
                   if(use_cre_z_gp)  'pred_gp_cre_z',
                   if(use_cre_t_gp)  'pred_gp_cre_t',
                   if(use_cre_tz_gp)  'pred_gp_cre_tz')
  
  message(names(preds))
  
  names(preds)[!is.na(names(preds))]
  # if there is more than one z, then return a list of length zl cell_preds
  pred_gp_list <- list()
  for (i in 1:length(names(preds)[!is.na(names(preds))])){
    if(names(preds)[i]== 'pred_gp_int_1'){
      if('age' %in% int_gp_1_effs) l_z <- length(zl) else l_z <- 1
        chunklength <- dim(preds[[i]][[1]])[1]/l_z
    }else if(names(preds)[i]== 'pred_gp_s'){
      chunklength <- dim(preds[[i]][[1]])[1]
    }else if(names(preds)[i]== 'pred_gp_sz'){
      chunklength <- dim(preds[[i]][[1]])[1]
    } else{
      chunklength <- dim(preds[[i]])[1]
    }
    pred_gp_list[[i]]<- list()
    if(names(preds)[i]== 'pred_gp_int_1'){
        for(z in 1:l_z) {
          pred_gp_list[[i]][[z]] <- preds[[i]][[1]][((z-1)*chunklength+1):(chunklength*z),1:samples]
          message(paste0('done with ', names(preds)[i], ' for agebin ', z))
          # pred_gp <- pred_gp_list[[i]]
        }
    } else {
      if(names(preds)[i]== 'pred_gp_s'|names(preds)[i]== 'pred_gp_sz'){
        pred_gp_list[[i]] <- preds[[i]][[1]][1:chunklength,1:samples]
      } else {
        pred_gp_list[[i]] <- preds[[i]][1:chunklength,1:samples]
      }
    }
  }
  names(pred_gp_list) <- names(preds)[!is.na(names(preds))]
  # return the predicted cell_pred object
  return(pred_gp_list)
  
}
