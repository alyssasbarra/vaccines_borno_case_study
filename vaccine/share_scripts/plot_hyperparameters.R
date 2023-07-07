plot_hyperparameters <- function(indicator,
                                 indicator_title,
                                 indicator_group,
                                 run_date,
                                 age,
                                 holdout, 
                                 plot_dir) {
  
  # get regions
  outputdir <- paste0('FILEPATH')
  config <- fread(paste0(outputdir, 'config.csv'))
  regions <- eval(parse(text = config[V1 == "Regions", V2]))
  
  # function to grab priors & posteriors from INLA object
  # draws heavily on INLA:::plot.inla() and related functions
  message("Loading posteriors...")
  get_posterior <- function(x) {
    
    all.hyper = INLA:::inla.all.hyper.postprocess(x$all.hyper)
    
    hyper = x$marginals.hyperpar
    hyper_names <- INLA:::inla.nameunfix(names(hyper))
    
    prior_list <- lapply(hyper_names, function(hn) {
      
      # Get the hyperparameter values
      hh <- hyper[[hn]]
      m = INLA:::inla.smarginal(hh)
      
      return(m)
    })
    
    names(prior_list) <- hyper_names
    return(prior_list)
  }

  reg_posterior_list <- lapply(regions, function(reg) {
  
    message(paste0(" ", reg))  
  
    # load inla object
    inla_obj <- readRDS(paste0("FILEPATH"))
    
    return(get_posterior(inla_obj))
  })
  
  names(reg_posterior_list) <- regions
    
  
  # Now grab priors from one region (all the same)
  message("Loading priors...")
  reg <- "cssa"
  inla_obj <- readRDS(paste0("FILEPATH"))
  
  get_prior <- function(x) {
    
    all.hyper = INLA:::inla.all.hyper.postprocess(x$all.hyper)
    
    hyper = x$marginals.hyperpar
    hyper_names <- INLA:::inla.nameunfix(names(hyper))
    
    prior_list <- lapply(hyper_names, function(hn) {
      
      if (grepl("Theta. for", hn)) the_range <- c(-5, 5)
      if (grepl("GroupRho for", hn)) the_range <- c(-0.999, 0.999)
      if (grepl("Group PACF. for", hn)) the_range <- c(-0.999, 0.999)
      if (grepl("Precision for", hn)) the_range <- c(0, 100)
      
      # Get the hyperparameter values
      id = unlist(strsplit(attr(hyper[[i]], "hyperid"), "\\|"))
      xy = (INLA:::inla.get.prior.xy(section = tolower(id[2]), hyperid = id[1],
                              all.hyper = all.hyper, range = the_range, intern = FALSE,
                              debug = debug))
      return(xy)
    })
    
    names(prior_list) <- hyper_names
    return(prior_list)
  }
  
  prior_list <- get_prior(inla_obj)
  
  # Plot
  message("Plotting...")
  plot_hyperparam <- function(hn) {
    
    # Grab posteriors by region
    reg_post_xy <- lapply(regions, function(rr) {
      reg_list <- reg_posterior_list[[rr]][[hn]]
      return(data.table(x = reg_list$x, y = reg_list$y, reg = rr))
    })
    
    df_post <- rbindlist(reg_post_xy)
    df_post[reg == "cssa", reg_title := "Central Sub-Saharan Africa"]
    df_post[reg == "essa", reg_title := "Eastern Sub-Saharan Africa"]
    df_post[reg == "sssa", reg_title := "Southern Sub-Saharan Africa"]
    df_post[reg == "wssa", reg_title := "Western Sub-Saharan Africa"]
    df_post[reg == "name", reg_title := "North Africa"]
    
    df_prior <- prior_list[[hn]]
    df_prior <- data.table(x = df_prior$x, y = df_prior$y)
    
    gg <- ggplot() +
      geom_line(data = df_post, aes(x=x, y=y, color = reg_title)) +
      geom_line(data = df_prior, aes(x=x, y=y), linetype = "dashed", color = "black") +
      theme_bw() +
      labs(x = "Value", y = "Density", title = hn, color = "Region", subtitle = indicator_title)
    
    return(gg)
    
  }
  
  for (h in names(prior_list)) {
    fn <- tolower(h)
    fn <- gsub(" |\\.", "_", fn)
    
    png(file = paste0(plot_dir, indicator, "_hyperparam_", fn, ".png"),
        width = 8,
        height = 3, 
        units = "in", 
        res = 200)
    
    gg <- plot_hyperparam(h)
    print(gg)
    dev.off()
  }

}

plot_hyperparameters(indicator = "dpt3_cov",
                     indicator_title = "DPT3 Coverage",
                     indicator_group = "vaccine", 
                     run_date = "RUN_DATE",
                     age = 0, holdout = 0,
                     plot_dir = "FILEPATH")


plot_hyperparameters(indicator = "dpt1_cond",
                     indicator_title = "P(doses = 1 | doses <= 1)",
                     indicator_group = "vaccine", 
                     run_date = "RUN_DATE",
                     age = 0, holdout = 0,
                     plot_dir = "FILEPATH")


plot_hyperparameters(indicator = "dpt2_cond",
                     indicator_title = "P(doses = 2 | doses <= 2)",
                     indicator_group = "vaccine", 
                     run_date = "RUN_DATE",
                     age = 0, holdout = 0,
                     plot_dir = "FILEPATH")
