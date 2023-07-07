# This script will demonstrate how to use the new aggregation of results tools
# to make plotting of adminstrative results and estimate vs data a bit
# more painless and allow for greater flexibility. You will notice that there is
# no new plotting functions in this scrip rather, the new functions aggregate
# the data in a format that allows you to use ggplot to make the plots most
# usefull for you and your team

.libPaths(c(.libPaths(), "FILEPATH"))
# installed via the following command on DATE against FILEPATH
# could not install via install.packages as no version was published for R 3.5.0
# devtools::install_github("tgouhier/biwavelet", upgrade_dependencies = FALSE)
library(biwavelet)
library(ggrepel)

## ## Set core_repo location and indicator group
## user              <- Sys.info()["user"]
## core_repo         <- sprintf("FILEPATH")
## remote            <- "origin"
## branch            <- "develop"
## indicator_group   <- "ort"
## indicator         <- "had_diarrhea"
## run_date          <- "RUN_DATE"

## indicator_group   <- 'training'
## indicator         <- 'tr_had_diarrhea'
## run_date          <- 'RUN_DATE'

sharedir <- paste('FILEPATH')
commondir <- paste('FILEPATH')
package_list <- c(
  t(read.csv(paste('FILEPATH'), header = FALSE)))

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0('FILEPATH'))
mbg_setup(package_list = package_list, repos = core_repo)
model_dir <- get_model_output_dir(indicator_group, indicator, run_date)
region <- get_output_regions(model_dir) ## this is ALL regions in model_dir

# this will tell the script later whether we want to make time series plots and
# maps of different administrative level data. A requirement of this script and
# the tools we will be using is that fractional raking has been run because we
# use those files to build the data structures used in the plotting below.
admLVLS <- 0:2

# Grab the shapefile version associated with your model run or define a new one
shapefile_version <- paste0(model_dir, "/config.csv") %>%
  read.csv(stringsAsFactors = FALSE) %>%
  filter(V1 == "raking_shapefile_version") %>%
  pull(V2)

# the plotlist will hold all of your ggplot objects AND ONLY ggplot objects
# it can be a nested list but the convience function at the end of this script
# will expect a ggplot object to be at the endpoints of this list in order to
# convert to the file format of your choice.
plotList <- list()

# This is the first new function which acts as the new get model pv table
# function. It grabs the full input data for a run and links it to your models
# estimates The default is for this function to grab results from holdout 0,
# however, you may change that iun the function arguments. The results is the
# original input data now with a number of draw columns attached to the
# data.table object which represent the estimates for that location, time,
# age group, and in the future sex.
estDF <- get_data_estimate_df(
      indicator_group,
      indicator,
      run_date,
      region,
      shapefile_version = shapefile_version)

# In addition to the point level data and how it compares to our estimates we
# also often want to plot higher level adminstrative data. While technically
# no input data for true geospatial models is passed in as admin data we can
# combine external admin data estimates either from surveys or gbd via left
# left joins here.
admDF0 <- load_admin_draws(
  indicator_group, indicator, run_date, region, level=0)
admDF1 <- load_admin_draws(
  indicator_group, indicator, run_date, region, level=1)
admDF2 <- load_admin_draws(
  indicator_group, indicator, run_date, region, level=2)

# Pull in the shape files for the SF objects so we can make plotting maps
# a lot easier using geom_sf we are also going to grab the centroids for
# ggrepel later.
adm0SF <- read_sf(get_admin_shapefile(0, version = shapefile_version)) %>%
  cbind(st_coordinates(st_centroid(.)))
adm1SF <- read_sf(get_admin_shapefile(1, version = shapefile_version)) %>%
  cbind(st_coordinates(st_centroid(.)))
adm2SF <- read_sf(get_admin_shapefile(2, version = shapefile_version)) %>%
  cbind(st_coordinates(st_centroid(.)))

# The last object that we are going to pull in is the resulting inla model for a
# run using the new function `read_inla_model` and the we wil extract the
# marginal distributions for the priors and posteriors with the function
# `extract_hyperpar_prior_posteriors` at the end of the script we show how to
# use the results to make comparisons of the priors and posteriors.

hyperDF <- as.data.table(bind_rows(lapply(region, function(r){
  read_inla_model(indicator_group, indicator, run_date, r, holdout=0) %>%
    extract_hyperpar_prior_posteriors() %>%
    mutate(region=r)
})))


# now we can actually strart making plots. We will also use the convenience
# function summarize draws which take confidence intervals of draws as well as
# a summary stat, such as median or mean. This function is much quicker than
# etheir data.table row wise quantile or apply functions as it is optimized
# for matrix applications and has its internals written in C. Note that
# in addition to giving a single value for width we also may give a vector of
# confidence intervals to compute.

(plotList$dataCoverageAdmin0 <- estDF %>%
  # Summarize all of the draw level information to the 95% CI
  summarize_draws(.width=.95) %>%
  # calculate the observed probability from data
  mutate(pobs = .[[indicator]]/N) %>%
  # Only look at five year age periods
  filter((year%%5) == 0) %>%
  # group up to admin zero and year also retaining the region column
  group_by(ADM_CODE, year, region) %>%
  # Summarize the total observed and estimated information based on N and weight
  dplyr::summarize(
    pobs=sum(weight*N*pobs)/sum(weight*N),
    .value=sum(weight*N*.value)/sum(weight*N),
    .lwr=sum(weight*N*.lwr)/sum(weight*N),
    .upr=sum(weight*N*.upr)/sum(weight*N),
    N=sum(N*weight)) %>%
  # plot the coverage
  ggplot(aes(x=pobs, y=.value, ymin=.lwr, ymax=.upr, color=region)) +
  geom_point() +
  geom_errorbar(size=.3) +
  theme_classic() +
  geom_abline(color="red") +
  facet_wrap(~year) +
  labs(x="Data Observed", y="Data Estimate") +
  ggtitle("Data Observed Vs 95% Coverage of Data Estimates: ADM0 Aggregation"))


(plotList$dataCovExpObs <- estDF %>%
  # Take draws from the binomial draws
  mutate_at(
    starts_with("draw", vars=names(.)),
    function(x) rbinom(nrow(.), round(.$N), x)) %>%
  as.data.table %>%
  # sumnmarize the draws at multiple levels returns a long data object
  summarize_draws(.width=c(.25, .5, .75, .95)) %>%
  # calculate coverage
  mutate(covered=(.[[indicator]] >= .lwr) & (.[[indicator]] <= .upr)) %>%
  # group by the different CI widths and regions
  group_by(.width, region) %>%
  # calculate overall coverage across regions
  dplyr::summarize(coverage=mean(covered, na.rm=T)) %>%
  ggplot(aes(x=.width, y=coverage, color=region, group=region)) +
  geom_line() +
  geom_point() +
  theme_classic() +
  geom_abline(color="red") +
  labs(x="Expected Coverage", y="Observed Coverage", color="", group=""))

# Make a mosser like plot for level 1 locations within level 0 location
# should we make these for all locations?
# right now the level zero plots are grouped by

if(0 %in% admLVLS){
  plotList$adm0TS <- lapply(region, function(r){
    admDF0 %>%
      filter(region==r) %>%
      summarize_draws() %>%
      ggplot(aes(x=year, y=.value, ymin=.lwr, ymax=.upr)) +
      theme_classic() +
      geom_line() +
      geom_ribbon(alpha=.2) +
      labs(x="Year", y="Estimate") +
      facet_wrap(~ADM0_NAME) +
      ggtitle(paste0(r, " ADM0 Time Series of Stunting"))
  })
  names(plotList$adm0TS) <- region

  plotList$adm0MAP <- lapply(region, function(r){
    admDF0 %>%
      filter(region==r & year == 2017) %>%
      summarize_draws() %>%
      {right_join(adm0SF, .)} %>%
      ggplot() +
      geom_sf(aes(fill=.value)) +
      scale_fill_distiller(palette="RdYlBu") +
      geom_label_repel(aes(x=X,
                           y=Y,
                           label = ADM0_NAME),
                       point.padding = unit(0.01, "lines"),
                       box.padding = unit(1.5, "lines"),
                       min.segment.length = unit(0, "lines"),
                       segment.alpha = 0.5) +
      theme_void() +
      labs(fill=indicator)
  })
  names(plotList$adm0MAP) <- region
}

if(1 %in% admLVLS){
  adm0 <- na.omit(unique(admDF1$ADM0_NAME))

  plotList$adm1TS <- lapply(adm0, function(r){
    admDF1 %>%
      filter(ADM0_NAME==r) %>%
      summarize_draws() %>%
      ggplot(aes(x=year, y=.value, ymin=.lwr, ymax=.upr)) +
      theme_classic() +
      geom_line() +
      geom_ribbon(alpha=.2) +
      labs(x="Year", y="Estimate") +
      facet_wrap(~ADM1_NAME) +
      ggtitle(paste0(r, " ADM1 Time Series of Stunting"))
  })
  names(plotList$adm1TS) <- adm0

  plotList$adm1MAP <- lapply(adm0, function(r){
    admDF1 %>%
      filter(ADM0_NAME==r & year == 2017) %>%
      summarize_draws() %>%
      {right_join(adm1SF, .)} %>%
      ggplot() +
      geom_sf(aes(fill=.value)) +
      scale_fill_distiller(palette="RdYlBu") +
      geom_label_repel(aes(x=X,
                           y=Y,
                           label = ADM1_NAME),
                       point.padding = unit(0.01, "lines"),
                       box.padding = unit(1.5, "lines"),
                       min.segment.length = unit(0, "lines"),
                       segment.alpha = 0.5) +
      theme_void() +
      labs(fill=indicator)
  })
  names(plotList$adm1MAP) <- adm0
}

if(2 %in% admLVLS){
  adm0 <- na.omit(unique(admDF1$ADM0_NAME))
  adm1 <- na.omit(unique(select(admDF2, ADM0_NAME, ADM1_NAME)))

  plotList$adm2TS <- list()
  plotList$adm2MAP <- list()

  for(z in adm0){
    print(z)
    adm1z <- pull(filter(adm1, ADM0_NAME == z), ADM1_NAME)

    plotList$adm2TS[[z]] <- lapply(adm1z, function(r){
      admDF2 %>%
        filter(ADM1_NAME==r) %>%
        summarize_draws() %>%
        ggplot(aes(x=year, y=.value, ymin=.lwr, ymax=.upr)) +
        theme_classic() +
        geom_line() +
        geom_ribbon(alpha=.2) +
        labs(x="Year", y="Estimate") +
        facet_wrap(~ADM2_NAME) +
        ggtitle(paste0(r, " ADM2 Time Series of Stunting"))
    })
    names(plotList$adm2TS[[z]]) <- adm1z

    plotList$adm2MAP[[z]] <- lapply(adm1z, function(r){
      admDF2 %>%
        filter(ADM1_NAME==r & year == 2017) %>%
        summarize_draws() %>%
        {right_join(adm2SF, .)} %>%
        ggplot() +
        geom_sf(aes(fill=.value)) +
        scale_fill_distiller(palette="RdYlBu") +
        geom_label_repel(aes(x=X,
                             y=Y,
                             label = ADM2_NAME),
                         point.padding = unit(0.01, "lines"),
                         box.padding = unit(1.5, "lines"),
                         min.segment.length = unit(0, "lines"),
                         segment.alpha = 0.5) +
        theme_void() +
        labs(fill=indicator)
    })
    names(plotList$adm2MAP[[z]]) <- adm1z

  }
}

plotList$hyperpriors <- hyperDF %>%
  mutate(label=paste0(label, "\nPrior: ", prior)) %>%
  # scale the densities to themselves since only the relative values matter
  # If we were not to do this then we couldnt effectively use facet grid
  group_by(region, label) %>%
  mutate(y=y/max(y)) %>%
  # plot the posterior and prior
  ggplot(aes(x, y, color=region, linetype=Distribution)) +
  geom_line() +
  facet_grid(region~label, scales = "free") +
  theme_classic() +
  labs(y="") +
  theme(
    strip.background.y = element_blank(), strip.text.y = element_blank(),
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.y = element_blank(), axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle=90, hjust=1))


# We can save the plot list using this function which will save not only the
# RDS file but also all of the plots as images with an extension `ext`
# specified by the user. `save_image_format` will need to be set to `TRUE`

save_plot_list(
  plotList, indicator_group, indicator, run_date, save_image_format=T, ext="png")

## plotList <- readRDS(paste0("FILEPATH))
