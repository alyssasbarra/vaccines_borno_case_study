indicator_family,binomial
check_cov_pixelcount,FALSE
coefs_sum1,FALSE
constrain_children_0_inf,FALSE
countries_not_to_rake,ESH+GUF
countries_not_to_subnat_rake,PHL
ctry_re_prior,"list(prior = 'pc.prec', param = c(5, 0.05))"
datatag,_resampled_with_age_cohort_with_svy_year
gbd_date,DATE
gbm_bf,0.75
gbm_bounded_0_1,TRUE
gbm_lr,0.005
gbm_tc,3
ho_mb,10
ho_ts,500
holdout_rundate,
holdout_strategy,nid
individual_countries,FALSE
inla_cores,10
intercept_prior,0
interval_mo,12
intstrat,eb
jn,sprintf('')
keep_inla_files,FALSE
lat_col,latitude
load_other_holdouts,FALSE
long_col,longitude
makeholdouts,FALSE
mesh_s_max_edge,"c(0.5, 5)"
mesh_s_offset,"c(1, 5)"
mesh_t_knots,"c(1:20)"
modeling_shapefile_version,DATE
shapefile_version,DATE
n_ho_folds,5
n_stack_folds,5
nugget_prior,"list(prior = 'loggamma', param = c(1, 5e-4))"
other_weight,
pop_measure,a0004t
rake_transform,logit
raking_shapefile_version,DATE
Regions,"c('ETH','IDN','IND','NGA','COD','PAK')"
resume_broken,TRUE
samples,100
scale_gaussian_variance_N,TRUE
skipinla,FALSE
skiptoinla,FALSE
skiptoinla_from_rundate,
slots,16
ss_col,weighted_n
st_targ,0.8
stacked_fixed_effects,gam + lasso + gbm
stackers_in_transform_space,TRUE
subnational_raking,TRUE
summstats,"c('mean', 'cirange', 'upper', 'lower', 'cfb')"
target_type,greater
test,FALSE
test_pct,5
time_stamp,TRUE
transform,inverse-logit
use_child_country_fes,FALSE
use_inla_country_fes,FALSE
use_nid_res,FALSE
use_raw_covs,FALSE
use_s2_mesh,TRUE
use_share,FALSE
use_stacking_covs,TRUE
validation_metrics_only,FALSE
withtag,TRUE
year_list,c(2000:2019)
yearload,annual
yr_col,year
z_list,c(1:5)
zcol,age_bin
z_map_file,FILEPATH
zcol_ag,age_bin_agg
use_subnat_res,FALSE
use_geos_nodes,TRUE
spat_strat,qt
subnat_country_to_get,
temp_strat,prop
rho_prior,"list(prior = 'normal', param = c(0, 1/(2.58^2)))"
run_time,16:00:00:00
s2_mesh_params,"c(40, 500, 500)"
s2_mesh_params_int,"c(40, 500, 500)"
skip.inla,0
skip.stacking,0
memory,10
metric_space,current
queue,long.q
rake_countries,TRUE
fit_with_tmb,TRUE
fixed_effects,
fixed_effects_measures,'mean + total'
gbd_fixed_effects,
gbd_fixed_effects_age,"c(2, 3, 4, 5)"
gbd_fixed_effects_measures,
use_adm2_res,FALSE
pop_release,DATE
no_nugget_predict,FALSE
spde_prior_s,"list(type = 'pc', prior = list(sigma = c(5, 0.05), range = c(0.01, 0.05)))"
spde_prior_int,"list(type = 'pc', prior = list(sigma = c(5, 0.05), range = c(0.01, 0.05)))"
error_iid_prior, "list(prior = 'pc.prec', param = c(3, 0.05))"
t_sigma_prior,"list(prior = 'pc.prec', param = c(5, 0.05))"
z_sigma_prior,"list(prior = 'pc.prec', param = c(5, 0.05))"
x_sigma_prior,"list(prior = 'pc.prec', param = c(7, 0.05))"
trho_st_prior,"list(prior = 'normal', param = c(3.5, 1/(1.44)))"
trho_tz_prior,"list(prior = 'normal', param = c(2, 1/(1.44)))"
trho_me_prior,"list(prior = 'normal', param = c(2, 1/(1.44)))"
zrho_prior,"list(prior = 'normal', param = c(2, 1/(1.44)))"
interacting_gp_1_effects,"c('space','year','age')"
poly_ag,FALSE
drop_low_weight,FALSE
low_weight_multiplier,2
mesh_z_knots,c(1:5)
use_gp,TRUE
use_space_only_gp,FALSE
use_time_only_gmrf,FALSE
use_age_only_gmrf,TRUE
use_error_iid_re,FALSE
use_sz_gp,FALSE
use_tz_gp,FALSE
use_covs_df,TRUE
use_country_res,TRUE
use_cre_z_gp,FALSE
use_z_fes,FALSE
use_inla_nugget,FALSE
use_time_spline,FALSE