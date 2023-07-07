library(data.table)
Regions <- c('vax_ansa','vax_caeu','vax_crbn','vax_cssa','vax_ctam','vax_eaas','vax_essa','vax_name','vax_seas','vax_soas','vax_sssa','vax_trsa','vax_wssa')
all_covars <- rbindlist(lapply(Regions, function(r){
  rbindlist(lapply(c('bcg','mcv1','polio3','dpt3'), function(v){
    v1 <- gsub('3','',v)
  df<-fread(paste0('FILEPATH'))
  return(df)
}))
}))
all_covars <- unique(all_covars)
as.data.table(table(all_covars$covariate))
orig<-copy(all_covars)
all_covars[,release:=NA]

all_covars[covariate %like% 'bias_corr_vacc', release :='DATE']
all_covars[covariate == 'access2', release :='DATE']
all_covars[covariate %like% 'cruts', release :='DATE']
all_covars[covariate == 'distriverslakes', release :='DATE']
all_covars[covariate == 'edu_mean_stage2pl', release :='DATE']
all_covars[covariate == 'elevation', release :='DATE']
all_covars[covariate == 'evi_v6', release :='DATE']
all_covars[covariate == 'ghslurbanicity', release :='DATE']
all_covars[covariate == 'growingseason', release :='DATE']
all_covars[covariate == 'irrigation', release :='DATE']
all_covars[covariate == 'ldi_pc', release :=NA] #gbd
all_covars[covariate == 'lst_day_v6', release :='DATE']
all_covars[covariate == 'lst_diurnal_diff_v6', release :='DATE']
all_covars[covariate == 'lst_night_v6', release :='DATE']
all_covars[covariate == 'mswep', release :='DATE']
all_covars[covariate == 'mx_warterror_10years', release :=NA]
all_covars[covariate == 'nexndvi', release :='DATE']
all_covars[covariate == 'tcb_v6', release :='DATE']
all_covars[covariate == 'tcw_v6', release :='DATE']
all_covars[covariate == 'worldpop_raked', release :='DATE']

all_covars<-unique(all_covars)
all_covars



for(r in Regions){
  for(v in c('bcg','mcv1','polio3','dpt3')){
    v1 <- gsub('3','',v)
    df<-fread(paste0('FILEPATH'))
    for(c in unique(df$covariate)){
      df[covariate==c,release:=all_covars[covariate==c]$release]
    }
    fwrite(df, paste0('FILEPATH'))
  }
}
