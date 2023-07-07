# Change source/polygon to factor and shorten survey names to fit on plot legend
input_clean <- function(df) {
  df %>%
    rowwise() %>% 
    mutate(source = ifelse(source == "MACRO_DHS", "DHS", source)) %>%
    mutate(source = ifelse(source == "MACRO_AIS", "AIS", source)) %>%
    mutate(source = ifelse(source == "UNICEF_MICS", "MICS", source)) %>%
    mutate(source = ifelse(source == "COUNTRY_SPECIFIC", "CS", source)) %>%
    mutate(source = ifelse(source == "WB_CWIQ", "CWIQ", source)) %>%
    mutate(source = ifelse(source == "WB_CWIQ", "CWIQ", source)) %>%
    mutate(source = ifelse(source == "WB_LSMS", "LSMS", source)) %>%
    mutate(source = ifelse(source == "WB_LSMS_ISA", "ISA", source)) %>%
    mutate(source = ifelse(source == "WB_PRIORITY_SURVEY", "PRI_S", source)) %>%
    mutate(source = ifelse(source == "ARAB_LEAGUE_PAPFAM", "PAPFAM", source)) %>%
    mutate(source = ifelse(source == "JHSPH_PERFORMANCE_MONITORING_ACCOUNTABILITY_SURVEY_PMA2020", "PMA", source)) %>% 
    mutate(source = ifelse(nchar(source) > 6, str_trunc(source, 6, ellipsis = ""), source)) %>% #truncate source if it is too long and not specified above
    ungroup() %>% 
    mutate(source = as.factor(source)) %>%
    mutate(point = as.factor(point)) %>%
    data.table()
}