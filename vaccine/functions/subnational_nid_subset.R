subnational_nid_subset <- function(df) {
  df %>% 
    mutate(point = as.numeric(levels(point))[point]) %>% 
    rowwise() %>% 
    mutate(point = ifelse(svy_id %in% subnational_nids, point + 2, point)) %>% 
    ungroup() %>% 
    mutate(point = as.factor(point)) %>% 
    data.table()
}