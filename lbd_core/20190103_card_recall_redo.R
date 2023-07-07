


##############################################################################
##### Card v recall: exploratory analysis with new data
##### USERNAME
##### DATE
##############################################################################

#---------- load in data ----------#
pacman::p_load(data.table, magrittr)
processed_logs <- list.files("FILEPATH", full.names=TRUE)
df <- lapply(processed_logs, fread) %>% rbindlist(., fill=TRUE)


#---------- check data was brought in correctly ----------#
dim(df)
colnames(df)

#---------- restrict to age range we actually care about (1-5 year olds),
#----------             kids who for sure have DPT doses                   ----------#
unique(df$age_year) # already age range captured
df <- df[which(!is.na(df$age_year)),]
df <- df[which(!is.na(df$dpt_dose)),]
dim(df)

#---------- how many kids have dpt card? ----------# 
frac <- ftable(df$dpt_card) 
frac[2] / (frac[1]+frac[2]) # overall card fraction
# if you use the has_vacc_card, the card fraction would be [1] 0.7652887 for has card (seen OR not seen)
# if you use the has_vacc_card, the card fraction would be [1] 0.3458046 for has card (seen ONLY)
#frac <- ftable(df$has_vacc_card) 

ftable(df$dpt1)
miss1 <- length(which(is.na(df$dpt1)))
ftable(df$dpt2)
miss2 <- length(which(is.na(df$dpt2)))
ftable(df$dpt3)
miss3 <- length(which(is.na(df$dpt3)))



#--------------- subset data to just kids who have information on all three doses
df <- df[which(!is.na(df$dpt1)),]
df <- df[which(!is.na(df$dpt2)),]
df <- df[which(!is.na(df$dpt3)),]


df <- df[which(df$year_start>1999),]

df$age_cohort <- ifelse( df$age_year<1, 0, 
                         ifelse(df$age_year<2, 1,
                                ifelse(df$age_year<3,2,
                                       ifelse(df$age_year<4,3,
                                              ifelse(df$age_year<5,4,5))))       )

df <- df[which(as.numeric(df$age_cohort)<5),]
df <- df[which(as.numeric(df$age_cohort)>0),]


df$card <- df$dpt_card
df$denom <- rep(1, dim(df)[1])

df$card1 <- ifelse(df$dpt_dose_from_card>=1,1,0)
df$card1_denom <- ifelse(!is.na(df$dpt_dose_from_card),1,0)

df$recall1 <- ifelse(df$dpt_dose_from_recall>=1,1,0)
df$recall1_denom <- ifelse(!is.na(df$dpt_dose_from_recall),1,0)

df$card1miss <- ifelse(is.na(df$card1), 999,df$card1)
df$recall1miss <- ifelse(is.na(df$recall1), 999,df$recall1)
df$combo1 <- df$card1miss+df$recall1miss
df$either1 <- ifelse(df$combo %in% c(1,2,1000),1,0)
df$either1_denom <- ifelse(df$combo %in% c(0,1,2,999,1000),1,0)


df$card3 <- ifelse(df$dpt_dose_from_card==3,1,0)
df$card3_denom <- ifelse(!is.na(df$dpt_dose_from_card),1,0)

df$recall3 <- ifelse(df$dpt_dose_from_recall>=3,1,0)
df$recall3_denom <- ifelse(!is.na(df$dpt_dose_from_recall),1,0)

df$card3miss <- ifelse(is.na(df$card3), 999,df$card3)
df$recall3miss <- ifelse(is.na(df$recall3), 999,df$recall3)
df$combo3 <- df$card3miss+df$recall3miss
df$either3 <- ifelse(df$combo3 %in% c(1,2,1000),1,0)
df$either3_denom <- ifelse(df$combo3 %in% c(0,1,2,999,1000),1,0)

df$card <- df$dpt_card
df$denom <-1

keep <- cbind(df$nid, df$age_cohort, df$card, df$denom, df$card1, df$card1_denom, df$either1, df$either1_denom, df$card3, df$card3_denom, df$either3, df$either3_denom)

surv_yr <- aggregate(keep, by=list(df$nid, df$age_cohort), FUN=sum, na.rm=T)
head(surv_yr)

colnames(surv_yr) <- c("nid","age","DROP","DROP1", "card","denom","card1","card1_denom","either1","either1_denom", "card3","card3_denom", "either3", "either3_denom")

surv_yr <- surv_yr[,c(1,2,5,6,7,8,9,10,11,12,13,14)]

surv_yr$card.fraction <- surv_yr$card / surv_yr$denom

surv_yr$card1_cov <- surv_yr$card1/surv_yr$card1_denom
surv_yr$either1_cov <- surv_yr$either1/surv_yr$either1_denom
surv_yr$card3_cov <- surv_yr$card3/surv_yr$card3_denom
surv_yr$either3_cov <- surv_yr$either3/surv_yr$either3_denom

a <- surv_yr$card3_cov
b <- surv_yr$either1_cov
c <-  surv_yr$card1_cov
d <- a*(b/c)
e <- surv_yr$either3_cov
surv_yr$ratio <- d/e

surv_yr$adjusted <- d
surv_yr$unadjusted <- e



library(ggplot2)

ggplot(surv_yr, aes(x=ratio, y=card.fraction, color=as.factor(age))) +
  geom_point(size=3, alpha=0.4, shape=16) + theme_minimal() +
  coord_equal() + xlab("WUENIC ratio") + ylab("Card fraction")

ggplot(surv_yr, aes(x=unadjusted, y=adjusted, color=as.factor(age))) + 
  geom_point(size=surv_yr$card.fraction*3, alpha=0.3, shape=16) + theme_minimal() +
  coord_equal() + xlab("crude IHME processed coverage") + ylab("coverage with WUENIC adjustment applied")










ftable(df$dpt_dose_from_card~df$dpt_dose_from_recall)








load('FILEPATH')




adm0 <-  readOGR(get_admin_shapefile(admin_level = 0))



ad0 <- sf::st_read(get_admin_shapefile(admin_level = 0, version = 'current'))

ad1 <- sf::st_read(get_admin_shapefile(admin_level = 1, version = 'current'))

ad2 <- sf::st_read(get_admin_shapefile(admin_level = 2, version = 'current'))



get_location_code_mapping()


location_metadata <- get_location_code_mapping(shapefile_version = 'current')
colnames(location_metadata) <- c( 'drop1', 'ADM0_NAME', 'drop3', "ISO3", 'drop4', 'drop5')


test0 <- merge(admin_0, ad0, by="ADM0_CODE")
final0 <- test0[,c(2,503,504)]
final0a <- merge(final0, location_metadata, by="ADM0_NAME")

test1 <- merge(admin_1, ad1, by="ADM1_CODE")
final1 <- test1[,c(2,503,504, 506)]
final1a <- merge(final1, location_metadata, by="ADM0_NAME")

test2 <- merge(admin_2, ad2, by="ADM2_CODE")
final2 <- test2[,c(2,503, 511, 513, 515)]
final2a <- merge(final2, location_metadata, by="ADM0_NAME", allow.cartesian=TRUE)

write.csv(final0a, 'FILEPATH')
write.csv(final1a, 'FILEPATH')
write.csv(final2a, 'FILEPATH')



ad0test <- merge(ad1, test1, by="ADM1_CODE")


test <- test1[,c(1,2,503,504, 506, 507)]







wuenic_ratio <- function(wf){
  
  ###### a ######
  a <- rep(NA, dim(wf)[1])
  for (k in 1:dim(wf)[1]){
    a[k] <- ifelse(wf$dpt_dose_from_card[k]==3, 1, 0)
  }
  
  a <- sum(c, na.rm=T) 
  
  ###### b ######
  b <- rep(NA, dim(wf)[1])
  b1 <- rep(NA, dim(wf)[1])
  b2 <- rep(NA, dim(wf)[1])
  b3 <- rep(NA, dim(wf)[1])
  
  for (k in 1:dim(wf)[1]){
    b1[k] <- ifelse(wf$dpt_dose_from_card[k]>=1,1,0)
    b2[k] <- ifelse(wf$dpt_dose_from_recall[k]>=1,1,0)
    b1[k] <- ifelse(is.na(wf$dpt_dose_from_card[k]),0,b1[k])
    b2[k] <- ifelse(is.na(wf$dpt_dose_from_recall[k]),0,b2[k])
    b3[k] <- b1[k] + b2[k]
    b[k] <- ifelse(b3[k]>=1,1,0)
  }
  
  b <- sum(b, na.rm=T)
  
  ###### c ######
  c<- rep(NA, dim(wf)[1])
  for (k in 1:dim(wf)[1]){
    c[k] <- ifelse(wf$dpt_dose_from_card[k]>=1, 1, 0)
  }
  
  c <- sum(c, na.rm=T) 
  
  
  ###### d ######
  d = a*(b/c)
  
  ###### e ######
  e <- rep(NA, dim(wf)[1])
  e1 <- rep(NA, dim(wf)[1])
  e2 <- rep(NA, dim(wf)[1])
  e3 <- rep(NA, dim(wf)[1])
  
  for (k in 1:dim(wf)[1]){
    e1[k] <- ifelse(wf$dpt_dose_from_card[k]==3,1,0)
    e2[k] <- ifelse(wf$dpt_dose_from_recall[k]>=3,1,0)
    e1[k] <- ifelse(is.na(wf$dpt_dose_from_card[k]),0,e1[k])
    e2[k] <- ifelse(is.na(wf$dpt_dose_from_recall[k]),0,e2[k])
    e3[k] <- e1[k] + e2[k]
    e[k] <- ifelse(e3[k]>=1,1,0)
  }
  
  e <- sum(e, na.rm=T)
  
  ###### ratio ######
  ratio <- d/e
  return(ratio)
  
}




DHS <- df[which(df$survey_name=="MACRO_DHS"),]
length(unique(DHS$nid))
dhs_ratio <- wuenic_ratio(wf=DHS)




MICS <- df[which(df$survey_name=="UNICEF_MICS"),]
length(unique(MICS$nid))
mics_ratio <- wuenic_ratio(wf=MICS)
