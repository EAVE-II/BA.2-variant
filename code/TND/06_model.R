
##############################################################################
# Name of file: 06_model.R
# Data release (if applicable):
# Original author(s): Chris Robertson chris.robertson@nhs.scot; eleftheria.vasileiou@ed.ac.uk
# Original date: 20 Jan 2021
# Latest update author (if not using version control) - eleftheria.vasileiou@ed.ac.uk
# Latest update date (if not using version control) - 
# Latest update description (if not using version control)
# Type of script: Descriptive stats
# Written/run onL: R Studio SERVER
# Version of R that the script was most recently run on: R 3.6.1
# Description of content: regression model
# Approximate run time: Unknown
##############################################################################

#load library for model
library(mgcv)

Location <- "/conf/"  # Server
#Location <- "//isdsf00d03/"  # Desktop

setwd("/conf/EAVE/GPanalysis/progs/EV/TND_VOC")
project_path <- paste0(Location,"EAVE/GPanalysis/progs/EV/TND_VOC")

#load TND dataset
df.tnd <- readRDS("/conf/EAVE/GPanalysis/progs/EV/TND_VOC/data/df.tnd.rds")

#remove rows with non-Scottish HB
df.tnd <- df.tnd %>% filter(is.na(HB_residence) | HB_residence != "ENGLAND/WALES/NORTHERN IRELAND")

#order vacc vrbs
df.tnd$vs_type2 <- as.factor(df.tnd$vs_type2)
df.tnd$vs_type2 <- relevel(df.tnd$vs_type2, "uv")
levels(df.tnd$vs_type2)
#"uv"         "AZ_v1"      "AZ_v2_0:13" "AZ_v2_14+"  "Mo_v1"      "Mo_v2_0:13" "Mo_v2_14+"  "PB_v1"      "PB_v2_0:13" "PB_v2_14+" 

saveRDS(df.tnd, paste0(project_path, "/data/df.tnd.rds"))

##########################
#Delta plus variant AY.4.2
##########################

#subset here 
delta.plus <- df.tnd %>% filter(outcome %in% c("negative", "delta_plus"))

#convert vrbs into numeric class
delta.plus$days <- as.numeric(delta.plus$days)
delta.plus$age <- as.numeric(delta.plus$age)
delta.plus <- delta.plus %>% mutate(result = if_else(test_result=="POSITIVE", 1, 0))

#model1
m1 <- gam(result ~ s(days) + s(age) + vacc_status2 + subject_sex + simd + n_risk_gps,
                                                    family=binomial, data=delta.plus)

summary(m1)

#alternative way to subset data within the model directly
# m1 <- gam(result ~ s(days) + s(age) + vacc_status2 + subject_sex + simd2020_sc_quintile + n_risk_gps
#             , family=binomial, data=df, subset=outcome %in% c("negative", "delta_plus"))

#get ORs & 95%CIs
round(exp(cbind("Odds Ratio" = coef(m1), confint.default(m1, level=0.95), digits=3)), digits = 3)


#model 2
m2 <- gam(result ~ s(days) + s(age) + vs_type2 + subject_sex + simd + n_risk_gps,
                                                           family=binomial, data=delta.plus)
                                                                             
summary(m2)
                                                                             
                                                                             
#get ORs & 95%CIs
round(exp(cbind("Odds Ratio" = coef(m2), confint.default(m2, level=0.95), digits=3)), digits = 3)

saveRDS(delta.plus, paste0(project_path, "/data/delta.plus.rds"))

#remove redundant objects
rm(m1, m2, delta.plus)

##########################
#Delta variant
##########################


#subset here 
delta <- df.tnd %>% filter(outcome %in% c("negative", "delta"))

#convert vrbs into numeric class
delta$days <- as.numeric(delta$days)
delta$age <- as.numeric(delta$age)
delta <- delta %>% mutate(result = if_else(test_result=="POSITIVE", 1, 0))

#model1
m1 <- gam(result ~ s(days) + s(age) + vacc_status2 + subject_sex + simd + n_risk_gps,
                                                         family=binomial, data=delta)

summary(m1)


#get ORs & 95%CIs
round(exp(cbind("Odds Ratio" = coef(m1), confint.default(m1, level=0.95), digits=3)), digits = 3)


#model 2
m2 <- gam(result ~ s(days) + s(age) + vs_type2 + subject_sex + simd + n_risk_gps,
                                                     family=binomial, data=delta)

summary(m2)


#get ORs & 95%CIs
round(exp(cbind("Odds Ratio" = coef(m2), confint.default(m2, level=0.95), digits=3)), digits = 3)

saveRDS(delta, paste0(project_path, "/data/delta.rds"))

#remove redundant objects
rm(m1, m2, delta)

###############################
#Pos. tests but not sequenced 
###############################

#load TND dataset
df.tnd <- readRDS("/conf/EAVE/GPanalysis/progs/EV/TND_VOC/data/df.tnd.rds")

#subset here 
no.sequence <- df.tnd %>% filter(outcome %in% c("negative", "ns_other"))

#convert vrbs into numeric class
no.sequence$days <- as.numeric(no.sequence$days)
no.sequence$age <- as.numeric(no.sequence$age)
no.sequence <- no.sequence %>% mutate(result = if_else(test_result=="POSITIVE", 1, 0))

#model1
m1 <- gam(result ~ s(days) + s(age) + vacc_status2 + subject_sex + simd + n_risk_gps,
          family=binomial, data=no.sequence)

summary(m1)

#get ORs & 95%CIs
round(exp(cbind("Odds Ratio" = coef(m1), confint.default(m1, level=0.95), digits=3)), digits = 3)


#model2
m2 <- gam(result ~ s(days) + s(age) + vs_type2 + subject_sex + simd + n_risk_gps,
                                               family=binomial, data=no.sequence)

summary(m2)

#get ORs & 95%CIs
round(exp(cbind("Odds Ratio" = coef(m2), confint.default(m2, level=0.95), digits=3)), digits = 3)


saveRDS(no.sequence, paste0(project_path, "/data/no.sequence.rds"))

