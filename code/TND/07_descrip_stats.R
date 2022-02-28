
##############################################################################
# Name of file: 07_descrip_stats.R
# Data release (if applicable):
# Original author(s): Chris Robertson chris.robertson@nhs.scot; eleftheria.vasileiou@ed.ac.uk
# Original date: 20 Jan 2021
# Latest update author (if not using version control) - eleftheria.vasileiou@ed.ac.uk
# Latest update date (if not using version control) - 
# Latest update description (if not using version control)
# Type of script: Descriptive stats
# Written/run onL: R Studio SERVER
# Version of R that the script was most recently run on: R 3.6.1
# Description of content: descriptive stats 
# Approximate run time: Unknown
##############################################################################

library(tidyverse)
library(finalfit)

Location <- "/conf/"  # Server
#Location <- "//isdsf00d03/"  # Desktop

setwd("/conf/EAVE/GPanalysis/progs/EV/TND_VOC")
project_path <- paste0(Location,"EAVE/GPanalysis/progs/EV/TND_VOC")

#load data
delta <- readRDS("/conf/EAVE/GPanalysis/progs/EV/TND_VOC/data/delta.rds")
delta.plus <- readRDS("/conf/EAVE/GPanalysis/progs/EV/TND_VOC/data/delta.plus.rds")
no.sequence <- readRDS("/conf/EAVE/GPanalysis/progs/EV/TND_VOC/data/no.sequence.rds")

#number of tests and pos. tests for Delta in each vaccine group
table(delta$vacc_status2, useNA = "always")
table(delta$outcome, delta$vacc_status2, useNA = "always")
table(delta$vs_type2, useNA = "always")
table(delta$outcome, delta$vs_type2, useNA = "always")

#number of tests and pos. tests for Delta Plus in each vaccine group
table(delta.plus$vacc_status2, useNA = "always")
table(delta.plus$outcome, delta.plus$vacc_status2, useNA = "always")
table(delta.plus$vs_type2, useNA = "always")
table(delta.plus$outcome, delta.plus$vs_type2, useNA = "always")

#number of tests and pos. tests but not sequenced in each vaccine group
table(delta.plus$vacc_status2, useNA = "always")
table(delta.plus$outcome, delta.plus$vacc_status2, useNA = "always")
table(delta.plus$vs_type2, useNA = "always")
table(delta.plus$outcome, delta.plus$vs_type2, useNA = "always")

#Population characteristics by test for Delta and Delta Plus datasets
dependent <- "outcome"
explanatory <- c("age", "subject_sex", "simd", "EAVE_Smoke", "EAVE_BP",
                 "Q_DIAG_AF", "Q_DIAG_ASTHMA","Q_DIAG_BLOOD_CANCER","Q_DIAG_CCF","Q_DIAG_CEREBRALPALSY",        
                 "Q_DIAG_CHD","Q_DIAG_CIRRHOSIS","Q_DIAG_CONGEN_HD","Q_DIAG_COPD","Q_DIAG_DEMENTIA",         
                 "Q_DIAG_DIABETES_1", "Q_DIAG_DIABETES_2","Q_DIAG_EPILEPSY","Q_DIAG_FRACTURE",           
                 "Q_DIAG_NEURO","Q_DIAG_PARKINSONS","Q_DIAG_PULM_HYPER", "Q_DIAG_PULM_RARE",            
                 "Q_DIAG_PVD","Q_DIAG_RA_SLE", "Q_DIAG_RESP_CANCER","Q_DIAG_SEV_MENT_ILL",               
                 "Q_DIAG_SICKLE_CELL", "Q_DIAG_STROKE","Q_DIAG_VTE","Q_HOME_CAT",             
                 "Q_LEARN_CAT", "Q_DIAG_CKD_LEVEL", "n_risk_gps","bmi_impute","HB_residence")


delta %>% summary_factorlist(dependent, explanatory, p=TRUE, add_dependent_label=TRUE)

delta.plus %>% summary_factorlist(dependent, explanatory, p=TRUE, add_dependent_label=TRUE)

no.sequence %>% summary_factorlist(dependent, explanatory, p=TRUE, add_dependent_label = TRUE)

