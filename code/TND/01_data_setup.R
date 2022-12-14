
##############################################################################
# Name of file: 02_lab_data.R
# Data release (if applicable):
# Original author(s): Chris Robertson chris.robertson@nhs.scot; eleftheria.vasileiou@ed.ac.uk
# Original date: 20 Jan 2021
# Latest update author (if not using version control) - eleftheria.vasileiou@ed.ac.uk
# Latest update date (if not using version control) - 
# Latest update description (if not using version control)
# Type of script: Descriptive stats
# Written/run onL: R Studio SERVER
# Version of R that the script was most recently run on: R 3.6.1
# Description of content: prepares data for TND
# Approximate run time: Unknown
##############################################################################

# 0 Setup ####
#Libraries
library(plyr)
library(tidyverse)
library(lubridate)

Location = '/conf/'

setwd(paste0(Location, 'EAVE/GPanalysis/analyses/BA.2-variant'))

#####################################################################
#Get all lab data from ECOSS
#####################################################################

#Start date for study period
a_begin = as.Date("2021-11-01")

#all laboratory data from ecoss 
lab.data <- readRDS("/conf/EAVE/GPanalysis/data/CDW_full.rds")

#convert specimen date to class date 
lab.data$date_ecoss_specimen <- as.Date(lab.data$date_ecoss_specimen)

#Number of pcr, pos/neg tests before a_begin
Tests <- lab.data %>% dplyr::select(EAVE_LINKNO, date_ecoss_specimen, test_result) %>%
  arrange(EAVE_LINKNO, date_ecoss_specimen) %>%
  filter(!duplicated(paste(EAVE_LINKNO, date_ecoss_specimen))) %>%
  filter(date_ecoss_specimen < a_begin) %>%
  group_by(EAVE_LINKNO) %>%
  dplyr::summarise(n_tests=n(), n_pos = sum(test_result=="POSITIVE"), n_neg = sum(test_result=="NEGATIVE")) %>%
  ungroup() %>% as.data.frame()

lab.data <- lab.data %>% left_join(Tests, by="EAVE_LINKNO") %>% 
  mutate(n_tests = if_else(is.na(n_tests), 0L, n_tests),
         n_pos = if_else(is.na(n_pos), 0L, n_pos),
         n_neg = if_else(is.na(n_neg), 0L, n_neg))

#####################################################################
#Get the positives and select first positive result
#####################################################################

#Select all positive tests from patients with covid symptoms and tested at lighthouse labs/community
z_pos <- lab.data %>% filter(test_result=="POSITIVE") %>%
  arrange(EAVE_LINKNO, date_ecoss_specimen) %>%
  filter(!duplicated(paste(EAVE_LINKNO, date_ecoss_specimen))) %>%
  mutate(lab = if_else(test_result_record_source == "ECOSS", "nhs","lh")) %>%
  mutate(date_onset_of_symptoms = as.Date(date_onset_of_symptoms)) %>%
  mutate(flag_covid_symptomatic = if_else(!is.na(flag_covid_symptomatic) & flag_covid_symptomatic=="true", 1L, 0L))

#keep only records with date of symptoms & tests from lighthouse lab
z_pos <- z_pos %>% filter(flag_covid_symptomatic==1 & lab == "lh")

#Remove all rows/records of people with a positive test before a_begin
z_pos <- z_pos %>% mutate(prev_test = if_else(date_ecoss_specimen < a_begin, 1, 0))

test_before_begin <- z_pos %>% 
  filter(prev_test=="1") %>%
  select(EAVE_LINKNO) %>%
  filter(!duplicated(EAVE_LINKNO))

z_pos <- z_pos %>% anti_join(test_before_begin, by = "EAVE_LINKNO") 

#Select the first positive test
z_pos <- z_pos %>%
  arrange(EAVE_LINKNO, date_ecoss_specimen) %>%
  filter(!duplicated(EAVE_LINKNO))

#####################################################################
#Get the negatives and select one at random
#####################################################################

#set a specific seed no. so whenever you re-run the random index you will get the same records
set.seed(12345678)

#Create a random ID number
lab.data <- lab.data %>% 
  mutate(random_index = runif(nrow(lab.data)))

#Get the negative tests
z_neg <- lab.data %>% filter(test_result=="NEGATIVE")

#Select all negative tests from patients with covid symptoms and tested at lighthouse labs/community
z_neg <- z_neg %>% mutate(lab = if_else(test_result_record_source == "ECOSS", "nhs","lh")) %>%
  mutate(date_onset_of_symptoms = as.Date(date_onset_of_symptoms)) %>%
  mutate(flag_covid_symptomatic = if_else(!is.na(flag_covid_symptomatic) & flag_covid_symptomatic=="true", 1L, 0L))

#keep only records with date of symptoms & tests from lighthouse lab
z_neg <- z_neg %>% filter(flag_covid_symptomatic==1 & lab == "lh") 

#find EAVE_LINKNOs with a neg result and who also have a positive one
z_sel <- z_neg$EAVE_LINKNO %in% z_pos$EAVE_LINKNO

#remove any patients in the negative dataset that also have a positive test
z_neg <- z_neg %>% filter(!z_sel) 

#remove any negative tests before a_begin
z_neg <- z_neg %>% mutate(prev_test = if_else(date_ecoss_specimen < a_begin, 1, 0))

#Remove all rows/records of people with a negative test before a_begin,
# and select one at random.
test_before_begin <- z_neg %>% 
  filter(prev_test=="1") %>%
  select(EAVE_LINKNO) %>%
  filter(!duplicated(EAVE_LINKNO))

z_neg <- z_neg %>% anti_join(test_before_begin, by = "EAVE_LINKNO") %>%
  arrange(EAVE_LINKNO, random_index) %>% 
  filter(!duplicated(EAVE_LINKNO))

#drop redundant datasets and variables
z_neg <- z_neg %>% select(-prev_test)

#####################################################################
#Bind positives and negatives
#####################################################################

#Bind pos. and neg. tests together 
lab.tests <- bind_rows(z_pos,z_neg) %>% 
  dplyr::select(-random_index)


# Import vaccination data
source('/conf/EAVE/GPanalysis/progs/Data_Cleaning/00_Read_DV_Vaccinations.R')

#link lab data with vaccine data 
df <- left_join(lab.tests, Vaccinations, by="EAVE_LINKNO") #%>%
#  filter(flag_incon==0) # %>% 
#  dplyr::select(-vacc_type_2, -flag_incon)

#Number of vaccinated before June 8, 2021 (yes/no)
df <- df %>% mutate(vacc_prev = case_when(date_vacc_1 < a_begin ~ 1,
                                          TRUE ~ 0))

# Creates vector of time period labels for vaccine status
create_time_period_labels <- function(min_days, max_days, interval){
  
  num_intervals = floor((max_days - min_days)/interval)
  
  list = c()
  
  for (i in 1:(num_intervals)){
       
       item = paste0((i-1)*interval, ':', i*interval -1)
       
       list = c(list, item)
  }
  
  list = c(list, paste0(max_days, '+'))
  
}

#create days since test 
df <- df %>% mutate(days = as.numeric(date_ecoss_specimen - a_begin),
              days_vacc_1_test = as.numeric(date_ecoss_specimen - date_vacc_1),
              days_vacc_2_test = as.numeric(date_ecoss_specimen - date_vacc_2),
              days_vacc_3_test = as.numeric(date_ecoss_specimen - date_vacc_3))

df <- df %>% 
  mutate(vacc_1_status = as.character( cut(days_vacc_1_test, breaks = c(0, 27, Inf), 
              labels = c("v1_0:27", "v1_28+") )),
         vacc_2_status = as.character( cut(days_vacc_2_test, breaks =  c(seq(0, 105, 7), Inf), 
              labels = paste0('v2_', create_time_period_labels(0, 105, 7)) ) ),
         vacc_3_status = as.character( cut(days_vacc_3_test, breaks =  c(seq(0, 105, 7), Inf), 
              labels = paste0('v3_', create_time_period_labels(0, 105, 7)) )) ) %>%
  mutate(vacc_status = case_when( is.na(date_vacc_1) | days_vacc_1_test < 0 ~ 'uv',
                        is.na(date_vacc_2) | days_vacc_2_test < 0 ~ vacc_1_status,
                        is.na(date_vacc_3) | days_vacc_3_test < 0 ~ vacc_2_status,
                        TRUE ~ vacc_3_status) ) %>%
  select(-vacc_1_status, -vacc_2_status, -vacc_2_status)

#Add vaccine type 
# df <- df %>% mutate(vs_type = case_when(vacc_status == 'uv' ~ vacc_status,
#                                         TRUE ~ paste(vacc_type, vacc_status, sep="_")) )
                    
df <- df %>% 
  mutate(vacc_2_status = as.character( cut(days_vacc_2_test, breaks =  c(seq(0, 14, 14), Inf), 
                                           labels = paste0('v2_', create_time_period_labels(0, 14, 14)) ) ),
         vacc_3_status = as.character( cut(days_vacc_3_test, breaks =  c(seq(0, 14, 14), Inf), 
                                           labels = paste0('v3_', create_time_period_labels(0, 14, 14)) )) ) %>%
  mutate(vacc_status_2 = case_when( is.na(date_vacc_1) | days_vacc_1_test < 0 ~ 'uv',
                                  is.na(date_vacc_2) | days_vacc_2_test < 0 ~ 'v1',
                                  is.na(date_vacc_3) | days_vacc_3_test < 0 ~ vacc_2_status,
                                  TRUE ~ vacc_3_status) ) %>%
  select(-vacc_2_status, -vacc_3_status)

#Add vaccine type 
# df <- df %>% mutate(vs_type_2 = case_when(vacc_status_2 == 'uv' ~ vacc_status_2,
#                                         TRUE ~ paste(vacc_type, vacc_status_2, sep="_")) )

# Make vaccine status a factor and relevel so baseline is uv
df <- mutate(df,
                 vacc_status = as.factor(vacc_status),
                 vacc_status_2 = as.factor(vacc_status_2)) %>%
  mutate(vacc_status = fct_relevel(vacc_status, "uv"),
         vacc_status = fct_relevel(vacc_status, "uv"))



EAVE_cohort <- readRDS(paste0(Location,"EAVE/GPanalysis/outputs/temp/Cohort_Demog_Endpoints_Dates2021-07-28.rds")) %>%
  filter(!duplicated(EAVE_LINKNO)) %>%
  #remove all who have died before the beginning
  filter(is.na(NRS.Date.Death) | (!is.na(NRS.Date.Death) & NRS.Date.Death > a_begin))

EAVE_Weights <- readRDS(paste0(Location,"EAVE/GPanalysis/outputs/temp/CR_Cohort_Weights.rds"))

EAVE_cohort  <- EAVE_cohort %>% left_join(EAVE_Weights, by="EAVE_LINKNO")
EAVE_cohort$eave_weight[is.na(EAVE_cohort$eave_weight)] <- mean(EAVE_cohort$eave_weight, na.rm=T)

#get QCOVID risk groups
qcovid_rg <- readRDS("/conf/EAVE/GPanalysis/progs/CR/Vaccine/output/temp/Qcovid_all.rds") %>% 
  dplyr::select(-(Sex:ur6_2016_name), -Q_BMI) %>%
  filter(!duplicated(EAVE_LINKNO))


#get BP and smoking groups
eave_rg_bpsmoke <- readRDS(paste0(Location,"EAVE/GPanalysis/outputs/temp/CR_Cohort_RG_EAVE_BP_Smoke.rds")) %>%
  filter(!duplicated(EAVE_LINKNO))%>% 
  dplyr::select(EAVE_LINKNO, EAVE_Smoking_Status_Worst, EAVE_BP) %>% 
  dplyr::rename(EAVE_Smoke = EAVE_Smoking_Status_Worst)

#add in household information
Cohort_Household <- readRDS("/conf/EAVE/GPanalysis/outputs/temp/Cohort_Household.rds") %>%
  mutate(n_hh_gp = cut(n_hh, breaks=c(0,1,2,5,10,30,100,max(n_hh)),
                       labels=c("1", "2", "3-5", "6-10", "11-30", "31-100", "101+")))%>% 
  mutate(ave_hh_age=if_else(is.na(ave_hh_age), mean(ave_hh_age, na.rm=T), ave_hh_age))

Cohort_Household <- Cohort_Household %>% 
  select(EAVE_LINKNO, n_hh_gp, ave_hh_age, care_home_elderly)

#Link all datasets with demographics
cohort <- EAVE_cohort %>% 
  #Select EAVE_cohort variables
  dplyr::select(EAVE_LINKNO:simd2020_sc_quintile) %>% 
  #Link to RG and BPsmoke
  left_join(eave_rg_bpsmoke, by="EAVE_LINKNO") %>% 
  #Link QCOVID risk groups
  left_join(qcovid_rg %>% 
              select(EAVE_LINKNO:Q_DIAG_CKD_LEVEL, n_risk_gps, bmi_impute), by="EAVE_LINKNO") %>%
  #Link to household size
  left_join(Cohort_Household, by="EAVE_LINKNO")

#Add demographics to main TND dataset
df_tnd <- df %>% left_join(cohort, by="EAVE_LINKNO")

#remove adults <16years old
df_tnd <- df_tnd %>% filter(age >= 16)


# Variant data
wgs <- readRDS(paste0(Location,"EAVE/GPanalysis/data/WGS_latest.rds")) %>% 
  mutate_at(c("Collection_Date","Sequencing_Date","Alignment_Date"), ~ as.Date(. , format="%d/%m/%Y")) %>%
  dplyr::rename(specimen_date = Collection_Date) %>%
  arrange(EAVE_LINKNO, specimen_date) %>% 
  filter(!duplicated(EAVE_LINKNO))

# S-gene status
sgene <- readRDS(paste0(Location,"/EAVE/GPanalysis/data/omicron_ctvals.rds"))

#create variant column
wgs <- wgs %>% mutate( variant = case_when(is.na(VariantofInterest) ~ "not_sequenced",
                                           VariantofInterest=="VOC-21APR-02" ~ "delta",
                                           TRUE ~ "other"))   %>%
      mutate( variant = ifelse(VariantShorthand %in% c('B.1.1.529', 'BA.2'), VariantShorthand, variant))

#Add VOCs in the TND study 
df_tnd <- df_tnd %>% left_join(wgs %>% select(EAVE_LINKNO, specimen_date, variant ), 
                            by=c( "EAVE_LINKNO" = "EAVE_LINKNO",
                              "date_ecoss_specimen" = "specimen_date") )

df_tnd <- df_tnd %>% mutate(variant = ifelse(is.na(variant), "not_sequenced", variant),
                            result = if_else(test_result=="POSITIVE", 1, 0)) %>%
                     mutate(outcome = case_when(result==0 ~ "negative",
                                                result==1 ~ variant,
                                                TRUE~"unknown"))

df_tnd <- left_join(df_tnd, sgene %>% select(EAVE_LINKNO, test_id, sgene_classification))

#remove rows with non-Scottish HB
df_tnd <- df_tnd %>% filter(is.na(HB_residence) | HB_residence != "ENGLAND/WALES/NORTHERN IRELAND")

saveRDS(df_tnd, './data/df_tnd.rds')
