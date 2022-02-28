##########################################################
# Name of file: 07a_s_gene_dropout_description.R
# Data release (if applicable):
# Original author(s): Chris Robertson chrisobertson@nhs.net
#                     Steven Kerrr steven.ker@ed.ac.uk
# Original date: 06 August 2020
# Latest update author (if not using version control) - Steven Kerrr steven.ker@ed.ac.uk
# Latest update date (if not using version control) - 
# Latest update description (if not using version control)
# Type of script: Descriptive stats
# Written/run onL: R Studio SERVER
# Version of R that the script was most recently run on: R 3.6.1
# Description of content: reads in the cohort and merges in the risk groups 
#                         selects out the positive cases only and fits a prediction model
#                         reads in new positives and forecasts 28 day deaths + Hospitalisations
# Approximate run time: Unknown
##########################################################

# 0 Setup ####
#Libraries
library(plyr)
library(tidyverse)
library(survival)
library(lubridate)

Location = '/conf/'

setwd(paste0(Location, 'EAVE/GPanalysis/analyses/AY4.2-variant'))

########################### 1 Load data  #######################################

# Variant data
wgs <- readRDS(paste0(Location,"EAVE/GPanalysis/data/WGS_latest.rds")) %>% 
  mutate_at(c("Collection_Date","Sequencing_Date","Alignment_Date"), ~ as.Date(. , format="%d/%m/%Y")) %>%
  rename(specimen_date = Collection_Date) %>%
  arrange(EAVE_LINKNO, specimen_date) %>% 
  filter(!duplicated(EAVE_LINKNO))

# Start date is the first AY4.2 case
a_begin <- filter(wgs, lineage == 'AY.4.2') %>%
           pull(specimen_date) %>%
           min()
#a_begin <- as.Date("2020-12-08")

EAVE_cohort <- readRDS(paste0(Location,"EAVE/GPanalysis/outputs/temp/Cohort_Demog_Endpoints_Times2021-10-21.rds")) %>%
               filter(!duplicated(EAVE_LINKNO)) %>%
               #remove all who have died before the beginning
               filter(is.na(NRS.Date.Death) | (!is.na(NRS.Date.Death) & NRS.Date.Death > a_begin))

EAVE_Weights <- readRDS(paste0(Location,"EAVE/GPanalysis/outputs/temp/CR_Cohort_Weights.rds"))

EAVE_cohort  <- EAVE_cohort %>% left_join(EAVE_Weights, by="EAVE_LINKNO")
EAVE_cohort$eave_weight[is.na(EAVE_cohort$eave_weight)] <- mean(EAVE_cohort$eave_weight, na.rm=T)

EAVE_cohort <- EAVE_cohort %>% mutate(death_covid = case_when(death_covid==1 & is.na(NRS.Date.Death) ~ 0,
                                                              TRUE ~ death_covid),
                                      icu_death = case_when(icu_death==1 & is.na(date_icu_death) ~ 0,
                                                            TRUE ~ icu_death ) )

# Get EAVE_BP and EAVE_Smoke risk groups
z <- readRDS(paste0(Location,"EAVE/GPanalysis/outputs/temp/CR_Cohort_RG_EAVE_BP_Smoke.rds")) %>%
     filter(!duplicated(EAVE_LINKNO)) %>%
     dplyr::select(EAVE_LINKNO, EAVE_Smoking_Status_Worst, EAVE_BP) %>% 
     dplyr::rename(EAVE_Smoke = EAVE_Smoking_Status_Worst)

# Get QCovid risk groups
qcovid_rg <- readRDS(paste0(Location, "EAVE/GPanalysis/progs/CR/Vaccine/output/temp/Qcovid_all.rds")) %>% 
  filter(!duplicated(EAVE_LINKNO)) %>%
  # Categorise BMI into groups
  mutate(bmi_cat = cut(bmi_impute, breaks = c(-Inf, 18.5,24.9,29.9,Inf),
                       labels=c("Underweight","Normal weight","Overweight","Obese")))

# Combine risk groups
rg <- qcovid_rg %>% 
      full_join(z, by="EAVE_LINKNO") %>%
      filter(!duplicated(EAVE_LINKNO)) 

# Previous Tests data
cdw_full  <- readRDS(paste0(Location,"EAVE/GPanalysis/data/CDW_deduped.rds")) %>%
             rename(specimen_date = SpecimenDate) %>%
             mutate(specimen_date = as_date(specimen_date)) %>%
             arrange(EAVE_LINKNO, specimen_date, desc(result)) %>%
             #get one test per person per day - preferentially positive test 
             filter(!duplicated(paste(EAVE_LINKNO, specimen_date)))  

z <- cdw_full %>% 
     filter(specimen_date < a_begin) %>% 
     group_by(EAVE_LINKNO) %>% dplyr::summarise(n_tests=n()) %>% 
     dplyr::select(EAVE_LINKNO, n_tests)

rg <- rg %>% 
      left_join(z, by="EAVE_LINKNO") %>%
      mutate(n_tests = if_else(is.na(n_tests), 0L, n_tests))

# Import vaccination data
source("/conf/EAVE/GPanalysis/progs/SRK/00_Read_GP_Vaccinations.R")
#source("/conf/EAVE/GPanalysis/progs/CR/00_Read_GP_Vaccinations.R")

Positive_Tests <- cdw_full %>% 
                  filter(result==1) %>% 
                  dplyr::select(EAVE_LINKNO, specimen_date, test_result_record_source) %>%
                  mutate(days = as.numeric(specimen_date - a_begin)) %>% 
                  mutate(test_before_begin = cut(days, breaks = c((min(days)-1), -28, -21, -14, -7, 0, max(days)),
                      labels=c("1+m", "4w","3w","2w","0-6d","post-start")))

# Covid hospitalisations
# covid_hospitalisations  <- EAVE_cohort %>% 
#   dplyr::select(EAVE_LINKNO, SpecimenDate, hosp_covid, date_hosp_covid, NRS.Date.Death) %>% 
#   filter(hosp_covid==1) %>% 
#   filter(date_hosp_covid > a_begin) %>% 
#   dplyr::rename(admission_date = date_hosp_covid) %>% 
#   mutate(admission_date = if_else(is.na(NRS.Date.Death) | !is.na(NRS.Date.Death)&(admission_date <= NRS.Date.Death), admission_date, NRS.Date.Death)) %>% 
#   dplyr::select(-hosp_covid, -NRS.Date.Death)

#update with more timely data
# z <- readRDS(paste0(Location,"EAVE/GPanalysis/data/covid_hospitalisations.rds"))
# 
# z1 <- z %>% filter(admission_date >= a_begin) %>% 
#   filter(!(EAVE_LINKNO %in% covid_hospitalisations$EAVE_LINKNO)) %>% # find id's not in current list
#   left_join(Positive_Tests, by="EAVE_LINKNO") %>% 
#   dplyr::select("EAVE_LINKNO","SpecimenDate","admission_date")
# 
# covid_hospitalisations <- bind_rows(covid_hospitalisations, z1) %>% 
#   filter(EAVE_LINKNO %in% EAVE_cohort$EAVE_LINKNO)
# 
# z <- covid_hospitalisations %>% group_by(admission_date) %>% dplyr::summarise(N=n())
# 
# z %>% ggplot(aes(x=admission_date, y=N)) + geom_point()

#Covid icu deaths
covid_icu_death <- EAVE_cohort %>%
  dplyr::select(EAVE_LINKNO, SpecimenDate, icu_death, date_icu_death, NRS.Date.Death) %>%
  filter(icu_death==1) %>%
  filter(date_icu_death > a_begin) %>%
  dplyr::rename(admission_date = date_icu_death) %>%
  mutate(admission_date = if_else(is.na(NRS.Date.Death) | !is.na(NRS.Date.Death)&(admission_date <= NRS.Date.Death), admission_date, NRS.Date.Death)) %>%
  dplyr::select(-icu_death, -NRS.Date.Death)

#Covid deaths
# covid_death <- EAVE_cohort %>% dplyr::select(EAVE_LINKNO, SpecimenDate, death_covid, NRS.Date.Death) %>% 
#   filter(death_covid==1) %>% 
#   filter(NRS.Date.Death > a_begin) %>% 
#   dplyr::rename(admission_date = NRS.Date.Death) %>% 
#   dplyr::select(-death_covid)

# All deaths
all_deaths  <- readRDS(paste0(Location,"EAVE/GPanalysis/data/all_deaths.rds"))

# Deaths with covid on death certificate
all_deaths_covid_dth_cert <- all_deaths %>%  
  mutate(across(UNDERLYING_CAUSE_OF_DEATH:CAUSE_OF_DEATH_CODE_9, ~if_else(. %in% c("U071","U072"), 1,0))) %>%
  rowwise() %>% 
  mutate(rowsum = sum(c_across(UNDERLYING_CAUSE_OF_DEATH:CAUSE_OF_DEATH_CODE_9))) %>% 
  mutate(covid_death_cert = if_else(rowsum>=1,1,0)) %>% 
  dplyr::select(EAVE_LINKNO, NRS.Date.Death, covid_death_cert)

# All hospitalisations
all_hospitalisations  <- readRDS(paste0(Location,"EAVE/GPanalysis/data/any_hospitalisation_post_01022020.rds"))


###################### 2 Create cohort ################################################

df_cohort <- EAVE_cohort %>% 
             dplyr::select(EAVE_LINKNO:ur6_2016, ageYear, eave_weight) %>% 
             left_join( select(rg, EAVE_LINKNO, DataZone, starts_with('Q'), EAVE_Smoke,
                               EAVE_BP, n_risk_gps, bmi_impute, bmi_cat, n_tests), by="EAVE_LINKNO") %>%
             mutate(age_grp = case_when(ageYear < 18 ~ "<18",
                                        ageYear < 65 ~"18-64", 
                                        ageYear < 80 ~"65-79",
                                        TRUE ~ "80+"),
                    ur6_2016 = replace_na(ur6_2016, "Unknown")) 
        

df_cohort <- select(Positive_Tests, EAVE_LINKNO, test_before_begin, days, specimen_date,
                    test_result_record_source) %>%
       # Keep latest positive test. Our analysis will be for covid hosp/death based on 
       # most recent positive test in cohort time period.
       arrange(desc(days)) %>%
       filter(!duplicated(EAVE_LINKNO)) %>%
       right_join(df_cohort) %>% 
       mutate(test_before_begin = as.character(test_before_begin)) %>% 
       mutate(test_before_begin = if_else(is.na(test_before_begin), "no pos test",test_before_begin)) %>%
       left_join(Vaccinations, by="EAVE_LINKNO") %>%
       mutate(flag_incon = if_else(is.na(flag_incon), 0,flag_incon)) %>%
       # Add most recent vaccination status
       mutate(vacc_type_comb = case_when(vacc_type_2  == 'AZ' ~ 'AZ_v2',
                                        vacc_type_2  == 'PB' ~ 'PB_v2',
                                        vacc_type  == 'AZ' ~ 'AZ_v1',
                                        vacc_type  == 'PB' ~ 'PB_v1',
                                        vacc_type == 'UNK' ~ NA_character_,
                                        flag_incon == 1 ~ NA_character_,
                                        TRUE ~ 'uv')) 

## Add in household information
# Household information (from Sept 2020)
Cohort_Household <- readRDS("/conf/EAVE/GPanalysis/outputs/temp/Cohort_Household.rds") %>%
  mutate(n_hh_gp = cut(n_hh, breaks=c(0,1,2,5,10,30,100,max(n_hh)),
                       labels=c("1", "2", "3-5", "6-10", "11-30", "31-100", "101+")))%>% 
  mutate(ave_hh_age=if_else(is.na(ave_hh_age), mean(ave_hh_age, na.rm=T), ave_hh_age) )

df_cohort <- df_cohort %>%
  left_join(select(Cohort_Household, EAVE_LINKNO, n_hh_gp, ave_hh_age)) %>%
  mutate(simd = as.character(simd)) %>%
  mutate(simd = case_when(simd == 'NA' ~ NA_character_,
                          TRUE ~ simd )) %>%
  mutate(simd = as.factor(simd)) %>%
  mutate(n_tests = as.factor(n_tests))

# Add in all deaths
df_cohort <- df_cohort %>%
              left_join(all_deaths_covid_dth_cert,  by="EAVE_LINKNO") 

z_ids <- c(Vaccinations$EAVE_LINKNO, all_deaths$EAVE_LINKNO, 
           filter(EAVE_cohort, tested==1)$EAVE_LINKNO, all_hospitalisations$EAVE_LINKNO) %>% 
         unique()

# Re-weighting taking into account contact with healthcare
z_N <- round(sum(df_cohort$eave_weight) )
z_k <- sum(df_cohort$EAVE_LINKNO %in% z_ids)
z_m <- round(sum(filter(df_cohort, (EAVE_LINKNO %in% z_ids))$eave_weight))

df_cohort <- df_cohort %>% 
  mutate(eave_weight = if_else(EAVE_LINKNO %in% z_ids, 1, eave_weight*(z_N - z_k)/(z_N - z_m)) )

########## 3 Create study datatset of all who tested positive  ########################

# Use all positive tests
df_pos <- Positive_Tests %>%
          filter(test_before_begin == 'post-start') %>%
          left_join(df_cohort) %>%
          filter(ageYear >= 18 & flag_incon == 0) %>%
          mutate(lab = if_else(test_result_record_source == "ECOSS", "nhs","lh")) 

# Add all hospitalisations. This is important because some will be in hospital at time of test
z <- df_pos %>% left_join( all_hospitalisations %>% 
                             dplyr::select(-validchi) %>% 
                             filter( !(!is.na(discharge_date) & discharge_date <= a_begin))) %>%
  mutate(admission_date = case_when( discharge_date <= specimen_date ~ as.Date(NA),
                                     TRUE ~ admission_date),
         emergency = case_when( discharge_date <= specimen_date ~ NA,
                                TRUE ~ emergency),
         discharge_date = case_when( discharge_date <= specimen_date ~ as.Date(NA),
                                     TRUE ~  discharge_date))

# Separate multiple admissions and single/zero admissions
z_id <- z %>% filter(duplicated(EAVE_LINKNO)) %>% pull(EAVE_LINKNO) %>% unique()

z_01 <- z %>% filter(!(EAVE_LINKNO %in% z_id)) # 0 or 1 admission
z_m <- z %>% filter((EAVE_LINKNO %in% z_id))  #multiple admissions

# Select first post specimen hospital admission for those with multiple admissions
z_m <- z_m %>%  
  arrange(EAVE_LINKNO, admission_date) %>% 
  filter(!duplicated(EAVE_LINKNO))  

df_pos <- bind_rows(z_01, z_m)

# End date is maximum of admission date and discharge date
a_end_pos <- max(c(max(df_pos$admission_date, na.rm=T), max(df_pos$discharge_date, na.rm=T))) 


# Add variables  
df_pos <- df_pos %>%        
  mutate(
    ageYear = ifelse(ageYear >= 100, 100, ageYear),
    days = as.numeric(specimen_date -min(specimen_date, na.rm = TRUE)),
    death = if_else(is.na(NRS.Date.Death), 0L, 1L),
    vs1 = case_when(is.na(date_vacc_1) | date_vacc_1 > specimen_date ~ "uv",
                    date_vacc_1 <= specimen_date &  date_vacc_1 > specimen_date - 28 ~ "v1_0:27",
                    TRUE ~ "v1_28+"),
    vs2 = case_when(is.na(date_vacc_2) | date_vacc_2 > specimen_date ~ "uv",
                    date_vacc_2 <= specimen_date &  date_vacc_2 > specimen_date - 28 ~ "v2_0:27",
                    TRUE ~ "v2_28+"),
    vs = if_else(vs2=="uv", vs1, vs2), 
    time_to_hosp = if_else(!is.na(admission_date), 
                           as.numeric(admission_date - specimen_date), 
                           as.numeric(a_end_pos - specimen_date)) ,
    time_to_death = if_else(!is.na(NRS.Date.Death), 
                            as.numeric(NRS.Date.Death - specimen_date), 
                            as.numeric(a_end_pos - specimen_date))      ) %>% 
  mutate(
    vs_type = case_when(vs %in% c("v1_0:27", "v1_28+") ~ paste0(vacc_type, '_', vs) ,
                        vs %in% c("v2_0:27", "v2_28+") ~ paste0(vacc_type_2, '_', vs) ,
                        TRUE ~ 'uv'),
    time_to_hosp_death = ifelse(time_to_hosp < time_to_death, time_to_hosp, time_to_death),
    age_gp = cut(ageYear, breaks = c(-1, 15, 39, 59, 74, 120),
                 labels=c("0-15","16-39","40-59","60-74","75+")),
    # People who were in hospital at test get excluded from hospital admissions analysis
    # People who were discharged before collection date have already been excluded
    in_hosp_at_test = if_else(time_to_hosp <= -3,  1L, 0L ),
    # if time to hosp is negative, they were in hospital at specimen collection date
    hosp_covid = case_when( (!is.na(admission_date) & (admission_date - specimen_date <= 14)) ~1,
                            TRUE ~ 0),
    death_covid = case_when(covid_death_cert == 1 | 
                              (!is.na(NRS.Date.Death) & (NRS.Date.Death - specimen_date <= 28)) ~1,
                            TRUE ~ 0),
    vacc = case_when(vs %in% c("v1_28+","v2_0:27","v2_28+") ~ "v1_28+v2",
                     TRUE ~ vs))  %>%
  mutate(time_to_hosp = if_else(!is.na(NRS.Date.Death) & NRS.Date.Death < a_end_pos & is.na(admission_date), 
                                as.numeric(NRS.Date.Death - specimen_date),
                                time_to_hosp ),
         hosp_covid_emerg = if_else(!is.na(emergency) & emergency ,hosp_covid, 0),
         vt = case_when(vacc %in% c("v1_0:27","v1_28+v2") ~ paste(vacc_type, vacc, sep="_"),
                        TRUE ~ 'uv')) %>%
  # time_to_hosp is capped at 15 days because otherwise it doesn't count as a covid hospitalisation
  # time_to_death is capped at 29 days
  mutate(time_to_hosp = case_when(time_to_hosp >=15 ~ 15,
                                  time_to_hosp <0 ~ 0,
                                  TRUE ~ time_to_hosp),
         time_to_death = case_when(time_to_death >=29 ~ 29,
                                   time_to_death < 0 ~ 0,
                                   TRUE ~ time_to_death),
         time_to_hosp_death = case_when(time_to_hosp_death >=29 ~ 29,
                                        time_to_hosp_death < 0 ~ 0,
                                        TRUE ~ time_to_hosp_death),
         hosp_death_covid = hosp_covid_emerg | death_covid,
         # Make uv the baseline category for vt
         vt = fct_relevel(vt, "uv")) %>%
  # Assume hospitalisation/death is 0.1 days after test if they happen on same day
  mutate(time_to_hosp = if_else(time_to_hosp==0,0.1,time_to_hosp),
         time_to_death = if_else(time_to_death==0,0.1,time_to_death),
         time_to_hosp_death = if_else(time_to_hosp_death==0,0.1,time_to_hosp_death)) %>%
  relocate(EAVE_LINKNO)


# Add in ICU deaths
df_pos <- df_pos %>% left_join(covid_icu_death %>% dplyr::select(-SpecimenDate) %>% 
                                 dplyr::rename(date_icu_death = admission_date)) %>% 
  mutate(covid_icu_death = if_else(is.na(date_icu_death), 0L, 1L),
         time_to_icu_death = if_else(is.na(date_icu_death), as.numeric(a_end_pos - specimen_date),
                                     as.numeric(date_icu_death - specimen_date))) %>% 
  mutate(time_to_icu_death = if_else(time_to_icu_death < 0, 0, time_to_icu_death)) %>%
  mutate(time_to_icu_death = if_else(time_to_icu_death==0,0.1,time_to_icu_death))


########## 4 Create study datatset of all who were virally sequenced ########################

# Gather data on all people who have been virally sequenced.
df_seq <- wgs %>% 
      filter(specimen_date >= a_begin) %>%
      left_join( cdw_full %>% select(EAVE_LINKNO, specimen_date, test_result_record_source)) %>%
      distinct() %>%
      # Create columns for lab type, and variants
      mutate(lab = if_else(test_result_record_source == "ECOSS", "nhs","lh"),
             variant = case_when(is.na(VariantofInterest) ~ "not_sequenced",
                             VariantofInterest=="VOC-21APR-02" ~ "delta",
                             TRUE ~ "other")) %>%
      select(EAVE_LINKNO, specimen_date, lineage, variant, lab) %>%
      # Combine with cohort data
      left_join( df_cohort %>% select(-specimen_date), by="EAVE_LINKNO") %>%
      filter(ageYear >= 18 & flag_incon == 0)

# Add all hospitalisations. This is important because some will be in hospital at time of test
z <- df_seq %>% left_join( all_hospitalisations %>% 
                               dplyr::select(-validchi) %>% 
                               filter( !(!is.na(discharge_date) & discharge_date <= a_begin))) %>%
          mutate(admission_date = case_when( discharge_date <= specimen_date ~ as.Date(NA),
                                     TRUE ~ admission_date),
                 emergency = case_when( discharge_date <= specimen_date ~ NA,
                                     TRUE ~ emergency),
                 discharge_date = case_when( discharge_date <= specimen_date ~ as.Date(NA),
                                     TRUE ~  discharge_date))

# Separate multiple admissions and single/zero admissions
z_id <- z %>% filter(duplicated(EAVE_LINKNO)) %>% pull(EAVE_LINKNO) %>% unique()

z_01 <- z %>% filter(!(EAVE_LINKNO %in% z_id)) # 0 or 1 admission
z_m <- z %>% filter((EAVE_LINKNO %in% z_id))  #multiple admissions

# Select first post specimen hospital admission for those with multiple admissions
z_m <- z_m %>%  
       arrange(EAVE_LINKNO, admission_date) %>% 
       filter(!duplicated(EAVE_LINKNO))  

df_seq <- bind_rows(z_01, z_m)

# End date is maximum of admission date and discharge date
a_end <- max(c(max(df_seq$admission_date, na.rm=T), max(df_seq$discharge_date, na.rm=T)))  

# Add variables  
df_seq <- df_seq %>%        
      mutate(
        variant = ifelse(lineage == 'AY.4.2', 'AY.4.2', variant),
        ageYear = ifelse(ageYear >= 100, 100, ageYear),
        days = as.numeric(specimen_date -min(specimen_date, na.rm = TRUE)),
        death = if_else(is.na(NRS.Date.Death), 0L, 1L),
        vs1 = case_when(is.na(date_vacc_1) | date_vacc_1 > specimen_date ~ "uv",
                       date_vacc_1 <= specimen_date &  date_vacc_1 > specimen_date - 28 ~ "v1_0:27",
                       TRUE ~ "v1_28+"),
        vs2 = case_when(is.na(date_vacc_2) | date_vacc_2 > specimen_date ~ "uv",
                       date_vacc_2 <= specimen_date &  date_vacc_2 > specimen_date - 28 ~ "v2_0:27",
                       TRUE ~ "v2_28+"),
        vs = if_else(vs2=="uv", vs1, vs2), 
        time_to_hosp = if_else(!is.na(admission_date), 
                              as.numeric(admission_date - specimen_date), 
                              as.numeric(a_end - specimen_date)) ,
        time_to_death = if_else(!is.na(NRS.Date.Death), 
                               as.numeric(NRS.Date.Death - specimen_date), 
                               as.numeric(a_end - specimen_date))      ) %>% 
      mutate(
        vs_type = case_when(vs %in% c("v1_0:27", "v1_28+") ~ paste0(vacc_type, '_', vs ),
                            vs %in% c("v2_0:27", "v2_28+") ~ paste0(vacc_type_2, '_', vs ),
                            TRUE ~ 'uv'),
        variant = as.factor(variant),
        time_to_hosp_death = ifelse(time_to_hosp < time_to_death, time_to_hosp, time_to_death),
        age_gp = cut(ageYear, breaks = c(-1, 15, 39, 59, 74, 120),
                                             labels=c("0-15","16-39","40-59","60-74","75+")),
        # People who were in hospital at test get excluded from hospital admissions analysis
        # People who were discharged before collection date have already been excluded
        in_hosp_at_test = if_else(time_to_hosp <= -3,  1L, 0L ),
        # if time to hosp is negative, they were in hospital at specimen collection date
        hosp_covid = case_when( (!is.na(admission_date) & (admission_date - specimen_date <= 14)) ~1,
                                        TRUE ~ 0),
        death_covid = case_when(covid_death_cert == 1 | 
                (!is.na(NRS.Date.Death) & (NRS.Date.Death - specimen_date <= 28)) ~1,
                                                                TRUE ~ 0),
        vacc = case_when(vs %in% c("v1_28+","v2_0:27","v2_28+") ~ "v1_28+v2",
                         TRUE ~ vs))  %>%
      mutate(time_to_hosp = if_else(!is.na(NRS.Date.Death) & NRS.Date.Death < a_end & is.na(admission_date), 
                        as.numeric(NRS.Date.Death - specimen_date),
                        time_to_hosp ),
             hosp_covid_emerg = if_else(!is.na(emergency) & emergency ,hosp_covid, 0),
      vt = case_when(vacc %in% c("v1_0:27","v1_28+v2") ~ paste(vacc_type, vacc, sep="_"),
                            TRUE ~ 'uv')) %>%
      # time_to_hosp is capped at 15 days because otherwise it doesn't count as a covid hospitalisation
      # time_to_death is capped at 29 days
      mutate(time_to_hosp = case_when(time_to_hosp >=15 ~ 15,
                              time_to_hosp <0 ~ 0,
                              TRUE ~ time_to_hosp),
             time_to_death = case_when(time_to_death >=29 ~ 29,
                              time_to_death < 0 ~ 0,
                              TRUE ~ time_to_death),
             time_to_hosp_death = case_when(time_to_hosp_death >=29 ~ 29,
                                            time_to_hosp_death < 0 ~ 0,
                                            TRUE ~ time_to_hosp_death),
             hosp_death_covid = hosp_covid_emerg | death_covid,
      # Make uv the baseline category for vt
             vt = fct_relevel(vt, "uv"),
      # Make delta the baseline variant cateogry
             variant = fct_relevel(variant, "delta")) %>%
      # Assume hospitalisation/death is 0.1 days after test if they happen on same day
      mutate(time_to_hosp = if_else(time_to_hosp==0,0.1,time_to_hosp),
            time_to_death = if_else(time_to_death==0,0.1,time_to_death),
            time_to_hosp_death = if_else(time_to_hosp_death==0,0.1,time_to_hosp_death)) %>%
      relocate(EAVE_LINKNO)


# Add in ICU deaths
df_seq <- df_seq %>% left_join(covid_icu_death %>% dplyr::select(-SpecimenDate) %>% 
                              dplyr::rename(date_icu_death = admission_date)) %>% 
  mutate(covid_icu_death = if_else(is.na(date_icu_death), 0L, 1L),
         time_to_icu_death = if_else(is.na(date_icu_death), as.numeric(a_end - specimen_date),
                                     as.numeric(date_icu_death - specimen_date))) %>% 
  mutate(time_to_icu_death = if_else(time_to_icu_death < 0, 0, time_to_icu_death)) %>%
  mutate(time_to_icu_death = if_else(time_to_icu_death==0,0.1,time_to_icu_death))

df_pos <- mutate(df_pos, sequenced = ifelse( EAVE_LINKNO %in% df_seq$EAVE_LINKNO, 1, 0))


####################### 5 Save ###############################################################

saveRDS(df_cohort, './data/df_cohort.rds')

saveRDS(df_pos, './data/df_pos.rds')

saveRDS(df_seq, './data/df_seq.rds')

# rm(all_deaths, all_deaths_covid_dth_cert, all_hospitalisations, cdw_full,
#    Cohort_Household, covid_death, covid_hospitalisations, covid_icu_death,
#    datazones, EAVE_cohort, EAVE_Weights, Positive_Tests, qcovid_rg, rg,
#    Vaccinations, wgs, z, z_01, z_m, z1)