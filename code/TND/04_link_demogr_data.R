
##############################################################################
# Name of file: 04_link_demogr_data.R
# Data release (if applicable):
# Original author(s): Chris Robertson chris.robertson@nhs.scot; eleftheria.vasileiou@ed.ac.uk
# Original date: 20 Jan 2021
# Latest update author (if not using version control) - eleftheria.vasileiou@ed.ac.uk
# Latest update date (if not using version control) - 
# Latest update description (if not using version control)
# Type of script: Descriptive stats
# Written/run onL: R Studio SERVER
# Version of R that the script was most recently run on: R 3.6.1
# Description of content: add demographic data to the linked dataset
# Approximate run time: Unknown
##############################################################################

#Load/read in linked data
df <- readRDS("/conf/EAVE/GPanalysis/progs/EV/TND_VOC/data/df.rds")

#Load demographic data 
EAVE_cohort <- readRDS(paste0(Location,"EAVE/GPanalysis/outputs/temp/Cohort_Demog_Endpoints_Times2021-10-21.rds"))
EAVE_cohort <- filter(EAVE_cohort, !duplicated(EAVE_LINKNO))


#remove all who have died before the beginning
EAVE_cohort <- filter(EAVE_cohort, is.na(NRS.Date.Death) | (!is.na(NRS.Date.Death) & NRS.Date.Death > a_begin))
EAVE_cohort <- mutate(EAVE_cohort, ageYear = ageYear +1) # add 1 year to get age at march 2021
#remove under 16s
EAVE_cohort <- filter(EAVE_cohort, ageYear >= 18)
#EAVE_cohort <- EAVE_cohort %>% mutate(age_gp = cut(ageYear, breaks=c(-1,seq(19,89,by=5),120), 
#                                                   labels=c(paste(c(16,seq(20,85, by=5)),seq(19,89, by=5),sep="-"),"90+")))

#get cohort weights
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
          dplyr::select(EAVE_LINKNO:simd) %>% 
          #Link to RG and BPsmoke
          left_join(eave_rg_bpsmoke, by="EAVE_LINKNO") %>% 
          #Link QCOVID risk groups
          left_join(qcovid_rg %>% 
          select(EAVE_LINKNO:Q_DIAG_CKD_LEVEL, n_risk_gps, bmi_impute), by="EAVE_LINKNO") %>%
          #Link to household size
          left_join(Cohort_Household, by="EAVE_LINKNO")

nrow(cohort)

#remove redundant datasets
rm(EAVE_Weights, EAVE_cohort, eave_rg_bpsmoke, qcovid_rg, Cohort_Household)

saveRDS(cohort, paste0(project_path, "/data/cohort.rds"))

#Add demographics to main TND dataset
df.tnd <- df %>% left_join(cohort, by="EAVE_LINKNO")


#lots of NAs in 'ageYear' and 'Sex' vrbs from cohort
#so, use 'age' and 'subject_sex' vrbs from lab data instead which have no NAs

#remove adults <16years old
df.tnd <- df.tnd %>% filter(age >= 18)

nrow(df.tnd.sym)

saveRDS(df.tnd, paste0(project_path, "/data/df.tnd.rds"))

#Last date of tested 
max(df.tnd$date_ecoss_specimen)
#"2021-10-24"

#remove redundant datasets
rm(df, cohort)




