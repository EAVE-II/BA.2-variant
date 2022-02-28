
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
# Description of content: reads in and prepares laboratory data
# Approximate run time: Unknown
##############################################################################


#####################################################################
#Get all lab data from ECOSS
#####################################################################

#all laboratory data from ecoss 
lab.data <- readRDS("/conf/EAVE/GPanalysis/data/CDW_full.rds")

table(lab.data$test_result, useNA = "always")

#convert specimen date to class date 
lab.data$date_ecoss_specimen <- as.Date(lab.data$date_ecoss_specimen)

#Start date for study period
a_begin = "2021-06-08"
a_begin = as.Date(a_begin)


#Number of pcr, pos/neg tests before a_begin date June 08, 2021 (First Delta Plus VOC in our genomic dataset)
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


#saveRDS(Tests, paste0(project_path, "/data/Tests.rds"))
saveRDS(lab.data, paste0(project_path, "/data/lab.data.rds"))

#remove redundant datasets
rm(Tests)

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

#Remove all rows/records of people with a positive test before June 8,2021 
z_pos <- z_pos %>% mutate(prev_test = if_else(date_ecoss_specimen < a_begin, 1, 0))

table(z_pos$prev_test, useNA = "always")


test_before_8jun <- z_pos %>% 
                    filter(prev_test=="1") %>%
                    select(EAVE_LINKNO) %>%
                    filter(!duplicated(EAVE_LINKNO))

z_pos <- z_pos %>% anti_join(test_before_8jun, by = "EAVE_LINKNO") 

#Select the first positive test
z_pos <- z_pos %>%
         arrange(EAVE_LINKNO, date_ecoss_specimen) %>%
         filter(!duplicated(EAVE_LINKNO))

#drop redundant datasets and variables
z_pos <- z_pos %>% select(-prev_test)
rm(test_before_8jun)

nrow(z_pos)


#####################################################################
#Get the negatives and select one at random
#####################################################################

#set a specific seed no. so whenever you re-run the random index you will get the same records
set.seed(12345678)

#Create a random ID number
lab.data <- lab.data %>% 
            mutate(random_index = runif(nrow(lab.data)))


#Get the negative tests and select one at random (n=2.8m)
z_neg <- lab.data %>% filter(test_result=="NEGATIVE") %>% 
                       arrange(EAVE_LINKNO, random_index) %>% 
                       filter(!duplicated(EAVE_LINKNO))


#Select all negative tests from patients with covid symptoms and tested at lighthouse labs/community
z_neg <- z_neg %>% mutate(lab = if_else(test_result_record_source == "ECOSS", "nhs","lh")) %>%
                   mutate(date_onset_of_symptoms = as.Date(date_onset_of_symptoms)) %>%
                   mutate(flag_covid_symptomatic = if_else(!is.na(flag_covid_symptomatic) & flag_covid_symptomatic=="true", 1L, 0L))
                      

#keep only records with date of symptoms & tests from lighthouse lab
z_neg <- z_neg %>% filter(flag_covid_symptomatic==1 & lab == "lh") 

#find EAVE_LINKNOs with a neg result and who also have a positive one
#remove the negative one
z_sel <- z_neg$EAVE_LINKNO %in% z_pos$EAVE_LINKNO

#checks
table(z_sel)


#remove any patients in the negative dataset that also have a positive test
z_neg <- z_neg %>% filter(!z_sel) 


#remove any negative tests before June 8, 2021
z_neg <- z_neg %>% mutate(prev_test = if_else(date_ecoss_specimen < a_begin, 1, 0))

table(z_neg$prev_test, useNA = "always")


#Remove all rows/records of people with a negative test before May 19,2021
test_before_8jun <- z_neg %>% 
                    filter(prev_test=="1") %>%
                    select(EAVE_LINKNO) %>%
                    filter(!duplicated(EAVE_LINKNO))

z_neg <- z_neg %>% anti_join(test_before_8jun, by = "EAVE_LINKNO")

#drop redundant datasets and variables
z_neg <- z_neg %>% select(-prev_test)
rm(test_before_8jun)

#####################################################################
#Bind positives and negatives
#####################################################################

#Bind pos. and neg. tests together 
lab.tests <- bind_rows(z_pos,z_neg) %>% 
             dplyr::select(-random_index)

#checks
table(lab.tests$test_result, useNA = "always")

#remove redundant datasets
rm(lab.data, z_pos, z_neg, z_sel)

#save datasets 
saveRDS(lab.tests, paste0(project_path, "/data/lab.tests.rds"))



