
##############################################################################
# Name of file: 03_link_lab_vacc_data.R
# Data release (if applicable):
# Original author(s): Chris Robertson chris.robertson@nhs.scot; eleftheria.vasileiou@ed.ac.uk
# Original date: 20 Jan 2021
# Latest update author (if not using version control) - eleftheria.vasileiou@ed.ac.uk
# Latest update date (if not using version control) - 
# Latest update description (if not using version control)
# Type of script: Descriptive stats
# Written/run onL: R Studio SERVER
# Version of R that the script was most recently run on: R 3.6.1
# Description of content: link laboratory data with vaccination data 
# Approximate run time: Unknown
##############################################################################

#Load/read in lab and vaccine data
lab.tests <- readRDS("/conf/EAVE/GPanalysis/progs/EV/TND_waning/data/lab.tests.rds")
Vaccinations <- readRDS("/conf/EAVE/GPanalysis/progs/EV/TND_waning/data/Vaccinations.rds")

#exclude people with inconsistent vacine types for dose 1 & 2
#all inconsistencies will be omitted so vacc_type_2 is the same as vacc_type

#link lab data with vaccine data 
df <- left_join(lab.tests, Vaccinations, by="EAVE_LINKNO") %>%
                filter(flag_incon==0) %>% 
                dplyr::select(-vacc_type_2, -flag_incon)

#convert from string to date class
library(lubridate)
a_begin <- ymd(a_begin)
class(a_begin) #Date

#Number of vaccinated before June 8, 2021 (yes/no)
df <- df %>% mutate(vacc_prev = case_when(date_vacc_1 < a_begin ~ 1,
                                          TRUE ~ 0))

#create days since test 
df <- df %>% mutate(days = as.numeric(date_ecoss_specimen - a_begin))
df <- df %>% mutate(days_vacc_1_test = as.numeric(date_ecoss_specimen - date_vacc_1))
df <- df %>% mutate(days_vacc_2_test = as.numeric(date_ecoss_specimen - date_vacc_2))

#create vacc.status for every week1,2,3,4,5,6,7,8,9,10,11,12,13+
df <- df %>% mutate(vacc_status = case_when(is.na(date_vacc_2) &  days_vacc_1_test <0 ~ "uv",
                                            is.na(date_vacc_2) & days_vacc_1_test >=0 & days_vacc_1_test <= 27 ~ "v1_0:27",
                                            is.na(date_vacc_2) & days_vacc_1_test >=28 ~ "v1_28+",
                                            !is.na(date_vacc_2) & days_vacc_1_test <0 ~ "test_before_vaccination1",
                                            !is.na(date_vacc_2) & days_vacc_2_test <0 ~ "test_before_vaccination2",
                                            !is.na(date_vacc_2) & days_vacc_2_test >=0 & days_vacc_2_test <=6 ~ "v2_0:6",
                                            !is.na(date_vacc_2) & days_vacc_2_test >=7 & days_vacc_2_test <=13 ~ "v2_7:13",
                                            !is.na(date_vacc_2) & days_vacc_2_test >=14 & days_vacc_2_test <=20 ~ "v2_14:20",
                                            !is.na(date_vacc_2) & days_vacc_2_test >=21 & days_vacc_2_test <=27 ~ "v2_21:27",
                                            !is.na(date_vacc_2) & days_vacc_2_test >=28 & days_vacc_2_test <=34 ~ "v2_28:34",
                                            !is.na(date_vacc_2) & days_vacc_2_test >=35 & days_vacc_2_test <=41 ~ "v2_35:41",
                                            !is.na(date_vacc_2) & days_vacc_2_test >=42 & days_vacc_2_test <=48 ~ "v2_42:48",
                                            !is.na(date_vacc_2) & days_vacc_2_test >=49 & days_vacc_2_test <=55 ~ "v2_49:55",
                                            !is.na(date_vacc_2) & days_vacc_2_test >=56 & days_vacc_2_test <=62 ~ "v2_56:62",
                                            !is.na(date_vacc_2) & days_vacc_2_test >=63 & days_vacc_2_test <=69 ~ "v2_63:69",
                                            !is.na(date_vacc_2) & days_vacc_2_test >=70 & days_vacc_2_test <=76 ~ "v2_70:76",
                                            !is.na(date_vacc_2) & days_vacc_2_test >=77 & days_vacc_2_test <=83 ~ "v2_77:83",
                                            !is.na(date_vacc_2) & days_vacc_2_test >=84 & days_vacc_2_test <=90 ~ "v2_84:90",
                                            !is.na(date_vacc_2) & days_vacc_2_test >=91 & days_vacc_2_test <=97 ~ "v2_91:97",
                                            !is.na(date_vacc_2) & days_vacc_2_test >=98 & days_vacc_2_test <=104 ~ "v2_98:104",
                                            !is.na(date_vacc_2) & days_vacc_2_test >=105  ~ "v2_105+"))

#checks
table(df$vacc_status, useNA = "always")

#test_before_vaccination1 refers to those that at the time of being tested had not taken their first/second doses
#so they were vaccinated after they were tested - assign them as unvaccinated
#test_before_vaccination2 refers to those that at the time of being tested had taken first dose but not second dose
#so assign them as v1_0:27 or v1_28+ depending on days between vaccine and test date
df <- df %>% mutate(vacc_status = case_when(vacc_status == "test_before_vaccination1" ~ "uv",
                                            vacc_status == "test_before_vaccination2" & days_vacc_1_test >= 0 & days_vacc_1_test <=27 ~ "v1_0:27",
                                            vacc_status == "test_before_vaccination2" & days_vacc_1_test >=28 ~ "v1_28+",
                                             TRUE ~ vacc_status))
                                             

#checks
table(df$vacc_status, useNA = "always")

#Assign NAs (if there are any) in vacc_type to uv 
table(df$vacc_type, useNA = "always")

#Add vaccine type 
df <- df %>% mutate(vs_type = paste(vacc_type, vacc_status, sep="_"))

#checks
table(df$vs_type, useNA = "always")

#Assign the AZ_uv, Mo_uv, PB_uv and uv_uv to unvaccinated
df <- df %>% mutate(vs_type = case_when(vs_type == "AZ_uv" ~ "uv",
                                        vs_type == "Mo_uv" ~ "uv",
                                        vs_type == "PB_uv" ~ "uv",
                                        TRUE ~ vs_type))
                                    
                                      
#checks
table(df$vs_type, useNA = "always")

#collapse vaccine groups into smaller groups cause don't have lots of numbers for sequence data
df$vacc_status <- as.character(df$vacc_status)
df$vs_type <- as.character(df$vs_type)


df <- df %>% mutate(vacc_status2 = case_when(vacc_status %in% c("v1_0:27", "v1_28+") ~ "v1",
                                                   vacc_status %in% c("v2_0:6","v2_7:13") ~ "v2_0:13",
                                                   vacc_status %in% c("v2_14:20", "v2_21:27", "v2_28:34",  
                                                                      "v2_35:41", "v2_42:48", "v2_49:55",
                                                                      "v2_56:62", "v2_63:69","v2_70:76",
                                                                      "v2_77:83", "v2_84:90","v2_91:97", 
                                                                      "v2_98:104","v2_105+") ~ "v2_14+", 
                                                                      TRUE~"uv"))
                                                   
#checks
table(df$vacc_status2, useNA = "always")

#Add vaccine type
df <- df %>% mutate(vs_type2 = paste(vacc_type, vacc_status2, sep="_"))

table(df$vs_type2, useNA = "always")

#Assign the AZ_uv, Mo_uv, PB_uv and uv_uv to unvaccinated
df <- df %>% mutate(vs_type2 = case_when(vs_type2 == "AZ_uv" ~ "uv",
                                        vs_type2 == "Mo_uv" ~ "uv",
                                        vs_type2 == "PB_uv" ~ "uv",
                                        TRUE ~ vs_type2))
#checks
table(df$vs_type2, useNA = "always")


saveRDS(df, paste0(project_path, "/data/df.rds"))


#remove redundant datasets
rm(Vaccinations, lab.tests)




