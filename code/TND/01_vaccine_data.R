##########################################################
# Name of file: 01_vaccine_data.R
# Data release (if applicable):
# Original author(s): Chris Robertson chris.robertson@phs.scot
# Original date: 27 Aug 2021
# Latest update author (if not using version control) - Chris Robertson chris.robertson@phs.scot
# Latest update date (if not using version control) - 
# Latest update description (if not using version control)
# Type of script: Descriptive stats
# Written/run onL: R Studio SERVER
# Version of R that the script was most recently run on: R 3.6.1
# Description of content: reads in the vaccination data and tidies
# Approximate run time: Unknown
##########################################################

library(plyr)
library(tidyverse)
library(survival)
library(lubridate)
library(dplyr)

Location <- "/conf/"  # Server
#Location <- "//isdsf00d03/"  # Desktop

setwd("/conf/EAVE/GPanalysis/progs/EV/TND_VOC")
project_path <- paste0(Location,"EAVE/GPanalysis/progs/EV/TND_VOC")


#read in the GP vaccination data
z  <- readRDS("/conf/EAVE/GPanalysis/data/cleaned_data/C19vaccine.rds")
#z_dv  <- readRDS(paste0(Location,"EAVE/GPanalysis/data/cleaned_data/C19vaccine_dvprod.rds"))
#print(table(z_dv$vacc_type_code, z_dv$dose_number, exclude=NULL))

#print(unique(z$type))
#print(unique(z$stage))
print(table(z$type, z$stage, exclude=NULL))

#exclude the flu vaccinations
z_flu_vacc <- c("27114211000001105", "36509011000001106","38973211000001108")

Vaccinations <- z %>% mutate(Date = as.Date(occurrence_time)) %>%
  filter( !(type %in% z_flu_vacc)) %>% 
  #  filter(stage != 0) %>% #virtually all of the stage 0 are BAD
  mutate(vacc_type = case_when(type == "39114911000001105" ~ "AZ",
                               type == "39115611000001103" ~ "PB",
                               type == "39326911000001101" ~ "Mo",
                               type == "39326911000001100" ~ "Mo",
                               TRUE ~ "UNK") , 
                               dose_number = stage )  # few stage 3 and 0

v0 <- filter(Vaccinations, dose_number==0) %>% 
  dplyr::select(EAVE_LINKNO, Date, vacc_type, dose_number) %>% 
  arrange(EAVE_LINKNO, Date) %>% 
  filter(!duplicated(EAVE_LINKNO))
v1 <- filter(Vaccinations, dose_number==1) %>% 
  dplyr::select(EAVE_LINKNO, Date, vacc_type, dose_number) %>% 
  arrange(EAVE_LINKNO, Date, vacc_type)
#find duplicates
z <- filter(v1, duplicated(EAVE_LINKNO)) %>% pull(EAVE_LINKNO)
z <- filter(v1, EAVE_LINKNO %in% z)
z <- z %>% group_by(EAVE_LINKNO) %>% mutate(diff = as.numeric(Date - first(Date)) ) %>% ungroup() 
z <- z %>% mutate(first_record = if_else(EAVE_LINKNO == lag(EAVE_LINKNO), 0,1)) %>% 
  mutate(first_record = if_else(is.na(first_record), 1, first_record))
z1 <- z %>%  filter(first_record==1 | first_record==0 & vacc_type != "UNK") # omit the unknowns in second place
z1 <- z1 %>% filter(first_record==1 | first_record==0 & diff > 18) # omit the short gaps
z_first <- z1 %>% filter(first_record==1) %>% dplyr::select(-diff, -first_record)
z_second <- z1 %>% filter(first_record==0) %>% dplyr::select(-diff, -first_record) %>% mutate(dose_number=2) # some duplicates in here
#from individuals with >= 1 first dose take the first one and add it back into v1
v1 <- v1 %>% filter(!(EAVE_LINKNO %in% z_first$EAVE_LINKNO) ) %>% bind_rows(z_first)

v2 <- filter(Vaccinations, dose_number==2) %>% 
  dplyr::select(EAVE_LINKNO, Date, vacc_type, dose_number) %>% 
  bind_rows(z_second) %>% 
  arrange(EAVE_LINKNO, Date, vacc_type) 

#find duplicates
z <- filter(v2, duplicated(EAVE_LINKNO)) %>% pull(EAVE_LINKNO)
z <- filter(v2, EAVE_LINKNO %in% z)
z <- z %>% group_by(EAVE_LINKNO) %>% mutate(diff = as.numeric(Date - first(Date)) ) %>% ungroup() 
z <- z %>% mutate(first_record = if_else(EAVE_LINKNO == lag(EAVE_LINKNO), 0,1)) %>% 
  mutate(first_record = if_else(is.na(first_record), 1, first_record))
z1 <- z %>%  filter(first_record==1 | first_record==0 & vacc_type != "UNK") # omit the unknowns in second place
z1 <- z1 %>% filter(first_record==1 | first_record==0 & diff > 18) # omit the short gaps
z_first <- z1 %>% filter(first_record==1) %>% dplyr::select(-diff, -first_record)
z_second <- z1 %>% filter(first_record==0) %>% dplyr::select(-diff, -first_record) %>% mutate(dose_number=3) # some duplicates in here
#from individuals with >= 1 first dose take the first one and add it back into v1
v2 <- v2 %>% filter(!(EAVE_LINKNO %in% z_first$EAVE_LINKNO) ) %>% bind_rows(z_first)


v3 <- filter(Vaccinations, dose_number==3) %>% 
  dplyr::select(EAVE_LINKNO, Date, vacc_type, dose_number) %>% 
  bind_rows(z_second) %>% 
  arrange(EAVE_LINKNO, Date) %>% 
  filter(!duplicated(EAVE_LINKNO))

#find duplicates - no duplicates

v999 <- filter(Vaccinations, dose_number==999) %>% 
  dplyr::select(EAVE_LINKNO, Date, vacc_type, dose_number) %>% 
  mutate(dose_number=3)

v3 <- v3 %>% bind_rows(v999)

#find duplicates
z <- filter(v3, duplicated(EAVE_LINKNO)) %>% pull(EAVE_LINKNO)
z <- filter(v3, EAVE_LINKNO %in% z)
z <- z %>% arrange(EAVE_LINKNO, Date) %>% group_by(EAVE_LINKNO) %>% mutate(diff = as.numeric(Date - first(Date)) ) %>% ungroup() 
z <- z %>% mutate(first_record = if_else(EAVE_LINKNO == lag(EAVE_LINKNO), 0,1)) %>% 
  mutate(first_record = if_else(is.na(first_record), 1, first_record))
z1 <- z %>%  filter(first_record==1 | first_record==0 & vacc_type != "UNK") # omit the unknowns in second place
z1 <- z1 %>% filter(first_record==1 | first_record==0 & diff > 18) # omit the short gaps
z_first <- z1 %>% filter(first_record==1) %>% dplyr::select(-diff, -first_record)
z_second <- z1 %>% filter(first_record==0) %>% dplyr::select(-diff, -first_record) %>% mutate(dose_number=4) # some duplicates in here
#from individuals with >= 1 first dose take the first one and add it back into v1
v3 <- v3 %>% filter(!(EAVE_LINKNO %in% z_first$EAVE_LINKNO) ) %>% bind_rows(z_first)

Vaccinations <- full_join(v1,v2, by="EAVE_LINKNO") %>% 
  mutate(date_vacc_1 = as.Date(Date.x), 
         date_vacc_2 = as.Date(Date.y) ) %>% 
  dplyr::rename(vacc_type=vacc_type.x,
                vacc_type_2=vacc_type.y) %>% 
  dplyr::select(-dose_number.x, -dose_number.y, -Date.x, -Date.y) %>% 
  mutate(z_sel = is.na(date_vacc_1)) %>% 
  mutate(vacc_type = if_else(z_sel, vacc_type_2, vacc_type),
         date_vacc_1 = if_else(z_sel, date_vacc_2, date_vacc_1),
         vacc_type_2 = if_else(z_sel, NA_character_, vacc_type_2),
         date_vacc_2 = if_else(z_sel, as.Date(NA), date_vacc_2)) %>%  # move vacc2 to vacc1 for those with no vacc 1 and change vacc2 to missing
  dplyr::select(-z_sel)

#add in booster
Vaccinations
z <- full_join(Vaccinations,v3, by="EAVE_LINKNO") %>% 
  mutate(date_vacc_3 = as.Date(Date)) %>% 
  dplyr::rename(vacc_type_3 = vacc_type.y, vacc_type=vacc_type.x) %>% 
  dplyr::select(-dose_number, -Date) %>% 
  # move vacc3 to vacc1 for those with no vacc 1 (and hence no vacc_2) and change vacc2/3 to missing
  mutate(z_sel = is.na(date_vacc_1)) %>% 
  mutate(vacc_type = if_else(z_sel, vacc_type_3, vacc_type),
         date_vacc_1 = if_else(z_sel, date_vacc_3, date_vacc_1),
         vacc_type_2 = if_else(z_sel, NA_character_, vacc_type_2),
         date_vacc_2 = if_else(z_sel, as.Date(NA), date_vacc_2),
         vacc_type_3 = if_else(z_sel, NA_character_, vacc_type_3),
         date_vacc_3 = if_else(z_sel, as.Date(NA), date_vacc_3)) %>% 
  # move vacc3 to vacc1 for those with no vacc 1 (and hence no vacc_2) and change vacc2/3 to missing
  mutate(z_sel = is.na(date_vacc_2)) %>% 
  mutate(vacc_type_2 = if_else(z_sel, vacc_type_3, vacc_type),
         date_vacc_2 = if_else(z_sel, date_vacc_3, date_vacc_2),
         vacc_type_3 = if_else(z_sel, NA_character_, vacc_type_3),
         date_vacc_3 = if_else(z_sel, as.Date(NA), date_vacc_3)) %>% 
  dplyr::select(-z_sel)        
#z <- filter(Vaccinations, is.na(date_vacc_1))
#table(z$vacc_type, exclude=NULL)

#investigate the errors
v03 <- full_join(v0,z_second, by="EAVE_LINKNO") %>% 
  mutate(date_vacc_0 = as.Date(Date.x), 
         date_vacc_4 = as.Date(Date.y) ) %>% 
  dplyr::rename(vacc_type_0=vacc_type.x,
                vacc_type_4=vacc_type.y) %>% 
  dplyr::select(-dose_number.x, -dose_number.y, -Date.x, -Date.y)

#remove any that are in Vaccinations
z_sel <- v03$EAVE_LINKNO %in% Vaccinations$EAVE_LINKNO
print(table(z_sel)) #all in there already so no need to modify
#v03 <- v03 %>% left_join(Vaccinations, by="EAVE_LINKNO")
rm(z,v0, v1,v2,v3, v999, v03, z_sel, z1 ,z_first, z_second)

print(table(Vaccinations$vacc_type, Vaccinations$vacc_type_2, exclude=NULL))
#flag inconsistent records
#Vaccinations <- Vaccinations %>% filter(vacc_type %in% c("AZ","PB", "Mo")) %>% 
#  filter(vacc_type_2 %in% c("AZ","PB", "Mo") | is.na(vacc_type_2)) %>% 
#  filter( !(!is.na(vacc_type_2) & (vacc_type_2 != vacc_type)))
Vaccinations <- Vaccinations %>%  
  mutate(flag_incon = case_when(vacc_type %in% c("AZ","PB", "Mo") & is.na(vacc_type_2) ~ 0,
                                vacc_type %in% c("AZ","PB", "Mo") & !is.na(vacc_type_2) & vacc_type==vacc_type_2 ~ 0,
                                TRUE ~ 1 ))
print(table(Vaccinations$vacc_type, Vaccinations$vacc_type_2, Vaccinations$flag_incon, exclude=NULL))

#second on same day as first - make one dose
Vaccinations <- Vaccinations %>% 
  mutate(vacc_type_2 = if_else(!is.na(date_vacc_2) & (date_vacc_2 == date_vacc_1), NA_character_, vacc_type_2 ) ) %>% 
  mutate(date_vacc_2 = as.Date(ifelse(!is.na(date_vacc_2) & (date_vacc_2 == date_vacc_1), NA, date_vacc_2 ), origin=as.Date("1970-01-01")) )
#omit records with second dose too close to first
Vaccinations <- Vaccinations %>% 
  mutate(flag_incon = if_else(!is.na(date_vacc_2)&(date_vacc_2 <= date_vacc_1 + 18), 1, flag_incon))
#z <- Vaccinations %>% filter(date_vacc_1 < as.Date("2020-12-07"))
Vaccinations <- Vaccinations %>% 
  mutate(flag_incon = if_else(date_vacc_1 <= as.Date("2020-12-08"), 1, flag_incon),
         flag_incon = if_else(!is.na(date_vacc_2)&date_vacc_2 <= as.Date("2020-12-08")+18, 1, flag_incon))


#remove records with vaccination date before 8 Dec. 2020
a_begin <- "2020-12-08"
a_begin <- as.Date(a_begin)

table(Vaccinations$date_vacc_1 < a_begin, useNA = "always") 
table(Vaccinations$date_vacc_2 < a_begin, useNA = "always") 

#remove vacc. records before Dec. 8, 2021 but keep make sure to keeps NAs/unvacc in the vacc data
Vaccinations <- Vaccinations %>% filter(is.na(date_vacc_1) | date_vacc_1 >= a_begin)
Vaccinations <- Vaccinations %>% filter(is.na(date_vacc_2) | date_vacc_2 >= a_begin)

#remove redundant objects
rm(z_flu_vacc)

nrow(Vaccinations)

saveRDS(Vaccinations, paste0(project_path, "/data/Vaccinations.rds"))

rm(Vaccinations)



