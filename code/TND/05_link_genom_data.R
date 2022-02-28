
##############################################################################
# Name of file: 05_link_genom_data.R
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
df.tnd <- readRDS("/conf/EAVE/GPanalysis/progs/EV/TND_VOC/data/df.tnd.rds")

#Variant/sequencing data 
wgs <- readRDS("/conf/EAVE/GPanalysis/data/WGS_latest.rds")

#create column with alpha and delta variants
wgs <- wgs %>% mutate(variant = case_when(is.na(VariantofInterest) ~ "not_sequenced",
                                                VariantofInterest=="VOC-20DEC-01" ~ "alpha",
                                                VariantofInterest=="VOC-21APR-02" ~ "delta",
                                                TRUE ~ "other"))     

#checks
table(wgs$variant, useNA = "always")

wgs <- wgs %>% mutate(variant = ifelse(lineage == 'AY.4.2', 'AY.4.2', variant))

#checks
table(wgs$variant, useNA = "always") 

#Add VOCs in the TND study 
wgs <- wgs %>% select(EAVE_LINKNO, variant, VariantofInterest)
wgs <- wgs %>% filter(!duplicated(EAVE_LINKNO)) # remove duplicates 
df.tnd <- df.tnd %>% left_join(wgs, by="EAVE_LINKNO")

#checks
table(df.tnd$variant, useNA = "always")


df.tnd <- df.tnd %>% mutate(variant = ifelse(is.na(variant), "not_sequenced", variant))

#checks
table(df.tnd$variant, useNA = "always")

df.tnd <- df.tnd %>% mutate(result = if_else(test_result=="POSITIVE", 1, 0))

df.tnd <- df.tnd %>% mutate(outcome = case_when(result==0 ~ "negative",
                                        result==1 & variant=="alpha" ~ "alpha",
                                        result==1 & variant=="delta" ~ "delta",
                                        result==1 & variant=="AY.4.2" ~ "delta_plus",
                                        result==1 & variant %in% c("not_sequenced","other") ~ "ns_other",
                                        TRUE~"unknown"))

#checks
table(df.tnd$outcome, useNA = "always")
table(df.tnd$result, df.tnd$variant, useNA = "always")

saveRDS(df.tnd, paste0(project_path, "/data/df.tnd.rds"))

#remove redundant objects
rm(wgs)



