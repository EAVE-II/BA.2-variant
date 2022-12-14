##########################################################
## Title: 1st dose COVID-19 vaccine waning
## Code author(s): Steven Kerr steven.kerr@ed.ac.uk
## Description: 02_descriptive - Descriptive analyses on the
##              baseline cohort
##########################################################

# Libraries
library("finalfit")
library(tidyverse)
library(survival)
library(lubridate)

Location = '/conf/'

setwd(paste0(Location, 'EAVE/GPanalysis/analyses/BA.2-variant'))

source('./code/hospitalisation_death/00_functions.R')

df_cohort = readRDS('./data/df_cohort.rds')

df_pos = readRDS('./data/df_pos.rds')

df_seq = readRDS('./data/df_seq.rds')

############## 0 Functions #######################

fun.extract <- function(z.fit) {
  #takes a coxph filt using penalised splines and drops off the ps terms
  #make sure no variable begins ps
  z <- summary(z.fit)
  z <- data.frame(z$conf.int)
  z <- z %>% mutate(names = row.names(z)) %>% 
    filter(!(grepl("^ps", names))) %>% 
    dplyr::relocate(names, .before=1) %>% 
    dplyr::select(-exp..coef.)
  names(z) <- c("names","HR","LCL","UCL")
  z
}


plot_HR <- function(model_fit, term){
  # plots hazard ratios for a single term in a fitted model
  
  hr <- termplot(model_fit, term = term, se = T, plot = F)
  
  var <- names(hr)
  
  hr <- hr[[var]]
  
  hr <- mutate(hr, ucl = y + 1.96*se,
               lcl = y - 1.96*se) %>%
    mutate_at(c('y', 'ucl', 'lcl'), exp)
  
  hr <- do.call(data.frame,lapply(hr, function(x) replace(x, is.infinite(x),NA)))
  
  output <- ggplot(data=hr, aes(x=x, y=y)) + geom_line() +
    geom_ribbon(aes(ymin=lcl, ymax=ucl), linetype=2, alpha=0.1, fill = 'steelblue')  + 
    ylab("Hazard Ratio")
  
  if (var == 'ageYear'){
    output <- output + xlab("Age")
  } else if (var == 'days'){
    output <- output + xlab("Days since first specimen collection date")
  }
}

a_begin <- as.Date("2021-11-01")
# End date is maximum of admission date and discharge date
# a_end <- max(c(max(df_seq$admission_date, na.rm=T), max(df_seq$discharge_date, na.rm=T)))
a_end = as.Date("2022-03-20")

# Create output directory
# output_dir = paste0("./output/", a_begin, "-", a_end)
output_dir = paste0("./output/", a_begin, "-", a_end, '_rerun')
if (!dir.exists(output_dir)) {dir.create(output_dir)}

####################### 0 Exploration ################################################

# Investigations on whether to use s gene positve/negative as a proxy for
# BA.2, BA.1
# The issue is that s gene positive is also delta.
cutoff_date = as.Date('2022-01-06')

df_seq_post = filter(df_seq, specimen_date >= cutoff_date)

df_pos_post = filter(df_pos, specimen_date >= cutoff_date)

# How good a proxy is s gene for BA.2, BA.1
sgene_variant_table =  as.data.frame.matrix( table(df_seq_post$sgene_classification, df_seq_post$variant) )

write.csv(sgene_variant_table, paste0(output_dir, "/sgene_variant_table_death_hosp.csv"))

# How much of each variant do we have in the cutoff sequencing data
table(df_seq_post %>% pull(variant)) 

# How many events in the cutoff sequencing data
events_table = as.data.frame.matrix(table(df_seq_post$hosp_death_covid, df_seq_post$variant))

write.csv(events_table , paste0(output_dir, "/events_table .csv"))

# How many events in the cutoff positive test data
table(df_pos_post $hosp_death_covid, df_pos_post $sgene_classification)

# How many events in the non cutoff sequencing data
table(df_seq$hosp_death_covid, df_seq$variant)


##### 1 Descriptive tables ####

# Uses function summary_factorlist_wt from 00_functions.R

df_cohort$Total = 'Total'

## Full cohort summary tables
rgs <- colnames(df_cohort)[startsWith(colnames(df_cohort), "Q")]

explanatory <- c("Total",
                 "Sex", 
                 "ageYear", 
                 "age_grp",  
                 "vacc_type_comb",
                 "simd2020_sc_quintile", 
                 "ur6_2016_name", 
                 "n_risk_gps",
                 "n_tests", 
                 "ave_hh_age", 
                 "n_hh_gp", 
                 "bmi_cat", 
                 rgs,
                 'EAVE_Smoke',
                 'EAVE_BP')

summary_tbl_wt_chrt <- summary_factorlist_wt(df_cohort, "Total", explanatory = explanatory) 

names(summary_tbl_wt_chrt) <- c('Characteristic', 'Levels', 'Total')

summary_tbl_wt_chrt$Characteristic[duplicated(summary_tbl_wt_chrt$Characteristic)] <- ''

summary_tbl_wt_chrt[1, 'Levels'] <- ''

write.csv(summary_tbl_wt_chrt , paste0(output_dir, "/summary_table_weights_cohort.csv"), row.names = F)


## Summary table for all who tested positive

df_pos$Total = 'Total'

convert_to_factor_cols <- c('in_hosp_at_test', 'lab', 'hosp_covid', 'hosp_covid_emerg',
                            'death_covid', 'death', 
                            setdiff(rgs, 'Q_BMI') )

df_pos <- df_pos %>% mutate_at(convert_to_factor_cols, as.factor)

explanatory <- c("Sex", 
                 "ageYear", 
                 "age_grp",
                 "vs",
                 "vs_type",
                 "in_hosp_at_test",
                 "lab",
                 "hosp_covid",
                 "hosp_covid_emerg",
                 "death_covid",
                 "death",
                 "simd2020_sc_quintile", 
                 "ur6_2016_name", 
                 "n_risk_gps",
                 "n_tests", 
                 "ave_hh_age", 
                 "n_hh_gp", 
                 "bmi_cat", 
                 rgs,
                 'EAVE_Smoke',
                 'EAVE_BP')

summary_tbl_wt_pos <- summary_factorlist_wt(df_pos, "Total", explanatory = explanatory) 

names(summary_tbl_wt_pos) <- c('Characteristic', 'Levels', 'Total')

summary_tbl_wt_pos$Characteristic[duplicated(summary_tbl_wt_pos$Characteristic)] <- ''

summary_tbl_wt_pos[1, 'Levels'] <- ''

write.csv(summary_tbl_wt_pos , paste0(output_dir, "/summary_table_weights_positive.csv"), row.names = F)



# Summary table for those who were sequenced

summary_tbl_seq <- summary_factorlist(df_seq %>% mutate(Q_DIAG_CKD_LEVEL = as.factor(Q_DIAG_CKD_LEVEL)), 
                                      "variant", explanatory = explanatory, add_col_totals = TRUE) %>%
                  rename(Characteristic = label, Levels = levels)

write.csv(summary_tbl_seq, paste0(output_dir, "/summary_table_weights_seq.csv"), row.names = F)



# Combine those who were positive and those who were sequenced
summary_tbl_wt_pos[1, 'Levels'] <- 'F'

# comb_table <- summary_tbl_wt_pos %>%
#                       mutate(Levels =  gsub("mean.sd","Mean (SD)", Levels),
#                              Characteristic = ifelse(Characteristic == '', NA, Characteristic)) %>% 
#                       fill(Characteristic, .direction = 'down' ) %>%
#                       right_join(summary_tbl_seq %>% 
#                         mutate(Characteristic = ifelse(Characteristic == '', NA, Characteristic)) %>%
#                         fill(Characteristic, .direction = 'down'), na_matches = 'na') %>%
#                       mutate(Characteristic = ifelse(!duplicated(Characteristic), Characteristic, ''  )    )
#    
# 
# comb_table <- comb_table[c(139, 1:138), ]
#        
# write.csv(comb_table, paste0(output_dir, "/comb_table.csv"), row.names = F)


##################### 2 Descriptive graphs ########################

# Positive tests by day

z <- wgs %>% select(EAVE_LINKNO, specimen_date) %>%
  mutate(sequenced = 1) %>%
  full_join(Positive_Tests %>% select(EAVE_LINKNO, specimen_date)) %>%
  arrange(EAVE_LINKNO, specimen_date) %>%
  #get one test per person per day
  filter(!duplicated(paste(EAVE_LINKNO, specimen_date)))  %>%
  replace_na(list(sequenced = 0)) %>%
  group_by(specimen_date, sequenced) %>%
  summarise(N=n()) %>%
  mutate(sequenced = case_when(sequenced == 1 ~ 'Sequenced',
                               sequenced == 0 ~ 'Not sequenced'))


grid <- expand.grid(seq(a_begin, a_end, by="days"), c('Not sequenced', 'Sequenced'))

names(grid) <- c('specimen_date', 'sequenced')  

grid <- grid %>%
  left_join(z) %>%
  replace_na( list(N = 0)) 

# For reasons unknown, geom_smooth doesn't like it when I give it a sequence of dates
# as the value for zseq. Need to give it numeric value of days starting at 1970-01-01
# Variant data is reliable until 14 days before end date, so truncate smooth
smooth_start = as.numeric(a_begin - as.Date('1970-01-01'))

smooth_end = smooth_start + as.numeric(a_end - a_begin) -14

grid %>%  ggplot(aes(x=specimen_date, y=N, colour = sequenced)) + geom_point() +
  # Delay of ~2 weeks in sequencing data, so truncate smooth 2 weeks before end
  geom_smooth(xseq = smooth_start:smooth_end) +
  labs(x="Specimen date",y ="Number", colour="Sequenced", title="Positive tests by day") + 
  theme(legend.title = element_blank()) 

ggsave(paste0(output_dir, "/pos_tests_by_day.png"), width=14, height=10, unit="cm")


# Cases by day and variant

z <- df_seq %>% 
  group_by(specimen_date, variant) %>% 
  dplyr::summarise(N=n()) 

grid <- expand.grid(seq(a_begin, a_end, by="days"), c('delta', 'BA.1', 'BA.2', 'other'))

names(grid) <- c('specimen_date', 'variant')  

grid <- grid %>%
  left_join(z) %>%
  replace_na( list(N = 0) )

# For reasons unknown, geom_smooth doesn't like it when I give it a sequence of dates
# as the value for zseq. Need to give it numeric value of days starting at 1970-01-01
# Variant data is reliable until two weeks before end date, so truncate smooth
smooth_start = as.numeric(a_begin - as.Date('1970-01-01'))

smooth_end = smooth_start + as.numeric(a_end - a_begin) - 14

grid %>%  ggplot(aes(x=specimen_date, y=N, colour = variant)) + geom_point() +
  # Delay of ~3 days in case data, so truncate smooth
  geom_smooth(xseq = smooth_start:smooth_end) +
  labs(x="Specimen Date",y ="Number", title="Number of cases per day", colour="Variant") +
  scale_y_log10()

ggsave(paste0(output_dir, "/cases_by_day.png"), width=14, height=10, unit="cm")



# Emergency covid hospitalisation or covid deaths by day

z <- df_seq %>% filter(hosp_death_covid==1 & in_hosp_at_test == 0 & lab == 'lh') %>%
  mutate(event_date = pmin(NRS.Date.Death, admission_date, na.rm = TRUE) ) %>%
  group_by(event_date, variant) %>% 
  dplyr::summarise(N=n())  


grid <- expand.grid(seq(a_begin, a_end, by="days"), c('delta', 'BA.1', 'BA.2', 'other'))

names(grid) <- c('event_date', 'variant')  

grid <- grid %>%
  left_join(z) %>%
  replace_na( list(N = 0) )

grid %>%  ggplot(aes(x=event_date, y=N, colour = variant)) + geom_point() +
  scale_y_continuous(limits = c(0, 12), breaks = 0:max(z$N + 1)) +
  geom_smooth(xseq = smooth_start:smooth_end) +
  labs(x="Event date",y ="Number", colour="Variant", title="Emergency covid hospital admissions or covid deaths by day") +
  facet_wrap(~variant)

ggsave(paste0(output_dir, "/hosp_death_by_variant_day.png"), width=14, height=10, unit="cm")




# Person years to event by variant
z.rv <- "hosp_death_covid" 
z.rv.time <- "time_to_hosp_death" 

fmla.plot <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~  variant "))

z.tab <- pyears(fmla.plot, data=df_seq, in_hosp_at_test == 0 & lab == 'lh', data.frame=TRUE)$data

write.csv(z.tab, paste0(output_dir, "/pyears_by_variant.csv"))



# Person years to event by vaccination status
z.rv <- "hosp_death_covid" 
z.rv.time <- "time_to_hosp_death" 

fmla.plot <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~  vs "))

z.tab <- pyears(fmla.plot, data=df_seq, in_hosp_at_test == 0 & lab == 'lh', data.frame=TRUE)$data

write.csv(z.tab, paste0(output_dir, "/pyears_by_vaccine_status.csv"))



# Person years to event by variant and vaccination status

# Create a combine variant, vaccination status column for pyears calculations
df_seq <- mutate(df_seq, variant_vs = as.factor(paste(variant, vs, sep = '_')))

fmla.plot <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~  variant_vs "))

z.tab <- pyears(fmla.plot, data=df_seq, in_hosp_at_test == 0 & lab == 'lh', data.frame=TRUE)$data

write.csv(z.tab, paste0(output_dir, "/pyears_by_variant_vaccine_status.csv"))



# Person years to event by variant and vaccination status by sex, simd2020_sc_quintile and number 
# of risk groups

fmla.plot <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~ Sex "))

z.tab1 <- pyears(fmla.plot, data=df_seq, in_hosp_at_test == 0 & lab == 'lh', data.frame=TRUE)$data %>%
          rename(Variable = Sex)

fmla.plot <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~  simd2020_sc_quintile  "))

z.tab2 <- pyears(fmla.plot, data=df_seq, in_hosp_at_test == 0 & lab == 'lh', data.frame=TRUE)$data %>%
  rename(Variable = simd2020_sc_quintile)

fmla.plot <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~  n_risk_gps "))

z.tab3 <- pyears(fmla.plot, data=df_seq, in_hosp_at_test == 0 & lab == 'lh', data.frame=TRUE)$data %>%
  rename(Variable = n_risk_gps)

z.tab <- bind_rows(z.tab1, z.tab2, z.tab3) %>%
          mutate(pyears = round(pyears))

write.csv(z.tab, paste0(output_dir, "/pyears_by_sex_simd2020_sc_quintile_n_risk_gps.csv"))


################### 3 Analysis ##########################

# Hazard Ratios for emergency covid hospitalisation or covid death from community
fmla.final <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~  
      pspline(ageYear) + pspline(days) + Sex + simd2020_sc_quintile + n_risk_gps + variant +vs "))

z.fit <- coxph(fmla.final , data=df_seq, subset = in_hosp_at_test == 0 
               & lab == 'lh')

z <- fun.extract(z.fit)

write.csv(z, paste0(output_dir, "/hosp_death_HR.csv"))

# Plot HRs for spline terms
plot_HR(z.fit, 1)

ggsave(paste0(output_dir, "/hosp_death_HR_age.png"), width=14, height=10, unit="cm")

plot_HR(z.fit, 2)

ggsave(paste0(output_dir, "/hosp_death_HR_days.png"), width=14, height=10, unit="cm")





# Hazard Ratios for emergency covid hospitalisation or covid death from community testing, 
# with interaction between variant and vaccine status

fmla.final <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   
    pspline(ageYear) + pspline(days)  + Sex  + simd2020_sc_quintile +  n_risk_gps + variant + variant:vs "))

z.fit <- coxph(fmla.final , data=df_seq, subset = in_hosp_at_test == 0 & lab == 'lh')
z <- fun.extract(z.fit)



write.csv(z, paste0(output_dir, "/hosp_death_HR_int.csv"))

# Plot HRs for spline terms
plot_HR(z.fit, 1)

ggsave(paste0(output_dir, "/hosp_death_HR_int_age.png"), width=14, height=10, unit="cm")

plot_HR(z.fit, 2)

ggsave(paste0(output_dir, "/hosp_death_HR_int_days.png"), width=14, height=10, unit="cm")









# Investigation of interaction nterm

# Same as before, except only include BA.1 and BA.2

# No interaction
fmla.final <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~  
                               pspline(ageYear) + pspline(days) + Sex + simd2020_sc_quintile + n_risk_gps + variant +vs "))

z.fit1 <- coxph(fmla.final , data=df_seq, subset = in_hosp_at_test == 0 
               & lab == 'lh' & variant %in% c("BA.1", "BA.2"))

z <- fun.extract(z.fit1)

write.csv(z, paste0(output_dir, "/hosp_death_HR_BA_only.csv"))

# With interaction

fmla.final <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   
    pspline(ageYear) + pspline(days)  + Sex  + simd2020_sc_quintile +  n_risk_gps + variant + variant:vs "))

z.fit2 <- coxph(fmla.final , data=df_seq, subset = in_hosp_at_test == 0 & lab == 'lh' 
                & variant %in% c("BA.1", "BA.2"))

z <- fun.extract(z.fit2)

write.csv(z, paste0(output_dir, "/hosp_death_HR_int_BA_only.csv"))


anova(z.fit1, z.fit2)





# Simpler vaccination status

df_seq = mutate(df_seq, vs3 = case_when( !is.na(date_vacc_3) ~  'v3',
                                         !is.na(date_vacc_2) ~  'v2',
                                         !is.na(date_vacc_1) ~  'v1',
                                         TRUE ~ 'uv'))


# No interaction
fmla.final <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~  
                               pspline(ageYear) + pspline(days) + Sex + simd2020_sc_quintile + n_risk_gps + variant +vs3 "))

z.fit1 <- coxph(fmla.final , data=df_seq, subset = in_hosp_at_test == 0 
                & lab == 'lh' & variant %in% c("BA.1", "BA.2"))

z <- fun.extract(z.fit1)

write.csv(z, paste0(output_dir, "/hosp_death_HR_BA_only_vs3.csv"))

# With interaction

fmla.final <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   
    pspline(ageYear) + pspline(days)  + Sex  + simd2020_sc_quintile +  n_risk_gps + variant + variant:vs3 "))

z.fit2 <- coxph(fmla.final , data=df_seq, subset = in_hosp_at_test == 0 & lab == 'lh' 
                & variant %in% c("BA.1", "BA.2"))

z <- fun.extract(z.fit2)

write.csv(z, paste0(output_dir, "/hosp_death_HR_int_BA_only_vs3.csv"))


anova(z.fit1, z.fit2)






