##########################################################
## Title: 1st dose COVID-19 vaccine waning
## Code author(s): Rachel Mulholland <rachel.mulholland@ed.ac.uk> 
##                 Chris Robertson <chrisrobertson@nhs.net>
## Description: 00_functions - Unique functions for analysis
##########################################################

#### Libraries ####
library("spatstat")

# Table defaults ---------------------------------------------------------------------------
# This makes table resize or continue over multiple pages in all output types
# PDF powered by kableExtra, Word by flextable
mytable = function(x, caption = "", row.names = FALSE, longtable = TRUE,
                   latex_options = c("hold_position"), font_size = 7.0, ...){
  
  # if not latex or html then else is Word
  if (is_latex_output()) {
    knitr::kable(x, row.names = row.names, align = c("l", "l", "r", "r", "r", "r", "r", "r", "r", "r", "r"),
                 booktabs = TRUE, caption = caption, #longtable = longtable,
                 linesep = "", ...) %>%
      kableExtra::kable_styling(font_size = font_size,
                                latex_options = latex_options)
  } else if(is_html_output()) {
    knitr::kable(x, row.names = row.names, align = c("l", "l", "r", "r", "r", "r", "r", "r", "r", "r", "r"),
                 booktabs = TRUE, caption = caption, longtable = longtable,
                 linesep = "", ...) %>%
      kableExtra::kable_styling(latex_options = c("scale_down", "hold_position"))
  } else {
    flextable::flextable(x) %>%
      flextable::autofit() %>%
      flextable::width(j = 1, width = 1.5) %>%
      flextable::height(i = 1, height = 0.5, part = "header")
  }
}


#### Summary table using weights ####
# Creates a table of cohort summaries using weights
# For categorical variables, the sum of the weights and the % is calculated
# For numerical variables, the weighed mean and weighted sd are calculated, as well as
# the weighted median and weighted IQR

# Input:
# - data = the dataset (must have weights named eave_weight)
# - dependent = a character of the dependent variables name
# - explanatory = a string of characters of the explanatory variables

# Output:
# A table with each explanatory variable as a row (multiple rows for each category if categorical)
# with two columns of the weighted summaries for the levels in the dependent variable

summary_factorlist_wt <- function(data, dependent, explanatory){
  # Create list to put in summaries into each element
  summary_tbl_list <- list()
  
  for(i in 1:length(explanatory)){
    
    # Extract variable
    n <- data %>%
      pull(!!sym(explanatory[i]))
    
    # If numeric then make weighted mean
    if(is.numeric(n)) {
      z_mean <- data %>%
        group_by(!!sym(dependent)) %>%
        summarise(mean = round(weighted.mean(!!sym(explanatory[i]), w = eave_weight, na.rm = TRUE),1),
                  sd = round(sqrt(spatstat.geom::weighted.var(!!sym(explanatory[i]), w = eave_weight)),1)) %>%
        mutate(mean.sd = paste0(mean, " (",sd,")")) %>%
        select(-mean, -sd) %>%
        mutate("characteristic" = explanatory[i]) %>%
        pivot_wider(names_from = !!sym(dependent), values_from = mean.sd) %>%
        relocate(characteristic) %>%
        mutate(levels = "mean.sd")
      
      
      z_median <- data %>%
        group_by(!!sym(dependent)) %>%
        summarise(median = spatstat.geom::weighted.median(!!sym(explanatory[i]), w = eave_weight),
                  q1 = spatstat.geom::weighted.quantile(!!sym(explanatory[i]), w = eave_weight, probs = 0.25),
                  q3 = spatstat.geom::weighted.quantile(!!sym(explanatory[i]), w = eave_weight, probs = 0.75)) %>%
        mutate("characteristic" = explanatory[i]) %>%
        mutate(iqr = q3 -q1) %>%
        mutate(median.iqr = paste0(median, " (",iqr,")")) %>%
        select(-q1, -q3, -median, -iqr) %>%
        pivot_wider(names_from = !!sym(dependent), values_from = median.iqr) %>%
        relocate(characteristic) %>%
        mutate(levels = "median.iqr")
      
      # Combine!!
      summary_tbl_list[[i]] <- full_join(z_mean, z_median)
      
      
      # Else get sum of weights of each level
    } else if (length(unique(data %>% pull(!!sym(dependent) ) ) ) ==1) {
      
      # This is for when there is only one level in the dependent variable
      summary_tbl_list[[i]] <- data %>%
        group_by(!!sym(explanatory[i])) %>%
        summarise(n = sum(eave_weight)) %>%
        ungroup() %>%
        mutate(perc = sprintf("%.1f",round(n/sum(n)*100,1))) %>%
        mutate_if(is.numeric, ~formatC(round(.,0), format = "f", big.mark = ",", drop0trailing = TRUE)) %>%
        mutate(n_perc := paste0(n, " (", perc,"%)")) %>%
        select(-n, -perc) %>%
        dplyr::rename("levels"=explanatory[i], !!dependent := n_perc) %>%
        mutate("characteristic" = explanatory[i]) %>%
        relocate(characteristic)
      
    } else {

      summary_tbl_list[[i]] <- data %>%
        group_by(!!sym(explanatory[i]), !!sym(dependent)) %>%
        summarise(n = sum(eave_weight)) %>%
        ungroup() %>%
        group_by(!!sym(dependent)) %>%
        mutate(perc = sprintf("%.1f",round(n/sum(n)*100,1))) %>%
        mutate_if(is.numeric, ~formatC(round(.,0), format = "f", big.mark = ",", drop0trailing = TRUE)) %>%
        mutate(n_perc := paste0(n, " (", perc,"%)")) %>%
        select(-n, -perc) %>%
        pivot_wider(names_from = !!sym(dependent), values_from = n_perc) %>%
        dplyr::rename("levels"=explanatory[i]) %>%
        mutate("characteristic" = explanatory[i]) %>%
        relocate(characteristic)
      
    }
  }
  
  # Combine list together to make dataset
  summary_tbl_wt <- summary_tbl_list %>%
    reduce(full_join)
  
  summary_tbl_wt
}

##### Table of events and event rates by vaccine category
# Calculates the standardised mean differences (smd) between the uv and vacc for each of 
# the categorical explanatory variables for a vaccine type.

# Input:
# - data = the cohort descriptive dataset - z_chrt_desc

# Output:
# A table with weighted millions of person years spent with in each vaccination category in
# the cohort, and weighted count of events and event rates per million person years for
# hospitaliation, death, and hospitalisation or death post vaccination, and more than 14
# days post vaccination

# Table output to be used to plot comparisons between the matched and overall population 
# (before matching - crude)

event_summary_wt <- function(data){
  
  summary_tbl_list <- list()
  
  # First row is person years spent with each vaccination status in cohort
  first_row <-  t(select(data, starts_with('days'))) %*% pull(data, eave_weight)/(365.21 * 1000) 
  
  dependent <- grep('vacc_at', names(data), value = TRUE)
  
  for (i in 1:length(dependent)){
    summary_tbl_list[[i]] <- data %>%
      group_by(!!sym(dependent[i]) ) %>%
      summarise(n = sum(eave_weight)) %>%
      na.omit() %>%
      mutate(rate = sprintf('%.2f',n/first_row))  %>% 
      # format numbers with commas every 3 digits,  
      mutate_if(is.numeric, ~formatC(., format = "f", big.mark = ",", drop0trailing = TRUE)) %>%
      mutate(n = paste0( n, ' (', rate, ')') )  %>%
      select(-rate) %>%
      pivot_wider(names_from = !!sym(dependent[i]), values_from = n) %>%
      mutate(Event = dependent[i])  
  }
  
  # Combine list together to make dataset
  summary_tbl_wt <- summary_tbl_list %>% 
    reduce(full_join) %>%
    mutate(Event = c('Hospitalisation',
                     '14 days prior to hospitalisation',
                     'Death',
                     '14 days prior to death',
                     'Hospitliation or death',
                     '14 days prior to hospitalisation or death')) 
  
  first_row <- formatC(sprintf('%.2f',first_row), format = "f", big.mark = ",", drop0trailing = TRUE)
  
  names(first_row) <- names(summary_tbl_wt)[1:5]
  
  first_row <- data.frame(as.list(first_row), stringsAsFactors = FALSE) %>% 
    mutate(Event = 'Person years (thousands)') %>%
    relocate(Event)
  
  summary_tbl_wt <- bind_rows(first_row, summary_tbl_wt) %>% relocate(uv, .after = Event)
  
  names(summary_tbl_wt) <- c('Event', 'Unvaccinated', 'First dose ChAdOx1', 
                             'Second dose ChAdOx1', 'First dose BNT162b2', 
                             'Second dose BNT162b2')
  
  summary_tbl_wt
}


