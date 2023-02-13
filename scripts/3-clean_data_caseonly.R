## ---------------------------
##
## Script name: 3-clean_data_caseonly
##
## Purpose of script: To clean case data
##
## Author: Yufan Gong
##
## Date Created: 2023-02-11
##
## Copyright (c) Yufan Gong, 2023
## Email: ivangong@ucla.edu
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

# clean duration data -----------------------------------------------------

#1. Get the clean peg_cov & peg_keyvar data

{
  heavy_metal <- chem_class %>% 
    filter(str_detect(`chem class (pan)`, "(?i)metal"))
  
  keyvars <- quote_all(pegid, smoker, age_diag, county, sex, race, 
                       minority, a1_schyrs, date_main_interview_collected)
  
  peg_keyvars_case <- list(peg1_keyvar, peg2_keyvar) %>% 
    map(function(data){
      data %>% 
        select(all_of(keyvars)) %>% 
        mutate_all(~replace(.x, .x == -6, NA)) %>% 
        filter(nchar(pegid)>5) %>% 
        mutate(interview_date = year(date_main_interview_collected))
    }) %>% 
    bind_rows() %>% 
    distinct() %>% 
    rename(smokers = smoker) %>% 
    set_variable_labels(
      pegid = "PEGID",
      sex = "Gender",
      minority = "Ethnicity",
      smokers = "Smoking status",
      a1_schyrs = "School years",
      county = "County",
      race = "Race",
    ) %>% 
    set_value_labels(
      minority = c("Non-White"=1, "White"=0),
      race = c("White" = 1, "Black" = 2, "Latino" = 3, 
               "Asian" = 4, "Native American" = 5, "Other" = 9),
      sex = c("Male"=1, "Female"=2),
      smokers = c("Non-smoker"=0, "Former smoker"=1, "Current smoker"=2),
      county = c("Fresno"=1, "Kern"=2, "Tulare"=3, "Other" =4)
    )

  case_ids <- datSamplePEG %>% 
    select(pegid, pdstudyparkinsonsdisease, 
           pdstudydate_main_interview_colle, pdstudydatediagnosed, 
           pdstudyphys_apptdate) %>% 
    filter(pdstudyparkinsonsdisease == 1)
  
  pest_case <- peg_keyvars_case %>% 
    full_join(pest_cov_more %>% select(pegid, pd), by = "pegid") %>% 
    full_join(datSamplePEG %>% 
                select(pegid, pdstudyparkinsonsdisease, age), by = "pegid") %>% 
    left_join(indexyr %>% 
                select(pegid, index_date), by = "pegid") %>% 
    rows_update(tibble(pegid = "80178CW37", a1_schyrs = 16, county = 3,
                       smokers = 1, race = 1, minority = 1)) %>%
    rows_update(tibble(pegid = "82413PR36", a1_schyrs = 9, county =2,
                       smokers = 1, race = 1, minority = 1, sex = 1)) %>%
    rows_update(tibble(pegid = "83620KP43", a1_schyrs = 12, county=1,
                       smokers = 2, race = 1, minority = 1, sex = 2)) %>% 
    mutate(study = if_else(grepl('^1', pegid), "PEG 1", "PEG 2"),
           race_new = case_when(race == 1 ~ "White",
                                race == 3 ~ "Hispanic",
                                TRUE ~ "Other"),
           index_date = if_else(is.na(index_date), interview_date, index_date), 
           age_diag = if_else(is.na(age_diag), age, age_diag),
           pd_new = case_when(pd == 1 | pdstudyparkinsonsdisease == 1 ~ "With PD",
                              pd == 0 | pdstudyparkinsonsdisease == 0 ~ "Without PD",
                              TRUE ~ as.character(pd)),
           indexyr = index_date,
           indexyr5 = indexyr-5,
           indexyr10 = indexyr-10) %>% 

    filter(pd_new == "With PD") %>% 
    modify_if(is.labelled, to_factor) %>% 
    mutate_at(vars(pd_new, race_new), fct_rev)

  
  pest_case_methylation <- pest_case %>% 
    filter(pegid %in% case_ids$pegid)
  summary(pest_case_methylation$sex)
}

{
  # list(c_dur_nolag_dr,r_dur_nolag_dr) %>% 
  #   map(function(data){
  #     data %>% 
  #       select_if(colSums(.>0)>=25) %>%  
  #       select(-pd) %>% 
  #       filter(pegid %in% pest_cov2$pegid) %>% 
  #       left_join(pest_cov2, by="pegid") %>% 
  #       set_variable_labels(agenew = "Age")
  #   }) %>% 
  #   set_names("c_dur_nolag_dr_clean",
  #             "r_dur_nolag_dr_clean") %>% 
  #   list2env(.,envir = .GlobalEnv)
  
  chemlist_new <- chemlist %>% 
    mutate(chemcode=glue("chem{chemcode}")) %>% 
    distinct()
}

#3. clean GRAPES data

{
  #limit to 1974-indexyr and remove -9999
  #and sum across chems for same year and person
  list(c_grape_out,r_grape_out) %>% 
    map(function(data){
      data %>% 
        filter(year > 1973 & year < 2008 & chempound >= 0 
               & pegid %in% case_ids$pegid) %>% 
        mutate(chemcode = paste0("chem",chemcode)) %>% 
        group_by(pegid, year, chemcode) %>%
        summarise(sum_total_lbs = sum(chempound),
                  .groups = "keep") %>%
        ungroup()
    }) %>% 
    set_names("c_grape_case_agg","r_grape_case_agg") %>% 
    list2env(.,envir = .GlobalEnv)
  
  #pull in input to get unexposed and limit to 1974-indexyr
  # merge each dataframe from two lists and add zero for unexposed
  list(
    list(c_grape_in,r_grape_in) %>% 
      map(function(data){
        data %>% 
          filter(year > 1973 & year < 2008 & pegid %in% case_ids$pegid) %>% 
          select(pegid, year)
      }), 
    list(c_grape_case_agg,r_grape_case_agg)) %>% 
    pmap(left_join, by = c("pegid","year")) %>% 
    map(function(data){
      data %>% 
        mutate(sum_total_lbs = ifelse(is.na(sum_total_lbs),0,sum_total_lbs))
    }) %>% 
    set_names("c_agg_yr_case_all","r_agg_yr_case_all") %>% 
    list2env(.,envir = .GlobalEnv)
  
  list(list(c_grape_in,r_grape_in),
       c("c","r")) %>% 
    pmap(function(data1,data2){
      data1 %>% 
        select(pegid,year) %>% 
        distinct() %>% 
        inner_join(pest_case_methylation %>% 
                     select(pegid, indexyr, indexyr5, indexyr10, pd_new, study), 
                   by = "pegid") %>% 
        filter(year > 1973) %>% 
        mutate(exp_yrs_lag_dr = ifelse(year<=indexyr,1,NA),
               exp_yrs_lag_5 = ifelse(year<=indexyr5,1,NA),
               exp_yrs_lag_10 = ifelse(year<=indexyr10,1,NA)) %>% 
        group_by(pegid) %>% 
        summarise_at(vars(starts_with("exp")),~sum(.x,na.rm=T),
                     .groups = "keep") %>% 
        mutate(location = data2)
    }) %>% 
    set_names("exp_yrs_lag_case_c", "exp_yrs_lag_case_r") %>% 
    list2env(.GlobalEnv)
  
  exp_yrs_lag_case <- list(exp_yrs_lag_case_c,exp_yrs_lag_case_r) %>% 
    rbindlist() %>% 
    pivot_wider(id_cols = pegid,
                names_from = location,
                values_from = exp_yrs_lag_dr:exp_yrs_lag_10) %>% 
    mutate_all(~replace(.,is.na(.),0))
  
}


# clean time window ------------------------------------------------------

{
  #number of years of exposure history
  list(c_agg_yr_case_all,r_agg_yr_case_all) %>% 
    map(function(data){
      data %>% 
        right_join(pest_case_methylation %>% 
                     select(pegid, index_date, pd_new, study), 
                   by = "pegid") %>% 
        mutate(window=cut(year,breaks = c(1973, 1989, Inf),
                          include.lowest = FALSE,
                          labels = c("1974-1989","1990-index"))) %>% 
        filter(year <= index_date) %>%
        distinct() 
    }) %>% 
    set_names("exp_window_address_case_c", "exp_window_address_case_r") %>% 
    list2env(.GlobalEnv)
  
  list(list(exp_window_address_case_c,exp_window_address_case_r),
       list(exp_yrs_lag_case_c,exp_yrs_lag_case_r))%>% 
    pmap(function(data1,data2){
      data1 %>% 
        mutate(ind = 1) %>% 
        #mutate_at(vars(sum_total_lbs), dec_out) %>%
        group_by(pegid, chemcode, window) %>% 
        summarise(duration = sum(ind),
                  chemuse = sum(sum_total_lbs),
                  .groups = "keep") %>% 
        #inner_join(pest_cov2 %>% select(pegid, pd, study), by = "pegid") %>% 
        full_join(data2 %>% select(-location), by = "pegid") %>% 
        right_join(pest_case_methylation, by = "pegid") %>% 
        ungroup() %>% 
        filter(!is.na(pd) & !is.na(chemuse)) %>% 
        mutate(chemuse_wt_10 = chemuse/exp_yrs_lag_10,
               chemuse_wt_5 = chemuse/exp_yrs_lag_5,
               chemuse_wt_dr = chemuse/exp_yrs_lag_dr)
    }) %>% 
    set_names("exp_lb_window_case_wt_c", "exp_lb_window_case_wt_r") %>% 
    list2env(.GlobalEnv)
  
  list(list(exp_window_address_case_c,exp_window_address_case_r),
       list(exp_yrs_lag_case_c,exp_yrs_lag_case_r))%>% 
    pmap(function(data1,data2){
      data1 %>% 
        mutate(ind = 1) %>% 
        #mutate_at(vars(sum_total_lbs), dec_out) %>%
        group_by(pegid, chemcode) %>% 
        summarise(duration = sum(ind),
                  chemuse = sum(sum_total_lbs),
                  .groups = "keep") %>% 
        full_join(data2 %>% select(-location), by = "pegid") %>% 
        right_join(pest_case_methylation, 
                   by = "pegid") %>% 
        ungroup() %>% 
        filter(!is.na(pd) & !is.na(chemuse)) %>% 
        mutate(chemuse_wt_10 = chemuse/exp_yrs_lag_10,
               chemuse_wt_5 = chemuse/exp_yrs_lag_5,
               chemuse_wt_dr = chemuse/exp_yrs_lag_dr)
    }) %>% 
    set_names("exp_lb_case_wt_c", "exp_lb_case_wt_r") %>% 
    list2env(.GlobalEnv)
  
}

# clean methylation data --------------------------------------------------

{
  myvar1 <- quote_all(sampleid,pegid,female,meanmethbysample,
                      aim_self_ethnicity,pdstudynumberyearswithpdatbloodd,
                      gds1,gds1_6,gds1_7,gds1_8,gds_cat,gds_cat1,
                      rfvotecaucasian, c11_depression,c11_depressionage)
  
  myvar2 <- quote_all(sampleid,meanxchromosome,aahoadjcellcounts,
                      bioage4hastaticadjage,ageaccelerationresidual,
                      cd8t,cd4t,nk,bcell,mono,gran,pdstudystudy,pdstudydiseasestatus)
  
  
  list(exp_lb_case_wt_c, exp_lb_case_wt_r) %>% 
    map(function(data){
      data %>% 
        select(pegid,chemcode,chemuse_wt_10) %>% 
        pivot_wider(
          id_cols = pegid,
          names_from = chemcode,
          values_from = chemuse_wt_10
        ) %>% 
        select(pegid, all_of(heavy_metal$chemcode)) %>% 
        full_join(datSamplePEG %>%
                    select(all_of(myvar1)), by = "pegid") %>%
        left_join(datSampleSteve %>%
                    select(all_of(myvar2)),
                  by = "sampleid") %>% 
        right_join(pest_case_methylation, by = "pegid") %>% 
        mutate(rfvotecaucasian2 = case_when((is.na(rfvotecaucasian) & race == 1) | 
                                              (is.na(rfvotecaucasian) & 
                                                 aim_self_ethnicity == "Caucasian") ~ 1,
                                            is.na(rfvotecaucasian) ~ 0.5,
                                            TRUE ~ rfvotecaucasian),
               ethnicity = case_when(is.na(aim_self_ethnicity) | 
                                       aim_self_ethnicity == "Asian" ~ "Other",
                                     TRUE ~ aim_self_ethnicity),
               caucasian = ifelse(ethnicity == "Caucasian",0,1),
               gds1_5 = ifelse(gds1<5,0,1),
               c11_3 = case_when(c11_depression == 0 ~ 0,
                                 c11_depression == 1 & (c11_depressionage < age_diag-4) ~ 2,
                                 TRUE ~ 1),
               gds5_pd = case_when(pd_new == "With PD" ~ 1 + gds1_5,
                                   pd_new == "Without PD" ~ 0),
               c11_3m = case_when(pd_new == "With PD" ~ 1+c11_3,
                                  pd_new == "Without PD" ~ 0),
               c11_2 = case_when(pd_new == "With PD" ~ 1 + c11_depression,
                                 pd_new == "Without PD" ~ 0),
               nlr = gran/(cd8t+cd4t+nk+bcell+mono)) %>% 
        set_value_labels(
          female = c("Male"=0, "Female"=1)) %>% 
        modify_if(is.labelled, to_factor) %>% 
        mutate_at(vars(county), fct_drop) 
    }) %>% 
    set_names("c_lb_sd_case_wt_10", "r_lb_sd_case_wt_10") %>% 
    list2env(.GlobalEnv)
  
  
  id_remove_c <- c("10259MK15", "10818RH27", "11022FC23", 
                   "11887JA27", "12158FS27", "12313LH31", "85491MA39") 
  id_remove_r <- c("11887JA27")
  
  #drop pegids which are not in grape in/out data & winsorize
  c_lb_sd_case_wt_10_new <- c_lb_sd_case_wt_10 %>%
    filter(pegid %notin% id_remove_c &
             !is.na(pegid)) %>% 
    mutate_at(vars(starts_with("chem")), ~dec_out(.x)) %>% 
    mutate_at(vars(starts_with("chem")), ~scale(replace(., is.na(.)|!is.finite(.), 0), center = FALSE)
              %>% as.vector()) 
  
  r_lb_sd_case_wt_10_new <- r_lb_sd_case_wt_10 %>%
    filter(pegid %notin% id_remove_r &
             !is.na(pegid)) %>% 
    mutate_at(vars(starts_with("chem")), ~dec_out(.x)) %>% 
    mutate_at(vars(starts_with("chem")), ~scale(replace(., is.na(.)|!is.finite(.), 0), center = FALSE)
              %>% as.vector()) 
}


list(c_lb_sd_case_wt_10_new, r_lb_sd_case_wt_10_new) %>% 
  map(function(data){
    data %>% 
      pull(sampleid) 
  }) %>% 
  set_names("sampleid_pd_c","sampleid_pd_r") %>% 
  list2env(.,envir = .GlobalEnv)

names(combined_resid_filter_pd_c) <- sampleid_pd_c
names(combined_resid_filter_pd_r) <- sampleid_pd_r
