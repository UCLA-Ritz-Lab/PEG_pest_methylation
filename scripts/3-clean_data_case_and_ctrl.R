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
    filter(str_detect(`chem class (pan)`, "(?i)metal|copper")) %>% 
    filter(chemcode %notin% c("chem1751", "chem153"))
  
  keyvars <- quote_all(pegid, smoker, age_diag, county, sex, race, 
                       minority, a1_schyrs, date_main_interview_collected)
  
  peg_keyvars <- list(peg1_keyvar, peg2_keyvar) %>% 
    map(function(data){
      data %>% 
        select(all_of(keyvars)) %>% 
        mutate_all(~replace(.x, .x == -6, NA)) %>% 
        filter(str_length(pegid)>5) %>% 
        mutate(interview_date = year(date_main_interview_collected))
    }) %>% 
    bind_rows() %>% 
    distinct() %>% 
    dplyr::rename(smokers = smoker) %>% 
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
  
datSamplePEG %>% 
  select(pegid, pdstudyparkinsonsdisease, 
           pdstudydate_main_interview_colle, pdstudydatediagnosed, 
           pdstudyphys_apptdate) %>% 
  filter(!is.na(pdstudyparkinsonsdisease)) %>% 
  group_by(pdstudyparkinsonsdisease) %>% 
  group_split() %>% 
  set_names("control_ids", "case_ids") %>% 
  list2env(.GlobalEnv)
  
  
peg_keyvars %>% 
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
  filter(!is.na(pd_new)) %>% 
  group_by(pd_new) %>% 
  modify_if(is.labelled, to_factor) %>% 
  mutate_at(vars(pd_new, race_new), fct_rev) %>% 
  group_split() %>% 
  set_names("pest_case", "pest_control") %>% 
  list2env(.GlobalEnv)

list(
  list(pest_case, pest_control),
  list(case_ids, control_ids)
) %>% 
  pmap(function(data1, data2){
    data1 %>% 
      filter(pegid %in% data2$pegid)
  }) %>% 
  set_names("pest_case_methylation", "pest_control_methylation") %>% 
  list2env(.GlobalEnv)

  
  # summary(pest_case_methylation$sex)
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
      list(pest_case_methylation, pest_control_methylation) %>% 
        map(function(df){
          data %>% 
            inner_join(df %>% 
                         select(pegid, indexyr), by = "pegid") %>% 
            group_by(pegid) %>% 
            filter(year > 1973 & year <= indexyr & chempound >= 0) %>% 
            ungroup() %>% 
            mutate(chemcode = paste0("chem",chemcode)) %>% 
            group_by(pegid, year, chemcode) %>%
            summarise(sum_total_lbs = sum(chempound),
                      .groups = "keep") %>%
            ungroup()
        })
    }) %>% 
    set_names("c_grape_agg","r_grape_agg") %>% 
    list2env(.,envir = .GlobalEnv)
  
  c_grape_agg %>% 
    set_names("c_grape_case_agg", "c_grape_control_agg") %>% 
    list2env(.,envir = .GlobalEnv)
  
  r_grape_agg %>% 
    set_names("r_grape_case_agg", "r_grape_control_agg") %>% 
    list2env(.,envir = .GlobalEnv)
  
  #pull in input to get unexposed and limit to 1974-indexyr
  # merge each dataframe from two lists and add zero for unexposed
  list(
    list(pest_case_methylation,pest_control_methylation),
    list(c_grape_case_agg,c_grape_control_agg),
    list(r_grape_case_agg, r_grape_control_agg)
  ) %>% 
    pmap(function(data1, data2, data3){
      list(
        list(c_grape_in,r_grape_in) %>% 
          map(function(data){
            data %>% 
              inner_join(data1 %>% 
                           select(pegid, indexyr), by = "pegid") %>% 
              group_by(pegid) %>% 
              filter(year > 1973 & year <= indexyr) %>% 
              ungroup() %>% 
              select(pegid, year)
          }), 
        list(data2, data3)) %>% 
        pmap(left_join, by = c("pegid","year")) %>% 
        map(function(df){
          df %>% 
            mutate(sum_total_lbs = ifelse(is.na(sum_total_lbs), 0, sum_total_lbs))
        })
    }) %>% 
    set_names("case_agg_yr_all","control_agg_yr_all") %>% 
    list2env(.,envir = .GlobalEnv)
  
  case_agg_yr_all %>% 
    set_names("c_agg_yr_case_all","r_agg_yr_case_all") %>% 
    list2env(.,envir = .GlobalEnv)
  
  control_agg_yr_all %>% 
    set_names("c_agg_yr_control_all","r_agg_yr_control_all") %>% 
    list2env(.,envir = .GlobalEnv)


  
  # Calculate lagged years
  
  list(pest_case_methylation, pest_control_methylation) %>% 
    map(function(df){
      list(list(c_grape_in,r_grape_in),
           c("c","r")) %>% 
        pmap(function(data1,data2){
          data1 %>% 
            select(pegid,year) %>% 
            distinct() %>% 
            inner_join(df %>% 
                         select(pegid, indexyr, indexyr5, 
                                indexyr10, pd_new, study), 
                       by = "pegid") %>% 
            filter(year > 1973) %>% 
            mutate(exp_yrs_lag_dr = ifelse(year <= indexyr,1,NA),
                   exp_yrs_lag_5 = ifelse(year <= indexyr5,1,NA),
                   exp_yrs_lag_10 = ifelse(year <= indexyr10,1,NA)) %>% 
            group_by(pegid) %>% 
            summarise_at(vars(starts_with("exp")),~sum(.x,na.rm=T),
                         .groups = "keep") %>% 
            mutate(location = data2)
        }) 
    }) %>% 
    set_names("exp_yrs_lag_case_list", "exp_yrs_lag_control_list") %>% 
    list2env(.GlobalEnv)
  
  exp_yrs_lag_case_list %>% 
    set_names("exp_yrs_lag_case_c", "exp_yrs_lag_case_r") %>% 
    list2env(.GlobalEnv)
  
  exp_yrs_lag_control_list %>% 
    set_names("exp_yrs_lag_control_c", "exp_yrs_lag_control_r") %>% 
    list2env(.GlobalEnv)
  
  
  exp_yrs_lag_case <- list(exp_yrs_lag_case_c,exp_yrs_lag_case_r) %>% 
    rbindlist() %>% 
    pivot_wider(id_cols = pegid,
                names_from = location,
                values_from = exp_yrs_lag_dr:exp_yrs_lag_10) %>% 
    mutate_all(~replace(.,is.na(.),0))
  
  exp_yrs_lag_control <- list(exp_yrs_lag_control_c, exp_yrs_lag_control_r) %>% 
    rbindlist() %>% 
    pivot_wider(id_cols = pegid,
                names_from = location,
                values_from = exp_yrs_lag_dr:exp_yrs_lag_10) %>% 
    mutate_all(~replace(.,is.na(.),0))
  
}


# clean time window ------------------------------------------------------

{
  #number of years of exposure history
  list(
    list(c_agg_yr_case_all,c_agg_yr_control_all),
    list(r_agg_yr_case_all, r_agg_yr_control_all),
    list(pest_case_methylation, pest_control_methylation)
  ) %>% 
    pmap(function(data1, data2, data3){
      list(data1,data2) %>% 
        map(function(data){
          data %>% 
            inner_join(data3 %>% 
                         select(pegid, indexyr, pd_new, study), 
                       by = "pegid") %>% 
            mutate(window=cut(year,breaks = c(1973, 1989, Inf),
                              include.lowest = FALSE,
                              labels = c("1974-1989","1990-index"))) %>% 
            # filter(year <= indexyr) %>%
            distinct()     
          }) 
    }) %>% 
    set_names("exp_window_address_case", "exp_window_address_control") %>% 
    list2env(.GlobalEnv)
  
  exp_window_address_case %>% 
    set_names("exp_window_address_case_c", "exp_window_address_case_r") %>% 
    list2env(.GlobalEnv)
  
  exp_window_address_control %>% 
    set_names("exp_window_address_control_c", "exp_window_address_control_r") %>% 
    list2env(.GlobalEnv)

  
  list(
    list(exp_window_address_case_c, exp_window_address_control_c),
    list(exp_window_address_case_r, exp_window_address_control_r),
    list(exp_yrs_lag_case_c, exp_yrs_lag_control_c),
    list(exp_yrs_lag_case_r, exp_yrs_lag_control_r),
    list(pest_case_methylation, pest_control_methylation)
  ) %>% 
    pmap(function(df1, df2, df3, df4, df5){
      list(list(df1, df2),
           list(df3, df4))%>% 
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
            right_join(df5, by = "pegid") %>% 
            ungroup() %>% 
            filter(!is.na(pd) & !is.na(chemuse)) %>% 
            mutate(chemuse_wt_10 = chemuse/exp_yrs_lag_10,
                   chemuse_wt_5 = chemuse/exp_yrs_lag_5,
                   chemuse_wt_dr = chemuse/exp_yrs_lag_dr)
        }) 
    }) %>% 
    set_names("exp_lb_window_case_wt", "exp_lb_window_control_wt") %>% 
    list2env(.GlobalEnv)
  
  exp_lb_window_case_wt %>% 
    set_names("exp_lb_window_case_wt_c", "exp_lb_window_case_wt_r") %>% 
    list2env(.GlobalEnv)
  
  exp_lb_window_control_wt %>% 
    set_names("exp_lb_window_control_wt_c", "exp_lb_window_control_wt_r") %>% 
    list2env(.GlobalEnv)
  
  list(
    list(exp_window_address_case_c, exp_window_address_control_c),
    list(exp_window_address_case_r, exp_window_address_control_r),
    list(exp_yrs_lag_case_c, exp_yrs_lag_control_c),
    list(exp_yrs_lag_case_r, exp_yrs_lag_control_r),
    list(pest_case_methylation, pest_control_methylation)
  ) %>% 
    pmap(function(df1, df2, df3, df4, df5){
      list(list(df1, df2),
           list(df3, df4))%>% 
        pmap(function(data1,data2){
          data1 %>% 
            mutate(ind = 1) %>% 
            #mutate_at(vars(sum_total_lbs), dec_out) %>%
            group_by(pegid, chemcode) %>% 
            summarise(duration = sum(ind),
                      chemuse = sum(sum_total_lbs),
                      .groups = "keep") %>% 
            full_join(data2 %>% select(-location), by = "pegid") %>% 
            right_join(df5, 
                       by = "pegid") %>% 
            ungroup() %>% 
            filter(!is.na(pd) & !is.na(chemuse)) %>% 
            mutate(chemuse_wt_10 = chemuse/exp_yrs_lag_10,
                   chemuse_wt_5 = chemuse/exp_yrs_lag_5,
                   chemuse_wt_dr = chemuse/exp_yrs_lag_dr)
        }) 
    }) %>% 
    set_names("exp_lb_case_wt", "exp_lb_control_wt") %>% 
    list2env(.GlobalEnv)
  
  exp_lb_case_wt %>% 
    set_names("exp_lb_case_wt_c", "exp_lb_case_wt_r") %>% 
    list2env(.GlobalEnv)
  
  exp_lb_control_wt %>% 
    set_names("exp_lb_control_wt_c", "exp_lb_control_wt_r") %>% 
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
  
  list(
    list(exp_lb_case_wt_c, exp_lb_control_wt_c),
    list(exp_lb_case_wt_r, exp_lb_control_wt_r),
    list(pest_case_methylation, pest_control_methylation)
  ) %>% 
    pmap(function(df1, df2, df3){
      list(df1, df2) %>% 
        map(function(data){
          data %>% 
            select(pegid,chemcode,chemuse_wt_10) %>% 
            pivot_wider(
              id_cols = pegid,
              names_from = chemcode,
              values_from = chemuse_wt_10
            ) %>% 
            select(pegid, any_of(heavy_metal$chemcode)) %>% 
            full_join(datSamplePEG %>%
                        select(all_of(myvar1)), by = "pegid") %>%
            left_join(datSampleSteve %>%
                        select(all_of(myvar2)),
                      by = "sampleid") %>% 
            right_join(df3, by = "pegid") %>% 
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
        }) 
    }) %>% 
    set_names("lb_sd_case_wt_10", "lb_sd_control_wt_10") %>% 
    list2env(.GlobalEnv)
  
  lb_sd_case_wt_10 %>% 
    set_names("c_lb_sd_case_wt_10", "r_lb_sd_case_wt_10") %>% 
    list2env(.GlobalEnv)
  
  lb_sd_control_wt_10 %>% 
    set_names("c_lb_sd_control_wt_10", "r_lb_sd_control_wt_10") %>% 
    list2env(.GlobalEnv)
  
  #calculate the n and % of exposure for each chemical: combine
  list(
    list(
      list(c_lb_sd_case_wt_10, r_lb_sd_case_wt_10),
      list(c_lb_sd_control_wt_10, r_lb_sd_control_wt_10)
    ),
    list(569, 237)
  ) %>% 
    pmap(function(data1,data2){
      data1 %>% 
        map(function(data){
          data %>% 
            dplyr::select(pegid,starts_with("chem")) 
        }) %>% 
        bind_rows() %>% 
        distinct() %>% 
        group_by(pegid) %>% 
        summarise_all(sum, na.rm=T) %>% 
        pivot_longer(
          cols = starts_with("chem"),
          names_to = c("chemcode"),
          values_to = "duration") %>% 
        mutate(dur_ind=if_else(duration>0,1,0)) %>% 
        group_by(chemcode) %>% 
        summarise(n_exp = sum(dur_ind),
                  .groups = "keep") %>% 
        mutate(pct = n_exp/data2) %>% 
        left_join(chemlist_new, by="chemcode") %>% 
        left_join(chem_class, by = c("chemname","chemcode")) %>% 
        relocate(chemname,chemcode,`chem class (pan)`)
    }) %>% 
    set_names("dur_long_combine_case", "dur_long_combine_control") %>% 
    list2env(.GlobalEnv)

  setdiff(dur_long_combine_case$chemcode, dur_long_combine_control$chemcode)
  
  dur_long_combine <- list(dur_long_combine_case, dur_long_combine_control) %>% 
    map(function(data){
      data %>% 
        select(-c(pct, `chem class (pan)`))
    }) %>% 
    reduce(inner_join, by = c("chemname", "chemcode")) %>% 
    rename(case = n_exp.x,
           ctrl = n_exp.y) %>% 
    mutate(n_exp = case + ctrl,
           pct = n_exp/806) %>% 
    left_join(chem_class, by = c("chemname","chemcode")) %>% 
    relocate(chemname,chemcode,`chem class (pan)`)
  
  
  metal_filter_case <- dur_long_combine %>% 
    filter(n_exp >= 30) %>% 
    filter(chemcode %notin% c("chem1638", "chem164", "chem283",
                              "chem353", "chem354")) %>%
    pull(chemcode) 
  
  metal_filter_ctrl <- dur_long_combine %>% 
    filter(n_exp >= 30) %>% 
    filter(chemcode %notin% c("chem1638", "chem164", "chem283",
                              "chem353", "chem354")) %>%
    pull(chemcode) 
  
  # metal_name <- str_sort(metal_filter) %>% 
  #   map(function(chem){
  #     heavy_metal %>% 
  #       filter(chemcode %in% chem) %>% 
  #       filter(chemcode %notin% c("chem1638", "chem164", "chem283", 
  #                                 "chem353", "chem354"))
  #   }) %>% 
  #   rbindlist()

  chem_todrop_case <- setdiff(heavy_metal$chemcode, metal_filter_case)
  chem_todrop_ctrl <- setdiff(heavy_metal$chemcode, metal_filter_ctrl)
  
  # setdiff(names(r_lb_sd_case_wt_10), names(c_lb_sd_case_wt_10))
  
  #drop pegids which are not in grape in/out data & z-transform
  
  id_remove_c <- c("10259MK15", "10818RH27", "11022FC23", 
                   "11887JA27", "12158FS27", "12313LH31", "85491MA39") 
  id_remove_r <- c("11887JA27")
  
  
  list(
    list(c_lb_sd_case_wt_10, r_lb_sd_case_wt_10),
    list(id_remove_c, id_remove_r)
  ) %>% 
    purrr::pmap(function(data1, data2){
      data1 %>% 
        filter(pegid %notin% data2 & !is.na(pegid)) %>%
        mutate_at(vars(starts_with("chem")), 
                  ~extreme_remove_percentile_win(.x)) %>%
        mutate_at(vars(starts_with("chem")),
                  ~replace(., is.na(.)|!is.finite(.), 0) %>% as.vector()) %>% 
        mutate_at(vars(starts_with("chem")), 
                  ~if_else(.x > quantile(.x, c(0.25, 0.75), na.rm = TRUE)[2], 
                           "high", "low"))
    }) %>% 
    set_names("c_lb_sd_case_wt_10_quantile", "r_lb_sd_case_wt_10_quantile") %>% 
    list2env(.GlobalEnv)
  
  
  list(
    list(c_lb_sd_case_wt_10, r_lb_sd_case_wt_10),
    list(id_remove_c, id_remove_r)
  ) %>% 
    purrr::pmap(function(data1, data2){
      data1 %>% 
        filter(pegid %notin% data2 & !is.na(pegid)) %>%
        mutate_at(vars(starts_with("chem")), 
                  ~extreme_remove_percentile_win(.x)) %>%
        mutate_at(vars(starts_with("chem")),
                  ~replace(., is.na(.)|!is.finite(.), 0) %>% as.vector()) %>% 
        mutate_at(vars(starts_with("chem")), ~if_else(.x > 0, "ever", "never"))
    }) %>% 
    set_names("c_lb_sd_case_wt_10_evernever", 
              "r_lb_sd_case_wt_10_evernever") %>% 
    list2env(.GlobalEnv)
  
  list(
    list(c_lb_sd_case_wt_10, r_lb_sd_case_wt_10),
    list(id_remove_c, id_remove_r)
  ) %>% 
    purrr::pmap(function(data1, data2){
      data1 %>% 
        filter(pegid %notin% data2 & !is.na(pegid)) %>%
        mutate_at(vars(starts_with("chem")), ~extreme_remove_percentile_win(.x)) %>%
        mutate_at(vars(starts_with("chem")),
                  ~scale(replace(., is.na(.)|!is.finite(.), 0), center = TRUE)
                  %>% as.vector())
    }) %>% 
    set_names("c_lb_sd_case_wt_10_ztrans", "r_lb_sd_case_wt_10_ztrans") %>% 
    list2env(.GlobalEnv)
  

  list(
    list(c_lb_sd_case_wt_10, c_lb_sd_control_wt_10),
    list(r_lb_sd_case_wt_10, r_lb_sd_control_wt_10),
    list(chem_todrop_case, chem_todrop_ctrl)
  ) %>% 
    pmap(function(df1, df2, df3){
      list(
        list(df1, df2),
        list(id_remove_c, id_remove_r)
      ) %>% 
        purrr::pmap(function(data1, data2){
          data1 %>% 
            select(-any_of(df3)) %>%
            filter(pegid %notin% data2 & !is.na(pegid)) %>%
            mutate_at(vars(starts_with("chem")), 
                      ~extreme_remove_percentile_win(.x)) %>%
            mutate_at(vars(starts_with("chem")),
                      ~replace(., is.na(.)|!is.finite(.), 0) %>% as.vector()) %>% 
            mutate(metal_count = rowSums(select(., starts_with("chem")) > 0)) %>% 
            select(-starts_with("chem")) %>% 
            relocate(pegid, metal_count)
        }) 
    }) %>% 
    set_names("lb_sd_case_wt_10_count", 
              "lb_sd_control_wt_10_count") %>% 
    list2env(.GlobalEnv)
  
  
  lb_sd_case_wt_10_count %>% 
    set_names("c_lb_sd_case_wt_10_count", 
              "r_lb_sd_case_wt_10_count") %>% 
    list2env(.GlobalEnv)
  
  lb_sd_control_wt_10_count %>% 
    set_names("c_lb_sd_control_wt_10_count", 
              "r_lb_sd_control_wt_10_count") %>% 
    list2env(.GlobalEnv)

 
}

list(
  list(c_lb_sd_case_wt_10_count, r_lb_sd_case_wt_10_count,
       c_lb_sd_control_wt_10_count, r_lb_sd_control_wt_10_count),
  list(id_remove_c, id_remove_r, id_remove_c, id_remove_r)
) %>% 
  pmap(function(data1, data2){
    data1 %>% 
      filter(pegid %notin% data2 & !is.na(pegid)) %>% 
      pull(sampleid) 
  }) %>% 
  set_names("sampleid_pd_c","sampleid_pd_r",
            "sampleid_ctrl_c","sampleid_ctrl_r") %>% 
  list2env(.,envir = .GlobalEnv)

peg_noob_nors_win_total <- list(PEG_NOOB_nors_win_filter_ctrl_r, 
                                PEG_NOOB_nors_win_filter_pd_r) %>% 
  bind_cols()



# names(peg_noob_nors_win_total)

names(combined_resid_filter_pd_c) <- sampleid_pd_c
names(combined_resid_filter_pd_r) <- sampleid_pd_r
names(combined_resid_filter_ctrl_c) <- sampleid_ctrl_c
names(combined_resid_filter_ctrl_r) <- sampleid_ctrl_r
