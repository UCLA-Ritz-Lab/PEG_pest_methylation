## ---------------------------
##
## Script name: 3-clean_data_case_and_ctrl.R
##
## Purpose of script: To clean the data for case and control
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
  list("(?i)metal|copper", "(?i)copper", "(?i)Organophosphorus") %>% 
    map(function(data){
      chem_class %>% 
        filter(str_detect(`chem class (pan)`, data)) %>% 
        filter(chemcode %notin% c("chem1751", "chem153"))
    }) %>% 
    set_names("heavy_metal", "chem_copper", "chem_op") %>% 
    list2env(.GlobalEnv)
  
  
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
           pd_new = case_when(
             pd == 1 | pdstudyparkinsonsdisease == 1 ~ "With PD",
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
  
  pest_methylation_clean <- list(
    list(pest_case, pest_control),
    list(case_ids, control_ids)
  ) %>% 
    pmap(function(data1, data2){
      data1 %>% 
        filter(pegid %in% data2$pegid) %>% 
        # select(-c(date_main_interview_collected, interview_date)) %>% 
        mutate(pd = if_else(pdstudyparkinsonsdisease == 1, 1, 0)) 
      # %>%
      #   mice(m=5, maxit=50, method="pmm", seed=305301666) %>% # impute missing data
      #   complete(1)
    }) %>% 
    set_names("case", "control")
  
  pest_methylation_covar <- pest_methylation_clean %>% 
    reduce(rbind)
  
  create_report(pest_methylation_covar, output_dir = here("reports"), 
                output_file = "covar_report.html")
  

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
  #limit to 1974-2008 and remove -9999
  #and sum across chems for same year and person
  
  list(c_grape_out,r_grape_out) %>% 
    map(function(data){
      pest_methylation_clean %>% 
        map(function(df){
          data %>% 
            inner_join(df %>%
                         select(pegid, indexyr), by = "pegid") %>%
            group_by(pegid) %>%
            filter(year > 1973 & year < 2008 & chempound >= 0) %>% 
            ungroup() %>%
            mutate(chemcode = paste0("chem",chemcode)) %>% 
            group_by(pegid, year, chemcode) %>%
            summarise(sum_total_lbs = sum(chempound),
                      .groups = "keep") %>%
            ungroup() %>% 
            distinct()
        }) %>% 
        set_names("case", "control")
    }) %>% 
    set_names("c_grape_agg", "r_grape_agg") %>% 
    list2env(.,envir = .GlobalEnv)
  
  
  #pull in input to get unexposed and limit to 1974-2008
  # merge each dataframe from two lists and add zero for unexposed
  list(
    pest_methylation_clean,
    c_grape_agg,
    r_grape_agg
  ) %>% 
    pmap(function(data1, data2, data3){
      list(
        list(c_grape_in, r_grape_in) %>% 
          map(function(data){
            data %>% 
              inner_join(data1 %>%
                           select(pegid, indexyr), by = "pegid") %>%
              group_by(pegid) %>%
              filter(year > 1973 & year < 2008) %>% 
              ungroup() %>%
              select(pegid, year)
          }), 
        list(data2, data3)) %>% 
        pmap(left_join, by = c("pegid","year")) %>% 
        map(function(df){
          df %>% 
            mutate(sum_total_lbs = if_else(
              is.na(sum_total_lbs), 0, sum_total_lbs)) %>% 
            distinct()
        }) %>% 
        set_names("occupational", "residential")
    }) %>% 
    set_names("case_agg_yr_all", "control_agg_yr_all") %>% 
    list2env(.,envir = .GlobalEnv)

  
  # Calculate lagged years
  
  pest_methylation_clean %>% 
    map(function(df){
      list(list(c_grape_in, r_grape_in),
           c("c","r")) %>% 
        pmap(function(data1, data2){
          data1 %>% 
            select(pegid, year) %>% 
            distinct() %>% 
            inner_join(df %>% 
                         select(pegid, indexyr, indexyr5, 
                                indexyr10, pd_new, study), 
                       by = "pegid") %>% 
            filter(year > 1973) %>% 
            mutate(exp_yrs_lag_dr = ifelse(year <= indexyr, 1, NA),
                   exp_yrs_lag_5 = ifelse(year <= indexyr5, 1, NA),
                   exp_yrs_lag_10 = ifelse(year <= indexyr10, 1, NA)) %>% 
            group_by(pegid) %>% 
            summarise_at(vars(starts_with("exp")),~sum(.x, na.rm=T),
                         .groups = "keep") %>% 
            mutate(location = data2)
        }) %>% 
        set_names("occupational", "residential")
    }) %>% 
    set_names("exp_yrs_lag_case_list", "exp_yrs_lag_control_list") %>% 
    list2env(.GlobalEnv)
  
  
  list(
    exp_yrs_lag_case_list,
    exp_yrs_lag_control_list
  ) %>% 
    map(function(datalist){
      datalist %>% 
        rbindlist() %>% 
        pivot_wider(id_cols = pegid,
                    names_from = location,
                    values_from = exp_yrs_lag_dr:exp_yrs_lag_10) %>% 
        mutate_all(~replace(.,is.na(.),0))
    }) %>% 
    set_names("exp_yrs_lag_case", "exp_yrs_lag_control") %>% 
    list2env(.GlobalEnv)
}



# clean exposure data -----------------------------------------------------

{
  #exposure matrix
  list(
    list(case_agg_yr_all, control_agg_yr_all),
    pest_methylation_clean
  ) %>% 
    pmap(function(data, covar){
      data %>% 
        map(function(df){
          df %>% 
            inner_join(covar %>% 
                         select(pegid, indexyr, pd_new, study), 
                       by = "pegid") %>%
            group_by(pegid) %>%
            filter(year <= indexyr) %>%
            ungroup() %>% 
            distinct() 
        }) %>% 
        set_names("occupational", "residential")
    }) %>% 
    set_names("exp_address_case", "exp_address_control") %>% 
    list2env(.GlobalEnv)
  
  list(list(exp_address_case, exp_address_control),
       list(exp_yrs_lag_case_list, exp_yrs_lag_control_list),
       pest_methylation_clean )%>%
    pmap(function(data, lag, covar){
      list(data, lag) %>% 
        pmap(function(data1, data2){
          data1 %>%
            group_by(pegid) %>%
            filter(year <= indexyr - 10) %>% 
            ungroup() %>% 
            mutate(ind = 1) %>%
            replace_na(list(sum_total_lbs = 0)) %>% 
            group_by(pegid, chemcode) %>%
            summarise(duration = sum(ind),
                      chemuse = sum(sum_total_lbs),
                      .groups = "keep") %>% 
            ungroup() %>% 
            left_join(data2 %>% select(-location), by = "pegid") %>% 
            left_join(covar,
                      by = "pegid") %>% 
            filter(!is.na(pd) & !is.na(duration)) %>% 
            mutate(dur_wt_10 = duration/exp_yrs_lag_10,
                   lb_wt_10 = chemuse/exp_yrs_lag_10,
                   dur_lb_10 = duration * lb_wt_10)
        }) %>% 
        set_names("occupational", "residential")
    }) %>%
    set_names("exp_wt_case", "exp_wt_control") %>%
    list2env(.GlobalEnv)
  
  test <- exp_wt_case[["occupational"]]

  # list(
  #   list(case_agg_yr_all[[1]], control_agg_yr_all[[1]]),
  #   list(case_agg_yr_all[[2]], control_agg_yr_all[[2]]),
  #   pest_methylation_clean
  # ) %>% 
  #   pmap(function(data1, data2, data3){
  #     list(data1,data2) %>% 
  #       map(function(data){
  #         data %>% 
  #           inner_join(data3 %>% 
  #                        select(pegid, indexyr, pd_new, study), 
  #                      by = "pegid") %>% 
  #           mutate(window=cut(year,breaks = c(1973, 1989, Inf),
  #                             include.lowest = FALSE,
  #                             labels = c("1974-1989","1990-index"))) %>% 
  #           # filter(year <= indexyr) %>%
  #           distinct()     
  #         }) 
  #   }) %>% 
  #   set_names("exp_window_address_case", "exp_window_address_control") %>% 
  #   list2env(.GlobalEnv)
  
  # exp_window_address_all_list <- list(exp_window_address_case, 
  #                                     exp_window_address_control) %>% 
  #   pmap(function(data1, data2){
  #     rbind(data1, data2)
  #   })
  # 
  # list(
  #   list(exp_window_address_case[[1]], exp_window_address_control[[1]]),
  #   list(exp_window_address_case[[2]], exp_window_address_control[[2]]),
  #   list(exp_yrs_lag_case_list[[1]], exp_yrs_lag_control_list[[1]]),
  #   list(exp_yrs_lag_case_list[[2]], exp_yrs_lag_control_list[[2]]),
  #   pest_methylation_clean
  # ) %>% 
  #   pmap(function(df1, df2, df3, df4, df5){
  #     list(list(df1, df2),
  #          list(df3, df4))%>% 
  #       pmap(function(data1,data2){
  #         data1 %>% 
  #           mutate(ind = 1) %>% 
  #           #mutate_at(vars(sum_total_lbs), dec_out) %>%
  #           group_by(pegid, chemcode, window) %>% 
  #           summarise(duration = sum(ind),
  #                     chemuse = sum(sum_total_lbs),
  #                     .groups = "keep") %>% 
  #           # inner_join(pest_cov2 %>% 
  #           # select(pegid, pd, study), by = "pegid") %>% 
  #           full_join(data2 %>% select(-location), by = "pegid") %>% 
  #           right_join(df5, by = "pegid") %>% 
  #           ungroup() %>% 
  #           filter(!is.na(pd) & !is.na(chemuse)) %>% 
  #           mutate(chemuse_wt_10 = chemuse/exp_yrs_lag_10,
  #                  chemuse_wt_5 = chemuse/exp_yrs_lag_5,
  #                  chemuse_wt_dr = chemuse/exp_yrs_lag_dr)
  #       }) 
  #   }) %>% 
  #   set_names("exp_lb_window_case_wt", "exp_lb_window_control_wt") %>% 
  #   list2env(.GlobalEnv)
  
  # list(
  #   list(exp_window_address_case[[1]], exp_window_address_control[[1]]),
  #   list(exp_window_address_case[[2]], exp_window_address_control[[2]]),
  #   list(exp_yrs_lag_case_list[[1]], exp_yrs_lag_control_list[[1]]),
  #   list(exp_yrs_lag_case_list[[2]], exp_yrs_lag_control_list[[2]]),
  #   pest_methylation_clean
  # ) %>% 
  #   pmap(function(df1, df2, df3, df4, df5){
  #     list(list(df1, df2),
  #          list(df3, df4))%>% 
  #       pmap(function(data1,data2){
  #         data1 %>% 
  #           mutate(ind = 1) %>% 
  #           #mutate_at(vars(sum_total_lbs), dec_out) %>%
  #           group_by(pegid, chemcode) %>% 
  #           summarise(duration = sum(ind),
  #                     chemuse = sum(sum_total_lbs),
  #                     .groups = "keep") %>% 
  #           full_join(data2 %>% select(-location), by = "pegid") %>% 
  #           right_join(df5, 
  #                      by = "pegid") %>% 
  #           ungroup() %>% 
  #           filter(!is.na(pd) & !is.na(chemuse)) %>% 
  #           mutate(chemuse_wt_10 = chemuse/exp_yrs_lag_10,
  #                  chemuse_wt_5 = chemuse/exp_yrs_lag_5,
  #                  chemuse_wt_dr = chemuse/exp_yrs_lag_dr)
  #       }) 
  #   }) %>% 
  #   set_names("exp_lb_case_wt", "exp_lb_control_wt") %>% 
  #   list2env(.GlobalEnv)

  
}

# clean methylation data --------------------------------------------------

{
  myvar1 <- quote_all(sampleid,pegid,female,meanmethbysample,
                      aim_self_ethnicity,pdstudynumberyearswithpdatbloodd,
                      gds1,gds1_6,gds1_7,gds1_8,gds_cat,gds_cat1,
                      rfvotecaucasian, c11_depression,c11_depressionage)
  
  myvar2 <- quote_all(sampleid,meanxchromosome,aahoadjcellcounts,
                      bioage4hastaticadjage,ageaccelerationresidual,
                      cd8t,cd4t,nk,bcell,mono,gran,pdstudystudy,
                      pdstudydiseasestatus)
  
  list(heavy_metal, chem_copper, chem_op) %>% 
    map(function(chemlist){
      list(
        list(exp_wt_case, exp_wt_control),
        pest_methylation_clean
      ) %>%
        pmap(function(data1, data2){
          data1 %>% 
            map(function(data){
              list("dur_wt_10", "lb_wt_10", "dur_lb_10") %>% 
                map(function(var){
                  data %>% 
                    select(pegid, chemcode, var) %>% 
                    pivot_wider(
                      id_cols = pegid,
                      names_from = chemcode,
                      values_from = var
                    ) %>% 
                    select(pegid, any_of(chemlist$chemcode)) %>% 
                    full_join(datSamplePEG %>%
                                select(all_of(myvar1)), by = "pegid") %>%
                    left_join(datSampleSteve %>%
                                select(all_of(myvar2)),
                              by = "sampleid") %>% 
                    right_join(data2, by = "pegid") %>% 
                    mutate(
                      rfvotecaucasian2 = case_when(
                        (is.na(rfvotecaucasian) & race == 1) | 
                          (is.na(rfvotecaucasian) & 
                             aim_self_ethnicity == "Caucasian") ~ 1,
                        is.na(rfvotecaucasian) ~ 0.5,
                        TRUE ~ rfvotecaucasian),
                      ethnicity = case_when(
                        is.na(aim_self_ethnicity) | 
                          aim_self_ethnicity == "Asian" ~ "Other",
                        TRUE ~ aim_self_ethnicity),
                      caucasian = if_else(ethnicity == "Caucasian", 0, 1),
                      gds1_5 = if_else(gds1 < 5, 0, 1),
                      c11_3 = case_when(
                        c11_depression == 0 ~ 0,
                        c11_depression == 1 & (c11_depressionage < age_diag-4) ~ 2,
                        TRUE ~ 1),
                      gds5_pd = if_else(
                        pd_new == "With PD", 1 + gds1_5, 0),
                      c11_3m = if_else(
                        pd_new == "With PD", 1 + c11_3, 0),
                      c11_2 = if_else(
                        pd_new == "With PD", 1 + c11_depression, 0),
                      nlr = gran/(cd8t + cd4t + nk + bcell + mono)
                    ) %>% 
                    set_value_labels(
                      female = c("Male" = 0, "Female" = 1)) %>% 
                    modify_if(is.labelled, to_factor) %>% 
                    mutate_at(vars(county), fct_drop) 
                }) %>% 
                set_names("duration", "intensity", "pound-year")
            }) %>% 
            set_names("occupational", "residential")
        }) %>% 
        set_names("case", "control")
    }) %>% 
    set_names("metal_wt_list", "copper_wt_list", "op_wt_list") %>% 
    list2env(.GlobalEnv)
  
  # list(heavy_metal, chem_copper, chem_op) %>% 
  #   map(function(chemlist){
  #     list(
  #       list(exp_wt_case[[1]], exp_wt_control[[1]]),
  #       list(exp_wt_case[[2]], exp_wt_control[[2]]),
  #       pest_methylation_clean
  #     ) %>% 
  #       pmap(function(df1, df2, df3){
  #         list(df1, df2) %>% 
  #           map(function(data){
  #             data %>% 
  #               select(pegid, chemcode, lb_wt_10) %>% 
  #               pivot_wider(
  #                 id_cols = pegid,
  #                 names_from = chemcode,
  #                 values_from = lb_wt_10
  #               ) %>% 
  #               select(pegid, any_of(chemlist$chemcode)) %>% 
  #               full_join(datSamplePEG %>%
  #                           select(all_of(myvar1)), by = "pegid") %>%
  #               left_join(datSampleSteve %>%
  #                           select(all_of(myvar2)),
  #                         by = "sampleid") %>% 
  #               right_join(df3, by = "pegid") %>% 
  #               mutate(
  #                 rfvotecaucasian2 = case_when(
  #                   (is.na(rfvotecaucasian) & race == 1) | 
  #                     (is.na(rfvotecaucasian) & 
  #                        aim_self_ethnicity == "Caucasian") ~ 1,
  #                   is.na(rfvotecaucasian) ~ 0.5,
  #                   TRUE ~ rfvotecaucasian),
  #                 ethnicity = case_when(
  #                   is.na(aim_self_ethnicity) | 
  #                     aim_self_ethnicity == "Asian" ~ "Other",
  #                   TRUE ~ aim_self_ethnicity),
  #                 caucasian = if_else(ethnicity == "Caucasian",0,1),
  #                 gds1_5 = if_else(gds1<5,0,1),
  #                 c11_3 = case_when(
  #                   c11_depression == 0 ~ 0,
  #                   c11_depression == 1 & (c11_depressionage < age_diag-4) ~ 2,
  #                   TRUE ~ 1),
  #                 gds5_pd = if_else(
  #                   pd_new == "With PD", 1 + gds1_5, 0),
  #                 c11_3m = if_else(
  #                   pd_new == "With PD", 1 + c11_3, 0),
  #                 c11_2 = if_else(
  #                   pd_new == "With PD", 1 + c11_depression, 0),
  #                 nlr = gran/(cd8t + cd4t + nk + bcell + mono)
  #               ) %>% 
  #               set_value_labels(
  #                 female = c("Male"=0, "Female"=1)) %>% 
  #               modify_if(is.labelled, to_factor) %>% 
  #               mutate_at(vars(county), fct_drop) 
  #           }) 
  #       }) 
  #   }) %>% 
  #   set_names("metal_wt_list", "copper_wt_list", "op_wt_list") %>% 
  #   list2env(.GlobalEnv)
  
  
  #calculate the n and % of exposure for each chemical: combine
  
  list(heavy_metal, chem_copper, chem_op) %>% 
    map(function(chemlist){
      list(
        list(exp_address_case, exp_address_control),
        pest_methylation_clean
      ) %>% 
        pmap(function(data1, data2){
          data1 %>% 
            map(function(data){
              data %>% 
                filter(pegid %in% data2$pegid & 
                         chemcode %in% chemlist$chemcode) %>%
                group_by(pegid) %>%
                filter(year <= indexyr - 10) %>%
                ungroup()
            }) %>% 
            bind_rows() %>% 
            group_by(pegid, chemcode) %>% 
            summarise(sum_total_lbs = sum(sum_total_lbs),
                      .groups = "keep") %>% 
            left_join(data2, by = "pegid") %>% 
            mutate(exp_ind = if_else(sum_total_lbs > 0, 1, 0)) %>%
            group_by(chemcode, pd_new) %>% 
            summarise(n_exp = sum(exp_ind),
                      .groups = "keep") %>% 
            ungroup() %>% 
            pivot_wider(
              id_cols = chemcode,
              names_from = pd_new,
              values_from = n_exp
            ) %>% 
            # rename(case = `With PD`,
            #        control = `Without PD`) %>% 
            # mutate(n_exp = case + control,
            #        pct = n_exp/806) %>% 
            left_join(chemlist_new, by="chemcode") %>% 
            left_join(chem_class, by = c("chemname","chemcode")) %>% 
            relocate(chemname,chemcode,`chem class (pan)`)
        }) %>% 
        set_names("case", "control")
    }) %>% 
    set_names("metal_long_combine_list", 
              "copper_long_combine_list",
              "op_long_combine_list") %>% 
    list2env(.GlobalEnv)

  
  list(metal_long_combine_list, copper_long_combine_list, 
       op_long_combine_list) %>% 
    map(function(df){
      df %>% 
        reduce(inner_join, by = c("chemname", "chemcode", "chem class (pan)")) %>% 
        rename(case = `With PD`,
               control = `Without PD`) %>%
        mutate(n_exp = case + control,
               pct = n_exp/806) %>% 
        filter(n_exp >= 50 
               # & !is.na(`chem class (pan)`) 
               & !is.na(chemname)
               # & `chem class (pan)` %notin% c("n/a")
        )
    }) %>% 
    set_names("metal_long_combine",
              "copper_long_combine",
              "op_long_combine") %>% 
    list2env(.GlobalEnv)
  
  
  list(metal_long_combine, copper_long_combine, op_long_combine) %>% 
    map(function(df){
      df %>% 
        # filter(n_exp >= 30) %>%
        # filter(chemcode %notin% c("chem1638", "chem164", "chem283",
        #                           "chem353", "chem354")) %>%
        pull(chemcode) 
    }) %>% 
    set_names("metal_filter", "copper_filter", "op_filter") %>% 
    list2env(.GlobalEnv)
  
  
  list(
    list(metal_filter, copper_filter, op_filter),
    list(heavy_metal, chem_copper, chem_op)
  ) %>% 
    pmap(function(data1, data2){
      data2 %>% 
        filter(chemcode %in% data1) %>% 
        # filter(chemcode %notin% c("chem1638", "chem164", "chem283",
        #                           "chem353", "chem354")) %>% 
        arrange(chemcode)
    }) %>%
    set_names("metal_filter_name", "copper_filter_name", "op_filter_name") %>%
    list2env(.GlobalEnv)

  
  metal_todrop <- setdiff(heavy_metal$chemcode, metal_filter)
  copper_todrop <- setdiff(chem_copper$chemcode, copper_filter)
  op_todrop <- setdiff(chem_op$chemcode, op_filter)
  
  # setdiff(names(r_lb_sd_case_wt_10), names(c_lb_sd_case_wt_10))
  
  #drop pegids which are not in grape in/out data & z-transform
  
  # id_remove_c <- c("10259MK15", "10818RH27", "11022FC23", 
  #                  "11887JA27", "12158FS27", "12313LH31", "85491MA39") 
  # id_remove_r <- c("11887JA27")
  
  
  list(
    list(metal_wt_list, copper_wt_list, op_wt_list),
    list(metal_todrop, copper_todrop, op_todrop)
  ) %>% 
    pmap(function(datalist, chem_todrop){
      datalist %>% 
        map(function(dflist){
          dflist %>% 
            map(function(datals){
              list(
                datals,
                list(`as.numeric`, `extreme_remove_percentile_win`, 
                     `extreme_remove_percentile_win`)
              ) %>% 
                pmap(function(df, fun){
                  list(create_quantile, create_evernever, 
                       create_ztrans, `create_counts`) %>% 
                    map(function(method){
                      df %>% 
                        select(-any_of(chem_todrop)) %>%
                        mutate_at(
                          vars(starts_with("chem")),
                          ~replace(.x, is.na(.x)|!is.finite(.x), 0)) %>% 
                        mutate_at(vars(starts_with("chem")), fun) %>% 
                        method
                    }) %>% 
                    set_names("quantile", "evernever", "ztrans", "count")
                })
            })
        }) 
    }) %>% 
    map(function(datalist_processed){
      list(
        datalist_processed,
        list(datalist_processed[["control"]], datalist_processed[["control"]])
      ) %>% 
        pmap(function(dflist1_processed, dflist2_processed){
          list(dflist1_processed, dflist2_processed) %>% 
            pmap(function(datals1_processed, datals2_processed){
              list(datals1_processed, datals2_processed) %>% 
                pmap(function(df1_processed, df2_processed){
                  list(df1_processed[["quantile"]],
                       df1_processed[["evernever"]],
                       df1_processed[["ztrans"]],
                       df1_processed[["count"]] %>% 
                         mutate(
                           count = rowSums(
                             across(starts_with("chem"), 
                                    ~. > median(df2_processed[["count"]][[cur_column()]])))) %>% 
                         select(-starts_with("chem")) %>% 
                         relocate(pegid, count)
                       ) %>% 
                    set_names("quantile", "evernever", "ztrans", "count")
                })
            })
        })
    }) %>% 
    set_names("metal_wt_list_processed", 
              "copper_wt_list_processed",
              "op_wt_list_processed") %>% 
    list2env(.GlobalEnv)
  
  # list(
  #   list(metal_wt_list, copper_wt_list, op_wt_list),
  #   list(metal_todrop, copper_todrop, op_todrop)
  # ) %>% 
  #   pmap(function(df, chem){
  #     list(create_quantile, create_evernever, 
  #          create_ztrans, create_counts) %>% 
  #       map(function(method){
  #         list(
  #           df %>% 
  #             map(function(data){
  #               data[[1]]
  #             }),
  #           df %>% 
  #             map(function(data){
  #               data[[2]]
  #             }),
  #           list(chem, chem)
  #         ) %>% 
  #           pmap(function(df1, df2, df3){
  #             list(
  #               list(df1, df2),
  #               list(id_remove_c, id_remove_r)
  #             ) %>% 
  #               pmap(function(data1, data2){
  #                 data1 %>% 
  #                   select(-any_of(df3)) %>%
  #                   filter(pegid %notin% data2 & !is.na(pegid)) %>%
  #                   mutate_at(vars(starts_with("chem")), 
  #                             ~extreme_remove_percentile_win(.x)) %>%
  #                   mutate_at(vars(starts_with("chem")),
  #                             ~replace(., is.na(.)|!is.finite(.), 0) 
  #                             %>% as.vector()) %>% 
  #                   method
  #               }) 
  #           }) 
  #       }) 
  #   }) %>% 
  #   set_names("metal_wt_list_processed", 
  #             "copper_wt_list_processed",
  #             "op_wt_list_processed") %>% 
  #   list2env(.GlobalEnv)
  
  # lb_sd_metal_wt_10_processed %>% 
  #   set_names("lb_sd_wt_10_quantile_metal", 
  #             "lb_sd_wt_10_evernever_metal",
  #             "lb_sd_wt_10_ztrans_metal",
  #             "lb_sd_wt_10_count_win_metal") %>% 
  #   list2env(.GlobalEnv)
  

  
  # further process with count data
  
  # test <- list(
  #   list(
  #     metal_wt_list_processed[["case"]][["occupational"]][["intensity"]][["count"]],
  #     metal_wt_list_processed[["control"]][["occupational"]][["intensity"]][["count"]]
  #   ), 
  #   list(
  #     metal_wt_list_processed[["control"]][["occupational"]][["intensity"]][["count"]],
  #     metal_wt_list_processed[["control"]][["occupational"]][["intensity"]][["count"]]
  #   ) 
  # ) %>% 
  #   pmap(function(df1, df2){
  #     df1 %>% 
  #       mutate(
  #         count = rowSums(
  #           across(starts_with("chem"), 
  #                  ~. > median(df2[[cur_column()]])))
  #         # ,
  #         # copper_count = rowSums(
  #         #   across(matches(chem_copper), 
  #         #          ~. > median(df2[[cur_column()]])))
  #       ) %>% 
  #       select(-starts_with("chem")) %>% 
  #       relocate(pegid, count)
  #   })
  
  # list(metal_wt_list_processed, 
  #      copper_wt_list_processed, 
  #      op_wt_list_processed) %>% 
  #   map(function(datalist){
  #     list(
  #       datalist,
  #       list(datalist[["control"]], datalist[["control"]])
  #     ) %>% 
  #       pmap(function(dflist1, dflist2){
  #         list(dflist1, dflist2) %>% 
  #           pmap(function(datals1, datals2){
  #             list(datals1, datals2) %>% 
  #               pmap(function(df1, df2){
  #                 df1[["count"]] %>% 
  #                   mutate(
  #                     count = rowSums(
  #                       across(starts_with("chem"), 
  #                              ~. > median(df2[["count"]][[cur_column()]])))
  #                     # ,
  #                     # copper_count = rowSums(
  #                     #   across(matches(chem_copper), 
  #                     #          ~. > median(df2[[cur_column()]])))
  #                   ) %>% 
  #                   select(-starts_with("chem")) %>% 
  #                   relocate(pegid, count)
  #               })
  #           })
  #       })
  #   }) %>% 
  #   set_names("metal_wt_count_list_processed", 
  #             "copper_wt_count_list_processed",
  #             "op_wt_count_list_processed") %>% 
  #   list2env(.GlobalEnv)
  


  
# list(lb_sd_metal_wt_10_processed, lb_sd_copper_wt_10_processed, 
#      lb_sd_op_wt_10_processed) %>% 
#     map(function(df){
#       list(
#         df[[4]],
#         list(df[[4]][[2]], 
#              df[[4]][[2]])
#       ) %>% 
#         map(function(dflist){
#           dflist %>% 
#             pmap(function(df1, df2){
#               df1 %>% 
#                 mutate(
#                   count = rowSums(
#                     across(starts_with("chem"), 
#                            ~. > median(df2[[cur_column()]])))
#                   # ,
#                   # copper_count = rowSums(
#                   #   across(matches(chem_copper), 
#                   #          ~. > median(df2[[cur_column()]])))
#                 ) %>% 
#                 select(-starts_with("chem")) %>% 
#                 relocate(pegid, count)
#             })
#         }) 
#     }) %>% 
#     set_names("lb_sd_metal_wt_10_count", 
#               "lb_sd_copper_wt_10_count",
#               "lb_sd_op_wt_10_count") %>% 
#     list2env(.GlobalEnv)
  

    
  
  # list(
  #   c(lb_sd_metal_wt_10_count[[1]], lb_sd_metal_wt_10_count[[2]]),
  #   list(id_remove_c, id_remove_r, id_remove_c, id_remove_r)
  # ) %>% 
  #   pmap(function(data1, data2){
  #     data1 %>% 
  #       filter(pegid %notin% data2 & !is.na(pegid)) %>% 
  #       pull(sampleid) 
  #   }) %>% 
  #   set_names("sampleid_pd_c", "sampleid_pd_r",
  #             "sampleid_ctrl_c", "sampleid_ctrl_r") %>% 
  #   list2env(.,envir = .GlobalEnv)
  
  peg_noob_nors_win_total <- list(PEG_NOOB_nors_win_filter_ctrl_r, 
                                  PEG_NOOB_nors_win_filter_pd_r) %>% 
    bind_cols()
  
  combined_resid_total <- list(combined_resid_filter_ctrl_r, 
                               combined_resid_filter_pd_r) %>% 
    bind_cols()
  
}



#--------------------------------End of the code--------------------------------
