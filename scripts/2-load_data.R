#------------------------------------------------------------------------------#
#-------------------------------Created by Yufan-------------------------------#
#------------------------------Date: 01/07/2022--------------------------------#
#-------------------------------To load datasets-------------------------------#
#------------------------------------------------------------------------------#


# Load data ---------------------------------------------------------------

{
  grape_name <- quote_all(c_grape_in,c_walk, c_grape_out, 
                          r_grape_out, r_grape_in, r_walk)
  list.dirs(here(),recursive = FALSE) %>% 
    list.files("\\.txt$", full.names = TRUE, recursive = T) %>%
    grep("(?i)XWALK|(?i)GRAPES|input2",., value=TRUE, ignore.case = TRUE) %>% 
    map(read_delim) %>%
    map(~rename_all(.x, str_to_lower)) %>% 
    set_names(grape_name) %>%
    list2env(.,envir = .GlobalEnv)
  
  exposure <- quote_all(exp_c,c_dur_nolag_dr,c_wt_long,chem_class,
                        chemlist,indexyr,serum,exp_r,r_dur_nolag_dr,r_wt_long)
  list.dirs(here(),recursive = FALSE) %>% 
    list.files("\\.xlsx$", full.names = TRUE, recursive = T) %>%
    grep("nolag|wt_long|chemlist|dr_long|SERUM|Indexyr2|chem_class",., 
         value=TRUE, ignore.case = TRUE) %>% 
    map(read_xlsx) %>% 
    map(~rename_all(.x, str_to_lower)) %>% 
    set_names(exposure) %>%
    list2env(.,envir = .GlobalEnv)
  
  
  pest_name <- quote_all(peg_c,cgep_intv,pest_cov_more,pest_cov,peg_r)
  list.dirs(here(),recursive = FALSE) %>%
    list.files("\\.csv$", full.names = TRUE, recursive = T) %>%
    grep("avg_lag10|Pest_cov|CGEP_intv",., value=TRUE, ignore.case = TRUE) %>%
    keep(~!str_detect(.x,"Summary")) %>%
    map(read_csv) %>%
    map(~rename_all(.x, str_to_lower)) %>%
    set_names(pest_name) %>%
    list2env(.,envir = .GlobalEnv)
  
  #import methylation data
  
  list.dirs(here(),recursive = FALSE) %>%
    list.files("\\.RData$", full.names = TRUE, recursive = T) %>%
    grep("filter_pd",., 
         value=TRUE, ignore.case = TRUE) %>%
    keep(~!str_detect(.x,"win")) %>%
    map(.,load,.GlobalEnv)

  
  
  methl_name <- quote_all(datSampleSteve,datSamplePEG,cpg)
  list.dirs(here(),recursive = FALSE) %>%
    list.files("\\.csv$", full.names = TRUE, recursive = T) %>%
    grep("normalize|DNAmage|SampleAnnotation",., 
         value=TRUE, ignore.case = TRUE) %>%
    map(read_csv) %>% 
    map(~rename_all(.x, str_to_lower)) %>%
    set_names(methl_name) %>%
    list2env(.,envir = .GlobalEnv)
  
  keyvar <- quote_all(cgep_keyvar,peg1_keyvar, peg2_keyvar)
  list.dirs(here(),recursive = FALSE) %>%
    list.files("\\.sas7bdat$", full.names = TRUE, recursive = T) %>%
    grep("keyvar|key_var",., value=TRUE, ignore.case = TRUE) %>% 
    keep(~!str_detect(.x,"2017|2018|(?i)raw|(?i)archive|(?i)grapes")) %>% 
    map(read_sas) %>% 
    map(~rename_all(.x, str_to_lower)) %>%
    set_names(keyvar) %>%
    list2env(.,envir = .GlobalEnv)
  # memory.size()
  # gc()
  # memory.limit(size=56000)
  # listname <- quote_all(chem_wd,chem_wd_new,
  #                       ewas.summary_list_all_c,ewas.summary_list_all_r)
  # list.dirs(here(),recursive = FALSE) %>%
  #   list.files("\\.rds$", full.names = TRUE, recursive = T) %>%
  #   map(read_rds) %>%
  #   #map(~rename_all(.x, str_to_lower)) %>%
  #   set_names(listname) %>%
  #   list2env(.,envir = .GlobalEnv)
  # chem_wd_new <- read_rds(here("data","chem_wd_new.rds"))
  
  # matrix_name <- quote_all(cpg_matrix_c,cpg_matrix_r, pest_label_c, pest_label_r, pest_matrix_c, pest_matrix_r)
  # list.dirs(here(),recursive = FALSE) %>% 
  #   list.files("\\.txt$", full.names = TRUE, recursive = T) %>%
  #   grep("(?i)pest_matrix|cpg_matrix|pest_label",., value=TRUE, ignore.case = TRUE) %>% 
  #   map(read_delim) %>%
  #   map(~rename_all(.x, str_to_lower)) %>% 
  #   set_names(matrix_name) %>%
  #   list2env(.,envir = .GlobalEnv)
  
  # load("meffil_list_all_yg_win_c.RData")
  # load("meffil_list_all_yg_win_r.RData")
  # load("combined_resid_new_c.RData")
  # load("combined_resid_new_r.RData")
  # load(here("data", "methylation", "processed", "champ_dmplist_all_kp_win_c.RData"))
  # load(here("data", "methylation", "processed", "champ_dmplist_all_kp_win_r.RData"))
}

  
#--------------------------------End of the code--------------------------------



