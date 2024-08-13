## ---------------------------
##
## Script name: get_residuals
##
## Purpose of script: To get residuals from the linear regression model
##
## Author: Yufan Gong
##
## Date Created: 2022-11-22
##
## Copyright (c) Yufan Gong, 2022
## Email: ivangong@ucla.edu
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

# DNA Methylation residuals -----------------------------------------------

## Y ~ covariates

myvars_covar <- quote_all(cd8t, cd4t, nk, mono, bcell, gran, age, female, 
                          smokers, rfvotecaucasian2, study)
#remove pdyears, a1_schyrs
#stratify by pd status??

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)

ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450kSub <- ann450k[, c(1:4,12:19,24:ncol(ann450k))]
ann450kSubT<-as.data.frame(ann450kSub[,])

#remove cpgs from X & Y chromosomes
ann450kRS_xy<-ann450kSubT %>% 
  filter(chr %in% c('chrX','chrY'))

#remove non cpg loci
PEG_NOOB_nors_filter <- PEG_NOOB_nors_new %>% 
  filter(!str_detect(row.names(.), "ch|rs")) %>% 
  filter(rownames(.) %notin% rownames(ann450kRS_xy))

save(PEG_NOOB_nors_filter, file = "PEG_NOOB_nors_filter.RData")

#limit to PD case only
list(
  list(c_lb_sd_case_wt_10, r_lb_sd_case_wt_10),
  list(id_remove_c, id_remove_r)
) %>% 
  pmap(function(data1, data2){
    data1 %>% 
      filter(pegid %notin% data2 & !is.na(pegid)) %>% 
      pull(sampleid) 
  }) %>% 
  set_names("sampleid_pd_c","sampleid_pd_r") %>% 
  list2env(.,envir = .GlobalEnv)

#limit to control only
list(
  list(c_lb_sd_control_wt_10_count, r_lb_sd_control_wt_10_count),
  list(id_remove_c, id_remove_r)
) %>% 
  pmap(function(data1, data2){
    data1 %>% 
      filter(pegid %notin% data2 & !is.na(pegid)) %>% 
      pull(sampleid) 
  }) %>% 
  set_names("sampleid_ctrl_c","sampleid_ctrl_r") %>% 
  list2env(.,envir = .GlobalEnv)


#select covariates
list(
  list(c_lb_sd_case_wt_10, r_lb_sd_case_wt_10),
  list(id_remove_c, id_remove_r)
) %>% 
  pmap(function(data1, data2){
    data1 %>% 
      filter(pegid %notin% data2 & !is.na(pegid)) %>% 
      filter(sampleid %in% sampleid_pd_r) %>% 
      select(all_of(myvars_covar)) %>% 
      na.omit()
  }) %>% 
  set_names("covar_c","covar_r") %>% 
  list2env(.,envir = .GlobalEnv)


list(
  list(c_lb_sd_control_wt_10_count, r_lb_sd_control_wt_10_count),
  list(id_remove_c, id_remove_r)
) %>% 
  pmap(function(data1, data2){
    data1 %>% 
      filter(pegid %notin% data2 & !is.na(pegid)) %>% 
      filter(sampleid %in% sampleid_ctrl_r) %>% 
      select(all_of(myvars_covar)) %>% 
      select(-study) %>% 
      na.omit()
  }) %>% 
  set_names("covar_ctrl_c","covar_ctrl_r") %>% 
  list2env(.,envir = .GlobalEnv)

#winsorize the beta matrix for occupational and residential, respectively
list(
  list(c_lb_sd_case_wt_10, r_lb_sd_case_wt_10),
  list(id_remove_c, id_remove_r)
) %>% 
  pmap(function(data1, data2){
    PEG_NOOB_nors_filter %>% 
      select(all_of(data1 %>% 
                      filter(pegid %notin% data2 & !is.na(pegid)) %>% 
                      filter(sampleid %in% sampleid_pd) %>% 
                      pull(sampleid)
      )) %>%
      # apply(MARGIN = 1, dec_out) %>%
      t() %>%
      as.data.frame()
  }) %>% 
  set_names("PEG_NOOB_nors_win_filter_pd_c","PEG_NOOB_nors_win_filter_pd_r") %>% 
  list2env(.,envir = .GlobalEnv)


list(
  list(c_lb_sd_control_wt_10_count, r_lb_sd_control_wt_10_count),
  list(id_remove_c, id_remove_r)
) %>% 
  pmap(function(data1, data2){
    PEG_NOOB_nors_filter %>% 
      select(all_of(data1 %>% 
                      filter(pegid %notin% data2 & !is.na(pegid)) %>% 
                      filter(sampleid %in% sampleid_ctrl_r) %>% 
                      pull(sampleid)
      )) %>%
      apply(MARGIN = 1, dec_out) %>%
      t() %>%
      as.data.frame()
  }) %>% 
  set_names("PEG_NOOB_nors_win_filter_ctrl_c","PEG_NOOB_nors_win_filter_ctrl_r") %>% 
  list2env(.,envir = .GlobalEnv)

save(PEG_NOOB_nors_win_filter_ctrl_c, file = "PEG_NOOB_nors_win_filter_ctrl_c.RData")
save(PEG_NOOB_nors_win_filter_ctrl_r, file = "PEG_NOOB_nors_win_filter_ctrl_r.RData")


#check rownames
table(rownames(PEG_NOOB_nors_win_filter_ctrl_c)==rownames(PEG_NOOB_nors_filter))

save(PEG_NOOB_nors_win_filter_pd_c, file = "PEG_NOOB_nors_win_filter_pd_c.RData")
save(PEG_NOOB_nors_win_filter_pd_r, file = "PEG_NOOB_nors_win_filter_pd_r.RData")

# load winsorized beta-matrix
# list.dirs(here(),recursive = FALSE) %>%
#   list.files("\\.RData$", full.names = TRUE, recursive = T) %>%
#   grep("filter_pd",., 
#        value=TRUE, ignore.case = TRUE) %>%
#   keep(~!str_detect(.x,"combine|dmplist|old|more|90")) %>%
#   map(.,load,.GlobalEnv)

library(WGCNA)
p_count <- list(1:100000, 100001:200000, 200001:300000, 300001:349789)

# separate the beta matrix to smaller datasets
list(list(p_count,p_count),
     list(PEG_NOOB_nors_win_filter_pd_c,PEG_NOOB_nors_win_filter_pd_r)) %>% 
  pmap(function(data1,data2){
    data1 %>% 
      map(function(list){
        data2 %>% 
          dplyr::slice(list)
      })
  }) %>% 
  set_names("methyllist_win_c","methyllist_win_r") %>% 
  list2env(.,envir = .GlobalEnv)

list(list(p_count,p_count),
     list(PEG_NOOB_nors_win_filter_ctrl_c,PEG_NOOB_nors_win_filter_ctrl_r)) %>% 
  pmap(function(data1,data2){
    data1 %>% 
      map(function(list){
        data2 %>% 
          dplyr::slice(list)
      })
  }) %>% 
  set_names("methyllist_win_c","methyllist_win_r") %>% 
  list2env(.,envir = .GlobalEnv)


# run empirical Bayes linear model
#for cases
list(list(methyllist_win_c,methyllist_win_r),
     list(covar_c,covar_r)) %>% 
  pmap(function(data1,data2){
    data1 %>% 
      map(function(df){
        empiricalBayesLM(data=t(df),
                         removedCovariates = data2,
                         automaticWeights = "bicov",
                         aw.maxPOutliers = 0.01)
      })
  }) %>% 
  set_names("methyllist_resid_c","methyllist_resid_r") %>% 
  list2env(.,envir = .GlobalEnv)

#for controls
list(list(methyllist_win_c,methyllist_win_r),
     list(covar_ctrl_c,covar_ctrl_r)) %>% 
  pmap(function(data1,data2){
    data1 %>% 
      map(function(df){
        empiricalBayesLM(data=t(df),
                         removedCovariates = data2,
                         automaticWeights = "bicov",
                         aw.maxPOutliers = 0.01)
      })
  }) %>% 
  set_names("methyllist_resid_c","methyllist_resid_r") %>% 
  list2env(.,envir = .GlobalEnv)

list(methyllist_resid_c,methyllist_resid_r) %>% 
  map(function(data){
    data %>% 
      map(function(df){
        as.data.frame(df$adjustedData)
      })
  }) %>% 
  set_names("adjusteddf_resid_c","adjusteddf_resid_r") %>% 
  list2env(.,envir = .GlobalEnv)

# list(adjusteddf_resid_c,adjusteddf_resid_r) %>% 
#   map(function(data){
#     data %>% 
#       bind_cols() %>% 
#       t() %>% 
#       as.data.frame()
#   }) %>% 
#   set_names("combined_resid_filter_pd_c","combined_resid_filter_pd_r") %>% 
#   list2env(.,envir = .GlobalEnv)

#for controls
list(adjusteddf_resid_c,adjusteddf_resid_r) %>% 
  map(function(data){
    data %>% 
      bind_cols() %>% 
      t() %>% 
      as.data.frame()
  }) %>% 
  set_names("combined_resid_filter_ctrl_c","combined_resid_filter_ctrl_r") %>% 
  list2env(.,envir = .GlobalEnv)

save(combined_resid_filter_ctrl_c, file = "combined_resid_filter_ctrl_95percentile_out0.01_c.RData")
save(combined_resid_filter_ctrl_r, file = "combined_resid_filter_ctrl_95percentile_out0.01_r.RData")

save(combined_resid_filter_pd_c, file = "combined_resid_filter_pd_95percentile_out0.01_c.RData")
save(combined_resid_filter_pd_r, file = "combined_resid_filter_pd_95percentile_out0.01_r.RData")



#--------------------------------End of the code--------------------------------


















