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

# library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# library(IlluminaHumanMethylation450kmanifest)
# 
# ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# ann450kSub <- ann450k[, c(1:4,12:19,24:ncol(ann450k))]
# ann450kSubT<-as.data.frame(ann450kSub[,])
# 
# #remove cpgs from X & Y chromosomes
# ann450kRS_xy<-ann450kSubT %>% 
#   filter(chr %in% c('chrX','chrY'))
# 
# #remove non cpg loci
# PEG_NOOB_nors_filter <- PEG_NOOB_nors_new %>% 
#   filter(!str_detect(row.names(.), "ch|rs")) %>% 
#   filter(rownames(.) %notin% rownames(ann450kRS_xy))
# 
# save(PEG_NOOB_nors_filter, file = "PEG_NOOB_nors_filter.RData")

# Load filtered methylation matrix data
load(here("data", "methylation", "processed", "PEG_NOOB_nors_filter.RData"))


#select covariates for cases and controls

list(lb_sd_copper_wt_count, lb_sd_op_wt_count) %>%
  pmap(function(df1, df2){
    list(df1, df2) %>%
      pmap(function(data1, data2){
        data1 %>%
          dplyr::select(pegid, all_of(myvars_covar)) %>%
          left_join(data2 %>%
                      dplyr::select(pegid, count), by = "pegid") %>%
          na.omit()
      })
  }) %>%
  set_names("covar_case_process", "covar_ctrl_process") %>%
  list2env(.,envir = .GlobalEnv)


covar_case_combind <- covar_case_process[[2]] %>% 
  left_join(count_combine_op[[1]], by = "pegid") %>% 
  dplyr::select(-c(pegid, count)) %>% 
  dplyr::rename(count = total)

covar_ctrl_combind <- covar_ctrl_process[[2]] %>% 
  left_join(count_combine_op[[2]], by = "pegid") %>% 
  dplyr::select(-c(pegid, count, study)) %>% 
  dplyr::rename(count = total)

#Generate the beta matrix for cases and controls, respectively

list(lb_sd_copper_wt_count[[1]][[2]], lb_sd_copper_wt_count[[2]][[2]]) %>% 
  map(function(data){
    PEG_NOOB_nors_filter %>% 
      select(all_of(data %>% 
                      pull(sampleid))) %>%
      apply(MARGIN = 1, dec_out) %>%
      t() %>%
      as.data.frame()
  }) %>% 
  set_names("PEG_NOOB_nors_win_filter_case","PEG_NOOB_nors_win_filter_ctrl") %>%
  list2env(.,envir = .GlobalEnv)

save(PEG_NOOB_nors_win_filter_case, file = "PEG_NOOB_nors_win_filter_case.RData")
save(PEG_NOOB_nors_win_filter_ctrl, file = "PEG_NOOB_nors_win_filter_ctrl.RData")

#check rownames
table(rownames(PEG_NOOB_nors_win_filter_case)==rownames(PEG_NOOB_nors_filter))
table(rownames(PEG_NOOB_nors_win_filter_ctrl)==rownames(PEG_NOOB_nors_filter))

####run empirical Bayes linear model####
library(WGCNA)
p_count <- list(1:100000, 100001:200000, 200001:300000, 300001:349789)

# separate the beta matrix to smaller datasets
list(list(p_count, p_count),
     list(PEG_NOOB_nors_win_filter_case, PEG_NOOB_nors_win_filter_ctrl)) %>% 
  pmap(function(data1,data2){
    data1 %>% 
      map(function(list){
        data2 %>% 
          dplyr::slice(list)
      })
  }) %>% 
  set_names("methyllist_case", "methyllist_ctrl") %>% 
  list2env(.,envir = .GlobalEnv)

# run empirical Bayes linear model

list(list(methyllist_case, methyllist_ctrl),
     list(covar_case_combind, covar_ctrl_combind)) %>% 
  pmap(function(data1,data2){
    data1 %>% 
      map(function(df){
        empiricalBayesLM(data=t(df),
                         removedCovariates = data2,
                         automaticWeights = "bicov",
                         aw.maxPOutliers = 0.01)
      })
  }) %>% 
  set_names("methyllist_resid_case","methyllist_resid_ctrl") %>% 
  list2env(.,envir = .GlobalEnv)


list(methyllist_resid_case, methyllist_resid_ctrl) %>% 
  map(function(data){
    data %>% 
      map(function(df){
        as.data.frame(df$adjustedData)
      })
  }) %>% 
  set_names("adjusteddf_resid_case", "adjusteddf_resid_ctrl") %>% 
  list2env(.,envir = .GlobalEnv)


list(adjusteddf_resid_case, adjusteddf_resid_ctrl) %>% 
  map(function(data){
    data %>% 
      bind_cols() %>% 
      t() %>% 
      as.data.frame()
  }) %>% 
  set_names("combined_resid_win_filter_case", "combined_resid_win_filter_ctrl") %>% 
  list2env(.,envir = .GlobalEnv)

save(combined_resid_win_filter_case, file = "combined_resid_win_filter_case.RData")
save(combined_resid_win_filter_ctrl, file = "combined_resid_win_filter_ctrl.RData")



#--------------------------------End of the code--------------------------------


















