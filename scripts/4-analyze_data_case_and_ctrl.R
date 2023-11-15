## ---------------------------
##
## Script name: 
##
## Purpose of script: 
##
## Author: Yufan Gong
##
## Date Created: 2023-02-12
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


#total population: stratified by study wave
myvar <- quote_all(pegid,sampleid,age, pd_new,female,
                   smokers, a1_schyrs, ethnicity, county,
                   pdstudystudy,meanmethbysample)


list(c_lb_sd_case_wt_10_count,
     r_lb_sd_case_wt_10_count,
     c_lb_sd_control_wt_10_count,
     r_lb_sd_control_wt_10_count) %>% 
  map(function(data){
    data %>% 
      dplyr::select(all_of(myvar))
  }) %>% 
  rbindlist() %>% 
  set_variable_labels(
    pegid = "PEGID",
    age = "Age",
    pd_new = "PD status",
    female = "Sex",
    #race_new = "Race/Ethnicity",
    ethnicity = "Ethnicity",
    smokers = "Smoking status",
    pdstudystudy = "PEG study",
    meanmethbysample = "Mean methylation",
    a1_schyrs = "School years",
    county = "County"
  ) %>% 
  distinct() %>% 
  #filter(pegid %in% datSampleSteve$externaldnacode) %>% 
  dplyr::select(-c(sampleid,pegid, pd_new)) %>%
  tbl_summary(missing = "no",
              #type=list(c(female) ~ "categorical"),
              statistic = list(all_continuous() ~ "{mean} ({sd})",
                               all_categorical() ~ "{n} ({p}%)"),
              digits = list(all_continuous() ~ 1,
                            all_categorical() ~ 1)) %>%
  modify_header(label = "**Characteristics**") %>%
  # modify_spanning_header(starts_with("stat_") ~ "**PD status**") %>%
  modify_caption("Table 1. Demographic characteristics of PEG cases (N = {N})") %>%
  bold_labels() %>% 
  table1()


#heavy metal chemical use in total

list(exp_window_address_case_c, exp_window_address_case_r) %>% 
  map(function(data){
    data %>% 
      filter(chemcode %in% metal_name$chemcode) %>%
      mutate_at(vars(sum_total_lbs), extreme_remove_percentile_win) %>% 
      group_by(pegid,year) %>% 
      summarise(sum_total_lbs = sum(sum_total_lbs),.groups = "keep") %>% 
      inner_join(pest_case_methylation, by = "pegid") %>%
      replace_na(list(sum_total_lbs = 0))
  }) %>% 
  set_names("outcheck_c","outcheck_r") %>% 
  list2env(.,envir = .GlobalEnv)

list(list(outcheck_c,outcheck_r),
     c("Occupational","Residential"))%>% 
  pmap(function(data1,data2){
    data1 %>% 
      drop_na() %>% 
      mutate(location = data2)
  }) %>% 
  rbindlist() %>% 
  filter(year <= 2008) %>% 
  group_by(year, location) %>% 
  summarise(y=mean(sum_total_lbs, na.rm=T),.groups="keep") %>% 
  #change to average per person, including unexposed people
  ggplot(aes(x=year, y=y, group=location, color = factor(location))) + 
  annotate("rect", fill = "azure2", alpha = 0.5,
           xmin = 1989, xmax = 2010,
           ymin = -Inf, ymax = Inf) +
  labs(x = "Years",
       y = "Chemical use (lbs)",
       color = "Location") +
  #geom_point()+
  geom_line(linewidth=1) +
  scale_color_colorblind() +
  # scale_color_ordinal(labels=c("Without PD", "With PD")) +
  geom_vline(xintercept = 1989, lty=2) +
  theme_classic()+
  theme( plot.title = element_text(hjust = 0.5, size=20,face="bold"),
         axis.text = element_text(size = 15),
         axis.title=element_text(size=20,face="bold"),
         legend.title = element_text(size = 14),
         legend.text = element_text(size = 14),
         strip.text = element_text(size=13)) 

#heavy metal chemical use: specific


list(exp_window_address_case_c, exp_window_address_case_r) %>% 
  map(function(data){
    data %>% 
      filter(chemcode %in% metal_name$chemcode) %>%
      # mutate_at(vars(sum_total_lbs), extreme_remove_percentile_win) %>% 
      group_by(chemcode) %>% 
      group_split() %>% 
      map(function(df){
        df %>% 
          mutate_at(vars(sum_total_lbs), extreme_remove_percentile_win) %>% 
          group_by(pegid,year) %>% 
          summarise(sum_total_lbs = sum(sum_total_lbs),.groups = "keep") %>% 
          inner_join(pest_case_methylation, by = "pegid") %>%
          replace_na(list(sum_total_lbs = 0))
      })
  }) %>% 
  set_names("outcheck_list_c","outcheck_list_r") %>% 
  list2env(.,envir = .GlobalEnv)

plotlist <- list(outcheck_list_c, outcheck_list_r) %>% 
  pmap(function(list1, list2){
    list(
      list(list1, list2),
      c("Occupational","Residential")
    ) %>% 
      pmap(function(data1,data2){
        data1 %>% 
          drop_na() %>% 
          mutate(location = data2)
      }) %>% 
      rbindlist() %>% 
      filter(year <= 2008) %>% 
      group_by(year, location) %>% 
      summarise(y=mean(sum_total_lbs, na.rm=T),.groups="keep") %>% 
      #change to average per person, including unexposed people
      ggplot(aes(x=year, y=y, group=location, color = factor(location))) + 
      annotate("rect", fill = "azure2", alpha = 0.5,
               xmin = 1989, xmax = 2010,
               ymin = -Inf, ymax = Inf) +
      labs(x = "Years",
           y = "Chemical use (lbs)",
           color = "Location") +
      #geom_point()+
      geom_line(linewidth=1) +
      scale_color_colorblind() +
      # scale_color_ordinal(labels=c("Without PD", "With PD")) +
      geom_vline(xintercept = 1989, lty=2) +
      theme_classic()+
      theme( plot.title = element_text(hjust = 0.5, size=20,face="bold"),
             axis.text = element_text(size = 15),
             axis.title=element_text(size=20,face="bold"),
             legend.title = element_text(size = 14),
             legend.text = element_text(size = 14),
             strip.text = element_text(size=13)) 
  })




png(file=here("figures","heavy metal exposure trend_fitered.png"), 
    width = 1920, height = 1920)
ggarrange(plotlist = plotlist,
          labels = metal_name$chemname,
          ncol = 3, nrow = 6)
dev.off()




# DMP --------------------------------------------------------------------

#Differentially methylated positions (DMP)
# Load meffil and set how many cores to use for parallelization

# source(here("scripts", "mclapply.hack.R"))
source(here::here("scripts", "meffil_fixed.R"))
source(here::here("scripts", "meffil_report_fixed.R"))


#Method 1: ChAMP
library(ChAMP)
library(limma)

# list(list(c_lb_sd_case_wt_10_new, r_lb_sd_case_wt_10_new),
#      list(combined_resid_filter_pd_c, combined_resid_filter_pd_r)) %>% 
#   pmap(function(data1,data2){
#     data1 %>% 
#       select(all_of(metal_name$chemcode)) %>%
#       # select(chem155, chem156, chem60) %>% 
#       # select(-c(chem151, chem158, chem161)) %>% 
#       map(~champ.DMP(beta = data2, pheno = .x, 
#                      compare.group = NULL, 
#                      adjPVal = 0.05, adjust.method = "BH", 
#                      arraytype = "450K")[[1]])
#   }) %>% 
#   set_names("champ_dmplist_case_heavymetal_signew_c",
#             "champ_dmplist_case_heavymetal_signew_r") %>% 
#   list2env(.,envir = .GlobalEnv)

# using the quartile pesticide dataset
dmp_quartile_c <- c_lb_sd_case_wt_10_quantile %>% 
  # select(all_of(metal_name$chemcode)) %>% 
  select(chem155, chem283, chem714) %>% 
  # select(-c(chem151, chem158, chem161)) %>% 
  map(~champ.DMP(beta = combined_resid_filter_pd_c, pheno = .x, 
                 compare.group = NULL, 
                 adjPVal = 0.05, adjust.method = "BH", 
                 arraytype = "450K")[[1]])


# using the evernever pesticide dataset
dmp_evernever_c <- c_lb_sd_case_wt_10_evernever %>% 
  # select(all_of(metal_name$chemcode)) %>% 
  select(chem155, chem283) %>% 
  # select(-c(chem151, chem158, chem161)) %>% 
  map(~champ.DMP(beta = combined_resid_filter_pd_c, pheno = .x, 
                 compare.group = NULL, 
                 adjPVal = 0.05, adjust.method = "BH", 
                 arraytype = "450K")[[1]])


# using the z-transformed pesticide dataset
dmp_ztrans_c <- c_lb_sd_case_wt_10_ztrans %>% 
  # select(all_of(metal_name$chemcode)) %>% 
  select(chem283) %>% 
  # select(-c(chem151, chem158, chem161)) %>% 
  map(~champ.DMP(beta = combined_resid_filter_pd_c, pheno = .x, 
                 compare.group = NULL, 
                 adjPVal = 0.05, adjust.method = "BH", 
                 arraytype = "450K")[[1]])

dmp_ztrans_r <- r_lb_sd_case_wt_10_ztrans %>% 
  # select(all_of(metal_name$chemcode)) %>% 
  select(chem1876, chem34) %>% 
  # select(-c(chem151, chem158, chem161)) %>% 
  map(~champ.DMP(beta = combined_resid_filter_pd_r, pheno = .x, 
                 compare.group = NULL, 
                 adjPVal = 0.05, adjust.method = "BH", 
                 arraytype = "450K")[[1]])

#using count data
list(
  list(c_lb_sd_case_wt_10_count, r_lb_sd_case_wt_10_count),
  list(c_lb_sd_control_wt_10_count, r_lb_sd_control_wt_10_count)
) %>% 
  map(function(datalist){
    datalist %>% 
      map(function(data){
        data %>% 
          dplyr::select(pegid, metal_count, copper_count)
      }) %>% 
      purrr::reduce(full_join, by = "pegid") %>% 
      mutate_all(~replace(., is.na(.), 0)) %>% 
      transmute(pegid = pegid,
                  total = metal_count.x + metal_count.y,
                copper = copper_count.x + copper_count.y)
  }) %>% 
  set_names("count_combine_case", "count_combine_ctrl") %>% 
  list2env(.GlobalEnv)

count_combine_total <- list(count_combine_ctrl, count_combine_case) %>% 
  bind_rows()

hist(count_combine_case$copper)
hist(count_combine_ctrl$total)
hist(count_combine_total$total)
skim(count_combine_case$copper)
skim(count_combine_ctrl$total)
skim(count_combine_total$total)

# plotlist <- c_lb_sd_case_wt_10_quantile %>% 
#   select(starts_with("chem")) %>% 
#   map(function(data){
#     data %>% 
#     ggplot(aes(x = ., fill = .)) +
#       scale_fill_colorblind() +
#       geom_bar()
#   })

list(
  list(count_combine_case, count_combine_ctrl),
  list(combined_resid_filter_pd_r, combined_resid_filter_ctrl_r)
) %>% 
  pmap(function(data1, data2){
    data1 %>% 
      select(total) %>% 
      map(~champ.DMP(beta = data2, pheno = .x, 
                     compare.group = NULL, 
                     adjPVal = 0.05, adjust.method = "BH", 
                     arraytype = "450K")[[1]])
  }) %>% 
  set_names("dmp_count_case", "dmp_count_ctrl") %>% 
  list2env(.GlobalEnv)


save(dmp_count_case, file = "dmp_count_case.RData")
save(dmp_count_ctrl, file = "dmp_count_ctrl.RData")

save(dmp_quartile_c, 
     file = "dmp_quartile_c.RData")
save(dmp_evernever_c, 
     file = "dmp_evernever_c.RData")
save(dmp_ztrans_c, 
     file = "dmp_ztrans_c.RData")
save(dmp_ztrans_r, 
     file = "dmp_ztrans_r.RData")

DMP.GUI(DMP = dmp_count_case[[1]], 
        beta = as.matrix(combined_resid_filter_pd_r), 
        pheno = count_combine_case$total,
        cutgroupnumber = 2)

DMP.GUI()
DMP.GUI(DMP = dmp_count_ctrl[[1]], 
        beta = as.matrix(combined_resid_filter_ctrl_r), 
        pheno = count_combine_ctrl$total,
        cutgroupnumber = 2)

# champ.GSEA(beta=combined_resid_filter_pd_r,
#            DMP=dmp_r_new[[1]],
#            DMR=dmr_r[[1]],
#            CpGlist=NULL,
#            Genelist=NULL,
#            pheno=r_lb_sd_case_wt_10_new$chem1876,
#            method="fisher",
#            arraytype="450K",
#            Rplot=TRUE,
#            adjPval=0.05,
#            cores=1)

#Method 2: meffil
library(meffil)
myvars_ewas <- quote_all(cd8t, cd4t, nk, mono, bcell, gran, age, female,smokers,
                         rfvotecaucasian2, study)

# myvars_ewas <- quote_all(cd4t, gran, age, female,
#                          pdstudynumberyearswithpdatbloodd,
#                          rfvotecaucasian2, pd_new)


#select covariates
list(
  list(c_lb_sd_case_wt_10_count, r_lb_sd_case_wt_10_count),
  list(c_lb_sd_control_wt_10_count, r_lb_sd_control_wt_10_count)
) %>% 
  map(function(datalist){
    datalist %>% 
      map(function(data){
        data %>% 
          dplyr::select(all_of(myvars_ewas)) %>% 
          # rename(pdyears = pdstudynumberyearswithpdatbloodd) %>% 
          # mutate(pdyears = ifelse(is.na(pdyears), 
          #                         mean(pdyears, na.rm=TRUE), pdyears)) %>% 
          # replace_na(list(pdyears = 0)) %>% 
          na.omit()
        #replace_na(list(c11_depression = 0,gds1_5 = 0)) %>% 
      })
  }) %>% 
  set_names("covar_case","covar_ctrl") %>% 
  list2env(.,envir = .GlobalEnv)


covar_case %>% 
  set_names("covar_case_c", "covar_case_r") %>% 
  list2env(.,envir = .GlobalEnv)

covar_ctrl %>% 
  map(function(data){
    data %>% 
      dplyr::select(-study)
  }) %>% 
  set_names("covar_ctrl_c", "covar_ctrl_r") %>% 
  list2env(.,envir = .GlobalEnv)


covar_total <- list(
  list(covar_ctrl[[2]], covar_case_r),
  c("ctrl", "case")
) %>% 
  pmap(function(data1, data2){
    data1 %>% 
      mutate(pd = data2)
  }) %>% 
  bind_rows()

ewas.parameters <- meffil.ewas.parameters(
  sig.threshold=1e-6,  ## EWAS p-value threshold
  max.plots=20, ## plot at most 20 CpG sites
  qq.inflation.method="median",  ## measure inflation using median
  model="all") ## select default EWAS model; 


list(list(count_combine_case, count_combine_ctrl, count_combine_total),
     list(PEG_NOOB_nors_win_filter_pd_r, PEG_NOOB_nors_win_filter_ctrl_r, 
          peg_noob_nors_win_total),
     list(covar_case_r, covar_ctrl_r, covar_total)) %>% 
  pmap(function(data1,data2,data3){
    data1 %>% 
      dplyr::select(total, copper) %>% 
      map(~meffil.ewas( beta=as.matrix(data2), 
                        variable=.x, covariates = data3, 
                        batch = NULL, weights = NULL,  cell.counts = NULL,
                        isva = F, sva = F, smartsva = F, n.sv = NULL, 
                        winsorize.pct = NA, robust = TRUE,
                        rlm = FALSE, outlier.iqr.factor = NA,  featureset = NA, 
                        random.seed = 20230922, lmfit.safer = F, verbose = T))
  }) %>% 
  set_names("meffil_count_case","meffil_count_ctrl", "meffil_count_total") %>% 
  list2env(.,envir = .GlobalEnv)

test <- meffil_count_total$copper$analyses$all$table %>% 
  filter(-log10(p.value) > 6)

save(meffil_count_case, file = "meffil_count_case.RData")
save(meffil_count_ctrl, file = "meffil_count_ctrl.RData")
save(meffil_count_total, file = "meffil_count_total.RData")

list(list(meffil_count_case, meffil_count_ctrl, meffil_count_total),
     list(PEG_NOOB_nors_win_filter_pd_r, PEG_NOOB_nors_win_filter_ctrl_r, 
          peg_noob_nors_win_total)) %>% 
  pmap(function(meffil_list,data){
    meffil_list %>% 
      map(~ meffil.ewas.summary_fix(.x,data,
                                    parameters = ewas.parameters))
  }) %>% 
  set_names("ewas.summary_count_case", "ewas.summary_count_ctrl", 
            "ewas.summary_count_total") %>% 
  list2env(.,envir = .GlobalEnv)

list(list(ewas.summary_count_case,
          names(ewas.summary_count_case), "case"),
     list(ewas.summary_count_ctrl,
          names(ewas.summary_count_ctrl), "ctrl"),
     list(ewas.summary_count_total,
          names(ewas.summary_count_total), "total")) %>%
  map(function(datalist){
    datalist %>%
      pmap(function(data1, data2, data3){
        meffil.report(
          data1,
          output.file = here::here(
            "reports", data3, paste(data2, data3, ".html", sep = "_")))
      })
  })

# DMR analysis ------------------------------------------------------------


# using the quartile pesticide dataset
dmr_champ_chem155_quartile_c <- c_lb_sd_case_wt_10_quantile %>% 
  select(chem155) %>% 
  map(~champ.DMR(beta = as.matrix(combined_resid_filter_pd_c), 
                 pheno = .x, compare.group = NULL, 
                 arraytype = "450K", method = "Bumphunter", minProbes = 0, 
                 adjPvalDmr = 0.05, cores = 5, maxGap = 300, cutoff = NULL, 
                 pickCutoff = TRUE, smooth = TRUE, 
                 smoothFunction = loessByCluster, useWeights = FALSE, 
                 permutations = NULL, B = 250, nullMethod = "bootstrap", 
                 meanLassoRadius = 375, minDmrSep = 1000, minDmrSize = 50, 
                 adjPvalProbe = 0.05, Rplot = T, PDFplot = T, 
                 resultsDir = here("dmr"), 
                 rmSNPCH = T, fdr = 0.05, dist = 2, 
                 mafcut = 0.05, lambda = 1000, C = 2))


# using the evernever pesticide dataset
dmr_champ_chem155_evernever_c <- c_lb_sd_case_wt_10_evernever %>% 
  select(chem155) %>% 
  map(~champ.DMR(beta = as.matrix(combined_resid_filter_pd_c), 
                 pheno = .x, compare.group = NULL, 
                 arraytype = "450K", method = "Bumphunter", minProbes = 0, 
                 adjPvalDmr = 0.05, cores = 5, maxGap = 300, cutoff = NULL, 
                 pickCutoff = TRUE, smooth = TRUE, 
                 smoothFunction = loessByCluster, useWeights = FALSE, 
                 permutations = NULL, B = 250, nullMethod = "bootstrap", 
                 meanLassoRadius = 375, minDmrSep = 1000, minDmrSize = 50, 
                 adjPvalProbe = 0.05, Rplot = T, PDFplot = T, 
                 resultsDir = here("dmr"), 
                 rmSNPCH = T, fdr = 0.05, dist = 2, 
                 mafcut = 0.05, lambda = 1000, C = 2))

  
save(dmr_champ_chem155_evernever_c, file = "dmr_champ_chem155_evernever_c.RData")
save(dmr_champ_chem283_evernever_c, file = "dmr_champ_chem283_evernever_c.RData")
save(dmr_champ_chem714_evernever_c, file = "dmr_champ_chem714_evernever_c.RData")

dmr_evernever_c <- list(chem155 = dmr_champ_chem155_evernever_c[[1]], 
                        chem283 = dmr_champ_chem283_evernever_c[[1]], 
                        chem714 = dmr_champ_chem714_evernever_c[[1]])

dmr_highlow_c <- list(chem155 = dmr_champ_chem155_c[[1]], 
                      chem283 = dmr_champ_chem283_c[1], 
                      chem714 = dmr_champ_chem714_c[[1]])

save(dmr_evernever_c, file = "dmr_evernever_c.RData")
save(dmr_highlow_c, file = "dmr_highlow_c.RData")

DMR.GUI(DMR  = dmr_evernever_c[[3]], 
        beta = as.matrix(combined_resid_filter_pd_c), 
        pheno = c_lb_sd_case_wt_10_evernever$chem714)


#Method 2 DRMcate

# list(list(count_combine_case, count_combine_ctrl),
#      list(PEG_NOOB_nors_win_filter_pd_r, PEG_NOOB_nors_win_filter_ctrl_r)) %>% 
#   pmap(function(data1,data2){
#     data1 %>% 
#       dplyr::select(total) %>% 
#       map(~
#             cpg.annotate(datatype="array",fdr = 0.05, 
#                          as.matrix(data2),
#                          what="Beta",arraytype = "450K",
#                          design=model.matrix(~ .x, data1),
#                          coef=2, 
#                          analysis.type="differential",
#                          annotation=c(array = "IlluminaHumanMethylation450k"
#                                       , annotation = "ilmn12.hg19")) %>% 
#             dmrcate(lambda=1000, C=2) %>% 
#             extractRanges(genome = "hg19")
#       )
#   }) %>% 
#   set_names("dmrlist_case","dmrlist_ctrl") %>% 
#   list2env(.,envir = .GlobalEnv)




# GESA analysis -----------------------------------------------------------

list(meffil_count_case, meffil_count_ctrl, meffil_count_total) %>% 
  map(function(meffil_list){
    meffil_list %>% 
      map(function(data){
        data$analyses$all$table %>% 
          # slice_min(p.value, n=2500) %>% 
          mutate(cpgs = row.names(.))
      })
  }) %>% 
  set_names("metal_meffil_case", "metal_meffil_control", 
            "metal_meffil_total") %>% 
  list2env(.,envir = .GlobalEnv)

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)

ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

list(metal_meffil_case, metal_meffil_control, metal_meffil_total) %>% 
  map(function(data){
    data %>% 
      map(function(df){
        as.data.frame(ann450k[match(df$cpgs,ann450k$Name), 
                              c(1:4,12:19,24:ncol(ann450k))])
      })
  }) %>% 
  set_names("ann450ksubt_case","ann450ksubt_ctrl", "ann450ksubt_total") %>% 
  list2env(.,envir = .GlobalEnv)

list(list(metal_meffil_case, ann450ksubt_case),
     list(metal_meffil_control, ann450ksubt_ctrl),
     list(metal_meffil_total, ann450ksubt_total)) %>% 
  map(function(data){
    data %>% 
      pmap(function(data1,data2){
        cbind(data1,data2)
      })
  }) %>% 
  set_names("ewas_annot_case", "ewas_annot_ctrl", "ewas_annot_total") %>% 
  list2env(.,envir = .GlobalEnv)


list(ewas_annot_case, ewas_annot_ctrl, ewas_annot_total) %>% 
  map(function(plist){
    plist %>% 
      map(~filter(.x, fdr < 0.05))
  }) %>% 
  set_names("ewas_annot_new_case","ewas_annot_new_ctrl", 
            "ewas_annot_new_total") %>% 
  list2env(.,envir = .GlobalEnv)

library(methylGSA)

meta.cpg_case <- ewas_annot_new_case %>% 
  map(function(data){
    data %>% 
      pull(p.value) %>% 
      set_names(data$cpgs)
  })

meta.cpg_total <- ewas_annot_new_total %>% 
  map(function(data){
    data %>% 
      pull(p.value) %>% 
      set_names(data$cpgs)
  })

meta.cpg_case %>% 
  set_names("meta.cpg_case_total", "meta.cpg_case_copper") %>% 
  list2env(.,envir = .GlobalEnv)

meta.cpg_total %>% 
  set_names("meta.cpg_total_total", "meta.cpg_total_copper") %>% 
  list2env(.,envir = .GlobalEnv)


#function 2: methylRRA
gsea_case_total <- methylRRA(cpg.pval = meta.cpg_case_total, method = "GSEA")
gsea_case_copper <- methylRRA(cpg.pval = meta.cpg_case_copper, method = "GSEA")
gsea_total_copper <- methylRRA(cpg.pval = meta.cpg_total_copper, method = "GSEA")

barplot(gsea_case_total, num = 10, colorby = "pvalue")
barplot(gsea_case_copper, num = 10, colorby = "pvalue")
barplot(gsea_total_copper, num = 10, colorby = "pvalue")

# Pathway analysis --------------------------------------------------------

library(PANTHER.db)
pthOrganisms(PANTHER.db) <- "HUMAN"
PANTHER.db


go_ids <- keys(PANTHER.db,keytype="PATHWAY_ID")
cols <- "ENTREZ"
panther <- mapIds(PANTHER.db, keys=go_ids, column=cols, 
                  as.numeric(cols), keytype="PATHWAY_ID", multiVals="list")
lengths(panther)

#ALL CpG sites
#GSEA for panther pathways
methylglm_pan_case <- methylglm(cpg.pval = meta.cpg_case_copper, 
                                GS.list = panther, GS.idtype = "ENTREZID")

methylglm_pan_total <- methylglm(cpg.pval = meta.cpg_total_copper, 
                                GS.list = panther, GS.idtype = "ENTREZID")

#get pathway names
cols <- "PATHWAY_TERM"
res.p <- mapIds(PANTHER.db, keys=go_ids, column=cols, 
                as.numeric(cols), keytype="PATHWAY_ID", multiVals="list")

panther.names <- res.p %>% 
  bind_rows() %>% 
  gather() %>% 
  transmute(
    ID = key,
    pathway = value
  )

methylglm_path_case <- methylglm_pan_case %>% 
  left_join(panther.names, by = "ID")

methylglm_path_total <- methylglm_pan_total %>% 
  left_join(panther.names, by = "ID")

# Visulization ------------------------------------------------------------

library(ggtext)

list(meffil_count_case, meffil_count_ctrl, meffil_count_total) %>% 
  map(function(meffil_list){
    meffil_list %>% 
      #keep(., map_lgl(., ~ nrow(.x)>=300)) %>% 
      map(function(data){
        datanew <- data$analyses$all$table %>% 
          dplyr::rename(chr = chromosome)
        
        chromosomes <- paste("chr", c(1:22, "X","Y"), sep="")
        chromosomes <- intersect(chromosomes, datanew$chr)
        chromosome.lengths <- sapply(chromosomes, function(chromosome)
          max(datanew$position[which(datanew$chr == chromosome)]))
        
        
        chromosome.lengths <- as.numeric(chromosome.lengths)
        chromosome.starts <- c(1,cumsum(chromosome.lengths)+1)
        names(chromosome.starts) <- c(chromosomes, "NA")
        datanew$global <- datanew$position + 
          chromosome.starts[datanew$chr] - 1
        datanew
      })
  }) %>% 
  set_names("metal_list_case","metal_list_ctrl", "metal_list_total") %>% 
  list2env(.,envir = .GlobalEnv)

list(metal_list_case, metal_list_ctrl, metal_list_total) %>% 
  map(function(metal_list){
    metal_list %>% 
      map(function(data){
        data %>% 
          group_by(chr) %>% 
          dplyr::summarize(center = mean(global)) %>% 
          ungroup()
      })
  }) %>% 
  set_names("axis_set_case", "axis_set_ctrl", "axis_set_total") %>% 
  list2env(.,envir = .GlobalEnv)


list(metal_list_case, metal_list_ctrl, metal_list_total) %>% 
  map(function(metal_list){
    metal_list %>% 
      map(function(data){
        data %>% 
          filter(p.value == min(p.value)) %>% 
          mutate(ylim = abs(floor(log10(p.value))) + 1) %>% 
          pull(ylim)
      })
  }) %>% 
  set_names("ylim_case","ylim_ctrl", "ylim_total") %>% 
  list2env(.,envir = .GlobalEnv)


list(list(metal_list_case, axis_set_case, ylim_case),
     list(metal_list_ctrl, axis_set_ctrl, ylim_ctrl),
     list(metal_list_total, axis_set_total, ylim_total)) %>% 
  map(function(datalist){
    datalist %>% 
      pmap(function(pest, axis, ylim){
        pest %>% 
          mutate(sig = if_else(-log10(p.value) > 6, 1, 0)) %>% 
          ggplot(aes(x = global, y = -log10(p.value), 
                     color = as_factor(sig), size = -log10(p.value))) +
          scale_color_manual(values = c("1" = "red", "0" = "black")) +
          geom_hline(yintercept = -log10(10e-7), color = "grey40", linetype = "dashed") + 
          geom_point(alpha = 0.75) +
          scale_x_continuous(label = axis$chr, breaks = axis$center) +
          scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
          #scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis$chr)))) +
          scale_size_continuous(range = c(0.5,3)) +
          labs(x = "chromosome", 
               y = "-log<sub>10</sub>(p)",
               color = "p-value") + 
          theme_minimal() +
          theme( 
            legend.position = "none",
            panel.border = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            axis.title.y = element_markdown(),
            axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
          )
      })
  }) %>% 
  set_names("manhattanlist_case","manhattanlist_ctrl","manhattanlist_total") %>% 
  list2env(.,envir = .GlobalEnv)

plotlist <- list(manhattanlist_case[[2]], 
                 manhattanlist_ctrl[[2]], 
                 manhattanlist_total[[2]])

plot_label <- c("case", "control", "total")

#print all plots in one figure
png(file=here::here("figures","manhattan plots for cases and controls_copper.png"), 
    width = 1920, height = 1080)
ggarrange(plotlist = plotlist,
          labels = plot_label,
          ncol = 3, nrow = 1)
dev.off()


#scatter plot

list(quote_all(total, copper),
     quote_all(cg11208322, cg11208322))%>% 
  pmap(function(chem,cpg){
    test <- meffil_count_case[[chem]]$analyses$all$table
    test2 <- t(PEG_NOOB_nors_win_filter_pd_r %>% 
                 filter(rownames(.) %in% rownames(test))) %>% 
      as.data.frame()
    test3 <- cbind(count_combine_case %>% 
                     dplyr::select(chem), test2)
    test3 %>% 
      dplyr::select(chem,cpg) %>% 
      ggplot(aes(x = test3[,1], y = test3[,2])) + 
      geom_point() +
      geom_smooth(method = "lm") +
      stat_cor(label.y = 0.96)+
      labs(x = chem,
           y = cpg)
  })

##### Legacy code ----------

# GSEA analysis -----------------------------------------------------------



# gesa_chem283_c <- champ.GSEA(beta=combined_resid_filter_pd_c,
#            DMP=dmp_evernever_c[[2]],
#            DMR=dmr_evernever_c[[2]],
#            CpGlist=NULL,
#            Genelist=NULL,
#            pheno=c_lb_sd_case_wt_10_evernever$chem283,
#            method="fisher",
#            arraytype="450K",
#            Rplot=TRUE,
#            adjPval=0.05,
#            cores=1)


# Annotation and GSEA -----------------------------------------------------

list(dmp_c_new, 
     dmp_r_new) %>% 
  map(function(champ_list){
    champ_list %>% 
      map(function(data){
        data %>% 
          # slice_min(p.value, n=2500) %>% 
          mutate(cpgs = row.names(.))
      })
  }) %>% 
  set_names("pest_champ_c","pest_champ_r") %>% 
  list2env(.,envir = .GlobalEnv)

total_champ <- dmp_count %>% 
  map(function(data){
    data %>% 
      mutate(cpgs = row.names(.)) %>% 
      filter(P.Value < 10e-7)
  })

ann450ksubt_count <- total_champ %>% 
  map(function(df){
  as.data.frame(ann450k[match(df$cpgs,ann450k$Name), 
                        c(1:4,12:19,24:ncol(ann450k))]) %>% 
    select(-c(DHS, Phantom, Enhancer))
    })

ewas_annot_count <- list(total_champ, ann450ksubt_count) %>% 
  pmap(function(data1,data2){
    cbind(data1,data2)
  })

meta.cpg_count <- ewas_annot_count %>% 
  map(function(df){
    df %>% 
      pull(P.Value) %>% 
      set_names(df$cpgs)
  })


library(methylGSA)
gesa_count <- methylRRA(cpg.pval = meta.cpg_count[[1]], method = "ORA")


library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)

ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

list(pest_champ_c,pest_champ_r) %>% 
  map(function(data){
    data %>% 
      map(function(df){
        as.data.frame(ann450k[match(df$cpgs,ann450k$Name), 
                              c(1:4,12:19,24:ncol(ann450k))]) %>% 
          select(-c(DHS, Phantom, Enhancer))
      })
  }) %>% 
  set_names("ann450ksubt_c","ann450ksubt_r") %>% 
  list2env(.,envir = .GlobalEnv)

list(list(pest_champ_c,ann450ksubt_c),
     list(pest_champ_r,ann450ksubt_r)) %>% 
  map(function(data){
    data %>% 
      pmap(function(data1,data2){
        cbind(data1,data2)
      })
  }) %>% 
  set_names("ewas_annot_c","ewas_annot_r") %>% 
  list2env(.,envir = .GlobalEnv)


list(list(ewas_annot_c,names(ewas_annot_c)),
     list(ewas_annot_r,names(ewas_annot_r))) %>% 
  map(function(plist){
    plist %>% 
      pmap(function(data, i){
        data %>% 
          mutate(chemcode = i)
      }) %>% 
      bind_rows() %>% 
      relocate(chemcode,cpgs)
  }) %>% 
  set_names("ewas_annot_new_c","ewas_annot_new_r") %>% 
  list2env(.,envir = .GlobalEnv)

library(methylGSA)

#combined
list(ewas_annot_new_c,ewas_annot_new_r) %>% 
  map(function(data){
    data %>% 
      pull(P.Value) %>% 
      set_names(data$cpgs)
  }) %>% 
  set_names("meta.cpg_c","meta.cpg_r") %>% 
  list2env(.,envir = .GlobalEnv)

#iteration
list(ewas_annot_c, ewas_annot_r) %>% 
  map(function(data){
    data %>% 
      map(function(df){
        df %>% 
          pull(P.Value) %>% 
          set_names(df$cpgs)
      })
  }) %>% 
  set_names("meta.cpg_list_c","meta.cpg_list_r") %>% 
  list2env(.,envir = .GlobalEnv)

#function 1: methylglm
list(meta.cpg_c, meta.cpg_r) %>% 
  map(function(data){
    c("KEGG","GO","Reactome") %>% 
      map(function(method){
        methylglm(cpg.pval = data, GS.type = method)
      }) 
  }) %>% 
  set_names("methylglm_c","methylglm_r") %>% 
  list2env(.,envir = .GlobalEnv)

#function 2: methylRRA
list(meta.cpg_c, meta.cpg_r) %>% 
  map(function(data){
    methylRRA(cpg.pval = data, method = "GSEA")
  }) %>% 
  set_names("methylRRA_c","methylRRA_r") %>% 
  list2env(.,envir = .GlobalEnv)

barplot(methylRRA_c, num = 10, colorby = "pvalue")
barplot(methylRRA_r, num = 10, colorby = "pvalue")

list(meta.cpg_list_c, meta.cpg_list_r) %>% 
  map(function(meta_list){
    meta_list %>% 
      map(function(data){
        methylRRA(cpg.pval = data, method = "ORA")
      })
  }) %>% 
  set_names("methylRRA_list_c","methylRRA_list_r") %>% 
  list2env(.,envir = .GlobalEnv)


methylRRA_list_r <- meta.cpg_list_r %>% 
  map(function(data){
    methylRRA(cpg.pval = data, method = "ORA")
  })

methylRRA_chem354_c <- methylRRA(cpg.pval = meta.cpg_list_c[["chem354"]], method = "GSEA")
methylRRA_chem251_c <- methylRRA(cpg.pval = meta.cpg_list_c[["chem251"]], method = "GSEA")
methylRRA_chem599_c <- methylRRA(cpg.pval = meta.cpg_list_c[["chem599"]], method = "GSEA")

barplot(methylRRA_chem354_c, num = 10, colorby = "pvalue")
barplot(methylRRA_chem251_c, num = 10, colorby = "pvalue")
barplot(methylRRA_chem599_c, num = 10, colorby = "pvalue")

methylRRA_chem354_r <- methylRRA(cpg.pval = meta.cpg_list_r[["chem354"]], method = "GSEA")
methylRRA_chem251_r <- methylRRA(cpg.pval = meta.cpg_list_r[["chem251"]], method = "GSEA")
methylRRA_chem599_r <- methylRRA(cpg.pval = meta.cpg_list_r[["chem599"]], method = "GSEA")

barplot(methylRRA_list_r[[2]], num = 10, colorby = "pvalue")
barplot(methylRRA_chem251_r, num = 10, colorby = "pvalue")
barplot(methylRRA_chem599_r, num = 10, colorby = "pvalue")


methylRRA_chem1638_r <- methylRRA(cpg.pval = meta.cpg_list_r[["chem1638"]], method = "GSEA")

# Pathway analysis --------------------------------------------------------

library(PANTHER.db)
pthOrganisms(PANTHER.db) <- "HUMAN"
PANTHER.db


go_ids <- keys(PANTHER.db,keytype="PATHWAY_ID")
cols <- "ENTREZ"
panther <- mapIds(PANTHER.db, keys=go_ids, column=cols, 
                  as.numeric(cols), keytype="PATHWAY_ID", multiVals="list")
lengths(panther)

#ALL CpG sites
#GSEA for panther pathways
list(meta.cpg_c, meta.cpg_r) %>% 
  map(function(data){
    methylglm(cpg.pval = data, GS.list = panther, GS.idtype = "ENTREZID")
  }) %>% 
  set_names("methylglm_pan_c","methylglm_pan_r") %>% 
  list2env(.,envir = .GlobalEnv)

methylglm_pan_count <- methylglm(cpg.pval = meta.cpg_count[[1]], 
                                 GS.list = panther, GS.idtype = "ENTREZID")
methylglm_path_count <- methylglm_pan_count %>% 
  left_join(panther.names, by = "ID")

#get pathway names
cols <- "PATHWAY_TERM"
res.p <- mapIds(PANTHER.db, keys=go_ids, column=cols, 
                as.numeric(cols), keytype="PATHWAY_ID", multiVals="list")

panther.names <- res.p %>% 
  bind_rows() %>% 
  gather() %>% 
  transmute(
    ID = key,
    pathway = value
  )

list(methylglm_pan_c, methylglm_pan_r) %>% 
  map(function(data){
    data %>% 
      left_join(panther.names, by = "ID")
  }) %>% 
  set_names("methylglm_path_c","methylglm_path_r") %>% 
  list2env(.,envir = .GlobalEnv)


# Visualization -----------------------------------------------------------


#occupational
dmp_count_rownames <- names(dmp_count) %>% 
  map(function(data){
    rownames(dmp_count[[data]])[1]
  }) %>% 
  as_vector()

test <- combined_resid_filter_pd_r["cg02052410",]
test2 <- t(test) %>% as_tibble()
hist(t(test))

plotlist  <-list(names(dmp_count),
     dmp_count_rownames)%>% 
  pmap(function(chem,cpg){
    test <- dmp_count[[chem]]
    test2 <- t(combined_resid_filter_pd_r %>% 
                 filter(rownames(.) %in% rownames(test))) %>% 
      as.data.frame()
    test3 <- cbind(count_combine %>% 
                     dplyr::select(all_of(chem)), test2)
    test3 %>% 
      dplyr::select(chem,cpg) %>% 
      ggplot(aes(x = test3[,1], y = test3[,2])) + 
      geom_point() +
      geom_smooth(method = "lm") +
      stat_cor(label.y = 0.99)+
      labs(x = chem,
           y = cpg)+
      theme_classic()
  })
plotlist

png(file=here("figures","scatter plots for chemicals and the corresponding top hit CpGs in residential setting_new.png"), 
    width = 1920, height = 1920)
ggarrange(plotlist = plotlist,
          labels = names(dmp_r_new),
          ncol = 2, nrow = 1)
dev.off()

#occupational
list(quote_all(chem1876,chem1673,chem1638),
     quote_all(cg12017057, cg15569292, cg10397063))%>% 
  pmap(function(chem,cpg){
    test <- champ_dmplist_case_heavymetal_sig_c[[chem]]
    test2 <- t(combined_resid_filter_pd_c %>% 
                 filter(rownames(.) %in% rownames(test))) %>% 
      as.data.frame()
    test3 <- cbind(c_lb_sd_case_wt_10_new %>% 
                     dplyr::select(chem), test2)
    test3 %>% 
      dplyr::select(chem,cpg) %>% 
      ggplot(aes(x = test3[,1], y = test3[,2])) + 
      geom_point() +
      geom_smooth(method = "loess") +
      labs(x = chem,
           y = cpg)
  })

#residential
list(quote_all(chem151,chem155,chem156),
     quote_all(cg22190861, cg09741070, cg26569590))%>% 
  pmap(function(chem,cpg){
    test <- champ_dmplist_case_heavymetal_sig_r[[chem]]
    test2 <- t(combined_resid_filter_pd_r %>% 
                 filter(rownames(.) %in% rownames(test))) %>% 
      as.data.frame()
    test3 <- cbind(r_lb_sd_case_wt_10_new %>% 
                     dplyr::select(all_of(chem)), test2)
    test3 %>% 
      dplyr::select(chem,cpg) %>% 
      ggplot(aes(x = test3[,1], y = test3[,2])) + 
      geom_point() +
      geom_smooth(method = "loess") +
      labs(x = chem,
           y = cpg)
  })

#residential
list(quote_all(chem1876,chem1673,chem1638),
     quote_all(cg04334751, cg15569292, cg21103170))%>% 
  pmap(function(chem,cpg){
    test <- champ_dmplist_case_heavymetal_sig_r[[chem]]
    test2 <- t(combined_resid_filter_pd_r %>% 
                 filter(rownames(.) %in% rownames(test))) %>% 
      as.data.frame()
    test3 <- cbind(r_lb_sd_case_wt_10_new %>% 
                     dplyr::select(chem), test2)
    test3 %>% 
      dplyr::select(chem,cpg) %>% 
      ggplot(aes(x = test3[,1], y = test3[,2])) + 
      geom_point() +
      geom_smooth(method = "loess") +
      labs(x = chem,
           y = cpg)
  })



# manhattan plot ----------------------------------------------------------



datanew <- dmp_count$total %>% 
  as_tibble() %>% 
  mutate(chr = paste("chr", CHR, sep=""))

chromosomes <- paste("chr", c(1:22, "X","Y"), sep="")
chromosomes <- intersect(chromosomes, datanew$chr)
chromosome.lengths <- sapply(chromosomes, function(chromosome)
  max(datanew$MAPINFO[which(datanew$chr == chromosome)]))


chromosome.lengths <- as.numeric(chromosome.lengths)
chromosome.starts <- c(1,cumsum(chromosome.lengths)+1)
names(chromosome.starts) <- c(chromosomes, "NA")
datanew$global <- datanew$MAPINFO + 
  chromosome.starts[datanew$chr] - 1


axis_set <- datanew %>% 
  group_by(chr) %>% 
  summarize(center = mean(global))


ylim <- datanew %>% 
  filter(P.Value == min(P.Value)) %>% 
  mutate(ylim = abs(floor(log10(P.Value))) + 2) %>% 
  pull(ylim)

library(ggtext)

datanew %>% 
  ggplot(aes(x = global, y = -log10(P.Value), 
             color = as_factor(chr), size = -log10(P.Value))) +
  geom_hline(yintercept = -log10(10e-7), color = "grey40", linetype = "dashed") + 
  geom_point(alpha = 0.75) +
  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(2.5, ylim)) +
  #scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis$chr)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, 
       y = "-log<sub>10</sub>(p)") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5))
