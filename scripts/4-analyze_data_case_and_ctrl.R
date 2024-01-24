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
myvar <- quote_all(pegid, sampleid, age, pd_new,female,
                   smokers, a1_schyrs, ethnicity, county,
                   pdstudystudy, meanmethbysample)


#using count data
list(lb_sd_metal_wt_10_count, lb_sd_copper_wt_10_count, 
     lb_sd_op_wt_10_count) %>% 
  map(function(df){
    list(
      df[[1]], df[[2]]
    ) %>% 
      map(function(datalist){
        datalist %>% 
          map(function(data){
            data %>% 
              select(pegid, count)
          }) %>% 
          reduce(full_join, by = "pegid") %>% 
          mutate_all(~replace(., is.na(.), 0)) %>% 
          transmute(pegid = pegid,
                    total = count.x + count.y)
      })
  }) %>% 
  set_names("count_combine_metal", 
            "count_combine_copper",
            "count_combine_op") %>% 
  list2env(.GlobalEnv)

list(count_combine_metal, count_combine_copper, count_combine_op) %>% 
  map(function(data){
    data %>% 
      bind_rows()
  }) %>% 
  set_names("count_combine_metal_total", "count_combine_copper_total", 
            "count_combine_op_total") %>% 
  list2env(.GlobalEnv)


hist(count_combine_metal_total$total)
skim(count_combine_metal_total$total)

hist(count_combine_copper_total$total)
skim(count_combine_copper_total$total)

hist(count_combine_op_total$total)
skim(count_combine_op_total$total)


c(lb_sd_metal_wt_10_count[[1]], 
  lb_sd_metal_wt_10_count[[2]]) %>% 
  map(function(data){
    data %>% 
      dplyr::select(all_of(myvar))
  }) %>% 
  rbindlist() %>% 
  left_join(count_combine_copper_total, by = "pegid") %>% 
  left_join(count_combine_op_total %>% 
              select(pegid, total), by = "pegid") %>%
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
    total.x = "Copper count",
    total.y = "OP count",
    county = "County"
  ) %>% 
  distinct() %>% 
  #filter(pegid %in% datSampleSteve$externaldnacode) %>% 
  dplyr::select(-c(sampleid,pegid, pdstudystudy, county, a1_schyrs)) %>%
  tbl_summary(missing = "no",
              by = pd_new,
              #type=list(c(female) ~ "categorical"),
              statistic = list(all_continuous() ~ "{mean} ({sd})",
                               all_categorical() ~ "{n} ({p}%)"),
              digits = list(all_continuous() ~ 1,
                            all_categorical() ~ 1)) %>%
  add_overall() %>% 
  modify_header(label = "**Characteristics**") %>%
  # modify_spanning_header(starts_with("stat_") ~ "**PD status**") %>%
  modify_caption("Table 1. Demographic characteristics of PEG participants (N = {N})") %>%
  bold_labels() %>% 
  table1()


#heavy metal chemical use in total
list(metal_filter, copper_filter, op_filter) %>% 
  map(function(chem){
    exp_window_address_all_list %>% 
      map(function(data){
        data %>% 
          filter(chemcode %in% chem) %>%
          mutate_at(vars(sum_total_lbs), extreme_remove_percentile_win) %>% 
          group_by(pegid,year) %>% 
          summarise(sum_total_lbs = sum(sum_total_lbs),.groups = "keep") %>% 
          inner_join(pest_methylation_clean[[1]], by = "pegid") %>%
          replace_na(list(sum_total_lbs = 0))
      }) 
  }) %>% 
  set_names("outcheck_metal_all", "outcheck_copper_all", "outcheck_op_all") %>% 
  list2env(.,envir = .GlobalEnv)

chem_usage_plot <- list(
  list(outcheck_metal_all, outcheck_copper_all, outcheck_op_all),
  list("Average heavy metal usage over time", "Average copper usage over time",
       "Average OP usage over time")
) %>% 
  pmap(function(data, x){
    list(data,
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
           color = "Location",
           title = x) +
      #geom_point()+
      geom_line(linewidth=1) +
      scale_color_colorblind() +
      # scale_color_ordinal(labels=c("Without PD", "With PD")) +
      geom_vline(xintercept = 1989, lty=2) +
      theme_classic()
  })


ggarrange(plotlist = list(chem_usage_plot[[2]], 
                          chem_usage_plot[[3]]),
          labels = c("A", "B"),
          ncol = 2, nrow = 1,
          common.legend = T,
          legend = "bottom")
#heavy metal & op chemical use: specific

list(metal_filter, copper_filter, op_filter) %>% 
  map(function(chem){
    exp_window_address_all_list %>% 
      map(function(data){
        data %>% 
          filter(chemcode %in% chem) %>%
          # mutate_at(vars(sum_total_lbs), extreme_remove_percentile_win) %>% 
          group_by(chemcode) %>% 
          group_split() %>% 
          map(function(df){
            df %>% 
              mutate_at(vars(sum_total_lbs), extreme_remove_percentile_win) %>% 
              group_by(pegid,year) %>% 
              summarise(sum_total_lbs = sum(sum_total_lbs),.groups = "keep") %>% 
              inner_join(pest_methylation_clean[[1]], by = "pegid") %>%
              replace_na(list(sum_total_lbs = 0))
          })
      }) 
  }) %>% 
  set_names("outcheck_metal_specific", "outcheck_copper_specific",
            "outcheck_op_specific") %>% 
  list2env(.,envir = .GlobalEnv)

list(outcheck_metal_specific, outcheck_copper_specific, 
     outcheck_op_specific) %>% 
  map(function(data){
    data %>% 
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
          theme_classic()
      })
  }) %>% 
  set_names("plotlist_metal", "plotlist_copper", "plotlist_op") %>% 
  list2env(.GlobalEnv)


list(
  list(plotlist_metal, plotlist_copper, plotlist_op),
  # list("metal_usage.png", "copper_usage.png", "op_usage.png"),
  list(metal_filter_name$chemname, copper_filter_name$chemname, 
       op_filter_name$chemname)
) %>% 
  pmap(function(df1, df2, chem){
    png(file=here::here("figures", df2), 
        width = 1920, height = 1080)
    ggarrange(plotlist = df1,
              labels = chem,
              common.legend = T,
              legend = "bottom")
    dev.off()
  })

png(file=here::here("figures", "metal_usage.png"), 
    width = 1920, height = 1080, res = 125)
ggarrange(plotlist = plotlist_metal,
          labels = metal_filter_name$chemname,
          common.legend = T,
          legend = "bottom")
dev.off()

png(file=here::here("figures", "copper_usage.png"), 
    width = 1920, height = 1080, res = 125)
ggarrange(plotlist = plotlist_copper,
          labels = copper_filter_name$chemname,
          common.legend = T,
          legend = "bottom")
dev.off()

png(file=here::here("figures", "op_usage.png"), 
    width = 1920, height = 1080, res = 125)
ggarrange(plotlist = plotlist_op,
          labels = op_filter_name$chemname,
          common.legend = T,
          legend = "bottom")
dev.off()

# DMP --------------------------------------------------------------------

#Differentially methylated positions (DMP)


#Method 1: ChAMP
library(ChAMP)

list(
  count_combine_metal,
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


DMP.GUI(DMP = dmp_count_case[[1]], 
        beta = as.matrix(combined_resid_filter_pd_r), 
        pheno = count_combine_case$total,
        cutgroupnumber = 2)

DMP.GUI(DMP = dmp_count_ctrl[[1]], 
        beta = as.matrix(combined_resid_filter_ctrl_r), 
        pheno = count_combine_ctrl$total,
        cutgroupnumber = 2)



#Method 2: meffil
# Load meffil and set how many cores to use for parallelization
source(here("scripts", "meffil_fixed.R"))
source(here("scripts", "meffil_report_fixed.R"))
library(meffil)
myvars_ewas <- quote_all(cd8t, cd4t, nk, mono, bcell, gran, 
                         age, female, smokers, rfvotecaucasian2, study)

# myvars_ewas <- quote_all(cd4t, gran, age, female,
#                          pdstudynumberyearswithpdatbloodd,
#                          rfvotecaucasian2, pd_new)


#select covariates
list(lb_sd_copper_wt_10_count, lb_sd_op_wt_10_count) %>%
  pmap(function(df1, df2){
    list(df1, df2) %>%
      pmap(function(data1, data2){
        data1 %>%
          dplyr::select(pegid, all_of(myvars_ewas)) %>%
          left_join(data2 %>%
                      dplyr::select(pegid, count), by = "pegid") %>%
          # dplyr::select(-pegid) %>%
          # rename(pdyears = pdstudynumberyearswithpdatbloodd) %>%
          # mutate(pdyears = ifelse(is.na(pdyears),
          #                         mean(pdyears, na.rm=TRUE), pdyears)) %>%
          # replace_na(list(pdyears = 0)) %>%
          na.omit()
        #replace_na(list(c11_depression = 0,gds1_5 = 0)) %>%
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


covar_total <- list(
  list(covar_ctrl_combind, covar_case_combind),
  c("ctrl", "case")
) %>% 
  pmap(function(data1, data2){
    data1 %>% 
      mutate(pd = data2) 
  }) %>% 
  bind_rows() %>% 
  mutate(study = if_else(is.na(study), "PEG 1", study))

ewas.parameters <- meffil.ewas.parameters(
  sig.threshold=1e-6,  ## EWAS p-value threshold
  max.plots=20, ## plot at most 20 CpG sites
  qq.inflation.method="median",  ## measure inflation using median
  model="all") ## select default EWAS model; 

list(list(count_combine_copper[[1]], count_combine_copper[[2]], 
          count_combine_copper_total, 
          count_combine_copper[[1]], count_combine_copper[[2]], 
          count_combine_copper_total),
     list(PEG_NOOB_nors_win_filter_pd_r, PEG_NOOB_nors_win_filter_ctrl_r, 
          peg_noob_nors_win_total,
          PEG_NOOB_nors_win_filter_pd_r, PEG_NOOB_nors_win_filter_ctrl_r, 
          peg_noob_nors_win_total),
     list(covar_case_combind %>% dplyr::select(-count), 
          covar_ctrl_combind %>% dplyr::select(-count), 
          covar_total %>% dplyr::select(-count),
          covar_case_combind, covar_ctrl_combind, covar_total)) %>% 
  pmap(function(data1,data2,data3){
    data1 %>% 
      dplyr::select(total) %>% 
      map(~meffil.ewas( beta=as.matrix(data2), 
                        variable=.x, 
                        covariates = data3, 
                        batch = NULL, weights = NULL,  cell.counts = NULL,
                        isva = F, sva = F, smartsva = F, n.sv = NULL, 
                        winsorize.pct = NA, robust = TRUE,
                        rlm = FALSE, outlier.iqr.factor = NA,  featureset = NA, 
                        random.seed = 20230922, lmfit.safer = F, verbose = T))
  }) %>% 
  set_names("meffil_count_noop_case", "meffil_count_noop_ctrl", 
            "meffil_count_noop_total",
            "meffil_count_op_case", "meffil_count_op_ctrl", 
            "meffil_count_op_total") %>% 
  list2env(.,envir = .GlobalEnv)


test <- meffil_count_noop_case$total$analyses$all$table %>% 
  filter(-log10(p.value) > 6)


save(meffil_count_op_case, file = "meffil_count_op_case.RData")
save(meffil_count_op_ctrl, file = "meffil_count_op_ctrl.RData")
save(meffil_count_op_total, file = "meffil_count_op_total.RData")

save(meffil_count_noop_case, file = "meffil_count_noop_case.RData")
save(meffil_count_noop_ctrl, file = "meffil_count_noop_ctrl.RData")
save(meffil_count_noop_total, file = "meffil_count_noop_total.RData")


list(list(meffil_count_op_case, 
          meffil_count_op_ctrl, meffil_count_op_total,
          meffil_count_noop_case, 
          meffil_count_noop_ctrl, meffil_count_noop_total),
     list(PEG_NOOB_nors_win_filter_pd_r, 
          PEG_NOOB_nors_win_filter_ctrl_r, 
          peg_noob_nors_win_total,
          PEG_NOOB_nors_win_filter_pd_r, 
          PEG_NOOB_nors_win_filter_ctrl_r, 
          peg_noob_nors_win_total)) %>% 
  pmap(function(meffil_list,data){
    meffil_list %>% 
      map(~ meffil.ewas.summary_fix(.x,data,
                                    parameters = ewas.parameters))
  }) %>% 
  set_names("ewas.summary_count_op_case", "ewas.summary_count_op_ctrl", 
            "ewas.summary_count_op_total",
            "ewas.summary_count_noop_case", "ewas.summary_count_noop_ctrl", 
            "ewas.summary_count_noop_total") %>% 
  list2env(.,envir = .GlobalEnv)

list(list(ewas.summary_count_op_case,
          names(ewas.summary_count_op_case), "case_op"),
     list(ewas.summary_count_op_ctrl,
          names(ewas.summary_count_op_ctrl), "ctrl_op"),
     list(ewas.summary_count_op_total,
          names(ewas.summary_count_op_total), "total_op"),
     list(ewas.summary_count_noop_case,
          names(ewas.summary_count_noop_case), "case_noop"),
     list(ewas.summary_count_noop_ctrl,
          names(ewas.summary_count_noop_ctrl), "ctrl_noop"),
     list(ewas.summary_count_noop_total,
          names(ewas.summary_count_noop_total), "total_noop")) %>%
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



# GESA analysis -----------------------------------------------------------

list(meffil_count_op_case, meffil_count_op_ctrl, meffil_count_op_total,
     meffil_count_noop_case, meffil_count_noop_ctrl, meffil_count_noop_total) %>% 
  map(function(meffil_list){
    meffil_list %>% 
      map(function(data){
        data$analyses$all$table %>% 
          # slice_min(p.value, n=2500) %>% 
          mutate(cpgs = row.names(.))
      })
  }) %>% 
  set_names("metal_meffil_op_case", "metal_meffil_op_ctrl", 
            "metal_meffil_op_total", "metal_meffil_noop_case", 
            "metal_meffil_noop_ctrl", "metal_meffil_noop_total") %>% 
  list2env(.,envir = .GlobalEnv)

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)

ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

list(metal_meffil_op_case, metal_meffil_op_ctrl, metal_meffil_op_total,
     metal_meffil_noop_case, metal_meffil_noop_ctrl, metal_meffil_noop_total) %>% 
  map(function(data){
    data %>% 
      map(function(df){
        as.data.frame(ann450k[match(df$cpgs,ann450k$Name), 
                              c(1:4,12:19,24:ncol(ann450k))])
      })
  }) %>% 
  set_names("ann450ksubt_op_case","ann450ksubt_op_ctrl", "ann450ksubt_op_total",
            "ann450ksubt_noop_case","ann450ksubt_noop_ctrl", 
            "ann450ksubt_noop_total") %>% 
  list2env(.,envir = .GlobalEnv)

list(
  list(meffil_count_op_case, meffil_count_op_ctrl, meffil_count_op_total,
       meffil_count_noop_case, meffil_count_noop_ctrl, meffil_count_noop_total),
  list(ann450ksubt_op_case, ann450ksubt_op_ctrl, ann450ksubt_op_total,
       ann450ksubt_noop_case, ann450ksubt_noop_ctrl, ann450ksubt_noop_total)
) %>% 
  pmap(function(data1, data2){
    data1$total$analyses$all$table %>% 
      filter(-log10(p.value) > 6) %>% 
      rownames_to_column("cpg") %>% 
      mutate(chr = chromosome, 
             pos = position) %>% 
      dplyr::select(-c(chromosome, position)) %>% 
      left_join(data2$total %>% 
                  dplyr::select(chr, pos, Probe_rs, Probe_maf, UCSC_RefGene_Name, 
                                UCSC_RefGene_Group), by = c("chr", "pos")) %>% 
      dplyr::rename(gene = UCSC_RefGene_Name,
                    gene_region = UCSC_RefGene_Group) %>% 
      relocate(cpg, chr, pos, Probe_rs, Probe_maf, gene, 
               gene_region)
  }) %>% 
  set_names("annote_table_op_case", "annote_table_op_ctrl", 
            "annote_table_op_total",
            "annote_table_noop_case", "annote_table_noop_ctrl", 
            "annote_table_noop_total") %>% 
  list2env(.GlobalEnv)


# list(meffil_count_op_case, meffil_count_op_ctrl, meffil_count_op_total,
#      meffil_count_noop_case, meffil_count_noop_ctrl, meffil_count_noop_total) %>% 
#   map(function(data){
#     data$copper$analyses$all$table %>% 
#       filter(-log10(p.value) > 6) %>% 
#       rownames_to_column("cpg") %>% 
#       mutate(chr = chromosome, 
#              pos = position) %>% 
#       dplyr::select(-c(chromosome, position)) %>% 
#       left_join(ann450ksubt_case$total %>% 
#                   dplyr::select(chr, pos, Probe_rs, Probe_maf, UCSC_RefGene_Name, 
#                                 UCSC_RefGene_Group), by = c("chr", "pos")) %>% 
#       dplyr::rename(gene = UCSC_RefGene_Name,
#                     gene_region = UCSC_RefGene_Group) %>% 
#       relocate(cpg, chr, pos, Probe_rs, Probe_maf, gene, 
#                gene_region)
#   }) %>% 
#   set_names("annote_table_copper_case", "annote_table_copper_ctrl", 
#             "annote_table_copper_total", "annote_table_copper_noop_case", 
#             "annote_table_copper_noop_ctrl", 
#             "annote_table_copper_noop_total") %>% 
#   list2env(.GlobalEnv)

table_names <- list(annote_table_op_case = annote_table_op_case,
                    annote_table_noop_case = annote_table_noop_case)

list(table_names,
     names(table_names)) %>% 
  pmap(function(data1,data2){
    write_csv(data1, 
              here("tables", paste(data2, ".csv", sep = "")))
  })



list(list(metal_meffil_op_case, ann450ksubt_op_case),
     list(metal_meffil_op_ctrl, ann450ksubt_op_ctrl),
     list(metal_meffil_op_total, ann450ksubt_op_total),
     list(metal_meffil_noop_case, ann450ksubt_noop_case),
     list(metal_meffil_noop_ctrl, ann450ksubt_noop_ctrl),
     list(metal_meffil_noop_total, ann450ksubt_noop_total)) %>% 
  map(function(data){
    data %>% 
      pmap(function(data1,data2){
        cbind(data1,data2)
      })
  }) %>% 
  set_names("ewas_annot_op_case", "ewas_annot_op_ctrl", "ewas_annot_op_total",
            "ewas_annot_noop_case", "ewas_annot_noop_ctrl", 
            "ewas_annot_noop_total") %>% 
  list2env(.,envir = .GlobalEnv)



list(ewas_annot_op_case, ewas_annot_op_ctrl, ewas_annot_op_total,
     ewas_annot_noop_case, ewas_annot_noop_ctrl, ewas_annot_noop_total) %>% 
  map(function(plist){
    plist %>% 
      map(~filter(.x, -log10(p.value) > 6))
  }) %>% 
  set_names("ewas_annot_new_op_case","ewas_annot_new_op_ctrl", 
            "ewas_annot_new_op_total", "ewas_annot_new_noop_case",
            "ewas_annot_new_noop_ctrl", "ewas_annot_new_noop_total") %>% 
  list2env(.,envir = .GlobalEnv)

library(methylGSA)

list(ewas_annot_op_case, ewas_annot_noop_case) %>% 
  map(function(datalist){
    datalist %>% 
      map(function(data){
        data %>% 
          pull(p.value) %>% 
          set_names(data$cpgs)
      })
  }) %>% 
  set_names("meta.cpg_op_case", "meta.cpg_noop_case") %>% 
  list2env(.GlobalEnv)

?methylRRA()
#function 2: methylRRA
# gsea_case_total <- methylRRA(cpg.pval = meta.cpg_case_total, method = "GSEA")
gsea_case_copper <- methylRRA(cpg.pval = meta.cpg_op_case[[1]], 
                              array.type = "450K", group = "all", sig.cut = 0.05, 
                              method = "GSEA", GS.idtype = "SYMBOL", 
                              GS.type = "GO", minsize = 100, maxsize = 500)


gsea_case_noop_copper <- methylRRA(cpg.pval = meta.cpg_noop_case[[1]], method = "ORA")

?methylRRA()

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
methylglm_pan_case_op <- methylglm(cpg.pval = meta.cpg_op_case[[1]], 
                                GS.list = panther, GS.idtype = "ENTREZID")

methylglm_pan_case_noop <- methylglm(cpg.pval = meta.cpg_noop_case[[1]], 
                                GS.list = panther, GS.idtype = "ENTREZID")
?methylglm()
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

methylglm_path_case_total <- methylglm_pan_case_total %>% 
  left_join(panther.names, by = "ID")

methylglm_path_case_copper <- methylglm_pan_case_copper %>% 
  left_join(panther.names, by = "ID")

# Visulization ------------------------------------------------------------

library(ggtext)

list(meffil_count_case, meffil_count_ctrl, meffil_count_total,
     meffil_count_noop_case, meffil_count_noop_ctrl, meffil_count_noop_total) %>% 
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
  set_names("metal_list_case","metal_list_ctrl", "metal_list_total",
            "metal_list_noop_case","metal_list_noop_ctrl", "metal_list_noop_total") %>% 
  list2env(.,envir = .GlobalEnv)

list(metal_list_case, metal_list_ctrl, metal_list_total,
     metal_list_noop_case, metal_list_noop_ctrl, metal_list_noop_total) %>% 
  map(function(metal_list){
    metal_list %>% 
      map(function(data){
        data %>% 
          group_by(chr) %>% 
          dplyr::summarize(center = mean(global)) %>% 
          ungroup()
      })
  }) %>% 
  set_names("axis_set_case", "axis_set_ctrl", "axis_set_total",
            "axis_set_noop_case", "axis_set_noop_ctrl", "axis_set_noop_total") %>% 
  list2env(.,envir = .GlobalEnv)


list(metal_list_case, metal_list_ctrl, metal_list_total,
     metal_list_noop_case, metal_list_noop_ctrl, metal_list_noop_total) %>% 
  map(function(metal_list){
    metal_list %>% 
      map(function(data){
        data %>% 
          filter(p.value == min(p.value)) %>% 
          mutate(ylim = abs(floor(log10(p.value))) + 1) %>% 
          pull(ylim)
      })
  }) %>% 
  set_names("ylim_case","ylim_ctrl", "ylim_total",
            "ylim_noop_case","ylim_noop_ctrl", "ylim_noop_total") %>% 
  list2env(.,envir = .GlobalEnv)


plotlist <- list(
  list(metal_list_case, axis_set_case, ylim_case),
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
          geom_hline(yintercept = -log10(10e-7), color = "red", linetype = "dashed") + 
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

list(
  list(metal_list_noop_case, axis_set_noop_case, ylim_noop_case),
  list(metal_list_noop_ctrl, axis_set_noop_ctrl, ylim_noop_ctrl),
  list(metal_list_noop_total, axis_set_noop_total, ylim_noop_total)) %>% 
  map(function(datalist){
    datalist %>% 
      pmap(function(pest, axis, ylim){
        pest %>% 
          mutate(sig = if_else(-log10(p.value) > 6, 1, 0)) %>% 
          ggplot(aes(x = global, y = -log10(p.value), 
                     color = as_factor(sig), size = -log10(p.value))) +
          scale_color_manual(values = c("1" = "red", "0" = "black")) +
          geom_hline(yintercept = -log10(10e-7), color = "red", linetype = "dashed") + 
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
            axis.text = element_text(size = 10,face="bold"),
            axis.title = element_text(size=20,face="bold"),
            axis.title.y = element_markdown(),
            axis.text.x = element_text(angle = 60, size = 10, vjust = 0.5)
          )
      })
  }) %>% 
  set_names("manhattanlist_noop_case","manhattanlist_noop_ctrl",
            "manhattanlist_noop_total") %>% 
  list2env(.,envir = .GlobalEnv)



plotlist <- list(manhattanlist_noop_total[[2]],
                 manhattanlist_noop_case[[2]], 
                 manhattanlist_noop_ctrl[[2]]
                 )


#print all plots in one figure
png(file=here::here("figures","manhattan plots for cases and controls_copper_adjop_122423.png"), 
    width = 1920, height = 1080)
ggarrange(plotlist = plotlist,
          labels = c("A", "B", "C"),
          ncol = 2, nrow = 2)
dev.off()


#scatter plot

mean_methylist <- list(
  list(combined_resid_total, combined_resid_filter_pd_r, 
       combined_resid_filter_ctrl_r),
  list(count_combine_metal_total, count_combine_metal[[1]], 
       count_combine_metal[[2]])
) %>% 
  pmap(function(df1, df2){
    df1 %>% 
      map(~mean(.x)) %>% 
      bind_rows() %>% 
      pivot_longer(
        cols = starts_with("X"),
        names_to = "sampleid",
        values_to = "mean_methyl") %>% 
      left_join(datSamplePEG %>% 
                  select(sampleid, pegid), 
                by = "sampleid") %>% 
      left_join(df2 %>% 
                  select(pegid, copper), by = "pegid")
  }) 

plotlist <- mean_methylist %>% 
  map(function(data){
    data %>% 
      ggplot(aes(x = copper, y = mean_methyl)) + 
      geom_point() +
      geom_smooth(method = "lm") +
      # geom_jitter()+
      stat_cor(label.y = 0.54)+
      theme_classic()+
      labs(x = "Copper count",
           y = "Mean methylation residual")
  })

ggarrange(plotlist = plotlist,
          labels = c("A", "B", "C"),
          ncol = 2, nrow = 2)

list(combined_resid_filter_ctrl_r, combined_resid_filter_pd_r) %>% 
  map(function(data){
    data %>% 
      t() %>% 
      as.data.frame() %>% 
      rownames_to_column("sampleid") %>% 
      left_join(datSamplePEG %>% 
                  select(pegid, sampleid), by ="sampleid") %>% 
      select(-sampleid) %>% 
      column_to_rownames("pegid") %>% 
      transmute(row_mean = rowMeans(across(everything())))
  }) %>% 
  set_names("residual_ctrl_mean", "residual_case_mean") %>% 
  list2env(.GlobalEnv)

list(
  list(residual_case_mean, residual_ctrl_mean),
  count_combine_metal
) %>% 
  pmap(function(df1, df2){
    df1 %>% 
      rownames_to_column("pegid") %>% 
      left_join(df2, by = "pegid") %>% 
      ggplot(aes(x = total, y = row_mean)) + 
      geom_point() +
      geom_smooth(method = "lm") +
      stat_cor(label.y = 0.54)+
      labs(x = "Metal counts",
           y = "Mean methylation residual")
  })



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

#--------------------------------End of the code--------------------------------