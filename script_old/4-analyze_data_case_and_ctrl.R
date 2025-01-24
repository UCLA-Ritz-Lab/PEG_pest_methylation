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
list(lb_sd_metal_wt_count, lb_sd_copper_wt_count, 
     lb_sd_op_wt_count) %>% 
  map(function(df){
    list(
      df[[1]], df[[2]]
    ) %>% 
      map(function(datalist){
        list(datalist, 
             list(exp_yrs_lag_c, exp_yrs_lag_r)) %>%
          pmap(function(data, exp_yrs){
            data %>% 
              dplyr::select(pegid, count) %>% 
              left_join(exp_yrs, by = "pegid")
              # mutate(count = count/exp_yrs_lag_dr)
          }) %>% 
          purrr::reduce(full_join, by = "pegid") %>% 
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
count_combine_copper[[2]] %>% 
  ggplot(aes(x = log(total+1))) + 
  geom_histogram(binwidth=1)
skim(count_combine_copper[[2]]$total)

hist(count_combine_op_total$total)
count_combine_op_total %>% 
  ggplot(aes(x = total)) + 
  geom_histogram(binwidth=1)

skim(count_combine_op_total$total)

#using duration data
list(heavy_metal$chemcode, chem_copper$chemcode, chem_op$chemcode) %>% 
  map(function(chemlist){
    list(
      list(exp_window_address_case, exp_window_address_control),
      count_combine_copper
    ) %>%
      pmap(function(data1, data2){
        data1 %>%
          map(function(df){
            df %>% 
              group_by(pegid) %>%
              filter(year <= indexyr) %>% 
              filter(chemcode %in% chemlist) %>% 
              ungroup() %>% 
              replace_na(list(sum_total_lbs = 0))
          }) %>% 
          bind_rows() %>% 
          group_by(pegid, year, chemcode) %>% 
          summarise(sum_total_lbs = sum(sum_total_lbs),
                    .groups = "keep") %>%
          ungroup() %>% 
          mutate(ind = if_else(sum_total_lbs > 0, 1, 0)) %>%
          group_by(pegid, chemcode) %>%
          summarise(duration = sum(ind),
                    .groups = "keep") %>% 
          ungroup() %>% 
          right_join(data2, by = "pegid") %>% 
          mutate_at(vars(duration), ~replace_na(., 0)) %>% 
          group_by(pegid) %>% 
          summarise(exp_duration = sum(duration),
                    .groups = "keep") %>%
          ungroup()
      })
  }) %>%
  set_names("exp_duration_metal", "exp_duration_copper", "exp_duration_op") %>%
  list2env(.GlobalEnv)

list(exp_duration_metal, exp_duration_copper, exp_duration_op) %>% 
  map(function(data){
    data %>% 
      bind_rows()
  }) %>% 
  set_names("exp_duration_metal_total", "exp_duration_copper_total", 
            "exp_duration_op_total") %>% 
  list2env(.GlobalEnv)

exp_duration_copper[[1]] %>%
  select(pegid) %>% 
  distinct() %>% 
  nrow()
  
exp_duration_copper_total %>% 
  ggplot(aes(x = exp_duration)) + 
  geom_histogram()
skim(exp_duration_copper$exp_duration)

exp_duration_op %>% 
  ggplot(aes(x = exp_duration)) + 
  geom_histogram()
skim(exp_duration_op$exp_duration)

## create table 1
c(lb_sd_metal_wt_count[[1]], 
  lb_sd_metal_wt_count[[2]]) %>% 
  map(function(data){
    data %>% 
      dplyr::select(all_of(myvar))
  }) %>% 
  rbindlist() %>% 
  left_join(count_combine_copper_total, by = "pegid") %>% 
  left_join(count_combine_op_total %>% 
              dplyr::select(pegid, total), by = "pegid") %>%
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
                            all_categorical() ~ 0)) %>%
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
              right_join(pest_methylation_clean_total, by = "pegid") %>%
              replace_na(list(sum_total_lbs = 0)) %>% 
              ungroup()
          })
      }) 
  }) %>% 
  set_names("outcheck_metal_specific", "outcheck_copper_specific",
            "outcheck_op_specific") %>% 
  list2env(.,envir = .GlobalEnv)

outcheck_copper_specific[[1]][[2]] %>% 
  select(pegid) %>% 
  distinct() %>% 
  nrow()

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
          group_by(year, pd_new) %>% 
          summarise(y=mean(sum_total_lbs, na.rm=T),.groups="keep") %>% 
          #change to average per person, including unexposed people
          ggplot(aes(x=year, y=y, group=pd_new, color = factor(pd_new))) + 
          annotate("rect", fill = "azure2", alpha = 0.5,
                   xmin = 1989, xmax = 2010,
                   ymin = -Inf, ymax = Inf) +
          labs(x = "Years",
               y = "Chemical use (lbs)",
               color = "PD") +
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
    width = 2600, height = 1600, res = 125)
ggarrange(plotlist = plotlist_metal,
          labels = metal_filter_name$chemname,
          common.legend = T,
          legend = "bottom")
dev.off()

png(file=here::here("figures", "copper_usage.png"), 
    width = 2600, height = 1600, res = 125)
ggarrange(plotlist = plotlist_copper,
          labels = copper_filter_name$chemname,
          common.legend = T,
          legend = "bottom")
dev.off()

png(file=here::here("figures", "op_usage.png"), 
    width = 2600, height = 1600, res = 125)
ggarrange(plotlist = plotlist_op,
          labels = op_filter_name$chemname,
          common.legend = T,
          legend = "bottom")
dev.off()

# DMP --------------------------------------------------------------------

#Differentially methylated positions (DMP)


#Method 1: ChAMP
library(ChAMP)

# list(
#   count_combine_copper,
#   list(combined_resid_filter_pd_r, combined_resid_filter_ctrl_r)
# ) %>% 
#   pmap(function(data1, data2){
#     data1 %>% 
#       dplyr::select(total) %>% 
#       map(~champ.DMP(beta = data2, pheno = .x, 
#                      compare.group = NULL, 
#                      adjPVal = 0.05, adjust.method = "BH", 
#                      arraytype = "450K")[[1]])
#   }) %>% 
#   set_names("dmp_count_case", "dmp_count_ctrl") %>% 
#   list2env(.GlobalEnv)


dmp_count_case <- count_combine_copper[[1]] %>%
  dplyr::select(total) %>%
  purrr::map(~champ.DMP(beta = combined_resid_win_filter_case, 
                 pheno = .x,
                 compare.group = NULL,
                 adjPVal = 1, adjust.method = "BH",
                 arraytype = "450K")[[1]])


dmp_count_ctrl <- count_combine_copper[[2]] %>%
  dplyr::select(total) %>%
  purrr::map(~champ.DMP(beta = combined_resid_win_filter_ctrl, 
                        pheno = .x,
                        compare.group = NULL,
                        adjPVal = 1, adjust.method = "BH",
                        arraytype = "450K")[[1]])

dmp_count_total <- count_combine_copper_total %>%
  dplyr::select(total) %>%
  purrr::map(~champ.DMP(beta = combined_resid_win_filter_total, 
                        pheno = .x,
                        compare.group = NULL,
                        adjPVal = 1, adjust.method = "BH",
                        arraytype = "450K")[[1]])

dmp_count_ctrlpd_total <- count_combine_copper_total %>%
  dplyr::select(total) %>%
  purrr::map(~champ.DMP(beta = combined_resid_nostudy_ctrlpd_win_filter_total, 
                        pheno = .x,
                        compare.group = NULL,
                        adjPVal = 0.05, adjust.method = "BH",
                        arraytype = "450K")[[1]])

save(dmp_count_case, file = "dmp_count_case.RData")
save(dmp_count_ctrl, file = "dmp_count_ctrl.RData")
save(dmp_count_total, file = "dmp_count_total.RData")
save(dmp_count_ctrlpd_total, file = "dmp_count_ctrlpd_total.RData")

test_case <- dmp_count_case$total %>% 
  rownames_to_column("cpg") %>% 
  filter(cpg %in% rownames(test_total))

test_control <- dmp_count_ctrl$total %>% 
  rownames_to_column("cpg") %>% 
  filter(cpg %in% rownames(test_total))

test_total <- dmp_count_total$total %>% 
  filter(-log10(P.Value) > 6)

test_ctrlpd_total <- dmp_count_ctrlpd_total$total %>% 
  filter(-log10(P.Value) > 6)

setdiff(rownames(test_case), rownames(test_total))

# dmp_count_total <- count_combine_copper_total %>%
#   dplyr::select(total) %>%
#   map(~champ.DMP(beta = combined_resid_total, 
#                  pheno = .x,
#                  compare.group = NULL,
#                  adjPVal = 0.05, adjust.method = "BH",
#                  arraytype = "450K")[[1]])

DMP.GUI(DMP = dmp_count_case[[1]], 
        beta = as.matrix(combined_resid_filter_pd_r), 
        pheno = count_combine_copper[[1]]$total)
# 
# DMP.GUI(DMP = dmp_count_ctrl[[1]], 
#         beta = as.matrix(combined_resid_filter_ctrl_r), 
#         pheno = count_combine_ctrl$total,
#         cutgroupnumber = 2)

## pearson correlation

test_total

#method 2: limma
library(limma)
case_count_combined <- count_combine_copper[[1]] %>% 
  left_join(covar_case_process[[2]], by = "pegid") %>% 
  left_join(count_combine_op[[1]], by = "pegid") %>% 
  dplyr::select(-c(pegid, count)) %>% 
  dplyr::rename(copper_count = total.x,
                op_count = total.y)

design_total <- model.matrix(~ total, data = count_combine_copper_total)
fit_total <- lmFit(combined_resid_win_filter_total, design)
fit2_total <- eBayes(fit_total)

design_case <- model.matrix(~ total, data = count_combine_copper[[1]])
fit_case <- lmFit(combined_resid_win_filter_case, design_case)
fit2_case <- eBayes(fit_case)

design_ctrl <- model.matrix(~ total, data = count_combine_copper[[2]])
fit_ctrl <- lmFit(combined_resid_win_filter_ctrl, design_ctrl)
fit2_ctrl <- eBayes(fit_ctrl)

list(
  list(fit2_total, fit2_case, fit2_ctrl),
  list(design_total, design_case, design_ctrl),
  list(combined_resid_win_filter_total, combined_resid_win_filter_case, 
       combined_resid_win_filter_ctrl),
  list(B_total, B_case, B_ctrl),
  list("total", "case", "ctrl")
) %>% 
  pmap(function(data1, data2, data3, data4, data5){
    topTable(data1, coef=ncol(data2), sort.by="p",
             number = nrow(data3), 
             adjust.method = "BY") %>% 
      rename_all(~str_c(., "_", data5)) %>%
      rownames_to_column("cpg") %>% 
      filter(cpg %in% bicor_total$ID) %>% 
      # sort rows by the cpg order in bicor_total$ID
      left_join(data4 %>% 
                  select(total) %>% 
                  rownames_to_column("cpg"), by = "cpg") %>% 
      rename(coefficient = total)
  }) %>% 
  set_names("sig_cpg_limma_total", "sig_cpg_limma_case", 
            "sig_cpg_limma_ctrl") %>%
  list2env(.GlobalEnv)

sig_cpg_limma_final <- list(sig_cpg_limma_total, 
                            sig_cpg_limma_case, 
                            sig_cpg_limma_ctrl) %>% 
  reduce(left_join, by = "cpg")

write_csv(sig_cpg_limma_final, here("tables", "Supplement Table 1.csv"))


par(mfrow=c(2,5))
sapply(rownames(sig_cpg_limma)[1:10], function(cpg){
  plotCpg(combined_resid_win_filter_total, cpg=cpg, 
          pheno=count_combine_copper_total$total, 
          ylab = "Beta values")
})
dev.off()

dat <- data.frame(foldchange = fit2[["coefficients"]][,2], 
                  logPvalue =  -log10(fit2[["p.value"]][,2]))
dat$threshold <- as.factor(abs(dat$foldchange) < 0.0005)

#Visualization
cols <- c("TRUE" = "grey", "FALSE" = "blue")
dat %>% 
  # mutate(threshold = if_else(logFC <= 0.0005, "sig", "notsig")) %>%
  ggplot(aes(x=foldchange, y = logPvalue, color=threshold)) +
  geom_point(alpha=.6, size=1.2) +
  scale_colour_manual(values = cols) +
  geom_vline(xintercept = 0.0005, colour="#990000", linetype="dashed") + 
  geom_vline(xintercept = - 0.0005, colour="#990000", linetype="dashed") +
  theme(legend.position="none") +
  xlab("Fold Change") +
  ylab("-log10 p value") +
  theme_bw() +
  theme(legend.position = "none")


#Method 3: meffil
# Load meffil and set how many cores to use for parallelization
source(here::here("scripts", "meffil_fixed.R"))
source(here::here("scripts", "meffil_report_fixed.R"))
library(meffil)
myvars_ewas <- quote_all(cd8t, cd4t, nk, mono, bcell, gran, 
                         age, female, smokers, rfvotecaucasian2, study)

# myvars_ewas <- quote_all(cd4t, gran, age, female,
#                          pdstudynumberyearswithpdatbloodd,
#                          rfvotecaucasian2, pd_new)


#select covariates
list(lb_sd_copper_wt_count, lb_sd_op_wt_count) %>%
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
  list(covar_case_combind, covar_ctrl_combind),
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
     list(PEG_NOOB_nors_win_filter_case, PEG_NOOB_nors_win_filter_ctrl, 
          PEG_NOOB_nors_win_filter_total,
          PEG_NOOB_nors_win_filter_case, PEG_NOOB_nors_win_filter_ctrl, 
          PEG_NOOB_nors_win_filter_total),
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
                        random.seed = 20240916, lmfit.safer = F, verbose = T))
  }) %>% 
  set_names("meffil_count_noop_case", "meffil_count_noop_ctrl", 
            "meffil_count_noop_total",
            "meffil_count_op_case", "meffil_count_op_ctrl", 
            "meffil_count_op_total") %>% 
  list2env(.,envir = .GlobalEnv)


test_meffil <- meffil_count_op_case$total$analyses$all$table %>% 
  filter(-log10(p.value) > 6) %>% 
  rownames_to_column("cpg") %>%
  arrange(p.value)

test_champ <- dmp_count_case$total %>% 
  filter(-log10(P.Value) > 6) %>% 
  rownames_to_column("cpg") %>%
  arrange(P.Value)
# write_csv(test, "copper_op_count_sig.csv")

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
     list(PEG_NOOB_nors_win_filter_case, PEG_NOOB_nors_win_filter_ctrl, 
          PEG_NOOB_nors_win_filter_total,
          PEG_NOOB_nors_win_filter_case, PEG_NOOB_nors_win_filter_ctrl, 
          PEG_NOOB_nors_win_filter_total)) %>% 
  pmap(function(meffil_list, data){
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


# Association test --------------------------------------------------------

### Pearson correlation
library(WGCNA)


bicor_total <- standardScreeningNumericTrait(
  datExpr = t(combined_resid_win_filter_total %>% 
                filter(rownames(.) %in% dmp_sig_total_op$cpg)),
  yNumeric = count_combine_copper_total$total,
  corFnc = "bicor",
  alternative = "two.sided"
)

bicor_case <- standardScreeningNumericTrait(
  datExpr = t(combined_resid_win_filter_case %>% 
                filter(rownames(.) %in% dmp_sig_total_op$cpg)),
  yNumeric = count_combine_copper[[1]]$total,
  corFnc = "bicor",
  alternative = "two.sided"
)

bicor_ctrl <- standardScreeningNumericTrait(
  datExpr = t(combined_resid_win_filter_ctrl %>% 
                filter(rownames(.) %in% dmp_sig_total_op$cpg)),
  yNumeric = count_combine_copper[[2]]$total,
  corFnc = "bicor",
  alternative = "two.sided"
)


save(bicor_total, file = "bicor_total.RData")
save(bicor_case, file = "bicor_case.RData")
save(bicor_ctrl, file = "bicor_ctrl.RData")

# DMR analysis ------------------------------------------------------------

# method 1: ipdmr

library(ENmix)
test <- meffil_count_op_case$total$analyses$all$table %>% 
  # filter(-log10(p.value) > 3) %>%
  dplyr::select(chromosome, position, p.value) %>% 
  dplyr::rename(end = position, 
                chr = chromosome,
                p = p.value) %>%
  mutate(start = end - 1,
         chr = str_remove(chr, "chr")) %>%
  rownames_to_column("probe") %>% 
  relocate(chr, start, end, p, probe)
?ipdmr()
ipdmr(data=test, include.all.sig.sites=TRUE,dist.cutoff=1000,bin.size=50,
      seed=0.0001,region_plot=TRUE,mht_plot=FALSE,verbose=TRUE)
combp(data = test, bin.size=50, nCores = 1)

combp
#method 2: DMRcate
library(DMRcate)

case_count_combined <- count_combine_copper[[1]] %>% 
  left_join(covar_case_process[[2]], by = "pegid") %>% 
  left_join(count_combine_op[[1]], by = "pegid") %>% 
  dplyr::select(-c(pegid, count)) %>% 
  dplyr::rename(copper_count = total.x,
                op_count = total.y)

ctrl_count_combined <- count_combine_copper[[2]] %>% 
  left_join(covar_ctrl_process[[2]], by = "pegid") %>% 
  left_join(count_combine_op[[2]], by = "pegid") %>% 
  dplyr::select(-c(pegid, count, study)) %>% 
  dplyr::rename(copper_count = total.x,
                op_count = total.y)

total_count_combined <- list(case_count_combined, ctrl_count_combined) %>% 
  bind_rows()




design_case <- model.matrix(~ total, data = count_combine_copper[[1]])
design_ctrl <- model.matrix(~ ., data = ctrl_count_combined)


test <- count_combine_copper_total %>% 
  mutate(copper_level = if_else(total > mean(total), "high", "low"))

design_total <- model.matrix(~ total, data = count_combine_copper_total)
# # Setting some annotation
mvalue_total <- logit2(as.matrix(combined_resid_win_filter_total))
myAnnotation_total <- cpg.annotate(object = mvalue_total, 
                                   datatype = "array",
                                   what = "Beta",
                                   analysis.type = "differential",
                                   design = design_total,
                                   contrasts = FALSE,
                                   coef = 2,
                                   arraytype = "450K",
                                   fdr = 0.05)

dmr_dmrcate_total <- dmrcate(myAnnotation_total, lambda=1000, C=2)
results.ranges_total <- extractRanges(dmr_dmrcate_total)
results.ranges_total

test_result_total <- results.ranges_total %>% 
  as.data.frame() %>% 
  filter(no.cpgs >= 7)


# ?cpg.annotate()

mvalue_case <- logit2(as.matrix(combined_resid_win_filter_case))
myAnnotation_case <- cpg.annotate(object = mvalue_case, 
                             datatype = "array",
                             what = "M",
                             analysis.type = "differential",
                             design = design_case,
                             coef = ncol(design_case),
                             arraytype = "450K",
                             annotation = c(array = "IlluminaHumanMethylation450k", 
                                            annotation = "ilmn12.hg19"),
                             fdr = 0.05)

myAnnotation_ctrl <- cpg.annotate(object = as.matrix(PEG_NOOB_nors_filter_ctrl), 
                                  datatype = "array",
                                  what = "Beta",
                                  analysis.type = "differential",
                                  design = design_ctrl,
                                  contrasts = FALSE,
                                  coef = 2,
                                  arraytype = "450K",
                                  fdr = 0.05)


str(myAnnotation_case)
# 
# # DMR analysis
dmr_dmrcate_case <- dmrcate(myAnnotation_case, lambda=1000, C=2)

save(dmr_dmrcate_case, file = "dmr_dmrcate_case.RData")
?dmrcate()
?extractRanges()
results.ranges_case <- extractRanges(dmr_dmrcate_case)
results.ranges_case

test <- results.ranges_case %>% 
  as.data.frame() %>% 
  filter(no.cpgs >= 7)

# cols <- c(2,4)[count_combine_copper[[1]]$total]
# names(cols) <-group
# 
# par(mfrow=c(1,1))
# DMR.plot(ranges=results.ranges, dmr=1, 
#          CpGs=as.matrix(PEG_NOOB_nors_win_filter_pd_r), 
#          phen.col=as.numeric(count_combine_copper[[1]]$total),
#          what="Beta", arraytype="450K", genome="hg19")

# method3: champ package
list.dirs(here::here(),recursive = FALSE) %>%
  list.files("\\.RData$", full.names = TRUE, recursive = T) %>%
  grep("bumphunter",.,
       value=TRUE, ignore.case = TRUE) %>%
  # keep(~!str_detect(.x,"win|dmplist|old")) %>%
  map(.,load,.GlobalEnv)


list(
  list(count_combine_copper[[1]], count_combine_copper[[2]], 
       count_combine_copper_total),
  list(combined_resid_win_filter_case, combined_resid_win_filter_ctrl, 
       combined_resid_win_filter_total)
) %>% 
  pmap(function(data1, data2){
    data1 %>% 
      dplyr::select(total) %>% 
      map(~champ.DMR(beta = as.matrix(data2), 
                     pheno = .x, compare.group = NULL, 
                     arraytype = "450K", method = "Bumphunter", minProbes = 5, 
                     adjPvalDmr = 0.05, cores = 10, maxGap=300,
                     cutoff=NULL,
                     pickCutoff=TRUE,
                     smooth=TRUE,
                     smoothFunction=loessByCluster,
                     useWeights=FALSE,
                     permutations=NULL,
                     B=250,
                     nullMethod="bootstrap"))
  }) %>%
  set_names("dmr_bumphunter_champ_copper_case", "dmr_bumphunter_champ_copper_ctrl", 
            "dmr_bumphunter_champ_copper_total") %>%
  list2env(.,envir = .GlobalEnv)

list(
  list(count_combine_copper[[1]], count_combine_copper[[2]], 
       count_combine_copper_total),
  list(combined_resid_win_filter_case, combined_resid_win_filter_ctrl, 
       combined_resid_win_filter_total)
) %>% 
  pmap(function(data1, data2){
    data1 %>% 
      dplyr::select(total) %>% 
      map(~champ.DMR(beta = as.matrix(data2), 
                     pheno = .x, compare.group = NULL, 
                     arraytype = "450K", method = "DMRcate", minProbes = 5, 
                     adjPvalDmr = 0.05, cores = 10,
                     rmSNPCH=T,
                     fdr=0.05,
                     dist=2,
                     mafcut=0.05,
                     lambda=1000,
                     C=2))
  }) %>%
  set_names("dmr_dmrcate_champ_copper_case", "dmr_dmrcate_champ_copper_ctrl", 
            "dmr_dmrcate_champ_copper_total") %>%
  list2env(.,envir = .GlobalEnv)





save(dmr_bumphunter_champ_copper_case, file = "dmr_bumphunter_champ_copper_case.RData")
save(dmr_bumphunter_champ_copper_ctrl, file = "dmr_bumphunter_champ_copper_ctrl.RData")
save(dmr_bumphunter_champ_copper_total, file = "dmr_bumphunter_champ_copper_total.RData")

DMR.GUI(DMR  = dmr_champ_copper_case[[1]], 
        beta = as.matrix(combined_resid_case), 
        pheno = count_combine_copper[[1]]$total)

test <- dmr_bumphunter_champ_copper_total$total$BumphunterDMR
test_dmrcate <- results.ranges_case %>% 
  as.data.frame() %>% 
  filter(no.cpgs >= 5)
# GSEA analysis -----------------------------------------------------------

### method 1: ChAMP package (gometh)

list(
  list(combined_resid_win_filter_case, combined_resid_win_filter_total),
  list(dmp_count_case, dmp_count_total),
  list(dmr_bumphunter_champ_copper_case, dmr_bumphunter_champ_copper_total),
  list(count_combine_copper[[1]],  
       count_combine_copper_total)
) %>% 
  pmap(function(data1, data2, data3, data4){
      champ.GSEA(beta = as.matrix(data1),
                  DMP = data2[[1]],
                  DMR = data3[[1]],
                  CpGlist=NULL,
                  Genelist=NULL,
                  pheno=data4$total,
                  method="gometh",
                  arraytype="450K",
                  Rplot=TRUE,
                  adjPval=0.05,
                  cores = 10)
  }) %>%
  set_names("gsea_champ_copper_case", "gsea_champ_copper_total") %>%
  list2env(.,envir = .GlobalEnv)


save(gsea_champ_copper_case, file = "gsea_champ_copper_case.RData")
save(gsea_champ_copper_total, file = "gsea_champ_copper_total.RData")

### method 2: ChAMP package (empirical bayes)
list(
  list(combined_resid_win_filter_case, 
       combined_resid_win_filter_ctrl, 
       combined_resid_win_filter_total),
  list(count_combine_copper[[1]],
       count_combine_copper[[2]],
       count_combine_copper_total)
) %>% 
  pmap(function(data1, data2){
    champ.GSEA(beta = as.matrix(data1),
               DMP = NULL,
               DMR = NULL,
               CpGlist = NULL,
               Genelist = NULL,
               pheno=data2$total,
               method="ebayes",
               arraytype="450K",
               Rplot=TRUE,
               adjPval=0.05,
               cores = 10)
  }) %>%
  set_names("gsea_champ_copper_ebayes_case", "gsea_champ_copper_ebayes_ctrl",
            "gsea_champ_copper_ebayes_total") %>%
  list2env(.,envir = .GlobalEnv)

save(gsea_champ_copper_ebayes_case, file = "gsea_champ_copper_ebayes_case.RData")
save(gsea_champ_copper_ebayes_ctrl, file = "gsea_champ_copper_ebayes_ctrl.RData")
save(gsea_champ_copper_ebayes_total, file = "gsea_champ_copper_ebayes_total.RData")

test <- gsea_champ_copper_total$DMP %>% 
  rownames_to_column("GO_term")




### method 2: missMethyl package

library(missMethyl)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

topCpGs_case <- meffil_count_op_case$total$analyses$all$table %>% 
  filter(-log10(p.value) > 6) %>% 
  rownames_to_column("cpg") %>%
  arrange(p.value)

sigCpGs <- topCpGs_case %>% 
  pull(cpg)

gst_go <- gometh(sig.cpg=sigCpGs, all.cpg=rownames(PEG_NOOB_nors_filter_case), 
                 collection="GO", plot.bias=TRUE)

topGSA(gst_go, n=10)

gst_kegg <- gometh(sig.cpg=sigCpGs, all.cpg=rownames(PEG_NOOB_nors_filter_case), 
                   collection="KEGG", plot.bias=TRUE)

topGSA(gst_kegg, n=10)

# ebgsea_champ_copper <- champ.ebGSEA(beta = as.matrix(combined_resid_filter_pd_r),
#                                     pheno = count_combine_copper[[1]]$total,
#                                     arraytype = "450K")


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
                  rownames_to_column("cpg")%>% 
                  dplyr::select(cpg, Probe_rs, Probe_maf, UCSC_RefGene_Name, 
                                UCSC_RefGene_Group), by = "cpg") %>% 
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


list(annote_table_op_case$cpg, annote_table_noop_case$cpg) %>% 
  ggVennDiagram(
    category.names = c("A", "B"),
    set_color = c("blue", "red"),
    set_size = 6,
    label = "both",
    label_geom = "label",
    label_alpha = 0,
    label_color = "firebrick"
  ) +
  scale_fill_gradient(low = "grey90",high = "grey60") +
  scale_color_manual(values = c("grey10","grey10"))

test2 <- annote_table_op_case %>% 
  select(cpg, gene, gene_region) %>% 
  filter(cpg %notin% annote_table_noop_case$cpg)

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
              here::here("tables", paste(data2, ".csv", sep = "")))
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
      # map(~filter(.x, p.value < 1E-06 
      #             # & !str_detect(UCSC_RefGene_Name, ";")
      #             )) %>% 
      map(~arrange(.x, p.value) %>%
            slice_head(n = 1000))
  }) %>%
  set_names("ewas_annot_filter_op_case","ewas_annot_filter_op_ctrl",
            "ewas_annot_filter_op_total", "ewas_annot_filter_noop_case",
            "ewas_annot_filter_noop_ctrl", "ewas_annot_filter_noop_total") %>%
  list2env(.,envir = .GlobalEnv)

library(methylGSA)

list(ewas_annot_filter_op_case, ewas_annot_filter_op_ctrl, 
     ewas_annot_filter_op_total, 
     ewas_annot_filter_noop_case, ewas_annot_filter_noop_ctrl, 
     ewas_annot_filter_noop_total) %>% 
  map(function(datalist){
    datalist %>% 
      map(function(data){
        data %>% 
          pull(p.value) %>% 
          set_names(data$cpgs)
      })
  }) %>% 
  set_names("meta.cpg_op_case", "meta.cpg_op_ctrl", "meta.cpg_op_total", 
            "meta.cpg_noop_case", "meta.cpg_noop_ctrl", 
            "meta.cpg_noop_total") %>% 
  list2env(.GlobalEnv)


#function 2: methylRRA


list(meta.cpg_op_case, meta.cpg_op_ctrl, meta.cpg_op_total,
     meta.cpg_noop_case, meta.cpg_noop_ctrl, meta.cpg_noop_total) %>% 
  map(function(cpglist){
    list("GO", "KEGG", "Reactome") %>% 
      map(function(data){
        methylRRA(cpg.pval = cpglist[[1]], 
                  sig.cut = 0.05, 
                  method = "ORA", GS.type = data)
      }) 
  }) %>%
  set_names("gsea_case_op_copper_list", "gsea_ctrl_op_copper_list", 
            "gsea_total_op_copper_list", 
            "gsea_case_noop_copper_list", "gsea_ctrl_noop_copper_list", 
            "gsea_total_noop_copper_list") %>%
  list2env(.GlobalEnv)


gsea_case_op_copper_list <- list("GO", "KEGG", "Reactome") %>% 
  map(function(data){
    methylRRA(cpg.pval = meta.cpg_op_case[[1]], 
              sig.cut = 0.05, 
              method = "GSEA", GS.type = data)
  })
test <- gsea_case_op_copper_list[[1]]
seq(1:3) %>% 
  map(~barplot(gsea_case_op_copper_list[[.x]], num = 10, colorby = "pvalue"))

barplot(gsea_case_op_copper_list[[1]], num = 10, colorby = "pvalue")
barplot(gsea_case_noop_copper_list[[1]], num = 10, colorby = "pvalue")

# Pathway analysis --------------------------------------------------------

library(PANTHER.db)
pthOrganisms(PANTHER.db) <- "HUMAN"
PANTHER.db
columns(PANTHER.db)
keytypes(PANTHER.db)

go_ids <- keys(PANTHER.db,keytype="PATHWAY_ID")
cols <- "GOSLIM_ID"
panther <- mapIds(PANTHER.db, keys=go_ids, column=cols, 
                  as.numeric(cols), keytype="PATHWAY_ID", multiVals="list")
lengths(panther)
panther_inner <- select(PANTHER.db, keys=go_ids, column=cols, 
                        keytype="PATHWAY_ID", jointype = "left")
#ALL CpG sites
#GSEA for panther pathways
# library(org.Hs.eg.db)
# 
# test <- org.Hs.egGO
# 
# mapped_genes <- mappedkeys(test)
# 
# test_list <- as.list(test[mapped_genes])
# 
# if(length(test_list) > 0){
#   got <- test_list[[1]]
#   got[[1]][["GOID"]]
#   got[[1]][["ONTOLOGY"]]
#   got[[1]][["Evidence"]]
# }


library(methylGSA)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
sig_cpg <- dmp_count_total$total %>% 
  filter(-log10(P.Value) > 6) %>% 
  arrange(P.Value) %>% 
  rownames_to_column("cpg")

methylglm_pan_total_op <- methylglm(cpg.pval = sig_cpg %>% 
                                     pull(P.Value) %>% 
                                     set_names(sig_cpg$cpg),
                                GS.list = panther, GS.idtype = "ENTREZID")

test <- gsea_champ_copper_total$DMP %>% 
  arrange(`P.DE`) %>% 
  rownames_to_column("GOSLIM_ID") %>% 
  left_join(panther_inner %>% 
              dplyr::select(PATHWAY_ID, GOSLIM_ID), by = "GOSLIM_ID") %>% 
  distinct()
# ?methylglm()
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

methylglm_path_total_op <- methylglm_pan_total_op %>% 
  left_join(panther.names, by = "ID")


# Visulization ------------------------------------------------------------

library(ggtext)

list(meffil_count_op_case, meffil_count_op_ctrl, meffil_count_op_total,
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
  set_names("metal_list_op_case","metal_list_op_ctrl", "metal_list_op_total",
            "metal_list_noop_case","metal_list_noop_ctrl", "metal_list_noop_total") %>% 
  list2env(.,envir = .GlobalEnv)

list(dmp_count_case, dmp_count_total) %>% 
  map(function(dmp_list){
    dmp_list %>% 
      #keep(., map_lgl(., ~ nrow(.x)>=300)) %>% 
      map(function(data){
        datanew <- data %>% 
          dplyr::rename(chr = CHR,
                        p.value = P.Value,
                        position = MAPINFO) %>% 
          mutate(chr = str_c("chr", chr),
                 cpg = rownames(.))
        
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
  set_names("dmp_list_op_case", "dmp_list_op_total") %>% 
  list2env(.,envir = .GlobalEnv)

list(dmp_list_op_case, dmp_list_op_total) %>% 
  map(function(dmp_list){
    dmp_list %>% 
      map(function(data){
        data %>% 
          group_by(chr) %>% 
          dplyr::summarize(center = mean(global)) %>% 
          ungroup()
      })
  }) %>% 
  set_names("axis_set_op_case",  "axis_set_op_total") %>% 
  list2env(.,envir = .GlobalEnv)


list(dmp_list_op_case, dmp_list_op_total) %>% 
  map(function(dmp_list){
    dmp_list %>% 
      map(function(data){
        data %>% 
          filter(p.value == min(p.value)) %>% 
          mutate(ylim = abs(floor(log10(p.value))) + 1) %>% 
          pull(ylim)
      })
  }) %>% 
  set_names("ylim_op_case", "ylim_op_total") %>% 
  list2env(.,envir = .GlobalEnv)

library(ggrepel)
list(
  list(dmp_list_op_case, axis_set_op_case, ylim_op_case),
  list(dmp_list_op_total, axis_set_op_total, ylim_op_total)
  ) %>% 
  map(function(datalist){
    datalist %>% 
      pmap(function(pest, axis, ylim){
        pest %>% 
          arrange(p.value) %>% 
          mutate(sig = if_else(-log10(p.value) > 6, 1, 0),

                 label = if_else(row_number() <= 10, cpg, "")) %>% 
          ggplot(aes(x = global, y = -log10(p.value), 
                     color = as_factor(sig), size = -log10(p.value))) +
          scale_color_manual(values = c("1" = "red", "0" = "black")) +
          geom_hline(yintercept = -log10(10e-7), color = "red", linetype = "dashed") + 
          geom_point(alpha = 0.75) +
          geom_label_repel(aes(label = label), 
                           box.padding = 1,
                           nudge_x = 0.25,
                           nudge_y = 0.25) +
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
            axis.title.y = element_markdown(face = "bold", size = 15),
            axis.title.x = element_text(face = "bold", size = 15),
            axis.text.y = element_text(size = 15), 
            axis.text.x = element_text(angle = 60, size = 15, vjust = 0.5)
          )
      })
  }) %>% 
  set_names("manhattanlist_op_case", "manhattanlist_op_total") %>% 
  list2env(.,envir = .GlobalEnv)

test <- metal_list_op_case[[1]] %>% 
  mutate(sig = if_else(-log10(p.value) > 6, 1, 0),
         label = if_else(sig == 1, rownames(.), "")) %>% 
  filter(sig == 1)


test2 <- metal_list_op_total[[1]] %>% 
  filter(-log10(p.value) > 6)
  
  mutate(sig = if_else(-log10(p.value) > 6, 1, 0),
         label = if_else(sig == 1, rownames(.), "")) %>% 
  filter(sig == 1)

setdiff(test$label, test2$label)

plotlist_noop <- list(manhattanlist_noop_total[[1]],
                 manhattanlist_noop_case[[1]], 
                 manhattanlist_noop_ctrl[[1]]
                 )

plotlist_op <- list(manhattanlist_op_total[[1]],
                      manhattanlist_op_case[[1]], 
                      manhattanlist_op_ctrl[[1]]
)


#print all plots in one figure
png(file=here::here("figures","manhattan plots for cases and controls_copper_adjop_012824.png"), 
    width = 1920, height = 1080)
ggarrange(plotlist = plotlist_op,
          labels = c("A", "B", "C"),
          font.label = list(size = 20, color = "black"),
          ncol = 2, nrow = 2)
dev.off()


png(file=here::here("figures","manhattan plots for cases and controls_copper_noop_012824.png"), 
    width = 1920, height = 1080)
ggarrange(plotlist = plotlist_noop,
          labels = c("A", "B", "C"),
          font.label = list(size = 20, color = "black"),
          ncol = 2, nrow = 2)
dev.off()


#scatter plot

# plot the top 10 most significantly differentially methylated CpGs in cases

dmp_sig_total_op <- dmp_count_total$total %>% 
  filter(-log10(P.Value) > 6) %>% 
  rownames_to_column("cpg") %>%
  arrange(P.Value)

dmp_sig_case_op <- dmp_count_case$total %>% 
  filter(-log10(P.Value) > 6) %>% 
  rownames_to_column("cpg") %>%
  arrange(P.Value)

dmp_sole_case_op <- dmp_count_case$total %>% 
  filter(-log10(P.Value) > 6) %>% 
  rownames_to_column("cpg") %>%
  arrange(P.Value) %>% 
  filter(cpg %notin% dmp_sig_total_op$cpg)

mean_methylist_copper_new <- list(
  list(PEG_NOOB_nors_win_filter_total, PEG_NOOB_nors_win_filter_case, 
       PEG_NOOB_nors_win_filter_ctrl),
  list(count_combine_copper_total, count_combine_copper[[1]], 
       count_combine_copper[[2]])
) %>% 
  pmap(function(df1, df2){
    dmp_sole_case_op$cpg[1:10] %>% 
      map(function(cpg){
        df1 %>% 
          filter(rownames(.) %in% cpg) %>% 
          pivot_longer(
            cols = starts_with("X"),
            names_to = "sampleid",
            values_to = "beta_val") %>% 
          left_join(datSamplePEG %>% 
                      dplyr::select(sampleid, pegid), 
                    by = "sampleid") %>% 
          left_join(df2, by = "pegid") %>% 
          mutate(pd = if_else(pegid %in% case_ids$pegid, "case", "control"))
      })
  }) %>% 
  set_names("total", "case", "control")


mean_methylist_copper_res_new <- list(
  list(combined_resid_win_filter_total, combined_resid_win_filter_case, 
       combined_resid_win_filter_ctrl),
  list(count_combine_copper_total, count_combine_copper[[1]], 
       count_combine_copper[[2]])
) %>% 
  pmap(function(df1, df2){
    dmp_sole_case_op$cpg[1:10] %>% 
      map(function(cpg){
        df1 %>% 
          filter(rownames(.) %in% cpg) %>% 
          pivot_longer(
            cols = starts_with("X"),
            names_to = "sampleid",
            values_to = "beta_val") %>% 
          left_join(datSamplePEG %>% 
                      dplyr::select(sampleid, pegid), 
                    by = "sampleid") %>% 
          left_join(df2, by = "pegid") %>% 
          mutate(pd = if_else(pegid %in% case_ids$pegid, "case", "control"))
      })
  }) %>% 
  set_names("total", "case", "control")

plotlist_raw_copper_new <- mean_methylist_copper_new %>% 
  map(function(data){
    list(data, dmp_sole_case_op$cpg[1:10]) %>% 
      pmap(function(df, cpg){
        df %>% 
          ggplot(aes(x = total, y = beta_val, color = pd)) + 
          scale_color_colorblind()+
          geom_point() +
          geom_smooth(method = "lm") +
          # geom_jitter()+
          # stat_cor(label.y = 0.54)+
          theme_classic()+
          labs(title = cpg, 
               x = "Copper count",
               y = "Beta values",
               legend = "PD status")
      })
  })


plotlist_res_copper_new <- mean_methylist_copper_res_new %>% 
  map(function(data){
    list(data, dmp_sole_case_op$cpg[1:10]) %>% 
      pmap(function(df, cpg){
        df %>% 
          ggplot(aes(x = total, y = beta_val, color = pd)) + 
          scale_color_colorblind()+
          geom_point() +
          geom_smooth(method = "lm") +
          # geom_jitter()+
          # stat_cor(label.y = 0.54)+
          theme_classic()+
          labs(title = cpg, 
               x = "Copper count",
               y = "Adjusted beta values",
               color = "PD status")
      })
  })


ggarrange(plotlist = plotlist_raw_copper_new[[1]],
          # labels = c("A", "B", "C"),
          ncol = 5, nrow = 2)

library(patchwork)
patchwork::wrap_plots(plotlist_res_copper_new[[1]], ncol = 5,
                      nrow = 2, guides = "collect")

## scatter plot to check the association between case and control methylation levels
# x-axis: dmp beta-value for controls
# y-axis: dmp beta-value for cases
library(ggtext)
# dmp_beta_ctrl <- dmp_count_ctrl$total %>% 
#   rownames_to_column("cpg") %>%
#   select(cpg, B) %>% 
#   rename(B_ctrl = B)
# 
# dmp_beta_case <- dmp_count_case$total %>%
#   rownames_to_column("cpg") %>%
#   select(cpg, B) %>%
#   rename(B_case = B)
B_ctrl <- fit2_ctrl$coefficients %>% 
  as.data.frame()
B_case <- fit2_case$coefficients %>% 
  as.data.frame()

B_total <- fit2_total$coefficients %>% 
  as.data.frame()

dmp_beta_combined <- tibble(
  B_ctrl = fit2_ctrl$coefficients %>% 
    as.data.frame() %>% 
    filter(rownames(.) %in% bicor_total$ID) %>%
    pull(total),
  B_case = fit2_case$coefficients %>% 
    as.data.frame() %>%
    filter(rownames(.) %in% bicor_total$ID) %>%
    pull(total)
)


  

mean_cpg_res_beta_ctrl <- combined_resid_win_filter_ctrl %>% 
  mutate(mean_beta_ctrl = base::rowMeans(dplyr::select(., starts_with("X")))) %>%
  dplyr::select(mean_beta_ctrl)

mean_cpg_beta_ctrl <- PEG_NOOB_nors_win_filter_ctrl %>% 
  mutate(mean_beta_ctrl = base::rowMeans(dplyr::select(., starts_with("X")))) %>%
  dplyr::select(mean_beta_ctrl)

mean_cpg_res_beta_case <- combined_resid_win_filter_case %>%
  mutate(mean_beta_case = base::rowMeans(dplyr::select(., starts_with("X")))) %>%
  dplyr::select(mean_beta_case)

mean_cpg_beta_case <- PEG_NOOB_nors_win_filter_case %>%
  mutate(mean_beta_case = base::rowMeans(dplyr::select(., starts_with("X")))) %>%
  dplyr::select(mean_beta_case)

mean_cpg_beta <- bind_cols(mean_cpg_beta_ctrl, mean_cpg_beta_case)
mean_cpg_res_beta <- bind_cols(mean_cpg_res_beta_ctrl, mean_cpg_res_beta_case)

ggplot(dmp_beta_combined, aes(x = B_ctrl, y = B_case)) +
  geom_point(alpha = 0.75) +
  geom_smooth(method = "lm") +
  stat_cor(label.y = 0.003, size = 6) +
  labs(
    title = "",
    x = "Beta coefficients of EWAS among controls",
    y = "Beta coefficients of EWAS among cases")+
  theme_classic() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(face = "bold", size = 15),
    axis.title.x = element_text(face = "bold", size = 15),
    axis.text.y = element_text(size = 15), 
    axis.text.x = element_text(angle = 60, size = 15, vjust = 0.5)
  )

ggplot(mean_cpg_beta, aes(x = mean_beta_ctrl, y = mean_beta_case)) +
  geom_point(alpha = 0.75) +
  geom_smooth(method = "lm") +
  stat_cor(label.y = 0.96, size = 12) +
  labs(x = "Mean original beta value for controls",
       y = "Mean original beta value for cases")+
  theme_classic() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(face = "bold", size = 15),
    axis.title.x = element_text(face = "bold", size = 15),
    axis.text.y = element_text(size = 15), 
    axis.text.x = element_text(angle = 60, size = 15, vjust = 0.5)
  )

# Making a volcano plot
#making dataset
volcano_dat <- data.frame(foldchange = meffil_count_op_case$total$analyses$all$table$coefficient, 
                  logPvalue =  -log10(meffil_count_op_case$total$analyses$all$table$p.value))
volcano_dat$threshold <- as.factor(abs(volcano_dat$foldchange) < 0.001)


#Visualization
cols <- c("TRUE" = "grey", "FALSE" = "blue")
ggplot(data=volcano_dat, aes(x=foldchange, y = logPvalue, color = threshold)) +
  geom_point(alpha=.6, size=1.2) +
  scale_colour_manual(values = cols) +
  geom_vline(xintercept = 0.001, colour="#990000", linetype="dashed") +
  geom_vline(xintercept = -0.001, colour="#990000", linetype="dashed") +
  theme(legend.position="none") +
  xlab("Methylation difference") +
  ylab("-log10 p value") +
  theme_bw() +
  theme(legend.position = "none")


mean_methylist_copper <- list(
  list(combined_resid_total, combined_resid_filter_pd_r, 
       combined_resid_filter_ctrl_r),
  list(count_combine_copper_total, count_combine_copper[[1]], 
       count_combine_copper[[2]])
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
                  dplyr::select(sampleid, pegid), 
                by = "sampleid") %>% 
      left_join(df2, by = "pegid")
  }) 

mean_methylist_raw_copper <- list(
  list(peg_noob_nors_win_total, PEG_NOOB_nors_win_filter_pd_r, 
       PEG_NOOB_nors_win_filter_ctrl_r),
  list(count_combine_copper_total, count_combine_copper[[1]], 
       count_combine_copper[[2]])
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
                  dplyr::select(sampleid, pegid), 
                by = "sampleid") %>% 
      left_join(df2, by = "pegid")
  }) 

test <- mean_methylist_copper[[1]]

mean_methylist_op <- list(
  list(combined_resid_total, combined_resid_filter_pd_r, 
       combined_resid_filter_ctrl_r),
  list(count_combine_op_total, count_combine_op[[1]], 
       count_combine_op[[2]])
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
      left_join(df2, by = "pegid")
  }) 

mean_methylist_raw_op <- list(
  list(peg_noob_nors_win_total, PEG_NOOB_nors_win_filter_pd_r, 
       PEG_NOOB_nors_win_filter_ctrl_r),
  list(count_combine_op_total, count_combine_op[[1]], 
       count_combine_op[[2]])
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
                  dplyr::select(sampleid, pegid), 
                by = "sampleid") %>% 
      left_join(df2, by = "pegid")
  }) 

plotlist_copper <- mean_methylist_copper %>% 
  map(function(data){
    data %>% 
      ggplot(aes(x = total, y = mean_methyl)) + 
      geom_point() +
      geom_smooth(method = "lm") +
      # geom_jitter()+
      stat_cor(label.y = 0.54)+
      theme_classic()+
      labs(x = "Copper count",
           y = "Mean methylation residual")
  })

plotlist_raw_copper <- mean_methylist_raw_copper %>% 
  map(function(data){
    data %>% 
      ggplot(aes(x = total, y = mean_methyl)) + 
      geom_point() +
      geom_smooth(method = "lm") +
      # geom_jitter()+
      stat_cor(label.y = 0.54)+
      theme_classic()+
      labs(x = "Copper count",
           y = "Mean methylation level")
  })

plotlist_op <- mean_methylist_op %>% 
  map(function(data){
    data %>% 
      ggplot(aes(x = total, y = mean_methyl)) + 
      geom_point() +
      geom_smooth(method = "lm") +
      # geom_jitter()+
      stat_cor(label.y = 0.54)+
      theme_classic()+
      labs(x = "OP count",
           y = "Mean methylation residual")
  })

plotlist_raw_op <- mean_methylist_raw_op %>% 
  map(function(data){
    data %>% 
      ggplot(aes(x = total, y = mean_methyl)) + 
      geom_point() +
      geom_smooth(method = "lm") +
      # geom_jitter()+
      stat_cor(label.y = 0.54)+
      theme_classic()+
      labs(x = "OP count",
           y = "Mean methylation level")
  })

ggarrange(plotlist = plotlist_raw_op,
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