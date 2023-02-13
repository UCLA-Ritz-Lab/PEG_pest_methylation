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

list(c_lb_sd_case_wt_10_new,
     r_lb_sd_case_wt_10_new) %>% 
  map(function(data){
    data %>% 
      select(all_of(myvar))
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
  select(-c(sampleid,pegid)) %>%
  tbl_strata(
    strata = pdstudystudy,
    .tbl_fun =
      ~ .x %>%
      tbl_summary(by = pd_new, missing = "no",
                  #type=list(c(female) ~ "categorical"),
                  statistic = list(all_continuous() ~ "{mean} ({sd})",
                                   all_categorical() ~ "{n} ({p}%)"),
                  digits = list(all_continuous() ~ 1,
                                all_categorical() ~ 1)) %>%
      modify_header(label = "**Characteristics**") %>%
      modify_spanning_header(starts_with("stat_") ~ "**PD status**") %>%
      modify_caption("Table 1. Demographic characteristics of PEG cases (N = 569)") %>%
      bold_labels() 
  ) %>% 
  table1()

#calculate the n and % of exposure for each chemical: combine
dur_long_combine <- list(c_lb_sd_case_wt_10_new, r_lb_sd_case_wt_10_new) %>% 
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
  mutate(pct = n_exp/569) %>% 
  left_join(chemlist_new, by="chemcode") %>% 
  left_join(chem_class, by = c("chemname","chemcode")) %>% 
  relocate(chemname,chemcode,`chem class (pan)`)


#heavy metal chemical use

list(exp_window_address_case_c, exp_window_address_case_r) %>% 
  map(function(data){
    data %>% 
      filter(chemcode %in% heavy_metal$chemcode) %>%
      # mutate_at(vars(sum_total_lbs), dec_out_na) %>% 
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
  group_by(year, location) %>% 
  summarise(y=mean(sum_total_lbs, na.rm=T),.groups="keep") %>% #change to average per person, including unexposed people
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
         legend.title = element_text(size = 18),
         legend.text = element_text(size = 15),
         strip.text = element_text(size=13)) 

# DMP --------------------------------------------------------------------

#Differentially methylated positions (DMP)
# Load meffil and set how many cores to use for parallelization

source(here("scripts", "mclapply.hack.R"))
source(here("scripts", "meffil_fixed.R"))
source(here("scripts", "meffil_report_fixed.R"))

#Method 1: ChAMP
library(ChAMP)

list(list(c_lb_sd_case_wt_10_new,r_lb_sd_case_wt_10_new),
     list(combined_resid_filter_pd_c,combined_resid_filter_pd_r)) %>% 
  pmap(function(data1,data2){
    data1 %>% 
      select(starts_with("chem")) %>% 
      map(~champ.DMP(beta = data2, pheno = .x, 
                     compare.group = NULL, 
                     adjPVal = 0.05, adjust.method = "BH", 
                     arraytype = "450K")[[1]])
  }) %>% 
  set_names("champ_dmplist_case_heavymetal_sig_c","champ_dmplist_case_heavymetal_sig_r") %>% 
  list2env(.,envir = .GlobalEnv)

save(champ_dmplist_case_heavymetal_sig_c, file = "champ_dmplist_case_heavymetal_sig_c.RData")
save(champ_dmplist_case_heavymetal_sig_r, file = "champ_dmplist_case_heavymetal_sig_r.RData")


# Annotation and GSEA -----------------------------------------------------

list(champ_dmplist_case_heavymetal_sig_c,champ_dmplist_case_heavymetal_sig_r) %>% 
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


list(list(ewas_annot_c,heavy_metal$chemcode),
     list(ewas_annot_r,heavy_metal$chemcode)) %>% 
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
        methylRRA(cpg.pval = data, method = "GSEA")
      })
  }) %>% 
  set_names("methylRRA_list_c","methylRRA_list_r") %>% 
  list2env(.,envir = .GlobalEnv)

methylRRA_chem354_c <- methylRRA(cpg.pval = meta.cpg_list_c[["chem354"]], method = "GSEA")
methylRRA_chem251_c <- methylRRA(cpg.pval = meta.cpg_list_c[["chem251"]], method = "GSEA")
methylRRA_chem599_c <- methylRRA(cpg.pval = meta.cpg_list_c[["chem599"]], method = "GSEA")

barplot(methylRRA_chem354_c, num = 10, colorby = "pvalue")
barplot(methylRRA_chem251_c, num = 10, colorby = "pvalue")
barplot(methylRRA_chem599_c, num = 10, colorby = "pvalue")

methylRRA_chem354_r <- methylRRA(cpg.pval = meta.cpg_list_r[["chem354"]], method = "GSEA")
methylRRA_chem251_r <- methylRRA(cpg.pval = meta.cpg_list_r[["chem251"]], method = "GSEA")
methylRRA_chem599_r <- methylRRA(cpg.pval = meta.cpg_list_r[["chem599"]], method = "GSEA")

barplot(methylRRA_chem354_r, num = 10, colorby = "pvalue")
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
list(quote_all(chem354,chem251,chem599),
     quote_all(cg07699440, cg16282242, cg24896860))%>% 
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
list(quote_all(chem354,chem251,chem599),
     quote_all(cg00540866, cg26842720, cg09010802))%>% 
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

