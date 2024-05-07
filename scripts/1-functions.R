## ---------------------------
##
## Script name: 1-functions.R
##
## Purpose of script: To create functions that will be used in the analysis
##
## Author: Yufan Gong
##
## Date Created: 2024-05-07
##
## Copyright (c) Yufan Gong, 2024
## Email: ivangong@ucla.edu
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------


# 1. Check the directory --------------------------------------------------

#please open all the R scripts in PEG_pest_methylation.Rproj

{
  getwd()
  options(mc.cores = 9)
}


# 2. Loading packages -----------------------------------------------------

{
  pacman::p_load(
    #For creating tables
    "kableExtra",  #create amazing tables: kbl()
    "skimr",       #summary statistics: skim()
    #"dataxray",    #explore data: report_xray()
    "arsenal",     #create tables: tableby()
    "expss",       #create contingency tables: calc_cro_cpct()  
    "huxtable",    #as_hux_table()
    "flextable",   #as_flex_table()
    "gtsummary",   #create amazing tables: tbl_summary()
    
    #For loading data
    "readxl",     #read in excel data: read_xlsx
    "haven",      #read in sas data: read_sas()
    "here",       #setting the directory in the project: here()
    
    #For manipulating data
    "rlang",      #for Non-standard evaluation: eval(), expr(), ensym(), caller_env(), exec(), !!
    "magrittr",   #for the pipe operator: %>% and %<>%
    "lubridate",  #for manipulating dates: intervals(), durations()
    "labelled",   #labelleling the data: set_variable_labels(), set_value_labels()
    
    # Enhancing plots
    "scales",      #makes easy to format percent, dollars, comas: percent()
    "ggalt",       #makes easy splines: geom_xsplines()
    "ggeasy",      #applies labels among other things: easy_labs()
    "gridExtra",   #combining plots and tables on plots: grid.arrange(), tableGrob()
    "ggpubr",      #combines plots: ggarrange()
    "ggthemes",     #blind colors
    "ggVennDiagram", #venn diagram
    "Amelia",      #check missing pattern: missmap()
    "pheatmap",    #create pretty heatmap: pheatmap(),
    "ggnewscale",  #create multiple scales: new_scale()
    "ggrepel",     #avoid overlapping labels: geom_text_repel()
    "DataExplorer",#create a report of the data: create_report()
    "patchwork",   #combines plots: wrap_plots(), plot_layout()
    "plotly",      #create interactive plots: ggplotly()
    
    # Other great packages
    "glue",        #replaces paste: glue()
    "Hmisc",       #explore the data: describe()
    "mice",        #imput missing data: mice()
    "gmodels",     #create contigency table: CrossTable()
    "meta",        #meta models: metagen()
    "codebook",    #amazing package to set labels: dict_to_list()
    "foreach",     #executing R code repeatedly: foreach(), %do%, %dopar%
    "doParallel",  #Provides a parallel backend for the %dopar% function: registerDoParallel(),
    "doFuture",    #Provides a parallel backend for the %dofuture% function: registerDoFuture()
    "future",      #Provides a future backend for the %dofuture% function: plan()
    "furrr",       #Provides a purrr like syntax for future: future_map(), future_pmap()
    
    # For analysis
    #"lme4",        #linear mixed-effects model: lmer()
    #"nlme",        #linear and nonlinear mixed effects:
    # "minfi",       #getBeta(), getSex()
    # "meffil",      #normalization: meffil.normalize.dataset(), ewas: meffil.ewas()
    # "DMRcate",      #For DMR analysis: dmrcate()
    # "WGCNA",       #weighted correlation network analysis: GOenrichmentAnalysis()
    # "ChAMP",       #champ.DMP(), champ.DMR()
    #"mixOmics",    #block.plsda()
    
    #For data cleaning
    "tidyverse"   #data manipulation and visualization:select(), mutate()
  )
}



# 3. Create sub folders ---------------------------------------------------


# c("script","rmd","tables","figures","data") %>%
#   map(dir.create)


# 4. clean env ------------------------------------------------------------

rm(list = ls())

# 5. create functions -----------------------------------------------------

{
  quote_all <- function(...){
    args<-rlang::ensyms(...)
    paste(purrr::map(args,as_string),sep = "")
  }
  
  `%notin%` <- Negate(`%in%`)
  
  # create table1
  table1 <- function(table) {
    
    table %>% 
      as_hux_table() -> hux
    
    table %>% 
      as_flex_table() -> flex
    
    return(list(hux=hux, flex=flex))
    
  }
  
  #create functions to remove outliers
  
  ## 1. use IQR, an alternative of 1
  extreme_remove_iqr <- function(x) {
    Q = quantile(x, c(0.25, 0.75), na.rm = TRUE)
    iqr = IQR(x, na.rm = TRUE)
    x = if_else(x < Q[1] - 1.5 * iqr | x > Q[2] + 1.5 * iqr, NA, x)
  }
  
  ## 2. use percentile/thousands: remove
  extreme_remove_percentile <- function(x) {
    Q = quantile(x, c(0.001, 0.999), na.rm = TRUE)
    x = if_else(x < Q[1] | x > Q[2], NA, x)
  }
  
  
  ## 3. use percentile/thousands: winsorization
  extreme_remove_percentile_win <- function(x) {
    Q = quantile(x, c(0.1, 0.90), na.rm = TRUE)
    x = case_when(x < Q[1] ~ Q[1],
                  x > Q[2] ~ Q[2],
                  TRUE ~ x)
  }
  
  # an equivalent as the third function
  dec_out <- function(x, na.rm = TRUE) {
    Q = quantile(x, c(0.05, 0.95), na.rm = na.rm)
    id_p = x > Q[2]
    id_n = x < Q[1]
    x[id_p] = Q[2]
    x[id_n] = Q[1]
    return(x)
  }
  
  # different data processing methods
  create_quantile <- function(data){
    data %>% 
      mutate_at(vars(starts_with("chem")), 
                ~if_else(.x > quantile(.x, 
                                       c(0.25, 0.75), na.rm = TRUE)[2], 
                         "high", "low"))
  }
  
  create_evernever <- function(data){
    data %>% 
      mutate_at(vars(starts_with("chem")), 
                ~if_else(.x > 0, "ever", "never"))
  }
  
  
  create_ztrans <- function(data){
    data %>% 
      mutate_at(vars(starts_with("chem")),
                ~scale(., center = TRUE)
                %>% as.vector())
  }
  
  create_counts <- function(data){
    data %>% 
      as_tibble()
  }
}



#--------------------------------End of the code--------------------------------