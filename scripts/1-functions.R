#------------------------------------------------------------------------------#
#-------------------------------Created by Yufan-------------------------------#
#------------------------------Date: 01/07/2022--------------------------------#
#-------------------------------To create functions----------------------------#
#------------------------------------------------------------------------------#


# 1. Check the directory --------------------------------------------------

#please open all the R scripts in PEG_Methylation.Rproj

# 2. Loading packages -----------------------------------------------------
# library(usethis) 
# usethis::edit_r_environ()
options(mc.cores = 9)
pacman::p_load(
  #For creating tables
  "kableExtra",  #create amazing tables: kbl()
  "skimr",       #summary statistics: skim()
  "dataxray",    #explore data: report_xray()
  "arsenal",     #create tables: tableby()
  "expss",       #create contingency tables: calc_cro_cpct()  
  "huxtable",    #as_hux_table()
  "flextable",   #as_flex_table()
  "gtsummary",   #create amazing tables: tbl_summary()
  
  #For manipulating data
  "rlang",       #for Non-standard evaluation: eval(), expr(), ensym(), caller_env(), exec(), !!
  "magrittr",    #for the pipe operator: %>% and %<>%
  "broom",       #for tidying up the results of a regression: tidy()
  "lubridate",   #for manipulating dates: intervals(), durations()
  "labelled",    #labelleling the data: set_variable_labels(), set_value_labels()
  
  # Enhancing plots
  #"scales",      #makes easy to format percent, dollars, comas: percent()
  #"ggalt",       #makes easy splines: geom_xsplines()
  #"ggeasy",      #applies labels among other things: easy_labs()
  "gridExtra",   #combining plots and tables on plots: grid.arrange(), tableGrob()
  "ggpubr",      #combines plots: ggarrange()
  "ggthemes",     #blind colors
  #"Amelia",      #check missing pattern: missmap()
  #"pheatmap",    #create pretty heatmap: pheatmap()
  #"qqman",       #create Manhattan plots: manhattan()
  
  # Other great packages
  "glue",        #replaces paste: glue()
  "Hmisc",       #explore the data: describe()
  "mise",        #clear environment space: mise()
  "gmodels",     #create contigency table: CrossTable()
  "meta",        #meta models: metagen()
  "codebook",    #amazing package to set labels: dict_to_list()
  
  # For analysis
  #"lme4",        #linear mixed-effects model: lmer()
  #"nlme",        #linear and nonlinear mixed effects:
  # "minfi",       #getBeta(), getSex()
  "meffil",      #normalization: meffil.normalize.dataset(), ewas: meffil.ewas()
  # "DMRcate",      #For DMR analysis: dmrcate()
  # "WGCNA",       #weighted correlation network analysis: GOenrichmentAnalysis()
  "ChAMP",       #champ.DMP(), champ.DMR()
  #"mixOmics",    #block.plsda()
  
  #For loading data
  "readxl",      #read in excel data: read_xlsx
  "haven",       #read in sas data: read_sas()
  "here",        #setting the directory in the project: here()
  
  #For data manipulation
  "tidyverse"   #data manipulation and visualization:select(), mutate()
)
#invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE))
# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
# BiocManager::install("PANTHER.db")
# install.packages("devtools")
# library(devtools)
# install_github("perishky/meffil",force = TRUE)
# remotes::install_github("perishky/meffil")
# library(minfi)
# library(minfiData)
# library(sva)
# library(meffil)
# library(mixOmics)

# 3. Create sub folders ---------------------------------------------------


# c("script","rmd","tables","figures","data") %>%
#   map(dir.create)


# 4. clean env ------------------------------------------------------------

mise()


# 5. create functions -----------------------------------------------------


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
  Q = quantile(x, c(0.001, 0.999), na.rm = TRUE)
  x = case_when(x < Q[1] ~ Q[1],
                x > Q[2] ~ Q[2],
                TRUE ~ x)
}
#--------------------------------End of the code--------------------------------