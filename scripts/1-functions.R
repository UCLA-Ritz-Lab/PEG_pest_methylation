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
  # "meffil",      #normalization: meffil.normalize.dataset(), ewas: meffil.ewas()
  # "DMRcate",      #For DMR analysis: dmrcate()
  # "WGCNA",       #weighted correlation network analysis: GOenrichmentAnalysis()
  #"ChAMP",       #champ.DMP(), champ.DMR()
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

isnt_out_z <- function(x, thres = 3, na.rm = TRUE) {
  abs(x - mean(x, na.rm = na.rm)) <= thres * sd(x, na.rm = na.rm)
}
isnt_out_mad <- function(x, thres = 3, na.rm = TRUE) {
  abs(x - median(x, na.rm = na.rm)) <= thres * mad(x, na.rm = na.rm)
}
isnt_out_tukey <- function(x, k = 1.5, na.rm = TRUE) {
  quar <- quantile(x, probs = c(0.25, 0.75), na.rm = na.rm)
  iqr <- diff(quar)
  
  (quar[1] - k * iqr <= x) & (x <= quar[2] + k * iqr)
}

isnt_out_funs <- list(
  z2 = isnt_out_z,
  mad = isnt_out_mad,
  tukey = isnt_out_tukey
)

dec_out_na <- function(x, thres = 3, na.rm = TRUE) {
  sd <- sd(x, na.rm = na.rm)
  mean <- mean(x, na.rm = na.rm)
  id_out = x > mean + thres * sd | x < mean - thres * sd
  x[id_out] = NA
  return(x)
}


dec_out <- function(x, thres = 3, na.rm = TRUE) {
  sd <- sd(x, na.rm = na.rm)
  mean <- mean(x, na.rm = na.rm)
  id_p = x > mean + thres * sd
  id_n = x < mean - thres * sd
  x[id_p] = mean + thres * sd
  x[id_n] = mean - thres * sd
  return(x)
}

#--------------------------------End of the code--------------------------------