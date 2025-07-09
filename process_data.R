# Most, if not all, datasets come from https://github.com/rguo12/awesome-causality-data
rm(list = ls())
require(dplyr)
require(data.table)
require(ggplot2)
setwd("~/projects/causal_sub_datasets/")
#mkdir if not exists
if (!dir.exists("clean_data")) {
  dir.create("clean_data")
}

process_ihdp = function(){
  d = read.csv("./IHDP/ihdp_data.csv")
  d = d[,c(4:ncol(d), 1, 2, 3)]
  d$mu0 = NULL
  d$mu1 = NULL
  d$treatment = as.numeric(as.factor(d$treatment)) - 1
  write.csv(d, file = "clean_data/IHDP.csv", row.names = FALSE, quote = FALSE)
}

process_IBM_scaling_single = function(){
  d = read.csv("./ibm-causal-data/LBIDD/x.csv")
  d_factual = read.csv("./ibm-causal-data/LBIDD/scaling/ff695a5ff5464ee1b7deda74ba5425d7.csv")
  d_cf = read.csv("./ibm-causal-data/LBIDD/scaling/ff695a5ff5464ee1b7deda74ba5425d7_cf.csv")
  # re-order the rows by sample_id
  d = d[order(d$sample_id), ]
  d_factual = d_factual[order(d_factual$sample_id), ]
  
  dd = left_join(d_factual, d_cf, by = "sample_id")
  
  # match the rows by sample_id, and keep only the rows that are in both dataframes
  dd <- left_join(dd, d, by = "sample_id")
  dd$sample_id = NULL
  dd = dd[,c(5:ncol(dd), 1:4)]
  # write.csv(dd, file = "clean_data/IBM_scaling_single.csv", row.names = FALSE, quote = FALSE)
}

process_IBM_scaling_all = function(){
  # Semi-synthetic datasets based on Linked Births and Infant Deaths Database (LBIDD),
  #.  widely referred to as Twins. 
  d = read.csv("./ibm-causal-data/LBIDD/x.csv")
  # list files
  dir = "./ibm-causal-data/LBIDD/scaling/"
  files = list.files(dir, pattern = "*.csv", full.names = TRUE)
  for(i in 1:length(files)){
    if(i %% 2 == 0) {
      next
    }
    cat("Processing file ", i, " of ", length(files), "\n")
    
    cf_file = files[i]
    factual_file = files[i + 1]
    
    # check if the file name matches
    stopifnot(cf_file == paste0(strsplit(factual_file, ".csv")[1], "_cf.csv"))
    
    d_factual = read.csv(factual_file)
    d_cf = read.csv(cf_file)
    
    # re-order the rows by sample_id
    d = d[order(d$sample_id), ]
    d_factual = d_factual[order(d_factual$sample_id), ]
    
    dd = left_join(d_factual, d_cf, by = "sample_id")
    
    # match the rows by sample_id, and keep only the rows that are in both dataframes
    dd <- left_join(dd, d, by = "sample_id")
    dd$sample_id = NULL
    dd = dd[,c(5:ncol(dd), 1:4)]
    
    dd$treatment = dd$z
    dd$y_factual = ifelse(dd$treatment == 1, dd$y1, dd$y0)
    dd$y_cfactual = ifelse(dd$treatment == 1, dd$y0, dd$y1)
    dd$y0 = NULL
    dd$y1 = NULL
    dd$z = NULL
    dd$y = NULL
    write.csv(dd, file = paste0("clean_data/IBM_scaling_",i, ".csv"),
              row.names = FALSE, quote = FALSE)
  }
}



process_acic2016 = function(){
  # if (require("remotes", quietly = TRUE) == FALSE) {
  #   install.packages("remotes")
  #   require("remotes")
  # }
  # remotes::install_github("vdorie/aciccomp/2016")
  require(aciccomp2016)
  for(i in 1:nrow(parameters_2016)){
    cat("Processing ACIC2016 dataset ", i, " of ", nrow(parameters_2016), "\n")
    acic2016 = dgp_2016(input_2016, parameters_2016[i,], random.seed = 1)
    dd2 = input_2016
    # one-hot encode dd$x2 and dd$x_24
    library(fastDummies)
    
    dd2$x_2 = as.factor(dd2$x_2)
    dd2$x_24 = as.factor(dd2$x_24)
    dd <- dummy_cols(dd2, select_columns = c("x_2", "x_21" ,"x_24"), remove_first_dummy = F)
    # remove the original columns as their one-hot encoded version have been included.
    dd$x_2 = NULL
    dd$x_24 = NULL
    dd$x_21 = NULL
    
    dd$treatment = acic2016$z
    dd$y0 = acic2016$y.0
    dd$y1 = acic2016$y.1
    
    dd$y_factual = ifelse(dd$treatment == 1, dd$y1, dd$y0)
    dd$y_cfactual = ifelse(dd$treatment == 1, dd$y0, dd$y1)
    
    dd$y0 = NULL
    dd$y1 = NULL
    
    write.csv(dd, file = paste0("clean_data/ACIC2016_", i, ".csv"), 
              row.names = FALSE, quote = FALSE)
  }
}

process_acic2016()

# process_IBM_scaling_all()
# process_ihdp()
