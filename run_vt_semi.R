rm(list=ls())
options(warn = 2)

library(rpart)
library(rpart.plot)

source("./utils.R")

collect_exp_res_VT = function(vt.sbgrps, gt_te_test, d_test){
  # What to collect: 
  # - The accuracy of the max subgroup and all other subgroups, in terms of individual treatment effects
  # - The accuracy of the max subgroup and all other subgroups, in terms of the subgroup treatment effect
  # - The variance, estimated from the gt_te of the test set, of the treatment effect of each subgroup
  
  subgroups = vt.sbgrps$Subgroup
  treatment_effects = rep(0, length(subgroups))
  gt_subgroup_effects = rep(0, length(subgroups))
  max_te = -Inf
  subgroup_gt_treatment_effects_var = rep(0, length(subgroups))
  
  vt_treatment_effects = rep(0, length(subgroups))
  vt_treatment_effects_var = rep(0, length(subgroups))
  diff_average = rep(0, length(subgroups)) # the difference between the subgroup te and the average gt_tex
  diff_individual <- vector("list", length(subgroups)) # the individual difference of the subgroup_te and the gt_te
  
  
  for(i in 1:length(subgroups)){
    sg = subgroups[i]
    items = strsplit(sg, " & ")[[1]]
    bool_test = rep(T, nrow(d_test))
    for(j in 1:length(items)){
      item = items[j]
      if(grepl(">", item)){
        feature_value = strsplit(item, ">=")[[1]]
        operator = ">="
      } else {
        feature_value = strsplit(item, "< ")[[1]]
        operator = "<"
      }
      feature = feature_value[1]
      value = as.numeric(feature_value[2])
      if(operator == ">="){
        bool_test = bool_test & (d_test[[feature]] >= value)
      } else {
        bool_test = bool_test & (d_test[[feature]] < value)
      }
    }
    subgroup_te = mean(d_test$Y[bool_test & d_test$T == 1]) - 
      mean(d_test$Y[bool_test & d_test$T == 0])
    treatment_effects[i] = subgroup_te
    subgroup_gt_treatment_effects_var[i] = var(gt_te_test[bool_test])
    diff_average[i] = te - mean(gt_te_test[bool_test])
    diff_individual[[i]] = te - gt_te_test[bool_test] 
    
    if(subgroup_te > max_te){
      max_te = subgroup_te
      abs_error_te = abs(max_te - gt_te)
      max_sg_bool = bool_test
    }
  }
  # collect the results
  which_max = which.max(treatment_effects)
  return(list(
    subgroup_treatment_effects = treatment_effects,
    subgroup_gt_treatment_effects_var = subgroup_gt_treatment_effects_var,
    subgroup_diff_average = diff_average,
    subgroup_diff_individual = diff_individual,
    max_subgroup_treatment_effect = tree_treatment_effects[which_max],
    which_max = which_max,
    max_subgroup_gt_treatment_effct = mean(gt_te_test[max_sg_bool])
  ))
}

fit_vt_model = function(d_train){
  tree_type = ifelse(length(unique(d_train$Y)) > 2, "reg", "class")
  vt.obj <- vt.data(dataset         = d_train,
                    outcome.field   = "Y",
                    treatment.field = "treatment",
                    interactions    = TRUE)
  # Run the VT method
  ## First step : create random forest model
  vt.for <- vt.forest(forest.type  = "one",
                      vt.data      = vt.obj,
                      interactions = TRUE,
                      ntree        = 500)
  ## Second step : find rules in data 
  vt.trees <- vt.tree(tree.type = tree_type,
                      vt.difft  = vt.for, 
                      threshold = quantile(vt.for$difft, seq(.5,.8,.1)),
                      maxdepth  = ncol(d) - 2)
  ## Print results
  vt.sbgrps <- vt.subgroups(vt.trees)
  return(vt.sbgrps)
}

semi_datasets_folder = "../causal_sub_datasets/clean_data/"
files = list.files(semi_datasets_folder, pattern = "*.csv", full.names = TRUE)

res_df = data.frame(
  dataset = character(),
  subgroup_treatment_effects = I(list()),
  subgroup_gt_treatment_effects_var = I(list()),
  subgroup_diff_average = I(list()),
  subgroup_diff_individual = I(list()),
  max_subgroup_treatment_effect = numeric(),
  which_max = integer(),
  max_subgroup_gt_treatment_effct = numeric()
)

res_df_summary = data.frame(
  dataset = character(),
  max_subgroup_treatment_effect = numeric(),
  subgroup_diff_average_weighted_mean = numeric(),
  subgroup_diff_average_unweighted_mean = numeric(),
  max_subgroup_diff_average = numeric(),
  max_subgroup_gt_treatment_effct = numeric()
)

for(i in 1:length(files)){
  # if(files[i] == "../causal_sub_datasets/clean_data//IBM_scaling_41.csv"){
  #   # browser()
  # } else{
  #   next
  # }
  data_file = files[i]
  cat("running VT on semi-synthetic dataset: ", data_file, "\n")
  
  read_and_process_data_res = read_and_process_semi_data(data_file)
  d = read_and_process_data_res$d
  gt_te = read_and_process_data_res$gt_te
  
  # train/test split
  set.seed(1)
  train_indices = sample(1:nrow(d), size = 0.5 * nrow(d), replace = FALSE)
  
  d_train = d[train_indices, ]
  d_test = d[-train_indices, ]
  gt_te_test = gt_te[-train_indices]
  
  # fit a VT model
  vt.sbgrps = fit_vt_model(d_train)
  
  
  
  exp_res = collect_exp_res_VT(vt.sbgrps, gt_te_test, d_test)
  
  
  res_df = rbind(res_df, data.frame(
    dataset = basename(data_file),
    subgroup_treatment_effects = I(list(exp_res$subgroup_treatment_effects)),
    subgroup_gt_treatment_effects_var = I(list(exp_res$subgroup_gt_treatment_effects_var)),
    subgroup_diff_average = I(list(exp_res$subgroup_diff_average)),
    subgroup_diff_individual = I(list(exp_res$subgroup_diff_individual)),
    max_subgroup_treatment_effect = exp_res$max_subgroup_treatment_effect,
    which_max = exp_res$which_max, 
    max_subgroup_gt_treatment_effct = exp_res$max_subgroup_gt_treatment_effct
  ))
  
  subgroup_coverage = sapply(exp_res$subgroup_diff_individual, length)
  
  res_df_summary = rbind(res_df_summary, data.frame(
    dataset = basename(data_file),
    max_subgroup_treatment_effect = exp_res$max_subgroup_treatment_effect,
    subgroup_diff_average_weighted_mean = weighted.mean(exp_res$subgroup_diff_average, 
                                                        w = subgroup_coverage, 
                                                        na.rm = TRUE),
    subgroup_diff_average_unweighted_mean = mean(exp_res$subgroup_diff_average, na.rm = TRUE),
    max_subgroup_diff_average = 
      exp_res$subgroup_diff_average[exp_res$which_max],
    max_subgroup_gt_treatment_effct = exp_res$max_subgroup_gt_treatment_effct
  ))
  if(i %% 10 == 0){
    write.csv(res_df, 
              file = paste0("./semi_VT_results_full.csv"), 
              row.names = FALSE, 
              quote = TRUE)
    write.csv(res_df_summary,
              file = paste0("./semi_VT_results.csv"), 
              row.names = FALSE, 
              quote = TRUE)
  }
}

require(data.table)
exp_df = fread("./semi_VT_results.csv", 
               stringsAsFactors = FALSE, 
               na.strings = c("", "NA", "NULL"))
# round up the numbers



