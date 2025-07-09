rm(list = ls())
library(rpart)
library(rpart.plot)

set.seed(1)

source("./util_synthetic.R")

require(rsides)
col_type <- function(x) {
  u <- unique(x)
  
  if (length(u) == 2) {                # exactly two values → binary
    "nominal"
  } else {                             # everything else
    "numeric"
  }
}


get_subgroup_bool = function(dd, rsides_results, gt_max_bool){
  # dd can be d_train or d_test
  # rsides_results is the output of SubgroupSearch of the RSIDES package
  # gt_max_bool is a boolean vector indicating the rows of dd that are in the ground_truth subgroup
  subgroups = rsides_results$patient_subgroups$Subgroups$subgroup
  treatment_effects = rep(0, length(subgroups))
  max_te = -Inf
  for(i in 1:length(subgroups)){
    sg = subgroups[[i]]
    if(sg == "Overall population"){
      bool = rep(TRUE, nrow(dd))
    } else {
      bool = rep(TRUE, nrow(dd))
      rules = strsplit(sg, ":")[[1]]
      for(j in 1:length(rules)){
        rule = rules[j]
        # get the operator from ">", "<", ">=", "<=", "="
        operator <- regmatches(rule, regexpr("(>=|<=|>|<|=)", rule))
        
        feature_value = strsplit(rule, operator)[[1]]
        feature = feature_value[1]
        value = as.numeric(feature_value[2])
        bool = bool & switch(
          operator,
          ">" = dd[[feature]] > value,
          "<" = dd[[feature]] < value,
          ">=" = dd[[feature]] >= value,
          "<=" = dd[[feature]] <= value,
          "=" = dd[[feature]] == value
        )
      }
    }
    # get the treatment effects
    bool_t1 = bool & (dd$T == 1)
    bool_t0 = bool & (dd$T == 0)
    if(sum(bool_t1) == 0 | sum(bool_t0) == 0){
      te = -Inf
    } else{
      te = mean(dd$Y[bool_t1]) - mean(dd$Y[bool_t0])
      
      
      if(max_te < te){
        max_te = te
        max_sg_bool = bool
      }  
    }
    
    treatment_effects[i] = te
  }
  
  which_max = which.max(treatment_effects)
  # jaccard for that subgroup
  jaccard_similarity = sum(gt_max_bool & max_sg_bool) / sum(gt_max_bool | max_sg_bool)
  return(list(
    jaccard_similarity = jaccard_similarity,
    treatment_effect = treatment_effects[which_max]
  )
  )
}

simulator_names = c(
  "simulate1",
  "simulate_imbalance_treatment",
  "simulate_long_rule"
)
# simulate low signal data
num_simulators = length(simulator_names)
ns = seq(1000,5000,500)
iters = 1:50
num_exp_res = length(ns) * length(iters)
gt_te = rep(0, num_exp_res)
tree_max_te = rep(0, num_exp_res)
tree_jaccard = rep(0, num_exp_res)
rsides_jaccard = rep(0, num_exp_res)
rsides_max_te = rep(0, num_exp_res)
counter = 0

full_results = data.frame(
  simulator_name = rep("", num_exp_res),
  n = rep(0, num_exp_res),
  iter = rep(0, num_exp_res),
  rsides_jaccard = rsides_jaccard,
  rsides_max_te = rsides_max_te,
  tree_jaccard = tree_jaccard,
  tree_max_te = tree_max_te,
  gt_te = gt_te,
  tree_diff_te = rep(0, num_exp_res),
  rsides_diff_te = rep(0, num_exp_res)
)


# already specified in "./util_synthetic.R"
# gt_rule = list(
#   simulator1 = "X1 > 1",
#   simulator_imbalance_treatment = "X1 > 1",
#   simulator_long_rule = "X1 > -1 & X2 > -1 & X3 > -1 & X4 > -1"
# )

for(simulator_name in simulator_names){
  # if(simulator_name != "simulate_long_rule") next()
  for(n in ns){
    for(iter_ in iters){
      counter = counter + 1
      set.seed(counter)
      
      file_path = paste0("./new_simulation/", simulator_name, "/n_", n, "_iter_", iter_, ".csv")
      d = read.csv(file_path)
      
      train_indices = sample(1:n, n * 0.5, replace = FALSE)
      d_train = d[train_indices, ]
      d_test = d[-train_indices, ]
      
      
      analysis_method = "Z-test for proportions"
      
      d_feature = d[, -c(ncol(d)-1, ncol(d))]  # exclude T and Y
      col_types = sapply(d_feature, col_type)
      
      class_cols = names(col_types[col_types == "nominal"])
      cont_cols = names(col_types[col_types == "numeric"])
      endpoint_parameters = list(outcome_variable = "Y",
                                 type = "binary",
                                 label = "Y",
                                 analysis_method = analysis_method,
                                 cont_covariates = paste(cont_cols, collapse = ", "),
                                 class_covariates = paste(class_cols, collapse = ", "),
                                 direction = 1)
      
      biomarker_names = colnames(d)[1:(ncol(d)-2)]
      biomarker_types = sapply(d[, biomarker_names], col_type)
      
      data_set_parameters = list(data_set = d_train,
                                 treatment_variable_name = "T",
                                 treatment_variable_control_value = "0",
                                 biomarker_names = biomarker_names,
                                 biomarker_types = biomarker_types)
      
      algorithm_parameters = list(
        n_perms_mult_adjust = 10,
        min_subgroup_size = 10,
        subgroup_search_algorithm = "SIDES procedure",
        ncores = 1,
        random_seed = 42,
        depth = ncol(d)-2, 
        width = 2)
      
      parameters = list(endpoint_parameters = endpoint_parameters,
                        data_set_parameters = data_set_parameters,
                        algorithm_parameters = algorithm_parameters)
      results = SubgroupSearch(parameters)
      # results
      
      # run the code
      gt_max_bool = get_gt_bool(d_test, simulator_name)
      res_list = get_subgroup_bool(d_test, results, gt_max_bool)
      gt_te[counter] =
        mean(d_test$Y[gt_max_bool & (d_test$T == 1)]) - mean(d_test$Y[gt_max_bool & (d_test$T == 0)])
  
      
      
      
      # benchmark against RPART
      
      # Fit a decision tree model
      tree_model <- rpart(
        Y ~ X1 + X2 + T,
        data = d,
        method = "class"
      )
      # rpart.plot(tree_model, type = 3, extra = 101)
      
      # get the rules
      tree_rules = rpart.rules(tree_model, cover=F)
      num_leaves = sum(tree_model$frame$var =="<leaf>")
      tree_treatment_effects = rep(0, nrow(tree_rules))
      tree_max_te = -Inf
      for(i in 1:num_leaves){
        
        rs = do.call(paste, tree_rules[i,])
        rs = strsplit(rs, "when")[[1]][2]
        if(is.na(rs)){
          bool_test = rep(T, nrow(d_test))
        } else{
          rs_vec = strsplit(rs, split="&")[[1]]
          bool_test = rep(T, nrow(d_test))
          for(j in 1:length(rs_vec)){
            r = rs_vec[j]
            
            operator <- regmatches(r, regexpr("(>=|<=|>|<|=)", r))
            if(length(operator) == 0){
              operator = "in"  
              feature = strsplit(r, "is")[[1]][1]
              feature = gsub(" ", "", feature)
              
              values = strsplit(r, "is")[[1]][2]
              values = as.numeric(strsplit(values, "to")[[1]])
              if(feature == "T"){
                next
              }
              bool_test = bool_test & (d_test[[feature]] >= values[1] & d_test[[feature]] <= values[2])
            } else{
              feature_value = strsplit(r, operator)[[1]]
              feature = feature_value[1]   
              feature = gsub(" ", "", feature)
              if(feature == "T"){
                next
              }
              value = as.numeric(feature_value[2])
              bool_test = bool_test & switch(
                operator,
                ">" = d_test[[feature]] > value,
                "<" = d_test[[feature]] < value,
                ">=" = d_test[[feature]] >= value,
                "<=" = d_test[[feature]] <= value,
                "=" = d_test[[feature]] == value
              )
            }
          }
        }
        if((sum(bool_test & (d_test$T == 1))==0) | 
           (sum(bool_test & (d_test$T == 0))==0)){
          next
        }
        te = mean(d_test$Y[bool_test & (d_test$T == 1)]) - mean(d_test$Y[bool_test & (d_test$T == 0)])
        if(tree_max_te < te){
          tree_max_te = te
          tree_max_sg_bool = bool_test
        }
      }
      which_max = which.max(tree_treatment_effects)
      # jaccard for that subgroup
      jaccard_similarity = sum(gt_max_bool & tree_max_sg_bool) / sum(gt_max_bool | tree_max_sg_bool)
      # cat("Jaccard similarity (RPART):", jaccard_similarity, "\n")
      # cat("Treatment effect (RPART):", tree_max_te, "\n")
      full_results[counter, ] = list(simulator_name, n, iter_, res_list$jaccard_similarity, 
                                  res_list$treatment_effect, jaccard_similarity, tree_max_te, gt_te[counter],
                                  tree_diff_te = abs(tree_max_te - gt_te[counter]), 
                                  rsides_diff_te = abs(res_list$treatment_effect - gt_te[counter])
                                  )
      print(full_results[counter, ,drop=F])
    }
    # round up and save the results
    full_results$rsides_jaccard = round(full_results$rsides_jaccard, 3)
    full_results$rsides_max_te = round(full_results$rsides_max_te, 3)
    full_results$tree_jaccard = round(full_results$tree_jaccard, 3)
    full_results$tree_max_te = round(full_results$tree_max_te, 3)
    full_results$gt_te = round(full_results$gt_te, 3)
    full_results$tree_diff_te = round(full_results$tree_diff_te, 3)
    full_results$rsides_diff_te = round(full_results$rsides_diff_te, 3)
    
    write.csv(full_results, "res_SIDES_and_CART.csv", row.names = FALSE)
  }
}
require(data.table)
full_res_dt = data.table(full_results)
# summarize as mean and sd with different n
summary_res = full_res_dt[, .(
  mean_rsides_jaccard = mean(rsides_jaccard),
  mean_tree_jaccard = mean(tree_jaccard),
  mean_rsides_max_te = mean(rsides_max_te),
  mean_tree_max_te = mean(tree_max_te),
  mean_gt_te = mean(gt_te),
  mean_rsides_diff_te = mean(rsides_diff_te),
  mean_tree_diff_te = mean(tree_diff_te)
), by = .(n, simulator_name)]
summary_res

write.csv(summary_res, "summary_res_SIDES_and_CART.csv", row.names = FALSE)
