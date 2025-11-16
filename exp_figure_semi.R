rm(list=ls())

require(data.table)

# d_cart = fread("./semi_cart_results.csv")
# d_rsides = fread("./semi_rsides_results.csv")
# d_quint = fread("./semi_quint_results.csv")
# 
# d_ct = fread("./competitors/causal_tree/semi_CausalTree_results.csv")
# d_it = fread("./competitors/interaction_tree/R code for IT_NOV2024/semi_IT_results.csv")
# d_distill = fread("./competitors/distill_tree/semi_distill_results.csv")
rm(list=ls())

require(data.table)
require(ggplot2)
require(dplyr)
d_cart = fread("./semi_cart_results.csv")
d_rsides = fread("./FINAL_semi_rsides_results.csv")
d_quint = fread("./semi_quint_results.csv")

d_ct = fread("./competitors/causal_tree/NEW_semi_CausalTree_results.csv")
d_it = fread("./competitors/interaction_tree/R code for IT_NOV2024/NEW_semi_IT_results.csv")
d_distill = fread("./competitors/distill_tree/semi_distill_results.csv")
csv_files = list.files("./competitors/curls/res_server/semi_synthetic/combined/", full.names = TRUE, pattern = "*.csv")
data_info = basename(csv_files)

d_curls = do.call(rbind, lapply(csv_files, read.csv))
dataset_name = sapply(data_info, function(a){
  strsplit(a, split = "_curls")[[1]][1]
})
dataset_name = unname(dataset_name) %>% paste0(".csv")
d_curls$dataset = dataset_name
d_curls = as.data.table(d_curls)

d_curls2 = d_curls[,c("dataset", "treatment_effect_estimated",
                      "diff_gt_estimated", "treatment_effect_gt"
)]
colnames(d_curls2)[3] = "max_subgroup_diff_average"
colnames(d_curls2)[2] = "max_subgroup_treatment_effect"
colnames(d_curls2)[4] = "max_subgroup_gt_treatment_effct"

d_curls2$subgroup_diff_average_weighted_mean = NA
d_curls2$subgroup_diff_average_unweighted_mean = NA
d_curls = d_curls2

require(ggplot2)
require(dplyr)
# # combine d_cart and d_rsides and add a column named "algorithm"
# d_cart[, algorithm := "CART"]
# d_rsides[, algorithm := "RSIDES"]
# d_cart_rsides = rbind(d_cart, d_rsides)
# d_cart_rsides$semi_simulator = sapply(d_cart_rsides$dataset, strsplit, split = "_[0-9]") %>%
#   sapply(function(x) x[[1]])
# 
# # make the scatter plot with x-axis max_subgroup_treatment_effect, 
# #.  y-axis as max_subgroup_gt_treatment_effct
# # together with a reference line y = x
# ggplot(d_cart_rsides, aes(x = max_subgroup_treatment_effect, 
#                           y = max_subgroup_gt_treatment_effct, 
#                           color = algorithm, 
#                           shape = semi_simulator)) + 
#   geom_point(size = 3) + 
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black")


make_figure = function(d_cart, d_other, algorithm_name, simulator_name = NULL){
  if(!is.null(simulator_name)){
    # filter the datasets by simulator_name
    d_cart$"simulator" = sapply(d_cart$dataset, function(a){
      if(grepl("ACIC2016", a)) {
        return("ACIC2016")
      } else if(grepl("IBM_scaling", a)) {
        return("IBM_scaling")
      } else if(grepl("IHDP", a)) {
        return("IHDP")
      } else {
        return(a)
      }
    })
    d_other$"simulator" = sapply(d_other$dataset, function(a){
      if(grepl("ACIC2016", a)) {
        return("ACIC2016")
      } else if(grepl("IBM_scaling", a)) {
        return("IBM_scaling")
      } else if(grepl("IHDP", a)) {
        return("IHDP")
      } else {
        return(a)
      }
    })
    d_cart = d_cart[d_cart$simulator == simulator_name, ]
    d_other = d_other[d_other$simulator == simulator_name, ]
  }
  d_cart[, algorithm := "CART"]
  d_other[, algorithm := algorithm_name]
  
  # for quint, some results are -Inf, report how many and remove from it
  if (algorithm_name == "QUINT") {
    n_inf = sum(d_other$max_subgroup_treatment_effect == -Inf)
    cat("Number of -Inf in QUINT: ", n_inf, "\n")
    d_other = d_other[d_other$max_subgroup_treatment_effect != -Inf, ]
  }
  
  d_combined = rbind(d_cart, d_other)
  d_combined$semi_simulator = sapply(d_combined$dataset, strsplit, split = "_[0-9]") %>%
    sapply(function(x) x[[1]])
  
  ggplot(d_combined, aes(x = max_subgroup_treatment_effect, 
                         y = max_subgroup_gt_treatment_effct, 
                         color = algorithm, 
                         shape = semi_simulator)) + 
    geom_point(size = 3, alpha=0.5) + 
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    labs(title = paste("CART vs", algorithm_name),
         x = "Max Subgroup Treatment Effect (Estimated)",
         y = "Ground Truth Treatment Effect \nfor the Learned Subgroup") +
    theme_minimal()
}

make_figure(d_cart, d_distill, "Distill-Tree")
make_figure(d_cart, d_rsides, "SIDES")
make_figure(d_cart, d_quint, "QUINT")
make_figure(d_cart, d_ct, "Causal Tree")
make_figure(d_cart, d_it, "Interaction Tree")
make_figure(d_cart, d_curls, "CURLS")

make_figure(d_cart, d_distill, "Distill-Tree", "ACIC2016")

make_figure(d_cart, d_rsides, "SIDES", "ACIC2016")
make_figure(d_cart, d_quint, "QUINT", "ACIC2016")
make_figure(d_cart, d_ct, "Causal Tree", "ACIC2016")
make_figure(d_cart, d_it, "Interaction Tree", "ACIC2016")
make_figure(d_cart, d_curls, "CURLS", "ACIC2016")

# 
# csv_files = list.files("./competitors/curls/res_server/semi_synthetic/combined/", full.names = TRUE, pattern = "*.csv")
# data_info = basename(csv_files)
# 
# d_curls = do.call(rbind, lapply(csv_files, read.csv))
# dataset_name = sapply(data_info, function(a){
#   strsplit(a, split = "_curls")[[1]][1]
# })
# dataset_name = unname(dataset_name) %>% paste0(".csv")
# d_curls$dataset = dataset_name
# d_curls = as.data.table(d_curls)
# 
# 
# # select the columns from d_curls only: dataset, treatment_effect_estimated, 
# #.  treatment_effect_gt
# d_curls = d_curls[, .(dataset, max_subgroup_treatment_effect = treatment_effect_estimated, 
#                        max_subgroup_gt_treatment_effct = treatment_effect_gt)]
# # select d_cart_sub as well
# d_cart_sub = d_cart[, .(dataset, max_subgroup_treatment_effect, 
#                          max_subgroup_gt_treatment_effct)]
# make_figure(d_cart_sub, d_curls, "CURLS")
# 
# 
# 
# 
# 
# 
# 
# 
