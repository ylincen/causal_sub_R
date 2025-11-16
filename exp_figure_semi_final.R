rm(list=ls())

require(data.table)
require(ggplot2)
require(dplyr)
options(warn=0)
# get_info_acic2016 = function(){
#   require(aciccomp2016)
#   nsamples = rep(0, nrow(parameters_2016))
#   ncols = rep(0, nrow(parameters_2016))
#   variance_te = rep(0, nrow(parameters_2016))
#   for(i in 1:nrow(parameters_2016)){
#     cat("Processing ACIC2016 dataset ", i, " of ", nrow(parameters_2016), "\n")
#     acic2016 = dgp_2016(input_2016, parameters_2016[i,], random.seed = 1)
#     dd2 = input_2016
#     # one-hot encode dd$x2 and dd$x_24
#     library(fastDummies)
#     
#     dd2$x_2 = as.factor(dd2$x_2)
#     dd2$x_24 = as.factor(dd2$x_24)
#     dd <- dummy_cols(dd2, select_columns = c("x_2", "x_21" ,"x_24"), remove_first_dummy = F)
#     # remove the original columns as their one-hot encoded version have been included.
#     dd$x_2 = NULL
#     dd$x_24 = NULL
#     dd$x_21 = NULL
#     
#     dd$treatment = acic2016$z
#     dd$y0 = acic2016$y.0
#     dd$y1 = acic2016$y.1
#     
#     variance_te[i] = var(dd$y1 - dd$y0, na.rm = TRUE)
#     
#     dd$y_factual = ifelse(dd$treatment == 1, dd$y1, dd$y0)
#     dd$y_cfactual = ifelse(dd$treatment == 1, dd$y0, dd$y1)
#     
#     dd$y0 = NULL
#     dd$y1 = NULL
#     
#     nsamples[i] = nrow(dd)
#     ncols[i] = ncol(dd) - 4 # remove treatment, y_factual, y_cfactual, and y0/y1
#     
#   }
#   return(list(nsamples = nsamples, ncols = ncols, variance_te=variance_te))
# }
# 
# acic2016_data_info = get_info_acic2016()
# nsamples = acic2016_data_info$nsamples
# ndims = acic2016_data_info$ncols
# variance_te = acic2016_data_info$variance_te

# d_cart = fread("./semi_cart_results.csv")
d_cart = fread("./FINAL_FINAL_semi_cart_results.csv")
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

# rbind all the datasets together, and add the algorithm
d_cart$algorithm = "Ours"
d_rsides$algorithm = "SIDES"
d_quint$algorithm = "QUINT"
d_ct$algorithm = "CausalTree"
d_it$algorithm = "InteractionTree"
d_distill$algorithm = "DistillTree"
d_curls$algorithm = "CURLS"
d_combine = rbind(d_cart, d_rsides, d_quint, d_ct, d_it, d_distill, d_curls)
d_combine$max_subgroup_diff_average = abs(d_combine$max_subgroup_diff_average)

# match the dataset with three simulator "ACIC2016", "IBM_scaling", "IHDP", and
#. create a new column called semi-synethtic-simulator
d_combine$simulator = sapply(d_combine$dataset, function(a){
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

d_combine = d_combine[d_combine$simulator == "ACIC2016", c("dataset", "algorithm",  "max_subgroup_treatment_effect", "max_subgroup_gt_treatment_effct")]
d_combine$abs_diff = abs(d_combine$max_subgroup_treatment_effect - d_combine$max_subgroup_gt_treatment_effct)
d_combine$diff = (d_combine$max_subgroup_treatment_effect - d_combine$max_subgroup_gt_treatment_effct)
d_wide = dcast(d_combine, 
        dataset ~ algorithm, 
        value.var = c("max_subgroup_gt_treatment_effct"))


ggplot(d_combine, aes(x=algorithm, y=max_subgroup_gt_treatment_effct)) + 
  geom_boxplot()

diff_plot = ggplot(d_combine, aes(x=algorithm, y=diff)) + 
  geom_boxplot(outlier.alpha = 0.5) + 
  theme_bw() + 
  xlab("") +
  scale_y_continuous(trans = "asinh") +
  ylab("Difference between estimated and \nground truth subgroup treatment effects\n (asinh scale for better visualization)")
  
diff_plot
ggsave(diff_plot, filename = "./semi_final_diff_boxplot.png", 
       width = 6, height = 3, units = "in", dpi = 300, device = "png")


tmp1 = apply(d_wide[, 2:ncol(d_wide)], 2, mean, na.rm=T)
tmp2 = apply(-d_wide[, 2:ncol(d_wide)], 1, rank, na.last="keep") %>% apply(1, function(a){
  mean(a<=1, na.rm = T)
})
final_table = data.frame(cbind(tmp1, tmp2))
colnames(final_table) = c("Mean Groud Truth Subgroup Treatment Effect",
                            "% of rank 1st among all 77 datasets")
require(xtable)
print(xtable(final_table, digits = 3), 
      include.rownames = TRUE, 
      include.colnames = TRUE, 
      floating = FALSE, 
      table.placement = "H", 
      sanitize.text.function = identity)

# apply(d_wide[, 2:ncol(d_wide)], 1, rank, na.last="keep") %>% apply(1, sd, na.rm=T)
# 
# d_wide_acc = dcast(d_combine, 
#         dataset ~ algorithm, 
#         value.var = c("diff"))
# d_wide_acc$QUINT[is.infinite(d_wide_acc$QUINT)] = NA
# apply(d_wide_acc[, 2:ncol(d_wide_acc)], 2, mean, na.rm=T)
# apply(d_wide_acc[, 2:ncol(d_wide_acc)], 2, sd, na.rm=T)
# 
# hist(d_wide_acc$Ours)
# hist(d_wide_acc$CausalTree, add=T, col="red")
